# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/20
content:    Compare our results with other datasets.
'''
import os
import sys
import glob
import gzip
import subprocess as sp
import numpy as np
import pandas as pd
import loompy

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet import Dataset, concatenate


fdn_data = '../../data/sequencing/datasets/immune_final/'
fns_loom = {
        'ours': fdn_data+'immune.loom',
        'tabulamuris': '../../data/tabulamuris/FACS_alltissues/immune.loom',
        }


# Add northstar
sys.path.insert(0, '/home/fabio/university/postdoc/northstar/build/lib')
import northstar



def create_loom_from_figshare(fn_loom):
    '''Create a loom file from our FigShare files'''
    print('Convert TSV into a loom file for convenience')
    fn_counts = fdn_data+'table_counts_lungimmune.tsv.gz'
    fn_meta = fdn_data+'table_cellmetadata_lungimmune.tsv.gz'

    print('Load metadata')
    samplesheet = pd.read_csv(
            fn_meta, sep='\t', index_col=0,
            compression='gzip',
            )

    print('Load counts')
    counts = pd.read_csv(
            fn_counts, sep='\t', index_col=0,
            compression='gzip',
            ).astype(np.float32)

    print('Normalize by coverage')
    counts = (1e6 * counts / samplesheet['Coverage']).astype(np.float32)

    col_attrs = {col: samplesheet[col].values for col in samplesheet.columns}
    col_attrs['CellID'] = samplesheet.index.values
    row_attrs = {'GeneName': counts.index.values}
    loompy.create(
        fn_loom,
        layers={'': counts.values},
        col_attrs=col_attrs,
        row_attrs=row_attrs,
        )


def create_loom_from_tabulamurisfacs(fn_loom):
    samplesheet = pd.read_csv(
        '../../data/tabulamuris/FACS_alltissues/annotations_facs.csv',
        sep=',',
        index_col='cell',
        low_memory=False,
        ).iloc[:, 2:]
    samplesheet.index.name = 'CellID'

    immune_types = [
        'B cell',
        'DN1 thymic pro-T cell',
        'regulatory T cell',
        'basophil',
        'pre-natural killer cell',
        'immature T cell',
        'myeloid cell',
        'T cell',
        'granulocyte',
        'naive B cell',
        'leukocyte',
        'precursor B cell',
        'macrophage',
        'immature B cell',
        'monocyte',
        'late pro-B cell',
        'natural killer cell',
        'granulocyte monocyte progenitor cell',
        'classical monocyte', 'lymphocyte',
        'professional antigen presenting cell',
        'mature natural killer cell',
        'immature NK T cell',
        'immature natural killer cell',
    ]

    samplesheet = samplesheet.loc[samplesheet['cell_ontology_class'].isin(
        immune_types)]

    print('Tabula Muris has a total of {:} immune cells of {:} types'.format(
        samplesheet.shape[0], len(immune_types)))

    cnames_unsort = samplesheet.index

    cellnames = []
    genes = []
    counts = []
    fns = glob.glob('../../data/tabulamuris/FACS_alltissues/FACS/FACS/*.loom')
    for ifn, fn in enumerate(fns):
        tissue = os.path.basename(fn)[:-len('-counts.loom')]
        print('Mining {:} ({:}/{:})'.format(tissue, ifn+1, len(fns)))
        with loompy.connect(fn) as dsl:
            cnsus = dsl.ca['CellID']
            idx = pd.Index(cnsus).isin(cnames_unsort).nonzero()[0]
            cns = cnsus[idx]
            cos = dsl[:, idx]

            cellnames.append(cns)
            counts.append(cos)
            genes.append(dsl.ra['GeneName'])

    # Check that they all have the same genes
    if len(set([tuple(x) for x in genes])) > 1:
        print('WARNING: not all tissues have the same genes')
        return {'ss': samplesheet, 'counts': counts, 'cns': cellnames, 'genes': genes}

    print('Merging into single loom file')
    cellnames = np.concatenate(cellnames)
    counts = np.hstack(counts)
    genes = genes[0]
    samplesheet = samplesheet.loc[cellnames]

    print('Writing loom file')
    col_attrs = {col: samplesheet[col].values for col in samplesheet.columns}
    col_attrs['CellID'] = samplesheet.index.values
    row_attrs = {'GeneName': genes}
    loompy.create(
        fn_loom,
        layers={'': counts},
        col_attrs=col_attrs,
        row_attrs=row_attrs,
        )


if __name__ == '__main__':

    if not os.path.isfile(fns_loom['ours']):
        create_loom_from_figshare(fns_loom['ours'])

    if not os.path.isfile(fns_loom['tabulamuris']):
        ss = create_loom_from_tabulamurisfacs(fns_loom['tabulamuris'])

    print('Load loom file')
    ds = Dataset(
            dataset={
                'path': fns_loom['ours'],
                'index_samples': 'CellID',
                'index_features': 'GeneName',
                },
        )
    ds.samplesheet['Tissue'] = 'lung'

    print('Load Tabula Muris FACS immune')
    ds_tm = Dataset(
            dataset={
                'path': fns_loom['tabulamuris'],
                'index_samples': 'CellID',
                'index_features': 'GeneName',
                },
        )
    ds_tm.samplesheet['Tissue'] = ds_tm.samplesheet['tissue']
    ds_tm.samplesheet['Cell Subtype'] = ds_tm.samplesheet['cell_ontology_class']
    ds_tm.samplesheet['Cell Subtype'].replace(
        {
            'natural killer cell': 'NK cell',
            'mature natural killer cell': 'NK cell',
            'pre-natural killer cell': 'NK cell',
            'immature natural killer cell': 'NK cell',
            'immature NK T cell': 'NKT cell',
            'immature T cell': 'T cell',
            'regulatory T cell': 'T cell',
            'DN1 thymic pro-T cell': 'T cell',
            'naive B cell': 'B cell',
            'immature B cell': 'B cell',
            'precursor B cell': 'B cell',
            'late pro-B cell': 'B cell',
            'classical monocyte': 'monocyte',
        },
        inplace=True)
    # Some cell types from Tabula Muris are very unclear
    skip = [
        'professional antigen presenting cell',
        'lymphocyte',
        'leukocyte',
        'myeloid cell',
        ]
    ds_tm.query_samples_by_metadata(
        'cell_ontology_class not in @skip',
        local_dict=locals(),
        inplace=True)

    print('Merge etc based on northstar')
    ns = northstar.Subsample(
            atlas={
                'cell_types': ds_tm.samplesheet['Cell Subtype'],
                'counts': ds_tm.counts,
                },
            join='intersection',
            n_pcs=35,
            resolution_parameter=0.001,
        )

    ns.new_data = ds.counts
    ns._check_init_arguments()
    ns.fetch_atlas_if_needed()
    ns.compute_feature_intersection()
    ns._check_feature_intersection()
    ns.prepare_feature_selection()
    ns.select_features()
    ns._check_feature_selection()
    ns.merge_atlas_newdata()

    print('Make PCA and graph')
    ns.compute_pca()
    ns.compute_similarity_graph()

    print('Cluster graph')
    ns.cluster_graph()

    print('Compute embedding')
    vs = ns.embed('umap')

    print('Make dataset with merged')
    genes = np.intersect1d(ds.featurenames, ds_tm.featurenames)
    ds.query_features_by_name(genes, inplace=True)
    ds.samplesheet['Data source'] = 'new_data'
    ds_tm.query_features_by_name(genes, inplace=True)
    ds_tm.samplesheet['Data source'] = 'Tabula Muris'
    dsme = concatenate([ds_tm, ds])
    dsme.samplesheet['northstar assignment'] = np.concatenate(
            [ds_tm.samplesheet['Cell Subtype'].values, ns.membership],
            )
    new_clusters = [x for x in np.unique(ns.membership) if x.isdigit()]

    print('Plot embedding')
    genes = ['Data source', 'Cell Subtype', 'northstar assignment', 'Tissue']
    cmaps = {
        'Data source': {'Tabula Muris': 'darkred', 'new_data': 'steelblue'},
        'Cell Subtype': {
            'Mac IV': 'darkolivegreen',
            'Mac II': 'lime',
            'Mac I': 'seagreen',
            'Mac III': 'greenyellow',
            'Mac V': 'turquoise',
            'DC I': 'darkviolet',
            'DC II': 'fuchsia',
            'DC III': 'violet',
            'B cell': 'lightcoral',
            'T cell': 'tomato',
            'NKT cell': 'saddlebrown',
            'NK cell': 'chocolate',
            'IL cell': 'tan',
            'mast cell': 'gold',
            'basophil': 'grey',
            'neutrophil': 'black',
            'monocyte': 'steelblue',
            'granulocyte': 'forestgreen',
            'granulocyte monocyte progenitor cell': 'seagreen',
            'macrophage': 'indigo',
        },
        'northstar assignment': {
            'B cell': 'lightcoral',
            'T cell': 'tomato',
            'NKT cell': 'saddlebrown',
            'NK cell': 'chocolate',
            'IL cell': 'tan',
            'basophil': 'grey',
            'neutrophil': 'black',
            'monocyte': 'steelblue',
            'granulocyte': 'forestgreen',
            'granulocyte monocyte progenitor cell': 'seagreen',
            'macrophage': 'indigo',
        },
    }
    cluster_additional_colors = [
            'darkolivegreen',
            'lime',
            'greenyellow',
            'turquoise',
            'darkviolet',
            'fuchsia',
            'violet',
            ]
    for i, newclu in enumerate(new_clusters):
        cmaps['northstar assignment'][newclu] = cluster_additional_colors[i]

    fig, axs = plt.subplots(1, 4, figsize=(19, 6.8), sharex=True, sharey=True)
    for i in range(len(axs)):
        gene = genes[i]
        ax = axs[i]
        cmap = cmaps.get(gene, 'viridis')
        dsme.plot.scatter_reduced(
                vs,
                color_by=gene,
                color_log=False,
                cmap=cmap,
                ax=ax,
                alpha=0.2 - (0.1 * (gene == 'tissue')),
                s=15,
                )
        #ax.set_axis_off()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(gene)
        handles, labels = [], []
        for key, color in ax._singlet_cmap.items():
            handles.append(ax.scatter([], [], color=color))
            labels.append(key)
        ax.legend(
                handles, labels, loc='upper center',
                fontsize=8, ncol=2,
                bbox_to_anchor=(0.5, -0.05), bbox_transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()

    plt.ion()
    plt.show()
