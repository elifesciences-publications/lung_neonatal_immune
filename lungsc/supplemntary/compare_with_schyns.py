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
        'schyns': '../../data/Schyns/loom/schyns.loom',
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
    ds.samplesheet['cellSubtype'] = ds.samplesheet['Cell Subtype']
    ds.query_samples_by_metadata(
        'cellSubtype in ("Mac I", "Mac II", "Mac III", "Mac IV", "Mac V")',
        inplace=True)

    print('Load Schyns et al loom file')
    ds_sc = Dataset(
            dataset={
                'path': fns_loom['schyns'],
                'index_samples': 'CellID',
                'index_features': 'Gene',
                'bit_precision': 32,
                },
        )
    ds_sc.samplesheet['cellSubtype'] = ds_sc.samplesheet['ClusterName']

    print('Merge the data by hand')
    genes = np.intersect1d(ds.featurenames, ds_sc.featurenames)
    ds.query_features_by_name(genes, inplace=True)
    ds.samplesheet['Data source'] = 'Domingo-Gonzales\net al. (this paper)'
    ds_sc.query_features_by_name(genes, inplace=True)
    ds_sc.samplesheet['Data source'] = 'Schyns et al.'
    dsme = concatenate([ds, ds_sc])
    dsme.counts.normalize('counts_per_million', inplace=True)

    print('Feature selection')
    features = dsme.feature_selection.overdispersed_within_groups(
            'Data source', inplace=False)
    dsf = dsme.query_features_by_name(features)

    print('PCA')
    dsc = dsf.dimensionality.pca(n_dims=30, robust=False, return_dataset='samples')

    print('Embedding')
    vs = dsc.dimensionality.umap()
    dsme.samplesheet['umap1'] = vs.iloc[:, 0]
    dsme.samplesheet['umap2'] = vs.iloc[:, 1]

    print('Plot embedding')
    fig, axs = plt.subplots(2, 5, figsize=(10, 4.5), sharex=True, sharey=True)
    axs = list(axs.ravel())
    fig2, axs2 = plt.subplots(1, 2, figsize=(8, 5.5), sharex=True, sharey=True)
    axs.extend(list(axs2))
    marker_genes = [
            #('Pecam1', 'Cd31'),
            #('Ptprc', 'Cd45'),
            #'Col6a2',
            #'Cd3e',
            #'Gzma',
            #'Arg1',
            #('Ms4a1', 'Cd20'),
            'Cd68',
            'Gal',
            'Itgax',
            'Car4',
            'C1qa',
            'Plac8',
            #'H2',
            'Cd72',
            'H2-Eb1',
            'Folr2',
            #'Batf3',
            #'Itgae',
            #'Mreg',
            #'Cd209c',
            #'Stfa2',
            #'Cpa3',
            #'Mcpt4',
            'Mki67',
            ]
    markers = [
            #'SortType',
            ] + marker_genes
    markers += [
            #('number_of_genes_1plusreads', 'n genes'),
            #'Timepoint',
            #'Gender',
            #'Mousename',
            #'cellType',
            #'Treatment',
            'cellSubtype',
            'Data source',
            ]
    mgs = [x if isinstance(x, str) else x[0] for x in marker_genes]
    for ipl, (gene, ax) in enumerate(zip(markers, axs)):
        print('Plotting gene {:} of {:}'.format(ipl+1, len(markers)))
        if isinstance(gene, str):
            gene, title = gene, gene
        else:
            gene, title = gene
        if gene == 'cellSubtype':
            cmap = {
                'Mac IV': 'darkolivegreen',
                'Mac II': 'lime',
                'Mac I': 'seagreen',
                'Mac III': 'greenyellow',
                'Mac V': 'turquoise',
                'AM': 'darkgrey',
                'CD16.2 MI': 'navy',
                'CD206 MI': 'darkorange',
                'MHCII MI': 'dodgerblue',
            }
        elif gene == 'cellType':
            cmap = sns.color_palette('Set1', n_colors=6)
        elif gene == 'Timepoint':
            cmap = {
                'E18.5': 'navy',
                'P1': 'gold',
                'P7': 'tomato',
                'P21': 'firebrick',
            }
        else:
            cmap = 'viridis'
        dsme.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=13 + 5 * (ax in axs[-2:]),
                alpha=0.25,
                color_by=gene,
                color_log=(gene in mgs + ['number_of_genes_1plusreads']),
                cmap=cmap,
                )
        ax.grid(False)
        ax.set_title(title)

        if gene in ('Treatment', 'Timepoint', 'cellSubtype', 'Data source'):
            import matplotlib.lines as mlines
            d = ax._singlet_cmap
            handles = []
            labels = []
            if gene == 'cellSubtype':
                keys = ['Mac I', 'Mac II', 'Mac III', 'Mac IV', 'Mac V',
                        'AM', 'CD16.2 MI', 'CD206 MI', 'MHCII MI']
            else:
                keys = d.keys()
            for key in keys:
                color = d[key]
                h = mlines.Line2D(
                    [], [], color=color, marker='o', lw=0,
                    markersize=5,
                    )
                handles.append(h)
                if gene == 'Treatment':
                    key = key[0]
                    loc = 'upper left'
                else:
                    loc = 'lower left'
                labels.append(key.upper())
            if gene == 'Timepoint':
                labels_old = list(labels)
                labels = ['E18.5', 'P1', 'P7', 'P21']
                handles = [handles[labels_old.index(li)] for li in labels]
            ax.legend(
                handles, labels, loc='upper center',
                fontsize=9,
                ncol=2,
                bbox_to_anchor=(0.5, -0.15),
                bbox_transform=ax.transAxes,
                )

    fig.tight_layout()
    fig2.tight_layout()
    fig_fdn = '../../figures/endomese_share/endo_paper_supplementary_figs/'
    fig.savefig(fig_fdn+'comparison_with_Schyns_genes_lowres.png', dpi=300)
    fig2.savefig(fig_fdn+'comparison_with_Schyns_metadata_lowres.png', dpi=300)

    plt.ion()
    plt.show()
