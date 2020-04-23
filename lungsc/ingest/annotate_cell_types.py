# vim: fdm=indent
'''
author:     Fabio Zanini
date:       04/12/19
content:    Annotate cell types
'''
import os
import sys
import glob
import gzip
import argparse
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

# Ensure leidenalg is correct
sys.path.insert(0, os.path.abspath('../../packages'))
from lungsc.ingest.load_dataset import versions, DatasetLung


# Manual annotation functions
def figure_ct_cst(ds, vs):
    x0, x1 = plt.gcf().get_axes()[0].get_xlim()
    y0, y1 = plt.gcf().get_axes()[0].get_ylim()
    print(x0, x1, y0, y1)
    ind = (vs['dimension 1'] >= x0) & (vs['dimension 1'] < x1)
    ind &= (vs['dimension 2'] >= y0) & (vs['dimension 2'] < y1)
    print(ds.samplesheet.loc[ind, 'Treatment'].value_counts())
    print(ds.samplesheet.loc[ind, 'cellType'].value_counts())
    print(ds.samplesheet.loc[ind, 'cellSubtype'].value_counts())


def assign_tsne_fig_ct(ds, vs, ct, cst, ds0, overwrite=False):
    x0, x1 = plt.gcf().get_axes()[0].get_xlim()
    y0, y1 = plt.gcf().get_axes()[0].get_ylim()
    print(x0, x1, y0, y1)
    ind = (vs['dimension 1'] >= x0) & (vs['dimension 1'] < x1)
    ind &= (vs['dimension 2'] >= y0) & (vs['dimension 2'] < y1)
    ind = vs.index[ind]
    ds0.samplesheet.loc[ind, 'cellType'] = ct
    ds0.samplesheet.loc[ind, 'cellSubtype'] = cst
    ds.samplesheet.loc[ind, 'cellType'] = ct
    ds.samplesheet.loc[ind, 'cellSubtype'] = cst
    for fign in plt.get_fignums():
        fig = plt.figure(fign)
        for ax in fig.get_axes():
            gene = ax.get_title()
            if gene == 'cellSubtype':
                ax.clear()
                ind = np.arange(40)
                np.random.shuffle(ind)
                cmap = sns.color_palette('husl', n_colors=40)
                cmap = [cmap[i] for i in ind]
                csts = ds.samplesheet['cellSubtype'].unique()
                cmap = dict(zip(csts, cmap))
                cmap[''] = (0, 0, 0)
            elif gene == 'cellType':
                cmap = sns.color_palette('Set1', n_colors=6)
            else:
                continue

            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=10,
                    alpha=0.3,
                    color_by=gene,
                    color_log=False,
                    cmap=cmap,
                    )
            ax.grid(False)
            ax.set_title('cellSubtype')

            break


def tag_doublets(ds, vs, ds0):
    x0, x1 = plt.gcf().get_axes()[0].get_xlim()
    y0, y1 = plt.gcf().get_axes()[0].get_ylim()
    print(x0, x1, y0, y1)
    ind = (vs['dimension 1'] >= x0) & (vs['dimension 1'] < x1)
    ind &= (vs['dimension 2'] >= y0) & (vs['dimension 2'] < y1)
    ind = vs.index[ind]
    ds0.samplesheet.loc[ind, 'doublet'] = 1
    ds.samplesheet.loc[ind, 'doublet'] = 1
    ds0.samplesheet.loc[ind, 'cellType'] = ''
    ds.samplesheet.loc[ind, 'cellType'] = ''
    ds0.samplesheet.loc[ind, 'cellSubtype'] = ''
    ds.samplesheet.loc[ind, 'cellSubtype'] = ''



if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument(
            '--cellType',
            choices=['immune', 'general'],
            required=True,
            )
    pa.add_argument(
            '--includeNoCellType',
            action='store_true',
            )

    args = pa.parse_args()

    version = versions[-1]
    ds = DatasetLung.load(version=version, include_doublets=True)
    ds.samplesheet.loc[ds.samplesheet['doublet'] == 1, 'cellSubtype'] = 'doublet?'

    if args.cellType == 'general':
        print('Plot dimensionality reduction of dataset')
        print('Feature selection')
        features = ds.feature_selection.overdispersed_within_groups(
                'Mousename',
                inplace=False,
                )
        dsf = ds.query_features_by_name(features)

        print('PCA')
        dsc = dsf.dimensionality.pca(
                n_dims=30,
                robust=False,
                return_dataset='samples',
                )

        print('tSNE')
        vs = dsc.dimensionality.tsne(perplexity=20)

        print('knn graph')
        neighbors, _, _ = dsc.graph.knn(
                axis='samples',
                n_neighbors=10,
                return_sparse=False,
                )

        print('Leiden clustering')
        edges = set()
        for nei in neighbors:
            for n in nei:
                edges.add(frozenset(n))
        edges = [tuple(e) for e in edges]

        ds.samplesheet['cluster'] = ds.cluster.leiden(
                axis='samples',
                edges=edges,
                resolution_parameter=0.0005,
                )

        print('Initial guess of cell type for each cluster')
        di = {'CD45': 'immune', 'CD31': 'endothelial', 'mesench': 'mesenchymal'}
        ds.samplesheet['cellType'] = ''
        for cl in ds.samplesheet['cluster'].unique():
            ind = ds.samplesheet['cluster'] == cl
            st = ds.samplesheet.loc[ind, 'SortType']
            ds.samplesheet.loc[ind, 'cellType'] = di[st.value_counts().idxmax()]

        def assign_ct(cluster_numbers, cell_type):
            if np.isscalar(cluster_numbers):
                cluster_numbers = [cluster_numbers]
            ind = ds.samplesheet['cluster'].isin(cluster_numbers)
            ds.samplesheet.loc[ind, 'cellType'] = cell_type

        fig, axs = plt.subplots(3, 5, figsize=(13, 7), sharex=True, sharey=True)
        axs = axs.ravel()
        marker_genes = [
                #('Pecam1', 'Cd31'),
                ('Ptprc', 'Cd45'),
                #'Col6a2',
                'Cd3e',
                'Gzma',
                ('Ms4a1', 'Cd20'),
                'Cpa3',
                'Mcpt4',
                'Plac8',
                'Cd68',
                'Batf3',
                'Cd209c',
                'Stfa2',
                #'Gja5',
                #'Maf',
                #'Car8',
                #'Car4',
                #'Tgfbi',
                #'Mcam',
                #'Pdgfra',
                #'Hhip',
                #'Wnt2',
                #'Trpc6',
                'Mki67',
                ]
        markers = [
                'SortType',
                ] + marker_genes[:-1]
        if 'Treatment' in ds.samplesheet.columns:
            markers.append('Treatment')
        else:
            markers.append(marker_genes[-1])
        markers += [
                #('number_of_genes_1plusreads', 'n genes'),
                #'Timepoint',
                #'Gender',
                #'Mousename',
                'cellType',
                'cellSubtype',
                ]
        mgs = [x if isinstance(x, str) else x[0] for x in marker_genes]
        for ipl, (gene, ax) in enumerate(zip(markers, axs)):
            print('Plotting gene {:} of {:}'.format(ipl+1, len(markers)))
            if isinstance(gene, str):
                gene, title = gene, gene
            else:
                gene, title = gene
            if gene == 'cellSubtype':
                ncol = len(ds.samplesheet['cellSubtype'].unique())
                ind = np.arange(ncol)
                np.random.shuffle(ind)
                cmap = sns.color_palette('Set1', n_colors=ncol)
                cmap = [cmap[i] for i in ind]
            elif gene == 'cellType':
                cmap = sns.color_palette('Set1', n_colors=6)
            else:
                cmap = 'viridis'
            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=10,
                    alpha=0.04 + 0.1 * (gene not in ['annotated', 'cluster', 'cellType', 'Mousename', 'Treatment', 'Timepoint']),
                    color_by=gene,
                    color_log=(gene in mgs + ['number_of_genes_1plusreads']),
                    cmap=cmap,
                    )
            ax.grid(False)
            ax.set_title(title)

            if gene == 'cluster':
                for com in ds.samplesheet['cluster'].unique():
                    vsc = vs.loc[ds.samplesheet[gene] == com]
                    xc, yc = vsc.values.mean(axis=0)
                    ax.scatter([xc], [yc], s=10, facecolor='none', edgecolor='red', marker='^')
                    ax.text(xc, yc, str(com), fontsize=8, ha='center', va='bottom')

            if gene in ('Treatment', 'Timepoint'):
                import matplotlib.lines as mlines
                d = ax._singlet_cmap
                handles = []
                labels = []
                for key, color in d.items():
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
                    ncol = 2
                else:
                    ncol = 1
                ax.legend(handles, labels, loc=loc, fontsize=6, ncol=ncol)

        fig.tight_layout()

        print('Plot big panel with all cell subtypes')
        import random
        import matplotlib.lines as mlines
        fig, ax = plt.subplots(figsize=(13.5, 9.8))
        gene = 'cellSubtype'
        cts = np.sort(ds.samplesheet['cellSubtype'].unique())
        colors = sns.color_palette('husl', n_colors=len(cts))
        random.shuffle(colors)
        cmap = dict(zip(cts, colors))
        markerd = {'immune': 'o', 'endothelial': 's', 'mesenchymal': '*'}
        handles = []
        labels = []
        for ct in ['immune', 'endothelial', 'mesenchymal']:
            dsi = ds.query_samples_by_metadata('(cellType == @ct) & (cellSubtype not in ("", "nan"))', local_dict=locals())
            vsi = vs.loc[dsi.samplenames]
            dsi.plot.scatter_reduced_samples(
                    vsi,
                    ax=ax,
                    s=25,
                    alpha=0.4,
                    marker=markerd[ct],
                    color_by=gene,
                    color_log=False,
                    cmap=cmap,
                    )
            for key, color in ax._singlet_cmap.items():
                if (dsi.samplesheet['cellSubtype'] == key).sum() == 0:
                    continue
                h = mlines.Line2D(
                    [], [], color=color, marker=markerd[ct], lw=0,
                    markersize=5,
                    )
                handles.append(h)
                labels.append(key.upper())
        ax.legend(
            handles, labels,
            loc='upper left',
            fontsize=10,
            ncol=1,
            title='Cell subtype:',
            bbox_to_anchor=(1.01, 1.0),
            bbox_transform=ax.transAxes,
            )
        ax.grid(False)
        ax.set_axis_off()
        fig.tight_layout()

    if args.cellType == 'immune':
        print('Focus on subtypes within immune cells')
        if args.includeNoCellType:
            dsi = ds.query_samples_by_metadata('cellType in ("immune", "")')
        else:
            dsi = ds.query_samples_by_metadata('cellType == "immune"')

        if False:
            print('Feature selection')
            features = dsi.feature_selection.overdispersed_within_groups('Mousename', inplace=False)
            dsf = dsi.query_features_by_name(features)

            print('PCA')
            dsc = dsf.dimensionality.pca(n_dims=30, robust=False, return_dataset='samples')

            print('tSNE')
            vsi = dsc.dimensionality.tsne(perplexity=20)
            vsi.to_csv('../../data/sequencing/datasets/all_{:}/tsne_immune.tsv'.format(version), sep='\t', index=True)
        else:
            vsi = pd.read_csv('../../data/sequencing/datasets/all_{:}/tsne_immune.tsv'.format(version), sep='\t', index_col=0)

        print('Plot dimensionality reduction of dataset')
        vs = vsi
        #fig, axs = plt.subplots(3, 11, figsize=(17.3, 5.5), sharex=True, sharey=True)
        fig, axs = plt.subplots(3, 7, figsize=(16, 7), sharex=True, sharey=True)
        axs = axs.ravel()
        marker_genes = [
                #('Pecam1', 'Cd31'),
                #('Ptprc', 'Cd45'),
                #'Col6a2',
                'Cd3e',
                'Gzma',
                'Arg1',
                ('Ms4a1', 'Cd20'),
                'Cd68',
                'Gal',
                'Itgax',
                'Car4',
                'C1qa',
                'Plac8',
                'Batf3',
                'Itgae',
                'Mreg',
                'Cd209c',
                'Stfa2',
                'Cpa3',
                'Mcpt4',
                'Mki67',
                ]
        markers = [
                #'SortType',
                ] + marker_genes
        markers += [
                #('number_of_genes_1plusreads', 'n genes'),
                'Timepoint',
                #'Gender',
                #'Mousename',
                #'cellType',
                'Treatment',
                'cellSubtype',
                ]
        mgs = [x if isinstance(x, str) else x[0] for x in marker_genes]
        for ipl, (gene, ax) in enumerate(zip(markers, axs)):
            print('Plotting gene {:} of {:}'.format(ipl+1, len(markers)))
            if isinstance(gene, str):
                gene, title = gene, gene
            else:
                gene, title = gene
            if gene == 'cellSubtype':
                ind = np.arange(40)
                np.random.shuffle(ind)
                cmap = sns.color_palette('husl', n_colors=40)
                cmap = [cmap[i] for i in ind]
                csts = dsi.samplesheet['cellSubtype'].unique()
                cmap = dict(zip(csts, cmap))
                cmap[''] = (0, 0, 0)
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
            dsi.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=10,
                    alpha=0.25,
                    color_by=gene,
                    color_log=(gene in mgs + ['number_of_genes_1plusreads']),
                    cmap=cmap,
                    )
            ax.grid(False)
            ax.set_title(title)

            if gene == 'cluster':
                for com in dsi.samplesheet['cluster'].unique():
                    vsc = vs.loc[dsi.samplesheet[gene] == com]
                    xc, yc = vsc.values.mean(axis=0)
                    ax.scatter([xc], [yc], s=10, facecolor='none', edgecolor='red', marker='^')
                    ax.text(xc, yc, str(com), fontsize=8, ha='center', va='bottom')

            if gene in ('Treatment', 'Timepoint'):
                import matplotlib.lines as mlines
                d = ax._singlet_cmap
                handles = []
                labels = []
                for key, color in d.items():
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
                    ncol = 2
                else:
                    ncol = 1
                ax.legend(handles, labels, loc=loc, fontsize=6, ncol=ncol)

        fig.tight_layout()

        # The block below is triggered manually once the annotation is
        # considered satisfactory
        if False:
            print('Add subtype annotation to loom file')
            if 'cellSubtype' not in ds.samplesheet.columns:
                ds.samplesheet['cellSubtype'] = ''

            version = versions[-1]
            import loompy
            fn_good = '../../data/sequencing/datasets/all_{:}/good.loom'.format(version)
            with loompy.connect(fn_good) as dsl:
                dsl.ca['cellSubtype'] = ds.samplesheet['cellSubtype'].values
                dsl.ca['cellType'] = ds.samplesheet['cellType'].values
