# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/04/20
content:    Pathway analysis for marker genes for Mac I-V.
'''
import os
import sys
import glob
import gzip
import subprocess as sp
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

from lungsc.ingest.load_dataset import DatasetLung, versions
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet import Dataset, CountsTable, FeatureSheet, CountsTableSparse, SampleSheet

fig_fdn = '../../figures/immune_paper_figs/immune_paper_supplementary_figures/'


if __name__ == '__main__':

    ds = DatasetLung.load(preprocess=True, version=versions[-2])
    ds.query_samples_by_metadata(
            '(cellType == "immune") & (doublet == 0) & (Treatment == "normal")',
            inplace=True)

    genes_all = 'Fgf2, Fgf9, Fgf10, Fgf18, Bmp4, Bmp5, Gdf2, Nog, Wnt2, Wnt5a, Wnt7b, Sema3c, Sema5a, Ntn1, Ntn4, Hhip, Tgfb1, Tgfb2, Tgfb3, Shh'.split(', ')

    print('Load tSNE from file')
    vs = pd.read_csv(
        '../../data/sequencing/datasets/all_20190828/tsne_immune.tsv',
        sep='\t',
        index_col=0,
        )
    vs = vs.loc[ds.samplenames]

    print('Plot known signaling genes')
    genes = [g for g in genes_all if (g in ds.featurenames)]
    markers = list(genes) + [
            'cellSubtype',
    ]

    fig, axs = plt.subplots(4, 5, figsize=(12, 8), sharex=True, sharey=True)
    axs = axs.ravel()
    for i, (gene, ax) in enumerate(zip(markers, axs)):
        if gene == 'cellSubtype':
            cmap = {
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
                'NK cell': 'chocolate',
                'IL cell': 'tan',
                'mast cell': 'gold',
                'basophil': 'grey',
                'neutrophil': 'black',
            }
        else:
            cmap = 'viridis'

        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=15,
                alpha=0.50,
                cmap=cmap,
                color_by=gene,
                color_log=(gene != 'cellSubtype'),
                )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(gene)

        if gene == 'cellSubtype':
            d = ax._singlet_cmap
            handles = []
            labels = []
            for key, color in d.items():
                h = mlines.Line2D(
                    [], [], color=color, marker='o', lw=0,
                    markersize=5,
                    )
                handles.append(h)
                labels.append(key.upper())
            ncol = 2
            fontsize = 6
            ax.legend(
                handles, labels, loc='lower left', fontsize=fontsize, ncol=ncol,
                bbox_to_anchor=(1.01, 0.01),
                bbox_transform=ax.transAxes,
                )

    fig.tight_layout()
    fig.savefig(fig_fdn+'tsne_signaling_epithelial.png')

    print('Feature selection')
    features = ds.feature_selection.overdispersed_within_groups('Mousename', inplace=False)
    dsf = ds.query_features_by_name(features)

    print('PCA')
    dsc = dsf.dimensionality.pca(n_dims=30, robust=False, return_dataset='samples')

    print('Smoothen along knn')
    edges = dsc.graph.knn('samples', 20, return_kind='edges')
    # Smoothen twice (reaching second-order neighbors)
    counts_smooth = ds.counts.smoothen_neighbors(edges, n_iterations=2)
    ds_smooth = Dataset(
            counts_table=counts_smooth,
            samplesheet=ds.samplesheet.copy(),
        )

    print('Plot known signaling genes, smoothened')
    genes = [g for g in genes_all if (g in ds.featurenames)]
    markers = list(genes) + [
            'cellSubtype',
    ]

    fig, axs = plt.subplots(4, 5, figsize=(12, 8), sharex=True, sharey=True)
    axs = axs.ravel()
    for i, (gene, ax) in enumerate(zip(markers, axs)):
        if gene == 'cellSubtype':
            cmap = {
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
                'NK cell': 'chocolate',
                'IL cell': 'tan',
                'mast cell': 'gold',
                'basophil': 'grey',
                'neutrophil': 'black',
            }
        else:
            cmap = 'viridis'

        ds_smooth.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=15,
                alpha=0.50,
                cmap=cmap,
                color_by=gene,
                color_log=(gene != 'cellSubtype'),
                )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(gene)

        if gene == 'cellSubtype':
            d = ax._singlet_cmap
            handles = []
            labels = []
            for key, color in d.items():
                h = mlines.Line2D(
                    [], [], color=color, marker='o', lw=0,
                    markersize=5,
                    )
                handles.append(h)
                labels.append(key.upper())
            ncol = 2
            fontsize = 6
            ax.legend(
                handles, labels, loc='lower left', fontsize=fontsize, ncol=ncol,
                bbox_to_anchor=(1.01, 0.01),
                bbox_transform=ax.transAxes,
                )

    fig.tight_layout()
    fig.savefig(fig_fdn+'tsne_signaling_epithelial_smoothed_twice.png')
