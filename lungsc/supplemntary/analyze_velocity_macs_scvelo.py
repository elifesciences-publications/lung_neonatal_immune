# vim: fdm=indent
'''
author:     Fabio Zanini
date:       30/04/20
content:    Merge loom files for RNA velocity from Mac I-III cells
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

from lungsc.ingest.load_dataset import DatasetLung, versions


fig_fdn = '../../figures/immune_paper_figs/immune_paper_supplementary_figures/'


if __name__ == '__main__':

    ds = DatasetLung.load(preprocess=True, version=versions[-2])
    ds.query_samples_by_metadata('(cellType == "immune") & (doublet == 0) & (Treatment == "normal") & (cellSubtype in ("Mac I", "Mac II", "Mac III"))', inplace=True)

    print('Feature selection')
    features = ds.feature_selection.overdispersed_within_groups('Mousename', inplace=False)
    dsf = ds.query_features_by_name(features)

    print('PCA')
    dsc = dsf.dimensionality.pca(n_dims=25, robust=False, return_dataset='samples')

    print('Knn')
    knn = dsc.graph.knn('samples', n_neighbors=10, return_kind='sparse')

    print('Load tSNE from file')
    vs = pd.read_csv(
        '../../data/sequencing/datasets/all_20190828/tsne_immune.tsv',
        sep='\t',
        index_col=0,
        )
    vs = vs.loc[ds.samplenames]
    ds.samplesheet['tSNE1'] = vs['dimension 1']
    ds.samplesheet['tSNE2'] = vs['dimension 2']

    print('Load combined velocity file')
    fn_combined = '../../data/sequencing/datasets/all_{:}/velocity_MacI-III.loom'.format(versions[-2])
    import scvelo as scv
    adata = scv.read(fn_combined, cache=True)
    adata.var_names_make_unique()

    print('Follow tutorial')
    # show proportions of spliced/unspliced abundances
    #scv.utils.show_proportions(adata)

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=25, n_neighbors=10)

    scv.tl.velocity(adata)

    scv.tl.velocity_graph(adata)

    # Use embedding from the paper
    adata.obsm['X_tsne'] = vs.loc[adata.obs_names].values

    fig, ax = plt.subplots(figsize=(4.2, 4))
    ds.plot.scatter_reduced(
        ('tSNE1', 'tSNE2'),
        color_by='cellSubtype',
        cmap={
            'Mac I': 'seagreen',
            'Mac II': 'lime',
            'Mac III': 'greenyellow',
            },
        color_log=False,
        ax=ax,
        s=30,
        alpha=0.7,
        )
    scv.pl.velocity_embedding_stream(
        adata,
        basis='tsne',
        size=0,
        ax=ax,
        alpha=0.5,
        )
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(fig_fdn+'velocity_stochastic_MacI-III.png')

    plt.ion()
    plt.show()
