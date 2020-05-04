# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/07/19
content:    Plot panels for Fig 6.
'''
import os
import sys
import glob
import gzip
import subprocess as sp
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from lungsc.ingest.load_dataset import DatasetLung, versions

fig_fdn = '../../figures/immune_paper_figs/immune_paper_figure_6/'


if __name__ == '__main__':

    macrophages = ['Mac I', 'Mac II', 'Mac III', 'Mac IV', 'Mac V']
    ds0 = DatasetLung.load(preprocess=True, version=versions[-2])
    ds = ds0.query_samples_by_metadata(
        '(cellType == "immune") & (doublet == 0) & (Treatment == "normal") & (cellSubtype in @macrophages)',
        local_dict=locals(),
    )

    if False:
        print('Feature selection')
        features = ds.feature_selection.overdispersed_within_groups('Mousename', inplace=False)
        dsf = ds.query_features_by_name(features)

        print('PCA')
        dsc = dsf.dimensionality.pca(n_dims=25, robust=False, return_dataset='samples')

        print('tSNE')
        vs = dsc.dimensionality.tsne(perplexity=30)

    print('Load tSNE from file')
    vs = pd.read_csv(
        '../../data/sequencing/datasets/all_20190828/tsne_immune.tsv',
        sep='\t',
        index_col=0,
        )
    vs = vs.loc[ds.samplenames]

    if True:
        print('Plot tSNEs')
        fns = ds.featurenames
        markers = [
            'Fcgr1', # CD64
            'Fcgr2b', # CD32
            'Fcgr3',  # also CD32
            'Fcgr4',  # CD16a
            'Fcer1a',
            'Fcer1g',
            'Fcer2a',  # CD23
            'Fcgrt',  # Fcrn
            ]
        markers += ['Timepoint', 'cellSubtype', 'Mousename']
        ds.samplesheet['before_birth'] = ds.samplesheet['Timepoint'] == 'E18.5'
        fig, axs = plt.subplots(3, 4, figsize=(14, 8))
        axs = axs.ravel()
        for gene, ax in zip(markers, axs):
            if gene == 'cellSubtype':
                cmap = {
                    'Mac IV': 'darkolivegreen',
                    'Mac II': 'lime',
                    'Mac I': 'seagreen',
                    'Mac III': 'greenyellow',
                    'Mac V': 'turquoise',
                }
            elif gene == 'Timepoint':
                cmap = {
                    'E18.5': 'navy',
                    'P1': 'gold',
                    'P7': 'tomato',
                    'P21': 'firebrick',
                    }
            elif gene == 'Mousename':
                cmap = sns.color_palette('Set1', n_colors=10)
            else:
                cmap = 'viridis'

            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=12,
                    alpha=0.20,
                    cmap=cmap,
                    color_by=gene,
                    color_log=gene in fns,
                    )

            if gene == 'cellSubtype':
                cst = ds.samplesheet['cellSubtype']
                cstu = ['Mac I', 'Mac II', 'Mac III', 'Mac IV', 'Mac V']
                for cs in cstu:
                    xm, ym = vs.loc[cst == cs].mean(axis=0)
                    ax.text(xm, ym, cs, ha='center', va='center',
                            fontsize=12)

            ax.grid(False)
            ax.set_axis_off()
            ax.set_title(gene)
        fig.tight_layout()

        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_tsne_gene_prepostbirth_{:}.{:}'.format(
                    gene, ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_gene_prepostbirth_{:}.{:}'.format(
                    gene, ext))
        plt.close(fig)

    plt.ion()
    plt.show()
