# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/09/19
content:    Plot panels for Fig 5.
'''
import os
import sys
import glob
import gzip
import pickle
import subprocess as sp
import numpy as np
import pandas as pd
import anndata
import scanpy

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

from lungsc.pilots.load_dataset import DatasetLung, versions

fig_fdn = '../../figures/immune_paper_figs/immune_paper_figure_5/'


if __name__ == '__main__':

    ds0 = DatasetLung.load(preprocess=True, version=versions[-1])
    ds0.query_samples_by_metadata(
        '(cellType == "immune") & (doublet == 0)', inplace=True)

    print('Dendritic cells')
    ds = ds0.query_samples_by_metadata(
        'cellSubtype in ("DC I", "DC II", "DC III")',
        )

    print('Load tSNE from file')
    vs = pd.read_csv(
        '../../data/sequencing/datasets/all_20190828/tsne_immune.tsv',
        sep='\t',
        index_col=0,
        )
    vs = vs.loc[ds.samplenames]

    if True:
        print('Plot tSNE with clusters')
        gene, title = ('cellSubtype', 'subtype')
        ln = ds.samplesheet[gene].str.len().max()

        fig, ax = plt.subplots(figsize=(2 + 0.07 * ln, 2))
        cmap = {
            'DC I': 'darkviolet',
            'DC II': 'violet',
            'DC III': 'fuchsia',
        }
        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=12,
                alpha=0.40,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )

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
            labels.append(key.upper())
        ncol = 1
        fontsize = 8
        ax.legend(
            handles, labels, loc='upper left', fontsize=fontsize, ncol=ncol,
            bbox_to_anchor=(1.01, 1.01),
            bbox_transform=ax.transAxes,
            )

        ax.grid(False)
        ax.set_axis_off()
        fig.tight_layout()

    if False:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_tsne_clusters.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_clusters.{:}'.format(
                    ext))

    if True:
        print('Plot timepoints')
        gene, title = 'Timepoint', 'timepoint'
        ln = ds.samplesheet[gene].str.len().max()

        fig, ax = plt.subplots(figsize=(2 + 0.07 * ln, 2))
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            }
        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=12,
                alpha=0.40,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )

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
            labels.append(key.upper())
        labels_old = list(labels)
        labels = ['E18.5', 'P1', 'P7', 'P21']
        handles = [handles[labels_old.index(li)] for li in labels]
        ncol = 1
        fontsize = 8
        ax.legend(
            handles, labels, loc='upper left', fontsize=fontsize, ncol=ncol,
            bbox_to_anchor=(1.01, 1.01),
            bbox_transform=ax.transAxes,
            )

        ax.grid(False)
        ax.set_axis_off()
        fig.tight_layout()

    if False:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_tsne_timepoint.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_timepoint.{:}'.format(
                    ext))

    if True:
        print('Plot tsne with genes')
        genes = ['Itgae', 'Mreg', 'Cd209a', 'Zbtb46', 'Flt3', 'Itgax',
                 'Cacnb3', 'Fscn1', 'Ccl5', 'Ccr7']
        fig, axs = plt.subplots(2, 5, figsize=(8.5, 4))
        axs = axs.ravel()
        for ax, gene in zip(axs, genes):
            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=12,
                    alpha=0.50,
                    cmap='viridis',
                    color_by=gene,
                    color_log=True,
                    )
            ax.grid(False)
            ax.set_title(gene)
            ax.set_axis_off()
        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_tsne_genes_DCs.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_genes_DCs.{:}'.format(
                    ext))

    if True:
        print('Plot cell type abundances')
        ctms = ['DC I', 'DC II', 'DC III']
        df0 = ds0.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df0['c'] = 1
        frac0 = df0.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac0 *= 100 / frac0.sum(axis=0)
        frac0 = frac0.loc[ctms]

        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1
        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)

        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1
        numb = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.int64).T

        colors = sns.color_palette('muted', n_colors=5)
        cmap = {
            'DC I': colors[0],
            'DC II': colors[1],
            'DC III': colors[2],
            }

        print('Plot fraction trajectories')
        from scipy import interpolate
        fig = plt.figure(figsize=(11, 3))
        axs = []
        axs.append(fig.add_subplot(1, 3, 1))
        axs.append(fig.add_subplot(1, 3, 2, sharex=axs[0], sharey=axs[0]))
        axs.append(fig.add_subplot(1, 3, 3))
        sts = frac.index.tolist()
        if '' in sts:
            sts.remove('')
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        for iax, (ax, fraci) in enumerate(zip(axs, [frac0, frac, numb])):
            for ist, st in enumerate(sts):
                y = fraci.loc[st, xorder]

                outx = np.linspace(0, 3, 100)
                outy = 10**interpolate.pchip_interpolate(x, np.log10(0.1 + y), outx) - 0.1
                out = np.vstack([outx, outy])
                m = ['o', 's', '^', 'd', 'v'][ist]

                ax.scatter(
                        x, y + 0.1,
                        marker=m,
                        lw=2, alpha=0.8,
                        edgecolor=cmap[st],
                        facecolor='none',
                        label=st,
                        zorder=10 - ist,
                        )
                ax.plot(out[0], out[1] + 0.1, lw=2,
                        alpha=0.4,
                        color=cmap[st],
                        zorder=10 - ist,
                        )
            ax.grid(True)
            ax.set_xticks(x)
            ax.set_xticklabels(xorder)
            ax.set_ylim(0.1, 101)
            ax.set_yscale('log')
            ax.set_yticks([0.1, 1, 10, 100])
            if iax == 0:
                ax.set_title('Percentage of\nimmune cells')
            elif iax == 1:
                ax.set_title('Percentage of\nDCs')
            else:
                ax.set_title('Number of\nDCs')
                ax.legend(title='Subtype:', bbox_to_anchor=(1.04, 1.02), loc="upper left", fontsize=9)
        axs[0].set_yticklabels(['0.1%', '1%', '10%', '100%'])
        axs[-1].set_yticklabels(['0', '1', '10', '100'])

        fig.tight_layout()

    if False:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'DC_population_abundances.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'DC_population_abundances.{:}'.format(
                    ext))

    print('Granulocytes')
    ds = ds0.query_samples_by_metadata(
        'cellSubtype in ("basophil", "mast cell", "neutrophil")',
        )

    print('Load tSNE from file')
    vs = pd.read_csv(
        '../../data/sequencing/datasets/all_20190828/tsne_immune.tsv',
        sep='\t',
        index_col=0,
        )
    vs = vs.loc[ds.samplenames]

    if True:
        print('Plot tSNE with clusters')
        gene, title = ('cellSubtype', 'subtype')
        ln = ds.samplesheet[gene].str.len().max()

        fig, ax = plt.subplots(figsize=(2 + 0.07 * ln, 2))
        cmap = {
            'mast cell': 'gold',
            'basophil': 'grey',
            'neutrophil': 'black',
        }
        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=12,
                alpha=0.40,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )

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
            labels.append(key.upper())
        ncol = 1
        fontsize = 8
        ax.legend(
            handles, labels, loc='upper left', fontsize=fontsize, ncol=ncol,
            bbox_to_anchor=(1.01, 1.01),
            bbox_transform=ax.transAxes,
            )

        ax.grid(False)
        ax.set_axis_off()
        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_tsne_clusters_granulocytes.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_clusters_granulocytes.{:}'.format(
                    ext))


    if True:
        print('Plot tSNE with genes')
        genes = [
            'Mcpt4', 'Mcpt8', 'Il6', 'Ccl4', 'Hgf', 'Osm',
            'Stfa1', 'Stfa2', 'Retnlg', 'S100a9', 'S100a8',
            ]
        for gene in genes:
            fig, ax = plt.subplots(figsize=(2, 2))
            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=12,
                    alpha=0.40,
                    cmap='viridis',
                    color_by=gene,
                    color_log=True,
                    )
            ax.grid(False)
            ax.set_title(gene)
            ax.set_axis_off()
            fig.tight_layout()

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'immune_tsne_clusters_granulocytes_{:}.{:}'.format(
                            gene, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'immune_tsne_clusters_granulocytes{:}.{:}'.format(
                            gene, ext))

    plt.ion()
    plt.show()
