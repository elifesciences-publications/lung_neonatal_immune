# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/07/19
content:    Plot panels for fig 2.
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

from lungsc.pilots.load_dataset import DatasetLung, versions


ctms = ['Mac I', 'Mac II', 'Mac III', 'Mac IV', 'Mac V']

fig_fdn = '../../figures/immune_paper_figs/immune_paper_figure_2/'


if __name__ == '__main__':

    ds0 = DatasetLung.load(preprocess=True, version=versions[-2])
    ds0.query_samples_by_metadata(
        '(cellType == "immune") & (doublet == 0)', inplace=True)

    ds = ds0.query_samples_by_metadata('cellSubtype in @ctms', local_dict=locals())

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
        '../../data/sequencing/datasets/all_{:}/tsne_immune.tsv'.format(versions[-2]),
        sep='\t',
        index_col=0,
        )
    vs = vs.loc[ds.samplenames]

    if True:
        print('Plot tSNE with clusters')
        gene, title = ('cellSubtype', 'subtype')
        ln = ds.samplesheet[gene].str.len().max()

        fig, ax = plt.subplots(figsize=(4.8 + 0.07 * ln, 4.2))
        cmap = {
            'Mac IV': 'darkolivegreen',
            'Mac II': 'lime',
            'Mac I': 'seagreen',
            'Mac III': 'greenyellow',
            'Mac V': 'turquoise',
            'DC I': 'darkviolet',
            'DC II': 'violet',
            'DC III': 'fuchsia',
            'B cell': 'lightcoral',
            'T cell': 'tomato',
            'NK cell': 'chocolate',
            'IL cell': 'tan',
            'mast cell': 'gold',
            'neutrophil': 'black',
        }

        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=12,
                alpha=0.30,
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
                fig.savefig(fig_fdn+'immune_tsne_metadata_{:}.{:}'.format(
                    gene, ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_metadata_{:}.{:}'.format(
                    gene, ext))

    if True:
        print('Plot tSNE with Cd68')
        gene = 'Cd68'
        title = gene

        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        ds.plot.scatter_reduced_samples(
                vs,
                ax=ax,
                s=12,
                alpha=0.30,
                cmap='viridis',
                color_by=gene,
                color_log=True,
                )
        ax.grid(False)
        ax.set_axis_off()
        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_tsne_{:}.{:}'.format(
                    gene, ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_{:}.{:}'.format(
                    gene, ext))

    if True:
        print('Plot tSNEs with Lgals3 and Itgam')
        genes = ['Lgals3', 'Itgam']
        for gene in genes:
            title = gene

            fig, ax = plt.subplots(figsize=(4.8, 4.2))
            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=12,
                    alpha=0.30,
                    cmap='viridis',
                    color_by=gene,
                    color_log=True,
                    )
            ax.grid(False)
            ax.set_title(title)
            ax.set_axis_off()
            fig.tight_layout()

            if True:
                for ext in ['svg', 'pdf', ['png', 600]]:
                    if isinstance(ext, list):
                        ext, dpi = ext
                        fig.savefig(fig_fdn+'immune_tsne_{:}.{:}'.format(
                            gene, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'immune_tsne_{:}.{:}'.format(
                            gene, ext))


    if True:
        print('Plot a few genes before/after birth separately')
        dsbb = ds.query_samples_by_metadata('Timepoint == "E18.5"')
        dsab = ds.query_samples_by_metadata('Timepoint != "E18.5"')

        genes = [
            'Timepoint',
            'Axl',
            'Dab2',
            'Ifitm6',
            'Plac8',
            ]

        fig, axs = plt.subplots(len(genes), 2, figsize=(4.5, 9.2), sharex=True, sharey=True)
        for ir, ax_row in enumerate(axs):
            gene = genes[ir]

            if gene == 'Timepoint':
                cmap = {
                    'E18.5': 'navy',
                    'P1': 'gold',
                    'P7': 'tomato',
                    'P21': 'firebrick',
                    }
            else:
                cmap = 'viridis'

            # Before birth
            ax = ax_row[0]
            dsbb.plot.scatter_reduced_samples(
                    vs.loc[dsbb.samplenames],
                    ax=ax,
                    s=8,
                    alpha=0.30,
                    cmap=cmap,
                    color_by=gene,
                    color_log=(gene != 'Timepoint'),
                    )
            ax.grid(False)
            #ax.set_axis_off()
            ax.set_xticklabels([])
            ax.set_yticklabels([])

            # After birth
            ax = ax_row[1]
            dsab.plot.scatter_reduced_samples(
                    vs.loc[dsab.samplenames],
                    ax=ax,
                    s=8,
                    alpha=0.30,
                    cmap=cmap,
                    color_by=gene,
                    color_log=(gene != 'Timepoint'),
                    )
            ax.grid(False)
            #ax.set_axis_off()
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax_row[0].set_ylabel(gene)
            ax_row[0].set_xlabel('')
            ax_row[1].set_xlabel('')
            if ir == 0:
                ax_row[0].set_title('Before birth')
                ax_row[1].set_title('After birth')

        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_tsne_metadata_{:}_split_only.{:}'.format(
                    gene, ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_metadata_{:}_split_only.{:}'.format(
                    gene, ext))

    if True:
        print('Plot cell type abundances')
        df0 = ds0.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df0['c'] = 1
        frac0 = df0.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac0 *= 100 / frac0.sum(axis=0)
        frac0 = frac0.loc[ctms]

        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1
        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)

        colors = sns.color_palette('muted', n_colors=5)
        cmap = {
            'Mac I': colors[0],
            'Mac II': colors[1],
            'Mac III': colors[2],
            'Mac IV': colors[3],
            'Mac V': colors[4],
            }

        print('Plot fraction trajectories')
        from scipy import interpolate
        fig, axs = plt.subplots(1, 2, figsize=(7, 3), sharex=True, sharey=True)
        sts = frac.index.tolist()
        if '' in sts:
            sts.remove('')
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        for iax, (ax, fraci) in enumerate(zip(axs, [frac0, frac])):
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
                ax.set_title('% of\nimmune cells')
            else:
                ax.set_title('% of\nmacrophages')
                ax.legend(title='Subtype:', bbox_to_anchor=(1.04, 1.02), loc="upper left", fontsize=9)
        axs[0].set_yticklabels(['0.1%', '1%', '10%', '100%'])

        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'macrophage_population_abundances.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'macrohpage_population_abundances.{:}'.format(
                    ext))

    if True:
        print('Print fraction of proliferating macrophages')

        df = ds.samplesheet[['Timepoint']].copy()
        df['Proliferating'] = ds.counts.loc['Mki67'] > 10
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'Proliferating']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)

        print('Plot fraction trajectories')
        from scipy import interpolate
        fig, axs = plt.subplots(1, 2, figsize=(6, 2.5), gridspec_kw={'width_ratios': [3, 1]})
        ax = axs[0]
        sts = [True]
        if '' in sts:
            sts.remove('')
        colors = sns.color_palette('muted', n_colors=len(sts))
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        for ist, st in enumerate(sts):
            y = np.maximum(0.5, frac.loc[st, xorder])

            outx = np.linspace(0, 3, 100)
            outy = interpolate.pchip_interpolate(x, y, outx)
            out = np.vstack([outx, outy])
            m = ['o', 's', '^', 'd'][ist]

            ax.scatter(
                    x, y,
                    marker=m,
                    lw=2, alpha=0.8,
                    edgecolor=colors[ist],
                    label=st,
                    facecolor='none',
                    zorder=10 - ist,
                    )
            ax.plot(out[0], out[1], lw=2,
                    alpha=0.4,
                    color=colors[ist],
                    zorder=10 - ist,
                    )
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(xorder)
        #ax.legend(
        #    title='Proliferating:',
        #    bbox_to_anchor=(1.04, 1.02),
        #    loc="upper left",
        #    fontsize=9)
        ax.set_ylabel('% of macrophages')
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax.set_ylim(0, 100)

        df = ds.samplesheet[['Timepoint', 'cellSubtype']].copy()
        df['Proliferating'] = ds.counts.loc['Mki67'] > 10
        df['c'] = 1
        df = df.loc[df['Timepoint'] == 'E18.5']
        fracc = df.groupby(['cellSubtype', 'Proliferating']).count()['c'].unstack().fillna(0).astype(np.float64).T
        fracc *= 100 / fracc.sum(axis=0)
        fracc = fracc.loc[True]
        cmap = {
            'Mac IV': 'darkolivegreen',
            'Mac II': 'lime',
            'Mac I': 'seagreen',
            'Mac III': 'greenyellow',
            'Mac V': 'turquoise',
        }

        ax = axs[1]
        order = ['Mac I', 'Mac II', 'Mac III', 'Mac IV', 'Mac V']
        for iy, ct in enumerate(order):
            ax.barh(4 - iy, fracc[ct], 0.65, color=cmap[ct], zorder=5)
        ax.set_ylim(-0.5, 4.5)
        ax.set_yticks(np.arange(5))
        ax.set_yticklabels(order[::-1])
        ax.set_xlim(-0.01, 101)
        ax.set_xticks([0, 50, 100])
        ax.set_xlabel('% proliferating\nat E18.5', ha='center')
        ax.grid(True)

        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'proliferating_population_abundance.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'proliferating_population_abundance.{:}'.format(
                    ext))

    if True:
        print('Plot proliferation genes')
        genes = [
            'Mki67',
            'Mcm5',
            'Mcm7',
            ]

        fig, axs = plt.subplots(1, len(genes), figsize=(6, 2.3), sharex=True, sharey=True)
        for ir, (gene, ax) in enumerate(zip(genes, axs)):
            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=8,
                    alpha=0.30,
                    cmap='viridis',
                    color_by=gene,
                    color_log=True,
                    )
            ax.grid(False)
            #ax.set_axis_off()
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_title(gene)

        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_tsne_proliferation.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_proliferation.{:}'.format(
                    ext))

    plt.ion()
    plt.show()
