# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/07/19
content:    Plot panels for fig 1.
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


fig_fdn = '../../figures/immune_paper_figs/immune_paper_figure_1/'


if __name__ == '__main__':

    ds = DatasetLung.load(preprocess=True, version=versions[-2])
    ds.query_samples_by_metadata('(cellType == "immune") & (doublet == 0) & (Treatment == "normal")', inplace=True)

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
        print('Plot timepoints')
        gene, title = 'Timepoint', 'timepoint'
        ln = ds.samplesheet[gene].str.len().max()

        fig, ax = plt.subplots(figsize=(4.8 + 0.07 * ln, 4.2))
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
                alpha=0.20,
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
        print('Plot before/after birth separately')
        dsbb = ds.query_samples_by_metadata('Timepoint == "E18.5"')
        dsab = ds.query_samples_by_metadata('Timepoint != "E18.5"')

        fig, axs = plt.subplots(1, 2, figsize=(8.5 + 0.07 * ln, 4.2))
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            }

        # Before birth
        ax = axs[0]
        dsbb.plot.scatter_reduced_samples(
                vs.loc[dsbb.samplenames],
                ax=ax,
                s=8,
                alpha=0.70,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )
        ax.grid(False)
        ax.set_axis_off()

        # After birth
        ax = axs[1]
        dsab.plot.scatter_reduced_samples(
                vs.loc[dsab.samplenames],
                ax=ax,
                s=8,
                alpha=0.70,
                cmap=cmap,
                color_by=gene,
                color_log=False,
                )
        ax.grid(False)
        ax.set_axis_off()

        d = cmap
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
        if gene == 'Timepoint':
            labels_old = list(labels)
            labels = ['E18.5', 'P1', 'P7', 'P21']
            handles = [handles[labels_old.index(li)] for li in labels]
            ncol = 1
            fontsize = 8
        else:
            ncol = 1
            fontsize = 8
        ax.legend(
            handles, labels, loc='upper left', fontsize=fontsize, ncol=ncol,
            bbox_to_anchor=(1.01, 1.01),
            bbox_transform=ax.transAxes,
            )
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
        print('Plot rough cell type abundances')
        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)

        fracr = frac.loc[[]]
        fracr.loc['Mac I'] = frac.loc[
            'Mac I'
            ]
        fracr.loc['Mac II-V'] = frac.loc[[
            'Mac II',
            'Mac IV',
            'Mac III',
            'Mac V',
            ]].sum(axis=0)
        fracr.loc['Granulocytes and DCs'] = frac.loc[[
            'mast cell',
            'basophil',
            'neutrophil',
            'DC I',
            'DC II',
            'DC III',
            ]].sum(axis=0)
        fracr.loc['Lymphocytes'] = frac.loc[[
            'B cell',
            'T cell',
            'NK cell',
            'IL cell',
            ]].sum(axis=0)

        print('Plot fraction trajectories')
        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(5.5, 2.5))
        sts = fracr.index.tolist()
        if '' in sts:
            sts.remove('')
        colors = sns.color_palette('muted', n_colors=len(sts))
        x = np.arange(4)
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        for ist, st in enumerate(sts):
            y = np.maximum(0.5, fracr.loc[st, xorder])

            outx = np.linspace(0, 3, 100)
            outy = 10**interpolate.pchip_interpolate(x, np.log10(0.1 + y), outx) - 0.1
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
        ax.legend(title='Cell subtype:', bbox_to_anchor=(1.04, 1.02), loc="upper left", fontsize=9)
        ax.set_ylim(0.4, 101)
        ax.set_yscale('log')
        ax.set_ylabel('Percentage of immune cells')
        ax.set_yticks([0.5, 1, 10, 100])
        ax.set_yticklabels(['0%', '1%', '10%', '100%'])

        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_population_abundances_coarse.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_population_abundances_coarse.{:}'.format(
                    ext))

    if True:
        print('Plot dotplot of selected genes')
        genes_dotplot = [
            #'Cd14',
            #'Fcgr3',
            'Cd68',
            'Gal',
            'Itgax',
            'Car4',
            #'Ccl12',
            'C1qa',
            'Plac8',
            'Batf3',
            #'Ifi205',
            #'Zmynd15',
            'Itgae',
            'Cd209a',
            'Mreg',
            'Mcpt4',
            'Mcpt8',
            #'Stfa1',
            'Retnlg',
            #'Ifitm6',
            'Ms4a1',
            'Gzma',
            'Cd3e',
            'Areg',
            ]
        group_order = [
            ('Mac I', 'Mac I'),
            ('Mac II', 'Mac II'),
            ('Mac III', 'Mac III'),
            ('Mac IV', 'Mac IV'),
            ('Mac V', 'Mac V'),
            ('DC I', 'DC I'),
            ('DC II', 'DC II'),
            ('DC III', 'DC III'),
            ('mast cell', 'mast cell'),
            ('basophil', 'basophil'),
            ('neutrophil', 'neutrophil'),
            ('B cell', 'B cell'),
            ('NK cell', 'NK cell'),
            ('T cell', 'T cell'),
            ('IL cell', 'IL cell'),
            ]
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
        fig, ax = plt.subplots(figsize=(5, 7))
        ds.plot.dot_plot(
            axis='samples',
            group_by='cellSubtype',
            group_order=[x[0] for x in group_order],
            plot_list=genes_dotplot[::-1],
            ax=ax,
            layout='vertical',
            cmap='viridis',
            tight_layout=False,
            vmax='max_single',
            vmin='min_single',
            min_size=0,
            )
        ax.set_xticklabels([x[1] for x in group_order], rotation=90)
        for ict, (ct, _) in enumerate(group_order):
            r = plt.Rectangle(
                (ict - 0.5, -0.5), 1, len(genes_dotplot),
                facecolor=cmap[ct],
                edgecolor='none',
                alpha=0.1,
                zorder=0.5,
                )
            ax.add_patch(r)
        # Macro
        r = plt.Rectangle(
            (-0.5, -0.5 + 11), 5, 6,
            facecolor='none',
            edgecolor='seagreen',
            alpha=0.6,
            zorder=10,
            lw=2,
            clip_on=False,
            )
        ax.add_patch(r)
        # DC
        r = plt.Rectangle(
            (-0.5 + 5, -0.5 + 7), 3, 4,
            facecolor='none',
            edgecolor='darkmagenta',
            alpha=0.6,
            zorder=10,
            lw=2,
            clip_on=False,
            )
        ax.add_patch(r)
        # Granulocytes
        r = plt.Rectangle(
            (-0.5 + 8, -0.5 + 4), 3, 3,
            facecolor='none',
            edgecolor='dimgrey',
            alpha=0.9,
            zorder=10,
            lw=2,
            clip_on=False,
            )
        ax.add_patch(r)
        # Lympho
        r = plt.Rectangle(
            (-0.5 + 11, -0.5), 4, 4,
            facecolor='none',
            edgecolor='tomato',
            alpha=1.0,
            zorder=10,
            lw=2,
            clip_on=False,
            )
        ax.add_patch(r)

        fig.text(0.26, 0.95, 'M', ha='center', va='bottom')
        fig.text(0.405, 0.95, 'DC', ha='center', va='bottom')
        fig.text(0.52, 0.95, 'G', ha='center', va='bottom')
        fig.text(0.65, 0.95, 'L', ha='center', va='bottom')

        # Legends
        sfun = ax._singlet_dotmap['fraction_size_map']
        handles1 = [
            ax.scatter([], [], color='grey', s=sfun(0)),
            ax.scatter([], [], color='grey', s=sfun(0.1)),
            ax.scatter([], [], color='grey', s=sfun(0.25)),
            ax.scatter([], [], color='grey', s=sfun(0.5)),
            ax.scatter([], [], color='grey', s=sfun(0.75)),
            ax.scatter([], [], color='grey', s=sfun(1)),
            ]
        labels1 = ['0%', '10%', '25%', '50%', '75%', '100%']
        leg1 = ax.legend(
                handles1, labels1,
                title='Fraction of\npopulation:    ',
                bbox_to_anchor=(1.01, 1.0),
                loc='upper left',
                )
        cfun = ax._singlet_dotmap['level_color_map']

        handles2 = [
            ax.scatter([], [], color=cfun(0.0), s=sfun(1)),
            ax.scatter([], [], color=cfun(0.3), s=sfun(1)),
            ax.scatter([], [], color=cfun(0.6), s=sfun(1)),
            ax.scatter([], [], color=cfun(1.0), s=sfun(1)),
            ]
        labels2 = ['None', 'Low', 'Medium', 'High']
        leg2 = ax.legend(
                handles2, labels2,
                title='Expression\ncompared to\nhighest\npopulation:',
                bbox_to_anchor=(1.01, 0.68),
                loc='upper left',
                )
        ax.add_artist(leg1)

        fig.tight_layout(rect=(0, 0, 1, 0.96))

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_dotplot.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_dotplot.{:}'.format(
                    ext))

    plt.ion()
    plt.show()
