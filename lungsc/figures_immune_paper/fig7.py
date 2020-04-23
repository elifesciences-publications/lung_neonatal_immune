# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/12/19
content:    Make fig7 for the paper.
'''
import os
import sys
import glob
import gzip
import subprocess as sp
import numpy as np
import pandas as pd
from scipy import interpolate

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

from lungsc.ingest.load_dataset import DatasetLung, versions


fig_fdn = '../../figures/immune_paper_figs/immune_paper_figure_7/'


if __name__ == '__main__':

    ds = DatasetLung.load(preprocess=True, version=versions[-1], include_hyperoxia=True)
    ds.query_samples_by_metadata('(cellType == "immune") & (doublet == 0)', inplace=True)
    ds7 = ds.query_samples_by_metadata('(Timepoint == "P7") & (cellSubtype != "")')

    if False:
        print('Population abundance changes')
        df = ds.samplesheet[['Timepoint', 'cellSubtype', 'Gender', 'Treatment']].copy()
        for sn, row in df.iterrows():
            if row['Treatment'] == 'hyperoxia':
                df.at[sn, 'Timepoint'] = 'hyperoxia'
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)

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
            'mast cell II': 'orange',
        }
        sts = frac.index.tolist()
        if '' in sts:
            sts.remove('')
        xorder = ['P7', 'hyperoxia']
        x = np.arange(len(xorder))
        groups = [
            ('Granulocytes', ('neutrophil', 'mast cell', 'basophil', 'mast cell II')),
            ('Lymphocytes', ('T cell', 'B cell', 'NK cell', 'IL cell')),
            ('Macrophages', ('Mac I', 'Mac II', 'Mac III', 'Mac IV', 'Mac V')),
            ('DCs', ('DC I', 'DC II', 'DC III')),
            ]
        fig, axs = plt.subplots(4, 1, figsize=(4, 7), sharex=True, sharey=True)
        axs = axs.ravel()
        for ist, (gname, group) in enumerate(groups):
            ax = axs[ist]
            ax.set_title(gname)

            for st in group:
                y = np.log10(0.1 + frac.loc[st, xorder])
                ls = '-'
                m = 'o'
                ax.scatter(
                        x, y,
                        marker=m,
                        lw=2, alpha=0.8,
                        edgecolor=cmap[st],
                        facecolor='none',
                        zorder=10 - ist,
                        label=st,
                        )

                x0, x1, y0, y1 = x[-2], x[-1], y[-2], y[-1]
                # Cut the first and last bits
                th = 0.05
                xs = x0 + th * (x1 - x0)
                ys = y0 + th * (y1 - y0)
                xe = x0 + (1 - th) * (x1 - x0)
                ye = y0 + (1 - th) * (y1 - y0)
                dx = xe - xs
                dy = ye - ys
                ax.arrow(
                    xs, ys, dx, dy,
                    facecolor=cmap[st],
                    edgecolor='none',
                    width=0.045,
                    length_includes_head=True,
                    alpha=0.9,
                    )

            ax.legend(
                loc='upper left',
                bbox_to_anchor=(1.01, 0.99),
                bbox_transform=ax.transAxes,
                )
            ax.grid(True)
            ax.set_xticks(x)
            ax.set_xticklabels(xorder[:-1] + ['HO'])
            ax.set_ylim(-1, 2.01)
            ax.set_yticks([-1, 0, 1, 2])
            ax.set_yticklabels(['0.1%', '1%', '10%', '100%'])
            ax.set_xlim(-0.4, 1.4)

        fig.text(0.06, 0.5, 'Percentage of immune cells', rotation=90, ha='center', va='center')
        fig.tight_layout(rect=(0.08, 0, 1, 1))

    if False:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'immune_population_changes_hyperoxia.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_population_changes_hyperoxia.{:}'.format(
                    ext))

    if False:
        print('Plot tSNE at P7 with hyperoxia')

        if False:
            print('Feature selection')
            features = ds7.feature_selection.overdispersed_within_groups('Mousename', inplace=False)
            dsf = ds7.query_features_by_name(features)

            print('PCA')
            dsc = dsf.dimensionality.pca(n_dims=30, robust=False, return_dataset='samples')

            print('tSNE')
            vs = dsc.dimensionality.tsne(perplexity=20)
            vs.to_csv('../../data/sequencing/datasets/all_20191204/tsne_P7_with_hyperoxia_immune.tsv', sep='\t', index=True)
        else:
            vs = pd.read_csv('../../data/sequencing/datasets/all_20191204/tsne_P7_with_hyperoxia_immune.tsv', sep='\t', index_col=0)

        vs = vs.loc[ds7.samplenames]

        print('Plot dimensionality reduction of dataset')
        fig = plt.figure(figsize=(7, 7))
        gs = fig.add_gridspec(4, 4)
        axs = []
        axs.append(fig.add_subplot(gs[0, 0]))
        axs.append(fig.add_subplot(gs[1, 0], sharex=axs[0], sharey=axs[0]))
        axs.append(fig.add_subplot(gs[2, 0], sharex=axs[0], sharey=axs[0]))
        axs.append(fig.add_subplot(gs[3, 0], sharex=axs[0], sharey=axs[0]))
        axs.append(fig.add_subplot(gs[0, 1], sharex=axs[0], sharey=axs[0]))
        axs.append(fig.add_subplot(gs[1, 1], sharex=axs[0], sharey=axs[0]))
        axs.append(fig.add_subplot(gs[2, 1], sharex=axs[0], sharey=axs[0]))
        axs.append(fig.add_subplot(gs[3, 1], sharex=axs[0], sharey=axs[0]))
        axs.append(fig.add_subplot(gs[0:2, 2:4], sharex=axs[0], sharey=axs[0]))
        axs.append(fig.add_subplot(gs[2:4, 2:4], sharex=axs[0], sharey=axs[0]))
        marker_genes = [
                #'Cd3e',
                #'Gzma',
                #'Arg1',
                #'Ms4a1',
                #'Cd68',
                #'Gal',
                #'Itgax',
                #'Car4',
                #'C1qa',
                #'Plac8',
                #'Batf3',
                #'Itgae',
                #'Mreg',
                #'Cd209c',
                #'Stfa2',
                #'Cpa3',
                #'Mcpt4',
                #'Mki67',
                'Notch2',
                'Slc39a1',
                'Fgf16',
                'Sod1',
                'Fli1',
                'Atp6ap1',
                ]
        markers = marker_genes + [
                #('number_of_genes_1plusreads', 'n genes'),
                'Gender',
                'Mousename',
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
                    'mast cell II': 'orange',
                }
            elif gene == 'Timepoint':
                cmap = {
                    'E18.5': 'navy',
                    'P1': 'gold',
                    'P7': 'tomato',
                    'P21': 'firebrick',
                }
            else:
                cmap = 'viridis'
            ds7.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=10,
                    alpha=0.2,
                    color_by=gene,
                    color_log=(gene in mgs + ['number_of_genes_1plusreads']),
                    cmap=cmap,
                    )
            ax.grid(False)
            ax.set_title(title)
            ax.set_axis_off()

            if gene == 'cellSubtype':
                dmoves = {
                    'mast cell II': {'y': -7},
                    'Mac I': {'y': +2},
                    'NK cell': {'x': +9},
                    'Mac IV': {'x': +4},
                    'T cell': {'x': +3},
                    'DC I': {'x': -4},
                    'basophil': {'x': +2},
                    }
                for com in ds7.samplesheet['cellSubtype'].unique():
                    vsc = vs.loc[ds7.samplesheet[gene] == com]
                    xc, yc = vsc.values.mean(axis=0)
                    xt, yt = xc, yc + 3
                    if com in dmoves:
                        dmove = dmoves[com]
                        if 'x' in dmove:
                            xt += dmove['x']
                        if 'y' in dmove:
                            yt += dmove['y']
                    ax.scatter([xc], [yc], s=10, facecolor='none', edgecolor='red', marker='^')
                    ax.text(xt, yt, str(com), fontsize=10, ha='center', va='bottom')

            if gene == 'Treatment':
                d = ax._singlet_cmap
                d2 = {'hyperoxia': 'HO', 'normal': 'NO'}
                handles = []
                labels = []
                for key in ['normal', 'hyperoxia']:
                    color = d[key]
                    h = mlines.Line2D(
                        [], [], color=color, marker='o', lw=0,
                        markersize=5,
                        )
                    handles.append(h)
                    loc = 'upper left'
                    labels.append(d2[key])
                if gene == 'Timepoint':
                    labels_old = list(labels)
                    labels = ['E18.5', 'P1', 'P7', 'P21']
                    handles = [handles[labels_old.index(li)] for li in labels]
                    ncol = 2
                else:
                    ncol = 1
                ax.legend(handles, labels, loc=loc, fontsize=10, ncol=ncol)

        fig.tight_layout()

    if False:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'tsne_P7_with_hyperoxia.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'tsne_P7_with_hyperoxia.{:}'.format(
                    ext))

    if True:
        print('Plot distributions of expression')
        dsu = ds7.split('Treatment')
        genes = [
            'Fgf16',
            'Sod1',
            'Notch2',
            'Slc39a1',
            'Fli1',
            'Atp6ap1',
        ]
        fig, axs = plt.subplots(6, 1, figsize=(2.25, 7), sharex=True)
        axs = axs.ravel()
        colors = {
            'normal': 'gold',
            'hyperoxia': plt.cm.viridis(0.0),
            }
        for i in range(len(axs)):
            ax = axs[i]
            gene = genes[i]

            x0 = np.log10(0.1 + dsu['normal'].counts.loc[gene].values)
            x1 = np.log10(0.1 + dsu['hyperoxia'].counts.loc[gene].values)

            sns.kdeplot(x0, ax=ax, color=colors['normal'], bw=0.3, clip=(-1.1, 6), label='Normoxia')
            sns.kdeplot(x1, ax=ax, color=colors['hyperoxia'], bw=0.3, clip=(-1.1, 6), label='Hyperoxia')
            if i == 0:
                ax.legend(
                    loc='lower center',
                    bbox_to_anchor=(0.5, 1.01),
                    bbox_transform=ax.transAxes,
                    ncol=1,
                    )
            else:
                h = ax.legend()
                h.remove()

            ax.grid(True)
            ax.set_xlim(-1.1, 5)
            ax.set_xticks([-1, 0, 2, 4])
            ax.set_xticklabels(['$0$', '$1$', '$10^2$', '$10^4$'])
            ax.set_yticklabels([])
            ax.set_ylabel(gene, rotation=0, ha='right', va='center')
            if i == len(axs) - 1:
                ax.set_xlabel('Gene expression\n[cpm]', ha='center')

        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'distributions_P7_with_hyperoxia.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'distributions_P7_with_hyperoxia.{:}'.format(
                    ext))

    if True:
        print('Plot distributions of expression by cell type')
        dsu = ds7.split(['Treatment', 'cellSubtype'])
        csts = ['basophil', 'mast cell', 'IL cell', 'B cell', 'Mac II', 'Mac III', 'Mac IV', 'Mac V']
        genes = [
            'Fgf16',
            'Sod1',
            'Notch2',
            'Slc39a1',
            'Fli1',
            'Atp6ap1',
        ]
        colors = {
            'normal': 'gold',
            'hyperoxia': plt.cm.viridis(0.0),
            }
        combos = {
            'Fgf16': csts,
            'Fli1': ['basophil', 'IL cell'],
            'Notch2': ['IL cell'],
            'Slc39a1': ['IL cell', 'Mac II', 'Mac III', 'Mac IV', 'Mac V', 'B cell'],
            'Sod1': ['mast cell', 'Mac II', 'B cell'],
            'Atp6ap1': ['basophil'],
            }
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
            'mast cell II': 'orange',
        }
        fig, axs = plt.subplots(len(genes), len(csts), figsize=(10, 7), sharex=True)
        for icol in range(len(csts)):
            for irow in range(len(genes)):
                ax = axs[irow, icol]
                cst = csts[icol]
                gene = genes[irow]

                x0 = np.log10(0.1 + dsu[('normal', cst)].counts.loc[gene].values)
                x1 = np.log10(0.1 + dsu[('hyperoxia', cst)].counts.loc[gene].values)

                alpha = 0.4
                if (gene in combos) and (cst in combos[gene]):
                    alpha += 0.5
                sns.kdeplot(x0, ax=ax, alpha=alpha, color=colors['normal'], bw=0.3, clip=(-1.1, 6), label='Normoxia')
                sns.kdeplot(x1, ax=ax, alpha=alpha, color=colors['hyperoxia'], bw=0.3, clip=(-1.1, 6), label='Hyperoxia')
                if (icol == 0) and False:
                    ax.legend(
                        loc='lower center',
                        bbox_to_anchor=(0.5, 1.01),
                        bbox_transform=ax.transAxes,
                        ncol=1,
                        )
                else:
                    h = ax.legend()
                    h.remove()

                ax.grid(True)
                ax.set_xlim(-1.1, 5)
                ax.set_xticks([-1, 0, 2, 4])
                ax.set_xticklabels(['$0$', '$1$', '$10^2$', '$10^4$'])
                ax.set_yticklabels([])
                if icol == 0:
                    ax.set_ylabel(gene, rotation=0, ha='right', va='center')
                if irow == 0:
                    ax.set_title(cst)
                    r = plt.Rectangle(
                        (0, 1.03), 1, 0.28,
                        facecolor='none',
                        edgecolor=cmap[cst],
                        lw=3,
                        transform=ax.transAxes,
                        clip_on=False,
                        )
                    ax.add_artist(r)

        fig.text(0.55, 0.03, 'Gene expression\n[cpm]', ha='center')
        fig.tight_layout(rect=(0, 0.05, 1, 1))

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'distributions_bycelltype_P7_with_hyperoxia.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'distributions_bycelltype_P7_with_hyperoxia.{:}'.format(
                    ext))

    plt.ion()
    plt.show()
