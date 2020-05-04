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
    ds = ds0.query_samples_by_metadata(
        'cellSubtype == "Mac V"',
        )

    # See other script for pseudotiming
    fn_pt = '../../data/sequencing/datasets/all_{:}/macV_pseudotime.tsv'.format(versions[-1])
    fn_tsne = '../../data/sequencing/datasets/all_{:}/macV_tsne.tsv'.format(versions[-1])
    if not os.path.isfile(fn_pt):
        print('Feature selection')
        features = ds.feature_selection.overdispersed_within_groups(
                'Mousename',
                n_features=500,
                inplace=False,
                )
        dsf = ds.query_features_by_name(features)

        print('PCA')
        dsc = dsf.dimensionality.pca(
            n_dims=25,
            robust=False,
            return_dataset='samples',
            )

        print('tSNE')
        vs = dsc.dimensionality.tsne(perplexity=20)

        print('Convert to AnnData')
        adata = anndata.AnnData(
            dsc.counts.T,
            obs=dsc.samplesheet,
            var=dsc.featuresheet,
            )

        print('Flip the tSNE if necessary')
        if vs.loc[ds.samplesheet['Timepoint'] == 'E18.5'].mean(axis=0)['dimension 1'] > 0:
            vs['dimension 1'] = -vs['dimension 1']
        if vs.loc[ds.samplesheet['Timepoint'] == 'E18.5'].mean(axis=0)['dimension 2'] < 0:
            vs['dimension 2'] = -vs['dimension 2']

        print('Find initial cell, e.g. highest cell cycle')
        stem = vs.loc[(ds.samplesheet['Timepoint'] == 'E18.5') & (vs['dimension 2'] > 2), 'dimension 1'].idxmin()
        stem_idx = ds.samplenames.tolist().index(stem)

        print('Perform a bunch of operations before pseudotime')
        adata.uns['iroot'] = stem_idx
        scanpy.pp.neighbors(
                adata,
                n_neighbors=15,
                n_pcs=0,
                use_rep='X',
                knn=True,
                random_state=0,
                method='umap',
                metric='correlation',
                metric_kwds={},
                copy=False,
                )
        scanpy.tl.diffmap(
                adata,
                n_comps=15,
                copy=False,
                )

        print('Compute pseudotime')
        scanpy.tl.dpt(
                adata,
                n_dcs=10,
                n_branchings=0,
                min_group_size=0.01,
                allow_kendall_tau_shift=True,
                copy=False,
                )
        ds.samplesheet['pseudotime'] = adata.obs['dpt_pseudotime']

        ds.samplesheet[['pseudotime']].to_csv(fn_pt, sep='\t', index=True, header=True)
        vs.to_csv(fn_tsne, sep='\t', index=True, header=True)

    vs = pd.read_csv(fn_tsne, sep='\t', index_col=0)
    pt = pd.read_csv(fn_pt, sep='\t', index_col=0, squeeze=True)
    pt = pt.loc[ds.samplenames]
    ds.samplesheet['pseudotime'] = pt
    stem = vs.loc[(ds.samplesheet['Timepoint'] == 'E18.5') & (vs['dimension 2'] > 2), 'dimension 1'].idxmin()
    stem_idx = ds.samplenames.tolist().index(stem)

    if True:
        print('Plot tSNE with pseudotime')
        fig, axs = plt.subplots(1, 7, figsize=(17, 2.5))
        for ax, gene in zip(axs, ['Ccr1', 'Ccr2', 'Treml4', 'Ace', 'Slc12a2', 'Timepoint', 'pseudotime']):
            if gene == 'pseudotime':
                cmap = 'viridis'
            elif gene == 'Timepoint':
                cmap = {
                    'E18.5': 'navy',
                    'P1': 'gold',
                    'P7': 'tomato',
                    'P21': 'firebrick',
                    }
            else:
                cmap = 'viridis'
            ds.plot.scatter_reduced_samples(
                    vs,
                    ax=ax,
                    s=12,
                    alpha=0.30,
                    cmap=cmap,
                    color_by=gene,
                    color_log=(gene not in ('pseudotime', 'Timepoint')),
                    )
            ax.scatter(
                [vs.at[stem, 'dimension 1']], [vs.at[stem, 'dimension 2']],
                s=200,
                marker='*',
                edgecolor='k',
                facecolor='none',
                lw=2,
                )
            ax.grid(False)
            ax.set_axis_off()
            ax.set_title(gene.capitalize())

        pt = ds.samplesheet['pseudotime']
        times = np.linspace(0, 1, 6)
        traj = []
        for it in range(len(times) - 1):
            ind = (pt >= times[it])
            if it != len(times) - 1:
                ind &= (pt < times[it + 1])
            if ind.sum() == 0:
                continue
            xm, ym = vs.loc[ind].mean(axis=0)
            tm = 0.5 * (times[it] + times[it + 1])
            traj.append({
                'pseudotime': tm,
                'x': xm,
                'y': ym,
                })
        traj = pd.DataFrame(traj)
        traj.sort_values('pseudotime', inplace=True)
        if len(traj) > 2:
            ax.plot(
                traj['x'].values[:-1], traj['y'].values[:-1],
                lw=2, color='k',
                )
        if len(traj) >= 2:
            xi = traj['x'].values[-2]
            yi = traj['y'].values[-2]
            dx = traj['x'].values[-1] - xi
            dy = traj['y'].values[-1] - yi
            ax.arrow(
                xi, yi, dx, dy,
                color='k',
                lw=2,
                length_includes_head=True,
                head_width=4,
                )

        # Legend for pseudotime
        labels = ['0', '0.25', '0.5', '0.75', '1']
        sfun = plt.cm.get_cmap('viridis')
        handles = [
            axs[-1].scatter([], [], marker='s', s=50, color=sfun(0)),
            axs[-1].scatter([], [], marker='s', s=50, color=sfun(0.25)),
            axs[-1].scatter([], [], marker='s', s=50, color=sfun(0.50)),
            axs[-1].scatter([], [], marker='s', s=50, color=sfun(0.75)),
            axs[-1].scatter([], [], marker='s', s=50, color=sfun(1.0)),
            ]
        leg = axs[-1].legend(
                handles, labels,
                title='Pseudotime:',
                bbox_to_anchor=(1.08, 0.99),
                loc='upper left',
                )

        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'MacV_tsne_pseudotime.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'MacV_tsne_pseudotime.{:}'.format(
                    ext))

    if True:
        print('Plot heatmap with single genes and pseudotime bins')
        pt_bins = [
            [0, 0.225],
            [0.175, 0.425],
            [0.375, 0.625],
            [0.575, 0.825],
            [0.775, 1.01],
            ]

        genes = [
            'Ear2',
            'Treml4',
            'Cd9',
            'Lair1',
            'Ace',
            'Eno3',
            'Cd36',
            'Ccr1',
            'Sell',
            'Vcan',
            'Gas7',
            'Ly6c2',
            'Fn1',
            'Ccr2',
            'Ifitm6',
            'Cd68',
            'Plac8',
            'Ptprc',
            ]
        data = pd.DataFrame([], index=genes)
        for (bl, br) in pt_bins:
            dsi = ds.query_samples_by_metadata(
                '(pseudotime >= @bl) & (pseudotime < @br)',
                local_dict=locals(),
                )
            mat = np.log10(0.1 + dsi.counts.loc[genes]).mean(axis=1)
            data[(bl, br)] = mat

        # Normalize by max expression of that gene
        data += 1
        data = (data.T / data.max(axis=1)).T

        fig, ax = plt.subplots(figsize=(8, 3))
        sns.heatmap(
            data.T,
            ax=ax,
            cmap='plasma',
            vmin=0,
            vmax=1,
            fmt='.1f',
            xticklabels=True,
            yticklabels=True,
            cbar=False,
            )
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.set_ylabel('Pseudotime', rotation=270, labelpad=30)
        ax.arrow(
                1.023, 0.8, 0, -0.6, color='k', lw=1.5,
                head_width=0.02,
                head_length=0.08,
                clip_on=False,
                transform=ax.transAxes,
                )
        for tk in ax.get_xticklabels():
            tk.set_rotation(90)
        ax.set_xlim(0, len(genes))
        ax.set_ylim(len(pt_bins), 0)
        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'MacV_heatmap_singlegenes_pseudotime.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'MacV_heatmap_singlegenes_pseudotime.{:}'.format(
                    ext))

    if True:
        print('Plot heatmap with pathways and pseudotime bins')
        pt_bins = [
            [0, 0.225],
            [0.175, 0.425],
            [0.375, 0.625],
            [0.575, 0.825],
            [0.775, 1.01],
            ]

        pathways = [
            ('cytoskeleton', ['Actb', 'Vim', 'Dstn', 'Rap2b', 'Ckap4']),
            ('matrix\nremodeling', ['S100a10', 'Fn1', 'F13a1', 'Vcan', 'Mmp8']),
            ('reduction\noxidation', ['Sod1', 'Prdx6', 'Glrx']),
            ('chemotaxis', ['Ccr2', 'Ly6c1', 'Ccr1']),
            ('macrophage\ndifferentiation', ['Csf1r', 'Spn', 'Cd300e']),
            ('alternative\nactivation', ['Fn1', 'Tgm2', 'Cd36']),
            ('neg reg of\ninflammation', ['Ceacam1', 'Dusp5', 'Ear2', 'Cd274']),
            ('innate\nimmunity', ['Pglyrp1', 'Cd300ld', 'Treml4', 'Slc12a2', 'Ace']),
            ]
        genes = sum((x[1] for x in pathways), [])

        data = pd.DataFrame([], index=genes)
        for (bl, br) in pt_bins:
            dsi = ds.query_samples_by_metadata(
                '(pseudotime >= @bl) & (pseudotime < @br)',
                local_dict=locals(),
                )
            mat = np.log10(0.1 + dsi.counts.loc[genes]).mean(axis=1)
            data[(bl, br)] = mat

        # Normalize by max expression of that gene
        data += 1
        data = (data.T / data.max(axis=1)).T

        fig, axs = plt.subplots(
                2, 1, figsize=(8, 4), sharex=True,
                gridspec_kw={'height_ratios': [1, 20]})
        sns.heatmap(
            data.T,
            ax=axs[1],
            cmap='plasma',
            vmin=0,
            vmax=1,
            fmt='.1f',
            xticklabels=True,
            yticklabels=True,
            cbar=False,
            )
        i = 0
        for ipw, (pw, gns) in enumerate(pathways):
            if i != 0:
                axs[1].plot([i] * 2, [0, len(genes)], lw=2, color='lightgrey', alpha=0.9)
            i += len(gns)

        axs[1].yaxis.set_label_position("right")
        axs[1].yaxis.tick_right()
        axs[1].set_yticklabels([])
        axs[1].set_yticks([])
        axs[1].set_ylabel('Pseudotime', rotation=270, labelpad=28)
        axs[1].arrow(
                len(genes) + 0.8, 1, 0, 0.6 * len(pt_bins), color='k', lw=1.5,
                head_width=0.4,
                head_length=0.6,
                clip_on=False,
                )
        for tk in axs[1].get_xticklabels():
            tk.set_rotation(90)
        axs[1].set_xlim(0, len(genes))
        axs[1].set_ylim(len(pt_bins), 0)

        # Legend
        labels = ['none', 'low', 'mid', 'high']
        sfun = plt.cm.get_cmap('plasma')
        handles = [
            axs[1].scatter([], [], marker='s', s=50, color=sfun(0)),
            axs[1].scatter([], [], marker='s', s=50, color=sfun(0.33)),
            axs[1].scatter([], [], marker='s', s=50, color=sfun(0.67)),
            axs[1].scatter([], [], marker='s', s=50, color=sfun(1.0)),
            ]
        leg = axs[1].legend(
                handles, labels,
                title='Gene\nexpression:',
                bbox_to_anchor=(1.08, 0.99),
                loc='upper left',
                )

        axs[0].set_ylim(0, 1)
        axs[0].set_xlim(0, len(genes))
        color_d = dict(zip(
            (x[0] for x in pathways),
            sns.color_palette('muted', n_colors=len(pathways)),
            ))
        i = 0
        for ipw, (pw, gns) in enumerate(pathways):
            w = len(gns)
            rect = plt.Rectangle(
                    (i, 0), w, 1,
                    facecolor=color_d[pw],
                    edgecolor='none',
                    lw=0,
                    )
            axs[0].add_artist(rect)

            wt = i + 0.5 * w + 0.15 * w * (pw == '')
            ht = 2 + 1.5 * (ipw % 2) + 1.9 * (pw in ['chemotaxis', 'alternative\nactivation'])
            axs[0].text(
                    wt, ht, pw, ha='center', va='bottom',
                    fontsize=10,
                    clip_on=False,
                    )

            if ipw % 2:
                axs[0].plot(
                        [wt] * 2, [ht - 0.2, 1.2], lw=1, color='k',
                        clip_on=False,
                        )

            i += w
        axs[0].set_axis_off()
        fig.tight_layout(h_pad=0.01)

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'MacV_heatmap_pathways_pseudotime.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'MacV_heatmap_pathways_pseudotime.{:}'.format(
                    ext))



    plt.ion()
    plt.show()
