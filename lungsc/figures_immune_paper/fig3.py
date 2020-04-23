# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/07/19
content:    Plot panels for Fig 3.
'''
import os
import sys
import glob
import gzip
import pickle
import subprocess as sp
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

from lungsc.pilots.load_dataset import DatasetLung, versions


ctms = ['Mac I', 'Mac II', 'Mac III', 'Mac IV', 'Mac V']

fig_fdn = '../../figures/immune_paper_figs/immune_paper_figure_3/'


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
        print('Plot tSNE with Gal, Car4, Itgax, C1qa, Plac8, Ifitm6')

        genes = ['Gal', 'Car4', 'Itgax', 'C1qa', 'Plac8', 'Ifitm6']
        for gene in genes:
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
        print('Make table with top DE genes within macrophages')
        fn_comp = '../../data/gene_lists/immune_DEGs_macros.pkl'
        if not os.path.isfile(fn_comp):
            comps = {}
            for cst in ctms:
                print('DE for {:}'.format(cst))

                ds.samplesheet['is_focal'] = ds.samplesheet['cellSubtype'] == cst
                dsp = ds.split('is_focal')

                # Subsample
                for key in dsp:
                    if dsp[key].n_samples > 300:
                        dsp[key].subsample(300, inplace=True)

                comp = dsp[True].compare(dsp[False])
                comp['log2_fc'] = np.log2(dsp[True].counts.mean(axis=1) + 0.1) - np.log2(dsp[False].counts.mean(axis=1) + 0.1)
                comp.name = cst
                comps[cst] = comp
            del ds.samplesheet['is_focal']
            with open(fn_comp, 'wb') as f:
                pickle.dump(comps, f)

        else:
            with open(fn_comp, 'rb') as f:
                comps = pickle.load(f)

        if False:
            print('Save tables to file')
            tops = {}
            for cst in ctms:
                fn_comp_tsv = '../../data/gene_lists/immune_DEGs_{:}.tsv'.format(mp)
                comp = comps[cst]
                top = comp.loc[comp['log2_fc'] > 0].nlargest(50, 'statistic')
                tops[cst] = top
                top.to_csv(fn_comp_tsv, sep='\t', index=True)
            top_sum = pd.DataFrame([], index=np.arange(50))
            for mp, top in tops.items():
                top_sum[mp] = top.index
            fn_comp_tsv_sum = '../../data/gene_lists/immune_DEGs_macros_summary.tsv'
            top_sum.to_csv(fn_comp_tsv_sum, sep='\t', index=False)

    if True:
        print('Plot heatmap with single top DE genes')
        tops = {}
        for cst, comp in comps.items():
            tops[cst] = comp.loc[comp['log2_fc'] > 0].nlargest(5, 'statistic').index.tolist()

        genes = sum(tops.values(), [])

        genes = [
            # Common
            'Ptprc',
            'Cd68',
            'Axl',
            'Dab2',
            # Mac I
            'Gal',
            'Mcm5',
            'Mcm2',
            'Mcm3',
            'Mcm4',
            'Mcm6',
            'Bub1',
            'Plk1',
            'Top2a',
            'Mki67',
            # Mac II,
            'Car4',
            'Atp6v0d2',
            'Mgll',
            'Krt19',
            'Slc39a2',
            'Coro6',
            'Marco',
            'Siglecf',
            'Gpnmb',
            'Ear1',
            'Cd200r4',
            'Ctsk',
            'Ly75',
            'Bhlhe41',
            'Slc7a2',
            'Cdh1',
            'Pex11a',
            # Mac III
            'Itgax',
            'Adgrl3',
            # Mac IV
            'Fcrls',
            'Pf4',
            'C1qa',
            'C1qb',
            'C1qc',
            'C3ar1',
            'Tmem176b',
            'Cxcl12',
            'Ccl12',
            'Cxcl16',
            'Stab1',
            'Ms4a7',
            'Ms4a4a',
            'Igfbp4',
            'Apoe',
            'Lgmn',
            'Maf',
            # Mac V
            'Pla2g7',
            'Ifitm2',
            'Ifitm3',
            'Ifitm6',
            'Plac8',
            'Pglyrp1',
            'Serpinb10',
            'Adgre4',
            'Adgre5',
            'Napsa',
            'Rnase6',
            'Fyb',
            'Clec4a1',
            'Itga4',
            'Samhd1',

            ]
        data = pd.DataFrame([], index=genes)
        for cst in ctms:
            dsi = ds.query_samples_by_metadata(
                'cellSubtype == @cst',
                local_dict=locals(),
                )
            mat = np.log10(0.1 + dsi.counts.loc[genes]).mean(axis=1)
            data[cst] = mat

        # Normalize by max expression of that gene
        data += 1
        data = (data.T / data.max(axis=1)).T

        fig, ax = plt.subplots(figsize=(3, 10.5))
        sns.heatmap(
            data,
            ax=ax,
            cmap='plasma',
            vmin=0,
            vmax=1,
            fmt='.1f',
            xticklabels=True,
            yticklabels=True,
            cbar=False,
            )
        for tk in ax.get_yticklabels():
            tk.set_rotation(0)
        for tk in ax.get_xticklabels():
            tk.set_rotation(90)
        ax.set_xlim(0, 5)
        ax.set_ylim(len(genes), 0)
        fig.tight_layout()

    if True:
        fig.savefig(fig_fdn+'heatmap_single_genes_full.png')

    if True:
        print('Plot heatmap with pathways')

        pathways = [
            ('cell cycle', ['Ncapd2', 'Mcm5', 'Mcm7', 'Cdca8', 'Smc2']),
            ('glycolysis', ['Cluh', 'Dbi', 'Eno1', 'Ldha', 'Pkm']),
            ('lipid metabolism', ['Lpl', 'Lipa', 'Abcg1', 'Sdc4', 'Abca9', 'Abca1']),
            ('matrix\nremodelling', ['Crispld2', 'Spint1']),
            ('angiogenesis', ['Fn1', 'Il18', 'Axl', 'Gas6', 'Pf4', 'Apoe']),
            ('alveolar', ['Adgrl3', 'Tgm2', 'Clec4n', 'Pparg', 'Ear2', 'Itgax', 'Car4', 'Bhlhe41', 'Trim29']),
            ('cell migration', ['Ccr2', 'Ccr5', 'Cx3cr1', 'Cxcl16', 'Cxcl2']),
            ('antibacterial', ['Acp5', 'Mpeg1', 'Plac8', 'Rnase6', 'Lyz2']),
            ('complement', ['C1qc', 'C1qa', 'C1qb', 'C3ar1']),
            ('alternative\nactivation', ['Ms4a8a', 'Axl', 'Il18', 'Maf', 'Lgmn']),
            ('type-I IFN', ['Adgrl3', 'Ifitm6', 'Ifitm3', 'Ifi27l2a', 'Ifitm2']),
            ('neg reg of\ninflammation', ['Cd200r4', 'Gpnmb', 'Il1rn', 'Dapk1', 'Dok2', 'Cd300a', 'Nr4a1', 'Lst1']),
            ]
        genes = sum((x[1] for x in pathways), [])

        data = pd.DataFrame([], index=genes)
        for cst in ctms:
            dsi = ds.query_samples_by_metadata(
                'cellSubtype == @cst',
                local_dict=locals(),
                )
            mat = np.log10(0.1 + dsi.counts.loc[genes]).mean(axis=1)
            data[cst] = mat

        # Normalize by max expression of that gene
        data += 1
        data = (data.T / data.max(axis=1)).T

        fig, axs = plt.subplots(
                2, 1, figsize=(11, 4), sharex=True,
                gridspec_kw={'height_ratios': [1, 20]})
        sns.heatmap(
                data.iloc[:, :5].T,
                ax=axs[1],
                cmap='plasma',
                vmin=0,
                vmax=1,
                fmt='.1f',
                xticklabels=True,
                yticklabels=True,
                cbar=False,
            )
        for tk in axs[1].get_yticklabels():
            tk.set_rotation(0)
        for tk in axs[1].get_xticklabels():
            tk.set_rotation(90)
            tk.set_fontsize(8)
        axs[1].set_ylim(5, 0)
        axs[1].set_xlim(0, len(genes))
        i = 0
        for ipw, (pw, gns) in enumerate(pathways):
            if i != 0:
                axs[1].plot([i] * 2, [0, len(genes)], lw=2, color='lightgrey', alpha=0.9)
            i += len(gns)

        # Legend
        labels = ['none', 'low', 'mid', 'high']
        sfun = plt.cm.plasma
        handles = [
            axs[1].scatter([], [], marker='s', s=50, color=sfun(0)),
            axs[1].scatter([], [], marker='s', s=50, color=sfun(0.33)),
            axs[1].scatter([], [], marker='s', s=50, color=sfun(0.67)),
            axs[1].scatter([], [], marker='s', s=50, color=sfun(1.0)),
            ]
        leg = axs[1].legend(
                handles, labels,
                title='Gene\nexpression:',
                bbox_to_anchor=(1.01, 0.99),
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

            wt = i + 0.5 * w - 0.15 * w * (pw == 'lipid metabolism')
            ht = 2 + 1.5 * (ipw % 2)
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
        fig.savefig(fig_fdn+'heatmap_pathways.png')

    plt.ion()
    plt.show()
