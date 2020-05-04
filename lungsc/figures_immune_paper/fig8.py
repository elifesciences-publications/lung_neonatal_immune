# vim: fdm=indent
'''
author:     Fabio Zanini
date:       18/09/19
content:    Plot panels for Fig 8.
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

fig_fdn = '../../figures/immune_paper_figs/immune_paper_figure_8/'


if __name__ == '__main__':

    ds0 = DatasetLung.load(preprocess=True, version=versions[-2])
    ds0.query_samples_by_metadata(
        '(cellType == "immune") & (doublet == 0)', inplace=True)

    if True:
        print('Plot abundances of B/T cells')
        df = ds0.samplesheet[['Timepoint', 'cellSubtype', 'Gender']].copy()
        df['c'] = 1

        frac = df.groupby(['Timepoint', 'cellSubtype']).count()['c'].unstack().fillna(0).astype(np.float64).T
        frac *= 100 / frac.sum(axis=0)

        fracr = frac.loc[['B cell', 'T cell']]

        print('Plot fraction trajectories')
        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(3, 2.5))
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
        ax.legend(
            title='Cell subtype:',
            loc='lower right',
            #bbox_to_anchor=(1.04, 1.02), loc="upper left",
            fontsize=9)
        ax.set_ylim(0.4, 101)
        ax.set_yscale('log')
        ax.set_ylabel('Percentage of immune cells')
        ax.set_yticks([0.5, 1, 10, 100])
        ax.set_yticklabels(['0%', '1%', '10%', '100%'])

        fig.tight_layout()

    if False:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'lymphocyte_population_abundances.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'lymphocyte_population_abundances.{:}'.format(
                    ext))

    ds = ds0.query_samples_by_metadata(
        'cellSubtype == "T cell"',
        )

    #if True:
    #    print('Percentage of CD4/CD8 T cells of all immune cells')
    #    ni = ds0.samplesheet.groupby(['cellType', 'Timepoint', 'Treatment', 'doublet']).count().iloc[:, 0].unstack()[0].unstack()['normal'].unstack().loc['immune']

    #    df = ds.counts.loc[['Cd4']].T > 0
    #    df['Cd8'] = (ds.counts.loc[['Cd8a', 'Cd8b1']] > 0).any()
    #    df['Timepoint'] = ds.samplesheet['Timepoint']
    #    df['c'] = 1

    #    tab = df.groupby(['Cd4', 'Cd8', 'Timepoint']).count().iloc[:, 0].unstack().fillna(0).astype(int)

    #    fig, axs = plt.subplots(2, 1, figsize=(5, 7), sharex=True)
    #    ax = axs[0]
    #    x = np.arange(4)
    #    xl = ['E18.5', 'P1', 'P7', 'P21']
    #    colors = sns.color_palette('muted', n_colors=4)
    #    for ir, ((cd4i, cd8i), row) in enumerate(tab.iterrows()):
    #        y = row.loc[xl].values
    #        label = 'Cd4' + ('+' if cd4i else '-') + ' Cd8' + ('+' if cd8i else '-')
    #        ax.plot(x, y, 'o-', lw=2, color=colors[ir], alpha=0.5, label=label)
    #    ax.set_xticks(x)
    #    ax.grid(True)
    #    ax.set_ylabel('Number of T cells')
    #    ax.legend()

    #    ax = axs[1]
    #    tabi = 100 * tab / ni
    #    sums = pd.Series(np.zeros(tabi.shape[1]), index=tabi.columns)
    #    for ir, ((cd4i, cd8i), row) in enumerate(tabi.iterrows()):
    #        y = row.loc[xl]
    #        y0 = sums.loc[xl]
    #        y1 = y0 + y
    #        sums += y
    #        #ax.plot(x, y, 'o-', lw=2, color=colors[ir], alpha=0.5)
    #        ax.fill_between(x, y0, y1, color=colors[ir], alpha=0.5)
    #    ax.plot(x, tabi.sum(axis=0).loc[xl], 'o-', lw=2, color='grey', alpha=0.5)
    #    ax.set_xticks(x)
    #    ax.set_xticklabels(xl)
    #    ax.set_yticks([0, 5, 10, 15])
    #    ax.set_yticklabels(['0%', '5%', '10%', '15%'])
    #    ax.grid(True)
    #    ax.set_ylabel('Percentage of all immune cells')

    #    fig.tight_layout()

    #if True:
    #    print('Plot violins of expression of CD4 and CD8')
    #    df = np.log10(0.1 + ds.counts.loc[['Cd4', 'Cd8a', 'Cd8b1']])
    #    df = df.stack().to_frame()
    #    df.rename(columns={0: 'exp'}, inplace=True)

    #    df['rec'] = df.index.get_level_values(0)
    #    df['c'] = 1

    #    fig, ax = plt.subplots(figsize=(4, 3))
    #    sns.violinplot(
    #        data=df,
    #        x='c',
    #        hue='rec',
    #        y='exp',
    #        palette='muted',
    #        hue_order=['Cd4', 'Cd8a', 'Cd8b1'],
    #        scale='width',
    #        )
    #    ax.legend(
    #        loc='upper left',
    #        bbox_to_anchor=(1.01, 0.99),
    #        bbox_transform=ax.transAxes,
    #        )
    #    ax.set_xticks([])
    #    ax.set_ylim(0.09, 5)
    #    ax.set_yticks([-1, 0, 1, 2, 3, 4, 5])
    #    ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
    #    ax.set_xlabel('')
    #    ax.set_ylabel('Gene expression [cpm]')

    #    fig.tight_layout()

    #if True:
    #    for ext in ['svg', 'pdf', ['png', 600]]:
    #        if isinstance(ext, list):
    #            ext, dpi = ext
    #            fig.savefig(fig_fdn+'T_cell_Cd4_8_expression.{:}'.format(
    #                ext),
    #                dpi=dpi)
    #        else:
    #            fig.savefig(fig_fdn+'T_cell_Cd4_8_expression.{:}'.format(
    #                ext))

    #if True:
    #    print('Percentage of CD4/CD8  within T cells')
    #    ni = ds0.samplesheet.groupby(['cellType', 'Timepoint', 'Treatment', 'doublet']).count().iloc[:, 0].unstack()[0].unstack()['normal'].unstack().loc['immune']

    #    df = ds.counts.loc[['Cd4']].T > 10
    #    df['Cd8'] = (ds.counts.loc[['Cd8a', 'Cd8b1']] > 10).any()
    #    df['Timepoint'] = ds.samplesheet['Timepoint']
    #    df['c'] = 1

    #    xl = ['P1', 'P7', 'P21']
    #    tab = df.groupby(['Cd4', 'Cd8', 'Timepoint']).count().iloc[:, 0].unstack().fillna(0).astype(int)
    #    tab = tab[xl]

    #    fig, ax = plt.subplots(figsize=(4, 3))
    #    x = np.arange(3)
    #    colors = sns.color_palette('muted', n_colors=4)
    #    y = np.zeros(len(x))
    #    for ir, ((cd4i, cd8i), row) in enumerate(tab.iterrows()):
    #        dy = 100.0 * row.loc[xl].values / tab.sum(axis=0)
    #        label = 'Cd4' + ('+' if cd4i else '-') + ' Cd8' + ('+' if cd8i else '-')
    #        ax.fill_between(
    #                x, y, y + dy,
    #                lw=2,
    #                color=colors[ir],
    #                alpha=0.5,
    #                label=label)
    #        y += dy
    #    ax.set_xticks(x)
    #    ax.grid(True)
    #    ax.set_ylabel('Percentage of T cells')
    #    ax.legend(
    #        loc='upper left',
    #        bbox_to_anchor=(1.01, 0.99),
    #        bbox_transform=ax.transAxes,
    #        )
    #    ax.set_xticklabels(xl)

    #    fig.tight_layout()

    #if True:
    #    for ext in ['svg', 'pdf', ['png', 600]]:
    #        if isinstance(ext, list):
    #            ext, dpi = ext
    #            fig.savefig(fig_fdn+'T_cell_Cd4_8_fractions.{:}'.format(
    #                ext),
    #                dpi=dpi)
    #        else:
    #            fig.savefig(fig_fdn+'T_cell_Cd4_8_fractions.{:}'.format(
    #                ext))

    #if True:
    #    print('Percentage of CD4/CD8  within T cells, with bars')
    #    ni = ds0.samplesheet.groupby(['cellType', 'Timepoint', 'Treatment', 'doublet']).count().iloc[:, 0].unstack()[0].unstack()['normal'].unstack().loc['immune']

    #    df = ds.counts.loc[['Cd4']].T > 10
    #    df['Cd8'] = (ds.counts.loc[['Cd8a', 'Cd8b1']] > 50).any()
    #    df['Timepoint'] = ds.samplesheet['Timepoint']
    #    df['c'] = 1

    #    xl = ['P1', 'P7', 'P21']
    #    tab = df.groupby(['Cd4', 'Cd8', 'Timepoint']).count().iloc[:, 0].unstack().fillna(0).astype(int)
    #    tab = tab[xl]

    #    fig, ax = plt.subplots(figsize=(4, 3))
    #    x = np.arange(3)
    #    colors = sns.color_palette('muted', n_colors=4)
    #    y = np.zeros(len(x))
    #    for ir, ((cd4i, cd8i), row) in enumerate(tab.iterrows()):
    #        dy = 100.0 * row.loc[xl].values / tab.sum(axis=0)
    #        label = 'Cd4' + ('+' if cd4i else '-') + ' Cd8' + ('+' if cd8i else '-')
    #        ax.bar(
    #                x + 0.1, dy, width=0.8, bottom=y,
    #                lw=2,
    #                color=colors[ir],
    #                alpha=0.8,

    #                label=label)
    #        y += dy
    #    ax.set_xticks(x)
    #    ax.grid(True)
    #    ax.set_ylabel('Percentage of T cells')
    #    ax.legend(
    #        loc='upper left',
    #        bbox_to_anchor=(1.01, 0.99),
    #        bbox_transform=ax.transAxes,
    #        )
    #    ax.set_xticklabels(xl)

    #    fig.tight_layout()

    #if True:
    #    for ext in ['svg', 'pdf', ['png', 600]]:
    #        if isinstance(ext, list):
    #            ext, dpi = ext
    #            fig.savefig(fig_fdn+'T_cell_Cd4_8_fractions_bars.{:}'.format(
    #                ext),
    #                dpi=dpi)
    #        else:
    #            fig.savefig(fig_fdn+'T_cell_Cd4_8_fractions_bars.{:}'.format(
    #                ext))

    if True:
        print('Plot tSNEs of T cells')
        vs = pd.read_csv(
            '../../data/sequencing/datasets/all_{:}/tsne_immune.tsv'.format(versions[-2]),
            sep='\t',
            index_col=0,
            )
        vs = vs.loc[ds.samplenames]

        genes = ['Trac', 'Tcrg-C4', 'Tcrg-C2', 'Tcrg-C1']
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
                        fig.savefig(fig_fdn+'immune_tsne_Tcells_{:}.{:}'.format(
                            gene, ext),
                            dpi=dpi)
                    else:
                        fig.savefig(fig_fdn+'immune_tsne_Tcells_{:}.{:}'.format(
                            gene, ext))

    plt.ion(); plt.show(); sys.exit()

    print('Look at B cells')
    dsB = ds0.query_samples_by_metadata(
        'cellSubtype == "B cell"',
        )

    print('Load tSNE from file')
    vsB = pd.read_csv(
        '../../data/sequencing/datasets/all_{:}/tsne_immune.tsv'.format(versions[-2]),
        sep='\t',
        index_col=0,
        )
    vsB = vsB.loc[dsB.samplenames]

    if True:
        print('Plot B cell maturation genes')

        genes = [
            'Ighm',
            'Ighd',
            'Mki67',
            'Tbx21',
            'Aicda',
            'Prdm1',
            'Ebi3',
            'coverage',
            ]

        fig, axs = plt.subplots(2, 4, figsize=(8.2, 4), sharex=True, sharey=True)
        axs = axs.ravel()
        for ir, (gene, ax) in enumerate(zip(genes, axs)):
            dsB.plot.scatter_reduced_samples(
                    vsB,
                    ax=ax,
                    s=12,
                    alpha=0.50,
                    cmap='viridis',
                    color_by=gene,
                    color_log=True,
                    #high_on_top=True,
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
                fig.savefig(fig_fdn+'immune_tsne_Bcell_maturation.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'immune_tsne_Bcell_maturation.{:}'.format(
                    ext))

    if True:
        print('Study time points of proliferative cells')
        df = dsB.counts.loc[['Mki67']].T > 0
        df['Timepoint'] = dsB.samplesheet['Timepoint']
        xorder = ['E18.5', 'P1', 'P7', 'P21']
        frac = 100 * df.groupby('Timepoint').mean().loc[xorder, 'Mki67']

        from scipy import interpolate
        fig, ax = plt.subplots(figsize=(4, 2.5))
        color = 'darkslateblue'
        x = np.arange(4)
        y = np.maximum(0.5, frac.loc[xorder])

        outx = np.linspace(0, 3, 100)
        outy = interpolate.pchip_interpolate(x, y, outx)
        out = np.vstack([outx, outy])

        ax.scatter(
                x, y,
                marker='o',
                lw=2, alpha=0.8,
                edgecolor=color,
                facecolor='none',
                zorder=10,
                )
        ax.plot(out[0], out[1], lw=2,
                alpha=0.4,
                color=color,
                zorder=10,
                )
        ax.grid(True)
        ax.set_xticks(x)
        ax.set_xticklabels(xorder)
        ax.set_ylabel('Mki67+ (% of B cells)')
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax.set_ylim(0, 100)

        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'proliferating_population_abundance_Bcells.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'proliferating_population_abundance_Bcells.{:}'.format(
                    ext))



    if True:
        print('Analyze BCRs and hypermutation')
        fn = '../../data/sequencing/datasets/all_20190828/filtered_BCR_summary/IMGT_gapped_db.tab'
        bcr = pd.read_csv(fn, sep='\t', index_col=0)
        bcrh = bcr.loc[(bcr['LOCUS'] == 'H') & (bcr['IN_FRAME'] == 'T')]
        bcrl = bcr.loc[bcr['LOCUS'].isin(['K', 'J']) & (bcr['IN_FRAME'] == 'T')]
        vmuth = bcrh['V_IDENTITY']
        vmutl = bcrl['V_IDENTITY']
        jmuth = bcrh['J_IDENTITY']
        jmutl = bcrl['J_IDENTITY']

        fig, ax = plt.subplots(figsize=(3, 2.5))
        pdata = [
            {'chain': 'Heavy', 'locus': 'V', 'datum': vmuth},
            {'chain': 'Light', 'locus': 'V', 'datum': vmutl},
            {'chain': 'Heavy', 'locus': 'J', 'datum': jmuth},
            {'chain': 'Light', 'locus': 'J', 'datum': jmutl},
            ]
        for datum in pdata:
            vm = datum['datum']
            title = datum['chain'] + ' ' + datum['locus']
            y = np.sort(vm.values)
            x = np.linspace(0, 1, len(y))
            ax.plot(x, y, lw=2, label=title)

        ax.set_yticks([0.9, 0.95, 1.0])
        ax.set_xticks([0.0, 0.25, 0.50, 0.75])
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(0.85, 1.005)
        ax.legend(loc='lower right')
        ax.set_ylabel('V/J gene identity')
        ax.set_xlabel('Cumulative cell fraction')
        ax.grid(True)

        fig.tight_layout()

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'Bcell_mutdistance.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'Bcell_mutdistance.{:}'.format(
                    ext))

    if True:
        print('Analyze BCRs and hypermutation by time')
        fn = '../../data/sequencing/datasets/all_20190828/filtered_BCR_summary/IMGT_gapped_db.tab'
        bcr = pd.read_csv(fn, sep='\t', index_col=0)

        tmp = dsB.samplesheet['fastq'].str.replace('/', '_')
        bcells_all = pd.Series(data=tmp.index, index=tmp.values)

        pdata = {}
        chains = [('Heavy', ['H']), ('Light', ['L', 'K'])]
        for chain, chs in chains:
            bcri = bcr.loc[bcr['LOCUS'].isin(chs) & (bcr['IN_FRAME'] == 'T')]

            bcri_names = bcri.index.str.slice(7)
            tp = dsB.samplesheet.loc[bcells_all.loc[bcri_names].values, 'Timepoint']

            for locus in ['V', 'J']:
                for t in ['E18.5', 'P1', 'P7', 'P21']:
                    ind = (tp == t).values.nonzero()[0]
                    mut = bcri.iloc[ind]['{:}_IDENTITY'.format(locus)]
                    pdata[(chain, locus, t)] = mut.values

        fig, axs = plt.subplots(2, 2, figsize=(5.2, 4.5), sharex=True, sharey=True)
        cmap = {
            'E18.5': 'navy',
            'P1': 'gold',
            'P7': 'tomato',
            'P21': 'firebrick',
            }
        for chain, ax_row in zip(chains, axs):
            chname = chain[0]
            for locus, ax in zip(['V', 'J'], ax_row):
                title = chname + ' ' + locus
                for t in ['E18.5', 'P1', 'P7', 'P21']:
                    y = np.sort(pdata[(chname, locus, t)])
                    x = np.linspace(0, 1, len(y))
                    ax.plot(x, y, lw=2, color=cmap[t], label=t)
                ax.set_yticks([0.9, 0.95, 1.0])
                ax.set_xticks([0.0, 0.1, 0.25, 0.50, 0.75])
                ax.set_xlim(-0.01, 1.01)
                ax.set_ylim(0.85, 1.005)
                ax.grid(True)
            axs[0][0].legend(loc='lower right')

            fig.text(0.02, 0.5, 'V/J gene identity', rotation=90, ha='center', va='center')
            fig.text(0.52, 0.02, 'Cumulative cell fraction', ha='center', va='center')

        ax.set_ylim(0.9, 1.005)
        ax.set_xlim(-0.01, 0.2)
        fig.tight_layout(rect=(0.02, 0.02, 1, 1))

    if True:
        for ext in ['svg', 'pdf', ['png', 600]]:
            if isinstance(ext, list):
                ext, dpi = ext
                fig.savefig(fig_fdn+'Bcell_mutdistance_bytime.{:}'.format(
                    ext),
                    dpi=dpi)
            else:
                fig.savefig(fig_fdn+'Bcell_mutdistance_bytime.{:}'.format(
                    ext))

    plt.ion()
    plt.show()
