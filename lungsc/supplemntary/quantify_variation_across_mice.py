# vim: fdm=indent
'''
author:     Fabio Zanini
date:       17/01/20
content:    Quantify variation across mice.
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


if __name__ == '__main__':

    ds0 = DatasetLung.load(preprocess=True, version=versions[-2])
    ds0.query_samples_by_metadata(
        '(cellType == "immune") & (doublet == 0)', inplace=True)

    print('Load tSNE from file')
    vs = pd.read_csv(
        '../../data/sequencing/datasets/all_{:}/tsne_immune.tsv'.format(versions[-2]),
        sep='\t',
        index_col=0,
        )

    csts = ds0.samplesheet['cellSubtype'].unique()

    data_plot = {}
    for cst in csts:
        dsi = ds0.query_samples_by_metadata('cellSubtype == @cst', local_dict=locals())

        vsi = vs.loc[dsi.samplenames]

        print('Aggregate info for cluster {:}'.format(cst))
        df = vsi.copy()
        for col in ['Timepoint', 'Gender', 'Mousename']:
            df[col] = dsi.samplesheet[col]

        from scipy.spatial.distance import cdist
        for tp, dfi in df.groupby('Timepoint'):
            mns = dfi['Mousename'].unique()
            if len(mns) < 2:
                continue

            # First mouse is always female
            if mns[1].startswith('F'):
                mns = mns[::-1]
            dfid = {mn: dfi.loc[dfi['Mousename'] == mn] for mn in mns}

            m1 = dfid[mns[0]][['dimension 1', 'dimension 2']].values
            m2 = dfid[mns[1]][['dimension 1', 'dimension 2']].values

            if (len(m1) < 20) or (len(m2) < 20):
                continue

            dm11 = cdist(m1, m1)
            dm22 = cdist(m2, m2)
            dm12 = cdist(m1, m2)

            # Get random pairs
            n = 100
            dself1 = dm11[
                    np.random.randint(len(m1), size=n),
                    np.random.randint(len(m1), size=n),
                    ]
            dself2 = dm22[
                    np.random.randint(len(m2), size=n),
                    np.random.randint(len(m2), size=n),
                    ]
            dcross = dm12[
                    np.random.randint(len(m1), size=n),
                    np.random.randint(len(m2), size=n),
                    ]
            data_plot[(cst, tp, '11')] = dself1
            data_plot[(cst, tp, '22')] = dself2
            data_plot[(cst, tp, '12')] = dcross

    fig, axs = plt.subplots(4, 4, figsize=(12, 10), sharex=True, sharey=True)
    axs = axs.ravel()
    iax = 0
    for cst in csts:
        for tp in ['E18.5', 'P1', 'P7', 'P21']:
            if (cst, tp, '11')not in data_plot:
                continue
            dself1 = data_plot[(cst, tp, '11')]
            dself2 = data_plot[(cst, tp, '22')]
            dcross = data_plot[(cst, tp, '12')]

            ax = axs[iax]
            iax += 1
            ax.set_title('{:}, {:}'.format(cst, tp))
            x11 = np.sort(dself1)
            x22 = np.sort(dself2)
            x12 = np.sort(dcross)
            y11 = 1.0 - np.linspace(0, 1, len(x11))
            y22 = 1.0 - np.linspace(0, 1, len(x22))
            y12 = 1.0 - np.linspace(0, 1, len(x12))

            ax.plot(x11, y11, lw=2, color='steelblue', label='W/in female')
            ax.plot(x22, y22, lw=2, color='olive', label='W/in male')
            ax.plot(x12, y12, lw=2, color='darkred', label='Between')

            ax.grid(True)
            ax.set_ylim(-0.01, 1.01)

    axs[0].legend(loc='upper right')

    fig.text(0.52, 0.02, 'Euclidean distance in t-SNE', ha='center')
    fig.text(0.02, 0.52, 'Fraction of pairs < distance x', va='center', rotation=90)
    fig.tight_layout(rect=(0.04, 0.04, 1, 0.96))
    fig.savefig('../../figures/immune_crossmouse_variation.png')

    plt.ion()
    plt.show()
