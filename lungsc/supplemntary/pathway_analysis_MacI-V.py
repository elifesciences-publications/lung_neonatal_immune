# vim: fdm=indent
'''
author:     Fabio Zanini
date:       29/04/20
content:    Pathway analysis for marker genes for Mac I-V.
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



if __name__ == '__main__':

    go_long = pd.read_csv(
            '../../data/GO_annotations/gene_association.mgi.gz',
            compression='gzip', sep='\t', comment='!',
            header=None,
            )

    go_terms = pd.read_csv(
            '../../data/GO_annotations/go_terms.mgi',
            sep='\t',
            header=None,
            ).set_index(1)
    go_terms = go_terms.loc[go_terms[0] == 'Biological Process'][2]

    def fun(df):
        gt = go_terms.loc[df[4]]
        idx = gt.str.contains('cell cycle')
        #idx |= ...
        return idx.sum() == 0

    gby = go_long.loc[go_long[4].isin(go_terms.index)].groupby(2)
    go_genes = {key for key, val in gby if fun(val)}
    go_genes = sorted(go_genes)

    ds = DatasetLung.load(preprocess=True, version=versions[-2])
    ds.query_samples_by_metadata(
            '(cellType == "immune") & (doublet == 0) & (Treatment == "normal") & (cellSubtype in ("Mac I", "Mac II", "Mac III", "Mac IV", "Mac V"))',
            inplace=True)


    print('Find markers')
    csts = ['Mac I', 'Mac II', 'Mac III', 'Mac IV', 'Mac V']
    comps = {}
    for cst in csts:
        print(cst)
        ds.samplesheet['is_{:}'.format(cst)] = ds.samplesheet['cellSubtype'] == cst
        dsp = ds.split('is_{:}'.format(cst))
        dsp[True].subsample(100, inplace=True)
        dsp[False].subsample(100, inplace=True)
        comp = dsp[True].compare(dsp[False])
        comp['log2_fold_change_abs'] = np.abs(comp['log2_fold_change'])
        comps[cst] = comp
        break

    sys.exit()

    for cst, comp in comps.items():
        idx = comp['statistic'] >= 0.30
        idx &= comp['log2_fold_change'] > 0
        idx &= comp.index.isin(go_genes)
        degs = comp.loc[idx].nlargest(500, 'log2_fold_change')
        print(','.join(degs.index))
