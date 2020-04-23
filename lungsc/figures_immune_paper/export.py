# vim: fdm=indent
'''
author:     Fabio Zanini
date:       21/09/19
content:    Export samplesheet for immune paper.
'''
import os
import sys
import glob
import gzip
import subprocess as sp
import numpy as np
import pandas as pd

from lungsc.pilots.load_dataset import DatasetLung, versions


if __name__ == '__main__':

    ds = DatasetLung.load(preprocess=True, version=versions[-1])
    ds.query_samples_by_metadata('(cellType == "immune") & (doublet == 0) & (Treatment == "normal")', inplace=True)

    print('Export sampleseet')
    cols = [
            'Gender',
            'Timepoint',
            'Mousename',
            'DC',
            'Well',
            'fastq',
            'cellType',
            'cellSubtype',
        ]
    fn_ss = '../../data/sequencing/datasets/all_{:}/samplesheet_immune_export.tsv'.format(versions[-1])
    ds.samplesheet[cols].to_csv(
        fn_ss,
        sep='\t',
        index=True,
        header=True,
        )

    print('Copy samplesheet to remote server')
    sp.run(
        'scp {:} sk2.0-1:{:}'.format(
            fn_ss,
            '/oak/stanford/groups/quake/fzanini/sequencing_data/lung_development/export/immune/',
            ),
        shell=True,
        check=True,
        )
