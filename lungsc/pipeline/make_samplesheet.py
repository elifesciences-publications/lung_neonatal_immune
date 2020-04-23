# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/19
content:    Make samplesheet from the various sequencing runs
'''
import os
import glob
import numpy as np
import pandas as pd


root_fdn = '/oak/stanford/groups/quake/fzanini/sequencing_data/lung_development/'

# Each DC is a sequencing run with various numbers of plates
# DC29 is a repeat and has swapped samplenames, let's not touch it
dcs = [
    'test',
    'DC20', 'DC21', 'DC23', 'DC25',
    'DC28', 'DC30', 'DC33', 'DC34',
    'DC35', 'DC36',
    ]


def make_samplename_from_fastq(dc, fdn):
    if 'test' in fdn:
        return 'test'

    p2 = fdn.split('_')[1]

    gender = p2[0]

    if 'yperoxia' in fdn:
        treatment = 'hyperoxia'
        timepoint = 'P7'
        p2_rest = p2[1+len('hyperoxia'):]
    else:
        treatment = 'normal'
        if 'E18' in p2:
            timepoint = 'E18.5'
            p2_rest = p2[4:]
        elif ('P1' in p2) or ('P7' in p2):
            timepoint = p2[1:3]
            p2_rest = p2[3:]
        elif 'P21' in p2:
            timepoint = 'P21'
            p2_rest = p2[4:]
        else:
            raise ValueError('Time point not found in foldername: {:}'.format(fdn))

    if 'mesench1' in p2_rest:
        sort = 'mesench'
        plate = '1'
    elif 'mesench2' in p2_rest:
        sort = 'mesench'
        plate = '2'
    elif 'CD452' in p2_rest:
        sort = 'CD45'
        plate = '2'
    elif 'CD45' in p2_rest:
        sort = 'CD45'
        plate = '1'
    elif 'CD31' in p2_rest:
        sort = 'CD31'
        plate = '1'
    else:
        raise ValueError(
            'Sort type not found in foldername: {:}: p2_rest = {:}'.format(
                fdn, p2_rest))

    well = fdn.split('_')[-1]

    sn = '_'.join([gender, timepoint, treatment, sort, plate, dc, well])

    return sn


if __name__ == '__main__':

    # Start from fastqs
    fastqs = []
    for x in dcs:
        xx = [(x, y, x + '/' + y) for y in os.listdir(root_fdn+'fastq/'+x) if 'merge_lanes' not in y]
        fastqs.extend(xx)
    sampletable = pd.DataFrame(
        fastqs, columns=['DC', 'samplename_fastq', 'fastq'],
        )
    sampletable['bam'] = ''
    sampletable['samplename_bam'] = ''

    # Add samplename
    sampletable['samplename'] = ''
    for idx, row in sampletable.iterrows():
        dc = row['DC']
        fdn = row['samplename_fastq']
        sampletable.at[idx, 'samplename'] = make_samplename_from_fastq(dc, fdn)
    sampletable.set_index('samplename', drop=True, inplace=True)

    # Integrate bam files
    for dc in dcs:
        print('DC: {:}'.format(dc))
        subtable = sampletable.loc[sampletable['DC'] == dc]
        if subtable.shape[0] == 0:
            print('MISSING FASTQ FILES!')
            continue

        bam_fdns = os.listdir(root_fdn+'bamfiles/{:}'.format(dc))

        # Add the bamfile foldernames
        if not len(bam_fdns):
            print('{:}: missing BAM folders'.format(dc))
            for idx, row in subtable.iterrows():
                fqname = row['samplename_fastq']
                sampletable.at[idx, 'samplename_bam'] = fqname
                sampletable.at[idx, 'bam'] = dc+'/'+fqname
            continue

        # We only get here for DCs that have both a fastq and a BAM
        # These are DC34 and DC35 only (!)
        # 1. Look for a fastq folder with the same name as the BAM folder
        snlist = list(subtable['samplename_fastq'])
        found = {key: False for key in snlist}
        for fdn in bam_fdns:
            if fdn in snlist:
                idx = subtable.index[snlist.index(fdn)]
                sampletable.at[idx, 'samplename_bam'] = fdn
                sampletable.at[idx, 'bam'] = dc+'/'+fdn
                found[idx] = True
                continue

            print("Fastq filename not found for BAM folder: {:}".format(fdn))
            break

    idx = sampletable.index
    sampletable['Gender'] = ''
    sampletable['Timepoint'] = ''
    sampletable['Treatment'] = ''
    sampletable['SortType'] = ''
    sampletable['Plate'] = ''
    sampletable['Well'] = ''
    mat = data = np.array(
        list(sampletable.iloc[1:].index.str.split('_', expand=False).values))
    sampletable.loc[idx[1:], 'Gender'] = mat[:, 0]
    sampletable.loc[idx[1:], 'Timepoint'] = mat[:, 1]
    sampletable.loc[idx[1:], 'Treatment'] = mat[:, 2]
    sampletable.loc[idx[1:], 'SortType'] = mat[:, 3]
    sampletable.loc[idx[1:], 'Plate'] = mat[:, 4]
    sampletable.loc[idx[1:], 'Well'] = mat[:, 6]

    sampletable.to_csv('../../data/samplesheet.tsv', sep='\t', index=True)
