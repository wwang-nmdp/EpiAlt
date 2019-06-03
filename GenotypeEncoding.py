#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur May 12 2019

@author: wwang

For testing the code of genotypeExtraction.
"""
import sys

# sys.path.append('/home/hhuang/WGS2018/utils/')
from utils import coreFunctions as cf
import numpy as np
import pandas as pd

from pandas_plink import read_plink

import gc

import pickle

mode = 'DR_genotype_encoding'
## EC2
plink_fp = "~/efs/GWASH_IMPUTED_DATA/ImputeQC/"
metadata_fp = "~/data/GWAS/available_cases.csv"
output_fp = '~/data/GWAS/MismatchEncoded/' + mode + '/'

chrList = range(1, 23)  # chr 1:22

metadata_avail_cases = cf.readCaseInfo(metadata_fp)

num_case = len(metadata_avail_cases['BMTcase'])

metadata_pd = pd.read_csv(metadata_fp, index_col=2)
metadata_pd.index = metadata_pd.index.map(str)

for chrom in chrList:

    try:
        (bim, fam, bed) = read_plink(plink_fp + str(chrom) + 'bed')
        # bim - pandas.DataFrame – Alleles.
        # fam - pandas.DataFrame – Samples.
        # G - numpy.ndarray – Genotype.

        DR_gt_encoding = np.empty([num_case, bim.shape[0]], dtype='float64') - 1

        # fam.iid.head()
        # G.head()

        for case_index in range(num_case):
            bmt_fams_d = fam.query("iid in ['" + metadata_avail_cases.iloc[case_index]['did'] + "']")

            bmt_fams_r = fam.query("iid in ['" + metadata_avail_cases.iloc[case_index]['rid'] + "']")

            gt_d = bed[:, bmt_fams_d.i.values].compute()
            gt_r = bed[:, bmt_fams_r.i.values].compute()

            bmt_gt = gt_d * 3 + gt_r  # donor * 3 + recipient = gt code

            DR_gt_encoding[case_index, :] = bmt_gt

        # Remove all duplicated columns (SNPs)
        if len(bim['snp'][bim['snp'].duplicated(False)]) > 0:
            print('>>>> Column(s) {0} have duplicates! Dropping all duplicated columns (SNPs)'.format(list(set(bim['snp'][bim['snp'].duplicated(False)]))))

        BMT_pdMtx = pd.DataFrame(data=DR_gt_encoding[:, ~(bim['snp'].duplicated(False))],
                                 index=metadata_avail_cases['BMTcase'],
                                 columns=bim['snp'][~(bim['snp'].duplicated(False))])

        BMT_pdMtx.to_hdf(output_fp+'EncodedMatrix/chr'+str(chrom)+'_EncodedMatrix_original_' + mode + '.h5',
                         key='chr_'+str(chrom), complib='blosc', complevel=9)

        # filter 95% call rate
        filtered_index = (1 - BMT_pdMtx.isnull().sum()/BMT_pdMtx.shape[0]) >= 0.95
        BMT_pdMtx_filtered = BMT_pdMtx.loc[:, filtered_index]
        print('Chromosome {0} removed {1} ({2:.2f} %) out of {3} variants (95 % call rate)'.format(chrom,
                                                                                                   (filtered_index == False).sum(),
                                                                                                   (filtered_index == False).sum()/filtered_index.count()*100,
                                                                                                   filtered_index.count()))
        BMT_pdMtx_filtered.to_hdf(output_fp + 'EncodedMatrix/chr' + str(chrom) + '_EncodedMatrix_95filtered_' + mode + '.h5',
                           key='chr_' + str(chrom), complib='blosc', complevel=9)


        gc.collect()

    except FileNotFoundError as e:
        print(e)

