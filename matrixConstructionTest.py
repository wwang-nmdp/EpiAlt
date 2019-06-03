#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur May 12 2019

@author: wwang

For testing the code of matrixConstruction.
"""

import pandas as pd
import numpy as np


mode = 'DR_genotype_encoding'
meta_fp = '~/GWASH/unimputed/cleaned_BMT_cases/GWAS_cleaned_BMT_cases_info.csv'
encodedMat_fp = '~/GWASH/unimputed/cleaned_BMT_cases/'+ mode + '/'

meta_info = pd.read_csv(meta_fp, index_col = 'bmt_case')

Encoded_mat = pd.read_hdf(encodedMat_fp + 'EncodedMatrix_95filtered_' + mode + '.h5')

agvhd24 = meta_info.loc[Encoded_mat.index.astype(int),'agvhi24']
agvhd34 = meta_info.loc[Encoded_mat.index.astype(int),'agvhi34']

disease = meta_info.loc[Encoded_mat.index.astype(int),'disease']

#######
# Donors Only
#######
donors_only_mat = (Encoded_mat/3).apply(np.floor)
donors_only_mat.insert(0, 'agvhi24', agvhd24.values)
donors_only_mat.insert(0, 'agvhi34', agvhd34.values)
donors_only_mat.insert(0, 'disease', disease.values)

donors_only_mat = donors_only_mat[donors_only_mat['agvhi24'].notna()]
# Encoded_mat34 = Encoded_mat[Encoded_mat['agvhi34'].notna()]
donors_only_mat = donors_only_mat.fillna(-9)
donors_only_mat = donors_only_mat.astype(int)

feather.write_dataframe(donors_only_mat, 'donors_only_cleaned_SNP_mat_outcome_disease.feather')

#######
# Patients Only
#######

patients_only_mat = Encoded_mat.mod(3)

patients_only_mat.insert(0, 'agvhi24', agvhd24.values)
patients_only_mat.insert(0, 'agvhi34', agvhd34.values)
patients_only_mat.insert(0, 'disease', disease.values)

patients_only_mat = patients_only_mat[patients_only_mat['agvhi24'].notna()]
# Encoded_mat34 = Encoded_mat[Encoded_mat['agvhi34'].notna()]
patients_only_mat = patients_only_mat.fillna(-9)
patients_only_mat = patients_only_mat.astype(int)

feather.write_dataframe(patients_only_mat, 'patient_only_cleaned_SNP_mat_outcome_disease.feather')

########
# Autosome only
########

meta_fp = '/mnt/cloudbiodata_nfs_2/users/hhuang/GWASH/unimputed/cleaned_BMT_cases/GWAS_cleaned_BMT_cases_info.csv'
encodedMat_fp = '/mnt/cloudbiodata_nfs_2/users/hhuang/GWASH/unimputed/cleaned_BMT_cases/'+ mode + '/'
autosome_avail_snps_fp = '/mnt/cloudbiodata_nfs_2/users/hhuang/GWASH/unimputed/cleaned_BMT_cases/patient_only/autosome_avail_snp_map.csv'

#meta_fp = '~/data/GWAS/unimputed/GWAS_cleaned_BMT_cases_info.csv'
#encodedMat_fp = '/home/hhuang/data/GWAS/unimputed/cleaned_BMT_cases/MismatchEncoded/' +mode+'/'
#autosome_avail_snps_fp = '/home/hhuang/data/GWAS/unimputed/cleaned_BMT_cases/patients_only/autosome_avail_snp_map.csv'

Encoded_mat = pd.read_hdf(encodedMat_fp + 'EncodedMatrix_95filtered_' + mode + '.h5')

agvhd24 = meta_info.loc[Encoded_mat.index.astype(int),'agvhi24']
agvhd34 = meta_info.loc[Encoded_mat.index.astype(int),'agvhi34']

disease = meta_info.loc[Encoded_mat.index.astype(int),'disease']

autosome_avail_snps = pd.read_csv(autosome_avail_snps_fp)
autosome_encoded = Encoded_mat[autosome_avail_snps['RS']]
autosome_encoded['BMT_case'] = autosome_encoded.index

autosome_encoded_index = pd.DataFrame(autosome_encoded.index.values, columns=['BMT_case'])
feather.write_dataframe(autosome_encoded_index, 'DR_genotype_encoding_autsome_only_index.feather')
feather.write_dataframe(autosome_encoded, 'DR_genotype_encoding_autsome_only.feather')

del Encoded_mat
gc.collect()

########
# AML case-control genotypes in one matrix
########

autosome_encoded = feather.read_dataframe('DR_genotype_encoding_autsome_only.feather')
autosome_encoded = autosome_encoded.set_index('BMT_case')
gc.collect()

disease = meta_info.loc[autosome_encoded.index.astype(int),'disease']

donors_only_mat = (autosome_encoded/3).apply(np.floor)
donors_only_mat.insert(0, 'disease', disease.values)

patients_only_mat = autosome_encoded.mod(3)
patients_only_mat.insert(0, 'disease', disease.values)

del autosome_encoded
gc.collect()

AML_donor = donors_only_mat[donors_only_mat['disease'] == 10]
del donors_only_mat
gc.collect()

AML_patient = patients_only_mat[patients_only_mat['disease'] == 10]
del patients_only_mat
gc.collect()

# AML case - control
AML_donor['disease'] = 0
AML_patient['disease'] = 1

genetic_data = pd.concat([AML_donor, AML_patient], axis=0)

del AML_donor, AML_patient
gc.collect()

BMT_cases=pd.DataFrame(data=genetic_data.index.values,columns=['BMT_case'])

feather.write_dataframe(genetic_data, 'AML_case_control_autosome_OneMatrix.feather')
feather.write_dataframe(BMT_cases, 'AML_case_control_autosome_OneMatrix_BMT_cases.feather')

#######
# acute GVHD:  donor+ AML recipient in one matrix
#######


patients_only_mat = feather.read_dataframe('patient_only_autosome_SNP_mat_outcome_disease.feather')
patients_only_mat = patients_only_mat.set_index('BMT_case')
gc.collect()

AML_patient = patients_only_mat[patients_only_mat['disease'] == 10]
del patients_only_mat
gc.collect()

donors_only_mat = feather.read_dataframe('donors_only_autosome_SNP_mat_outcome_disease.feather')
donors_only_mat = donors_only_mat.set_index('BMT_case')
gc.collect()

AML_donor = donors_only_mat[donors_only_mat['disease'] == 10]
del donors_only_mat
gc.collect()

agvhd24 = AML_patient['agvhi24']
agvhd34 = AML_patient['agvhi34']

AML_patient_mat = AML_patient.drop(['disease', 'agvhi24', 'agvhi34'], axis=1)
AML_patient_mat = AML_patient_mat.add_prefix('r_')
del AML_patient
gc.collect()

AML_donor_mat = AML_donor.drop(['disease', 'agvhi24', 'agvhi34'], axis=1)
AML_donor_mat = AML_donor_mat.add_prefix('d_')
del AML_donor
gc.collect()

patient_donor_mat = pd.concat([AML_patient_mat, AML_donor_mat], axis=1, sort=False)
del AML_patient_mat
del AML_donor_mat
gc.collect()

patient_donor_mat = patient_donor_mat.astype(int)

patient_donor_mat['BMT_case'] = patient_donor_mat.index

feather.write_dataframe(patient_donor_mat, 'AML_patient_donor_autosome_SNP_mat.feather')

patient_donor_mat.insert(0, 'agvhi34', agvhd34.values)
patient_donor_mat.insert(0, 'agvhi24', agvhd24.values)
feather.write_dataframe(patient_donor_mat, 'AML_patient_donor_autosome_SNP_mat_TXoutcome.feather')

