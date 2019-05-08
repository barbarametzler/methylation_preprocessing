#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 16:46:57 2019

@author: nicolasagrotis
"""
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyreadr
from math import sqrt
import timeit
#import quadprog
#from CH3.python.preprocess import preprocess 


## load dependencies
import pandas as pd
import numpy as np
import sys
import os
import pyreadr


from CH3.python.quality_control_1 import visualisation_plots, remove_unreliable_samples, k_mean_sex_infer, infer_sex, call_snps, identify_replicates, compare_sex


'''
control_beads = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/hm450_controls.Rds')
control_beads = control_beads[None]


# In functions its dnam
cpgs = pd.read_csv("/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData2/cpgs.csv",index_col='Unnamed: 0')
#1min 25s ± 3.06 s per loop (mean ± std. dev. of 7 runs, 1 loop each)



#had to be downloaded from R
coefs = pd.read_csv("/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/coefs.csv",index_col='Unnamed: 0')

controls = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/eira_controls.Rds')
controls = controls[None]

#snps = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData2/eira_snps.Rds')
#snps = snps[None]

snps = pd.read_csv("/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData2/snps.csv",index_col='Unnamed: 0')

samples = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/eira_samples.Rds')
samples = samples[None]
samples.set_index('sample.id',inplace=True)

# Load individual characteristics data

covars = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/Covariates.Rds')
covars = covars[None]
covars.set_index('gsm',inplace=True)

samples_sheet = pyreadr.read_r('/Users/nicolasagrotis/Desktop/illuminAlysis/illumiData/Sample_sheet.Rds')
samples_sheet = samples_sheet[None]
samples_sheet.set_index('sample.name',inplace=True)

common= samples_sheet.index.intersection(covars.index)

covars=covars.loc[common]
covars['sample_id']=samples_sheet['sample.id']
covars.set_index('sample_id',inplace=True)


'''
##### added Barbs


## load dependencies
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
import timeit
import seaborn as sns
import matplotlib.pyplot as plt
from CH3.python.illuminaio import list_idat
from CH3.python.preprocess import preprocess
from CH3.python.quality_control_1 import visualisation_plots


probes_file = 'CH3/python/testing/probes_1.csv'
controls_file = 'CH3/python/testing/control_beads.csv'
idat_files_folder = 'idat'


#control_beads = spyreadr.read_r('/Users/metzlerabarbara/Documents/GitHub/methylation_preprocessing/Python_package/CH3/illumina_manifests/hm450_controls.rds')
#control_beads = control_beads[None]

# load covars
covars = pyreadr.read_r('Covariates.Rds')
covars = covars[None]
covars.set_index('gsm',inplace=True)


# load sample sheet
samples_sheet = pyreadr.read_r('Sample_sheet.Rds')
samples_sheet = samples_sheet[None]



##Preprocessing

#for id, df in enumerate(data_list):
samples, cpgs, snps, intensities_A, intensities_B, controls_red, controls_grn = preprocess(probes_file, controls_file,
    idat_files_folder, min_beads=3, detection=0.05, return_intensities=True)

#cpgs = pd.read_csv("/Users/metzlerabarbara/Downloads/cpgs.csv",index_col='Unnamed: 0', engine='c')


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Visualisation of data 
# i is the ith sample that will be viewed for the last graph showing the signature of the methylation beta values
# all graphs functioning correctly

visualisation_plots(snps,3,samples,cpgs)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Remove variables based on percentage of missing variables and quality check of 
# bc1.grn,bc1.red and bc2
# functioning correctly


samples,cpgs,covars=remove_unreliable_samples(samples,0.1,cpgs,covars)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Apply a K-means classification model to investigate the clustering of males and 
# females based on the missing Y chromosome and the median X chromosome

#k_mean_sex_infer(samples)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Use the print below in the workflow document
# plt.scatter(x=samples['median_chrX'],y=samples['missing_chrY'])
samples=infer_sex(samples,threshold_chrX=0.37,threshold_chrY=0.39)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# snps distribution plot for ith snp
# The allele at a SNP locus can be inferred from SNP intensities measured on the BeadChip.
# Chose the sample number for which the SNPS will be visualised

call_snps(snps,4)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Need samples['sex'] for this function
# Need the infer_sex function
identify_replicates(snps,42,samples)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Compare the infered sex with the true sex available throught the covars dataset
# NEED TO REVERSE THE PRINTING RECCOMMENDATION
samples=compare_sex(covars,samples)






