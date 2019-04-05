## load dependencies
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
import timeit
from CH3.python.illuminaio import list_idat
from CH3.python.preprocess import preprocess
from CH3.python.quality_control_1 import snps_distribution_box, remove_unreliable_samples, k_mean_sex_infer, infer_sex, snps_distribution, identify_replicates, compare_sex, estimate_leukocytes

start_time = timeit.default_timer()

'''
##load data for testing
beads5 = '/Users/metzlerabarbara/Library/Mobile Documents/com~apple~CloudDocs/dnam/R05C01_beads.csv'
beads4 = '/Users/metzlerabarbara/Library/Mobile Documents/com~apple~CloudDocs/dnam/R04C01_beads.csv'
beads3 = '/Users/metzlerabarbara/Library/Mobile Documents/com~apple~CloudDocs/dnam/R03C01_beads.csv'
beads2 = '/Users/metzlerabarbara/Library/Mobile Documents/com~apple~CloudDocs/dnam/R02C01_beads.csv'
beads1 = '/Users/metzlerabarbara/Library/Mobile Documents/com~apple~CloudDocs/dnam/R01C01_beads.csv'

data1 = load_data(beads1)
data2 = load_data(beads2)
data3 = load_data(beads3)
data4 = load_data(beads4)
data5 = load_data(beads5)

data_list = [data1, data2, data3, data4, data5]
# shape (622399, 7)
'''

#read manifests
#probes_file = 'DNAm/python/illumina_manifests/hm450_probes.rds'
#controls_file = 'DNAm/python/illumina_manifests/hm450_controls.rds'

#quick fix -csv file
probes_file = 'CH3/python/testing/probes_1.csv'
controls_file = 'CH3/python/testing/control_beads.csv'
idat_files_folder = 'idat'



##Preprocessing

#for id, df in enumerate(data_list):
samples, cpgs, snps, intensities_A, intensities_B, controls_red, controls_grn = preprocess(probes_file, controls_file,
    idat_files_folder, min_beads=3, detection=0.05, return_intensities=True)
print ('----------------------')
print (samples)
print ('----------------------')

print(timeit.default_timer() - start_time)


## Quality control

snps_distribution_box(snps, i, samples, cpgs)


#, remove_unreliable_samples, k_mean_sex_infer, infer_sex, snps_distribution, identify_replicates, compare_sex, estimate_leukocytes



