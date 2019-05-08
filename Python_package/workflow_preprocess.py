## load dependencies
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
import timeit
import seaborn as sns
import matplotlib.pyplot as plt
from tabulate import tabulate
from CH3.python.illuminaio import list_idat
from CH3.python.preprocess import preprocess
from CH3.python.quality_control import visualisation_plots

start_time = timeit.default_timer()

#read manifests
#probes_file = 'DNAm/python/illumina_manifests/hm450_probes.rds'
#controls_file = 'DNAm/python/illumina_manifests/hm450_controls.rds'

#quick fix -csv file
probes_file = 'CH3/python/testing/probes_1.csv'
controls_file = 'CH3/python/testing/control_beads.csv'
idat_files_folder = 'idat'


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

print(tabulate(samples, headers='keys', tablefmt='psql', numalign="right"))


print(timeit.default_timer() - start_time)


#print (type(samples)) #[5 rows x 23 columns]
#print (type(cpgs)) #[5 rows x 485577 columns]n
#print (type(snps)) # [5 rows x 65 columns]
#print (type(intensities_A)) # [485577 rows x 5 columns]
#print (type(intensities_B)) # [485577 rows x 5 columns]
#print (type(controls_grn)) #[835 rows x 5 columns]
#print (type(controls_red)) #[835 rows x 5 columns]




