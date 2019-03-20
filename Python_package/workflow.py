## load dependencies
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
from DNAm.python.illuminaio import list_idat
from DNAm.python.preprocess import load_data, read_manifests, preprocess

##load data for testing
data_file = '/Users/metzlerabarbara/Library/Mobile Documents/com~apple~CloudDocs/dnam/R05C01_beads.csv'
data = load_data(data_file)
# shape (622399, 7)
#data = data[1:200000]


#read manifests
#probes_file = 'DNAm/python/illumina_manifests/hm450_probes.rds'
#controls_file = 'DNAm/python/illumina_manifests/hm450_controls.rds'

#quick fix -csv file
probes_file = '/Users/metzlerabarbara/OneDrive - Imperial College London/IMPERIAL/CE/Week 1/Practical1/Data/preprocessing/probes_1.csv'
controls_file = '/Users/metzlerabarbara/OneDrive - Imperial College London/IMPERIAL/TDS/control_beads.csv'
idat_files_folder = '/Users/metzlerabarbara/OneDrive - Imperial College London/IMPERIAL/CE/Week 1/Practical1/Data/preprocessing/idat/'


probes, controls = read_manifests(probes_file, controls_file)

##Preprocessing
samples, cpgs, snps, intensities_A, intensities_B, controls_red, controls_grn = preprocess(data, probes, controls,
    idat_files_folder, arg_beads=3, arg_detection=0.05, return_intensities=True)

print (samples)