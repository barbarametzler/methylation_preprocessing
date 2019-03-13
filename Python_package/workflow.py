## load dependencies
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
from read_idat import list_idat
import preprocess

##load data for testing
# define paths
csv_file = '/Users/metzlerabarbara/Library/Mobile Documents/com~apple~CloudDocs/dnam/R05C01_beads.csv'
data = load_data(csv_file)


#read manifests
probes_file = '/Users/metzlerabarbara/OneDrive - Imperial College London/IMPERIAL/CE/Week 1/Practical1/Data/preprocessing/illumina_methylation/manifests/hm450_probes.rds'
controls_file = '/Users/metzlerabarbara/OneDrive - Imperial College London/IMPERIAL/CE/Week 1/Practical1/Data/preprocessing/illumina_methylation/manifests/hm450_controls.rds'

probes, controls = read_manifests(probes_file, controls_file)


##Preparation of the outputs
idat_files_folder = '/Users/metzlerabarbara/OneDrive - Imperial College London/IMPERIAL/CE/Week 1/Practical1/Data/preprocessing/idat/'

preperation_outputs(probes, idat_files_folder)

##Create intensities
intensities, controls = create_intensities(probes, controls, idat_files, arg_beads=3, data)
