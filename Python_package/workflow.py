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
csv_file= 'DNAm/python/dnam/R05C01_beads.csv'
data = load_data(csv_file)


#read manifests
probes_file = 'DNAm/python/illumina_manifests/hm450_probes.rds'
#controls_file = 'DNAm/python/illumina_manifests/hm450_controls.rds'
#quick fix -csv file
#controls_file = 

probes, controls = read_manifests(probes_file, controls_file)


##Preparation of the outputs
idat_files_folder = 'idat/'

preperation_outputs(probes, idat_files_folder)

##Create intensities
intensities, controls = create_intensities(data, probes, controls, idat_files, arg_beads=3, )
