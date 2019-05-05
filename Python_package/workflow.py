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
from CH3.python.quality_control_1 import visualisation_plots

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

print(tabulate(samples, headers='keys', tablefmt='psql'))

#print (samples.info)
#print (samples["7800246024_R05C01"].round(4))


print(timeit.default_timer() - start_time)


#print (type(samples)) #[5 rows x 23 columns]
#print (type(cpgs)) #[5 rows x 485577 columns]n
#print (type(snps)) # [5 rows x 65 columns]
#print (type(intensities_A)) # [485577 rows x 5 columns]
#print (type(intensities_B)) # [485577 rows x 5 columns]
#print (type(controls_grn)) #[835 rows x 5 columns]
#print (type(controls_red)) #[835 rows x 5 columns]

#samples.to_csv('samples_pre.csv')
#cpgs.to_csv('cpgs_pre.csv')
#snps.to_csv('snps_pre.csv')


## Quality control
# create 3 plots

'''
# Boxplot of the 'bc1.grn','bc1.red','bc2' for sample i
#samples = samples.T
df = samples[['bc1_grn','bc1_red','bc2']]
sns.boxplot(x="variable", y="value", data=pd.melt(df)).set_title("Boxplot")
plt.show()

# Plot of the histgram of the missing sample per row with the threshold
samples_prop_na=samples['missing']
thresh = 0.1
plt.hist(samples_prop_na)
plt.title('Distribution of missing variables per row')
plt.axvline(x=thresh, color='r', linestyle='dashed', linewidth=2)
plt.show()

# Plot of the histgram of the missing sample per column with the threshold
samples_col_missing = samples.isnull().mean(axis=0)
thresh = 0.1
plt.hist(samples_col_missing,bins=100, edgecolor="none")
plt.title('Distribution of missing variables per column')
plt.axvline(x=thresh, color='r', linestyle='dashed', linewidth=2)
plt.show()

# Plot of the methylation β-value-distribution, expecting a bimodal distribution
common_1=cpgs.index.intersection(samples.index)
cpgs=cpgs.loc[common_1]
cpgs_1=cpgs.iloc[i,:]
sns.distplot(cpgs_1, hist=True, kde=True, 
bins=int(180/5), color = 'darkblue', 
hist_kws={'edgecolor':'black'},
kde_kws={'linewidth': 4}).set_title("Methylation β-value-distribution")
plt.show()
'''

