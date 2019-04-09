#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:30:07 2019

@author: nicolasagrotis
"""

# Quality Control Function architecture:
            
# 1) remove unreliable functions 'remove_unreliable_samples'
            #- missing
            #-outliers
            #- Remove corresponding samples in other list entries

#2) infer sex 'infer_sex':
            # - get the F&M column
            # - maybe try and also disply the results
            # - ask for threshold
            
#3) call SNPS 'call_snps' :
            # - identify the snps using the thresholds
            # - plot? ASK tim
            
#4) snps distance matrix 'snp_distance':
            # - plot the heatmap
            # - display the matrix
            
#5) identify duplicates 'identify_replicates'
            # - print duplicates/replicates

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

#-----------------------------------------------------------------------------------------#

# Data Loading
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


#-----------------------------------------------------------------------------------------#
# visualise the distribution of one of the snps as a boxplot
      
def visualisation_plots(snps,i,samples,cpgs):
    
        # Boxplot of the 'bc1.grn','bc1.red','bc2' for sample i
        df = samples[['bc1.grn','bc1.red','bc2']]
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
        
        
        

#-----------------------------------------------------------------------------------------#    

# Remove variables based on percentage of missing variables and quality check of 
# bc1.grn,bc1.red and bc2
            
def remove_unreliable_samples(samples,threshold,cpgs,covars):
        
    # Apply the threshold for the missing data per row for samples
    samples=samples.loc[samples['missing']<threshold]
    
    # Apply the threshold for the missing data per column for cpgs
    cpgs=cpgs[cpgs.columns[cpgs.isnull().mean() < threshold]]   
    
    # Remove unreliable reading of less than 1
    samples=samples.loc[(samples['bc1.grn']>1)&(samples['bc1.red']>1)&(samples['bc2']>1)]

    # Subset the data based on the previous conditions
    cpgs=cpgs.loc[samples.index.intersection(cpgs.index)]
    
    
    return samples,cpgs,covars




#-----------------------------------------------------------------------------------------#    


# Apply a K-means classification model to investigate the clustering of males and 
# females based on the missing Y chromosome and the median X chromosome

def k_mean_sex_infer(samples):

    from sklearn import datasets
    from sklearn.cluster import KMeans
    import sklearn.metrics as sm
    

    x=samples[['median.chrX','missing.chrY']]
    
    # Apply the Kmeans clustering method with predefined clusters
    model = KMeans(n_clusters=2, init=np.array(((0.25,0.1),(0.5,0.7))),algorithm="elkan")
    model.fit(x)
    centroids = model.cluster_centers_
    
    # View the results
    # Set the size of the plot
    plt.figure(figsize=(14,7))
     
    # Create a colormap
    colormap = np.array(['red', 'lime'])
    
    # Plot the Original Classifications
    plt.subplot(1, 2, 1)
    plt.scatter(x['median.chrX'], x['missing.chrY'], s=40)
    plt.title('Original Plot Without Classification')
    plt.show()
     
    # Plot the Model Classifications
    plt.subplot(1, 2, 2)
    plt.scatter(x['median.chrX'], x['missing.chrY'], c=colormap[model.labels_], s=40)
    plt.title('K Mean Classification')
    plt.show()
    
    # Print the final centroids
    print('The K-Means classification centroids are: ',centroids)
    
    samples['sex_Kmeans']=model.labels_
    
    samples.loc[(samples['sex_Kmeans']==1),'sex_Kmeans']='F'
    samples.loc[(samples['sex_Kmeans']==0),'sex_Kmeans']='M'






#-----------------------------------------------------------------------------------------#    

# Use the print below in the workflow document
# plt.scatter(x=samples['median_chrX'],y=samples['missing_chrY'])


def infer_sex(samples,threshold_chrX=0.37,threshold_chrY=0.39):
    
    
    # From bibliography set hard boundaries to descriminate between males and females
    # hard boundaries: if median_chrX' < 0.37 and missing Y chromosome is smaller than 0.39 then M
                      #if median_chrX' > 0.37 and missing Y chromosome is bigger than 0.39 then F
                      # Otherwise set the value to NaN
                      
    samples.loc[(samples['median.chrX'] < threshold_chrX) & (samples['missing.chrY'] < threshold_chrY), 'sex'] = 'M'
    samples.loc[(samples['median.chrX'] > threshold_chrX) & (samples['missing.chrY'] > threshold_chrY), 'sex'] = 'F'
    samples.loc[(samples['median.chrX'] < threshold_chrX) & (samples['missing.chrY'] > threshold_chrY), 'sex'] = np.nan
    samples.loc[(samples['median.chrX'] > threshold_chrX) & (samples['missing.chrY'] < threshold_chrY), 'sex'] = np.nan    
    
    # Count the number of males and females
    num_males=samples.loc[samples.sex == 'M', 'sex'].count()
    num_females=samples.loc[samples.sex == 'F', 'sex'].count()
    print("Number of Males:",num_males)
    print("Number of Females:",num_females)
    
    genders=['F','M']
    fg = sns.FacetGrid(data=samples, hue='sex', hue_order=genders, aspect=1.61)
    fg.map(plt.scatter, 'median.chrX', 'missing.chrY')
    plt.legend(loc='upper left')
    plt.show()
 
    return samples


#-----------------------------------------------------------------------------------------#    

# snps distribution plot for ith snp
# The allele at a SNP locus can be inferred from SNP intensities measured on the BeadChip.

def call_snps(snps,i):
    
    
        snps=snps.T
        snps_i=snps.iloc[:,1]
        snps_i=pd.DataFrame(snps_i)
        snps_i.columns=['sample']
        
        snps_i.loc[(snps_i['sample'] <= 0.2),'allele']=0
        snps_i.loc[(snps_i['sample'] >= 0.8),'allele']=2        
        snps_i.loc[(snps_i['sample'] > 0.2)&(snps_i['sample'] < 0.8),'allele']=1

        ax=sns.boxplot(x="allele", y="sample", data=snps_i)
        ax.set(xlabel='Carrier Status',ylabel='Theta intensities')
        ax.set_title('Distribution of the 65 SNPS')
        ax.axhline(0.2, ls='--')
        ax.axhline(0.8, ls='--')
        plt.show()
       
      
#-----------------------------------------------------------------------------------------#    


# Need samples['sex'] for this function
# Need the infer_sex function
def identify_replicates(snps,threshold,samples):
    
    #snps.set_index('Unnamed: 0',inplace=True)
    snps.index.name='sample_id'
    
    #set thresholds for snps classifciation
    snps[snps<=0.2]=0
    snps[snps>=0.8]=2
    snps[(snps>0.2 ) & (snps<0.8)]=1
    
    # Create a distance matrix
    dist_matrix = np.empty((snps.shape[0], snps.shape[0]))
    dist_matrix[:,:] = np.nan
    
    # Populate the distance matric with the values of the classification
    for i in range(0,snps.shape[0]):
        for j in range(i+1,snps.shape[0]):
            dist_matrix[j, i] = abs(snps.iloc[i,:]-snps.iloc[j,:]).sum()
            dist_m=pd.DataFrame(dist_matrix)
            dist_m.index=snps.index
            dist_m.columns=snps.index
    
    # Plot a heatmat of the distance matrix       
    ax = sns.heatmap(dist_m, annot=True)
    ax.set_title('Distance Matrix for SNPSs')
    
    dist_m=dist_m < threshold
    
    rows=[]
    columns=[]
    
    for i in range(0,dist_m.shape[0]):
        for j in range(0,dist_m.shape[0]):
            
            if dist_m.iloc[i,j] == True:
                
                rows.append(dist_m.index[i])
                columns.append(dist_m.index[j])
                
    
    sex_ident=samples['sex']
    
    plt.show()

    # Print the identified replicates
    for n in range (0,len(rows)):
            
        if sex_ident[rows[n]] == sex_ident[columns[n]]:
                
            print('Replicate detected:',rows[n],columns[n]) 


#-----------------------------------------------------------------------------------------#

# Compare the infered sex with the true sex available throught the covars dataset

def compare_sex(covars,samples):
    
    
    # Create new columns in samples with the covars gender
    samples['true_sex']=covars['gender']

    # Modify the samples to enable the comparisson
    samples.loc[(samples['true_sex']=='f'),'True_sex']='F'
    samples.loc[(samples['true_sex']=='m'),'True_sex']='M'
    samples.drop(columns='true_sex',inplace=True)
    samples.loc[(samples['sex']== samples['True_sex']),'final_sex']=samples['sex']
    samples.loc[(samples['sex']!= samples['True_sex']),'final_sex']='U'

    # Plot the graph for sex classification and identify the mismmaches
    genders=['F','M','U']
    fg = sns.FacetGrid(data=samples, hue='sex', hue_order=genders, aspect=1.61)
    fg.map(plt.scatter, 'median.chrX', 'missing.chrY').add_legend()
    plt.legend(loc='upper left')
    plt.show()
    
    # Print the IDs of the mismached samples
    mismatch_sex=samples.loc[samples['final_sex']=='U'].index
    print('The mismatched samples are: ',mismatch_sex)
    
    # Print the % of mismatched samples
    sum_sex=samples['sex'].loc[(samples['sex'])==(samples['True_sex'])].count()
    perc_sex=(sum_sex/samples['sex'].count())*100
    print('The percentage of correctly infered sex samples is: ',perc_sex)
    
    return samples
#-----------------------------------------------------------------------------------------#

def estimate_leukocytes(coefs,cpgs):
  
    cpgs=cpgs.T
    
    common=coefs.index.intersection(cpgs.index)
    
    cpgs=cpgs.loc[common]
    
    
    betas=cpgs
    
    cpgs_in=len(betas.columns)
    coefs_col=len(coefs.columns)
    
    wbc_predictions = np.zeros([cpgs_in,coefs_col])
    
    A = np.identity(len(coefs.columns))
    
    b = np.repeat(0,len(coefs.columns))
    
    #D_1=coefs.loc[idx,:].values
    
    for i in range(wbc_predictions.shape[0]):
        idx=betas[betas.iloc[:,i].notna()].index
        betas=betas.loc[idx]
        beta=betas.iloc[:,i]
        D=np.matmul(coefs.loc[idx,:].T.values,coefs.loc[idx,:].values)
        d=np.matmul(coefs.loc[idx,:].T.values,betas.values)
        wbc_predictions[i,:]=quadprog_solve_qp(D,-d,A,b)
    
    #cpgs.to_csv('coef_clean.csv')

#$ sudo pip install quadprog

def quadprog_solve_qp(P, q, G, A):
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if A is not None:
        qp_C = -np.vstack([A, G]).T
        #qp_b = -np.hstack([b, h])
        meq = A.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_C, meq)[0]

# experiment with sign of d(-d/d)
# b=x
# D=P
# d=-q
# A=G
# b_0=h
    
