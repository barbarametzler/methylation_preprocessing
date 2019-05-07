# -*- coding: utf-8 -*-
# For license information, see LICENSE.TXT


import sys
import os
import pyreadr
import pandas as pd
import numpy as np
from scipy.stats import norm
from sklearn import preprocessing

from CH3.python.illuminaio import list_idat


def preprocess(probes_file, controls_file, idat_files_folder, min_beads=3, detection=0.05, return_intensities=False, return_snps_r=False):
    """
    Preprocesses Illumina Infinium DNA methylation bead chips

    Parses .idat files and performs probe censoring, background subtraction and dye-bias correction. 
    SNPs and summary statistics including information from control beads are also provided.

    Parameters
    -----------
        data (): 
        probes():
        idat_files_folder (path): path to folder containing .idat files (2 per sample)
        min_beads (int, optional): probes with less beads will be censored (default 3)
        detection (float, optional): p-value for probe-detection, probes that aren't significantly different from negative control beads are censored (default 0.05)
        return_intensities (bool, optional): returns four (large) matrices containing preprocessed intensities: intensities_A, intensities_B and controls_red, controls_grn
        return_snps_r (bool, optional): returns matrix containing SNP r-coordinate in polar coordinate system
        verbose (bool, optional): prints timestamp per sample and overall time taken

    Returns
    --------
        Returns list with at least three elements
        
        dataframe 
            containing sample metadata and summary statistics
        dataframe
            containing CpG beta-values (= methylation proportions)
        dataframe
            containing theta-values of SNPs in polar coordinate system
        dataframe  
            containing r-values of SNPs in polar coordinate system (optional))
        dataframe 
            containing unmethylated-intensities (A-beads for Illumina I, red channel for Illumina II, optional))
        dataframe
            containing methylated-intensities (B-beads for Illumina I, green channel for Illumina II, optional))
        dataframe
            containing control-bead-intensities on red channel (optional))
        dataframe
            containing control-bead-intensities on green channel (optional))

        {
            'idat_files_folder': idat_files_folder
            'min_beads': min_beads
            'detection': detection
            'return_intensities': return_intensities
            'return_snps_r': return_snps_r
            'verbose': verbose

        }
    """

    ## check argument values
    assert detection > 0 
    assert detection <= 1
    assert min_beads > 0


    def load_data(csv_file):
        data = pd.read_csv(csv_file)
        data.columns = ['probe_address', 'grn_n', 'grn_mean', 'grn_sd', 'red_n', 'red_mean', 'red_sd']
        data.set_index(['probe_address'], inplace=True)
        return data

    def read_manifests(probes_file, controls_file):
        #result_p = pyreadr.read_r(probes_file)
        #probes = pd.DataFrame((result_p[None])) #columns=['chr', 'pos', 'type', 'address_a', 'address_b'])

        #result_c = pyreadr.read_r(controls_file)
        #controls = pd.DataFrame((result_c[None])) #columns=['type', 'color', 'description', 'comment'])
        
        ## temporary solution (pyreader does not recognise the index labels)
        controls = pd.read_csv(controls_file, low_memory=True)
        controls.set_index(['Unnamed: 0'], inplace=True)
        controls.index.names = ['sample_id']

        probes = pd.read_csv(probes_file, low_memory=True)
        probes.set_index(['Unnamed: 0'], inplace=True)
        probes.index.names = ['probe_address']
        
        return probes, controls


    ### quick solution - needs to be removed in final package

    beads5 = 'CH3/python/dnam/R05C01_beads.csv'
    beads4 = 'CH3/python/dnam/R04C01_beads.csv'
    beads3 = 'CH3/python/dnam/R03C01_beads.csv'
    beads2 = 'CH3/python/dnam/R02C01_beads.csv'
    beads1 = 'CH3/python/dnam/R01C01_beads.csv'


    data1 = load_data(beads1)
    data2 = load_data(beads2)
    data3 = load_data(beads3)
    data4 = load_data(beads4)
    data5 = load_data(beads5)

    data_list = [data1, data2, data3, data4, data5]

    probes, controls = read_manifests(probes_file, controls_file)


    ## preparation of outputs
    inf1grn = probes[probes['type'] == "I-Grn"]
    inf1red = probes[probes['type'] == "I-Red"]
    inf2 = probes[probes['type'] == "II"]
    idat_files = list_idat(idat_files_folder)

    ## create intensities
    ## create empty dataframes to append to
    intensities_A = pd.DataFrame(np.nan, index=probes.index, columns=np.sort(pd.unique(idat_files['sample.id'])))
    intensities_B = pd.DataFrame(np.nan, index=probes.index, columns=np.sort(pd.unique(idat_files['sample.id'])))

    controls_grn = pd.DataFrame(np.nan, index=controls.index, columns=np.sort(pd.unique(idat_files['sample.id'])))
    controls_red = pd.DataFrame(np.nan, index=controls.index, columns=np.sort(pd.unique(idat_files['sample.id'])))


    ##Separation of unmethylated and methylated intensities
    # Separate Grn/Red intensities into A (unmethylated) and B (methylated) intensities
    # depending on their type

    ad_a_grn = (probes['address.a'].loc[probes['type'] == 'I-Grn'])
    ad_b_grn = (probes['address.b'].loc[probes['type'] == 'I-Grn'])
    ad_a_red = (probes['address.b'].loc[probes['type'] == 'I-Red'])
    ad_b_red = (probes['address.b'].loc[probes['type'] == 'I-Red'])
    ad_a_inf = (probes['address.a'].loc[probes['type'] == 'II'])

    inf1grn = probes[probes['type'] == "I-Grn"].index
    inf1red = probes[probes['type'] == "I-Red"].index
    inf2 = probes[probes['type'] == "II"].index

    con_ind = controls.index
    sample_ids = np.sort(pd.unique(idat_files['sample.id']))


    #loop over sample id index and fill out rows based on if value is in data

    for column, data in zip(intensities_A, data_list):  
        intensities_A[column].loc[inf1grn] = np.where(data.loc[ad_a_grn, 'grn_n'] >= min_beads, data.loc[ad_a_grn, 'grn_mean'], np.nan)
        intensities_A[column].loc[inf1red] = np.where(data.loc[ad_a_red, 'red_n'] >= min_beads, data.loc[ad_a_red, 'red_mean'], np.nan)
        intensities_A[column].loc[inf2] = np.where(data.loc[ad_a_inf, 'grn_n'] >= min_beads, data.loc[ad_a_inf, 'grn_mean'], np.nan)

    

        ##checked - the same

    for column, data in zip(intensities_B, data_list):
        intensities_B[column].loc[inf1grn] = np.where(data.loc[ad_b_grn, 'grn_n'] >= min_beads, data.loc[ad_b_grn, 'grn_mean'], np.nan)
        intensities_B[column].loc[inf1red] = np.where(data.loc[ad_b_red, 'red_n'] >= min_beads, data.loc[ad_b_red, 'red_mean'], np.nan)
        intensities_B[column].loc[inf2] = np.where(data.loc[ad_a_inf, 'red_n'] >= min_beads, data.loc[ad_a_inf, 'grn_mean'], np.nan)

    #print (intensities_B.head())
    #print (intensities_A.describe())
    #print (intensities_B.describe())

    z = norm.ppf(1 - detection)

    neg_means_grn_list = []
    neg_sds_grn_list = []
    threshold_inf1grn_list = []

    for column, data in zip(controls_grn, data_list):
        controls_grn[column].loc[con_ind] = np.where(data.loc[con_ind, 'grn_n'] > 0,
                                                      data.loc[con_ind, 'grn_mean'],
                                                      np.nan)
        #checked- the same

        neg_beads_grn = controls[controls['type'] == "NEGATIVE"].index
        neg_means_grn = (controls_grn[column].loc[neg_beads_grn]).mean()
        neg_sds_grn = (controls_grn[column].loc[neg_beads_grn]).std(axis=0)
        threshold_inf1grn = 2 * neg_means_grn + z * np.sqrt(2) * neg_sds_grn
        neg_means_grn_list.append(neg_means_grn)
        neg_sds_grn_list.append(neg_sds_grn)
        threshold_inf1grn_list.append(threshold_inf1grn)


    neg_means_red_list = []
    neg_sds_red_list = []
    threshold_inf1red_list = []

    for column, data in zip(controls_red, data_list):
        controls_red[column].loc[con_ind] = np.where(data.loc[con_ind, 'red_n'] > 0,
                                                      data.loc[con_ind, 'red_mean'],
                                                      np.nan)

        neg_beads_red = controls[controls['type'] == "NEGATIVE"].index
        neg_means_red = (controls_red[column].loc[neg_beads_red]).mean()
        neg_sds_red = (controls_red[column].loc[neg_beads_red]).std(axis=0)
        threshold_inf1red = 2 * neg_means_red + z * np.sqrt(2) * neg_sds_red

        neg_means_red_list.append(neg_means_red)
        neg_sds_red_list.append(neg_sds_red)
        threshold_inf1red_list.append(threshold_inf1red)


    threshold_inf2_list =[]

    ## Defining a threshold of detection
    for x, y, j, i in zip(neg_means_grn_list, neg_means_red_list, neg_sds_grn_list, neg_sds_red_list):
        neg_means_ = (x + y)
        threshold_inf2 = neg_means_ + z * np.sqrt(j**2 + i**2 )
        threshold_inf2_list.append(threshold_inf2)

    ## thresholds are the same

    # Censoring of values below the detection limit and background subtraction
    # Background subtraction


    I = pd.concat([intensities_A, intensities_B], axis=1)

    for column in intensities_A:
        I[column] = I[column].sum(axis=1) 

    I = I.loc[:,~I.columns.duplicated()]
    print (I)

    for column, x, y, z, a, b in zip(intensities_A, threshold_inf1grn_list, threshold_inf1red_list, threshold_inf2_list, neg_means_grn_list, neg_means_red_list):

        #print ((intensities_A[column].loc[inf1grn]).describe()) this is still the same

        #print ((intensities_A[column].loc[inf1grn] > a).describe()) # this is correct 
        #print ((intensities_A[column].loc[inf1grn] > x).describe()) this is not the same
        #print (I[column])


        intensities_A[column].loc[inf1grn] = (np.where((intensities_A[column].loc[inf1grn] > a) & (I[column].loc[inf1grn] > x),
                                                    (intensities_A[column].loc[inf1grn] - a),
                                                    np.nan))

        #print ((intensities_A[column].loc[inf1grn] > a).sum())

        #print ((intensities_A[column].loc[inf1grn] > 301.7267).describe())
        #print ((intensities_A[column].loc[inf1grn] - a).describe())

        #print ((intensities_A[column].loc[inf1grn]).describe()) not the same anymore

        intensities_A[column].loc[inf1red] = (np.where((intensities_A[column].loc[inf1red] > b) & (I[column].loc[inf1red] > y),
                                                    (intensities_A[column].loc[inf1red] - b),
                                                    np.nan))   

        intensities_A[column].loc[inf2] = (np.where((intensities_A[column].loc[inf2] > b) & (I[column].loc[inf2] > z),
                                                    (intensities_A[column].loc[inf2] - b),
                                                    np.nan))
    

    for column, x, y, z, a, b in zip(intensities_B, threshold_inf1grn_list, threshold_inf1red_list, threshold_inf2_list, neg_means_grn_list, neg_means_red_list):


        intensities_B[column].loc[inf1grn] = (np.where((intensities_B[column].loc[inf1grn] > a) & (I[column].loc[inf1grn] > x),
                                                    (intensities_B[column].loc[inf1grn] - a),
                                                    np.nan))

        intensities_B[column].loc[inf1red] = (np.where((intensities_B[column].loc[inf1red] > b) & (I[column].loc[inf1red] > y),
                                                    (intensities_B[column].loc[inf1red] - b),
                                                    np.nan))   

        intensities_B[column].loc[inf2] = (np.where((intensities_B[column].loc[inf2] > a) & (I[column].loc[inf2] > z),
                                                    (intensities_B[column].loc[inf2] - a),
                                                    np.nan))


    
    # Extract normalization probes for Grn and Red, and form the dye bias correction constant
    norm_grn_beads = controls[controls['type'].isin(['NORM_C', 'NORM_G'])].index
    norm_red_beads1 = controls['description'].loc[norm_grn_beads].str.translate(str.maketrans('CG', 'TA')) 
    norm_red_beads = controls[controls['description'].isin(norm_red_beads1)].index

    corrections_grn_list = []
    corrections_red_list = []


    for column1, column2 in zip(controls_grn, controls_red):
        grn = controls_grn[column1].loc[norm_grn_beads].values
        red = controls_red[column2].loc[norm_red_beads].values
        
        norm_data = pd.DataFrame([grn, red]).T

        corrections_grn = (np.mean(norm_data,  axis=1)/ grn).mean(axis=0)
        corrections_red = (np.mean(norm_data, axis=1)/ red).mean(axis=0)

        corrections_grn_list.append(corrections_grn)
        corrections_red_list.append(corrections_red)


    #print (corrections_red_list)
    ## Apply dye bias correction
    for column1, column2, x, y in zip(intensities_A, intensities_B, corrections_red_list, corrections_grn_list):
        
        #print ((intensities_A[column1].loc[inf2]).head())

        intensities_A[column1].loc[inf2] = intensities_A[column1].loc[inf2] * x
        intensities_B[column2].loc[inf2] = intensities_B[column2].loc[inf2] * y

        #print (intensities_A[column].isna().sum())
        print ("****")
        #print (intensities_B[column].isna().sum())
 
    #print (intensities_A.isna().sum()) 
    #print (intensities_B.shape) #(485577, 5)


    ## Computing DNA methylation ratios (β values)
    #Some of the probes are SNPs (N=65), these can be identified because they start with the prefix “rs”

    # Create DNAm ratios as B (methylated) over total

    idx = probes[probes.index.str.startswith('rs')].index ##'rs10796216', 'rs715359',..

    intensities = intensities_A.add(intensities_B, fill_value=0)
    intensities.set_index(probes.index, inplace=True)

    intensities_A_wo = intensities_A.index.isin(idx)
    intensities_B_wo = intensities_B.index.isin(idx)
    intensities_wo = intensities.index.isin(idx)
    
    #print (intensities_B[~intensities_B_wo]) #485512, 5
 
    dnam = intensities_B[~intensities_B_wo] / np.sum(intensities[~intensities_wo])

    '''
    x = dnam1.values #returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    dnamm = pd.DataFrame(x_scaled)

    dnam = dnamm.set_index(dnam1.index)
    dnam.columns=dnam1.columns.values

    #print (dnam.head(50))
    #print (dnam.info)
    #print (dnam.iloc[:,4].mean())
    
    '''

    ## SNPS

    snps_theta = pd.DataFrame(np.nan, index=probes.loc[idx].index, columns=np.sort(pd.unique(idat_files['sample.id'])))
    snps_r = pd.DataFrame(np.nan, index=probes.loc[idx].index, columns=np.sort(pd.unique(idat_files['sample.id'])))
    
    idx = probes[probes.index.str.contains('rs')].index
    snps_theta[column] = np.arctan2(intensities_B[intensities_B_wo][column], intensities_A[intensities_A_wo][column])/ (np.pi/2)
    snps_r[column] = np.sqrt(np.sum(intensities.loc[idx]**2))

    ### matching 
    # Extract all control probes data, and add summary statistics to samples table

    summary = pd.DataFrame(np.nan, columns=np.sort(pd.unique(idat_files['sample.id'])), 
        index=['bc1_grn', 'bc1_red', 'bc2', 'ext_a', 'ext_c', 'ext_g', 'ext_t',
                'hyp_low', 'hyp_med', 'hyp_high', 'np_a', 'np_c',
                'np_g', 'np_t', 'spec1_grn', 'spec1_red', 'spec2', 'st_grn', 'st_red',
                'tr', 'missing', 'median_chrX', 'missing_chrY'])

    for column in summary:

        # match 1

        bg = ['BS Conversion I-U1', 'BS Conversion I-U2','BS Conversion I-U3']
        idx_bg = (controls[controls['description'].str.contains('|'.join(bg))]).index
        match_ = controls['description'].loc[idx_bg].str.translate(str.maketrans('U', 'C'))
        idx_signal = controls['description'].isin(match_)
        new_ind = controls_grn[idx_signal].index
        summary[column].loc['bc1_grn'] = (np.nanmean(controls_grn[column].loc[new_ind])/np.mean(controls_grn[column].loc[idx_bg]))

        bg = ['BS Conversion I-U4', 'BS Conversion I-U5','BS Conversion I-U6']
        idx_bg = (controls[controls['description'].str.contains('|'.join(bg))]).index
        match_ = controls['description'].loc[idx_bg].str.translate(str.maketrans('U', 'C'))
        idx_signal = controls['description'].isin(match_)
        new_ind = controls_red[idx_signal].index
        summary[column].loc['bc1_red'] = (np.nanmean(controls_red[column].loc[new_ind])/np.mean(controls_red[column].loc[idx_bg]))

        idx = controls[controls['type'] == 'BISULFITE CONVERSION II'].index
        summary[column].loc['bc2'] = np.nanmean(controls_red[column].loc[idx]/np.nanmean(controls_grn[column].loc[idx]))

        idd = controls[controls['description'] == 'Extension (A)'].index
        summary[column].loc['ext_a'] = controls_red[column].loc[idd].values
        
        idd = controls[controls['description'] == "Extension (C)"].index
        summary[column].loc['ext_c'] = controls_grn[column].loc[idd].values 

        idd = controls[controls['description'] == 'Extension (G)'].index
        summary[column].loc['ext_g'] = controls_grn[column].loc[idd].values

        idd = controls[controls['description'] == 'Extension (T)'].index
        summary[column].loc['ext_t'] = controls_red[column].loc[idd].values
        
        idd = controls[controls['description'] == 'Hyb (Low)'].index
        summary[column].loc['hyp_low'] = controls_grn[column].loc[idd].values 

        idd = controls[controls['description'] == 'Hyb (Medium)'].index
        summary[column].loc['hyp_med'] = controls_grn[column].loc[idd].values 

        idd = controls[controls['description'] == 'Hyb (High)'].index
        summary[column].loc['hyp_high'] = controls_grn[column].loc[idd].values 

        idd = controls[controls['description'] == 'NP (A)'].index
        summary[column].loc['np_a'] = controls_red[column].loc[idd].values 

        idd = controls[controls['description'] == 'NP (C)'].index
        summary[column].loc['np_c'] = controls_grn[column].loc[idd].values 

        idd = controls[controls['description'] == 'NP (G)'].index
        summary[column].loc['np_g'] = controls_grn[column].loc[idd].values 

        idd = controls[controls['description'] == 'NP (T)'].index
        summary[column].loc['np_t'] = controls_red[column].loc[idd].values 

        # match 2

        controls[controls['description'] == 'GT Mismatch 1 (MM)'] = 'gt_mismatch_1_mm'
        controls[controls['description'] == 'GT Mismatch 2 (MM)'] = 'gt_mismatch_2_mm'
        controls[controls['description'] == 'GT Mismatch 3 (MM)'] = 'gt_mismatch_3_mm'
        bg = ['gt_mismatch_1_mm', 'gt_mismatch_2_mm', 'gt_mismatch_3_mm']
        idx_bg = controls[controls['description'].str.contains('|'.join(bg))].index

        controls[controls['description'] == 'GT Mismatch 1 (PM)'] = 'gt_mismatch_1_pm'
        controls[controls['description'] == 'GT Mismatch 2 (PM)'] = 'gt_mismatch_2_pm'
        controls[controls['description'] == 'GT Mismatch 3 (PM)'] = 'gt_mismatch_3_pm'
        signal = ['gt_mismatch_1_pm', 'gt_mismatch_2_pm', 'gt_mismatch_3_pm']
        idx_signal = controls[controls['description'].str.contains('|'.join(signal))].index

        summary[column].loc['spec1_grn'] = (np.nanmean(controls_grn[column].loc[idx_signal])/np.mean(controls_grn[column].loc[idx_bg]))

    
        controls[controls['description'] == 'GT Mismatch 4 (MM)'] = 'gt_mismatch_4_mm'
        controls[controls['description'] == 'GT Mismatch 5 (MM)'] = 'gt_mismatch_5_mm'
        controls[controls['description'] == 'GT Mismatch 6 (MM)'] = 'gt_mismatch_6_mm'
        bg = ['gt_mismatch_4_', 'gt_mismatch_5_', 'gt_mismatch_6_']
        idx_bg = controls[controls['description'].str.contains('|'.join(bg))].index

        controls[controls['description'] == 'GT Mismatch 4 (PM)'] = 'gt_mismatch_4_pm'
        controls[controls['description'] == 'GT Mismatch 5 (PM)'] = 'gt_mismatch_5_pm'
        controls[controls['description'] == 'GT Mismatch 6 (PM)'] = 'gt_mismatch_6_pm'
        signal = ['gt_mismatch_4_pm', 'gt_mismatch_5_pm', 'gt_mismatch_6_pm']
        idx_signal = controls[controls['description'].str.contains('|'.join(signal))].index

        summary[column].loc['spec1_red'] = (np.nanmean(controls_red[column].loc[idx_signal])/np.mean(controls_red[column].loc[idx_bg]))


        idx = controls[controls['type'] == 'SPECIFICITY II'].index
        summary[column].loc['spec2'] = (np.nanmean(controls_red[column].loc[idx])/np.mean(controls_grn[column].loc[idx]))

        # match 3 
        idx_bg = controls[controls['description'] == ('Biotin (Bkg)')].index
        idx_signal = controls[controls['description'] == ('Biotin (High)')].index
        summary[column].loc['st_grn'] = controls_grn[column].loc[idx_signal].values / controls_grn[column].loc[idx_bg].values

        # match 4
        idx_bg = controls[controls['description'] == ('DNP (Bkg)')].index
        idx_signal = controls[controls['description'] == ('DNP (High)')].index
        summary[column].loc['st_red'] = controls_red[column].loc[idx_signal].values / controls_red[column].loc[idx_bg].values
        
        # match 5
        idx = controls[controls['type'] == ('TARGET REMOVAL')].index
        summary[column].loc['tr'] = np.nanmax(controls_grn[column].loc[idx], axis=0)
        summary[column].loc['missing'] = (dnam[column].isna().sum()) / len(dnam)

        # match 6
        idx = probes[probes['chr'] == 'X'].index
        summary[column].loc['median_chrX'] = np.nanmedian(dnam[column].loc[idx])
        #print (np.nanmedian(dnam[column].loc[idx]).round(4))

        idy = probes[probes['chr'] == 'Y'].index
        ## less missing values as in R code? 
        print ("missing values dnam per column")
        print (dnam[column].loc[idy].isna().sum())
        summary[column].loc['missing_chrY'] = dnam[column].loc[idy].isna().mean()


        ## 
    #cpgs = pd.read_csv("CH3/python/testing/cpgs.csv",index_col='Unnamed: 0', engine='c')


    samples = summary.T
    cpgs = dnam.T
    snps = snps_theta.T


    if return_intensities == True:
        return samples, cpgs, snps, intensities_A, intensities_B, controls_red, controls_grn

    elif return_snps_r == True:

        snps_r = np.sqrt(np.sum(intensities[idx]**2))
        
        return snps_r

    else:
        return samples, cpgs, snps




