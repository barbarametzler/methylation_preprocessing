# -*- coding: utf-8 -*-
# For license information, see LICENSE.TXT

import pandas as pd
import numpy as np
import sys
import os
import pyreadr
from read_idat import list_idat


def load_data(csv_file):
    data = pd.read_csv(csv_file)
    data.columns = ['sample_id', 'grn_n', 'grn_mean', 'grn_sd', 'red_n', 'red_mean', 'red_sd']
    return data

def read_manifests(probes_file, controls_file):
    result_p = pyreadr.read_r(probes_file)
    probes = pd.DataFrame((result_p[None])) #columns=['chr', 'pos', 'type', 'address_a', 'address_b'])

    result_c = pyreadr.read_r(controls_file)
    #controls = pd.DataFrame((result_c[None])) #columns=['type', 'color', 'description', 'comment'])
    
    ## temporary solution (pyreader does not recognise the index labels)
    controls = pd.read_csv(controls_file)
    controls.set_index(['Unnamed: 0'], inplace=True)
    controls.index.names = ['sample_id']
    return probes, controls

def preperation_outputs(probes, idat_files_folder):
    inf1grn = probes[probes['type'] == "I-Grn"]
    inf1red = probes[probes['type'] == "I-Red"]
    inf2 = probes[probes['type'] == "II"]
    idat_files = list_idat(idat_files_folder)
    return inf1grn, inf1red, inf2, idat_files


def create_intensities(data, probes, controls, idat_files, arg_beads=3, arg_detection = 0.05):
    ## create empty dataframes to append to
    intensities_A = pd.DataFrame(np.nan, index=probes.index, columns=pd.unique(idat_files['sample.id']))
    intensities_B = pd.DataFrame(np.nan, index=probes.index, columns=pd.unique(idat_files['sample.id']))

    controls_grn = pd.DataFrame(np.nan, index=controls.index, columns=pd.unique(idat_files['sample.id']))
    controls_red = pd.DataFrame(np.nan, index=controls.index, columns=pd.unique(idat_files['sample.id']))


    ##Separation of unmethylated and methylated intensities
    # Separate Grn/Red intensities into A (unmethylated) and B (methylated) intensities
    # depending on their type


    ad_a_grn = (probes['address.a'].loc[probes['type'] == 'I-Grn']).index
    ad_b_grn = (probes['address.b'].loc[probes['type'] == 'I-Grn']).index
    ad_a_red = (probes['address.b'].loc[probes['type'] == 'I-Red']).index
    ad_b_red = (probes['address.b'].loc[probes['type'] == 'I-Red']).index
    ad_a_inf = (probes['address.a'].loc[probes['type'] == 'inf2']).index

    inf1grn = probes[probes['type'] == "I-Grn"].index
    inf1red = probes[probes['type'] == "I-Red"].index
    inf2 = probes[probes['type'] == "II"].index

    con_ind = controls.index
    sample_ids = pd.unique(idat_files['sample.id'])

    ## temporary solution, since data for R05C01 is only available
    intensities_AA = pd.DataFrame(intensities_A['7800246024_R05C01'])
    intensities_BB = pd.DataFrame(intensities_B['7800246024_R05C01'])

    #loop over sample id index and fill out rows based on if value is in data

    for column in intensities_AA:
        #if (data['grn_n'].isin(ad_a_grn).any()) >= arg_beads:
        intensities_AA[column].loc[inf1grn] = (np.where(data.loc[ad_a_grn, 'grn_n'] >= arg_beads, data.loc[ad_a_grn, 'grn_mean'], np.nan))
        intensities_AA[column].loc[inf1red] = (np.where(data.loc[ad_a_red, 'red_n'] >= arg_beads, data.loc[ad_a_red, 'red_mean'], np.nan))
        intensities_AA[column].loc[inf2] = (np.where(data.loc[ad_a_inf, 'grn_n'] >= arg_beads, data.loc[ad_a_inf, 'grn_mean'], np.nan))


    for column in intensities_BB:
        intensities_BB[column].loc[inf1grn] = (np.where(data.loc[ad_b_grn, 'grn_n'] >= arg_beads, data.loc[ad_b_grn, 'grn_mean'], np.nan))
        intensities_BB[column].loc[inf1red] = (np.where(data.loc[ad_b_red, 'red_n'] >= arg_beads, data.loc[ad_b_red, 'red_mean'], np.nan))
        intensities_BB[column].loc[inf2] = (np.where(data.loc[ad_a_inf, 'red_n'] >= arg_beads, data.loc[ad_a_inf, 'grn_mean'], np.nan))

     for column in controls_grn:
        dataa = data.set_index('sample_id')
        controls_grn[column].loc[con_ind] = np.where(dataa.loc[con_ind, 'grn_n'] > 0,
                                                      dataa.loc[con_ind, 'grn_mean'],
                                                      np.nan)

        neg_beads_grn = controls[controls['type'] == "NEGATIVE"].index
        neg_means_grn = (controls_grn[column].loc[neg_beads_grn]).mean()
        neg_sds_grn = (controls_grn[column].loc[neg_beads_grn]).std(axis=0)

    for column in controls_red:
        controls_red[column].loc[con_ind] = np.where(dataa.loc[con_ind, 'red_n'] > 0,
                                                      dataa.loc[con_ind, 'red_mean'],
                                                      np.nan)

        neg_beads_red = controls[controls['type'] == "NEGATIVE"].index
        neg_means_red = (controls_red[column].loc[neg_beads_red]).mean()
        neg_sds_red = (controls_red[column].loc[neg_beads_red]).std(axis=0)


    ## Defining a threshold of detection
    neg_means_ = np.mean([neg_means_grn, neg_means_red])
    neg_sds_ = np.std([neg_sds_grn, neg_sds_red])

    z = norm.ppf(1 - arg_detection)

    threshold_inf1grn = 2 * neg_means_grn + z * np.sqrt(2) * neg_sds_grn
    threshold_inf1red = 2 * neg_means_red + z * np.sqrt(2) * neg_sds_red
    threshold_inf2 = np.sum(neg_means_) + z * np.sqrt(np.sum(neg_sds_ ** 2))

    # Censoring of values below the detection limit and background subtraction
    # Background subtraction

    
    

    for column in intensities_AA:
        I_A = (intensities_AA[column]).sum(axis=1)
        intensities_AA[column].loc[inf1grn] = (np.where((data.loc[inf1grn] > neg_means_grn &
                                                        I_A.loc[inf1grn] > threshold_inf1grn),
                                                        intensities_AA[column].loc[inf1grn] - neg_means_grn,
                                                        np.nan))


        intensities_AA[column].loc[inf1red] = (np.where((data.loc[inf1red] > neg_means_red &
                                                        I_B.loc[inf1red] > threshold_inf1red),
                                                        intensities_AA[column].loc[inf1red] - neg_means_red,
                                                        np.nan))

        intensities_AA[column].loc[inf2] = (np.where(data.loc[inf2] > neg_means_ &
                                                        I_A.loc[inf2] > threshold_inf2,
                                                        intensities_AA[column].loc[inf2] - neg_means_,
                                                        np.nan))


    for column in intensities_BB:
        I_B = (intensities_B[column]).sum(axis=1)
        intensities_BB[column].loc[inf1grn] = (np.where(data.loc[inf1grn] > neg_means_grn &
                                                        I_B.loc[inf1grn] > threshold_inf1grn,
                                                        intensities_BB[column].loc[inf1grn] - neg_means_grn,
                                                        np.nan))

        intensities_BB[column].loc[inf1red] = (np.where(data.loc[inf1red] > neg_means_red &
                                                        I_B.loc[inf1red] > threshold_inf1red,
                                                        intensities_BB[column].loc[inf1red] - neg_means_red,
                                                        np.nan))

        intensities_BB[column].loc[inf2] = (np.where(data.loc[inf2] > neg_means_ &
                                                        I_B.loc[inf2] > threshold_inf2,
                                                        intensities_AA[column].loc[inf2] - neg_means_,
                                                        np.nan))




# Extract normalization probes for Grn and Red, and form the dye bias correction constant
    norm_grn_beads = control['type'].isin(['NORM_C', 'NORM_G'])




