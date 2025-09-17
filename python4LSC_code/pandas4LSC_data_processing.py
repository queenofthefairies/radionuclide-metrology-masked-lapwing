# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 14:55:29 2017

@author: siobhant
"""
import os
import pandas as pd
import numpy as np
import datetime

#==============================================================================
# Basic stats functions
#==============================================================================
def weighted_mean_calc(unc_temp, n0_temp, weight_mean):
    weight = 1/unc_temp**2
    if weight_mean == 1:
        wght_x_n0 = (n0_temp)*weight
        wm_out = np.nansum(wght_x_n0)/np.nansum(weight)
    else:
        wm_out = np.nanmean(n0_temp)    
    return wm_out, weight 

def stdev_mean_theor_calc(weight):
    stdev_mean_theor = np.sqrt(1/np.nansum(weight))    
    return stdev_mean_theor

#==============================================================================
# Get data and put it into a dataframe
#==============================================================================
def get_data(file_dir,file_list, file_type = 'excel'):

    for i in range(len(file_list)):
        file_name = file_list[i][0]
        file_threshold_mV = file_list[i][1]
        file_path="{0}/{1}".format(file_dir,file_name)

        # na_values = 0 means it treats zero counts as corrupted files / ignore them
        # TO READ CSV FILES
        if file_type == 'csv':
            thresh_df_i = pd.read_csv(file_path, sep=',', skiprows=5, na_values=0, index_col=0,dayfirst=True, parse_dates=[20], infer_datetime_format=True)
        
        # TO READ EXCEL FILES
        elif file_type == 'excel':
            thresh_df_i = pd.read_excel(file_path, skiprows=5, na_values=0, index_col=0)
        
        # Deletes empty rows in dataframe
        thresh_df_i = thresh_df_i.dropna(how='all')

        thresh_df_i['Threshold voltage (mV)'] = file_threshold_mV

        # if first data file, initialise threshold dataframe
        if i == 0:
            thresh_df = thresh_df_i
        # else add data from this threshold to existing data frame
        else:
            thresh_df = pd.concat([thresh_df, thresh_df_i])

        print('Read file {0}'.format(file_path))

    return thresh_df

#==============================================================================
# Calculate decay factor for threshold data
#==============================================================================
def decay_factor(df, halflifeseconds, refdatetime, mass, dilution):
    #calculate time difference between start time and reference time
    df['Real as s']=pd.to_timedelta(df[' Real'],unit='s')
    #df['Finished g']=pd.to_datetime(df[' Finished'])
    df['timedif']=refdatetime-(df[' Finished']-df['Real as s'])
    #change to seconds
    df['timedif(s)']=df['timedif']/pd.Timedelta(seconds=1)
    #calculate decay factor
    df['Decay Factor']=(np.log(2)*df[' Real']/halflifeseconds/
             (1-np.exp(-np.log(2)*df[' Real']/halflifeseconds))*
             0.5**(df['timedif(s)']/halflifeseconds))

    #Constants to include in dataframe
    df['Reference datetime']=refdatetime
    df['Source mass mg']=mass
    df['Source dilution factor']=dilution
    print('calculated decay factor')
    return df

#==============================================================================
# Accidental coincidences
#==============================================================================
def accidental_coincidences(df, rt):
    # rt = resolving time IN NANOSECONDS
    # All the bits of the Venn diagram for calculating accidental coincidences

    ################# ASSUMING 3 CHANNELS
    df['p3A'] = (df[' A']-df[' AB']-df[' AC']+df[' ABC'])/df[' Live']
    df['p3B'] = (df[' B']-df[' AB']-df[' BC']+df[' ABC'])/df[' Live']
    df['p3C'] = (df[' C']-df[' BC']-df[' AC']+df[' ABC'])/df[' Live']
    
    df['p3S'] = df['p3A']+df['p3B']+df['p3C']

    df['p3AB'] = (df[' AB']-df[' ABC'])/df[' Live']
    df['p3AC'] = (df[' AC']-df[' ABC'])/df[' Live']
    df['p3BC'] = (df[' BC']-df[' ABC'])/df[' Live']

    df['p3D'] = df['p3AB'] +df['p3AC'] +df['p3BC']

    df['p3T'] = df[' ABC']/df[' Live']

    # Accidental coincidences calculation for 3 CHANNELS
    df['aAB'] = (2*(df['p3A']*df['p3B'] +df['p3A']*df['p3BC'] +df['p3B']*df['p3AC']
                    +df['p3AC']*df['p3BC']) 
                 +(df['p3S'] +df['p3D'] -df['p3AB'])*(df['p3AB'] +df['p3T']))*rt*10**(-9)

    df['aBC'] = (2*(df['p3B']*df['p3C'] +df['p3B']*df['p3AC'] +df['p3C']*df['p3AB']
                    +df['p3AB']*df['p3AC']) 
                 +(df['p3S'] +df['p3D'] -df['p3BC'])*(df['p3BC'] +df['p3T']))*rt*10**(-9)

    df['aAC'] = (2*(df['p3A']*df['p3C'] +df['p3A']*df['p3BC'] +df['p3C']*df['p3AB']
                    +df['p3AB']*df['p3BC']) 
                 +(df['p3S'] +df['p3D'] -df['p3AC'])*(df['p3AC'] +df['p3T']))*rt*10**(-9)

    df['aD'] = (2*(df['p3A']*df['p3B'] +df['p3B']*df['p3C'] +df['p3A']*df['p3C']) 
                 +df['p3S']*(df['p3D'] +df['p3T']))*rt*10**(-9)
    
    df['aT'] = (2*(df['p3A']*df['p3BC'] +df['p3B']*df['p3AC'] +df['p3C']*df['p3AB'])
                 +(df['p3S'] +df['p3D'])*df['p3T'] 
                 +2*(df['p3AB']*df['p3AC'] +df['p3AB']*df['p3BC'] + +df['p3AC']*df['p3BC']))*rt*10**(-9)

                    
    ############## ASSUMING 4 CHANNELS
    df['p4A'] = (df[' A']-df[' AB']-df[' AC']-df[' AX']
                +df[' ABX']+df[' ACX']+df[' ABC']-df[' ABCX'])/df[' Live']
    df['p4B'] = (df[' B']-df[' AB']-df[' BC']-df[' BX']
                +df[' ABX']+df[' BCX']+df[' ABC']-df[' ABCX'])/df[' Live']
    df['p4C'] = (df[' C']-df[' BC']-df[' AC']-df[' CX']
                +df[' BCX']+df[' ACX']+df[' ABC']-df[' ABCX'])/df[' Live']
    df['p4X'] = (df[' X']-df[' AX']-df[' BX']-df[' CX']
                +df[' ABX']+df[' ACX']+df[' BCX']-df[' ABCX'])/df[' Live']
    
    df['p4S'] = df['p4A']+df['p4B']+df['p4C']+df['p4X']

    df['p4AB'] = (df[' AB']-df[' ABX']-df[' ABC']+df[' ABCX'])/df[' Live']
    df['p4AC'] = (df[' AC']-df[' ACX']-df[' ABC']+df[' ABCX'])/df[' Live']
    df['p4BC'] = (df[' BC']-df[' BCX']-df[' ABC']+df[' ABCX'])/df[' Live']

    df['p4AX'] = (df[' AX']-df[' ABX']-df[' ACX']+df[' ABCX'])/df[' Live']
    df['p4BX'] = (df[' BX']-df[' ABX']-df[' BCX']+df[' ABCX'])/df[' Live']
    df['p4CX'] = (df[' CX']-df[' ACX']-df[' BCX']+df[' ABCX'])/df[' Live']

    df['p4D'] = df['p4AB'] +df['p4AC'] +df['p4BC'] +df['p4AX'] +df['p4BX'] +df['p4CX'] 

    df['p4ABX'] = (df[' ABX']-df[' ABCX'])/df[' Live']
    df['p4ACX'] = (df[' ACX']-df[' ABCX'])/df[' Live']
    df['p4BCX'] = (df[' BCX']-df[' ABCX'])/df[' Live']

    df['p4ABC'] = (df[' ABC']-df[' ABCX'])/df[' Live']

    df['p4ABCX'] = df[' ABCX']/df[' Live']

    # Accidental coincidences calculation for LDX types 4 CHANNELS
    df['aLDX Type1'] = ((df['p4X']+df['p4AX']+df['p4BX']+df['p4CX'])*(df['p4ABC']+df['p4AB']+df['p4AC']+df['p4BC'])
                         +df['p4A']*(df['p4BX']+df['p4CX']) +df['p4B']*(df['p4AX']+df['p4CX']) +df['p4C']*(df['p4AX']+df['p4BX'])
                         +df['p4AX']*df['p4BX'] +df['p4AX']*df['p4CX'] +df['p4BX']*df['p4CX'])*2*rt*10**(-9)

    df['aLDX Type2'] = (df['p4ABCX']*(df['p4S']+df['p4D']+df['p4ABC'])
                        +df['p4ABX']*(df['p4S']-df['p4C']+df['p4AB']+df['p4AX']+df['p4BX'])
                        +df['p4ACX']*(df['p4S']-df['p4B']+df['p4AC']+df['p4AX']+df['p4CX'])
                        +df['p4BCX']*(df['p4S']-df['p4A']+df['p4BC']+df['p4BX']+df['p4CX']))*rt*10**(-9)
    
    df['aLDX'] = df['aLDX Type1'] + df['aLDX Type2']
    print('Done accidental coincidences correction for threshold data')
    return df

#==============================================================================
# Calculate corrected rates for logical doubles (LD, LDX, X)
#==============================================================================
def doubles_rates_corrected(df,back_df, accidental_coincidence_corr = 0):
    # get the threshold voltages in the data file
    unique_thresholds = df['Threshold voltage (mV)'].unique()

    # initialise background dictionary
    background_dict = {}
    # go through background data one threshold at a time
    for threshold_i in unique_thresholds:
        # only look at backgrounds for that threshold
        back_df_i = back_df.loc[back_df['threshold'] == threshold_i]
        try:
            # add dictionary entry for that threshold background
            background_dict[str(threshold_i)] = [back_df_i['backLDr8'].iloc[0], 
                                                back_df_i['unc backLDr8'].iloc[0], 
                                                back_df_i['backXr8'].iloc[0], 
                                                back_df_i['unc backXr8'].iloc[0],
                                                back_df_i['backLDXr8'].iloc[0],
                                                back_df_i['unc backLDXr8'].iloc[0]]
        except KeyError:
            print('\n there is no background for the {0} mV threshold voltage in the background file \n'.format(threshold_i))
            break

    # then assign those backgrounds to the corresponding threshold data 
    # in the main dataframe 
    df["backLDr8"] = df["Threshold voltage (mV)"].apply(lambda x: background_dict.get(str(x))[0])
    df["unc_backLDr8"] = df["Threshold voltage (mV)"].apply(lambda x: background_dict.get(str(x))[1])
    df["backXr8"] = df["Threshold voltage (mV)"].apply(lambda x: background_dict.get(str(x))[2])
    df["unc_backXr8"] = df["Threshold voltage (mV)"].apply(lambda x: background_dict.get(str(x))[3])
    df["backLDXr8"] = df["Threshold voltage (mV)"].apply(lambda x: background_dict.get(str(x))[4])
    df["unc_backLDXr8"] = df["Threshold voltage (mV)"].apply(lambda x: background_dict.get(str(x))[5])

    if accidental_coincidence_corr == 1:
        # LD Rate, decay, accidental coinc and background corrected
        df['LDr8']=(df[' LD']/df[' Live'] -df['backLDr8'] -df['aD'])*df['Decay Factor']
        # X Rate, decay and background corrected
        df['Xr8']=(df[' X']/df[' Live'] -df['backXr8'])*df['Decay Factor']
        # LDX Rate, accidental coinc, decay and background corrected
        df['LDXr8']=(df[' LDX']/df[' Live'] -df['backLDXr8'] -df['aLDX'])*df['Decay Factor']
        print('done background, accidental coincidence and decay corrections')

    else: # don't correct for accidental coincidences
        # LD Rate, decay and background corrected
        df['LDr8']=(df[' LD']/df[' Live'] -df['backLDr8'])*df['Decay Factor']
        # X Rate, decay and background corrected
        df['Xr8']=(df[' X']/df[' Live'] -df['backXr8'])*df['Decay Factor']
        # LDX Rate, decay and background corrected
        df['LDXr8']=(df[' LDX']/df[' Live'] -df['backLDXr8'])*df['Decay Factor']
        print('done background and decay corrections')

    return df
  
#==============================================================================
# Calculate BG/C and associated stats
#============================================================================== 
def linearise_threshold_data(df, dilution, mass):
    # linearise the threshold data
    df['BG/C'] = df['LDr8']*df['Xr8']/df['LDXr8']*dilution/mass
    df['B'] = df['LDr8']*(dilution/mass)
    df['G/C-1'] = df['Xr8']/df['LDXr8'] - 1
    df['1-C/G'] = 1 - df['LDXr8']/df['Xr8']

    # uncertainty equations
    df['uncBG/C'] = df['BG/C']*np.sqrt((2*df[' LDX'])/(df[' LD']*df[' X']) - 1/df[' X'] - 1/df[' LD'] + 1/df[' LDX'])
    df['uncB'] = df['B']*np.sqrt((2*df[' LDX'])/(df[' LD']*df[' X']) - 1/df[' X'] - 1/df[' LD'] + 1/df[' LDX'])
    df['uncG/C-1'] = df['G/C-1']*np.sqrt(df[' X']/(df[' LDX']*(df[' X'] - df[' LDX'])))
    df['unc1-C/G'] = df['1-C/G']*np.sqrt(df[' LDX']/(df[' X']*(df[' X'] - df[' LDX'])))

    print('Linearised threshold data')
    return df

#==============================================================================
# STATSGET 
# Generates values for use in regressions, and calculates a variety of 
# statistics...
# For G/C-1 and 1-C/G ('independent variables'): weighted mean, 
# theoretical standard deviation of weighted mean.
# 
# For BG/C and B ('Dependent variables'): weighted mean, 
# theoretical standard deviation of weighted mean, 
# theoretical standard deviation of weighted mean as a relative %, 
# observed variance in the weighted mean, 
# observed standard deviation in weighted mean, 
# and the observed standard deviation of the weighted mean as a relative % 
# mean
# 
#==============================================================================

def stats_get(thresh_df, weight_mean_bool):
    # this function produces the regression data 
    # i.e. takes weighted means and calculates st devs for data from each threshold, 
    # puts stats from all thresholds into a dataframe called reg_df

    # weight_mean_bool is a boolean 
    # i.e. if weight_mean_bool == 1, calculate weighted mean

    # get the threshold voltages in the data file
    unique_thresholds = thresh_df['Threshold voltage (mV)'].unique()

    # names of columns in the eventual regression dataframe
    regcolumns=['G/C-1 WM', 'Stdev of Mean G/C-1 (Theor)', 'Stdev of Mean G/C-1 (Obs)',
                'BG/C WM', 'Stdev of Mean BGC (Theor)', 'Stdev of Mean BGC (Obs)',
                '1-C/G WM', 'Stdev of Mean 1-C/G (Theor)', 'Stdev of Mean 1-C/G (Obs)', 
                'B WM', 'Stdev of Mean B (Theor)', 
                'Stdev of Mean B (Obs)', 'Stdev of Mean B (Obs)(ODR only)']

    # go through data one threshold at a time
    for i in range(len(unique_thresholds)):
        threshold_i = unique_thresholds[i]
        # only look at data for that threshold
        df = thresh_df.loc[thresh_df['Threshold voltage (mV)'] == threshold_i]

        # BG/C Freda's equations
        wmBGC_CT, weightBGC_CT = weighted_mean_calc(df['uncBG/C'], df['BG/C'], weight_mean_bool) 
        stdev_theorBGC = stdev_mean_theor_calc(weightBGC_CT)
        stdev_obsBGC = wmBGC_CT*np.sqrt(np.nanvar(df['LDr8'], ddof=1)/np.nanmean(df['LDr8'])**2 + np.nanvar(df['Xr8'], ddof=1)/np.nanmean(df['Xr8'])**2 + np.nanvar(df['LDXr8'], ddof=1)/np.nanmean(df['LDXr8'])**2 
                                        + (2*np.cov(df['LDr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDr8'])*np.nanmean(df['Xr8'])) 
                                            - (2*np.cov(df['LDXr8'],df['LDr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['LDr8']))  
                                                - (2*np.cov(df['LDXr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['Xr8'])))/np.sqrt(np.count_nonzero(~np.isnan(df['LDr8'])) - 1)

        # B Freda's equations
        wmB_CT_Freda, weightB_CT_Freda = weighted_mean_calc(df['uncB'], df['B'], weight_mean_bool) 
        stdev_theorB_Freda = stdev_mean_theor_calc(weightB_CT_Freda)
        stdev_obsB_Freda = wmB_CT_Freda*np.sqrt(np.nanvar(df['LDr8'], ddof=1)/np.nanmean(df['LDr8'])**2 + np.nanvar(df['Xr8'], ddof=1)/np.nanmean(df['Xr8'])**2 + np.nanvar(df['LDXr8'], ddof=1)/np.nanmean(df['LDXr8'])**2 
                                        + (2*np.cov(df['LDr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDr8'])*np.nanmean(df['Xr8'])) 
                                            - (2*np.cov(df['LDXr8'],df['LDr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['LDr8']))  
                                                - (2*np.cov(df['LDXr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['Xr8'])))/np.sqrt(np.count_nonzero(~np.isnan(df['LDr8'])) - 1)
        stdev_obsB_ODRonly = wmB_CT_Freda*np.sqrt(np.nanvar(df['LDr8'], ddof=1)/np.nanmean(df['LDr8'])**2)/np.sqrt(np.count_nonzero(~np.isnan(df['LDr8'])) - 1)

        # G/C - 1 Freda's equations
        wmGC_1_CT, weightGC_1_CT = weighted_mean_calc(df['uncG/C-1'], df['G/C-1'], weight_mean_bool) 
        stdev_theorGC_1 = stdev_mean_theor_calc(weightGC_1_CT)
        stdev_obsGC_1 = (wmGC_1_CT+1)*np.sqrt(np.nanvar(df['LDXr8'], ddof=1)/np.nanmean(df['LDXr8'])**2 + np.nanvar(df['Xr8'], ddof=1)/np.nanmean(df['Xr8'])**2
                        - (2*np.cov(df['LDXr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['Xr8'])))/np.sqrt(np.count_nonzero(~np.isnan(df['LDXr8'])) - 1)  

        # 1 - C/G Freda's equations
        wm1_CG_CT, weight1_CG_CT = weighted_mean_calc(df['unc1-C/G'], df['1-C/G'], weight_mean_bool) 
        stdev_theor1_CG = stdev_mean_theor_calc(weight1_CG_CT)
        stdev_obs1_CG = (1-wm1_CG_CT)*np.sqrt(np.nanvar(df['LDXr8'], ddof=1)/np.nanmean(df['LDXr8'])**2 + np.nanvar(df['Xr8'], ddof=1)/np.nanmean(df['Xr8'])**2
                        - (2*np.cov(df['LDXr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['Xr8'])))/np.sqrt(np.count_nonzero(~np.isnan(df['LDXr8'])) - 1)    

        threshold_data_stats = np.array([[wmGC_1_CT, stdev_theorGC_1, stdev_obsGC_1, wmBGC_CT, stdev_theorBGC, stdev_obsBGC, 
                                         wm1_CG_CT, stdev_theor1_CG, stdev_obs1_CG, wmB_CT_Freda, stdev_theorB_Freda, stdev_obsB_Freda, stdev_obsB_ODRonly]])
        reg_df_i = pd.DataFrame(data=threshold_data_stats,index=[i],columns=regcolumns)

        if i == 0:
            # initialise regression dataframe with calculated statistics from first threshold data
            reg_df = reg_df_i
        else:
            # append calculated statistics from next threshold data to regression dataframe
            reg_df = pd.concat([reg_df, reg_df_i])

    print('Finished calculating data for regressions. Calculated weighted means and st devs for each threshold')
    return reg_df
             