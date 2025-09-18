# -*- coding: utf-8 -*-
"""
Siobhan Test February 2018 on Cu64 D3 LS3 data

@author: siobhant
"""
import os
import pandas as pd
import numpy as np
import datetime
from scipy.optimize import curve_fit
import scipy.odr.odrpack as odear
import sys
import pandas4LSC_data_processing as LSC_data_processing
import pandas4LSC_data_analysis as LSC_data_analysis

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" YOUR INPUT IS REQUIRED BELOW """
#==============================================================================
""" CONSTANTS """
#==============================================================================
# set reference date yyyy,mm,dd,hh,MM,SS,fractions of second. leave out
# leading zeroes, e.g. June is 6 not 06
refdatetime=datetime.datetime(2024,12,15,23,0,0,0)
# Half life in seconds
halflifeseconds=1133*365.24219878*24*60*60
# Dilution factor
dilution=4.43253556514143
# Active solution mass in mg
mass=55.4836969556208 # LS2
# Branch ratio
branchratio=1
# What do you want the output files to be called? 
# Source name is generally good
outputfilename='Ho166m-LS2-win2_siobhan_AC_test'
# Resolving time in ns
rt=250
# Dead time in us
dt=50
# Gamma shift
gs=-5050
#==============================================================================
""" OPTIONS """
#==============================================================================
# Doubles or triples? Select one at a time
doubles_data = 1 # 0=NO 1=YES
triples_data = 0 # 0=NO 1=YES
# Manually specified background?
specified_background = 0 # 0=NO 1=YES
# Accidental coincidences correction?
accidental_coinc = 1 # 0=NO 1=YES
# weighting of points? Select one at a time
standard_deviation = 'OSD' # observed standard deviation or 'TSD' for theoretical standard deviation
# Use weighted mean?
weight_mean = 1 # 0=NO 1=YES

#==============================================================================
""" DATA """
#==============================================================================
# LIST OF DATA FILES
# DIRECTORY, where are the data files located?
data_dir="csvs_LS2_win2"
# format is [filename, threshold_voltage_in_mV]
data_filename_list = [["LS2_win2_20mV.xlsx", 20],
                      ["LS2_win2_50mV.xlsx", 50],
                      ["LS2_win2_100mV.xlsx", 100],
                      ["LS2_win2_200mV.xlsx", 200],
                      ["LS2_win2_300mV.xlsx", 300],
                      ["LS2_win2_400mV.xlsx", 400],
                      ["LS2_win2_500mV.xlsx", 500],
                      ["LS2_win2_600mV.xlsx", 600],
                      ["LS2_win2_700mV.xlsx", 700],
                      ["LS2_win2_800mV.xlsx", 800],
                      ["LS2_win2_900mV.xlsx", 900]
                     ]

#==============================================================================
""" If you are specifying backgrounds with a single excel spreadsheet """
# i.e. specfied_background = 1 
""" enter the name of the spreadsheet here... """
# DIRECTORY, where is the specified background file located?
specified_background_dir="csvs_LS2_win2"
# BACKGROUND EXCEL FILE
specified_background_excel="backgroundall_Doubles.xlsx"
#==============================================================================
""" Or if you are just analysing csv files to get the backgrounds for
each threshold, put those files here... """
# i.e. specfied_background = 0
# BACKGROUND XLSX FILES
# DIRECTORY, where are the data files located?
background_dir="csvs_BKG_win2"
# format is [filename, threshold_voltage_in_mV]
background_filename_list = [["BKG-2_20mV_com.xlsx", 20],
                            ["BKG-2_50mV_com.xlsx", 50],
                            ["BKG-2_100mV_com.xlsx", 100],
                            ["BKG-2_200mV_com.xlsx", 200],
                            ["BKG-2_300mV_com.xlsx", 300],
                            ["BKG-2_400mV_com.xlsx", 400],
                            ["BKG-2_500mV_com.xlsx", 500],
                            ["BKG-2_600mV.xlsx", 600],
                            ["BKG-2_700mV.xlsx", 700],
                            ["BKG-2_800mV.xlsx", 800],
                            ["BKG-2_900mV.xlsx", 900]
                            ]

#==============================================================================
#~~~~~~~~~~~~~~~~ NO MORE USER INPUT REQUIRED UNLESS PROMPTED ~~~~~~~~~~~~~~~~~
if accidental_coinc == 1:
    ACcorr = '_ACcorr'
else:
    ACcorr = ''
############### ANALYSING BACKGROUND DATA
# Load background data
if specified_background == 1:
    back_file_path = specified_background_dir + '/' + specified_background_excel
    back_df = pd.read_excel(back_file_path, skiprows=0, na_values=0, index_col=0)
# or load raw background data and process
elif specified_background == 0:
    raw_background_df = LSC_data_processing.get_data(background_dir, background_filename_list, 
                                                     file_type = 'excel')
    if accidental_coinc == 1:
        raw_background_df = LSC_data_processing.accidental_coincidences(raw_background_df, rt)
    
    raw_background_df = LSC_data_processing.background_doubles_rates_corrected(raw_background_df, 
                                                                               accidental_coincidence_corr = accidental_coinc)
    all_background_data_filename = "Background_Thresh_Data_rt{0}dt{1}_{2}{3}.xlsx".format(rt,dt,'doubles',ACcorr)
    raw_background_df.to_excel(all_background_data_filename)

    back_df = LSC_data_processing.background_average(raw_background_df)
    averaged_background_data_filename = "Background_Averaged_rt{0}dt{1}_{2}{3}.xlsx".format(rt,dt,'doubles',ACcorr)
    back_df.to_excel(averaged_background_data_filename)
   
    print('\n Processed background data and saved averaged background to {0}'.format(averaged_background_data_filename))
    print(back_df)
    print()

############### ANALYSING THRESHOLD DATA
# Load threshold data
thresh_df = LSC_data_processing.get_data(data_dir, data_filename_list, 
                                         file_type = 'excel')
# calculate decay factor
thresh_df = LSC_data_processing.decay_factor(thresh_df, halflifeseconds, 
                                             refdatetime, mass, dilution)
# calculate accidental coincidences
if accidental_coinc == 1:
    threshdf = LSC_data_processing.accidental_coincidences(thresh_df, rt)

# apply corrections (background, decay factor, accidental coinc if requested)
thresh_df = LSC_data_processing.doubles_rates_corrected(thresh_df,back_df,
                                                        accidental_coincidence_corr = accidental_coinc)

# linearise threshold data
thresh_df = LSC_data_processing.linearise_threshold_data(thresh_df, dilution, mass)

# Set up filename suffix
if doubles_data==1:
    DorT='doubles'
    if weight_mean == 1:
        WM="_WM"
    else:
        WM="_unWM"
    if specified_background==1:
        Sb='_SB'
    else:
        Sb=''

# save threshold dataframe to excel spreadsheet
threshold_data_filename = "{0}_rt{1}dt{2}_ThreshData_{3}{4}{5}{6}_siobhan2025.xlsx".format(outputfilename,rt,dt,DorT,Sb,WM,ACcorr)
threshdf.to_excel(threshold_data_filename)
print('written threshold data to {0}'.format(threshold_data_filename))

# calculate stats from each threshold to produce regression dataframe
regression_df = LSC_data_processing.stats_get(thresh_df, weight_mean)
reg_filename="{0}_rt{1}dt{2}_RegData_{3}{4}{5}{6}_newunceqns.xlsx".format(outputfilename,rt,dt,DorT,Sb,WM,ACcorr)
regression_df.to_excel(reg_filename)
print('written regression data to {0}'.format(reg_filename))
        
# # Or, directly load in the data to use for the for regression
# regdf=pd.read_excel(my_regression_data_filename,header=0,index_col=0)

# # Which regression code to run???
# do regression on regression dataframe
LSC_data_analysis.regression(regression_df, outputfilename, branchratio, 
                             refdatetime, halflifeseconds, dilution, mass, rt, 
                             dt, gs, Sb, DorT, ACcorr, WM,
                             StandDev = standard_deviation)

 

# print()
# print()
# if specified_background==1:
#     print("Background specified by user")
# print("{0} data used. Points weighted by {1}".format(DorT,StandDev))
print("Finished python4LSC. See you later :) ")
