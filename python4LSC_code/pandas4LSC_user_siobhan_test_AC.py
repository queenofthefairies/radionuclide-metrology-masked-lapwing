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
import pandas4LSC_threshold_data_processing as LSC_data_processing
#import pandas4LSC_regression_new as LSC_regression
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
# Doubles or triples? Select one at a time
doubles_data = 1 # 0=NO 1=YES
triples_data = 0 # 0=NO 1=YES
# Manually specified background?
specified_background = 1 # 0=NO 1=YES
# Accidental coincidences correction?
accidental_coinc = 1 # 0=NO 1=YES
# weighting of points? Select one at a time
theoretical_sd = 1 # 0=NO 1=YES
observed_sd = 0 # 0=NO 1=YES
# Use weighted mean?
weight_mean = 1 # 0=NO 1=YES
# DIRECTORY, where are the data files located?
data_dir="csvs_LS2_win2"

#==============================================================================
# FILES
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
# BACKGROUND EXCEL FILE
backgroundexcel="backgroundall_Doubles.xlsx"
#==============================================================================
""" Or if you are just analysing csv files to get the backgrounds for
each threshold, put those files here... """
# i.e. specfied_background = 0
# BACKGROUND CSV FILES
# back1A must be the background for thresh1A etc
back1A="BKG-2_20mV_com.xlsx"
back1B="BKG-2_50mV_com.xlsx"
back1C="BKG-2_100mV_com.xlsx"
back1D="BKG-2_200mV_com.xlsx"
back1E="BKG-2_300mV_com.xlsx"
back1F="BKG-2_400mV_com.xlsx"
back1G="BKG-2_500mV_com.xlsx"
back1H="BKG-2_600mV.xlsx"
back1I="BKG-2_700mV.xlsx"
back1J="BKG-2_800mV.xlsx"
back1K="BKG-2_900mV.xlsx"
back1L=0
back1M=0
back1N=0
back1O=0
back1P=0
back1Q=0

#==============================================================================
#~~~~~~~~~~~~~~~~ NO MORE USER INPUT REQUIRED UNLESS PROMPTED ~~~~~~~~~~~~~~~~~

# Filename suffix
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

# Load threshold data
thresh_df = LSC_data_processing.get_data(data_dir, data_filename_list, file_type = 'excel')
# calculate decay factor
thresh_df = LSC_data_processing.decay_factor(thresh_df, halflifeseconds, refdatetime, mass, dilution)
        
print(thresh_df)
        
# # Data for regression
# filename="{0}_rt{1}dt{2}_RegData_{3}{4}{5}_newunceqns.xlsx".format(outputfilename,rt,dt,DorT,Sb,WM)
# regdf=pd.read_excel(filename,header=0,index_col=0)

# # Which regression code to run???

# if theoretical_sd==1:
#     StandDev='TSD'
#     print("Calling pandas4LSC_regression")
#     import pandas4LSC_regression  
    
# if observed_sd==1:
#     StandDev='OSD'
#     print("Calling pandas4LSC_regression")
#     import pandas4LSC_regression    

 

# print()
# print()
# if specified_background==1:
#     print("Background specified by user")
# print("{0} data used. Points weighted by {1}".format(DorT,StandDev))
print("Finished python4LSC. See you later :) ")
