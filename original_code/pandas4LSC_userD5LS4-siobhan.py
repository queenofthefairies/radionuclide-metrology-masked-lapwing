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
mass=65.3471827893177
# Branch ratio
branchratio=1
# What do you want the output files to be called? 
# Source name is generally good
outputfilename='Ho166m-LS7-win2'
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
# weighting of points? Select one at a time
theoretical_sd = 1 # 0=NO 1=YES
observed_sd = 0 # 0=NO 1=YES
# Use weighted mean?
weight_mean = 1 # 0=NO 1=YES
# DIRECTORY, where are the csv files located?
rtdt1dir="csvs_LS7_win2"
#==============================================================================
""" DO YOU WANT TO USE THE SPEEDY FILE GRABBER? 0=NO, 1=YES """
speedyfilegrabber=0
bkgprefix=0 #for using speedy file grabber section
dataprefix=0 #for using speedy file grabber section
#==============================================================================
""" ONLY USE THIS SECTION IF NOT USING SPEEDY FILE GRABBER """
# FILES
thresh1A="LS7_20mV.xlsx"
thresh1B="LS7_50mV.xlsx"
thresh1C="LS7_100mV.xlsx"
thresh1D="LS7_200mV.xlsx"
thresh1E="LS7_300mV.xlsx"
thresh1F="LS7_400mV.xlsx"
thresh1G="LS7_500mV.xlsx"
thresh1H="LS7_600mV.xlsx"
thresh1I="LS7_700mV.xlsx"
thresh1J="LS7_800mV.xlsx"
thresh1K="LS7_900mV.xlsx"
thresh1L="LS7_20mV_r.xlsx"
thresh1M="LS7_50mV_r.xlsx"
thresh1N="LS7_100mV_r.xlsx"
thresh1O="LS7_200mV_r.xlsx"
thresh1P="LS7_300mV_r.xlsx"
thresh1Q="LS7_400mV_r.xlsx"
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

# Which main program to run???
if doubles_data==1:
    DorT='doubles'
    if weight_mean == 1:
        WM="_WM"
    else:
        WM="_unWM"
    if specified_background==1:
        Sb='_SB'
        print("Calling pandas4LSCdoubles_main_specifiedbkg") 
        import pandas4LSCdoubles_main_specifiedbkg
    else:
        Sb=''
        print("Calling pandas4LSCdoubles_main") 
        import pandas4LSCdoubles_main
        
if triples_data==1:
    DorT='triples'
    if weight_mean == 1:
        WM="_WM"
    else:
        WM="_unWM"
    if specified_background==1:
        Sb='_SB'
        print("Calling pandas4LSCtriples_main_specifiedbkg")
        import pandas4LSCtriples_main_specifiedbkg
    else:
        Sb=''
        print("Calling pandas4LSCtriples_main") 
        import pandas4LSCtriples_main
        
# Data for regression
filename="{0}_rt{1}dt{2}_RegData_{3}{4}{5}_newunceqns.xlsx".format(outputfilename,rt,dt,DorT,Sb,WM)
regdf=pd.read_excel(filename,header=0,index_col=0)

# Which regression code to run???

if theoretical_sd==1:
    StandDev='TSD'
    print("Calling pandas4LSC_regression")
    import pandas4LSC_regression  
    
if observed_sd==1:
    StandDev='OSD'
    print("Calling pandas4LSC_regression")
    import pandas4LSC_regression    

 

print()
print()
if specified_background==1:
    print("Background specified by user")
print("{0} data used. Points weighted by {1}".format(DorT,StandDev))
print("Finished python4LSC. See you later :) ")
