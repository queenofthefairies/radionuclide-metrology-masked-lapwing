# -*- coding: utf-8 -*-
"""
PLOTSTUFF4LSC USER

This script takes user input:
- regression data file 
- all fits file

@author: siobhant
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import scipy.odr.odrpack as odear
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

import PlotStuff4LSC

#==============================================================================
""" Start of user input """
#==============================================================================

regression_data_filename='Ho166m-LS2-win2_siobhan_AC_test_rt250dt50_RegData_doubles_SB_WM_ACcorr_newunceqns.xlsx'
fits_filename="Ho166m-LS2-win2_siobhan_AC_test_rt250dt50_AllFits_doubles_SB_OSD_WM_ACcorr_newunceqns.xlsx"
sourcename="LS2_doubles"
branchingratio=1
DorT='Doubles'
SD='TSD' # or 'TSD' for weighting by theoretical standard deviation

#CHANGE DOMAIN (plot horizontal axis maximum)
xMaxBGC=.1 #5
#CHANGE RANGE (plot vertical axis maximum)
yMinBGC=1030
yMaxBGC=1050 #50000
#change domain of resid plot
xMinRBGC=0
xMaxRBGC=.1 #5
#change range resid plot
yMinRBGC=-2
yMaxRBGC=2

#CHANGE DOMAIN (plot horizontal axis maximum)
xMaxB=0.1
#CHANGE RANGE (plot vertical axis maximum)
yMinB=950 #6000
yMaxB=1040
#change domain of resid plot
xMinRB=0
xMaxRB=.1 #0.9
#change range resid plot
yMinRB=-2
yMaxRB=2

# FONT SIZE (not for residual plots, only for extrapolation plots)
fsize=10 # 10 recommended for papers, 12 recommended for powerpoints
# PLOT SIZE (not for residual plots, only for extrapolation plots)
pwidth=6 # 6 recommended for papers, 8 recommended for powerpoints
pheight=4 # 4 recommended for papers and powerpoints

# POSITION OF LEGEND AND TEXT ON PLOTS
# 'upper right'  : 1
#'upper left'   : 2
#'lower left'   : 3
#'lower right'  : 4
# BGC vs GC-1
BGC_legpos = 4 # 4 is default
BGC_txtpos = 2 # 2 is default
# B vs 1-CG
B_legpos = 3 # 3 is default
B_txtpos = 1 # 1 is default


type_of_fit = 'linear'
type_of_reg = 'LS'
stand_dev = 'OSD'
what_to_plot = 'BGC'
#==============================================================================
""" end of user input """
#==============================================================================
#datfileref = regression_data_filename[:-5]
print()

#reg_df = pd.read_excel(regression_data_filename, index_col=0)
#reg_df.head()

#datatofit = reg_df.sort_values('1-C/G WM')

PlotStuff4LSC.plotter(regression_data_filename, fits_filename, what_to_plot, 
                      type_of_fit, type_of_reg, stand_dev)


  
    
# print()
# print("Plotting the least squares fits of B vs 1-C/G")
# import PlotStuff4LSC_1


# print()
# print("Plotting the least squares fits of BG/C vs G/C-1")
# import PlotStuff4LSC_2


# print()
# print("Plotting the ODR fits of B vs 1-C/G")
# import PlotStuff4LSC_3

# print()
# print("Plotting the ODR fits of BG/C vs G/C-1")
# import PlotStuff4LSC_4

# print()
# print("Finished plotting stuff!")

# plt.close("all")