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

# FONT SIZE (not for residual plots, only for extrapolation plots)
fsize=10 # 10 recommended for papers, 12 recommended for powerpoints
# PLOT SIZE (not for residual plots, only for extrapolation plots)
pwidth=6 # 6 recommended for papers, 8 recommended for powerpoints
pheight=4 # 4 recommended for papers and powerpoints

type_of_fit = 'cubic' # 'linear' OR 'cubic'
type_of_reg = 'ODR' # 'LS' OR 'ODR'
what_to_plot = 'B' # 'BGC' for BG/C vs G/C-1, OR 'B' for B vs 1-C/G 
x_max = 0.1 # x axis maximum for plot
#==============================================================================
""" end of user input """
#==============================================================================

PlotStuff4LSC.plotter(regression_data_filename, fits_filename, what_to_plot, 
                      type_of_fit, type_of_reg, x_max, fsize, pwidth, pheight)

print('finished plotting stuff')