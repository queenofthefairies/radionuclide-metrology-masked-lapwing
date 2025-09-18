# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 20:49:22 2017

@author: hepvis
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import scipy.odr.odrpack as odear
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText


def f1(x, m, c):
    """ The linear function y= m*x + c """
    return m*x + c

def f3(x, a, c, d):
    """ The cubic function y= a*x^3 + c*x + d """
    return a*x**3  + c*x + d

def tick_function(X, what_to_plot):
    if what_to_plot == 'B':
        V = (1-X)*100
    elif what_to_plot == 'BGC':
        V = 1/(1+X)*100
    return ["%.1f" % z for z in V]

def plotter(reg_data_filename, fits_data_filename, what_to_plot, 
            type_of_fit, type_of_reg, x_max, fsize, pwidth, pheight):

    # get summary about source, stand dev etc
    summary_df = pd.read_excel(fits_data_filename, 
                               sheet_name='Summary')
    source_name_df = summary_df.loc[(summary_df['Function'] == 'Source name')]
    source_name = source_name_df['Regression data'][9]
    std_dev_df = summary_df.loc[(summary_df['Function'] == 'Standard deviation used as weights')]
    stand_dev = std_dev_df['Regression data'][18]
    backg_df = summary_df.loc[(summary_df['Function'] == 'Background')]
    backg = backg_df['Regression data'][19]
    if backg != '_SB':
        backg = ''


    # read in regression data
    reg_df = pd.read_excel(reg_data_filename, index_col=0)

    if what_to_plot == 'BGC':
        x_data = np.array(reg_df['G/C-1 WM'])
        y_data = np.array(reg_df['BG/C WM'])
        if stand_dev == 'OSD':
            x_unc = np.array(reg_df['Stdev of Mean G/C-1 (Obs)'])
            y_unc = np.array(reg_df['Stdev of Mean BGC (Obs)'])
        elif stand_dev == 'TSD':
            x_unc = np.array(reg_df['Stdev of Mean G/C-1 (Theor)'])
            y_unc = np.array(reg_df['Stdev of Mean BGC (Theor)'])

    elif what_to_plot == 'B':
        x_data = np.array(reg_df['1-C/G WM'])
        y_data = np.array(reg_df['B WM'])
        if stand_dev == 'OSD':
            x_unc = np.array(reg_df['Stdev of Mean 1-C/G (Obs)'])
            y_unc = np.array(reg_df['Stdev of Mean B (Obs)'])
            if type_of_fit == 'ODR':
                y_unc = np.array(reg_df['Stdev of Mean B (Obs)(ODR only)'])
        elif stand_dev == 'TSD':
            x_unc = np.array(reg_df['Stdev of Mean 1-C/G (Theor)'])
            y_unc = np.array(reg_df['Stdev of Mean B (Theor)'])
    print('read in regression data from {0}'.format(reg_data_filename))
    
    # read in fits data
    sheet_dictionary = {'linear' : 'LinearFits', 'cubic': 'CubicFits'}
    cols_range_dictionary = {'linear' : range(1,9), 'cubic': range(1,11)}
    fits_df = pd.read_excel(fits_data_filename, index_col=None, usecols=cols_range_dictionary[type_of_fit], 
                            sheet_name=sheet_dictionary[type_of_fit])
    fitted_data_dictionary = {'BGC':'BG/C vs G/C-1', 'B': 'B vs 1-C/G'}
    fit_method_dictionary = {'LS':'Least Squares', 'ODR':'Orthogonal Distance'}
    specific_fit_df = fits_df.loc[(fits_df['Regression data'] == fitted_data_dictionary[what_to_plot]) & (fits_df['Regression method'] == fit_method_dictionary[type_of_reg])]
    print('read in fit information from {0}'.format(fits_data_filename))
    print()
    print('Plotting this specific {0} fit'.format(type_of_fit))
    print(specific_fit_df)
    print()
    
    fit_params = []
    if type_of_fit == 'linear':
        fit_gradient = specific_fit_df['gradient'].to_numpy()[0]
        fit_intercept = specific_fit_df['intercept'].to_numpy()[0]
        fit_params = [fit_gradient, fit_intercept]
    elif type_of_fit == 'cubic':
        fit_a = specific_fit_df['a'].to_numpy()[0]
        fit_c = specific_fit_df['c'].to_numpy()[0]
        fit_intercept = specific_fit_df['intercept'].to_numpy()[0]
        fit_params = [fit_a,fit_c,fit_intercept]
    activity_conc = specific_fit_df['Activity conc'].to_numpy()[0]
    activity_conc_rel_unc = specific_fit_df['rel unc intercept %'].to_numpy()[0]
    activity_conc_abs_unc = activity_conc*activity_conc_rel_unc/100

    fit_x_to_plot = np.linspace(0, x_max, 200)
    if type_of_fit == 'linear':
        fit_y_to_plot = f1(fit_x_to_plot, *fit_params)
        residuals_to_plot = y_data - f1(x_data, *fit_params)
    elif type_of_fit == 'cubic':
        fit_y_to_plot = f3(fit_x_to_plot, *fit_params)
        residuals_to_plot = y_data - f3(x_data, *fit_params)
    
    if what_to_plot == 'B':
        x_label = r"$1 - N_C /N_\gamma$"
        y_label = r"$N_\beta $ (s$^{-1}$ mg $^{-1}$)"
        leg_pos = 3 # legend position
        txt_pos = 1 # corner text position
        plot_fig_filename = '{0}_B_vs_ineff_{1}_{2}_PLOT_{2}{3}.png'.format(source_name, type_of_fit, type_of_reg, stand_dev,backg)
        resid_plot_fig_filename = '{0}_B_vs_ineff_{1}_{2}_RESIDUALS_PLOT_{3}{4}.png'.format(source_name,type_of_fit,type_of_reg,stand_dev,backg)

    elif what_to_plot == 'BGC':
        x_label = r"$N_\gamma / N_C - 1$"
        y_label = r"$N_\beta N_\gamma / N_C $ (s$^{-1}$ mg $^{-1}$)"
        leg_pos = 4 # legend position
        txt_pos = 2 # corner text position
        plot_fig_filename = '{0}_BGC_vs_GoC-1_{1}_{2}_PLOT_{3}{4}.png'.format(source_name,type_of_fit,type_of_reg,stand_dev,backg)
        resid_plot_fig_filename = '{0}_BGC_vs_GoC-1_{1}_{2}_RESIDUALS_PLOT_{3}{4}.png'.format(source_name,type_of_fit,type_of_reg,stand_dev,backg)
    
    cornertext="Activity from {0} intercept: \n".format(type_of_fit) +r"{0:.4f} $\pm$ {1:.4f} MBq/g ({2:.4f}%)".format(activity_conc,
                                                                                                                       activity_conc_abs_unc,
                                                                                                                       activity_conc_rel_unc)
    # MAKE PLOT   
    f, ax = plt.subplots(1,1)
    # plot regression data
    ax.errorbar(x_data, y_data, xerr=x_unc, yerr=y_unc, fmt= 'o', capsize=2, label='{0}'.format(source_name))
    # plot fit
    ax.plot(fit_x_to_plot, fit_y_to_plot, label='{0} fit {1}'.format(type_of_fit,type_of_reg))
    # tweak plot layout
    plt.locator_params(axis='y', nbins=6)
    plt.locator_params(axis='x', nbins=6)
    plt.xticks(fontsize=fsize-1, rotation=0)
    plt.yticks(fontsize=fsize-1, rotation=0)
    plt.xlim(0,x_max)
    ax.legend(loc=leg_pos, fontsize=fsize-1)
    # add corner text
    ax.add_artist(AnchoredText(cornertext,prop=dict(size=fsize-1),loc=txt_pos,frameon=False,borderpad=1))
    # add axes labels
    ax.set_xlabel(x_label,fontsize=fsize)
    ax.set_ylabel(y_label,fontsize=fsize)

    # plot secondary x-axis to show beta efficiency
    ax2 = ax.twiny()
    new_tick_locations = np.linspace(0, x_max, 6)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations,what_to_plot))
    ax2.set_xlabel(r"$\beta$ efficiency %",fontsize=fsize)
    plt.locator_params(axis='y', nbins=6)
    plt.locator_params(axis='x', nbins=6)
    plt.xticks(fontsize=fsize-1, rotation=0)
    plt.yticks(fontsize=fsize, rotation=0)
    plt.xlim(0,x_max)

    # final plot
    plt.tight_layout()
    plt.gcf().set_size_inches(pwidth,pheight)
    plt.savefig(plot_fig_filename,dpi=240)
    plt.show()

    print('done plotting data. plot saved to {0}'.format(plot_fig_filename))
    print()

    # MAKE RESIDUAL PLOT
    f, ax = plt.subplots(1,1)
    # plot zero line
    ax.axhline(y=0,linewidth=1,color='C1')
    # plot residuals
    ax.errorbar(x_data, residuals_to_plot, xerr=x_unc, yerr=y_unc, fmt= 'o', capsize=2, label='{0}'.format(source_name))
    # tweak plot layout
    plt.locator_params(axis='y', nbins=6)
    plt.locator_params(axis='x', nbins=6)
    plt.xticks(fontsize=fsize-1, rotation=0)
    plt.yticks(fontsize=fsize-1, rotation=0)
    plt.xlim(0,x_max)
    ax.legend(loc=leg_pos, fontsize=fsize-1)
    # add axes labels
    ax.set_xlabel(x_label,fontsize=fsize)
    ax.set_ylabel('Residuals ' + y_label,fontsize=fsize)
    # add plot title
    ax.set_title('{0} {1} {2} residuals'.format(source_name, type_of_fit, type_of_reg))

    # final plot
    plt.tight_layout()
    plt.gcf().set_size_inches(pwidth,pheight)
    plt.savefig(resid_plot_fig_filename,dpi=240)
    plt.show()


    # f, (ax1,ax3) = plt.subplots(2, sharex=True, sharey=True)
    # ax1.errorbar(x_dataB, residLT, xerr=0, yerr=y_uncB, fmt= 'o', capsize=2, color='#ff7f0e', label='Linear fit residuals')
    # ax1.axhline(y=0,color='#1f77b4',linewidth=1)
    # ax1.legend(loc='lower right', fontsize=10)
    # ax1.set_title('{0} LS residuals'.format(sourcename))
    # ax1.set_ylabel(r"$N_\beta $ (s$^{-1}$ mg $^{-1}$)")
    # ax3.errorbar(x_dataB, residCT, xerr=0, yerr=y_uncB, fmt= 'o', capsize=2, color='#2ca02c',label='Cubic fit residuals')
    # ax3.axhline(y=0,color='#1f77b4',linewidth=1)
    # ax3.legend(loc='lower right',fontsize=10)
    # ax3.set_xlabel(r"$1 - N_C / N_\gamma $")
    # ax3.set_ylabel(r"$N_\beta $ (s$^{-1}$ mg $^{-1}$)")
    # f.subplots_adjust(hspace=0)
    # plt.xlim(xMinRB,xMaxRB)
    # plt.ylim(yMinRB,yMaxRB)
    # plt.tight_layout()
    # plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    # plt.gcf().set_size_inches(6,4.5)
    # plt.savefig('{0}Resid_B_vs_ineff_LeastSq{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
    # plt.show()

    # f, (ax1,ax3) = plt.subplots(2, sharex=True, sharey=True)
    # ax1.errorbar(x_dataB, residLP, xerr=0, yerr=percy_unc, fmt= 'o', capsize=2, color='#ff7f0e', label='Linear fit residuals')
    # ax1.axhline(y=0,color='#1f77b4',linewidth=1)
    # ax1.legend(loc='lower right', fontsize=10)
    # ax1.set_title('{0} LS residuals'.format(sourcename))
    # ax1.set_ylabel("%")
    # ax3.errorbar(x_dataB, residCP, xerr=0, yerr=percy_unc, fmt= 'o', capsize=2, color='#2ca02c',label='Cubic fit residuals')
    # ax3.axhline(y=0,color='#1f77b4',linewidth=1)
    # ax3.legend(loc='lower right',fontsize=10)
    # ax3.set_xlabel(r"$1 - N_C / N_\gamma $")
    # ax3.set_ylabel("%")
    # f.subplots_adjust(hspace=0)
    # plt.xlim(xMinRB,xMaxRB)
    # plt.tight_layout()
    # plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    # plt.gcf().set_size_inches(6,4.5)
    # plt.savefig('{0}ResidP_B_vs_ineff_LeastSq{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
    # plt.show()
    

    return

