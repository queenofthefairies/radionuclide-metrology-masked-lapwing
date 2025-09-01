# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 20:49:22 2017

@author: hepvis
"""

from __main__ import *

def f1(x, m, c):
    """ The linear function y= m*x + c """
    return m*x + c

def f3(x, a, c, d):
    """ The cubic function y= a*x^3 + c*x + d """
    return a*x**3  + c*x + d

# LINEAR BG/C vs G/C-1
poptLW, pcovLW = curve_fit(f1, xdataBGC, ydataBGC, p0=None, sigma=yuncBGC, absolute_sigma=True)
xtp = np.linspace(0, xMaxBGC, 200)
ylinfitW = f1(xtp, *poptLW)
beta_lin2=[poptLW[0],poptLW[1]]
# CUBIC BG/C vs G/C-1
poptcW, pcovcW = curve_fit(f3, xdataBGC, ydataBGC, p0=None, sigma=yuncBGC, absolute_sigma=True)
ycfitW = f3(xtp, *poptcW)
beta_cub2=[poptcW[0],poptcW[1]]
#==============================================================================
# BGC vs G/C-1 PLOT MAKER
#==============================================================================
residLT = ydataBGC - f1(xdataBGC, *poptLW)
residCT = ydataBGC - f3(xdataBGC, *poptcW)
residLP = (ydataBGC - f1(xdataBGC, *poptLW))/ydataBGC*100
residCP = (ydataBGC - f3(xdataBGC, *poptcW))/ydataBGC*100
percyunc = yuncBGC / ydataBGC *100
unc_intL=np.sqrt(pcovLW[1][1])
unc_intC=np.sqrt(pcovcW[2][2])
reluncL=unc_intL/poptLW[1]*100
reluncC=unc_intC/poptcW[2]*100
actfromL=poptLW[1]/(branchingratio*1000)
actfromC=poptcW[2]/(branchingratio*1000)
uncL=unc_intL/poptLW[1]*actfromL
uncC=unc_intC/poptcW[2]*actfromC
cornertext="Activity from linear intercept: \n"+r"{0:.3f} $\pm$ {1:.3f} MBq/g ({2:.1f}%)".format(actfromL,uncL,reluncL)+"\n"+"Activity from cubic intercept: \n"+r"{0:.3f}  $\pm$ {1:.3f} MBq/g ({2:.1f}%)".format(actfromC,uncC,reluncC)

#PLOT
f, ax = plt.subplots(1,1)
ax.errorbar(xdataBGC, ydataBGC, xerr=0, yerr=yuncBGC, fmt= 'o', capsize=2, label='{0}'.format(sourcename))
ax.plot(xtp, ylinfitW, label='Linear fit LS')
ax.plot(xtp, ycfitW, label='Cubic fit LS')
plt.locator_params(axis='y', nbins=6)
plt.locator_params(axis='x', nbins=6)
plt.xticks(fontsize=fsize-1, rotation=0)
plt.yticks(fontsize=fsize-1, rotation=0)
plt.xlim(0,xMaxBGC)
plt.ylim(yMinBGC,yMaxBGC)
ax.legend(loc=BGC_legpos, fontsize=fsize-1)
ax.add_artist(AnchoredText(cornertext,prop=dict(size=fsize-1),loc=BGC_txtpos,frameon=False,borderpad=1))
ax.set_xlabel(r"$N_\gamma / N_C - 1$",fontsize=fsize)
ax.set_ylabel(r"$N_\beta N_\gamma / N_C $ (s$^{-1}$ mg $^{-1}$)",fontsize=fsize)

ax2 = ax.twiny()
new_tick_locations = np.linspace(0, xMaxBGC, 6)
def tick_function(X):
    V = 1/(1+X)*100
    return ["%.1f" % z for z in V]

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r"$\beta$ efficiency %",fontsize=fsize)
plt.locator_params(axis='y', nbins=6)
plt.xticks(fontsize=fsize-1, rotation=0)
plt.yticks(fontsize=fsize, rotation=0)
plt.ylim(yMinBGC,yMaxBGC)
plt.tight_layout()
plt.gcf().set_size_inches(pwidth,pheight)
plt.savefig('{0}_BGC_vs_GoC-1_LeastSq_PLOT_{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
plt.show()
    
#RESIDUAL PLOTS
f, (ax1,ax3) = plt.subplots(2, sharex=True, sharey=True)
ax1.errorbar(xdataBGC, residLT, xerr=0, yerr=yuncBGC, fmt= 'o', capsize=2, color='#ff7f0e', label='Linear fit residuals')
ax1.axhline(y=0,color='#1f77b4',linewidth=1)
ax1.legend(loc='lower right', fontsize=10)
ax1.set_title('{0} LS residuals'.format(sourcename))
ax1.set_ylabel(r"$N_\beta N_\gamma / N_C $ (s$^{-1}$ mg $^{-1}$)")
ax3.errorbar(xdataBGC, residCT, xerr=0, yerr=yuncBGC, fmt= 'o', capsize=2, color='#2ca02c',label='Cubic fit residuals')
ax3.axhline(y=0,color='#1f77b4',linewidth=1)
ax3.legend(loc='lower right',fontsize=10)
ax3.set_xlabel(r"$N_\gamma / N_C - 1$")
ax3.set_ylabel(r"$N_\beta N_\gamma / N_C $ (s$^{-1}$ mg $^{-1}$)")
f.subplots_adjust(hspace=0)
plt.xlim(xMinRBGC,xMaxRBGC)
plt.ylim(yMinRBGC,yMaxRBGC)
plt.tight_layout()
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.gcf().set_size_inches(6,4.5)
plt.savefig('{0}Resid_BGC_vs_GoC-1_LeastSq{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
plt.show()

f, (ax1,ax3) = plt.subplots(2, sharex=True, sharey=True)
ax1.errorbar(xdataBGC, residLP, xerr=0, yerr=percyunc, fmt= 'o', capsize=2, color='#ff7f0e', label='Linear fit residuals')
ax1.axhline(y=0,color='#1f77b4',linewidth=1)
ax1.legend(loc='lower right', fontsize=10)
ax1.set_title('{0} LS residuals'.format(sourcename))
ax1.set_ylabel("%")
ax3.errorbar(xdataBGC, residCP, xerr=0, yerr=percyunc, fmt= 'o', capsize=2, color='#2ca02c',label='Cubic fit residuals')
ax3.axhline(y=0,color='#1f77b4',linewidth=1)
ax3.legend(loc='lower right',fontsize=10)
ax3.set_xlabel(r"$N_\gamma / N_C - 1$")
ax3.set_ylabel("%")
f.subplots_adjust(hspace=0)
plt.xlim(xMinRBGC,xMaxRBGC)
plt.tight_layout()
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.gcf().set_size_inches(6,4.5)
plt.savefig('{0}ResidP_BGC_vs_GoC-1_LeastSq{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
plt.show()
