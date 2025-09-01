# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 20:49:22 2017

@author: hepvis
"""
from __main__ import *

def f1(p,x):
    """ The linear function y= m*x + c """
    m,c=p
    return m*x + c

def f3(p,x):
    """ The cubic function y= a*x^3 + c*x + d """
    a,c,d=p
    return a*x**3  + c*x + d

# GET ESTIMATES FROM FITS
estfilename=fitsfilename
estlin_df = pd.read_excel(estfilename,sheet_name="LinearFits",index_col=0)
estcub_df = pd.read_excel(estfilename,sheet_name="CubicFits",index_col=0)
estlin = np.array(estlin_df)
estcub = np.array(estcub_df)
beta_lin2=[estlin[0,2],estlin[0,4]]
beta_cub2=[estcub[0,2],estcub[0,4],estcub[0,6]]

#linear weighted

lin_model=odear.Model(f1)
lin_reg_data=odear.RealData(xdataBGC,ydataBGC,sx=xuncBGC,sy=yuncBGC)
lin_ODR=odear.ODR(lin_reg_data, lin_model, beta0=beta_lin2)
lin_out= lin_ODR.run()

xtp = np.linspace(0, xMaxBGC, 500)
ylinfitW = f1(lin_out.beta,xtp)

#CUBIC
cub_model=odear.Model(f3)
cub_reg_data=odear.RealData(xdataBGC,ydataBGC,sx=xuncBGC,sy=yuncBGC)
cub_ODR=odear.ODR(cub_reg_data, cub_model, beta0=beta_cub2)
cub_out= cub_ODR.run()

ycfitW = f3(cub_out.beta,xtp)    


residLT = ydataBGC - f1(lin_out.beta,xdataBGC)
residCT = ydataBGC - f3(cub_out.beta,xdataBGC)
residLP = (ydataBGC - f1(lin_out.beta,xdataBGC))/ydataBGC*100
residCP = (ydataBGC - f3(cub_out.beta,xdataBGC))/ydataBGC*100

#==============================================================================
# BG/C vs G/C-1 PLOT MAKER
#==============================================================================

percyunc = yuncBGC / ydataBGC *100
unc_intL=lin_out.sd_beta[1]
unc_intC=cub_out.sd_beta[2]
reluncL=unc_intL/lin_out.beta[1]*100
reluncC=unc_intC/cub_out.beta[2]*100
actfromL=lin_out.beta[1]/(branchingratio*1000)
actfromC=cub_out.beta[2]/(branchingratio*1000)
uncL=unc_intL/lin_out.beta[1]*actfromL
uncC=unc_intC/cub_out.beta[2]*actfromC
cornertext="Activity from linear intercept: \n"+r"{0:.3f} $\pm$ {1:.3f} MBq/g ({2:.1f}%)".format(actfromL,uncL,reluncL)+"\n"+"Activity from cubic intercept: \n"+r"{0:.3f}  $\pm$ {1:.3f} MBq/g ({2:.1f}%)".format(actfromC,uncC,reluncC)

#PLOT
f, ax = plt.subplots(1,1)
ax.errorbar(xdataBGC, ydataBGC, xerr=xuncBGC, yerr=yuncBGC, fmt= 'o', capsize=2, label='{0}'.format(sourcename))
ax.plot(xtp, ylinfitW, label='Linear fit ODR')
ax.plot(xtp, ycfitW, label='Cubic fit ODR')
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
plt.locator_params(axis='x', nbins=6)
plt.xticks(fontsize=fsize-1, rotation=0)
plt.yticks(fontsize=fsize, rotation=0)
plt.ylim(yMinBGC,yMaxBGC)
plt.tight_layout()
plt.gcf().set_size_inches(pwidth,pheight)
plt.savefig('{0}_BGC_vs_GoC-1_ODR_PLOT_{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
plt.show()
    
#RESIDUAL PLOTS
f, (ax1,ax3) = plt.subplots(2, sharex=True, sharey=True)
ax1.errorbar(xdataBGC, residLT, xerr=xuncBGC, yerr=yuncBGC, fmt= 'o', capsize=2, color='#ff7f0e', label='Linear fit residuals')
ax1.axhline(y=0,color='#1f77b4',linewidth=1)
ax1.legend(loc='lower right', fontsize=10)
ax1.set_title('{0} ODR residuals'.format(sourcename))
ax1.set_ylabel(r"$N_\beta N_\gamma / N_C $ (s$^{-1}$ mg $^{-1}$)")
ax3.errorbar(xdataBGC, residCT, xerr=xuncBGC, yerr=yuncBGC, fmt= 'o', capsize=2, color='#2ca02c',label='Cubic fit residuals')
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
plt.savefig('{0}Resid_BGC_vs_GoC-1_ODR{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
plt.show()

f, (ax1,ax3) = plt.subplots(2, sharex=True, sharey=True)
ax1.errorbar(xdataBGC, residLP, xerr=xuncBGC, yerr=percyunc, fmt= 'o', capsize=2, color='#ff7f0e', label='Linear fit residuals')
ax1.axhline(y=0,color='#1f77b4',linewidth=1)
ax1.legend(loc='lower right', fontsize=10)
ax1.set_title('{0} ODR residuals'.format(sourcename))
ax1.set_ylabel("%")
ax3.errorbar(xdataBGC, residCP, xerr=xuncBGC, yerr=percyunc, fmt= 'o', capsize=2, color='#2ca02c',label='Cubic fit residuals')
ax3.axhline(y=0,color='#1f77b4',linewidth=1)
ax3.legend(loc='lower right',fontsize=10)
ax3.set_xlabel(r"$N_\gamma / N_C - 1$")
ax3.set_ylabel("%")
f.subplots_adjust(hspace=0)
plt.xlim(xMinRBGC,xMaxRBGC)
plt.tight_layout()
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.gcf().set_size_inches(6,4.5)
plt.savefig('{0}ResidP_BGC_vs_GoC-1_ODR{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
plt.show()
