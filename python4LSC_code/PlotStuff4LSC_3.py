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
beta_lin1=[estlin[2,2],estlin[2,4]]
beta_cub1=[estcub[2,2],estcub[2,4],estcub[2,6]]

#linear weighted

lin_model=odear.Model(f1)
lin_reg_data=odear.RealData(xdataB,ydataB,sx=xuncB,sy=yuncB)
lin_ODR=odear.ODR(lin_reg_data, lin_model, beta0=beta_lin1)
lin_out= lin_ODR.run()

xtp = np.linspace(0, xMaxB, 500)
ylinfitW = f1(lin_out.beta,xtp)

#CUBIC
cub_model=odear.Model(f3)
cub_reg_data=odear.RealData(xdataB,ydataB,sx=xuncB,sy=yuncB)
cub_ODR=odear.ODR(cub_reg_data, cub_model, beta0=beta_cub1)
cub_out= cub_ODR.run()

ycfitW = f3(cub_out.beta,xtp)  

#==============================================================================
# B vs 1-C/G PLOT MAKER
#==============================================================================
residLT = ydataB - f1(lin_out.beta,xdataB)
residCT = ydataB - f3(cub_out.beta,xdataB)
residLP = (ydataB - f1(lin_out.beta,xdataB))/ydataB*100
residCP = (ydataB - f3(cub_out.beta,xdataB))/ydataB*100
percyunc = yuncB / ydataB *100
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
ax.errorbar(xdataB, ydataB, xerr=xuncB, yerr=yuncB, fmt= 'o', capsize=2, label='{0}'.format(sourcename))
ax.plot(xtp, ylinfitW, label='Linear fit ODR')
ax.plot(xtp, ycfitW, label='Cubic fit ODR')
plt.locator_params(axis='y', nbins=6)
plt.locator_params(axis='x', nbins=6)
plt.xticks(fontsize=fsize-1, rotation=0)
plt.yticks(fontsize=fsize-1, rotation=0)
plt.xlim(0,xMaxB)
plt.ylim(yMinB,yMaxB)
ax.legend(loc=B_legpos, fontsize=fsize-1)
ax.add_artist(AnchoredText(cornertext,prop=dict(size=fsize-1),loc=B_txtpos,frameon=False,borderpad=1))
ax.set_xlabel(r"$1 - N_C /N_\gamma$",fontsize=fsize)
ax.set_ylabel(r"$N_\beta $ (s$^{-1}$ mg $^{-1}$)",fontsize=fsize)

ax2 = ax.twiny()
new_tick_locations = np.linspace(0, xMaxB, 6)
def tick_function(X):
    V = (1-X)*100
    return ["%.1f" % z for z in V]

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r"$\beta$ efficiency %",fontsize=fsize)
plt.locator_params(axis='y', nbins=6)
plt.locator_params(axis='x', nbins=6)
plt.xticks(fontsize=fsize-1, rotation=0)
plt.yticks(fontsize=fsize, rotation=0)
plt.ylim(yMinB,yMaxB)
plt.tight_layout()
plt.gcf().set_size_inches(pwidth,pheight)
plt.savefig('{0}_B_vs_ineff_ODR_PLOT_{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
plt.show()
    
#RESIDUAL PLOTS
f, (ax1,ax3) = plt.subplots(2, sharex=True, sharey=True)
ax1.errorbar(xdataB, residLT, xerr=xuncB, yerr=yuncB, fmt= 'o', capsize=2, color='#ff7f0e', label='Linear fit residuals')
ax1.axhline(y=0,color='#1f77b4',linewidth=1)
ax1.legend(loc='lower right', fontsize=10)
ax1.set_title('{0} ODR residuals'.format(sourcename))
ax1.set_ylabel(r"$N_\beta $ (s$^{-1}$ mg $^{-1}$)")
ax3.errorbar(xdataB, residCT, xerr=xuncB, yerr=yuncB, fmt= 'o', capsize=2, color='#2ca02c',label='Cubic fit residuals')
ax3.axhline(y=0,color='#1f77b4',linewidth=1)
ax3.legend(loc='lower right',fontsize=10)
ax3.set_xlabel(r"$1 - N_C / N_\gamma $")
ax3.set_ylabel(r"$N_\beta $ (s$^{-1}$ mg $^{-1}$)")
f.subplots_adjust(hspace=0)
plt.xlim(xMinRB,xMaxRB)
plt.ylim(yMinRB,yMaxRB)
plt.tight_layout()
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.gcf().set_size_inches(6,4.5)
plt.savefig('{0}Resid_B_vs_ineff_ODR{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
plt.show()

f, (ax1,ax3) = plt.subplots(2, sharex=True, sharey=True)
ax1.errorbar(xdataB, residLP, xerr=xuncB, yerr=percyunc, fmt= 'o', capsize=2, color='#ff7f0e', label='Linear fit residuals')
ax1.axhline(y=0,color='#1f77b4',linewidth=1)
ax1.legend(loc='lower right', fontsize=10)
ax1.set_title('{0} ODR residuals'.format(sourcename))
ax1.set_ylabel("%")
ax3.errorbar(xdataB, residCP, xerr=xuncB, yerr=percyunc, fmt= 'o', capsize=2, color='#2ca02c',label='Cubic fit residuals')
ax3.axhline(y=0,color='#1f77b4',linewidth=1)
ax3.legend(loc='lower right',fontsize=10)
ax3.set_xlabel(r"$1 - N_C / N_\gamma $")
ax3.set_ylabel("%")
f.subplots_adjust(hspace=0)
plt.xlim(xMinRB,xMaxRB)
plt.tight_layout()
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.gcf().set_size_inches(6,4.5)
plt.savefig('{0}ResidP_B_vs_ineff_ODR{1}{2}.png'.format(sourcename,DorT,SD),dpi=240)
plt.show()
