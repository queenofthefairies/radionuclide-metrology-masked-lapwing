# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 17:55:30 2017

@author: hepvis
"""

from __main__ import *
#==============================================================================
# Least Squares regression functions
#==============================================================================
def f1(x, m, c):
    """ The linear function y= m*x + c """
    return m*x + c

def fit(xdata, ydata, yunc, branchratio):
    poptLW, pcovLW = curve_fit(f1, xdata, ydata, p0=None, sigma=yunc, absolute_sigma=True)
    stdevm=np.sqrt(pcovLW[0][0])
    stdevc=np.sqrt(pcovLW[1][1])
    relunc=stdevc/poptLW[1]*100
    actconc=poptLW[1]/(1000*branchratio)
    fitparams=[(branchratio, poptLW[0],stdevm,poptLW[1],stdevc,actconc,relunc)]
    fitcols=['Branch ratio','gradient','unc gradient', 'intercept','unc intercept','Activity conc', 'rel unc intercept %']
    transferdf=pd.DataFrame(data=fitparams,columns=fitcols)
    return transferdf

def f3(x, a, c, d):
    """ The cubic function y= a*x**3 + c*x + d """
    return a*x**3 + c*x + d

def fit3(xdata, ydata, yunc, branchratio):
    poptLW, pcovLW = curve_fit(f3, xdata, ydata, p0=None, sigma=yunc, absolute_sigma=True)
    stdeva=np.sqrt(pcovLW[0][0])
    stdevc=np.sqrt(pcovLW[1][1])
    stdevd=np.sqrt(pcovLW[2][2])
    relunc=stdevd/poptLW[2]*100
    actconc=poptLW[2]/(1000*branchratio)
    fit3params=[(branchratio,poptLW[0],stdeva,poptLW[1],stdevc,poptLW[2],stdevd,actconc,relunc)]
    fit3cols=['Branch ratio','a','unc a', 'c', 'unc c', 'intercept','unc intercept','Activity conc', 'rel unc intercept %']
    transfer3df=pd.DataFrame(data=fit3params,columns=fit3cols)
    return transfer3df
    
#==============================================================================
# Orthogonal Distance Regression functions
#==============================================================================
def f1ODR(p,x):
    """ The linear function y= m*x + c """
    m,c=p
    return m*x + c

def f3ODR(p,x):
    """ The cubic function y= a*x**3 + c*x + d """
    a,c,d=p
    return a*x**3 + c*x + d

def fitODR(xdata, ydata, xunc, yunc, branchratio): 
    lin_model=odear.Model(f1ODR)
    lin_reg_data=odear.RealData(xdata,ydata,sx=xunc,sy=yunc)
    lin_ODR=odear.ODR(lin_reg_data, lin_model, beta0=beta_lin)
    lin_out= lin_ODR.run()
    fitcols=['Branch ratio','gradient','unc gradient', 'intercept','unc intercept','Activity conc', 'rel unc intercept %']
    linreg_info=[(branchratio, lin_out.beta[0],lin_out.sd_beta[0],lin_out.beta[1],lin_out.sd_beta[1],lin_out.beta[1]/(1000*branchratio),lin_out.sd_beta[1]/lin_out.beta[1]*100)]
    linregdf=pd.DataFrame(data=linreg_info,columns=(fitcols))
    return linregdf

def fit3ODR(xdata, ydata, xunc, yunc, branchratio):
    cub_model=odear.Model(f3ODR)
    cub_reg_data=odear.RealData(xdata,ydata,sx=xunc,sy=yunc)
    cub_ODR=odear.ODR(cub_reg_data, cub_model, beta0=beta_cub)
    cub_out= cub_ODR.run()
    fit3cols=['Branch ratio','a','unc a', 'c', 'unc c', 'intercept','unc intercept','Activity conc', 'rel unc intercept %']
    cubreg_info=[(branchratio, cub_out.beta[0],cub_out.sd_beta[0],cub_out.beta[1],
                  cub_out.sd_beta[1],cub_out.beta[2],cub_out.sd_beta[2],cub_out.beta[2]/(1000*branchratio),cub_out.sd_beta[2]/cub_out.beta[2]*100)]
    cubregdf=pd.DataFrame(data=cubreg_info,columns=(fit3cols))
    return cubregdf

#==============================================================================
# Select data from data frame to perform regressions on
#==============================================================================
    
if StandDev=='OSD':    
    xdataBGC = np.array(regdf)[:,0]
    ydataBGC = np.array(regdf)[:,3]
    xuncBGC = np.array(regdf)[:,2]
    yuncBGC = np.array(regdf)[:,5]
    
    xdataB = np.array(regdf)[:,7]
    ydataB = np.array(regdf)[:,10]
    xuncB = np.array(regdf)[:,9]
    yuncB = np.array(regdf)[:,12]
    yuncB_ODRonly = np.array(regdf)[:,13]
    
if StandDev=='TSD':    
    xdataBGC = np.array(regdf)[:,0]
    ydataBGC = np.array(regdf)[:,3]
    xuncBGC = np.array(regdf)[:,1]
    yuncBGC = np.array(regdf)[:,4]
    
    xdataB = np.array(regdf)[:,7]
    ydataB = np.array(regdf)[:,10]
    xuncB = np.array(regdf)[:,8]
    yuncB = np.array(regdf)[:,11]



regdf.drop(regdf.index, inplace=True)

#==============================================================================
# Data frames that will contain fit parameters etc.
#==============================================================================
fitparamscolumns=[
'Regression data', 'Regression method',
'gradient','unc gradient','intercept','unc intercept',
'Activity conc','rel unc intercept %']

fitparamsdf=pd.DataFrame(data=[],columns=(fitparamscolumns))

fit3paramscolumns=['Regression data', 'Regression method',
'a','unc a', 'c', 'unc c', 'intercept',
'unc intercept','Activity conc', 'rel unc intercept %']

fit3paramsdf=pd.DataFrame(data=[],columns=(fit3paramscolumns))

fitsummarycols=['Function', 'Regression data', 'Regression method',
'Activity conc', 'rel unc intercept %']

fitsummarydf=pd.DataFrame(data=[],columns=(fitsummarycols))

infodat=[('','','','',''),
         ('Source name',outputfilename,'','',''),
         ('Branch ratio',branchratio,'','',''),
         ('Reference date time',refdatetime,'','',''),
         ('Half life s',halflifeseconds,'','',''), 
         ('Dilution factor',dilution,'','',''),
         ('Mass mg',mass,'','',''),
         ('Resolving time ns',rt,'','',''),
         ('Dead time us',dt,'','',''),
         ('Gamma shift',gs,'','',''),
         ('Standard deviation used as weights',StandDev,'','',''),
         ('Background',Sb,'','',''),
         ('Data',DorT,'','','')]
         
infodf=pd.DataFrame(data=infodat,columns=(fitsummarycols))
#==============================================================================
# REGRESSION TIME
#==============================================================================
# LINEAR regression. Data = BG/C vs G/C-1 Method = Least Squares
linregBGC_LS=fit(xdataBGC,ydataBGC,yuncBGC,branchratio)
BGC_LS=pd.DataFrame(data=[('Linear','BG/C vs G/C-1','Least Squares')],columns=['Function','Regression data','Regression method'])
linregBGC_LS=pd.concat([linregBGC_LS,BGC_LS],axis=1)
#fitparamsdf=fitparamsdf.append(linregBGC_LS)


# LINEAR regression. Data = BG/C vs G/C-1 Method = ODR
# Least squares becomes initial estimate for ODR
beta_lin=[np.array(linregBGC_LS)[0,3],np.array(linregBGC_LS)[0,5]]
linregBGC_ODR=fitODR(xdataBGC, ydataBGC, xuncBGC, yuncBGC,branchratio)
linBGC_ODR=pd.DataFrame(data=[('Linear','BG/C vs G/C-1','Orthogonal Distance')],columns=['Function','Regression data','Regression method'])
linregBGC_ODR=pd.concat([linregBGC_ODR,linBGC_ODR],axis=1)
#fitparamsdf=fitparamsdf.append(linregBGC_ODR)

# LINEAR regression. Data = B vs 1-C/G Method = Least Squares
linregB_LS=fit(xdataB,ydataB,yuncB,branchratio)
linB_LS=pd.DataFrame(data=[('Linear','B vs 1-C/G','Least Squares')],columns=['Function','Regression data','Regression method'])
linregB_LS=pd.concat([linregB_LS,linB_LS],axis=1)
#fitparamsdf=fitparamsdf.append(linregB_LS)

#LINEAR regression. Data = B vs 1-C/G Method = ODR
# Least squares becomes initial estimate for ODR
beta_lin=[np.array(linregB_LS)[0,3],np.array(linregB_LS)[0,5]]
linregB_ODR=fitODR(xdataB, ydataB, xuncB, yuncB, branchratio)
linB_ODR=pd.DataFrame(data=[('Linear','B vs 1-C/G','Orthogonal Distance')],columns=['Function','Regression data','Regression method'])
linregB_ODR=pd.concat([linregB_ODR,linB_ODR],axis=1)
#fitparamsdf=fitparamsdf.append(linregB_ODR)

# Bring all LINEAR regression results into one data frame
fitparamsdf = pd.concat([linregBGC_LS,linregBGC_ODR,linregB_LS,linregB_ODR])

# CUBIC regression. Data = BG/C vs G/C-1 Method = Least Squares
cubregBGC_LS=fit3(xdataBGC,ydataBGC,yuncBGC,branchratio)
cubBGC_LS=pd.DataFrame(data=[('Cubic','BG/C vs G/C-1','Least Squares')],columns=['Function','Regression data','Regression method'])
cubregBGC_LS=pd.concat([cubregBGC_LS,cubBGC_LS],axis=1)
#fit3paramsdf=fit3paramsdf.append(cubregBGC_LS)

# CUBIC regression. Data = BG/C vs G/C-1 Method = ODR
# Least squares becomes initial estimate for ODR
beta_cub=[np.array(cubregBGC_LS)[0,3],np.array(cubregBGC_LS)[0,5],np.array(cubregBGC_LS)[0,7]]
cubregBGC_ODR=fit3ODR(xdataBGC, ydataBGC, xuncBGC, yuncBGC, branchratio)
cubBGC_ODR=pd.DataFrame(data=[('Cubic','BG/C vs G/C-1','Orthogonal Distance')],columns=['Function','Regression data','Regression method'])
cubregBGC_ODR=pd.concat([cubregBGC_ODR,cubBGC_ODR],axis=1)
#fit3paramsdf=fit3paramsdf.append(cubregBGC_ODR)

# CUBIC regression. Data = B vs 1-C/G Method = Least Squares
cubregB_LS=fit3(xdataB,ydataB,yuncB,branchratio)
cubB_LS=pd.DataFrame(data=[('Cubic','B vs 1-C/G','Least Squares')],columns=['Function','Regression data','Regression method'])
cubregB_LS=pd.concat([cubregB_LS,cubB_LS],axis=1)
#fit3paramsdf=fit3paramsdf.append(cubregB_LS)

# CUBIC regression. Data = B vs 1-C/G Method = ODR
# Least squares becomes initial estimate for ODR
beta_cub=[np.array(cubregB_LS)[0,3],np.array(cubregB_LS)[0,5],np.array(cubregB_LS)[0,7]]
cubregB_ODR=fit3ODR(xdataB, ydataB, xuncB, yuncB, branchratio)
cubB_ODR=pd.DataFrame(data=[('Cubic','B vs 1-C/G','Orthogonal Distance')],columns=['Function','Regression data','Regression method'])
cubregB_ODR=pd.concat([cubregB_ODR,cubB_ODR],axis=1)
#fit3paramsdf=fit3paramsdf.append(cubregB_ODR)

# Bring all CUBIC regression results into one data frame
fit3paramsdf = pd.concat([cubregBGC_LS,cubregBGC_ODR,cubregB_LS,cubregB_ODR])

# CLEAR DATA FRAMES
linregBGC_LS.drop(linregBGC_LS.index, inplace=True)
linregBGC_ODR.drop(linregBGC_ODR.index, inplace=True)
linregB_LS.drop(linregB_LS.index, inplace=True)
linregB_ODR.drop(linregB_ODR.index, inplace=True)
cubregBGC_LS.drop(cubregBGC_LS.index, inplace=True)
cubregBGC_ODR.drop(cubregBGC_ODR.index, inplace=True)
cubregB_LS.drop(cubregB_LS.index, inplace=True)
cubregB_ODR.drop(cubregB_ODR.index, inplace=True)

#==============================================================================
# Saving fit information
#==============================================================================
print()
fitsummarydf=pd.concat([fitsummarydf,fit3paramsdf,fitparamsdf,infodf])
#fitsummarydf=fitsummarydf,fitparamsdf)
#fitsummarydf=fitsummarydfinfodf)

fit1df=pd.DataFrame(data=[('','','','','','','',''),('y=m*x+c','','','','','','','')],columns=(fitparamscolumns))
fitparamsdf=pd.concat([fitparamsdf,fit1df])
fit3df=pd.DataFrame(data=[('','','','','','','','','',''),('y=a*x^3 + c*x + d','','','','','','','','','')],columns=(fit3paramscolumns))
fit3paramsdf=pd.concat([fit3paramsdf,fit3df])

print("Finished with data")
fitparamsdf=fitparamsdf[fitparamscolumns]
fit3paramsdf=fit3paramsdf[fit3paramscolumns]
fitsummarydf=fitsummarydf[fitsummarycols]
output_excel_filename = "{0}_rt{1}dt{2}_AllFits_{3}{4}_{5}{6}_newunceqns.xlsx".format(outputfilename,rt,dt,DorT,Sb,StandDev,WM)
with pd.ExcelWriter(output_excel_filename) as writer:
    fitsummarydf.to_excel(writer,sheet_name="Summary")
    fitparamsdf.to_excel(writer,sheet_name="LinearFits")
    fit3paramsdf.to_excel(writer,sheet_name="CubicFits")
#fits_writer.save()
print("All fit info saved in {0}".format(output_excel_filename))
fitparamsdf.drop(fitparamsdf.index, inplace=True)
fit3paramsdf.drop(fit3paramsdf.index, inplace=True)
