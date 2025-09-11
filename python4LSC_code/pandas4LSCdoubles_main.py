# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 14:55:29 2017

@author: siobhant
"""
from __main__ import *

#==============================================================================
# BACKGROUND stats functions
#==============================================================================
def bckgrndX(group):
    """ Mean X background """
    d = group['Xr8']
    c1 = group.count(axis=0)[0]
    return d.sum()/c1

def bckgrndLD(group):
    """ Mean LD background """
    d = group['LDr8']
    c1 = group.count(axis=0)[0]
    return d.sum()/c1

def bckgrndLDX(group):
    """ Mean LDX background """
    d = group['LDXr8']
    c1 = group.count(axis=0)[0]
    return d.sum()/c1
    
def bckgrndX_var(group):
    """ variance X background """
    d = group['Xr8']
    m = bckgrndX(group)
    c1 = group.count(axis=0)[0]
    return ((d-m)**2).sum()/(c1-1)
    
def bckgrndLD_var(group):
    """ variance X background """
    d = group['LDr8']
    m = bckgrndLD(group)
    c1 = group.count(axis=0)[0]
    return ((d-m)**2).sum()/(c1-1)

def bckgrndLDX_var(group):
    """ variance X background """
    d = group['LDXr8']
    m = bckgrndLDX(group)
    c1 = group.count(axis=0)[0]
    return ((d-m)**2).sum()/(c1-1)

#==============================================================================
# Data extraction functions
#==============================================================================
def background(df):
    df.head()
    # background ABC Rate
    df['ABCr8']=df[' ABC']/df[' Live']
    # background X Rate
    df['Xr8']=df[' X']/df[' Live']
    # background ABCX Rate
    df['ABCXr8']=df[' ABCX']/df[' Live']
        # background LD Rate
    df['LDr8']=df[' LD']/df[' Live']
        # background LDX Rate
    df['LDXr8']=df[' LDX']/df[' Live']
    return df



def weighted_mean_calc(unc_temp, n0_temp):
    weight = 1/unc_temp**2
    if weight_mean == 1:
        wght_x_n0 = (n0_temp)*weight
        wm_out = np.nansum(wght_x_n0)/np.nansum(weight)
    else:
        wm_out = np.nanmean(n0_temp)    
    return wm_out, weight 


def stdev_mean_theor_calc(weight):
    stdev_mean_theor = np.sqrt(1/np.nansum(weight))    
    return stdev_mean_theor




def threshold(df,backdf):
    backdf.head()
    backdf=background(backdf)
    df.head()
    df['backLDr8']=bckgrndLD(backdf)
    df['backXr8']=bckgrndX(backdf)
    df['backLDXr8']=bckgrndLDX(backdf)
    #calculate time difference between start time and reference time
    df['Real as s']=pd.to_timedelta(df[' Real'],unit='s')
    #df['Finished g']=pd.to_datetime(df[' Finished'])
    df['timedif']=refdatetime-(df[' Finished']-df['Real as s'])
    #change to seconds
    df['timedif(s)']=df['timedif']/pd.Timedelta(seconds=1)
    #calculate decay factor
    df['Decay Factor']=(np.log(2)*df[' Real']/halflifeseconds/
             (1-np.exp(-np.log(2)*df[' Real']/halflifeseconds))*
             0.5**(df['timedif(s)']/halflifeseconds))
    # LD Rate, decay and background corrected
    df['LDr8']=(df[' LD']/df[' Live']-df['backLDr8'])*df['Decay Factor']
    # X Rate, decay and background corrected
    df['Xr8']=(df[' X']/df[' Live']-df['backXr8'])*df['Decay Factor']
    # LDX Rate, decay and background corrected
    df['LDXr8']=(df[' LDX']/df[' Live']-df['backLDXr8'])*df['Decay Factor']
    
  
    
    ### FREDAS EQN MATCHING
    df['uncBG/C'] = (df['LDr8']*df['Xr8']/df['LDXr8'])*(dilution/mass)*np.sqrt((2*df[' LDX'])/(df[' LD']*df[' X']) - 1/df[' X'] - 1/df[' LD'] + 1/df[' LDX'])
    wmBGC_CT, weightBGC_CT = weighted_mean_calc(df['uncBG/C'], (df['LDr8']*df['Xr8']/df['LDXr8'])*(dilution/mass)) 
    stdev_theorBGC = stdev_mean_theor_calc(weightBGC_CT)
    stdev_obsBGC = wmBGC_CT*np.sqrt(np.nanvar(df['LDr8'], ddof=1)/np.nanmean(df['LDr8'])**2 + np.nanvar(df['Xr8'], ddof=1)/np.nanmean(df['Xr8'])**2 + np.nanvar(df['LDXr8'], ddof=1)/np.nanmean(df['LDXr8'])**2 
                                       + (2*np.cov(df['LDr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDr8'])*np.nanmean(df['Xr8'])) 
                                          - (2*np.cov(df['LDXr8'],df['LDr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['LDr8']))  
                                              - (2*np.cov(df['LDXr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['Xr8'])))/np.sqrt(np.count_nonzero(~np.isnan(df['LDr8'])) - 1)
    
    
    df['uncB'] = df['LDr8']*(dilution/mass)*np.sqrt((2*df[' LDX'])/(df[' LD']*df[' X']) - 1/df[' X'] - 1/df[' LD'] + 1/df[' LDX'])
    wmB_CT_Freda, weightB_CT_Freda = weighted_mean_calc(df['uncB'], df['LDr8']*(dilution/mass)) 
    stdev_theorB_Freda = stdev_mean_theor_calc(weightB_CT_Freda)
    stdev_obsB_Freda = wmB_CT_Freda*np.sqrt(np.nanvar(df['LDr8'], ddof=1)/np.nanmean(df['LDr8'])**2 + np.nanvar(df['Xr8'], ddof=1)/np.nanmean(df['Xr8'])**2 + np.nanvar(df['LDXr8'], ddof=1)/np.nanmean(df['LDXr8'])**2 
                                       + (2*np.cov(df['LDr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDr8'])*np.nanmean(df['Xr8'])) 
                                          - (2*np.cov(df['LDXr8'],df['LDr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['LDr8']))  
                                              - (2*np.cov(df['LDXr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['Xr8'])))/np.sqrt(np.count_nonzero(~np.isnan(df['LDr8'])) - 1)
    stdev_obsB_ODRonly = wmB_CT_Freda*np.sqrt(np.nanvar(df['LDr8'], ddof=1)/np.nanmean(df['LDr8'])**2)/np.sqrt(np.count_nonzero(~np.isnan(df['LDr8'])) - 1)
    
    
    df['uncG/C-1'] = (df['Xr8']/df['LDXr8'] - 1)*np.sqrt(df[' X']/(df[' LDX']*(df[' X'] - df[' LDX'])))
    wmGC_1_CT, weightGC_1_CT = weighted_mean_calc(df['uncG/C-1'], (df['Xr8']/df['LDXr8'] - 1)) 
    stdev_theorGC_1 = stdev_mean_theor_calc(weightGC_1_CT)
    stdev_obsGC_1 = (wmGC_1_CT+1)*np.sqrt(np.nanvar(df['LDXr8'], ddof=1)/np.nanmean(df['LDXr8'])**2 + np.nanvar(df['Xr8'], ddof=1)/np.nanmean(df['Xr8'])**2
                     - (2*np.cov(df['LDXr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['Xr8'])))/np.sqrt(np.count_nonzero(~np.isnan(df['LDXr8'])) - 1)    
    
    
    
    df['unc1-C/G'] = (1-df['LDXr8']/df['Xr8'])*np.sqrt(df[' LDX']/(df[' X']*(df[' X'] - df[' LDX'])))
    wm1_CG_CT, weight1_CG_CT = weighted_mean_calc(df['unc1-C/G'], (1-df['LDXr8']/df['Xr8'])) 
    stdev_theor1_CG = stdev_mean_theor_calc(weight1_CG_CT)
    stdev_obs1_CG = (1-wm1_CG_CT)*np.sqrt(np.nanvar(df['LDXr8'], ddof=1)/np.nanmean(df['LDXr8'])**2 + np.nanvar(df['Xr8'], ddof=1)/np.nanmean(df['Xr8'])**2
                     - (2*np.cov(df['LDXr8'],df['Xr8'],rowvar=False)[0,1])/(np.nanmean(df['LDXr8'])*np.nanmean(df['Xr8'])))/np.sqrt(np.count_nonzero(~np.isnan(df['LDXr8'])) - 1)    
  
    
   
    
    unc_params = np.array([wmGC_1_CT, stdev_theorGC_1, stdev_obsGC_1, wmBGC_CT, stdev_theorBGC, stdev_obsBGC, 
                           wm1_CG_CT, stdev_theor1_CG, stdev_obs1_CG, wmB_CT_Freda, stdev_theorB_Freda, stdev_obsB_Freda, stdev_obsB_ODRonly])
    
    
       
    #Constants to include in dataframe
    df['Reference datetime']=refdatetime
    df['Source mass mg']=mass
    df['Source dilution factor']=dilution
    return df, unc_params 


def dataget(dirname,filename,backfilename):
    dfpath="{0}/{1}".format(dirname,filename)
    print()
    backdfpath="{0}/{1}".format(dirname,backfilename)
    #Treats zero counts as corrupted files with na_values=0
    
    #UNCOMMENT TO READ CSV FILES
#    thresh1df = pd.read_csv(dfpath, sep=',', skiprows=5, na_values=0, index_col=0,dayfirst=True, parse_dates=[20], infer_datetime_format=True)
#    back1df = pd.read_csv(backdfpath, sep=',', skiprows=5, na_values=0, index_col=0,dayfirst=True, parse_dates=[20], infer_datetime_format=True)
    
    #UNCOMMENT TO READ EXCEL FILES
    thresh1df = pd.read_excel(dfpath, skiprows=5, na_values=0, index_col=0)
    back1df = pd.read_excel(backdfpath, skiprows=5, na_values=0, index_col=0)
    
    
    #Deletes NaN row that exists whenever read_excel is used (NaN row seems to be automatically skipped over when using read_csv)
    for ind_name in thresh1df.index:
        thresh1df_nan_test = thresh1df.loc[ind_name]
        if np.isnan(thresh1df_nan_test.iloc[0]) and np.isnan(thresh1df_nan_test.iloc[1]) and np.isnan(thresh1df_nan_test.iloc[2]):
            thresh1df = thresh1df.drop([ind_name])
    
    
    thresh1df, unc_params = threshold(thresh1df,back1df)
    print('done {0}'.format(dfpath[:-5]))
    return thresh1df, unc_params

#==============================================================================
# STATSGET 
# Generates values for use in regressions, and calculates a variety of 
# statistics...
# For G/C-1 and 1-C/G ('independent variables'): weighted mean, 
# theoretical standard deviation of weighted mean.
# 
# For BG/C and B ('Dependent variables'): weighted mean, 
# theoretical standard deviation of weighted mean, 
# theoretical standard deviation of weighted mean as a relative %, 
# observed variance in the weighted mean, 
# observed standard deviation in weighted mean, 
# and the observed standard deviation of the weighted mean as a relative % 
# mean
# 
#==============================================================================

def statsget(unc_params):
    return [(unc_params[0], unc_params[1], unc_params[2],
             unc_params[3], unc_params[4], unc_params[5], '',
             unc_params[6], unc_params[7], unc_params[8],
             unc_params[9], unc_params[10], unc_params[11], unc_params[12])]
    
#==============================================================================
# For Speedy File Grabber
#==============================================================================

def startstring(char, stringlist):
    newlist = []
    for string in stringlist:
        if string.startswith(char):
            newlist.append(string)
    newlist.sort()
    return newlist

#==============================================================================
# Speedy file grabber
#==============================================================================
if speedyfilegrabber:
    path  = os.getcwd()
    path1 = "{0}/{1}".format(path,rtdt1dir)
    filenames1 = os.listdir(path1)
    
    bkgs1 = startstring(bkgprefix,filenames1)
    files1 = startstring(dataprefix, filenames1)
    
    
    print()
    print('Running pandas4LSC')
    print('Checking files')
    LB=len(bkgs1)
    LD=len(files1)
    if LB == LD:
        pass
    else:
        print()
        print('!!!ATTENTION!!! Different number of background files and data files')
        print('Since there is a mismatch with background and data files')
        print('Take another look at your files and rename or delete as necessary, then run code again')
        print()
        sys.exit('')
    
    print("Check background files have corresponding threshold data files:")
    i=0
    if LD > LB:
        l=LB
    else:
        l=LD
    while i < l:
        print("{0} --- {1}".format(bkgs1[i],files1[i]))
        print()
        i=i+1
    g2g=input('Does every entry in the list match its background file? type Y or N:  ')
    
    if g2g=='N':
        print()
        print('Since there is a mismatch with background and data files')
        print('Take another look at your files and rename or delete as necessary, then run code again')
        print()
        sys.exit('')
    if g2g=='Y':
        print('Data analysis happening now!')
else:
    pass
#==============================================================================
# Organising all those files and directories
#=========================================="170914_bkg_30mV.csv"====================================

rtdtdirectories=(rtdt1dir)
rtdtinfo=(rt,dt)

if speedyfilegrabber:
    bckgrnds1=bkgs1 #for speedy file grabber
    threshfiles1=files1 #for speedy file grabber
else:
    threshfiles1=(thresh1A,thresh1B,thresh1C,thresh1D,thresh1E,thresh1F,thresh1G,
              thresh1H,thresh1I,thresh1J,thresh1K,thresh1L,thresh1M,thresh1N,
              thresh1O,thresh1P,thresh1Q)
    bckgrnds1=(back1A,back1B,back1C,back1D,back1E,back1F,back1G,back1H,back1I,
           back1J,back1K,back1L,back1M,back1N,back1O,back1P,back1Q)
##==============================================================================
## Building Data Frames
##==============================================================================

threshcolumns= [
' A', ' B', ' C',' AB', ' AC', ' BC', ' ABC', ' X', ' ABX', ' ACX', ' BCX', 
' ABCX', ' Real', ' Live', ' LD', ' LDX', ' Started', 
' Finished', 'Source mass mg', 'Source dilution factor', 'Reference datetime', 
'timedif(s)','Decay Factor', 'backLDr8','LDr8','backXr8', 'Xr8','backLDXr8', 'LDXr8', 
'uncBG/C', 'uncB', 'uncG/C-1', 'unc1-C/G']

threshdf=pd.DataFrame(data=[], columns=(threshcolumns))


regcolumns=['G/C-1 WM', 'Stdev of Mean G/C-1 (Theor)', 'Stdev of Mean G/C-1 (Obs)',
'BG/C WM', 'Stdev of Mean BGC (Theor)', 'Stdev of Mean BGC (Obs)', '',
'1-C/G WM', 'Stdev of Mean 1-C/G (Theor)', 'Stdev of Mean 1-C/G (Obs)', 
'B WM', 'Stdev of Mean B (Theor)', 'Stdev of Mean B (Obs)', 'Stdev of Mean B (Obs)(ODR only)']

regdf=pd.DataFrame(data=[],columns=(regcolumns))

#==============================================================================
# ANALYSIS
#==============================================================================

dirname=rtdt1dir
i=0
while i < len(threshfiles1):
    if threshfiles1[i]:
        thresh1df, unc_params = dataget(dirname,threshfiles1[i],bckgrnds1[i])
        stats1=statsget(unc_params)
        reg1df=pd.DataFrame(data=stats1,index=[i],columns=regcolumns)
        threshdf=threshdf.append(thresh1df)
        regdf=regdf.append(reg1df)
        thresh1df.drop(thresh1df.index, inplace=True)
        reg1df.drop(reg1df.index, inplace=True)
    i=i+1

regdf=regdf[regcolumns]
regdf=regdf.sort_values('G/C-1 WM')
regdf.to_excel("{0}_rt{1}dt{2}_RegData_{3}{4}{5}_newunceqns.xlsx".format(outputfilename,rt,dt,DorT,Sb,WM))
regdf.drop(regdf.index, inplace=True)
threshdf=threshdf[threshcolumns]
threshdf.to_excel("{0}_rt{1}dt{2}_ThreshData_{3}{4}{5}_newunceqns.xlsx".format(outputfilename,rt,dt,DorT,Sb,WM))
threshdf.drop(threshdf.index, inplace=True)


print()
print('Finished with different thresholds in {0}'.format(dirname))

    
   ############ END ############
