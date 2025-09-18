# Python4LSC code
This contains code used to analyse data from liquid scintillation coincidence counting experiments in radionuclide metrology. It also includes a test dataset (Ho166m). The original version of the code was written in 2018; the code was updated in 2025 to account for accidental coincidences, and updates to various python packages. The 2025 scripts have a modular structure consisting of:     
 - `pandas4LSC_user` script, which calls functions in    
    -  `pandas4LSC_data_processing` to deal with raw data taken at different threshold voltages,     
    - plus `pandas4LSC_data_analysis` to do the regression on the data. 

Scripts to plot the results are still completely separate.    
- - `PlotStuff4LSC_user` script, which calls functions in    
    -  `PlotStuff4LSC` to plot the linearised threshold data, along with the cubic or linear fit, and a separate residuals plot






 

