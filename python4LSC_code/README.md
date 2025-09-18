# Python4LSC code
This contains code used to analyse data from liquid scintillation coincidence counting experiments in radionuclide metrology. It also includes a test dataset (Ho166m). The original version of the code was written in 2018; the code was updated in 2025 to account for accidental coincidences, and updates to various python packages. The 2025 scripts have a modular structure consisting of:     
 - `pandas4LSC_user` script, which calls functions in    
    -  `pandas4LSC_data_processing` to deal with raw data taken at different threshold voltages,     
    - plus `pandas4LSC_data_analysis` to do the regression on the data. 

## Analysis of threshold data and regression
Brief outline of the `pandas4LSC_user` script.

### User input
- user needs to input constants such as reference date time, half life, dilution factor, mass, branch ratio, source name, resolving time, dead time and gamma shift    
- user needs to choose options:
    - doubles or triples
    - specified background (i.e. background xlsx file already prepared with background values for each threshold), or process raw background files
    - accidental coincidences correction
    - which standard deviation to use to weight the points (observed standard deviation, OSD, or theoretical standard deviation, TSD)
    - whether to use the weighted mean or not (default is yes)
- user needs to supply threshold data files, and background threshold data files _or_ specified background file

### Functions called
**The module `pandas4LSC_data_processing` is imported as `LSC_data_processing`:**
- `LSC_data_processing.get_data` reads in the raw threshold data and returns a dataframe of that data 
    - `thresh_df` refers to threshold dataframe (measurements on source)
    - `raw_back_df` refers to background dataframe (raw background measurements)
- `LSC_data_processing.accidental_coincidences` calculates accidental coincidence corrections and adds these to the dataframe
- `LSC_data_processing.background_doubles_rates_corrected` calculates counting rates for the background data with optional accidental coinc correction (if not using specified background)
- `LSC_data_processing.background_average` calculates average background for each threshold and returns a new background dataframe 
    - `back_df` is the background dataframe either taken from the specified background file, or returned by the `background_average` function
- `LSC_data_processing.decay_factor` calculates the decay factor and adds this to the `thresh_df`
- `LSC_data_processing.doubles_rates_corrected` calculates the counting rates, corrected for background, decay factor, and optional accidental coinc correction for `thresh_df`
- `LSC_data_processing.linearise_threshold_data` calculates quantities such as _BG/C_ etc. for `thresh_df`
- `LSC_data_processing.stats_get` calculates the weighted means, std devs for the linearised quantities and puts these in a new dataframe `reg_df`
    - `reg_df` is the dataframe that contains the data that is used for regression

**The module `pandas4LSC_data_analysis is imported as `LSC_data_analysis`:**
- `LSC_data_analysis.regression` performs regression on `reg_df` and saves the results to the `_AllFits.xlsx` file
    - by default, this script runs every kind of regression (cubic/linear, least squares/orthogonal distance) on _BG/C_ vs _G/C-1_, and _B_ vs _1-C/G_ data

## Plotting
Scripts to plot the results are still completely separate.    
- `PlotStuff4LSC_user` script, which calls the `plotter` function in    
    -  `PlotStuff4LSC` to plot the linearised threshold data, along with the cubic or linear fit, and a separate residuals plot
    - only one fit is shown on the plot at a time






 

