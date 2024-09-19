# WL_RestEEG

This repo contains the Matlab scripts used to harmonise the EEG pre-processing across WL data collection sites.
It also contains the script used by all sites to generate the spectral power estimates using the pre-processed data.

## Organisation

The scripts are organised as follows:
1. Files specific to each site (i.e., processing script and accompanying electrode location files) are allocated into the relevant folder, these include:
    1. Helsinki
    2. McLean
    3. Stanford
    4. UCSD
2. Files used to generate parameters for analysis are in the "ParameterDerivation" Folder.
    1. Currently this only includes the spectral power estimates. 
