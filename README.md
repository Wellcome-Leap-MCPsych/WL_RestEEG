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

## Notes:
1. The test data are in various forms based on what was shared to me. There are some weird things about it from memory.
   1. Stanford was unable to share a trimmed raw file, so shared an existing .SET file (i.e., output from EEGLAB). This meant I had to work with them to get their import functions correct without directly being able to test it.
   2. McLean has two different test sets; not because they are different internally but due to file naming conventions. The example data has been renamed using a different naming convention and so might the naming of output files might be a bit strange.
3. The .ced files are electrode location files and contain 3-D co-ordinates that allow the electrodes to be mapped on the head.
    1. These are site specific and based on how the data are collected.     
