

## Genrating TEC maps

1. Use "TEC_maps.ipynb" to download the relevant ionex files into "ionex..." folders
2. Download the goes15 1m EUVs .nc files for 2013 and 2014 from : https://satdat.ngdc.noaa.gov/sem/goes/data/euvs/netcdf/goes15/ and place it inside "EUVs" folder
3. Run "DMDc_TEC_EUV.ipync" file in order to generate DMD TEC predictions (at code block #10 change the "
"days=3" variable in order to utilize different "control input" lenght)


## Generating NEU errors
1. Install "gLab" on you pc and make sure you can run it from terminal
2. Use "gLab_executer.ipynb" in order to execute the possition estimation using gLAB tool
3. To do so you will have to place the nessesary file, such as sattelite observation(.14o), orbits(.sp3) and clocks(.clk_30s) 
corsponding to a choosen date and place it inside ./gLAB/files/disturbance (or quiet)/<date>
###see 27_10_2014 date example.
4. The generated files will be in the output directory. 
5.Use "gLab_output.ipynb" file inorder to generate the NEU graphs.