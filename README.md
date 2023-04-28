### NOPP_Utilities  

Helper codes for tranferring data to/from COAWST grids and analyzing data / model output associated with the NOPP Hurricane Coastal Impacts project.  

Trying to keep the `NOPP_environment.yml` file current...this lists the packages needed by these notebooks.  

#### Examples of how to do stuff
##### Landcover
`Demo_CCAP2COAWST` - Demonstrates lookup table to convert C-CAP coastal landcover classes to ROMS vegetation parameters.  
##### Water levels
`Demo_hwm_Ian` - Demonstrates how to read high-water mark and storm tide gauge files, find those locations in ROMS output files, and make (crude) maps (no masking of non-water points in this version) and make scatterplots with multiple linear regressions. TODO - Add masking, add other statistics. Demonstarates with Ian data.  
`Demo_hwm_Michael` - Similar to `Demo_hwm_Ian`. No masking or fance map. Scatterplot has discrete color bar.  
`Demo_plot_RDG_Ian` - This reads time series of water-level stations during Ian and makes plots. No model data included.  
`Demo_plot_WL_Michael` - Similar to `Demo_plot_RDG_Ian`...reads water-level data and plots. No model data.  
`Plot_Michael_unfiltered` - More complete demo of reading storm tide instruments and comparing with time series of modeled peak water levels. Includes examples of masking, rudimentary `plotly` map. Has a demo of using `hvplot` at the end.  
##### Waves
`Ian_L1_wave_comp` - Demonstrates loading NDBC buoy data with xarray and plots time series of Hs. Compares with UFL L1 model results. No statistics.  
`commpare_model_drifter` - Demonstrates reading spotter buoy data from pickle files, and reading SWAN model output using `wavespectra` and comparing time series.  Calculates RMSE...need skill statistic.  
`compare_model_ndbc` - Demonstrates reading NDBC data with `wavespectra`. (Can update this now that `ndbc` is part of `wavespectra` distro). Needs to be run again and tested.  
`example_model_buoy_comparison` - From *Isabelle Houghton*. Demos reading spotter buoys and building a 2dspec. Also, estimating wave components (a1, b1, ...) for model data.  
`spotter_pickle2wavespec` - Modified from *Jake Davis* - Demonstrates wrangling spotter buoy data into `wavespectra` format...but only works for 1dspec so far. Need to get it working for 2dspec.  
`test_model_buoy_comparison` - Rambling demo including reading via `xarray` and `wavespectra, based largely on `example_model_buoy_comparison`
`to_wavespec_estimator` - Derived from `spotter_pickle2wavespec` - Work in progress trying to get spotter 2dspec into `wavespectra` format.
`wave_analysis_nopp_ian_ppt_reggie` - From *Maitane*. Example includes nice maps of best tracks with spotter trajectories and comparison with SWAN output.  
`wave_plot_funcs.py` - Collection of functions to do model-data comparison plots.  
##### Input for wave models
`make_spec_input` - Short script to make `SWAN` input files requesting output at locations of NDBC buoys.  
`make_swan_input` - Reads spotter, microswift, and dwsd  buoy locations, plot them, and write `SWAN` input files to request 2dspec at these times and locations.  
##### Other
`yanda_wget_script` - Example for constructing `wget` commands for downloading NDBC buoy data, from *Yanda*.  
`test_labels_to_netcdf` - Test of converting Doodled labels to `.nc` files.  
`load_doodler_class` - Reads the Doodler `.nc` files.  












