Follow instructions for installing IOOS environment.  
Download environment.yml file
Edit .yml file:
rename environment to NOPP
Ensure the following are in the .yml file (for XROMS) 
cf_xarray
dask
jupyter
jupyterlab
netcdf4
numba >= 0.49
numpy
xarray
xgcm

cartopy >= 0.18
cmocean
datashader >= 0.11
geoviews
holoviews
hvplot
zarr

wavespectra
 
Install environment  


`conda activate NOPP`
`conda install -c conda-forge jupyterlab`

#### get xroms code
`git clone git@github.com:xoceanmodel/xroms.git`
`cd ~crs/xroms`
`$ pip install -e ` 