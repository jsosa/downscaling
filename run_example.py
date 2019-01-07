#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

from downscaling import downscaling

wd_hresf  = './example/wd_hres.tif'
wsl_lresf = './example/034_1533.tif'
dem_hresf = './example/034_dem.tif'
nodata    = -9999 # NODATA value in both data sets

downscaling(dem_hresf=dem_hresf, wsl_lresf=wsl_lresf, wd_hresf=wd_hresf, thr_val=0.00416*3, thr_dpt=0, nodata=nodata)
