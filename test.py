#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

from downscaling import main

wd_hresf  = './test/wd_hres.tif'
wsl_lresf = './test/wsl_lres.elev'
dem_hresf = './test/dem.tif'
nodata    = -9999 # NODATA value in both data sets

main(dem_hresf=dem_hresf,wsl_lresf=wsl_lresf,wd_hresf=wd_hresf,nodata=nodata)
