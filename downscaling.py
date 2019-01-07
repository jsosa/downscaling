#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import sys
import numpy as np
import gdalutils as gu

# TODO:
# thr_dpt should be applied to water depth not water surface elevation!

def downscaling(dem_hresf,wsl_lresf,wd_hresf,thr_val,thr_dpt,nodata):

    """
    Downscale LFP coarse water surface elevation results to higher reslution water depth
    based on the high resolution DEM. It receives a threshold indicating the window in the DEM
    where the water will be spread.
    """

    # Reading high resolution DEM
    dem_hres = gu.get_data(dem_hresf)

    # Reading coarse resolution water surface elevation
    wsl_lres = gu.get_data(wsl_lresf)
    wsl_lres_msk = np.ma.masked_values(wsl_lres,nodata)
    
    # Reading geo information for both
    geo_hres = gu.get_geo(dem_hresf)
    geo_lres = gu.get_geo(wsl_lresf)

    # Selecting water surface elevatoin > thr_dpt, IT SHOULD BE WATER DEPTH (COARSE RESOLUTION
    iy,ix = np.where(wsl_lres_msk>=thr_dpt)
    x_lres = geo_lres[8][ix]
    y_lres = geo_lres[9][iy]

    # Create a empty array, the output array
    wd_hres = np.zeros_like(dem_hres)

    # Iterate over every coarse water surface elevation pixel
    for i in range(len(x_lres)):

        # Searching window over the pixel
        xmin = x_lres[i] - thr_val
        ymin = y_lres[i] - thr_val
        xmax = x_lres[i] + thr_val
        ymax = y_lres[i] + thr_val

        # Clipping the high resolution DEM for previous window
        dat_dem, geo_dem = gu.clip_raster(dem_hresf, xmin, ymin, xmax, ymax)
        dat_dem_msk = np.ma.masked_values(dat_dem,nodata)

        # Getting water surface elevation value
        wsl_val = wsl_lres_msk[iy[i],ix[i]]

        # Substracting water surface eleavation to values in DEM window
        # Only positive values are valid cells
        wd = wsl_val-dat_dem_msk
        wd = np.ma.masked_where(wd<0,wd)
        wd = np.ma.masked_where(wd.filled(0)<0,wd)

        # Getting indexes from high DEM
        yind0,yind1,xind0,xind1 = get_index_geo(geo_dem,geo_hres)
        
        # Reading high resolution water depth values
        win = wd_hres[yind0:yind1,xind0:xind1]

        # If window is with values and mean is greater than 0
        if np.mean(win)>0:
            
            # Average previous values
            c = np.ones([win.shape[0],win.shape[1],2]) * np.nan
            c[:,:,0] = win
            c[:,:,1] = wd.filled(0)
            wd_hres[yind0:yind1,xind0:xind1] = np.mean(c,axis=2)
        
        # Otherwise substitute by values
        else:
            wd_hres[yind0:yind1,xind0:xind1] = wd.filled(0)
    
    # Write final high resolution water depth map
    gu.write_raster(wd_hres,wd_hresf,geo_hres,'Float64',0)

def get_index_geo(geo,geo2):
	
    xmin = geo[8].min()
    ymin = geo[9].min()
    xmax = geo[8].max()
    ymax = geo[9].max()
    xx = geo2[8]
    yy = geo2[9]
    xind0 = abs(xx-xmin).argmin()
    xind1 = abs(xx-xmax).argmin()+1  # inclusive slicing
    yind0 = abs(yy-ymax).argmin()
    yind1 = abs(yy-ymin).argmin()+1  # inclusive slicing
    return yind0,yind1,xind0,xind1
