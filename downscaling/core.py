#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os
import sys
import shutil
import numpy as np
import gdalutils as gu
from subprocess import call

def downscaling(dem_hresf,msk_hresf,wsl_lresf,wd_lresf,thr_val,thr_dpt,thr_decay,nodata,wd_hresf):

    """
    Downscale flood maps from coarse to a high resolution. Water
    surface elevation from the coarser grid is substracted from high
    resolution DEM, it will wet cells with elevation lower than water surface
    elvation. Water surface elevation can be threshold by `thr_dpt` applied
    to the coarse resolution water depth. `thr_val` is a threshold of spreadness.
    A decay in resulting high resolution water depths can be applied by `thr_decay`
    where thr_decay=0 produces no decay and thr_val>0 produces a decay following
    the Gaussian equation D=C**thr_val, where D is water depth and C is the proximity.
    Proximity C is given by proximity (distance) of mask specified by `msk_hresf`,
    `msk_hresf` can be location of main river and tributaries in the study area.
    """

    # Reading high resolution DEM
    dem_hres = gu.get_data(dem_hresf)

    # Reading coarse resolution water surface elevation
    wsl_lres = gu.get_data(wsl_lresf)
    wsl_lres_msk = np.ma.masked_values(wsl_lres,nodata)
    
    # Reading coarse resolution water depth
    wd_lres = gu.get_data(wd_lresf)
    wd_lres_msk = np.ma.masked_values(wd_lres,0)

    # Reading geo information for both
    geo_hres = gu.get_geo(dem_hresf)
    geo_lres = gu.get_geo(wsl_lresf)

    # Selecting water depth > thr_dpt
    iy,ix = np.where(wd_lres_msk>thr_dpt)
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

    # Create a temp folder
    outfolder = os.path.dirname(wd_hresf) + '/tmp/'
    try:
        os.mkdir(outfolder)
    except FileExistsError:
        pass

    # Calc proximity around `msk_hresf`
    call(['gdal_proximity.py','-distunits','GEO',
                              '-co','COMPRESS=LZW',
                              '-nodata','0',
                              msk_hresf,outfolder+'buffer_dist.tif'])

    mask = gu.get_data(outfolder+'buffer_dist.tif')
    weig = (1-mask)**thr_decay
    wd_final = wd_hres * weig

    # Write final high resolution water depth map
    gu.write_raster(wd_final,wd_hresf,geo_hres,'Float64',0)

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
