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

    # Reading coarse resolution water surface elevation
    wsl_lres = gu.get_data(wsl_lresf)
    
    # Reading coarse resolution water depth
    wd_lres = gu.get_data(wd_lresf)

    # Reading geo information for both
    geo_hres = gu.get_geo(dem_hresf)
    geo_lres = gu.get_geo(wsl_lresf)

    # Selecting water depth > thr_dpt
    iy,ix = np.where(wd_lres>thr_dpt)
    x_lres = geo_lres[8][ix]
    y_lres = geo_lres[9][iy]

    # Create a empty array, the output array
    wd_hres = np.zeros((len(geo_hres[9]),len(geo_hres[8])))

    # Iterate over every coarse water surface elevation pixel
    for i in range(len(x_lres)):

        # Searching window over the pixel
        xmin = x_lres[i] - thr_val
        ymin = y_lres[i] - thr_val
        xmax = x_lres[i] + thr_val
        ymax = y_lres[i] + thr_val

        # Clipping the high resolution DEM for previous window
        dat_dem, geo_dem = gu.clip_raster(dem_hresf, xmin, ymin, xmax, ymax)

        # Getting indexes from high resolution DEM
        yind0,yind1,xind0,xind1 = get_index_geo(geo_dem[8],geo_dem[9],geo_hres[8],geo_hres[9])
        
        # Reading 'from previous iteration' high resolution water depth values
        win = wd_hres[yind0:yind1,xind0:xind1]

        # Getting water surface elevation value
        # Substracting water surface eleavation to values in DEM window
        # Only positive values are valid cells
        wd = wsl_lres[iy[i],ix[i]]-dat_dem
        wd[wd<0]=0        

        # Average previous values
        wd_hres[yind0:yind1,xind0:xind1] = (win+wd)/2 # np.mean([win,wd],axis=0)        

    # Apply a exonential decay over the inundated area
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
                              msk_hresf,
                              outfolder+'buffer_dist.tif'])

    # Multiply by decay function
    mask = gu.get_data(outfolder+'buffer_dist.tif')
    wd_hres = wd_hres * (1-mask)**thr_decay

    # Write final high resolution water depth map
    gu.write_raster(wd_hres,wd_hresf,geo_hres,'Float64',0)

def get_index_geo(x1,y1,x2,y2):
    
    xmin  = x1.min()
    ymin  = y1.min()
    xmax  = x1.max()
    ymax  = y1.max()
    xind0 = abs(x2-xmin).argmin()
    xind1 = abs(x2-xmax).argmin()+1  # inclusive slicing
    yind0 = abs(y2-ymax).argmin()
    yind1 = abs(y2-ymin).argmin()+1  # inclusive slicing
    return yind0,yind1,xind0,xind1
