#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os
import rasterio as rio
import xarray as xr
import numpy as np
import multiprocessing as mp
from subprocess import call

def downscaling(dem_hresf,msk_hresf,wsl_lresf,wd_lresf,thr_val,thr_dpt,thr_decay,wd_hresf,gtiff=True):

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
    wsl_lres = xr.open_rasterio(wsl_lresf).sel(band=1)

    # Reading coarse resolution water depth
    wd_lres = xr.open_rasterio(wd_lresf).sel(band=1)

    # Reading geo information for both
    dem_hres = xr.open_rasterio(dem_hresf).sel(band=1)

    # Selecting water depth > thr_dpt
    wsl_lres_msk = wsl_lres.where(wd_lres>thr_dpt) #thr_dpt
    wsl_lres_msk = wsl_lres_msk.to_series().to_frame().dropna().reset_index(level=['x','y'])
    wsl_lres_msk['xmin'] = wsl_lres_msk['x'] - thr_val
    wsl_lres_msk['ymin'] = wsl_lres_msk['y'] - thr_val
    wsl_lres_msk['xmax'] = wsl_lres_msk['x'] + thr_val
    wsl_lres_msk['ymax'] = wsl_lres_msk['y'] + thr_val

    # Create a empty array, the output array
    wd_hres = xr.full_like(dem_hres,fill_value=0)
    x_hres = wd_hres.x
    y_hres = wd_hres.y

    for xmin,ymin,xmax,ymax,wsl in zip(wsl_lres_msk['xmin'],wsl_lres_msk['ymin'],
                                       wsl_lres_msk['xmax'],wsl_lres_msk['ymax'],
                                       wsl_lres_msk[0]):

        wd = wsl-dem_hres.sel(x=slice(xmin,xmax),y=slice(ymax,ymin))
        wd = wd.where(wd>0,0)
        win = wd_hres.sel(x=slice(xmin,xmax),y=slice(ymax,ymin))
        wd_hres.loc[dict(x=x_hres[(x_hres>xmin)&(x_hres<xmax)], y=y_hres[(y_hres>ymin)&(y_hres<ymax)])] = ((wd+win)/2).values

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
    mask = xr.open_rasterio(outfolder+'buffer_dist.tif').sel(band=1)
    wd_hres[:,:] = wd_hres[:,:] * ((1-mask)**thr_decay).data

    # Write final high resolution water depth map
    xr.Dataset({'wd':wd_hres}).to_netcdf(wd_hresf,encoding={'wd':{'zlib':True}})

    # Output GTIFF file
    if gtiff:
        name = os.path.dirname(wd_hresf) + '/' + os.path.basename(wd_hresf).split('.')[0] + '.tif'
        with rio.open(name, 'w',
                driver ='GTiff',
                height = len(wd_hres.y),
                width = len(wd_hres.x),
                count=1,
                dtype=str(wd_hres.dtype),
                crs=wd_hres.crs,
                compress='lzw',
                transform=wd_hres.transform[0:6]) as dst:
            dst.write(wd_hres.data,1)

def downscaling_mp(dem_hresf,msk_hresf,wsl_lresf,wd_lresf,thr_val,thr_dpt,thr_decay,wd_hresf,gtiff=True):

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
    wsl_lres = xr.open_rasterio(wsl_lresf).sel(band=1)

    # Reading coarse resolution water depth
    wd_lres = xr.open_rasterio(wd_lresf).sel(band=1)

    # Reading geo information for both
    dem_hres = xr.open_rasterio(dem_hresf).sel(band=1)

    # Selecting water depth > thr_dpt
    wsl_lres_msk = wsl_lres.where(wd_lres>thr_dpt) #thr_dpt
    wsl_lres_msk = wsl_lres_msk.to_series().to_frame().dropna().reset_index(level=['x','y'])
    wsl_lres_msk['xmin'] = wsl_lres_msk['x'] - thr_val
    wsl_lres_msk['ymin'] = wsl_lres_msk['y'] - thr_val
    wsl_lres_msk['xmax'] = wsl_lres_msk['x'] + thr_val
    wsl_lres_msk['ymax'] = wsl_lres_msk['y'] + thr_val

    # Create a empty array, the output array
    wd_hres = xr.full_like(dem_hres,fill_value=0)

    def downscaling_mp_func(output,dem_hres,wd_hres,wsl_lres_msk):
        
        x_hres = dem_hres.x
        y_hres = dem_hres.y
        
        for xmin,ymin,xmax,ymax,wsl in zip(wsl_lres_msk['xmin'],wsl_lres_msk['ymin'],
                                        wsl_lres_msk['xmax'],wsl_lres_msk['ymax'],
                                        wsl_lres_msk[0]):

            wd = wsl-dem_hres.sel(x=slice(xmin,xmax),y=slice(ymax,ymin))
            wd = wd.where(wd>0,0)
            win = wd_hres.sel(x=slice(xmin,xmax),y=slice(ymax,ymin))
            wd_hres.loc[dict(x=x_hres[(x_hres>xmin)&(x_hres<xmax)], y=y_hres[(y_hres>ymin)&(y_hres<ymax)])] = ((wd+win)/2).values
        
        output.put(wd_hres)

    # Define processors
    nproc = mp.cpu_count()

    # Split x and y in nproc parts
    split_df = [df for g, df in wsl_lres_msk.groupby(np.arange(len(wsl_lres_msk)) // (len(wsl_lres_msk)/nproc))]

    # Define an output queue
    output = mp.Queue()

    # Setup a list of processes that we want to run
    processes = [mp.Process(target=downscaling_mp_func, args=(output,dem_hres,wd_hres,split_df[i])) for i in range(len(split_df))]

    # Run processes
    for p in processes:
        p.start()

    # Get process results from the output queue
    results = [output.get() for p in processes]

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
    mask = xr.open_rasterio(outfolder+'buffer_dist.tif').sel(band=1)
    wd_hres[:,:] = (sum(results)).data * ((1-mask)**thr_decay).data

    # Write final high resolution water depth map
    xr.Dataset({'wd':wd_hres}).to_netcdf(wd_hresf,encoding={'wd':{'zlib':True}})

    # Output GTIFF file
    if gtiff:
        name = os.path.dirname(wd_hresf) + '/' + os.path.basename(wd_hresf).split('.')[0] + '.tif'
        with rio.open(name, 'w',
                driver ='GTiff',
                height = len(wd_hres.y),
                width = len(wd_hres.x),
                count=1,
                dtype=str(wd_hres.dtype),
                crs=wd_hres.crs,
                compress='lzw',
                transform=wd_hres.transform[0:6]) as dst:
            dst.write(wd_hres.data,1)
