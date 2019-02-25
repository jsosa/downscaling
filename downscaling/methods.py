#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import pandas as pd
import xarray as xr
from subprocess import call


def read_static_inputs(dem_hresf, msk_hresf):

    # Reading geo information for both
    dem_hres = xr.open_rasterio(dem_hresf).sel(band=1)

    # Reading high resolution rivers, channels and lakes mask
    # Replace all values >0 with ones
    call(['gdal_proximity.py', '-distunits', 'GEO',
          msk_hresf, './msk_hres_prox.tif'])
    msk_hres = xr.open_rasterio('./msk_hres_prox.tif').sel(band=1)

    return dem_hres, msk_hres


def read_dynamic_inputs(wsl_lresf, wd_lresf, net_lresf):

    # Reading coarse resolution water surface elevation
    wsl_lres = xr.open_rasterio(wsl_lresf).sel(band=1)

    # Reading coarse resolution water depth
    wd_lres = xr.open_rasterio(wd_lresf).sel(band=1)

    # Reading coarser river network
    call(['gdal_proximity.py', '-distunits', 'GEO',
          net_lresf, './net_lres_prox.tif'])
    net_lres_prox = xr.open_rasterio('./net_lres_prox.tif').sel(band=1)
    net_lres_prox_wd = net_lres_prox[:, :] * (wd_lres > 0).astype(int).data

    # Creating dataframes
    wsl_lres_df = wsl_lres.to_dataframe('elev').reset_index(level=['x', 'y'])
    wd_lres_df = wd_lres.to_dataframe('wd').reset_index(level=['x', 'y'])
    net_lres_prox_wd_df = net_lres_prox_wd.to_dataframe(
        'prox').reset_index(level=['x', 'y'])
    del wsl_lres_df['band']
    del wd_lres_df['band']
    del net_lres_prox_wd_df['band']

    wsl_lres_df = wsl_lres_df[wd_lres_df['wd'] > 0]
    net_lres_prox_wd_df = net_lres_prox_wd_df[wd_lres_df['wd'] > 0]
    wd_lres_df = wd_lres_df[wd_lres_df['wd'] > 0]

    wsl_lres_df['wd'] = wd_lres_df['wd']
    wsl_lres_df['vol'] = wd_lres_df['wd'] * 30/3600 * 30/3600
    wsl_lres_df['prox'] = net_lres_prox_wd_df['prox']

    return wsl_lres_df


def _process_df(df, dem, msk, wd_hres, x_hres, y_hres):

    # Read proximity high-res map slice
    # Read previous wd_hres high-res slice
    msk_slice = msk.sel(x=slice(df.xmin, df.xmax), y=slice(df.ymax, df.ymin))
    win_before = wd_hres.sel(x=slice(df.xmin, df.xmax),
                             y=slice(df.ymax, df.ymin))

    # Substract low-res water surface elevation from high-res dem slice
    wd = df.elev-dem.sel(x=slice(df.xmin, df.xmax), y=slice(df.ymax, df.ymin))

    # Replace by zero negative values and values where previous wd_hres high-res is zero
    wd = wd.where((wd > 0) & (win_before == 0), 0)

    # Convert to dataframes
    wd_df = wd.to_dataframe('wd')
    msk_df = msk_slice.to_dataframe('prox')

    # Sorting water depth df by proximity
    # Check when values are lower than low-res VOL, otherwise assign zero
    # Sort values back to its original index
    wd_df['prox'] = msk_df['prox']
    wd_df['vol'] = wd_df['wd'] * 3/3600 * 3/3600
    wd_df['idx'] = range(len(wd_df))
    wd_df.sort_values(by='prox', inplace=True)
    wd_df['cumvol'] = wd_df['vol'].cumsum()
    wd_df.wd.where(wd_df['cumvol'] <= df.vol, 0)
    wd_df.sort_values(by='idx', inplace=True)

    # Water depth dataframe to xr.DataArray
    # Adding values from previous water depth slice
    # Assign new values to big array, wd_hres
    wd = wd_df.wd.to_xarray().sel(y=slice(None, None, -1))
    wd = wd + win_before
    wd_hres.loc[dict(x=x_hres[(x_hres > df.xmin) & (x_hres < df.xmax)], y=y_hres[(
        y_hres > df.ymin) & (y_hres < df.ymax)])] = wd.values


def prepare_df(df, thr_val=(30/3600)*4):

    df = df.sort_values(by=['elev','prox','wd'],ascending=[False,True,True])
    df = df[df['prox'] >= 0]
    df['xmin'] = df['x'] - thr_val
    df['ymin'] = df['y'] - thr_val
    df['xmax'] = df['x'] + thr_val
    df['ymax'] = df['y'] + thr_val

    return df


def run_algorithm(df, dem, msk):

    # Array to be filled
    wd_hres = xr.full_like(dem, fill_value=0)
    x_hres = wd_hres.x
    y_hres = wd_hres.y

    # Run algorithm for every dataframe row
    df.apply(_process_df, dem=dem, msk=msk, wd_hres=wd_hres, x_hres=x_hres, y_hres=y_hres, axis=1)

    return wd_hres


def apply_decay(wd_hres, msk_hres, val=100):

    wd_hres[:, :] = wd_hres[:, :] * ((1-msk_hres)**val).data
    
    return wd_hres


def save_dataset(fname, wd_hres):

    xr.Dataset({'wd': wd_hres}).to_netcdf(
        fname, encoding={'wd': {'zlib': True}})
