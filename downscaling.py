#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import sys
import numpy as np
import gdalutils as gu

def main(dem_hresf,wsl_lresf,wd_hresf,nodata):

    dem_hres = gu.get_data(dem_hresf)
    wsl_lres = gu.get_data(wsl_lresf)
    wsl_lres_msk = np.ma.masked_values(wsl_lres,nodata)
    geo_hres = gu.get_geo(dem_hresf)
    geo_lres = gu.get_geo(wsl_lresf)
    iy,ix = np.where(wsl_lres>nodata)
    x_lres = geo_lres[8][ix]
    y_lres = geo_lres[9][iy]

    thresh = 0.00416
    wd_hres = np.zeros_like(dem_hres)
    for i in range(len(x_lres)):
        xmin = x_lres[i] - thresh
        ymin = y_lres[i] - thresh
        xmax = x_lres[i] + thresh
        ymax = y_lres[i] + thresh
        dat_dem, geo_dem = gu.clip_raster(dem_hresf, xmin, ymin, xmax, ymax)
        dat_dem_msk = np.ma.masked_values(dat_dem,nodata)
        dem_krnl = wsl_lres_msk[iy[i],ix[i]]
        wd = dem_krnl-dat_dem_msk
        wd = np.ma.masked_where(wd<0,wd)
        wd = np.ma.masked_where(wd.filled(0)<0,wd)
        yind0,yind1,xind0,xind1 = get_index_geo(geo_dem,geo_hres)
        wd_hres[yind0:yind1,xind0:xind1] = wd.filled(0)
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

wd_hresf  = './test/wd_hres.tif'
wsl_lresf = './test/wsl_lres.elev'
dem_hresf = './test/dem.tif'
nodata    = -9999 # NODATA value in both data sets

main(dem_hresf=dem_hresf,wsl_lresf=wsl_lresf,wd_hresf=wd_hresf,nodata=nodata)
