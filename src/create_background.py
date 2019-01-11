#!/usr/bin/env python

""" Create background field for optimal interpolation.

"""
import numpy as np
from scipy.spatial import cKDTree
import xarray as xr

from tools import read_grid, ll2xyz

import matplotlib.pyplot as plt


def main():

    # read model grid
    grd = read_grid('data/raw/GoMSS/Configuration/coordinates.nc')

    # domain boundaries - hard coded for now
    lonmin = grd.glamt.min() - .1
    lonmax = grd.glamt.max() + .1
    latmin = grd.gphit.min() - .1
    latmax = grd.gphit.max() + .1

    ds = read_gebco('data/external/GEBCO/GEBCO_2014_2D.nc',
                    lonmin=lonmin, lonmax=lonmax,
                    latmin=latmin, latmax=latmax)

    glon, glat = np.meshgrid(ds.lon.values, ds.lat.values)

    # convert lon/lat to Cartesian coordinates
    xs, ys, zs = ll2xyz(glon.ravel(), glat.ravel())
    xt, yt, zt = ll2xyz(grd.glamt.values.ravel(), grd.gphit.values.ravel())

    # create cKDTree object to represent source grid
    tree = cKDTree(list(zip(xs, ys, zs)))

    # find indices of the nearest points in the flattened array
    dist, inds = tree.query(list(zip(xt, yt, zt)),
                           k=4,
                           n_jobs=-1)

    # calculate weights based on iverse distance
    wghts = 1.0 / dist**2
    ele = np.sum(wghts * ds.elevation.values.ravel()[inds], axis=1) / np.sum(wghts, axis=1)
    ele = ele.reshape((grd.dims['y'], grd.dims['x']))

    # create output dataset
    dsout = xr.Dataset(data_vars={'ele': (['x', 'y'], ele)},
                       coords={'lon': (['x', 'y'], grd.glamt),
                               'lat': (['x', 'y'], grd.gphit)})

    # write to file
    dsout.to_netcdf('data/interim/GEBCO_background.nc')


def read_gebco(fname, lonmin=-180., lonmax=180., latmin=-90., latmax=90.):
    '''Read and subset GEBCO data.'''

    # load GEBCO file
    ds = xr.open_dataset(fname)

    # subset data
    ds = ds.sel(lon=slice(lonmin, lonmax), lat=slice(latmin, latmax))

    return ds

if __name__ == '__main__':

    # run main method
    main()
