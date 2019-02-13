#!/usr/bin/env python

""" Set minimum and maximum depth.

"""
import numpy as np
from scipy.spatial import cKDTree
import xarray as xr

from tools import ll2xyz, read_webtide

from dask.distributed import Client

import matplotlib.pyplot as plt

def main():

    client = Client()

    wtdset = 'HRglobal'
    hmin = 5.
    hmax = 4000.
    ampfac = 1.2

    # read masked bathymetry
    ana = xr.open_dataset('data/interim/oi_bathymetry_masked.nc')['ana'].stack(xy=('x', 'y'))
    ana = ana.where(ana != 0.)

    # read WebTide data
    df = read_webtide(wtdset, ampmax=True)

    # source coordinates: WebTide
    xs, ys, zs = ll2xyz(df['lon'].values, df['lat'].values)

    # target coordinates: observations
    xt, yt, zt = ll2xyz(ana['lon'].values, ana['lat'].values)

    # create cKDTree object to represent source grid
    kdt = cKDTree(list(zip(xs, ys, zs)))

    # find indices of the nearest npt grid points in the flattened array
    dist, inds = kdt.query(list(zip(xt, yt, zt)),
                           k=1,
                           n_jobs=-1)

    mindepth = np.maximum(hmin, ampfac * df.ampmax[inds].values)

    # replace shallow grid points with minimum depth
    ana = np.maximum(ana, mindepth).unstack()

    # prepare output dataset
    dsout = ana.to_dataset(name='ana')

    dsout['ana'].plot()
    plt.show()

    # write to file
    # dsout.to_netcdf('data/interim/oi_bathymetry_masked.nc')


if __name__ == '__main__':

    # run main method
    main()
