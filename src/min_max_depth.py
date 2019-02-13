#!/usr/bin/env python

""" Set minimum and maximum depth.

"""
import numpy as np
from scipy.spatial import cKDTree
import xarray as xr

from tools import ll2xyz, read_webtide


def main():

    wtdset = 'HRglobal'
    hmin = 5.
    hmax = 4000.
    ampfac = 1.2

    # read masked bathymetry
    ana = xr.open_dataset('data/interim/oi_bathymetry_masked.nc')['ana'].stack(xy=('x', 'y'))
    mask = ana.where(ana == 0., other=1.)

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

    # compute the minimum depth according to Maraldi et al. 2013
    mindepth = np.maximum(hmin, ampfac * df.ampmax[inds].values)

    # replace shallow grid points with minimum depth
    ana = np.maximum(ana, mindepth)
    ana[mask == 0.] = 0.

    # set maximum depth
    ana[ana > hmax] = hmax

    # prepare output dataset
    dsout = ana.unstack().to_dataset(name='ana')

    # write to file
    dsout.to_netcdf('data/interim/oi_bathymetry_depth_corrected.nc')


if __name__ == '__main__':

    # run main method
    main()
