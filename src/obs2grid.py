#!/usr/bin/env python

""" Interpolate observations to model grid.

"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import xarray as xr

from tools import read_grid, ll2xyz


def main():

    # read pandas.DataFrames with extracted data and concatenate datasets
    msldf = pd.read_hdf('data/interim/NWATL21_subset_msl.h5')
    msldf = msldf.rename(columns={'Med_depth': 'depth'})
    llwltdf = pd.read_hdf('data/interim/NWATL21_subset_llwlt_corrected.h5')
    df = pd.concat([msldf, llwltdf], ignore_index=True, sort=False)

    # free up some memory
    del msldf, llwltdf

    # read model grid
    grd = read_grid('data/raw/GoMSS/Configuration/coordinates.nc')

    # maximum distance from grid points within which to find observations
    dmax = np.sqrt(grd.e1t.values**2 + grd.e2t.values**2).max()

    # convert lon/lat to Cartesian coordinates
    xm, ym, zm = ll2xyz(grd.glamt.values.ravel(), grd.gphit.values.ravel())
    xo, yo, zo = ll2xyz(df.Lon.values, df.Lat.values)

    # create cKDTree object to represent source grid
    tree = cKDTree(list(zip(xm, ym, zm)))

    # find indices of the nearest points in the flattened array
    dist, ind = tree.query(list(zip(xo, yo, zo)),
                           k=1,
                           distance_upper_bound=dmax,
                           n_jobs=-1)

    # add index as column to data frame, this way we can use pandas groupby
    # functionality to easily calculate the median for each model grid point
    df['ind'] = ind
    df = df.loc[~np.isinf(dist)] # this removes all observations w/o grid point

    # calculate median of observations at each grid point
    df = df.groupby('ind')['depth'].median()

    # create grid
    depth = np.zeros((grd.dims['y'], grd.dims['x'])).ravel()
    depth[df.index] = df
    depth = depth.reshape((grd.dims['y'], grd.dims['x']))

    # create output dataset
    dsout = xr.Dataset(data_vars={'depth': (['x', 'y'], depth)},
                       coords={'lon': (['x', 'y'], grd.glamt),
                               'lat': (['x', 'y'], grd.gphit)})

    # write to file
    dsout.to_netcdf('data/interim/NWATL21_subset_gridded.nc')

if __name__ == '__main__':

    # run main method
    main()
