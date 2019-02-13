#!/usr/bin/env python

""" Correct the vertical datum of observation to MSL.

"""
import numpy as np
import os
import pandas as pd
from scipy.spatial import cKDTree

from tools import ll2xyz, read_webtide


def main():

    # specify WebTide dataset
    wtdset = 'HRglobal'

    # read dask.DataFrame with extracted data
    obsdf = pd.read_hdf('data/interim/NWATL21_subset_llwlt.h5')

    # read WebTide data
    wtdf = read_webtide(wtdset, ampmax=True)

    # source coordinates: WebTide
    xs, ys, zs = ll2xyz(wtdf['lon'].values, wtdf['lat'].values)

    # target coordinates: observations
    xt, yt, zt = ll2xyz(obsdf['Lon'].values, obsdf['Lat'].values)

    # create cKDTree object to represent source grid
    kdt = cKDTree(list(zip(xs, ys, zs)))

    # find indices of the nearest npt grid points in the flattened array
    dist, inds = kdt.query(list(zip(xt, yt, zt)), k=1)

    # initialize
    outdf = obsdf.drop('Med_depth', axis=1)
    outdf['depth'] = np.zeros_like(outdf.index)

    # compute correction factor (from bathycor.c by D. Greenberg)
    ffact = .85 + .15 * (wtdf.ampmax[inds].values - 2.) / 3.
    ffact[wtdf.ampmax[inds].values <= 2.] = .85
    ffact[wtdf.ampmax[inds].values >= 5.] = 1.

    # correction step - it is assumed that the input depths are negative down
    outdf.depth = obsdf.Med_depth - ffact * wtdf.ampmax[inds].values

    # write output files
    store = pd.HDFStore('data/interim/NWATL21_subset_llwlt_corrected.h5')
    store.put('df', outdf, data_columns=outdf.columns)
    store.close()


if __name__ == '__main__':

    # run main method
    main()
