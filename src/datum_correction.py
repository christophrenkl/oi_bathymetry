#!/usr/bin/env python

""" Correct the vertical datum of observation to MSL.

"""
import numpy as np
import os
import pandas as pd
from scipy.spatial import cKDTree

from tools import ll2xyz


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


def read_webtide(wtdset, constituents=['M2', 'S2', 'N2', 'K1', 'O1'],
                 ampmax=False):
    '''Read WebTide data and compute maximum amplitude if desired.'''

    # directory with WebTide data
    wtdir = 'data/external/WebTide/data'

    # name of coordinate file depends on WebTide set
    if wtdset == 'HRglobal':
        cfile = 'HRgloballl.nod'.format(wtdset)
    else:
        cfile = '{}_ll.nod'.format(wtdset)

    # load file with WebTide coordinates
    df = pd.read_csv(os.path.join(wtdir, wtdset, cfile),
                     header=None,
                     delimiter='\s+',
                     usecols=[1, 2],
                     names=['lon', 'lat'],
                     engine='python')

    # load file with WebTide bathymetry and append to coordinates
    df['depth'] = pd.read_csv(os.path.join(wtdir, wtdset,
                                           '{}.bat'.format(wtdset)),
                              header=None,
                              delimiter='\s+',
                              usecols=[1],
                              names=['depth'],
                              engine='python')

    # initialize column for maximum amplitude if desired
    if ampmax:
        df['ampmax'] = np.zeros(len(df.index))

    # loop through constituents
    for name in constituents:

        # load file with WebTide coordinates
        tmp = pd.read_csv(os.path.join(wtdir, wtdset,
                                       '{}.barotropic.s2c'.format(name)),
                          skiprows=[0, 1, 2],
                          header=None,
                          delimiter='\s+',
                          usecols=[1, 2],
                          names=['amp', 'pha'],
                          engine='python')

        # write amplitude to output dictionary
        df[name+'_amp'] = tmp['amp']

        # write phase to output dictionary and map to +/- 180 degrees
        tmp.loc[tmp['pha'] > 180., 'pha'] = tmp['pha'] - 360.
        df[name+'_pha'] = tmp['pha']

        # compute maximum amplitude if desired
        if ampmax:
            df['ampmax'] = df['ampmax'] + tmp['amp']

    return df


if __name__ == '__main__':

    # run main method
    main()
