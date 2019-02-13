#!/usr/bin/env python

""" tools

A collection of functions to create a gridded bathymetry.

"""
import numpy as np
import os
import pandas as pd
import xarray as xr

def read_grid(fname):
    """Read output grid"""

    grd = xr.open_dataset(fname)

    return grd


def ll2xyz(lon, lat, R=6371000):
    """
    Calculates coordinates of a point on a sphere with radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)

    return x,y,z


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
