#!/usr/bin/env python

""" tools

A collection of functions to create a gridded bathymetry.

"""

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
