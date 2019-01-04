#!/usr/bin/env python

""" tools

A collection of functions to create a gridded bathymetry.

"""

import xarray as xr

def read_grid(fname):
    """Read output grid"""

    grd = xr.open_dataset(fname)

    return grd
