#!/usr/bin/env python

""" Create land-sea mask.

    https://stackoverflow.com/questions/53322952/creating-a-land-ocean-mask-in-cartopy
"""

import cartopy.io.shapereader as shpreader
import numpy as np
from shapely.geometry import box
import shapely.vectorized
import xarray as xr

def main():

    # read OI bathymetry
    ana = xr.open_dataset('data/interim/oi_bathymetry.nc')['ana'].stack(xy=('x', 'y'))

    # get coordinates
    lons = ana.lon.values
    lats = ana.lat.values

    # create bounding box of model domain
    bbox = box(lons.min(), lats.min(), lons.max(), lats.max())

    # read shapefile
    shpfilename = shpreader.natural_earth(resolution='50m',
                                          category='physical',
                                          name='land')
    reader = shpreader.Reader(shpfilename)

    # select all polygons that fall within model bounding box
    rec = list(reader.geometries())
    polys = [land for land in rec if land.intersects(bbox)]

    # get indices of grid points which are within each polygon
    inds = [np.where(shapely.vectorized.contains(poly, lons, lats)) for poly in polys]
    inds = np.unique(np.concatenate(inds, axis=1))

    # set land points to zero and multiply by -1 to make depth values positive
    ana[inds] = 0.
    ana = -1 * ana.unstack()

    # prepare output dataset
    dsout = ana.to_dataset(name='ana')

    # write to file
    dsout.to_netcdf('data/interim/oi_bathymetry_masked.nc')


if __name__ == '__main__':

    # run main method
    main()
