#!/usr/bin/env python

""" Create land-sea mask.

    https://stackoverflow.com/questions/53322952/creating-a-land-ocean-mask-in-cartopy
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import numpy as np
from shapely.geometry import Point
from shapely.prepared import prep
import xarray as xr

import matplotlib.pyplot as plt

def main():

    ana = xr.open_dataset('data/interim/oi_bathymetry.nc')['ana'].stack(xy=('x', 'y'))

    lons = ana.lon.values
    lats = ana.lat.values

    # load the shapefile geometries
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '50m')
    land_polygons = list(land_10m.geometries())

    # turn the model coordinates into shapely Points
    points = [Point(point) for point in zip(lons, lats)]

    # prepare land polygons for speed-up
    land_polygons_prep = [prep(land_polygon) for land_polygon in land_polygons]

    # determine the model grid points that fall on land
    land = []
    for land_polygon in land_polygons_prep:
        land.extend([tuple(point.coords)[0] for point in filter(land_polygon.covers, points)])

    # find indeces of land points
    lon_land, lat_land = zip(*land)
    _, ind, _ = np.intersect1d(lons, lon_land, return_indices=True)

    # set land points to zero and multiply by -1 to make depth values positive
    ana[ind] = 0.
    ana = -1 * ana.unstack()

    # prepare output dataset
    dsout = ana.to_dataset(name='ana')

    # write to file
    dsout.to_netcdf('data/interim/oi_bathymetry_masked.nc')


if __name__ == '__main__':

    # run main method
    main()
