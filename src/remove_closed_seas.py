#!/usr/bin/env python

""" Remove closed seas.

"""
import numpy as np
import xarray as xr

def main():

    # read masked and depth-corrected bathymetry
    ana = xr.open_dataset('data/interim/oi_bathymetry_depth_corrected.nc')['ana']

    nj, ni = ana.shape

    # create dummy array
    btest = np.zeros((nj, ni))

    # treat domain borders
    btest[ 0,  :][ana[ 0,  :] > 0] = 1
    btest[-1,  :][ana[-1,  :] > 0] = 1
    btest[ :,  0][ana[ :,  0] > 0] = 1
    btest[ :, -1][ana[ :, -1] > 0] = 1

    # set flag
    nbadd = 1

    while nbadd != 0:

        nbadd = 0

        # iterate throug model domain
        it = np.nditer(ana[1:-1, 1:-1], flags=['multi_index'])
        while not it.finished:

            # get multi_index of iterator
            jj, ii = it.multi_index
            jj = jj + 1
            ii = ii + 1

            # if bathymetry has a wet point
            if it[0] > 0.:

                # check if any of the surrounding grid points are land
                if np.maximum.reduce([btest[jj, ii+1],
                                      btest[jj, ii-1],
                                      btest[jj+1, ii],
                                      btest[jj-1, ii]]) == 1:

                       if btest[jj, ii] != 1:
                           nbadd = nbadd + 1

                       # mark grid point as land
                       btest[jj, ii] = 1

            # move to next grid point
            it.iternext()

    # apply new land mask
    ana.values[btest == 0] = 0.

    # prepare output dataset and write to file
    dsout = ana.to_dataset(name='ana')
    dsout.to_netcdf('data/interim/oi_bathymetry_final.nc')


if __name__ == '__main__':

    # run main method
    main()
