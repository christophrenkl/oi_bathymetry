#!/usr/bin/env python

""" Prepare output grid.

"""
import tools

def main():

    grd = tools.read_grid('data/raw/coordinates.nc')

    print('lonmin = ', grd.glamt.min())
    print('lonmax = ', grd.glamt.max())
    print('latmin = ', grd.gphit.min())
    print('latmax = ', grd.gphit.max())



if __name__ == '__main__':

    # run main method
    main()
