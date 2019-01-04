#!/usr/bin/env python

""" Extract observations.

"""
import dask.dataframe as dd
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

def main():

    # domain boundaries - hard coded for now
    lonmin = -71.923204
    lonmax = -51.113073
    latmin = 37.567801
    latmax = 48.848009

    # data types of columns
    dtypes = {'Year': float,
              'Month': float,
              'Day': float,
              'Hr': float,
              'Min': float,
              'Sec': float,
              'Lon': float,
              'Lon': float,
              'Lat': float,
              'Obj_type': str,
              'Med_depth': float,
              'Agency': str,
              'Depth_type': str,
              'Vertical_ref': str,
              'Group_Id': str}

    # read file with observations as Dask.DataFrame
    df = dd.read_csv('/media/chrenkl/external/NWATL21_out.txt',
                     sep='\s+',
                     na_values = {'Day' :      'NORTH2008',
                                  'Hr':        'NORTH2008',
                                  'Min':       'NORTH2008',
                                  'Lon':       '-157.6566-1436.630000',
                                  'Lat':       'UNH',
                                  'Med_depth': ['MSL:2006',
                                                '-2400.01TIBEAM',
                                                'CHS']},
                     dtype=dtypes,
                     assume_missing=True,
                     error_bad_lines=False)
    # Note that the last argument drops all lines which have more columns than
    # suggested by the header. A warning will be printed to the terminal.

    # preprocessing ---

    # remove incomplete data
    df = df.dropna()

    # extract data within given domain
    df = df.loc[(df['Lon'] >= lonmin) &
                (df['Lon'] <= lonmax) &
                (df['Lat'] >= latmin) &
                (df['Lat'] <= latmax)]

    # remove data with Med_depth < -10000
    df = df.loc[df['Med_depth'] > -1e4]

    # remove data from gridded and interpolated products
    df = df.loc[(df['Agency'] != 'GMT')    |
                (df['Agency'] != 'NETCDF') |
                (df['Agency'] != 'NGDC')   |
                (df['Agency'] != 'IBCAO')  |
                (df['Agency'] != 'ABC')    |
                (df['Agency'] != 'CHSQUE') |
                (df['Agency'] != 'NSTOPO') |
                (df['Agency'] != 'INTERP')]

    # remove data with certain Group_Id which have been identified to have
    # inconsistent depth values
    df = df.loc[(df['Group_Id'] != 'ect18-38')      |
                (df['Group_Id'] != 'ch036l01')      |
                (df['Group_Id'] != 'c2207')         |
                (df['Group_Id'] != 'p885ns')        |
                (df['Group_Id'] != 'a2091l01')      |
                (df['Group_Id'] != 'kn151l4')       |
                (df['Group_Id'] != 'KJACK2006')     |
                (df['Group_Id'] != 'BROWNSBANK1996')]


    # remove data with inappropriate reference level
    # df = df.loc[(df[''])]
    df = df.compute()

    plt.scatter(df['Lon'], df['Lat'],
                c=df['Med_depth'])

    plt.show()

    # print(df.head())



if __name__ == '__main__':

    # run main method
    main()
