#!/usr/bin/env python

""" Extract observations.

"""
import dask.dataframe as dd

def main():

    # domain boundaries - hard coded for now
    lonmin = -71.923204
    lonmax = -51.113073
    latmin = 37.567801
    latmax = 48.848009

<<<<<<< HEAD
=======
    # fname = '/media/chrenkl/external/NWATL21_out.txt'
>>>>>>> 58109c156795b42742054199248e0cd1da53b549
    fname = 'data/raw/NWATL21_out.txt'

    # client = Client(processes=False)

    # data types of columns
    dtypes = {'Lon': float,
              'Lat': float,
              'Obj_type': str,
              'Med_depth': float,
              'Agency': str,
              'Depth_type': str,
              'Vertical_ref': str,
              'Group_Id': str}

    # read file with observations as Dask.DataFrame
    df = dd.read_csv(fname,
                     sep='\s+',
                     usecols=['Lon', 'Lat', 'Med_depth',
                              'Agency', 'Group_Id', 'Vertical_ref'],
                     na_values = {'Lon':       '-157.6566-1436.630000',
                                  'Lat':       'UNH',
                                  'Med_depth': ['MSL:2006',
                                                '-2400.01TIBEAM',
                                                'CHS']},
                     dtype=dtypes,
                     assume_missing=True,
                     error_bad_lines=False)
    # Note that the last argument drops all lines which have more columns than
    # suggested by the header. A warning will be printed to the terminal.

    # data extraction ----------------------------------------------------------

    # remove incomplete data
    df = df.dropna()

    # extract data within given domain
    df = df.loc[(df['Lon'] >= lonmin) &
                (df['Lon'] <= lonmax) &
                (df['Lat'] >= latmin) &
                (df['Lat'] <= latmax)]

    # remove data from gridded and interpolated products
    df = df.loc[(df['Agency'] != 'GMT')    &
                (df['Agency'] != 'NETCDF') &
                (df['Agency'] != 'NGDC')   &
                (df['Agency'] != 'IBCAO')  &
                (df['Agency'] != 'ABC')    &
                (df['Agency'] != 'CHSQUE') &
                (df['Agency'] != 'NSTOPO') &
                (df['Agency'] != 'INTERP')]

    # drop Agency column
    df = df.drop('Agency', axis=1)

    # write output file
    df.to_parquet('data/interim/NWATL21_subset')

if __name__ == '__main__':

    # run main method
    main()
