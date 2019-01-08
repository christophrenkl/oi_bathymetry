#!/usr/bin/env python

""" Extract observations.

"""
import dask.dataframe as dd
import numpy as np
import pandas as pd

from dask.distributed import Client

def main():

    client = Client(processes=False)

    # read dask.DataFrame with extracted data
    df = dd.read_parquet('data/interim/NWATL21_subset')

    # remove data with Med_depth < -10000
    df = df.loc[df['Med_depth'] > -1e4]

    # remove data with certain Group_Id which have been identified to have
    # inconsistent depth values
    df = df.loc[(df['Group_Id'] != 'ect18-38')      &
                (df['Group_Id'] != 'ch036l01')      &
                (df['Group_Id'] != 'c2207')         &
                (df['Group_Id'] != 'p885ns')        &
                (df['Group_Id'] != 'a2091l01')      &
                (df['Group_Id'] != 'kn151l4')       &
                (df['Group_Id'] != 'KJACK2006')     &
                (df['Group_Id'] != 'BROWNSBANK1996')]

    # make sure all depths are negative
    # df['Med_depth'] = -1 * df['Med_depth'].abs()

    # compute up to here, keep results in memory
    df = client.persist(df)

    # Some observations are referenced to LLWLT and need to be corrected.
    # Therefore, the dataset is split up into two subsets for further
    # processing. It is further assumed that the tidal correction does not
    # affect the accuracy for observations with Med_depth < -200.

    # only select coordinates and Med_depth
    outcols = ['Lon', 'Lat', 'Med_depth']

    # observations relative to MSL or MWL - no correction needed.
    dfmsl = df.loc[(df['Vertical_ref'] == 'MSL:2005')         |
                   (df['Vertical_ref'] == 'MSL:2006')         |
                   (df['Vertical_ref'] == 'MWL:2006')         |
                   (((df['Vertical_ref'] == 'LLWLT:2005')     |
                     (df['Vertical_ref'] == 'LLWLT:2006')     |
                     (df['Vertical_ref'] == 'VER_DAT:LLWLT')) &
                     (df['Med_depth'] < -200)), outcols]

    # observations relative to LLWLT - correction needed.
    dfllwlt = df.loc[((df['Vertical_ref'] == 'LLWLT:2005')     |
                      (df['Vertical_ref'] == 'LLWLT:2006')     |
                      (df['Vertical_ref'] == 'VER_DAT:LLWLT')) &
                     (df['Med_depth'] >= -200), outcols]

    # write output files
    dfmsl.to_parquet('data/interim/NWATL21_subset_msl')
    dfllwlt.to_parquet('data/interim/NWATL21_subset_llwlt')


if __name__ == '__main__':

    # run main method
    main()
