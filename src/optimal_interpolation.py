#!/usr/bin/env python

""" Optimal interplation of background and observations.

"""
import dask.array as da
import numpy as np
from scipy.spatial import cKDTree
import xarray as xr

import dask

from tools import ll2xyz

import matplotlib.pyplot as plt


def main():

    # read background field and gridded observations
    bg = xr.open_dataset('data/interim/GEBCO_background.nc')['ele']
    obs = xr.open_dataset('data/interim/NWATL21_subset_gridded.nc')['depth']

    # reshape 2D fields to vectors
    bg = bg.stack(z=('x', 'y'))
    obs = obs.stack(z=('x', 'y'))

    # drop all grid points without observations
    obs =  obs.where(obs != 0., drop=True)

    # number of grid points with observations
    nobs = len(obs.z)

    # convert lon/lat to Cartesian coordinates
    xbg, ybg, zbg = ll2xyz(bg.lon.values, bg.lat.values)
    xo, yo, zo = ll2xyz(obs.lon.values, obs.lat.values)

    # create cKDTree object to represent source grid
    tree = cKDTree(list(zip(xbg, ybg, zbg)))

    # find indices of the nearest neighbor
    _, ind = tree.query(list(zip(xo, yo, zo)), k=1, n_jobs=-1)

    # create H-matrix with dimension [nobs, len(bg.z)]
    H = h_matrix(ind, len(bg.z))

    print(H)


def h_matrix(ind, ncols):

    create_row_lazy = dask.delayed(create_row, pure=True)

    lazy_rows = [create_row_lazy(ncols, ii) for ii in ind]

    sample = lazy_rows[0].compute()  # load the first image (assume rest are same shape/dtype)

    rows = [da.from_delayed(lazy_row,           # Construct a small Dask array
                            dtype=sample.dtype,   # for every lazy value
                            shape=sample.shape) for lazy_row in lazy_rows]

    return da.stack(rows, axis=0).rechunk(chunks=(1000, 1000))


def create_row(length, ind):
    row = np.zeros((length))
    row[ind] = 1
    return row




    # # compute innovation (y - Hx)
    # yHx = (yo - xb)# .fillna(value=0.)
    #
    # print(yHx[np.where((yHx != 0.))])
    #
    # # yHx.plot(cmap='viridis')
    # # plt.show()
    # #

def oi_gaspari_cohn(x1m, x2m, x1o, x2o, ymHx, c, nsobs, sigobs, nstar):
    '''Optimal interpolation using Gaspari-Cohn correlation function.'''

    # number of observations
    nobs = len(x1o)

    # R-matrix: observation error covariance matrix
    R = np.diag(sigobs / nsobs + sigobs / nstar)

    # initialize H-matrix
    H = np.zeros((nobs, len(x1m)))

    # fill H-matrix at indices of nearest neighbour
    for irow in np.arange(nobs):

        z = x1o[irow] - x1m + 1j * (x2o[irow] - x2m)
        ind = np.abs(z).argmin()
        H[irow, ind] = 1.

    # compute intra-grid distances
    z = x1m + 1j * x2m
    Z = np.array([z,] * len(z))
    igdist = np.abs(Z - Z.T)

    # Q-matrix: background error
    Q = 100000. * np.ones(len(x1m))
    #Q[bckgrnd[imod] < 10.] = 2.
    Q = np.diag(Q)

    # B-matrix: covariance matrix
    B = np.zeros_like(igdist)
    i1 = np.where(igdist < c)
    i2 = np.where((igdist >= c) & (igdist < 2 * c))
    z = igdist[i1] / c

    B[i1] = -z**5 / 4 + z**4 /2 + 5/8 * z**3 - 5/3 * z**2 + 1

    z = igdist[i2] / c

    B[i2] = z**5 / 12 - z**4 /2 + 5/8 * z**3 + 5/3 * z**2 - 5*z \
            + 4 - 2/3 * 1/z

    # weighting with background error
    B = Q.dot(B).dot(Q)

    # optimal interpolation
    I = np.eye(len(x1m))
    W = B.dot(H.T).dot(np.linalg.inv((H.dot(B).dot(H.T) + R)))

    pred = W.dot(ymHx)
    Pa = np.diag((I - W.dot(H)).dot(B)) / np.diag(B)

    return pred, Pa

if __name__ == '__main__':

    # run main method
    main()
