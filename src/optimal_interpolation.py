#!/usr/bin/env python

""" Optimal interplation of background and observations.

"""
import numpy as np
from numpy.linalg import inv
import xarray as xr

from dask.distributed import Client

def main():

    client = Client()

    # parameters - hard-coded for now
    c = 1e4 # [m]
    sigobs = 1.
    N = 60

    # read background field and gridded observations
    bg = xr.open_dataset('data/interim/GEBCO_background.nc')['ele']
    obs = xr.open_dataset('data/interim/NWATL21_subset_gridded.nc')['depth']

    # create data array for optimally interpolated bathymetry
    ana = xr.DataArray(np.ones_like(bg) * np.nan)

    # Since working on the whole model domain can be computationally expensive,
    # the procedure is to select sliding sub-domains to apply the OI scheme to.
    # For the Gaspari-Cohn localization, information outside the sub-domains is
    # needed, so we have do extend the subdomains by enough grid points to have
    # all observations within a radius of two times the length scale c included.

    # calulate grid spacing in x- and y-direction, choose the maximum value
    dx = haversine(bg.lon.isel({'x': 0, 'y': 0}), bg.lat.isel({'x': 0, 'y': 0}),
                   bg.lon.isel({'x': 1, 'y': 0}), bg.lat.isel({'x': 1, 'y': 0}))
    dy = haversine(bg.lon.isel({'x': 0, 'y': 0}), bg.lat.isel({'x': 0, 'y': 0}),
                   bg.lon.isel({'x': 0, 'y': 1}), bg.lat.isel({'x': 0, 'y': 1}))
    delta = np.max((dx, dy))

    # number of ghost cells for each sub-domain extension
    ng = (c // delta + 1).astype(int)

    # get the size of the whole model domain to determine number of sub-domains
    # which should contain NxN grid points
    nx, ny = bg.shape
    nx = nx // N + 1
    ny = ny // N + 1

    # subdomain counter
    dcount = 0

    # itereate through sub-domains
    for ix in range(nx):
        # start and end indices of sub-domain in i-direction
        i1 = ix * N
        i2 = min((ix+1) * N,  bg.shape[0])

        # start and end indices of extended sub-domain in i-direction
        i1e = max(0, i1-ng)
        i2e = min(bg.shape[0], i2+ng)

        for iy in range(ny):

            dcount = dcount + 1
            print('Process sub-domain {}/{}'.format(dcount, nx*ny))

            # start and end indices of sub-domain in j-direction
            j1 = iy * N
            j2 = min((iy+1) * N,  bg.shape[1])

            # start and end indices of extended sub-domain in j-direction
            j1e = max(0, j1-ng)
            j2e = min(bg.shape[1], j2+ng)

            # select sub-domain and reshape 2D fields to vectors
            xb = bg.isel(x=slice(i1e, i2e), y=slice(j1e, j2e)).stack(xy=('x', 'y'))
            xo = obs.isel(x=slice(i1e, i2e), y=slice(j1e, j2e)).stack(xy=('x', 'y'))

            # drop all grid points without observations
            xo = xo.where(xo != 0., drop=True)

            # observation error covariance matrix
            R = np.eye(len(xo.xy))

            # Gaspari-Cohn localization matrix
            B = gaspari_cohn_matrix(xb.lon.values, xb.lat.values, c)

            # H operator - maps model and observation grid points
            H = h_operator(xb, xo)

            # compute weight matrix
            S = H.dot(B).dot(H.T) + R
            W = B.dot(H.T).dot(inv(S))

            # optimal interpolation
            xa = xb + W.dot(xo - H.dot(xb))
            xa = xa.unstack()

            # write to output dataset
            ana[slice(i1, i2), slice(j1, j2)] = xa[i1-i1e:i1-i1e+N, j1-j1e:j1-j1e+N]

    # create output dataset
    dsout = xr.Dataset(data_vars={'ana': (['x', 'y'], ana)},
                       coords={'lon': (['x', 'y'], bg.lon),
                               'lat': (['x', 'y'], bg.lat)})

    # write to file
    dsout.to_netcdf('data/interim/oi_bathymetry.nc')


def h_operator(xb, xo):

    # indeces of grid points where observations exist
    _, ind, _ = np.intersect1d(xb.xy, xo.xy, return_indices=True)

    # number of grid points in background and observations
    nxb = len(xb.xy)
    nxo = len(xo.xy)

    # create H operator
    H = np.zeros((nxo, nxb))
    H[range(nxo), ind]  = 1.

    return H


def gaspari_cohn_matrix(xblon, xblat, c):

    # create matrix with all combinations of grid points
    tmp = xblon + 1j * xblat
    tmp = np.array([tmp,] * len(tmp))
    a1 = tmp.flatten()
    a2 = tmp.T.flatten()

    # compute intra-grid distances
    d = haversine(np.real(a1), np.imag(a1), np.real(a2), np.imag(a2))
    d = d.reshape(tmp.shape)

    # clean up some memory
    del tmp, a1, a2

    B = np.zeros_like(d)
    z1 = d[d < c] / c
    B[d < c] = -z1**5 / 4 + z1**4 /2 + 5/8 * z1**3 - 5/3 * z1**2 + 1

    z2 = d[(d >= c) & (d < 2 * c)] / c
    B[(d >= c) & (d < 2 * c)] = z2**5 / 12 - z2**4 /2 + 5/8 * z2**3 \
                                + 5/3 * z2**2 - 5*z2 + 4 - 2/3 * 1/z2

    return B


def haversine(lon1, lat1, lon2, lat2):
    """ """

    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    d = 6371e3 * c

    return d


if __name__ == '__main__':

    # run main method
    main()
