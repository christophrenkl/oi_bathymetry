#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import modules
import os
import xarray as xr
import pandas as pd
import numpy as np
import mpl_toolkits.basemap
from mpl_toolkits.basemap import Basemap
from matplotlib import path

import nemo_tools as nt

"""
Module create_bathymetry

Date: 2017-09-01

Description: Module to create a bathymetry file for a NEMO ocean model using an
             optimal interpolation technique. The code is mostly based on the
             original MATLAB code described in Lei et al. (2014).
             
             Ji Lei et al. "A procedure for processing in situ bathymetric
             observations for modeling Canadian shelf seas", DFO Technical
             Report, 2014.
"""

__author__   = 'Christoph Renkl'
__email__    = 'christoph.renkl@dal.ca'

# Define functions below. -----------------------------------------------------

def main():
    ''' '''

    # project directory
    projdir = '/home/chrenkl/Projects/nemo_bathymetry'
    
    # NEMO configuration
    conf = 'GoMSS'
    
    # WebTide data set
    wtset = 'nwatl'
    
    # tidal constituents for datum correction and minimum depth
    constituents = ['M2', 'S2', 'N2', 'K1', 'O1']
    
    # increment for limited area of domain on which OI is performed    
    delta = 150e3
    
    # OI parameters
    sigobs = 1.
    c = 10000. # [m]
    
    # parameters for minimum depth
    hmin = 5.
    ampfac = 1.2
    
    # NO MODIFICATIONS BELOW THIS LINE ########################################
    
    # STEP 0: read data =======================================================
    print('STEP 0: Data preparation. ========================================')
    
    # read NEMO grid
    NEMO = nt.nemo_grid(os.path.join(projdir, 'data/raw', conf,
                                     'Configuration',
                                     'coordinates.nc'))
    
    # read ETOPO data
    ETOPO = read_etopo(os.path.join(projdir, 'data/external/ETOPO1',
                                    'ETOPO1_Bed_g_gmt4.grd'),
                       lonmin=NEMO['lonmin'], lonmax=NEMO['lonmax'],
                       latmin=NEMO['latmin'], latmax=NEMO['latmax'])
    
    # read in-situ data
    # LLWLT, MSL = read_insitu(os.path.join(projdir, 'data/interim'))
    
    # read WebTide data
    WebTide = read_webtide(os.path.join(projdir, 'data/external/WebTide'),
                           wtset,
                           constituents=constituents,
                           ampmax=True)
            
    # create map projection for interpolations
    proj = Basemap(projection='merc',
                   llcrnrlat=NEMO['latmin'], urcrnrlat=NEMO['latmax'],
                   llcrnrlon=NEMO['lonmin'], urcrnrlon=NEMO['lonmax'],
                   resolution='h')
        
    # STEP 1: correct verical datum of in-situ data ===========================
    print('STEP 1: Vertical datum correction. ===============================')
    
    # file name of cleaned and corrected insitu data
    insitufile = os.path.join(projdir, 'data/interim/in-situ_corrected.h5')

    try:
        
        InSitu = pd.read_hdf(insitufile, 'InSitu')
        print('Found in-situ_corrected.h5. Skip vertical datum correction.')
        
    except:
        
        print('Correct vertical datum of in-situ data using WebTide.')
        print('NOT IMPLEMENTED YET!!!')
    
        ### IMPLEMENT DATUM CORRECTION HERE !!! ###
        
    # STEP 2: interpolate in-situ data to NEMO bathymetry =====================
    print('STEP 2: Gridding of in-situ data. ================================')
    
    # file name of cleaned and corrected insitu data
    ginsitufile = os.path.join(projdir, 'data/interim/gridded_in-situ.npz')
    
    try:
        
        gridInSitu = np.load(ginsitufile)
        print('Found gridded_in-situ.npz. Skip gridding.')
        
    except:
        
        print('Interpolate in-situ data to NEMO grid using nearest neighbour.')
        print('NOT IMPLEMENTED YET!!!')
    
        ### IMPLEMENT GRIDDING HERE !!! ###
           
    # STEP 3: create background field for OI (ETOPO on NEMO grid) =============    
    print('STEP 3: Create background field for OI. ==========================')
    
    # interpolate ETOPO to NEMO grid
    ETOPOint= etopo2nemo(NEMO, ETOPO)
           
    # STEP 4: optimal interpolation ===========================================    
    print('STEP 4: Optimal Interpolation. ===================================')
    
    ana = optimal_interpolation(ETOPOint, gridInSitu, proj, delta, c, sigobs)                                
                
    # STEP 5: create land-sea mask ============================================    
    print('STEP 5: Create land-sea mask. ====================================')
    
    # create land-sea mask
    mask = land_sea_mask(NEMO, proj)
    
    # apply land-sea mask and set land values to zero
    ana[(ana < 0) | (mask == 1.)] = 0.
               
    # STEP 6: set minimum depth (Maraldi et al., 2013) ========================  
    print('STEP 6: Set minimum depth. =======================================')
    
    ana = minimum_depth(ana, NEMO, WebTide, proj, hmin, ampfac)
        
    # STEP 7: remove closed seas ==============================================
    print('STEP 7: Remove closed seas. ======================================')
    
    
    ana = remove_closed_seas(ana)
    
    
def remove_closed_seas(ana):
    '''Remove closed seas in model domain. Copied from AGRIF nesting tools.'''
    
    print('Remove closed seas as in AGRIF nesting tools')
    
    # domain size
    nj, ni = np.shape(ana)
    
    # create dummy array
    btest = np.zeros((nj, ni))
    
    # treat domain borders
    btest[ 0,  :][ana[ 0,  :] > 0] = 1
    btest[-1,  :][ana[-1,  :] > 0] = 1
    btest[ :,  0][ana[ :,  0] > 0] = 1
    btest[ :, -1][ana[ :, -1] > 0] = 1
    
    nbadd = 1
    
    while nbadd != 0:
        
        nbadd = 0
        
        for jj in np.arange(1, nj-1):
            for ii in np.arange(1, ni-1):

                if ana[jj, ii] > 0.:
                    
                    if np.maximum.reduce([btest[jj, ii+1], 
                                          btest[jj, ii-1],
                                          btest[jj+1, ii],
                                          btest[jj-1, ii]]) == 1:
                           
                           if btest[jj, ii] != 1:
                               
                               nbadd = nbadd + 1
                               
                           btest[jj, ii] = 1
    
    ana[btest == 0] = 0.

    return ana
    
    
def minimum_depth(ana, NEMO, WebTide, proj, hmin=5., ampfac=1.2):
    '''Set minimum depth depending on according to Maraldi et al. (2013) using
       tidal amplitudes from WebTide.'''
    
    print('Set minimum depth depending on according to Maraldi et al. (2013).')
        
    # NEMO grid dimension
    ndata = NEMO.dims['x'] * NEMO.dims['y']
    
    # WebTide has wetting and drying, don't use grid points shallower than 10m    
    WebTide = WebTide[WebTide.depth >= 10]
    
    # target coordinates: NEMO
    xm, ym = proj(NEMO.glamt.values.ravel(), NEMO.gphit.values.ravel())
    
    # source coordinates: WebTide
    xwt, ywt = proj(WebTide.lon.values, WebTide.lat.values)

    # initialization
    tmp = np.zeros(ndata)
    
    # loop through nemo grid
    for ipt in np.arange(ndata):
        
        # only correct wet grid points
        if ana.ravel()[ipt] > 0.:
        
            # index of minimum distance to each WebTide node
            mind = WebTide.index[np.abs(xm[ipt] - xwt + 1j * (ym[ipt] - ywt)).argmin()]
            
            # compute median of in-situ observations closest to grid point
            tmp[ipt] = np.maximum(hmin, ampfac * WebTide.ampmax[mind])
                
    # reshape
    mindepth = tmp.reshape((NEMO.dims['y'], NEMO.dims['x']))
    
    # replace shallow grid points with minimum depth
    ana = np.maximum(ana, mindepth)
    
    return ana
    
    
def land_sea_mask(NEMO, proj):
    '''Create land-sea mask based on coastline coordinates.'''
    
    # NEMO grid coordinates
    xgrd, ygrd = proj(NEMO.glamt.values.ravel(), NEMO.gphit.values.ravel())
    
    # get coastline coordinates
    cl = proj.drawcoastlines()
    coords = cl.get_segments()
    
    mask = np.zeros(NEMO.dims['y'] * NEMO.dims['x'])    
    
    for ipath in np.arange(len(coords)):
            
        # we need coordinates as list of tuples
        tmp = list(map(tuple, coords[ipath]))
        
        # add a corner to avoid funny artefacts
        # tmp.append((tmp[0][0], tmp[-1][-1]))
        tmp.append((coords[ipath][:,0].min(), coords[ipath][:,1].max()))
        
        
        codes = [path.Path.MOVETO] + [path.Path.LINETO]*(len(tmp)-1)
    
        # create a path object of coastline
        cp = path.Path(tmp, codes)
        
        test = np.hstack((xgrd[:, np.newaxis], ygrd[:, np.newaxis]))
    
        # find all NEMO grid points within path
        mask = mask + cp.contains_points(test)
    
    
    mask = mask.reshape((NEMO.dims['y'], NEMO.dims['x']))

    return mask


def optimal_interpolation(bckgrnd, obs, proj, delta, c, sigobs):
    '''Optimal interpolation routine.'''
    
    # initialize analysis field with background
    ana = np.copy(bckgrnd['depth'])   # analysis

    # innovation (observations minus background)
    innov = obs['depth'] - bckgrnd['depth']
    innov[innov == -bckgrnd['depth']] = 0.

    # projected coordinates: in-situ grid
    xx, yy = proj(bckgrnd['lon'], bckgrnd['lon'])
    
    # loop through limited areas of domain
    for x1 in np.arange(xx.min()-c, xx.max()+delta+c, delta):
        for y1 in np.arange(yy.min()-c, yy.max()+delta+c, delta):
                
            # upper bounds of limited area
            x2 = x1 + delta + 2 * c
            y2 = y1 + delta + 2 * c

            # initialise temporary analysis field which will include overlap
            anatmp = np.zeros_like(ana)
            
            # select background grid indices in bounding box
            imod = np.where((xx >= x1) & (xx <  x2) & (yy >= y1) & (yy <  y2))
            
            # select innovation grid indices in bounding box
            iobs = np.where((xx >= x1) & (xx <  x2) & (yy >= y1) & (yy <  y2) &
                            (innov != 0.))
            
            # select domain grid indices in bounding box without overlap 
            imoda = np.where((xx >= x1+c) & (xx <  x2-c) &
                             (yy >= y1+c) & (yy <  y2-c))
        
            # select grid points within selected area (including overlap)
            x1m = xx[imod]; x2m = yy[imod]
            x1o = xx[iobs]; x2o = yy[iobs]

            # select innovation (y - Hx) within selected area
            ymHx = innov[iobs]
            
            # number of observations
            nobs = len(ymHx)
            
            # ony apply optimal interpolation if observations exist
            if nobs > 0 and len(x1m) > 0:
              
                # weighting: number of observations that went into super obs
                nsobs = obs['nobs'][iobs]
                
#                # optimal interpolations using Gaspari-Cohn correlation
#                pred = oi_gaspari_cohn(x1m, x2m, x1o, x2o, ymHx,
#                                       c, nsobs, sigobs)
#        
#                # temporary analysis including overlap
#                anatmp[imod] = bckgrnd['depth'][imod] + pred
                
                # analysis without overlap
                ana[imoda] = anatmp[imoda]    

    return ana


def oi_gaspari_cohn(x1m, x2m, x1o, x2o, ymHx, c, nsobs, sigobs):
    '''Optimal interpolation using Gaspari-Cohn correlation function.'''
    
    print('Optimal interpolation using Gaspari-Cohn correlation function.')
    
    # number of observations
    nobs = len(x1o)
    
    # R-matrix: observation error covariance matrix
    R = np.diag(sigobs / nsobs)
        
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
    Q = np.ones(len(x1m))
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
    tmp = np.linalg.solve((H.dot(B).dot(H.T) + R), ymHx)
    pred = B.dot(H.T).dot(tmp)
    
    return pred


def etopo2nemo(NEMO, ETOPO):
    '''Interpolate ETOPO data to NEMO grid.'''
    
    print('Interpolate ETOPO data to NEMO grid.')
    
    # initialization
    ETOPOint = {'lon': NEMO.glamt.values,
                'lat': NEMO.gphit.values}
    
    ETOPOint['depth'] = mpl_toolkits.basemap.interp(-ETOPO.z.values,
                                                    ETOPO.x.values,
                                                    ETOPO.y.values,
                                                    NEMO.glamt.values,
                                                    NEMO.gphit.values,
                                                    order=1)
    
    return ETOPOint


def read_webtide(indir, wtset, constituents=['M2', 'S2', 'N2', 'K1', 'O1'],
                 ampmax=False):
    '''Read WebTide data and compute maximum amplitude if desired.'''
    
    print('Read WebTide data.')
    
    # load file with WebTide coordinates
    WebTide = pd.read_csv(os.path.join(indir, wtset, wtset+'_ll.nod'),
                                       header=None,
                                       delimiter='\s+',
                                       usecols=[1, 2],
                                       names=['lon', 'lat'],
                                       engine='python')
    
    # load file with WebTide bathymetry and append to coordinates
    WebTide['depth'] = pd.read_csv(os.path.join(indir, wtset, wtset+'.bat'),
                                     header=None,
                                     delimiter='\s+',
                                     usecols=[1],
                                     names=['depth'],
                                     engine='python')
    
    # initialize column for maximum amplitude if desired
    if ampmax:
        WebTide['ampmax'] = np.zeros(len(WebTide.index))
    
    # loop through constituents
    for cst in constituents:

        # load file with WebTide coordinates
        tmp = pd.read_csv(os.path.join(indir, wtset, cst+'.barotropic.s2c'),
                                       skiprows=[0, 1, 2],
                                       header=None,
                                       delimiter='\s+',
                                       usecols=[1, 2],
                                       names=['amp', 'pha'],
                                       engine='python')
        
        # write amplitude to output dictionary
        WebTide[cst+'_amp'] = tmp['amp']
        
        # write phase to output dictionary and map to +/- 180 degrees
        tmp['pha'][tmp['pha'] > 180.] = tmp['pha'][tmp['pha'] > 180.] - 360.
        WebTide[cst+'_pha'] = tmp['pha']
        
        # compute maximum amplitude if desired
        if ampmax:
            WebTide['ampmax'] = WebTide['ampmax'] + tmp['amp']
    
    return WebTide


def read_insitu(indir):
    '''Read in-situ data and sort by vertical datum.'''
    
    print('Read in-situ data.')
    
    # in-situ observations with respect to lower low water large tides
    LLWLT = pd.read_csv(os.path.join(indir, 'LLWLTData.lld'),
                                     header=None,
                                     delimiter='\s+',
                                     names=['lon', 'lat', 'depth'],
                                     engine='python')
    
    # in-situ observations with respect to mean sea level
    MSL = pd.read_csv(os.path.join(indir, 'MSLData.lld'),
                                   header=None,
                                   delimiter='\s+',
                                   names=['lon', 'lat', 'depth'],
                                   engine='python')
    
    return (LLWLT, MSL)
    
    
def read_etopo(fname, lonmin=-180., lonmax=180., latmin=-90., latmax=90.):
    '''Read ETOPO data.'''
    
    print('Read ETOPO data.')
    
    # load ETOPO file
    ETOPO = xr.open_dataset(fname)
    
    # subset data
    ETOPO = ETOPO.sel(x=slice(lonmin, lonmax), y=slice(latmin, latmax))
    
    return ETOPO
    

if __name__ == '__main__':
    
    # run main method
    main()
