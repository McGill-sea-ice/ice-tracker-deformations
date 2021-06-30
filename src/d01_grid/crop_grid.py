'''
Author: Beatrice Duval (bdu002)

--------------------------------------------------------------------
Code for creation of a netcdf grid, adapted from original RIOPS grid
--------------------------------------------------------------------

Crops the RIOPS grid to region of interest. 
Creates a new netcdf cropped grid containing tracer points (lat,lon and 
    x,y coordinates in the aeqd transform), speed points, distance between 
    tracer points and land, and a sea-land mask.
'''

import os

import numpy as np
import pandas as pd
from haversine import haversine
from netCDF4 import Dataset
from pyproj import Proj
from math import nan
from src.d00_utils.grid_coord_system import find_nearestGridTracerPt

################################### MAIN PROGRAM #################################

#-------------------------Trim NEMO grid-----------------------------------------

# Load netcdf NEMO dataset
gridPath = '/home/socn000/env_ubuntu-18.04-skylake-64/eccc-ppp4/datafiles/constants/oce/repository/master/CONCEPTS/creg012pe/grid/coordinates_CREG12_ext.nc'
grid_ds = Dataset(gridPath)

# Set values at which the grid is to be cropped
y1 = 1075
y2 = 1750
x1 = 375
x2 = 1200

# Create a cropped matrix for grid latitudes 
# and longitudes in the region of interest
LAT = grid_ds['nav_lat'][y1:y2, x1:x2] # 'Point traceur' t
LON = grid_ds['nav_lon'][y1:y2, x1:x2]

fLAT =grid_ds['gphif'][y1:y2, x1:x2]   # 'Point vitesse f'
fLON = grid_ds['glamf'][y1:y2, x1:x2]

# Retrieve new size for grid matrices
ydim, xdim = LAT.shape

#-------------Create 'points traceurs' matrix in x,y coordinates-----------------

# Convert grid points from lat/lon to x,y coordinates following 
# the Azimuthal Equidistant transform
p = Proj(proj='aeqd', ellps='WGS84', preserve_units=False)
X, Y = p(LON, LAT)

#---------------------Create land-sea mask -------------------------------------

# Load bathymetry dataset
bathyPath = '/home/socn000/env_ubuntu-18.04-skylake-64/eccc-ppp4/datafiles/constants/oce/repository/master/CONCEPTS/creg012pe/grid/bathy_creg12_ext_mask.nc'
bathy_ds = Dataset(bathyPath)

# Create a cropped matrix from bathymetry 
# of same size as LAT and LON
MASK = bathy_ds['Bathymetry'][y1:y2, x1:x2]

# Set all land values (>0) to 0 and ocean values (=0)
# to 1 
MASK[ MASK[:,:] > 0 ] = 1
MASK[ MASK[:,:] == 0 ] = 0


#-----------------Create distance to coastline matrix --------------------------
    
# Initialize distances matrix to a zero matrix
DIST = np.zeros( (ydim, xdim) )

# Create a masked 'points traceurs' matrix in x,y coordinates
# where all ocean points are set to nan
maskedX = X.copy()
maskedY = Y.copy()

maskedX[ MASK[:,:] == 1 ] = nan
maskedY[ MASK[:,:] == 1 ] = nan

# Iterate through all grid points
for y in range(ydim) :
    for x in range(xdim):
        # If the point is located on sea, find its minimal 
        # distance to land
        # Else, the distance to land is 0
        if MASK[y,x] == 1:
            j, i = find_nearestGridTracerPt(maskedX, maskedY, (X[y,x], Y[y,x]))
            # Add the minimal distance found to the 
            # distances to coastline matrix
            DIST[y,x] = haversine( (LAT[y,x], LON[y,x]), (LAT[j,i], LON[j,i]) )

#---------------------Create new netcdf grid-------------------------------------

# Find absolute path in which cropped_grid.nc is to be stored
projPath = os.path.dirname(os.path.realpath(__file__))
cropGridPath = projPath + '/../../data/00_grid/cropped_grid.nc'

# Create cropped_grid.nc netcdf file and associated crop_ds dataset
crop_ds = Dataset(cropGridPath, 'w', format = 'NETCDF4')

# Create (y, x) matrix to store LAT, LON, Y, X, DIST, MASK and fLAT, fLON
y = crop_ds.createDimension('y', ydim)
x = crop_ds.createDimension('x', xdim)

# Create variables for netcdf data set
nav_lat = crop_ds.createVariable('nav_lat', 'f8', ('y', 'x')) # 'Point traceur' (lat)
nav_lon = crop_ds.createVariable('nav_lon', 'f8', ('y', 'x')) # 'Point traceur' (lon)
nav_y   = crop_ds.createVariable('nav_y', 'f8', ('y', 'x'))   # 'Point traceur' (y/lat)
nav_x   = crop_ds.createVariable('nav_x', 'f8', ('y', 'x'))   # 'Point traceur' (x/lon)
gphif   = crop_ds.createVariable('gphif', 'f8', ('y', 'x'))   # 'Point vitesse' (lat)
glamf   = crop_ds.createVariable('glamf', 'f8', ('y', 'x'))   # 'Point vitesse' (lon)
mask    = crop_ds.createVariable('mask', 'i1', ('y', 'x'))    # Land-ocean mask
dist    = crop_ds.createVariable('dist', 'f8', ('y', 'x'))    # Ocean-coastline distances

# Specify units
nav_lat.units = 'degrees'
nav_lon.units = 'degrees'
gphif.units   = 'degrees'
glamf.units   = 'degrees'
dist.units    = 'meters'

# Attribute data matrices to variables
nav_lat[:, :] = LAT
nav_lon[:,:]  = LON
nav_y[:, :]   = Y
nav_x[:,:]    = X
gphif[:, :]   = fLAT
glamf[:,:]    = fLON
mask[:, :]    = MASK
dist[:,:]     = DIST

# Close dataset
crop_ds.close()
