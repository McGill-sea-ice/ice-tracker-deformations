'''
Author: Beatrice Duval (bdu002)

--------------------------------------------------------------------
Code for creation of a netcdf grid, adapted from original RIOPS grid
--------------------------------------------------------------------

Crops the RIOPS grid to region of interest. 
Creates a new netcdf cropped grid containing central 
    grid points ('points traceurs'), f speed points 
    ('points de vitesse'), distance between cental 
    points and land, and a sea-land mask.
'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from haversine import haversine
import os

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
LAT = grid_ds['nav_lat'][y1:y2, x1:x2] # 'Point traceur t'
LON = grid_ds['nav_lon'][y1:y2, x1:x2]

fLAT =grid_ds['gphif'][y1:y2, x1:x2]   # 'Point vitesse f'
fLON = grid_ds['glamf'][y1:y2, x1:x2]

# Retrieve new size for grid matrices
ydim, xdim = LAT.shape

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


#-----------------Create list of coastline (lat, lon) points---------------------

# Initialize coastline list
COAST = []

# Iterate through all grid points and add coastline points to the list
for y in range(ydim):
    for x in range(xdim):
        
        # A land point (0) is on the coastline if an ocean
        # point (1) is somewhere around it, for example:
        #
        #               1    1    1
        #               1    0    1   
        #               1    1    1
        
        isCoast = False
        
        for i in [-1, 0, 1]: 
            for j in [-1, 0, 1]: 
                # Check if location is within boundary
                if 0 <= x + i < xdim and 0 <= y + j < ydim:
                    isCoast = isCoast or MASK[y+j,x+i] == 1

        if MASK[y,x] == 0 and isCoast:
            COAST.append((LAT[y,x], LON[y,x]))


#-----------------Create distance to coastline matrix --------------------------
    

# Initialize distances matrix to a zero matrix
DIST = np.zeros( (ydim, xdim) )

# Iterate through all grid points
for y in range(ydim):
    for x in range(xdim):

        # If the point is located on sea, find its minimal 
        # distance to land
        # Else, the distance to land is 0
        if MASK[y,x] == 1:
            # Initialize a minimum distance to land variable
            min_dist = np.Inf

            # Iterate through all coast points
            for index, coastPoint in enumerate(COAST):
                    # Find the distance between the current coast point and the sea point
                    dist = haversine(( LAT[y,x], LON[y,x]), coastPoint)
                    
                    # If the distance found is less than the current
                    # minimal distance, update the minimal distance
                    if ( dist < min_dist ):
                        min_dist = dist
            
            # Add the minimal distance found to the 
            # distances to coastline matrix
            DIST[y,x] = min_dist



#---------------------Create new netcdf grid-------------------------------------

# Find absolute path in which cropped_grid.nc is to be stored
projPath = os.path.dirname(os.path.realpath(__file__))
cropGridPath = projPath + '/cropped_grid.nc'

# Create cropped_grid.nc netcdf file and associated crop_ds dataset
crop_ds = Dataset(cropGridPath, 'w', format = 'NETCDF4')

# Create (y, x) matrix to store LAT, LON, DIST, MASK and fLAT, fLON
y = crop_ds.createDimension('y', ydim)
x = crop_ds.createDimension('x', xdim)

# Create variables for netcdf data set
nav_lat = crop_ds.createVariable('nav_lat', 'f8', ('y', 'x')) # 'Point traceur' (lat)
nav_lon = crop_ds.createVariable('nav_lon', 'f8', ('y', 'x')) # 'Point traceur' (lon)
gphif   = crop_ds.createVariable('gphif', 'f8', ('y', 'x'))   # f (lat)
glamf   = crop_ds.createVariable('glamf', 'f8', ('y', 'x'))   # f (lon)
mask    = crop_ds.createVariable('mask', 'i1', ('y', 'x'))    # Land-ocean mask
dist    = crop_ds.createVariable('dist', 'f8', ('y', 'x'))    # Ocean-coastline distances

# Specify units
nav_lat.units = 'degrees'
nav_lon.units = 'degrees'
gphif.units = 'degrees'
glamf.units = 'degrees'
dist.units = 'meters'

# Attribute LAT and LON data to nav_lat and nav_lon
nav_lat[:, :] = LAT
nav_lon[:,:]  = LON
gphif[:, :]   = fLAT
glamf[:,:]    = fLON
mask[:, :]    = MASK
dist[:,:]     = DIST

# Close dataset
crop_ds.close()
