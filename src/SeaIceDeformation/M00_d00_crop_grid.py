'''
Author: Beatrice Duval (bdu002)

--------------------------------------------------------------------
Crop grid
--------------------------------------------------------------------

Module for creation of a grid (netcdf format), adapted from original RIOPS grid.

    - Crops the RIOPS grid to region of interest
    - Creates a new netcdf cropped grid containing tracer points (lat,lon and 
      x,y coordinates in the aeqd transform), speed points (lat,lon), distance 
      between tracer points and land, and a sea-land mask
    - Prints the % of the work done when creating the distance to coastline matrix
'''

from math import nan

import numpy as np
from haversine import haversine
from netCDF4 import Dataset
from pyproj import Proj

import config
from utils_grid_coord_system import find_nearestGridTracerPt

################################### MAIN PROGRAM #################################

# Set indices at which the grid is to be cropped, i.e. set the range of interest
y1 = 1075
y2 = 1750
x1 = 375
x2 = 1200


#-------------------------Load RIOPS grid-----------------------------------------

# Load netcdf RIOPS grid and create a dataset
gridPath = config.config['Grid']['riops_grid']
grid_ds = Dataset(gridPath)

# Create a matrix for grid latitudes and longitudes
LAT = grid_ds['nav_lat'][:, :] # Tracer point t
LON = grid_ds['nav_lon'][:, :]

fLAT =grid_ds['gphif'][:, :]   # Speed point f
fLON = grid_ds['glamf'][:, :]


#-------------Create a tracer point matrix in x,y coordinates-----------------

# Convert grid points from lat/lon to x,y coordinates following 
# the Azimuthal Equidistant projection
p = Proj(proj='aeqd', ellps='WGS84', preserve_units=False)
X, Y = p(LON, LAT)

#---------------------Create land-sea mask -------------------------------------

# Load bathymetry dataset
bathyPath = config.config['Grid']['bathymetry']
bathy_ds = Dataset(bathyPath)

# Create a matrix from bathymetry
MASK = bathy_ds['Bathymetry'][:, :]

# Set all land values (>0) to 0 and ocean values (=0)
# to 1 
MASK[ MASK[:,:] > 0 ] = 1
MASK[ MASK[:,:] == 0 ] = 0

#-----------------Create distance to coastline matrix --------------------------
    
# Initialize distances matrix to a zero matrix
DIST = np.zeros( LAT.shape )

# Create a masked tracer point matrix in x,y coordinates
# where all ocean points are set to nan
maskedX = X.copy()
maskedY = Y.copy()

maskedX[ MASK[:,:] == 1 ] = nan
maskedY[ MASK[:,:] == 1 ] = nan

# Crop the masked tracer point matrix (to reduce run time) but make it cover a 
# large enough area so that the nearest land point to every point of interest can be found
maskedX = maskedX[1000:1795, 300:1300]
maskedY = maskedY[1000:1795, 300:1300]

# Find number of data points to process and set number of processed data points to 0
# This will be used to print the % of the job done
num = (y2-y1) * (x2-x1)
count = 0

print('--- Populating a distance to coastline matrix ---')

# Iterate through the grid points of interest (from y1 to y2, and x1 to x2)
for y in range(y1, y2):
    for x in range(x1, x2):
        
        # Increase number of processed data points by 1 and print the % of work done if we are at a multiple of 10,000
        count += 1
        if count % 10000 == 0:
            print( str( count/num * 100 ) + '% done.' )

        # If the point is located on sea, find its minimal 
        # distance to land
        # Else, the distance to land is 0
        if MASK[y,x] == 1:
            j, i = find_nearestGridTracerPt(maskedX, maskedY, (X[y,x], Y[y,x]))

            # Adjust the cropped indices we found to the original (uncropped) indices
            j, i = j+1000, i+300 

            # Add the minimal distance found to the distances to coastline matrix
            DIST[y,x] = haversine( (LAT[y,x], LON[y,x]), (LAT[j,i], LON[j,i]) )


print('--- Distance to coastline matrix is complete ---')


#-------------------------Trim RIOPS grid-----------------------------------------

# Crop all grid matrices using the range of interest defined at the beginning of the module
LAT  = LAT[y1:y2, x1:x2]
LON  = LON[y1:y2, x1:x2]
Y    = Y[y1:y2, x1:x2]
X    = X[y1:y2, x1:x2]
fLAT = fLAT[y1:y2, x1:x2]
fLON = fLON[y1:y2, x1:x2]
MASK = MASK[y1:y2, x1:x2]
DIST = DIST[y1:y2, x1:x2]

# Retrieve new size for grid matrices
ydim, xdim = LAT.shape

#---------------------Create new netcdf grid-------------------------------------

# Find absolute path in which cropped_grid.nc is to be stored
cropGridPath = config.config['Grid']['cropped_grid']

# Create cropped_grid.nc netcdf file and associated crop_ds dataset
crop_ds = Dataset(cropGridPath, 'w', format = 'NETCDF4')

# Create (y, x) matrix to store LAT, LON, Y, X, DIST, MASK and fLAT, fLON
y = crop_ds.createDimension('y', ydim)
x = crop_ds.createDimension('x', xdim)

# Create variables for netcdf data set
nav_lat = crop_ds.createVariable('nav_lat', 'f8', ('y', 'x')) # 'Point traceur' (lat)
nav_lon = crop_ds.createVariable('nav_lon', 'f8', ('y', 'x')) # 'Point traceur' (lon)
nav_y   = crop_ds.createVariable('nav_y', 'f8', ('y', 'x'))   # 'Point traceur' (y in aeqd transf.)
nav_x   = crop_ds.createVariable('nav_x', 'f8', ('y', 'x'))   # 'Point traceur' (x in aeqd transf.)
gphif   = crop_ds.createVariable('gphif', 'f8', ('y', 'x'))   # 'Point vitesse' (lat)
glamf   = crop_ds.createVariable('glamf', 'f8', ('y', 'x'))   # 'Point vitesse' (lon)
mask    = crop_ds.createVariable('mask', 'i1', ('y', 'x'))    # Land-ocean mask
dist    = crop_ds.createVariable('dist', 'f8', ('y', 'x'))    # Ocean-coastline distances

# Specify units
nav_lat.units = 'degrees'
nav_lon.units = 'degrees'
gphif.units   = 'degrees'
glamf.units   = 'degrees'
dist.units    = 'kilometers'

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
