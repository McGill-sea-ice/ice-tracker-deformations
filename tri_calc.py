'''
Author: Beatrice Duval (bdu002)

---------------------------------------
Code that computes sea ice deformation.
---------------------------------------

Gathers the data points into triangular data cells.
Associates each triangle center point to a grid cell using the Azimuthal 
    Equidistant (aeqd) transform.
Uses these grid cells as local cartesian coordinate systems to calculate 
    sea ice deformations.
'''

import os
from math import sqrt

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from haversine import haversine
from netCDF4 import Dataset
from pyproj import Proj
from scipy.spatial import Delaunay

################################ FUNCTIONS #################################

def nearestGrid(gridX, gridY, xyPt):
    '''  (float(j,i), float(j,i), float, float) -> (int, int)

    Finds the nearest grid box for an input point. 
    All points have (x,y) coordinates that follow the Azimuthal 
    Equidistant (aeqd) transform.
    Returns the (j,i) index of the nearest grid box 'point traceur'.
    
    Keyword arguments:
    gridX -- x coordinates for a grid matrix of 'points traceurs' (jxi)
    gridY -- y coordinates for a grid matrix of 'points traceurs' (jxi)
    xPt   -- x coordinate for the input point
    yPt   -- y coordinate for the input point
    
    source: https://kbkb-wx-python.blogspot.com/2016/08/find-nearest-latitude-and-longitude.html
    '''
    
    xPt, yPt = xyPt

    # Compute x and y differences between all x,y grid 
    # points and the input x,y data point
    absDeltaX = np.abs(gridX - xPt)
    absDeltaY = np.abs(gridY - yPt)

    # Find local grid delta maximums
    maxDeltas = np.maximum(absDeltaX, absDeltaY)

    # The minimum local maximum is at the nearest grid location
    nearest_idx = np.argmin(maxDeltas)

    # Retrieve the nearest grid point's (j,i) index and return it
    nearest_j, nearest_i = np.unravel_index(nearest_idx, gridX.shape)

    return nearest_j, nearest_i

def aeqdTriCenter(tri_simplice, xs, ys):
    ''' (array, array, array) ->  tuple

    Finds the center of a triangle in Azimuthal Equidistant 
    (aeqd) transform coordinates. 
    Returns the aqed (x, y) coordinates.

    Keyword arguments:
    tri_simplice -- Array of indices (size: 3) of the data points forming the vertices of a triangle
    xs -- Array of x coordinates of data points in the aeqd transform 
    ys -- Array of y coordinates of data points in the aeqd transform 
    '''
    # Retrieve the x, y coords of all triangle vertices
    y1 = ys[tri_simplice[0]]
    x1 = xs[tri_simplice[0]]

    y2 = ys[tri_simplice[1]]
    x2 = xs[tri_simplice[1]]

    y3 = ys[tri_simplice[2]]
    x3 = xs[tri_simplice[2]]

    # Average the x and the y coords (respectively) to find the center
    center_y = (y1 + y2 + y3) / 3
    center_x = (x1 + x2 + x3) / 3

    return center_x, center_y


def gridCS(fLAT, fLON, ji):
    '''  (float(j,i), float(j,i), tuple) -> (tuple, tuple, tuple, tuple, float, float)
    
    Finds the data needed to define and use the local 
    grid coordinate system (CS), i.e.:

            B(0,b)     A

            C(0,0)   D(d,0)
    
    Returns the (lat, lon) coords of B/C/D and the value of b and d.

    Keyword arguments:
    fLAT -- latitudes for a grid matrix of 'points vitesse' (jxi)
    fLON -- longitudes for a grid matrix of 'points vitesse' (jxi)
    ji   -- (j, i) tuple of the matrix index at which the reference 
            grid 'point traceur' is located
    
    '''
    # Unpack the input 'point traceur' indices
    j, i = ji
    
    # Find the lat/lon of the 'point traceur''s needed 'points vitesses'
    B = (fLAT[j, i-1], fLON[j, i-1])
    C = (fLAT[j-1, i-1], fLON[j-1, i-1])
    D = (fLAT[j-1, i], fLON[j-1, i])  

    # Find the distance between the origin (C) and B/D (respectively)
    b =  haversine(B, C)
    d =  haversine(D, C)

    return B, C, D, b, d

def xy_gridCS(gridCS, latlonPt ):
    ''' (tuple, tuple) -> (tuple)
    
    Returns the x,y coordinates of a lat,lon point in the local 
    grid CS defined by the return value of the gridCS function.
    
    Keyword arguments:
    gridCS   -- return value of gridCS function
    latlonPt -- (lat, lon) tuple  
    '''

    B, C, D, b, d = gridCS

    a1 = haversine(B, latlonPt)
    a2 = haversine(C, latlonPt)
    a3 = haversine(D, latlonPt)
    
    x = ( a2**2 - a3**2 + d**2 ) / (2 * d)
    y = ( a2**2 - a1**2 + b**2 ) / (2 * b)

    return x, y

def dudx():
    '''
    '''

def dvdy():
    '''
    '''

def dudy():
    '''
    '''

def dvdx():
    '''
    '''

def A():
    '''
    '''



def eps_I():
    '''
    '''

    return dudx + dvdy

def eps_II():
    '''
    '''

    return sqrt( (dudx - dvdy )**2 + ( dudy + dvdx )**2 )

def eps_tot(eps_I, eps_II):
    '''
    '''

    return sqrt( eps_I**2 + eps_II**2 )

################################### MAIN PROGRAM ###########################


#-------------------------Load grid-----------------------------------------

# Find absolute path in which cropped_grid.nc is stored
projPath = os.path.dirname(os.path.realpath(__file__))
gridPath = projPath + '/cropped_grid.nc'

# Load netcdf cropped NEMO dataset
grid_ds = Dataset(gridPath)

# Retrieve grid 'points traceurs' grid 
LAT = grid_ds['nav_lat'][:,:] # (lon, lat)
LON = grid_ds['nav_lon'][:,:]
Y = grid_ds['nav_y'][:,:]     #  (y, x)
X = grid_ds['nav_x'][:,:]


#-------------------------Load data set-------------------------------------

# Set a csv file path and create a data frame
path = '/home/bdu002/2021_SeaIceDeformation/data/2020_MarApr_RCM/pairs_20200428190223_20200501192643_1.csv'
df = pd.read_csv(path)

# Load lon, lat data point arrays
lon = df['sLon']
lat = df['sLat']

# Convert data points from lon,lat to x,y coordinates following 
# the Azimuthal Equidistant transform
p = Proj(proj='aeqd', ellps='WGS84', preserve_units=False)
x, y = p(lon, lat)



#-------------------------Delaunay triangulation------------------------------

# Generate a Delaunay triangulation
xy = list(zip(x, y))
tri = Delaunay(xy)


'''
#-------------------Crop grid to cover only data set---------------------------

# Find max and min latitudes and longitudes in the data set
max_lon = lon.max()
min_lon = lon.min()
max_lat = lat.max()
min_lat = lat.min()

# Retrieve grid indices that are within that range of latitudes and longitudes
j, i = np.where((min_lat <= LAT) & (LAT <= max_lat) &  (min_lon <= LON) & (LON <= max_lon))

# Find max j,i and min j,i and slice the grid matrices
max_i =  
max_j =
min_i = 
min_j = 
'''

