'''
Author: Beatrice Duval (bdu002)

----------------------------------------------
Testing code for tri_calc python script.
----------------------------------------------

Tests the find_nearestGrid function.
'''

import os
from datetime import datetime
from math import asin, degrees, sqrt, tan

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from haversine import haversine
from netCDF4 import Dataset
from pyproj import Proj
from scipy.interpolate import griddata
from scipy.spatial import Delaunay

import os,sys,inspect
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir) 
from tri_calc import *


def show_tri(lats, lons, tri_simplices):
    ''' (array, array, Delaunay) -> None
    
    Plots lat, lon data points and a triangular mesh.

    Keyword arguments:
    lats -- array of data point latitudes
    lons -- array of data point longitudes
    tri  -- Delaunay object made using input lat,lon data points
    '''

    # Initialize figure
    fig = plt.figure()

    # Initialize subplot
    ax = fig.add_subplot(111,
                            projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                                central_latitude=90))

    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()

    # Plot the starting data points and the mesh triangles
    ax.scatter(lons, lats, transform=ccrs.Geodetic())
    ax.triplot(lons, lats, tri.simplices, transform=ccrs.Geodetic())

    plt.show()


################################### MAIN PROGRAM ###############################

#-------------------------Load grid-----------------------------------------

# Find absolute path in which cropped_grid.nc is stored
projPath = os.path.dirname(os.path.realpath(__file__))
gridPath = projPath + '/../cropped_grid.nc'

# Load netcdf cropped NEMO dataset
grid_ds = Dataset(gridPath)

# Retrieve grid 'points traceurs' and 'points vitesse' 
LAT  = grid_ds['nav_lat'][:,:]   # 'Points traceurs'(lon, lat)
LON  = grid_ds['nav_lon'][:,:]
Y    = grid_ds['nav_y'][:,:]     # 'Points traceurs'(y, x)
X    = grid_ds['nav_x'][:,:]
fLAT = grid_ds['gphif'][:,:]     # 'Points vitesse'(lon, lat)
fLON = grid_ds['glamf'][:,:]

#-------------------------Load data set-------------------------------------

# Set a csv file path and create a data frame
path = '/home/bdu002/2021_SeaIceDeformation/data/2020_MarApr_S1/pairs_20200303194831_20200304194135_1.csv'
df = pd.read_csv(path)

# Load lon, lat data point arrays
lons = df['sLon']
lats = df['sLat']

# Convert data points from lon,lat to x,y coordinates following 
# the Azimuthal Equidistant transform
p = Proj(proj='aeqd', ellps='WGS84', preserve_units=False)
xs, ys = p(lons, lats)


#-----------------------TEST 1 - find_nearestGrid-------------------------------

'''
We test the function that finds the nearest grid box (i.e. 'point traceur') 
for some (lat, lon) point. 

We produce scatter plots of distances between data points and the nearest grid 
point for a single csv file. Since grid 'points traceurs' are about 4-8 km away 
from each other, we expect the distance between data points and their nearest grid
point to be less than 4-8 km.
'''

# Initialize an array that will store the distances between data 
# points and the nearest grid point
distFromGrid = []

# Iterate through all data points
for idx in range(len(lons)):
    # Find the nearest grid point and compute the distance
    # between that and the data point
    j, i = nearestGrid(X,Y, (xs[idx],ys[idx]))
    dist=haversine((lats[idx], lons[idx]),(LAT[j,i], LON[j, i]))
    
    # Append that distance to the distances list
    distFromGrid.append(dist)

# Create a list of indices for plotting purposes
indices = list(range(0, len(distFromGrid)))

# Initialize a figure and add a title
fig = plt.figure()
plt.suptitle('Location of Data Points (right)\nand Distance Between Data Points and Their Nearest Grid Point (left)',fontsize=12, y=1)

#-----------
# Initialize a first subplot for a scatter plot of distances
ax1 = fig.add_subplot(121)

# Plot distances and corresponding indices
ax1.scatter(indices, distFromGrid, marker=".")

# Add x and y labels
ax1.set(xlabel='Data point index', ylabel='Distance (km)')

#-----------
# Initialize a second subplot to visualize where the data points are located 
ax2 = fig.add_subplot(122,
                        projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                            central_latitude=90))
    
# Set the map extent in order to see the entire region of interest
ax2.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

# Add coastlines and gridlines
ax2.coastlines()
ax2.gridlines()

# Plot the data points on the map
ax2.scatter(lons, lats, transform=ccrs.Geodetic())
    
plt.show()


#-----------------------TEST 2 - find_nearestGrid & aeqdTriCenter-------------------------------

'''
We simultaneously test the find_nearestGrid (again) and aeqdTriCenter 
functions by producing a plot of a triangle mesh, all of the triangles's
center points and all of the center points' nearest grid 'point traceur'.
'''

# Generate a Delaunay triangulation
xy = list(zip(xs, ys))
tri = Delaunay(xy)

# Initialize a list for triangle center lat,lon points
center_lons = []
center_lats = []

# Initialize a list for center points' nearest 'points traceurs'
center_gridlons = []
center_gridlats = []

# Iterate through all mesh triangles
for i in range(len(tri.simplices)):
    # Find center point for the triangle in aeqd x,y coordinates
    center_y, center_x = aeqdTriCenter(tri.simplices[i], ys, xs)

    # Convert center point back to lat,lon 
    # and add it to the list
    center_lon, center_lat = p(center_x, center_y, inverse = True)
    center_lons.append(center_lon)
    center_lats.append(center_lat)

    # Find the nearest grid 'point traceur' relative to the center point
    # and append it to the list
    j, i = nearestGrid(X,Y, (center_x, center_y))
    center_gridlons.append( LON[j, i])
    center_gridlats.append( LAT[j, i])

#-----------
# Plot center points, the triangle mesh and the nearest grid 'points traceurs' 
# relative to the center points
#-----------

# Initialize a figure and add a title
fig = plt.figure()
ax = fig.add_subplot(111,
                        projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                            central_latitude=90))
    
# Set the map extent in order to see the entire region of interest
ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

# Add coastlines and gridlines
ax.coastlines()
ax.gridlines()

# Plot the triangle mesh, center points and nearest grid 'points traceurs'
ax.triplot(lons, lats, tri.simplices, transform=ccrs.Geodetic())
ax.scatter(center_lons, center_lats, transform=ccrs.Geodetic(), color='green', marker = '*', label="Triangle centers")
ax.scatter(center_gridlons, center_gridlats, transform=ccrs.Geodetic(), color='red', marker = '.', label = "'Points traceurs'")

# Add a legend and show the plot
plt.legend(loc='upper right', fontsize='xx-small')
plt.show()



#-----------------------TEST 3 - gridCS & xy_gridCS-----------------------------

'''
We investigate the accuracy of the conversion from lat,lon to x,y 
coordinates in a local cartesian grid coordinate system (CS) using the
gridCS and xy_gridCS functions.

We choose some arbitrary grid 'point traceur' (i,j) and draw its 
'points de vitesse' (A(i,j), B(i-1,j), C(i-1, j-1) and D(i,j-1)) 
and some adjacent 'points traceurs'. 

If the gridCS and xy_gridCS functions are accurate, C should be at 
(0,0) in the local coordinate system and we should observe the following 
positionning:

                (i-1, j+1)               (i,j+1)                  (i+1, j+1)


                            B(i-1, j)                 A(i,j)


                (i-1, j)                  (i,j)                    (i+1, j)


                            C(i-1, j-1)              D(i, j-1)


                (i-1, j-1)               (i, j-1)                 (i+1, j-1)  


'''

# Initialize lists for 'points de vitesse' x, y coordinates in the local grid CS
xABCDgridCS = []
yABCDgridCS = []

CS = gridCS(fLAT, fLON, (200, 200))
for i in [-1, 0]:
    for j in [-1, 0]:
        xgridCS, ygridCS = xy_gridCS( CS, (fLAT[200+j, 200+i], fLON[200+j, 200+i]) ) 
        xABCDgridCS.append(xgridCS)
        yABCDgridCS.append(ygridCS)

# Initialize lists for 'points traceurs' x, y coordinates in the local grid CS
xsgridCS = []
ysgridCS = []

for i in [-1, 0, 1]: 
    for j in [-1, 0, 1]:
        xgridCS, ygridCS = xy_gridCS( CS, (LAT[200+j, 200+i], LON[200+j, 200+i]) ) 
        xsgridCS.append(xgridCS)
        ysgridCS.append(ygridCS)

# Initialize a figure to plot the grid points in the local grid CS
fig = plt.figure()
ax1 = fig.add_subplot(111)

# Add a title
plt.suptitle("Points Plotted in the Cartesian Grid Coordinate\n System Relative to the (i,j) Grid Cell (or 'Point traceur')",fontsize=12, y=1) 

# Plot all 'points traceurs' and 'points de vitesse' in the local grid CS
ax1.scatter(xsgridCS, ysgridCS, marker = '.', color = 'red', label = "'Points traceurs'")
ax1.scatter(xABCDgridCS, yABCDgridCS, marker = 'v', color = 'blue', label = "'Points de vitesse'")

# Annotate the grid points with their position relative to the 
# reference grid 'point traceur'(i,j)
ax1.annotate("C (i-1, j-1)", (xABCDgridCS[0], yABCDgridCS[0]))
ax1.annotate("B (i-1, j)", (xABCDgridCS[1], yABCDgridCS[1]))
ax1.annotate("D (i, j-1)", (xABCDgridCS[2], yABCDgridCS[2]))
ax1.annotate("A (i, j)", (xABCDgridCS[3], yABCDgridCS[3]))
ax1.annotate("(i-1, j-1)", (xsgridCS[0], ysgridCS[0]))
ax1.annotate("(i-1, j)", (xsgridCS[1], ysgridCS[1]))
ax1.annotate("(i-1, j+1)", (xsgridCS[2], ysgridCS[2]))
ax1.annotate("(i, j-1)", (xsgridCS[3], ysgridCS[3]))
ax1.annotate("(i, j)", (xsgridCS[4], ysgridCS[4]))
ax1.annotate("(i, j+1)", (xsgridCS[5], ysgridCS[5]))
ax1.annotate("(i+1, j-1)", (xsgridCS[6], ysgridCS[6]))
ax1.annotate("(i+1, j)", (xsgridCS[7], ysgridCS[7]))
ax1.annotate("(i+1, j+1)", (xsgridCS[8], ysgridCS[8]))

# Add a legend and show the plot
plt.legend(loc='best', fontsize='xx-small')
plt.show()