'''
Author: Beatrice Duval (bdu002)

----------------------------------------------
Testing code for tri_calc python script.
----------------------------------------------

Tests the find_nearestGrid function.
'''


import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from haversine import haversine
from pyproj import Proj
from scipy.spatial import Delaunay
from src.d00_utils.grid_coord_system import *
from src.d01_data.load01_raw_csv import sLat, sLon
from src.d01_data.load00_grid import LAT, LON, X_aeqd, Y_aeqd, fLAT, fLON


# Convert data points from lon,lat to x,y coordinates following 
# the Azimuthal Equidistant transform
p = Proj(proj='aeqd', ellps='WGS84', preserve_units=False)
xs, ys = p(sLon, sLat)


#-----------------------TEST 1 - find_nearestGridTracerPt-------------------------------

'''
We test the function that finds the nearest grid box (i.e. tracer point) 
for some (lat, lon) point. 

We produce scatter plots of distances between data points and the nearest grid 
point for a single csv file. Since grid tracer points are about 4-8 km away 
from each other, we expect the distance between data points and their nearest grid
point to be less than 4-8 km.
'''

# Initialize an array that will store the distances between data 
# points and the nearest grid point
distFromGrid = []

# Iterate through all data points
for idx in range(len(sLon)):
    # Find the nearest grid point and compute the distance
    # between that and the data point
    j, i = find_nearestGridTracerPt(X_aeqd,Y_aeqd, (xs[idx],ys[idx]))
    dist=haversine((sLat[idx], sLon[idx]),(LAT[j,i], LON[j, i]))
    
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
ax2.scatter(sLon, sLat, transform=ccrs.Geodetic())
    
plt.show()


#-----------------------TEST 2 - find_nearestGridTracerPt & find_aeqdTriCenter-------------------------------

'''
We simultaneously test the find_nearestGrid (again) and aeqdTriCenter 
functions by producing a plot of a triangle mesh, all of the triangles's
center points and all of the center points' nearest grid tracer point.
'''

# Generate a Delaunay triangulation
xy = list(zip(xs, ys))
tri = Delaunay(xy)

# Initialize a list for triangle center lat,lon points
center_lons = []
center_lats = []

# Initialize a list for center points' nearest tracer points
center_gridlons = []
center_gridlats = []


# Iterate through all mesh triangles
for n in range(len(tri.simplices)):
    # Find all x,y aeqd coordinates of the triangle's vertices
    x1_aeqd = xs[tri.simplices[n][0]]
    x2_aeqd = xs[tri.simplices[n][1]]
    x3_aeqd = xs[tri.simplices[n][2]]
    y1_aeqd = ys[tri.simplices[n][0]]
    y2_aeqd = ys[tri.simplices[n][1]]
    y3_aeqd = ys[tri.simplices[n][2]]

    # Find center point for the triangle in aeqd x,y coordinates
    center_x, center_y = find_aeqdTriCenter(x1_aeqd, x2_aeqd, x3_aeqd, y1_aeqd, y2_aeqd, y3_aeqd)

    # Convert center point back to lat,lon 
    # and add it to the list
    center_lon, center_lat = p(center_x, center_y, inverse = True)
    center_lons.append(center_lon)
    center_lats.append(center_lat)

    # Find the nearest grid tracer point relative to the center point
    # and append it to the list
    j, i = find_nearestGridTracerPt(X_aeqd,Y_aeqd, (center_x, center_y))
    center_gridlons.append( LON[j, i])
    center_gridlats.append( LAT[j, i])

#-----------
# Plot center points, the triangle mesh and the nearest grid tracer points 
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

# Plot the triangle mesh, center points and nearest grid tracer points
ax.triplot(sLon, sLat, tri.simplices, transform=ccrs.Geodetic())
ax.scatter(center_lons, center_lats, transform=ccrs.Geodetic(), color='green', marker = '*', label="Triangle centers")
ax.scatter(center_gridlons, center_gridlats, transform=ccrs.Geodetic(), color='red', marker = '.', label = "tracer points")

# Add a legend and show the plot
plt.legend(loc='upper right')
plt.show()



#-----------------------TEST 3 - define_gridCS & xy_gridCS-----------------------------

'''
We investigate the accuracy of the conversion from lat,lon to x,y 
coordinates in a local cartesian grid coordinate system (CS) using the
define_gridCS and xy_gridCS functions.

We choose some arbitrary grid tracer point (i,j) and draw its 
speed points (A(i,j), B(i-1,j), C(i-1, j-1) and D(i,j-1)) 
and some adjacent tracer points in the local cartesian grid 
CS defined by A-D. 

If the define_gridCS and xy_gridCS functions are accurate, C should be at 
(0,0) in the local coordinate system and we should observe the following 
positionning:

                (i-1, j+1)               (i,j+1)                  (i+1, j+1)


                            B(i-1, j)                 A(i,j)


                (i-1, j)                  (i,j)                    (i+1, j)


                            C(i-1, j-1)              D(i, j-1)


                (i-1, j-1)               (i, j-1)                 (i+1, j-1)  


'''

# Initialize lists for speed points x, y coordinates in the local grid CS
xABCDgridCS = []
yABCDgridCS = []

CS = define_gridCS(fLAT, fLON, (200, 200))

for i in [-1, 0]:
    for j in [-1, 0]:
        xgridCS, ygridCS = get_xy_gridCS( CS, (fLAT[200+j, 200+i], fLON[200+j, 200+i]) ) 
        xABCDgridCS.append(xgridCS)
        yABCDgridCS.append(ygridCS)

# Initialize lists for tracer points x, y coordinates in the local grid CS
xsgridCS = []
ysgridCS = []

for i in [-1, 0, 1]: 
    for j in [-1, 0, 1]:
        xgridCS, ygridCS = get_xy_gridCS( CS, (LAT[200+j, 200+i], LON[200+j, 200+i]) ) 
        xsgridCS.append(xgridCS)
        ysgridCS.append(ygridCS)

# Initialize a figure to plot the grid points in the local grid CS
fig = plt.figure()
ax1 = fig.add_subplot(111)

# Add a title
plt.suptitle("Points Plotted in the Cartesian Grid Coordinate\n System Relative to the (i,j) Grid Cell (or tracer point)",fontsize=12, y=1) 

# Plot all tracer points and speed points in the local grid CS
ax1.scatter(xsgridCS, ysgridCS, marker = '.', color = 'red', label = "tracer points")
ax1.scatter(xABCDgridCS, yABCDgridCS, marker = 'v', color = 'blue', label = "speed points")

# Annotate the grid points with their position relative to the 
# reference grid tracer point(i,j)
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