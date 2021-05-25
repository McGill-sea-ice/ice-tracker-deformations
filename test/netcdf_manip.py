'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------
Testing code for netcdf grid data visualization.
-----------------------------------------------

Produces a plot of grid points and colors them as a
    function of their distance from land.
'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os

################################### MAIN PROGRAM #################################

# Find absolute path in which cropped_grid.nc is stored
projPath = os.path.dirname(os.path.realpath(__file__))
gridPath = projPath + '/cropped_grid.nc'

# Load netcdf cropped NEMO dataset
grid_ds = Dataset(gridPath)

# Retrieve grid latitudes and longitudes
LAT = grid_ds['nav_lat'][:,:]
LON = grid_ds['nav_lon'][:,:]

# Retrieve distances from coast
DIST = grid_ds['dist'][:,:]

# Initialize figure
fig = plt.figure()
ax = fig.add_subplot(121,
                        projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                            central_latitude=90))
# Plot the grid points and color them as a function of
# their distance from land
ax.scatter(LON, LAT, c=DIST, transform=ccrs.Geodetic())

 # Add coastlines and gridlines
ax.coastlines()
ax.gridlines()

# Show plot
plt.show()
