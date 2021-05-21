'''
Author: Beatrice Duval (bdu002)

Crops the NEMO grid to region of interest. 
Plots the resulting grid points. 
Creates a new netcdf cropped grid.

'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

################################### MAIN PROGRAM #################################

#-------------------------Trim NEMO grid-----------------------------------------

# Load netcdf NEMO dataset
gridPath = '/home/socn000/env_ubuntu-18.04-skylake-64/eccc-ppp4/datafiles/constants/oce/repository/master/CONCEPTS/creg012pe/grid/coordinates_CREG12_ext.nc'
ds = Dataset(gridPath)

# Create a cropped matrix for grid latitudes 
# and longitudes in the region of interest
LAT = ds['nav_lat'][1075:1750, 375:1200]
LON = ds['nav_lon'][1075:1750, 375:1200]

# Remove all points for which the latitude is
# less than 60 degrees
LAT[ LAT[:,:] < 60 ] = np.nan
LON[ LAT[:,:] == np.nan ] = np.nan

# TO-DO
# Remove all points for which the longitude is
# between -60 and -120 and latitude is less than
# 80 degrees (i.e. remove part of CAN archipelago)
#np.where( ( LAT < 80 ) & ( LON < -60 ) & ( LON > -120 ), np.nan, LAT)

# TO-DO
# Remove all points that are close to the coast

#---------------------Plot modified NEMO grid -------------------------------------

# Initialize figure
fig = plt.figure()
ax = fig.add_subplot(111,
                     projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                          central_latitude=90))
# Plot the cropped and filtered grid points
ax.scatter(LON, LAT, transform=ccrs.Geodetic())

# Add coastlines and gridlines
ax.coastlines()
ax.gridlines()

# Show plot
plt.show()


#---------------------Create new netcd4 grid------------------------------------

# Create cropped_grid.nc netcd4 file and associated crop_ds dataset
cropGridPath = '/home/bdu002/2021_SeaIceDeformation/cropped_grid.nc'
crop_ds = Dataset(cropGridPath, 'w', format = 'NETCDF4')

# Create (y, x) matrix to store cropped nav_lat and nav_lon
# whose size is the same as LAT and LON matrices
ydim, xdim = LAT.shape
y = crop_ds.createDimension('y', ydim)
x = crop_ds.createDimension('x', xdim)

# Create nav_lat and nav_lon variables for netcd4 data set
nav_lat = crop_ds.createVariable('nav_lat', 'f8', ('y', 'x'))
nav_lon = crop_ds.createVariable('nav_lon', 'f8', ('y', 'x'))

# Specify units for nav_lat and nav_lon
nav_lat.units = 'degrees'
nav_lon.units = 'degrees'

# Attribute LAT and LON data to nav_lat and nav_lon
nav_lat[:, :] = LAT
nav_lon[:,:] = LON

# Close dataset
crop_ds.close()
