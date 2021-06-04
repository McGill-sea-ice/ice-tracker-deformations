'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------
Testing code for netcdf grid data visualization.
-----------------------------------------------

Provides visualization tools for distance from land matrice,
    'points traceurs' and 'points vitesse' positionning, and
    grid cropping.
'''

import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from haversine import haversine
from netCDF4 import Dataset
from pyproj import Proj

################################ FUNCTIONS #################################

def show_distCoast(LAT, LON, DIST):
    ''' (float(y,x), float(y,x), float(y,x)) -> none
    
    Produces a plot of netcdf grid points and colors
    them as a functon of their distance from land

    Keyword arguments:
    LAT  -- 'Points traceurs' latitudes matrix (nxm)
    LON  -- 'Points traceurs' longitdes matrix (nxm)
    DIST --  Distances from land matrix (nxm)
    '''
    
    # Initialize figure
    fig = plt.figure()
    ax = fig.add_subplot(111,
                            projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                                central_latitude=90))

    # Set the map extent in order to see the 
    # entire region of interest
    ax.set_extent((-3000000, 4000000, 8500000, 11900000), ccrs.AzimuthalEquidistant())

    # Plot the grid points and color them as a function of
    # their distance from land
    cb = ax.scatter(LON, LAT, c=DIST, transform=ccrs.Geodetic())

    # Add a colorbar
    plt.colorbar(cb, orientation='vertical',ticklocation='auto')
    
    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()

    # Add a title
    fig.suptitle('Grid Points Colored as a Function of Distance (km) from Coast', fontsize=14, y=1)

    # Show plot
    plt.show()


def show_TraceurVitessePt(LAT, LON, fLAT, fLON, y, x):
    ''' (float(y,x), float(y,x), float(y,x), float(y,x), int, int) -> none
    
    Produces a plot of 'points traceurs' and 'points vitesse'
    that are annotated relatively to a (y,x) 'point traceur'
    
    Keyword arguments:
    LAT  -- 'Points traceurs' latitudes matrix (nxm)
    LON  -- 'Points traceurs' longitdes matrix (nxm)
    fLAT -- 'Points vitesse' latitudes matrix (nxm)
    fLON -- 'Points vitesse' longitudes matrix (nxm)
    y    -- Matrix single 'point traceur' index (column)
    x    -- Matrix single 'point traceur' index (row)
    '''

    # Initialize figure
    fig = plt.figure()
    ax = fig.add_subplot(111,
                            projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                                    central_latitude=90))


    # Create a matplotlib transform from the cartopy coordinate system
    crs = ccrs.Geodetic()
    transform = crs._as_mpl_transform(ax)

    # Plot 'points traceurs' and their neighboring 'points vitesse' (f points) 
    # and annotate them
    for i in [-1, 0, 1]: 
        for j in [-1, 0, 1]:
            if i==0 and j==0:
                # Create a label only once for 'points traceurs' and 'points vitesse'
                ax.scatter(LON[y+j,x+i], LAT[y+j,x+i], transform=ccrs.Geodetic(), color='red', label="'Points traceurs'")
                ax.scatter(fLON[y+j,x+i], fLAT[y+j,x+i], transform=ccrs.Geodetic(), color='blue', label="'Points vitesse'") 
            
            else:    
                ax.scatter(fLON[y+j,x+i], fLAT[y+j,x+i], transform=ccrs.Geodetic(), color='blue')
                ax.scatter(LON[y+j,x+i], LAT[y+j,x+i], transform=ccrs.Geodetic(), color='red')


            # Create text variable for annotation
            if i == 0:
                xtxt=""
            if i == 1:
                xtxt="+1"
            if i == -1:
                xtxt="-1"
            if j == 0:
                ytxt=""
            if j == 1:
                ytxt="+1"
            if j == -1:
                ytxt="-1" 
            
            # Indicate the relative (y,x) positionning 
            ax.annotate("(y" + ytxt + ", x" + xtxt + ")", (fLON[y+j,x+i], fLAT[y+j,x+i]), xytext=(1, 1), 
                                                                    xycoords=transform,
                                                                    textcoords='offset points')
            
            ax.annotate("(y" + ytxt + ", x" + xtxt + ")", (LON[y+j,x+i], LAT[y+j,x+i]), xytext=(1, 1), 
                                                                    xycoords=transform,
                                                                    textcoords='offset points')
            
    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()
    plt.legend(loc='lower left', fontsize='x-small')

    # Add a title
    plt.suptitle("Visualization of 'Points Traceurs' \nand Neighboring 'Points Vitesse'", fontsize=14, y=1)

    # Show plot
    plt.show()


def show_grid(LAT, LON, y1, y2, x1, x2):
    ''' (float(y,x), float(y,x), float(y,x), float(y,x), int, int) -> none
    
    Produces a plot of 'points traceurs' from a grid that is 
    cropped y-wise from y1 to y2 and x-wise from x1 to x2

    Keyword arguments:
    LAT  -- 'Points traceurs's latitudes matrix (yxx)
    LON  -- 'Points traceurs' longitdes matrix (yxx)
    y1   -- Starting column index
    y2   -- Ending column index (not included)
    x1   -- Starting row index
    x2   -- Ending row index (not included)
    '''

    # Crop LAT/LON matrices 
    LAT = LAT[y1:y2, x1:x2] 
    LON = LON[y1:y2, x1:x2]

    # Initialize figure
    fig = plt.figure()
    ax = fig.add_subplot(111,
                            projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                                    central_latitude=90))
    # Set the map extent in order to see the 
    # entire region of interest
    ax.set_extent((-3000000, 4000000, 8500000, 11900000), ccrs.AzimuthalEquidistant())

    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()

    # Plot LAT/LON grid points and data points
    ax.scatter(LON, LAT, transform=ccrs.Geodetic(), color="red", marker='.')

    plt.show()



################################### MAIN PROGRAM ###############################

if __name__ == '__main__':

    #-----------------------Retrieve cropped RIOPS data-------------------------------

    # Find absolute path in which cropped_grid.nc is stored
    projPath = os.path.dirname(os.path.realpath(__file__))
    gridPath = projPath + '/../cropped_grid.nc'

    # Load netcdf cropped NEMO dataset
    cropGrid_ds = Dataset(gridPath)

    # Retrieve grid latitudes and longitudes
    cropLAT = cropGrid_ds['nav_lat'][:,:]
    cropLON = cropGrid_ds['nav_lon'][:,:]

    # Retrieve distances from coast
    cropDIST = cropGrid_ds['dist'][:,:]

    #-----------------------Retrieve original RIOPS data-------------------------------

    # Load netcdf NEMO dataset
    grid_ds = Dataset('/home/socn000/env_ubuntu-18.04-skylake-64/eccc-ppp4/datafiles/constants/oce/repository/master/CONCEPTS/creg012pe/grid/coordinates_CREG12_ext.nc')

    LAT = grid_ds['nav_lat'][:,:] # 'Point traceur'
    LON = grid_ds['nav_lon'][:,:]

    # Retrieve f grid latitudes and longitudes
    fLAT = grid_ds['gphif'][:,:]
    fLON = grid_ds['glamf'][:,:]






