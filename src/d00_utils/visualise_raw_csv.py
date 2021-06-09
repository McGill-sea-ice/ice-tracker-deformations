'''
Author: Beatrice Duval (bdu002)

-------------------------------------------
Visualisation tools for raw csv data files. 
-------------------------------------------

'''

import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
from haversine import haversine

from .datetime_raw_csv import *

def show_SpatialCoverage( dir ):
    ''' (string) -> none

    Produces a plot of data points of all csv files under 
    the given directory on a single artic map
    
    Keyword arguments:
    dir -- path to directory containing csv files
    '''

    # Initialize figure and map
    fig = plt.figure()
    ax = fig.add_subplot(111,
                        projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                             central_latitude=90))
    
    # Set the map extent in order to see the 
    # entire region of interest
    ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

    # Iterate through files in dir
    for filename in os.listdir(dir):
        # Check if the file is a csv file
        if filename.endswith(".csv"):
            df = pd.read_csv(os.path.join(dir, filename))
            # Check if the file contains data
            if df.shape[0] > 1 :
                # Plot starting data points on the map
                ax.scatter(list(df['sLon']), list(df['sLat']), transform=ccrs.Geodetic())

    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()

    # Add a title
    fig.suptitle('Spatial coverage of data stored in ' + os.path.basename(os.path.normpath(dir)) + ' directory', fontsize=14, y=1)

    plt.show()


def show_ConsecutiveDists(path):
    ''' (string) -> none

    Produces a scatter plot of distances between consecutive 
    data points in csv file

    Keyword arguments:
    path -- csv file path 
    '''

    # Create a data frame from a csv file
    df = pd.read_csv(path)
    
    # Create datetime objects for the starting and ending
    # times of the data
    start, end = dataDatetimes(path)

    # Initialize scatter plot arrays
    distances = []
    indices = list(range(1, len(df.index)))

    # Iterate through rows of the csv file
    for index, row in df.iterrows():
        if index != 0: # Indices count starts at 0
            # Compute the distance between the current point
            # and the previous one and add it to the distances
            # array
            distances.append(haversine((row['sLat'], row['sLon']), (prevLat, prevLon)))

        # Keep the current latitude and longitude for the next distance 
        # calculation
        prevLat = row['sLat']
        prevLon = row['sLon']

    # Plot distances and indices
    plt.scatter(indices, distances, marker=".")

    # Add x and y labels, a title, subtitle and show the plot
    plt.xlabel("Data point index")
    plt.ylabel("Distance (km)")
    plt.suptitle('Distance between consecutive starting data points - S1',fontsize=14, y=1)
    plt.title(start.strftime("%b %d %Y %H:%M:%S"), fontsize=10)
    
    plt.show()

def show_numberedDP(path):
    ''' (string) -> none

    Plots numbered data points from a csv file on a map.

    Keyword arguments:
    path --  csv file path (absolute or relative) 
    '''

    # Create a data frame for a csv file
    df = pd.read_csv(path)
    
    # Create datetime objects for the starting and ending
    # times of the data
    start, _ = dataDatetimes(path)

    # Initialize figure and plot
    fig = plt.figure()
    ax = fig.add_subplot(111,
                        projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                            central_latitude=90))
    
    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()

    # Plot the starting data points
    ax.scatter(list(df['sLon']), list(df['sLat']), transform=ccrs.Geodetic())

    # Create a matplotlib transform from the cartopy coordinate system
    crs = ccrs.Geodetic()
    transform = crs._as_mpl_transform(ax)

    # Number each data point by its row number in the csv file
    for i, (slon, slat) in enumerate(zip(list(df['sLon']), list(df['sLat'])), start=0):
        ax.annotate(str(i), (slon, slat), xytext=(1, 1), 
                                           xycoords=transform,
                                           textcoords='offset points')
    
    # Add a title
    fig.suptitle('Numbered data points on ' + start.strftime("%b %d %Y %H:%M:%S") + ' - S1', fontsize=14, y=1)

    plt.show()


def show_DP(path):
    ''' (array, array) -> none

    Plots data points from a csv file on a map of the entire 
    arctic ocean.

    Keyword arguments:
    path --  csv file path (absolute or relative) 
    '''

    # Create a data frame
    df = pd.read_csv(path)

    # Retrieve data points
    sLon = df['sLon']       # Starting longitude
    sLat = df['sLat']       # Starting latitude
    
    # Create datetime objects for the starting and ending
    # times of the data
    start, end = dataDatetimes(path)

    # Initialize figure
    fig = plt.figure()
    
    # Initialize first subplot
    ax = fig.add_subplot(111,
                        projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                            central_latitude=90))
    
    # Set the map extent in order to see the entire region of interest
    ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()

    # Plot the starting data points
    #ax.scatter(list(df['sLosn']), list(df['sLat']), transform=ccrs.Geodetic())
    ax.scatter(sLon, sLat, transform=ccrs.Geodetic())

    # Add a title
    fig.suptitle('Data points on ' + start.strftime("%b %d %Y %H:%M:%S") + ' - S1', fontsize=14, y=1)

    plt.show()



