'''
Author: Beatrice Duval (bdu002)

----------------------------------------------
Testing code for csv data file visualization. 
----------------------------------------------

Provides tools that:
    Display S1 data points on a map (numbered) for individual 
        csv files.  
    Plot the distances between consecutive data points for 
        individual csv files. 
    Display S1 data points on a map from all csv files under 
        a directory.
'''

import os
from datetime import datetime
from math import asin, degrees, sqrt, tan

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Proj
from haversine import haversine
from scipy.spatial import Delaunay
from scipy.interpolate import griddata

################################ FUNCTIONS #################################


def find_gridLatDeg(dist):
    ''' (float) -> float 
    
    Takes as input desired spacing (km) between parallels to be mapped
    and returns the angle (degrees) between parallels to be displayed on 
    a map

    Keyword arguments:
    dist -- Spacing between parallels (km)
    
    >>> find_gridDeg(0)
    0.0
    >>> find_gridDeg(50)
    0.44966080296032374
    '''
    # Set the approx. value of Earth radius (km)
    R = 6371

    # Use the 'haversine' formula backwards to find dLat
    # for a given input distance
    c = dist / R
    a = (1 + tan(c / 2) ** 2) ** (-1)

    # Set dLon to zero, and dLat to be the value
    # we want to compute
    dLat = 2 * asin(sqrt(a))

    return 180 - degrees(dLat)
    

def dataDate(p):
    ''' (str) -> datetime, datetime
    
    Takes as input RCM or S1 csv data file path
    and returns Start and end times as datetime objects
    
    Keyword arguments:
    p -- csv data file path (absolute or relative)

    >>> start, end = dataDate('2020_MarApr_S1/pairs_20200301033644_20200303032021_1.csv')
    >>> print(start.strftime("%b %d %Y %H:%M:%S"))
    Mar 01 2020 03:36:44
    >>> print(end.strftime("%b %d %Y %H:%M:%S"))
    Mar 03 2020 03:20:21
    '''

    # Create datetime objects for data starting and ending times
    # using the times specified in the csv file names
    start = datetime(int(p[-35:-31]), # Year
                     int(p[-31:-29]), # Month
                     int(p[-29:-27]), # Day
                     int(p[-27:-25]), # Hour
                     int(p[-25:-23]), # Minute
                     int(p[-23:-21])) # Second

    end   = datetime(int(p[-20:-16]), # Year
                     int(p[-16:-14]), # Month
                     int(p[-14:-12]), # Day
                     int(p[-12:-10]), # Hour
                     int(p[-10:-8]),  # Minute
                     int(p[-8:-6]))   # Second

    return start, end

def dT(s, e):
    ''' (datetime, datetime) -> float
    
    Takes as input start and end times as datetime objects
    and returns delta time (s)

    Keyword arguments:
    s -- start time (datetime)
    e -- end time (datetime)

    >>> start, end = dataDate('2020_MarApr_S1/pairs_20200301033644_20200303032021_1.csv')
    >>> print(dT(start, end))
    171817.0
    '''
    return (e-s).total_seconds()

def show_SpatialCoverage( dir ):
    ''' (string) -> none

    Produces a plot of data points of all csv files under 
    the given directory on a single artic map
    
    Keyword arguments:
    dir -- path to directory containing csv files (absolute or relative) 
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
    path -- csv file path (absolute or relative) 
    '''

    # Create a data frame from a csv file
    df = pd.read_csv(path)
    
    # Create datetime objects for the starting and ending
    # times of the data
    start, end = dataDate(path)

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
    start, end = dataDate(path)

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
    
    fig.suptitle('Numbered data points on ' + start.strftime("%b %d %Y %H:%M:%S") + ' - S1', fontsize=14, y=1)

    
    plt.show()


def show_DP(path):
    ''' (string) -> none

    Plots data points from a csv file on a map of the entire 
    arctic ocean.

    Keyword arguments:
    path --  csv file path (absolute or relative) 
    '''

    # Retrieve data from a csv file
    lat, lon, _ , _ = get_trackerDP(path)

    # Create datetime objects for the starting and ending
    # times of the data
    start, end = dataDate(path)

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
    ax.scatter(lon, lat, transform=ccrs.Geodetic())

    # Add a title
    fig.suptitle('Data points on ' + start.strftime("%b %d %Y %H:%M:%S") + ' - S1', fontsize=14, y=1)

    plt.show()


def get_trackerDP(path):
    ''' (string) -> array, array, array, array

    Retrieves data from csv file.

    Keyword arguments:
    path --  csv file path (absolute or relative)
    '''

    # Create a data frame
    df = pd.read_csv(path)

    # Retrieve data points
    slon = df['sLon']       # Starting longitude
    slat = df['sLat']       # Starting latitude
    elon = df['eLon']       # Ending longitude
    elat = df['eLat']       # Ending latitude

    return slat, slon, elat, elon


################################### MAIN PROGRAM #################################

if __name__ == '__main__':
    print("Main")


