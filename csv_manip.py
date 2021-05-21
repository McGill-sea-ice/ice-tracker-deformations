'''
Author: Beatrice Duval (bdu002)

Preliminary code. 

Displays S1 data points on a map (numbered) from individual 
    csv files.  
Plots the distances between consecutive data points. 
Displays S1 data points on a map from all csv files in 
    a subdirectory.
'''

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from math import sin, asin, cos, sqrt, atan2, tan, radians, degrees
from datetime import datetime
import os
from haversine import haversine

################################ FUNCTIONS #################################

def dist(sLat, sLon, eLat, eLon):
    ''' (float, float, float, float) -> float 
    Finds the distance between two points.
    Input: starting and ending latitudes and longitudes (degrees)
    Output: returns the distance (km) between the starting and ending points
    >>> print( dist(50, 60, 50, 60) )
    0.0
    >>> print( dist(45.5016889, -73.567256, 46.829853, -71.254028) ) # Montreal-QC distance
    231.38209668942255
    '''
    # Set the approx. value of Earth radius (km)
    R = 6371

    # Convert latitudes and longitudes to radians
    sLat = radians(sLat)
    sLon = radians(sLon)
    eLat = radians(eLat)
    eLon = radians(eLon)

    # Compute variation in latitude and longitude
    dLat = eLat - sLat
    dLon = eLon - sLon

    # Compute the distance between both points using the 'haversine' formula
    a = sin(dLat / 2)**2 + cos(sLat) * cos(eLat) * sin(dLon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    dist = R * c

    return dist


def find_gridLatDeg(dist):
    ''' (float) -> float 
    Input: desired spacing (km) between parallels to be mapped
    Output: returns the angle (degrees) between parallels to 
            be displayed on a map
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
    Input: RCM or S1 csv data file path (absolute or relative)
    Output: Start and end times as datetime objects
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
    Input: Start and end times as datetime objects
    Output: Delta time (s)
    >>> start, end = dataDate('2020_MarApr_S1/pairs_20200301033644_20200303032021_1.csv')
    >>> print(dT(start, end))
    171817.0
    '''
    return (e-s).total_seconds()

################################### MAIN PROGRAM #################################

# List of S1 data paths

paths = ['/home/bdu002/2021_SeaIceDeformation/data/2020_MarApr_S1/pairs_20200301033644_20200303032021_1.csv', 
         '/home/bdu002/2021_SeaIceDeformation/data/2020_MarApr_RCM/pairs_20200302024700_20200305023923_1.csv', 
         '/home/bdu002/2021_SeaIceDeformation/data/2020_MarApr_RCM/pairs_20200305023923_20200305181645_1.csv'
        ]

# Iterate through data paths
for path in paths:

    # Create a data frame from a csv file
    df = pd.read_csv(path)
    
    # Create datetime objects for the starting and ending
    # times of the data
    start, end = dataDate(path)


    #-----------------------Plot data points on artic map-------------------------

    fig, ax = plt.subplots()

    # Setup Arctic basemap
    m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l', ax=ax)

    # Add a coastline and fill continents
    m.drawcoastlines()
    m.fillcontinents(color='gray')

    # Draw parallels and meridians
    #m.drawparallels(np.arange(-90.,90.,0.0899321605947705),labels=[0,0,0,0])  # 10 km between each parallel
    m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1])

    # Plot data points (start and end) on the basemap
    sLon, sLat = m(list(df['sLon']), list(df['sLat']))
    #eLon, eLat = m(list(df['eLon']), list(df['eLat'])) # end points
    ax.scatter( sLon, sLat, marker=".", color="red", label='Starting positions')
    #m.scatter( eLon, eLat,marker=".", color="blue", label='Ending positions') # end points

    # Add a title, a subtitle and a legend, and show the plot
    plt.suptitle('S1 Data Points',fontsize=14, y=1)
    plt.title(start.strftime("%b %d %Y %H:%M:%S"), fontsize=10)
    plt.legend(fontsize=8)

    plt.show()

    #------------------Plot numbered data points on artic map-------------------------
    
    fig, ax = plt.subplots()

    # Setup Arctic basemap
    m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l', ax=ax)

    # Add a coastline and fill continents
    m.drawcoastlines()
    m.fillcontinents(color='gray')

    # Draw parallels and meridians
    m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1])

    # Plot data points (start) on the basemap
    sLon, sLat = m(list(df['sLon']), list(df['sLat']))
    ax.scatter( sLon, sLat,marker=".", color="red", label='Starting positions')
    
    # Number each data point by its row number in the csv file
    for i, (slon, slat) in enumerate(zip(sLon, sLat), start=0):
        ax.annotate(str(i), (slon, slat), xytext=(1, 1), textcoords='offset points')

    # Add a title, a subtitle and a legend, and show the plot
    plt.suptitle('S1 Data Points - Numbered',fontsize=14, y=1)
    plt.title(start.strftime("%b %d %Y %H:%M:%S") + ' (Shown is the upper corner of data)', fontsize=10)
    plt.legend(fontsize=8)
    
    plt.show()

    #---------Scatter plot of distances between consecutive data points---------

    # Initialize scatter plot arrays
    distances = []
    indices = list(range(1, len(df.index)))

    # Iterate through rows of the csv file
    for index, row in df.iterrows():
        if index != 0: # Indices count starts at 0
            # Compute the distance between the current point
            # and the previous one and add it to the distances
            # array
            distances.append(dist(row['sLat'], row['sLon'], prevLat, prevLon))
            #indices.append(index)
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

   
#---------------Plot data points of all csv files on single artic map----------------------

# Setup Arctic basemap
fig, ax = plt.subplots()
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l', ax=ax)

# Add a coastline and fill continents
m.drawcoastlines()

# Draw parallels and meridians
m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1])

# Add a title, a subtitle and a legend, and show the plot
plt.suptitle('S1 Data Points-All files combined',fontsize=14, y=1)

directory = '/home/bdu002/2021_SeaIceDeformation/data/2020_MarApr_S1'

for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        df = pd.read_csv(os.path.join(directory, filename))
        if df.shape[0] > 1 :
            # Plot data points (start and end) on the basemap
            sLon, sLat = m(list(df['sLon']), list(df['sLat']))
            ax.scatter( sLon, sLat, marker=".", color="red", label='Starting positions')
    else:
        continue


plt.show()

