'''
Author: Beatrice Duval (bdu002)

Preliminary code. Displays S1 data points on a map.
'''

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from math import sin, asin, cos, sqrt, atan2, tan, radians, degrees

########################### FUNCTIONS ####################################

def dist(sLat, sLon, eLat, eLon):
    ''' (float, float, float, float) -> float 
    Finds the distance between two points.
    Input: starting and ending latitudes and longitudes (degrees)
    Output: returns the distance (km) between the starting and ending points
    >>> find_dist(50, 60, 50, 60)
    0.0
    >>> find_dist(68.5, 68.6, 68.7, 68.65)
    22.331316376877414
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



########################## MAIN PROGRAM #################################


#-----------------------Plot points on artic map-------------------------

# Create a data frame from a csv file
df = pd.read_csv('/home/bdu002/2021_SeaIceDeformation/data/2020_MarApr_S1/pairs_20200301015600_20200301200351_1.csv')

# Setup Arctic basemap
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')

# Add a coastline and fill continents
m.drawcoastlines()
m.fillcontinents(color='gray')

# Draw parallels so that there is 10 km between each
m.drawparallels(np.arange(-90.,90.,0.0899321605947705),labels=[0,0,0,0])
#m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1])

# Plot data points (start and end) on the basemap
sLon, sLat = m(list(df['sLon']), list(df['sLat']))
eLon, eLat = m(list(df['eLon']), list(df['eLat']))
m.scatter( sLon, sLat,marker=".", color="red", label='Starting positions')
m.scatter( eLon, eLat,marker=".", color="blue", label='Ending positions')

# Add a title, a subtitle and a legend, and show the plot
plt.suptitle('S1 Data Points on March 01 2020',fontsize=14, y=1)
plt.title('Distance between parallels: 10 km',fontsize=10)
plt.legend(fontsize=8)
plt.show()

#---------Scatter plot of distances between consecutive data points---------

# Initialize scatter plot arrays
distances = []
indices = []

# Iterate through rows of the csv file
for index, row in df.iterrows():
    if index != 0: # Indices count starts at 0
        # Compute the distance between the current point
        # and the previous one and add it to the distances
        # array
        distances.append(dist(row['sLat'], row['sLon'], prevLat, prevLon))
        indices.append(index)
    # Keep the current latitude and longitude for the next distance 
    # calculation
    prevLat = row['sLat']
    prevLon = row['sLon']

# Plot distances and indices
plt.scatter(indices, distances, marker=".")

# Add a x and y label, a title and show the plot
plt.xlabel("Data point index")
plt.ylabel("Distance (km)")
plt.title('Distance between consecutive starting data points - S1, March 01 2020')
plt.show()

# Test the function dist
print("Quebec-Montreal distance: " + str(dist(45.5016889, -73.567256, 46.829853, -71.254028)) + " km")