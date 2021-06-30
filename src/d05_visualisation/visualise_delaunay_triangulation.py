'''
Author: Beatrice Duval (bdu002)

--------------------------------------------------------------------------
Code that plots the Delaunay triangulation using the processed csv file.
--------------------------------------------------------------------------

'''

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import src.config
from src.d00_utils.load_csv import load_processed_csv, load_raw_csv

'''
_________________________________________________________________________________________
LOAD DATA SET
'''

# Initialize the config global variables (i.e. .csv file paths for all stages of data processing)
src.config.init()

# Retrieve a single data set (n is the index of an element in the list of .csv file paths)
n = 0
raw_csv_path = src.config.raw_csv_paths[n]
processed_csv_path = src.config.processed_csv_paths[n]

# Load a raw dataset
sLat, sLon, eLat, eLon = load_raw_csv( raw_csv_path )

# Load a processed data set
_, _, _, _, _, _, _, vertice_idx1, vertice_idx2, vertice_idx3 = load_processed_csv( processed_csv_path )


'''
_________________________________________________________________________________________
PLOT
'''

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

# Plot all triangles
for n in range(len(vertice_idx1)):
    ax.plot([sLon[vertice_idx1[n]], sLon[vertice_idx2[n]], sLon[vertice_idx3[n]], sLon[vertice_idx1[n]]],
            [sLat[vertice_idx1[n]], sLat[vertice_idx2[n]], sLat[vertice_idx3[n]], sLat[vertice_idx1[n]]], 
            color = 'xkcd:sky blue', 
            transform=ccrs.Geodetic())

# Add a title
fig.suptitle('Delaunay Triangulation', fontsize=14, y=1)

plt.show()

