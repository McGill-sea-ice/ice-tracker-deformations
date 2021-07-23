'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------------------------------
Visualise raw data files from config and identity them by the filename.
-----------------------------------------------------------------------

'''

import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import config
import utils_load_data as load_data

# Initialize figure and map
fig = plt.figure()
ax = fig.add_subplot(111, projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                             central_latitude=90))

# Set the map extent in order to see the 
# entire region of interest
ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

# Iterate though all datasets
for raw_path in config.data_paths['raw']:
    # Load the raw data
    try:
        raw_data = load_data.load_raw(raw_path)

    except load_data.DataFileError as dfe:
            print(str(dfe) + 'It will not be plotted.')
            continue

    # Plot the data points
    ax.scatter(list(raw_data['sLon']), list(raw_data['sLat']), transform=ccrs.Geodetic(), label=os.path.basename(raw_path))

# Add a legend and a title
ax.legend()
ax.set_title('Raw datasets Positionning')

# Add coastlines and gridlines
ax.coastlines()
ax.gridlines()

plt.show()
