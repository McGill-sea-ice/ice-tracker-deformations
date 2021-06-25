'''
Author: Beatrice Duval (bdu002)

-------------------------------------------
Code that plots total sea-ice deformations.
-------------------------------------------

'''
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from src.d00_utils.datetime_raw_csv import dataDatetimes
from src.d01_data.get_data_paths import calculations_csv_path
from src.d01_data.load01_raw_csv import sLat, sLon
from src.d01_data.load02_processed_csv import (vertice_idx1, vertice_idx2,
                                               vertice_idx3)
from src.d01_data.load04_calculations_csv import eps_tot

# Get start and end times as datetime objects
start, end = dataDatetimes(calculations_csv_path)

# Stack the 3 vertices index arrays
triangles = np.stack((vertice_idx1, vertice_idx2, vertice_idx3), axis=-1)

# Initialize figure and add subplots
fig = plt.figure()

# Set the matplotlib projection and transform
proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
trans = ccrs.Geodetic()

# Initialize a subplot for a linear representation of total deformation
ax = fig.add_subplot(111, projection=proj)

# Plot the grid points and color them as a function of
# their total deformation rate
cb_epsTot_lin = ax.tripcolor( sLon, sLat, triangles, facecolors=eps_tot, transform=trans, cmap = 'plasma', vmin=0, vmax=0.1 )

# Hide deformations over land
ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

# Set the map extent in order to see the entire region of interest
ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

# Add a colorbar and include the displacements plot in order 
# to have it in the same size as the rest of the plots
plt.colorbar(cb_epsTot_lin, ax=ax)

# Add coastlines and gridlines
ax.coastlines()
ax.gridlines()

# Add a title and subtitle
plt.suptitle( 'Total Deformation Rate ($days^{-1}$)', fontsize=14, y=1)
plt.title(start.strftime("%b %d %Y %H:%M:%S") + " - " + end.strftime("%b %d %Y %H:%M:%S"), fontsize=10, y=1)

plt.show()
