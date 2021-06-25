'''
Author: Beatrice Duval (bdu002)

--------------------------------------------------------------------------
Code that plots the Delaunay triangulation using the processed csv file.
--------------------------------------------------------------------------

'''

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from src.d01_data.load02_processed_csv import sLat1, sLat2, sLat3, sLon1, sLon2, sLon3, eLat1, eLat2, eLat3, eLon1, eLon2, eLon3

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
for i in range(len(sLat1)):
    ax.plot([sLon1[i], sLon2[i], sLon3[i], sLon1[i]],
            [sLat1[i], sLat2[i], sLat3[i], sLat1[i]], 
            color = 'xkcd:sky blue', 
            transform=ccrs.Geodetic())

# Add a title
fig.suptitle('Delaunay Triangulation', fontsize=14, y=1)

plt.show()

