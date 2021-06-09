'''
Author: Beatrice Duval (bdu002)

-------------------------------------------------------------------------------------
Code that plots individual triangles in the local grid CS and in lat, lon coordinates
with tracer point and x/y axis shown
-------------------------------------------------------------------------------------

'''

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from src.d01_data.load_processed_csv import sLat1, sLat2, sLat3, sLon1, sLon2, sLon3
from src.d01_data.load_converted_csv import sX1, sX2, sX3, sY1, sY2, sY3, j, i
from src.d01_data.load_grid import LAT, LON, fLAT, fLON
from src.d00_utils.grid_coord_system import define_gridCS

# Select which triangle to visualise
n = 100

# Initialize figure
fig = plt.figure()
    
# Initialize first subplot which will plot a triangle in lat,lon coordinates
ax_latlon = fig.add_subplot(211,
                        projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                            central_latitude=90))


# Plot the triangle
ax_latlon.plot([sLon1[n], sLon2[n], sLon3[n], sLon1[n]],
            [sLat1[n], sLat2[n], sLat3[n], sLat1[n]], 
            color = 'xkcd:sky blue', 
            transform=ccrs.Geodetic())

# Plot the tracer point
ax_latlon.scatter(LON[j[n],i[n]], LAT[j[n],i[n]], 
            color = 'red', 
            transform=ccrs.Geodetic(), label = "Tracer point")

# Get the triangle's B,C and D points which define the local grid CS
(BLat, BLon), (CLat, CLon), (DLat, DLon), _, _ = define_gridCS(fLAT, fLON, (j[n],i[n]))


 # Create a matplotlib transform from the cartopy coordinate system
crs = ccrs.Geodetic()
transform = crs._as_mpl_transform(ax_latlon)

# Draw the x and y axis of reference for the local grid CS using B, C and D
ax_latlon.plot([BLon, CLon, DLon], [BLat, CLat, DLat], 
            color = 'green', 
            transform=ccrs.Geodetic(), label = "Axis of reference in the local grid CS")

# Add labels to the drawn axis of reference
ax_latlon.annotate("y", (BLon, BLat), xytext=(1, 1), xycoords=transform, textcoords='offset points', color="green") #
ax_latlon.annotate("x", (DLon, DLat), xytext=(1, 1), xycoords=transform, textcoords='offset points', color="green")

# Add a legend 
ax_latlon.legend(loc='upper right', fontsize='x-small')

# Initialize second subplot which will plot a triangle in x,y coordinates in its local grid CS
ax_gridCS =  fig.add_subplot(212)

# Draw the triangle
ax_gridCS.plot([sX1[n], sX2[n], sX3[n], sX1[n]], [sY1[n], sY2[n], sY3[n], sY1[n]],
            color = 'xkcd:sky blue')

# Add labels to the x and y axis
ax_gridCS.set_ylabel('y')
ax_gridCS.set_xlabel('x')

# Add a title
fig.suptitle('Individual Data Point Triangle Drawn Using (lat, lon) Coordinates (top)\n and (x, y) Coordinates in a Local Grid Coordinate System (CS) (bottom) \n (triangle no. ' + str(n) + ')', fontsize=10, y=1)



plt.show()

