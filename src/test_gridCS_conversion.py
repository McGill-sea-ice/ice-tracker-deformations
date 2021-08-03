'''
Author: Beatrice Duval (bdu002)
-------------------------------------------------------------------------------------
Code that plots individual triangles in the local grid CS and in lat, lon coordinates
with tracer point and x/y axis shown
-------------------------------------------------------------------------------------
'''
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import src.config
from utils_grid_coord_system import define_gridCS
import utils_load_data as load_data
import utils_load_grid as load_grid

'''
_________________________________________________________________________________________
LOAD DATA SET AND GRID
'''
# Initialize the config global variables (i.e. .csv file path lists for all stages of data processing)
src.config.init()
# Retrieve a single data set (n is the index of an element in the list of .csv file paths)
n = 0
raw_csv_path = src.config.raw_csv_paths[n]
processed_csv_path = src.config.processed_csv_paths[n]
converted_csv_path = src.config.converted_csv_paths[n]
# Load a raw dataset
sLat, sLon, eLat, eLon = load_raw_csv( raw_csv_path )
# Load a processed dataset
_, _, _, _, _, _, _, vertice_idx1, vertice_idx2, vertice_idx3 = load_processed_csv( processed_csv_path )
# Load a converted dataset
_, sX1,sX2, sX3, sY1, sY2, sY3, _, _, _, _, _, _, j_tracers, i_tracers = load_converted_csv( converted_csv_path )


# Load the RIOPS grid
grid = load_grid.load_grid()
fLAT = grid['fLAT']
fLON = grid['fLON']
LAT = grid['LAT']
LON = grid['LON']
X_aeqd = grid['X_aeqd']
Y_aeqd = grid['Y_aeqd']

'''
_________________________________________________________________________________________
PLOT
'''
# Select which triangle to visualise
n = 100
# Initialize figure
fig = plt.figure()
    
# Initialize first subplot which will plot a triangle in lat,lon coordinates
ax_latlon = fig.add_subplot(211,
                        projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                            central_latitude=90))
# Plot the triangle
ax_latlon.plot([sLon[vertice_idx1[n]], sLon[vertice_idx2[n]], sLon[vertice_idx3[n]], sLon[vertice_idx1[n]]],
            [sLat[vertice_idx1[n]], sLat[vertice_idx2[n]], sLat[vertice_idx3[n]], sLat[vertice_idx1[n]]], 
            color = 'xkcd:sky blue', 
            transform=ccrs.Geodetic())
# Plot the tracer point
ax_latlon.scatter(LON[ int(j_tracers[n]), int(i_tracers[n])], LAT[int(j_tracers[n]), int(i_tracers[n])], 
            color = 'red', 
            transform=ccrs.Geodetic(), label = "Tracer point")
# Get the triangle's B,C and D points which define the local grid CS
(BLat, BLon), (CLat, CLon), (DLat, DLon), _, _ = define_gridCS(fLAT, fLON, (int(j_tracers[n]),int(i_tracers[n])))
 # Create a matplotlib transform from the cartopy coordinate system
crs = ccrs.Geodetic()
transform = crs._as_mpl_transform(ax_latlon)
# Draw the x and y axis of reference for the local grid CS using B, C and D
ax_latlon.plot([BLon, CLon, DLon], [BLat, CLat, DLat], 
            color = 'green', 
            transform=ccrs.Geodetic(), label = "Axis of reference in the local grid CS")
# Add labels to the drawn axis of reference
ax_latlon.annotate("y", (BLon, BLat), xytext=(1, 1), xycoords=transform, textcoords='offset points', color="green") 
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