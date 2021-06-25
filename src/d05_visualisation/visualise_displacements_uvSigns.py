'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------------------------------------------
Code that plots data points' displacement and the sign of their velocity components.
-----------------------------------------------------------------------------------

'''

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from src.d01_data.load01_raw_csv import sLat, sLon, eLat, eLon
from src.d00_utils.datetime_raw_csv import dataDatetimes, dT
from src.d01_data.get_data_paths import calculations_csv_path
from src.d00_utils.grid_coord_system import (define_gridCS,
                                             find_nearestGridTracerPt,
                                             get_xy_gridCS)
from src.d00_utils.deformation_computations import (calculate_uv_lists)
from pyproj import Proj
from src.d01_data.load00_grid import X_aeqd, Y_aeqd, fLAT, fLON
from src.d01_data.load03_converted_csv import i_tracers, j_tracers

# Get start and end times as datetime objects
start, end = dataDatetimes(calculations_csv_path)

# Compute the time interval (days)
dt = dT( (start, end) )

# Set the matplotlib projection and transform
proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
trans = ccrs.Geodetic()

# Create a projection object in the aeqd transform
p = Proj(proj='aeqd', ellps='WGS84', preserve_units=False)

# Initialize a list to store starting and ending x and y positions
sx_list = []
sy_list = []
ex_list = []
ey_list = []

# Compute the speed of every data point in a local grid CS
for n in range(len(sLat)):
    
    # Retrieve the data point starting and ending positions
    nsLat, nsLon, neLat, neLon =  sLat[n], sLon[n], eLat[n], eLon[n]

    # Express these coordinates in the aeqd transform
    nXYLat_aeqd = p(nsLon, nsLat)

    # Express these coordinates in a local cartesian grid CS
    # Find the nearest grid tracer point
    ji = find_nearestGridTracerPt(X_aeqd, Y_aeqd, nXYLat_aeqd)
    j, i = ji

    # Get the x, y coords in the local grid CS
    sx, sy = get_xy_gridCS( define_gridCS(fLAT, fLON,  ji), (nsLat, nsLon))
    ex, ey = get_xy_gridCS( define_gridCS(fLAT, fLON,  ji), (neLat, neLon))
 
    # Add starting and ending x and y positions to the lists
    sx_list.append(sx)
    sy_list.append(sy)
    ex_list.append(ex)
    ey_list.append(ey)

    
# Compute the u and v speed components
u_list, v_list = calculate_uv_lists( sx_list, sy_list, sy_list, ey_list, dt)    

# Initialize figure
fig = plt.figure()

# Initialize subplots for u and v signs
ax_u = fig.add_subplot(121, projection=proj)
ax_v = fig.add_subplot(122, projection=proj)

for ax, speed in zip([ax_u, ax_v], [u_list, v_list]):
    # Create a matplotlib transform from the cartopy coordinate system
    as_mpl_transform = trans._as_mpl_transform(ax)

    # Create arrays that contain only data points that have a positive, 
    # negative or 0 velocity
    positiveLat = sLat.copy()                       # Positive velocities data points
    positiveLon = sLon.copy()
    positiveLat[np.array(speed) <= 0] = np.nan
    positiveLon[np.array(speed) <= 0] = np.nan

    negativeLat = sLat.copy()                       # Negative velocities data points
    negativeLon = sLon.copy()
    negativeLat[np.array(speed) >= 0] = np.nan
    negativeLon[np.array(speed) >= 0] = np.nan

    zeroLat = sLat.copy()                           # Zero velocities data points
    zeroLon = sLon.copy()
    zeroLat[np.array(speed) != 0] = np.nan
    zeroLon[np.array(speed) != 0] = np.nan

    # Plot the sign of the speed component
    ax.scatter( negativeLon, negativeLat, color='blue', transform=trans )
    ax.scatter( zeroLon, zeroLat, color='gray', transform=trans )
    ax.scatter( positiveLon, positiveLat, color='red', transform=trans)


    # Plot the data points displacements
    for n in range(len(sLat)):
        ax.annotate('', xy = (eLon[n], eLat[n]),  xycoords = as_mpl_transform, \
        xytext = (sLon[n], sLat[n]), textcoords = as_mpl_transform, fontsize = 7, \
        color = '#303030', arrowprops=dict(edgecolor='black', arrowstyle = '->'))

    # Add a reference xy axis
    # Find the max i,j grid point
    max_j, max_i = int(max(j_tracers)), int(max(i_tracers))

    # Plot the y axis
    ax.annotate('y', xy = (fLON[max_j-1, max_i-1], fLAT[max_j-1, max_i-1]),  xycoords = as_mpl_transform, \
                            xytext = (fLON[max_j+10, max_i-1], fLAT[max_j+10, max_i-1]), textcoords = as_mpl_transform, fontsize = 7, \
                            color = 'red', arrowprops=dict(edgecolor='red', arrowstyle = '<-'))

    # Plot the x axis
    ax.annotate('x', xy = (fLON[max_j-1, max_i-1], fLAT[max_j-1, max_i-1]),  xycoords = as_mpl_transform, \
                            xytext = (fLON[max_j-1, max_i+10], fLAT[max_j-1, max_i+10]), textcoords = as_mpl_transform, fontsize = 7, \
                            color = 'red', arrowprops=dict(edgecolor='red', arrowstyle = '<-'))


# Add a title and subtitles
plt.suptitle('Speed components colored as a function of their sign \n(red: -, blue: +, gray: 0)\n' + start.strftime("%b %d %Y %H:%M:%S") + " - " + end.strftime("%b %d %Y %H:%M:%S"), fontsize=14, y=1)
ax_u.set_title('u-component', fontsize=8, y=1)
ax_v.set_title('v-component', fontsize=8, y=1)

# Show plot
plt.show()