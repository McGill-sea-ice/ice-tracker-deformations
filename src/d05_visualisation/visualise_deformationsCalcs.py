'''
Author: Beatrice Duval (bdu002)

----------------------------------------------------------------------
Code that plots results from sea-ice deformations calculations 
(strain rates, total deformation, triangle area, and displacement).
----------------------------------------------------------------------

'''
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from src.d00_utils.datetime_raw_csv import dataDatetimes
from src.d00_utils.load_csv import load_raw_csv, load_processed_csv, load_converted_csv, load_calculations_csv
from src.d01_data.load00_grid import fLAT, fLON
import src.config

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
converted_csv_path = src.config.converted_csv_paths[n]
calculations_csv_path = src.config.calculations_csv_paths[n]

# Load a raw dataset
sLat, sLon, eLat, eLon = load_raw_csv( raw_csv_path )

# Load a processed dataset
_, _, _, _, _, _, _, vertice_idx1, vertice_idx2, vertice_idx3 = load_processed_csv( processed_csv_path )

# Load a converted dataset
_, _, _, _, _, _, _, _, _, _, _, _, _, j_tracers, i_tracers = load_converted_csv( converted_csv_path )

# Load a calculations dataset
_, dudx, dudy, dvdx, dvdy, eps_I, eps_II, eps_tot = load_calculations_csv( calculations_csv_path )

# Get start and end times as datetime objects
start, end = dataDatetimes(raw_csv_path)


'''
_________________________________________________________________________________________
PLOT
'''

# Initialize figure and add subplots
fig = plt.figure()

# If plot_obs_pt is set to True, an observation point will be circled in red on all plots
plot_obs_pt = False

# Set an observation point
#obs_lon, obs_lat = -172.068658, 72.776678   # Low deformation point -> displacement arrows same direction
#obs_lon, obs_lat = -174.838872, 73.905923   # High deformation point (~0.27 days^-1) -> displacement arrows have different directions
#obs_lon, obs_lat = -171.456486, 74.080925   # Highest deformation point (~16 days^-1) -> area ~ 0, is on the side
#obs_lon, obs_lat = -168.537858, 73.781656   # Highest (2) deformation point (~7 days^-1) -> very narrow triangle, area ~ 0, is on the side
#obs_lon, obs_lat = -175.581566, 71.731767    # High deformation point (~0.65 days^-1) -> displacement arrows have different directions
obs_lon, obs_lat = -177.709838, 73.513994    # High deformation point (~0.85days^-1) -> area ~ 0, is on the side, vertices are far appart

# Set the matplotlib projection and transform
proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
trans = ccrs.Geodetic()

#------------------ ROW 1 ----------------------------------------

#---- COL 1 - Displacements and triangulation

# Initialize the subplot
ax_displacement = fig.add_subplot(241, projection = proj)

# Add coastlines and gridlines
ax_displacement.coastlines()
ax_displacement.gridlines()

# Create a matplotlib transform from the cartopy coordinate system
as_mpl_transform = trans._as_mpl_transform(ax_displacement)

for n in range(len(sLat)):

    # Plot the data points displacements
    ax_displacement.annotate('', xy = (eLon[n], eLat[n]),  xycoords = as_mpl_transform, \
    xytext = (sLon[n], sLat[n]), textcoords = as_mpl_transform, fontsize = 7, \
    color = '#303030', arrowprops=dict(edgecolor='black', arrowstyle = '->'))
    
for n in range( len(vertice_idx1)):
    # Plot the triangles
    ax_displacement.plot([sLon[vertice_idx1[n]], sLon[vertice_idx1[n]], sLon[vertice_idx1[n]], sLon[vertice_idx1[n]]],
                         [sLat[vertice_idx1[n]], sLat[vertice_idx1[n]], sLat[vertice_idx1[n]], sLat[vertice_idx1[n]]], 
                         color = 'xkcd:sky blue', 
                         transform=trans)

# Plot observation point
if plot_obs_pt:
    ax_displacement.scatter(obs_lon, obs_lat, transform=trans, facecolors='none', edgecolors='r')

# Add a reference xy axis
# Find the max i,j grid point
max_j, max_i = int(max(j_tracers)), int(max(i_tracers))

# Plot the y axis
ax_displacement.annotate('y', xy = (fLON[max_j-1, max_i-1], fLAT[max_j-1, max_i-1]),  xycoords = as_mpl_transform, \
                        xytext = (fLON[max_j+10, max_i-1], fLAT[max_j+10, max_i-1]), textcoords = as_mpl_transform, fontsize = 7, \
                        color = 'red', arrowprops=dict(edgecolor='red', arrowstyle = '<-'))

# Plot the x axis
ax_displacement.annotate('x', xy = (fLON[max_j-1, max_i-1], fLAT[max_j-1, max_i-1]),  xycoords = as_mpl_transform, \
                        xytext = (fLON[max_j-1, max_i+10], fLAT[max_j-1, max_i+10]), textcoords = as_mpl_transform, fontsize = 7, \
                        color = 'red', arrowprops=dict(edgecolor='red', arrowstyle = '<-'))

# Add a title
ax_displacement.set_title('Delaunay Triangulation and \nData Points Displacements', fontsize=8, y=1)


#---- COL 2 - Total deformation

# Stack the 3 vertices index arrays
triangles = np.stack((vertice_idx1, vertice_idx2, vertice_idx3), axis=-1)

# Initialize a subplot for a linear representation of total deformation
ax_epsTot = fig.add_subplot(242, projection=proj)

# Set the map extent in order to see the entire region of interest
ax_epsTot.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

# Plot the grid points and color them as a function of
# their total deformation rate
cb_epsTot = ax_epsTot.tripcolor( sLon, sLat, triangles, facecolors=eps_tot, transform=trans, cmap = 'plasma', vmin=0, vmax=0.1 )
ax_epsTot.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

# Add coastlines and gridlines
ax_epsTot.coastlines()
ax_epsTot.gridlines()

# Plot observation point
if plot_obs_pt:
    ax_epsTot.scatter(obs_lon, obs_lat, transform=ccrs.Geodetic(), facecolors='none', edgecolors='r')

# Add a title
ax_epsTot.set_title('Total Deformation Rates ($days^{-1}$)', fontsize=8, y=1)


#---- COL 3 - Divergence (eps_I)

# Initialize a subplot for a linear representation of total deformation
ax_epsI = fig.add_subplot(243, projection=proj)

# Set the map extent in order to see the entire region of interest
ax_epsI.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

# Plot the grid points and color them as a function of
# their total deformation rate
cb_epsI = ax_epsI.tripcolor( sLon, sLat, triangles, facecolors=eps_I, transform=trans, cmap = 'plasma', vmin=0, vmax=0.1 )

# Add coastlines and gridlines
ax_epsI.coastlines()
ax_epsI.gridlines()

# Plot observation point
if plot_obs_pt:
    ax_epsI.scatter(obs_lon, obs_lat, transform=ccrs.Geodetic(), facecolors='none', edgecolors='r')

# Add a title
ax_epsI.set_title('Divergence Strain Rates ($days^{-1}$)', fontsize=8, y=1)


#---- COL 4 - Maximum shear (eps_II)

# Initialize a subplot for a linear representation of total deformation
ax_epsII = fig.add_subplot(244, projection=proj)

# Set the map extent in order to see the entire region of interest
ax_epsII.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

# Plot the grid points and color them as a function of
# their total deformation rate
cb_epsII = ax_epsII.tripcolor( sLon, sLat, triangles, facecolors=eps_II, transform=trans, cmap = 'plasma', vmin=0, vmax=0.1 )

# Add coastlines and gridlines
ax_epsII.coastlines()
ax_epsII.gridlines()

# Add a colorbar and include the displacements plot in order 
# to have it in the same size as the rest of the plots
plt.colorbar(cb_epsII, ax=[ax_epsTot, ax_epsI, ax_epsII])

# Plot observation point
if plot_obs_pt:
    ax_epsII.scatter(obs_lon, obs_lat, transform=ccrs.Geodetic(), facecolors='none', edgecolors='r')

# Add a title
ax_epsII.set_title('Maximum Shear Strain Rates ($days^{-1}$)', fontsize=8, y=1)



#------------------ ROW 2 ----------------------------------------

# Initialize strain rate subplots
ax_dudx = fig.add_subplot(245, projection = proj)
ax_dudy = fig.add_subplot(246, projection = proj)
ax_dvdx = fig.add_subplot(247, projection = proj)
ax_dvdy = fig.add_subplot(248, projection = proj)

# Create a list of subplot axes and strain rates to plot
strain_axes = [ax_dudx, ax_dudy, ax_dvdx, ax_dvdy]
strain_rates_list = [dudx, dudy, dvdx, dvdy]
strain_strings = ['du/dx', 'du/dy', 'dv/dx', 'dv/dy']

# Iterate through strain rates plot
for (ax, strain_rates, strain_string) in zip(strain_axes, strain_rates_list, strain_strings):

    # Plot the grid points and color them as a function of the current strain rate
    cb = ax.tripcolor( sLon, sLat, triangles, facecolors=strain_rates, transform=trans, cmap = 'plasma', vmin=-0.2, vmax=0.2 )

    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()
    
    # Set the map extent in order to see the entire region of interest
    ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

    # Plot observation point
    if plot_obs_pt:
        ax.scatter(obs_lon, obs_lat, transform=ccrs.Geodetic(), facecolors='none', edgecolors='r')

    # Add a title
    ax.set_title('Strain Rates ($days^{-1}$) \n - ' + strain_string + ' -', fontsize=8, y=1)

# Add a colorbar for all strain rates
plt.colorbar(cb, ax=strain_axes)


#------------------ Add main titles and show the plot ----------------------------------------

# Add main title indicating starting and ending times
plt.suptitle( start.strftime("%b %d %Y %H:%M:%S") + " - " + end.strftime("%b %d %Y %H:%M:%S"), fontsize=14, y=1)

plt.show()