'''
Author: Beatrice Duval (bdu002)

------------------------------------------------------------------------------
Code that plots total sea-ice deformations, divergence and shear strain rates.
------------------------------------------------------------------------------

'''

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

import config
import utils_datetime
import utils_load_data as load_data

'''
_________________________________________________________________________________________
INITIALIZE PLOTS
'''

# Set the matplotlib projection and transform
proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
trans = ccrs.Geodetic()

# Create a colormap for the divergence ax in order that 0 is 
# associated to the same color as the other plots
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","#110788","yellow"])

# Initialize figures for total deformation (tot), divergence (I) and shear (II)
fig_tot = plt.figure()
fig_I = plt.figure()
fig_II = plt.figure()

# Initialize subplots
ax_tot = fig_tot.add_subplot(111, projection=proj)
ax_I = fig_I.add_subplot(111, projection=proj)
ax_II = fig_II.add_subplot(111, projection=proj)

# Create a list of axes to be iterated over
ax_list = [ax_tot, ax_I, ax_II]

for ax in ax_list:
    # Set the map extent in order to see the entire region of interest
    ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

# Iterate through data files in config
for raw_path, triangulated_path, calculated_path in zip(config.data_paths['raw'], config.data_paths['triangulated'], config.data_paths['calculations']):

    try:
        '''
        _________________________________________________________________________________________
        LOAD DATA
        '''
        # Load the raw data
        raw_data    = load_data.load_raw(raw_path)
        sLat        = raw_data['sLat']       # Starting latitudes
        sLon        = raw_data['sLon']       # Starting longitudes

        # Load triangulated data
        triangulated_data = load_data.load_triangulated(triangulated_path)
        vertice_idx1      = triangulated_data['vertice_idx1']   # Vertex indices in raw csv file
        vertice_idx2      = triangulated_data['vertice_idx2']
        vertice_idx3      = triangulated_data['vertice_idx3']

        # Load calculated data
        calculated_data = load_data.load_calculations(calculated_path)
        eps_I           = calculated_data['eps_I']    # Divergence rate
        eps_II          = calculated_data['eps_II']   # Maximum shear strain rate
        eps_tot         = calculated_data['eps_tot']  # Total deformation rate
    
    except:
        continue
    
    # Get start and end times as datetime objects
    start, end = utils_datetime.dataDatetimes(raw_path)

    '''
    _________________________________________________________________________________________
    PLOT DEFORMATIONS
    '''

    # Stack the 3 vertices index arrays
    triangles = np.stack((vertice_idx1, vertice_idx2, vertice_idx3), axis=-1)

    # Plot total deformations, divergence and shear
    cb_tot = ax_tot.tripcolor( sLon, sLat, triangles, facecolors=eps_tot, transform=trans, cmap='plasma', vmin=0, vmax=0.1 )
    cb_I = ax_I.tripcolor( sLon, sLat, triangles, facecolors=eps_I, transform=trans, cmap=cmap, vmin=-0.04, vmax=0.04 )
    cb_II = ax_II.tripcolor( sLon, sLat, triangles, facecolors=eps_II, transform=trans, cmap='plasma', vmin=0, vmax=0.1 )


'''
_________________________________________________________________________________________
ADD FEATURES TO THE PLOTS
'''

# Create a list of colorbars and titles to be iterated over
cb_list = [cb_tot, cb_I, cb_II]
title_list = ['Total Deformation Rate $(Days^{-1})$', 'Divergence Rate $(Days^{-1})$', 'Shear Rate $(Days^{-1})$']

# Iterate through all axes
for ax, title, cb in zip(ax_list, title_list, cb_list):
    # Add a colorbar
    plt.colorbar(cb, ax=ax)

    # Add a title
    ax.set_title(title)

    # Add gridlines
    ax.gridlines()

    # Hide deformations over land
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')


fig_tot.savefig('/home/bdu002/2021_SeaIceDeformation/figs/tot.png', bbox_inches='tight')
fig_I.savefig('/home/bdu002/2021_SeaIceDeformation/figs/I.png', bbox_inches='tight')
fig_II.savefig('/home/bdu002/2021_SeaIceDeformation/figs/II.png', bbox_inches='tight')
