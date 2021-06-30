'''
Author: Beatrice Duval (bdu002)

-------------------------------------------
Code that plots total sea-ice deformations.
-------------------------------------------

Parameters:
draw_file_index -- (True, False) If set to True, the index (in the config.py lists) 
                   of each plotted sets of data will be drawn on the map
draw_single_csv -- (True, False) If set to True, each set of data in the config.py 
                   will be plotted individually and shown one at the time. If set to 
                   False, the sets of data in the config.py file will be plotted 
                   simultaneously on the same map.

'''
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import src.config
from src.d00_utils.datetime_raw_csv import dataDatetimes
from src.d00_utils.load_csv import (load_calculations_csv, load_converted_csv,
                                    load_processed_csv, load_raw_csv)

'''
_________________________________________________________________________________________
PARAMETERS
'''

draw_file_index = False
draw_single_csv = True

'''
_________________________________________________________________________________________
PLOT
'''

# Initialize the config global variables (i.e. .csv file paths for all stages of data processing)
src.config.init()

if not draw_single_csv:
    # Initialize figure and add subplots
    fig = plt.figure()

    # Set the matplotlib projection and transform
    proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
    trans = ccrs.Geodetic()

    # Initialize a subplot for a linear representation of total deformation
    ax = fig.add_subplot(111, projection=proj)

# Iterate through all raw and processed .csv files listed in config.py
for raw_csv_path, processed_csv_path, converted_csv_path, calculations_csv_path, n in zip(src.config.raw_csv_paths, src.config.processed_csv_paths, src.config.converted_csv_paths, src.config.calculations_csv_paths, range(len(src.config.calculations_csv_paths))):

    try:
                
        '''
        _________________________________________________________________________________________
        LOAD DATA SET
        '''

        # Load a raw dataset
        sLat, sLon, _, _ = load_raw_csv( raw_csv_path )

        # Load a processed dataset
        _, _, _, _, _, _, _, vertice_idx1, vertice_idx2, vertice_idx3 = load_processed_csv( processed_csv_path )

        # Load a converted dataset
        _, _, _, _, _, _, _, _, _, _, _, _, _, j_tracers, i_tracers = load_converted_csv( converted_csv_path )

        # Load a calculations dataset
        _, _, _, _, _, _, _, eps_tot = load_calculations_csv( calculations_csv_path )

    except:
        continue

    # Retrieve start and end times
    start, end = dataDatetimes(raw_csv_path)

    

    if draw_single_csv:
        # Initialize figure and add subplots
        fig = plt.figure()

        # Set the matplotlib projection and transform
        proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
        trans = ccrs.Geodetic()

        # Initialize a subplot for a linear representation of total deformation
        ax = fig.add_subplot(111, projection=proj)

    # Stack the 3 vertices index arrays
    triangles = np.stack((vertice_idx1, vertice_idx2, vertice_idx3), axis=-1)

    # Plot the grid points and color them as a function of
    # their total deformation rate
    cb_epsTot_lin = ax.tripcolor( sLon, sLat, triangles, facecolors=eps_tot, transform=trans, cmap = 'plasma', vmin=0, vmax=0.1 )
    
    # Create a mpl transform
    crs = ccrs.Geodetic()
    transform = crs._as_mpl_transform(ax)

    if draw_file_index:
        # Write the index of the data set in the .csv file paths
        ax.annotate(str(n), (sLon[int(len(sLon)/2)], sLat[int(len(sLon)/2)]), xytext=(1, 1), xycoords=transform, textcoords='offset points', fontsize=16) 

    if draw_single_csv:
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
        plt.title(start.strftime("%b %d %Y %H:%M:%S") + " - " + end.strftime("%b %d %Y %H:%M:%S") + '\n(' + os.path.basename(raw_csv_path) + ')', fontsize=10, y=1)

        plt.show()

if not draw_single_csv:
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

    plt.show()
