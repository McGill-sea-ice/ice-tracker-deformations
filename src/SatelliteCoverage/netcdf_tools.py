"""
Author: Lekima Yakuden
GitHub: LekiYak

--------------------------------------------------------------------------------
Tools for analysing and processing netCDF files
--------------------------------------------------------------------------------

This file contains functions for analysing and processing netCDF files.

"""

from lib2to3.pytree import convert
from time import strftime

import cartopy.crs as ccrs
import matplotlib.tri as tri
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pyproj
import math
from netCDF4 import Dataset
from src.SatelliteCoverage.utils import convert_to_grid
from shapely.ops import transform
import os

from tqdm import tqdm

# IGNORING WARNINGS, COMMENT IF YOU WANT TO SEE THEM
import warnings
warnings.filterwarnings("ignore")

def plot_start_end_points(path:str):
    """
    Plots the start and end points of the triangle vertices on a map.
    Start points are blue, end points are red.

    INPUTS:
    path -- Path to netCDF {str}

    OUTPUTS:
    None -- Saves plot to the output directory

    """
    print('--- Plotting start and end points ---')

    # Load data from netCDF file
    data = load_netcdf(path)

    satellite = data['icetracker']
    start_lats = data['start_lats']
    start_lons = data['start_lons']
    end_lats = data['end_lats']
    end_lons = data['end_lons']

    """
    Plot Preamble
    """

    # Setting projection as North Pole Stereographic
    proj = ccrs.NorthPolarStereo(central_longitude=0)

    # Initialize figure
    fig = plt.figure(figsize=(5.5, 5.5), )
    ax = fig.add_subplot(projection = proj, frameon=False)

    # Upper (u) and lower (l) extents of x, y (metres)
    lxextent = -4400000
    uxextent = 2500000
    uyextent = 3500000
    lyextent = -2500000

    ax.set_extent((lxextent, uxextent, uyextent, lyextent), ccrs.NorthPolarStereo())

    # Show lat/lon grid and coastlines
    ax.gridlines(draw_labels=True)
    ax.coastlines()

    # Hide datapoints over land
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

    """
    Plotting
    """

    # Iterating plotting thrice to plot all three points of the triangles
    for i in range(2):

        # Converting start/end lat/lons to x/y (North pole stereographic)
        start_x, start_y = convert_to_grid(start_lons[i], start_lats[i])
        end_x, end_y = convert_to_grid(end_lons[i], end_lats[i])

        # Plotting start points (Blue)
        ax.scatter(start_x, start_y, color = 'blue', s = 0.01)

        # Plotting end points (Red)
        ax.scatter(end_x, end_y, color = 'red', s = 0.01)

    # Set title
    ax.set_title(f'Start and End points, {start_year}-{start_month}-{start_day}, {end_year}-{end_month}-{end_day}')

    # Set a directory to store figures
    figsPath =  output_folder + '/' + '/figs/'

    # Create directory if it doesn't exist
    os.makedirs(figsPath, exist_ok=True)

    # Set prefix for filename
    prefix = satellite + '_' + start_year + start_month + start_day + '_' + end_year + end_month + end_day + '_dt' + str(timestep) + '_tol' + str(tolerance)

    # Full path of figure
    fig_path = figsPath + prefix + '_start_end_points.png'

    # Check if the path exists and overwriting if it does
    if os.path.isfile(fig_path):
            os.remove(fig_path)

    # Saving figure
    plt.savefig(fig_path, bbox_inches='tight', dpi=600)

def recreate_coordinates(start_lat1, start_lat2, start_lat3, start_lon1, start_lon2, start_lon3, start_id1, start_id2, start_id3):
    """
    This function takes in a list of start lats/lons and corresponding start lat/lon IDs (triangle vertices)
    and outputs corresponding, larger lists of start lats/lons (with 0s in some indices) that reflect the
    coordinates' placement in the original data files for use with ax.tripcolor

    INPUTS:
    start_lat1,2,3 -- Arrays of starting latitudes {np.array, list}
    start_lons1,2,3 -- Arrays of starting longitudes {np.array, list}
    start_id1,2,3 -- Array of starting IDs corresponding to start_lats1,2,3 and start_lons1,2,3 {np.array, list}

    OUTPUTS:
    new_lat -- Array of latitude values at the positions they were orignally in, in the data file
    new_lon -- Array of longitude values at the positions they were originally in, in the data file

    """

    # Combined list of start IDs
    start_ids = np.hstack((start_id1, start_id2, start_id3))

    # Skipping blank data files
    if len(start_ids) == 0:
        return 0, 0

    # Initializing new lists of coordinates
    new_lon, new_lat = ([np.nan] * (max(start_ids) + 1) for i in range(2))

    for i in range(len(start_lat1)):
        new_lat[start_id1[i]] = (start_lat1[i])
        new_lat[start_id2[i]] = (start_lat2[i])
        new_lat[start_id3[i]] = (start_lat3[i])

    for i in range(len(start_lon1)):
        new_lon[start_id1[i]] = start_lon1[i]
        new_lon[start_id2[i]] = start_lon2[i]
        new_lon[start_id3[i]] = start_lon3[i]

    return new_lat, new_lon

def plot_deformations(path=None, data_in=None, config=None):
    """
    This function plots deformations from a netCDF file using matplotlib's ax.tripcolor.
    The function assumes that the netCDF was generated from src/SeaIceDeformation's M01_d03_compute_deformations.py.

    INUPTS:
    path -- Path to netCDF {str}

    OUTPUTS:
    None -- Saves divergence, shear, and vorticity plots to the output directory
    """

    # Loading data from netcdf as a dictionary
    if path != None and data_in == None :
        data = load_netcdf(path)
    elif path == None and data_in != None :
        data = data_in
    elif path == None and data_in == None or path != None and data != None :
        print('You need to give at least one path or data to plot and not both ')

    """
    Preamble
    """
    # Set the matplotlib projection and transform
    proj = ccrs.NorthPolarStereo(central_longitude=0)
    trans = ccrs.Geodetic()

    # Initialize figures for total deformation (tot), divergence (I) and shear (II)
    fig_div = plt.figure()
    fig_shr = plt.figure()
    fig_vrt = plt.figure()

    # Initialize subplots
    ax_div = fig_div.add_subplot(111, projection=proj)
    ax_shr = fig_shr.add_subplot(111, projection=proj)
    ax_vrt = fig_vrt.add_subplot(111, projection=proj)

    # Create a list of axes to be iterated over
    ax_list = [ax_div, ax_shr, ax_vrt]

    for ax in ax_list:
        # Set the map extent in order to see the entire region of interest
        ax.set_extent((-4400000, 2500000, 3500000, -2500000), ccrs.NorthPolarStereo())

    """
    Plotting
    """

    print('--- Creating sea-ice deformation figures ---')

    # Iterating over all files (Unique triangulations)
    for i in tqdm(np.unique(data['no'])):

        # Obtaining number of rows corresponding to triangles in given file that will be iterated over
        file_length = np.count_nonzero(data['no'] == i)

        # Setting maximum index
        min_index = np.where(data['no'] == i)[0][0]
        max_index = np.where(data['no'] == i)[0][-1]+1

        # Arranging triangle vertices in array for use in ax.tripcolor
        triangles = np.stack((data['idx1'][min_index:max_index], data['idx2'][min_index:max_index],
                    data['idx3'][min_index:max_index]), axis=-1)

        # Filtering data range to that of the current "file"
        start_lat1_temp, start_lat2_temp, start_lat3_temp = data['start_lat1'][min_index:max_index], \
            data['start_lat2'][min_index:max_index], data['start_lat3'][min_index:max_index]

        start_lon1_temp, start_lon2_temp, start_lon3_temp = data['start_lon1'][min_index:max_index], \
            data['start_lon2'][min_index:max_index], data['start_lon3'][min_index:max_index]

        start_lat, start_lon = recreate_coordinates(start_lat1_temp, start_lat2_temp, start_lat3_temp,
                                                        start_lon1_temp, start_lon2_temp, start_lon3_temp,
                                                        data['idx1'][min_index:max_index],
                                                        data['idx2'][min_index:max_index],
                                                        data['idx3'][min_index:max_index])

        # Extracting deformation data
        div_colours = data['div'][min_index:max_index]
        shr_colours = data['shr'][min_index:max_index]
        vrt_colours = data['vrt'][min_index:max_index]

        # tranform the coordinates already to improve the plot efficiency
        new_coords = ax_div.projection.transform_points(trans, np.array(start_lon), np.array(start_lat))

        # create one triangulation object
        tria = tri.Triangulation(new_coords[:,0], new_coords[:,1], triangles=triangles)

        if len(triangles) != 0:
        # Plotting
            cb_div = ax_div.tripcolor(tria, facecolors=div_colours, cmap='coolwarm', vmin=-0.04, vmax=0.04)
            cb_shr = ax_shr.tripcolor(tria, facecolors=shr_colours, cmap='plasma', vmin=0, vmax=0.1)
            cb_vrt = ax_vrt.tripcolor(tria, facecolors=vrt_colours, cmap='coolwarm', vmin=-0.1, vmax=0.1)


    # Create a list of colorbars and titles to be iterated over
    cb_list = [cb_div, cb_shr, cb_vrt]
    title_list = ['Divergence Rate $(Days^{-1})$', 'Shear Rate $(Days^{-1})$', 'Vorticity $(Days^{-1})$']

    Date_options = config['Date_options']
    start_year   = Date_options['start_year']
    end_year   = Date_options['end_year']
    start_month  = Date_options['start_month']
    end_month  = Date_options['end_month']
    start_day    = Date_options['start_day']
    end_day    = Date_options['end_day']
    timestep     = Date_options['timestep']

    IO            = config['IO']
    output_folder = IO['output_folder']
    exp           = IO['exp']

    # Iterate through all axes
    for ax, title, cb in zip(ax_list, title_list, cb_list):
        # Add a colorbar
        plt.colorbar(cb, ax=ax)

        # Add a title
        ax.set_title(title + '\n' + start_year + '-' + start_month + '-' + start_day + ' to ' +
                        end_year + '-' + end_month + '-' + end_day + ', ' + timestep + 'hr')

        # Add gridlines
        ax.gridlines()

        # Hide deformations over land
        ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

        '''
        _________________________________________________________________________________________
        SAVE PLOTS
        '''

    # Set a directory to store figures for the current experiment
    figsPath =  output_folder + '/' + exp + '/figs/'

    # Create the directory if it does not exist already
    os.makedirs(figsPath, exist_ok=True)

    itr = data.icetracker
    tol = data.tolerance
    # Create a prefix for the figure filenames
    prefix = itr + '_' + start_year + start_month + start_day + '_' + end_year + end_month + end_day + '_dt' + str(timestep) + '_tol' + str(tol)

    # Create the figure filenames
    div_path   = figsPath + prefix + '_div.png'
    shr_path = figsPath + prefix + '_shr.png'
    rot_path   = figsPath + prefix + '_vrt.png'


    for fig, fig_path in zip([fig_div, fig_shr, fig_vrt], [div_path, shr_path, rot_path]):

        # Check if the figures already exist. If they do, delete them.
        if os.path.isfile(fig_path):
            os.remove(fig_path)

        # Save the new figures
        fig.savefig(fig_path, bbox_inches='tight', dpi=600)

    return None

if __name__ == '__main__':
    # Reading config
    config = read_config()

    # Initializing more specific ConfigParser objects
    IO = config['IO']
    Date_options = config['Date_options']
    options = config['options']
    meta = config['Metadata']
    netcdf_tools = config['netcdf_tools']

    path = IO['netcdf_path']
    output_folder = IO['output_folder']

    start_year = Date_options['start_year']
    start_month = Date_options['start_month']
    start_day = Date_options['start_day']

    end_year = Date_options['end_year']
    end_month = Date_options['end_month']
    end_day = Date_options['end_day']

    timestep = Date_options['timestep']
    tolerance = Date_options['tolerance']

    area_filter = options['area_filter']
    centre_lat = options['centre_lat']
    centre_lon = options['centre_lon']
    radius = options['radius']

    if netcdf_tools['plot_start_end_points'] == 'True':
        plot_start_end_points(path)

    if netcdf_tools['plot_deformation'] == 'True':
        plot_deformations(path)
