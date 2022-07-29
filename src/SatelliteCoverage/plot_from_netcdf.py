from lib2to3.pytree import convert
from time import strftime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import math
from netCDF4 import Dataset
from config import *
from coverage_frequency_map import convert_to_grid
from shapely.ops import transform
import os

import warnings
warnings.filterwarnings("ignore")

def seconds_to_date(path:str):

    ds = Dataset(path, mode='r')

    reftime = ds.getncattr('referenceTime')

    reftime = reftime[0:10]
    reftime = datetime.strptime(reftime, '%Y-%m-%d')
    
    ds.close()

    start_date = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
    end_date = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

    start_ref_dt = (start_date - reftime).total_seconds()
    end_ref_dt = (end_date - reftime).total_seconds()

    return start_ref_dt, end_ref_dt

def load_netcdf(path:str):
    # Load netCDF as Dataset from *path*
    ds = Dataset(path, mode='r')

    # Start and end times
    start_date_str = start_year + start_month + start_day
    end_date_str = end_year + end_month + end_day
    start_date = datetime.strptime(start_date_str, '%Y%m%d')
    end_date = datetime.strptime(end_date_str, '%Y%m%d')

    # Get start / end times from user
    start_time_s, end_time_s = seconds_to_date(path)

    # Extracting time variables
    start_time = ds.variables['start_time'][:]
    end_time = ds.variables['end_time'][:]

    # Indices of data in desired time frame
    time_indices = np.where( (start_time > start_time_s) & (start_time < end_time_s) )[0]

    # Extracting data
    start_lat1 = (ds.variables['start_lat1'][:])[time_indices]
    start_lat2 = (ds.variables['start_lat2'][:])[time_indices]
    start_lat3 = (ds.variables['start_lat3'][:])[time_indices]

    start_lats = np.array([start_lat1, start_lat2, start_lat3])

    start_lon1 = (ds.variables['start_lon1'][:])[time_indices]
    start_lon2 = (ds.variables['start_lon2'][:])[time_indices]
    start_lon3 = (ds.variables['start_lon3'][:])[time_indices]

    start_lons = np.array([start_lon1, start_lon2, start_lon3])

    end_lat1 = (ds.variables['end_lat1'][:])[time_indices]
    end_lat2 = (ds.variables['end_lat2'][:])[time_indices]
    end_lat3 = (ds.variables['end_lat3'][:])[time_indices]

    end_lats = np.array([end_lat1, end_lat2, end_lat3])

    end_lon1 = (ds.variables['end_lon1'][:])[time_indices]
    end_lon2 = (ds.variables['end_lon2'][:])[time_indices]
    end_lon3 = (ds.variables['end_lon3'][:])[time_indices]

    end_lons = np.array([end_lon1, end_lon2, end_lon3])

    div = (ds.variables['div'][:])[time_indices]
    shear = (ds.variables['shear'][:])[time_indices]
    vrt = (ds.variables['vrt'][:])[time_indices]

    idx1 = (ds.variables['idx1'][:])[time_indices]
    idx2 = (ds.variables['idx2'][:])[time_indices]
    idx3 = (ds.variables['idx3'][:])[time_indices]
    no = (ds.variables['no'][:])[time_indices]

    id_start_lat1 = (ds.variables['id_start_lat1'][:])[time_indices]
    id_start_lat2 = (ds.variables['id_start_lat2'][:])[time_indices]
    id_start_lat3 = (ds.variables['id_start_lat3'][:])[time_indices]

    reftime = ds.getncattr('referenceTime')

    # Closing dataset
    ds.close()

    return {'start_lats': start_lats, 'start_lons': start_lons, 'end_lats': end_lats, 'end_lons': end_lons, 'div': div, 'shear': shear, 'vrt': vrt, 'start_time': start_time, 'end_time': end_time, 'time_indices': time_indices, 'reftime': reftime, 'idx1': idx1, 'idx2': idx2, 'idx3': idx3, 'no': no, 'start_id1': id_start_lat1, 'start_id2': id_start_lat2, 'start_id3': id_start_lat3}

def plot_start_end_points(path:str):

    start_lats = load_netcdf(path)['start_lats']
    start_lons = load_netcdf(path)['start_lons']
    end_lats = load_netcdf(path)['end_lats']
    end_lons = load_netcdf(path)['end_lons']

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
    plt.suptitle('Start and end points, 2021-02-01 to 2021-02-28 \nBlue points: Start points. Red points: End points.')

    # Saving figure
    plt.savefig('20202021_tripoints_month_72.png')

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
    
    # Loading start lats/lons
    # start_lat1 = list(start_lat1)
    # start_lat2 = list(start_lat2)
    # start_lat3 = list(start_lat3)

    # start_lon1 = list(start_lon1)
    # start_lon2 = list(start_lon2)
    # start_lon3 = list(start_lon3)

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
    
def plot_deformations(path:str):
    """
    This function plots deformations from a netCDF file using matplotlib's ax.tripcolor. 
    The function assumes that the netCDF was generated from src/SeaIceDeformation's M01_d03_compute_deformations.py.

    INUPTS:
    path -- Path to netCDF {str}

    OUTPUTS:
    none
    Saves divergence, shear, and vorticity plots.
    """

    # Loading data from netcdf as a dictionary
    data = load_netcdf(path)

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

    # Obtaining start index of data
    min_no = np.nanmin(data['no']) # Minimum (smallest) ID number
    min_index = np.nanmin(np.where(data['no'] == min_no)) # This value will be updated for each triangle
    
    # Iterating over all files (Unique triangulations)
    for i in range(len(set(data['no']))):
        # Here the file ID is zero-indexed, so file 1 has the ID 0

        # Obtaining number of rows to iterate over
        file_length = np.count_nonzero(data['no'] == (i))

        # Setting maximum index
        max_index = min_index + file_length

        # Arranging triangle vertices in array for use in ax.tripcolor
        triangles = np.stack((data['idx1'][min_index:max_index], data['idx2'][min_index:max_index], 
                    data['idx3'][min_index:max_index]), axis=-1)

        # Filtering data range to that of the current "file"

        start_lat1_temp, start_lat2_temp, start_lat3_temp = data['start_lats'][0][min_index:max_index], \
            data['start_lats'][1][min_index:max_index], data['start_lats'][2][min_index:max_index]
        
        start_lon1_temp, start_lon2_temp, start_lon3_temp = data['start_lons'][0][min_index:max_index], \
            data['start_lons'][1][min_index:max_index], data['start_lons'][2][min_index:max_index]

        start_lat, start_lon = recreate_coordinates(start_lat1_temp, start_lat2_temp, start_lat3_temp, 
                                                        start_lon1_temp, start_lon2_temp, start_lon3_temp, 
                                                        data['start_id1'][min_index:max_index],
                                                        data['start_id2'][min_index:max_index], 
                                                        data['start_id3'][min_index:max_index])

        # Extracting deformation data
        div_colours = data['div'][min_index:max_index]
        shr_colours = data['shear'][min_index:max_index]
        vrt_colours = data['vrt'][min_index:max_index]

        if len(triangles) != 0:
        # Plotting
            cb_div = ax_div.tripcolor(start_lon, start_lat, triangles, transform=trans, facecolors=div_colours, cmap='coolwarm', vmin=-0.04, vmax=0.04)
            cb_shr = ax_shr.tripcolor(start_lon, start_lat, triangles, transform=trans, facecolors=shr_colours, cmap='plasma', vmin=0, vmax=0.1)
            cb_vrt = ax_vrt.tripcolor(start_lon, start_lat, triangles, transform=trans, facecolors=vrt_colours, cmap='coolwarm', vmin=-0.1, vmax=0.1)

        # Updating minimum index
        min_index = max_index

    # Create a list of colorbars and titles to be iterated over
    cb_list = [cb_div, cb_shr, cb_vrt]
    title_list = ['Divergence Rate $(Days^{-1})$', 'Shear Rate $(Days^{-1})$', 'Rotation Rate $(Days^{-1})$']

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
    figsPath =  output_folder + '/' + '/figs/'
        
        # Create the directory if it does not exist already
    os.makedirs(figsPath, exist_ok=True)

        # Create a prefix for the figure filenames
    prefix = start_year + start_month + start_day + '_' + str(timestep) 
        
        # Create the figure filenames
    div_path   = figsPath + prefix + '_div.png'
    shear_path = figsPath + prefix + '_shear.png'
    rot_path   = figsPath + prefix + '_rot.png'

        
    for fig, fig_path in zip([fig_div, fig_shr, fig_vrt], [div_path, shear_path, rot_path]):
            
        # Check if the figures already exist. If they do, delete them.
        if os.path.isfile(fig_path):
            os.remove(fig_path)

        # Save the new figures
        fig.savefig(fig_path, bbox_inches='tight', dpi=300)

if __name__ == '__main__':
    # Reading config
    config = read_config()

    # Initializing more specific ConfigParser objects
    IO = config['IO']
    options = config['options']

    path = IO['netcdf_folder']
    output_folder = IO['output_folder']

    start_year = options['start_year']
    start_month = options['start_month']
    start_day = options['start_day']

    end_year = options['end_year']
    end_month = options['end_month']
    end_day = options['end_day']

    timestep = options['timestep']

    #plot_start_end_points(path)
    plot_deformations(path)
