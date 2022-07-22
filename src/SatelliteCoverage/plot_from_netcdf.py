from lib2to3.pytree import convert
from time import strftime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon
from matplotlib.collections import PathCollection
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import math
from netCDF4 import Dataset
from config import *
from coverage_frequency_map import convert_to_grid
import heapq
import xarray

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

    reftime = ds.getncattr('referenceTime')

    # Closing dataset
    ds.close()

    return {'start_lats': start_lats, 'start_lons': start_lons, 'end_lats': end_lats, 'end_lons': end_lons, 'div': div, 'shear': shear, 'vrt': vrt, 'start_time': start_time, 'end_time': end_time, 'time_indices': time_indices, 'reftime': reftime, 'idx1': idx1, 'idx2': idx2, 'idx3': idx3, 'no': no}

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

def plot_deformations_polygons(path:str):

    # Loading data from netcdf as a dictionary
    data = load_netcdf(path)

    # Initializing patches list
    patches = []

    # Initializing colormaps
    cmap_div = matplotlib.cm.get_cmap('coolwarm')
    cmap_shr = matplotlib.cm.get_cmap('plasma')
    cmap_rot = matplotlib.cm.get_cmap('coolwarm')

    # Initializing vertices list
    vertices = []

    for i in range(len(data['div'])):

        # Initializing temporary vertices list
        temp_vertices = []

        # Iterating over all vertices in triangle
        for j in range(3):
            temp_vertices.append((data['start_lats'][j][i], data['start_lons'][j][i]))

        # Appending triangle vertices to main list
        vertices.append(temp_vertices)
    
    polygon = Polygon(vertices[0])
    
    norm_div = matplotlib.colors.Normalize(vmin=-0.04, vmax=0.04)

def plot_deformations(path:str):
    # Loading data from netcdf as a dictionary
    data = load_netcdf(path)



    """
    Preamble
    """
    # Set the matplotlib projection and transform
    proj = ccrs.NorthPolarStereo(central_longitude=0)

    # Initialize figures for total deformation (tot), divergence (I) and shear (II)
    fig_div = plt.figure()
    fig_shr = plt.figure()
    fig_rot = plt.figure()

    # Initialize subplots
    ax_div = fig_div.add_subplot(111, projection=proj)
    ax_shr = fig_shr.add_subplot(111, projection=proj)
    ax_rot = fig_rot.add_subplot(111, projection=proj)

    # Create a list of axes to be iterated over
    ax_list = [ax_div, ax_shr, ax_rot]

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
        # Here i + 1 is the ID number of the file

        # Obtaining number of rows to iterate over
        file_length = np.count_nonzero(data['no'] == (min_no + i))

        # Setting maximum index
        max_index = min_index + file_length

        # Arranging triangle vertices in array for use in ax.tripcolor
        triangles = np.stack((data['idx1'][min_index:max_index], data['idx2'][min_index:max_index], 
                    data['idx3'][min_index:max_index]), axis=-1)
        
        # Initializing start lat/lon lists
        start_lats = []
        start_lons = []

        # Extracting lat/lons associated with the i-th triangulation
        for j in range(3):
            start_lats.append(data['start_lats'][j][min_index:max_index])
            start_lons.append(data['start_lons'][j][min_index:max_index])

        # Stacking vertices
        start_lats = np.stack((start_lats[0], start_lats[1], start_lats[2]), axis=-1)
        start_lons = np.stack((start_lons[0], start_lons[1], start_lons[2]), axis=-1)

        # Flattening lists
        start_lats = [item for sublist in start_lats for item in sublist]
        start_lons = [item for sublist in start_lons for item in sublist]

        # Converting to polar stereographic coordinates
        start_x, start_y = convert_to_grid(start_lons, start_lats)

        # Extracting div data
        colours = data['div'][min_index:max_index]

        # Printing debugging info
        print(f'Indices: {min_index} ~ {max_index}')
        print(f'Triangles: {len(triangles)}')
        print(f'Colours: {len(colours)}')
        print(triangles)

        # Plotting
        cb_div = ax_div.tripcolor(start_x, start_y, triangles, facecolor=colours, cmap='coolwarm', vmin=-0.04, vmax=0.04)

        # Updating minimum index
        min_index = max_index
    
    for fig in [fig_div, fig_shr, fig_rot]:
        plt.savefig(f"{fig}PLOTDEFORMATOINTEST.png")

def test_function(path:str):

    # Loading data
    data = load_netcdf(path)

    cmap_div = matplotlib.cm.get_cmap('coolwarm')
    cmap_shr = matplotlib.cm.get_cmap('plasma')
    cmap_rot = matplotlib.cm.get_cmap('coolwarm')

def tri_area(n):
    # Iterating over every triangle
    for i in n:
        
        #Initializing list of cartesian vertices
        vertices_xy = []

        # Iterating over each vertex
        for j in [0, 1, 2]:
            x, y = convert_to_grid(vertices[j][0][i], vertices[j][1][i])

            # Appending cartesian coordinates to list
            vertices_xy.append(x) 
            vertices_xy.append(y)

        # Now, vertices_xy = [x1, y1, x2, y2, x3, y3]  
        x1 = vertices_xy[0]
        y1 = vertices_xy[1]
        x2 = vertices_xy[2]
        y2 = vertices_xy[3]
        x3 = vertices_xy[4]
        y3 = vertices_xy[5]

        # Calculating area
        area = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

        print(f"{i} / {len(start_lat1) // 1000}")

        # Appending to main list
        areas.append(area)

        if i == (len(start_lat1) // 1000 - 1):
        
            max_areas = heapq.nlargest(10, areas)

            print(max_areas)

    """
    At this point the list tri_xy should look like this:
    [[(startx1, starty1), (startx2, starty2), (startx3, starty3)], [(startx1, starty1), (startx2, starty2), (startx3, starty3)],  ... ] where each list (element n in the list) contains tuples of the cartesian coordinates of each vertex
    in the triangle. e.g. tri_xy[n] = [(startx1, starty1), (startx2, starty2), (startx3, starty3)]
    """


if __name__ == '__main__':
    # Reading config
    config = read_config()

    # Initializing more specific ConfigParser objects
    IO = config['IO']
    options = config['options']

    path = IO['netcdf_folder']

    start_year = options['start_year']
    start_month = options['start_month']
    start_day = options['start_day']

    end_year = options['end_year']
    end_month = options['end_month']
    end_day = options['end_day']

    #plot_start_end_points(path)
    #test_function(path)
    #plot_deformations(path)
    #tri_area(path)[0]
    #seconds_to_date(path)
    plot_deformations(path)

    """
    Loading Data
    ------------------------------------------------------------------------------------------
    """

    # # Load netCDF as Dataset from *path*
    # ds = Dataset(path, mode='r')

    # # Extracting start/end points
    # start_lat1 = ds.variables['start_lat1'][:]
    # start_lat2 = ds.variables['start_lat2'][:]
    # start_lat3 = ds.variables['start_lat3'][:]

    # start_lon1 = ds.variables['start_lon1'][:]
    # start_lon2 = ds.variables['start_lon2'][:]
    # start_lon3 = ds.variables['start_lon3'][:]

    # # Closing Dataset
    # ds.close()

    # # Listing vertices
    # vertices = [(start_lon1, start_lat1), (start_lon2, start_lat2), (start_lon3, start_lat3)]

    # # Initialising area list
    # areas = []

    # # Multiprocessing process
    # process1 = multiprocessing.Process(target=tri_area, args=(range(len(start_lat1) // 1000), ))

    # process1.start()

    # process1.join()

    # print(f"AREA CALCULATIONS DONE, SORTING...")

    # print(max(areas))
    # max_areas = heapq.nlargest(10, areas)

    # print(max_areas)

    