import os
from time import strftime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import utils_load_grid as grid
import math

from config import *

def compile_data(raw_paths):

    """ 
    Compiles all data points in the list of data files (raw_paths) into a dataframe (df)

    Returns the dataframe with all datapoints' starting lat and lon in columns.

    INPUTS:
    raw_paths -- List of paths to data files {List}

    OUTPUTS: 
    df -- Pandas dataframe with the following columns: Index  lat  lon {Dataframe}
    """

    df = pd.DataFrame()

    # Initialising progression counters
    num_files = len(raw_paths)
    i = 0

    # Appending each file's datapoints to the dataframe
    for filepath in raw_paths:
        
        # Initialize temporary dataframe
        temp_df = pd.DataFrame()

        # Reading datapoints into temporary dataframe
        temp_df = pd.read_csv(filepath, sep='\s\s+', engine='python', usecols = ['sLat','sLon'])

        temp_df.rename(columns = {'sLat':'lat', 'sLon':'lon'}, inplace = True)

        df = df.append(temp_df)

        # Updating counter
        i += 1
        print(f'{i} / {num_files}')

    return df

def visualise_coverage_histogram2d(xy, max_date, min_date, delta_t):
    print('Plotting heat map')
    """
    Preamble
    """
    # Reading config
    config = read_config()

    IO = config['IO']
    options = config['options']
    meta = config['meta']

    delta_t = options['delta_t']
    tolerance = options['tolerance']
    resolution = float(options['resolution'])
    strres = options['resolution']
    output_path = IO['output_folder']
    tracker = meta['ice_tracker']

    proj = ccrs.NorthPolarStereo(central_longitude=0)

    fig = plt.figure(figsize=(6.5, 5.5), )
    ax = fig.add_subplot(projection = proj, frameon=False)

    # Timestep calculation
    date_delta = max_date - min_date
    delta_hours = date_delta.total_seconds() // 3600
    length = delta_hours // int(delta_t)

    """
    Terrain
    """
    # Upper (u) and lower (l) extents of x, y (metres)
    lxextent = -3100000
    uxextent = 2500000
    uyextent = 2500000
    lyextent = -1900000
    
    # Show lat/lon grid
    ax.gridlines(draw_labels=True)

    # Hide datapoints over land
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
    
    """
    Data
    """
    # Grid resolution calculations
    xscale = uxextent - lxextent
    yscale = uyextent - lyextent

    xscale = math.floor(xscale / (1000 * resolution))
    yscale = math.floor(yscale / (1000 * resolution))

    # Extracting x and y coordinates of datapoints (Numpy arrays)
    xi, yj = xy
    
    # Plotting histogram (cmin=1 to make values = 0 transparent)
    norm = plt.Normalize(0, int(length))
    hh = ax.hist2d(xi, yj, bins=(xscale, yscale), cmap='plasma', cmin=1, norm=norm)

    # Adding colourbar (hh[3] normalizes the colourbar)
    cbar = fig.colorbar(hh[3], ax=ax)
    cbar.set_label(f'% of total period tile has data')
    cbar.set_ticklabels(np.arange(0, 110, 10))
    cbar.set_ticks(np.linspace(0, int(length), num=11))

    """
    Extraneous
    """
    # Max and min dates for title and file name (_str)
    max_date_str = max_date.strftime("%Y%m%d")
    min_date_str = min_date.strftime("%Y%m%d")
    max_date = max_date.strftime("%x")
    min_date = min_date.strftime("%x")

    # if/elif for title creation, for grammatical correctness
    if delta_t != '0':
        ax.set_title(f'{tracker}, {min_date} to {max_date}, {delta_t} \u00B1 {tolerance} hours, {resolution} km resolution')
    
    elif delta_t == '0':
        ax.set_title(f'{tracker}, {min_date} to {max_date} encompassing all time intervals')

    # Saving figure as YYYYMMDD_YYYYMMDD_deltat_tolerance_resolution_'res'_tracker_freq.png
    prefix = min_date_str + '_' + max_date_str + '_' + delta_t + '_' + tolerance + '_' + 'res' + str(int(resolution)) + '_' + tracker
    plt.savefig(output + '/' + prefix + '_' + 'freq.png')

    print(f'Saved as {prefix}_freq.png')

    return hh[1], hh[2]

def coverage_histogram2d(xy, xbins, ybins):
    """
    Returns a 2D numpy array representing grid cells with or without data
    (1 or 0).

    INPUTS:
    xy -- Array of tuples containing x and y coordinates of data {numpy array}

    OUTPUTS:
    H -- 2D numpy array representing polar stereographic grid with 1s and 0s,
         representing the presence of lack of data. {numpy array}
    """

    """
    Preamble
    """
    # Reading config & preamble
    config = read_config()

    IO = config['IO']
    options = config['options']
    meta = config['meta']

    resolution = float(options['resolution'])

    proj = ccrs.NorthPolarStereo(central_longitude=0)

    fig = plt.figure(figsize=(6.5, 5.5), )
    ax = fig.add_subplot(projection = proj, frameon=False)

    """
    Data
    """
    # Upper (u) and lower (l) extents of histogram grid (metres)
    lxextent = -3100000
    uxextent = 2500000
    uyextent = 2500000
    lyextent = -1900000
 
    # Grid resolution calculations (x,yscale final values in *resolution* km)
    xscale = uxextent - lxextent
    yscale = uyextent - lyextent

    xscale = math.floor(xscale / (1000 * resolution))
    yscale = math.floor(yscale / (1000 * resolution))

    # Extracting x and y coordinates of datapoints (numpy arrays)
    xi, yj = xy

    # Plotting histogram (H) and converting bin values to 0 or 1 range=[[lxextent,uxextent], [lyextent,uyextent]]
    H, xbins, ybins = np.histogram2d(xi, yj, bins=(xbins, ybins))
    H[H>1] = 1

    return H

def convert_to_grid(lon, lat):
    """
    WARNING: INPUT IS LON, LAT (X, Y), AND THE OUTPUT IS ALSO X, Y

    Takes in a point in EPSG 4326 (Lat/Lon) and transforms said point to
    a polar stereographic projection (EPSG 3413).

    lon, lat {float} -> x, y {float}

    INPUTS:
    lat, lon -- Float values representing a point in WGS 84

    OUTPUTS:
    x, y -- Float values representing the transformed point in EPSG 3413
    """

    # Input projection (lat/lon)
    in_proj = pyproj.Proj(init='epsg:4326')

    # Output projection (EPSG 3413, Polar Stereographic)
    out_proj = pyproj.Proj('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ', preserve_units=True)

    # Transform function
    x, y = np.array(pyproj.transform(in_proj, out_proj, lon, lat))
    return x, y

def coverage_timeseries(interval_list, resolution, date_pairs):
    """
    Plots a time series of the area coverage (in % of the Arctic ocean) for a given list of lists containing
    data file paths [interval_list], where each list of files defines a user-set interval (i.e. interval of 72hrs)

    INPUTS:
    interval_list -- List of lists, each sublist containing data file paths and each sublist (index n) representing 
                     the data files which share temporal overlap with the n th interval. {list}

    resolution -- Resolution of grid to be used, in km. Read from config. {str}

    date_pairs -- List of tuples containing the start and end dates of the n th interval, whose data contents can be
                  found at the same index in *interval_list* {list}

    OUTPUTS:
    Plot of % of area covered as a function of time.

    """

    # Setting constants (Adjusting ocean area to units of histogram grid)
    arctic_ocean_area = 15558000 # Square kilometres
    arctic_ocean_area = arctic_ocean_area / int(resolution) ** 2

    # Initialising dataframe to store interval data
    df = pd.DataFrame(columns=['percentage', 'start_date', 'end_date'])

    # Iterating over each interval
    for i in range(len(interval_list)):
        # Loads data and converts to x/y for each interval
        interval_df = compile_data(interval_list[i])

            # Skips empty lists
        try:
            xy = convert_to_grid(interval_df['lon'], interval_df['lat'])
        except KeyError:
            continue

        # Generates histogram (2d np array)
        histogram = coverage_histogram2d(xy, xbins, ybins)

        # Computing area of arctic ocean covered
        covered_area = len((np.flatnonzero(histogram)))
        covered_percentage = (covered_area / arctic_ocean_area) * 100

        # Extracting dates
        start_date = date_pairs[i][0]
        end_date = date_pairs[i][1]

        # Appending to main dataframe
        df.loc[len(df.index)] = [covered_percentage, start_date, end_date]

    df.plot(x='start_date', y='percentage', kind='line')

    plt.savefig('coverage_areaRSCMS1202011102021060172hrs.png')

    print(df.sort_values(by=['percentage']))

    return 

def interval_frequency_histogram2d(interval_list):
    """
    
    """
    # Initializing empty numpy array (2D histogram)
    H = np.array([])

    # Iterating over each interval
    for i in range(len(interval_list)):
        # Loads data and converts to x/y for each interval
        interval_df = compile_data(interval_list[i])

            # Skips empty lists
        try:
            xy = convert_to_grid(interval_df['lon'], interval_df['lat'])
        except KeyError:
            xy = (0,0)
            continue

        # Generates histogram (2D numpy array)
        histogram = coverage_histogram2d(xy, xbins, ybins)

        # Changing size of total histogram (only on first run)
        if i == 0:
            H.resize(histogram.shape)

        # Adding interval-specific histogram to total histogram
        H = H + histogram

    # Transposing histogram for plotting
    H = H.T

    """
    Land and projection
    """

    proj = ccrs.NorthPolarStereo(central_longitude=0)

    fig = plt.figure(figsize=(6.5, 5.5), )
    ax = fig.add_subplot(projection = proj, frameon=False)

    """
    Terrain
    """
    # Upper (u) and lower (l) extents of x, y (metres)
    lxextent = -3100000
    uxextent = 2500000
    uyextent = 2500000
    lyextent = -1900000
    
    # Show lat/lon grid
    ax.gridlines(draw_labels=True)

    # Hide datapoints over land
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
    
    """
    Data
    """
    # Grid resolution calculations
    xscale = uxextent - lxextent
    yscale = uyextent - lyextent

    xscale = math.floor(xscale / (1000 * int(resolution)))
    yscale = math.floor(yscale / (1000 * int(resolution)))

    length = len(interval_list)

    # Plotting histogram using axis.contourf function
    norm = plt.Normalize(0, length) # Normalizing colours (values)
    levels = np.linspace(1, length, 11) # Dividing contours into 10 levels (10%, 20% etc)
    h = ax.contourf(np.delete(xbins, len(xbins) - 1), np.delete(ybins, len(ybins) - 1), H, cmap='plasma', norm=norm, levels=levels)

    # Colourbar
    cbar = fig.colorbar(h, ax=ax)
    cbar.set_label(f'% of total period tile has data')
    cbar.set_ticklabels(np.arange(0, 110, 10))

    # Saving the file
    plt.axis('scaled')

    max_date_title = str(max_date.strftime('%Y')) + '-' + str(max_date.strftime('%m')) + '-' + str(max_date.strftime('%d'))
    min_date_title = str(min_date.strftime('%Y')) + '-' + str(min_date.strftime('%m')) + '-' + str(min_date.strftime('%d'))
    max_date_str = max_date.strftime("%Y%m%d")
    min_date_str = min_date.strftime("%Y%m%d") 
  

    # if/elif for title creation, for grammatical correctness
    if delta_t != '0':
        ax.set_title(f'{tracker}, {min_date_title} to {max_date_title}, {delta_t} \u00B1 {tolerance} hrs, {resolution} km, {interval} hr intervals')
    
    elif delta_t == '0':
        ax.set_title(f'{tracker}, {min_date} to {max_date}, all timesteps, {resolution} km, {interval} hr intervals')

    # Saving figure as YYYYMMDD_YYYYMMDD_deltat_tolerance_resolution_'res'_tracker_freq.png
    prefix = min_date_str + '_' + max_date_str + '_' + delta_t + '_' + tolerance + '_' + 'res' + resolution + '_' + tracker + '_' + interval
    plt.savefig(output + '/' + prefix + '_' + 'intervalfreq.png')

    print(f'Saved as {prefix}_freq.png')
    
    return h

config = read_config()

# Initializing more specific ConfigParser objects
IO = config['IO']
options = config['options']
meta = config['meta']

data_path = IO['data_folder']
output = IO['output_folder']
tracker = meta['ice_tracker']

sYear = options['start_year']
sMonth = options['start_month']
sDay = options['start_day']

eYear = options['end_year']
eMonth = options['end_month']
eDay = options['end_day']

delta_t = options['delta_t']
tolerance = options['tolerance']

resolution = options['resolution']

interval = options['interval']

# Fetching filter information
raw_list = filter_data(sYear, sMonth, sDay, eYear, eMonth, eDay, delta_t, tolerance, data_path)['raw_paths']
max_date = filter_data(sYear, sMonth, sDay, eYear, eMonth, eDay, delta_t, tolerance, data_path)['max_date']
min_date = filter_data(sYear, sMonth, sDay, eYear, eMonth, eDay, delta_t, tolerance, data_path)['min_date']

interval_list = divide_intervals(raw_list, max_date, min_date, interval)['interval_list']
date_pairs = divide_intervals(raw_list, max_date, min_date, interval)['date_pairs']

# Compiling master dataframe
df = compile_data(raw_list)

# Converting points from lat/lon to EPSG 3413
xy = convert_to_grid(df['lon'], df['lat'])

# Plotting heat map
xbins, ybins = visualise_coverage_histogram2d(xy, max_date, min_date, delta_t)

#xbins = np.delete(xbins, len(xbins)-3)
#ybins = np.delete(ybins, len(ybins)-3)

# Plotting time series
#coverage_timeseries(interval_list, resolution, date_pairs)

# Plotting interval heat map
print(interval_frequency_histogram2d(interval_list))