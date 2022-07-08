import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
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
    
    df -- Pandas dataframe with the following columns: Index  sLat  sLon {Dataframe}
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

        df = df.append(temp_df)

        # Updating counter
        i += 1
        print(f'{i} / {num_files}')

    return df

def count_duplicates(df):
    """
    Takes dataframe with sLat and sLon values and counts duplicates, returns dataframe
    with coordinates and frequency count.

    INPUT:
    df -- Pandas dataframe with sLat and sLong in two respective columns, representing data points

    OUTPUT:
    df -- Input dataframe with frequency column added
    
    """

    # Pivot df onto its side to count duplicate Lat & Lon entires
    df = df.pivot_table(columns=['sLat', 'sLon'], aggfunc='size')

    # Add size Multi index as column named 'frequency'
    df = df.reset_index()  
    
    df.rename(columns = {0:'frequency'}, inplace = True)

    return df

def visualise_coverage_old(df):

    print("Plotting frequency map")

    # Set the matplotlib projection and transform
    #proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
    proj=ccrs.NorthPolarStereo(central_longitude=0)
    trans = ccrs.Geodetic()

    # Initialize figure
    fig = plt.figure()

    # Add projection
    ax = fig.add_subplot(111, projection = proj)
    #ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())
    ax.set_extent((-3100000, 2500000, -1900000, 2500000), ccrs.NorthPolarStereo())
    
    # Setting coordinates
    ax.scatter(x=df['sLon'], y=df['sLat'], transform= trans, s=0.1)

    # Set grid
    #ax.gridlines()

    # Hide datapoints over land
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

    # Set title
    ax.set_title('Heatmap of datapoint frequency over the 2019-2020 Winter (S1)')
    
    # Save figure
    plt.savefig('TEST_ALL2.png')

    return

def visualise_coverage_histogram2d(xy, max_date, min_date):
    print('Plotting heat map')
    """
    Preamble
    """
    # Reading config
    config = read_config()

    IO = config['IO']
    options = config['options']
    meta = config['meta']

    sYear = options['start_year']
    sMonth = options['start_month']
    sDay = options['start_day']
    eYear = options['end_year']
    eMonth = options['end_month']
    eDay = options['end_day']
    delta_t = options['delta_t']
    tolerance = options['tolerance']
    resolution = float(options['resolution'])
    strres = options['resolution']
    output_path = IO['output_folder']
    tracker = meta['ice_tracker']

    proj=ccrs.NorthPolarStereo(central_longitude=0)

    fig = plt.figure(figsize=(6.5, 5.5), )
    ax = fig.add_subplot(projection = proj, frameon=False)

    """
    Terrain
    """
    # Left, right (x) and up, down (y) extents of figure (metres)
    lxextent = -3100000
    rxextent = 2500000
    uyextent = 2500000
    dyextent = -1900000
    
    # Set grid
    ax.gridlines(draw_labels=True)

    # Hide datapoints over land
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')
    
    """
    Data
    """
    # Resolution calculations
    xscale = rxextent - lxextent
    yscale = uyextent - dyextent

    xscale = math.floor(xscale / (1000 * resolution))
    yscale = math.floor(yscale / (1000 * resolution))

    # Extracting x and y coordinates (Numpy arrays)
    xi, yj = xy

    # Plotting histogram (cmin=1 to make values = 0 transparent)
    hh = ax.hist2d(xi, yj, bins=(xscale, yscale), cmap='plasma', cmin=1)

    # Adding colourbar (hh[3] normalizes the colourbar)
    fig.colorbar(hh[3], ax=ax)

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
        ax.set_title(f'{min_date} to {max_date}, {delta_t} \u00B1 {tolerance} hours, {resolution} km resolution')
    
    elif delta_t == '0':
        ax.set_title(f'{min_date} to {max_date} encompassing all time intervals')

    # Saving figure as YYYYMMDD_YYYYMMDD_deltat_tolerance_resolution_'res'_tracker_freq.png
    prefix = min_date_str + '_' + max_date_str + '_' + delta_t + '_' + tolerance + '_' + 'res' + strres + '_' + tracker
    plt.savefig(prefix + '_' + 'freq.png')

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






config = read_config()

# Initializing more specific ConfigParser objects
IO = config['IO']
options = config['options']

data_path = IO['data_folder']

sYear = options['start_year']
sMonth = options['start_month']
sDay = options['start_day']

eYear = options['end_year']
eMonth = options['end_month']
eDay = options['end_day']

delta_t = options['delta_t']
tolerance = options['tolerance']

resolution = options['resolution']

# Fetching filter information
raw_list = filter_data(sYear, sMonth, sDay, eYear, eMonth, eDay, delta_t, tolerance, data_path)['raw_paths']
max_date = filter_data(sYear, sMonth, sDay, eYear, eMonth, eDay, delta_t, tolerance, data_path)['max_date']
min_date = filter_data(sYear, sMonth, sDay, eYear, eMonth, eDay, delta_t, tolerance, data_path)['min_date']

# Compiling master dataframe
df = compile_data(raw_list)

# Converting points from lat/lon to EPSG 3413
xy = convert_to_grid(df['sLon'], df['sLat'])

# Plotting heat map
visualise_coverage_histogram2d(xy, max_date, min_date)

