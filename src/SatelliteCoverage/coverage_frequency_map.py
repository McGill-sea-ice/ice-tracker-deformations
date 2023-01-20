"""
Author: Lekima Yakuden
GitHub: LekiYak

--------------------------------------------------------------------------------
Tools for analysing raw data files
--------------------------------------------------------------------------------

This file contains functions for analysing raw data files' spatial and temporal coverage.
"""

import os
from time import strftime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import Normalize
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

from mpl_toolkits.axes_grid1 import make_axes_locatable
from config import *

# IGNORING WARNINGS, COMMENT IF YOU WANT TO SEE THEM
import warnings
warnings.filterwarnings("ignore")

# Visualises spatio-temporal coverage on a basemap
def visualise_coverage_histogram2d(xy, max_date, min_date, timestep):
    print('--- Plotting heat map ---')
    """
    Preamble
    """
    proj = ccrs.NorthPolarStereo(central_longitude=0)

    fig = plt.figure(figsize=(6.5, 5.5), )
    ax = fig.add_subplot(projection = proj, frameon=False)

    # Timestep calculation
    date_delta = max_date - min_date
    delta_hours = date_delta.total_seconds() // 3600
    length = delta_hours // int(timestep)

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

    xscale = math.floor(xscale / (1000 * float(resolution)))
    yscale = math.floor(yscale / (1000 * float(resolution)))

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
    if timestep != '0':
        ax.set_title(f'{tracker}, {min_date} to {max_date}, {timestep} \u00B1 {tolerance} hours, {resolution} km resolution')

    elif timestep == '0':
        ax.set_title(f'{tracker}, {min_date} to {max_date} encompassing all time intervals')

    # Saving figure as YYYYMMDD_YYYYMMDD_deltat_tolerance_resolution_'res'_tracker_freq.png
    prefix = tracker + '_' + min_date_str + '_' + max_date_str + '_dt' + timestep + '_tol' + tolerance + '_res' + str(int(resolution)) 

    # Set a directory to store figures
    figsPath =  output + '/' + '/figs/'

    # Create directory if it doesn't exist
    os.makedirs(figsPath, exist_ok=True)

    plt.savefig(figsPath + prefix + '_' + 'freq.png')

    print(f'Saved as {prefix}_freq.png')

    return hh[1], hh[2]

# Returns histogram of coverage representing the arctic ocean
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
    Date_options = config['Date_options']
    options = config['options']
    meta = config['Metadata']

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

# Plots timeseries of spatial coverage
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
    None -- Plot of % of area covered as a function of time.

    """
    print('--- Plotting coverage time series ---')

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

    # Plotting timeseries
    df.plot(x='start_date', y='percentage', kind='line')

    # Set a directory to store figures
    figsPath =  output + '/' + '/figs/'

    # Create directory if it doesn't exist
    os.makedirs(figsPath, exist_ok=True)

    # Title
    title = tracker + '_' + start_year + start_month + start_day + '_' +end_year + end_month + end_day + '_dt'+ timestep + '_tol' + tolerance + '_res' + resolution  + '_coverage_area_timeseries'+ '.png'

    # Saving figure
    plt.savefig(figsPath + title, bbox_inches='tight')

    print('Saved as ' + figsPath)

# Visualises coverage as a heatmap, split between user-set intervals
def interval_frequency_histogram2d(interval_list):
    """
    Plots a heatmap showing data availaibility (in % of time intervals covered) in a specified range of times.

    INPUTS:
    interval_list -- List of lists, with each sub-list (n-th) containing the paths to data files contained within
                     the n-th interval.

    OUTPUTS:
    Heatmap image (.png) in user-set output folder, with file name 'TRACKER_STARTDATE_to_ENDDATE_DELTAt+-TOLERANCE_
    hrs_RESOLUTION_km_INTERVAL_hrs"

    """
    print('--- Plotting interval frequency histogram ---')

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
        histogram[histogram>0.0] = 1.0
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
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],projection = proj)

    out_proj = pyproj.Proj(init='epsg:4326')
    in_proj = pyproj.Proj('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ', preserve_units=True)

    xx, yy = np.meshgrid(xbins, ybins)

    binslon,binslat = pyproj.transform(in_proj,out_proj,xx,yy)
    H[H==0.0] = np.nan
    H = H*100.0
    print(binslon.shape,binslat.shape,H.shape)

    cmap1 = mpl.colormaps['plasma']
    cmap1.set_bad('w')

    im = ax.pcolormesh(binslon,binslat,H/len(interval_list),
                           transform=ccrs.PlateCarree(),
                           vmin=0,vmax=100.0,
                           cmap=cmap1)

    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%",pad=0.2)
    clb = plt.colorbar(im)#,shrink=0.5
    """
    Terrain
    """
    # Upper (u) and lower (l) extents of x, y (metres)
    lxextent = -3100000
    uxextent = 2500000
    uyextent = 2500000
    lyextent = -1900000

    # Show lat/lon grid
    ax.gridlines()

    # Hide datapoints over land
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

    """
    Data
    """
    # Grid resolution calculations
#    xscale = uxextent - lxextent
#    yscale = uyextent - lyextent

#    xscale = math.floor(xscale / (1000 * int(resolution)))
#    yscale = math.floor(yscale / (1000 * int(resolution)))

#    length = len(interval_list)

    # Converting dates for title and file name purposes
    max_date_title = str(max_date.strftime('%Y')) + '-' + str(max_date.strftime('%m')) + '-' + str(max_date.strftime('%d'))
    min_date_title = str(min_date.strftime('%Y')) + '-' + str(min_date.strftime('%m')) + '-' + str(min_date.strftime('%d'))
    max_date_str = max_date.strftime("%Y%m%d")
    min_date_str = min_date.strftime("%Y%m%d")

    # Set a directory to store figures
    figsPath =  output + '/' + '/figs/'

    # Create directory if it doesn't exist
    os.makedirs(figsPath, exist_ok=True)

    # if/elif for title creation, for grammatical correctness
    if timestep != '0':
        ax.set_title(f'Percent coverage ({interval}h intervals), {tracker}, \n {min_date_title} - {max_date_title}, {timestep} \u00B1 {tolerance} h pairs')

    elif timestep == '0':
        ax.set_title(f'{tracker}, {min_date} to {max_date}, all timesteps, {resolution} km, {interval} hr intervals')

    # Saving figure as YYYYMMDD_YYYYMMDD_timestep_tolerance_resolution_'res'_tracker_freq.png
    prefix = tracker + '_'+ min_date_str + '_' + max_date_str + '_dt' + timestep + '_tol' + tolerance + '_res' + resolution + '_' + interval
    plt.savefig(figsPath + prefix + '_' + 'intervalfreq.png')

    print(f'Saved as {prefix}_intervalfreq.png')




if __name__ == '__main__':
    config = read_config()

    # Initializing more specific ConfigParser objects
    IO = config['IO']
    Date_options = config['Date_options']
    options = config['options']
    meta = config['Metadata']
    coverage_frequency = config['coverage_frequency']

    data_path = IO['data_folder']
    output = IO['output_folder']
    tracker = meta['ice_tracker']

    start_year = Date_options['start_year']
    start_month = Date_options['start_month']
    start_day = Date_options['start_day']

    end_year = Date_options['end_year']
    end_month = Date_options['end_month']
    end_day = Date_options['end_day']

    timestep = Date_options['timestep']
    tolerance = Date_options['tolerance']

    resolution = options['resolution']
    interval = options['interval']

    # Fetching filter information
    raw_list, max_date, min_date = filter_data(Date_options = Date_options, IO = IO, Metadata = meta)

    # Dividing data into intervals if the user desires
    if coverage_frequency['visualise_timeseries'] == 'True' or coverage_frequency['visualise_interval'] == True:
        interval_list, date_pairs = divide_intervals(raw_list, max_date, min_date, interval)

    # Compiling master dataframe
    df = compile_data(raw_list)

    # Converting points from lat/lon to EPSG 3413
    xy = convert_to_grid(df['lon'], df['lat'])

    # Plotting coverage heat map
    xbins, ybins = visualise_coverage_histogram2d(xy, max_date, min_date, timestep)

    if coverage_frequency['visualise_timeseries'] == 'True':

        # Plotting time series
        coverage_timeseries(interval_list, resolution, date_pairs)

    if coverage_frequency['visualise_interval'] == 'True':

        # Plotting interval heat map
        interval_frequency_histogram2d(interval_list)
