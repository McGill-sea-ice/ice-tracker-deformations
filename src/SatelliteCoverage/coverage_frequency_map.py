"""
Author: Lekima Yakuden
GitHub: LekiYak

--------------------------------------------------------------------------------
Tools for analysing raw data files
--------------------------------------------------------------------------------

This file contains functions for analysing raw data files' spatial and temporal coverage.
"""

# Loading from default packages
import os
import sys
parent = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,parent)
import time
import math
from tqdm import tqdm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib as mpl
import cartopy.feature as cfeature

# Loading from other files
from SatelliteCoverage.config import *
from SatelliteCoverage.utils import *

# Generate map x/y bins that will be used to compute frequency at each cell on map
def get_map_bins(xy, config=None):

    resolution = float(config['options']['resolution'])

    # Upper (u) and lower (l) extents of map_x, map_y (metres)
    lxextent = -3100000
    uxextent = 2500000
    uyextent = 2500000
    lyextent = -1900000

    # Grid resolution calculations
    xscale = uxextent - lxextent
    yscale = uyextent - lyextent

    xscale = math.floor(xscale / (1000 * float(resolution)))
    yscale = math.floor(yscale / (1000 * float(resolution)))

    # Extracting x and y coordinates of datapoints (Numpy arrays)
    xi, yj = xy

    # Make bins vectors
    dxi = (np.max(xi)-np.min(xi)) / xscale
    dyj = (np.max(yj)-np.min(yj)) / yscale

    xbins_out = np.arange(np.min(xi),np.max(xi)+dxi,dxi)
    ybins_out = np.arange(np.min(yj),np.max(yj)+dyj,dyj)

    return xbins_out, ybins_out

# Returns histogram of coverage representing the arctic ocean
def coverage_histogram2d(xy, xbins_map, ybins_map):
    """
    Returns a 2D numpy array representing grid cells with or without data
    (1 or 0).

    INPUTS:
    xy -- Array of tuples containing x and y coordinates of data {numpy array}
    xbins_map, ybins_map -- Map bins for frequency histogram

    OUTPUTS:
    H -- 2D numpy array representing polar stereographic grid with 1s and 0s,
         representing the presence of lack of data. {numpy array}
    """

    # Extracting x and y coordinates of datapoints (numpy arrays)
    xi, yj = xy

    # Plotting histogram (H) and converting bin values to 0 or 1 range=[[lxextent,uxextent], [lyextent,uyextent]]
    H, _, _ = np.histogram2d(xi, yj, bins=(xbins_map, ybins_map))
    H[H>1] = 1

    return H


# Plots timeseries of spatial coverage
def coverage_timeseries(interval_list, date_pairs, xbins_map, ybins_map, config=None):
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

    resolution = config['options']['resolution']
    interval = config['options']['interval']

    # Setting constants (Adjusting ocean area to units of histogram grid)
    arctic_ocean_area = 15558000 # Square kilometres
    arctic_ocean_area = arctic_ocean_area / int(resolution) ** 2

    # Initialising dataframe to store interval data
    df = pd.DataFrame(columns=['percentage', 'start_date', 'end_date'])

    # Iterating over each interval
    for i in tqdm(range(len(interval_list)), position=0, leave=True):
        # Loads data and converts to x/y for each interval
        interval_df = compile_data(raw_paths=interval_list[i])

        # Skips empty lists
        try:
            xy = convert_to_grid(interval_df['lon'], interval_df['lat'])
        except KeyError:
            continue

        # Generates histogram (2d np array)
        histogram = coverage_histogram2d(xy, xbins_map, ybins_map)

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

    output_folder = config['IO']['output_folder']
    exp = config['IO']['exp']
    icetracker = config['Metadata']['icetracker']
    Date_options = config['Date_options']
    # Converting dates for title and file name purposes
    start_year  = str(Date_options['start_year'])
    start_month = str(Date_options['start_month'])
    start_day   = str(Date_options['start_day'])
    end_year    = str(Date_options['end_year'])
    end_month   = str(Date_options['end_month'])
    end_day     = str(Date_options['end_day'])
    timestep    = str(Date_options['timestep'])
    tolerance   = str(Date_options['tolerance'])

    # Set a directory to store figures
    figsPath =  output_folder + '/' + exp + '/figs/'

    # Create directory if it doesn't exist
    os.makedirs(figsPath, exist_ok=True)

    # prefix
    prefix = get_prefix(config=config)

    # figname
    figname = prefix + '_res' + resolution  + '_int' + interval + '_coverage_area_timeseries'

    # Saving figure
    print('Saving coverage timeserie figure at ' + figsPath + figname + '.png')
    plt.savefig(figsPath + figname + '.png', bbox_inches='tight', dpi=600)

    # save the time serie in pickle format
    print('Saving coverage timeserie data at ' + figsPath + figname + '.pkl')
    df.to_pickle(figsPath + figname + '.pkl')


# Visualises coverage as a heatmap, split between user-set intervals
def interval_frequency_histogram2d(interval_list, xbins_map, ybins_map, config=None):
    """
    Plots a heatmap showing data availaibility (in % of time intervals covered) in a specified range of times.

    INPUTS:
    interval_list -- List of lists, with each sub-list (n-th) containing the paths to data files contained within
                     the n-th interval.
    Date_options --

    OUTPUTS:
    Heatmap image (.png) in user-set output folder, with file name 'STARTDATE_to_ENDDATE_DELTAt+-TOLERANCE_
    hrs_RESOLUTION_km_INTERVAL_hrs"

    """
    print('--- Plotting interval frequency histogram ---')

    Date_options = config['Date_options']

    # Initializing empty numpy array (2D histogram)
    H = np.array([])

    # Iterating over each interval
    for i in tqdm(range(len(interval_list)), position=0, leave=True):
        # Loads data and converts to x/y for each interval
        interval_df = compile_data(interval_list[i])

        # Skips empty lists
        try:
            xy = convert_to_grid(interval_df['lon'], interval_df['lat'])
        except KeyError:
            xy = (0,0)
            continue

        # Generates histogram (2D numpy array)
        histogram = coverage_histogram2d(xy, xbins_map, ybins_map)
        histogram[histogram > 0.0] = 1.0
        # Changing size of total histogram (only on first run)
        if i == 0:
            H.resize(histogram.shape)


        # Adding interval-specific histogram to total histogram
        H = H + histogram

    # Transposing histogram for plotting
    H = H.T
    H[H == 0] = np.nan

    """
    Land and projection
    """
    proj = ccrs.NorthPolarStereo(central_longitude=0)
    trans = ccrs.Geodetic()

    fig = plt.figure(figsize=(6.5, 5.5), )
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],projection = proj)

    # out_proj = pyproj.CRS('epsg:4326')
    # in_proj = pyproj.CRS('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ', preserve_units=True)

    xx, yy = np.meshgrid(xbins_map, ybins_map)

    binslat, binslon = convert_from_grid(xx, yy)

    H[H==0.0] = np.nan
    H = H*100.0

    """
    Plot Data
    """
    cmap1 = mpl.colormaps['plasma']
    cmap1.set_bad('w')

    new_coords = ax.projection.transform_points(trans, binslon, binslat)

    im = ax.pcolormesh(new_coords[:,:,0], new_coords[:,:,1], H/len(interval_list), vmin=0, vmax=100.0, cmap=cmap1)
    # im = ax.pcolormesh(xx, yy, H/len(interval_list), vmin=0, vmax=100.0, cmap=cmap1)

    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%",pad=0.2)
    clb = plt.colorbar(im)#,shrink=0.5
    clb.set_label('% of total period tile has data')
    clb.set_ticks(np.arange(0, 110, 10))
    clb.set_ticklabels(np.arange(0, 110, 10))
    # plt.axis('scaled')

    # Show lat/lon grid
    ax.gridlines()

    # Hide datapoints over land
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

    # Converting dates for title and file name purposes
    start_year  = str(Date_options['start_year'])
    start_month = str(Date_options['start_month'])
    start_day   = str(Date_options['start_day'])
    end_year    = str(Date_options['end_year'])
    end_month   = str(Date_options['end_month'])
    end_day     = str(Date_options['end_day'])
    timestep    = str(Date_options['timestep'])
    tolerance   = str(Date_options['tolerance'])
    resolution  = str(config['options']['resolution'])
    interval    = str(config['options']['interval'])
    icetracker  = str(config['Metadata']['icetracker'])

    # Concatenate start and end dates
    sDate = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
    eDate = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

    eDate_title = end_year + '-' + end_month + '-' + end_day
    sDate_title = start_year + '-' + start_month + '-' + start_day
    eDate_str = eDate.strftime("%Y%m%d")
    sDate_str = sDate.strftime("%Y%m%d")

    output_folder = config['IO']['output_folder']
    exp = config['IO']['exp']
    # Set a directory to store figures
    figsPath =  output_folder + '/' + exp + '/figs/'

    # Create directory if it doesn't exist
    os.makedirs(figsPath, exist_ok=True)

    # if/elif for title creation, for grammatical correctness
    if timestep != '0':
        ax.set_title(f'Percent coverage ({interval}h intervals), {icetracker}, \n {sDate_title} - {eDate_title}, {timestep} \u00B1 {tolerance} h pairs')

    elif timestep == '0':
        ax.set_title(f'{tracker}, {sDate_title} to {eDate_title}, all timesteps, {resolution} km, {interval} hr intervals')

    # Saving figure as YYYYMMDD_YYYYMMDD_timestep_tolerance_resolution_'res'_tracker_freq.png
    prefix = get_prefix(config=config)

    fig_name = figsPath + prefix + '_res' + resolution  + '_int' + interval + 'coverage_area_map.png'
    print('Saving coverage 2D histogram figure at ' + fig_name)
    plt.savefig(fig_name, dpi=600, bbox_inches='tight')


if __name__ == '__main__':

    # Retrieve the starting time
    start_time = time.time()

    config = read_config()

    # load options from config
    viz_ts = config['coverage_frequency']['visualise_timeseries']
    viz_cf = config['coverage_frequency']['visualise_interval']

    # Fetching filter information
    raw_list = filter_data(config=config)
    config['raw_list'] = raw_list

    # Dividing data into intervals if the user desires
    if viz_ts or viz_cf:
        interval_list, date_pairs = divide_intervals(config=config)

        # Compiling master dataframe
        df = compile_data(raw_paths=config['raw_list'])

        # Converting points from lat/lon to EPSG 3413
        xy = convert_to_grid(df['lon'], df['lat'])

        # Plotting coverage heat map
        xbins, ybins = get_map_bins(xy, config=config)

        if viz_ts:
            # Plotting time series of coverage in % of total Arctic Ocean area
            coverage_timeseries(interval_list, date_pairs, xbins, ybins, config=config)

        if viz_cf:
            # Plotting coverage heat map in % of intervals with data
            interval_frequency_histogram2d(interval_list, xbins, ybins, config=config)

    # Display the run time
    print("--- %s seconds ---" % (time.time() - start_time))