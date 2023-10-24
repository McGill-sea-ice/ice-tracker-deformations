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
from matplotlib import dates as mdates
import cartopy.crs as ccrs
import matplotlib as mpl
from shapely.geometry import Point
from shapely.prepared import prep
import pyproj
import cartopy.feature as cfeature

# Loading from other files
from SatelliteCoverage.config import *
from SatelliteCoverage.utils import *

# Generate map x/y bins that will be used to compute frequency at each cell on map
def get_map_bins(config=None):

    resolution = float(config['options']['resolution'])

    # Upper (u) and lower (l) extents of map_x, map_y (metres)
    lxextent = -4400000
    uxextent =  2600000
    uyextent =  4000000
    lyextent = -2600000

    # Make bins vectors
    dxi = (1000*float(resolution))
    dyj = (1000*float(resolution))

    xbins_out = np.arange(lxextent,uxextent+dxi,dxi)
    ybins_out = np.arange(lyextent,uyextent+dyj,dyj)

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

    percentage = stb(config['options']['percentage'])
    ref_lat = int(config['options']['ref_lat'])
    resolution = config['options']['resolution']
    interval = config['options']['interval']

    # Default to plotting area if min_lat is lower than 50N.
    if ref_lat < 50 :
        ref_lat = 0
        percentage = False

    # Get land points and map cells center lat/lon
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m')
    land_polygons = list(land_10m.geometries())

    dxi = (1000*float(resolution))
    dyj = (1000*float(resolution))
    xx_center, yy_center = np.meshgrid(xbins_map[0:-1]+0.5*dxi, ybins_map[0:-1]+0.5*dyj)
    cell_center_lat,cell_center_lon = convert_from_grid(xx_center, yy_center)

    # Check which map point is land
    points = [Point(point) for point in zip(cell_center_lon.ravel(), cell_center_lat.ravel())]
    land_polygons_prep = [prep(land_polygon) for land_polygon in land_polygons]
    land_lon = []
    land_lat = []
    for land_polygon in land_polygons_prep:
        land_lon.extend([point.coords[0][0] for point in filter(land_polygon.covers, points)])
        land_lat.extend([point.coords[0][1] for point in filter(land_polygon.covers, points)])
    land_lat = np.array(land_lat)
    land_lon = np.array(land_lon)

    # Compute reference area for percentage
    n_pts_land = len(land_lat[land_lat >= ref_lat])
    n_pts_tot = len(cell_center_lat[cell_center_lat >= ref_lat])
    ref_area = (n_pts_tot-n_pts_land)*(int(resolution)** 2)

    # Initialising dataframe to store interval data
    df = pd.DataFrame(columns=['start_date', 'end_date', 'area','percentage'])

    # Iterating over each interval
    for i in tqdm(range(len(interval_list)), position=0, leave=True):
        # Loads data and converts to x/y for each interval
        interval_df = compile_data(raw_paths=interval_list[i])

        # Extracting dates
        start_date = date_pairs[i][0]
        end_date = date_pairs[i][1]

        # Skips empty lists
        try:
            xy = convert_to_grid(interval_df['lon'], interval_df['lat'])
        except KeyError:
            df.loc[len(df.index)] = [start_date, end_date, 10**30, 10**30]
            continue

        # Generates histogram (2d np array)
        histogram = coverage_histogram2d(xy, xbins_map, ybins_map)
        histogram = histogram.T
        histogram[histogram > 0] = 1
        histogram[histogram < 1] = np.nan

        # Computing area of arctic ocean covered
        covered_area = (np.nansum(histogram[cell_center_lat >= ref_lat]))*(int(resolution)** 2)
        covered_percentage = (covered_area / ref_area) * 100

        # Appending to main dataframe
        df.loc[len(df.index)] = [start_date, end_date, covered_area, covered_percentage]

    # Place nans where there is no coverage
    df.loc[df['area'] >= 10**30,'area'] = np.nan
    df.loc[df['percentage'] >= 10**30, 'percentage'] = np.nan
    # Plotting timeseries
    f, ax = plt.subplots()
    # ax.set_xlabel('Date')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    f.autofmt_xdate()
    if percentage:
        ax.plot(df['start_date'],df['percentage'],'-')
        ax.set_ylabel('Coverage ( % of ocean points above ' + str(ref_lat) + '$^{\circ}$N)')
    else:
        ax.plot(df['start_date'],df['area']/(1e6),'-')
        ax.set_ylabel('Coverage ' + ('above '+ str(ref_lat) + '$^{\circ}$N ')*(ref_lat > 0) + '(10$^{6}~$km$^2$)')

    output_folder = config['IO']['output_folder']
    exp = config['IO']['exp']

    # Set a directory to store figures
    figsPath =  output_folder + '/' + exp + '/figs/'

    # Create directory if it doesn't exist
    os.makedirs(figsPath, exist_ok=True)

    # prefix
    prefix = get_prefix(config=config)

    # figname
    figname = prefix + '_res' + resolution  + '_int' + interval + (ref_lat > 0)*('_reflat' + str(ref_lat)) + '_coverage'+ percentage*'_percent' +'_area_timeseries'

    # Saving figure
    print('Saving coverage timeserie figure at ' + figsPath + figname + '.png')
    plt.savefig(figsPath + figname + '.png', bbox_inches='tight', dpi=600)

    # save the time serie in pickle format
    print('Saving coverage timeserie data at ' + figsPath + figname + '.pkl')
    df.to_pickle(figsPath + figname + '.pkl')

    return None


# Visualises coverage as a heatmap, split between user-set intervals
def interval_frequency_histogram2d(interval_list, date_pairs, xbins_map, ybins_map, config=None):
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

    # List of summer months (Won't include data from these months)
    summer_months = [6, 7, 8, 9, 10]
    # Iterating over each interval
    for i in tqdm(range(len(interval_list)), position=0, leave=True):
        if date_pairs[i][0].month not in summer_months:
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
            if H.shape == (0,):
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
    n_int = int(np.sum([date_pairs[i][0].month not in summer_months for i in range(len(date_pairs))]))

    """
    Plot Data
    """
    cmap1 = mpl.colormaps['plasma']
    cmap1.set_bad('w')

    new_coords = ax.projection.transform_points(trans, binslon, binslat)

    im = ax.pcolormesh(new_coords[:,:,0], new_coords[:,:,1], H/n_int, vmin=0, vmax=100.0, cmap=cmap1)
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

    fig_name = figsPath + prefix + '_res' + resolution  + '_int' + interval + '_coverage_area_map.png'
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

        # Define map bins for coverage analysis
        xbins, ybins = get_map_bins(config=config)

        if viz_ts:
            # Plotting time series of coverage in % of total ocean area above
            # a given 'ref_lat', or of covered area in km^2 above a 'ref_lat'
            # Note: if 'ref_lat' is equal to zero, it takes all available data.
            coverage_timeseries(interval_list, date_pairs, xbins, ybins, config=config)

        if viz_cf:
            # Plotting coverage heat map in % of intervals with data
            interval_frequency_histogram2d(interval_list, date_pairs, xbins, ybins, config=config)


    # Display the run time
    print("--- %s seconds ---" % (time.time() - start_time))
