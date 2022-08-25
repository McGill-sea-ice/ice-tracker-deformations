"""
Author: Lekima Yakuden 
GitHub: LekiYak

--------------------------------------------------------------------------------
Configuration file for data / netCDF analysis tools
--------------------------------------------------------------------------------

This file contains functions for loading and processing user options, raw data, and netCDF files.
"""

import configparser
import os
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset
import haversine as hs
import pandas as pd
import pyproj

# Loading config file
def read_config():
    """
    Imports user options from options.ini

    INPUTS:
    None

    OUTPUTS:
    ConfigParser object with user preferences {ConfigParser}
    """
    # Path of current directory (where options.ini should be located)
    cwd = os.path.dirname(os.path.realpath(__file__))

    # Reading options.ini
    config = configparser.ConfigParser()
    config.read(cwd + '/options.ini')

    # Return a ConfigParser object
    return config

"""
Raw Data 
"""
# Compiles raw data into pandas dataframe
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

# Converts lat/lon to the north pole stereographic projection (EPSG 3413)
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

# Divides raw data into intervals specified by the user
def divide_intervals(raw_paths, max_date, min_date, interval):
    """
    Divides delta-t filtered data into chunks (intervals) of *interval* hours for
    processing.

    INPUTS:
    raw_paths -- List of data file paths with the desired delta t {list}
    max_date -- Upper limit of selected date range {datetime object}
    min_date -- Lower limit of selected date range {datetime object}
    interval -- Desired interval length in hours {str}

    OUTPUTS:
    interval_list -- List of n lists, where n is the number of full intervals which
                     fit in the min and max dates of coverage. Each nested list (n th)
                     contains the paths to files which share a date range (even partially)
                     with the n th interval.
    date_pairs -- List of tuples, each containing the start and end dates of the n th interval (datetime
                  objects in tuples, all in a list)
    
    """

    # Converting date range to timedelta object (hours)
    date_range = (max_date - min_date)
    date_range = date_range.days * 24

    # Counting number of full intervals over date range {int}
    interval_count = date_range // int(interval)

    # Allowing *min_date* to be updated and initializing time difference object
    min_date = min_date
    dtime = timedelta(hours=int(interval))

    # Creating list containing tuples of date ranges [(dt1, dt2), (dt2, dt3) ...]
    date_pairs = []
    for i in range(interval_count):
        date_pairs.append((min_date, min_date + dtime))
        min_date = min_date + dtime

    # Sorting files into intervals
    interval_list = []
    for pair in date_pairs:
        
        # Initializing temporary list to store dates in interval
        temp_list = []

        # Checking if date range of file overlaps with interval (if true, append)
        for filepath in raw_paths:
            initial_date = datetime.strptime(filepath[-35:-21], '%Y%m%d%H%M%S')
            final_date = datetime.strptime(filepath[-20:-6], '%Y%m%d%H%M%S')            

            if (initial_date <= pair[1]) and (final_date >= pair[0]):
                temp_list.append(filepath)

        # Appending list of interval-contained files        
        interval_list.append(temp_list)      

    return {'interval_list': interval_list, 'date_pairs': date_pairs}

# Filters raw data based on user set parameters
def filter_data(start_year, start_month, start_day, end_year, end_month, end_day, timestep, tolerance, data_path):
    """
    Filters through the data files located in 'data_path' using the user 
    options in 'options.ini'. Outputs a list of paths to data files which
    satisfy the user's criteria.

    Automatically changes date range to match data availability. i.e. if the user specifies
    a date range between 01-11-2020 and 01-06-2021, but data is only available from 
    05-11-2020 and 24-05-2021, dates to be processed will be set to the latter, and
    the user will be notified (Line 164)

    INPUTS:
    start_year -- Starting year YYYY {str}
    start_month -- Starting month MM {str}
    start_day -- Starting day DD {str}

    end_year -- Ending year YYYY {str}
    end_month -- Ending month MM {str}
    end_day -- Ending day DD {str} 

    timestep -- Desired timestep in hours {str}
    tolerance -- Number of hours around timestep that will be filtered through {str}

    data_path -- Path to directory containing data files {str}

    OUTPUTS:
    dictionary object -- raw_paths: List of file paths, max_date: Latest date in file list, 
                         min_date: Earliest date in file list
    """

    # Concatenate start and end dates
    start_date = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
    end_date = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

    # Set delta t tolerance 
    upper_timestep = timedelta(hours=(int(timestep) + int(tolerance)))
    lower_timestep = timedelta(hours=(int(timestep) - int(tolerance)))

    # Initializing file list and date count variables
    raw_paths = []
    min_date = datetime(3000, 12, 25)
    max_date = datetime(1000, 12, 25)

    # Filtering data files by date
    for filename in os.listdir(data_path):
        
        # Extracting initial and final dates from data file names
        initial_date = datetime.strptime(filename[6:20], '%Y%m%d%H%M%S')
        final_date = datetime.strptime(filename[21:35], '%Y%m%d%H%M%S')

        # Checking if all files from initial_date to final_date will be loaded (timestep == '0')
        if timestep != '0':
            # Filtering by date range and delta t and appending to the file list
            if start_date.date() <= initial_date.date() <= end_date.date() and start_date.date() <= final_date.date() <= end_date.date() and lower_timestep <= (final_date-initial_date) <= upper_timestep: 
                raw_paths.append(data_path + '/' + filename)

                # Updating date tracker
                if initial_date < min_date:
                    min_date = initial_date
                if final_date > max_date:
                    max_date = final_date    
        
        elif timestep == '0':
            # Filtering by date range only
            if start_date.date() <= initial_date.date() <= end_date.date() and start_date.date() <= final_date.date() <= end_date.date(): 
                raw_paths.append(data_path + '/' + filename)
                
                # Updating date tracker
                if initial_date < min_date:
                    min_date = initial_date
                if final_date > max_date:
                    max_date = final_date

    # Notifying user of date range change
    if start_date != min_date or end_date != max_date:
        print(f'Start and end dates of data updated to {min_date} and {max_date}')
            
    return {'raw_paths': raw_paths, 'max_date': max_date, 'min_date': min_date}


"""
netCDF Analysis
"""
# Converts desired time range to seconds for netCDF analysis
def seconds_to_date(path:str, start_year, start_month, start_day, end_year, end_month, end_day):
    """
    This function takes in a netCDF's reference time and start and end times of each triangle
    and calculates the seconds ellapsed between the former and latter. This is done to filter
    the data loaded from the netCDF by time, as the start and end times in the netCDF are stored
    as "seconds after the reference time".

    INPUTS:
    path -- String of path to the netCDF file the data will be read from {str}

    start_year -- Start year of the user-specified period (Hereinafter the "period") {str}
    start_month -- Start month of the period {str}
    start_day -- Start day of the period {str}

    end_year -- End year of the period {str}
    end_month -- End month of the period {str}
    end_day -- End day of the period {str}
    """

    # Open dataset
    ds = Dataset(path, mode='r')

    # Fetch reference time (start timestamp)
    reftime = ds.getncattr('referenceTime')

    # Closing dataset
    ds.close()

    # Converting reference time from string to datetime object
    reftime = reftime[0:10]
    reftime = datetime.strptime(reftime, '%Y-%m-%d')

    # Converting start and end date strings to datetime objects
    start_date = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
    end_date = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

    # Taking the difference between the start and end times of each triangle
    # and converting it to seconds
    start_ref_dt = (start_date - reftime).total_seconds()
    end_ref_dt = (end_date - reftime).total_seconds()

    return start_ref_dt, end_ref_dt

# Filters netCDF data by area
def filter_area(centre_lat, centre_lon, radius, start_lats, start_lons):
    """
    This function filters data read from a netCDF (Stored as NumPy arrays) by area.
    The user specifies a centre point in WGS 84 (Lat/Lon) and a radius in km, which are then
    combined to create a mask, leaving only data *radius* kilometres from the centre point.
    Only triangles with all three vertices within this circle will be loaded.

    INPUTS:
    centre_lat -- Latitude of the centre point {float}
    centre_lon -- Longitude of the centre point {float}
    radius -- Radius of the circle (mask) in kilometres {float}
    start_lats -- Starting latitudes of the triangles in the format 
                  [[Latitudes (1)], [Latitudes (2)], [Latitudes (3)]] {NumPy array}
    start_lons -- Starting latitudes of the triangles in the format 
                  [[Longitudes (1)], [Longitudes (2)], [Longitudes (3)]] {NumPy array}

    OUTPUTS:
    indices -- List of indices of triangles located within the area {list}
    """

    print('--- Filtering area ---')

    # Turning inputs into floats if they aren't already
    centre_lat, centre_lon, radius = float(centre_lat), float(centre_lon), float(radius)
    centre_point = (centre_lat, centre_lon)

    # Loading start lat/lon values
    start_lat1, start_lat2, start_lat3 = start_lats[0], start_lats[1], start_lats[2]
    start_lon1, start_lon2, start_lon3 = start_lons[0], start_lons[1], start_lons[2]

    # List of 1s and 0s, 1 if the point is within the radius and 0 else
    tf_list1, tf_list2, tf_list3 = ([] for i in range(3))

    # Iterating over all lists
    for coordinate in zip(start_lat1, start_lon1):
        if hs.haversine(coordinate, centre_point) <= radius: # hs.haversine calculates the 
            tf_list1.append(1)                               # distance between two lats/lons
        else:
            tf_list1.append(0)
    
    for coordinate in zip(start_lat2, start_lon2):
        if hs.haversine(coordinate, centre_point) <= radius:
            tf_list2.append(1)
        else:
            tf_list2.append(0)

    for coordinate in zip(start_lat3, start_lon3):
        if hs.haversine(coordinate, centre_point) <= radius:
            tf_list3.append(1)
        else:
            tf_list3.append(0)

    # Filtering indices where all three points (entire triangle) is within the radius
    indices = [i for i in range(len(tf_list1)) if tf_list1[i] == tf_list2[i] == tf_list3[i] == 1]

    return indices

# Loads netCDF data
def load_netcdf(path:str):
    """
    This function reads and loads data from a netCDF file in the same format as those output by the deformation
    calculation script in src/SeaIceDeformation. The user is able to filter the data by time and area,
    allowing for analytical tools to be applied to selected snippets of data.
    
    INPUTS:
    path -- String of path to the netCDF file the data will be read from {str}

    OUTPUTS:
    dictionary object -- Dictionary storing filtered data as NumPy arrays {dict}
    """

    print('--- Loading data ---')
    
    # Reading config
    config = read_config()
    options = config['options']

    # Reading user options
    start_year, start_month, start_day = options['start_year'], options['start_month'], options['start_day']
    end_year, end_month, end_day = options['end_year'], options['end_month'], options['end_day']
    area_filter, centre_lat, centre_lon, radius = options['area_filter'], options['centre_lat'], options['centre_lon'], options['radius']

    # Load netCDF as Dataset from *path*
    ds = Dataset(path, mode='r')

    # Get start / end times from user
    start_time_s, end_time_s = seconds_to_date(path, start_year, start_month, start_day, end_year, end_month, end_day)

    # Extracting time variables
    start_time = ds.variables['start_time'][:]
    end_time = ds.variables['end_time'][:]

    # Indices of data in desired time frame
    time_indices = np.where( (start_time > start_time_s) & (start_time < end_time_s) )[0]

    # Extracting data (Filtered by time only)

    start_time = start_time[time_indices]
    end_time = end_time[time_indices]

    start_lat1 = (ds.variables['start_lat1'][:])[time_indices]
    start_lat2 = (ds.variables['start_lat2'][:])[time_indices]
    start_lat3 = (ds.variables['start_lat3'][:])[time_indices]

    start_lon1 = (ds.variables['start_lon1'][:])[time_indices]
    start_lon2 = (ds.variables['start_lon2'][:])[time_indices]
    start_lon3 = (ds.variables['start_lon3'][:])[time_indices]

    end_lat1 = (ds.variables['end_lat1'][:])[time_indices]
    end_lat2 = (ds.variables['end_lat2'][:])[time_indices]
    end_lat3 = (ds.variables['end_lat3'][:])[time_indices]

    end_lon1 = (ds.variables['end_lon1'][:])[time_indices]
    end_lon2 = (ds.variables['end_lon2'][:])[time_indices]
    end_lon3 = (ds.variables['end_lon3'][:])[time_indices]

    div = (ds.variables['div'][:])[time_indices]
    shr = (ds.variables['shr'][:])[time_indices]
    vrt = (ds.variables['vrt'][:])[time_indices]

    idx1 = (ds.variables['idx1'][:])[time_indices]
    idx2 = (ds.variables['idx2'][:])[time_indices]
    idx3 = (ds.variables['idx3'][:])[time_indices]
    no   = (ds.variables['no'][:])[time_indices]

    id_start_lat1 = (ds.variables['id_start_lat1'][:])[time_indices]
    id_start_lat2 = (ds.variables['id_start_lat2'][:])[time_indices]
    id_start_lat3 = (ds.variables['id_start_lat3'][:])[time_indices]

    dudx = (ds.variables['dudx'][:])[time_indices]
    dudy = (ds.variables['dudy'][:])[time_indices]
    dvdx = (ds.variables['dvdx'][:])[time_indices]
    dvdy = (ds.variables['dvdy'][:])[time_indices]

    # Compressing coordinates into arrays of arrays
    start_lats = np.array([start_lat1, start_lat2, start_lat3])
    start_lons = np.array([start_lon1, start_lon2, start_lon3])
    end_lats = np.array([end_lat1, end_lat2, end_lat3])
    end_lons = np.array([end_lon1, end_lon2, end_lon3])

    reftime = ds.getncattr('referenceTime')
    icetracker = ds.getncattr('iceTracker')
    trackingerror = ds.getncattr('trackingError')

    if area_filter == 'True':
        print('--- Filtering by area ---')
        # Filtering by area
        area_indices = filter_area(centre_lat, centre_lon, radius, start_lats, start_lons)

        # Applying filter
        start_time = start_time[area_indices]
        end_time = end_time[area_indices]

        start_lats = np.array([start_lat1[area_indices], start_lat2[area_indices], start_lat3[area_indices]])
        start_lons = np.array([start_lon1[area_indices], start_lon2[area_indices], start_lon3[area_indices]])
        end_lats = np.array([end_lat1[area_indices], end_lat2[area_indices], end_lat3[area_indices]])
        end_lons = np.array([end_lon1[area_indices], start_lon2[area_indices], end_lon3[area_indices]])

        div = div[area_indices]
        shr = shr[area_indices]
        vrt = vrt[area_indices]

        idx1 = idx1[area_indices]
        idx2 = idx2[area_indices]
        idx3 = idx3[area_indices]
        no   = no[area_indices]

        id_start_lat1 = id_start_lat1[area_indices]
        id_start_lat2 = id_start_lat2[area_indices]
        id_start_lat3 = id_start_lat3[area_indices]

        dudx = dudx[area_indices]
        dudy = dudy[area_indices]
        dvdx = dvdx[area_indices]
        dvdy = dvdy[area_indices]

    # Closing dataset
    ds.close()

    return {'start_lats': start_lats, 'start_lons': start_lons, 'end_lats': end_lats, 'end_lons': end_lons, 
            'div': div, 'shr': shr, 'vrt': vrt, 
            'start_time': start_time, 'end_time': end_time, 'time_indices': time_indices, 
            'reftime': reftime, 'icetracker': icetracker, 'trackingerror': trackingerror, 
            'idx1': idx1, 'idx2': idx2, 'idx3': idx3, 'no': no, 
            'start_id1': id_start_lat1, 'start_id2': id_start_lat2, 'start_id3': id_start_lat3, 
            'dudx': dudx, 'dudy': dudy, 'dvdx': dvdx, 'dvdy': dvdy}
