"""
Author: Lekima Yakuden
GitHub: LekiYak

--------------------------------------------------------------------------------
Configuration file for data / netCDF analysis tools
--------------------------------------------------------------------------------

This file contains functions for loading and processing user options, raw data, and netCDF files.
"""

# import from default packages
import os
import configparser
from datetime import datetime, timedelta
import pandas as pd
from tqdm import tqdm

# Imports from other files
from SatelliteCoverage.utils import stb

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

    # Choose the default or user file
    def_fname = '/options.def'
    usr_fname = '/options.ini'
    if os.path.exists(cwd + usr_fname):
        fname = usr_fname
    elif os.path.exists(cwd + def_fname):
        print('--- Using default parameters options.def ---')
        print('--- Create the options.ini file to define user parameters ---')
        fname = def_fname
    else:
        print('/!/ No config file found! /!/')

    #opening the file
    config.read(cwd + fname)

    # Return a ConfigParser object
    config_dict = {sect: dict(config.items(sect)) for sect in config.sections()}

    # Iterate through the dictionnary to change strings to bools
    for key in config_dict:
        if isinstance(config_dict[key], dict):
            for keyy in config_dict[key]:
                if isinstance(config_dict[key][keyy], str):
                    config_dict[key][keyy] = stb(config_dict[key][keyy])
        elif isinstance(config_dict[key], str):
            config_dict[key] = stb(config_dict[key])

    return config_dict

"""
Raw Data
"""
# Compiles raw data into pandas dataframe
def compile_data(raw_paths=None):

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
    satellite = 'RCM'
    # Appending each file's datapoints to the dataframe
    for filepath in tqdm(raw_paths, position=1, leave=False):

        if satellite == 'RCM_new':
            # Initialize temporary dataframe
            temp_df = pd.DataFrame()
            with open(filepath) as fd:
                headers = [ next(fd) for i in range(7) ]
                temp_df = pd.read_csv(fd,engine='python',sep = '\s\s+', usecols = ['lat_beg','lon_beg'])
            my_list = list(temp_df)
            temp_df.drop([0], axis=0, inplace=True)
            temp_df.rename(columns = {'lat_beg':'lat','lon_beg':'lon'}, inplace = True)

        else:
            # Reading datapoints into temporary dataframe
            temp_df = pd.read_csv(filepath, sep='\s\s+', engine='python', usecols = ['sLat','sLon'])
            temp_df.rename(columns = {'sLat':'lat', 'sLon':'lon'}, inplace = True)

        df = pd.concat([df,temp_df])

        # Updating counter
        i += 1
        # print(f'{i} / {num_files}')

    return df

# Divides raw data into intervals specified by the user
# def divide_intervals(raw_paths, max_date, min_date, interval):
def divide_intervals(config=None):
    """
    Divides delta-t filtered data into chunks (intervals) of *interval* hours for
    processing.

    INPUTS:
    raw_paths -- List of data file paths with the desired delta t {list}
    Date_options -- ConfigParser object to determine start/end dates of intervals as specififed by user
    Options -- ConfigParser object to determine the interval length as specififed by user

    OUTPUTS:
    interval_list -- List of n lists, where n is the number of full intervals which
                     fit in the min and max dates of coverage. Each nested list (n th)
                     contains the paths to files which share a date range (even partially)
                     with the n th interval.
    date_pairs -- List of tuples, each containing the start and end dates of the n th interval (datetime
                  objects in tuples, all in a list)

    """

    raw_paths = config['raw_list']
    Date_options = config['Date_options']
    options = config['options']
    satellite = config['Metadata']['icetracker']

    start_year  = str(Date_options['start_year'])
    start_month = str(Date_options['start_month'])
    start_day   = str(Date_options['start_day'])
    end_year    = str(Date_options['end_year'])
    end_month   = str(Date_options['end_month'])
    end_day     = str(Date_options['end_day'])
    interval    = options['interval']


    # Concatenate start and end dates
    sDate = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
    eDate = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

    # Converting date range to timedelta object (hours)
    date_range = (eDate - sDate)
    date_range = date_range.days * 24

    # Counting number of full intervals over date range {int}
    interval_count = date_range // int(interval)

    # Allowing *min_date* to be updated and initializing time difference object
    min_date_it = sDate
    dtime = timedelta(hours=int(interval))

    # Creating list containing tuples of date ranges [(dt1, dt2), (dt2, dt3) ...]
    date_pairs = []
    for i in range(interval_count):
        date_pairs.append((min_date_it, min_date_it + dtime))
        min_date_it = min_date_it + dtime

    # Sorting files into intervals
    interval_list = []
    for pair in date_pairs:

        # Initializing temporary list to store dates in interval
        temp_list = []

        # Checking if date range of file overlaps with interval (if true, append)
        for filepath in raw_paths:
            if satellite == 'RCM_new':
                txt = filepath.split('_')
                initial_date = datetime.strptime("%s%s" % (txt[9],txt[10]), '%Y%m%d%H%M%S')
                final_date = datetime.strptime("%s%s" % (txt[19],txt[20]), '%Y%m%d%H%M%S')
            else:
                initial_date = datetime.strptime(filepath[-35:-21], '%Y%m%d%H%M%S')
                final_date = datetime.strptime(filepath[-20:-6], '%Y%m%d%H%M%S')

            if (initial_date <= pair[1]) and (final_date >= pair[0]):
                temp_list.append(filepath)

        # Appending list of interval-contained files
        interval_list.append(temp_list)

    return interval_list, date_pairs


# Filters raw data based on user set parameters
def filter_data(config=None):
    """
    Filters through the data files located in 'data_path' using the user
    options in 'options.ini'. Outputs a list of paths to data files which
    satisfy the user's criteria.

    Note: If the available data does not span the full start/end dates interval
    specified by the user in 'options.ini', the user will be informed that the
    the files to be processed have a different time span the the one requested.

    INPUTS:
    Date_options -- ConfigParser object to determine start/end dates, timestep, and tolerance as specififed by user.
    IO -- ConfigParser object to determine the path to directory containing the data files
    Metadata -- ConfigParser object to determine which satellite(s) is being used.

    OUTPUTS:
    raw_paths -- List of file paths {list}
    """

    IO = config['IO']
    Date_options = config['Date_options']
    Metadata = config['Metadata']

    start_year  = str(Date_options['start_year'])
    start_month = str(Date_options['start_month'])
    start_day   = str(Date_options['start_day'])
    end_year    = str(Date_options['end_year'])
    end_month   = str(Date_options['end_month'])
    end_day     = str(Date_options['end_day'])
    timestep    = int(Date_options['timestep'])
    tolerance   = int(Date_options['tolerance'])

    # Concatenate start and end dates
    sDate = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
    eDate = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

    # Set delta t tolerance
    upper_timestep = timedelta(hours=(int(timestep) + int(tolerance)))
    lower_timestep = timedelta(hours=(int(timestep) - int(tolerance)))

    # Initializing file list and date count variables
    raw_paths = []
    min_date = datetime(3000, 12, 25)
    max_date = datetime(1000, 12, 25)

    # List of summer months (Won't include data from these months)
    summer_months = [6, 7, 8, 9, 10]

    date_path = IO['data_folder']
    satellite = Metadata['icetracker']

    if satellite == 'RCMS1':
        sat_list = ['rcm/','s1/']
    elif satellite == 'S1':
        sat_list = ['s1/']
    elif satellite == 'RCM':
        sat_list = ['rcm/']
    elif satellite == 'RCM_new':
        sat_list = ['rcm_new/']
    else:
        raise Exception("Oh, original, but satellite data other than RCM or S1 is not defined!!")

    for sat_type in sat_list:
        for year in range(int(start_year), int(end_year)+1):

            if not(os.path.exists(date_path + sat_type + str(year) + '/')):
                print('No data for '+ sat_type + ' in ' + str(year) )
            else:
                data_path = date_path + sat_type + str(year) + '/'

                #--------------------------------------------------------------------------
                # listing the files in the folder and adding the pairs if in the right dates
                #     This could be made a function if used when a IO class is made
                #---------------------------------------------------------------------------

                # Filtering data files by date
                for filename in os.listdir(data_path):

                    # Extracting initial and final dates from data file names
                    if satellite == 'RCM_new':
                        txt = filename.split('_')
                        iDate = datetime.strptime("%s%s" % (txt[5],txt[6]), '%Y%m%d%H%M%S')
                        fDate = datetime.strptime("%s%s" % (txt[15],txt[16]), '%Y%m%d%H%M%S')
                    else:
                        iDate = datetime.strptime(filename[6:20], '%Y%m%d%H%M%S')
                        fDate = datetime.strptime(filename[21:35], '%Y%m%d%H%M%S')

                    # Checking if all files from iDate to fDate will be loaded (timestep == '0')
                    if timestep != '0':
                        # Filtering by date range and delta t and appending to the file list
                        if iDate < eDate and fDate > sDate and lower_timestep <= (fDate-iDate) <= upper_timestep and iDate.month not in summer_months:
                            raw_paths.append(data_path + '/' + filename)

                            # Updating date tracker
                            if iDate < min_date:
                                min_date = iDate
                            if fDate > max_date:
                                max_date = fDate

                    elif timestep == '0':
                        # Filtering by date range only

                        if iDate < eDate and fDate > sDate and iDate.month not in summer_months:

                            raw_paths.append(data_path + '/' + filename)

                            # Updating date tracker
                            if iDate < min_date:
                                min_date = iDate
                            if fDate > max_date:
                                max_date = fDate


    # Notifying user of data start/end times of files corresponding to the specified interval
    print(f"Valid start/end dates of coverage map: {sDate} to {eDate}")
    print(f"Earliest/latest start/end dates of data included in coverage map: {min_date} to {max_date}")

    return raw_paths

