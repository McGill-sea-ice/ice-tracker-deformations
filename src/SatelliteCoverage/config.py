import configparser
import os
from datetime import datetime, timedelta
import numpy as np

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

# Reading options
config = read_config()

# Initializing more specific ConfigParser objects
IO = config['IO']
options = config['options']

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
            iDate = datetime.strptime(filepath[-35:-21], '%Y%m%d%H%M%S')
            fDate = datetime.strptime(filepath[-20:-6], '%Y%m%d%H%M%S')            

            if (iDate <= pair[1]) and (fDate >= pair[0]):
                temp_list.append(filepath)

        # Appending list of interval-contained files        
        interval_list.append(temp_list)      

    return {'interval_list': interval_list, 'date_pairs': date_pairs}

def filter_data(sYear, sMonth, sDay, eYear, eMonth, eDay, delta_t, tolerance, data_path):
    """
    Filters through the data files located in 'data_path' using the user 
    options in 'options.ini'. Outputs a list of paths to data files which
    satisfy the user's criteria.

    INPUTS:
    sYear -- Starting year YYYY {str}
    sMonth -- Starting month MM {str}
    sDay -- Starting day DD {str}

    eYear -- Ending year YYYY {str}
    eMonth -- Ending month MM {str}
    eDay -- Ending day DD {str} 

    delta_t -- Desired timestep in hours {str}
    tolerance -- Number of hours around timestep that will be filtered through {str}

    data_path -- Path to directory containing data files {str}

    OUTPUTS:
    dictionary object -- raw_paths: List of file paths, max_date: Highest date in file list, min_date: Lowest date in file list
    """

    # Concatenate start and end dates
    sDate = datetime.strptime(sYear + sMonth + sDay, '%Y%m%d')
    eDate = datetime.strptime(eYear + eMonth + eDay, '%Y%m%d')

    # Set delta t tolerance 
    upper_delta_t = timedelta(hours=(int(delta_t) + int(tolerance)))
    lower_delta_t = timedelta(hours=(int(delta_t) - int(tolerance)))

    # Initializing file list and date count variables
    raw_paths = []
    min_date = datetime(3000, 12, 25)
    max_date = datetime(1000, 12, 25)

    # Filtering data files by date
    for filename in os.listdir(data_path):
        
        # Extracting initial and final dates from data file names
        iDate = datetime.strptime(filename[6:20], '%Y%m%d%H%M%S')
        fDate = datetime.strptime(filename[21:35], '%Y%m%d%H%M%S')

        # Checking if all files from iDate to fDate will be loaded (delta_t == '0')
        if delta_t != '0':
            # Filtering by date range and delta t and appending to the file list
            if sDate.date() <= iDate.date() <= eDate.date() and sDate.date() <= fDate.date() <= eDate.date() and lower_delta_t <= (fDate-iDate) <= upper_delta_t: 
                raw_paths.append(data_path + '/' + filename)

                # Updating date tracker
                if iDate < min_date:
                    min_date = iDate
                if fDate > max_date:
                    max_date = fDate    
        
        elif delta_t == '0':
            # Filtering by date range only
            if sDate.date() <= iDate.date() <= eDate.date() and sDate.date() <= fDate.date() <= eDate.date(): 
                raw_paths.append(data_path + '/' + filename)
                
                # Updating date tracker
                if iDate < min_date:
                    min_date = iDate
                if fDate > max_date:
                    max_date = fDate

    return {'raw_paths': raw_paths, 'max_date': max_date, 'min_date': min_date}