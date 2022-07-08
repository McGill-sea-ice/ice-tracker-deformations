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
