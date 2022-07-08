import configparser
import os
from datetime import datetime, timedelta

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
    raw_paths -- List of paths to filtered data files {list}
    """

    # Concatenate start and end dates
    sDate = datetime.strptime(sYear + sMonth + sDay, '%Y%m%d')
    eDate = datetime.strptime(eYear + eMonth + eDay, '%Y%m%d')

    # Set delta t tolerance 
    upper_delta_t = timedelta(hours=(int(delta_t) + int(tolerance)))
    lower_delta_t = timedelta(hours=(int(delta_t) - int(tolerance)))

    # Initializing empty list to store desired file names
    raw_paths = []
    # Filtering data files by date
    for filename in os.listdir(data_path):
        
        # Extracting initial and final dates from data file names
        iDate = datetime.strptime(filename[6:20], '%Y%m%d%H%M%S')
        fDate = datetime.strptime(filename[21:35], '%Y%m%d%H%M%S')

        # Checking if all files from iDate to fDate will be loaded (delta_t = '0')
        if delta_t != '0':
            # Filtering by date range and delta t and appending to the file list
            if sDate.date() <= iDate.date() <= eDate.date() and sDate.date() <= fDate.date() <= eDate.date() and lower_delta_t <= (fDate-iDate) <= upper_delta_t: 
                raw_paths.append(data_path + '/' + filename)
        
        elif delta_t == '0':
            # Filtering by date range only
            if sDate.date() <= iDate.date() <= eDate.date() and sDate.date() <= fDate.date() <= eDate.date(): 
                raw_paths.append(data_path + '/' + filename)            

    return raw_paths



data_path = IO['data_folder']

sYear = options['start_year']
sMonth = options['start_month']
sDay = options['start_day']

eYear = options['end_year']
eMonth = options['end_month']
eDay = options['end_day']

delta_t = options['delta_t']
tolerance = options['tolerance']

print(len(filter_data(sYear, sMonth, sDay, eYear, eMonth, eDay, delta_t, tolerance, data_path)))
