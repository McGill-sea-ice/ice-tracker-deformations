'''
Author: Beatrice Duval (bdu002)

-----------------------------------------
Configuration script for data processing
-----------------------------------------

- Retrieves configuration arguments from namelist.ini
- Selects which raw datasets are to be processed using namelist.ini arguments
- Obtain paths to which data files of all stages of data processing are to be stored

'''

import configparser
import os
import re
import sys
from datetime import datetime, timedelta

import utils_datetime
import utils_get_data_paths as get_data_paths


'''
_______________________________________________________________________
DEFINE ERROR CLASS
'''


# Create a class of errors for datasets selection
class datasetSelectionError(Exception):
    pass


'''
_______________________________________________________________________
DEFINE CONFIG FUNCTIONS
'''

def get_config_args():
    ''' None -> ConfigParser

    Function that reads the namelist.ini file and returns a
    ConfigParser object, assuming the namelist.ini file is
    located under the current directory.

    '''

    # Retrieve the path of the src folder, i.e. the current directory
    srcPath = os.path.dirname(os.path.realpath(__file__))

    # Read the namelist.ini file using configparser
    config = configparser.ConfigParser()
    config.read(srcPath + '/namelist.ini')

    # Return a ConfigParser object
    return config


def filter_data(Date_options = None, IO = None, Metadata = None):
    """
    Filters through the data files located in 'data_path' using the user
    options in 'options.ini'. Outputs a list of paths to data files which
    satisfy the user's criteria.

    Automatically changes date range to match data availability. i.e. if the user specifies
    a date range between 01-11-2020 and 01-06-2021, but data is only available from
    05-11-2020 and 24-05-2021, dates to be processed will be set to the latter, and
    the user will be notified (Line

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
    raw_paths -- List of file paths {list}
    """

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
    satellite = Metadata['ice_tracker']

    if satellite == 'RCMS1':
        sat_list = ['rcm/','s1/']
    elif satellite == 'S1':
        sat_list = ['s1/']
    elif satellite == 'RCM':
        sat_list = ['rcm/']
    else:
        sys.exit("Oh, original, but satellite data other than RCM or S1 is not defined!!")

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
                    iDate = datetime.strptime(filename[6:20], '%Y%m%d%H%M%S')
                    fDate = datetime.strptime(filename[21:35], '%Y%m%d%H%M%S')

                    # Checking if all files from iDate to fDate will be loaded (timestep == '0')
                    if timestep != '0':
                        # Filtering by date range and delta t and appending to the file list
                        if sDate.date() <= iDate.date() <= eDate.date() and sDate.date() <= fDate.date() <= eDate.date() and lower_timestep <= (fDate-iDate) <= upper_timestep and iDate.month not in summer_months:
                            raw_paths.append(data_path + '/' + filename)

                            # Updating date tracker
                            if iDate < min_date:
                                min_date = iDate
                            if fDate > max_date:
                                max_date = fDate

                    elif timestep == '0':
                        # Filtering by date range only
                        if sDate.date() <= iDate.date() <= eDate.date() and sDate.date() <= fDate.date() <= eDate.date() and iDate.month not in summer_months:
                            raw_paths.append(data_path + '/' + filename)

                            # Updating date tracker
                            if iDate < min_date:
                                min_date = iDate
                            if fDate > max_date:
                                max_date = fDate

    # Notifying user of date range change
    if sDate != min_date or eDate != max_date:
        print(f"Start and end dates of data updated to {min_date} and {max_date}")

    return raw_paths


def get_datapaths(raw_paths, Date_options = None, IO = None, Metadata = None):
    # def get_datapaths(raw_paths, output_path, exp, start_year, start_month, start_day, ice_tracker):
    ''' (str, str, str, str) -> dict[str, Any]

    Function that creates lists that store data file paths for every dataset
    and for every stage of data processing (triangulation, conversion and
    calculations), as well as the output netcdf file path.

    Returns a dictionnary of the lists of .csv file paths (one per stage of data
    processing) and of the output netcdf file path.

    Every raw data file (e.g. pairs_20200320010317_20200401010317_1.dat) is
    associated to a triangulated (e.g. tri_20200320010317_20200401010317_1.csv),
    and calculated (e.g. calc_20200320010317_20200401010317_1.csv)
    .csv file.

    The output netcdf file combines the output data of all processed datasets
    (listed in the input *raw_paths*). Its name indicates the common starting times
    of all datasets (e.g. RCMS1SID_20200301_dx.nc).

    All output files are stored under output_path/exp, and under a folder
    associated to the stage of dataprocessing.

    For example, given a raw data file 'pair_A.dat', the output triangulation file
    will be stored under 'output_path/exp/02_triangulated/tri_A.csv', the output
    converted file under 'output_path/exp/03_converted/conv_A.csv', the output
    calculations files under 'output_path/exp/04_calculations/calc_A.csv', and the
    output netcdf file under 'output_path/exp/05_output/RCMS1SID_YYYYMMDD_dx.nc',
    where YYYY, MM and DD refer to the input *start_year*, *start_month* and *start_day*,
    respectively.

    Keyword arguments: \\
    raw_paths   -- list of raw datasets' absolute paths \\
    output_path -- path in which all output data files are to be stored \\
    exp         -- name of the experiment
    start_year  -- starting year of the raw datasets listed in raw_paths \\
    start_month -- starting month of the raw datasets listed in raw_paths \\
    start_day   -- starting day of the raw datasets listed in raw_paths \\
    ice_tracker -- sea-ice motion tracker used for calculations \\
    '''

    # Raise an error if the input raw paths list is empty
    if raw_paths ==  []:
        raise datasetSelectionError('The list of raw datasets to process is empty.')

    # Initialize lists of data paths for the subsequent stages of data processing
    triangulated_paths  = []
    calculations_paths  = []

    # Iterate through all input raw file paths
    for raw in raw_paths:

        # Retrive the raw filename
        raw_filename = os.path.basename(raw)

        # For each raw file path, find the appropriate path for each stage of data processing
        tri  = get_data_paths.get_triangulated_csv_path(raw_filename, IO['output_folder'], IO['exp'])  # path for triangulated data file
        cal  = get_data_paths.get_calculations_csv_path(raw_filename, IO['output_folder'], IO['exp'])  # path for calculated data file

        # Append the file paths to the file path lists
        triangulated_paths.append(tri)
        calculations_paths.append(cal)

    # Find the appropriate path for the output netcdf file
    nc_output_path = get_data_paths.get_output_nc_path(IO, Date_options, Metadata)

    # Create the output data paths dictionnary
    # we are using to process data
    data_paths =  { 'raw': raw_paths,
                    'triangulated': triangulated_paths,
                    'calculations': calculations_paths,
                    'nc_output': nc_output_path,
                    'output_folder': IO['output_folder']+IO['exp']}

    return data_paths



'''
_______________________________________________________________________
PERFORM CONFIGURATION
'''

# Retrieve configuration arguments from namelist.ini
config = get_config_args()

# Retrieve the IO and Date_options sections
IO = config['IO']
Date_options = config['Date_options']
Exp_options = config['Metadata']

raw_paths = filter_data( Date_options = Date_options,
                         IO           = IO,
                         Metadata     = Exp_options           )

# Get the paths to which data files of all stages of data processing will be stored
data_paths = get_datapaths( raw_paths,
                            Date_options = Date_options,
                            IO           = IO,
                            Metadata     = Exp_options)
