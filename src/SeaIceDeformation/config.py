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



def select_ds(folder, start_year, start_month, start_day, duration):
    ''' (str, str, str, str, float) -> list(str)
    
    Function that selects the raw datasets to process. 
    
    The selected raw datasets are all under *folder*, start at 
    *start_year*-*start_month*-*start_day*, and span over a period of 
    *duration* days, +/- 6 hours.

    Keyword arguments: \\
    folder      -- absolute path of the folder that contains the raw
                  data files over which the selection process will 
                  be performed  \\
    start_year  -- starting year of the raw datasets to select for data processing \\
    start_month -- starting month of the raw datasets to select for data processing \\
    start_day   -- starting day of the raw datasets to select for data processing \\
    duration    -- time span of the raw datasets to select for data processing

    '''
    # Initialize a string of raw data paths to select for processing
    raw_paths = []
    
    # Create constants for the beginning and the end of the regex for raw data file names
    BEG = '^pairs_'
    END = '[0-2][0-9][0-6][0-9][0-6][0-9]_[1-2][0-9][0-9][0-9][0-1][0-9][0-3][0-9][0-2][0-9][0-6][0-9][0-6][0-9]_[0-9].dat$'
    
    # Create a regex sequence for raw data file names using the input start time
    regex = re.compile('{0}{1}{2}{3}{4}'.format(BEG, start_year, start_month, start_day, END))

    # Iterate through all files in the input folder to find matches
    for filename in os.listdir(folder):
        
        # Check if the current file matches the regex
        if regex.match(filename):
            
            # Compute the time over which the dataset spans (days)
            dt = utils_datetime.dT( utils_datetime.dataDatetimes(filename) ) / 86400
            
            # Check if the current file spans over a time that is approx. 
            # equal to the input duration (± 6 hours, or ± 0.25 days). 
            if duration-0.25 < dt < duration+0.25:

                # The current raw file starts at the input date and spans over 
                # the input duration +/- 6 hours.
                # Append the raw data path to the raw_paths list
                raw_paths.append(folder + '/' + filename) 

    return raw_paths


def get_datapaths(raw_paths, start_year, start_month, start_day, method):
    ''' (str, str, str, str) -> dict[str, Any]

    Function that creates the lists that store data file paths for every dataset
    and for every stage of data processing (triangulation, conversion and 
    calculations), as well as the output netcdf file path.

    Returns a dictionnary of the lists of .csv file paths and of the output 
    netcdf file path.

    Every raw data file (e.g. pairs_20200320010317_20200401010317_1.dat) is 
    associated to a triangulated (e.g. tri_20200320010317_20200401010317_1.csv), 
    converted (e.g. tri_gridCS_20200320010317_20200401010317_1.csv) (for M00 
    method only) and calculated (e.g. calc_20200320010317_20200401010317_1.csv) 
    .csv file.

    The output netcdf file combines the output data of all processed datasets 
    (listed in the input *raw_paths*). Its name indicates the common starting times 
    of all datasets (e.g. RCMS1SID_20200301_dx.nc).

    The files are stored in the following tree structure, where 
    some_parent_folder can have any name, but will be the same for all
    data files (i.e. triangulated, converted, calculations, output) 
    that result from a common raw data file (which is originally stored 
    under some_parent_folder) :

    ├── 2021_SeaIceDeformations
    |    ...
    |    ├── data
    |    |   ├── 00_grid /etc
    |    |   ├── 01_raw
    |    |   |   └── some_parent_folder
    |    |   |       └── pairs_20200320010317_20200401010317_1.csv
    |    |   ├── 02_triangulated
    |    |   |   └── some_parent_folder
    |    |   |       └── tri_20200320010317_20200401010317_1.csv
    |    |   ├── 03_converted
    |    |   |   └── some_parent_folder
    |    |   |       └── tri_gridCS_20200320010317_20200401010317_1.csv
    |    |   ├── 04_calculations
    |    |   |   └── some_parent_folder
    |    |   |       └── calc_20200320010317_20200401010317_1.csv
    |    |   └──05_output
    |    |       └── some_parent_folder
    |    |           └── RCMS1SID_20200320010317_dx.nc
    ...
    
    Keyword arguments: \\
    raw_paths   -- list of raw datasets' absolute paths \\
    start_year  -- starting year of the raw datasets listed in raw_paths \\
    start_month -- starting month of the raw datasets listed in raw_paths \\
    start_day   -- starting day of the raw datasets listed in raw_paths \\
    method      -- method that will be used to process datasets
    '''

    # Raise an error if the input raw paths list is empty
    if raw_paths ==  []:
        raise datasetSelectionError('The list of raw datasets to process is empty.')
    
    # Initialize lists of data paths for the subsequent stages of data processing
    triangulated_paths  = []
    calculations_paths  = []

    # If we are processing data using method M00, initialize an additionnal list for 
    # converted data paths
    if method == 'M00':
        converted_paths     = []

    # Iterate through all input raw file paths
    for raw in raw_paths:
     
        # For each raw file path, find the appropriate path for each stage of data processing
        tri  = get_data_paths.get_processed_csv_path(raw)     # path for triangulated data file
        cal  = get_data_paths.get_calculations_csv_path(raw)  # path for calculated data file

        # Append the file paths to the file path lists
        triangulated_paths.append(tri)
        calculations_paths.append(cal)

        # If we are processing data using method M00, find the appropriate path for the 
        # converted file, and append it the converted file path list
        if method == 'M00':
            conv = get_data_paths.get_converted_csv_path(raw)     
            converted_paths.append(conv) 

    # Find the appropriate path for the output netcdf file
    output_path = get_data_paths.get_output_nc_path(raw, start_year, start_month, start_day)

    # Create the output data paths dictionnaty depending on the method 
    # we are using to process data
    if method == 'M00':
        data_paths =  { 'raw': raw_paths, 
                        'triangulated': triangulated_paths, 
                        'converted': converted_paths, 
                        'calculations': calculations_paths, 
                        'output': output_path }
    
    elif method == 'M01':
        data_paths =  { 'raw': raw_paths, 
                        'triangulated': triangulated_paths, 
                        'calculations': calculations_paths, 
                        'output': output_path }        

    return data_paths



'''
_______________________________________________________________________
PERFORM CONFIGURATION
'''

# Retrieve configuration arguments from namelist.ini 
config = get_config_args()

# Retrieve the Date_options section 
Date_options = config['Date_options']

# Select the raw files to be processed using config arguments
raw_paths = select_ds(  config['IO']['raw_data_folder'], 
                        Date_options['start_year'], 
                        Date_options['start_month'], 
                        Date_options['start_day'], 
                        float(Date_options['duration']) )

# Get the paths to which data files of all stages of data processing will be stored
data_paths = get_datapaths( raw_paths, 
                            Date_options['start_year'], 
                            Date_options['start_month'], 
                            Date_options['start_day'], 
                            config['Processing_options']['method'] )


