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



def select_ds(data_folder, start_year, start_month, start_day, duration):
    '''
    
    Function that selects the raw datasets to process. 
    
    The selected raw datasets are all under *data_folder*, start at 
    *start_year*-*start_month*-*start_day*, and span over a period of 
    *duration* days, +/- 6 hours.

    Keyword arguments: \\
    data_folder -- absolute path of the data_folder that contains the raw
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
    for filename in os.listdir(data_folder):
        
        # Check if the current file matches the regex
        if regex.match(filename):
            
            # Compute the time over which the dataset spans (days)
            dt = utils_datetime.dT( utils_datetime.dataDatetimes(filename) ) / 86400
            
            # Check if the current file spans over a time that is approx. 
            # equal to the input duration (+/- 6 hours, or +/- 0.25 days). 
            if duration-0.25 < dt < duration+0.25:

                # The current raw file starts at the input date and spans over 
                # the input duration +/- 6 hours.
                # Append the raw data path to the raw_paths list
                raw_paths.append(data_folder + '/' + filename) 

    return raw_paths


def get_datapaths(raw_paths, output_path, exp, start_year, start_month, start_day, method):
    ''' (str, str, str, str) -> dict[str, Any]

    Function that creates lists that store data file paths for every dataset
    and for every stage of data processing (triangulation, conversion and 
    calculations), as well as the output netcdf file path.

    Returns a dictionnary of the lists of .csv file paths (one per stage of data 
    processing) and of the output netcdf file path.

    Every raw data file (e.g. pairs_20200320010317_20200401010317_1.dat) is 
    associated to a triangulated (e.g. tri_20200320010317_20200401010317_1.csv), 
    converted (e.g. conv_20200320010317_20200401010317_1.csv) (for M00 
    method only) and calculated (e.g. calc_20200320010317_20200401010317_1.csv) 
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
     
        # Retrive the raw filename
        raw_filename = os.path.basename(raw)

        # For each raw file path, find the appropriate path for each stage of data processing
        tri  = get_data_paths.get_triangulated_csv_path(raw_filename, output_path, exp)     # path for triangulated data file
        cal  = get_data_paths.get_calculations_csv_path(raw_filename, output_path, exp)  # path for calculated data file

        # Append the file paths to the file path lists
        triangulated_paths.append(tri)
        calculations_paths.append(cal)

        # If we are processing data using method M00, find the appropriate path for the 
        # converted file, and append it the converted file path list
        if method == 'M00':
            conv = get_data_paths.get_converted_csv_path(raw_filename, output_path, exp)     
            converted_paths.append(conv) 

    # Find the appropriate path for the output netcdf file
    output_path = get_data_paths.get_output_nc_path(output_path, exp, start_year, start_month, start_day)

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

# Retrieve the IO and Date_options sections
IO = config['IO']
Date_options = config['Date_options']

# Select the raw files to be processed using config arguments
raw_paths = select_ds(  IO['data_folder'], 
                        Date_options['start_year'], 
                        Date_options['start_month'], 
                        Date_options['start_day'], 
                        float(Date_options['duration']) )

# Get the paths to which data files of all stages of data processing will be stored
data_paths = get_datapaths( raw_paths, 
                            IO['output_folder'],
                            IO['exp'],
                            Date_options['start_year'], 
                            Date_options['start_month'], 
                            Date_options['start_day'], 
                            config['Processing_options']['method'] )


