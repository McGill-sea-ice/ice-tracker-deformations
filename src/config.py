'''
Author: Beatrice Duval (bdu002)

-----------------------------------------
Configuration script for data processing
-----------------------------------------

'''

import configparser
import os
import re

import utils_datetime
import utils_get_data_paths as get_data_paths


# Create a class of errors for datasets selection
class datasetSelectionError(Exception):
    pass


def get_config_args():
    
    # Retrieve the path of the src folder
    srcPath = os.path.dirname(os.path.realpath(__file__))
    
    # Read the namelist.ini file using configparser
    config = configparser.ConfigParser()
    config.read(srcPath + '/namelist.ini')

    return config



def select_ds(foldername, txt, startYear, startMonth, startDay, duration):
    # Initialize a string of data paths to process
    raw_paths = []
    
    # Create constants for the beginning and the end of the regex for raw data file names
    BEG = '^pairs_'
    END = '[0-2][0-9][0-6][0-9][0-6][0-9]_[1-2][0-9][0-9][0-9][0-1][0-9][0-3][0-9][0-2][0-9][0-6][0-9][0-6][0-9]_[0-9].dat$'
    
    # Create a regex sequence for raw data file names using the specified start time in the command-line arguments
    regex = re.compile('{0}{1}{2}{3}{4}'.format(BEG, startYear, startMonth, startDay, END))

    # Iterate through all files in the given directory to find matches
    for filename in os.listdir(foldername):
        
        # Check if the current file matches the regex
        if regex.match(filename):
            
            # Compute the time over which the dataset spans
            dt = utils_datetime.dT( utils_datetime.dataDatetimes(filename) ) / 86400
            
            # Select the current file for processing if:
            #   1) the duration over which datasets must span is not specified in the command-line argument (is None);
            #   2) it spans over a time that is approx. equal to the specified duration (± 6 hours, or ± 0.25 days). 
            if duration is None or (duration-0.25 < dt < duration+0.25):
                # Append the data path to the datapaths string
                raw_paths.append(foldername + '/' + filename) 

    return raw_paths


def get_datapaths(raw_paths, start_year, start_month, start_day):
    ''' (string) -> dict[str, list]

    Function that initializes the lists that store .csv file paths for 
    every stage of data processing (triangulation, conversion and calculations),
    and the output netcdf file path, using the global list of raw .csv data paths.

    Returns a dictionnary of the lists of .csv file paths and of the output 
    netcdf file path.

    Every raw .csv file (e.g. pairs_20200320010317_20200401010317_1.csv) is 
    associated to a triangulated (e.g. tri_20200320010317_20200401010317_1.csv), 
    converted (e.g. tri_gridCS_20200320010317_20200401010317_1.csv) and 
    calculated (e.g. calc_20200320010317_20200401010317_1.csv) .csv file.

    The output netcdf file combines the output data of all processed datasets 
    listed in a .txt file. Its name indicates the common starting times 
    of all datasets.

    The files are stored in the following tree structure, where 
    some_parent_folder can have any name, but will be the same for all
    .csv files resulting from a common raw .csv file :

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
    
    '''

    # Raise an error if the raw paths list is empty
    if raw_paths ==  []:
        raise datasetSelectionError( 'The list of raw datasets to process is empty.')
    
    # Initialize lists of data paths for the subsequent stages of data processing
    triangulated_paths  = []
    converted_paths     = []
    calculations_paths  = []

    # Iterate through all raw .csv file paths
    for raw in raw_paths:
     
        # For each raw .csv file path, find the appropriate path for each subsequent stages of processing
        tri = get_data_paths.get_processed_csv_path(raw)     # path for triangulated data file
        conv = get_data_paths.get_converted_csv_path(raw)     # path for converted data file
        cal = get_data_paths.get_calculations_csv_path(raw)  # path for calculated data file

        # Add the .csv file paths to the file path lists
        triangulated_paths.append(tri)
        converted_paths.append(conv)
        calculations_paths.append(cal)
    
    # Find the appropriate path for the output netcdf file
    output_path = get_data_paths.get_output_nc_path(raw, start_year, start_month, start_day)

    return { 'raw': raw_paths, 
             'triangulated': triangulated_paths, 
             'converted': converted_paths, 
             'calculations': calculations_paths, 
             'output': output_path }


# Retrieve namelist.ini config arguments
config = get_config_args()
Date_options = config['Date_options']

# TEMPORARY
select_file = '/home/bdu002/2021_SeaIceDeformation/dataset_paths.txt'

# Select files to be processed and store their paths in a .txt file
raw_paths = select_ds(config['IO']['raw_data_folder'], select_file, Date_options['start_year'], Date_options['start_month'], Date_options['start_day'], float(Date_options['duration']))

# Load file paths for all stages of data processing
data_paths = get_datapaths(raw_paths, Date_options['start_year'], Date_options['start_month'], Date_options['start_day'])

if __name__ == '__main__':
    projPath = os.path.dirname(os.path.realpath(__file__))
    config = configparser.ConfigParser()
    print(config.read(projPath + '/namelist.ini'))
    print(config.sections())

