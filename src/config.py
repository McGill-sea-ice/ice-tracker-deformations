'''
Author: Beatrice Duval (bdu002)

-----------------------------------------
Configuration script for data processing.
-----------------------------------------

'''
import argparse
import os
import re

import utils_datetime
import utils_get_data_paths as get_data_paths


def load_config(filename):
    ''' (string) -> dict[str, list]

    Function that initializes the lists that store .csv file paths for 
    every stage of data processing (triangulation, conversion and calculations), 
    using the global list of raw .csv data paths.

    Returns a dictionnary of the lists of .csv file paths.

    Every raw .csv file (e.g. pairs_20200320010317_20200401010317_1.csv) is 
    associated to a triangulated (e.g. tri_20200320010317_20200401010317_1.csv), 
    converted (e.g. tri_gridCS_20200320010317_20200401010317_1.csv) and 
    calculated (e.g. calc_20200320010317_20200401010317_1.csv) .csv file.

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
    |    |   └── 04_calculations
    |    |       └── some_parent_folder
    |    |           └── calc_20200320010317_20200401010317_1.csv
    ...
    
    '''
    with open(filename, 'r') as f:
        lines = f.read().splitlines()

    triangulated_paths  = []
    converted_paths     = []
    calculations_paths  = []

    # Iterate through all raw .csv file paths
    for raw_path in lines:
     
        # For each raw .csv file path, find the appropriate path for each subsequent stages of processing
        pr = get_data_paths.get_processed_csv_path(raw_path)     # path for triangulated data file
        cv = get_data_paths.get_converted_csv_path(raw_path)     # path for converted data file
        cl = get_data_paths.get_calculations_csv_path(raw_path)  # path for calculated data file

        # Add the .csv file paths to the file path lists
        triangulated_paths.append(pr)
        converted_paths.append(cv)
        calculations_paths.append(cl)
    
    return {'raw': lines, 'triangulated': triangulated_paths, 'converted': converted_paths, 'calculations': calculations_paths}


def select_ds(foldername, txt, startYear, startMonth, startDay, duration):
    # Initialize a string of data paths to process
    data_paths = ''
    
    # Create constants for the beginning and the end of the regex for raw data file names
    BEG = '^pairs_'
    END = '[0-2][0-9][0-6][0-9][0-6][0-9]_[1-2][0-9][0-9][0-9][0-1][0-9][0-3][0-9][0-2][0-9][0-6][0-9][0-6][0-9]_[0-9].dat$'
    
    # Create a regex sequence for raw data file names using the specified start time in the command-line arguments
    regex=re.compile('{0}{1}{2}{3}{4}'.format(BEG, startYear, startMonth, startDay, END))

    # Iterate through all files in the given directory to find matches
    for filename in os.listdir(foldername):
        
        # Check if the current file matches the regex
        if regex.match(filename):
            
            # Compute the time over which the dataset spans
            dt = utils_datetime.dT( utils_datetime.dataDatetimes(filename) )
            
            # Select the current file for processing if:
            #   1) the duration over which datasets must span is not specified in the command-line argument (is None);
            #   2) it spans over a time that is approx. equal to the specified duration (± 6 hours, or ± 0.25 days). 
            if duration is None or (duration-0.25 < dt < duration+0.25):
                # Append the data path to the datapaths string
                data_paths += foldername + '/' + filename + '\n'

    # Remove the last escape sequence from the string of datasets to process
    if len(data_paths) != 0:
        data_paths = data_paths[0:-1]
    
    # Write the file paths to the .txt file
    with open(txt, "w") as myfile:
        myfile.write(data_paths)
    

def get_args():
    
    p = argparse.ArgumentParser(description='Compute deformations from Sentinel-1 and RCM datasets.')
    
    p.add_argument('--overwrite', action='store_true', help='When the overwrite argument is specified,           \
                                                            all raw data sets that have already been             \
                                                            processed will be re-processed and the resulting     \
                                                            .csv files for all stages of data processing         \
                                                            will be overwritten')
    
    p.add_argument('--method', default='M01', choices=['M00', 'M01'], help='Method to be used when processing  \
                                                                            datasets. M00 refers to the method       \
                                                                            that uses the RIOPS grid, while M01      \
                                                                            refers to the method that uses X/Y       \
                                                                            coordinates from the input datasets.')
    
    p.add_argument('--data_folder', default='/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_RCM', 
                                    help='Absolute path of the folder from which datasets are selected for processing.')
    
    p.add_argument('--select_file', default='/home/bdu002/2021_SeaIceDeformation/dataset_paths.txt', 
                                    help='Absolute path of a .txt file in which are stored the paths of the raw        \
                                        datasets that have been selected to be processed.')
    
    p.add_argument('--keep_select', action='store_true', help='When the keep_select argument is specified, no dataset \
                                                               selection is performed and the select_file is kept as is.')

    p.add_argument('--start_year', default='[1-2][0-9][0-9][0-9]', 
                                   help='Starting year (4 digits format) of the datasets to process. A regular expression \
                                        can be provided, e.g. 201[0-9].')

    p.add_argument('--start_month', default='[0-1][0-9]', 
                                    help='Starting month (2 digits format) of the datasets to process. A regular expression \
                                         can be provided, e.g. 0[16].')

    p.add_argument('--start_day', default='[0-3][0-9]', 
                                  help='Starting day (2 digits format) of the datasets to process. A regular expression \
                                       can be provided, e.g. [0-3]0.')

    p.add_argument('--duration', help='Time span (in days) of the datasets to process.', type=int)

    return p.parse_args()
    
# Retrieve command line arguments
args = get_args()

if not args.keep_select:
    # Select files to be processed and store their paths in a .txt file
    select_ds(args.data_folder, args.select_file, args.start_year, args.start_month, args.start_day, args.duration)

# Load file paths for all stages of data processing
data_paths = load_config(args.select_file)

if __name__ == '__main__':

    from pprint import pprint
    pprint(args)
