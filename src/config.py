'''
Author: Beatrice Duval (bdu002)

-----------------------------------------
Configuration script for data processing.
-----------------------------------------

Argparse Parameters:
overwrite -- (True, False) If overwrite is set to true, files that have already 
              been processed will be processed again and overwritten. 
              If overwrite is set to false, files that have aleady been processed 
              will not be processed during the program's execution.

'''

import argparse
from utils_get_data_paths import (get_calculations_csv_path,
                                              get_converted_csv_path,
                                              get_processed_csv_path)


def load_config(filename):
    ''' (string) -> dict[str, list]

    Function that initializes the lists that store .csv file paths for 
    every stage of data processing (processing, conversion and calculations), 
    using the global list of raw .csv data paths.

    Returns a dictionnary of the lists of .csv file paths.

    Every raw .csv file (e.g. pairs_20200320010317_20200401010317_1.csv) is 
    associated to a processed (e.g. tri_20200320010317_20200401010317_1.csv), 
    converted (e.g. tri_gridCS_20200320010317_20200401010317_1.csv) and 
    calculated (e.g. calc_20200320010317_20200401010317_1.csv) .csv file.

    The .csv files are stored in the following tree structure, where 
    some_parent_folder can have any name, but will be the same for all
    .csv files resulting from a common raw .csv file :

    ├── 2021_SeaIceDeformations
    |    ...
    |    ├── data
    |    |   ├── 00_grid /etc
    |    |   ├── 01_raw
    |    |   |   └── some_parent_folder
    |    |   |       └── pairs_20200320010317_20200401010317_1.csv
    |    |   ├── 02_processed
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

    processed_csv_paths     = []
    converted_csv_paths     = []
    calculations_csv_paths  = []

    # Iterate through all raw .csv file paths
    for raw_csv_path in lines:
     
        # For each raw .csv file path, find the appropriate path for each subsequent stages of processing
        processed_csv_path    = get_processed_csv_path(raw_csv_path)     # path for processed .csv data file
        converted_csv_path    = get_converted_csv_path(raw_csv_path)     # path for converted .csv data file
        calculations_csv_path = get_calculations_csv_path(raw_csv_path)  # path for calculated .csv data file

        # Add the .csv file paths to the file path lists
        processed_csv_paths.append(processed_csv_path)
        converted_csv_paths.append(converted_csv_path)
        calculations_csv_paths.append(calculations_csv_path)
    
    return {'raw': lines, 'processed': processed_csv_paths, 'converted': converted_csv_paths, 'calculations': calculations_csv_paths}


def get_args():
    p = argparse.ArgumentParser(description='Compute deformations from Sentinel-1 and RCM datasets.')
    p.add_argument('--overwrite', action='store_true', help='When the overwrite argument is specified,           \
                                                            all raw data sets that have already been             \
                                                            processed will be re-processed and the resulting     \
                                                            .csv files for all stages of data processing         \
                                                            will be overwritten')
    
    p.add_argument('--ds_path', default='/home/bdu002/2021_SeaIceDeformation/dataset_paths.txt', 
                                help='Absolute path of .txt file in which are stored the paths of the raw        \
                                        files to process.')

    return p.parse_args()

# Retrieve command line arguments
args = get_args()

# Load .csv file paths for all stages of data processing
csv_paths = load_config(args.ds_path)

if __name__ == '__main__':
    from pprint import pprint
    pprint(load_config('/home/bdu002/.config/SeaIceDeformation.txt'))

    args = get_args()
    print(args.overwrite)
