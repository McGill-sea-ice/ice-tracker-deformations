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



def get_args():
    p = argparse.ArgumentParser(description='Compute deformations from Sentinel-1 and RCM datasets.')
    p.add_argument('-o', '--overwrite', action='store_true', help='When the overwrite argument is specified,           \
                                                            all raw data sets that have already been             \
                                                            processed will be re-processed and the resulting     \
                                                            .csv files for all stages of data processing         \
                                                            will be overwritten')
    
    p.add_argument('-p', '--ds_path', default='/home/bdu002/2021_SeaIceDeformation/dataset_paths.txt', 
                                help='Absolute path of .txt file in which are stored the paths of the raw        \
                                        files to process.')

    p.add_argument('-m', '--method', default='M01', choices=['M00', 'M01'], help='Method to be used when processing  \
                                                                            datasets. M00 refers to the method       \
                                                                            that uses the RIOPS grid, while M01      \
                                                                            refers to the method that uses X/Y       \
                                                                            coordinates from the input datasets.')


    return p.parse_args()

# Retrieve command line arguments
args = get_args()

# Load file paths for all stages of data processing
data_paths = load_config(args.ds_path)

if __name__ == '__main__':
    from pprint import pprint
    pprint(load_config('/home/bdu002/.config/SeaIceDeformation.txt'))

    args = get_args()
    pprint(args)
