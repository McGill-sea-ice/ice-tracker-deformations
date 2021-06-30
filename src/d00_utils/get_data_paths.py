'''
Author: Beatrice Duval (bdu002)

------------------------------------------------------------------
Code that provides functions that produce .csv data file paths for   
each stage of data processing given an initial raw .csv file path.
------------------------------------------------------------------

'''

import os

def get_raw_csv_name_dir_parentfolder(raw_csv_path):
    ''' (string) -> tuple(string, string, string)

    Takes as input a raw .csv file's absolute path and returns the 
    name of the file, its directory and its parent folder name 
    (raw_filename, raw_directory, parent_folder).

    Keyword arguments:
    raw_csv_path -- raw .csv file absolute path

    >>> raw_csv_path = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301002108_20200313002108_1.csv'
    >>> print(get_raw_csv_name_dir_parentfolder(raw_csv_path))
    ('pairs_20200301002108_20200313002108_1.csv', '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1', '2020_MarApr_S1')
    '''

    # Get the raw .csv file's name from its path
    raw_filename = os.path.basename(raw_csv_path)
    
    # Get the raw .csv file's directory from its path
    raw_directory = os.path.dirname(raw_csv_path)
    
    # Get the raw .csv file's parent folder
    parent_folder = os.path.basename(raw_directory)

    # Return the raw .csv's name, directory, and parent folder
    return raw_filename, raw_directory, parent_folder


def get_processed_csv_path(raw_csv_path):
    ''' (string) -> (string, string)

    From a raw .csv file's absolute path, this function produces a processing 
    stage data .csv file absolute path. The processed file is to be stored under the
    project's 'data/02_processed' directory, and under a parent directory that 
    has the same name as the raw .csv file's parent directory (i.e. if the raw
    .csv file is stored under a '2020_MarApr_S1' directory, the processed .csv 
    file will be stored under a '2020_MarApr_S1' directory as well).

    Returns the processed file's path and name (processed_csv_path, 
    processed_filename).

    Keyword arguments:
    raw_csv_path -- raw .csv file absolute path

    >>> raw_csv_path = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301002108_20200313002108_1.csv'
    >>> print( get_processed_csv_path(raw_csv_path) )
    ('/home/bdu002/2021_SeaIceDeformation/data/02_processed/2020_MarApr_S1/tri_20200301002108_20200313002108_1.csv', 'tri_20200301002108_20200313002108_1.csv')
    '''

    # Get the raw .csv file's name, directory and parent folder name from its path
    raw_filename, raw_directory, parent_folder = get_raw_csv_name_dir_parentfolder(raw_csv_path)
    
    # Create the processing stage data .csv filename using the raw filename
    processed_filename = 'tri' + raw_filename[5:len(raw_filename)]

    # Get the directory in which the processed .csv file is to be stored
    processed_csv_path = raw_directory + '/../../02_processed/' + parent_folder + '/' +processed_filename

    # Return a normalized processed csv path
    return os.path.normpath(processed_csv_path), processed_filename


def get_converted_csv_path(raw_csv_path):
    ''' (string) -> (string, string)

    From a raw .csv file's absolute path, this function produces a conversion 
    stage data .csv file absolute path. The converted file is to be stored under the
    project's 'data/03_converted' directory, and under a parent directory that 
    has the same name as the raw .csv file's parent directory (i.e. if the raw
    .csv file is stored under a '2020_MarApr_S1' directory, the converted .csv 
    file will be stored under a '2020_MarApr_S1' directory as well).

    Returns the converted file's path and name (converted_csv_path, 
    converted_filename).

    Keyword arguments:
    raw_csv_path -- raw .csv file absolute path

    >>> raw_csv_path = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301002108_20200313002108_1.csv'
    >>> print( get_converted_csv_path(raw_csv_path) )
    ('/home/bdu002/2021_SeaIceDeformation/data/03_converted/2020_MarApr_S1/tri_gridCS_20200301002108_20200313002108_1.csv', 'tri_gridCS_20200301002108_20200313002108_1.csv')
    '''

    # Get the raw .csv file's name, directory and parent folder name from its path
    raw_filename, raw_directory, parent_folder = get_raw_csv_name_dir_parentfolder(raw_csv_path)

    # Create the conversion stage data .csv filename using the raw filename
    converted_filename = 'tri_gridCS' + raw_filename[5:len(raw_filename)]

    # Get the directory in which the converted .csv file is to be stored
    converted_csv_path = raw_directory + '/../../03_converted/' + parent_folder + '/' + converted_filename

    # Return a normalized converted csv path
    return os.path.normpath(converted_csv_path), converted_filename


def get_calculations_csv_path(raw_csv_path):
    ''' (string) -> (string, string)

    From a raw .csv file's absolute path, this function produces a calculations 
    stage data .csv file absolute path. The calculated file is to be stored under the
    project's 'data/04_calculations' directory, and under a parent directory that 
    has the same name as the raw .csv file's parent directory (i.e. if the raw
    .csv file is stored under a '2020_MarApr_S1' directory, the calculated .csv 
    file will be stored under a '2020_MarApr_S1' directory as well).

    Returns the calculations file's path and name (calculations_csv_path, 
    calculations_filename).

    Keyword arguments:
    raw_csv_path -- raw .csv file absolute path

    >>> raw_csv_path = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301002108_20200313002108_1.csv'
    >>> print( get_calculations_csv_path(raw_csv_path) )
    ('/home/bdu002/2021_SeaIceDeformation/data/04_calculations/2020_MarApr_S1/calc_20200301002108_20200313002108_1.csv', 'calc_20200301002108_20200313002108_1.csv')
    '''

    # Get the raw .csv file's name, directory and parent folder name from its path
    raw_filename, raw_directory, parent_folder = get_raw_csv_name_dir_parentfolder(raw_csv_path)

    # Create the calculations stage data .csv filename using the raw filename
    calculations_filename = 'calc' + raw_filename[5:len(raw_filename)]

    # Get the directory in which the calculated .csv file is to be stored
    calculations_csv_path = raw_directory + '/../../04_calculations/' + parent_folder + '/' + calculations_filename

    # Return a normalized calculated csv path
    return os.path.normpath(calculations_csv_path), calculations_filename

