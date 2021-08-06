'''
Author: Beatrice Duval (bdu002)

------------------------------------------------------------------
Utils - Get data paths
------------------------------------------------------------------

Code that provides functions that produce data file paths for   
each stage of data processing given an initial raw file path.

'''

import os

def get_raw_csv_name_dir_parentfolder(raw_path):
    ''' (str) -> tuple(str, str, str)

    Takes as input a raw file's absolute path and returns the 
    name of the file, its directory and its parent folder name 
    (raw_filename, raw_directory, parent_folder).

    Keyword arguments: \\
    raw_path -- raw file absolute path

    >>> raw_csv_path = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301002108_20200313002108_1.csv'
    >>> print(get_raw_csv_name_dir_parentfolder(raw_csv_path))
    'pairs_20200301002108_20200313002108_1.csv', '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1', '2020_MarApr_S1'
    '''

    # Get the raw file's name from its path
    raw_filename = os.path.basename(raw_path)
    
    # Get the raw file's directory from its path
    raw_directory = os.path.dirname(raw_path)
    
    # Get the raw file's parent folder
    parent_folder = os.path.basename(raw_directory)

    # Return the raw file's name, directory, and parent folder
    return raw_filename, raw_directory, parent_folder


def get_processed_csv_path(raw_path):
    ''' (str) -> str

    From a raw file's absolute path, this function produces a processing 
    stage data .csv file absolute path. The processed file is to be stored under the
    project's 'data/02_processed' directory, and under a parent directory that 
    has the same name as the raw file's parent directory (i.e. if the raw
    .csv file is stored under a '2020_MarApr_S1' directory, the processed .csv 
    file will be stored under a '2020_MarApr_S1' directory as well).

    Returns the processed file's path.

    Keyword arguments: \\
    raw_path -- raw file absolute path

    >>> raw_csv_path = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301002108_20200313002108_1.csv'
    >>> print( get_processed_csv_path(raw_csv_path) )
    '/home/bdu002/2021_SeaIceDeformation/data/02_processed/2020_MarApr_S1/tri_20200301002108_20200313002108_1.csv'
    '''

    # Get the raw file's name, directory and parent folder name from its path
    raw_filename, raw_directory, parent_folder = get_raw_csv_name_dir_parentfolder(raw_path)
    
    # Create the processing stage data .csv filename using the raw filename
    processed_filename = 'tri' + raw_filename[5:len(raw_filename)]

    # Get the directory in which the processed .csv file is to be stored
    processed_csv_path = raw_directory + '/../../02_triangulated/' + parent_folder + '/' +processed_filename

    # Return a normalized processed csv path
    return os.path.normpath(processed_csv_path)


def get_converted_csv_path(raw_path):
    ''' (str) -> str

    From a raw file's absolute path, this function produces a conversion 
    stage data .csv file absolute path. The converted file is to be stored under the
    project's 'data/03_converted' directory, and under a parent directory that 
    has the same name as the raw file's parent directory (i.e. if the raw
    .csv file is stored under a '2020_MarApr_S1' directory, the converted .csv 
    file will be stored under a '2020_MarApr_S1' directory as well).

    Returns the converted file's path.

    Keyword arguments: \\
    raw_path -- raw file absolute path

    >>> raw_csv_path = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301002108_20200313002108_1.csv'
    >>> print( get_converted_csv_path(raw_csv_path) )
    '/home/bdu002/2021_SeaIceDeformation/data/03_converted/2020_MarApr_S1/tri_gridCS_20200301002108_20200313002108_1.csv'
    '''

    # Get the raw file's name, directory and parent folder name from its path
    raw_filename, raw_directory, parent_folder = get_raw_csv_name_dir_parentfolder(raw_path)

    # Create the conversion stage data .csv filename using the raw filename
    converted_filename = 'tri_gridCS' + raw_filename[5:len(raw_filename)]

    # Get the directory in which the converted .csv file is to be stored
    converted_csv_path = raw_directory + '/../../03_converted/' + parent_folder + '/' + converted_filename

    # Return a normalized converted csv path
    return os.path.normpath(converted_csv_path)


def get_calculations_csv_path(raw_path):
    ''' (str) -> str

    From a raw file's absolute path, this function produces a calculations 
    stage data .csv file absolute path. The calculated file is to be stored under the
    project's 'data/04_calculations' directory, and under a parent directory that 
    has the same name as the raw file's parent directory (i.e. if the raw
    .csv file is stored under a '2020_MarApr_S1' directory, the calculated .csv 
    file will be stored under a '2020_MarApr_S1' directory as well).

    Returns the calculations file's path.

    Keyword arguments: \\
    raw_path -- raw file absolute path

    >>> raw_csv_path = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301002108_20200313002108_1.csv'
    >>> print( get_calculations_csv_path(raw_csv_path) )
    '/home/bdu002/2021_SeaIceDeformation/data/04_calculations/2020_MarApr_S1/calc_20200301002108_20200313002108_1.csv'
    '''

    # Get the raw file's name, directory and parent folder name from its path
    raw_filename, raw_directory, parent_folder = get_raw_csv_name_dir_parentfolder(raw_path)

    # Create the calculations stage data .csv filename using the raw filename
    calculations_filename = 'calc' + raw_filename[5:len(raw_filename)]

    # Get the directory in which the calculated .csv file is to be stored
    calculations_csv_path = raw_directory + '/../../04_calculations/' + parent_folder + '/' + calculations_filename

    # Return a normalized calculated csv path
    return os.path.normpath(calculations_csv_path)


def get_output_nc_path(raw_path, start_year, start_month, start_day):
    ''' (str, str, str, str) -> str

    From a raw file's absolute path and a start time, this function 
    produces an output netcdf file absolute path. The output file is to 
    be stored under the project's 'data/05_output' directory, and under a
    parent directory that has the same name as the raw file's parent 
    directory (i.e. if the raw files are stored under a '2020_MarApr_S1' 
    directory, the output netcdf file will be stored under a '2020_MarApr_S1' 
    directory as well).

    The output netcdf file combines the deformation results from all datasets 
    that have been processed simultaneously (listed in config). 

    Returns the ouptut netcdf file's path.

    Keyword arguments: \\
    raw_path    -- one of the raw file absolute paths \\
    start_year  -- Common starting year of all datasets that have been processed \\
    start_month -- Common starting month of all datasets that have been processed \\
    start_day   -- Common starting day of all datasets that have been processed

    >>> raw_csv_path = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1/pairs_20200301002108_20200313002108_1.csv'
    >>> print( get_output_nc_path(raw_csv_path, '2020', '03', '01') )
    '/home/bdu002/2021_SeaIceDeformation/data/05_output/2020_MarApr_S1/RCMS1SID_20200301_dx.nc'
    '''

    # Get the raw file's name, directory and parent folder name from its path
    _, raw_directory, parent_folder = get_raw_csv_name_dir_parentfolder(raw_path)

    # Create the calculations stage data .csv filename using the raw filename
    output_filename = 'RCMS1SID_' + start_year + start_month + start_day + '_dx.nc'

    # Get the directory in which the calculated .csv file is to be stored
    output_nc_path = raw_directory + '/../../05_output/' + parent_folder + '/' + output_filename

    # Return a normalized calculated csv path
    return os.path.normpath(output_nc_path)

