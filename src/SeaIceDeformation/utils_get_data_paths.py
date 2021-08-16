'''
Author: Beatrice Duval (bdu002)

------------------------------------------------------------------
Utils - Get data paths
------------------------------------------------------------------

Code that provides functions that produce data file paths for   
each stage of data processing given an initial raw file name, an 
output folder, an experiment name, and a starting date.

'''


def get_triangulated_csv_path(raw_filename, output_path, exp):
    ''' (str, str, str) -> str

    This function produces a triangulated stage data .csv file absolute path. 
    The triangulated file is to be stored under the output_path/exp/02_triangulated
    directory, and its name will be created from the raw_filename, i.e. we replace 
    the raw_filename's prefix by 'tri'

    Returns the triangulated file's path.

    Keyword arguments: \\
    raw_filename -- raw filename \\
    output_path  -- path in which all output data files are to be stored \\
    exp          -- name of the experiment

    >>> raw_filename = 'pairs_20200301002108_20200313002108_1.csv'
    >>> output_path = '/home/bdu002/outputs'
    >>> exp = '2020_MarApr_S1'
    >>> print( get_triangulated_csv_path(raw_filename, output_path, exp) )
    /home/bdu002/outputs/2020_MarApr_S1/02_triangulated/tri_20200301002108_20200313002108_1.csv
    '''

    # Create the triangulation stage data .csv filename using the raw filename
    tri_filename = 'tri' + raw_filename[-36:]

    # Get the directory in which the triangulated .csv file is to be stored
    tri_path = output_path + '/' + exp + '/02_triangulated/'  + tri_filename

    # Return the triangulated csv path
    return tri_path


def get_converted_csv_path(raw_filename, output_path, exp):
    ''' (str, str, str) -> str

    This function produces a converted stage data .csv file absolute path. 
    The converted file is to be stored under the output_path/exp/03_converted
    directory, and its name will be created from the raw_filename, i.e. we replace 
    the raw_filename's prefix by 'conv'

    Returns the converted file's path.

    Keyword arguments: \\
    raw_filename -- raw filename \\
    output_path  -- path in which all output data files are to be stored \\
    exp          -- name of the experiment

    >>> raw_filename = 'pairs_20200301002108_20200313002108_1.csv'
    >>> output_path = '/home/bdu002/outputs'
    >>> exp = '2020_MarApr_S1'
    >>> print( get_converted_csv_path(raw_filename, output_path, exp) )
    /home/bdu002/outputs/2020_MarApr_S1/03_converted/conv_20200301002108_20200313002108_1.csv
    '''

    # Create the conversion stage data .csv filename using the raw filename
    conv_filename = 'conv' + raw_filename[-36:]

    # Get the directory in which the converted .csv file is to be stored
    conv_path = output_path + '/' + exp + '/03_converted/'  + conv_filename

    # Return the converted csv path
    return conv_path


def get_calculations_csv_path(raw_filename, output_path, exp):
    ''' (str, str, str) -> str

    This function produces a calculations stage data .csv file absolute path. 
    The calculations file is to be stored under the output_path/exp/04_calculations
    directory, and its name will be created from the raw_filename, i.e. we replace 
    the raw_filename's prefix by 'calc'

    Returns the calculations file's path.

    Keyword arguments: \\
    raw_filename -- raw filename \\
    output_path  -- path in which all output data files are to be stored \\
    exp          -- name of the experiment

    >>> raw_filename = 'pairs_20200301002108_20200313002108_1.csv'
    >>> output_path = '/home/bdu002/outputs'
    >>> exp = '2020_MarApr_S1'
    >>> print( get_calculations_csv_path(raw_filename, output_path, exp) )
    /home/bdu002/outputs/2020_MarApr_S1/04_calculations/calc_20200301002108_20200313002108_1.csv
    '''

    # Create the calculations stage data .csv filename using the raw filename
    calc_filename = 'calc' + raw_filename[-36:]

    # Get the directory in which the calculated .csv file is to be stored
    calc_csv_path = output_path + '/' + exp + '/04_calculations/'  + calc_filename

    # Return a normalized calculated csv path
    return calc_csv_path


def get_output_nc_path(output_path, exp, start_year, start_month, start_day):
    ''' (str, str, str, str, str) -> str

    This function produces an output netcdf file absolute path. The output file is to 
    be stored under the output_path/exp/05_output directory, and its name will be created 
    from the input start date.

    The output netcdf file combines the deformation results from all datasets 
    that have been processed simultaneously (listed in config). 

    Returns the ouptut netcdf file's path.

    Keyword arguments: \\
    output_path -- path in which all output data files are to be stored \\
    exp         -- name of the experiment \\
    start_year  -- Common starting year of all datasets that have been processed \\
    start_month -- Common starting month of all datasets that have been processed \\
    start_day   -- Common starting day of all datasets that have been processed

    >>> output_path = '/home/bdu002/outputs'
    >>> exp = '2020_MarApr_S1'
    >>> print( get_output_nc_path(output_path, exp, '2020', '03', '01') )
    '/home/bdu002/outputs/2020_MarApr_S1/05_output/RCMS1SID_20200301_dx.nc'
    '''

    # Create the calculations stage data .csv filename using the raw filename
    output_filename = 'RCMS1SID_' + start_year + start_month + start_day + '_dx.nc'

    # Get the directory in which the calculated .csv file is to be stored
    output_nc_path = output_path + '/' + exp + '/05_output/' + output_filename

    # Return a normalized calculated csv path
    return output_nc_path

