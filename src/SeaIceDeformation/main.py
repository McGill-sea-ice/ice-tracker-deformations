'''
Author: Beatrice Duval (bdu002)

---------------------------------
Main script for data processing
--------------------------------

Script that executes all steps towards the calculation of sea-ice deformations and displays the execution time.

'''

# Loading from default packages
import os
import sys
parent = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,parent)
import time
import datetime
from tqdm import tqdm
# Loading from other files
from config import get_config_args, filter_data
from config import divide_intervals, get_datapaths
from utils_delaunay_triangulation import stb, delaunay_triangulation
from utils_compute_deformations import compute_deformations
from visualise_deformation import visualise_deformations
from SatelliteCoverage.netcdf_tools import plot_deformations

# Retrieve the starting time
start_time = time.time()

'''
0) Configuration extraction
_______________________________________________________________________
PERFORM CONFIGURATION
'''
 # Retrieve configuration arguments from namelist.ini
config = get_config_args()

raw_paths = filter_data(config=config)
#config['raw_paths'] = raw_paths
print(raw_paths)
date_pairs = divide_intervals(config=config)
print(date_pairs)
# Iterating over each interval
for i in tqdm(range(len(date_pairs)), position=0, leave=True):
    # Loads data and converts to x/y for each interval
    datepairs = date_pairs[i]
    config['Date_options']['start_year'] = datepairs[0].strftime("%Y")
    config['Date_options']['start_month'] = datepairs[0].strftime("%m")
    config['Date_options']['start_day'] = datepairs[0].strftime("%d")
    config['Date_options']['end_year'] = datepairs[1].strftime("%Y")
    config['Date_options']['end_month'] = datepairs[1].strftime("%m")
    config['Date_options']['end_day'] = datepairs[1].strftime("%d")
    raw_paths = filter_data(config=config)
    config['raw_paths'] = raw_paths
    # Get the paths to which data files of all stages of data processing will be stored
    data_paths = get_datapaths(config=config)
    config['data_paths'] = data_paths


    '''
    1) Triangulation

    Perform a Delaunay triangulation and store the results in a .csv file

    '''

    delaunay_triangulation(config=config)


    '''
    2) Calculations

    Compute sea-ice deformations rates using the X/Y triangulations results

    '''

    dataset = compute_deformations(config=config)


    '''
    3) Visualise Deformations
    '''

    if config['Processing_options']['visualise']:

        # Ploting using csv file
        # visualise_deformations(config=config)

        # Plotting using the netCDF dataset
        plot_deformations(data_in=dataset, config=config)


    # Close netCDF dataset
    dataset.close()

# Display the run time
print("--- %s seconds ---" % (time.time() - start_time))


