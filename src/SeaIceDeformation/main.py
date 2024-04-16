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
date_pairs = divide_intervals(config=config)

# Iterating over each (daily) interval
for i in tqdm(range(len(date_pairs)), position=0, leave=False):

    #Set the start and end time according to the given date
    datepairs = date_pairs[i]
    config['Date_options']['start_year'] = datepairs[0].strftime("%Y")
    config['Date_options']['start_month'] = datepairs[0].strftime("%m")
    config['Date_options']['start_day'] = datepairs[0].strftime("%d")
    config['Date_options']['end_year'] = datepairs[1].strftime("%Y")
    config['Date_options']['end_month'] = datepairs[1].strftime("%m")
    config['Date_options']['end_day'] = datepairs[1].strftime("%d")

    # Loads sea ice motion data for the given interval
    raw_paths = filter_data(config=config)
    config['raw_paths'] = raw_paths

    # Get the paths to which data files of all stages of data processing will be stored
    data_paths = get_datapaths(config=config)
    config['data_paths'] = data_paths

    #Triangulate and store results to temp .csv files (to be improved)
    delaunay_triangulation(config=config)

    #Compute sea-ice deformations rates and stack in daily netcdf output.
    dataset = compute_deformations(config=config)

    # Close netCDF dataset
    dataset.close()

# Display the run time
print("--- %s seconds ---" % (time.time() - start_time))


