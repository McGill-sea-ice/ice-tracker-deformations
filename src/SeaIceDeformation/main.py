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
# Loading from other files
from config import get_config
from utils_delaunay_triangulation import stb, delaunay_triangulation
from utils_compute_deformations import compute_deformations
from visualise_deformation import visualise_deformations
from SatelliteCoverage.netcdf_tools import plot_deformations

# Retrieve the starting time
start_time = time.time()

'''
0) Configuration extraction

Read the namelist.ini file

'''
config = get_config()

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


