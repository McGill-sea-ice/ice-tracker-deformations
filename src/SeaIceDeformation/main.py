'''
Author: Beatrice Duval (bdu002)

---------------------------------
Main script for data processing
--------------------------------

Script that executes all steps towards the calculation of sea-ice deformations and displays the execution time.

'''

import sys
sys.path.insert(0, '/aos/home/dringeisen/code/ice-tracker-deformations/')

import time

import config

import utils_delaunay_triangulation as delaunay_triangulation
import utils_compute_deformations as compute_deformations

import visualise_deformation

from src.SatelliteCoverage import netcdf_tools

# Retrieve the starting time
start_time = time.time()


'''
1) Triangulation

Perform a Delaunay triangulation and store the results in a .csv file

'''

print('--- Performing a Delaunay triangulation ---')

delaunay_triangulation.delaunay_triangulation()


'''
2) Calculations

Compute sea-ice deformations rates using the X/Y triangulations results

'''

print('--- Computing sea-ice deformations ---')
# Output the netcdf dataset
dataset = compute_deformations.compute_deformations()


'''
3) Visualise Deformations
'''

if config.config['Processing_options'].getboolean('visualise'):

    print('--- Creating sea-ice deformation figures ---')

    # Ploting using csv file
    # visualise_deformation.visualise_deformations()

    # Plotting using the netCDF dataset
    netcdf_tools.plot_deformations(data_in=dataset)


# Close netCDF dataset
dataset.close()

# Display the run time
print("--- %s seconds ---" % (time.time() - start_time))


