'''
Author: Beatrice Duval (bdu002)

---------------------------------------------------------------------------------
Script that executes all steps towards the calculation of sea-ice deformations.

All raw .csv files listed in config.py are processed. See config.py file to modify 
parameters.

The execution time is displayed.
---------------------------------------------------------------------------------

'''

import time
start_time = time.time()


'''
1) Initialize global variables 

Initialize .csv file path lists for all data processing stages (processed, 
converted and calculated .csv's).

'''

import src.config
src.config.init()


'''
2) Processing 

Perform a Delaunay triangulation and store the results in a csv file

'''

print('--- Performing a Delaunay triangulation ---')

import src.d02_processing.delaunay_triangulation


'''
3) Conversion 

Convert the triangulation results to a local cartesian grid coordinate system

'''

print('--- Converting to local grid coordinate systems ---')

import src.d03_conversion.to_grid_coord_system


'''
4) Calculation

Compute sea-ice deformations rates using the converted triangulations results

'''

print('--- Computing sea-ice deformations ---')

import src.d04_calculating.compute_deformations


# Display the run time
print("--- %s seconds ---" % (time.time() - start_time))