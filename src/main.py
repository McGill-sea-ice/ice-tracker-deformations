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

import d02_delaunay_triangulation 
import d03_to_grid_coord_system  
import d04_compute_deformations 


# Retrieve the starting time
start_time = time.time()


'''
2) Processing 

Perform a Delaunay triangulation and store the results in a csv file

'''

print('--- Performing a Delaunay triangulation ---')
d02_delaunay_triangulation.delaunay_triangulation()


'''
3) Conversion 

Convert the triangulation results to a local cartesian grid coordinate system

'''

print('--- Converting to local grid coordinate systems ---')
d03_to_grid_coord_system.to_grid_coord_syst()


'''
4) Calculation

Compute sea-ice deformations rates using the converted triangulations results

'''

print('--- Computing sea-ice deformations ---')
d04_compute_deformations.compute_deformations()


# Display the run time
print("--- %s seconds ---" % (time.time() - start_time))
