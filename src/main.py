'''
Author: Beatrice Duval (bdu002)

---------------------------------------------------------------------------------
Script that executes all steps towards the calculation of sea-ice deformations.
The execution time is displayed.
---------------------------------------------------------------------------------

'''

import time
start_time = time.time()


'''
1) Processing 

Perform a Delaunay triangulation and store the lat,lon results in a csv file

'''

import src.d02_processing.delaunay_triangulation

'''
2) Conversion 

Convert lat,lon triangulation results to x,y coordinates in a local cartesian 
grid coordinate system

'''

import src.d03_conversion.to_grid_coord_system


'''
3) Calculation

Compute sea-ice deformations rates using the converted triangulations results

'''


print("--- %s seconds ---" % (time.time() - start_time))