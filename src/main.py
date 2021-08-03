'''
Author: Beatrice Duval (bdu002)

---------------------------------
Main script for data processing
--------------------------------

Script that executes all steps towards the calculation of sea-ice deformations and displays the execution time.

'''

import time

import config

if config.args.method == 'M00':
    import M00_d01_delaunay_triangulation as delaunay_triangulation
    import M00_d02_to_grid_coord_system as to_grid_coord_system
    import M00_d03_compute_deformations as compute_deformations

elif config.args.method == 'M01':
    import M01_d01_delaunay_triangulation as delaunay_triangulation
    import M01_d03_compute_deformations as compute_deformations

import visualise_deformation

# Retrieve the starting time
start_time = time.time()


'''
1) Triangulation

Perform a Delaunay triangulation and store the results in a .csv file

'''

print('--- Performing a Delaunay triangulation ---')

delaunay_triangulation.delaunay_triangulation()


'''
2) Conversion (M00 method only)

Convert the triangulation results to a local cartesian grid coordinate system

'''

if config.args.method == 'M00':

    print('--- Converting to local grid coordinate systems ---')

    to_grid_coord_system.to_grid_coord_syst()


'''
3) Calculations

Compute sea-ice deformations rates using the X/Y triangulations results

'''

print('--- Computing sea-ice deformations ---')

compute_deformations.compute_deformations()


'''
4) Visualise Deformations


print('--- Creating sea-ice deformations figures ---')

visualise_deformation.visualise_deformations()
'''

# Display the run time
print("--- %s seconds ---" % (time.time() - start_time))


