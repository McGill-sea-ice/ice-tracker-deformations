'''
Author: Beatrice Duval (bdu002)

----------------------------------------------------------------------------
Code that converts the triangles' lat,lon coordinates from the processed 
csv file to x,y coordinates in a local cartesian grid coordinate system (CS). 

These results are then stored in a converted csv file.
----------------------------------------------------------------------------

Output csv file format:
    _________________________________________________________________________________________________
    no. | sX1 | sX2 | sX3 | sY1 | sY2 | sY3 | eX1 | eX2 | eX3 | eY1 | eY2 | eY3 | tracer_j | tracer_i

where:
    - no. is the triangle number;
    - sX1, sX2, sX3 are the starting x coordinates for the 1st, 2nd and 3rd vertices;
    - sY1, sY2, sY3 are the ending y coordinates for the 1st, 2nd and 3rd vertices;
    - eX1, eX2, eX3, eY1, eY2, eY3 are the ending x and y coordinates;
    - tracer_j and tracer_i are the indices of the tracer points defining local CS.

'''

import csv
import os
from math import nan
from src.d00_utils.grid_coord_system import (define_gridCS, find_aeqdTriCenter,
                                             find_nearestGridTracerPt,
                                             get_xy_gridCS, get_tri_angles)
from src.d00_utils.load_csv import load_raw_csv, load_processed_csv

from src.d01_grid.load00_grid import X_aeqd, Y_aeqd, fLAT, fLON
import src.config

# Iterate through all raw, processed and converted .csv files listed in config.py
for raw_csv_path, processed_csv_path, converted_csv_path in zip(src.config.raw_csv_paths, src.config.processed_csv_paths,  src.config.converted_csv_paths):

    # If the converted file already exists and overwrite (in config.py) is set to false,
    # go to the next iteration.
    # ELse, process the raw and processed files and write/overwrite the converted file.
    if os.path.exists(converted_csv_path) and not src.config.overwrite:
        continue
    
    if not( os.path.exists(raw_csv_path) and os.path.exists(processed_csv_path) ):
        continue

    # Load the raw and the processed data sets
    sLat, sLon, eLat, eLon = load_raw_csv( raw_csv_path )

    _, sX1_aeqd, sX2_aeqd, sX3_aeqd, \
       sY1_aeqd, sY2_aeqd, sY3_aeqd, \
       vertice_idx1, vertice_idx2, vertice_idx3 = load_processed_csv( processed_csv_path )

       
    # Create a header and a list of data rows that will be used to create the output csv file
    header = ['no.', 'sX1', 'sX2', 'sX3',             
                    'sY1', 'sY2', 'sY3', 
                    'eX1', 'eX2', 'eX3', 
                    'eY1', 'eY2', 'eY3', 
                    'tracer_j', 'tracer_i']
    
    row_list = [header]

    # Iterate through every triangle in the processed csv file 
    for n in range(len(sX1_aeqd)):

        # Compute the triangle's starting central point in the aeqd transform
        x_aeqdCenter, y_aeqdCenter = find_aeqdTriCenter( sX1_aeqd[n], sX2_aeqd[n], sX3_aeqd[n], sY1_aeqd[n], sY2_aeqd[n], sY3_aeqd[n])

        # Find the matrix indices of the nearest grid tracer point relative to 
        # the triangle center point
        j, i = find_nearestGridTracerPt( X_aeqd, Y_aeqd, (x_aeqdCenter, y_aeqdCenter) )

        # The grid coordinate system (CS) is defined by the speed points 
        # (A-D points shown below) associated to and around the tracer point:
        #
        #                       B(0,b)      A
        #
        #                            (j,i)
        #
        #                       C(0,0)    D(d,0)

        # Find the (lat, lon) position of B, C and D and the distances b and d,
        # and store these as a tuple (gridCS)
        gridCS = define_gridCS(fLAT, fLON, (j,i))    # = (B, C, D, b, d)

        # Retrieve each vertex's index in the raw csv file
        v1_index = vertice_idx1[n]
        v2_index = vertice_idx2[n]
        v3_index = vertice_idx3[n]

        # Find the x,y starting and ending positions of the triangle vertices in 
        # the local CS
        sX1, sY1 = get_xy_gridCS(gridCS, (sLat[v1_index], sLon[v1_index]) ) # Starting positions
        sX2, sY2 = get_xy_gridCS(gridCS, (sLat[v2_index], sLon[v2_index]) )
        sX3, sY3 = get_xy_gridCS(gridCS, (sLat[v3_index], sLon[v3_index]) )

        eX1, eY1 = get_xy_gridCS(gridCS, (eLat[v1_index], eLon[v1_index]) ) # Ending positions
        eX2, eY2 = get_xy_gridCS(gridCS, (eLat[v2_index], eLon[v2_index]) )
        eX3, eY3 = get_xy_gridCS(gridCS, (eLat[v3_index], eLon[v3_index]) )

        # Compute the initial triangle angles
        angles = get_tri_angles( (sX1, sY1), (sX2, sY2), (sX3, sY3))

        # If one of the angles is inferior to 10 degrees (= 0.175 rad), reject it
        if any(angle < 0.175 for angle in angles):
            sX1 = sX2 = sX3 = sY1 = sY2 = sY3 = eX1 = eX2 = eX3 = eY1 = eY2 = eY3 = j = i = nan

        # Add the data row corresponding to the current triangle to the list of data rows
        row_list.append( [ n, sX1, sX2, sX3, sY1, sY2, sY3, eX1, eX2, eX3, eY1, eY2, eY3, j, i] )


    #--------------------Write the results to a csv file---------------------------------

    # Write the results in the converted_csv_path file path
    with open(converted_csv_path, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)

        # Write the data rows to the csv file
        writer.writerows(row_list)
