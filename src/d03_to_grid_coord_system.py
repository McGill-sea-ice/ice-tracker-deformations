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
from math import nan, sqrt
import utils_grid_coord_system as grid_cs
import utils_load_grid as load_grid
import utils_load_csv as load_csv
import config

def to_grid_coord_syst():

    # Iterate through all raw, processed and converted .csv files listed in config.py
    for raw_csv_path, processed_csv_path, converted_csv_path in zip(config.csv_paths['raw'], config.csv_paths['processed'], config.csv_paths['converted']):

        # If the converted file already exists and overwrite (in config.py) is set to false,
        # go to the next iteration.
        # ELse, process the raw and processed files and write/overwrite the converted file.
        if os.path.exists(converted_csv_path) and not config.args.overwrite:
            continue
        
        # Check if the raw and processed .csv files exist
        if not( os.path.exists(raw_csv_path) or os.path.exists(processed_csv_path) ):
            continue
        
        # Load the raw data set. If an error is encountered (no or not enough data points), 
        # print the error message and go to the next raw .csv file.
        try:
            sLat, sLon, eLat, eLon = load_csv.load_raw_csv( raw_csv_path )
        except load_csv.DataFileError as dfe:
            print(dfe + 'It will not be processed.')
            continue

        # Load the processed data set
        _, sX1_aeqd, sX2_aeqd, sX3_aeqd, \
        sY1_aeqd, sY2_aeqd, sY3_aeqd, \
        vertice_idx1, vertice_idx2, vertice_idx3 = load_csv.load_processed_csv( processed_csv_path )

        # Load the RIOPS grid
        grid = load_grid.load_grid()

        # Create a header and a list of data rows that will be used to create the output csv file
        header = ['no.', 'sX1', 'sX2', 'sX3',             
                        'sY1', 'sY2', 'sY3', 
                        'eX1', 'eX2', 'eX3', 
                        'eY1', 'eY2', 'eY3', 
                        'tracer_j', 'tracer_i']
        
        row_list = [header]

        go_to_next_dataset = False

        # Iterate through every triangle in the processed csv file 
        for n in range(len(sX1_aeqd)):

            # Compute the triangle's starting central point in the aeqd transform
            x_aeqdCenter, y_aeqdCenter = grid_cs.find_aeqdTriCenter( sX1_aeqd[n], sX2_aeqd[n], sX3_aeqd[n], sY1_aeqd[n], sY2_aeqd[n], sY3_aeqd[n])

            # Find the matrix indices of the nearest grid tracer point relative to 
            # the triangle center point
            j, i = grid_cs.find_nearestGridTracerPt( grid['X_aeqd'], grid['Y_aeqd'], (x_aeqdCenter, y_aeqdCenter) )

            # Compute the distance (in meters) between the center point and the nearest grid tracer point
            dist = sqrt((grid['X_aeqd'][j,i] - x_aeqdCenter)**2 + (grid['Y_aeqd'][j,i] - y_aeqdCenter)**2)

            # If the distance is greater than 10 km, we are not in the region of interest.
            if dist > 10000:
                # Delete the processed .csv file and go to the next iteration, i.e. do not process further this data set
                os.remove(processed_csv_path)
                print( os.path.basename(raw_csv_path) + ' is not in the region of interest. It will not be processed.')
                go_to_next_dataset = True
                break
            
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
            gridCS = grid_cs.define_gridCS(grid['fLAT'], grid['fLON'], (j,i))    # = (B, C, D, b, d)

            # Retrieve each vertex's index in the raw csv file
            v1_index = vertice_idx1[n]
            v2_index = vertice_idx2[n]
            v3_index = vertice_idx3[n]

            # Find the x,y starting and ending positions of the triangle vertices in 
            # the local CS
            sX1, sY1 = grid_cs.get_xy_gridCS(gridCS, (sLat[v1_index], sLon[v1_index]) ) # Starting positions
            sX2, sY2 = grid_cs.get_xy_gridCS(gridCS, (sLat[v2_index], sLon[v2_index]) )
            sX3, sY3 = grid_cs.get_xy_gridCS(gridCS, (sLat[v3_index], sLon[v3_index]) )

            eX1, eY1 = grid_cs.get_xy_gridCS(gridCS, (eLat[v1_index], eLon[v1_index]) ) # Ending positions
            eX2, eY2 = grid_cs.get_xy_gridCS(gridCS, (eLat[v2_index], eLon[v2_index]) )
            eX3, eY3 = grid_cs.get_xy_gridCS(gridCS, (eLat[v3_index], eLon[v3_index]) )

            # Compute the initial triangle angles
            angles = grid_cs.get_tri_angles( (sX1, sY1), (sX2, sY2), (sX3, sY3))

            # If one of the angles is inferior to 10 degrees (= 0.175 rad), or
            # if the tracer point is on land (distance to land = 0), reject it
            if any(angle < 0.175 for angle in angles) or grid['DIST'][j,i] == 0.0:
                sX1 = sX2 = sX3 = sY1 = sY2 = sY3 = eX1 = eX2 = eX3 = eY1 = eY2 = eY3 = j = i = nan

            # Add the data row corresponding to the current triangle to the list of data rows
            row_list.append( [ n, sX1, sX2, sX3, sY1, sY2, sY3, eX1, eX2, eX3, eY1, eY2, eY3, j, i] )

        if go_to_next_dataset:
            continue

        #--------------------Write the results to a csv file---------------------------------
        
        # Create a directory to store the converted csv path if it does not exist already
        os.makedirs(os.path.dirname(converted_csv_path), exist_ok=True)

        # Write the results in the converted_csv_path file path
        with open(converted_csv_path, 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)

            # Write the data rows to the csv file
            writer.writerows(row_list)


if __name__ == '__main__':
    to_grid_coord_syst()