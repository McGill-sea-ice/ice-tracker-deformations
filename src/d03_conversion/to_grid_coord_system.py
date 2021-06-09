'''
Author: Beatrice Duval (bdu002)

----------------------------------------------------------------------------
Code that converts the triangles' lat,lon coordinates from the processed 
csv file to x,y coordinates in a local cartesian grid coordinate system (CS). 

These results are then stored in a converted csv file.
----------------------------------------------------------------------------

Output csv file format:
    no., sX1, sX2, sX3, sY1, sY2, sY3, eX1, eX2, eX3, eY1, eY2, eY3, tracer_j, tracer_i

where:
    - no. is the triangle number;
    - sX1, sX2, sX3 are the starting x coordinates for the 1st, 2nd and 3rd vertices;
    - sY1, sY2, sY3 are the ending y coordinates for the 1st, 2nd and 3rd vertices;
    - eX1, eX2, eX3, eY1, eY2, eY3 are the ending x and y coordinates;
    - tracer_j and tracer_i are the indices of the tracer points defining local CS.


'''

from src.d01_data.load_processed_csv import sLat1, sLat2, sLat3, sLon1, sLon2, sLon3, eLat1, eLat2, eLat3, eLon1, eLon2, eLon3, sX1_aeqd, sX2_aeqd, sX3_aeqd, sY1_aeqd, sY2_aeqd, sY3_aeqd
from src.d01_data.load_grid import X_aeqd, Y_aeqd, fLAT, fLON
from src.d00_utils.grid_coord_system import find_nearestGridTracerPt, find_aeqdTriCenter, define_gridCS, get_xy_gridCS
from src.d01_data.data_paths import converted_csv_path
import csv

# Create a header and a list of data rows that will be used to create the output csv file
header = ['no.', 'sX1', 'sX2', 'sX3',             
                 'sY1', 'sY2', 'sY3', 
                 'eX1', 'eX2', 'eX3', 
                 'eY1', 'eY2', 'eY3', 
                 'tracer_j', 'tracer_i']
row_list = [header]

# Iterate through every triangle in the processed csv file 
for n in range(len(sLat1)):

    # Compute the triangle's starting central point in the aeqd transform
    xy_aeqdCenter = find_aeqdTriCenter( sX1_aeqd[n], sX2_aeqd[n], sX3_aeqd[n], sY1_aeqd[n], sY2_aeqd[n], sY3_aeqd[n])

    # Find the matrix indices of the nearest grid tracer point relative to 
    # the triangle center point
    j, i = find_nearestGridTracerPt( X_aeqd, Y_aeqd, xy_aeqdCenter )

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

    # Find the x,y starting and ending positions of the triangle vertices in 
    # the local CS
    sX1, sY1 = get_xy_gridCS(gridCS, (sLat1[n], sLon1[n]) ) # Starting positions
    sX2, sY2 = get_xy_gridCS(gridCS, (sLat2[n], sLon2[n]) )
    sX3, sY3 = get_xy_gridCS(gridCS, (sLat3[n], sLon3[n]) )

    eX1, eY1 = get_xy_gridCS(gridCS, (eLat1[n], eLon1[n]) ) # Ending positions
    eX2, eY2 = get_xy_gridCS(gridCS, (eLat2[n], eLon2[n]) )
    eX3, eY3 = get_xy_gridCS(gridCS, (eLat3[n], eLon3[n]) )

    # Add the data row corresponding to the current triangle to the list of data rows
    row_list.append( [ n, sX1, sX2, sX3, sY1, sY2, sY3, eX1, eX2, eX3, eY1, eY2, eY3, j, i ] )


#--------------------Write the results to a csv file---------------------------------

# Write the results in the converted_csv_path file path
with open(converted_csv_path, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # Write the data rows to the csv file
    writer.writerows(row_list)