'''
Author: Beatrice Duval (bdu002)

----------------------
M01 - Delaunay triangulation
----------------------

Code that performs a Delaunay triangulation on a set of X/Y data points and stores the results in a csv file.

1. Output csv file
    
    Format:
        no. | vertice_idx1 | vertice_idx2 | vertice_idx3

    Variables:
    - no. is the triangle number
    - vertice_idx1, vertice_idx2, vertice_idx3 are the triangle vertices' indices in the raw file

'''

import csv
import os

from scipy.spatial import Delaunay

import config
import utils_load_data as load_data
import utils_grid_coord_system as grid_coord_system

def delaunay_triangulation():

    # Iterate through all raw data file paths listed in config
    for raw_path, triangulated_path in zip(config.data_paths['raw'], config.data_paths['triangulated']):
        # If the triangulated file already exists and overwrite (in config) is set to false,
        # go to the next iteration.
        # ELse, process the raw file and write/overwrite the triangulated file.
        if os.path.exists(triangulated_path) and not config.config['Processing_options'].getboolean('overwrite'):
            continue

        # Load the raw data set. If an error is encountered (no or not enough data points), 
        # print the error message and go to the next raw .csv file.
        try:
            raw_data = load_data.load_raw( raw_path )
            startX = raw_data['startX']
            startY = raw_data['startY']
        
        except load_data.DataFileError as dfe:
            print(str(dfe) + 'It will not be processed.')
            continue

        # Generate a Delaunay triangulation
        startXY = list(zip(startX, startY))
        tri = Delaunay(startXY)

        # Create a header and a list of data rows that will be used to create the output csv file
        header = ['no.', 'vertice_idx1', 'vertice_idx2',  'vertice_idx3']
        row_list = [header]

        # Iterate through all triangles
        for n in range(len(tri.simplices)):
                
            # Find the index of all 3 data points that form the current triangle
            vertice_idx1 = tri.simplices[n][0]
            vertice_idx2 = tri.simplices[n][1]
            vertice_idx3 = tri.simplices[n][2]

            # Find the starting X/Y position of each triangle vertex
            sX1 = startX[vertice_idx1]
            sY1 = startY[vertice_idx1]
            sX2 = startX[vertice_idx2]
            sY2 = startY[vertice_idx2]
            sX3 = startX[vertice_idx3]
            sY3 = startY[vertice_idx3]

            # Compute the current triangle's angles
            angles = grid_coord_system.get_tri_angles( (sX1, sY1), (sX2, sY2), (sX3, sY3))

            # Keep the triangle if and only if none of its angles are inferior to 10 degrees (= 0.175 rad)
            if not any(angle < 0.175 for angle in angles):
                
                # Add the data row corresponding to the current triangle to the list of data rows
                row_list.append([ n, vertice_idx1, vertice_idx2,  vertice_idx3] )


        #--------------------Write the results to a csv file---------------------------------

        # Create a directory to store the triangulated csv path if it does not exist already
        os.makedirs(os.path.dirname(triangulated_path), exist_ok=True)

        # Write the results in the triangulated file path
        with open(triangulated_path, 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)

            # Write the data rows to the csv file
            writer.writerows(row_list)


if __name__ == "__main__":
    delaunay_triangulation()

