'''
Author: Beatrice Duval (bdu002)

--------------------------------------------------------------------
Code that performs a Delaunay triangulation on a set of X/Y data 
points and stores the results in a csv file.
--------------------------------------------------------------------

Output csv file format:
    ________________________________________________
    no. | vertice_idx1 | vertice_idx2 | vertice_idx3

where:
    - no. is the triangle number;
    - vertice_idx1, vertice_idx2, vertice_idx3 are the triangle vertices' indices in the raw file

'''

import csv
import os

from scipy.spatial import Delaunay

import config
import utils_load_data as load_data


def delaunay_triangulation():

    # Iterate through all raw data file paths listed in config
    for raw_path, triangulated_path in zip(config.data_paths['raw'], config.data_paths['triangulated']):
        # If the triangulated file already exists and overwrite (in config) is set to false,
        # go to the next iteration.
        # ELse, process the raw file and write/overwrite the triangulated file.
        if os.path.exists(triangulated_path) and not config.args.overwrite:
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

