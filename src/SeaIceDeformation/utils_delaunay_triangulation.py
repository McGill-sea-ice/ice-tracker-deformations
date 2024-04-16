'''
Author: Beatrice Duval (bdu002)

----------------------
M01 - Delaunay triangulation
----------------------

Module that performs a Delaunay triangulation on a set of X/Y data points and stores the results in a csv file.

1. Output csv file

    Format:
        no. | vertice_idx1 | vertice_idx2 | vertice_idx3

    Variables:
    - no. is the triangle number
    - vertice_idx1, vertice_idx2, vertice_idx3 are the triangle vertices' indices in the raw file

'''

# Loading from default packages
import csv
import os
from scipy.spatial import Delaunay
from tqdm import tqdm

# Loading from other files
import SeaIceDeformation.utils_grid_coord_system as grid_coord_system
import SeaIceDeformation.utils_load_data as load_data
from config import stb

def delaunay_triangulation(config=None):

    # Retrieve data_paths from config arguments
    dp = config['data_paths']

    # Iterate through all raw and triangulated data file paths listed in config
    empty_files = list(['# Files without data in them'])
    nbfb = 0
    nbfg = 0
    nbpg = 0
    print('--- Performing a Delaunay triangulation ---')
    for raw_path, triangulated_path, calculations_path in zip(tqdm(dp['raw']), dp['triangulated'], dp['calculations']):
        '''
        _________________________________________________________________________________________
        LOAD DATA
        '''

        # If the triangulated file already exists and overwrite (in config) is set to 'no',
        # go to the next iteration.
        # Else, process the raw file and write the triangulated file.
        if os.path.exists(triangulated_path) and not stb(config['Processing_options']['overwrite']):
            continue

        # Load the raw data set. If an error is encountered (no or not enough data points),
        # print the error message and go to the next raw file.
        try:
            raw_data, nbfb, empty_files, nbfg, nbpg = load_data.load_raw( raw_path, nbfb, empty_files, nbfg, nbpg )
            startX = raw_data['startX']
            startY = raw_data['startY']

        except load_data.DataFileError as dfe:
            print(str(dfe) + 'It will not be processed.')
            continue

        '''
        _________________________________________________________________________________________
        PERFORM DELAUNAY TRIANGULATION
        '''

        # Generate a Delaunay triangulation
        startXY = list(zip(startX, startY))
        try:
            tri = Delaunay(startXY)

        except Exception as e:
            # print(e)
            continue

        '''
        _________________________________________________________________________________________
        INITIALIZE AND POPULATE LIST OF DATA ROWS FOR CSV FILE
        '''

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

        '''
        _________________________________________________________________________________________
        WRITE RESULTS
        '''

        # Check if we have results to write
        if len(row_list) > 1:

            # Create a directory to store the triangulated csv path if it does not exist already
            os.makedirs(os.path.dirname(triangulated_path), exist_ok=True)

            # Write the results in the triangulated file path
            with open(triangulated_path, 'w', encoding='UTF8', newline='') as f:
                writer = csv.writer(f)

                # Write the data rows to the csv file
                writer.writerows(row_list)

        # Else, and if we are overwriting, delete the existing triangulation file and subsequent calculations file
        elif stb(config['Processing_options']['overwrite']):

            try:
                os.remove(triangulated_path)
            except OSError:
                pass

            try:
                os.remove(calculations_path)
            except OSError:
                pass

    fname_b = dp["output_folder"]+"/empty_files.txt"
    with open(fname_b, 'w') as file:
        for row in empty_files:
            file.write(row+'\n')

    print("The number of read files is ", nbfg, " with a total of ", nbpg, "data points.")
    print("The number of files with no or not enough data is ", nbfb, ". These files were not processed")
    print("The list of not processed files is save in file ", fname_b)

    return None

if __name__ == "__main__":
    delaunay_triangulation()

