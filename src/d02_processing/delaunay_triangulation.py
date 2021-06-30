'''
Author: Beatrice Duval (bdu002)

--------------------------------------------------------------------
Code that performs a Delaunay triangulation on a set of data points 
and stores the results in a csv file.
--------------------------------------------------------------------

Output csv file format:
    ___________________________________________________________________________________________________
    no. | sX1_aeqd | sX2_aeqd | sX3_aeqd | sY1_aeqd | sY2_aeqd | sY3_aeqd | vertice_idx1 | vertice_idx2

where:
    - no. is the triangle number;
    - sX1_aeqd, sX2_aeqd, sX3_aeqd, sY1_aeqd, sY2_aeqd, sY3_aeqd are the x,y coordinates 
        of the starting vertices in the Azimuthal Equidistant (aeqd) transform.
    - vertice_idx1, vertice_idx2, vertice_idx3 are the triangle vertices' indices in the raw .csv file
'''

import csv
import os

import src.config
from pyproj import Proj
from scipy.spatial import Delaunay
from src.d00_utils.load_csv import load_raw_csv

# Create a Azimuthal Equidistant (aeqd) transform projection object
p = Proj(proj='aeqd', ellps='WGS84', preserve_units=False)

# Iterate through all raw and processed .csv file paths listed in config.py
for raw_csv_path, processed_csv_path in zip(src.config.raw_csv_paths, src.config.processed_csv_paths):

    # If the processed file already exists and overwrite (in config.py) is set to false,
    # go to the next iteration.
    # ELse, process the raw file and write/overwrite the processed file.
    if os.path.exists(processed_csv_path) and not src.config.overwrite:
        continue

    # Load the raw data set. If an error is encountered (no or not enough data points), 
    # go to the next raw .csv file.
    try:
        sLat, sLon, eLat, eLon = load_raw_csv( raw_csv_path )
    except:
        continue

    # Convert starting data points from lon,lat to x,y coordinates 
    # following the aeqd transform 
    sX_aeqd, sY_aqed = p(sLon, sLat)

    # Generate a Delaunay triangulation
    sXY_aeqd = list(zip(sX_aeqd, sY_aqed))
    tri = Delaunay(sXY_aeqd)

    # Create a header and a list of data rows that will be used to create the output csv file
    header = ['no.', 'sX1_aeqd', 'sX2_aeqd', 'sX3_aeqd', 
                    'sY1_aeqd', 'sY2_aeqd', 'sY3_aeqd',
                    'vertice_idx1', 'vertice_idx2',  'vertice_idx3'
                    ]
    row_list = [header]

    # Iterate through all triangles
    for n in range(len(tri.simplices)):
            
        # Find the index of all 3 data points that form the current triangle
        vertice_idx1 = tri.simplices[n][0]
        vertice_idx2 = tri.simplices[n][1]
        vertice_idx3 = tri.simplices[n][2]

        # Retrieve the starting X and Y coordinates in the aeqd transform
        sX1_aeqd = sX_aeqd[vertice_idx1]   # Starting latitudes
        sX2_aeqd = sX_aeqd[vertice_idx2]
        sX3_aeqd = sX_aeqd[vertice_idx3]

        sY1_aeqd = sY_aqed[vertice_idx1]   # Starting longitudes
        sY2_aeqd = sY_aqed[vertice_idx2]
        sY3_aeqd = sY_aqed[vertice_idx3]

        # Add the data row corresponding to the current triangle to the list of data rows
        row_list.append([ n, sX1_aeqd, sX2_aeqd, sX3_aeqd, 
                            sY1_aeqd, sY2_aeqd, sY3_aeqd,
                            vertice_idx1, vertice_idx2,  vertice_idx3] )


    #--------------------Write the results to a csv file---------------------------------

    # Write the results in the processed_csv_path file path
    with open(processed_csv_path, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)

        # Write the data rows to the csv file
        writer.writerows(row_list)

    