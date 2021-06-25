'''
Author: Beatrice Duval (bdu002)

--------------------------------------------------------------------
Code that performs a Delaunay triangulation on a set of data points 
and stores the results in a csv file.
--------------------------------------------------------------------

Output csv file format:

    no., sLat1, sLat2, sLat3, sLon1, sLon2, sLon3, eLat1, eLat2, eLat3, eLon1, eLon2, eLon3, sX1_aeqd, sX2_aeqd, sX3_aeqd, sY1_aeqd, sY2_aeqd, sY3_aeqd, vertice_idx1, vertice_idx2

where:

    - no. is the triangle number;
    - sLat1, sLat2, sLat3 are the starting latitudes of the 1st, 2nd and 3rd vertices;
    - sLon1, sLon2, sLon3 are the starting longitudes;
    - eLat1, eLat2, eLat3, eLon1, eLon2, eLon3 are the ending latitudes and latitudes;
    - sX1_aeqd, sX2_aeqd, sX3_aeqd, sY1_aeqd, sY2_aeqd, sY3_aeqd are the x,y coordinates 
        of the starting vertices in the Azimuthal Equidistant (aeqd) transform.
    - vertice_idx1, vertice_idx2, vertice_idx3 are the triangle vertices' indices in the raw csv file
'''

import csv

from pyproj import Proj
from scipy.spatial import Delaunay
from src.d01_data.get_data_paths import processed_csv_path
from src.d01_data.load01_raw_csv import sLat, sLon, eLat, eLon

# Convert starting data points from lon,lat to x,y coordinates 
# following the Azimuthal Equidistant (aeqd) transform 
p = Proj(proj='aeqd', ellps='WGS84', preserve_units=False)
sX_aeqd, sY_aqed = p(sLon, sLat)

# Generate a Delaunay triangulation
sXY_aeqd = list(zip(sX_aeqd, sY_aqed))
tri = Delaunay(sXY_aeqd)

# Create a header and a list of data rows that will be used to create the output csv file
header = ['no.', 'sLat1', 'sLat2', 'sLat3',             
                 'sLon1', 'sLon2', 'sLon3', 
                 'eLat1', 'eLat2', 'eLat3', 
                 'eLon1', 'eLon2', 'eLon3', 
                 'sX1_aeqd', 'sX2_aeqd', 'sX3_aeqd', 
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

    # Retrieve the starting and ending latitudes and longitudes of the data points
    sLat1 = sLat[vertice_idx1]   # Starting latitudes
    sLat2 = sLat[vertice_idx2]
    sLat3 = sLat[vertice_idx3]

    sLon1 = sLon[vertice_idx1]    # Starting longitudes
    sLon2 = sLon[vertice_idx2]
    sLon3 = sLon[vertice_idx3]

    eLat1 = eLat[vertice_idx1]    # Ending latitudes
    eLat2 = eLat[vertice_idx2]
    eLat3 = eLat[vertice_idx3]

    eLon1 = eLon[vertice_idx1]    # Ending longitudes
    eLon2 = eLon[vertice_idx2]
    eLon3 = eLon[vertice_idx3]  
        
    # Retrieve the starting X and Y coordinates in the aeqd transform
    sX1_aeqd = sX_aeqd[vertice_idx1]   # Starting latitudes
    sX2_aeqd = sX_aeqd[vertice_idx2]
    sX3_aeqd = sX_aeqd[vertice_idx3]

    sY1_aeqd = sY_aqed[vertice_idx1]   # Starting longitudes
    sY2_aeqd = sY_aqed[vertice_idx2]
    sY3_aeqd = sY_aqed[vertice_idx3]

    # Add the data row corresponding to the current triangle to the list of data rows
    row_list.append([ n, sLat1, sLat2, sLat3, 
                         sLon1, sLon2, sLon3, 
                         eLat1, eLat2, eLat3,
                         eLon1, eLon2, eLon3, 
                         sX1_aeqd, sX2_aeqd, sX3_aeqd, 
                         sY1_aeqd, sY2_aeqd, sY3_aeqd,
                         vertice_idx1, vertice_idx2,  vertice_idx3] )


#--------------------Write the results to a csv file---------------------------------

# Write the results in the processed_csv_path file path
with open(processed_csv_path, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # Write the data rows to the csv file
    writer.writerows(row_list)

    