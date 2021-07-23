'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------------------------------------------
Code that computes triangle cells' sea-ice deformations using a triangulated csv file. 

The results are then stored in a calculations csv file.
-----------------------------------------------------------------------------------

Output csv file format:
    __________________________________________________________
    no. | dudx | dudy | dvdx | dvdy | eps_I | eps_II | eps_tot

where:
    - no. is the triangle number;
    - dudx, dudy, dvdx, dvdy are strain rates;
    - eps_I is the divergence rate;
    - eps_II is the maximum shear strain rate;
    - eps_tot is the total sea-ice deformation rate;

'''

import csv
import os
from math import sqrt

import config
import utils_datetime
import utils_deformation_computations as deformation_comp
import utils_load_data as load_data


def compute_deformations():

    # Iterate through all raw and triangulated data files listed in config
    for raw_path, triangulated_path, calculations_path in zip(config.data_paths['raw'], config.data_paths['triangulated'], config.data_paths['calculations']):
        
        # If the calculations file already exists and overwrite (in args config) is set to false,
        # go to the next iteration.
        # ELse, process the triangulated file and write/overwrite the calculations file.
        if os.path.exists(calculations_path) and not config.args.overwrite:
            continue
        
        # Check if the triangulated data set exists. If it does not, go to the next iteration
        if not( os.path.exists(triangulated_path) ):
            continue

        # Load the raw data set. If an error is encountered (no or not enough data points), 
        # go to the next dataset.
        try:
            raw_data = load_data.load_raw( raw_path )
            startX  = raw_data['startX']  # Starting X positions (px)
            startY  = raw_data['startY']  # Starting Y positions (px)
            endX    = raw_data['endX']    # Ending X positions (px)
            endY    = raw_data['endY']    # Ending Y positions (px)

        except load_data.DataFileError as dfe:
            # The error has already been printed in the delaunay triangulation stage            
            continue
        
        # Load the triangulated dataset
        triangulated_data = load_data.load_triangulated( triangulated_path )
        vertice_idx1 = triangulated_data['vertice_idx1'] # Vertex indices in raw csv file
        vertice_idx2 = triangulated_data['vertice_idx2']
        vertice_idx3 = triangulated_data['vertice_idx3']

        # Create a header and a list of data rows that will be used to create the output csv file
        header = ['no.', 'dudx', 'dudy', 'dvdx', 'dvdy', 'eps_I', 'eps_II', 'eps_tot']
        row_list = [header]

        # Compute the time interval (days)
        dt = utils_datetime.dT( utils_datetime.dataDatetimes(raw_path) )
     
        # Iterate through every triangle in the converted csv file 
        for n in range(len(vertice_idx1)):
            
            # Retrieve the index of each triangle vertex in the raw dataset
            i1 = vertice_idx1[n]
            i2 = vertice_idx2[n]
            i3 = vertice_idx3[n]

            # Retrieve the starting and ending X/Y positions of each triangle vertices
            sX1 = startX[i1]
            sX2 = startX[i2]
            sX3 = startX[i3]

            sY1 = startY[i1]             
            sY2 = startY[i2]   
            sY3 = startY[i3]   

            eX1 = endX[i1]             
            eX2 = endX[i2]        
            eX3 = endX[i3]        

            eY1 = endY[i1]          
            eY2 = endY[i2] 
            eY3 = endY[i3] 

            # Create lists of triangular cell vertices' positions
            sx_list = [sX1, sX2, sX3]  # Starting x positions
            sy_list = [sY1, sY2, sY3]  # Starting y positions
            ex_list = [eX1, eX2, eX3]  # Ending x positions
            ey_list = [eY1, eY2, eY3]  # Ending y positions
            
            # Create a list of velocity components for each triangle vertex
            u_list, v_list = deformation_comp.calculate_uv_lists( sx_list, ex_list, sy_list, ey_list, dt)

            # Compute the strain rates
            dudx, dudy, dvdx, dvdy = deformation_comp.calculate_strainRates( u_list, v_list, sx_list, sy_list )

            # Compute the divergence rate
            eps_I = dudx + dvdy

            # Compute the maximum shear strain rate
            eps_II = sqrt(  (dudx - dvdy)**2 + (dudy + dvdx)**2  )

            # Compute the total sea-ice deformation rate
            eps_tot = sqrt( eps_I**2 + eps_II**2 )
            
            # Add the data row corresponding to the current triangle to the list of data rows
            row_list.append( [n, dudx, dudy, dvdx, dvdy, eps_I, eps_II, eps_tot] )
            

        #--------------------Write the results to a csv file---------------------------------

        # Create a directory to store the calculations csv path if it does not exist already
        os.makedirs(os.path.dirname(calculations_path), exist_ok=True)

        # Write the results in the calculations_csv_path file path
        with open(calculations_path, 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)

            # Write the data rows to the csv file
            writer.writerows(row_list)


if __name__ == '__main__':
    compute_deformations()