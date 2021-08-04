'''
Author: Beatrice Duval (bdu002)

------------------------------------------
M00 - Sea-ice deformations calculations
------------------------------------------

Code that computes triangle cells' sea-ice deformations using a converted csv file. 
The results are then stored in a calculations csv file.

1. Output csv file

    Format:
        no. | dudx | dudy | dvdx | dvdy | eps_I | eps_II | eps_tot

    Variables:
    - no.: triangle number
    - dudx, dudy, dvdx, dvdy: strain rates
    - eps_I: divergence rate
    - eps_II: maximum shear strain rate
    - eps_tot: total sea-ice deformation rate

'''

import csv
import os
from math import sqrt

import numpy as np
from numpy.core.fromnumeric import mean

import config
import utils_datetime
import utils_deformation_computations as deformation_comp
import utils_load_data as load_data


def compute_deformations():

    # Iterate through all raw, processed and converted .csv files listed in config.py
    for converted_path, calculations_path in zip(config.data_paths['converted'], config.data_paths['calculations']):
        
        # If the calculations file already exists and overwrite (in config.py) is set to false,
        # go to the next iteration.
        # ELse, process the converted file and write/overwrite the calculations file.
        if os.path.exists(calculations_path) and not config.config['Processing_options'].getboolean('overwrite'):
            continue
        
        # Check if the converted data set exists
        if not( os.path.exists(converted_path) ):
            continue

        # Load converted data set
        converted_data = load_data.load_converted(converted_path)

        sX1         = converted_data['sX1']             # Starting x coordinates in a local grid CS
        sX2         = converted_data['sX2']
        sX3         = converted_data['sX3']

        sY1         = converted_data['sY1']             # Starting y coordinates in a local grid CS
        sY2         = converted_data['sY2']
        sY3         = converted_data['sY3']

        eX1         = converted_data['eX1']             # Starting x coordinates in a local grid CS
        eX2         = converted_data['eX2']
        eX3         = converted_data['eX3']

        eY1         = converted_data['eY1']             # Starting y coordinates in a local grid CS
        eY2         = converted_data['eY2']
        eY3         = converted_data['eY3']

        # Create a header and a list of data rows that will be used to create the output csv file
        header = ['no.', 'dudx', 'dudy', 'dvdx', 'dvdy', 'eps_I', 'eps_II', 'eps_tot']
        row_list = [header]

        # Compute the time interval (days)
        dt = utils_datetime.dT( utils_datetime.dataDatetimes(calculations_path) ) / 86400
        num=0
        # Iterate through every triangle in the converted csv file 
        for n in range(len(sX1)):

            # Create lists of triangular cell vertices' positions
            sx_list = np.array([sX1[n], sX2[n], sX3[n]])  # Starting x positions
            sy_list = np.array([sY1[n], sY2[n], sY3[n]])  # Starting y positions
            ex_list = np.array([eX1[n], eX2[n], eX3[n]])  # Ending x positions
            ey_list = np.array([eY1[n], eY2[n], eY3[n]])  # Ending y positions
            
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
