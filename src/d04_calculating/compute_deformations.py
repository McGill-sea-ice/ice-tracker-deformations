'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------------------------------------------
Code that computes triangle cells' sea-ice deformations using a converted csv file. 

The results are then stored in a calculations csv file.
-----------------------------------------------------------------------------------

Output csv file format:
    no., dudx, dudy, dvdx, dvdy, eps_I, eps_II, eps_tot

where:
    - no. is the triangle number;
    - dudx, dudy, dvdx, dvdy are strain rates;
    - eps_I is the divergence rate;
    - eps_II is the maximum shear strain rate;
    - eps_tot is the total sea-ice deformation rate;

'''

import csv
from math import sqrt

from src.d00_utils.datetime_raw_csv import dataDatetimes, dT
from src.d00_utils.deformation_computations import (calculate_strainRates,
                                                    calculate_uv_lists)
from src.d01_data.get_data_paths import calculations_csv_path
from src.d01_data.load03_converted_csv import (eX1,eX2, eX3, eY1, eY2, eY3, 
                                               sX1, sX2, sX3, sY1, sY2, sY3)

# Create a header and a list of data rows that will be used to create the output csv file
header = ['no.', 'dudx', 'dudy', 'dvdx', 'dvdy', 'eps_I', 'eps_II', 'eps_tot']
row_list = [header]

# Compute the time interval (days)
dt = dT( dataDatetimes(calculations_csv_path) )

# Iterate through every triangle in the converted csv file 
for n in range(len(sX1)):

    # Create lists of triangular cell vertices' positions
    sx_list = [sX1[n], sX2[n], sX3[n]]  # Starting x positions
    sy_list = [sY1[n], sY2[n], sY3[n]]  # Starting y positions
    ex_list = [eX1[n], eX2[n], eX3[n]]  # Ending x positions
    ey_list = [eY1[n], eY2[n], eY3[n]]  # Ending y positions
    
    # Create a list of velocity components for each triangle vertex
    u_list, v_list = calculate_uv_lists( sx_list, ex_list, sy_list, ey_list, dt)

    # Compute the strain rates
    dudx, dudy, dvdx, dvdy = calculate_strainRates( u_list, v_list, sx_list, sy_list )

    # Compute the divergence rate
    eps_I = dudx + dvdy

    # Compute the maximum shear strain rate
    eps_II = sqrt(  (dudx - dvdy)**2 + (dudy + dvdx)**2  )

    # Compute the total sea-ice deformation rate
    eps_tot = sqrt( eps_I**2 + eps_II**2 )
    
    # Add the data row corresponding to the current triangle to the list of data rows
    row_list.append( [n, dudx, dudy, dvdx, dvdy, eps_I, eps_II, eps_tot] )
    

#--------------------Write the results to a csv file---------------------------------

# Write the results in the calculations_csv_path file path
with open(calculations_csv_path, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # Write the data rows to the csv file
    writer.writerows(row_list)
