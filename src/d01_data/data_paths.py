'''
Author: Beatrice Duval (bdu002)

-------------------------------------------------------------
Code that produces csv data file paths for all stages of data 
processing given an initial raw csv file and parent folder. 
-------------------------------------------------------------

'''

import os

# Input a raw csv file name and the parent folder name
raw_filename = 'pairs_20200301001804_20200301015705_1.csv'
parentfolder = '2020_MarApr_S1'

# Find absolute path in which the project is stored
curr_path = os.path.dirname(os.path.realpath(__file__))
proj_path = curr_path + '/../../'

# Find the path in which the raw csv file is stored
raw_csv_path = proj_path + 'data/01_raw/' + parentfolder + '/' + raw_filename

# Find the processing stage data csv filename and path
processed_filename = '/tri' + raw_filename[5:len(raw_filename)]
processed_csv_path = proj_path + 'data/02_processed/' + parentfolder + '/' + processed_filename

# Find the conversion stage data csv filename and path
converted_filename = '/tri_gridCS' + raw_filename[5:len(raw_filename)]
converted_csv_path = proj_path + 'data/03_converted/' + parentfolder + '/' + converted_filename

# Find the calculating stage data csv filename and path
calculations_filename = '/calc' + raw_filename[5:len(raw_filename)]
calculations_csv_path = proj_path + 'data/03_converted/' + parentfolder + '/' + calculations_filename