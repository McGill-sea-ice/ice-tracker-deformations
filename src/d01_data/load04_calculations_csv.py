'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------
Code that loads data from a calculations csv file.
-----------------------------------------------

'''

import pandas as pd
from src.d01_data.get_data_paths import calculations_csv_path

# Create a data frame
df = pd.read_csv(calculations_csv_path)

# Retrieve data
no          = df['no.']             # Triangle number

dudx        = df['dudx']            # Strain rates
dudy        = df['dudy']
dvdx        = df['dvdx']
dvdy        = df['dvdy']

eps_I       = df['eps_I']           # Divergence rate

eps_II      = df['eps_II']          # Strain rate

eps_tot     = df['eps_tot']         # Total deformation rate


