'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------
Code that loads data from a processed csv file.
-----------------------------------------------

'''

import pandas as pd
from src.d01_data.data_paths import processed_csv_path

# Create a data frame
df = pd.read_csv(processed_csv_path)

# Retrieve data
no       = df['no.']          # Triangle number

sLat1    = df['sLat1']        # Starting latitudes
sLat2    = df['sLat2']
sLat3    = df['sLat3']

sLon1    = df['sLon1']        # Starting longitudes
sLon2    = df['sLon2']
sLon3    = df['sLon3']

eLat1    = df['eLat1']        # Ending latitudes
eLat2    = df['eLat2']
eLat3    = df['eLat3']

eLon1    = df['eLon1']        # Ending longitudes
eLon2    = df['eLon2'] 
eLon3    = df['eLon3']  

sX1_aeqd = df['sX1_aeqd']     # Starting X coordinates in the aeqd transform
sX2_aeqd = df['sX2_aeqd']
sX3_aeqd = df['sX3_aeqd']

sY1_aeqd = df['sY1_aeqd']     # Starting Y coordinates in the aeqd transform
sY2_aeqd = df['sY2_aeqd']
sY3_aeqd = df['sY3_aeqd']

