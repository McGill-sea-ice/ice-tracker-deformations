'''
Author: Beatrice Duval (bdu002)

-----------------------------------------
Code that loads data from a raw csv file.
-----------------------------------------

'''

import pandas as pd
from src.d01_data.data_paths import raw_csv_path

# Create a data frame
df = pd.read_csv(raw_csv_path)

# Retrieve data points
sLon = df['sLon']       # Starting longitudes
sLat = df['sLat']       # Starting latitudes
eLon = df['eLon']       # Ending longitudes
eLat = df['eLat']       # Ending latitudes


