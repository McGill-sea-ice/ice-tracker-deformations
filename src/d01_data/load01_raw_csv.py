'''
Author: Beatrice Duval (bdu002)

-----------------------------------------
Code that loads data from a raw csv file.
-----------------------------------------

'''

import pandas as pd
from src.d01_data.get_data_paths import raw_csv_path, raw_filename

# Create a data frame
df = pd.read_csv(raw_csv_path)

# Retrieve data points
sLon = df['sLon']       # Starting longitudes
sLat = df['sLat']       # Starting latitudes
eLon = df['eLon']       # Ending longitudes
eLat = df['eLat']       # Ending latitudes

# Terminate the process if the input raw csv file is empty
if ( len(sLon) < 3 ):
    print("The input raw csv file " + raw_filename + " is empty or does not have enough data to perform a triangulation. The process has been terminated.")
    exit()

# Terminate the process if the data points in the input raw csv file are not in the region of interest
