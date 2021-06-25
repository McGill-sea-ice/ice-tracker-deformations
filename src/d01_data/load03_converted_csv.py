'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------
Code that loads data from a converted csv file.
-----------------------------------------------

'''

from numpy.core.numeric import array_equal
import pandas as pd
from src.d01_data.get_data_paths import converted_csv_path

# Create a data frame
df = pd.read_csv(converted_csv_path)

# Retrieve data
no          = df['no.']             # Triangle number

sX1         = df['sX1']             # Starting x coordinates in a local grid CS
sX2         = df['sX2']
sX3         = df['sX3']

sY1         = df['sY1']             # Starting y coordinates in a local grid CS
sY2         = df['sY2']
sY3         = df['sY3']

eX1         = df['eX1']             # Starting x coordinates in a local grid CS
eX2         = df['eX2']
eX3         = df['eX3']

eY1         = df['eY1']             # Starting y coordinates in a local grid CS
eY2         = df['eY2']
eY3         = df['eY3']

j_tracers   = df['tracer_j']        # Tracer point (defines the local CS) indices 
i_tracers   = df['tracer_i']

