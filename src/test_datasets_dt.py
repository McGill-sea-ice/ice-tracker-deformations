'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------------------------
Code that creates a histogram of delta times for all S1 datasets.
-----------------------------------------------------------------

'''

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

import utils_datetime
import utils_load_data as load_data

'''
_________________________________________________________________________________________
CREATE A DELTA TIMES LIST
'''

# Initialize a delta times list
dt_list = []

# Iterate through all files in the folder specified below
folder = '/home/bdu002/2021_SeaIceDeformation/data/01_raw/2020_MarApr_S1'
for filename in os.listdir(folder):

    # Add the delta time list only if the current file would be processed in the main script 
    # (i.e. is a .csv or .dat file and has more than 3 data points)
    try:
        load_data.load_raw(folder + '/' + filename)
        
        # Compute the time interval (hours)
        dt = utils_datetime.dT( utils_datetime.dataDatetimes(filename) ) * 24

        # Append the time interval to the delta times list
        dt_list.append(dt)

    except load_data.DataFileError as dfe:
            continue

print(len(dt_list))

'''
_________________________________________________________________________________________
CREATE A BINS LIST FOR THE HISTOGRAM
'''

# The first elements of the bin list are 0 and 6
bins = [0,6]

# Add elements to the list by increments of 12, until we have enough bins to cover the data
while bins[len(bins)-1] <= max(dt_list):
    lb = bins[len(bins)-1]
    bins.append(lb+12)

'''
_________________________________________________________________________________________
PLOT THE HISTOGRAM
'''

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Time Intervals of all RCM Data Sets for March and April 2020')
ax.hist(dt_list, bins=bins, weights=np.ones(len(dt_list)) / len(dt_list), edgecolor='k')
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))

plt.show()


