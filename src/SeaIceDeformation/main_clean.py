'''
Author: Beatrice Duval (bdu002)

---------------------------------
Main script for data processing
--------------------------------

Script that executes all steps towards the calculation of sea-ice deformations and displays the execution time.

'''

# Loading from default packages
import os
import sys
parent = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,parent)
import time
import datetime
from tqdm import tqdm
# Loading from other files
from config import get_config_args, filter_data
from config import divide_intervals, get_datapaths

from make_triangulated_array import make_into_tri_arrays
from CDFoutput import CDFoutput
from ASITS_dat_files import Load_ASITS_dat_file
from Deformations import compute_SIDRR

# Retrieve the starting time
start_time = time.time()
'''
0) Configuration extraction
_______________________________________________________________________
PERFORM CONFIGURATION
'''

# Retrieve configuration arguments from namelist.ini
config = get_config_args()
raw_paths = filter_data(config=config)
date_pairs = divide_intervals(config=config)

# Iterating over each (daily) interval
for i in tqdm(range(len(date_pairs)), position=0, leave=False):

    #Set the start and end time according to the given date
    datepairs = date_pairs[i]
    config['Date_options']['start_year'] = datepairs[0].strftime("%Y")
    config['Date_options']['start_month'] = datepairs[0].strftime("%m")
    config['Date_options']['start_day'] = datepairs[0].strftime("%d")
    config['Date_options']['end_year'] = datepairs[1].strftime("%Y")
    config['Date_options']['end_month'] = datepairs[1].strftime("%m")
    config['Date_options']['end_day'] = datepairs[1].strftime("%d")

    DateString = datepairs[0].strftime("%Y%m%d")
    # Loads sea ice motion data for the given interval
    raw_paths = filter_data(config=config)
    config['raw_paths'] = raw_paths

    #Recording files treated or discarded
    numtot = 0
    num_short = 0
    num_tri_fault = 0
    num_data = 0
    num_large = 0
    idpair = 0
    nbfb = 0
    nbfg = 0
    #Initializing the data arrays that will be stored in netcdf
    CDFdata = CDFoutput()

    for path in raw_paths:
        numtot += 1

        #Get data
        data = Load_ASITS_dat_file(FileName= path, config=config)
        if len(data.sLat) < 3:
            del data
            num_tri_fault += 1
            continue
        if data.dt < 0.5:
            del data
            num_short += 1
            print(path)
            continue
        #Triangulate data
        tri_array = make_into_tri_arrays(data = data,config=config)
        if tri_array.fault == 1:
            num_tri_fault += 1
            del data
            del tri_array
            continue

        #Compute sea-ice deformations rates and stack in daily netcdf output.
        idpair += 1
        SIDRRdata = compute_SIDRR( config = config,
                                   tri_array = tri_array,
                                   data = data,
                                   CDFout = CDFdata,
                                   idpair = idpair)
        num_data += SIDRRdata.num_tri
        num_large += SIDRRdata.num_large

        nbfb += SIDRRdata.nbfb
        nbfg += SIDRRdata.nbfg
        del SIDRRdata
        del tri_array
        del data
        num_ok = numtot-num_tri_fault-num_short
        print("Out of %s file processed, %s are ok, %s too short, %s could not triangulate." % (numtot, num_ok, num_short, num_tri_fault))
        print(nbfb+nbfg, nbfb)
        print("Printing netcdf for %s" % datepairs[0])
    CDFdata.write_CDF(config = config, DateString = DateString)
    del CDFdata

# Display the run time
print("--- %s seconds ---" % (time.time() - start_time))


