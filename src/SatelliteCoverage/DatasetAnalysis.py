"""
Author: Mathieu Plante, Lekima Yakuden
GitHub: mathieuslplante, LekiYak

--------------------------------------------------------------------------------
Tools for analysing and processing netCDF files
--------------------------------------------------------------------------------

This file contains functions for analysing and processing netCDF files.

"""

# Loading from default packages
import os
import sys
parent = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,parent)
from time import strftime
from datetime import datetime, date, time, timedelta
import time
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
from tqdm import tqdm
import haversine as hs

# Code from other files
from visualisation import visualisation
from LoadDataset   import SID_dataset
from SatelliteCoverage.config import read_config
from TimeUtil import TimeUtil
from SIDD_data_object import SIDD_Data
sys.path.append(r'/storage/mathieu/SIDD/src/')
from SIDD import SIDD

if __name__ == '__main__':

    # Start the clock
    start_time = time.time()
    Delta_days_init = 6
    # Reading config
    config = read_config()
    path = config['IO']['netcdf_path']
    time = TimeUtil(config = config['Date_options'])
    print(config['netcdf_tools'])
    if config['netcdf_tools']['compute_pdfs']:
        namelist_SIDD = config['SIDD_options']
        OutputFolder =  namelist_SIDD['outputfolder']+ config['IO']['exp'] + '/'
        analysis = SIDD(config = namelist_SIDD,
                        time=time)

    #--------------------------------------------------------------
    # Initialise with 6 earlier days to get all data for given day.
    #    This step is necessary now that the data is stored according to
    #    the start time. Here, we fetch data from earlier date that are
    #    valid for the first date of the required analysis period.
    #
    #    We are thus building a Memory term containing data to be carried
    #    over the next step
    #    - Data_Mem: Earlier data carrier
    #--------------------------------------------------------------

    refTime = time.ThisTime
    Data_Mem = None

    #Looping over earlier days
    for tic in range(0,Delta_days_init):

        ThisTime = time.ThisTime - timedelta(days=Delta_days_init-tic)
        NextTime = ThisTime + timedelta(seconds=time.tstep*60*60)
        ThisTime_str = ThisTime.strftime("%Y%m%d")
        NextTime_str = NextTime.strftime("%Y%m%d")

        # Head_start is the Delta-t associated with different reference date
        # in the earlier files
        Head_start = (tic)*60*60*24
        ref = (ThisTime-refTime).total_seconds()

        #Fetch the path to netcdf
        Sat = config['Metadata']['icetracker']
        ThisTimeFile = "%sSID_%s_%s_dt72_tol72_dx.nc" % (Sat,ThisTime_str, NextTime_str)
        filePath = path + ThisTimeFile

        #Load netcdf data
        Data = SID_dataset(FileName= filePath, config=config)
        indices = [i for i in range(len(Data.A)) if  Data.end_time[i] > -ref ]

        #Filter to keep data that are valid for the aim date
        if len(indices) == 0:
            continue
        Data.filter_data(indices = indices)
        Data.start_time = Data.start_time + Head_start
        Data.end_time = Data.end_time + Head_start

        #Stack in the Mem data carrier
        if tic == 0 or Data_Mem is None:
            Data_Mem = Data
        elif len(Data_Mem.A) == 0:
            Data_Mem = Data
        else:
            Data_Mem.Concatenate_data(Data2 = Data)
            print("Data_Mem length is now: ", len(Data_Mem.A[:]))
        Data_Mem.day_flag = Data_Mem.day_flag + 1

    #------------------------------------------------------------------
    # Do data analysis for each date in namelist time period
    #------------------------------------------------------------------

    # Head_start is the Delta-t associated with different reference date
    # in the data carried over from previous dates.
    Head_start = Delta_days_init*60*60*24

    # Iterating over each day
    for ThisTime in time.daterange():

        # Update time and data paths
        time.ThisTime = ThisTime
        time.NextTime = time.ThisTime + timedelta(seconds=time.tstep*60*60)
        ThisTime_str = time.ThisTime.strftime("%Y%m%d")
        NextTime_str = time.NextTime.strftime("%Y%m%d")
        Sat = config['Metadata']['icetracker']
        time.ThisTimeFile = "%sSID_%s_%s_dt72_tol72_dx.nc" % (Sat,ThisTime_str, NextTime_str)
        filePath = path + time.ThisTimeFile

        #Load visualisation tool
        visuals = visualisation()
        #Load new Data from SIDRR and stack to carrier
        Data = SID_dataset(FileName= filePath, config=config)
        Data.start_time = Data.start_time + Head_start
        Data.end_time = Data.end_time + Head_start
        Data.Concatenate_data(Data2 = Data_Mem)

        #Mask triangles with too large area. This should not be hardcoded!!
        indices = [i for i in range(len(Data.A)) if  Data.A[i] > 20000**2.0 ]
        Data.mask_data(indices = indices)

        if config['netcdf_tools']['plot_start_end_points']:
           # #Figure showing the start and end points in a file
           # visuals.plot_start_end_points(data = Data, config=config)

           # #Figure showing the stacked SAR image areas.
           # visuals.plot_pairs(data = Data, config = config)

            #Figure showing the triangulated data of specified SAR image pair ID.
            visuals.plot_triangles(data=Data, config=config, no = [2], show_ID = True)
            visuals.plot_triangles(data=Data, config=config, no = [34], show_ID = True)
            visuals.plot_triangles(data=Data, config=config, no = [1], show_ID = True)

        if config['netcdf_tools']['plot_deformation']:
            visuals.plot_deformations(data = Data, config=config, datestring = ThisTime_str)

        # Should be moved to Data analysis tool, not part of the published dataset.
        if config['netcdf_tools']['compute_pdfs']:
            SIDDdata = SIDD_Data(InputData = Data)
            analysis.add2pdf(data = SIDDdata, time = time)
            analysis.show_pdf(time = time,output_folder = OutputFolder  + "figs/PDFs/" + ThisTime_str)

        #--------------------------------------------------------------------------------
        #Update the data carrier for the next timestep
        #--------------------------------------------------------------------------------
        Data_Mem = Data
        indices = [i for i in range(len(Data_Mem.A)) if  Data_Mem.end_time[i] > Head_start + 24*60*60]
        Data_Mem.filter_data(indices = indices)
        print('Removing the data with earlier dates: ', len(Data.A),len(Data_Mem.A))
        Data_Mem.start_time = Data_Mem.start_time - 60*60*24
        Data_Mem.end_time = Data_Mem.end_time - 60*60*24
        Data_Mem.day_flag = Data_Mem.day_flag + 1
        del visuals
        del Data

    #Save SIDD netcdf analysis.
    if config['netcdf_tools']['compute_pdfs']:
        analysis.Save2netcdf(time = time, output_folder = OutputFolder + ThisTime_str)


    # Display the computation time
    print("--- %s seconds ---" % (time.time() - start_time))
