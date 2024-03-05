"""
Author: Mathieu Plante, Lekima Yakuden
GitHub: mathieuslplante, LekiYak

--------------------------------------------------------------------------------
Code to investigate SIDRR dataset format and caracteristics
--------------------------------------------------------------------------------

This code is used to load the SIDRR data from the daily netcdf files,
plot the included data and determine the spatio-temporal coverage.

SIDRR data computed from SAR image pairs are stored in daily netcdf files,
according to the acquistion time of the earliest image from the pair.
SIDRR data spanning a given date are found accross several netcdf files.

To load all data valid for a specific date, it is necessary to:

1. Span cdf files from days prior to the analysis initial time,
2. Keep track of the previously loaded data until past the
   acquisition time of the latest image from the pair.

"""

# Loading from default packages
import os
import sys
parent = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,parent)
sys.path.append(r'src/')
from time import strftime
import configparser
from datetime import datetime, date, time, timedelta
import time
import numpy as np
#from netCDF4 import Dataset

# Code from src files
from visualisation import visualisation
from LoadDataset   import SID_dataset
from TimeUtil import TimeUtil
from Statistics_objects import Coverage_map, Data_Distribution


if __name__ == '__main__':

    # Start the clock
    start_time = time.time()
    Delta_days_init = 6

    # Reading config
    config = configparser.ConfigParser()
    config.read('options.ini')

    path = config['IO']['netcdf_path']
    TimeTool = TimeUtil(config = config['Date_options'])

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

    refTime = TimeTool.ThisTime
    Data_Mem = None

    #Looping over earlier days
    for tic in range(0,Delta_days_init):

        ThisTime = TimeTool.ThisTime - timedelta(days=Delta_days_init-tic)
        NextTime = ThisTime + timedelta(seconds=TimeTool.tstep*60*60)
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
    # STARTING ANALYSIS
    #------------------------------------------------------------------

    # Head_start is the Delta-t associated with different reference date
    # in the data carried over from previous dates.
    Head_start = Delta_days_init*60*60*24

    #Initialize statistics objects
    if config['netcdf_tools']['show_spatial_coverage']:
        resolution = float(config['options']['resolution'])
        FreqMap = Coverage_map(resolution = resolution)

    #Initialise 1d histograms for the triangle areas.
    if config['netcdf_tools']['show_spatial_scale_dist']:
        A_dist_all = Data_Distribution(LeftBC= 0.0, RightBC = 20.0, nbins = 80, label = "all")
        A_dist_S1 = Data_Distribution(LeftBC= 0.0, RightBC = 20.0, nbins = 80, label = "S1")
        A_dist_RCM = Data_Distribution(LeftBC= 0.0, RightBC = 20.0, nbins = 80, label = "RCM")


    # Iterating over each day
    for ThisTime in TimeTool.daterange():

        #------------------------------------------------------------------
        # Fetch data
        #------------------------------------------------------------------

        # Update time and data paths
        TimeTool.ThisTime = ThisTime
        TimeTool.NextTime = TimeTool.ThisTime + timedelta(seconds=TimeTool.tstep*60*60)

        #Use to limit the analysis to specific months. Otherwise, ignore.
        if TimeTool.ThisTime.month > 12 and TimeTool.ThisTime.month < 1:
            continue
        ThisTime_str = TimeTool.ThisTime.strftime("%Y%m%d")
        NextTime_str = TimeTool.NextTime.strftime("%Y%m%d")
        Sat = config['Metadata']['icetracker']
        TimeTool.ThisTimeFile = "%sSID_%s_%s_dt72_tol72_dx.nc" % (Sat,ThisTime_str, NextTime_str)
        filePath = path + TimeTool.ThisTimeFile

        #Load visualisation tool
        visuals = visualisation(config=config)

        #Load new Data from SIDRR and stack to carrier
        Data = SID_dataset(FileName= filePath, config=config)
        Data.start_time = Data.start_time + Head_start
        Data.end_time = Data.end_time + Head_start
        Data.Concatenate_data(Data2 = Data_Mem)

        #------------------------------------------------------------------
        # Produce figures to visualise the dataset characteristics
        #------------------------------------------------------------------

        #Figure showing the start and end points in a file
        if config['netcdf_tools']['plot_start_end_points']:
           visuals.plot_start_end_points(data = Data, datestring  = ThisTime_str)

        #Figure showing the stacked SAR image areas.
        if config['netcdf_tools']['plot_stacked_pairs']:
           visuals.show_stacked_pairs(data = Data, datestring  = ThisTime_str)

        #Figure showing the triangulated data of specified SAR image pair ID.
        if config['netcdf_tools']['plot_triangulated_data']:
            visuals.plot_triangles(data=Data, no = int(2), triangle_zoom = True, datestring = ThisTime_str)

        #Figure showing normal, shear and rotation rates.
        if config['netcdf_tools']['plot_deformation']:
            visuals.plot_deformations(data = Data, datestring = ThisTime_str)
            visuals.show_tripcolor_field(data=Data, Field = Data.A,
                                          title = "Triangle Areas", label = "Area",
                                          datestring=ThisTime_str)
        #Figure showing the spatial coverage of the dataset
        if config['netcdf_tools']['show_spatial_coverage']:
            FreqMap.add2hist_2D(Data = Data)

        #Distribute the triangle areas in 1d histogram, and add to bins.
        if config['netcdf_tools']['show_spatial_scale_dist']:
            A_dist_all.add2hist_1D(Data = Data.A[:])
            A_dist_S1.add2hist_1D(Data = Data.A[Data.satellite[:]==1])
            A_dist_RCM.add2hist_1D(Data = Data.A[Data.satellite[:]==0])

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

    visuals = visualisation(config = config)

    # Make figure showing the spatio-temportal coverage of the SIDRR data in the analysed period
    if config['netcdf_tools']['show_spatial_coverage']:
        visuals.show_spatial_coverage(distribution_2D = FreqMap, datestring = TimeTool.StartDate_str + '_' + TimeTool.EndDate_str)

    # Make figure showing the distribution of data spatial scales
    if config['netcdf_tools']['show_spatial_scale_dist']:
        visuals.plot_area_dist(dist1 = A_dist_S1,
                               datestring = TimeTool.StartDate_str + '_' + TimeTool.EndDate_str)

    # Display the computation time
    print("--- %s seconds ---" % (time.time() - start_time))
