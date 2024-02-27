"""
Authors: Lekima Yakuden, Mathieu Plante, Amelie Bouchat, Damien Ringeisen
GitHub: LekiYak

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
from SatelliteCoverage.config import read_config
from SatelliteCoverage.utils import date_to_seconds, seconds_to_date, convert_to_grid, get_prefix




# Loads netCDF data
class SID_dataset:

    def __init__(self,FileName= None, config=None):
        """
        This function reads and loads data from a netCDF file in the same format as those output by the deformation
        calculation script in src/SeaIceDeformation. The user is able to filter the data by time and area,
        allowing for analytical tools to be applied to selected snippets of data.

        INPUTS:
        path -- String of path to the netCDF file the data will be read from {str}

        OUTPUTS:
        dictionary object -- Dictionary storing filtered data as NumPy arrays {dict}
        """

        print('--- Loading data (netcdf) ---')

        path = FileName
        print(path)

        # Reading config
        Date_options = config['Date_options']
        options = config['options']

        # Converting dates for title and file name purposes
        self.timestep    = Date_options['timestep']
        self.tolerance   = Date_options['tolerance']
        self.resolution  = config['options']['resolution']
        self.interval    = config['options']['interval']
        self.icetracker  = config['Metadata']['icetracker']
        area_filter = options['area_filter']
        # Reading user options
        start_year, start_month, start_day = Date_options['start_year'], Date_options['start_month'], Date_options['start_day']
        end_year, end_month, end_day = Date_options['end_year'], Date_options['end_month'], Date_options['end_day']

        self.area_filter= options['area_filter']

        # Load netCDF as Dataset from *path*
        ds = Dataset(path, mode='r')

        # Get start / end times from user
        self.start_time_s, self.end_time_s = date_to_seconds(path, start_year, start_month, start_day, end_year, end_month, end_day)

        # Extracting time variables
        self.start_time = ds.variables['start_time'][:]
        self.end_time = ds.variables['end_time'][:]

        # Indices of data in desired time frame (not needed in daily files)
        #time_indices = np.where( ((start_time < end_time_s) & (end_time > start_time_s)))[0]
        #start_time = start_time[time_indices]
        #end_time = end_time[time_indices]

        # Extracting data (Filtered by time only)
        self.start_lat1 = (ds.variables['start_lat1'][:])[:]
        self.start_lat2 = (ds.variables['start_lat2'][:])[:]
        self.start_lat3 = (ds.variables['start_lat3'][:])[:]

        self.start_lon1 = (ds.variables['start_lon1'][:])[:]
        self.start_lon2 = (ds.variables['start_lon2'][:])[:]
        self.start_lon3 = (ds.variables['start_lon3'][:])[:]

        self.end_lat1 = (ds.variables['end_lat1'][:])[:]
        self.end_lat2 = (ds.variables['end_lat2'][:])[:]
        self.end_lat3 = (ds.variables['end_lat3'][:])[:]

        self.end_lon1 = (ds.variables['end_lon1'][:])[:]
        self.end_lon2 = (ds.variables['end_lon2'][:])[:]
        self.end_lon3 = (ds.variables['end_lon3'][:])[:]

        self.A = (ds.variables['A'][:])[:]

        self.div = (ds.variables['div'][:])[:]
        self.shr = (ds.variables['shr'][:])[:]
        self.vrt = (ds.variables['vrt'][:])[:]

        self.idx1 = (ds.variables['idx1'][:])[:]
        self.idx2 = (ds.variables['idx2'][:])[:]
        self.idx3 = (ds.variables['idx3'][:])[:]
        self.no   = (ds.variables['no'][:])[:]

        self.dudx = (ds.variables['dudx'][:])[:]
        self.dudy = (ds.variables['dudy'][:])[:]
        self.dvdx = (ds.variables['dvdx'][:])[:]
        self.dvdy = (ds.variables['dvdy'][:])[:]

        #min_date = seconds_to_date(path,np.min(start_time))
        #max_date = seconds_to_date(path,np.max(end_time))
        #print(f"Earliest/latest start/end dates of data included in deformation map: {min_date} to {max_date}")

        reftime = ds.getncattr('referenceTime')
        icetracker = ds.getncattr('icetracker')
        timestep = ds.getncattr('timestep')
        tolerance = ds.getncattr('tolerance')
        trackingerror = ds.getncattr('trackingError')

        #closing the dataset
        ds.close()
        self.Mask = self.A.copy()
        self.Mask[:] = 1

        self.day_flag = self.no.copy()
        self.day_flag[:]= 1
#        if area_filter:
#            print(np.nanmean(self.A))
#            print(len(self.A))
            # Filtering indices where all three points (entire triangle) is within the radius
        indices = [i for i in range(len(self.A)) if  self.A[i] < 20000.0**2.0 ]
        self.Mask[indices] = 0
        self.filter_data(indices = indices)
        #    print(len(self.A))


    def filter_data(self, indices = None):
        # This function applies a filter by only keeping triangle
        # in the list of vector indices

        print('--- Filtering ---')
        # Applying filter
        print('before', np.nanmin(self.start_time),np.nanmax(self.start_time))
        self.start_time = self.start_time[indices]
        print('after',np.nanmin(self.start_time),np.nanmax(self.start_time))
        self.end_time = self.end_time[indices]
        self.start_lat1 = self.start_lat1[indices]
        self.start_lat2 = self.start_lat2[indices]
        self.start_lat3 = self.start_lat3[indices]
        self.start_lon1 = self.start_lon1[indices]
        self.start_lon2 = self.start_lon2[indices]
        self.start_lon3 = self.start_lon3[indices]
        self.end_lat1 = self.end_lat1[indices]
        self.end_lat2 = self.end_lat2[indices]
        self.end_lat3 = self.end_lat3[indices]
        self.end_lon1 = self.end_lon1[indices]
        self.end_lon2 = self.end_lon2[indices]
        self.end_lon3 = self.end_lon3[indices]
        self.A = self.A[indices]
        self.div = self.div[indices]
        self.shr = self.shr[indices]
        self.vrt = self.vrt[indices]
        self.idx1 = self.idx1[indices]
        self.idx2 = self.idx2[indices]
        self.idx3 = self.idx3[indices]
        self.no   = self.no[indices]
        self.dudx = self.dudx[indices]
        self.dudy = self.dudy[indices]
        self.dvdx = self.dvdx[indices]
        self.dvdy = self.dvdy[indices]
        self.Mask = self.Mask[indices]
        self.day_flag = self.day_flag[indices]

    def mask_data(self, indices = None):
        # This function applies a filter by only keeping triangle
        # in the list of vector indices

        print('--- Filtering ---')
        # Applying filter
#        self.start_time[indices] = np.nan
#        self.end_time[indices] = np.nan
#        self.start_lat1[indices] = np.nan
#        self.start_lat2[indices] = np.nan
#        self.start_lat3[indices] = np.nan
#        self.start_lon1[indices] = np.nan
#        self.start_lon2[indices] = np.nan
#        self.start_lon3[indices] = np.nan
        self.A[indices] = np.nan
        self.div[indices] = np.nan
        self.shr[indices] = np.nan
        self.vrt[indices] = np.nan
#        self.idx1[indices] = np.nan
#        self.idx2[indices] = np.nan
#        self.idx3[indices] = np.nan
#        self.no[indices] = np.nan
        self.dudx[indices] = np.nan
        self.dudy[indices] = np.nan
        self.dvdx[indices] = np.nan
        self.dvdy[indices] = np.nan
        self.Mask[indices] = np.nan

    def Concatenate_data(self, Data2 = None):
        # This function concatenate the new data to
        # the current object
        # in the list of vector indices

        print('--- Concatenate new data to history ---')
        print('    :: lengths are: %s, %s' % (len(self.A), len(Data2.A)))

        self.start_time = np.append(self.start_time,Data2.start_time)
        self.end_time = np.append(self.end_time,Data2.end_time)

        self.start_lat1 = np.append(self.start_lat1, Data2.start_lat1)
        self.start_lat2 = np.append(self.start_lat2,Data2.start_lat2)
        self.start_lat3 = np.append(self.start_lat3,Data2.start_lat3)
        self.start_lon1 = np.append(self.start_lon1,Data2.start_lon1)
        self.start_lon2 = np.append(self.start_lon2,Data2.start_lon2)
        self.start_lon3 = np.append(self.start_lon3,Data2.start_lon3)

        self.end_lat1 = np.append(self.end_lat1,Data2.end_lat1)
        self.end_lat2 = np.append(self.end_lat2,Data2.end_lat2)
        self.end_lat3 = np.append(self.end_lat3,Data2.end_lat3)
        self.end_lon1 = np.append(self.end_lon1,Data2.end_lon1)
        self.end_lon2 = np.append(self.end_lon2,Data2.end_lon2)
        self.end_lon3 = np.append(self.end_lon3,Data2.end_lon3)

        self.A = np.append(self.A,Data2.A)
        self.div = np.append(self.div,Data2.div)
        self.shr = np.append(self.shr,Data2.shr)
        self.vrt = np.append(self.vrt,Data2.vrt)

        self.idx1 = np.append(self.idx1,Data2.idx1)
        self.idx2 = np.append(self.idx2,Data2.idx2)
        self.idx3 = np.append(self.idx3,Data2.idx3)
        self.no = np.append(self.no,Data2.no)
        self.day_flag = np.append(self.day_flag,Data2.day_flag)

        self.dudx = np.append(self.dudx,Data2.dudx)
        self.dudy = np.append(self.dudy,Data2.dudy)
        self.dvdx = np.append(self.dvdx,Data2.dvdx)
        self.dvdy = np.append(self.dvdy,Data2.dvdy)

        self.Mask = np.append(self.Mask,Data2.Mask)
        #self.Mask = self.A.copy()
        #self.Mask[:] = 1
        # Filtering indices where all three points (entire triangle) is within the radius
       # indices = [i for i in range(len(self.A)) if  self.A[i] < 20000.0**2.0 ]
       # self.Mask[indices] = 0
