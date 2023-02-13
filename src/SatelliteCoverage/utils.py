"""
Author: Lekima Yakuden
GitHub: LekiYak

--------------------------------------------------------------------------------
Configuration file for data / netCDF analysis tools
--------------------------------------------------------------------------------

This file contains functions for loading and processing user options, raw data, and netCDF files.
"""

from netCDF4 import Dataset
from datetime import datetime, timedelta
import pyproj
import numpy as np

# Function to change strings to bools
def stb(s):
    if s in ['yes','Yes','true','True']:
         return True
    elif s in ['no','False','No','false']:
         return False
    else:
         return s

# Converts lat/lon to the north pole stereographic projection (EPSG 3413)
def convert_to_grid(lon, lat):
    """
    WARNING: INPUT IS LON, LAT (X, Y), AND THE OUTPUT IS ALSO X, Y

    Takes in a point in EPSG 4326 (Lat/Lon) and transforms said point to
    a polar stereographic projection (EPSG 3413).

    lon, lat {float} -> x, y {float}

    INPUTS:
    lat, lon -- Float values representing a point in WGS 84

    OUTPUTS:
    x, y -- Float values representing the transformed point in EPSG 3413
    """

    # Input projection (lat/lon)
    in_proj = pyproj.CRS('epsg:4326')

    # Output projection (EPSG 3413, Polar Stereographic)
    out_proj = pyproj.CRS('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ', preserve_units=True)

    # Transform function
    y, x = np.array(pyproj.Transformer.from_crs(in_proj, out_proj).transform(lat, lon))
    return y, x

# Converts lat/lon to the north pole stereographic projection (EPSG 3413)
def convert_from_grid(x, y):
    """
    WARNING: INPUT IS LON, LAT (X, Y), AND THE OUTPUT IS ALSO X, Y

    Takes in a point in polar stereographic  and transforms said point to
    a projection EPSG 4326 (Lat/Lon) .

    x, y {float} -> lon, lat {float}

    INPUTS:
    x, y -- Float values representing the transformed point in EPSG 3413

    OUTPUTS:
    lat, lon -- Float values representing a point in WGS 84

    """

    # output projection (lat/lon)
    out_proj = pyproj.CRS('epsg:4326')

    # input projection (EPSG 3413, Polar Stereographic)
    in_proj = pyproj.CRS('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ', preserve_units=True)

    # Transform function
    lat, lon = np.array(pyproj.Transformer.from_crs(in_proj, out_proj).transform(x, y))
    return lat, lon

"""
netCDF Analysis
"""
# Converts desired time range to seconds for netCDF analysis
def date_to_seconds(path:str, start_year, start_month, start_day, end_year, end_month, end_day):
    """
    This function takes in a netCDF's reference time and start and end times of each triangle
    and calculates the seconds elapsed between the former and latter. This is done to filter
    the data loaded from the netCDF by time, as the start and end times in the netCDF are stored
    as "seconds after the reference time".

    INPUTS:
    path -- String of path to the netCDF file the data will be read from {str}

    start_year -- Start year of the user-specified period (Hereinafter the "period") {str}
    start_month -- Start month of the period {str}
    start_day -- Start day of the period {str}

    end_year -- End year of the period {str}
    end_month -- End month of the period {str}
    end_day -- End day of the period {str}
    """

    # Open dataset
    ds = Dataset(path, mode='r')

    # Fetch reference time (start timestamp)
    reftime = ds.getncattr('referenceTime')

    # Closing dataset
    ds.close()

    # Converting reference time from string to datetime object
    reftime = reftime[0:10]
    reftime = datetime.strptime(reftime, '%Y-%m-%d')

    # Converting start and end date strings to datetime objects
    start_date = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
    end_date = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

    # Taking the difference between the start and end times of each triangle
    # and converting it to seconds
    start_ref_dt = (start_date - reftime).total_seconds()
    end_ref_dt = (end_date - reftime).total_seconds()

    return start_ref_dt, end_ref_dt

# Converts date in seconds from netCDF's reference time to datetime object
def seconds_to_date(path:str, date_sec_in):
    # Open dataset
    ds = Dataset(path, mode='r')

    # Fetch reference time (start timestamp)
    reftime = ds.getncattr('referenceTime')

    # Closing dataset
    ds.close()

    # Converting reference time from string to datetime object
    reftime = reftime[0:10]
    reftime = datetime.strptime(reftime, '%Y-%m-%d')

    date_out = reftime+timedelta(seconds=int(date_sec_in))

    return date_out

def get_prefix(config=None):

    Date_options = config['Date_options']
    sy = str(Date_options['start_year'])
    ey = str(Date_options['end_year'])
    sm = str(Date_options['start_month'])
    em = str(Date_options['end_month'])
    sd = str(Date_options['start_day'])
    ed = str(Date_options['end_day'])
    ts = str(Date_options['timestep'])
    to = str(Date_options['tolerance'])
    it = str(config['Metadata']['icetracker'])

    prefix = it + '_' + sy + sm + sd + '_' + ey + em + ed + '_dt' + ts + '_tol' + to

    if 'options' in config:
        if config['options']['area_filter'] :
            cla = str(config['options']['centre_lat'])
            clo = str(config['options']['centre_lon'])
            rad = str(config['options']['radius'])

            prefix = prefix + '_filt_lat'+ str(cla) + '_lon' + str(clo) + '_rad' + str(rad)

    return prefix
