"""
Author: Lekima Yakuden
GitHub: LekiYak

--------------------------------------------------------------------------------
Configuration file for data / netCDF analysis tools
--------------------------------------------------------------------------------

This file contains functions for loading and processing user options, raw data, and netCDF files.
"""

import configparser
import os
import sys
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset
import haversine as hs
import pandas as pd
import pyproj

def stb(s):
    if s in ['yes','Yes','true','True']:
         return True
    elif s in ['no','False','No','false']:
         return False
    else:
         raise ValueError

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
    in_proj = pyproj.Proj(init='epsg:4326')

    # Output projection (EPSG 3413, Polar Stereographic)
    out_proj = pyproj.Proj('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ', preserve_units=True)

    # Transform function
    x, y = np.array(pyproj.transform(in_proj, out_proj, lon, lat))
    return x, y
