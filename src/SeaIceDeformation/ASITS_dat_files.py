"""
Authors: Mathieu Plante, Amelie Bouchat, Damien Ringeisen, Lekima Yakuden, Beatrice Duval
GitHub: mathieuslplante

--------------------------------------------------------------------------------
Python object to load and process ECCC-ASITS output files
--------------------------------------------------------------------------------
This object contains information from motion vectors from the ASITS.
It includes functions to load a specific output file, get the motion data.

"""

# Loading from default packages
import os
import sys
parent = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,parent)
from netCDF4 import Dataset
import numpy as np
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import csv



# Loads netCDF data
class Load_ASITS_dat_file:

    def __init__(self,FileName= None, config=None):
        """
        This function reads and loads data from the ASITS output files.
        This creates an object containing the data fields.

        INPUTS:
        FileName = path to the .dat file with SIM data.
        config   = namelist object, from the python config parser object.

        OBJECT CARACTERISTICS:

        """

        path_raw = FileName
        # Create a data frame
        df = pd.read_csv(path_raw, sep='\s\s+', engine='python')
        self.sigx = 200
        # Retrieve data points
        self.sLat    = df['sLat']    # Starting latitudes
        self.sLon    = df['sLon']    # Starting longitudes
        self.eLat    = df['eLat']    # Ending latitudes
        self.eLon    = df['eLon']    # Ending longitudes
        self.startX  = df['startX']  # Starting X positions (px)
        self.startY  = df['startY']  # Starting Y positions (px)
        self.endX    = df['endX']    # Ending X positions (px)
        self.endY    = df['endY']    # Ending Y positions (px)
        self.dispX   = df['dispX']   # X displacement (px)
        self.dispY   = df['dispY']   # Y displacement (px)
        del df

        self.FileName = path_raw

        # Retrieve the starting and ending times and compute the time interval (days)
        self.start = datetime.datetime(int(self.FileName[-35:-31]), # Year
                                  int(self.FileName[-31:-29]), # Month
                                  int(self.FileName[-29:-27]), # Day
                                  int(self.FileName[-27:-25]), # Hour
                                  int(self.FileName[-25:-23]), # Minute
                                  int(self.FileName[-23:-21])) # Second

        self.end   = datetime.datetime(int(self.FileName[-20:-16]), # Year
                                  int(self.FileName[-16:-14]), # Month
                                  int(self.FileName[-14:-12]), # Day
                                  int(self.FileName[-12:-10]), # Hour
                                  int(self.FileName[-10:-8]),  # Minute
                                  int(self.FileName[-8:-6]))   # Second

        self.dt = (self.end-self.start).total_seconds()/86400.0


        # Get the start and end time with respect to 0h00 of the ref date (daily netcdf)
        Date_options = config['Date_options']
        YYYY = Date_options['start_year']
        MM = Date_options['start_month']
        DD = Date_options['start_day']
        self.year = str(YYYY)
        self.month = str(MM)
        self.day = str(DD)

        self.refTime = datetime.datetime(int(YYYY), # Year
                       int(MM),   # Month
                       int(DD),   # Day
                       0, 0, 0)   # Hour, minute, second

        self.sTime = (self.start- self.refTime).total_seconds()
        self.eTime = (self.end-self.refTime).total_seconds()

        self.sat = (0)*('rcm' in FileName)+(1)*('s1' in FileName)

