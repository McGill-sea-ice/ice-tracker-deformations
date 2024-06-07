'''
Author: Beatrice Duval (bdu002)

-----------------------------------------
Configuration script for data processing
-----------------------------------------

- Retrieves configuration arguments from namelist.ini
- Selects which raw datasets are to be processed using namelist.ini arguments
- Obtain paths to which data files of all stages of data processing are to be stored

'''

# Loading from default packages
import configparser
import os
import re
import sys
from datetime import datetime, timedelta
import pandas as pd
import matplotlib.pyplot as plt


'''
_______________________________________________________________________
DEFINE CONFIG FUNCTIONS
'''

class dataset_config:
    '''
    This object is used to make the dataset production configuration,
    get the namilist arguments and list the ASITS output files to be processed.
    '''

    def __init__(self):

        # Retrieve the path of the src folder, i.e. the current directory
        self.srcPath = os.path.dirname(os.path.realpath(__file__))

        # Choose the default or user file
        def_fname = '/namelist.def'
        usr_fname = '/namelist.ini'

        if os.path.exists(self.srcPath + usr_fname):
            self.fname = usr_fname
        elif os.path.exists(self.srcPath + def_fname):
            print('--- Using default parameters namelist.def ---')
            print('--- Create the namelits.ini file to define user parameters ---')
            self.fname = def_fname
        else:
            print('/!/ No config file found! /!/')

        #----------------
        # Get the namelist
        #----------------

        # Read the namelist.ini file using configparser
        config = configparser.ConfigParser()
        config.read(self.srcPath + self.fname)

        # Return a dictionnary object
        config_dict = {sect: dict(config.items(sect)) for sect in config.sections()}

        # Iterate through the dictionnary to change strings to bools
        for key in config_dict:
            if isinstance(config_dict[key], dict):
                for keyy in config_dict[key]:
                    if isinstance(config_dict[key][keyy], str):
                        config_dict[key][keyy] = self.stb(config_dict[key][keyy])
            elif isinstance(config_dict[key], str):
                config_dict[key] = self.stb(config_dict[key])

        self.namelist = config_dict


        #----------------
        # Make list of daily start-end time pairs (to itterate from)
        #----------------


        Date_options = self.namelist['Date_options']
        options = self.namelist['options']
        satellite = self.namelist['Metadata']['icetracker']

        #Get the start and end dates
        start_year  = str(Date_options['start_year'])
        start_month = str(Date_options['start_month'])
        start_day   = str(Date_options['start_day'])
        end_year    = str(Date_options['end_year'])
        end_month   = str(Date_options['end_month'])
        end_day     = str(Date_options['end_day'])

        sDate = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
        eDate = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

        # Get the number of days in the production period
        date_range = (eDate - sDate)
        date_range = date_range.days

        # Allowing *min_date* to be updated and initializing time difference object
        start_date = sDate
        dtime = timedelta(hours=24)

        # Creating list containing tuples of date ranges [(dt1, dt2), (dt2, dt3) ...]
        self.date_pairs = []
        for i in range(date_range):
            end_date = start_date + dtime
            self.date_pairs.append((start_date, end_date))
            start_date = end_date

    def get_daily_SARpair_data(self):
        """
        Filters through the data files located in the input data folder to make
        a list of ASITS outputs with first image acquisition time in the specified date.


        Output: raw_paths -- List of ASITS output files to be processed
        """

        # Retrieve the IO and Date_options sections
        IO = self.namelist['IO']
        Date_options = self.namelist['Date_options']
        Metadata = self.namelist['Metadata']

        # Get the start and end date of given day
        start_year  = str(Date_options['start_year'])
        start_month = str(Date_options['start_month'])
        start_day   = str(Date_options['start_day'])
        end_year    = str(Date_options['end_year'])
        end_month   = str(Date_options['end_month'])
        end_day     = str(Date_options['end_day'])
        timestep    = int(Date_options['timestep'])
        tolerance   = int(Date_options['tolerance'])

        sDate = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
        eDate = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

        # Initializing file list and date count variables
        raw_paths = []

        #Fetch the path to the SAR-derived sea ice motion files
        date_path = IO['data_folder']
        satellite = Metadata['icetracker']

        if satellite == 'RCMS1':
            sat_list = ['rcm/','s1/']
        elif satellite == 'S1':
            sat_list = ['s1/']
        elif satellite == 'RCM':
            sat_list = ['rcm/']
        else:
            sys.exit("Oh! Original, but satellite data %s is not defined!!" % satellite )


        # List the pair files that correspond to time interval
        for sat_type in sat_list:
            for year in range(int(start_year), int(end_year)+1):

                if not(os.path.exists(date_path + sat_type + str(year) + '/')):
                    print('No data for '+ sat_type + ' in ' + str(year) )
                else:
                    data_path = date_path + sat_type + str(year) + '/'

                    #--------------------------------------------------------------------------
                    # Seek pairs that are in the specific date
                    #---------------------------------------------------------------------------

                    # Filtering data files by date
                    listfiles = os.listdir(data_path)
                    for filename in sorted(listfiles):
                        # Extracting initial and final dates from data file names

                        iDate = datetime.strptime(filename[6:20], '%Y%m%d%H%M%S')
                        fDate = datetime.strptime(filename[21:35], '%Y%m%d%H%M%S')

                        if iDate <= eDate and iDate > sDate and (fDate-iDate) <= timedelta(days=6):
                            raw_paths.append(data_path + '/' + filename)

        return raw_paths


    def stb(self, s):
        if s in ['yes','Yes','true','True']:
             return True
        elif s in ['no','False','No','false']:
             return False
        else:
             return s
