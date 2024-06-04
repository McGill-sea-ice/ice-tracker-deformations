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

# Loading from other files
import SeaIceDeformation.utils_get_data_paths as get_data_paths
#from SatelliteCoverage.utils import stb

'''
_______________________________________________________________________
DEFINE ERROR CLASS
'''


# Create a class of errors for datasets selection
class datasetSelectionError(Exception):
    pass

'''
_______________________________________________________________________
DEFINE CONFIG FUNCTIONS
'''

def stb(s):
    if s in ['yes','Yes','true','True']:
         return True
    elif s in ['no','False','No','false']:
         return False
    else:
         return s


def get_config_args():
    ''' None -> ConfigParser

    Function that reads the namelist.ini file and returns a
    ConfigParser object, assuming the namelist.ini file is
    located under the current directory.

    '''

    # Retrieve the path of the src folder, i.e. the current directory
    srcPath = os.path.dirname(os.path.realpath(__file__))

    # Read the namelist.ini file using configparser
    config = configparser.ConfigParser()

    # Choose the default or user file
    def_fname = '/namelist.def'
    usr_fname = '/namelist.ini'
    if os.path.exists(srcPath + usr_fname):
        fname = usr_fname
    elif os.path.exists(srcPath + def_fname):
        print('--- Using default parameters namelist.def ---')
        print('--- Create the namelits.ini file to define user parameters ---')
        fname = def_fname
    else:
        print('/!/ No config file found! /!/')

    config.read(srcPath + fname)

    # Return a dictionnary object
    config_dict = {sect: dict(config.items(sect)) for sect in config.sections()}

    # Iterate through the dictionnary to change strings to bools
    for key in config_dict:
        if isinstance(config_dict[key], dict):
            for keyy in config_dict[key]:
                if isinstance(config_dict[key][keyy], str):
                    config_dict[key][keyy] = stb(config_dict[key][keyy])
        elif isinstance(config_dict[key], str):
            config_dict[key] = stb(config_dict[key])

    return config_dict



# Divides raw data into intervals specified by the user
# def divide_intervals(raw_paths, max_date, min_date, interval):
def divide_intervals(config=None):
    """
    Divides delta-t filtered data into chunks (intervals) of *interval* hours for
    processing.
    INPUTS:
    raw_paths -- List of data file paths with the desired delta t {list}
    Date_options -- ConfigParser object to determine start/end dates of intervals as specififed by user
    Options -- ConfigParser object to determine the interval length as specififed by user
    OUTPUTS:
    interval_list -- List of n lists, where n is the number of full intervals which
                     fit in the min and max dates of coverage. Each nested list (n th)
                     contains the paths to files which share a date range (even partially)
                     with the n th interval.
    date_pairs -- List of tuples, each containing the start and end dates of the n th interval (datetime
                  objects in tuples, all in a list)
    """

    #raw_paths = config['raw_paths']
    Date_options = config['Date_options']
    options = config['options']
    satellite = config['Metadata']['icetracker']

    start_year  = str(Date_options['start_year'])
    start_month = str(Date_options['start_month'])
    start_day   = str(Date_options['start_day'])
    end_year    = str(Date_options['end_year'])
    end_month   = str(Date_options['end_month'])
    end_day     = str(Date_options['end_day'])
    interval    = options['interval']
    print(start_year,start_month,start_day)

    # Concatenate start and end dates
    sDate = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
    eDate = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

    # Converting date range to timedelta object (hours)
    date_range = (eDate - sDate)
    date_range = date_range.days * 24

    # Counting number of full intervals over date range {int}
    interval_count = date_range // int(interval)

    # Allowing *min_date* to be updated and initializing time difference object
    min_date_it = sDate
    dtime = timedelta(hours=int(interval))

    # Creating list containing tuples of date ranges [(dt1, dt2), (dt2, dt3) ...]
    date_pairs = []
    for i in range(interval_count):
        date_pairs.append((min_date_it, min_date_it + dtime))
        min_date_it = min_date_it + dtime

    return date_pairs

def filter_data(config=None):
    """
    Filters through the data files located in 'data_path' using the user
    options in 'options.ini'. Outputs a list of paths to data files which
    satisfy the user's criteria.

    Automatically changes date range to match data availability. i.e. if the user specifies
    a date range between 01-11-2020 and 01-06-2021, but data is only available from
    05-11-2020 and 24-05-2021, dates to be processed will be set to the latter, and
    the user will be notified (Line

    INPUTS:
    start_year -- Starting year YYYY {str}
    start_month -- Starting month MM {str}
    start_day -- Starting day DD {str}

    end_year -- Ending year YYYY {str}
    end_month -- Ending month MM {str}
    end_day -- Ending day DD {str}

    timestep -- Desired timestep in hours {str}
    tolerance -- Number of hours around timestep that will be filtered through {str}

    data_path -- Path to directory containing data files {str}

    OUTPUTS:
    raw_paths -- List of file paths {list}
    """

    # Retrieve the IO and Date_options sections
    IO = config['IO']
    Date_options = config['Date_options']
    Metadata = config['Metadata']

    start_year  = str(Date_options['start_year'])
    start_month = str(Date_options['start_month'])
    start_day   = str(Date_options['start_day'])
    end_year    = str(Date_options['end_year'])
    end_month   = str(Date_options['end_month'])
    end_day     = str(Date_options['end_day'])
    timestep    = int(Date_options['timestep'])
    tolerance   = int(Date_options['tolerance'])
    print(start_year,start_month,start_day)

    # Concatenate start and end dates
    sDate = datetime.strptime(start_year + start_month + start_day, '%Y%m%d')
    eDate = datetime.strptime(end_year + end_month + end_day, '%Y%m%d')

    print(sDate,eDate)
    # Set delta t tolerance
    upper_timestep = timedelta(hours=(int(timestep) + int(tolerance)))
    lower_timestep = timedelta(hours=(int(timestep) - int(tolerance)))

    # Initializing file list and date count variables
    raw_paths = []
    min_date = datetime(3000, 12, 25)
    max_date = datetime(1000, 12, 25)

    # List of summer months (Won't include data from these months)
    #summer_months = [6, 7, 8, 9, 10]
    summer_months = [0]

    #Fetch the path to the SAR-derived sea ice motion files
    date_path = IO['data_folder']
    satellite = Metadata['icetracker']

    if satellite == 'RCMS1':
        sat_list = ['rcm/','s1/']
    elif satellite == 'S1':
        sat_list = ['s1/']
    elif satellite == 'RCM':
        sat_list = ['rcm/']
    elif satellite == 'RCM_new':
        sat_list = ['rcm_corrected/']
    else:
        sys.exit("Oh, original, but satellite data %s is not defined!!" % satellite )

    #--------------------------------------------------------
    # This below should be moved to data analysis tools.
    # Plotting the dt distribution for the images pair
    if config['Processing_options']['viz_tstp_dist']:
        # options (to have in the namelist?)
        dens = False
        log = False
        bin_step = 1
        ys = 2017
        ye = 2022
        alp = 0.5
        hsty = 'stepfilled' # 'bar', 'barstacked', 'step', 'stepfilled'
        bins = range(0,100+bin_step,bin_step) # bins initialization

        # Dictionnaries initialization
        interv={}
        interv['all'] =[]

        # figure initialization
        plt.figure()

        # Read the files names and create the data lists
        for sat in sat_list:
            interv[sat]=[]
            for year in range(int(ys), int(ye)+1): # all the years
                if os.path.exists(date_path + sat + str(year) + '/'):
                    data_path = date_path + sat + str(year) + '/'
                    listfiles = os.listdir(data_path)
                    for filename in sorted(listfiles):
                        # Extracting initial and final dates from data file names
                        iDate = datetime.strptime(filename[6:20], '%Y%m%d%H%M%S')
                        fDate = datetime.strptime(filename[21:35], '%Y%m%d%H%M%S')
                        diff = (fDate - iDate).total_seconds() / 3600.0 # timesteps in hour
                        interv[sat] += [float(diff)]

            interv['all'] += interv[sat]

        interv_list = [ interv[name] for name in sat_list ]
        label_list = sat_list
        if hsty != 'barstacked':
            interv_list += [interv['all']]
            label_list += ['all']
        hist = plt.hist(interv_list, bins=bins, density=dens, label=label_list, histtype=hsty, alpha = alp, log=log)

        # the figure bells and whistles
        plt.legend()
        if dens:
            plt.ylabel('PDF of the intervals')
        else :
            plt.ylabel('Number of intervals')
        plt.xlabel('Interval length [h]')
        plt.xticks(range(0,100,12))
        plt.title('Timestep between SAR images for S1 and RCM between {} and {} '.format(ys,ye))

        # saving the figure
        # Create the directory if it does not exist already
        output_folder = config['IO']['output_folder']  + '/' + config['IO']['exp'] + '/figs/'
        os.makedirs(output_folder, exist_ok=True)
        print('Saving timestep distribution figure at ' + output_folder + 'RCMS1_{}_{}_bin{}_dt_hist.png'.format(ys,ye,bin_step))
        plt.savefig(output_folder + 'RCMS1_{}_{}_bin{}_dt_hist.png'.format(ys,ye,bin_step), bbox_inches='tight', dpi=600)

        # Saving the data for after
        data={}
        data['bins_start'] = hist[1][:-1]
        data['bins_end'] = hist[1][1:]
        for i in range(len(sat_list)):
            name = label_list[i]
            data[name] = hist[0][i]
        df = pd.DataFrame(data, columns=['bins_start','bins_end'] + label_list)
        print('Saving timestep distribution data at ' + output_folder + 'RCMS1_{}_{}_bin{}_dt_hist.pkl'.format(ys, ye, bin_step))
        df.to_pickle(output_folder + 'RCMS1_{}_{}_bin{}_dt_hist.pkl'.format(ys, ye, bin_step))
        # df.to_csv(output_folder + 'RCMS1_{}_{}_bin{}_dt_hist.csv'.format(ys, ye, bin_step), index=False)
    #---------------------------------------------------------------------------------------------------


    # List the pair files that correspond to time interval
    for sat_type in sat_list:
        for year in range(int(start_year), int(end_year)+1):

            if not(os.path.exists(date_path + sat_type + str(year) + '/')):
                print('No data for '+ sat_type + ' in ' + str(year) )
            else:
                data_path = date_path + sat_type + str(year) + '/'

                #--------------------------------------------------------------------------
                # listing the files in the folder and adding the pairs if in the right dates
                #     This could be made a function if used when a IO class is made
                #---------------------------------------------------------------------------

                # Filtering data files by date
                listfiles = os.listdir(data_path)
                for filename in sorted(listfiles):
                    if satellite == 'RCM_new':
                        txt = filename.split('_')
                        iDate = datetime.strptime("%s%s" % (txt[5],txt[6]), '%Y%m%d%H%M%S')
                        fDate = datetime.strptime("%s%s" % (txt[15],txt[16]), '%Y%m%d%H%M%S')
                        mean_Date = iDate + (fDate - iDate)/2.0
                    else:
                        # Extracting initial and final dates from data file names
                        iDate = datetime.strptime(filename[6:20], '%Y%m%d%H%M%S')
                        fDate = datetime.strptime(filename[21:35], '%Y%m%d%H%M%S')
                        mean_Date = iDate + (fDate - iDate)/2.0


                    # Checking if all files from iDate to fDate will be loaded (timestep == '0')
                    if timestep != '0':

                        ##---- Use the part below to select data based on valid time-----
                        #if iDate < eDate and fDate > sDate and lower_timestep <= (fDate-iDate) <= upper_timestep and iDate.month not in summer_months:
                        ##---------------------------------------------------------------

                        ##---- Use the part below to select data based on start time-----
                        if iDate <= eDate and iDate > sDate and lower_timestep <= (fDate-iDate) <= upper_timestep and iDate.month not in summer_months:
                        ##---------------------------------------------------------------

                            raw_paths.append(data_path + '/' + filename)

                            # I don't think this is needed
                            # Updating date tracker
                            if iDate < min_date:
                                min_date = iDate
                            if fDate > max_date:
                                max_date = fDate

                    # This below never occurs
                    elif timestep == '0':
                        # Filtering by date range only
                        if iDate < eDate and fDate > sDate and iDate.month not in summer_months:
                            raw_paths.append(data_path + '/' + filename)

                            # Updating date tracker
                            if iDate < min_date:
                                min_date = iDate
                            if fDate > max_date:
                                max_date = fDate

    # Notifying user of date range change
    if sDate != min_date or eDate != max_date:
        print(f"Start and end dates of data updated to {min_date} and {max_date}")

    return raw_paths


def get_datapaths(config=None):
    ''' (str, str, str, str) -> dict[str, Any]

    Function that creates lists that store data file paths for every dataset
    and for every stage of data processing (triangulation, conversion and
    calculations), as well as the output netcdf file path.

    Returns a dictionnary of the lists of .csv file paths (one per stage of data
    processing) and of the output netcdf file path.

    Every raw data file (e.g. pairs_20200320010317_20200401010317_1.dat) is
    associated to a triangulated (e.g. tri_20200320010317_20200401010317_1.csv),
    and calculated (e.g. calc_20200320010317_20200401010317_1.csv)
    .csv file.

    The output netcdf file combines the output data of all processed datasets
    (listed in the input *raw_paths*). Its name indicates the common starting times
    of all datasets (e.g. RCMS1SID_20200301_dx.nc).

    All output files are stored under output_path/exp, and under a folder
    associated to the stage of dataprocessing.

    For example, given a raw data file 'pair_A.dat', the output triangulation file
    will be stored under 'output_path/exp/02_triangulated/tri_A.csv', the output
    converted file under 'output_path/exp/03_converted/conv_A.csv', the output
    calculations files under 'output_path/exp/04_calculations/calc_A.csv', and the
    output netcdf file under 'output_path/exp/05_output/RCMS1SID_YYYYMMDD_dx.nc',
    where YYYY, MM and DD refer to the input *start_year*, *start_month* and *start_day*,
    respectively.

    Keyword arguments: \\
    raw_paths   -- list of raw datasets' absolute paths \\
    output_path -- path in which all output data files are to be stored \\
    exp         -- name of the experiment
    start_year  -- starting year of the raw datasets listed in raw_paths \\
    start_month -- starting month of the raw datasets listed in raw_paths \\
    start_day   -- starting day of the raw datasets listed in raw_paths \\
    icetracker -- sea-ice motion tracker used for calculations \\
    '''

    # As this is cleaned, we would like to avoid saving data
    # every step. Some of these paths will not be needed.

    raw_paths = config['raw_paths']

    # Retrieve the IO and Date_options sections
    IO = config['IO']
    Date_options = config['Date_options']
    Metadata = config['Metadata']

    # Raise an error if the input raw paths list is empty
    if raw_paths ==  []:
        raise datasetSelectionError('The list of raw datasets to process is empty.')

    # Initialize lists of data paths for the subsequent stages of data processing
    triangulated_paths  = []
    calculations_paths  = []

    # Iterate through all input raw file paths
    for raw in raw_paths:

        # Retrive the raw filename
        raw_filename = os.path.basename(raw)

        # For each raw file path, find the appropriate path for each stage of data processing
        tri  = get_data_paths.get_triangulated_csv_path(raw_filename, IO['output_folder'], IO['exp'])  # path for triangulated data file
        cal  = get_data_paths.get_calculations_csv_path(raw_filename, IO['output_folder'], IO['exp'])  # path for calculated data file

        # Append the file paths to the file path lists
        triangulated_paths.append(tri)
        calculations_paths.append(cal)

    # Find the appropriate path for the output netcdf file
    nc_output_path = get_data_paths.get_output_nc_path(IO, Date_options, Metadata)

    # Create the output data paths dictionnary
    # we are using to process data
    data_paths =  { 'raw': raw_paths,
                    'triangulated': triangulated_paths,
                    'calculations': calculations_paths,
                    'nc_output': nc_output_path,
                    'output_folder': IO['output_folder']+IO['exp']}

    return data_paths

