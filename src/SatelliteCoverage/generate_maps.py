import os
import pandas as pd
from datetime import date
import numpy as np
import matplotlib.pyplot as plt


def list_dates(data_folder):

    # Grab names of all files in directory
    dirList = os.listdir(data_folder)

    # Dump filenames into dataframe and sort by ascending date
    filenameDf = pd.DataFrame(dirList, columns=['Filename'])
    filenameDf = filenameDf.sort_values(by=['Filename'])
    
    # Converting YYYYMMDDHHMMSS strings to dates and calcuating delta T in hours
    filenameDf['deltaT'] = (pd.to_datetime(filenameDf['Filename'].str[21:35], format='%Y%m%d%H%M%S') -
                            pd.to_datetime(filenameDf['Filename'].str[6:20], format='%Y%m%d%H%M%S'))

    filenameDf['deltaT'] = filenameDf['deltaT'] / np.timedelta64(1, 'h')

    return filenameDf

def plotting_histogram(ds):
    """
    Pandas Series -> histogram.png

    Plots histogram of series containing delta T's and their frequency

    REQUIRES INPUT OF filenameDf['deltaT'] from list_dates()
    """
    upper_lim = 100
    plt.hist(ds, bins=np.arange(0, upper_lim, 1))
    plt.xticks(np.arange(0, upper_lim, 6))
    plt.yticks(np.arange(0, 5000, 300))
    plt.savefig('deltaT_96hrs_6hrbins.png')

#plotting_histogram(list_dates())

def filter_df(df, duration):
    filtered_df = df.loc[(df['deltaT'] < duration + 3) & (df['deltaT'] > duration - 3)]

    filtered_df = filtered_df.sort_values(by=['deltaT'])

    raw_paths = filtered_df['Filename'].tolist()

    return raw_paths

#print(filter_df(list_dates('./data/02_coverage/S1/S1_combined'), 24)[0])
