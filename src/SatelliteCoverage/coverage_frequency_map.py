import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import utils_load_grid as grid

from config import *

def compile_data(raw_paths):
    """ 
    Compiles all data points in the list of data files (raw_paths) into a dataframe (df)

    Returns the dataframe with all datapoints' starting lat and lon in columns.

    INPUTS:
    raw_paths -- List of paths to data files {List}

    OUTPUTS: 
    
    df -- Pandas dataframe with the following columns: Index  sLat  sLong {Dataframe}
    """

    df = pd.DataFrame()

    # Appending each file's datapoints to the dataframe
    for filepath in raw_paths:
        
        # Initialize temporary dataframe
        temp_df = pd.DataFrame()

        # Reading datapoints into temporary dataframe
        temp_df = pd.read_csv(filepath, sep='\s\s+', engine='python', usecols = ['sLat','sLon'])

        df = df.append(temp_df)

    return df

def count_duplicates(df):
    """
    Takes dataframe with sLat and sLon values and counts duplicates, returns dataframe
    with coordinates and frequency count.

    INPUT:
    df -- Pandas dataframe with sLat and sLong in two respective columns, representing data points

    OUTPUT:
    df -- Input dataframe with frequency column added
    
    """

    # Pivot df onto its side to count duplicate Lat & Lon entires
    df = df.pivot_table(columns=['sLat', 'sLon'], aggfunc='size')

    # Add size Multi index as column named 'frequency'
    df = df.reset_index()  
    
    df.rename(columns = {0:'frequency'}, inplace = True)

    return df

def visualise_coverage_old(df):

    print("Plotting frequency map")

    # Set the matplotlib projection and transform
    #proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
    proj=ccrs.NorthPolarStereo(central_longitude=0)
    trans = ccrs.Geodetic()
    

    # Initialize figure
    fig = plt.figure()

    # Add projection
    ax = fig.add_subplot(111, projection = proj)
    #ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())
    ax.set_extent((-2500000, 2500000, -2500000, 2500000), ccrs.NorthPolarStereo())
    

    # Setting coordinates
    ax.scatter(x=df['sLon'], y=df['sLat'], transform= trans, s=0.1)

    # Set grid
    ax.gridlines()

    # Hide datapoints over land
    ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

    # Set title
    ax.set_title('Heatmap of datapoint frequency over the 2019-2020 Winter (S1)')
    
    # Save figure
    plt.savefig('TEST_ALL2.png')

    return

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
    out_proj = pyproj.Proj('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ', preserve_units=True)

    # Transform function
    x, y = pyproj.transform(in_proj, out_proj, lon, lat)
    return x, y

print(convert_to_grid(-76, 62))



config = read_config()

# Initializing more specific ConfigParser objects
IO = config['IO']
options = config['options']

data_path = IO['data_folder']

sYear = options['start_year']
sMonth = options['start_month']
sDay = options['start_day']

eYear = options['end_year']
eMonth = options['end_month']
eDay = options['end_day']

delta_t = options['delta_t']
tolerance = options['tolerance']

raw_list = filter_data(sYear, sMonth, sDay, eYear, eMonth, eDay, delta_t, tolerance, data_path)

df = compile_data(raw_list)

visualise_coverage_old(df)
