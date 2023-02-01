'''
Author: Beatrice Duval (bdu002)

------------------------------------------------------------------------------
Visualise deformations
------------------------------------------------------------------------------

Code that plots total sea-ice deformations, divergence and shear strain rates
using the datasets listed in a txt file (see config).

'''

import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np

import config
import utils_load_data as load_data
from tqdm import tqdm
from timeit import default_timer as dt

def visualise_deformations():

    '''
    _________________________________________________________________________________________
    INITIALIZE PLOTS
    '''

    # Set the matplotlib projection and transform
    proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
    trans = ccrs.Geodetic()

    # Initialize figures for total deformation (tot), divergence (I) and shear (II)
    fig_tot = plt.figure()
    fig_I = plt.figure()
    fig_II = plt.figure()
    fig_rot = plt.figure()

    # Initialize subplots
    ax_tot = fig_tot.add_subplot(111, projection=proj)
    ax_I = fig_I.add_subplot(111, projection=proj)
    ax_II = fig_II.add_subplot(111, projection=proj)
    ax_rot = fig_rot.add_subplot(111, projection=proj)

    # Create a list of axes to be iterated over
    ax_list = [ax_tot, ax_I, ax_II, ax_rot]

    for ax in ax_list:
        # Set the map extent in order to see the entire region of interest
        ax.set_extent((-3000000, 4000000, 8500000, 11500000), ccrs.AzimuthalEquidistant())

    # Iterate through data files in config
    for raw_path, triangulated_path, calculated_path in zip(config.data_paths['raw'], config.data_paths['triangulated'], config.data_paths['calculations']):

        try:
            '''
            _________________________________________________________________________________________
            LOAD DATA
            '''
            # Load the raw data
            raw_data    = load_data.load_raw(raw_path)
            sLat        = raw_data['sLat']       # Starting latitudes
            sLon        = raw_data['sLon']       # Starting longitudes

            # Load triangulated data
            triangulated_data = load_data.load_triangulated(triangulated_path)
            vertice_idx1      = triangulated_data['vertice_idx1']   # Vertex indices in raw csv file
            vertice_idx2      = triangulated_data['vertice_idx2']
            vertice_idx3      = triangulated_data['vertice_idx3']

            # Load calculated data
            calculated_data = load_data.load_calculations(calculated_path)
            eps_I           = calculated_data['eps_I']    # Divergence rate
            eps_II          = calculated_data['eps_II']   # Maximum shear strain rate
            eps_tot         = calculated_data['eps_tot']  # Total deformation rate
            dvdx            = calculated_data['dvdx']
            dudy            = calculated_data['dudy']
            rot             = dvdx - dudy

        except:
            print('Error in reading the file : ',raw_path,' continuing.')
            continue

        '''
        _________________________________________________________________________________________
        PLOT DEFORMATIONS
        '''

        # Stack the 3 vertices index arrays
        triangles = np.stack((vertice_idx1, vertice_idx2, vertice_idx3), axis=-1)

        # Plot total deformations, divergence and shear
        cb_tot = ax_tot.tripcolor( sLon, sLat, triangles, facecolors=eps_tot, transform=trans, cmap='plasma', vmin=0, vmax=0.1 )
        cb_I   = ax_I.tripcolor( sLon, sLat, triangles, facecolors=eps_I, transform=trans, cmap='coolwarm', vmin=-0.04, vmax=0.04 )
        cb_II  = ax_II.tripcolor( sLon, sLat, triangles, facecolors=eps_II, transform=trans, cmap='plasma', vmin=0, vmax=0.1 )
        cb_rot = ax_rot.tripcolor( sLon, sLat, triangles, facecolors=rot, transform=trans, cmap='coolwarm', vmin=-0.1, vmax=0.1 )

    '''
    _________________________________________________________________________________________
    ADD FEATURES TO THE PLOTS
    '''

    # Check if at least one dataset has been plotted,
    # ie if the figures have been plotted (cb_tot should
    # have been initialized in that case)
    try:
        cb_tot
    except NameError:
        cb_tot = None

    # Retrieve namelist arguments that define the dataset
    # that has been selected for processing
    Date_options = config.config['Date_options']
    start_year   = Date_options['start_year']
    start_month  = Date_options['start_month']
    start_day    = Date_options['start_day']
    duration     = Date_options['timestep']

    IO            = config.config['IO']
    output_folder = IO['output_folder']
    exp           = IO['exp']

    # All figures have been plotted
    if cb_tot is not None:

        # Create a list of colorbars and titles to be iterated over
        cb_list = [cb_tot, cb_I, cb_II, cb_rot]
        title_list = ['Total Deformation Rate $(Days^{-1})$', 'Divergence Rate $(Days^{-1})$', 'Shear Rate $(Days^{-1})$', 'Rotation Rate $(Days^{-1})$']

        # Iterate through all axes
        for ax, title, cb in zip(tqdm(ax_list), title_list, cb_list):
            # Add a colorbar
            plt.colorbar(cb, ax=ax)

            # Add a title
            ax.set_title(title + '\n' + start_year + '-' + start_month + '-' + start_day + ' -- over ' + str(duration) + ' day')

            # Add gridlines
            ax.gridlines()

            # Hide deformations over land
            ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

        '''
        _________________________________________________________________________________________
        SAVE PLOTS
        '''

        # Set a directory to store figures for the current experiment
        figsPath =  output_folder + '/' + exp + '/figs/'

        # Create the directory if it does not exist already
        os.makedirs(figsPath, exist_ok=True)

        # Create a prefix for the figure filenames
        prefix = start_year + start_month + start_day + '_' + str(duration)

        # Create the figure filenames
        tot_path   = figsPath + prefix + '_tot.png'
        div_path   = figsPath + prefix + '_div.png'
        shear_path = figsPath + prefix + '_shear.png'
        rot_path   = figsPath + prefix + '_rot.png'


        for fig, fig_path in zip([fig_tot, fig_I, fig_II, fig_rot], [tot_path, div_path, shear_path, rot_path]):

            # Check if the figures already exist. If they do, delete them.
            if os.path.isfile(fig_path):
                os.remove(fig_path)

            # Save the new figures
            fig.savefig(fig_path, bbox_inches='tight', dpi=600)

        plt.show()

    else:
        print('No figures have been created.')

def recreate_coordinates(start_lat1, start_lat2, start_lat3, start_lon1, start_lon2, start_lon3, start_id1, start_id2, start_id3):
    """
    This function takes in a list of start lats/lons and corresponding start lat/lon IDs (triangle vertices)
    and outputs corresponding, larger lists of start lats/lons (with 0s in some indices) that reflect the
    coordinates' placement in the original data files for use with ax.tripcolor

    INPUTS:
    start_lat1,2,3 -- Arrays of starting latitudes {np.array, list}
    start_lons1,2,3 -- Arrays of starting longitudes {np.array, list}
    start_id1,2,3 -- Array of starting IDs corresponding to start_lats1,2,3 and start_lons1,2,3 {np.array, list}

    OUTPUTS:
    new_lat -- Array of latitude values at the positions they were orignally in, in the data file
    new_lon -- Array of longitude values at the positions they were originally in, in the data file

    """

    # Combined list of start IDs
    start_ids = np.hstack((start_id1, start_id2, start_id3))

    # Skipping blank data files
    if len(start_ids) == 0:
        return 0, 0

    # Initializing new lists of coordinates
    new_lon, new_lat = ([np.nan] * (max(start_ids) + 1) for i in range(2))

    for i in range(len(start_lat1)):
        new_lat[start_id1[i]] = (start_lat1[i])
        new_lat[start_id2[i]] = (start_lat2[i])
        new_lat[start_id3[i]] = (start_lat3[i])

    for i in range(len(start_lon1)):
        new_lon[start_id1[i]] = start_lon1[i]
        new_lon[start_id2[i]] = start_lon2[i]
        new_lon[start_id3[i]] = start_lon3[i]

    return new_lat, new_lon


def plot_deformations_netdcf(data):
    """
    This function plots deformations from a netCDF file using matplotlib's ax.tripcolor.
    The function assumes that the netCDF was generated from src/SeaIceDeformation's M01_d03_compute_deformations.py.

    INUPTS:
    path -- Path to netCDF {str}

    OUTPUTS:
    None -- Saves divergence, shear, and vorticity plots to the output directory
    """

    """
    Preamble
    """

    # Set the matplotlib projection and transform
    proj = ccrs.NorthPolarStereo(central_longitude=0)
    trans = ccrs.Geodetic()

    # Initialize figures for total deformation (tot), divergence (I) and shear (II)
    fig_div = plt.figure()
    fig_shr = plt.figure()
    fig_vrt = plt.figure()

    # Initialize subplots
    ax_div = fig_div.add_subplot(111, projection=proj)
    ax_shr = fig_shr.add_subplot(111, projection=proj)
    ax_vrt = fig_vrt.add_subplot(111, projection=proj)

    # Create a list of axes to be iterated over
    ax_list = [ax_div, ax_shr, ax_vrt]

    for ax in ax_list:
        # Set the map extent in order to see the entire region of interest
        ax.set_extent((-4400000, 2500000, 3500000, -2500000), ccrs.NorthPolarStereo())

    timer_recreate = 0.
    timer_plot = 0.
    timer_tria = 0.
    timer_trans = 0.

    """
    Plotting
    """

    # Obtaining start index of data
    min_no = np.nanmin(data['no']) # Minimum (smallest) ID number
    min_index = np.nanmin(np.where(data['no'] == min_no)) # This value will be updated for each triangle

    # Iterating over all files (Unique triangulations)
    for i in tqdm(range(len(set(np.array(data['no']))))):

        # Here the file indice is zero-indexed, so the first file has the indice 0, regardless of its ID

        # Obtaining number of rows to iterate over
        file_length = np.count_nonzero(data['no'] == i + min_no)

        # Setting maximum index
        max_index = min_index + file_length

        # Arranging triangle vertices in array for use in ax.tripcolor
        triangles = np.stack((data['id_start_lat1'][min_index:max_index], data['id_start_lat2'][min_index:max_index],
                    data['id_start_lat3'][min_index:max_index]), axis=-1)

        # Filtering data range to that of the current "file"
        start_lat1_temp, start_lat2_temp, start_lat3_temp = data['start_lat1'][min_index:max_index], \
            data['start_lat2'][min_index:max_index], data['start_lat3'][min_index:max_index]

        start_lon1_temp, start_lon2_temp, start_lon3_temp = data['start_lon1'][min_index:max_index], \
            data['start_lon2'][min_index:max_index], data['start_lon3'][min_index:max_index]

        start_rec = dt()
        start_lat, start_lon = recreate_coordinates(start_lat1_temp, start_lat2_temp, start_lat3_temp,
                                                        start_lon1_temp, start_lon2_temp, start_lon3_temp,
                                                        data['id_start_lat1'][min_index:max_index],
                                                        data['id_start_lat2'][min_index:max_index],
                                                        data['id_start_lat3'][min_index:max_index])
        timer_recreate += dt() - start_rec

        # Extracting deformation data
        div_colours = data['div'][min_index:max_index]
        shr_colours = data['shr'][min_index:max_index]
        vrt_colours = data['vrt'][min_index:max_index]

        start_tria = dt()
        tria = tri.Triangulation(start_lon, start_lat, triangles=triangles)
        timer_tria = dt() - start_tria

        start_trans = dt()
        new_lon_lat = ax_div.projection.transform_points(trans, np.array(start_lon), np.array(start_lat))
        timer_trans = dt() - start_trans

        if len(triangles) != 0:
        # Plotting
            start_plt = dt()
            cb_div = ax_div.tripcolor(new_lon_lat[:,0], new_lon_lat[:,1], facecolors=div_colours, triangles = triangles, cmap='coolwarm', vmin=-0.04, vmax=0.04)
            cb_shr = ax_shr.tripcolor(new_lon_lat[:,0], new_lon_lat[:,1], facecolors=shr_colours, triangles = triangles, cmap='plasma', vmin=0, vmax=0.1)
            cb_vrt = ax_vrt.tripcolor(new_lon_lat[:,0], new_lon_lat[:,1], facecolors=vrt_colours, triangles = triangles, cmap='coolwarm', vmin=-0.1, vmax=0.1)
            timer_plot += dt() - start_plt
        # Updating minimum index
        min_index = max_index

    print('Recreate time = ',timer_recreate, 'seconds' )
    print('Triangulation time = ',timer_tria, 'seconds' )
    print('Transform time = ',timer_tria, 'seconds' )
    print('Ploting time = ',timer_plot, 'seconds' )

    # Create a list of colorbars and titles to be iterated over
    cb_list = [cb_div, cb_shr, cb_vrt]
    title_list = ['Divergence Rate $(Days^{-1})$', 'Shear Rate $(Days^{-1})$', 'Vorticity Rate $(Days^{-1})$']

    Date_options = config.config['Date_options']
    start_year   = Date_options['start_year']
    end_year   = Date_options['end_year']
    start_month  = Date_options['start_month']
    end_month  = Date_options['end_month']
    start_day    = Date_options['start_day']
    end_day    = Date_options['end_day']
    timestep     = Date_options['timestep']

    IO            = config.config['IO']
    output_folder = IO['output_folder']
    exp           = IO['exp']

    # Iterate through all axes
    for ax, title, cb in zip(ax_list, title_list, cb_list):
        # Add a colorbar
        plt.colorbar(cb, ax=ax)

        # Add a title
        ax.set_title(title + '\n' + start_year + '-' + start_month + '-' + start_day + ' to ' +
                        end_year + '-' + end_month + '-' + end_day + ', ' + timestep + 'hr')

        # Add gridlines
        ax.gridlines()

        # Hide deformations over land
        ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

        '''
        _________________________________________________________________________________________
        SAVE PLOTS
        '''

    # Set a directory to store figures for the current experiment
    figsPath =  output_folder + '/' + exp + '/figs/'

    # Create the directory if it does not exist already
    os.makedirs(figsPath, exist_ok=True)

    # Create a prefix for the figure filenames
    prefix = data.icetracker + '_' + start_year + start_month + start_day + '_' + end_year + end_month + end_day + '_dt' + str(timestep) + '_tol' + str(data.tolerance)

    # Create the figure filenames
    div_path   = figsPath + prefix + '_div.png'
    shr_path = figsPath + prefix + '_shr.png'
    rot_path   = figsPath + prefix + '_vrt.png'


    for fig, fig_path in zip([fig_div, fig_shr, fig_vrt], [div_path, shr_path, rot_path]):

        # Check if the figures already exist. If they do, delete them.
        if os.path.isfile(fig_path):
            os.remove(fig_path)

        # Save the new figures
        fig.savefig(fig_path, bbox_inches='tight', dpi=600)

if __name__ == '__main__':
    visualise_deformations()
