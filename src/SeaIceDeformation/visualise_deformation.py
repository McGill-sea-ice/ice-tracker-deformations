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

def visualise_deformations():

    '''
    _________________________________________________________________________________________
    INITIALIZE PLOTS
    '''

    # Set the matplotlib projection and transform
    proj = ccrs.NorthPolarStereo(central_longitude=0)
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
        ax.set_extent((-4400000, 2500000, 3500000, -2500000), ccrs.NorthPolarStereo())

    dp = config['data_paths']
    # Iterate through data files in config
    for raw_path, triangulated_path, calculated_path in zip(tqdm(dp['raw']), dp['triangulated'], dp['calculations']):

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

        # tranform the coordinates already to improve the plot efficiency
        new_coords = proj.transform_points(trans, np.array(s_lon), np.array(s_lat))

        # create one triangulation object
        tria = tri.Triangulation(new_coords[:,0], new_coords[:,0], triangles=triangles)

        # Plot total deformations, divergence and shear
        cb_tot = ax_tot.tripcolor( tria, facecolors=eps_tot, cmap='plasma', vmin=0, vmax=0.1 )
        cb_I   = ax_I.tripcolor( tria, facecolors=eps_I, cmap='coolwarm', vmin=-0.04, vmax=0.04 )
        cb_II  = ax_II.tripcolor( tria, facecolors=eps_II, cmap='plasma', vmin=0, vmax=0.1 )
        cb_rot = ax_rot.tripcolor( tria, facecolors=rot, cmap='coolwarm', vmin=-0.1, vmax=0.1 )

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

if __name__ == '__main__':
    visualise_deformations()
