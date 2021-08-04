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
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

import config
import utils_load_data as load_data


def visualise_deformations():

    '''
    _________________________________________________________________________________________
    INITIALIZE PLOTS
    '''

    # Set the matplotlib projection and transform
    proj = ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90)
    trans = ccrs.Geodetic()

    # Create a colormap for the divergence ax in order that 0 is 
    # associated to the same color as the other plots 
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","#110788","yellow"])
    #cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","#110788","#3e0b15"])
    
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
    for raw_path, triangulated_path, calculated_path, converted_path in zip(config.data_paths['raw'], config.data_paths['triangulated'], config.data_paths['calculations'], config.data_paths['converted']):

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

            # Load converted data
            converted_data = load_data.load_converted(converted_path)
            sX1         = converted_data['sX1']             # Starting x coordinates in a local grid CS
            sX2         = converted_data['sX2']
            sX3         = converted_data['sX3']

            sY1         = converted_data['sY1']             # Starting y coordinates in a local grid CS
            sY2         = converted_data['sY2']
            sY3         = converted_data['sY3']

            eX1         = converted_data['eX1']             # Starting x coordinates in a local grid CS
            eX2         = converted_data['eX2']
            eX3         = converted_data['eX3']

            eY1         = converted_data['eY1']             # Starting y coordinates in a local grid CS
            eY2         = converted_data['eY2']
            eY3         = converted_data['eY3']
            
            # Load calculated data
            calculated_data = load_data.load_calculations(calculated_path)
            eps_I           = calculated_data['eps_I']    # Divergence rate
            eps_II          = calculated_data['eps_II']   # Maximum shear strain rate
            eps_tot         = calculated_data['eps_tot']  # Total deformation rate
            dvdx            = calculated_data['dvdx']       
            dudy            = calculated_data['dudy']
            dvdy            = calculated_data['dvdy']       
            dudx            = calculated_data['dudx']
            rot             = dvdx - dudy


        except:
            continue

        '''
        _________________________________________________________________________________________
        PLOT DEFORMATIONS
        '''

        # Stack the 3 vertices index arrays
        triangles = np.stack((vertice_idx1, vertice_idx2, vertice_idx3), axis=-1)

        # Plot total deformations, divergence and shear
        cb_tot = ax_tot.tripcolor( sLon, sLat, triangles, facecolors=eps_tot, transform=trans, cmap='plasma', vmin=0, vmax=0.1 )
        cb_I = ax_I.tripcolor( sLon, sLat, triangles, facecolors=eps_I, transform=trans, cmap='coolwarm', vmin=-0.04, vmax=0.04 )
        cb_II = ax_II.tripcolor( sLon, sLat, triangles, facecolors=eps_II, transform=trans, cmap='plasma', vmin=0, vmax=0.1 )
        cb_rot = ax_rot.tripcolor( sLon, sLat, triangles, facecolors=rot, transform=trans, cmap='coolwarm', vmin=-0.1, vmax=0.1 )

    '''
    _________________________________________________________________________________________
    ADD FEATURES TO THE PLOTS
    '''

    # Check if at least one dataset has been plotted
    try: 
        cb_tot
    except NameError: 
        cb_tot = None
    
    # Retrieve the namelist arguments that define the dataset 
    # that has been selected for processing
    data_folder  = os.path.basename(config.config['IO']['raw_data_folder'])
    Date_options = config.config['Date_options']
    start_year   = Date_options['start_year']
    start_month  = Date_options['start_month']
    start_day    = Date_options['start_day']
    duration     = Date_options['duration']

    if cb_tot is not None:

        # Create a list of colorbars and titles to be iterated over
        cb_list = [cb_tot, cb_I, cb_II, cb_rot]
        title_list = ['Total Deformation Rate $(Days^{-1})$', 'Divergence Rate $(Days^{-1})$', 'Shear Rate $(Days^{-1})$', 'Rotation Rate $(Days^{-1})$']

        # Iterate through all axes
        for ax, title, cb in zip(ax_list, title_list, cb_list):
            # Add a colorbar
            plt.colorbar(cb, ax=ax)

            # Add a title
            ax.set_title(title + '\n' + start_year + '-' + start_month + '-' + start_day + ' -- over ' + str(duration) + ' day')

            # Add gridlines
            ax.gridlines()

            # Hide deformations over land
            ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

        

        # Set a directory to store figures
        currPath = os.path.dirname(os.path.realpath(__file__))
        figsPath = currPath + '/../figs/' + data_folder

        # Create the directory if it does not exist already
        os.makedirs(figsPath, exist_ok=True)

        # Create a prefix for the figure filenames
        prefix = start_year + start_month + start_day + '_' + str(duration) + '.png'
        
        # Create the figure filenames
        tot_path = figsPath + '/tot_' + prefix
        div_path = figsPath + '/div_'+ prefix
        shear_path = figsPath + '/shear_' + prefix
        rot_path = figsPath + '/rot_' + prefix

        # Check if the files already exist. If they do, delete them.
        for fig_path in [tot_path, div_path, shear_path, rot_path]:
            if os.path.isfile(fig_path):
                os.remove(fig_path)

        # Save all figures
        fig_tot.savefig(tot_path , bbox_inches='tight')
        fig_I.savefig( div_path, bbox_inches='tight')
        fig_II.savefig( shear_path, bbox_inches='tight')
        fig_rot.savefig( rot_path, bbox_inches='tight')

        plt.show()

    else:
        print('')

if __name__ == '__main__':
    visualise_deformations()