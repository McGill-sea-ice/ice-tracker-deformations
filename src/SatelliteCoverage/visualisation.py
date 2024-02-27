"""
Author: Mathieu Plante, Lekima Yakuden, Beatrice Duval

--------------------------------------------------------------------------------
Visualisation tool to explore the netcdf dataset
--------------------------------------------------------------------------------

This file contains functions for analysing the SID dataset from the netCDF files.

"""

# Loading from default packages
import os
import sys
parent = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,parent)
from time import strftime
import time
from netCDF4 import Dataset
import numpy as np
import pylab as p
import pyproj as pyproj
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
from tqdm import tqdm
import haversine as hs
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import make_interp_spline
from scipy.spatial import ConvexHull, convex_hull_plot_2d
#mpl.rcParams['text.usetex']=True

# Code from other files
from SatelliteCoverage.config import read_config
from SatelliteCoverage.utils import date_to_seconds, seconds_to_date, convert_to_grid, get_prefix


class visualisation:
    def __init__(self):

        print("Initializing plotting tools")

    def plot_start_end_points(self, data = None, config=None):
        """
        Plots the start and end points of the triangle vertices on a map.
        Start points are blue, end points are red.

        INPUTS:
        data -- data object from the dataset
        config -- Configuration object from the namelist

        OUTPUTS:
        None -- Saves plot to the output directory

        """
        print('--- Plotting start and end points ---')


        icetracker = data.icetracker

        # Reading user options
        Date_options = config['Date_options']
        start_year, start_month, start_day = Date_options['start_year'], Date_options['start_month'], Date_options['start_day']
        end_year, end_month, end_day = Date_options['end_year'], Date_options['end_month'], Date_options['end_day']
        timestep, tolerance = Date_options['timestep'], Date_options['tolerance']

        options = config['options']
        area_filter, centre_lat, centre_lon, radius = options['area_filter'], options['centre_lat'], options['centre_lon'], options['radius']

        start_lats = [data.start_lat1,data.start_lat2,data.start_lat3]
        start_lons = [data.start_lon1, data.start_lon2, data.start_lon3]
        end_lats = [data.end_lat1, data.end_lat2, data.end_lat3]
        end_lons = [data.end_lon1, data.end_lon2, data.end_lon3]

        IO            = config['IO']
        output_folder = IO['output_folder']
        exp           = IO['exp']

        """
        Plot Preamble
        """

        # Setting projection as North Pole Stereographic
        proj = ccrs.NorthPolarStereo(central_longitude=0)

        # Initialize figure
        fig = plt.figure(figsize=(5.5, 5.5), )
        ax = fig.add_subplot(projection = proj, frameon=False)

        # Upper (u) and lower (l) extents of x, y (metres)
        lxextent = -4400000
        uxextent = 2500000
        uyextent = 3500000
        lyextent = -2500000

        ax.set_extent((lxextent, uxextent, uyextent, lyextent), ccrs.NorthPolarStereo())

        # Show lat/lon grid and coastlines
        ax.gridlines(draw_labels=True)
        ax.coastlines()

        # Hide datapoints over land
        ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

        """
        Plotting
        """

        # Iterating plotting thrice to plot all three points of the triangles
        # This could be improved by having the lats lons in separate vectors.
        # Here, points are repeated in many triangles
        for i in range(3):

            # Converting start/end lat/lons to x/y (North pole stereographic)
            start_x, start_y = convert_to_grid(start_lons[i], start_lats[i])
            end_x, end_y = convert_to_grid(end_lons[i], end_lats[i])

            # Plotting start points (Blue)
            ax.scatter(start_x, start_y, color = 'blue', s = 0.1, marker='x')

            # Plotting end points (Red)
            ax.scatter(end_x, end_y, color = 'red', s = 0.1, marker='+')

        # Set title
        ax.set_title(f'Start and End points, {start_year}-{start_month}-{start_day}, {end_year}-{end_month}-{end_day}')

        # Set a directory to store figures
        figsPath =  output_folder + '/' + exp + '/figs/'

        # Create directory if it doesn't exist
        os.makedirs(figsPath, exist_ok=True)

        # Set prefix for filename
        prefix = get_prefix(config=config)

        # Full path of figure
        fig_path = figsPath + prefix + '_start_end_points.png'

        # Check if the path exists and overwriting if it does
        if os.path.isfile(fig_path):
            os.remove(fig_path)

        # Saving figure
        print('Saving start and end points figure at ',fig_path)
        plt.savefig(fig_path, bbox_inches='tight', dpi=600)
        plt.close(fig)

    def recreate_coordinates(self, start_lat1, start_lat2, start_lat3, start_lon1, start_lon2, start_lon3, start_id1, start_id2, start_id3):
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

    def plot_deformations(self, data=None, config=None, datestring=None):
        """
        This function plots deformations from a netCDF file using matplotlib's ax.tripcolor.
        The function assumes that the netCDF was generated from src/SeaIceDeformation's M01_d03_compute_deformations.py.

        INPUTS:
        data -- data object from the dataset
        config -- Configuration object from the namelist

        OUTPUTS:
        None -- Saves plot to the output directory
        """
        try:
            path = config['IO']['netcdf_path']
        except:
            path = None

        # Loading data from netcdf as a dictionary
        if path != None and data == None :
            data = load_data_netcdf(config=config)
        elif path == None and data == None:
            sys.exit('You need to give at least one path or data to plot and not both ')

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
        fig_are = plt.figure()
        fig_defs = plt.figure(figsize=(5, 9))

        # Initialize subplots
#        ax_div = fig_div.add_subplot(111, projection=proj)
#        ax_shr = fig_shr.add_subplot(111, projection=proj)
#        ax_vrt = fig_vrt.add_subplot(111, projection=proj)

        ax_div = fig_defs.add_axes([0.1, 0.62, 0.8, 0.25], projection=proj)
        ax_shr = fig_defs.add_axes([0.1, 0.34, 0.8, 0.25], projection=proj)
        ax_vrt = fig_defs.add_axes([0.1, 0.06, 0.8, 0.25], projection=proj)

        ax_are = fig_are.add_axes([0.1, 0.1, 0.8, 0.8], projection=proj)
        # Create a list of axes to be iterated over
        ax_list = [ax_div, ax_shr, ax_vrt, ax_are]

        for ax in ax_list:
            # Set the map extent in order to see the entire region of interest
            #ax.set_extent((-4400000, 2500000, 3500000, -2500000), ccrs.NorthPolarStereo())
            ax.set_extent((-3800000, 2300000, 3000000, -2500000), ccrs.NorthPolarStereo())
        """
        Plotting
        """

        print('--- Creating sea-ice deformation figures ---')

        # Iterating over all files (Unique triangulations)
        for j in tqdm(np.unique(data.day_flag)):
          no_day = data.no[np.where(data.day_flag==j)]
          for i in tqdm(np.unique(no_day)):
            #if i > 15:
            #   continue

            # Obtaining number of rows corresponding to triangles in given file that will be iterated over
            file_length = np.count_nonzero(data.no == i)

            # Setting maximum index
            print(len(data.no),len(data.day_flag))
            condi = (data.no[:] == i) & (data.day_flag[:] == j)
            min_index = np.where(condi)[0][0]
            max_index = np.where(condi)[0][-1]+1



            # Arranging triangle vertices in array for use in ax.tripcolor
            triangles = np.stack((data.idx1[min_index:max_index],
                                  data.idx2[min_index:max_index],
                                  data.idx3[min_index:max_index]), axis=-1)
            print(triangles.shape)
            Mask = triangles.copy()
            Mask[:,:] = 0
            Mask[data.Mask[min_index:max_index]== 0,:] = 1
            print(Mask.shape)
            # Filtering data range to that of the current "file"
            (start_lat1_temp,
             start_lat2_temp,
             start_lat3_temp) = (data.start_lat1[min_index:max_index],
                                 data.start_lat2[min_index:max_index],
                                 data.start_lat3[min_index:max_index])

            (start_lon1_temp,
             start_lon2_temp,
             start_lon3_temp) = (data.start_lon1[min_index:max_index],
                                 data.start_lon2[min_index:max_index],
                                 data.start_lon3[min_index:max_index])

            start_lat, start_lon = self.recreate_coordinates(start_lat1_temp, start_lat2_temp, start_lat3_temp,
                                                        start_lon1_temp, start_lon2_temp, start_lon3_temp,
                                                        data.idx1[min_index:max_index],
                                                        data.idx2[min_index:max_index],
                                                        data.idx3[min_index:max_index])

            # Extracting deformation data
            div_colours = data.div[min_index:max_index]
            shr_colours = data.shr[min_index:max_index]
            vrt_colours = data.vrt[min_index:max_index]
            are_colours = data.A[min_index:max_index]

            # tranform the coordinates already to improve the plot efficiency
            new_coords = proj.transform_points(trans, np.array(start_lon), np.array(start_lat))
            # create one triangulation object
            print(triangles.shape,Mask.shape)
            print(len(triangles),len(Mask))
            tria = tri.Triangulation(new_coords[:,0], new_coords[:,1], triangles=triangles) #s, mask = Mask[:,0])

            if len(triangles) != 0:
            # Plotting
                cb_div = ax_div.tripcolor(tria, facecolors=div_colours, cmap='coolwarm', vmin=-0.04, vmax=0.04)
                cb_shr = ax_shr.tripcolor(tria, facecolors=shr_colours, cmap='plasma', vmin=0, vmax=0.1)
                cb_vrt = ax_vrt.tripcolor(tria, facecolors=vrt_colours, cmap='coolwarm', vmin=-0.1, vmax=0.1)
                cb_are = ax_are.tripcolor(tria, facecolors=are_colours, cmap='plasma',vmin=0,vmax=20000.0**2.0)

        # Create a list of colorbars and titles to be iterated over
        test_str = "Triangulated data for $\dot{\epsilon}^{-1}$ pair 'no=%s'"
        cb_list = [cb_div, cb_shr, cb_vrt]
        eI_title ="$\dot{\epsilon}_{I}$ (day$^{-1}$)"
        eII_title = "$\dot{\epsilon}_{II}$ (day$^{-1}$)"
        zeta_title = "$\zeta$ (day$^{-1}$)"
        are_title = "A (m)"
        title_list = [eI_title, eII_title,zeta_title]
        ax_list = [ax_div, ax_shr, ax_vrt]

        Date_options  = config['Date_options']
        start_year    = Date_options['start_year']
        end_year      = Date_options['end_year']
        start_month   = Date_options['start_month']
        end_month     = Date_options['end_month']
        start_day     = Date_options['start_day']
        end_day       = Date_options['end_day']
        timestep      = Date_options['timestep']
        tolerance     = Date_options['tolerance']

        IO            = config['IO']
        output_folder = IO['output_folder']
        exp           = IO['exp']

        lbl_list = ["a)","b)","c)"]
        # Iterate through all axes
        for ax, title, cb, lbl in zip(ax_list, title_list, cb_list, lbl_list):
            # Add a colorbar
            clb = plt.colorbar(cb, ax = ax, shrink = 0.9,pad=.04)
            clb.ax.tick_params(labelsize=8)
            plt.text(1.12,1.02,title,ha='center', va='center', transform=ax.transAxes,fontsize=8)
            plt.text(-0.1,0.95,lbl,ha='center', va='center', transform=ax.transAxes,fontsize=12)
#            ax.text(0.5, -0.3,
#                title,
#                va='bottom', ha='center',
#                rotation='horizontal',
#                rotation_mode='anchor',
#                transform=ax.transAxes,fontsize=14)
            # Add gridlines
            ax.gridlines()

            # Hide deformations over land
            ax.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

            '''
            _________________________________________________________________________________________
            SAVE PLOTS
            '''
        if datestring is None:
            prefix = get_prefix(config = config)
        else:
            prefix = datestring
        # Set a directory to store figures for the current experiment
        figsPath =  output_folder + '/' + exp + '/figs/Deformations/'

        # Create the directory if it does not exist already
        os.makedirs(figsPath, exist_ok=True)

        itr = config['Metadata']['icetracker']
        # Create a prefix for the figure filenames
        prefix = get_prefix(config = config)
        if datestring is not None:
            prefix = datestring
        # Create the figure filenames
        div_path   = figsPath + prefix + '_div.png'
        shr_path   = figsPath + prefix + '_shr.png'
        rot_path   = figsPath + prefix + '_vrt.png'
        defs_path  = figsPath + prefix + '_defs.png'
        are_path   = figsPath + prefix + '_area.png'

#        for fig, fig_path in zip([fig_div, fig_shr, fig_vrt, fig_are], [div_path, shr_path, rot_path, are_path]):
        for fig, fig_path in zip([fig_defs, fig_are],[defs_path,are_path]):
            # Check if the figures already exist. If they do, delete them.
            if os.path.isfile(fig_path):
                os.remove(fig_path)

            # Save the new figures
            print('Saving deformation figure at ',fig_path)
            fig.savefig(fig_path, bbox_inches='tight', dpi=600)
            plt.close(fig)
        print(np.max(data.A[:]**0.5))

        return


    def plot_pairs(self, data=None, config=None):
        """
        This function plots deformations from a netCDF file using matplotlib's ax.tripcolor.
        The function assumes that the netCDF was generated from src/SeaIceDeformation's M01_d03_compute_deformations.py.

        INPUTS:
        data -- data object from the dataset
        config -- Configuration object from the namelist

        OUTPUTS:
        None -- Saves plot to the output directory
        """
        print("Drawing the contours of each pairs on a map")
        try:
            path = config['IO']['netcdf_path']
        except:
            path = None

        # Loading data from netcdf as a dictionary
        if path != None and data == None :
            data = load_data_netcdf(config=config)
        elif path == None and data == None:
            sys.exit('You need to give at least one path or data to plot and not both ')


        """
        Preamble
        """
        # Set the matplotlib projection and transform
        proj = ccrs.NorthPolarStereo(central_longitude=0)
        trans = ccrs.Geodetic()

        # Initialize figure
        fig_pairs = plt.figure()

        # Initialize subplots
        ax_pairs = fig_pairs.add_subplot(111, projection=proj)
        ax_pairs.set_extent((-4400000, 2500000, 3500000, -2500000), ccrs.NorthPolarStereo())

        """
        Plotting
        """

        print('--- Creating sea-ice pairs coverage figures ---')

        # Iterating over all files (Unique triangulations)
        for j in tqdm(np.unique(data.day_flag)):
          no_day = data.no[np.where(data.day_flag==j)]
          for i in tqdm(np.unique(no_day)):
            #if i > 15:
            #   continue

            # Obtaining number of rows corresponding to triangles in given file that will be iterated over
            file_length = np.count_nonzero(data.no == i)

            # Setting maximum index
            print(len(data.no),len(data.day_flag))
            condi = (data.no[:] == i) & (data.day_flag[:] == j)
            min_index = np.where(condi)[0][0]
            max_index = np.where(condi)[0][-1]+1

            # Arranging triangle vertices in array for use in ax.tripcolor
            triangles = np.stack((data.idx1[min_index:max_index],
                                  data.idx2[min_index:max_index],
                                  data.idx3[min_index:max_index]), axis=-1)

            # Filtering data range to that of the current "file"
            (start_lat1_temp,
             start_lat2_temp,
             start_lat3_temp) = (data.start_lat1[min_index:max_index],
                                 data.start_lat2[min_index:max_index],
                                 data.start_lat3[min_index:max_index])

            (start_lon1_temp,
             start_lon2_temp,
             start_lon3_temp) = (data.start_lon1[min_index:max_index],
                                 data.start_lon2[min_index:max_index],
                                 data.start_lon3[min_index:max_index])

            start_lons = np.concatenate((start_lon1_temp, start_lon2_temp, start_lon3_temp), axis=0, out=None, dtype=None, casting="same_kind")
            start_lats = np.concatenate((start_lat1_temp, start_lat2_temp, start_lat3_temp), axis=0, out=None, dtype=None, casting="same_kind")
            start_x, start_y = convert_to_grid(start_lons[:], start_lats[:])

            points = np.ndarray((len(start_x),2))
            points[:,0],points[:,1] = start_x,start_y

            # Make the hull around the points
            hull = ConvexHull(points)
#            hull_pts = self.convex_hull(points)
#            if hull_pts is not None:
#                #hull_pts = scale*hull_pts
#                print(hull_pts)
#                a,b = hull_pts.shape
#                x= np.zeros(a)
#                y= np.zeros(a)
#                for j in range(0,a):
#                    x[j] = hull_pts[j][0]
#                    y[j] = hull_pts[j][1]
#                print(x)
#                print(y)
#            nt = np.linspace(0, 1, 100)
#            t = np.zeros(x.shape)
#            t[1:] = np.sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2)
#            t = np.cumsum(t)
#            t /= t[-1]
#            spl = make_interp_spline(x, y)
#            x2, y2 = spl(nt).T
#make_interp_spline(t, y)
#                plt.plot(x, y,'r--',linewidth=2)
            for simplex in hull.simplices:
                plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

            x_mean = np.nanmean(np.squeeze(points[hull.simplices,0]))
            y_mean = np.nanmean(np.squeeze(points[hull.simplices,1]))
            ax_pairs.scatter(start_x, start_y, color = 'red', s = 0.1, marker='+')
            ax_pairs.text(x_mean, y_mean,str(i),fontsize = 3,
                          verticalalignment ='center',
                                  horizontalalignment = 'center')
        title = ['Matching pairs area']

        Date_options  = config['Date_options']
        start_year    = Date_options['start_year']
        end_year      = Date_options['end_year']
        start_month   = Date_options['start_month']
        end_month     = Date_options['end_month']
        start_day     = Date_options['start_day']
        end_day       = Date_options['end_day']
        timestep      = Date_options['timestep']
        tolerance     = Date_options['tolerance']

        IO            = config['IO']
        output_folder = IO['output_folder']
        exp           = IO['exp']


        # Add gridlines
        ax_pairs.gridlines()

        # Hide deformations over land
        ax_pairs.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

        '''
        _________________________________________________________________________________________
        SAVE PLOTS
        '''

        # Set a directory to store figures for the current experiment
        figsPath =  output_folder + '/' + exp + '/figs/'

        # Create the directory if it does not exist already
        os.makedirs(figsPath, exist_ok=True)

        itr = config['Metadata']['icetracker']
        # Create a prefix for the figure filenames
        prefix = get_prefix(config = config)

        # Create the figure filenames
        pair_path   = figsPath + prefix + '_pairs.png'

        # Check if the figures already exist. If they do, delete them.
        if os.path.isfile(pair_path):
            os.remove(pair_path)

        # Save the new figures
        print('Saving pairs area figure at ',pair_path)
        fig_pairs.savefig(pair_path, bbox_inches='tight', dpi=600)
        plt.close(fig_pairs)
        return



    def plot_triangles(self, data=None, config=None, no = None, show_ID = None):
        """
        This function plots deformations from a netCDF file using matplotlib's ax.tripcolor.
        The function assumes that the netCDF was generated from src/SeaIceDeformation's M01_d03_compute_deformations.py.

        INPUTS:
        data -- data object from the dataset
        config -- Configuration object from the namelist

        OUTPUTS:
        None -- Saves plot to the output directory
        """
        print("Drawing the contours of each pairs on a map")
        try:
            path = config['IO']['netcdf_path']
        except:
            path = None

        if show_ID is None:
            show_ID = False


        title = ['Triangulated data']
        Date_options  = config['Date_options']
        start_year    = Date_options['start_year']
        end_year      = Date_options['end_year']
        start_month   = Date_options['start_month']
        end_month     = Date_options['end_month']
        start_day     = Date_options['start_day']
        end_day       = Date_options['end_day']
        timestep      = Date_options['timestep']
        tolerance     = Date_options['tolerance']

        IO            = config['IO']
        output_folder = IO['output_folder']
        exp           = IO['exp']

        # Set a directory to store figures for the current experiment
        figsPath =  output_folder + '/' + exp + '/figs/'

        # Create the directory if it does not exist already
        os.makedirs(figsPath, exist_ok=True)

        itr = config['Metadata']['icetracker']
        # Create a prefix for the figure filenames
        prefix = get_prefix(config = config)





        # Loading data from netcdf as a dictionary
        if path != None and data == None :
            data = load_data_netcdf(config=config)
        elif path == None and data == None:
            sys.exit('You need to give at least one path or data to plot and not both ')


        """
        Preamble
        """
        # Set the matplotlib projection and transform
        proj = ccrs.NorthPolarStereo(central_longitude=0)
#        trans = ccrs.epsg(3995)
        trans = ccrs.Geodetic()
#        trans = ccrs.NorthPolarStereo(central_longitude=0)
#        transformer = pyproj.Transformer.from_crs("EPSG:4326",'+proj=stere +lat_0=90 +lat_ts=50 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ')
        transformer = pyproj.Transformer.from_crs("EPSG:4326",'EPSG:3995')
        trans_inv =  pyproj.Transformer.from_crs('EPSG:3995',"EPSG:4326")
        # Initialize figure
        fig_tris = plt.figure()
        ax_tris = fig_tris.add_subplot(111, projection=proj)
        ax_tris.set_extent([-4400000, 2500000, 3500000, -2500000], proj)

        """
        Plotting
        """

        print('--- Creating sea-ice pairs coverage figures ---')

        # Iterating over all files (Unique triangulations)

        if no is not None:
            pair_listing = np.array(no)
            list_string = str(no[0])
            print(pair_listing, list_string)
        else:
            pair_listing = tqdm(np.unique(data.no))
            list_string = 'all'

        for i in pair_listing:
            j = np.unique(data.day_flag[data.no == i])[0]
            condi = (data.no[:] == i) & (data.day_flag[:] == j)
            print("j is : ", j)
            # Obtaining number of rows corresponding to triangles in given file that will be iterated over
            file_length = np.count_nonzero(condi)

            # Setting maximum index
            min_index = np.where(condi)[0][0]
            max_index = np.where(condi)[0][-1]+1
            print("size of SAR scene :", file_length, min_index, max_index)

            # Arranging triangle vertices in array for use in ax.tripcolor
            triangles = np.stack((data.idx1[min_index:max_index],
                                  data.idx2[min_index:max_index],
                                  data.idx3[min_index:max_index]), axis=-1)

            # Filtering data range to that of the current "file"
            (start_lat1_temp,
             start_lat2_temp,
             start_lat3_temp) = (data.start_lat1[min_index:max_index],
                                 data.start_lat2[min_index:max_index],
                                 data.start_lat3[min_index:max_index])

            (start_lon1_temp,
             start_lon2_temp,
             start_lon3_temp) = (data.start_lon1[min_index:max_index],
                                 data.start_lon2[min_index:max_index],
                                 data.start_lon3[min_index:max_index])

            start_lons = np.concatenate((start_lon1_temp, start_lon2_temp, start_lon3_temp), axis=0, out=None, dtype=None, casting="same_kind")
            start_lats = np.concatenate((start_lat1_temp, start_lat2_temp, start_lat3_temp), axis=0, out=None, dtype=None, casting="same_kind")

            start_lat, start_lon = self.recreate_coordinates(start_lat1_temp, start_lat2_temp, start_lat3_temp,
                                                        start_lon1_temp, start_lon2_temp, start_lon3_temp,
                                                        data.idx1[min_index:max_index],
                                                        data.idx2[min_index:max_index],
                                                        data.idx3[min_index:max_index])

            # Extracting deformation data
            div_colours = data.div[min_index:max_index]


            # tranform the coordinates already to improve the plot efficiency
            new_coords = proj.transform_points(trans, np.array(start_lon), np.array(start_lat))

            if show_ID is True:

                # Repeating for the end points
                (end_lat1_temp,
                 end_lat2_temp,
                 end_lat3_temp) = (data.end_lat1[min_index:max_index],
                                 data.end_lat2[min_index:max_index],
                                 data.end_lat3[min_index:max_index])


                (end_lon1_temp,
                 end_lon2_temp,
                 end_lon3_temp) = (data.end_lon1[min_index:max_index],
                                 data.end_lon2[min_index:max_index],
                                 data.end_lon3[min_index:max_index])

                end_lons = np.concatenate((end_lon1_temp, end_lon2_temp, end_lon3_temp), axis=0, out=None, dtype=None, casting="same_kind")
                end_lats = np.concatenate((end_lat1_temp, end_lat2_temp, end_lat3_temp), axis=0, out=None, dtype=None, casting="same_kind")
                #start_x, start_y = convert_to_grid(start_lons[:], start_lats[:])
                end_lat, end_lon = self.recreate_coordinates(end_lat1_temp, end_lat2_temp, end_lat3_temp,
                                                        end_lon1_temp, end_lon2_temp, end_lon3_temp,
                                                        data.idx1[min_index:max_index],
                                                        data.idx2[min_index:max_index],
                                                        data.idx3[min_index:max_index])

                new_coords_end = proj.transform_points(trans, np.array(end_lon), np.array(end_lat))
                new_idx1 = proj.transform_points(trans, data.start_lon1[min_index:max_index], data.start_lat1[min_index:max_index])
                new_idx2 = proj.transform_points(trans, data.start_lon2[min_index:max_index], data.start_lat2[min_index:max_index])
                new_idx3 = proj.transform_points(trans, data.start_lon3[min_index:max_index], data.start_lat3[min_index:max_index])

                new_idx1_end = proj.transform_points(trans, data.end_lon1[min_index:max_index], data.end_lat1[min_index:max_index])
                new_idx2_end = proj.transform_points(trans, data.end_lon2[min_index:max_index], data.end_lat2[min_index:max_index])
                new_idx3_end = proj.transform_points(trans, data.end_lon3[min_index:max_index], data.end_lat3[min_index:max_index])

                IDs = np.concatenate((data.idx1[min_index:max_index], data.idx2[min_index:max_index], data.idx3[min_index:max_index]), axis=0)
                IDxs = np.concatenate((new_idx1[:,0], new_idx2[:,0], new_idx3[:,0]), axis=0)
                IDys = np.concatenate((new_idx1[:,1], new_idx2[:,1], new_idx3[:,1]), axis=0)
                IDxs_end = np.concatenate((new_idx1_end[:,0], new_idx2_end[:,0], new_idx3_end[:,0]), axis=0)
                IDys_end = np.concatenate((new_idx1_end[:,1], new_idx2_end[:,1], new_idx3_end[:,1]), axis=0)
                tria_end = tri.Triangulation(new_coords_end[:,0], new_coords_end[:,1], triangles=triangles)
                print(IDxs.shape)
                deltaX = IDxs_end[:] - IDxs[:]
                deltaY = IDys_end[:] - IDys[:]

            # create one triangulation object
            tria = tri.Triangulation(new_coords[:,0], new_coords[:,1], triangles=triangles)

            # tranform the coordinates already to improve the plot efficiency
            new_coords = proj.transform_points(trans, np.array(start_lon), np.array(start_lat))
            new_coords_ref = proj.transform_points(trans,np.array([0.0,-180.0, -90.0,90.0]),np.array([62.0,62.0,62.0,62.0]))
            x_ref, y_ref  = new_coords_ref[:,0], new_coords_ref[:,1]

            for i  in range(len(x_ref)):
                print(x_ref[i],y_ref[i])
            x,y = new_coords[:,0],new_coords[:,1]
            a1 = np.nanmin(x)
            a2 = np.nanmax(x)
            a3 = np.nanmin(y)
            a4 = np.nanmax(y)
            a = len(x)

            points = np.ndarray((a,2))
            points[:,0],points[:,1] = x,y
            points = points[~np.isnan(points)]
            try:
                ax_tris.set_extent((a1, a2, a3, a4),proj)
            except:
                sys.exit("error in extent : %s, %s, %s, %s" % (a1, a2, a3, a4))

            if len(triangles) != 0:
                # Plotting
                cb_div = ax_tris.tripcolor(tria, facecolors=div_colours,cmap='coolwarm', vmin=-0.04, vmax=0.04, edgecolors='k')
                # Make the hull around the points
                try:
                    hull = ConvexHull(points)
                    for simplex in hull.simplices:
                        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
                except:
                    show_ID = False
                    print("Not sure what is wrong with the hull... Continuing without it")
                #print the pair number
                #ax_tris.scatter(start_x, start_y, color = 'red', s = 0.1, marker='+')

                #ax_tris.set_extent((a1-50.0, a2 - 5.0, a1+50.0, a2+5.0))
                plt.title("Triangulated data for pair 'no=%s'" % no )


                # Add gridlines
                ax_tris.gridlines()

                # Hide deformations over land
                ax_tris.add_feature(cfeature.LAND, zorder=100, edgecolor='k')

                if show_ID is True:

                    # Initialize figure
                    fig_zoom = plt.figure()
                    ax_zoom = fig_zoom.add_subplot(111, projection=proj)
                    x_mean = np.nanmean(np.squeeze(points[hull.simplices,0]))
                    y_mean = np.nanmean(np.squeeze(points[hull.simplices,1]))

                    #Show triangles
                    cb_div_zoom = ax_zoom.tripcolor(tria, facecolors=div_colours*np.nan,cmap='coolwarm', vmin=-0.04, vmax=0.04, edgecolors='k')

                    #Scatter start and end positions
                    ax_zoom.scatter(new_coords[:,0], new_coords[:,1], color = 'k', s = 5, marker='*',label = 'start points')
                    ax_zoom.scatter(new_coords_end[:,0], new_coords_end[:,1], color = 'b', s = 5, marker='*',label = 'end points')

#                    cb_div_zoom_end = ax_zoom.tripcolor(tria_end, facecolors=div_colours*np.nan,cmap='coolwarm', vmin=-0.04, vmax=0.04, edgecolors='b')
                    a= len(IDs)
                    b= len(IDxs)
                    #print(a,b)
                    ks = np.arange(a)
                    for idk in tqdm(np.unique(IDs)):
                        # Setting maximum index
                        #print(idk)
                        indk = ks[IDs[:] == idk]
                        if len(indk)> 1: indk = indk[0]
                        #print(indk)
                        #print(IDxs[indk], IDys[indk], deltaX[indk], deltaY[indk])
                        t = ax_zoom.text(IDxs[indk], IDys[indk],str(IDs[indk]),fontsize = 4,clip_on=True)
                        t.clipbox = ax_zoom.bbox
                        ax_zoom.quiver(IDxs[indk], IDys[indk], deltaX[indk], deltaY[indk], angles='xy', scale_units='xy', scale=1, width = 0.005,color = 'b')
#                        ax_zoom.set_extent((a1, a2, a3, a4),proj)
                    ax_zoom.set_extent([x_mean - (a2-a1)/20, x_mean + (a2-a1)/20, y_mean - (a4-a3)/30, y_mean + (a4-a3)/30], proj)
                    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                                        ncol=2, mode="expand", borderaxespad=0.,fontsize=10.0 )


        '''
        _________________________________________________________________________________________
        SAVE PLOTS
        '''

        # Create the figure filenames
        tri_path   = figsPath + prefix + '_triangulated_no_%s.png' % list_string

        # Check if the figures already exist. If they do, delete them.
        if os.path.isfile(tri_path):
            os.remove(tri_path)

        # Save the new figures
        print('Saving triangulated data figure at ',tri_path)
        fig_tris.savefig(tri_path, dpi=600)
        plt.close(fig_tris)

        if show_ID is True:
            zoom_path   = figsPath + prefix + '_triangulated_no_%s_zoomed.png' % list_string
            if os.path.isfile(zoom_path):
                os.remove(zoom_path)
            print('Saving zoomed data figure at ',zoom_path)
            fig_zoom.savefig(zoom_path, dpi=600)
            plt.close(fig_zoom)
        return


    def plot_area_dist(config, data = None):

        A = np.sqrt(data.A[:]/1e3)
#######
#What is this in the dataset?
        satellite = np.array(d['satellite'][:])
#######

        dens = True
        log = False
        bin_step = 1
        ys = 2017
        ye = 2022
        alp = 0.5
        hsty = 'stepfilled' # 'bar', 'barstacked', 'step', 'stepfilled'
        bins = range(0,25+bin_step,bin_step) # bins initialization

        # Dictionnaries initialization
        y={}
        y['all'] =[]

        # figure initialization
        plt.figure()

        sat_list = ['rcm','s1']
        sat_nb = [0,1]
        # Read the files names and create the data lists
        for s,sat in enumerate(sat_nb):
            y[sat_list[s]] = A[np.where(satellite == sat)[0].tolist()].tolist()
            y['all'] += y[sat_list[s]]

        y_list = [ y[name] for name in sat_list ]
        label_list = sat_list
        if hsty != 'barstacked':
            y_list += [y['all']]
            label_list += ['all']
        hist = plt.hist(y_list, bins=bins, density=dens, label=label_list, histtype=hsty, alpha = alp, log=log)

        # the figure bells and whistles
        plt.legend()
        if dens:
            plt.ylabel('PDF')
        else :
            plt.ylabel('Number of triangles')
        plt.xlabel('Nominal size (km)')
        plt.xticks(range(0,25,5))
        plt.title(r'Triangle size ($\sqrt{A}$) in km for S1 and RCM between ' + str(ys) + ' and '+ str(ye))

        # saving the figure
        # Create the directory if it does not exist already
        output_folder = '/storage/amelie/RCMS1_analysis/outputs/'+'add_area'+'/figs/'
        os.makedirs(output_folder, exist_ok=True)
#        plt.savefig(output_folder + 'RCMS1_{}_{}_bin{}_A_hist.png'.format(ys,ye,bin_step), bbox_inches='tight', dpi=600)

       # Saving the data for after
#        data={}
#        data['bins_start'] = hist[1][:-1]
#        data['bins_end'] = hist[1][1:]
#        for i in range(len(sat_list)):
#            name = label_list[i]
#            data[name] = hist[0][i]
#        df = pd.DataFrame(data, columns=['bins_start','bins_end'] + label_list)
#        print('Saving size distribution data at ' + output_folder + 'RCMS1_{}_{}_bin{}_dt_hist.pkl'.format(ys, ye, bin_step$
#        df.to_pickle(output_folder + 'RCMS1_{}_{}_bin{}_A_hist.pkl'.format(ys, ye, bin_step))

        return None








