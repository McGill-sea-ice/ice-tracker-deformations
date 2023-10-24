'''
Author: Beatrice Duval (bdu002)

---------------------------------------
M01 - Sea-ice deformations calculations
---------------------------------------

Module that computes triangle cells' sea-ice deformations using a triangulated csv file.

The results are then stored in a calculations csv file for each dataset, and in a
single output netcdf4 file that combines all processed datasets.

1. Calculations csv file:

    Description: Each dataset to process results in a calculations .csv file (if possible).

    Format:

        no. | dudx | dudy | dvdx | dvdy | eps_I | eps_II | vrt | eps_tot

    Variables:

    - no.: the triangle number
    - dudx, dudy, dvdx, dvdy: strain rates (days^-1)
    - eps_I: divergence rate (days^-1)
    - eps_II: maximum shear strain rate (days^-1)
    - vrt: vorticity rate (days^-1)
    - eps_tot: total sea-ice deformation rate (days^-1)


2. Output netcdf4 file:

    Description: The deformation results of all processed datasets are stored
                 in a single output netcdf file. Each data row corresponds to
                 a triangle.


    Metadata:

    - iceTracker: RCM or S1 ice tracker
    - referenceTime: Reference time (time units are in seconds since the reference time)
    - trackingError: Ice tracker error (or resolution)

    Variables:

    - start_time, end_time: Triangle starting and ending times (in seconds since the reference time)
    - start_lat1, start_lat2, start_lat3: Starting latitudes of each triangle's vertices (degrees North)
    - start_lon1, start_lon2, start_lon3: Starting latitudes of each triangle's vertices (degrees East)
    - end_lat1, end_lat2, end_lat3: Ending latitudes of each triangle's vertices (degrees North)
    - end_lon1, end_lon2, end_lon3: Ending latitudes of each triangle's vertices (degrees East)
    - div: Divergence rate of each triangle (days^-1)
    - shear: Shear strain rate of each triangle (days^-1)
    - vrt: Vorticity rate of each triangle (days^-1)

'''

# Loading from default packages
import csv
import datetime
import os
from math import sqrt
import numpy as np
from netCDF4 import Dataset
from tqdm import tqdm
#from datetime import datetime, date, time, timedelta
# Loading from other files
import SeaIceDeformation.utils_datetime as utils_datetime
import SeaIceDeformation.utils_load_data as load_data


def compute_deformations(config=None):
    '''
    _________________________________________________________________________________________
    INITIALIZE OUTPUT NETCDF VARIABLES
    '''

    # Initialize arrays of data to be stored in a netcdf file.
    # These arrays will combine the output data from all datasets
    # that are being processed.
    sTime, eTime, \
        sLat1, sLat2, sLat3, sLon1, sLon2, sLon3, \
        eLat1, eLat2, eLat3, eLon1, eLon2, eLon3, \
        div, shr, vrt, idx1, idx2, idx3, no, idx_sLat1, idx_sLat2, idx_sLat3, A_l, dudx_l, dudy_l, dvdx_l, dvdy_l, sat_l \
        = ([] for i in range(30))

    # Retrieve the starting date common to all processed datasets (from namelist.ini)
    Date_options = config['Date_options']
    YYYY = Date_options['start_year']
    MM = Date_options['start_month']
    DD = Date_options['start_day']

    # Create a datetime object for the reference start time
    refTime = datetime.datetime(int(YYYY), # Year
                       int(MM),   # Month
                       int(DD),   # Day
                       0, 0, 0)   # Hour, minute, second

    # Counter of current file being processed
    file_num = 0

    # Iterate through all raw, triangulated and calculations data files listed in config
    empty_files_list = []
    nbfb = 0
    nbfg = 0
    nbpg = 0
    print('--- Computing sea-ice deformations ---')
    dp = config['data_paths']
    for raw_path, triangulated_path, calculations_path in zip(tqdm(dp['raw']), dp['triangulated'], dp['calculations']):
        '''
        _________________________________________________________________________________________
        LOAD DATA
        '''



        # Check if the triangulated data set exists. If it does not, go to the next iteration
        if not( os.path.exists(triangulated_path) ):
            continue

        # Load the raw data set. If an error is encountered (no or not enough data points),
        # go to the next dataset.
        try:
            raw_data, nbfb, empty_files_list, nbfg, nbpg = load_data.load_raw( raw_path, nbfb, empty_files_list, nbfg, nbpg )
            startX  = raw_data['startX']  # Starting X positions (px)
            startY  = raw_data['startY']  # Starting Y positions (px)
            endX    = raw_data['endX']    # Ending X positions (px)
            endY    = raw_data['endY']    # Ending Y positions (px)
            sLat    = raw_data['sLat']    # Starting latitudes
            sLon    = raw_data['sLon']    # Starting longitudes
            eLat    = raw_data['eLat']    # Ending latitudes
            eLon    = raw_data['eLon']    # Ending longitudes
            sat     = (0)*('rcm' in raw_path)+(1)*('s1' in raw_path)+(2)*('rcm_new' in raw_path) # If we have 'rcm_new' in the path, it will return both 0+2 = 2, so it works for now, but we may want change it eventually.

        except load_data.DataFileError as dfe:
            # The error has already been printed in the triangulation stage
            continue

        #print(raw_path)
        #print(triangulated_path)
        #print(len(startX))
        # Load the triangulated dataset
        triangulated_data = load_data.load_triangulated( triangulated_path )
        vertice_idx1 = triangulated_data['vertice_idx1'] # Vertex indices in raw data file
        vertice_idx2 = triangulated_data['vertice_idx2']
        vertice_idx3 = triangulated_data['vertice_idx3']
        #print(vertice_idx1)

        # Retrieve the starting and ending times and compute the time interval (days)
        start, end = utils_datetime.dataDatetimes(raw_path)
        dt = utils_datetime.dT( (start, end) ) / 86400

        # Create a header and a list of data rows that will be used to create the output csv file
        header = ['no.', 'Sat',  'A', 'dudx', 'dudy', 'dvdx', 'dvdy', 'eps_I', 'eps_II', 'vrt', 'eps_tot']
        row_list = [header]

        # Get the number of triangles in the current dataset
        num_tri = len(vertice_idx1)

        # Retrieve the starting and ending X/Y positions of each triangle vertices
        sX1 = np.array(startX)[vertice_idx1]
        sX2 = np.array(startX)[vertice_idx2]
        sX3 = np.array(startX)[vertice_idx3]

        sY1 = np.array(startY)[vertice_idx1]
        sY2 = np.array(startY)[vertice_idx2]
        sY3 = np.array(startY)[vertice_idx3]

        eX1 = np.array(endX)[vertice_idx1]
        eX2 = np.array(endX)[vertice_idx2]
        eX3 = np.array(endX)[vertice_idx3]

        eY1 = np.array(endY)[vertice_idx1]
        eY2 = np.array(endY)[vertice_idx2]
        eY3 = np.array(endY)[vertice_idx3]

        '''
        _________________________________________________________________________________________
        COMPUTE DEFORMATIONS FOR EACH TRIANGLE
        '''

        # Iterate through every triangle in the triangulated csv file
        for n in range(num_tri):

            # Create lists of triangular cell vertices' positions
            sx_list = np.array([sX1[n], sX2[n], sX3[n]])  # Starting x positions
            sy_list = np.array([sY1[n], sY2[n], sY3[n]])  # Starting y positions
            ex_list = np.array([eX1[n], eX2[n], eX3[n]] ) # Ending x positions
            ey_list = np.array([eY1[n], eY2[n], eY3[n]])  # Ending y positions

            # Write the vertices' IDs and triangle ID to lists
            idx1.append(vertice_idx1[n])
            idx2.append(vertice_idx2[n])
            idx3.append(vertice_idx3[n])
            no.append(file_num)

            # Create a list of velocity components for each triangle vertex
            u_list, v_list = calculate_uv_lists( sx_list, ex_list, sy_list, ey_list, dt)

            # Compute the strain rates
            dudx, dudy, dvdx, dvdy, A = calculate_strainrates( u_list, v_list, sx_list, sy_list )

            # Compute the divergence rate
            eps_I = dudx + dvdy

            # Compute the maximum shear strain rate
            eps_II = sqrt(  (dudx - dvdy)**2 + (dudy + dvdx)**2  )

            # Compute the vorticity
            rot = dvdx - dudy

            # Compute the total sea-ice deformation rate
            eps_tot = sqrt( eps_I**2 + eps_II**2 )


            '''
            _________________________________________________________________________________________
            POPULATE CSV ROW LIST
            '''

            # Add the data row corresponding to the current triangle to the list of data rows
            row_list.append( [n, sat, A, dudx, dudy, dvdx, dvdy, eps_I, eps_II, rot, eps_tot] )


            '''
            _________________________________________________________________________________________
            POPULATE NETCDF ARRAYS
            '''
            # Add the satellite to the netcdf lists
            sat_l.append(sat)

            # Add the divergence and the shear strain rates to the netcdf lists
            div.append(eps_I)
            shr.append(eps_II)
            vrt.append(rot)

            # Add the strain rates and area to the netcdf lists
            A_l.append(A)
            dudx_l.append(dudx)
            dudy_l.append(dudy)
            dvdx_l.append(dvdx)
            dvdy_l.append(dvdy)

        # Update file counter
        file_num += 1

        # Add the starting and ending times (in seconds since the reference time)
        # to the times list
        s = utils_datetime.dT((refTime, start))
        sTime.extend([s for i in range(num_tri)])

        e = utils_datetime.dT((refTime, end))
        eTime.extend([e for i in range(num_tri)])


        # Add the starting and ending Lat/Lon positions of each triangle vertices
        # to the netcdf list
        sLat1.extend(np.array(sLat)[vertice_idx1])
        sLat2.extend(np.array(sLat)[vertice_idx2])
        sLat3.extend(np.array(sLat)[vertice_idx3])

        sLon1.extend(np.array(sLon)[vertice_idx1])
        sLon2.extend(np.array(sLon)[vertice_idx2])
        sLon3.extend(np.array(sLon)[vertice_idx3])

        eLat1.extend(np.array(eLat)[vertice_idx1])
        eLat2.extend(np.array(eLat)[vertice_idx2])
        eLat3.extend(np.array(eLat)[vertice_idx3])

        eLon1.extend(np.array(eLon)[vertice_idx1])
        eLon2.extend(np.array(eLon)[vertice_idx2])
        eLon3.extend(np.array(eLon)[vertice_idx3])

        idx_sLat1.extend(vertice_idx1)
        idx_sLat2.extend(vertice_idx2)
        idx_sLat3.extend(vertice_idx3)

        '''
        _________________________________________________________________________________________
        WRITE THE CURRENT DATASET'S RESULTS TO A CSV FILE (DEPRECATED)
        '''
        # If the calculations file already exists and overwrite (in namelist.ini) is set to 'no',
        # do not write the calculations csv file
        # if not ( os.path.exists(calculations_path) and not config.config['Processing_options'].getboolean('overwrite')):
            # Create a directory to store the calculations csv path if it does not exist already
            # os.makedirs(os.path.dirname(calculations_path), exist_ok=True)

            # Write the results in the calculations_csv_path file path
            # with open(calculations_path, 'w', encoding='UTF8', newline='') as f:
                # writer = csv.writer(f)

                # Write the data rows to the csv file
                # writer.writerows(row_list)

    '''
    _________________________________________________________________________________________
    WRITE ALL RESULTS COMBINED TO A NETCDF FILE
    '''
    # Find absolute path in which the output netcdf file is to be stored
    nc_output_path = config['data_paths']['nc_output']
    # Create a directory to store the output netcdf file if it does not exist already
    os.makedirs(os.path.dirname(nc_output_path), exist_ok=True)

    print("Printing the deformations in the netcdf: ", nc_output_path)
    # Create an output netcdf file and dataset
    output_ds = Dataset(nc_output_path, 'w', format = 'NETCDF4')

    # Add metadata
    Metadata = config['Metadata']
    output_ds.icetracker = Metadata['icetracker']
    output_ds.referenceTime = YYYY + '-' + MM + '-' + DD + ' 00:00:00'
    output_ds.trackingError = Metadata['tracking_error'] + ' m'
    output_ds.timestep = Date_options['timestep'] + ' hours'
    output_ds.tolerance = Date_options['tolerance'] + ' hours'

    # Create x array to store output data
    x = output_ds.createDimension('x', len(sTime))

    # Create variables for netcdf data set
    start_time = output_ds.createVariable('start_time', 'u4', 'x') # Start and end times
    end_time   = output_ds.createVariable('end_time', 'u4', 'x')

    satellite  = output_ds.createVariable('satellite', 'u4', 'x') # Satellite

    start_lat1 = output_ds.createVariable('start_lat1', 'f8', 'x') # Starting Lat/Lon triangle vertices
    start_lat2 = output_ds.createVariable('start_lat2', 'f8', 'x')
    start_lat3 = output_ds.createVariable('start_lat3', 'f8', 'x')
    start_lon1 = output_ds.createVariable('start_lon1', 'f8', 'x')
    start_lon2 = output_ds.createVariable('start_lon2', 'f8', 'x')
    start_lon3 = output_ds.createVariable('start_lon3', 'f8', 'x')

    end_lat1   = output_ds.createVariable('end_lat1', 'f8', 'x') # Ending Lat/Lon triangle vertices
    end_lat2   = output_ds.createVariable('end_lat2', 'f8', 'x')
    end_lat3   = output_ds.createVariable('end_lat3', 'f8', 'x')
    end_lon1   = output_ds.createVariable('end_lon1', 'f8', 'x')
    end_lon2   = output_ds.createVariable('end_lon2', 'f8', 'x')
    end_lon3   = output_ds.createVariable('end_lon3', 'f8', 'x')

    d          = output_ds.createVariable('div', 'f8', 'x') # Divergence and shear strain and vorticity rates
    s          = output_ds.createVariable('shr', 'f8', 'x')
    v          = output_ds.createVariable('vrt', 'f8', 'x')

    id1        = output_ds.createVariable('idx1', 'u4', 'x') # Triangle vertices
    id2        = output_ds.createVariable('idx2', 'u4', 'x')
    id3        = output_ds.createVariable('idx3', 'u4', 'x')
    idtri      = output_ds.createVariable('no', 'u4', 'x')

    Aa         = output_ds.createVariable('A', 'f8', 'x') # Triangle area

    dux        = output_ds.createVariable('dudx', 'f8', 'x') # Strain rates
    duy        = output_ds.createVariable('dudy', 'f8', 'x')
    dvx        = output_ds.createVariable('dvdx', 'f8', 'x')
    dvy        = output_ds.createVariable('dvdy', 'f8', 'x')

    # Specify units for each variable
    start_time.units = 'seconds since the reference time'
    end_time.units   = 'seconds since the reference time'

    satellite.units   = '0: RCM; 1: S1'

    start_lat1.units = 'degrees North'
    start_lat2.units = 'degrees North'
    start_lat3.units = 'degrees North'
    start_lon1.units = 'degrees East'
    start_lon2.units = 'degrees East'
    start_lon3.units = 'degrees East'

    end_lat1.units   = 'degrees North'
    end_lat2.units   = 'degrees North'
    end_lat3.units   = 'degrees North'
    end_lon1.units   = 'degrees East'
    end_lon2.units   = 'degrees East'
    end_lon3.units   = 'degrees East'

    d.units          = '1/days'
    s.units          = '1/days'
    v.units          = '1/days'

    dux.units       = '1/days'
    duy.units       = '1/days'
    dvx.units       = '1/days'
    dvy.units       = '1/days'

    Aa.units         = 'square meters'

    # Attribute data arrays to each variable
    start_time[:] = sTime
    end_time[:]   = eTime

    satellite[:]  = sat_l

    start_lat1[:] = sLat1
    start_lat2[:] = sLat2
    start_lat3[:] = sLat3
    start_lon1[:] = sLon1
    start_lon2[:] = sLon2
    start_lon3[:] = sLon3

    end_lat1[:]   = eLat1
    end_lat2[:]   = eLat2
    end_lat3[:]   = eLat3
    end_lon1[:]   = eLon1
    end_lon2[:]   = eLon2
    end_lon3[:]   = eLon3

    d[:]          = div
    s[:]          = shr
    v[:]          = vrt

    id1[:]        = idx1
    id2[:]        = idx2
    id3[:]        = idx3
    idtri[:]      = no

    Aa[:]        = [Ai*(int(Metadata['tracking_error'])**2.) for Ai in A_l]

    dux[:]       = dudx_l
    duy[:]       = dudy_l
    dvx[:]       = dvdx_l
    dvy[:]       = dvdy_l

    return output_ds



def calculate_uv_lists( sx_list, ex_list, sy_list, ey_list, dT):
    ''' (list, list, float) -> float(list, list)

    Function that calculates u and v velocity components of
    each cell vertices and returns them as a list.

    Returns a tuple (u_list, v_list) of the lists of velocity components.

    Keyword arguments: \\
    sx_list -- list of starting x positions of each cell vertices \\
    ex_list -- list of ending x positions of each cell vertices \\
    sy_list -- list of starting x positions of each cell vertices \\
    ey_list -- list of ending x positions of each cell vertices \\
    dT      -- time interval
    '''

    # Compute the u and v velocity components at the current vertex
    u_list = (ex_list - sx_list) / dT
    v_list = (ey_list - sy_list) / dT

    return u_list, v_list


def calculate_strainrates( u_list, v_list, sx_list, sy_list ):
    ''' (list, list, list, list) -> tuple(float, float, float, float)

    Computes the strain rates (or velocity derivatives).

    Keyword arguments: \\
    u_list -- list of u component velocities for each cell vertices \\
    v_list -- list of v component velocities for each cell vertices \\
    sx_list -- list of starting x positions for each cell vertices \\
    sy_list -- list of starting y positions for each cell vertices \\
    '''

    # Find the number of cell vertices
    n = len(sx_list)

    #----- Compute the Lagrangian cell area A -----

    # Initialize the cell area to 0
    A = 0.0

    # Perform a summation to compute A (see Bouchat et al. (2020) eqn. 6)
    for i in range( n ):
        A += (1.0/2.0) * ( sx_list[i] * sy_list[((i+1) % n)] - sx_list[((i+1) % n)] *  sy_list[i])

    #----- Compute the strain rates ---------------

    # Initialize the strain rates to 0
    dudx = 0.0
    dudy = 0.0
    dvdx = 0.0
    dvdy = 0.0

    # Perform a summation to compute strain rates (see Bouchat et al. (2020) eqn. 5)
    for i in range( n ):

        dudx += 1.0/(2.0*A)  * ( u_list[((i+1) % n)] + u_list[i] ) * ( sy_list[((i+1) % n)] - sy_list[i] )

        dudy += -1.0/(2.0*A) * ( u_list[((i+1) % n)] + u_list[i] ) * ( sx_list[((i+1) % n)] - sx_list[i] )

        dvdx += 1.0/(2.0*A)  * ( v_list[((i+1) % n)] + v_list[i] ) * ( sy_list[((i+1) % n)] - sy_list[i] )

        dvdy += -1.0/(2.0*A) * ( v_list[((i+1) % n)] + v_list[i] ) * ( sx_list[((i+1) % n)] - sx_list[i] )

    return dudx, dudy, dvdx, dvdy, A


if __name__ == '__main__':
    compute_deformations()
