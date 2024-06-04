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
        div, shr, vrt, ids1, ids2, ids3, \
        idpair, ids_sLat1, ids_sLat2, ids_sLat3, A_l, \
        dudx_l, dudy_l, dvdx_l, dvdy_l, sat_l , \
        delA_l, delI_l, delII_l, delvrt_l, s2n \
        = ([] for i in range(35))
    sigx = 200.0

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
    file_tot = 0
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


        file_tot += 1
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


        # Retrieve the starting and ending times and compute the time interval (days)
        start, end = utils_datetime.dataDatetimes(raw_path)
        dt = utils_datetime.dT( (start, end) ) / 86400.0
        if dt < 0.5:
            print(dt)
            continue

        # Load the triangulated dataset
        triangulated_data = load_data.load_triangulated( triangulated_path )
        vertice_ids1 = triangulated_data['vertice_idx1'] # Vertex indices in raw data file
        vertice_ids2 = triangulated_data['vertice_idx2']
        vertice_ids3 = triangulated_data['vertice_idx3']
        #print(vertice_idx1)

        # Create a header and a list of data rows that will be used to create the output csv file
        header = ['no.', 'Sat',  'A', 'dudx', 'dudy', 'dvdx', 'dvdy', 'eps_I', 'eps_II', 'vrt', 'eps_tot']
        row_list = [header]

        # Get the number of triangles in the current dataset
        num_tri = len(vertice_ids1)

        # Retrieve the starting and ending X/Y positions of each triangle vertices
        sX1 = np.array(startX)[vertice_ids1]
        sX2 = np.array(startX)[vertice_ids2]
        sX3 = np.array(startX)[vertice_ids3]

        sY1 = np.array(startY)[vertice_ids1]
        sY2 = np.array(startY)[vertice_ids2]
        sY3 = np.array(startY)[vertice_ids3]

        eX1 = np.array(endX)[vertice_ids1]
        eX2 = np.array(endX)[vertice_ids2]
        eX3 = np.array(endX)[vertice_ids3]

        eY1 = np.array(endY)[vertice_ids1]
        eY2 = np.array(endY)[vertice_ids2]
        eY3 = np.array(endY)[vertice_ids3]

        '''
        _________________________________________________________________________________________
        COMPUTE DEFORMATIONS FOR EACH TRIANGLE
        '''

        # Iterate through every triangle in the triangulated csv file
        for n in range(num_tri):

            # Create lists of triangular cell vertices' positions
            sx_list = np.array([sX1[n], sX2[n], sX3[n]])*sigx  # Starting x positions
            sy_list = np.array([sY1[n], sY2[n], sY3[n]])*sigx  # Starting y positions
            ex_list = np.array([eX1[n], eX2[n], eX3[n]])*sigx # Ending x positions
            ey_list = np.array([eY1[n], eY2[n], eY3[n]])*sigx  # Ending y positions


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

            # Compute the propagation errors
            del11, del12, del21, del22, delA = calculate_trackerrors( u_list, v_list, sx_list, sy_list,dt)

            del11 = dudx*dudx*delA/A**2.0 + del11*(sigx**2.0)
            del22 = dvdy*dvdy*delA/A**2.0 + del22*(sigx**2.0)
            del12 = dudy*dudy*delA/A**2.0 + del12*(sigx**2.0)
            del21 = dvdx*dvdx*delA/A**2.0 + del21*(sigx**2.0)
            delI = del11 + del22
            delrot = del12 + del21
            term1 = (delI*(dudx - dvdy)**2.0 + (del12+del21)*(dudy + dvdx)**2.0) #this is eps_II^2 * delII^2
            term2 = eps_I**2.0 * delI
            epstot_deltot = term1 + term2

            if eps_II > 0.0:
                delII = (delI*(((dudx - dvdy)/eps_II)**2.0) + (del12+del21)*(((dudy + dvdx)/eps_II)**2.0))
                deltot = delII*((eps_II/eps_tot)**2.0) + delI*(( eps_I/eps_tot)**2.0)
                sig2noise = eps_tot**2.0 / epstot_deltot**0.5
            else:
                delII = 1000.0
                deltot = 1000.0
                sig2noise = 0.0

            '''
            _________________________________________________________________________________________
            Add the computed data to the output vectors
            '''
            if (A**0.5 < 20000.0):
                # Write the vertices' IDs and triangle ID to lists
                ids1.append(vertice_ids1[n])
                ids2.append(vertice_ids2[n])
                ids3.append(vertice_ids3[n])
                idpair.append(file_num)

                # Add the data row corresponding to the current triangle to the list of data rows
                row_list.append( [n, sat, A, dudx, dudy, dvdx, dvdy, eps_I, eps_II, rot, eps_tot] )

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

                # Add the strain rates and area to the netcdf lists
                delA_l.append(delA**0.5)
                delI_l.append(delI**0.5)
                delII_l.append(delII**0.5)
                delvrt_l.append(delrot**0.5)
                s2n.append(sig2noise)


                # Add the starting and ending Lat/Lon positions of each triangle vertices
                # to the netcdf list
                sLat1.append(np.array(sLat)[vertice_ids1[n]])
                sLat2.append(np.array(sLat)[vertice_ids2[n]])
                sLat3.append(np.array(sLat)[vertice_ids3[n]])

                sLon1.append(np.array(sLon)[vertice_ids1[n]])
                sLon2.append(np.array(sLon)[vertice_ids2[n]])
                sLon3.append(np.array(sLon)[vertice_ids3[n]])

                eLat1.append(np.array(eLat)[vertice_ids1[n]])
                eLat2.append(np.array(eLat)[vertice_ids2[n]])
                eLat3.append(np.array(eLat)[vertice_ids3[n]])

                eLon1.append(np.array(eLon)[vertice_ids1[n]])
                eLon2.append(np.array(eLon)[vertice_ids2[n]])
                eLon3.append(np.array(eLon)[vertice_ids3[n]])

                ids_sLat1.append(vertice_ids1[n])
                ids_sLat2.append(vertice_ids2[n])
                ids_sLat3.append(vertice_ids3[n])


            else:
                num_tri = num_tri - 1

        # Update file counter
        file_num += 1

        # Add the starting and ending times (in seconds since the reference time)
        # to the times list
        s = utils_datetime.dT((refTime, start))
        sTime.extend([s for i in range(num_tri)])

        e = utils_datetime.dT((refTime, end))
        eTime.extend([e for i in range(num_tri)])


    '''
    _________________________________________________________________________________________
    WRITE ALL RESULTS COMBINED TO A NETCDF FILE
    '''
    print(file_tot,file_num)
    # Find absolute path in which the output netcdf file is to be stored
    nc_output_path = config['data_paths']['nc_output']
    # Create a directory to store the output netcdf file if it does not exist already
    os.makedirs(os.path.dirname(nc_output_path), exist_ok=True)

    print("Printing the deformations in the netcdf: ", nc_output_path)

    # Create an output netcdf file and dataset
    output_ds = Dataset(nc_output_path, 'w', format = 'NETCDF4')

    # Add metadata as global attribute to the netcdf
    Metadata = config['Metadata']
    maxDeltat = int(Date_options['tolerance']) + int(Date_options['timestep'])
    if Metadata['icetracker'] == 'RCMS1':
        SARsource = "RCM and S1"
    else:
        SARsource = Metadata['icetracker']
    print(SARsource)
    Metadata = config['Metadata']
    output_ds.SatelliteSource = SARsource
    output_ds.referenceTime = YYYY + '-' + MM + '-' + DD + ' 00:00:00'
    output_ds.trackingError = Metadata['tracking_error'] + ' m'
    output_ds.Max_Deltat = '%s hours' % maxDeltat
    #output_ds.tolerance = Date_options['tolerance'] + ' hours'

    # Create x array to store output data
    n = output_ds.createDimension('n', len(sTime))

    # Create variables for netcdf data set
    start_time = output_ds.createVariable('start_time', 'u4', 'n') # Start and end times
    end_time   = output_ds.createVariable('end_time', 'u4', 'n')

    satellite  = output_ds.createVariable('satellite', 'u4', 'n') # Satellite

    start_lat1 = output_ds.createVariable('start_lat1', 'f8', 'n') # Starting Lat/Lon triangle vertices
    start_lat2 = output_ds.createVariable('start_lat2', 'f8', 'n')
    start_lat3 = output_ds.createVariable('start_lat3', 'f8', 'n')
    start_lon1 = output_ds.createVariable('start_lon1', 'f8', 'n')
    start_lon2 = output_ds.createVariable('start_lon2', 'f8', 'n')
    start_lon3 = output_ds.createVariable('start_lon3', 'f8', 'n')

    end_lat1   = output_ds.createVariable('end_lat1', 'f8', 'n') # Ending Lat/Lon triangle vertices
    end_lat2   = output_ds.createVariable('end_lat2', 'f8', 'n')
    end_lat3   = output_ds.createVariable('end_lat3', 'f8', 'n')
    end_lon1   = output_ds.createVariable('end_lon1', 'f8', 'n')
    end_lon2   = output_ds.createVariable('end_lon2', 'f8', 'n')
    end_lon3   = output_ds.createVariable('end_lon3', 'f8', 'n')

    d          = output_ds.createVariable('div', 'f8', 'n') # Divergence and shear strain and vorticity rates
    s          = output_ds.createVariable('shr', 'f8', 'n')
    v          = output_ds.createVariable('vrt', 'f8', 'n')

    id1        = output_ds.createVariable('ids1', 'u4', 'n') # Triangle vertices
    id2        = output_ds.createVariable('ids2', 'u4', 'n')
    id3        = output_ds.createVariable('ids3', 'u4', 'n')
    id_pair     = output_ds.createVariable('idpair', 'u4', 'n')

    Aa         = output_ds.createVariable('A', 'f8', 'n') # Triangle area

    dux        = output_ds.createVariable('dudx', 'f8', 'n') # Strain rates
    duy        = output_ds.createVariable('dudy', 'f8', 'n')
    dvx        = output_ds.createVariable('dvdx', 'f8', 'n')
    dvy        = output_ds.createVariable('dvdy', 'f8', 'n')


    sA         = output_ds.createVariable('errA', 'f8', 'n') # Triangle area
    sI         = output_ds.createVariable('err_div', 'f8', 'n') # Strain rates
    sII        = output_ds.createVariable('err_shr', 'f8', 'n')
    svrt        = output_ds.createVariable('err_vrt', 'f8', 'n')
    sig2n      = output_ds.createVariable('s2n', 'f8', 'n')


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

    id1.units        = 'ID vertex 1'
    id2.units        = 'ID vertex 2'
    id3.units        = 'ID vertex 3'
    id_pair.units     = 'ID SAR pair'

    d.units          = '1/days'
    s.units          = '1/days'
    v.units          = '1/days'

    dux.units        = '1/days'
    duy.units        = '1/days'
    dvx.units        = '1/days'
    dvy.units        = '1/days'

    Aa.units         = 'square meters'

    sI.units         = '1/days'
    sII.units        = '1/days'
    svrt.units       = '1/days'

    sA.units         = 'square meters'
    sig2n.units      = 'none'

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

    id1[:]        = ids1
    id2[:]        = ids2
    id3[:]        = ids3
    id_pair[:]      = idpair

    Aa[:]        = A_l #[Ai*(int(Metadata['tracking_error'])**2.) for Ai in A_l]

    dux[:]       = dudx_l
    duy[:]       = dudy_l
    dvx[:]       = dvdx_l
    dvy[:]       = dvdy_l

    sA[:]        = delA_l
    sI[:]        = delI_l
    sII[:]       = delII_l
    svrt[:]      = delvrt_l
    sig2n[:]     = s2n

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


def calculate_trackerrors( u_list, v_list, sx_list, sy_list, Deltat ):
    ''' (float, list, list) float

    Computes the strain rates (or velocity derivatives).

    Keyword arguments:
    u_list -- list of u component velocities for each cell vertices
    v_list -- list of v component velocities for each cell vertices
    sx_list -- list of starting x positions for each cell vertices
    sy_list -- list of starting y positions for each cell vertices
    '''

    # Find the number of cell vertices
    n = len(sx_list)

    #----- Compute the Lagrangian cell area A -----

    # Initialize the cell area to 0
    A = 0.0
    delA = 0.0

    # Perform a summation to compute A (see Bouchat et al. (2020) eqn. 6)
    for i in range( n ):
        A += (1.0/2.0) * ( sx_list[i] * sy_list[((i+1) % n)] - sx_list[((i+1) % n)] *  sy_list[i])

        delA += ((sx_list[i] - sx_list[((i+2) % n)])**2.0 + (sy_list[((i+2) % n)] - sy_list[i])**2.0)/4.0

    #----- Compute the strain rates ---------------

    # Initialize the strain rates to 0
    trackerror11 = 0.0
    trackerror12 = 0.0
    trackerror21 = 0.0
    trackerror22 = 0.0

    # Perform a summation to compute strain rates (see Bouchat et al. (2020) eqn. 5)
    for i in range( n ):

        trackerror11 += (sy_list[i] - sy_list[((i+2) % n)])**2.0/(2.0*A*A*Deltat*Deltat) + (u_list[i]-u_list[((i+2) % n)])**2.0/(4.0*A*A)

        trackerror22 += (sx_list[i] - sx_list[((i+2) % n)])**2.0/(2.0*A*A*Deltat*Deltat) + (v_list[i]-v_list[((i+2) % n)])**2.0/(4.0*A*A)

        trackerror12 += (sx_list[i] - sx_list[((i+2) % n)])**2.0/(2.0*A*A*Deltat*Deltat) + (u_list[i]-u_list[((i+2) % n)])**2.0/(4.0*A*A)

        trackerror21 += (sy_list[i] - sy_list[((i+2) % n)])**2.0/(2.0*A*A*Deltat*Deltat) + (v_list[i]-v_list[((i+2) % n)])**2.0/(4.0*A*A)


    return trackerror11, trackerror12, trackerror21, trackerror22, delA



if __name__ == '__main__':
    compute_deformations()
