'''
Author: Beatrice Duval (bdu002)

-----------------------------------------------------------------
Utils - Load data
-----------------------------------------------------------------

Helper code for loading files from all stages of data processing.
'''

# Loading from default packages
import os
import pandas as pd
import numpy as np

# Create a class of errors for csv loading and two subclasses
class DataFileError(Exception):
    pass

class FileFormatError(DataFileError):
    pass

class DataError(DataFileError):
    pass

def load_raw(path_raw, nbfb, empty_files, nbfg, nbpg):
    ''' (str) -> dict[str, list]

    Loads data from a raw .dat file.

    Returns a dictionnary containing sLat, sLon, eLat, eLon, startX, startY,
    endX, endY, dispX, and dispY.

    An error is raised if the file is not a .dat file, or if the file
    is empty or contains less than 3 data points (not enough to perform a triangulation).

    Keyword arguments: \\
    path_raw -- absolute path to raw .dat file
    '''

    # Determine if the raw file is a .dat file
    file_type = path_raw[-3:len(path_raw)]
    # If the file is not a .dat file, raise an error

    if file_type == 'dat':
        # Create a data frame
        df = pd.read_csv(path_raw, sep='\s\s+', engine='python')

        # Retrieve data points
        sLat    = df['sLat']    # Starting latitudes
        sLon    = df['sLon']    # Starting longitudes
        eLat    = df['eLat']    # Ending latitudes
        eLon    = df['eLon']    # Ending longitudes
        startX  = df['startX']  # Starting X positions (px)
        startY  = df['startY']  # Starting Y positions (px)
        endX    = df['endX']    # Ending X positions (px)
        endY    = df['endY']    # Ending Y positions (px)
        dispX   = df['dispX']   # X displacement (px)
        dispY   = df['dispY']   # Y displacement (px)

        # Raise an error if the input raw csv file is empty or if it contains less than 3 data points
        if ( len(sLon) < 4 ):

            # Retrieve the raw filename and raise an error
            filename_raw_csv  = os.path.basename( path_raw )

            # raise DataError(filename_raw_csv  + " is empty or does not have enough data to perform a triangulation. ")

            nbfb += 1
            empty_files.append(filename_raw_csv)


        else: # print the name of the file and the number of data points in there
            nbfg += 1
            nbpg += len(sLon)

        # Create a dictionnary of raw data and return it
        raw_data = {'sLat': sLat, 'sLon': sLon, 'eLat': eLat, 'eLon': eLon,
                        'startX': startX, 'startY': startY, 'endX': endX, 'endY': endY,
                        'dispX': dispX, 'dispY': dispY}

    elif file_type == 'trk':

        # Create a data frame
        df = pd.DataFrame()
        with open(path_raw) as fd:
            headers = [ next(fd) for i in range(9) ]
            df = pd.read_csv(fd,engine='python',sep = '\s\s+')
        df.drop([0], axis=0, inplace=True)
        nb = len(df.index)
        df.reset_index(inplace=True)

        # Retrieve data points
        sLat    = df['lat_beg']    # Starting latitudes
        sLon    = df['lon_beg']    # Starting longitudes
        eLat    = df['lat_end']    # Ending latitudes
        eLon    = df['lon_end']    # Ending longitudes
        startX  = df['pix']  # Starting X positions (px)
        startY  = df['lin'] # Starting Y positions (px)
        dispX   = df['dpix']   # X displacement (px)
        dispY   = df['dlin']   # Y displacement (px)
        endX    = startX.astype(float) + dispX.astype(float)    # Ending X positions (px)
        endY    = startY.astype(float) + dispY.astype(float)    # Ending Y positions (px)

        sLat = sLat.astype(float)
        sLon = sLon.astype(float)
        eLat = eLat.astype(float)
        eLon = eLon.astype(float)
        startX = startX.astype(float)
        startY = startY.astype(float)
        endX = endX.astype(float)
        endY = endY.astype(float)


        # Raise an error if the input raw csv file is empty or if it contains less than 3 data points
        if ( len(sLon) < 4 ):

            # Retrieve the raw filename and raise an error
            filename_raw_csv  = os.path.basename( path_raw )

            # raise DataError(filename_raw_csv  + " is empty or does not have enough data to perform a triangulation. ")

            nbfb += 1
            empty_files.append(filename_raw_csv)

        else: # print the name of the file and the number of data points in there
            nbfg += 1
            nbpg += len(sLon)

        # Create a dictionnary of raw data and return it
        raw_data = {'sLat': sLat, 'sLon': sLon, 'eLat': eLat, 'eLon': eLon,
                        'startX': startX, 'startY': startY, 'endX': endX, 'endY': endY,
                        'dispX': dispX, 'dispY': dispY}


    else:
        raise FileFormatError(os.path.basename(path_raw) + \
                        ' is not a .dat or .trk file.')

    return raw_data, nbfb, empty_files, nbfg, nbpg



def load_triangulated( path_triangulated ):
    ''' (str) -> dict[str, list]

    Loads data from a triangulated .csv file.

    Returns a dictionnary of arrays of triangle numbers, and the vertices'
        index in the raw file: (no, vertice_idx1, vertice_idx2, vertice_idx3).

    Keyword arguments: \\
    path_triangulated -- absolute path to processed .csv file

    '''
    # Create a data frame
    df = pd.read_csv(path_triangulated)

    # Retrieve data
    no           = df['no.']          # Triangle number

    # Create a return dictionnary of triangulated data

    vertice_idx1 = df['vertice_idx1'] # Vertex indices in raw file
    vertice_idx2 = df['vertice_idx2']
    vertice_idx3 = df['vertice_idx3']

    tri_data = {'no': no, 'vertice_idx1': vertice_idx1, 'vertice_idx2': vertice_idx2, 'vertice_idx3': vertice_idx3}

    return tri_data


def load_converted( path_converted ):
    ''' (str) -> dict[str, list]

    Loads data from a converted .csv file. Returns a dictionnary containing arrays
    of triangle numbers, starting and ending (x,y) positions of the vertices in a
    local cartesian coordinate system, and the (j, i) indices of the nearest grid
    tracer point that defines the local cartesian coordinate system (no, sX1,sX2,
    sX3, sY1, sY2, sY3, eX1, eX2, eX3, eY1, eY2, eY3, j_tracers, i_tracers).

    Keyword arguments: \\
    path_converted -- absolute path to converted .csv file
    '''
    # Create a data frame
    df = pd.read_csv( path_converted )

    # Retrieve data
    no          = df['no.']             # Triangle number

    sX1         = df['sX1']             # Starting x coordinates in a local grid CS
    sX2         = df['sX2']
    sX3         = df['sX3']

    sY1         = df['sY1']             # Starting y coordinates in a local grid CS
    sY2         = df['sY2']
    sY3         = df['sY3']

    eX1         = df['eX1']             # Ending x coordinates in a local grid CS
    eX2         = df['eX2']
    eX3         = df['eX3']

    eY1         = df['eY1']             # Ending y coordinates in a local grid CS
    eY2         = df['eY2']
    eY3         = df['eY3']

    j_tracers   = df['tracer_j']        # Tracer point (which defines the local CS) indices
    i_tracers   = df['tracer_i']


    converted_data = {'no': no, 'sX1': sX1, 'sX2': sX2, 'sX3': sX3, 'eX1': eX1, 'eX2': eX2, 'eX3': eX3,
                                'sY1': sY1, 'sY2': sY2, 'sY3': sY3, 'eY1': eY1, 'eY2': eY2, 'eY3': eY3,
                                'j_tracers': j_tracers, 'i_tracers': i_tracers}

    return converted_data


def load_calculations( path_calculations ):
    ''' (str) -> dict[str, list]

    Loads data from a calculations .csv file. Returns a dictionnary of arrays
    of triangle numbers, strain rates, divergence rate, maximum shear strain
    rate, and total deformation rate (no, dudx, dudy, dvdx, dvdy, eps_I, eps_II, eps_tot).

    Keyword arguments:
    path_calculations -- absolute path to calculations .csv file
    '''
    # Create a data frame
    df = pd.read_csv(path_calculations)

    # Retrieve data
    no          = df['no.']             # Triangle number

    dudx        = df['dudx']            # Strain rates
    dudy        = df['dudy']
    dvdx        = df['dvdx']
    dvdy        = df['dvdy']

    eps_I       = df['eps_I']           # Divergence rate

    eps_II      = df['eps_II']          # Maximum shear strain rate

    eps_tot     = df['eps_tot']         # Total deformation rate

    calculations_data = {'no': no, 'dudx': dudx, 'dudy': dudy, 'dvdx': dvdx, 'dvdy': dvdy,
                                   'eps_I': eps_I, 'eps_II': eps_II, 'eps_tot': eps_tot}

    return calculations_data
