'''
Author: Beatrice Duval (bdu002)

----------------------------------------------------------------------
Helper code for loading csv files from all stages of data processing.
----------------------------------------------------------------------

'''

import os
import pandas as pd

def load_raw_csv( path_raw_csv ):
    ''' (string) -> (array, array, array, array)

    Loads data from a raw .csv file. Returns starting and ending point
    positions in latitudes and longitudes: (sLat, sLon, eLat, eLon).

    Keyword arguments:
    path_raw_csv -- absolute path to raw .csv file
    '''
    # Create a data frame
    df = pd.read_csv(path_raw_csv)

    # Retrieve data points
    sLat = df['sLat']       # Starting latitudes
    sLon = df['sLon']       # Starting longitudes
    eLat = df['eLat']       # Ending latitudes
    eLon = df['eLon']       # Ending longitudes
    

    # Raise an error if the input raw csv file is empty
    if ( len(sLon) < 3 ):

        # Retrieve the raw csv filename
        filename_raw_csv  = os.path.basename( path_raw_csv )

        # Raise an error
        raise ValueError("The input raw csv file " + filename_raw_csv  + \
                         " is empty or does not have enough data to perform a triangulation." + \
                         "This file will not be processed.")

    return sLat, sLon, eLat, eLon



def load_processed_csv( path_processed_csv ):
    ''' (string) -> (array, array, array, array, array, array, array, array, array, array)

    Loads data from a processed .csv file. Returns an array of triangle numbers, 
    starting and ending positions of the vertices in the aeqd transform, 
    and the vertices' index in the raw csv file: (no, sX1_aeqd, sX2_aeqd, sX3_aeqd, 
    sY1_aeqd, sY2_aeqd, sY3_aeqd,  vertice_idx1, vertice_idx2, vertice_idx3).

    Keyword arguments:
    path_processed_csv -- absolute path to processed .csv file
    '''
    # Create a data frame
    df = pd.read_csv(path_processed_csv)

    # Retrieve data
    no           = df['no.']          # Triangle number

    sX1_aeqd     = df['sX1_aeqd']     # Starting X coordinates in the aeqd transform
    sX2_aeqd     = df['sX2_aeqd']
    sX3_aeqd     = df['sX3_aeqd']

    sY1_aeqd     = df['sY1_aeqd']     # Starting Y coordinates in the aeqd transform
    sY2_aeqd     = df['sY2_aeqd']
    sY3_aeqd     = df['sY3_aeqd']

    vertice_idx1 = df['vertice_idx1'] # Vertex indices in raw csv file
    vertice_idx2 = df['vertice_idx2']
    vertice_idx3 = df['vertice_idx3']

    return (no, sX1_aeqd, sX2_aeqd, sX3_aeqd, 
                sY1_aeqd, sY2_aeqd, sY3_aeqd, 
                vertice_idx1, vertice_idx2, vertice_idx3)


def load_converted_csv( path_converted_csv ):
    ''' (string) -> (array, array, array, array, array, array, array, array, array, array)

    Loads data from a converted .csv file. Returns an array of triangle numbers, 
    starting and ending (x,y) positions of the vertices in a local cartesian 
    coordinate system, and the (j, i) indices of the nearest grid tracer point 
    that defines the local cartesian coordinate system: (no, sX1,sX2, sX3, sY1, 
    sY2, sY3, eX1, eX2, eX3, eY1, eY2, eY3, j_tracers, i_tracers).

    Keyword arguments:
    path_converted_csv -- absolute path to converted .csv file
    '''
    # Create a data frame
    df = pd.read_csv( path_converted_csv )

    # Retrieve data
    no          = df['no.']             # Triangle number

    sX1         = df['sX1']             # Starting x coordinates in a local grid CS
    sX2         = df['sX2']
    sX3         = df['sX3']

    sY1         = df['sY1']             # Starting y coordinates in a local grid CS
    sY2         = df['sY2']
    sY3         = df['sY3']

    eX1         = df['eX1']             # Starting x coordinates in a local grid CS
    eX2         = df['eX2']
    eX3         = df['eX3']

    eY1         = df['eY1']             # Starting y coordinates in a local grid CS
    eY2         = df['eY2']
    eY3         = df['eY3']

    j_tracers   = df['tracer_j']        # Tracer point (defining the local CS) indices 
    i_tracers   = df['tracer_i']

    return no, sX1,sX2, sX3, sY1, sY2, sY3, eX1, eX2, eX3, eY1, eY2, eY3, j_tracers, i_tracers


def load_calculations_csv( path_calculations_csv ):
    ''' (string) -> (array, array, array, array, array, array, array, array, array, array)

    Loads data from a calculations .csv file. Returns an array of triangle numbers, 
    strain rates, divergence rate, maximum shear strain rate, and total deformation rate:
    (no, dudx, dudy, dvdx, dvdy, eps_I, eps_II, eps_tot).

    Keyword arguments:
    path_calculations_csv -- absolute path to calculations .csv file
    '''
    # Create a data frame
    df = pd.read_csv(path_calculations_csv)

    # Retrieve data
    no          = df['no.']             # Triangle number

    dudx        = df['dudx']            # Strain rates
    dudy        = df['dudy']
    dvdx        = df['dvdx']
    dvdy        = df['dvdy']

    eps_I       = df['eps_I']           # Divergence rate

    eps_II      = df['eps_II']          # Maximum shear strain rate

    eps_tot     = df['eps_tot']         # Total deformation rate

    return no, dudx, dudy, dvdx, dvdy, eps_I, eps_II, eps_tot