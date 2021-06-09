'''
Author: Beatrice Duval (bdu002)

----------------------------------------------------------------------
Helper code for the conversion of a triangle's lat,lon coordinates to 
x,y coordinates in a local cartesian grid coordinate system.
----------------------------------------------------------------------

'''

import numpy as np
from haversine import haversine

def find_aeqdTriCenter(x1_aeqd, x2_aeqd, x3_aeqd, y1_aeqd, y2_aeqd, y3_aeqd):
    ''' (float, float, float, float, float, float) ->  tuple

    Finds the center of a triangle in Azimuthal Equidistant 
    (aeqd) transform coordinates. 
    Returns the aqed (x, y) coordinates of the center.

    Keyword arguments:
    x1_aeqd -- x coordinate of 1st triangle vertex in the aeqd transform
    x2_aeqd -- x coordinate of 2nd triangle vertex in the aeqd transform 
    x3_aeqd -- x coordinate of 3rd triangle vertex in the aeqd transform 
    y1_aeqd -- y coordinate of 1st triangle vertex in the aeqd transform
    y2_aeqd -- y coordinate of 2nd triangle vertex in the aeqd transform
    y3_aeqd -- y coordinate of 3rd triangle vertex in the aeqd transform 
    '''

    # Average the x and the y coords (respectively) to find the center
    center_y = (y1_aeqd + y2_aeqd + y3_aeqd) / 3
    center_x = (x1_aeqd + x2_aeqd + x3_aeqd) / 3

    return center_x, center_y

def find_nearestGridTracerPt(gridX_aeqd, gridY_aeqd, xy_aeqd):
    '''  (float(j,i), float(j,i), tuple) -> (int, int)

    Finds the nearest grid box for an input point. 
    Returns the (j,i) index of the nearest grid box tracer point.
    
    Note: all (x,y) points follow the Azimuthal Equidistant (aeqd) transform.
    
    Keyword arguments:
    gridX_aeqd -- x coordinates for a grid matrix of tracer points (jxi)
    gridY_aeqd -- y coordinates for a grid matrix of tracer points (jxi)
    xy_aeqd    -- x,y coordinates for the input point
    
    source: https://kbkb-wx-python.blogspot.com/2016/08/find-nearest-latitude-and-longitude.html
    '''
    
    x_aeqd, y_aeqd = xy_aeqd

    # Compute x and y differences between all x,y grid 
    # points and the input x,y data point
    absDeltaX = np.abs(gridX_aeqd - x_aeqd)
    absDeltaY = np.abs(gridY_aeqd - y_aeqd)

    # Find local grid delta maximums
    maxDeltas = np.maximum(absDeltaX, absDeltaY)

    # The minimum local maximum is at the nearest grid tracer point
    nearest_idx = np.argmin(maxDeltas)

    # Retrieve the nearest grid tracer point's (j,i) index and return it
    nearest_j, nearest_i = np.unravel_index(nearest_idx, gridX_aeqd.shape)

    return nearest_j, nearest_i


def define_gridCS(fLAT, fLON, ji):
    ''' (float(j,i), float(j,i), tuple) -> (tuple, tuple, tuple, float, float)
    
    Finds the data needed to define and use the local 
    grid coordinate system (CS), i.e.:

            B(0,b)     A

            C(0,0)   D(d,0)
    
    Returns the (lat, lon) coords of B/C/D and the value of b and d.

    Keyword arguments:
    fLAT -- latitudes for a grid matrix of 'points vitesse' (jxi)
    fLON -- longitudes for a grid matrix of 'points vitesse' (jxi)
    ji   -- (j, i) tuple of the matrix index at which the reference 
            grid tracer point is located
    
    '''
    # Unpack the input tracer point indices
    j, i = ji
    
    # Find the lat/lon of the tracer point's associated speed points
    B = (fLAT[j, i-1], fLON[j, i-1])
    C = (fLAT[j-1, i-1], fLON[j-1, i-1])
    D = (fLAT[j-1, i], fLON[j-1, i])  

    # Compute the distance between the origin (C) and B/D (respectively)
    b =  haversine(B, C)
    d =  haversine(D, C)
    
    return B, C, D, b, d


def get_xy_gridCS(gridCS, latlonPt ):
    ''' (tuple, tuple) -> (tuple)
    
    Returns the x,y coordinates of a lat,lon point in the local 
    grid CS defined by the return value of the define_gridCS function.
    
    Keyword arguments:
    gridCS   -- return value of define_gridCS function
    latlonPt -- (lat, lon) tuple  
    '''

    # Unpack the gridCS tuple
    B, C, D, b, d = gridCS

    # Compute the distance between B, C, D and the input point (respectively)
    a1 = haversine(B, latlonPt)
    a2 = haversine(C, latlonPt)
    a3 = haversine(D, latlonPt)
    
    # Compute the x and y coordinates in the local CS
    '''
    Note: see RCM_deformation.pdf (section 3) and xy_gridCoords_derivation.pdf 
    in the docs folder for more information on this calculation.
    '''
    x = ( a2**2 - a3**2 + d**2 ) / (2 * d)
    y = ( a2**2 - a1**2 + b**2 ) / (2 * b)

    return x, y