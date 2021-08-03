'''
Author: Beatrice Duval (bdu002)

----------------------------------------------------------------
Utils - Load modified RIOPS grid
----------------------------------------------------------------

Code that provides a function that loads the project RIOPS grid.

'''

import os

from netCDF4 import Dataset

def load_grid():

    # Find absolute path in which the RIOPS grid is stored
    projPath = os.path.dirname(os.path.realpath(__file__))
    gridPath = projPath + '/../data/00_grid/cropped_grid.nc'

    # Load RIOPS grid dataset
    ds = Dataset(gridPath)

    # Retrieve grid matrices
    LAT    = ds['nav_lat'][:,:]    # Tracer points (latitudes)
    LON    = ds['nav_lon'][:,:]    # Tracer points (longitudes)

    Y_aeqd = ds['nav_y'][:,:]      # Tracer points (y - aeqd coord. syst.)
    X_aeqd = ds['nav_x'][:,:]      # Tracer points (x - aeqd coord. syst.)

    fLAT    = ds['gphif'][:,:]     # Speed points (latitudes)
    fLON    = ds['glamf'][:,:]     # Speed points (longitudes)

    MASK    = ds['mask'][:,:]      # Ocean-land mask (ocean = 1, land = 0)
    DIST    = ds['dist'][:,:]      # Distances from land matrix

    return {'LAT': LAT, 'LON': LON, 'Y_aeqd': Y_aeqd, 'X_aeqd': X_aeqd, 'fLAT': fLAT, 'fLON': fLON, 'MASK': MASK, 'DIST': DIST}