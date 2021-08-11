'''
Author: Beatrice Duval (bdu002)

----------------------------------------------------------------
Utils - Load modified RIOPS grid
----------------------------------------------------------------

Code that provides a function that loads the project RIOPS grid.

'''


from netCDF4 import Dataset


def load_grid(cropped_grid):
    ''' (str) -> dict[str, array]

    Function that loads the project RIOPS grid.

    Keyword arguments:
    cropped_grid -- absolute path of the netcdf file that stores the cropped RIOPS grid
    '''

    # Load RIOPS grid dataset
    ds = Dataset(cropped_grid)

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