'''
Author: Beatrice Duval (bdu002)

---------------------------------------------------------------------------
Code that plots the distances from land matrix from the cropped RIOPS grid.
---------------------------------------------------------------------------

'''

from src.d00_utils.visualise_grid import show_distCoast
from src.d01_data.load00_grid import LAT, LON, DIST

show_distCoast(LAT, LON, DIST)