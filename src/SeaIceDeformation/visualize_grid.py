import utils_visualise_grid
import utils_load_grid

grid_dict = utils_load_grid.load_grid('/home/leki/ice-tracker-deformations/data/00_grid/cropped_grid.nc')

LAT = grid_dict["LAT"]
LON = grid_dict['LON']

utils_visualise_grid.show_grid(LAT, LON, 1, 600, 1, 800)