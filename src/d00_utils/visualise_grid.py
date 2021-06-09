'''
Author: Beatrice Duval (bdu002)

------------------------------------------
Visualisation tools for netcdf grid data.
------------------------------------------

'''

import cartopy.crs as ccrs
import matplotlib.pyplot as plt


def show_distCoast(LAT, LON, DIST):
    ''' (float(j,i), float(j,i), float(j,i)) -> none
    
    Produces a plot of netcdf grid points and colors
    them as a functon of their distance from land

    Keyword arguments:
    LAT  -- Tracer point latitudes matrix (jxi)
    LON  -- Tracer point longitdes matrix (jxi)
    DIST -- Distances from land matrix (jxi)
    '''
    
    # Initialize figure
    fig = plt.figure()
    ax = fig.add_subplot(111,
                            projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                                central_latitude=90))

    # Set the map extent in order to see the 
    # entire region of interest
    ax.set_extent((-3000000, 4000000, 8500000, 11900000), ccrs.AzimuthalEquidistant())

    # Plot the grid points and color them as a function of
    # their distance from land
    cb = ax.scatter(LON, LAT, c=DIST, transform=ccrs.Geodetic())

    # Add a colorbar
    plt.colorbar(cb, orientation='vertical',ticklocation='auto')
    
    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()

    # Add a title
    fig.suptitle('Grid Points Colored as a Function of Distance (km) from Coast', fontsize=14, y=1)

    # Show plot
    plt.show()


def show_TraceurVitessePt(LAT, LON, fLAT, fLON, j, i):
    ''' (float(j,i), float(j,i), float(j,i), float(j,i), int, int) -> none
    
    Produces a plot of tracer points and speed points
    that are annotated relatively to a (j,i) tracer point
    
    Keyword arguments:
    LAT  -- Tracer point latitudes matrix (j,i)
    LON  -- Tracer point longitdes matrix (j,i)
    fLAT -- Speed point latitudes matrix (j,i)
    fLON -- Speed point longitudes matrix (j,i)
    j    -- Matrix single tracer point index
    i    -- Matrix single tracer point index
    '''

    # Initialize figure
    fig = plt.figure()
    ax = fig.add_subplot(111,
                            projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                                    central_latitude=90))


    # Create a matplotlib transform from the cartopy coordinate system
    crs = ccrs.Geodetic()
    transform = crs._as_mpl_transform(ax)

    # Plot 'points traceurs' and their neighboring 'points vitesse' (f points) 
    # and annotate them
    for m in [-1, 0, 1]: 
        for n in [-1, 0, 1]:
            if m==0 and n==0:
                # Create a label only once for 'points traceurs' and 'points vitesse'
                ax.scatter(LON[j+n,i+m], LAT[j+n,i+m], transform=ccrs.Geodetic(), color='red', label="Tracer points")
                ax.scatter(fLON[j+n,i+m], fLAT[j+n,i+m], transform=ccrs.Geodetic(), color='blue', label="Speed points") 
            
            else:    
                ax.scatter(fLON[j+n,i+m], fLAT[j+n,i+m], transform=ccrs.Geodetic(), color='blue')
                ax.scatter(LON[j+n,i+m], LAT[j+n,i+m], transform=ccrs.Geodetic(), color='red')


            # Create text variable for annotation
            if m == 0:
                xtxt=""
            if m == 1:
                xtxt="+1"
            if m == -1:
                xtxt="-1"
            if n == 0:
                ytxt=""
            if n == 1:
                ytxt="+1"
            if n == -1:
                ytxt="-1" 
            
            # Indicate the relative (y,x) positionning 
            ax.annotate("(y" + ytxt + ", x" + xtxt + ")", (fLON[j+j,i+i], fLAT[j+j,i+i]), xytext=(1, 1), 
                                                                    xycoords=transform,
                                                                    textcoords='offset points')
            
            ax.annotate("(y" + ytxt + ", x" + xtxt + ")", (LON[j+j,i+i], LAT[j+j,i+i]), xytext=(1, 1), 
                                                                    xycoords=transform,
                                                                    textcoords='offset points')
            
    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()
    plt.legend(loc='lower left', fontsize='x-small')

    # Add a title
    plt.suptitle("Visualization of Tracer Points \nand Neighboring Speed Points", fontsize=14, y=1)

    # Show plot
    plt.show()


def show_grid(LAT, LON, j1, j2, i1, i2):
    ''' (float(j,i), float(j,i), int, int, int, int) -> none
    
    Produces a plot of tracer points from a grid that is 
    cropped y-wise from y1 to y2 and x-wise from x1 to x2

    Keyword arguments:
    LAT  -- Tracer points latitudes matrix (jxi)
    LON  -- Tracer points longitdes matrix (jxi)
    y1   -- Starting column index
    y2   -- Ending column index (not included)
    x1   -- Starting row index
    x2   -- Ending row index (not included)
    '''

    # Crop LAT/LON matrices 
    LAT = LAT[j1:j2, i1:i2] 
    LON = LON[j1:j2, i1:i2]

    # Initialize figure
    fig = plt.figure()
    ax = fig.add_subplot(111,
                            projection=ccrs.AzimuthalEquidistant(central_longitude=0,
                                                                    central_latitude=90))
    # Set the map extent in order to see the 
    # entire region of interest
    ax.set_extent((-3000000, 4000000, 8500000, 11900000), ccrs.AzimuthalEquidistant())

    # Add coastlines and gridlines
    ax.coastlines()
    ax.gridlines()

    # Plot LAT/LON grid points and data points
    ax.scatter(LON, LAT, transform=ccrs.Geodetic(), color="red", marker='.')

    plt.show()