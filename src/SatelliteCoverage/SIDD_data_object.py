#================================================================
# DATA OBJECT FOR SIDD: Each data class defines an 
# object (e.g. data) that include, amound other variables,
# the divergence (data.eI) and shear (data.eII) on which the 
# the SIDDs will be computed.
#
# - class Eulerian Data: Creates a SIDD data object out of  
#   deformation model outputs, and recomputing the velocity
#   gradients either from Finite Difference on by Line Integral
#   (as specified by the user).
#
#   
# Author: Mathieu Plante, June 2022
#================================================================

import os
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, date, time, timedelta
import subprocess
from netCDF4 import Dataset 
import numpy.ma as ma
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature
import shapely.geometry as sgeom
import cmocean as cmo
from math import floor, pi


import pyproj
import matplotlib.path as mpath


##-----------------------------------------------------------
class SIDD_Data:

# Class specific to the use of gridded and/or masked eI and
# eII variables from model outputs (e.g. CICE, RIOPS).
#
# Optional Grid object is used to print visuals of the field
# variables
#
# Optional Mask object to mask low-concentration or
# close-to-coastline data.
#
#------------------------------------------------------------

    def __init__(self,
                 InputData = None,
                 meta= None,
                 Grid = None,
                 mask = None,
                 Method = None,
                 time = None):

        #if meta.VarNameDiv != None:
        self.eI = InputData.div.copy()
            #self.eII = self.eII*meta.DeformationScaleFactor
            #Mask missing values
            #self.eI = ma.array(self.eI,mask = self.eI >400.0)
        #if meta.VarNameShr != None:
        self.eII = InputData.shr.copy()
            #self.eII = self.eII*meta.DeformationScaleFactor
            #Mask missing values
            #self.eII = ma.array(self.eII,mask = self.eII >400.0)

        self.SIC = None
        #self.u = InputData.u.copy()
        #self.v = InputData.v.copy()
        self.h = None

        self.dudx = InputData.dudx.copy()
        self.dudy = InputData.dudy.copy()
        self.dvdx = InputData.dvdx.copy()
        self.dvdy = InputData.dvdy.copy()

        #self.del11 = GradData.del11.copy()
        #self.del22 = GradData.del22.copy()
        #self.del12 = GradData.del12.copy()
        #self.del21 = GradData.del21.copy()
        self.A = InputData.A.copy()

        self.del11 = None
        self.del22 = None
        self.del12 = None
        self.del21 = None

        self.lon = InputData.start_lon1.copy()
        self.lat = InputData.start_lat1.copy()

    def aggregate(self, data_n = None, n = None):

        if data_n == None:
            self.eI = self.eI.copy()/n
            self.eII = self.eII.copy()/n
            self.dudx = self.dudx.copy()/n
            self.dudy = self.dudy.copy()/n
            self.dvdx = self.dvdx.copy()/n
            self.dvdy = self.dvdy.copy()/n
        
            self.del11 = self.del11.copy()/n 
            self.del22 = self.del22.copy()/n
            self.del12 = self.del12.copy()/n
            self.del21 = self.del21.copy()/n 
        
            self.A = self.A.copy()/n
            self.u = self.u.copy()/n
            self.v = self.v.copy()/n
        else:
            self.eI = self.eI + data_n.eI/n
            self.eII = self.eII+ data_n.eII/n
	
            self.dudx = self.dudx+ data_n.dudx/n
            self.dudy = self.dudy+ data_n.dudy/n
            self.dvdx = self.dvdx+ data_n.dvdx/n
            self.dvdy = self.dvdy+ data_n.dvdy/n
        
            self.del11 = self.del11+ data_n.del11/n 
            self.del22 = self.del22+ data_n.del22/n
            self.del12 = self.del12+ data_n.del12/n
            self.del21 = self.del21+ data_n.del21/n 
        
            self.A = self.A+ data_n.A/n
            self.u = self.u+ data_n.u/n
            self.v = self.v+ data_n.v/n
            

    def show_deformations(self,
                          Grid=None,
                          time=None,
                          output_folder=None,
                          tracker = None):
    # Prints a map of the sea ice deformation data
    
    #Note: as of oct. 2021, the cartopy coastline features are shifted when used with a pcolormesh. We thus use a land mask with black edgecolor. 
    
        cmap = cmo.cm.balance
        cmap.set_bad('grey')
        
       
        Figure1 = plt.figure(figsize=[8.0,4.0])
        
        
        ax =  Figure1.add_axes([0.05, 0.2, 0.4, 0.7],projection=Grid.proj)
        
        im = ax.pcolormesh(Grid.lon,Grid.lat,self.eI,
                           transform=ccrs.PlateCarree(),
                           vmin=-0.1,vmax=0.1,
                           cmap=cmap)
        ax.set_extent(Grid.extent)

        ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')

        ax.text(0.5, -0.15, 
                'divergence (day$^{-1}$)', 
                va='bottom', ha='center',
                rotation='horizontal', 
                rotation_mode='anchor',
                transform=ax.transAxes,fontsize=14)
        plt.colorbar(im)
  
        ax2 =  Figure1.add_axes([0.55, 0.2, 0.4, 0.7],projection=Grid.proj)
        im = ax2.pcolormesh(Grid.lon,Grid.lat,self.eII,
                            transform=ccrs.PlateCarree(),
                            vmin=-0.1,vmax=0.1,cmap=cmap)
        ax2.set_extent(Grid.extent)
        ax2.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
        ax2.text(0.5, -0.15, 
                'shear (day$^{-1}$)',
                va='bottom', ha='center',
                rotation='horizontal', 
                rotation_mode='anchor',
                transform=ax2.transAxes,fontsize=14)
  
        plt.colorbar(im)
  
        Figure1.savefig("%sdeformations_%s.png" % (output_folder,time.ThisTimeFile))
        plt.close(Figure1)
          
    def show_SIC(self,
                 Grid=None,
                 time=None,
                 output_folder=None):
    # Prints a map of the sea ice deformation data
    
    #Note: as of oct. 2021, the cartopy coastline features are shifted when used with a pcolormesh. We thus use a land mask with black edgecolor. 
    
        cmap = cmo.cm.ice
        cmap.set_bad('grey')
        
       
        Figure1 = plt.figure(figsize=[8.0,4.0])
        
        
        ax =  Figure1.add_axes([0.05, 0.2, 0.4, 0.7],projection=Grid.proj)
        
        im = ax.pcolormesh(Grid.lon,Grid.lat,self.h,
                           transform=ccrs.PlateCarree(),
                           vmin=0.0,vmax=3.0,
                           cmap=cmap)
        ax.set_extent(Grid.extent)

        ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')

        ax.text(0.5, -0.15, 
                'ice thickness (m)', 
                va='bottom', ha='center',
                rotation='horizontal', 
                rotation_mode='anchor',
                transform=ax.transAxes,fontsize=14)
        plt.colorbar(im)
        plt.title(time.ThisTimeFile)
  
        ax2 =  Figure1.add_axes([0.55, 0.2, 0.4, 0.7],projection=Grid.proj)
        im = ax2.pcolormesh(Grid.lon,Grid.lat,self.SIC,
                            transform=ccrs.PlateCarree(),
                            vmin=0.0,vmax=1.0,cmap=cmap)
        ax2.set_extent(Grid.extent)
        ax2.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
        ax2.text(0.5, -0.15, 
                'SIC (%)',
                va='bottom', ha='center',
                rotation='horizontal', 
                rotation_mode='anchor',
                transform=ax2.transAxes,fontsize=14)
  
        plt.colorbar(im)
        plt.title(time.ThisTimeFile)
        Figure1.savefig("%sSIC_%s.png" % (output_folder,time.ThisTimeFile))
        plt.close(Figure1)   
               
      
    
