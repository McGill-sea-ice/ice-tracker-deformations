'''
Author: Mathieu Plante

----------------------
Python object opening and saving the netcdf outputs
----------------------

'''

# Loading from default packages
import os
import numpy as np
from netCDF4 import Dataset


# Load data and organize is into a triangular arrays
class CDFoutput:

    def __init__(self):


        #====================================================
        #
        # Create netcdf dataset variables
        #
        #====================================================

        self.sTime = []
        self.eTime =  []

        self.sat =  []

        self.sLat1 =  []
        self.sLat2 = []
        self.sLat3 =  []
        self.sLon1 =  []
        self.sLon2 =  []
        self.sLon3 =  []

        self.eLat1 =  []
        self.eLat2 =  []
        self.eLat3 = []
        self.eLon1 =  []
        self.eLon2 =  []
        self.eLon3 =  []

        self.div = []
        self.shr = []
        self.vrt = []

        self.ids1 = []
        self.ids2 = []
        self.ids3 = []
        self.idpair = []

        self.A = []

        self.dudx = []
        self.dudy = []
        self.dvdx = []
        self.dvdy = []

        self.delA = []
        self.delI = []
        self.delII = []
        self.delvrt = []
        self.s2n = []

    def append_data_from_pair(self, SIDRRdata= None):


        #====================================================
        #
        # Create netcdf dataset variables
        #
        #====================================================

        self.sTime = self.sTime.append(SIDRRdata.sTime)
        self.eTime = self.eTime.append(SIDRRdata.eTime)

        self.sat   = self.sat.append(SIDRRdata.sat)

        self.sLat1 = self.sLat1.append(SIDRRdata.sLat1)
        self.sLat2 = self.sLat2.append(SIDRRdata.sLat2)
        self.sLat3 = self.sLat3.append(SIDRRdata.sLat3)
        self.sLon1 = self.sLon1.append(SIDRRdata.sLon1)
        self.sLon2 = self.sLon2.append(SIDRRdata.sLon2)
        self.sLon3 = self.sLon3.append(SIDRRdata.sLon3)

        self.eLat1 = self.eLat1.append(SIDRRdata.eLat1)
        self.eLat2 = self.eLat2.append(SIDRRdata.eLat2)
        self.eLat3 = self.eLat3.append(SIDRRdata.eLat3)
        self.eLon1 = self.eLon1.append(SIDRRdata.eLon1)
        self.eLon2 = self.eLon2.append(SIDRRdata.eLon2)
        self.eLon3 = self.eLon3.append(SIDRRdata.eLon3)

        self.div   = self.div.append(SIDRRdata.div)
        self.shr   = self.shr.append(SIDRRdata.shr)
        self.vrt   = self.vrt.append(SIDRRdata.vrt)

        self.ids1  = self.ids1.append(SIDRRdata.ids1)
        self.ids2  = self.ids2.append(SIDRRdata.ids2)
        self.ids3  = self.ids3.append(SIDRRdata.ids3)
        self.idpair = self.idpair.append(SIDRRdata.idpair)

        self.A = self.A.append(SIDRRdata.A)

        self.dudx = self.dudx.append(SIDRRdata.dudx)
        self.dudy = self.dudy.append(SIDRRdata.dudy)
        self.dvdx = self.dvdx.append(SIDRRdata.dvdx)
        self.dvdy = self.dvdy.append(SIDRRdata.dvdy)

        self.delA = self.delA.append(SIDRRdata.delA)
        self.delI = self.delI.append(SIDRRdata.delI)
        self.delII = self.delII.append(SIDRRdata.delII)
        self.delvrt = self.delvrt.append(SIDRRdata.delvrt)
        self.s2n = self.s2n.append(SIDRRdata.s2n)

    def write_CDF(self, config = None,DateString = None):

        Date_options = config['Date_options']
        # Find absolute path in which the output netcdf file is to be stored
        output_path = config['IO']['output_folder']
        # Create a directory to store the output netcdf file if it does not exist already
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        output_filename = 'SIDRR_' + DateString + '.nc'
        # Get the directory in which the calculated .csv file is to be stored
        nc_output_path = output_path + '/' + output_filename

        # Create an output netcdf file and dataset
        output_ds = Dataset(nc_output_path, 'w', format = 'NETCDF4')

        #-------------------------------------------------
        # Add metadata to the netcdf
        #-------------------------------------------------
        Metadata = config['Metadata']
        maxDeltat = int(Date_options['tolerance']) + int(Date_options['timestep'])
        if Metadata['icetracker'] == 'RCMS1':
            SARsource = "RCM and S1"
        else:
            SARsource = Metadata['icetracker']

        Metadata = config['Metadata']
        output_ds.SatelliteSource = SARsource
        output_ds.referenceTime = DateString + ' 00:00:00'
        output_ds.trackingError = str(self.sigx) + ' m'
        output_ds.Max_Deltat = '%s hours' % maxDeltat


        #-------------------------------------------------
        # Create dimension variable
        #-------------------------------------------------

        n = output_ds.createDimension('n', len(self.sTime))

        #-------------------------------------------------
        # Create netcdf variables
        #-------------------------------------------------

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



        #-------------------------------------------------
        # Specify units for each variable
        #-------------------------------------------------

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

        #-------------------------------------------------
        # Populate the variables with the stored values
        #-------------------------------------------------

        start_time[:] = self.sTime
        end_time[:]   = self.eTime

        satellite[:]  = self.sat

        start_lat1[:] = self.sLat1
        start_lat2[:] = self.sLat2
        start_lat3[:] = self.sLat3
        start_lon1[:] = self.sLon1
        start_lon2[:] = self.sLon2
        start_lon3[:] = self.sLon3

        end_lat1[:]   = self.eLat1
        end_lat2[:]   = self.eLat2
        end_lat3[:]   = self.eLat3
        end_lon1[:]   = self.eLon1
        end_lon2[:]   = self.eLon2
        end_lon3[:]   = self.eLon3

        d[:]          = self.div
        s[:]          = self.shr
        v[:]          = self.vrt

        id1[:]        = self.ids1
        id2[:]        = self.ids2
        id3[:]        = self.ids3
        id_pair[:]    = self.idpair

        Aa[:]         = self.A

        dux[:]        = self.dudx
        duy[:]        = self.dudy
        dvx[:]        = self.dvdx
        dvy[:]        = self.dvdy

        sA[:]         = self.delA
        sI[:]         = self.delI
        sII[:]        = self.delII
        svrt[:]       = self.delvrt
        sig2n[:]      = self.s2n

