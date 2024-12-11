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

        self.sigx = 200

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

        self.sLat = []
        self.sLon = []
        self.eLat = []
        self.eLon = []
        self.tri_inv = []
        self.pts_idpair =[]

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
        self.sTime.extend(SIDRRdata.sTime)
        self.eTime.extend(SIDRRdata.eTime)

        self.sat.extend(SIDRRdata.sat)

        self.sLat1.extend(SIDRRdata.sLat1)
        self.sLat2.extend(SIDRRdata.sLat2)
        self.sLat3.extend(SIDRRdata.sLat3)
        self.sLon1.extend(SIDRRdata.sLon1)
        self.sLon2.extend(SIDRRdata.sLon2)
        self.sLon3.extend(SIDRRdata.sLon3)

        self.eLat1.extend(SIDRRdata.eLat1)
        self.eLat2.extend(SIDRRdata.eLat2)
        self.eLat3.extend(SIDRRdata.eLat3)
        self.eLon1.extend(SIDRRdata.eLon1)
        self.eLon2.extend(SIDRRdata.eLon2)
        self.eLon3.extend(SIDRRdata.eLon3)

        #This is Anton Korosov proposed modifications
        ids = np.hstack([SIDRRdata.ids1, SIDRRdata.ids2, SIDRRdata.ids3])
        sLat = np.hstack([SIDRRdata.sLat1, SIDRRdata.sLat2, SIDRRdata.sLat3])
        sLon = np.hstack([SIDRRdata.sLon1, SIDRRdata.sLon2, SIDRRdata.sLon3])
        eLat = np.hstack([SIDRRdata.eLat1, SIDRRdata.eLat2, SIDRRdata.eLat3])
        eLon = np.hstack([SIDRRdata.eLon1, SIDRRdata.eLon2, SIDRRdata.eLon3])
        _, unq_idx, unq_inv = np.unique(ids, return_index=True, return_inverse=True)


        start_lat = sLat[unq_idx]
        start_lon = sLon[unq_idx]
        end_lat = eLat[unq_idx]
        end_lon = eLon[unq_idx]

        pts_idpair = start_lat.copy()*0.0 + np.unique(SIDRRdata.idpair)

        tri = unq_inv.reshape(3, -1).T



        self.sLat.extend(start_lat)
        self.sLon.extend(start_lon)
        self.eLat.extend(end_lat)
        self.eLon.extend(end_lon)
        self.pts_idpair.extend(pts_idpair)

        self.div.extend(SIDRRdata.div)
        self.shr.extend(SIDRRdata.shr)
        self.vrt.extend(SIDRRdata.vrt)

        self.ids1.extend(tri[:,0])
        self.ids2.extend(tri[:,1])
        self.ids3.extend(tri[:,2])

        self.idpair.extend(SIDRRdata.idpair)
        self.A.extend(SIDRRdata.A)

        self.dudx.extend(SIDRRdata.dudx)
        self.dudy.extend(SIDRRdata.dudy)
        self.dvdx.extend(SIDRRdata.dvdx)
        self.dvdy.extend(SIDRRdata.dvdy)

        self.delA.extend(SIDRRdata.delA)
        self.delI.extend(SIDRRdata.delI)
        self.delII.extend(SIDRRdata.delII)
        self.delvrt.extend(SIDRRdata.delvrt)
        self.s2n.extend(SIDRRdata.s2n)

        self.sigx = SIDRRdata.sigx

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
        n = output_ds.createDimension('npts', len(self.sLat))

        #-------------------------------------------------
        # Create netcdf variables
        #-------------------------------------------------

        start_time = output_ds.createVariable('start_time', 'u4', 'n') # Start and end times
        end_time   = output_ds.createVariable('end_time', 'u4', 'n')

        satellite  = output_ds.createVariable('satellite', 'u4', 'n') # Satellite

        start_lat = output_ds.createVariable('start_lat', 'f8', 'npts') # tracked features start locations
        start_lon = output_ds.createVariable('start_lon', 'f8', 'npts')
        end_lat = output_ds.createVariable('end_lat', 'f8', 'npts') # tracked features end locations
        end_lon = output_ds.createVariable('end_lon', 'f8', 'npts')
        pts_idpair = output_ds.createVariable('pts_idpair','u4','npts')

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


        start_lat.units = 'degrees North'
        start_lon.units = 'degrees East'
        end_lat.units   = 'degrees North'
        end_lon.units   = 'degrees East'

        id1.units        = 'ID vertex 1'
        id2.units        = 'ID vertex 2'
        id3.units        = 'ID vertex 3'
        id_pair.units     = 'ID SAR pair'
        pts_idpair.units     = 'ID SAR pair'

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

        start_lat[:] = self.sLat
        start_lon[:] = self.sLon
        end_lat[:]   = self.eLat
        end_lon[:]   = self.eLon

        d[:]          = self.div
        s[:]          = self.shr
        v[:]          = self.vrt

        id1[:]        = self.ids1
        id2[:]        = self.ids2
        id3[:]        = self.ids3
        id_pair[:]    = self.idpair
        pts_idpair[:]    = self.pts_idpair

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

