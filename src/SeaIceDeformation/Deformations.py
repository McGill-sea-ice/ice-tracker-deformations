'''
Author: Mathieu Plante, Beatrice Duval (bdu002)

----------------------
Python object calculating the sea ice deformations and store the SIDRR information
from one specific pair of SAR images.
----------------------

'''

# Loading from default packages
import csv
import os
from scipy.spatial import Delaunay
from tqdm import tqdm
import datetime
import numpy as np
from math import sqrt


# Load data and organize is into a triangular arrays
class compute_SIDRR:

    def __init__(self, data = None, tri_array= None, CDFout = None,config=None, idpair = None):
        """
        This function takes the triangulated arrays and compute the
        SIDRR for each triangle.

        INPUTS:

        config   = namelist object, from the python config parser object.

        """

        #====================================================
        #
        # COMPUTE DEFORMATIONS FOR EACH TRIANGLE
        #
        #====================================================
        sigx = data.sigx
        CDFout.sigx = sigx
        dt = data.dt
        self.num_tri = len(tri_array.vertice_ids1)
        num_tri = int(self.num_tri)
        self.nbfg = 0
        self.nbfb = 0
        # Iterate through every triangle in the triangulated csv file
        for n in range(num_tri):

            # Create lists of triangular cell vertices' positions
            sx_list = np.array([tri_array.sX1[n], tri_array.sX2[n], tri_array.sX3[n]])*sigx  # Starting x positions
            sy_list = np.array([tri_array.sY1[n], tri_array.sY2[n], tri_array.sY3[n]])*sigx  # Starting y positions
            ex_list = np.array([tri_array.eX1[n], tri_array.eX2[n], tri_array.eX3[n]])*sigx  # Ending x positions
            ey_list = np.array([tri_array.eY1[n], tri_array.eY2[n], tri_array.eY3[n]])*sigx  # Ending y positions


            # Create a list of velocity components for each triangle vertex
            u_list, v_list = self.calculate_uv_lists( sx_list, ex_list, sy_list, ey_list, dt)

            # Compute the strain rates
            dudx, dudy, dvdx, dvdy, A = self.calculate_strainrates( u_list, v_list, sx_list, sy_list )

            # Compute the divergence rate
            eps_I = dudx + dvdy

            # Compute the maximum shear strain rate
            eps_II = sqrt(  (dudx - dvdy)**2 + (dudy + dvdx)**2  )

            # Compute the vorticity
            vrt = dvdx - dudy

            # Compute the total sea-ice deformation rate
            eps_tot = sqrt( eps_I**2 + eps_II**2 )

            # Compute the propagation errors
            del11, del12, del21, del22, delA = self.calculate_trackerrors( u_list, v_list, sx_list, sy_list,dt)

            del11 = dudx*dudx*delA/A**2.0 + del11*(sigx**2.0)
            del22 = dvdy*dvdy*delA/A**2.0 + del22*(sigx**2.0)
            del12 = dudy*dudy*delA/A**2.0 + del12*(sigx**2.0)
            del21 = dvdx*dvdx*delA/A**2.0 + del21*(sigx**2.0)
            delI = del11 + del22
            delvrt = del12 + del21
            term1 = (delI*(dudx - dvdy)**2.0 + (del12+del21)*(dudy + dvdx)**2.0) #this is eps_II^2 * delII^2
            term2 = eps_I**2.0 * delI
            epstot_deltot = term1 + term2

            if eps_II > 0.0:
                delII = (delI*(((dudx - dvdy)/eps_II)**2.0) + (del12+del21)*(((dudy + dvdx)/eps_II)**2.0))
                deltot = delII*((eps_II/eps_tot)**2.0) + delI*(( eps_I/eps_tot)**2.0)
                sig2noise = eps_tot**2.0 / epstot_deltot**0.5
            else:
                delII = 1000.0
                deltot = 1000.0
                sig2noise = 0.0

            #====================================================
            #
            # Add computed data to variables if area is ok
            #
            #====================================================

            if (A**0.5 < 20000.0):
                self.nbfg += 1
                # Write the vertices' IDs and triangle ID to lists
                CDFout.ids1.append(tri_array.vertice_ids1[n])
                CDFout.ids2.append(tri_array.vertice_ids2[n])
                CDFout.ids3.append(tri_array.vertice_ids3[n])
                if idpair is None:
                    idpair = 0
                CDFout.idpair.append(idpair)

                # Add the satellite to the netcdf lists
                CDFout.sat.append(int(data.sat))

                # Add the divergence and the shear strain rates to the netcdf lists
                CDFout.div.append(eps_I)
                CDFout.shr.append(eps_II)
                CDFout.vrt.append(vrt)

                # Add the strain rates and area to the netcdf lists
                CDFout.A.append(A)
                CDFout.dudx.append(dudx)
                CDFout.dudy.append(dudy)
                CDFout.dvdx.append(dvdx)
                CDFout.dvdy.append(dvdy)

                # Add the strain rates and area to the netcdf lists
                CDFout.delA.append(delA**0.5)
                CDFout.delI.append(delI**0.5)
                CDFout.delII.append(delII**0.5)
                CDFout.delvrt.append(delvrt**0.5)
                CDFout.s2n.append(sig2noise)


                # Add the starting and ending Lat/Lon positions of each triangle vertices
                # to the netcdf list
                CDFout.sLat1.append(np.array(data.sLat)[tri_array.vertice_ids1[n]])
                CDFout.sLat2.append(np.array(data.sLat)[tri_array.vertice_ids2[n]])
                CDFout.sLat3.append(np.array(data.sLat)[tri_array.vertice_ids3[n]])

                CDFout.sLon1.append(np.array(data.sLon)[tri_array.vertice_ids1[n]])
                CDFout.sLon2.append(np.array(data.sLon)[tri_array.vertice_ids2[n]])
                CDFout.sLon3.append(np.array(data.sLon)[tri_array.vertice_ids3[n]])

                CDFout.eLat1.append(np.array(data.eLat)[tri_array.vertice_ids1[n]])
                CDFout.eLat2.append(np.array(data.eLat)[tri_array.vertice_ids2[n]])
                CDFout.eLat3.append(np.array(data.eLat)[tri_array.vertice_ids3[n]])

                CDFout.eLon1.append(np.array(data.eLon)[tri_array.vertice_ids1[n]])
                CDFout.eLon2.append(np.array(data.eLon)[tri_array.vertice_ids2[n]])
                CDFout.eLon3.append(np.array(data.eLon)[tri_array.vertice_ids3[n]])

            else:
                self.nbfb += 1
                self.num_tri = self.num_tri - 1

        # Add the starting and ending times (in seconds since the reference time)
        # to the times list
        CDFout.sTime.extend([data.sTime for i in range(self.num_tri)])
        CDFout.eTime.extend([data.eTime for i in range(self.num_tri)])

        self.num_large = len(tri_array.vertice_ids1) - self.num_tri



    def calculate_uv_lists(self, sx_list, ex_list, sy_list, ey_list, dT):
        '''
        (list, list, float) -> float(list, list)

        Function that calculates u and v velocity components of
        each cell vertices and returns them as a list.

        Returns a tuple (u_list, v_list) of the lists of velocity components.

        Keyword arguments: \\
        sx_list -- list of starting x positions of each cell vertices \\
        ex_list -- list of ending x positions of each cell vertices \\
        sy_list -- list of starting x positions of each cell vertices \\
        ey_list -- list of ending x positions of each cell vertices \\
        dT      -- time interval
        '''

        # Compute the u and v velocity components at the current vertex
        u_list = (ex_list - sx_list) / dT
        v_list = (ey_list - sy_list) / dT

        return u_list, v_list


    def calculate_strainrates(self, u_list, v_list, sx_list, sy_list ):

        '''
        (list, list, list, list) -> tuple(float, float, float, float)
        Computes the strain rates (or velocity derivatives).

        Keyword arguments: \\
        u_list -- list of u component velocities for each cell vertices \\
        v_list -- list of v component velocities for each cell vertices \\
        sx_list -- list of starting x positions for each cell vertices \\
        sy_list -- list of starting y positions for each cell vertices \\
        '''

        # Find the number of cell vertices
        n = len(sx_list)

        #----- Compute the Lagrangian cell area A -----
        # Initialize the cell area to 0
        A = 0.0

        # Perform a summation to compute A (see Bouchat et al. (2020) eqn. 6)
        for i in range( n ):
            A += (1.0/2.0) * ( sx_list[i] * sy_list[((i+1) % n)] - sx_list[((i+1) % n)] *  sy_list[i])


        #----- Compute the strain rates ---------------

        # Initialize the strain rates to 0
        dudx = 0.0
        dudy = 0.0
        dvdx = 0.0
        dvdy = 0.0

        # Perform a summation to compute strain rates (see Bouchat et al. (2020) eqn. 5)
        for i in range( n ):

            dudx += 1.0/(2.0*A)  * ( u_list[((i+1) % n)] + u_list[i] ) * ( sy_list[((i+1) % n)] - sy_list[i] )

            dudy += -1.0/(2.0*A) * ( u_list[((i+1) % n)] + u_list[i] ) * ( sx_list[((i+1) % n)] - sx_list[i] )

            dvdx += 1.0/(2.0*A)  * ( v_list[((i+1) % n)] + v_list[i] ) * ( sy_list[((i+1) % n)] - sy_list[i] )

            dvdy += -1.0/(2.0*A) * ( v_list[((i+1) % n)] + v_list[i] ) * ( sx_list[((i+1) % n)] - sx_list[i] )

        return dudx, dudy, dvdx, dvdy, A


    def calculate_trackerrors(self, u_list, v_list, sx_list, sy_list, Deltat ):

        ''' (float, list, list) float
        Computes the strain rates (or velocity derivatives).

        Keyword arguments:
        u_list -- list of u component velocities for each cell vertices
        v_list -- list of v component velocities for each cell vertices
        sx_list -- list of starting x positions for each cell vertices
        sy_list -- list of starting y positions for each cell vertices
        Deltat  -- Time difference between the 2 images
        '''

        # Find the number of cell vertices
        n = len(sx_list)

        #----- Compute error on the Lagrangian cell area A -----

        A = 0.0
        delA = 0.0
        for i in range( n ):
            A += (1.0/2.0) * ( sx_list[i] * sy_list[((i+1) % n)] - sx_list[((i+1) % n)] *  sy_list[i])
            delA += ((sx_list[i] - sx_list[((i+2) % n)])**2.0 + (sy_list[((i+2) % n)] - sy_list[i])**2.0)/4.0

        #----- Compute the strain rate errors ---------------

        # Initialize the strain rates to 0
        trackerror11 = 0.0
        trackerror12 = 0.0
        trackerror21 = 0.0
        trackerror22 = 0.0

        # Perform a summation to compute strain rates (see Bouchat et al. (2020) eqn. 5)
        for i in range( n ):

            trackerror11 += (sy_list[i] - sy_list[((i+2) % n)])**2.0/(2.0*A*A*Deltat*Deltat) + (u_list[i]-u_list[((i+2) % n)])**2.0/(4.0*A*A)
            trackerror22 += (sx_list[i] - sx_list[((i+2) % n)])**2.0/(2.0*A*A*Deltat*Deltat) + (v_list[i]-v_list[((i+2) % n)])**2.0/(4.0*A*A)
            trackerror12 += (sx_list[i] - sx_list[((i+2) % n)])**2.0/(2.0*A*A*Deltat*Deltat) + (u_list[i]-u_list[((i+2) % n)])**2.0/(4.0*A*A)
            trackerror21 += (sy_list[i] - sy_list[((i+2) % n)])**2.0/(2.0*A*A*Deltat*Deltat) + (v_list[i]-v_list[((i+2) % n)])**2.0/(4.0*A*A)

        return trackerror11, trackerror12, trackerror21, trackerror22, delA


