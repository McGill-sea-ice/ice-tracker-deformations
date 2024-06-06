'''
Author: Mathieu Plante, Beatrice Duval (bdu002)

----------------------
Python object containing triangulated tracked features information
----------------------

This performs a Delaunay triangulation on the input data file, and contains the results.

    - no. is the triangle number
    - vertice_idx1, vertice_idx2, vertice_idx3 are the triangle vertices' indices in the raw file

'''

# Loading from default packages
import csv
import os
from scipy.spatial import Delaunay
import numpy as np
from math import acos

# Load data and organize is into a triangular arrays
class make_into_tri_arrays:

    def __init__(self, data= None, config=None):
        """
        This function reads and loads data from the ASITS output files.
        This creates an object containing the data fields.

        INPUTS:
        FileName = path to the .dat file with SIM data.
        config   = namelist object, from the python config parser object.

        OBJECT CARACTERISTICS:
        """

        #_________________________________________________________________________________________
        #PERFORM DELAUNAY TRIANGULATION
        self.fault = 0
        startXY = list(zip(data.startX, data.startY))

        try:
            tri = Delaunay(startXY)
        except:
            self.fault = 1
            print("Cant triangulate")
        if self.fault == 1:
            return

        #_________________________________________________________________________________________
        # List triangles and store vertex informations

        self.n   = []
        self.vertice_ids1 = []
        self.vertice_ids2 = []
        self.vertice_ids3 = []
        self.sX1 = []
        self.sX2 = []
        self.sX3 = []
        self.sY1 = []
        self.sY2 = []
        self.sY3 = []
        self.eX1 = []
        self.eX2 = []
        self.eX3 = []
        self.eY1 = []
        self.eY2 = []
        self.eY3 = []

        for n in range(len(tri.simplices)):

            # Find the index of all 3 data points that form the current triangle
            vertice_ids1 = tri.simplices[n][0]
            vertice_ids2 = tri.simplices[n][1]
            vertice_ids3 = tri.simplices[n][2]

            # Find the starting X/Y position of each triangle vertex
            sX1 = data.startX[vertice_ids1]
            sY1 = data.startY[vertice_ids1]
            sX2 = data.startX[vertice_ids2]
            sY2 = data.startY[vertice_ids2]
            sX3 = data.startX[vertice_ids3]
            sY3 = data.startY[vertice_ids3]
            # Find the end X/Y position of each triangle vertex
            eX1 = data.endX[vertice_ids1]
            eY1 = data.endY[vertice_ids1]
            eX2 = data.endX[vertice_ids2]
            eY2 = data.endY[vertice_ids2]
            eX3 = data.endX[vertice_ids3]
            eY3 = data.endY[vertice_ids3]

            # Compute the current triangle's angles
            angles = self.get_tri_angles( (sX1, sY1), (sX2, sY2), (sX3, sY3))

            # Keep the triangle if and only if none of its angles are inferior to 10 degrees (= 0.175 rad)
            if not any(angle < 0.175 for angle in angles):
                self.n.append(n)
                self.vertice_ids1.append(vertice_ids1)
                self.vertice_ids2.append(vertice_ids2)
                self.vertice_ids3.append(vertice_ids3)
                self.sX1.append(sX1)
                self.sX2.append(sX2)
                self.sX3.append(sX3)
                self.sY1.append(sY1)
                self.sY2.append(sY2)
                self.sY3.append(sY3)
                self.eX1.append(eX1)
                self.eX2.append(eX2)
                self.eX3.append(eX3)
                self.eY1.append(eY1)
                self.eY2.append(eY2)
                self.eY3.append(eY3)


        if len(self.vertice_ids1) < 1:
            self.fault = 1

    def get_tri_angles(self, xy1, xy2, xy3):
        ''' (tuple, tuple, tuple) -> tuple

        Returns the triangle's angle at each vertex in radians.

        Keyword arguments: \\
        xy1, xy2, xy3  -- (x, y) cartesian coordinates for the 1st, 2nd and 3rd triangle vertices
        '''

        # Unpack the input coordinates
        x1, y1 = xy1
        x2, y2 = xy2
        x3, y3 = xy3

        # Get the distance between each vertices
        d12 = ( (x2-x1)**2 + (y2-y1)**2 )**0.5
        d13 = ( (x3-x1)**2 + (y3-y1)**2 )**0.5
        d23 = ( (x3-x2)**2 + (y3-y2)**2 )**0.5

        # Compute the angle at each vertex (in radians)
        ang1 = acos( (d13**2 + d12**2 - d23**2) / (2*d13*d12) )
        ang2 = acos( (d12**2 + d23**2 - d13**2) / (2*d23*d12) )
        ang3 = acos( (d23**2 + d13**2 - d12**2) / (2*d13*d23) )

        return ang1, ang2, ang3
