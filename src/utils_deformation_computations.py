'''
Author: Beatrice Duval (bdu002)

-------------------------------------------------------
Utils - Deformation computations
-------------------------------------------------------

Helper code for the calculation of sea-ice deformations.
'''

def calculate_uv_lists( sx_list, ex_list, sy_list, ey_list, dT):
    ''' (list, list, float) -> float(list, list)

    Function that calculates u and v velocity components of 
    each cell vertices and returns them as a list each. 

    Returns a tuple (u_list, v_list) of the lists of velocity components.

    Keyword arguments:
    sx_list -- list of starting x positions of each cell vertices
    ex_list -- list of ending x positions of each cell vertices
    sy_list -- list of starting x positions of each cell vertices
    ey_list -- list of ending x positions of each cell vertices
    dT      -- time interval
    '''
    
    # Find the number of cell vertices
    n = len(sx_list)

    # Compute the u and v velocity components at the current vertex
    u_list = (ex_list - sx_list) / dT
    v_list = (ey_list - sy_list) / dT

    return u_list, v_list


def calculate_strainRates( u_list, v_list, sx_list, sy_list ):
    ''' (float, list, list) float

    Computes the strain rates (or velocity derivatives).

    Keyword arguments:
    u_list -- list of u component velocities for each cell vertices
    v_list -- list of v component velocities for each cell vertices
    sx_list -- list of starting x positions for each cell vertices
    sy_list -- list of starting y positions for each cell vertices
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

    return dudx, dudy, dvdx, dvdy
