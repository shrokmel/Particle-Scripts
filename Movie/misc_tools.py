
##############################################################################

import numpy as np
import math

##############################################################################

def nearest_neighbor(x1, x2, lx):
    """ compute vector to nearest neighbor"""
    
    dx1 = x1 - x2
    dx2 = x1 - x2 + lx
    dx3 = x1 - x2 - lx
    if dx1**2 < dx2**2 and dx1**2 < dx3**2:
        return dx1
    if dx2**2 < dx3**2:
        return dx2
    return dx3

##############################################################################

def compute_orientation(x, y, lx, ly, npol):
    """ compute orientation of all beads from bond vectores"""
    
    # number of molecules
    
    natoms = len(x)
    nmol = natoms/npol
    
    # allocate aray for results
    
    phi = np.zeros((natoms))
    
    # loop over all polymers
    
    k = 0
    
    for i in range(nmol):
        
        for j in range(npol):
            
            if j == 0:
                x1 = x[k]
                y1 = y[k]
                x2 = x[k+1]
                y2 = y[k+1]
            elif j == npol-1:
                x1 = x[k-1]
                y1 = y[k-1]
                x2 = x[k]
                y2 = y[k]
            else:
                x1 = x[k-1]
                y1 = y[k-1]
                x2 = x[k+1]
                y2 = y[k+1]

            # compute nearest neighbor
            
            dx = nearest_neighbor(x1, x2, lx)
            dy = nearest_neighbor(y1, y2, ly)
            
            # compute angle using atan2
            
            pi = math.atan2(dy, dx)
            phi[k] = pi

            # increment k
            
            k = k + 1
            
    return phi

##############################################################################

def compute_velocity(x1, y1, x3, y3, lx, ly, tstep1, tstep3, natoms):
    """ compute the velocity of each bead, consider pbc"""
    
    # define target arrays
    
    vx = np.zeros(natoms)
    vy = np.zeros(natoms)
    
    # loop over all beads
    
    for i in range(natoms):
        
        xi = x1[i]
        yi = y1[i]
        xj = x3[i]
        yj = y3[i]
        dx = nearest_neighbor(xj,xi,lx)
        dy = nearest_neighbor(yj,yi,ly)
        vx[i] = dx/(tstep3-tstep1)
        vy[i] = dy/(tstep3-tstep1)
        
    return vx, vy

##############################################################################