#!/usr/users/iff_th2/isele/Applications/Anaconda/bin/python2.7

import numpy as np
import codecs

################################################################
# collection of functions to read in the data of a char file
################################################################

def read_first(hfile):
    line = hfile.readline()
    line = line.split()
    nsteps = int(line[3])
    natoms = int(line[4])
    return natoms, nsteps

#################################################################

def read_snapshot(hfile, ifile):
    ### read in the information from the header
    # timestep
    hfile.readline()
    line = hfile.readline()
    line = line.split()
    tstep = int(line[0])
    # natoms
    hfile.readline()
    line = hfile.readline()
    line = line.split()
    natoms = int(line[0])
    # box dimensions
    hfile.readline()
    line = hfile.readline()
    line = line.split()
    xlo = float(line[0])
    xhi = float(line[1])
    lx = xhi - xlo
    line = hfile.readline()
    line = line.split()
    ylo = float(line[0])
    yhi = float(line[1])
    ly = yhi - ylo
    hfile.readline()
    hfile.readline()

    ### allocate memory
    xs = np.zeros((natoms), dtype = np.float64)
    ys = np.zeros((natoms), dtype = np.float64)

    ### read in the entire snapshot
    b = ifile.reader.read(natoms*4,natoms*4)
    k = 0
    for i in range(natoms):
        b1 = b[k]
        b2 = b[k+1]
        b3 = b[k+2]
        b4 = b[k+3]
        k = k + 4

        x1 = ord(b1)
        x2 = ord(b2)
        y1 = ord(b3)
        y2 = ord(b4)

        x = float(x1*256 + x2)/256**2
        y = float(y1*256 + y2)/256**2

        xs[i] = x
        ys[i] = y

    return xs,ys,lx,ly,tstep,natoms

################################################################
  
def skip_snapshots(hfile, ifile, nskip):
    """ skip some snapshots"""
    if nskip < 1:
        return
    # get number of atoms from the first header
    for i in range(3):
        hfile.readline()
    line = hfile.readline()
    line = line.split()
    natoms = int(line[0])
    for i in range(5):
        hfile.readline()
    # skip the remaining header lines
    for i in range(9*(nskip-1)):
        hfile.readline()
    # skip the body
    for i in range(nskip):
        ifile.reader.read(size = natoms*4, chars = natoms*4)
    return
