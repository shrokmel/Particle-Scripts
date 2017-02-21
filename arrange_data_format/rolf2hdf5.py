
##############################################################################


""" Turn from the data format of Rolf to hdf5 file format"""

###
# Rolf saves data in standard .dump format of LAMMPS,
# and then converts that into .char and .header files through his read_char.py script,
# which applies a half-homemade compression algorithm to the .dump LAMMPS data.
# Here the aim is to convert these .char and .header files through his read_char.py script,
# to a hdf5 file with a given data structure convention.

##### DATA STRUCTURE

# /positions/x   -> [nsteps, 2, nbeads]

# /info/box/x
# /info/box/y
# /info/dt
# /info/nsteps
# /info/nbeads   -> (total number of beads)
# /info/nsteps   
# /info/nsamp

# /param/kappa
# /param/fp
# /param/nbpf    -> (number of beads per filament)
# /param/bl      -> (bond length)
# /param/sigma
# /param/gamma
# /param/kT

#####

##############################################################################

## load needed libraries and necessary files

import argparse
import numpy as np
import h5py
import os
import codecs
import read_char

##############################################################################

def get_sim_info():
    """ read general simulation information through argument parsing"""
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--density", help="Packing fraction of the system")    
    parser.add_argument("-k", "--kappa", help="Bending rigidity")
    parser.add_argument("-f", "--fp", help="Propulsion force")
    parser.add_argument("-dt", "--timestep", help="Timestep of the simulation")
    parser.add_argument("-ns", "--nsamp", help="Sampling rate of data")
    parser.add_argument("-b", "--bl", help="Bond length of the simulation")    
    parser.add_argument("-s", "--sigma", help="Lennard Jones length")        
    parser.add_argument("-g", "--gamma", help="Friction coefficient")        
    parser.add_argument("-kt", "--kT", help="kT constant")        
    parser.add_argument("-nb", "--nbpf", help="Number of beads per filament")  
    parser.add_argument("-ch", "--charfile", help=".char file path")        
    parser.add_argument("-he", "--headerfile", help=".header file path")      
    args = parser.parse_args()
    
    return args

##############################################################################
    
def read_rolf_data(charfile, headerfile):
    """ read the data from the file Rolf saved and load them into arrays"""
    
    ### open files for reading
    
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    
    ### get general information from header file
    
    nbeads, nsteps = read_char.read_first(hfile)
 
    ### skip the initial snapshots
    
    nskip = 1000
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    
    if nsteps <= 0:
        print 'Error, nsteps < 0'
        exit()

    ### allocate arrays to store the data
    
    nsteps = 100
    pos = np.zeros((nsteps,2,nbeads), dtype=np.float64)
        
    ### loop over all steps and save data
   
    for step in range(nsteps):
        
        ### print stats
        
        print 'step / nsteps ', step, ' / ', nsteps
         
        ### read in the data
        
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
        
        pos[step][0] = x
        pos[step][1] = y
        
    return pos, lx, ly, nsteps, nbeads
    
##############################################################################

def save_data(pos, lx, ly, dt, nsteps, nbeads, nsamp, density, kappa, fp, nbpf, bl, sigma, gamma, kT):
    """ save the data in hdf5 file format"""
    
    ### open the output file
    
    ofilepath = "out.h5"
    ofile = h5py.File(ofilepath, 'w')
    
    ### save the positions
    
    positions = ofile.create_group('positions')
    print pos.shape
    positions.create_dataset('x', (nsteps,2,nbeads), data=pos, dtype=np.float64, compression='gzip')
    
    ### save general simulation info
    
    info = ofile.create_group('info')
    box = info.create_group('box')
    box.create_dataset('x', data=lx)
    box.create_dataset('y', data=ly)
    info.create_dataset('dt', data=dt)
    info.create_dataset('nsteps', data=nsteps)
    info.create_dataset('nbeads', data=nbeads)
    info.create_dataset('nsamp', data=nsamp)
    
    ### save simulation parameters

    param = ofile.create_group('param')
    param.create_dataset('density', data=density)
    param.create_dataset('kappa', data=kappa)
    param.create_dataset('fp', data=fp)
    param.create_dataset('nbpf', data=nbpf)    
    param.create_dataset('bl', data=bl)    
    param.create_dataset('sigma', data=sigma)    
    param.create_dataset('gamma', data=gamma)    
    param.create_dataset('kT', data=kT)    

    ofile.close()

    return
       
##############################################################################
    
def main():
    
    ### get the necessary general info about the simulation by argument parsing
    
    args = get_sim_info()

    ### load the data generated in the simulation 
    
    positions, lx, ly, nsteps, nbeads = read_rolf_data(args.charfile, args.headerfile)
    
    ### save the data 
    
    save_data(positions, lx, ly, args.timestep, nsteps, nbeads, args.nsamp, \
                  args.density, args.kappa, args.fp, args.nbpf, args.bl, \
                      args.sigma, args.gamma, args.kT) 
    
    return
    
##############################################################################    
    
if __name__ == '__main__':
    main()
    
##############################################################################
