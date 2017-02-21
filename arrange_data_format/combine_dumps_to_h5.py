
""" Combine multiple dump files from different restart instances into a single hdf5 file"""

### example command line arguments: 
###    -fl=/homea/ias2/duman/MASOUD/motor/ -t -k=1 -m=0
###         -dt=0.001 -ns=100000 -b=0.5 -s=1.0 -nf=6000

##############################################################################

import argparse
import numpy as np
import os
import h5py

##############################################################################

def read_contextual_info():
    """ read the contextual information provided by the user"""
    
    ### get the data folder and the last timestep info
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", help="Folder containing data")
    parser.add_argument("-t", "--tstep", nargs="?", const="10000000", \
                            type=int, help="The last time step that is being searched for")
    parser.add_argument("-d", "--density", type=str, help="Density")  
    parser.add_argument("-k", "--kappa", type=str, help="Bending rigidity")  
    parser.add_argument("-m", "--km", type=str, help="Motor spring strength")
    parser.add_argument("-p", "--pa", type=str, help="Probability of motors in case of antialignment")      
    parser.add_argument("-dt", "--timestep", type=float, help="Timestep of the simulation")
    parser.add_argument("-ns", "--nsamp", type=int, help="Sampling rate of data")
    parser.add_argument("-b", "--bl", type=float, help="Bond length of the simulation")    
    parser.add_argument("-s", "--sigma", type=float, help="Lennard Jones length")     
    parser.add_argument("-nf", "--nfils", type=int, help="Number of filaments")
    args = parser.parse_args()
    
    ### generate folder path

    path = 'density_' + args.density + '/kappa_' + args.kappa + '/km_' + args.km + '/panti_' + args.pa + '/' 
    path += 'out1.dump'
    fpath = args.folder + path
    assert os.path.exists(fpath), "\nOUT1.DUMP DOES NOT EXIST FOR: " + fpath 

    print fpath

    ### determine the total number of beads and the box size
    
    fl = open(fpath, 'r')
    fl.readline()
    fl.readline()
    fl.readline()
    line = fl.readline()
    line = line.split()
    nbeads = int(line[0])
    
    fl.readline()
    line = fl.readline()
    line = line.split()
    lx = float(line[1])
    line = fl.readline()
    line = line.split()
    ly = float(line[1])

    fl.close()
    
    ### total number of steps
    
    nsteps = args.tstep/args.nsamp
    nsteps -= 10
    nsteps += 1
            
    return args.folder, nbeads, nsteps, lx, ly, args

##############################################################################
    
def get_number_of_snaps(f, nbeads):
    """ determine the number of snapshots in the file"""
    
    os.system('wc -l ' + f + ' > tmp.txt')
    ifile = open('tmp.txt')
    line = ifile.readline()
    line = line.split()
    nlines = int(line[0])
    nsnaps = nlines/(nbeads+9)
    ifile.close()
    os.system('rm tmp.txt')
    
    return nsnaps

##############################################################################
    
def read_pos(fl, x, nbeads, nsnaps, lx, ly, checked, tstep_cnt, T):
    """ read the position data from the file and return the last timestep at the end"""

    ### read the positions unique per each tstep in a single dump file
    
    already_checked = False
    for snap in range(nsnaps):
        
        ### read the headers to check the uniqueness of current tstep
        
        fl.readline()
        line = fl.readline()
        line = line.split()
        tstep = int(line[0])
        
        ### finish if the last tstep is exceeded already
        
        if tstep > T:
            return tstep_cnt, tstep
        
        ### make sure the current tstep is unique
        
        if tstep not in checked:
            checked.append(tstep)
            tstep_cnt += 1
            already_checked = False
        else:
            print "ALREADY CHECKED " + str(tstep)             
            already_checked = True
    
        ### read the remaning part of the header part
        
        for j in range(7):
            fl.readline()
         
        ### read the positions per bead if the tstep is unique
        
        for j in range(nbeads):
            line = fl.readline()
            if already_checked:
                continue
            line = line.split()
            bid = int(line[0]) - 1
            x[tstep_cnt, 0, bid] = float(line[3])
            x[tstep_cnt, 1, bid] = float(line[4])

        if tstep_cnt%1000 == 0:
            print tstep_cnt, tstep

    return tstep_cnt, tstep

##############################################################################
    
def read_pos_from_dump_files(folder, nbeads, nsteps, T, lx, ly, args):
    """ read the position data of each dump file until the last tstep is reached"""

    ### generate file path and the total number of snapshots in the file
    
    current_file_number = 0
    x = np.zeros((nsteps, 2, nbeads), dtype=np.float32)
    tstep_cnt = -1
    checked = []
    tstep = 0
    
    while tstep < T:
        current_file_number += 1
        path = 'density_' + args.density + '/kappa_' + args.kappa + '/km_' + args.km + '/panti_' + args.pa + '/out' + \
            str(current_file_number) + '.dump'
        fpath = folder + path        
        assert os.path.exists(fpath), "out dump file does NOT exist for: " + fpath
        fl = open(fpath, 'r')
        nsnaps = get_number_of_snaps(fpath, nbeads)
        
        ### read the positions unique per each tstep in a single dump file
        
        tstep_cnt, tstep = read_pos(fl, x, nbeads, nsnaps, lx, ly, checked, tstep_cnt, T)
        fl.close()

    return x
    
##############################################################################
    
def write_h5_file(folder, x, nbeads, nsteps, lx, ly, args):
    """ write data to hdf5 file"""
    
    ### file path
    
    path = 'density_' + args.density + '/kappa_' + args.kappa + '/km_' + args.km + '/panti_' + args.pa + '/' 
    fpath = folder + path + 'out.h5'
    fl = h5py.File(fpath, 'w')
    
    ### positions of beads
    
    pos = fl.create_group('positions')
    pos.create_dataset('x', (nsteps, 2, nbeads), data=x, dtype=np.float32, compression='gzip') 
    
    ### simulation information
    
    info = fl.create_group('info')
    box = info.create_group('box')
    box.create_dataset('x', data=lx, dtype=np.float32)
    box.create_dataset('y', data=ly, dtype=np.float32)
    info.create_dataset('dt', data=args.timestep, dtype=np.float32)
    info.create_dataset('nsteps', data=nsteps, dtype=np.int32)
    info.create_dataset('nbeads', data=nbeads, dtype=np.int32)
    info.create_dataset('nsamp', data=args.nsamp, dtype=np.int32)
    info.create_dataset('nfils', data=args.nfils, dtype=np.int32)
    nbpf = 17
    info.create_dataset('nbpf', data=nbpf, dtype=np.int32)
    
    ### simulation parameters
    
    param = fl.create_group('param')
    param.create_dataset('density', data=float(args.density), dtype=np.float32)
    param.create_dataset('kappa', data=float(args.kappa), dtype=np.float32)
    param.create_dataset('km', data=float(args.km), dtype=np.float32)
    param.create_dataset('pa', data=float(args.pa), dtype=np.float32)    
    param.create_dataset('bl', data=args.bl, dtype=np.float32)
    param.create_dataset('sigma', data=args.sigma, dtype=np.float32)
    
    fl.close()
    
    return
       
##############################################################################
    
def main():

    folder, nbeads, nsteps, lx, ly, args = read_contextual_info()
    x = read_pos_from_dump_files(folder, nbeads, nsteps, args.tstep, lx, ly, args)
    #savefolder = '/homea/ias2/duman/MASOUD/motor/Data/'
    #savefolder = '/local/duman/SIMULATIONS/motor/Data/'
    #savefolder = '/homec/jiff26/jiff2610/Motor/Data/'
    write_h5_file(args.folder, x, nbeads, nsteps, lx, ly, args)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
