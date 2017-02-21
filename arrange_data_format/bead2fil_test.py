
""" Turn bead information to filament information and save the data to make things faster"""

##############################################################################

import argparse
import numpy as np
import os
import h5py
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math

##############################################################################

def neigh_min(dx, lx):
    """ calculate minimum image distance between neighboring beads"""
    
    dx1 = dx + lx
    dx2 = dx - lx
    if dx**2 < dx1**2 and dx**2 < dx2**2:
        return dx
    if dx1**2 < dx2**2:
        return dx1
    return dx2 

##############################################################################

def neigh_min_AR(dX, lx, ly):
    """ (numpy) calculate minimum image distance between neighbouring beads for an array of displacements"""
    dX = dX.T		# please bear with me, this is easier for Arvind)

    dX1 = dX + np.array([lx,ly])
    dX2 = dX - np.array([lx,ly])

    iter = zip(dX,dX1,dX2)

    X = [dx if (dx**2<dx1**2 and dx**2<dx2**2) else dx1 if (dx1**2<dx2**2) else dx2 for (dx,dy),(dx1,dy1),(dx2,dy2) in iter]
    Y = [dy if (dy**2<dy1**2 and dy**2<dy2**2) else dy1 if (dy1**2<dy2**2) else dy2 for (dx,dy),(dx1,dy1),(dx2,dy2) in iter]

    return np.vstack((X,Y))
 
def neigh_min_array(dx, lx):
    """ calculate minimum image distance between neighboring beads for an array of displacements"""
    
    for j in range(len(dx)):
        dx[j] = neigh_min(dx[j], lx)
    
    return dx
    
##############################################################################        
   
class Simulation:
    """ data structure for storing general simulation information"""

    def __init__(self, lx, ly, dt, nsteps, nbeads, nsamp, nfils, nbpf, density, kappa, km, pa, bl, sigma):
        
        self.lx = lx
        self.ly = ly
        self.dt = dt
        self.nsteps = nsteps
        self.nbeads = nbeads
        self.nsamp = nsamp
        self.nfils = nfils
        self.nbpf = nbpf
        self.kappa = kappa
        self.km = km
        self.bl = bl
        self.sigma = sigma
        self.density = density
	self.panti = pa
        
        ### normalize certain variables
        
#        self.lx /= self.bl
#        self.ly /= self.bl   
#        self.dt *= self.nsamp
        
        ### define more simulation parameters
        
#        self.kT = 1
#        self.gamma_n = 1
#        self.N_avg = np.average(self.nbpf)
#        self.r_avg = self.bl*self.N_avg/2/np.pi
#        self.tau_D = self.r_avg**2 * self.gamma_n * self.N_avg / self.kT
#        self.tau_A = 2 * self.r_avg * self.gamma_n / self.fp
        
        return
        
##############################################################################

class Beads:
    """ data structure for storing particle/bead information"""
    
    def __init__(self, x, sim):
        
        ### assign bead positions
        
        self.x = x
        
        ### assign filament indices to beads
        
        self.cid = np.zeros((sim.nbeads), dtype=np.int32)-1
        k = 0
        for j in range(sim.nfils):
            for n in range(sim.nbpf):
                self.cid[k] = j
                k += 1
        
        return

##############################################################################

class Cells:
    """ data structure for storing cell/filament information"""

    def __init__(self, beads, sim):
        
        ### calculate the centre of mass positions of cells 

        self.com = np.zeros((sim.nsteps, 2, sim.nfils), dtype=np.float32)
        for tstep in range(sim.nsteps):
            print tstep, sim.nsteps
            k = 0
            for j in range(sim.nfils):
                self.com[tstep, 0, j], self.com[tstep, 1, j] = \
                    self.compute_com_single(beads.x[tstep, 0, k:(k+sim.nbpf)], \
                        beads.x[tstep, 1, k:(k+sim.nbpf)], \
                            sim.lx, sim.ly, sim.nbpf)
                k += sim.nbpf
        
        return  

    def compute_filament_orientation(self, x, y, lx, ly, nbeads): #AR
        """ compute the average orientation of filament (only valid for filament) """

	### convention that capital letter variables contain both x and y axis
        ### correct pbcs - in a numpy way when possible (AR)
	X = np.vstack((x,y))
	dX = np.diff(X)
	X[:,1:] = X[:,:-1] + neigh_min_AR(dX,lx,ly)

	### compute dX again, corrected for pbcs
	dX = np.diff(X)

	### mean orientation of filament
	dXmean = dX.mean(axis=1)

	print dXmean
	return dXmean	# returns (dx_mean, dy_mean) as per Ozer's convention

    def compute_com_single_AR(self, x, y, lx, ly, nbeads): #AR
        """ compute the center of mass of filaments per timestep """  

	### convention that capital letter variables contain both x and y axis information
        ### correct pbcs - in a numpy way when possible (AR)
	X = np.vstack((x,y))
	dX = np.diff(X)
	X[:,1:] = X[:,:-1] + neigh_min_AR(dX,lx,ly)

	### compute comx and comy (AR)
	COM = X.mean(axis=1)		# COM is (2,) np array

	### put COM into simulation box
	COM /= np.vstack((lx,ly))
	COM -= np.floor(COM)
	COM *= np.vstack((lx,ly))

	return COM 	# return (comx, comy) as per Ozer's convention
 
    def compute_com_single(self, x, y, lx, ly, nbeads):
        """ compute the centre of mass of cells per timestep, correcting pbc on the fly"""
            
        ### correct pbcs
        
        for i in range(1, nbeads):
            dx = x[i] - x[i-1]
            dy = y[i] - y[i-1]
            x[i] = x[i-1] + neigh_min(dx,lx)
            y[i] = y[i-1] + neigh_min(dy,ly)
            
        ### compute comx and comy
        
        comx = np.average(x)
        comy = np.average(y)
        
        ### put comx and comy into simulation box
        
        comx /= lx
        comx = comx - math.floor(comx)
        comx *= lx
        comy /= ly
        comy = comy - math.floor(comy)
        comy *= ly
    
        return comx, comy
        
    def correct_for_pbc(self, sim):
        """ correct the positions for periodic boundary conditions"""
        
        self.comi = np.zeros((sim.nsteps, 2, sim.nfils), dtype=np.float32)
        self.comi[0, :, :] = self.com[0, :, :]
        for tstep in range(1, sim.nsteps):
            dx = self.com[tstep, 0, :] - self.com[tstep-1, 0, :]
            dy = self.com[tstep, 1, :] - self.com[tstep-1, 1, :]
            self.comi[tstep, 0, :] = self.comi[tstep-1, 0, :] + neigh_min_array(dx, sim.lx)
            self.comi[tstep, 1, :] = self.comi[tstep-1, 1, :] + neigh_min_array(dy, sim.ly)            

        return        
        
##############################################################################

def read_data(folder):
    """ read simulation data through hdf5 file"""

    ### access the file
    
    fpath = folder + '/out.h5'
    assert os.path.exists(fpath), "out.h5 does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')
    
    ### read in the positions of beads
    
    x = np.array(fl['/positions/x'], dtype=np.float32)

    ### read in the box info

    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]

    ### read in the general simulation info
    
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]

    ### read in the filament information
    
    nfils = fl['/info/nfils'][...]
    nbpf = fl['/info/nbpf'][...]

    ### read in the simulation parameters
    
    density = fl['/param/density'][...]
    kappa = fl['/param/kappa'][...]
    km = fl['/param/km'][...]
    bl = fl['/param/bl'][...]
    pa = fl['/param/pa'][...]
    sigma = fl['/param/sigma'][...]

    ### close the file

    fl.close()

    sim = Simulation(lx, ly, dt, nsteps, nbeads, nsamp, nfils, nbpf, density, kappa, km, pa, bl, sigma)
    beads = Beads(x, sim)
    fils = Cells(beads, sim)
    
    return beads, fils, sim
    
##############################################################################

def write_fil_data(fils, sim, folder):
    """ write filament data to hdf5 files in corresponding folders"""
    
    ### data file to write to
    
    fpath = folder + '/out_fil.h5'
    fl = h5py.File(fpath, 'w')
    
    ### positions of cells
    
    pos = fl.create_group('positions')
    pos.create_dataset('x', (sim.nsteps, 2, sim.nfils), data=fils.com, dtype=np.float64, compression='gzip') 
    pos.create_dataset('xi', (sim.nsteps, 2, sim.nfils), data=fils.comi, dtype=np.float64, compression='gzip') 
    
    ### simulation information
    
    info = fl.create_group('info')
    box = info.create_group('box')
    box.create_dataset('x', data=sim.lx)
    box.create_dataset('y', data=sim.ly)
    info.create_dataset('dt', data=sim.dt)
    info.create_dataset('nsteps', data=sim.nsteps)
    info.create_dataset('nbeads', data=sim.nbeads)
    info.create_dataset('nsamp', data=sim.nsamp)
    
    ### simulation parameters
    
    param = fl.create_group('param')
    param.create_dataset('kappa', data=sim.kappa)
    param.create_dataset('km', data=sim.km)
    param.create_dataset('bl', data=sim.bl)
    param.create_dataset('sigma', data=sim.sigma)
    param.create_dataset('density', data=sim.density)
    
    ### filament parameters
    
    info.create_dataset('nfils', data=sim.nfils)
    info.create_dataset('nbpf', data=sim.nbpf)
    
    fl.close()    
        
    return
        
##############################################################################

def main():

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", help="Folder containing data")
    args = parser.parse_args()
    
    ### read the data and general informaion from the folder
    
    beads, fils, sim = read_data(args.folder)
        
    ### write the cell data to hdf5 files
    
    fils.correct_for_pbc(sim)
    write_fil_data(fils, sim, args.folder)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
