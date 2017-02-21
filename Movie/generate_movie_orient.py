
""" Generate a movie of the beads of filaments"""

### run:
#  python generate_movie.py -fl=/local/duman/SIMULATIONS/motor/Data/kappa_128/km_0.1/ -ti -tf    

##############################################################################

import argparse
import numpy as np
import os
import h5py
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import misc_tools

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
        self.density = density
        self.kappa = kappa
        self.km = km
        self.panti = pa
        self.bl = bl
        self.sigma = sigma
        
        ### normalize certain variables
        
#        self.lx /= self.bl
#        self.ly /= self.bl   
#        self.dt *= self.nsamp
        
        ### define more simulation parameters
        
#        self.kT = 1
#        self.gamma_n = 1
#        self.N_avg = np.average(self.nbpc)
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
        
        ### assign cell indices to beads
        
        self.cid = np.zeros((sim.nbeads), dtype=np.int32)-1
        k = 0
        for j in range(sim.nfils):
            for n in range(sim.nbpf):
                self.cid[k] = j
                k += 1
                
        return
                
    def calc_orientation(self, step, sim):        
        """ calculate the bond orientations of beads"""
        
        self.phi = np.zeros((sim.nbeads))
        self.phi = misc_tools.compute_orientation(self.x[step,0,:],self.x[step,1,:],sim.lx,sim.ly,sim.nbpf)
        
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
    panti = fl['/param/pa'][...]
    bl = fl['/param/bl'][...]
    sigma = fl['/param/sigma'][...]

    ### close the file

    fl.close()
    
    sim = Simulation(lx, ly, dt, nsteps, nbeads, nsamp, nfils, nbpf, density, kappa, km, panti, bl, sigma)
    beads = Beads(x, sim)
    
    return beads, sim

##############################################################################
        
class Subplots:
    """ plot structure"""
    
    totcnt = -1             # Total number of subplots 
    
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
        return
        
    def addSubplot(self):
        """ add a subplot in the grid structure"""
        
        ### increase the number of subplots in the figure
        
        self.totcnt += 1
        
        ### get indices of the subplot in the figure
        
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])
        
##############################################################################

def plot_frames(beads, sim, ti, tf, save_eps):
    """ plot frames within the specified time window"""

    ### normalise variables
    
    lx = sim.lx/sim.bl
    ly = sim.ly/sim.bl
    
    ### set general plot properties

    savebase = './' #'/usr/users/iff_th2/duman/Motor/MOVIES/'
    savebase += 'density_' + str(sim.density)[0:4] + '_kappa_' + str(sim.kappa) + '_km_' + str(sim.km) + '_panti_' + str(sim.panti) + '/'
    os.system("mkdir -p " + savebase)
    quant_steps = 2056
    norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)  
    downlim = -2
    uplim = lx+2
    num_ticks = 5
    ax_len = 0.9                          # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.1                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    ### plot the frames
    
    for step in range(ti, tf):
        
        time = step*sim.dt
        print 'Step / Total : ', step, tf

        x = beads.x[step, 0, :]/sim.bl
        y = beads.x[step, 1, :]/sim.bl
        phi = beads.calc_orientation(step, sim)
        
        ### plot 

        subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
        ax0 = subp.addSubplot()
        
        line0 = ax0.scatter(x, y, s=1, c=beads.phi, \
                            cmap=plt.cm.get_cmap('hsv',quant_steps), \
                            edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, \
                            norm=norm, rasterized=True)
        
        ### title
        
#        ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/sim.tau_D) + \
#            ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=30)
        
        ### labels
            
        ax0.set_xlabel("$x/r_{0}$", fontsize=40)
        ax0.set_ylabel("$y/r_{0}$", fontsize=40)

        ### limits

        ax0.set_xlim((downlim, uplim))
        ax0.set_ylim((downlim, uplim))
        
        ### ticks
        
        ax0.xaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
        ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
        ax0.tick_params(axis='both', which='major', labelsize=30)
        
        cax0 = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg+ax_len/3, ax_len/4.6, ax_len/4.6], projection='polar')
        xval = np.arange(-np.pi, np.pi, 0.01)
        yval = np.ones_like(xval)
        cax0.scatter(xval, yval, c=xval, s=300, cmap=plt.cm.get_cmap('hsv',quant_steps), norm=norm, linewidths=0)
        cax0.set_xticks([])
        cax0.set_yticks([])
        cax0.set_title('$\\phi$',fontsize=20)
        cax0.set_rlim([-1,1])
        cax0.set_axis_off()        
        
        ### save

        savepath1 = savebase + "frame-" + "{0:05d}".format(int(step)) + ".png"
        if save_eps:
            savepath2 = savebase + "frame-" + "{0:05d}".format(int(step)) + ".eps"
            
        plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
        if save_eps:
            plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)        
        fig.clf()                

        
    return
        
##############################################################################

def main():
    

    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", help="Folder containing data")
    parser.add_argument("-ti","--init_time", nargs="?", const=10, type=int, \
                        help="First frame of the video, you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=2000, type=int, \
                        help="Last frame of the video, you can also leave it empty")
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### read the data and general informaion from the folder
    
    beads, sim = read_data(args.folder)
        
    ### plot the data in the given time window
    
    plot_frames(beads, sim, args.init_time, args.fin_time, args.save_eps)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
