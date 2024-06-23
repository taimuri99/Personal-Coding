

''' 
 * This program should carry out a simple Monte Carlo method
 * for molecules falling and diffusing on a surface (DDA-model).
 '''


import time
import sys
import numpy as np
from surface_objects import *   # explicitly import all the function defined in the module "surface_objects": 
                                # this means that everything defined there, is available in this namespace, thus avoiding the "module.object()" syntax 
from clusters_analysis import *
from math import *
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
    
print("Version 1.1 (Debugged and pythonized) by:\nMattia Fiorentini, January 2016\n\n")
print("(Orignal C version 0.1 by Elliot Gross, February 2012)\n\n")

# Initialization
random.seed()
n_pic = 0
max_occ = 0   
percent = 0
t_total = 0
inv_kb=11604.51911 # Inverse Boltzmann constant (K/eV)
# Duration
itt_max = 600 # number of Monte Carlo steps
fps = 100# output period [number of MC steps]
# Ambient
temperature = 300 # [K]
frequency = 0.01 # [fs^-1]
# Deposition
drop = 200# number of droplets that will be deposited on the surface 
dep_rate = 0.001 # deposition rate [fs^-1]
# Droplets type
delta_stack = 0.1 # energy barrier of particle diffusion [eV]
energy_ev_stack = 0.4 # energy barrier of aprticle desorbtion [eV]
energy_bind = 0.5 # binding energy between particles [eV]
delta_climb = 1# energy barrier of particle diffusing by climb on top of another one [eV]
# Substrate of diffusion
gx = 20 # number of sites along the X dir.
gy = 20 # number of sites along the y dir.
max_height = 100 # maximum number of dropets allowed to stack on each-other
delta_surface = 0.2 # energy barrier for surface difusion [eV]
energy_ev_surface = 0.9 # energy barrier for particle desorbsion from the surface [eV]

dep_proc = Deposition(drop,dep_rate,temperature) # defines the deposition procedure
particle = Droplet(delta_stack,energy_ev_stack,energy_bind,delta_climb) # defines the type of particles that will be deposited on the substrate
substrate = Surface([gx,gy,max_height],delta_surface,energy_ev_surface) # defines the substrate
reservoir = Surface([1,1,drop],0,0) # defines the tipe of the experimental dep_proc that will deposit the particles on the substrate
tip = reservoir.get_site([0,0]) # this is the site on the of the tip on which the particle that will be deposited will be stored momentarily
for i in range(drop):
    reservoir.add_part(tip,particle)
atmosphere = Surface([1,1,drop],0,0) # this is a ficticious surface object, the atmosphere, to which the desorbed particles will be stored
    
plt.ion()
fig, (axs, axh) = plt.subplots(nrows=1, ncols=2, figsize=(20, 12))
for itt in range(itt_max + 1):
    '''Fill "moves_list", with Move obj. It is a list containing all the possible moves'''
    # 1) Dropping of a new droplet on the surface
    tip = reservoir.get_site([0,0])
    if tip.occ == 0:
        moves_list = []
    else:
        drop_site = substrate.get_site([random.randrange(substrate.lx),random.randrange(substrate.ly)])
        while drop_site.occ == max_height:
            drop_site = substrate.get_site([random.randrange(substrate.lx),random.randrange(substrate.ly)])
        moves_list = [Move(dep_rate,reservoir,tip,substrate,drop_site)]
    # 2.0) Movement of each particle on the substrate
    i = 0
    for o_site in substrate:
        s_surf = substrate
        s_site = o_site
        i+=1
        first_neigh_sites = substrate.site_neighbours(o_site)
        proxy = 0 
        for neighbour in first_neigh_sites:
            if neighbour.occ == o_site.occ: proxy += 1
        j = 0  
        # 2.1) Diffusion
        for neighbour in first_neigh_sites: # diffusion  can happen ...
            j+=1
            if neighbour.occ < o_site.occ: # ... in plane/falling, or ...
                delta_total = neighbour.de_diff + proxy*particle.de_bind
            elif neighbour.occ == o_site.occ: # ... climbing.
                delta_total = neighbour.de_diff + proxy*particle.de_bind + particle.de_climb
            rate = frequency * exp(-delta_total*inv_kb/dep_proc.t)
            e_surf = substrate
            e_site = neighbour
            moves_list.append(Move(rate,s_surf,s_site,e_surf,e_site))
        # 2.2) Desorbation
        delta_total = o_site.de_desorb + proxy*particle.de_bind
        rate = frequency * exp(-delta_total*inv_kb/dep_proc.t) # rate for this moved
        e_surf = atmosphere
        e_site = atmosphere.get_site([0,0]) # final site,  to which particle will be added
        moves_list.append(Move(rate,s_surf,s_site,e_surf,e_site))
    
    '''Pick one of the Move obj using a "Monte Carlo" procedure'''
    # Total rate: sum each move's rate
    R=0                             
    for move in moves_list:
        #print("Rate ",move.rate)
        R += move.rate
    #print("R",R)
    # Pick a random number between (0,1)
    chi = 0.0
    while chi == 0.0:
        chi = random.random()
    #print("Chi ",chi)
    t = log(1/chi)/R     # time interval
    t_total+=t           # total time taken
    # Realize the corresponding move
    chi = random.random()
    norm_rsum=0                        
    for move in moves_list:
        norm_rsum += move.rate/R
        if norm_rsum >= chi:
            move.in_surf.remove_part(move.in_site)
            move.fi_surf.add_part(move.fi_site,particle)  
            break
   
    if 100*itt//itt_max > percent:
        percent = 100*itt//itt_max
        sys.stdout.write("\r%d%%\n" % percent )
        sys.stdout.flush() 

    # analysing and plotting a snapshot of the surface
    if itt%fps == 0 or 100*itt//itt_max == 100:
        fig.clf()
        axh = fig.add_subplot(1,2,2)
        N_clusters, heights, areas, volumes = surface_counts(substrate)
        multi_h = [heights, areas, volumes]
        n_bins = max(volumes)
        axh.hist(multi_h, n_bins, histtype='bar', label=("Height", "Area", "Volume"))
        axh.legend(fontsize = 'medium', loc= 1)
        axh.set_xticks([i+1 for i in range(n_bins)])
        axh.set_xticklabels([i+1 for i in range(n_bins)])

        '''
        N_clusters, height_distr, area_distr, vol_distr = surface_stats(substrate)
        fd_name = 'surf_test' + '{0:d}x{1:d}_{2:d}'.format(substrate.lx,substrate.ly,n_pic) + '.dat'
        fd = open(fd_name, 'w')
        fd.write('# Key, height distr., area distr., volume distr.\n')
        axh = 
        for i in range(max(vol_distr.keys())): # volume of a cluster is always larget than or equal to its height and area
            k = i+1
            pl =  '{0:d} {1:d} {2:d} {3:d}\n'.format(k,height_distr.get(k,0), area_distr.get(k,0), vol_distr.get(k,0))
            fd.write(pl)
        fd.close
        '''    
        n_pic += 1
        print("Snapshot n {0:d}".format(n_pic))
        tiling = np.zeros((substrate.lx,substrate.ly))
        for s in substrate:
            tiling[s.coord.x,s.coord.y] = s.occ
            if s.occ > max_occ: max_occ = s.occ
        scale = [el for el in range(max_occ+1)]
        
        axs = fig.add_subplot(1,2,1)
        im_content = axs.imshow(tiling, cmap=cm.afmhot, interpolation='nearest')
        cbar = fig.colorbar(im_content, ax=axs, ticks=scale)
        #im_content.set_data(tiling)
        #im_content.autoscale()
        #cbar.set_ticklabels(scale)
        #cbar.update_normal(im_content)
        plt.title(' Time = {0:.4e} ({1:d}), Clusters # = {2:d}'.format(t_total,itt,N_clusters))
        plt.draw()
        fig_name = 'surf_test' + '{0:d}x{1:d}_{2:d}'.format(substrate.lx,substrate.ly,n_pic) + '.png'
        plt.savefig(fig_name, bbox_inches='tight')
        
          

percent = 100
sys.stdout.write("\r%d%%\n" % percent )
sys.stdout.flush()

