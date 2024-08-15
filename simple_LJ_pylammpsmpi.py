## system: Simple LJ fluid with the pylammpsmpi package
## Author: Dr. Luke Davis
## Year: 2024
## E: ld731@cam.ac.uk

from __future__ import print_function
import numpy as np
import cupy as cp
import sys
import random
import math
import ctypes
import os
import time
import matplotlib.pyplot as plt
from pylammpsmpi import LammpsLibrary
from lammps import LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR 

# number of mpi threads for the LAMMPS engine
nprocs = int(sys.argv[1])

start_time = time.time()

# Basic simulation parameters
dimension = 2
bl = 100
hbl = bl*0.5
N = 2000
sigma = 1
epsilon = 10
cutoff = 2**(1./6.)*sigma

T =1
damp =1
dt =0.004

Nruns = 10_000
Nthermo = 1
Nxyz = 1

# seeds for random things in LAMMPS script
seed0 = random.randint(10000,1000000)
seed1 = random.randint(10000,1000000)
seed2 = random.randint(10000,1000000)
seed3 = random.randint(10000,1000000)

# instantiate lammps object with 4 processors
lmp = LammpsLibrary(cores=nprocs,mode='local',working_directory=".")
print(f"LAMMPS VERS: {lmp.version}")

# global settings
lmp.command(f"log log.simple_LJ_pylammpsmpi")
lmp.command(f"units lj")
if dimension < 3:
    lmp.command(f"dimension {dimension}")
lmp.command(f"atom_style atomic")
lmp.command(f"atom_modify map yes")
lmp.command(f"boundary p p p")
lmp.command(f"region box block -{hbl} {hbl} -{hbl} {hbl} -{cutoff} {cutoff} side in")
lmp.command(f"create_box 1 box")
lmp.command(f"create_atoms 1 random {N} {seed0} box overlap {cutoff} maxtry 100000")
print(f"Created {N} atoms")

# set the masses
lmp.command(f"mass 1 1")

# Pair interactions
lmp.command(f"pair_style lj/cut {cutoff}")
lmp.command(f"pair_coeff 1 1 {epsilon} {sigma} {cutoff}")

lmp.command(f"neighbor 0.3 bin")
lmp.command("neigh_modify every 1 delay 2 check yes")


# outputs
lmp.command(f"dump 1 all xyz {Nxyz} simple_LJ_pylammpsmpi.xyz")
lmp.command(f"thermo {Nthermo}")
lmp.command(f"thermo_style custom step temp pe etotal press")

# dynamics
lmp.command(f"fix f0 all nve")
lmp.command(f"fix f1 all langevin {T} {T} {damp} {seed3}")
lmp.command(f"fix f2 all enforce2d")
lmp.command(f"timestep {dt}")

# Check global number of atoms
Natomslmps =lmp.get_natoms()
print(f"Natoms: {N}, Natomslmps: {Natomslmps}")

# Test extract compute
lmp.command("compute 1 all reduce min x")
xmin = lmp.extract_compute("1",LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)
print(xmin)

# Run the dynamics!

n = 0
while (n <= Nruns):
    lmp.run(Nthermo)
    n+=Nthermo

# Check timings
print(f"Total Wall Time: {time.time() - start_time} seconds for {nprocs} MPI THREADS") 

# Save a figure of the particles
#== Visualize the final conditions == #                                                                 
#plt.clf()                                                                                               
#plt.cla()                                                                                               
coords = lmp.gather_atoms("x")
#print(coords)
                                                                                                       
typedias = [sigma]
typecols = ['blue']
                                                                                                        
pnts_inches = 72                                                                                        
pnts_whole_panel = bl*bl                                                                                
figx = 8                                                                                                
vis_sizes = [ (pnts_inches**2)*(figx**2)*((typedias[0]/2.0)**2*np.pi)/(pnts_whole_panel\
) for i in range(0,Natomslmps)]                                                                         
vis_cols = [ typecols[0] for i in range(0,Natomslmps)]                                  
                                                                                                        
plt.figure(figsize=(figx,figx))                                                                         
plt.scatter(coords[:, 0], coords[:, 1],s=vis_sizes, c=vis_cols, alpha=0.6)                              
plt.xlim(-hbl, hbl)                                                                                     
plt.ylim(-hbl, hbl)                                                                                     
plt.gca().set_aspect('equal', adjustable='box')                                                         
                                                                                                        
plt.savefig(f"simple_LJ_pylammpsmpi_finalconfig.pdf") 

# close lammps object
lmp.close()
