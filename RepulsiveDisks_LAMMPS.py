## This uses solely a Python (LAMMPS API) approach
## system: WCA [repulsive] Brownian Particles
## Author: Dr. Luke Davis
## Year: 2022
## E: luke.davis@ucl.ac.uk

from __future__ import print_function
import numpy as np
import sys
import random
import math
import ctypes
import os

wd = "./"
xyzdata_dir = f"{wd}xyzdata/"

try:
    os.mkdir(xyzdata_dir)
except FileExistsError:
    # directory already exists so ignore
    pass

## SIMULATION PARAMETERS
dimension = 2 # spatial dimension
bl = 80 # box length
hbl = 0.5*bl
pack_frac = 0.01 # requested packing fraction


T = 1.0 # reduced temperature (all units are reduced according to LAMMPS units LJ)
dt = 0.004 #timestep
Nthermo = 10 # display thermodynamic information every this number of steps
Nxyz = 1000 # dump xyz file every this number of steps
Nruns = 10000 # number of total timesteps
damp = 1.0 # inverse mobility
instance = 0 # label to distinguish different random runs
#active_types = [2]

# INTERACTION POTENTIAL PARAMETERS
epsilon = 10.0
sigma = 1.0
rad = sigma*0.5
HS_sigma = round(sigma*(2./( 1. + math.sqrt(T/epsilon)) )**(1./6.),3) # approximate HS diameter
HS_rad = 0.5*HS_sigma
cutoff = 2.**(1./6.)*sigma # this default value gives the WCA potential i.e. only repulsion (rc = 2^(1/6)*sigam ~ 1.122 when sigma = 1)
comm_cutoff = cutoff*1.2
vol_atom = 0.
if dimension == 2:
    vol_atom = math.pi*HS_rad**2
Natoms = math.floor((pack_frac*(bl**dimension))/vol_atom) # number of atoms 
dens = Natoms/(bl**dimension) # density
unique_sim_label = f"RABPS_LAMMPS_N{Natoms}_bl{bl}_HSsig{HS_sigma}_T{T}_damp{damp}_dt{dt}_Nruns{Nruns}_instance{instance}"

# seeds for random things in LAMMPS script
seed0 = random.randint(10000,1000000)
seed1 = random.randint(10000,1000000)
seed2 = random.randint(10000,1000000)
seed3 = random.randint(10000,1000000)


# BEGIN LAMMPS SIMULATION
#from lammps import lammps, PyLammps
## from mpi4py import MPI
from lammps import lammps
lmp = lammps()
#lmp.file(infile)
#L =PyLammps(ptr=lmp)

##comm = MPI.COMM_WORLD
##rank = comm.Get_rank()
##nprocs = comm.Get_size()


lmp.command(f"units lj")
if dimension < 3:
    lmp.command(f"dimension {dimension}")
lmp.command(f"atom_style atomic")
lmp.command(f"atom_modify map yes")
lmp.command(f"boundary p p p")
lmp.command(f"region box block -{hbl} {hbl} -{hbl} {hbl} -{rad} {rad} side in")
lmp.command(f"create_box 1 box")

# initialize configuration of non-overlapping hard disks
#lmp.command(f"create_atoms 1 random {Natoms} {seed0} box overlap {cutoff} maxtry 1000")

# create lattice of particles (then only place Natoms of those)

ncnt = 0
z = 0.
y = -hbl + cutoff*0.5
while y <= hbl - cutoff*0.5 and ncnt < Natoms:
    x = -hbl + cutoff*0.5
    while x <= hbl - cutoff*0.5 and ncnt < Natoms:
        x += cutoff
        lmp.command(f"create_atoms 1 single {x} {y} {z}")
        ncnt+=1
        #lattice_coords.append([x,y,z])
    y += cutoff


    
lmp.command(f"mass 1 1.")

# interactions
lmp.command(f"pair_style lj/smooth/linear {cutoff}")
lmp.command(f"pair_coeff * * {epsilon} {sigma}")

# outputs
lmp.command(f"dump 1 all xyz {Nxyz} {xyzdata_dir}{unique_sim_label}.xyz")
lmp.command(f"thermo {Nthermo}")
lmp.command(f"thermo_style custom step temp pe etotal press ")

#parallel stuff
#lmp.command(f"comm_style tiled")
#lmp.command(f"comm_modify cutoff {comm_cutoff}")

lmp.command(f"neighbor 0.3 bin")
lmp.command("neigh_modify every 1 delay 2 check yes")

# fixes
lmp.command(f"fix f0 all nve")
lmp.command(f"fix f1 all langevin {T} {T} {damp} {seed2}")
lmp.command(f"fix f2 all enforce2d")

# running params
lmp.command(f"timestep {dt}")

Natoms =lmp.extract_global("nlocal")
#if rank == 0:

print(f"\n\nPython Natoms = {Natoms}\n\n")

# implement dynamics
n = 0
while(n<Nruns-1):
    #L.run(1)
    #print(f"{L.runs[n].thermo.Step[0]} {L.runs[n].thermo.Temp[0]} ")
    lmp.command(f"run 1 pre no post no")
    n+=1
lmp.command(f"run 1 pre no post yes")

#print("Proc %d out of %d procs has" % (rank,nprocs),lmp)
#MPI.Finalize()
