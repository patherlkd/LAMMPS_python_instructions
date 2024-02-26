#!/usr/bin/env python
# coding: utf-8

# In[7]:
import subprocess
import glob as glb
import os
from os.path import exists
from os.path import basename
from pathlib import Path
import math
import matplotlib.pyplot as plt # matplotlib shorthand
import matplotlib.ticker as mtick
import sys # python system library
import numpy as np # numerical python library
import math # python math library
import glob
from matplotlib import cm
import random
#import scipy
#from scipy import ndimage
#from PIL import Image

#from IPython.display import display, HTML
#display(HTML("<style>.container { width:100% !important; }</style>"))

#sys.path.append('/clusternfs/ldavis/ovito-basic-3.7.4-x86_64/lib') # to find ovito
#sys.path.append('/clusternfs/ldavis/ovito-basic-3.7.4-x86_64/lib/ovito/lib') # to find ovito
sys.path.append('/home/ldavis/.localpython/lib/python3.9/site-packages') # to find ovito

#import ovito
from ovito.data import *
from ovito.io import import_file
from ovito.vis import Viewport
from ovito.vis import ParticlesVis
from ovito.vis import VectorVis
from ovito.vis import TachyonRenderer
from ovito.modifiers import ColorByTypeModifier
from ovito.modifiers import AffineTransformationModifier
from ovito.modifiers import ComputePropertyModifier
from ovito.modifiers import AssignColorModifier
from ovito.modifiers import ColorCodingModifier
from ovito.qt_compat import QtCore, QtGui
from ovito.vis import PythonViewportOverlay

from visualize_functions_package import *

import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False

from mpi4py import MPI


ff = 3
space         = ff*0.07
nb_lines      = ff*1
fig_width_pt  = 246.
inches_per_pt = 1./72.
golden_mean   = .66
fig_width     = fig_width_pt*inches_per_pt
fig_height    = (fig_width*golden_mean)+space
fig_size      = [ff*fig_width, ff*fig_height]
# Plot parameters
#plt.rcParams['figure.figsize'] = (7,5)
#plt.rcParams['figure.dpi'] = 150
#plt.rcParams["font.family"] = "Serif"
#plt.rcParams["font.size"] = 22
#params = {'mathtext.default': 'regular' }
params = {'font.family': "Serif",
        'legend.fontsize': ff*7,
          'axes.linewidth': ff*5e-1,
          'axes.labelsize': ff*10,
          #'text.fontsize': ff*4,
          'xtick.labelsize': ff*7,
          'ytick.labelsize': ff*7,
          'text.usetex': True,
          'text.latex.preamble':"\\usepackage{color}",
          'figure.figsize': fig_size}

plt.rcParams.update(params)


#trans_table = [0.95,0.95,0.95,0.0] # for transparency of particles
#trans_table = [0.9,0.9,0.9,0.,0.9]
#trans_table = [0.95,0.95,0.,0.95]
#trans_table = [0.95,0.95,0.95,0.0,0.0,0.0]
#trans_table = [0.0,0.0,0.0,0.0,0.0,0.0]
trans_table = [0.0,0.0,0.0,0.0,0.0,0.0]
ctable = [(1.0, 0.0, 0.0),  # red
    (.5, .5, .5),  # grey
          (0.0, 0.0, 1.0),  # blue
          (0.0, 0.5, 0.0),  # green
          (0.2,0.5,0.0),
          (0.0,0.5,0.2)
          ]

# for visualizing clusters
nclusttypes_max = 500
#trans_table1 = [ 0. for i in range(0,nclusttypes_max) ]
#ctable1 = [ [random.uniform(0.0,0.9),random.uniform(0.0,0.9),random.uniform(0.0,0.9)] for i in range(0,nclusttypes_max) ]

#trans_table += trans_table1
#ctable += ctable1

polylength = 5
polywidth = 1.
bl = '70.5'
viewmode = 'single'
axis = 'X'
cam_dir=(0.0,0.0,0.0)
cam_pos=(0.0,0.0,0.0)

RADIUS = '0.5'
CENTER=False
#SNAPFRAME = -1
SNAPFRAME = 0
#SNAPFRAME = 1
VMSNAP = 1
VMMOV = 0
EVERYNTH = 1
viewtype = 'ortho'
#viewtype = 'ortho'
#extralabelmedia = "_axis"+axis+"_vt-"+viewtype
#extralabelmedia = "_opaque_axis"+axis+"_vt-"+viewtype+"_frm"+str(SNAPFRAME)
extralabelmedia = "_axis"+axis+"_vt-"+viewtype+"_frm"+str(SNAPFRAME)
#extralabelmedia = "_fancyCloser_axis"+axis+"_vt-"+viewtype
#extralabelmedia = "_axis"+axis+"_vt-"+viewtype

#arg1 = sys.argv[1]


names = glb.glob("/local/ldavis/XMA/xyzdata/*.xyz") # grab relevant files

#names += glb.glob("/local/ldavis/Phase_Sep_Polymers/observables/cluster_analysis/*Ls100_*Lc20_*lspa5*Np50_*nruns1000*/*.xyz") # grab relevant files

# cluster coords
#names = glb.glob("/local/ldavis/Phase_Sep_Polymers/observables/cluster_analysis/cg_1apb_Ls45_Lc15_N120_lspa5_eta3.0_cohv2_Np50_bl75_rcoh1.0_sigma0.4_epscoh2.000_dieconst10.0_temp1.0_damp1_nruns10000000_dt0.004_instance0/*group*.xyz") # grab relevant files

names_len = len(names) #names = [names[0]]
names_len_og = len(names)

print("Number of xyz files: "+str(names_len))

# Initialize the MPI environment
MPI.Init()
comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
rank = comm.Get_rank()
names_len = math.floor(names_len/nprocs)*nprocs
assert(names_len%nprocs == 0)
rank_chunk = int((names_len)/nprocs)

rank_end = rank_chunk*(rank+1)
rank_begin = rank_chunk*rank

# ==================== Making of movies ==========================

ninn = 0
max_cut = 200 # hard wall on attempts
max_cut_iter = 0
for nin in range(rank_begin,rank_end):

    if nin > names_len_og:
        continue
    ninn = nin
    name = names[nin]
    bname = basename(name).replace(".xyz","")
    # line = subprocess.check_output(['tail', '-1', "/local/ldavis/Phase_Sep_Polymers/output/out_"+bname])
    # line = line.decode("utf-8")
    # if "Total wall time:" in line:
    #     print("SIMULATION EXISTS AND HAS FINISHED")
    # else:
    #     print("[!] Not finished or defected simulation")
    #     max_cut_iter+=1
    #     if max_cut_iter >= max_cut:
    #         break
    #     continue

    print(bname)
    image_width_fact = 1
    #xbl = get_var_from_name(name,'x',usetokenindex=True,index=4)
    #bl = get_var_from_name(name,'bl')

    #slabwidth = get_var_from_name(name,'slabwidth')

    print("bl: "+bl)
    #print("zybl: "+zybl)
    #print("slabwidth: "+slabwidth)
    FOV = float(bl)*0.5

    if viewtype == 'pers':
        FOV = float(30)
    if 'single' in viewmode:
        if axis == 'X':
            #image_width_fact = float(xbl)/float(zybl)
            cam_dir=(0.0,-1,0.0)
            cam_pos=(0.0,float(bl),0.0)

        if axis == 'Y':
            cam_dir=(0.0,0.0,-1)
            cam_pos=(0.0,0.0,float(bl))
        if axis == 'Z':
            cam_dir=(-1,0.0,0.0)
            cam_pos=(float(bl),0.0,0.0)
        if axis == 'XY':
            #image_width_fact = float(xbl)/float(zybl)
            cam_dir=(-1,-1,0.0)
            cam_pos=(float(bl),float(bl),0.0)


    visualize_simulation(name,label=''+extralabelmedia+'_HPC',plength=float(1),vlx=float(bl),vly=float(bl),vlz=float(bl),vmakesnap=VMSNAP,vviewmovie=0,vmakemovie=VMMOV,every_Nth=EVERYNTH,vradius=RADIUS,camera_dir=cam_dir,camera_pos=cam_pos,camera_fov=FOV
                         ,viewtype=viewtype
                         ,IMAGE_HEIGHT = 600
                         ,IMAGE_WIDTH_FACT = image_width_fact
                         ,ctable = ctable
                         ,drawarrows=False,addbox=False
                         ,trans_modify=True, ttable= trans_table, snapframe=SNAPFRAME
                         ,boxwidth=float(bl),boxheight=float(bl),boxleftpos=0.5-0.5*float(bl)/float(bl),boxtoppos=0.,box_line_width=0.0001,recenter_system=CENTER)


max_cut_iter = 1
## finish off remaining n < 10
ninn=names_len
if rank == 0:
    print(f"FINISHING OFF REMAINING n < {nprocs}")
    while ninn < names_len_og:

        name = names[ninn]
        bname = basename(name).replace(".xyz","")
#         line = subprocess.check_output(['tail', '-1', "/local/ldavis/Phase_Sep_Polymers/output/out_"+bname]\
# )
#         line = line.decode("utf-8")
#         if "Total wall time:" in line:
#             print("SIMULATION EXISTS AND HAS FINISHED")
#         else:
#             print("[!] Not finished or defected simulation")
#             max_cut_iter+=1
#             if max_cut_iter >= max_cut:
#                 break
#             continue

        print(bname)
        print(name)
        image_width_fact = 1
        #xbl = get_var_from_name(name,'x',usetokenindex=True,index=4)
        #bl = get_var_from_name(name,'bl')

        #slabwidth = get_var_from_name(name,'slabwidth')

        print("bl: "+bl)
        #print("zybl: "+zybl)
        #print("slabwidth: "+slabwidth)
        FOV = float(bl)*0.5

        if viewtype == 'pers':
            FOV = float(30)
        if 'single' in viewmode:
            if axis == 'X':
                #image_width_fact = float(xbl)/float(zybl)
                cam_dir=(0.0,-1,0.0)
                cam_pos=(0.0,float(bl),0.0)

            if axis == 'Y':
                cam_dir=(0.0,0.0,-1)
                cam_pos=(0.0,0.0,float(bl))
            if axis == 'Z':
                cam_dir=(-1,0.0,0.0)
                cam_pos=(float(bl),0.0,0.0)
            if axis == 'XY':
                #image_width_fact = float(xbl)/float(zybl)
                cam_dir=(-1,-1,0.0)
                cam_pos=(float(bl),float(bl),0.0)
        visualize_simulation(name,label=''+extralabelmedia+'_HPC',plength=float(1),vlx=float(bl),vly=float(bl),vlz=float(bl),vmakesnap=VMSNAP,vviewmovie=0,vmakemovie=VMMOV,every_Nth=EVERYNTH,vradius=RADIUS,camera_dir=cam_dir,camera_pos=cam_pos,camera_fov=FOV
                         ,viewtype=viewtype
                         ,IMAGE_HEIGHT = 600
                         ,IMAGE_WIDTH_FACT = image_width_fact
                         ,ctable = ctable
                         ,drawarrows=False,addbox=False
                         ,trans_modify=True, ttable= trans_table, snapframe=SNAPFRAME
                             ,boxwidth=float(bl),boxheight=float(bl),boxleftpos=0.5-0.5*float(bl)/float(bl),boxtoppos=0.,box_line_width=0.0001,recenter_system=CENTER)


        ninn+=1
MPI.Finalize()

















# axis = 'Z'
# extralabelmedia = "_axis"+axis

# for nin in range(rank_begin,rank_end):
#     name = names[nin]
#     image_width_fact = 1
#     xbl = get_var_from_name(name,'x',usetokenindex=True,index=4)
#     zybl = get_var_from_name(name,'zybl')
#     slabwidth = get_var_from_name(name,'slabwidth')

#     print("xbl: "+xbl)
#     print("zybl: "+zybl)
#     FOV = float(zybl)*0.5
#     if 'single' in viewmode:
#         if axis == 'X':
#             image_width_fact = float(xbl)/float(zybl)
#             cam_dir=(0.0,-1,0.0)
#             cam_pos=(0.0,float(zybl),0.0)
#         if axis == 'Y':
#             cam_dir=(0.0,0.0,-1)
#             cam_pos=(0.0,0.0,float(xbl))
#         if axis == 'Z':
#             image_width_fact = 1.0
#             cam_dir=(-1,0.0,0.0)
#             cam_pos=(float(xbl),0.0,0.0)


    # visualize_simulation(name,label=''+extralabelmedia+'_HPC',plength=float(1),vlx=float(xbl),vly=float(zybl),vlz=float(zybl),vmakesnap=VMSNAP,vviewmovie=0,vmakemovie=VMMOV,every_Nth=EVERYNTH,vradius=RADIUS,camera_dir=cam_dir,camera_pos=cam_pos,camera_fov=FOV
    #                      ,IMAGE_HEIGHT = 600
    #                      ,IMAGE_WIDTH_FACT = image_width_fact
    #                      ,drawarrows=False,addbox=False, boxwidth=float(zybl),boxheight=float(zybl),boxleftpos=0.0,boxtoppos=0.0)

# close the MPI environment
#MPI.Finalize()
