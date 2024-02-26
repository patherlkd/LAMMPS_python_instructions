(base) ldavis@dps124:/local/ldavis/XMA/visualisation$ cat visualize_functions_package.py
#update
import glob as glb
import os
from os.path import exists
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
from functools import partial

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

polylength = 5
polywidth = 1.

def get_var_from_name(name,var_string,usetokenindex=False,index=0):
    tokens = name.split("_")
    var_value_string = ''
    if not usetokenindex:

        for token in tokens:
            if var_string in token:
                var_value_string = token.replace(var_string,"")
                return var_value_string
    else:
        var_value_string = tokens[index].replace(var_string,"")
        return var_value_string

def add_ellipsoid_modify(frame, data,ellipsoid_ABC):
    N = data.particles.count
    ptypes = data.particles_.particle_types_
    ellip_index = data.particles_.add_particle((0.0,0.0,0.0))
    ptypes[ellip_index] = 5
    ptypes.types.append(ParticleType( id = 5, name = "E"))

    data.particles_.create_property('Aspherical Shape')
    p_shapes = data.particles_['Aspherical Shape_']
    p_shapes[ellip_index] = ellipsoid_ABC


def transparency_modify(frame,data,trans_table):
    tprop = data.particles_.particle_types_
    data.particles_.create_property('Transparency',data=np.zeros(data.particles.count))
    trans_list = data.particles_['Transparency_']
    for i in range(data.particles.count):
        trans_list[i] = trans_table[tprop[i]-1]
    
def color_modify(frame, data,color_table):


    tprop = data.particles_.particle_types_
    color = data.particles_.create_property('Color')
    for i in range(data.particles.count):
        color[i] = color_table[tprop[i]-1]

vector_vis = VectorVis(
    alignment = VectorVis.Alignment.Base,
    color = (1,0.54,0.01),
    width = polywidth,
   # scaling = 1.5, # for Nm = 5
    scaling = 1., # for Nm = 37
    flat_shading = False,
    transparency = 0.0
)

def recenter_system_modify(frame,data):

    rcom = np.array([0.0,0.0,0.0])
    i = 0

    while(i < data.particles.count):

        if i == 0:
            rcom[0] = data.particles.positions[i][0]
            rcom[1] = data.particles.positions[i][1]
            rcom[2] = data.particles.positions[i][2]
        else:
            rcom += data.particles.positions[i] 
        i+=1

    rcom /= data.particles.count
    # translate all particles such that rcom will be 0.0,0.0,0.0 again
    data.particles_.positions_[...] -= rcom
    
def display_rods_as_arrows_modify(frame,data):

    arrow_data = []
    hbl = 0.5*float(bl)
    i = 0
    tailcnt = 0
    headcnt = tailcnt + (polylength-1)
    while(i < data.particles.count):

        if i==tailcnt:
            dr = data.particles.positions[i+(polylength-1)] - data.particles.positions[i]
            ddr = math.sqrt(dr[0]*dr[0] + dr[1]*dr[1] +dr[2]*dr[2])
            if ddr <= polylength:
                arrow_data.append(dr)
            else:
                arrow_data.append([0.0,0.0,0.0]) # i.e. just ignore arrows that extend hbl or more

            tailcnt += polylength
        else:
            arrow_data.append([0.0,0.0,0.0])
        i += 1
    arrow_data = np.array(arrow_data)
    arrow = data.particles_.create_property('Arrow',data=arrow_data)
    arrow.vis = vector_vis

def visualize_simulation(name,label="",plength=5,vlx=100,vly=100,vlz=100,vradius="0.55",vmakemovie=1,vmakesnap=1,vviewmovie=0,every_Nth=1,camera_pos=(1,1,1),camera_dir=(1,-0.5,-1),camera_fov=math.radians(60.0),viewtype='ortho',
                         IMAGE_HEIGHT = 600,
                         IMAGE_WIDTH_FACT = 1.33,
                         ctable = [(.5, .5, .5),  # grey
                                   (1.0, 0.0, 0.0),  # red
                                   (0.0, 0.0, 1.0),  # blue
                                   (0.0, 0.5, 0.0),  # green
                                   ],
                         drawarrows=False
                         ,trans_modify = False, ttable = [ ]
                         , snapframe=-1
                         ,addellipsoid=False,ellipsoid_abc=[0.0,0.0,0.0]
                         ,addbox=False, boxwidth=1,boxheight=1,boxleftpos=0.0,boxtoppos=0.0,boxcolor=QtGui.QColor(0,0,0)
                         ,box_line_width=0.1
                         ,piplinestruct=["Particle Type", "Position.X", "Position.Y", "Position.Z"]
                         ,recenter_system=False
                         ,othernames=[]
                         ,othershapes=[] # "Sphere","Box", "Circle", "Square", "Cylinder", "Spherocylinder"
                         ,otherradii=[]):

    # == Generate the data pipeline ==
    pipeline = import_file(name, columns = piplinestruct)

    pipeline.add_to_scene()
    pipeline.source.num_frames

    otherpipelines = []
    othercnt = 0
    for othername in othernames:
        temp_pipeline = import_file(othername, columns = piplinestruct)
        temp_pipeline.add_to_scene()
        
        temp_pipeline.modifiers.append(ComputePropertyModifier(output_property="Radius", expressions=[otherradii[othercnt]]))
        vis_element = temp_pipeline.compute().particles.vis

        shape = othershapes[othercnt]

        if shape == "Sphere":
            vis_element.shape=ParticlesVis.Shape.Sphere
        if shape == "Box":
            vis_element.shape=ParticlesVis.Shape.Box
        if shape == "Circle":
            vis_element.shape=ParticlesVis.Shape.Circle
        if shape == "Square":
            vis_element.shape=ParticlesVis.Shape.Square
        if shape == "Cylinder":
            vis_element.shape=ParticlesVis.Shape.Cylinder
        if shape == "Spherocylinder":
            vis_element.shape=ParticlesVis.Shape.Spherocylinder
            
        otherpipelines.append(temp_pipeline)
        othercnt += 1

    totframes = pipeline.source.num_frames-1
    if snapframe == -1:
        snapframe = totframes
        
    pipeline.modifiers.append(ComputePropertyModifier(output_property="Radius", expressions=[vradius]))

    polylength = plength

    if recenter_system:
        pipeline.modifiers.append(recenter_system_modify)
    if drawarrows:
        pipeline.modifiers.append(display_rods_as_arrows_modify)
    if addellipsoid:
        pipeline.modifiers.append(partial(add_ellipsoid_modify,ellipsoid_ABC=ellipsoid_abc))

    pipeline.modifiers.append(partial(color_modify,color_table=ctable))
    if trans_modify:
        pipeline.modifiers.append(partial(transparency_modify,trans_table=ttable))
    
    
    # == Setup the simulation cell ==
    lx = vlx; ly = vly; lz = vlz


    def cell_modifier(frame,data):
        data.cell_[:,0] = (lx, 0, 0)
        data.cell_[:,1] = (0, ly, 0)
        data.cell_[:,2] = (0, 0, lz)
        data.cell_[:,3] = np.dot((-0.5, -0.5, -0.5), data.cell[:3,:3])
        data.cell_.pbc = (True, True, True)
        data.cell_.vis.rendering_color = (0.5, 0.5, .5)
        data.cell_.vis.line_width=box_line_width
        data.cell_.vis.enabled = True

    pipeline.modifiers.append(cell_modifier)
#data = pipeline.compute()
    #list(data.particles.keys())


    data = pipeline.compute()

    def render_box(args: PythonViewportOverlay.Arguments):
        if args.is_perspective:
            raise Exception("This only works in ortho mode")

        # convert box lengths to screen units
        s_boxwidth = args.project_size((0,0,0),boxwidth)
        s_boxheight = args.project_size((0,0,0),boxheight)

        s_boxleftpos = boxleftpos*args.painter.window().width()
        s_boxtoppos = boxtoppos*args.painter.window().height()
        box = QtCore.QRectF(s_boxleftpos,s_boxtoppos,s_boxwidth,s_boxheight)

        pen = QtGui.QPen(boxcolor)
        pen.setWidth(5)
        pen.setStyle(QtCore.Qt.DashDotLine)
        args.painter.setPen(pen)

        args.painter.drawRect(box)



    # == Setup the Viewport part ==
    vp= Viewport()
    if viewtype == 'ortho':
        vp.type = Viewport.Type.Ortho
    else:
        vp.type = Viewport.Type.Perspective
    vp.camera_pos = camera_pos
    vp.camera_dir = camera_dir
    vp.fov = camera_fov
    #vp.zoom_all()

    # == Add box to show slab == #
    if addbox:
        vp.overlays.append(PythonViewportOverlay(function=render_box))

    # == Tachyon renderer setup for movie ==
    tachymov=TachyonRenderer(shadows=False,ambient_occlusion=False,direct_light_intensity=1.1,depth_of_field=False)


    # == Tachyon renderer setup for movie ==
    tachysnap=TachyonRenderer(shadows=True,ambient_occlusion=True,direct_light_intensity=1.1,depth_of_field=False)

    name = Path(name).stem

    # == Make movies and/or snap ==

    if(vmakesnap):
        print("Making snap: "+name+"\n")
        vp.render_image(size=(int(IMAGE_WIDTH_FACT*IMAGE_HEIGHT),IMAGE_HEIGHT), filename=name+label+f"_totfrmes{totframes}.png", background=(1,1,1), frame=snapframe,renderer=tachysnap)

    if(vmakemovie):
        print("Making movie: "+name+"\n")
        vp.render_anim(size=(int(IMAGE_WIDTH_FACT*IMAGE_HEIGHT),IMAGE_HEIGHT), filename=name+label+f"_totfrmes{totframes}"+".mp4", fps=20,every_nth=every_Nth,background=(1,1,1),renderer=tachymov)


    # == View the movie ==
    if(vmakemovie and vviewmovie):
        os.system("vlc "+name+"_animation.mp4 > /dev/null 2>&1");
    pipeline.remove_from_scene()
