
#import modules
import sys
import os
import os.path as osp
from pathlib import Path
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import remesh
import utils as ut
import cutting
import volume_mesh
from capping import cap
import identification as id
from tkinter import Tk
from tkinter.filedialog import askdirectory
import mapping
import xml.etree.ElementTree as ET
import subprocess
import febioxml as feb
import qc
#----------------------------------------------------------------------------------------------------------------------------
# Setup
#----------------------------------------------------------------------------------------------------------------------------


file_dir = osp.dirname(osp.realpath(__file__))
temp_dir = osp.join(file_dir, r'temp')
output_dir = osp.join(file_dir, r'output')
input_dir = askdirectory(title='Select Input Folder') # shows dialog box and return the path 
vel_profile_dir = askdirectory(title='Select Velocity Profile Folder') # shows dialog box and return the path   

os.makedirs(temp_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

inlet_path = osp.join(input_dir, r'inlet.stl')
wall_path = osp.join(input_dir, r'wall.stl')
outlet_path = osp.join(input_dir, r'outlet.stl')

print('Created necessary files and directories')
#paths to the aorta geometry (for standard testing)
#inlet_path = "geometries\input\inlet.stl"
#wall_path = "geometries\input\wall.stl"
#outlet_path = "geometries\input\outlet.stl"

#Reading the files with pyvista
inlet = pv.read(inlet_path)
wall = pv.read(wall_path)
outlet = pv.read(outlet_path)

#Meshing parameters
mmg_parameters = {
    'mesh_density': '0.1', #hausdorf parameter of mmg, defines amount of added detail at curvature
    'sizing': '1', #forces similarly sized poly's, legacy option (used for stitching), can possibly be dropped
    'detection angle': '35'}

mmg3d_parameters = {
    'hausd': '0.1',
    'detection angle': '35'}

mmg3d_sol_parameters = {
    'bl_thickness': 1,
    'bl_edgelength': 1,
    'edgelength': 15}

tetgen_parameters = dict(
    order=1, 
    mindihedral=20, 
    minratio=1.5,
    nobisect=True,
    fixedvolume=True,
    maxvolume=1)

#Run qualifications (code will terminate if not met)
max_elements = 1000000
min_jacobian = 0.1
max_aspect   = 5

#Angle for identification of surfaces
id_angle = 35

#Mapping interpolation options
intp_options = {
    'zero_boundary_dist': 0.2,  # percentage of border with zero velocity (smooth damping at the border)
    'zero_backflow': False,     # set all backflow components to zero
    'kernel': 'linear',         # RBF interpolation kernel (linear is recommended)
    'smoothing': 0.5,           # interpolation smoothing, range recommended [0, 2]
    'degree': 0,                # degree of polynomial added to the RBF interpolation matrix
    'hard_noslip': False}       # check if no-slip condition on walls is met


#Plotting boolean, when True: code generates intermediate plots of workflow
show_plot = False

print('Setup and import done')

#--------------------------------------------------------------------------------------------------------------------------
# 3D-meshing algorithm
#--------------------------------------------------------------------------------------------------------------------------

#Cut the wall geometry after the aortic root
wall_cut = cutting.main_cutter(inlet, wall, plot=show_plot)
pv.save_meshio(osp.join(temp_dir, r'wall_cut.mesh'), wall_cut)

#Create caps
inlet_cap, outlet_cap = cap(wall_cut, plot=show_plot)
pv.save_meshio(osp.join(temp_dir, r'inlet_cap.mesh'), inlet_cap)
pv.save_meshio(osp.join(temp_dir, r'outlet_cap.mesh'), outlet_cap)

#Combine cutted wall and inlet/outlet caps
combined = (wall_cut + inlet_cap + outlet_cap).clean()
print('Meshes succesfully combined')

#Plot result of mesh combining. When edgetest gives an error, the structure has dublicate faces/nodes
if show_plot:
    plt = pv.Plotter()
    plt.add_mesh(combined, style='wireframe')
    edgetest = combined.extract_feature_edges(boundary_edges=True, non_manifold_edges=True, manifold_edges=False, feature_edges=False)
    plt.show()

pv.save_meshio(osp.join(temp_dir, r'combined_mesh.mesh'), combined)



#Run remesh (takes predetermined internally defined file path as input, DON'T CHANGE)
combined_remeshed = remesh.remesh_edge_detect(osp.join(temp_dir, r'combined_mesh.mesh'), osp.join(temp_dir, r'combined_mmg.mesh'), temp_dir, mmg_parameters, plot=show_plot)

#triangulation step to make sure Tetgen only gets triangles as input
combined_remeshed = combined_remeshed.extract_surface().triangulate()

#Report quality
qc.meshreport(combined_remeshed, 'Surface mesh quality report')

#Make an initial 3D mesh from the combined mesh using TetGen
combined_remeshed = combined_remeshed
tetmesh = volume_mesh.tetgen(combined_remeshed, tetgen_parameters, plot=show_plot)

#Plot bisection
if show_plot:
    qc.clip_plot(tetmesh, 'Initial 3D mesh')

#Report quality
qc.meshreport(tetmesh, 'Initial 3D mesh quality report')

#Create a .sol file for mmg3d
volume_mesh.write_sol(tetmesh, wall_cut, mmg3d_sol_parameters, osp.join(temp_dir, r'initial_volume_mesh.sol'), plot=show_plot)

#Save initial mesh
tetmesh.point_data.clear()
tetmesh.cell_data.clear()
pv.save_meshio(osp.join(temp_dir, r'initial_volume_mesh.mesh'), tetmesh)

#Refine 3D mesh with mmg3d
tetmesh = volume_mesh.mmg3d(osp.join(temp_dir, r'initial_volume_mesh.mesh'), osp.join(temp_dir, r'mmg3d_mesh.mesh'),
                                  temp_dir, mmg3d_parameters, plot=show_plot)

#Report quality
report = qc.meshreport(tetmesh, 'Refined 3D mesh quality report')

#Plot bad cells
jac = report['jac']
if show_plot:
    jac = report['jac']
    if np.any(jac<0.3):
        plt = pv.Plotter()
        plt.add_mesh(tetmesh, style='wireframe')
        plt.add_mesh(tetmesh.extract_cells(jac<0.3), color='red', show_edges=True)
        plt.add_text('Bad cells')
        plt.show()
    else: print('No bad cells')

#Plot bisection
if show_plot:
    qc.clip_plot(tetmesh, 'Refined 3D mesh')

#Save 3D mesh
tetmesh.save(osp.join(temp_dir, r'3D_output_mesh.vtk'))

#Run qualification
numcells = report['cells']
aspect = report['aspect']
if numcells > max_elements or any(aspect > max_aspect) or any(jac < min_jacobian): run = False
else: run=True

#Plot 3D_mesh
if show_plot:
    if run:text='Final 3D mesh'
    else:text='Final 3D mesh - insufficient quality or too many nodes'
    tetmesh.plot(show_edges = True, text=text)

if run==False:
    print('Terminating, 3D mesh insufficient quality or too many nodes')
    exit()

#----------------------------------------------------------------------------------------------------------------------------
# Identification
#----------------------------------------------------------------------------------------------------------------------------

#Create the seeds for the surface identification based on the center points. INLET FIRST!, OUTLET SECOND!
seeds = np.array([inlet_cap.points.mean(0),outlet_cap.points.mean(0)])

#Detect the surfaces of the 3D mesh whilst keeping the original ID's, if there is a seed delivered
if len(seeds) > 0:
    surface_identification = id.identify_surfaces(tetmesh, id_angle, seeds, show_plot)

#Seperate the indentified surfaces in inlet/outlet/wall
id_inlet = surface_identification[0]
id_outlet = surface_identification[1]
id_wall = surface_identification[2]

print('Surface identification done')
#Plot the identified surfaces for general overview
if show_plot:
    plt = pv.Plotter()
    plt.add_mesh(id_inlet, show_edges = True, color = 'red')
    plt.add_mesh(id_outlet, show_edges = True, color = 'blue')
    plt.add_mesh(id_wall, show_edges = True, color = 'green')
    plt.show()

#----------------------------------------------------------------------------------------------------------------------------
# Mapping
#----------------------------------------------------------------------------------------------------------------------------

#Perform the mapping of the velocity profiles on the inlet
#Output is a point cloud on every inlet node with the respective velocity data and the amount of mapped velocity profiles

velocity_mapped, n_maps = mapping.vel_mapping(vel_profile_dir, id_inlet, output_dir, intp_options, show_plot)
single_profile = velocity_mapped[5]
#----------------------------------------------------------------------------------------------------------------------------
# FEBio
#----------------------------------------------------------------------------------------------------------------------------

#Create a solver compatible file based on the 3D-mesh and meshing parameters
feb.xml_creator(tetmesh, id_inlet, id_outlet, id_wall, file_dir, output_dir)

#Run FEBio
#FEBio_path = r"C:/Program Files/bin/febio4.exe" #Path voor Yarran
FEBio_path = r"C:/Program Files/bin/febio4.exe" #Path voor normale mensen
#Use the current
FEBio_inputfile = osp.join(output_dir, r'simulation.feb')
subprocess.run([FEBio_path, FEBio_inputfile], check = True)


#Deletes all the temporary files in the temp folder (have to fix)
""" temp_files = glob(temp_dir)
for f in temp_files:
    os.remove(f) """

print('Done!')

