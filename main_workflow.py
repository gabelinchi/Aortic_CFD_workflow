
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
import volume_remesh_mmg
from capping import cap
import identification as id
from tkinter import Tk
from tkinter.filedialog import askdirectory
import mapping
import xml.etree.ElementTree as ET
import subprocess
import febioxml as feb
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
    'mesh_density': '0.1',
    'sizing': '1',
    'detection angle': '35'}

mmg3d_parameters = {
    'hausd': '0.1',
    'max_edgelength': '1',
    'detection angle': '35'}

tetgen_parameters = dict(
    nobisect=True)

#Angle for identification
id_angle = 30

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

#Make an initial 3D mesh from the combined mesh using TetGen
combined_remeshed = combined_remeshed
tetmesh = volume_mesh.volume_meshing(combined_remeshed, tetgen_parameters, plot=show_plot)

#Save initial mesh
pv.save_meshio(osp.join(temp_dir, r'initial_volume_mesh.mesh'), tetmesh)

#Refine 3D mesh with mmg3d
tetmesh = volume_remesh_mmg.mmg3d(osp.join(temp_dir, r'initial_volume_mesh.mesh'), osp.join(temp_dir, r'mmg3d_mesh.mesh'), temp_dir, mmg3d_parameters, plot=show_plot)

# Calculate mesh stats and quality
jac = tetmesh.compute_cell_quality(quality_measure='scaled_jacobian')['CellQuality']
aspect = tetmesh.compute_cell_quality(quality_measure='aspect_ratio')['CellQuality']
num_points = tetmesh.number_of_points
num_cells = tetmesh.number_of_cells

# Print report
print('--------------------------------------------------------------')
print('3D mesh report:')
print('num points:', num_points)
print('num_cells: ', num_cells)
print('3D mesh quality (mean scaled jacobian)')
print('min:', jac.min())
print('max:', jac.max())
print('avg:', jac.mean())
print('3D mesh quality (aspect ratio)')
print('min:', aspect.min())
print('max:', aspect.max())
print('avg:', aspect.mean())
print('--------------------------------------------------------------')

#Plot bad cells
if show_plot:
    plt = pv.Plotter()
    plt.add_mesh(tetmesh, style='wireframe')
    plt.add_mesh(tetmesh.extract_cells(jac<0.3), color='red', show_edges=True)
    plt.add_text('Bad cells')
    plt.show()


#Plot bisection------------------------------
if show_plot:
    # get cell centroids

    cells = tetmesh.cells.reshape(-1, 5)[:, 1:]
    cell_center = tetmesh.points[cells].mean(1)

    # extract cells below the 0 xy plane
    mask = cell_center[:, 2] < 0
    cell_ind = mask.nonzero()[0]
    subgrid = tetmesh.extract_cells(cell_ind)

    # advanced plotting
    plotter = pv.Plotter()
    plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
    plotter.add_mesh(combined_remeshed, 'r', 'wireframe')
    plotter.add_legend([[' Input Mesh ', 'r'], [' Tessellated Mesh ', 'black']])
    plotter.show()
#--------------------------------------------

#Save 3D mesh
tetmesh.save(osp.join(temp_dir, r'3D_output_mesh.vtk'))

#Plot 3D_mesh
if show_plot:
    tetmesh.plot(show_edges = True)

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
#Plot the identifies surfaces for general overview
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

velocity_map, n_maps = mapping.vel_mapping(vel_profile_dir, id_inlet, output_dir, intp_options, show_plot)

#----------------------------------------------------------------------------------------------------------------------------
# FEBio
#----------------------------------------------------------------------------------------------------------------------------

#Create a solver compatible file based on the 3D-mesh and meshing parameters
feb.xml_creator(tetmesh, id_inlet, id_outlet, id_wall, file_dir, output_dir)

#Run FEBio
FEBio_path = r"C:/Program Files/FEBioStudio2/bin/febio4.exe"
#Use the current
FEBio_inputfile = osp.join(temp_dir, r'simulation.feb')
subprocess.run([FEBio_path, FEBio_inputfile], check = True)


#Deletes all the temporary files in the temp folder (have to fix)
""" temp_files = glob(temp_dir)
for f in temp_files:
    os.remove(f) """

print('Done!')

