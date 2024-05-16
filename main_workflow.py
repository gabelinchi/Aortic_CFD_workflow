
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

tetgen_parameters = dict(
    nobisect=True,
    fixedvolume=True,
    quality=True,
    maxvolume=1,
    order=1, 
    mindihedral=35, 
    minratio=2)

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

#Make a 3D mesh from the combined mesh
combined_remeshed = combined_remeshed.delaunay_3d()
pv.save_meshio('combined_remeshed.mesh', combined_remeshed)

""" tetmesh = volume_mesh.volume_meshing(combined_remeshed, tetgen_parameters, plot=show_plot) """

import subprocess as sub
import meshio
#sub.run(f"{'py -m mmg3d -hausd 1 combined_remeshed.mesh combined_remeshed_mmg.mesh'}")
sub.run(f"{'py -m mmg3d -hausd 0.001 -hgrad 1.3 -hmax 0.05 combined_remeshed.mesh combined_remeshed_mmg.mesh'}")
meshio.write('combined_remeshed_mmg.vtk', meshio.read('combined_remeshed_mmg.mesh'))
tetmesh = pv.read('combined_remeshed_mmg.vtk')


#Report the quality of the mesh
jac = tetmesh.compute_cell_quality(quality_measure='scaled_jacobian')['CellQuality']
aspect = tetmesh.compute_cell_quality(quality_measure='aspect_ratio')['CellQuality']
print('3D mesh quality (mean scaled jacobian)')
print('min:', jac.min())
print('max:', jac.max())
print('avg:', jac.mean())
print('3D mesh quality (aspect ratio)')
print('min:', aspect.min())
print('max:', aspect.max())
print('avg:', aspect.mean())

#Plot bad cells
plt = pv.Plotter()
plt.add_mesh(tetmesh, style='wireframe')
plt.add_mesh(tetmesh.extract_cells(jac<0.1), color='red', show_edges=True)
plt.add_text('Bad cells')
plt.show()

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


#Deletes all the temporary files in the temp folder (have to fix)
""" temp_files = glob(temp_dir)
for f in temp_files:
    os.remove(f) """

print('Done!')

