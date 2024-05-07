
#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import remesh
import meshclosing
import utils as ut
import cutting
import volume_mesh

#paths to the aorta geometry

inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"
outlet_path = "geometries\input\outlet.stl"

#Filename of cutted wall
fn_wall_cut = "wall_cut.mesh"

#Reading the files with pyvista
inlet = pv.read(inlet_path)
wall = pv.read(wall_path)
outlet = pv.read(outlet_path)

#Meshing parameters
adjustment = np.array([0,0,20])
search_tolerance = 0.1 #Need to check this value and what it means. 0.1 works with the expected output nodes

mmg_parameters = {
    'mesh_density': '0.1',
    'sizing': '1'}

tetgen_parameters = dict(
    order=1, 
    mindihedral=20, 
    minratio=1.5)

#Cut the wall geometry after the aortic root
wall_cut = cutting.main_cutter(inlet, wall, plot=True)
pv.save_meshio(fn_wall_cut, wall_cut)

#Run remesh (takes the file Path !!! as input)
inlet_remeshed, wall_remeshed, outlet_remeshed = remesh.remesh(inlet_path, fn_wall_cut, outlet_path, mmg_parameters, plot=True)

#Close remeshed meshes

#Extracting the excess points
outlet_point_selection, wall_point_selection = meshclosing.point_selection(outlet_remeshed, wall_remeshed, adjustment, search_tolerance, 'outlet')
inlet_point_selection, wall_point_selection_final = meshclosing.point_selection(inlet_remeshed, wall_remeshed, adjustment, search_tolerance, 'inlet')

#Delete the excess points and perform a local remesh
wall_reduced = meshclosing.remove_points_and_fill(wall_remeshed, np.vstack((wall_point_selection, wall_point_selection_final)))
outlet_reduced = meshclosing.remove_points_and_fill(outlet_remeshed, outlet_point_selection)
inlet_reduced = meshclosing.remove_points_and_fill(inlet_remeshed, inlet_point_selection)

#Combine the closed meshes into one geometry
combined_mesh = (wall_reduced + outlet_reduced+inlet_reduced).clean(tolerance=search_tolerance)
combined_mesh.plot(show_edges=True)

#Make a 3D mesh from the combined mesh
tetmesh = volume_mesh.volume_meshing(combined_mesh, tetgen_parameters, False)

tetmesh.plot(show_edges = True)

#print(mesh)

#use an excessive amount of comments on everything.
#extra line of code changed
