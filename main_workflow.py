
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

#paths to the aorta geometry

inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"
outlet_path = "geometries\input\outlet.stl"

#Meshing parameters
adjustment = np.array([0,0,20])
search_tolerance = 0.1 #Need to check this value and what it means. 0.1 works with the expected output nodes

mmg_parameters = {
    'mesh_density': '0.1',
    'sizing': '1'}


#Run remesh
inlet_remeshed, wall_remeshed, outlet_remeshed = remesh.remesh(inlet_path, wall_path, outlet_path, mmg_parameters)

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

# create 3D tetmesh from surface mesh
tetmesh = tet.TetGen(combined_mesh)
tetmesh.tetrahedralize(order=1, mindihedral=20, minratio=1.5)
grid = tetmesh.grid
grid.plot(show_edges=True)

# Save the generated mesh as a .stl
grid.save('aorta_tetmesh.vtk', binary=False)

# get cell centroids
cells = grid.cells.reshape(-1, 5)[:, 1:]
cell_center = grid.points[cells].mean(1)

# extract cells below the 0 xy plane
mask = cell_center[:, 2] < 0
cell_ind = mask.nonzero()[0]
subgrid = grid.extract_cells(cell_ind)

# advanced plotting
plotter = pv.Plotter()
plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
plotter.add_mesh(combined_mesh, 'r', 'wireframe')
plotter.add_legend([[' Input Mesh ', 'r'],
                    [' Tessellated Mesh ', 'black']])
plotter.show()

#print(mesh)

#use an excessive amount of comments on everything.
#extra line of code changed
