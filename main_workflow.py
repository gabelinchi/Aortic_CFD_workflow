# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:16:40 2024

@author: lmorr
"""

#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import trimesh as tri

def extract_nodes_and_faces(unstrucuredgrid):
    """
    Converts a pyvista unstructured grid into an array of faces and an array of nodes
    Written by Yarran
    """
    nodes = unstrucuredgrid.points
    poly = inlet_mesh.extract_surface()
    faces = np.asarray(poly.faces).reshape((-1, 4))[:, 1:]
    return(nodes, faces)

wall_mesh = pv.read("wall.stl")
outlet_mesh = pv.read("outlet.stl")
inlet_mesh = pv.read("inlet.stl")
#wall = wall_mesh.plot()

inlet_nodes, inlet_faces = extract_nodes_and_faces(inlet_mesh)
#print(inlet_nodes)
#print(inlet_faces)

inlet_from_arrays = pv.PolyData.from_regular_faces(inlet_nodes, inlet_faces)
print(inlet_from_arrays)
plotter = pv.Plotter()
plotter.add_mesh(inlet_from_arrays, show_edges=True)
plotter.show()

inlet_remeshed_nodes, inlet_remeshed_faces = tri.remesh.subdivide_to_size(inlet_nodes, inlet_faces, 1)
inlet_remeshed = pv.PolyData.from_regular_faces(inlet_remeshed_nodes, inlet_remeshed_faces)
plotter = pv.Plotter()
plotter.add_mesh(inlet_remeshed, show_edges=True)
plotter.show()

""" # create 3D tetmesh from surface mesh
tetmesh = tet.TetGen(wall_and_io)
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
plotter.add_mesh(wall_and_io, 'r', 'wireframe')
plotter.add_legend([[' Input Mesh ', 'r'],
                    [' Tessellated Mesh ', 'black']])
plotter.show()

#print(mesh) """

#use an excessive amount of comments on everything.
#extra line of code changed
