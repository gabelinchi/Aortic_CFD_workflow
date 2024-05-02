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
import trimesh as tri

'''
The goal of this script is to calculate quality metrics for surface meshes, to begin I will use the wall
as a comparison point.
'''

# Import geometry
inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"
outlet_path = "geometries\input\outlet.stl"

wall_mesh = pv.read(wall_path)
outlet_mesh = pv.read(outlet_path)
inlet_mesh = pv.read(inlet_path)


# Remesh with mmg
"""
inlet_remeshed_mmg, wall_remeshed_mmg, outlet_remeshed_mmg = remesh.remesh(inlet_path, wall_path, outlet_path) """


# Remesh with trimesh
def extract_nodes_and_faces(unstrucuredgrid):
    """
    Converts a pyvista unstructured grid into an array of faces and an array of nodes
    Written by Yarran
    """
    nodes = unstrucuredgrid.points
    poly = unstrucuredgrid.extract_surface()
    faces = np.asarray(poly.faces).reshape((-1, 4))[:, 1:]
    return(nodes, faces)

inlet_nodes, inlet_faces = extract_nodes_and_faces(inlet_mesh)
outlet_nodes, outlet_faces = extract_nodes_and_faces(outlet_mesh)
wall_nodes, wall_faces = extract_nodes_and_faces(wall_mesh)

inlet_from_arrays = pv.PolyData.from_regular_faces(inlet_nodes, inlet_faces)
outlet_from_arrays = pv.PolyData.from_regular_faces(outlet_nodes, outlet_faces)
wall_from_arrays = pv.PolyData.from_regular_faces(wall_nodes, wall_faces)

plotter = pv.Plotter()
plotter.add_mesh(inlet_from_arrays, show_edges=True)
plotter.add_mesh(wall_from_arrays, show_edges=True)
plotter.show()

inlet_remeshed_nodes, inlet_remeshed_faces = tri.remesh.subdivide_to_size(inlet_nodes, inlet_faces, 2)
outlet_remeshed_nodes, outlet_remeshed_faces = tri.remesh.subdivide_to_size(outlet_nodes, outlet_faces, 2)
wall_remeshed_nodes, wall_remeshed_faces = tri.remesh.subdivide_to_size(wall_nodes, wall_faces, 2)
inlet_remeshed_trimesh = pv.PolyData.from_regular_faces(inlet_remeshed_nodes, inlet_remeshed_faces)
outlet_remeshed_trimesh = pv.PolyData.from_regular_faces(outlet_remeshed_nodes, outlet_remeshed_faces)
wall_remeshed_trimesh = pv.PolyData.from_regular_faces(wall_remeshed_nodes, wall_remeshed_faces)

plotter = pv.Plotter()
plotter.add_mesh(inlet_remeshed_trimesh, show_edges=True)
plotter.add_mesh(outlet_remeshed_trimesh, show_edges=True)
plotter.add_mesh(wall_remeshed_trimesh, show_edges=True)
plotter.show()

# Remesh with pyacvd

