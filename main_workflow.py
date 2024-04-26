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
import pyacvd as acvd

wall_mesh = pv.read("wall.stl")
outlet_mesh = pv.read("outlet.stl")
inlet_mesh = pv.read("inlet.stl")


def ac_remesh(mesh, subdivide, cluster, plots):
    """
    remeshes a surface using pyacvd
    """
    if cluster == 0:
        cluster = mesh.n_points * 2
    if plots == True:
        mesh.plot(show_edges=True)
    clus = acvd.Clustering(mesh)
    clus.subdivide(subdivide)
    clus.cluster(cluster)
    if plots == True:
        clus.plot()
    remesh = clus.create_mesh().clean(lines_to_points=False, polys_to_lines=False)
    if plots == True:
        remesh.plot(show_edges=True)
        remesh.plot_normals()
    return(remesh)

#Combine first
""" wall_and_inlet = wall_mesh.merge(inlet_mesh).clean(polys_to_lines=False, strips_to_polys=False) #combines two meshes and removes duplicate points
wall_and_io = wall_and_inlet.merge(outlet_mesh).clean(polys_to_lines=False, strips_to_polys=False)
wall_and_io.plot(show_edges=True)
wall_and_io_remesh = ac_remesh(wall_and_io, 3, 0, plots=True) """



#mesh first, combine later

inlet_remesh = ac_remesh(inlet_mesh, 4, 0, plots=True)
wall_remesh = ac_remesh(wall_mesh, 3, 0, plots=True)
wall = wall_mesh.plot()
wall_and_inlet = wall_remesh.merge(inlet_remesh).clean(polys_to_lines=False, strips_to_polys=False) #combines two meshes and removes duplicate points
wall_and_io = wall_and_inlet.merge(outlet_mesh).clean(polys_to_lines=False, strips_to_polys=False)
wall_and_io.plot(show_edges=True)

""" # create 3D tetmesh from surface mesh
tetmesh = tet.TetGen(wall_and_io_remesh)
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
plotter.add_mesh(wall_and_io_remesh, 'r', 'wireframe')
plotter.add_legend([[' Input Mesh ', 'r'],
                    [' Tessellated Mesh ', 'black']])
plotter.show() """

#print(mesh)

#use an excessive amount of comments on everything.
#extra line of code changed
