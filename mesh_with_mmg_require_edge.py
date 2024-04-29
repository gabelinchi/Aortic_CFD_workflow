#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import meshio
import subprocess as sub

sub.run('py -m mmgs -hausd 0.1 inlet_handedit.mesh inlet_handedit_mmg.mesh')
meshio.write('inlet_handedit_mmg.vtk', meshio.read('inlet_handedit_mmg.mesh'))
pv.read('inlet_handedit_mmg.vtk')

""" # Convert input file from stl to .mesh
meshio.write('inlet.mesh', meshio.read('geometries/input/inlet.stl'))
meshio.write('wall.mesh', meshio.read('geometries/input/wall.stl'))

#import .stl's to pyvista
pv_wall_mesh = pv.read("geometries\input\wall.stl")
pv_inlet_mesh = pv.read("geometries\input\inlet.stl")

# Extract boundaries
inlet_edge = pv_inlet_mesh.extract_feature_edges(manifold_edges=False, non_manifold_edges=False, feature_edges = False, boundary_edges=True)
inlet_edge.plot(show_edges=True)

# Save boundaries as .mesh
inlet_edge.save('geometries/temp/inlet_edge.ply')
meshio.write('geometries/temp/inlet_edge.mesh', meshio.read('geometries/temp/inlet_edge.ply'))
pv.read('geometries/temp/inlet_edge.ply').plot(show_edges=True) #plot inlet vertices

# Run mmg
sub.run('py -m mmgs -hausd 0.1 inlet.mesh inlet_mmg.mesh -hsiz 1')
sub.run('py -m mmgs -hausd 0.1 wall.mesh wall_mmg.mesh -hsiz 1')

# Convert back to .vtk and plot with pyvista
meshio.write('inlet_mmg.vtk', meshio.read('inlet_mmg.mesh'))
meshio.write('wall_mmg.vtk', meshio.read('wall_mmg.mesh'))
inlet_remeshed = pv.read('inlet_mmg.vtk')
wall_remeshed = pv.read('wall_mmg.vtk')
combined = inlet_remeshed + wall_remeshed
combined.plot(show_edges = True)

# Merge and plot with pv
wall_and_inlet = wall_remeshed.merge(inlet_remeshed).clean(tolerance=0.001) #combines two meshes and removes duplicate points
wall_and_inlet.plot(show_edges=True) """