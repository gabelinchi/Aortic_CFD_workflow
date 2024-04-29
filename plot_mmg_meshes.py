#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import meshio
import subprocess as sub

# Convert input file from stl to .mesh
meshio.write('inlet.mesh', meshio.read('geometries/input/inlet.stl'))
meshio.write('wall.mesh', meshio.read('geometries/input/wall.stl'))

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