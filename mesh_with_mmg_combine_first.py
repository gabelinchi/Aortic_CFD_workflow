#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import meshio
import subprocess as sub


#Combine .stl files with pv
wall_mesh = pv.read("geometries\input\wall.stl")
outlet_mesh = pv.read("geometries\input\outlet.stl")
inlet_mesh = pv.read("geometries\input\inlet.stl")

wall_and_inlet = wall_mesh.merge(inlet_mesh).clean(polys_to_lines=False, strips_to_polys=False) #combines two meshes and removes duplicate points
wall_and_io = wall_and_inlet.merge(outlet_mesh).clean(polys_to_lines=False, strips_to_polys=False)
wall_and_io.plot(show_edges=True)
wall_and_io.save('geometries/temp/wall_and_io.stl')

# Convert input file from .vtk to .mesh
meshio.write("geometries/temp/wall_and_io.mesh", meshio.read("geometries/temp/wall_and_io.stl"))

# Run mmg
sub.run('py -m mmgs -ar 30 -hausd 0.1 -hsiz geometries/temp/wall_and_io.mesh geometries/remeshed/wall_and_io_mmg.mesh')


# Convert back to .vtk and plot with pyvista
meshio.write('geometries/temp/wall_and_io_mmg.vtk', meshio.read('geometries/remeshed/wall_and_io_mmg.mesh'))

# Plot with pv
pv.read('geometries/temp/wall_and_io_mmg.vtk').plot(show_edges = True)