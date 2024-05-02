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

'''
The goal of this script is to calculate quality metrics for surface meshes, to begin I will use the wall
as a comparison point.
'''

# Remesh with mmg
inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"
outlet_path = "geometries\input\outlet.stl"

inlet_remeshed_mmg, wall_remeshed_mmg, outlet_remeshed_mmg = remesh.remesh(inlet_path, wall_path, outlet_path)

# Remesh with trimesh