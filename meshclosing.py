import pyvista as pv
import numpy as np
import meshio

meshio.write('inlet_mmg.vtk', meshio.read('inlet_mmg.mesh'))
meshio.write('wall_mmg.vtk', meshio.read('wall_mmg.mesh'))

inlet = pv.read('inlet_mmg.vtk')
wall = pv.read('wall_mmg.vtk')

inlet_edges = inlet.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
wall_edges = wall.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
