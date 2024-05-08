import numpy as np
import pyvista as pv

# Read 3D mesh
tetmesh = pv.read('aorta_tetmesh.vtk')

# Extract surface
surface = tetmesh.extract_surface()

# Extract edges
edges = surface.extract_feature_edges(50)
edges.plot()

# Delete edges from surface
split_surface = surface.remove_cells(edges.CellData["vtkOriginalCellIds"])
split_surface.plot(show_edges=True)
