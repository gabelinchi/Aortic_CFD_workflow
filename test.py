import pyvista as pv
import numpy as np
import meshio 
inlet_mesh = pv.read("inlet.stl")

inlet = pv.plot (inlet_mesh)
"""
Converts a pyvista unstructured grid into an array of faces and an array of nodes   """
nodes = inlet_mesh.points
poly = inlet_mesh.extract_surface()
faces = np.asarray(poly.faces).reshape((-1, 4))[:, 1:]
meshio.write_points_cells("inlet.stl", nodes, faces)



