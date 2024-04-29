import pyvista as pv
import numpy as np

inlet_remeshed = pv.read('inlet_mmg.vtk').extract_surface()
wall_remeshed = pv.read('wall_mmg.vtk').extract_surface()
#inlet_remeshed.plot(show_edges=True)
print(inlet_remeshed.points)
print(inlet_remeshed.faces.reshape(-1, 4)[:,[1,2,3]])
faces = inlet_remeshed.faces.reshape(-1, 4)[:,[1,2,3]]

point_to_remove = 1

#find faces containing point
rows_with_x = np.any(faces == point_to_remove, axis=1)

#extract neighbouring points
edge_points = []
for row in faces[rows_with_x]:
    edge_points.extend(row[row != point_to_remove]) 
edge_points = list(dict.fromkeys(edge_points))
print(edge_points)

#delete faces containing the point
faces_reduced = faces[~rows_with_x]

#next step is to create faces to connect the edge points, add those to faces_reduced, and rebuild the Polydata
new_faces = np.empty(3)
edge_points = np.array(edge_points)
for i in range(len(edge_points) - 1):
    print(i)
    print(edge_points)
    new_faces = np.vstack([new_faces, edge_points.take([i, i+1, i+2], mode='wrap')])
print(new_faces)



""" points_to_keep = ~np.isin(np.arange(inlet_remeshed.n_points), points_to_remove)
inlet_reduced = inlet_remeshed.extract_points(points_to_keep, adjacent_cells=False)
#inlet_reduced.plot(show_edges=True) """