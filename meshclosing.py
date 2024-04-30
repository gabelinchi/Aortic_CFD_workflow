import pyvista as pv
import numpy as np
import meshio

adjustment = np.array([0,0,20])
tolerance = 0.1 #Need to check this value and what it means

#prep .mesh for modification with pyvista
meshio.write('inlet_mmg.vtk', meshio.read('inlet_mmg.mesh'))
meshio.write('wall_mmg.vtk', meshio.read('wall_mmg.mesh'))

#import remeshed geometry
inlet = pv.read('inlet_mmg.vtk')
wall = pv.read('wall_mmg.vtk')

#Select edges
inlet_edges = inlet.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
wall_edges = wall.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
wall_inlet_edges = wall_edges.clip('-z', origin=(wall.center - adjustment))

#Get points
inlet_points = inlet_edges.points
wall_points = wall_inlet_edges.points


#Calculate the amount of rows (points) within the matrix
n_rows_inlet = inlet_points.shape[0]
n_rows_wall = wall_points.shape[0]

print(n_rows_inlet)
print(n_rows_wall)

#Create storage variables
inlet_excess_index = np.array([])
wall_excess_index = np.array([])
distance_inlet = np.zeros(n_rows_wall)
distance_wall = np.zeros(n_rows_inlet)


#Check for every inlet point the distance all the other wall points. If the inlet point doesn't have a neighbouring wall point
#within tolerance, the index of the point is added to inlet_excess_index. The same operation is done for the wall points.
if n_rows_inlet != n_rows_wall:
    for i in range(n_rows_inlet):
        for k in range(n_rows_wall):
            distance_inlet[k] = np.linalg.norm(inlet_points[i] - wall_points[k])
        print(distance_inlet)
        if not any(distance <= tolerance for distance in distance_inlet):
            inlet_excess_index = np.append(inlet_excess_index, i)
    
    for n in range(n_rows_wall):
        for t in range(n_rows_inlet):
            distance_wall[t] = np.linalg.norm(wall_points[n] - inlet_points[t])
        if not any(distance <= tolerance for distance in distance_wall):
            wall_excess_index = np.append(wall_excess_index, n)

print(inlet_excess_index)
print(wall_excess_index)
        


