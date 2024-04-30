import pyvista as pv
import numpy as np
import meshio

""" 
adjustment = np.array([0,0,20])
tolerance = 0.1 #Need to check this value and what it means. 0.1 works with the expected output nodes

#prep .mesh for modification with pyvista
meshio.write('inlet_mmg.vtk', meshio.read('inlet_mmg.mesh'))
meshio.write('wall_mmg.vtk', meshio.read('wall_mmg.mesh'))

#import remeshed geometry
inlet = pv.read('inlet_mmg.vtk').extract_surface()
wall = pv.read('wall_mmg.vtk').extract_surface()
 """

'''Function that selects the points on the boundary edge that don't have a neighbouring to connect with .clean. These points will be
removed and remeshed later in the workflow. Can be used for inlet as well as outlet (do still need to adjust code for that)'''
def point_selection(inlet, wall, adjustment, tolerance):
    
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

    #Create storage variables
    inlet_excess = np.empty((0,3))
    wall_excess = np.empty((0,3))
    distance_inlet = np.zeros(n_rows_wall)
    distance_wall = np.zeros(n_rows_inlet)


    #Check for every inlet point the distance all the other wall points. If the inlet point doesn't have a neighbouring wall point
    #within tolerance, the index of the point is added to inlet_excess_index. The same operation is done for the wall points.
    if n_rows_inlet != n_rows_wall:
        for i in range(n_rows_inlet):
            for k in range(n_rows_wall):
                distance_inlet[k] = np.linalg.norm(inlet_points[i] - wall_points[k])
            if not any(distance <= tolerance for distance in distance_inlet):
                inlet_excess = np.vstack([inlet_excess, inlet_points[i]])
        
        for n in range(n_rows_wall):
            for t in range(n_rows_inlet):
                distance_wall[t] = np.linalg.norm(wall_points[n] - inlet_points[t])
            if not any(distance <= tolerance for distance in distance_wall):
                wall_excess = np.vstack([wall_excess, wall_points[n]])

    return inlet_excess, wall_excess

inlet_excess, wall_excess = point_selection(inlet, wall, adjustment, tolerance)



        


