import pyvista as pv
import numpy as np
import meshio


adjustment = np.array([0,0,20])
tolerance = 0.1 #Need to check this value and what it means. 0.1 works with the expected output nodes

#prep .mesh for modification with pyvista
meshio.write('inlet_mmg.vtk', meshio.read('inlet_mmg.mesh'))
meshio.write('wall_mmg.vtk', meshio.read('wall_mmg.mesh'))

#import remeshed geometry
inlet = pv.read('inlet_mmg.vtk').extract_surface()
wall = pv.read('wall_mmg.vtk').extract_surface()


'''Function that selects the points on the boundary edge that don't have a neighbouring to connect with .clean. These points will be
removed and remeshed later in the workflow. Can be used for inlet as well as outlet'''
def point_selection(inlet_outlet, wall, adjustment, tolerance, config):
    
    #Selects clipping axis dependent if an input or output is inputted
    if config == 'inlet':
        clip_axis = '-z'
    else:
        clip_axis = 'z'

    #Select edges
    i_o_edges = inlet_outlet.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
    wall_edges = wall.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
    wall_i_o_edges = wall_edges.clip(clip_axis, origin=(wall.center - adjustment))

    #Get points
    i_o_points = i_o_edges.points
    wall_points = wall_i_o_edges.points


    #Calculate the amount of rows (points) within the matrix
    n_rows_i_o = i_o_points.shape[0]
    n_rows_wall = wall_points.shape[0]

    #Create storage variables
    i_o_excess = np.empty((0,3))
    wall_excess = np.empty((0,3))
    distance_inlet = np.zeros(n_rows_wall)
    distance_wall = np.zeros(n_rows_i_o)


    #Check for every input/output point the distance to all the other wall points. If the inlet/outlet point doesn't have a neighbouring wall point
    #within tolerance, the point is added to i_o_excess. The same operation is done for the wall points.
    if n_rows_i_o != n_rows_wall:
        for i in range(n_rows_i_o):
            for k in range(n_rows_wall):
                distance_inlet[k] = np.linalg.norm(i_o_points[i] - wall_points[k])
            if not any(distance <= tolerance for distance in distance_inlet):
                i_o_excess = np.vstack([i_o_excess, i_o_points[i]])
        
        for n in range(n_rows_wall):
            for t in range(n_rows_i_o):
                distance_wall[t] = np.linalg.norm(wall_points[n] - i_o_points[t])
            if not any(distance <= tolerance for distance in distance_wall):
                wall_excess = np.vstack([wall_excess, wall_points[n]])

    return i_o_excess, wall_excess


wall_excess = point_selection(inlet, wall, adjustment, tolerance, 'inlet')[1]

print(wall_excess)

        


