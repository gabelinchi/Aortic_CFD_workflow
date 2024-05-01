import pyvista as pv
import numpy as np
import meshio


'''Function that selects the points on the boundary edge that don't have a neighbouring to connect with .clean. These points will be
removed and remeshed later in the workflow. Can be used for inlet as well as outlet'''
def point_selection(inlet_outlet, wall, adjustment, tolerance, config):
    
    print('Start point selection')

    #Selects clipping axis dependent if an input or output is inserted in the function
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

    print('Finished point selection')
    return i_o_excess, wall_excess


# Function to pad an unpadded face array
def pad(faces):
    num_rows = faces.shape[0]
    return(np.hstack((np.full((num_rows, 1), 3),faces)))

def remove_points_and_fill(polydata, coords_to_remove, plot=False):
    '''
    Function that removes a point from a pyvistas polydata object and remeshes the hole left behind

    polydata: pyvista.PolyData object
    coords_to_remove: numpy array (3 columns) with coords of points to be removed
    plot: boolean, determines if process stepps are plotted, false by default

    Known bugs: (1) if plotting is set to True and coords_to_remove is empty, the function fails
                (2) if two neighbouring points are to be deleted, the resulting edgde contains a gap (in progress)
    '''
    print('Start point removal and remesh')
    # Convert input to polydata
    polydata = polydata.extract_surface()

    # Split polydata into np arrays
    faces = polydata.faces.reshape(-1, 4)[:,[1,2,3]]
    points = polydata.points

    # Obtain point indices that match input coords
    matches = (points == coords_to_remove[:, np.newaxis]).all(axis=2)
    points_to_remove = np.where(matches)[1]

    # Loop over points to populate fill_mesh with the pathches
    fill_mesh = pv.PolyData()
    for point in points_to_remove:
        # Find faces containing point
        rows_with_point = np.any(faces == point, axis=1)

        # Extract neighbouring points
        edge_points = []
        for row in faces[rows_with_point]:
            edge_points.extend(row[row != point]) 
        edge_points = np.array(list(dict.fromkeys(edge_points)))

        # Check if any of the neighbouring points is a point to be deleted (fix for bug 2)
        duplicates_bool = np.in1d(points_to_remove, edge_points)
        duplicates_list = points_to_remove[duplicates_bool]
        if len(duplicates_list) == 1: 
            print('Neighbouring point detected to point to be deleted:', duplicates_list)
            neighbour = duplicates_list[0]
            # Execute regular routine for the bordering point as well
            #--------------------------------------------------------
            # Find faces containing point
            rows_with_neighbour = np.any(faces == neighbour, axis=1)

            # Extract neighbouring points
            edge_points_neighbour = []
            for row in faces[rows_with_neighbour]:
                edge_points_neighbour.extend(row[row != neighbour]) 
            edge_points_neighbour = np.array(list(dict.fromkeys(edge_points_neighbour)))
            #--------------------------------------------------------
            #Add new faces and edge points together
            edge_points = np.array(list(dict.fromkeys(np.hstack((edge_points, edge_points_neighbour)))))
            rows_with_point = np.any([rows_with_point, rows_with_neighbour], axis=0)

            # Remove neighbour and orignal point from edge_points (ensures the touching edge remains straight)
            edge_points = np.delete(edge_points, np.where((edge_points == neighbour) | (edge_points == point)))

        # Throw an error if too many points in a row are to be removed
        elif len(duplicates_list) > 1:
            raise Exception("More than two neighbouring points to be deleted: mesh mismatch too large")
            

        # Delete faces containing point
        faces = faces[~rows_with_point]

        # Create a mesh to fill the gap and add to fill_mesh
        if len(edge_points) > 2: #Checks if there are points to mesh
            next_fill_mesh = pv.PolyData(points[edge_points]).delaunay_2d()
            if plot == True:
                next_fill_mesh.plot(show_edges=True)
            fill_mesh = fill_mesh.merge(next_fill_mesh)

    # Rebuild polydata of mesh to fill (now with the gaps)
    mesh_to_fill = pv.PolyData(points, pad(faces))

    # Combine with the patches
    reduced_surface = fill_mesh.merge(mesh_to_fill).clean()
    if plot==True:
        mesh_to_fill.plot(show_edges=True)
        reduced_surface.plot(show_edges=True)
    
    print('Finished point removal and remesh')
    return(reduced_surface) # Return pv.PolyData of the reconstructed surface


""" inlet_point_selection, wall_point_selection = point_selection(inlet, wall, adjustment, tolerance, 'inlet')
inlet_reduced = remove_points_and_fill(inlet, inlet_point_selection, plot=False)
wall_reduced = remove_points_and_fill(wall, wall_point_selection, plot=True) """

#wall_reduced.merge(inlet_reduced).clean(tolerance=0.1).plot(show_edges=True)