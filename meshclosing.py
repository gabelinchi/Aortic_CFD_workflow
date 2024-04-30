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

    Known bugs: if plotting is set to True and coords_to_remove is empty, the function fails
    '''
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
        rows_with_x = np.any(faces == point, axis=1)

        # Extract neighbouring points
        edge_points = []
        for row in faces[rows_with_x]:
            edge_points.extend(row[row != point]) 
        edge_points = np.array(list(dict.fromkeys(edge_points)))

        # Delete faces containing point
        faces = faces[~rows_with_x]

        # Create a mesh to fill the gap and add to fill_mesh
        fill_mesh = fill_mesh.merge(pv.PolyData(points[edge_points]).delaunay_2d())

    # Rebuild polydata of mesh to fill (now with the gaps)
    mesh_to_fill = pv.PolyData(points, pad(faces))

    # Combine with the patches
    reduced_surface = fill_mesh.merge(mesh_to_fill).clean()
    if plot==True:
        fill_mesh.plot(show_edges=True)
        mesh_to_fill.plot(show_edges=True)
        reduced_surface.plot(show_edges=True)
    return(reduced_surface) # Return pv.PolyData of the reconstructed surface


wall_reduced = remove_points_and_fill(wall, point_selection(inlet, wall, adjustment, tolerance)[1], plot=True)
inlet_reduced = remove_points_and_fill(inlet, point_selection(inlet, wall, adjustment, tolerance)[0])

wall_reduced.merge(inlet_reduced).clean(tolerance=0.1).plot(show_edges=True)