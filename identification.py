import numpy as np
import pyvista as pv

""" # Read 3D mesh
tetmesh = pv.read('aorta_tetmesh.vtk') """

#Function that removes points outside the target area for selection 
def remove_far_points(selection, center):

    #Create distance storage variable
    distances = np.array([])

    #Calculate the distances from the center of the inlet/outlet to every surface node
    for i in range(len(selection)):
        distance = np.linalg.norm(selection[i] - center)
        distances = np.append(distances, distance)
    
    #Calculate the average distance of all the distances
    avg_distance = np.average(distances)

    #Create new storage area for all nodes that belong to the inlet
    new_selection = np.empty((0,3))

    #Select all nodes that are NOT on a large distance from the center
    for i in range(len(distances)):
        if distances[i] < 4 * avg_distance:
            new_selection = np.vstack([new_selection, selection[i]])
        
    return new_selection

#Identificate inlet/outlet
def selection(wall, extraction_surface):

    print('Start inlet/outlet identification')
    #Extract the nodes on the surface of the 3D_mesh
    wall_surface = wall.extract_surface()
    wall_points = wall_surface.points

    #Create storage variable for nodes after cutting the surface
    first_selection = np.empty((0,3))

    #Calculate center and normal of the extraction_surface
    center = extraction_surface.points.mean(0)
    normal = extraction_surface.compute_normals()['Normals'].mean(0)

    #Calculate if nodes are in the plane of the cutting surface, if they are, they are added to first_selection
    for i in range(len(wall_points)):
        if -0.01 < np.dot(np.subtract(wall_points[i], center), normal) < 0.01:
            first_selection = np.vstack([first_selection, wall_points[i]])

    #Remove the points that are too far away from the center to belong to the inlet/outlet
    second_selection = remove_far_points(first_selection, center)

    #Create a closed surface mesh from the selected points
    second_selection = pv.PolyData(second_selection)
    second_selection = second_selection.delaunay_2d()

    print('Identification succesfull')

    return second_selection



""" # Extract surface
surface = tetmesh.extract_surface()

# Extract edges
edges = surface.extract_feature_edges(50)
edges.plot()

# Delete edges from surface
split_surface = surface.remove_cells(edges.CellData["vtkOriginalCellIds"])
split_surface.plot(show_edges=True) """
