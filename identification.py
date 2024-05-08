import numpy as np
import pyvista as pv

""" # Read 3D mesh
tetmesh = pv.read('aorta_tetmesh.vtk') """

def remove_far_points(selection, center):
    distances = np.array([])
    for i in range(len(selection)):
        distance = np.linalg.norm(selection[i] - center)
        distances = np.append(distances, distance)
    
    avg_distance = np.average(distances)
    new_selection = np.empty((0,3))

    for i in range(len(distances)):
        if distances[i] < 4 * avg_distance:
            new_selection = np.vstack([new_selection, selection[i]])
        

    return new_selection

#Identificate inlet/outlet
def selection(wall, extraction_surface):
    wall_surface = wall.extract_surface()
    wall_points = wall_surface.points
    first_selection = np.empty((0,3))
    center = extraction_surface.points.mean(0)
    normal = extraction_surface.compute_normals()['Normals'].mean(0)

    for i in range(len(wall_points)):
        if -0.01 < np.dot(np.subtract(wall_points[i], center), normal) < 0.01:
            first_selection = np.vstack([first_selection, wall_points[i]])

    second_selection = remove_far_points(first_selection, center)
    second_selection = pv.PolyData(second_selection)
    second_selection = second_selection.delaunay_2d()

    return second_selection



""" # Extract surface
surface = tetmesh.extract_surface()

# Extract edges
edges = surface.extract_feature_edges(50)
edges.plot()

# Delete edges from surface
split_surface = surface.remove_cells(edges.CellData["vtkOriginalCellIds"])
split_surface.plot(show_edges=True) """
