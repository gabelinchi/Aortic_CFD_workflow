import pyvista as pv
import numpy as np

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
    return() # Return pv.PolyData of the reconstructed surface

coords =np.array([[5.34289188,17.62569172, -7.65842196],[5.04848285, 19.24025361, -6.50630999],[7.61087718, 26.40826548,  4.68480583],[6.9490417,  26.22630622,  3.47412098],[5.15573387, 18.44954256, -7.09167542],[14.43523888, 10.74629191, -0.39696456],[16.97434086,  9.96784403,  2.8223432 ],[17.49438435,  9.89288418,  3.50644435]])
print(coords)
remove_points_and_fill(pv.read('wall_mmg.vtk'), coords, plot=True)