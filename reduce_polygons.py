import pyvista as pv
import numpy as np

# Function to pad an unpadded face array
def pad(faces):
    num_rows = faces.shape[0]
    return(np.hstack((np.full((num_rows, 1), 3),faces)))

def remove_point_and_fill(polydata, points_to_remove):
    # Split polydata into np arrays
    faces = polydata.faces.reshape(-1, 4)[:,[1,2,3]]
    points = polydata.points

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
        fill_mesh.plot(show_edges=True)

    # Rebuild polydata of mesh to fill (now with the gaps)
    mesh_to_fill = pv.PolyData(points, pad(faces))
    mesh_to_fill.plot(show_edges=True)

    # Combine with the patches
    reduced_surface = fill_mesh.merge(mesh_to_fill).clean()
    reduced_surface.plot(show_edges=True)
    return() # Return pv.PolyData of the reconstructed surface

#------------------------------------------------------------------------------------------------------
# Load mesh data as pv.Polydata and split into numpy arrays
inlet_remeshed = pv.read('inlet_mmg.vtk').extract_surface()
wall_remeshed = pv.read('wall_mmg.vtk').extract_surface()

remove_point_and_fill(wall_remeshed, [1, 2, 102])