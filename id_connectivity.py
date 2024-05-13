import numpy as np
import pyvista as pv

def boolplot (surf, boolarray, text=''):
    plt = pv.Plotter()
    plt.add_text(text)
    plt.add_mesh(surf, style='wireframe', color='black')
    plt.add_mesh(surf.extract_cells(boolarray), color = 'red')
    plt.show()
    return

def id_con (surf, seeds=None, plot=False):

    # Tag edgenodes
    edgenodes = surf.edge_mask(30)

    # Extract np arrays
    faces = surf.faces.reshape(-1, 4)[:,[1,2,3]]

    # Create array of edge faces
    edgenode_indices = np.where(edgenodes)[0] # Converts bool array to array of indices
    edgefaces = np.any(np.isin(faces, edgenode_indices), axis=1) # Creates bool array of rows containing edgepoints

    # Create surfaces array, first column is edges:
    surfaces = edgefaces.reshape(-1,1)

    # Start of identification loop
    num_surfaces = 1

    print('Start surface identification')

    while num_surfaces < 10:
        # Create a new surface
        surfaces = np.hstack((surfaces, np.full((len(faces), 1), False)))
        
        # Check if there are unlisted faces
        if np.any(~np.any(surfaces, axis=1)) == False:
            surfaces = np.delete(surfaces, num_surfaces, 1)
            print(num_surfaces-1, 'surfaces found in total')
            break
        print('Unidentified surface found')

        # Pick a nonlisted face
        front = np.argmax(~np.any(surfaces, axis=1)) # Front contains indices, not bools
        if len(seeds) != 0 and len(seeds) >= num_surfaces:
            front = surf.find_closest_cell(seeds[num_surfaces-1, :])

        # Add face to surface
        surfaces[front, num_surfaces] = True

        #boolplot(surf, surfaces[:, num_surfaces], text = 'Seed')

        counter = 0
        while counter < 10000:
            # Find neigbouring faces that arent part of a surface, add to surface
            frontnodes = faces[front].flatten() # Node indices
            neighbours = np.any(np.isin(faces, frontnodes), axis=1) # All neighbouring faces to front faces

            #boolplot(surf, neighbours, text = 'All neighbours')

            neighbours[front | np.any(surfaces[:, 1:], axis=1)] = False # Removes front faces and faces part of a surface

            #boolplot(surf, neighbours, text = 'Cleaned neighbours')

            surfaces[:,num_surfaces] = surfaces[:,num_surfaces] | neighbours # Adds neighbours to surface
            neighbours[edgefaces] = False # Removes edges
            front = neighbours # Creates new front


            # Check if the front is empty
            if np.any(front) == False:
                #boolplot(surf, surfaces[:, num_surfaces], text = 'Surface')
                break
            counter += 1
        num_surfaces += 1
    return surfaces

#-- for debugging & writing --#
#Load input files

tetmesh = pv.read('3D_output_mesh.vtk')

#paths to the aorta geometry
inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"
outlet_path = "geometries\input\outlet.stl"

#Reading the files with pyvista
inlet = pv.read(inlet_path)
wall = pv.read(wall_path)
outlet = pv.read(outlet_path)

seeds = np.array([inlet.points.mean(0), outlet.points.mean(0)])
print(seeds)

surf = tetmesh.extract_surface()

surfaces = id_con(surf, seeds=seeds)

inlet = surf.extract_cells(surfaces[:, 1])
outlet = surf.extract_cells(surfaces[:, 2])
wall = surf.extract_cells(surfaces[:, 3])


inlet_faces = inlet.cells.reshape(-1, 4)[:,[1,2,3]]
print(inlet['point_ind'])

with np.nditer(inlet_faces, op_flags=['readwrite']) as it:
        for x in it:
            x[...]= inlet['point_ind'][inlet_faces[...]]

print(inlet_faces)