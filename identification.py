import numpy as np
import pyvista as pv



def identify_surfaces (mesh, angle, seeds=[], plot=False, include_original_connectivity=True):
    '''
    This function identifies surfaces of a pyvista 3D or surface mesh using the pv.edge_mask filter to determine
    edges. It returns a pv MultiBlock containing the surfaces as separate blocks. Also included as point and cell
    data embedded in the MultiBlock are the original point indices (point_data['orig_point_indices']) and the faces of the identified surfaces
    defined by the original point indices (cell_data['original_connectivity']).

    bugs & limis:
    BUG If closest point to seed is on an edge, function misbehaves
    LIMIT Function only works for meshes where te surface consists of only triangles

    :mesh   : pyvista PolyData or Multiblock, mesh to be split up
    :angle  : int, sharp edge detection angle
    :seeds  : 2D numpy array containing xyz coords, can be used to order the surfaces if points near the surfaces
              are known. If empty, seeds are automatically generated 
              !WARNING! known bug: if closest point to seed is on an edge, function misbehaves
    :plot   : bool, show intermediate plots
    :include_original_connectivity: bool
    :returns: pv.MultiBlock containging identified surfaces
    '''

    # Plotter tool for better visualisation
    def boolplot (surf, boolarray, text=''):
        plt = pv.Plotter()
        plt.add_text(text)
        plt.add_mesh(surf, style='wireframe', color='black')
        plt.add_mesh(surf.extract_cells(boolarray), color = 'red')
        plt.show()
        return

    # Add original point data and extract surface
    mesh['orig_point_indices'] = np.arange(mesh.n_points, dtype=np.int32)
    surf = mesh.extract_surface()

    # Tag edgenodes
    edgenodes = surf.edge_mask(angle)

    # Extract np arrays
    faces = surf.faces.reshape(-1, 4)[:,[1,2,3]]

    # Create array of edge faces
    edgenode_indices = np.where(edgenodes)[0] # Converts bool array to array of indices
    edgefaces = np.any(np.isin(faces, edgenode_indices), axis=1) # Creates bool array of rows containing edgepoints

    if plot==True:
        boolplot(surf, edgefaces, text='Edges')

    # Create surfaces array, first column is edges:
    surfaces = edgefaces.reshape(-1,1)

    # Create PyVista PolyData block
    surf_block = pv.MultiBlock()

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
        if len(seeds) >= num_surfaces:
            front = surf.find_closest_cell(seeds[num_surfaces-1, :])

        # Add face to surface
        surfaces[front, num_surfaces] = True

        if plot==True:
            boolplot(surf, surfaces[:, num_surfaces], text = 'Seed')

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
                break
            counter += 1

        if plot==True:
            boolplot(surf, surfaces[:, num_surfaces], text = f'Identified surface no. {num_surfaces}')

        # Create new surface as Polydata
        new_surf = surf.extract_cells(surfaces[:, num_surfaces])

        # Get original node data and add as cell data
        if include_original_connectivity == True:
            new_surf_faces = new_surf.cells.reshape(-1, 4)[:,[1,2,3]]
            map = lambda x: new_surf['orig_point_indices'][x]
            new_surf.cell_data['original_connectivity'] = map(new_surf_faces)

        # Add new surface to the block
        surf_block.append(new_surf)

        num_surfaces += 1
    return surf_block

#-- for debugging & writing --#
""" #Load input files

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

surfaces = identify_surfaces(tetmesh, 30, seeds=seeds, include_original_connectivity=True)

tetmesh.extract_points(surfaces[0].cell_data['original_connectivity'].flatten()).plot()
surfaces[0].plot(text='inlet')
surfaces[1].plot(text='outlet')
surfaces[2].plot(text='wall') """
#-------------------------------#