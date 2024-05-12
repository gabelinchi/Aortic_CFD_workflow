import pyvista as pv
import tetgen as tet

def volume_meshing(combined_mesh, tetgen_parameters, plot=False):
    '''
    3D-meshing function that uses TetGen to create a 3D-mesh from a closed surface mesh. The function has the following inputs:
    combined_mesh: Closed surface mesh in pyvista format
    tetgen_parameters: meshing parameters as defined in the Setup part from the main_workflow
    returns: 3D-mesh

    The function has een advances plotting part to visually inspect the quality of the mesh
    '''

    print('Start 3D-meshing')
    # create 3D tetmesh from surface mesh
    tetmesh = tet.TetGen(combined_mesh)
    tetmesh.tetrahedralize(**tetgen_parameters)
    grid = tetmesh.grid
    

    print('3D meshing succesfull')
    if plot:
        grid.plot(show_edges=True)

        # Save the generated mesh as a .stl
        grid.save('aorta_tetmesh.vtk', binary=False)

        # get cell centroids
        cells = grid.cells.reshape(-1, 5)[:, 1:]
        cell_center = grid.points[cells].mean(1)

        # extract cells below the 0 xy plane
        mask = cell_center[:, 2] < 0
        cell_ind = mask.nonzero()[0]
        subgrid = grid.extract_cells(cell_ind)

        # advanced plotting
        plotter = pv.Plotter()
        plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
        plotter.add_mesh(combined_mesh, 'r', 'wireframe')
        plotter.add_legend([[' Input Mesh ', 'r'],
                            [' Tessellated Mesh ', 'black']])
        plotter.show()

    return grid