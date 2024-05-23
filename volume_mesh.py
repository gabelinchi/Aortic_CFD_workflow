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

        clipped = grid.clip('z', crinkle=True)
        # advanced plotting
        plotter = pv.Plotter()
        plotter.add_mesh(clipped, 'lightgrey', lighting=True, show_edges=True)
        plotter.add_mesh(combined_mesh, 'r', 'wireframe')
        plotter.add_legend([[' Input Mesh ', 'r'],
                            [' Tessellated Mesh ', 'black']])
        plotter.add_text('Initial 3D mesh')
        plotter.show()

    return grid