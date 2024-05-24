import pyvista as pv
import numpy as np

def write_sol(mesh, surf, parameters, dir, plot=False):
    '''
    Function that creates a .sol file for a mesh. This file can be used by mmg to specify local mesh density
    paremeter. Outputs the input mesh with an added array 'sol' containing the data written to the .sol file
    :mesh   : input mesh, must be Pyvista UnstructuredGrid volume mesh
    :surf   : wall of the geometry, can be any Pyvista surface
    :dist   : float, distance to surf that defines where the high density area ends
    :dlow   : far from wall density parameter for mmg
    :dhigh  : close to wall density parameter for mmg
    :dir    : location to sace .sol
    plot    : show plots
    '''
    print('start tagging')

    # Unpack parameters
    dist = parameters['bl_thickness']
    dlow = parameters['edgelength']
    dhigh = parameters['bl_edgelength']

    # Get closest point to wall for each point in mesh & calculate distance
    closest_points = surf.find_closest_cell(mesh.points, return_closest_point=True)[1]
    d_exact = np.linalg.norm(mesh.points - closest_points, axis=1)

    # Apply the density parameters
    density = np.where(d_exact > dist, dlow, dhigh)
    mesh["sol"] = density
    length = len(density)
    
    # Plot results
    if plot:
        clipped = mesh.clip('z', crinkle=True)
        plotter = pv.Plotter()
        plotter.add_mesh(clipped, 'lightgrey', lighting=True, show_edges=True, scalars='sol')
        plotter.add_mesh(surf, 'r', 'wireframe')
        plotter.add_legend([[' Input Mesh ', 'r'], [' Tessellated Mesh ', 'black']])
        plotter.add_text('Initial 3D mesh with local tags')
        plotter.show()

    # Write file
    header = 'MeshVersionFormatted \n2\n\nDimension 3\n\nSolAtVertices\n' + str(length) + '\n1 1\n'
    np.savetxt(dir, density.reshape(-1,1), header=header, footer='End', comments='', fmt='%f')

    return(mesh)