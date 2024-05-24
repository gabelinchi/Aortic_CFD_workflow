#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import meshio
import subprocess as sub

def tetgen(combined_mesh, tetgen_parameters, plot=False):
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
    return grid

def mmg3d(in_path, out_path, temp_path, parameters, plot=False):
    '''
    Remeshes volume mesh using mmg3d, saves a .mesh and returns pyvista PolyData/UnstructuredGrid
    :in_path : path to input file (.mesh)
    :out_path : path to output file (.mesh)
    :parameters : dict of mmg3d parameters
    :plot : bool, show intermediate plots
    :returns : PyvistaUnstructuredGrid of the remeshed geo
    '''
    print('Start 3D mesh refinement')

    hausd = parameters['hausd']
    angle = parameters['detection angle']
    ind = ' '

    # Run mmg
    sub.run(f"{'py -m mmg3d -hausd'}{ind}{hausd}{ind}{'-ar'}{ind}{angle}{ind}{in_path}{ind}{out_path}")

    # Convert back to .vtk
    meshio.write(osp.join(temp_path, r'temp_for_plot_3d_remesh.vtk'), meshio.read(out_path))
    remeshed = pv.read(osp.join(temp_path, r'temp_for_plot_3d_remesh.vtk'))

    # Remove non-tetrahedal elements (mmg3d outputs detected edges as line segments alongside the generated mesh)
    cell_types = remeshed.celltypes
    mask = cell_types == 10 # vtk type indication for tets
    remeshed = remeshed.extract_cells(mask)

    # Plot with pyvista
    if plot==True:
        remeshed.plot(show_edges = True, text='mmg3d remeshed')

    print('Succesfully remeshed geometry')
    return(remeshed)

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