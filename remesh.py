#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import meshio
import subprocess as sub

#General remesh function, currently not used in the workflow (MIGHT DELETE IN FINAL PRODUCT!)
def remesh(wall_path, temp_dir,parameters, plot=False):

    file_input_location = osp.join(temp_dir, r'wall.mesh')
    file_output_location = osp.join(temp_dir, r'wall_remeshed.mesh')

    # Convert input file from stl to .mesh
    meshio.write(file_input_location, meshio.read(wall_path))

    density = parameters['mesh_density']
    sizing = parameters['sizing']
    ind = ' '

    print('Start 2D remesh of wall')

    # Run mmg
    sub.run(f"{'py -m mmgs -hausd'}{ind}{density}{ind}{'-nr'}{ind}{file_input_location}{ind}{file_output_location}{ind}{'-hsiz'}{ind}{sizing}")

    # Convert back to .vtk and plot with pyvista
    meshio.write(osp.join(temp_dir, r'wall_remeshed.vtk'), meshio.read(file_output_location))

    if plot:
        wall_remeshed = pv.read(osp.join(temp_dir, r'wall_remeshed.vtk'))
        wall_remeshed.plot(show_edges = True, text='Remeshed surface')

    print('Succesfully remeshed wall')

    return

def remesh_edge_detect(in_path, out_path, temp_path, parameters, plot=False):
    '''
    Remeshes geometry using mmg with edge detection, saves a .mesh and returns pyvista PolyData/UnstructuredGrid
    :in_path : path to input file (.mesh)
    :out_path : path to output file (.mesh)
    :parameters : dict of mmg parameters
    :plot : bool, show intermediate plots
    :returns : PyvistaPolydata of the remeshed geo
    '''
    print('Start 2D remeshing')

    density = parameters['mesh_density']
    sizing = parameters['sizing']
    angle = parameters['detection angle']
    ind = ' '

    # Run mmg
    sub.run(f"{'py -m mmgs -ar'}{ind}{angle}{ind}{'-hausd'}{ind}{density}{ind}{in_path}{ind}{out_path}{ind}{'-hsiz'}{ind}{sizing}")

    # Convert back to .vtk and plot with pyvista
    meshio.write(osp.join(temp_path, r'temp_for_plot_remesh.vtk'), meshio.read(out_path))
    remeshed = pv.read(osp.join(temp_path, r'temp_for_plot_remesh.vtk'))

    if plot==True:
        remeshed.plot(show_edges = True, text='Remeshed surface')

    print('Succesfully remeshed geometry')
    return(remeshed)