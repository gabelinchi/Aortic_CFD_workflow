#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import meshio
import subprocess as sub

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