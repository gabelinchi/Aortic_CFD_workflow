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
def remesh(inlet_path, wall_path, output_path, parameters, plot=False):
    # Convert input file from stl to .mesh
    meshio.write('inlet.mesh', meshio.read(inlet_path))
    meshio.write('wall.mesh', meshio.read(wall_path))
    meshio.write('outlet.mesh', meshio.read(output_path))

    density = parameters['mesh_density']
    sizing = parameters['sizing']
    indentation = ' '

    #inlet_command = f"{'py -m mmgs -hausd'}{indentation}{density}{indentation}{'inlet.mesh inlet_mmg.mesh -hsiz'}{indentation}{sizing}"

    print('Start 2D remesh')

    # Run mmg
    sub.run(f"{'py -m mmgs -hausd'}{indentation}{density}{indentation}{'inlet.mesh inlet_mmg.mesh -hsiz'}{indentation}{sizing}")
    sub.run(f"{'py -m mmgs -hausd'}{indentation}{density}{indentation}{'wall.mesh wall_mmg.mesh -hsiz'}{indentation}{sizing}")
    sub.run(f"{'py -m mmgs -hausd'}{indentation}{density}{indentation}{'outlet.mesh outlet_mmg.mesh -hsiz'}{indentation}{sizing}")

    # Convert back to .vtk and plot with pyvista
    meshio.write('inlet_mmg.vtk', meshio.read('inlet_mmg.mesh'))
    meshio.write('wall_mmg.vtk', meshio.read('wall_mmg.mesh'))
    meshio.write('outlet_mmg.vtk', meshio.read('outlet_mmg.mesh'))

    if plot:
        inlet_remeshed = pv.read('inlet_mmg.vtk')
        wall_remeshed = pv.read('wall_mmg.vtk')
        outlet_remeshed = pv.read('outlet_mmg.vtk')

        inlet_remeshed.plot(show_edges = True)
        wall_remeshed.plot(show_edges = True)
        outlet_remeshed.plot(show_edges = True)

    print('Succesfully remeshed geometry')

    return inlet_remeshed, wall_remeshed, outlet_remeshed

def remesh_edge_detect(in_path, out_path, temp_path, parameters, plot=False):
    '''
    Remeshes geometry using mmg with edge detection, saves a .mesh and returns pyvista PolyData/UnstructuredGrid
    Detection angle and filenames are still hardcoded, should be fixed
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
        remeshed.plot(show_edges = True)

    print('Succesfully remeshed geometry')
    return(remeshed)