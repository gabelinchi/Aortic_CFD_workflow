#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import meshio
import subprocess as sub

def remesh(inlet_path, wall_path, output_path, parameters, plot=False):
    # Convert input file from stl to .mesh
    meshio.write('inlet.mesh', meshio.read(inlet_path))
    meshio.write('wall.mesh', meshio.read(wall_path))
    meshio.write('outlet.mesh', meshio.read(output_path))

    print('Succesfully imported geometry')

    if plot:
        #import .stl's to pyvista
        pv_inlet_mesh = pv.read(inlet_path)
        pv_wall_mesh = pv.read(wall_path)
        pv_outlet_mesh = pv.read(output_path)

        pv_inlet_mesh.plot(show_edges = True)
        pv_wall_mesh.plot(show_edges = True)
        pv_outlet_mesh.plot(show_edges = True)

    """ # Extract boundaries
    inlet_edge = pv_inlet_mesh.extract_feature_edges(manifold_edges=False, non_manifold_edges=False, feature_edges = False, boundary_edges=True)
    inlet_edge.plot(show_edges=True)

    # Save boundaries as .mesh
    inlet_edge.save('geometries/temp/inlet_edge.ply')
    meshio.write('geometries/temp/inlet_edge.mesh', meshio.read('geometries/temp/inlet_edge.ply'))
    pv.read('geometries/temp/inlet_edge.ply').plot(show_edges=True) #plot inlet vertices """
    density = parameters['mesh_density']
    sizing = parameters['sizing']
    indentation = ' '

    #inlet_command = f"{'py -m mmgs -hausd'}{indentation}{density}{indentation}{'inlet.mesh inlet_mmg.mesh -hsiz'}{indentation}{sizing}"


    # Run mmg
    sub.run(f"{'py -m mmgs -hausd'}{indentation}{density}{indentation}{'inlet.mesh inlet_mmg.mesh -hsiz'}{indentation}{sizing}")
    sub.run(f"{'py -m mmgs -hausd'}{indentation}{density}{indentation}{'wall.mesh wall_mmg.mesh -hsiz'}{indentation}{sizing}")
    sub.run(f"{'py -m mmgs -hausd'}{indentation}{density}{indentation}{'outlet.mesh outlet_mmg.mesh -hsiz'}{indentation}{sizing}")

    # Convert back to .vtk and plot with pyvista
    meshio.write('inlet_mmg.vtk', meshio.read('inlet_mmg.mesh'))
    meshio.write('wall_mmg.vtk', meshio.read('wall_mmg.mesh'))
    meshio.write('outlet_mmg.vtk', meshio.read('outlet_mmg.mesh'))

    inlet_remeshed = pv.read('inlet_mmg.vtk')
    wall_remeshed = pv.read('wall_mmg.vtk')
    outlet_remeshed = pv.read('outlet_mmg.vtk')

    inlet_remeshed.plot(show_edges = True)
    wall_remeshed.plot(show_edges = True)
    outlet_remeshed.plot(show_edges = True)

    print('Succesfully remeshed geometry')

    if plot:
        combined = inlet_remeshed + wall_remeshed + outlet_remeshed
        combined.plot(show_edges = True)

    """ # Merge and plot with pv
    wall_and_inlet = wall_remeshed.merge(inlet_remeshed).clean(tolerance=0.001) #combines two meshes and removes duplicate points
    wall_and_inlet.plot(show_edges=True) """
    return inlet_remeshed, wall_remeshed, outlet_remeshed

def remesh_edge_detect(in_path, out_path, parameters, plot=False):
    '''
    Remeshes geometry using mmg with edge detection, saves a .mesh and returns pyvista PolyData/UnstructuredGrid
    Detection angle and filenames are still hardcoded, should be fixed
    :in_path : path to input file (.mesh)
    :out_path : path to output file (.mesh)
    :parameters : dict of mmg parameters
    :plot : bool, show intermediate plots
    :returns : PyvistaPolydata of the remeshed geo
    '''
    density = parameters['mesh_density']
    sizing = parameters['sizing']
    angle = parameters['detection angle']
    ind = ' '

    # Run mmg
    sub.run(f"{'py -m mmgs -ar'}{ind}{angle}{ind}{'-hausd'}{ind}{density}{ind}{in_path}{ind}{out_path}{ind}{'-hsiz'}{ind}{sizing}")

    # Convert back to .vtk and plot with pyvista
    meshio.write('temp_for_plot_remesh.vtk', meshio.read(out_path))
    remeshed = pv.read('temp_for_plot_remesh.vtk')

    if plot==True:
        remeshed.plot(show_edges = True)

    print('Succesfully remeshed geometry')
    return(remeshed)