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

def remesh_edge_detect(mesh, parameters, plot=False):
    '''
    Remeshes geometry using mmg with edge detection
    Detection angle and filenames are still hardcoded, should be fixed
    mesh must be pyvista PolyData
    '''
    if plot==True:
        mesh.plot(show_edges = True)

    density = parameters['mesh_density']
    sizing = parameters['sizing']
    indentation = ' '

    # Convert to .mesh
    mesh.save('wall_and_io.ply')
    meshio.write('wall_and_io.mesh', meshio.read('wall_and_io.ply'))

    # Run mmg
    sub.run(f"{'py -m mmgs -ar 35 -hausd'}{indentation}{density}{indentation}{'wall_and_io.mesh wall_and_io_mmg.mesh -hsiz'}{indentation}{sizing}")

    # Convert back to .vtk and plot with pyvista
    meshio.write('wall_and_io_mmg.vtk', meshio.read('wall_and_io_mmg.mesh'))
    wall_and_io_remeshed = pv.read('wall_and_io_mmg.vtk')

    if plot==True:
        wall_and_io_remeshed.plot(show_edges = True)

    print('Succesfully remeshed geometry')
    return wall_and_io_remeshed