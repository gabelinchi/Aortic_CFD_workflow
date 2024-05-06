#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import remesh
import meshclosing
import utils as ut

#paths to the aorta geometry
inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"
outlet_path = "geometries\input\outlet.stl"

# Import geomtetry to pv

inlet = pv.read(inlet_path)
outlet = pv.read(outlet_path)
wall = pv.read(wall_path)

def cut(point, normal, wall):
    """
    Function that cuts a vessel geometry along a plane defined by a point & vector, and keeps all regions upstream
    of the point.
    :arg1 point: numpy array, xyz
    :arg2 normal: vector
    :arg3 wall: pyvista PolyData

    returns clipped geometry as PolyData
    """
    # Clip geometry
    reg1, reg2 = wall.clip(normal=normal, origin=point, return_clipped=True) #reg2 is in direction of vector

    # Extract geometry to keep from reg2
    reg2 = reg2.connectivity('closest', point)

    # Extract geometry to keep from reg1
    reg1 = reg1.connectivity('all')
    del_id = reg1.point_data['RegionId'][reg1.find_closest_point(point)]
    reg1 = reg1.extract_cells(np.where(reg1.cell_data['RegionId'] != del_id)[0])

    # Recombine regions
    clipped = reg1.merge(reg2).clean()

    return(clipped)

def get_clip_perimeter(point, normal, wall):
    '''
    Function that cuts a vessel geometry along a plane defined by a point & vector, and returns the resulting perimeter
    :arg1 point: numpy array, xyz
    :arg2 normal: vector
    :arg3 wall: pyvista PolyData

    returns pyvista PolyData of the perimeter resulting from the cut
    '''
    # Clip geometry
    clipped = wall.clip(normal=normal, origin=point) #reg2 is in direction of vector

    # Extract edges
    edges = clipped.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)

    # Extract edge closest to point
    edge = edges.connectivity('closest', point)

    return(edge)

