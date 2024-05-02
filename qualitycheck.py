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
import trimesh as tri
import pyacvd as acvd

'''
The goal of this script is to calculate quality metrics for surface meshes, to begin I will use the wall
as a comparison point.
'''

# Import geometry
inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"
outlet_path = "geometries\input\outlet.stl"

wall_mesh = pv.read(wall_path)
outlet_mesh = pv.read(outlet_path)
inlet_mesh = pv.read(inlet_path)

""" qual=wall_mesh.compute_cell_quality()
qual.plot()
print('qual')
print(wall_mesh.compute_cell_quality().cell_data['CellQuality'])
print(np.mean(wall_mesh.compute_cell_quality().cell_data['CellQuality'])) """


# Remesh with mmg
inlet_remeshed_mmg, wall_remeshed_mmg, outlet_remeshed_mmg = remesh.remesh(inlet_path, wall_path, outlet_path)


# Remesh with trimesh
def extract_nodes_and_faces(unstrucuredgrid):
    """
    Converts a pyvista unstructured grid into an array of faces and an array of nodes
    Written by Yarran
    """
    nodes = unstrucuredgrid.points
    poly = unstrucuredgrid.extract_surface()
    faces = np.asarray(poly.faces).reshape((-1, 4))[:, 1:]
    return(nodes, faces)

inlet_nodes, inlet_faces = extract_nodes_and_faces(inlet_mesh)
outlet_nodes, outlet_faces = extract_nodes_and_faces(outlet_mesh)
wall_nodes, wall_faces = extract_nodes_and_faces(wall_mesh)

inlet_from_arrays = pv.PolyData.from_regular_faces(inlet_nodes, inlet_faces)
outlet_from_arrays = pv.PolyData.from_regular_faces(outlet_nodes, outlet_faces)
wall_from_arrays = pv.PolyData.from_regular_faces(wall_nodes, wall_faces)

inlet_remeshed_nodes, inlet_remeshed_faces = tri.remesh.subdivide_to_size(inlet_nodes, inlet_faces, 1.97)
outlet_remeshed_nodes, outlet_remeshed_faces = tri.remesh.subdivide_to_size(outlet_nodes, outlet_faces, 1.97)
wall_remeshed_nodes, wall_remeshed_faces = tri.remesh.subdivide_to_size(wall_nodes, wall_faces, 1.97)

inlet_remeshed_trimesh = pv.PolyData.from_regular_faces(inlet_remeshed_nodes, inlet_remeshed_faces)
outlet_remeshed_trimesh = pv.PolyData.from_regular_faces(outlet_remeshed_nodes, outlet_remeshed_faces)
wall_remeshed_trimesh = pv.PolyData.from_regular_faces(wall_remeshed_nodes, wall_remeshed_faces)


# Remesh with pyacvd

def ac_remesh(mesh, subdivide, cluster, plots=False):
    """
    remeshes a surface using pyacvd
    """
    if cluster == 0:
        cluster = mesh.n_points * 7.5
    if plots == True:
        mesh.plot(show_edges=True)
    clus = acvd.Clustering(mesh)
    clus.subdivide(subdivide)
    clus.cluster(cluster)
    if plots == True:
        clus.plot()
    remesh = clus.create_mesh().clean(lines_to_points=False, polys_to_lines=False)
    if plots == True:
        remesh.plot(show_edges=True)
    return(remesh)

inlet_remeshed_acvd = ac_remesh(inlet_mesh, 4, 0)
wall_remeshed_acvd = ac_remesh(wall_mesh, 4, 0)
outlet_remeshed_acvd = ac_remesh(outlet_mesh, 4, 0)

(inlet_remeshed_mmg+outlet_remeshed_mmg+wall_remeshed_mmg).plot(show_edges=True)
(inlet_remeshed_trimesh+outlet_remeshed_trimesh+wall_remeshed_trimesh).plot(show_edges=True)
(inlet_remeshed_acvd+outlet_remeshed_acvd+wall_remeshed_acvd).plot(show_edges=True)

# Start of quality asessment, quality is determined on two points: similarity to the original mesh, and quality of the resulting
# mesh. The resulting mesh quality is based on the scaled jacobian method from pyvista
# To assess the similarity, two meshes are compared and the distance between any point on the original mesh and the new mesh is found

def mesh_quality(startmesh, remesh, meshname='', plot=False):
    '''
    Function to assess the quality of a surface mesh. Two factors are consdidered: similarity to the input in terms of shape
    and quality of the elements of the remesh.
    The first metric is based on the smallest distance between each node of the starting mesh and the surface of the remesh
    The second metric is based on the scaled jacobian of eacht cell element


    startmesh: (pyvista.PolyData)
    remesh: (pyvista.Polydata)
    meshname: (str) Name for the comparison, use to keep oversight if running the function multiple times
    plot: (bool) Show 3D plots of the quality and similarity

    Returns: (mean_similarity(#lower is better), mean_quality)
    '''
    closest_cells, closest_points = remesh.find_closest_cell(startmesh.points, return_closest_point=True)
    d_exact = np.linalg.norm(startmesh.points - closest_points, axis=1)
    startmesh.point_data['absdistance'] = d_exact
    quality = remesh.compute_cell_quality()
    mean_similarity = np.mean(d_exact)
    mean_quality = np.mean(quality.cell_data['CellQuality'])
    print('Mean similarity of',meshname,'(lower is better) =', mean_similarity)
    print('Mean quality of',meshname,' =', mean_quality)
    if plot==True:
        startmesh.plot(text=(meshname +' similarity to original mesh'))
        quality.plot(text=('Cell quality of ' + meshname))
    return(mean_similarity, mean_quality)

print(wall_remeshed_mmg)
print(wall_remeshed_acvd)
print(wall_remeshed_trimesh)
mesh_quality(wall_mesh, wall_remeshed_mmg, meshname='mmg remesh', plot=True)
mesh_quality(wall_mesh, wall_remeshed_acvd, meshname='acdv remesh', plot=True)
mesh_quality(wall_mesh, wall_remeshed_trimesh, meshname='trimesh remesh', plot=True)
