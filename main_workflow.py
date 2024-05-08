
#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import remesh
import utils as ut
import cutting
import volume_mesh
from capping import cap

#paths to the aorta geometry

inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"
outlet_path = "geometries\input\outlet.stl"

#Filename of cutted wall, inlet and outlet
fn_wall_cut = "wall_cut.mesh"
fn_inlet = "inlet.mesh"
fn_outlet = "outlet.mesh"

#Reading the files with pyvista
inlet = pv.read(inlet_path)
wall = pv.read(wall_path)
outlet = pv.read(outlet_path)

#Meshing parameters
adjustment = np.array([0,0,20])
search_tolerance = 0.1 #Need to check this value and what it means. 0.1 works with the expected output nodes

mmg_parameters = {
    'mesh_density': '0.1',
    'sizing': '1',
    'detection angle': '35'}

tetgen_parameters = dict(
    order=1, 
    mindihedral=20, 
    minratio=1.5)

#Cut the wall geometry after the aortic root
wall_cut = cutting.main_cutter(inlet, wall, plot=True)
pv.save_meshio(fn_wall_cut, wall_cut)

#Create caps
inlet_cap, outlet_cap = cap(wall_cut, plot=True)
pv.save_meshio('inlet_cap.mesh', inlet_cap)
pv.save_meshio('outlet_cap.mesh', outlet_cap)

#Combine parts
combined = (wall_cut+inlet_cap+outlet_cap).clean()
plt = pv.Plotter()
plt.add_mesh(combined, style='wireframe')
edgetest = combined.extract_feature_edges(boundary_edges=True, non_manifold_edges=True, manifold_edges=False, feature_edges=False)
plt.show()
pv.save_meshio('combined_mesh.mesh', combined)



#Run remesh (takes the file Path !!! as input)
combined_remeshed = remesh.remesh_edge_detect('combined_mesh.mesh', 'combined_mmg.mesh', mmg_parameters, plot=True) #HARDCODED FILENAMES!!
combined_remeshed = combined_remeshed.extract_surface().triangulate()


#Make a 3D mesh from the combined mesh
tetmesh = volume_mesh.volume_meshing(combined_remeshed, tetgen_parameters, False)

tetmesh.plot(show_edges = True)

def remove_far_points(selection, center):
    print(selection)
    distances = np.array([])
    for i in range(len(selection)):
        distance = np.linalg.norm(selection[i] - center)
        distances = np.append(distances, distance)
    
    avg_distance = np.average(distances)
    new_selection = np.empty((0,3))

    for i in range(len(distances)):
        if distances[i] < 4 * avg_distance:
            new_selection = np.vstack([new_selection, selection[:, i]])
        

    return new_selection

#Identificate inlet/outlet
def identification(wall, extraction_surface):
    wall_surface = wall.extract_surface()
    wall_edge = wall_surface.extract_feature_edges(45)
    wall_points = wall_surface.points
    first_selection = np.empty((0,3))
    center = extraction_surface.points.mean(0)
    normal = extraction_surface.compute_normals()['Normals'].mean(0)

    for i in range(len(wall_points)):
        if -0.01 < np.dot(np.subtract(wall_points[i], center), normal) < 0.01:
            first_selection = np.vstack([first_selection, wall_points[i]])

    second_selection = remove_far_points(first_selection, center)
    second_selection = pv.PolyData(second_selection)
    second_selection = second_selection.delaunay_2d()

    return second_selection, wall_edge

inlet_selected, edge = identification(tetmesh, inlet_cap)
plt = pv.Plotter()
plt.add_mesh(inlet_selected, show_edges = True, color = 'red')
plt.add_mesh(edge, color = 'blue')
plt.show()


#print(mesh)

#use an excessive amount of comments on everything.
#extra line of code changed
