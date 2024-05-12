
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
import identification as id

#paths to the aorta geometry

inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"
outlet_path = "geometries\input\outlet.stl"

#Reading the files with pyvista
inlet = pv.read(inlet_path)
wall = pv.read(wall_path)
outlet = pv.read(outlet_path)

#Meshing parameters
mmg_parameters = {
    'mesh_density': '0.1',
    'sizing': '1',
    'detection angle': '35'}

tetgen_parameters = dict(
    order=1, 
    mindihedral=20, 
    minratio=1.5)

#Plotting boolean, when True: code generates intermediate plots of workflow
show_plot = True


print('Setup and import done')

#Cut the wall geometry after the aortic root
wall_cut = cutting.main_cutter(inlet, wall, plot=show_plot)
pv.save_meshio("wall_cut.mesh", wall_cut)

#Create caps
inlet_cap, outlet_cap = cap(wall_cut, plot=show_plot)
pv.save_meshio('inlet_cap.mesh', inlet_cap)
pv.save_meshio('outlet_cap.mesh', outlet_cap)

#Combine parts
combined = (wall_cut+inlet_cap+outlet_cap).clean()

if show_plot:
    plt = pv.Plotter()
    plt.add_mesh(combined, style='wireframe')
    edgetest = combined.extract_feature_edges(boundary_edges=True, non_manifold_edges=True, manifold_edges=False, feature_edges=False)
    plt.show()

pv.save_meshio('combined_mesh.mesh', combined)



#Run remesh (takes the file Path !!! as input)
combined_remeshed = remesh.remesh_edge_detect('combined_mesh.mesh', 'combined_mmg.mesh', mmg_parameters, plot=show_plot) #HARDCODED FILENAMES!!
combined_remeshed = combined_remeshed.extract_surface().triangulate()


#Make a 3D mesh from the combined mesh
tetmesh = volume_mesh.volume_meshing(combined_remeshed, tetgen_parameters, False)

tetmesh.plot(show_edges = True)

inlet_selected = id.selection(tetmesh, inlet_cap)
outlet_selected = id.selection(tetmesh, outlet_cap)
plt = pv.Plotter()
plt.add_mesh(inlet_selected, show_edges = True, color = 'red')
plt.add_mesh(outlet_selected, show_edges = True, color = 'blue')
plt.show()


#print(mesh)

