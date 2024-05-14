
#import modules
import sys
import os
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
#----------------------------------------------------------------------------------------------------------------------------
# Setup
#----------------------------------------------------------------------------------------------------------------------------

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

#Angle for identification
id_angle = 30

#Plotting boolean, when True: code generates intermediate plots of workflow
show_plot = False


print('Setup and import done')

#--------------------------------------------------------------------------------------------------------------------------
# 3D-meshing algorithm
#--------------------------------------------------------------------------------------------------------------------------

#Cut the wall geometry after the aortic root
wall_cut = cutting.main_cutter(inlet, wall, plot=show_plot)
pv.save_meshio("temp\wall_cut.mesh", wall_cut)

#Create caps
inlet_cap, outlet_cap = cap(wall_cut, plot=show_plot)
pv.save_meshio('temp\inlet_cap.mesh', inlet_cap)
pv.save_meshio('temp\outlet_cap.mesh', outlet_cap)

#Combine cutted wall and inlet/outlet caps
combined = (wall_cut + inlet_cap + outlet_cap).clean()
print('Meshes succesfully combined')

#Plot result of mesh combining. When edgetest gives an error, the structure has dublicate faces/nodes
if show_plot:
    plt = pv.Plotter()
    plt.add_mesh(combined, style='wireframe')
    edgetest = combined.extract_feature_edges(boundary_edges=True, non_manifold_edges=True, manifold_edges=False, feature_edges=False)
    plt.show()

pv.save_meshio('temp\combined_mesh.mesh', combined)



#Run remesh (takes predetermined internally defined file path as input, DON'T CHANGE)
combined_remeshed = remesh.remesh_edge_detect('temp\combined_mesh.mesh', 'temp\combined_mmg.mesh', mmg_parameters, plot=show_plot)
#triangulation step to make sure Tetgen only gets triangels as input
combined_remeshed = combined_remeshed.extract_surface().triangulate()


#Make a 3D mesh from the combined mesh
tetmesh = volume_mesh.volume_meshing(combined_remeshed, tetgen_parameters, plot=show_plot)

#Save 3D mesh
tetmesh.save('3D_output_mesh.vtk')

#Plot 3D_mesh
if show_plot:
    tetmesh.plot(show_edges = True)

#----------------------------------------------------------------------------------------------------------------------------
# Identification
#----------------------------------------------------------------------------------------------------------------------------

#Create the seeds for the surface identification based on the center points. INLET FIRST!, OUTLET SECOND!
seeds = np.array([inlet.points.mean(0),outlet.points.mean(0)])

#Detect the surfaces of the 3D mesh whilst keeping the original ID's, if there is a seed delivered
if len(seeds) > 0:
    surface_identification = id.identify_surfaces(tetmesh, id_angle, seeds, show_plot)

#Seperate the indentified surfaces in inlet/outlet/wall
id_inlet = surface_identification[0]
id_outlet = surface_identification[1]
id_wall = surface_identification[2]

#Plot the identifies surfaces for general overview
if show_plot:
    plt = pv.Plotter()
    plt.add_mesh(id_inlet, show_edges = True, color = 'red')
    plt.add_mesh(id_outlet, show_edges = True, color = 'blue')
    plt.add_mesh(id_wall, show_edges = True, color = 'green')
    plt.show()

#----------------------------------------------------------------------------------------------------------------------------
# Mapping
#----------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------------
# FEBio
#----------------------------------------------------------------------------------------------------------------------------


#Deletes all the temporary files in the temp folder (have to fix)
""" temp_files = glob('temp')
for f in temp_files:
    os.remove(f) """

print('Done!')

