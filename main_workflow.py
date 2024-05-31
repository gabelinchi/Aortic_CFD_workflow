
#import modules
import sys
import os
import os.path as osp
from pathlib import Path
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
from tkinter import Tk
from tkinter.filedialog import askdirectory
import mapping
import xml.etree.ElementTree as ET
import subprocess
import febioxml as feb
import quality_control
import shutil
#----------------------------------------------------------------------------------------------------------------------------
# Setup
#----------------------------------------------------------------------------------------------------------------------------

#Meshing parameters
mmg_parameters = {
    'mesh_density': '0.1', #hausdorf parameter of mmg, defines amount of added detail at curvature
    'sizing': '1', #forces similarly sized poly's, legacy option (used for stitching), can possibly be dropped
    'detection angle': '45'}

mmg3d_parameters = {
    'hausd': '0.1',
    'detection angle': '45'}

mmg3d_sol_parameters = {
    'bl_thickness': 1,
    'bl_edgelength': 1,
    'edgelength': 15}

tetgen_parameters = dict(
    order=1, 
    mindihedral=20, 
    minratio=1.5,
    nobisect=True,
    fixedvolume=True,
    maxvolume=1)

#Run qualifications (code will terminate if not met)
max_elements = 1000000
min_jacobian = 0.1
max_aspect   = 5

#Angle for identification of surfaces
id_angle = 35

#Mapping interpolation options
intp_options = {
    'zero_boundary_dist': 0.2,  # percentage of border with zero velocity (smooth damping at the border)
    'zero_backflow': False,     # set all backflow components to zero
    'kernel': 'linear',         # RBF interpolation kernel (linear is recommended)
    'smoothing': 0.5,           # interpolation smoothing, range recommended [0, 2]
    'degree': 0,                # degree of polynomial added to the RBF interpolation matrix
    'hard_noslip': False}       # check if no-slip condition on walls is met


#Plotting boolean, when True: code generates intermediate plots of workflow
show_plot = False


#Create file environment before looping
#Get directory of main_workflow file
file_dir = osp.dirname(osp.realpath(__file__))  

#Ask user input, velocity profile and output directory. If output directory is not given it creates an output directory in file_dir
input_dir = askdirectory(title='Select Folder Containing Geometries') # shows dialog box and return the path 
vel_profile_dir = askdirectory(title='Select Velocity Profile Folder') # shows dialog box and return the path
output_dir = askdirectory(title='Select Output Folder')

if output_dir == '':
    output_dir = osp.join(file_dir, r'output')
    os.makedirs(output_dir, exist_ok=True)

#Create a temporary directory
temp_dir = osp.join(file_dir, r'temp')
os.makedirs(temp_dir, exist_ok=True)

#Select input folder names and count the amount of input geometries
input_list = os.listdir(input_dir)
n_geometries = len(input_list)


print('Setup done')

#-----------------------------------------Start automatic meshing workflow-------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# Environment creation
#--------------------------------------------------------------------------------------------------------------------------

for i in range(n_geometries):
    #Creates an output folder for specific case, based on the input folder name
    input_folder = osp.join(input_dir, input_list[i])
    output_name = f'0{i}_Result_{input_list[i]}'
    output_folder = osp.join(output_dir, output_name)


    #Grab the path of the geometry files
    inlet_path = osp.join(input_folder, osp.join(r'meshes', r'inlet.stl'))
    wall_path = osp.join(input_folder, osp.join(r'meshes', r'wall.stl'))
    outlet_path = osp.join(input_folder, osp.join(r'meshes', r'outlet.stl'))

    print('Created necessary files and directories')

    #Reading the files with pyvista
    inlet = pv.read(inlet_path)
    wall = pv.read(wall_path)
    outlet = pv.read(outlet_path)

    print('import of geometry done')



    #--------------------------------------------------------------------------------------------------------------------------
    # 3D-meshing algorithm
    #--------------------------------------------------------------------------------------------------------------------------

    #Cut the wall geometry after the aortic root
    wall_cut = cutting.main_cutter(inlet, wall, plot=show_plot)
    pv.save_meshio(osp.join(temp_dir, r'wall_cut.mesh'), wall_cut)

    #Create caps
    inlet_cap, outlet_cap = cap(wall_cut, plot=show_plot)
    pv.save_meshio(osp.join(temp_dir, r'inlet_cap.mesh'), inlet_cap)
    pv.save_meshio(osp.join(temp_dir, r'outlet_cap.mesh'), outlet_cap)

    #Combine cutted wall and inlet/outlet caps
    combined = (wall_cut + inlet_cap + outlet_cap).clean()
    print('Meshes succesfully combined')

    #Plot result of mesh combining. When edgetest gives an error, the structure has dublicate faces/nodes
    if show_plot:
        plt = pv.Plotter()
        plt.add_mesh(combined, style='wireframe')
        edgetest = combined.extract_feature_edges(boundary_edges=True, non_manifold_edges=True, manifold_edges=False, feature_edges=False)
        plt.show()

    pv.save_meshio(osp.join(temp_dir, r'combined_mesh.mesh'), combined)



    #Run remesh (takes predetermined internally defined file path as input, DON'T CHANGE)
    combined_remeshed = remesh.remesh_edge_detect(osp.join(temp_dir, r'combined_mesh.mesh'), osp.join(temp_dir, r'combined_mmg.mesh'), temp_dir, mmg_parameters, plot=show_plot)

    #triangulation step to make sure Tetgen only gets triangles as input
    combined_remeshed = combined_remeshed.extract_surface().triangulate()

    #Report quality
    quality_control.meshreport(combined_remeshed, 'Surface mesh quality report')

    #Make an initial 3D mesh from the combined mesh using TetGen
    combined_remeshed = combined_remeshed
    tetmesh = volume_mesh.tetgen(combined_remeshed, tetgen_parameters, plot=show_plot)

    #Plot bisection
    if show_plot:
        quality_control.clip_plot(tetmesh, 'Initial 3D mesh')

    #Report quality
    quality_control.meshreport(tetmesh, 'Initial 3D mesh quality report')

    #Create a .sol file for mmg3d
    volume_mesh.write_sol(tetmesh, wall_cut, mmg3d_sol_parameters, osp.join(temp_dir, r'initial_volume_mesh.sol'), plot=show_plot)

    #Save initial mesh
    tetmesh.point_data.clear()
    tetmesh.cell_data.clear()
    pv.save_meshio(osp.join(temp_dir, r'initial_volume_mesh.mesh'), tetmesh)

    #Refine 3D mesh with mmg3d
    tetmesh = volume_mesh.mmg3d(osp.join(temp_dir, r'initial_volume_mesh.mesh'), osp.join(temp_dir, r'mmg3d_mesh.mesh'), temp_dir, mmg3d_parameters, plot=show_plot)

    #Report quality
    report, report_text = quality_control.meshreport(tetmesh, f'Refined 3D mesh quality report {input_list[i]}')

    #Plot bad cells
    jac = report['jac']
    if show_plot:
        jac = report['jac']
        if np.any(jac<0.3):
            plt = pv.Plotter()
            plt.add_mesh(tetmesh, style='wireframe')
            plt.add_mesh(tetmesh.extract_cells(jac<0.3), color='red', show_edges=True)
            plt.add_text('Bad cells')
            plt.show()
        else: print('No bad cells')

    #Plot bisection
    if show_plot:
        quality_control.clip_plot(tetmesh, 'Refined 3D mesh')

    #Save 3D mesh
    tetmesh.save(osp.join(temp_dir, r'3D_output_mesh.vtk'))

    #Run qualification
    numcells = report['cells']
    aspect = report['aspect']
    if numcells > max_elements or any(aspect > max_aspect) or any(jac < min_jacobian): run = False
    else: run=True

    #Plot 3D_mesh
    if show_plot:
        if run:text='Final 3D mesh'
        else:text='Final 3D mesh - insufficient quality or too many nodes'
        tetmesh.plot(show_edges = True, text=text)

    #Write a log file of the mesh quality and continue to next geometry if quality is not sufficient.
    if run==False:
        print('Terminating, 3D mesh insufficient quality or too many nodes')
        print('See log files for quality rapport')
        log_dir = osp.join(file_dir, r'log\failed')
        ut.save_string_to_file(report_text, osp.join(log_dir, f'Qualityreport_failed_geometry_{input_list[i]}'))
        continue
    else:
        print('Quality is sufficient')
        log_dir = osp.join(file_dir, r'log')
        ut.save_string_to_file(report_text, osp.join(log_dir, f'Qualityreport_geometry_{input_list[i]}'))


    #----------------------------------------------------------------------------------------------------------------------------
    # Identification
    #----------------------------------------------------------------------------------------------------------------------------

    #Create the seeds for the surface identification based on the center points. INLET FIRST!, OUTLET SECOND!
    seeds = np.array([inlet_cap.points.mean(0),outlet_cap.points.mean(0)])

    #Detect the surfaces of the 3D mesh whilst keeping the original ID's, if there is a seed delivered
    if len(seeds) > 0:
        surface_identification = id.identify_surfaces(tetmesh, id_angle, seeds, show_plot)

    #Seperate the indentified surfaces in inlet/outlet/wall
    id_inlet = surface_identification[0]
    id_outlet = surface_identification[1]
    id_wall = surface_identification[2]

    print('Surface identification done')
    #Plot the identified surfaces for general overview
    if show_plot:
        plt = pv.Plotter()
        plt.add_mesh(id_inlet, show_edges = True, color = 'red')
        plt.add_mesh(id_outlet, show_edges = True, color = 'blue')
        plt.add_mesh(id_wall, show_edges = True, color = 'green')
        plt.show()

    #At this point meshing is succesfull and output files are written, geometry specific output folder is created.
    os.makedirs(output_folder, exist_ok=True)
    #----------------------------------------------------------------------------------------------------------------------------
    # Mapping
    #----------------------------------------------------------------------------------------------------------------------------

    #Perform the mapping of the velocity profiles on the inlet
    #Output is a point cloud on every inlet node with the respective velocity data and the amount of mapped velocity profiles

    velocity_mapped, n_maps = mapping.vel_mapping(vel_profile_dir, id_inlet, output_folder, intp_options, show_plot)
    #----------------------------------------------------------------------------------------------------------------------------
    # FEBio Creation
    #----------------------------------------------------------------------------------------------------------------------------

    #Create a solver compatible file based on the 3D-mesh and meshing parameters
    feb.xml_creator(tetmesh, id_inlet, id_outlet, id_wall, velocity_mapped, file_dir, output_folder)

#Deletes the temporary folder (might want to modify it to only delete the files)
shutil.rmtree(temp_dir)

#-----------------------------------------Start automatic simulation workflow-------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------
# FEBio Run
#----------------------------------------------------------------------------------------------------------------------------

#Grab the names of the output folders
output_list = os.listdir(output_dir)

#Path for FEBio solver executable
#FEBio_path = r"C:/Program Files/bin/febio4.exe" #Path voor Yarran
FEBio_path = r"C:/Program Files/FEBioStudio2/bin/febio4.exe" #Path voor normale mensen

#Run for every geometry a simulation
for sim in output_list:
    #Run FEBio
    sim_folder = osp.join(output_dir, sim)
    #Use the current
    FEBio_inputfile = osp.join(sim_folder, r'simulation.feb')
    subprocess.run([FEBio_path, FEBio_inputfile], check = True)

print('Done!')

