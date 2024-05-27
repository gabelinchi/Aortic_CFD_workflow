import sys
import os
import os.path as osp
from pathlib import Path
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import xml.etree.ElementTree as ET
import subprocess

def xml_creator(tetmesh, id_inlet, id_outlet, id_wall, velocity_profile, file_dir, output_dir):
    #constants:
    T = 4
    P = 0
    R = 8.31446
    Fc = 96485.3

    #fluidconstants
    materialtype = 'fluid'
    density = 1000
    k = 2200000                                     #bulkmodulus
    viscoustype = "Newtonian fluid"
    kappa = 1
    mu = 0.056

    #Boundarycondition zerofluidvelocity
    zerofluidvelocity_x = True
    zerofluidvelocity_y = True
    zerofluidvelocity_z = True

    #Loads
    velocity = -1                          #negative want moet aorta in, normals staan mesh uit (naar buiten)
    prescribe_nodal_velocities = True
    parabolic = True
    prescribe_rim_pressure= True

    #SimulationControl
    analysis = 'DYNAMIC'
    time_steps = 200                         #initial time steps, changes towards dtmax     
    step_size = 0.01                        #seconds
    plot_zero_state = 0
    plot_range = 0,-1
    plot_level = 'PLOT_MAJOR_ITRS'
    output_level = 'OUTPUT_MAJOR_ITRS'
    plot_stride = 1
    output_stride = 1
    adaptor_re_solve = 1
    time_stepper_type="default"
    max_retries = 20
    opt_iter = 50
    dtmin = 0
    dtmax = 0.1                           #max timestepsize
    aggressiveness = 0
    cutback = 0.5
    dtforce = 0 

    solvertype = 'fluid'
    symmetric_stiffness = 'non-symmetric'
    equation_scheme = 'staggered'
    equation_order = 'default'
    optimize_bw = 0
    lstol = 0.9
    lsmin = 0.01
    lsiter = 5
    max_refs = 5
    check_zero_diagonal = 0
    zero_diagonal_tol = 0
    force_partition = 0
    reform_each_time_step = 0
    reform_augment = 0
    diverge_reform = 0
    min_residual = 1e-20
    max_residual = 1e+20
    vtol = 0.001
    ftol = 0.001
    rhoi = 0
    etol = 0.01
    rtol = 0.001
    predictor = 0
    min_volume_ratio = 0
    qn_method_type ="Broyden"
    max_ups =50
    max_buffer_size = 0
    cycle_buffer = 1
    cmax = 100000



    #make nodes and elements arrays from mesh.vtk
    polydata = tetmesh 
    Meshnodes = polydata.points * 0.001
    ZeroElem = polydata.cells.reshape(-1,5)[:,[1,2,3,4]]
    MeshElem = []
    for num in ZeroElem:
        MeshElem.append(num + 1) 

    InletNodes = id_inlet.cell_data['original_connectivity']
    OutletNodes = id_outlet.cell_data['original_connectivity']
    WallNodes = id_wall.cell_data['original_connectivity']

    InletNodes = InletNodes + 1
    OutletNodes = OutletNodes + 1
    WallNodes = WallNodes + 1

    #get XML file
    tree = ET.parse(osp.join(file_dir, r'template_xml.feb'))
    root = tree.getroot()

    print('emptying .xml')
    #remove pre-existing nodes and elements
    Nodes_to_remove = []
    for element in root.iter('node'):               # Traverse the tree to find the element to edit
        Nodes_to_remove.append(element)             # Add Nodes to List
    Nodes = root.find("./Mesh/Nodes")
    for element in Nodes_to_remove:                 # Remove collected elements
        Nodes.remove(element)
    elements_to_remove = []
    for element in root.iter('elem'):               # Traverse the tree to find the specific element to edit
        elements_to_remove.append(element)          # Add elements to the list
    Elem = root.find("./Mesh/Elements")
    for element in elements_to_remove:
        Elem.remove(element)                       # Remove collected elements

    print('emptying surfaces')
    Surface = root.find("./Mesh/Surface[@name='FluidNormalVelocity1']")
    # If Surface element is found, remove all the tri3 elements under it
    if Surface is not None:
        for element in Surface.findall('tri3'):
            Surface.remove(element)
    Surface = root.find("./Mesh/Surface[@name='ZeroFluidVelocity1']")
    # If Surface element is found, remove all the tri3 elements under it
    if Surface is not None:
        for element in Surface.findall('tri3'):
            Surface.remove(element)
    Surface = root.find("./Mesh/Surface[@name='ZeroFluidDilatation3']")
    # If Surface element is found, remove all the tri3 elements under it
    if Surface is not None:
        for element in Surface.findall('tri3'):
            Surface.remove(element)
    print('emptied.xml')



    print('adding nodes and elements')
    # Add our own meshnodes to the file as XML elements
    nodes_list = []
    for i, array in enumerate(Meshnodes):
        element_name = "Node"
        element = ET.Element(element_name)
        element.text = ', '.join(map(str, array))   #add array as a string
        element.set('id', str(i+1))                 #set element id (python counts from 0, FEBio from 1)
        if i > 0:
            nodes_list[-1].tail = '\n\t\t\t'        #get the nodes in the right tree with tabs
        nodes_list.append(element)

    # Append XML elements to the 'Mesh/Nodes' element
    Element = root.find('.//Mesh/Nodes')
    if Element is not None:
        for node in nodes_list:
            Element.append(node)                    #add the nodes to the element <Nodes>
        nodes_list[-1].tail = '\n\t\t'              #get the closing </Nodes> in the right tree
    else:
        print("Element 'Mesh/Nodes' not found.") 


    # Add our own meshelements to the file as XML elements
    Elem_list = []
    for i, array in enumerate(MeshElem):
        element_name = "elem"
        element = ET.Element(element_name)
        element.text = ','.join(map(str, array))   #add array as a string
        element.set('id', str(i+1))                #set element id (python counts from 0, FEBio from 1)
        if i > 0:
            Elem_list[-1].tail = '\n\t\t\t'        #get the elements in the right tree with tabs
        Elem_list.append(element)

    # Append XML elements to the 'Mesh/Elements' element
    Element = root.find('.//Mesh/Elements')
    if Element is not None:
        for Elem in Elem_list:
            Element.append(Elem)                   #add the nodes to the element <Elements>
        Elem_list[-1].tail = '\n\t\t'              #get the closing </Elements> in the right tree
    else:
        print("Element 'Mesh/Elements' not found.") 


    # Add boundary conditions to the file as XML elements
    BC1_list = []
    for i, array in enumerate(WallNodes):
        element_name = "tri3"
        element = ET.Element(element_name)
        element.text = ','.join(map(str, array))   #add array as a string
        element.set('id', str(i+1))                #set element id (python counts from 0, FEBio from 1)
        if i > 0:
            BC1_list[-1].tail = '\n\t\t\t'        #get the elements in the right tree with tabs
        BC1_list.append(element)

    # Append XML elements to the 'Mesh/Surface' element
    Element = root.find('.//Mesh/Surface[@name = "ZeroFluidVelocity1"]')
    if Element is not None:
        for tri3 in BC1_list:
            Element.append(tri3)                   #add the nodes to the element <Elements>
        BC1_list[-1].tail = '\n\t\t'              #get the closing </Elements> in the right tree
    else:
        print("Element 'Mesh/Surface @name = ZeroFluidVelocity1' not found.") 

    # Add boundary conditions to the file as XML elements
    BC2_list = []
    for i, array in enumerate(OutletNodes):
        element_name = "tri3"
        element = ET.Element(element_name)
        element.text = ','.join(map(str, array))   #add array as a string
        element.set('id', str(i+1))                #set element id (python counts from 0, FEBio from 1)
        if i > 0:
            BC2_list[-1].tail = '\n\t\t\t'        #get the elements in the right tree with tabs
        BC2_list.append(element)

    # Append XML elements to the 'Mesh/Surface' element
    Element = root.find('.//Mesh/Surface[@name = "ZeroFluidDilatation1"]')
    if Element is not None:
        for tri3 in BC2_list:
            Element.append(tri3)                   #add the nodes to the element <Elements>
        BC2_list[-1].tail = '\n\t\t'              #get the closing </Elements> in the right tree
    else:
        print("Element 'Mesh/Surface @name = ZeroFluidDilatation' not found.") 

        
    # Add Fluidvelocity load to the file as XML elements
    load_list = []
    for i, array in enumerate(InletNodes):
        element_name = "tri3"
        element = ET.Element(element_name)
        element.text = ','.join(map(str, array))   #add array as a string
        element.set('id', str(i+1))                #set element id (python counts from 0, FEBio from 1)
        if i > 0:
            load_list[-1].tail = '\n\t\t\t'        #get the elements in the right tree with tabs
        load_list.append(element)

    # Append XML elements to the 'Mesh/Surface' element
    Element = root.find('.//Mesh/Surface[@name = "FluidNormalVelocity1"]')
    if Element is not None:
        for face in load_list:
            Element.append(face)                   #add the nodes to the element <Elements>
        load_list[-1].tail = '\n\t\t'              #get the closing </Elements> in the right tree
    else:
        print("Element load not found.") 
    print('added nodes')


    # Add velocity profile to the file as XML elements
    profile_list = []
    for i, array in enumerate(velocity_profile):
        element_name = "face"
        element = ET.Element(element_name)
        element.text = ','.join(map(str, array))   #add array as a string
        element.set('lid', str(i+1))                #set element id (python counts from 0, FEBio from 1)    
        if i > 0:
            profile_list[-1].tail = '\n\t\t\t'        #get the elements in the right tree with tabs
        profile_list.append(element)
    # Append XML elements to the 'MeshData/SurfaceData' element
    Element = root.find('.//MeshData/SurfaceData[@name = "velocityprofile"]')
    if Element is not None:
        for vel in profile_list:
            Element.append(vel)                   #add the nodes to the element <Elements>
        profile_list[-1].tail = '\n\t\t'              #get the closing </Elements> in the right tree
    else:
        print("Element SurfaceData[@name = 'velocityprofile'] not found.")
    
    print('added velocity profile')

    #Edit constants
    def Parent1(parent, variable, value):
        element_to_edit = root.find(f"./{parent}/{variable}") # Find the XML element based on the provided 'place' and 'variable'
        if element_to_edit is not None:                     # Check if the element is found
            element_to_edit.text = value              # Set the text content of the element to 'kwaliteit'
        else:
            """element_to_edit = root.find(f"./{parent}")   #weet niet zeker of dit gaat werken omdat de volgorde van variablen uit zou kunnen maken
            element_name = ET.Element(variable)
            element_to_edit.tail = '\n\t\t'
            element_to_edit.append(element_name)
            element_to_edit.text = value """
            
            print(f"Element '{variable}' under '{parent}' not found.")
        return

    def Parent2(pparent, parent, variable, value):
        element_to_edit = root.find(f"./{pparent}/{parent}/{variable}") # Find the XML element based on the provided 'parent' and 'variable'
        if element_to_edit is not None:                     # Check if the element is found
            element_to_edit.text = value 
        else:
                print(f"Element '{variable}' under '{parent}' not found.")
        return

    def Parent3(ppparent,pparent, parent, variable, value):
        element_to_edit = root.find(f"./{ppparent}/{pparent}/{parent}/{variable}") # Find the XML element based on the provided 'parent' and 'variable'
        if element_to_edit is not None:                     # Check if the element is found
            element_to_edit.text = value 
        else:
                print(f"Element '{variable}' under '{parent}' not found.")
        return

    def Parent4(pppparent, ppparent,pparent, parent, variable, value):
        element_to_edit = root.find(f"./{pppparent}/{ppparent}/{pparent}/{parent}/{variable}") # Find the XML element based on the provided 'parent' and 'variable'
        if element_to_edit is not None:                     # Check if the element is found
            element_to_edit.text = value 
        else:
                print(f"Element '{variable}' under '{parent}' not found.")
        return
    #single parent
    Parent1('Control', 'analysis', str(analysis))
    Parent1('Control', 'time_steps', str(time_steps))
    Parent1('Control', 'step_size', str(step_size))
    Parent1('Control', 'plot_zero_state', str(plot_zero_state))

    plot_range_str = ','.join(str(x) for x in plot_range)
    Parent1('Control', 'plot_range', plot_range_str)
    Parent1('Control', 'plot_level', str(plot_level))
    Parent1('Control', 'output_level', str(output_level))
    Parent1('Control', 'plot_stride', str(plot_stride))
    Parent1('Control', 'output_stride', str(output_stride))
    Parent1('Control', 'adaptor_re_solve', str(adaptor_re_solve))


    #double parent
    Parent2('Control', 'time_stepper', 'max_retries', str(max_retries))
    Parent2('Control', 'time_stepper', 'opt_iter', str(opt_iter))
    Parent2('Control', 'time_stepper', 'dtmin', str(dtmin))
    Parent2('Control', 'time_stepper', 'aggressiveness', str(aggressiveness))
    Parent2('Control', 'time_stepper', 'cutback', str(cutback))
    Parent2('Control', 'time_stepper', 'dtforce', str(dtforce))
    Parent2('Control', 'time_stepper', 'max_retries', str(max_retries))
    Parent2('Control', 'time_stepper', 'opt_iter', str(opt_iter))
    Parent2('Control', 'time_stepper', 'dtmin', str(dtmin))
    Parent2('Control', 'time_stepper', 'aggressiveness', str(aggressiveness))
    Parent2('Control', 'time_stepper', 'cutback', str(cutback))
    Parent2('Control', 'time_stepper', 'dtforce', str(dtforce))
    Parent2('Material', 'material', 'density', str(density))
    Parent2('Material', 'material', 'k', str(k))

    Parent2('Control','solver', 'symmetric_stiffness', str(symmetric_stiffness))
    Parent2('Control', 'solver','equation_scheme', str(equation_scheme))
    Parent2('Control','solver', 'equation_order', str(equation_order))
    Parent2('Control','solver', 'optimize_bw', str(optimize_bw))
    Parent2('Control', 'solver','lstol', str(lstol))
    Parent2('Control', 'solver','lsmin', str(lsmin))
    Parent2('Control','solver', 'lsiter', str(lsiter))
    Parent2('Control','solver', 'max_refs', str(max_refs))
    Parent2('Control','solver', 'check_zero_diagonal', str(check_zero_diagonal))
    Parent2('Control','solver', 'zero_diagonal_tol', str(zero_diagonal_tol))
    Parent2('Control','solver', 'force_partition', str(force_partition))
    Parent2('Control','solver', 'reform_each_time_step', str(reform_each_time_step))
    Parent2('Control','solver', 'reform_augment', str(reform_augment))
    Parent2('Control','solver', 'diverge_reform', str(diverge_reform))
    Parent2('Control','solver', 'min_residual', str(min_residual))
    Parent2('Control','solver', 'max_residual', str(max_residual))
    Parent2('Control','solver', 'vtol', str(vtol))
    Parent2('Control','solver', 'ftol', str(ftol))
    Parent2('Control','solver', 'rhoi', str(rhoi))
    Parent2('Control','solver', 'etol', str(etol))
    Parent2('Control','solver', 'rtol', str(rtol))
    Parent2('Control','solver', 'predictor', str(predictor))
    Parent2('Control','solver', 'min_volume_ratio', str(min_volume_ratio))
    Parent2('Control','solver', 'symmetric_stiffness', str(symmetric_stiffness))
    Parent2('Control', 'solver','equation_scheme', str(equation_scheme))
    Parent2('Control','solver', 'equation_order', str(equation_order))
    Parent2('Control','solver', 'optimize_bw', str(optimize_bw))
    Parent2('Control', 'solver','lstol', str(lstol))
    Parent2('Control', 'solver','lsmin', str(lsmin))
    Parent2('Control','solver', 'lsiter', str(lsiter))
    Parent2('Control','solver', 'max_refs', str(max_refs))
    Parent2('Control','solver', 'check_zero_diagonal', str(check_zero_diagonal))
    Parent2('Control','solver', 'zero_diagonal_tol', str(zero_diagonal_tol))
    Parent2('Control','solver', 'force_partition', str(force_partition))
    Parent2('Control','solver', 'reform_each_time_step', str(reform_each_time_step))
    Parent2('Control','solver', 'reform_augment', str(reform_augment))
    Parent2('Control','solver', 'diverge_reform', str(diverge_reform))
    Parent2('Control','solver', 'min_residual', str(min_residual))
    Parent2('Control','solver', 'max_residual', str(max_residual))
    Parent2('Control','solver', 'vtol', str(vtol))
    Parent2('Control','solver', 'ftol', str(ftol))
    Parent2('Control','solver', 'rhoi', str(rhoi))
    Parent2('Control','solver', 'etol', str(etol))
    Parent2('Control','solver', 'rtol', str(rtol))
    Parent2('Control','solver', 'predictor', str(predictor))
    Parent2('Control','solver', 'min_volume_ratio', str(min_volume_ratio))

    Parent2('Boundary', 'bc', 'wx_dof', str(int(zerofluidvelocity_x)))
    Parent2('Boundary', 'bc', 'wy_dof', str(int(zerofluidvelocity_y)))
    Parent2('Boundary', 'bc', 'wz_dof', str(int(zerofluidvelocity_z)))

    #Parent2('Loads', 'surface_load', 'velocity', str(velocity))
    #Parent2('Loads', 'surface_load', 'prescribe_nodal_velocities', str(int(prescribe_nodal_velocities)))
    #Parent2('Loads', 'surface_load', 'parabolic', str(int(parabolic)))
    #Parent2('Loads', 'surface_load', 'parabolic', str(int(prescribe_rim_pressure)))

    #triple parent
    Parent3('Material', 'material', 'viscous', 'kappa', str(kappa))
    Parent3('Material', 'material', 'viscous', 'mu', str(mu))

    """ timestepper  = root.find('./Control/time_stepper')      #dit moet nog onder step/step1
    timestepper.set('type',str(time_stepper_type))
    solver  = root.find('./Control/solver')
    solver.set('type',str(solvertype))
    qnmethod  = root.find('./Control/solver/qn_method')
    qnmethod.set('type',str(qn_method_type))
    material  = root.find('./Material/material')
    material.set('type',str(materialtype))
    viscous  = root.find('./Material/material/viscous')
    viscous.set('type',str(viscoustype))
    """
    #create FEBio file
    tree.write(osp.join(output_dir, r'simulation.feb'), encoding='ISO-8859-1', xml_declaration=True,)

    return
