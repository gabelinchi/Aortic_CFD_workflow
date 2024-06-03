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
    #-------------------------------------------------------------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------------------------------------------------------------
    
    #outputs
    output_displacement = False
    output_fluid_pressure = False
    output_nodal_fluid_velocity = False
    output_fluid_stress = True
    output_fluid_velocity = True
    output_fluid_acceleration = False
    output_fluid_vorticity = False
    output_fluid_rate_of_deformation = False
    output_fluid_dilatation = False
    output_fluid_volume_ratio = False

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
    time_steps = 400                         #initial time steps, changes towards dtmax     
    step_size = 0.005                        #seconds
    plot_zero_state = 0
    plot_range = 0,-1
    plot_level = 'PLOT_MAJOR_ITRS'
    output_level = 'OUTPUT_MAJOR_ITRS'
    plot_stride = 1
    output_stride = 1
    adaptor_re_solve = 1
    time_stepper_type="default"
    max_retries = 5                        # make smaller for quicker converge
    opt_iter = 25
    dtmin = 0
    dtmax = 0.05                           #max timestepsize
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

    #loadcontroller
    
    interpolate = 'LINEAR'                    #can be'STEP',  'LINEAR' or 'SMOOTH'
    extend = "CONSTANT"                     #sws niet nodig, overal gedefineerd

    #-------------------------------------------------------------------------------------------------------------------------
    # Adding nodes and elements
    # ------------------------------------------------------------------------------------------------------------------------
    #make nodes and elements arrays from mesh.vtk
    polydata = tetmesh 
    Meshnodes = polydata.points * 0.001                         #convert to meters
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
    Surface = root.find("./Mesh/Surface[@name='ZeroFluidVelocity2']")
    # If Surface element is found, remove all the tri3 elements under it
    if Surface is not None:
        for element in Surface.findall('tri3'):
            Surface.remove(element)
    Surface = root.find("./Mesh/Surface[@name='ZeroFluidDilatation1']")
    # If Surface element is found, remove all the tri3 elements under it
    if Surface is not None:
        for element in Surface.findall('tri3'):
            Surface.remove(element)
    loadcontroller = root.find("./LoadData/load_controller/points")
    # If loadcontroller element is found, remove all the loadpoints under it
    if loadcontroller is not None:
        for element in Surface.findall('pt'):
            loadcontroller.remove(element)

    
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

    

    #-------------------------------------------------------------------------------------------------------------------------
    # Add loadcurves
    # ------------------------------------------------------------------------------------------------------------------------
    #amplification factor of the load as loadcurve
    def loadcurve(interpolate,extend, time_steps, step_size,cycle, root, t_start, t_end, n):
        loadcurve_list = np.zeros((time_steps, 2))
        
        dt = t_end-t_start
        if n == 1:
            for i in range(time_steps):
                t = (i) * step_size                                       #timepoint at each step
                if (t % cycle) < (t_start-dt):                                         #amplification factor = 0 for time before t start
                    f = 0
                elif (t_start+cycle-dt)<= (t % cycle) < cycle:                              #amplification factor = 1 for time between t start and t end
                    f = ((t % cycle)-t_start-cycle+dt)/dt
                elif t_start<= (t % cycle) < (t_start+dt):                              #amplification factor = 1 for time between t start and t end
                    f = (-(t % cycle)+t_start+dt)/dt
                else:                                                   #amplification factor = 0 for time after t end
                    f = 0
                loadcurve_list[i] = np.array([t, f])                         #add function f(t) as points

        else: 
            for i in range(time_steps):
                t = (i) * step_size                                       #timepoint at each step
                if (t % cycle) < (t_start-dt):                                         #amplification factor = 0 for time before t start
                    f = 0
                elif (t_start-dt)<= (t % cycle) < t_start:                              #amplification factor = 1 for time between t start and t end
                    f = ((t % cycle)-t_start+dt)/dt
                elif t_start<= (t % cycle) < (t_start+dt):                              #amplification factor = 1 for time between t start and t end
                    f = (-(t % cycle)+t_start+dt)/dt
                else:                                                   #amplification factor = 0 for time after t end
                    f = 0
                loadcurve_list[i] = np.array([t, f])                         #add function f(t) as points

        load_main = root.find('.//LoadData')
        new_controller = ET.SubElement(load_main, 'load_controller')
        new_controller.set('id', f"{n}")
        new_controller.set('type', "loadcurve")
        new_controller.tail = '\n\t\t\t'

        interpolate_lc = ET.SubElement(new_controller, 'interpolate')
        interpolate_lc.text = str(interpolate)
        interpolate_lc.tail = '\n\t\t\t'

        extend_lc = ET.SubElement(new_controller, 'extend')
        extend_lc.text = str(extend)
        extend_lc.tail = '\n\t\t\t'

        points_lc = ET.SubElement(new_controller, 'points')
        points_lc.tail = '\n\t\t\t'

            
        load_controller = root.find(f'.//LoadData/load_controller[@id="{n}"]')
        points_element = load_controller.find('points')
        
        for array in loadcurve_list:                     # Create XML elements for each point in the load curve:
            element_name = "pt"
            element = ET.Element(element_name)              #make the string an .xml element
            t, f = array
            element.text = str(t) + ',' + str(f)
            element.tail= '\n\t\t\t\t'
            points_element.append(element)                  #add points to the xml

    #one loadcurve for each velocityprofile snapshot:
    hr = 60             #bpm
    cycle = hr/60       #time in seconds
    profiles = len(velocity_profile)
    for i in range(profiles):
        loadcurve(interpolate, extend, time_steps, step_size,cycle, root, (cycle/profiles)*i, (cycle/profiles)*(i+1), i +1)
    
    print('finished loadcontroller')

    #-------------------------------------------------------------------------------------------------------------------------
    # Boundary conditions
    # ------------------------------------------------------------------------------------------------------------------------

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

    # Append XML elements to the ZeroFluidVelocity (=wall) surface
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

    # Append XML elements to the ZeroFluidDilatation (=outlet) surface
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

    # Append XML elements to the FluidNormalVelocity surface (=inlet) element
    Element = root.find('.//Mesh/Surface[@name = "FluidNormalVelocity1"]')
    if Element is not None:
        for tri3 in load_list:
            Element.append(tri3)                   #add the nodes to the element <Elements>
        load_list[-1].tail = '\n\t\t'              #get the closing </Elements> in the right tree
    else:
        print("Element load not found.") 
    print('added nodes')

    #-------------------------------------------------------------------------------------------------------------------------
    # Velocity profile
    # ------------------------------------------------------------------------------------------------------------------------

    def vel_profile_adding(velocity_profile, frame_number):
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

        profile_main = root.find('.//MeshData')
        new_profile = ET.SubElement(profile_main, 'SurfaceData')
        profile_name = "velocityprofile" + str(frame_number)
        new_profile.set('name', profile_name)
        new_profile.set('data_type', "vec3")
        new_profile.set('surface', "FluidNormalVelocity1")
        new_profile.tail = '\n\t\t\t'


        Element = root.find(f'.//MeshData/SurfaceData[@name ="{profile_name}"]')
        if Element is not None:
            for vel in profile_list:
                Element.append(vel)                   #add the nodes to the element <Elements>
            profile_list[-1].tail = '\n\t\t'              #get the closing </Elements> in the right tree
        else:
            print("Element SurfaceData[@name = 'velocityprofile'] not found.")
    
    for i in range(len(velocity_profile)):
        vel_profile_adding(velocity_profile[i], i + 1)


    print('added velocity profile')

    #-------------------------------------------------------------------------------------------------------------------------
    # Add loads
    # ------------------------------------------------------------------------------------------------------------------------
    def add_loads(num_frames):
        load_main = root.find('.//Loads')
        new_load = ET.SubElement(load_main, 'surface_load')
        new_load.set('surface', "FluidNormalVelocity1")
        new_load.set('type', "fluid velocity")
        new_load.tail = '\n\t\t\t'

        velocity_l = ET.SubElement(new_load, 'velocity')
        velocity_l.set('type', 'map')
        velocity_l.set('lc', str(num_frames))
        velocity_l.text = 'velocityprofile' + str(num_frames)
        velocity_l.tail = '\n\t\t\t\t'

    for i in range(len(velocity_profile)):          #for the amount of velocity profiles, add different loads
        add_loads(i + 1)
    #-------------------------------------------------------------------------------------------------------------------------
    # Edit constants
    # ------------------------------------------------------------------------------------------------------------------------

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
    Parent2('Control', 'time_stepper', 'dtmax', str(dtmax))
    Parent2('Control', 'time_stepper', 'aggressiveness', str(aggressiveness))
    Parent2('Control', 'time_stepper', 'cutback', str(cutback))
    Parent2('Control', 'time_stepper', 'dtforce', str(dtforce))
    Parent2('Control', 'time_stepper', 'max_retries', str(max_retries))
    Parent2('Control', 'time_stepper', 'opt_iter', str(opt_iter))
    Parent2('Control', 'time_stepper', 'dtmin', str(dtmin))
    Parent2('Control', 'time_stepper', 'dtmax', str(dtmax))
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
    
    #Parent2('LoadData', 'load_controller', 'interpolate', str(interpolate))
    #Parent2('LoadData', 'load_controller', 'extend', str(extend))

    #triple parent
    Parent3('Material', 'material', 'viscous', 'kappa', str(kappa))
    Parent3('Material', 'material', 'viscous', 'mu', str(mu))
    Parent3('Control', 'solver', 'qn_method', 'max_ups', str(max_ups))
    Parent3('Control', 'solver', 'qn_method', 'max_buffer_size', str(max_buffer_size))
    Parent3('Control', 'solver', 'qn_method', 'cycle_buffer', str(cycle_buffer))
    Parent3('Control', 'solver', 'qn_method', 'cmax', str(cmax))


    timestepper  = root.find('./Control/time_stepper')      
    timestepper.set('type',str(time_stepper_type))
    solver  = root.find('./Control/solver')
    solver.set('type',str(solvertype))
    qnmethod  = root.find('./Control/solver/qn_method')
    qnmethod.set('type',str(qn_method_type))
    material  = root.find('./Material/material')
    material.set('type',str(materialtype))
    viscous  = root.find('./Material/material/viscous')
    viscous.set('type',str(viscoustype))
    

    #-------------------------------------------------------------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------------------------------------------------------------
    def outputs(output):
        output_main = root.find('.//Output/plotfile[@type="febio"]')
        new_output = ET.SubElement(output_main, 'var')
        new_output.set('type', f"{output}")
        new_output.tail = '\n\t\t'
        return

    if output_displacement==True:
        outputs('displacement')
    if output_fluid_pressure == True:
        outputs('fluid pressure')
    if output_nodal_fluid_velocity == True:
        outputs('nodal fluid velocity')
    if output_fluid_stress == True:
        outputs('fluid stress')
    if output_fluid_velocity == True:
        outputs('fluid velocity')
    if output_fluid_acceleration == True:
        outputs('fluid acceleration')
    if output_fluid_vorticity == True:
        outputs('fluid vorticity')
    if output_fluid_rate_of_deformation == True:
        outputs('fluid rate of deformation')
    if output_fluid_dilatation == True:
        outputs('fluid dilatation')
    if output_fluid_volume_ratio == True:
        outputs('fluid volume ratio')


    #create FEBio file
    tree.write(osp.join(output_dir, r'simulation.feb'), encoding='ISO-8859-1', xml_declaration=True,)

    return
