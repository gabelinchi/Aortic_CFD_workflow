import numpy as np
import pyvista as pv
import os
import os.path as osp
import identification as id
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
import cutting
import utils as ut

def wss_cut_and_hist(maxrange=50):
    '''
    Asks the user to select vtk files to process. For each file, geometry is cut where the artery is horizontal,
    wss data is added and the result is plotted. Also prints a histogram of wss for each .vtk

    :maxrange   : float, maximum value for the color plot
    :return     :pyvista MultiBlock with cut geometry
    '''
    # Read vtk
    vtk_path = askopenfilenames(title='Select VTK file(s)')
    block = pv.MultiBlock()
    filenames = []
    pos = 'xy'
    filecount = 1

    for filepath in vtk_path:
        vtk = pv.read(filepath).scale([1000,1000,1000]) #Scale the model back to mm
        filename = os.path.basename(filepath)
        filenames.append(filename)

        # Identify surfaces & select largest (wall)
        split_surf = id.identify_surfaces(vtk, 30)
        areas = np.empty(len(split_surf))
        for i in range(len(areas)):
            areas[i] = split_surf[i].area
        wall_index = np.argmax(areas, axis=0)
        wall = split_surf[int(wall_index)]

        # Select second largest (inlet)
        areas[wall_index] = 0
        inlet_index = np.argmax(areas, axis=0)
        inlet = split_surf[int(inlet_index)].extract_surface()

        # Cut at horizontal point
        wall = cutting.post_cutter(inlet, wall, flip_norm=True).extract_surface()

        """# Add array containing max. wss
        tensors = wall['fluid_stress'].reshape((wall['fluid_stress'].shape[0],3,3))
        eig = np.linalg.eigvals(tensors)
        min = eig.min(axis=1)
        max = eig.max(axis=1)
        wss = (max-min)/2
        wall['wss']=wss
        block.append(wall) """

        # Add array containing max. wss
        tensors = wall['fluid_stress'].reshape((wall['fluid_stress'].shape[0],3,3))
        wall = wall.compute_normals(point_normals=False)
        Tn = np.einsum('ijk,ik->ij', tensors, wall['Normals']) # Matrices * normals
        sigma_n = np.sum(Tn * wall['Normals'], axis=1) # Stress vectors * normals
        wss = abs((np.linalg.norm(Tn, axis=1)**2 - sigma_n**2))**0.5 # Extremely small values produce rounding errors
        wall['wss']=wss
        block.append(wall)
        print('File nr.', filecount,'of',len(vtk_path), 'processed')
        filecount += 1

    for i in range(len(filenames)):
        wss=block[i]['wss']

        #Create report text
        report_text = '''
--------------------------------------------------------------
Histogram of wss for: {0}
'(Percentile: Value[Pa])'
99: {1}
98: {2}
96: {3}
94: {4}
92: {5}
90: {6}
--------------------------------------------------------------
'''.format(filenames[i], np.percentile(wss,99), np.percentile(wss,98), np.percentile(wss,96), np.percentile(wss,94), np.percentile(wss,92), np.percentile(wss,90) )

        print(report_text)

        # Save to .txt
        ut.save_string_to_file(report_text, osp.join(osp.dirname(vtk_path[i]), f'WSS_{filenames[i][:-4]}.txt'))
    
        # Plot result
        #block[i].plot(scalars='wss', clim=[0, maxrange], text=filenames[i], scalar_bar_args={'title': 'wss [Pa]'},
        #              camera_position=pos, screenshot=osp.join(osp.dirname(vtk_path[i]), f'WSS_plot_{filenames[i][:-4]}.png'))

        # Plot and save camera position
        cam = pv.Plotter()
        cam.add_mesh(block[i], clim=[0, maxrange], scalar_bar_args={'title': 'wss [Pa]'}, scalars='wss')
        cam.add_text(filenames[i])
        cam.camera_position = pos
        cam.show()
        pos = cam.camera_position


        # Save screenshot
        plt = pv.Plotter(off_screen=True)
        plt.add_mesh(block[i], clim=[0, maxrange], scalar_bar_args={'title': 'wss [Pa]'}, scalars='wss')
        plt.add_text(filenames[i])
        plt.camera_position = pos
        plt.show(screenshot=osp.join(osp.dirname(vtk_path[i]), f'WSS_plot_{filenames[i][:-4]}.png'))
    return block

def wss_simple(maxrange=50):
    '''
    Asks the user to select vtk files to process. For each file, wss data is added and the result is plotted.

    :maxrange   : float, maximum value for the color plot
    :return     :pyvista MultiBlock with cut geometry
    '''
    # Read vtk
    vtk_path = askopenfilenames(title='Select VTK file(s)')
    block = pv.MultiBlock()
    filenames = []

    for filepath in vtk_path:
        vtk = pv.read(filepath).scale([1000,1000,1000])
        filename = os.path.basename(filepath)
        filenames.append(filename)

        # Identify surfaces & select largest (wall)
        split_surf = id.identify_surfaces(vtk, 30)
        lengths = np.empty(len(split_surf))
        for i in range(len(lengths)):
            lengths[i] = split_surf[i].number_of_points
        wall_index = np.argmax(lengths, axis=0)
        wall = split_surf[int(wall_index)].extract_surface()

        # Add array containing max. wss
        tensors = wall['fluid_stress'].reshape((wall['fluid_stress'].shape[0],3,3))
        wall = wall.compute_normals(point_normals=False)
        Tn = np.einsum('ijk,ik->ij', tensors, wall['Normals']) # Matrices * normals
        sigma_n = np.sum(Tn * wall['Normals'], axis=1) # Stress vectors * normals
        wss = abs((np.linalg.norm(Tn, axis=1)**2 - sigma_n**2))**0.5 # Extremely small values produce rounding errors
        wall['wss']=wss
        block.append(wall)
    
    # Plot results
    for i in range(len(block)):
        block[i].plot(scalars='wss', clim=[0,maxrange], text=filenames[i], scalar_bar_args={'title': 'wss [Pa]'})
    return

wss_cut_and_hist(100)