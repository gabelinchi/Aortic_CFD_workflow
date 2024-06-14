import numpy as np
import pyvista as pv
import os
import os.path as osp
import identification as id
from tkinter import Tk
from tkinter.filedialog import askopenfilenames
import cutting

def pp_wss():
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
        wall = split_surf[int(wall_index)]

        # Select second largest (inlet)
        inlet_index = np.argmax(np.delete(lengths, wall_index), axis=0)
        inlet = split_surf[int(inlet_index)].extract_surface()

        # Cut at horizontal point
        wall = cutting.post_cutter(inlet, wall, flip_norm=True)

        # Add array containing max. wss
        tensors = wall['fluid_stress'].reshape((wall['fluid_stress'].shape[0],3,3))
        eig = np.linalg.eigvals(tensors)
        min = eig.min(axis=1)
        max = eig.max(axis=1)
        wss = (max-min)/2
        wall['wss']=wss
        block.append(wall)

    for i in range(len(filenames)):
        wss=block[i]['wss']
        # Print statistics
        print('Histogram of wss for:', filenames[i])
        print('(Percentile: Value[Pa])') #Should be pascal (I think :) )
        print('99:', np.percentile(wss,99))
        print('98:', np.percentile(wss,98))
        print('96:', np.percentile(wss,96))
        print('94:', np.percentile(wss,94))
        print('92:', np.percentile(wss,92))
        print('90:', np.percentile(wss,90))
    
        # Plot result
        block[i].plot(scalars='wss', clim=[0,50], text=filenames[i])
    return

pp_wss(100)
