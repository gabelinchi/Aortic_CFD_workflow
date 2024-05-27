import numpy as np
import pyvista as pv
import os
import os.path as osp
import identification as id


file_dir = osp.dirname(osp.realpath(__file__))
temp_dir = osp.join(file_dir, r'temp')
output_dir = osp.join(file_dir, r'output')
vtk_path = osp.join(output_dir, r'simulation.2.vtk')

vtk = pv.read(vtk_path)

split_surf = id.identify_surfaces(vtk, 30)
print(split_surf)

lengths = np.empty(len(split_surf))
for i in range(len(lengths)):
    lengths[i] = split_surf[i].number_of_points
wall_index = np.argmax(lengths, axis=0)
wall = split_surf[int(wall_index)]

