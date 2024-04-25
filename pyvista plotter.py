import numpy as np
import sys 
import os 
import os.path as osp 
from glob import glob 
import pyvista as pv

aorta_geometry = r'C:\Users\lmorr\Documents\TU\23-24\BEP\case_0163\meshes\wall.stl'
target_profile = r'C:\Users\lmorr\Documents\TU\23-24\BEP\Target_profile\inlet.stl'
read_file_input = r'C:\Users\lmorr\Documents\TU\23-24\BEP\Velocity_profiles\001_00.vtp'
read_file_output = r'C:\Users\lmorr\Documents\TU\23-24\BEP\Mapping_output\Mapped_velocity_profile_00.vtp'
plot_file = pv.read(read_file_input)
plot_file.plot(show_edges = True)

