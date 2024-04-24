import numpy as np
import sys 
import os 
import os.path as osp 
from glob import glob 
import pyvista as pv

target_profile = r'C:\Users\lmorr\Documents\TU\23-24\BEP\Target_profile\inlet.stl'
read_file_input = r'C:\Users\lmorr\Documents\TU\23-24\BEP\Velocity_profiles\001_00.vtp'
read_file_output = r'C:\Users\lmorr\Documents\TU\23-24\BEP\Mapping_output\Mapped_velocity_profile_test_01.vtp'
plot_file = pv.read(target_profile)
plot_file.plot()

