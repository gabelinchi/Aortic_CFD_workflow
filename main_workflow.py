# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:16:40 2024

@author: lmorr
"""

#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv

wall_mesh = pv.read("wall.stl")
outlet_mesh = pv.read("outlet.stl")
inlet_mesh = pv.read("inlet.stl")
wall = wall_mesh.plot()

print(mesh)

#use an excessive amount of comments on everything.
#extra line of code changed
