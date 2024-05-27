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
import volume_remesh_mmg
from capping import cap
import identification as id
from tkinter import Tk
from tkinter.filedialog import askdirectory
import mapping
import xml.etree.ElementTree as ET
import subprocess
import febioxml as feb
from write_sol import write_sol

#Run FEBio
FEBio_path = r"C:\Program Files\FEBioStudio2\bin\febio4.exe"
#Use the current
FEBio_inputfile = osp.join(r'C:\Users\lmorr\Documents\TU\23-24\BEP\Git_repository\Aortic_CFD_workflow-3\output', r'simulation.feb')
subprocess.run([FEBio_path, FEBio_inputfile], check = True)