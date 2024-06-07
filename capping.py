#import modules
import sys
import os.path as osp
from glob import glob
import numpy as np
import pyvista as pv
import tetgen as tet
import utils as ut



def cap(wall, inlet_center, outlet_center, plot=False):
    '''
    Creates caps for a vessel with a single inlet & outlet
    Inlet & outlet are assigned based on the number of nodes, cap with the most nodes becomes the inlet by default
    :wall : Pyvista Polydata/UnstructuredGrid of the wall geometry
    :plot : bool, show intermediate steps in plots
    :returns : (inlet, outlet)
    '''
    print('Start cap generation from cutted wall')
    # Clear data arrays (somehow breaks the function if not done)
    wall.clear_data()
    
    # Extract edges to cap
    edges = wall.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
    
    #Plot extracted edges
    if plot==True:
        edges.plot(color='red', text='Perimeters to cap')
    
    # Split edges into separate data
    inlet_edges = edges.connectivity('closest', inlet_center)
    outlet_edges = edges.connectivity('closest', outlet_center)

    #Create a solid spiderweb-like surface mesh
    inlet = ut.make_spiderweb(pv.UnstructuredGrid(inlet_edges))
    outlet = ut.make_spiderweb(pv.UnstructuredGrid(outlet_edges))

    #Plot final result of capping function 
    if plot==True:
        plt = pv.Plotter()
        plt.add_mesh(inlet, label='Inlet', color='red', show_edges=True)
        plt.add_mesh(outlet, label='Outlet', color='blue', show_edges=True)
        plt.add_legend()
        plt.add_text('Generated caps')
        plt.show()

    print('Finished cap generation')

    return(inlet, outlet)