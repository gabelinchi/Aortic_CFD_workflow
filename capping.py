import pyvista as pv
import numpy as np
import utils as ut

#wall_path = "geometries\input\wall.stl"
#wall = pv.read(wall_path)



def cap(wall, plot=False):
    '''
    Creates caps for a vessel with a single inlet & outlet
    Inlet & outlet are assigned based on the number of nodes, cap with the most nodes becomes the inlet by default
    :wall : Pyvista Polydata/UnstructuredGrid of the wall geometry
    :plot : bool, show intermediate steps in plots
    :returns : (inlet, outlet)
    '''
    print('Start cap generation from cutted wall')

    # Extract edges to cap
    edges = wall.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
    
    #Plot extracted edges
    if plot==True:
        edges.plot()
    
    # Split edges into separate data
    edges = edges.connectivity('all')
    edges = edges.split_bodies()
    
    #Create a solid spiderweb-like surface mesh
    inlet = ut.make_spiderweb(edges[0])
    outlet = ut.make_spiderweb(edges[1])

    #Plot final result of capping function 
    if plot==True:
        plt = pv.Plotter()
        plt.add_mesh(inlet, label='Inlet', color='red', show_edges=True)
        plt.add_mesh(outlet, label='Outlet', color='blue', show_edges=True)
        plt.add_legend()
        plt.show()

    print('Finished cap generation')

    return(inlet, outlet)