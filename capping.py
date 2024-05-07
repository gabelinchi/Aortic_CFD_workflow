import numpy as np
import pyvista as pv
import cutting as cutting

inlet_path = "geometries\input\inlet.stl"
wall_path = "geometries\input\wall.stl"

inlet = pv.read(inlet_path)
wall = pv.read(wall_path)

wall_cut = cutting.main_cutter(inlet, wall, plot=True)

def cap(wall, plot=False):
    '''
    Creates caps for a vessel with a single inlet & outlet
    Inlet & outlet are assigned based on the number of nodes, cap with the most nodes becomes the inlet by default
    :wall : Pyvista Polydata/UnstructuredGrid of the wall geometry
    :plot : bool, show intermediate steps in plots
    '''
    # Extract edges to cap
    edges = wall.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
    if plot==True:
        edges.plot()
    
    # Split edges into separate data
    edges = edges.connectivity('all')
    edges = edges.split_bodies()
    
    inlet = edges[0].delaunay_2d()
    outlet = edges[1].delaunay_2d()
    if plot==True:
        plt = pv.Plotter()
        plt.add_mesh(inlet, label='Inlet', color='red', show_edges=True)
        plt.add_mesh(outlet, label='Outlet', color='blue', show_edges=True)
        plt.add_legend()
        plt.show()
    return(inlet, outlet)

cap(wall_cut, plot=True)
    
