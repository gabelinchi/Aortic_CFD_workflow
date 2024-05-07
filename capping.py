import pyvista as pv
import numpy as np

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
    # Extract edges to cap
    edges = wall.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
    if plot==True:
        edges.plot()
    
    # Split edges into separate data
    edges = edges.connectivity('all')
    edges = edges.split_bodies()
    points0 = edges[0].points
    cells0 = np.hstack((np.array([len(points0)]),np.arange(len(points0))))
    points1 = edges[1].points
    print('edges', edges[0].cells)
    cells1 = np.hstack((np.array([len(points1)]),np.arange(len(points1))))
    print(cells0)
    #inlet = pv.PolyData(points0, cells0)
    #outlet = pv.PolyData(points1, cells1)
    inlet = edges[0].delaunay_2d()
    outlet = edges[1].delaunay_2d()
    if plot==True:
        plt = pv.Plotter()
        plt.add_mesh(inlet, label='Inlet', color='red', show_edges=True)
        plt.add_mesh(outlet, label='Outlet', color='blue', show_edges=True)
        plt.add_legend()
        plt.show()
    return(inlet, outlet)
    
#cap(wall, plot=True)