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
    def pad(faces):
        num_rows = faces.shape[0]
        return(np.hstack((np.full((num_rows, 1), 3),faces)))
    
    # Extract edges to cap
    edges = wall.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False, feature_edges=False)
    if plot==True:
        edges.plot()
    
    # Split edges into separate data
    edges = edges.connectivity('all')
    edges = edges.split_bodies()
    #------------
    # One big poly

    # Extract points, cells
    points0 = edges[0].points
    points1 = edges[1].points
    cells0 = edges[0].cells
    cells1 = edges[1].cells
    def make_spiderweb(perimeter):
        points = perimeter.points
        cells = perimeter.cells
        cells = cells.reshape(-1, 3)[:,[1,2]]
        centerpoint = points.mean(0)
        points = np.vstack((points, centerpoint))
        cells = np.hstack((cells, np.full((len(cells), 1), len(points)-1)))
        cells = pad(cells)
        return pv.PolyData(points, cells)

    #------------
    inlet = make_spiderweb(edges[0])
    outlet = make_spiderweb(edges[1])
    if plot==True:
        plt = pv.Plotter()
        plt.add_mesh(inlet, label='Inlet', color='red', show_edges=True)
        plt.add_mesh(outlet, label='Outlet', color='blue', show_edges=True)
        plt.add_legend()
        plt.show()
    return(inlet, outlet)