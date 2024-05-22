import pyvista as pv
import numpy as np
import os
import os.path as osp



def get_bl_nodes(mesh, dist, dlow, dhigh, dir):
    print('start tagging')
    surf = mesh.extract_surface()
    closest_points = surf.find_closest_cell(mesh.points, return_closest_point=True)[1]
    d_exact = np.linalg.norm(mesh.points - closest_points, axis=1)
    density = np.where(d_exact > dist, dlow, dhigh)
    mesh["sol"] = density
    length = len(density)
    
    """ # get cell centroids
    plotmesh = mesh
    cells = plotmesh.cells.reshape(-1, 5)[:, 1:]
    cell_center = plotmesh.points[cells].mean(1)

    # extract cells below the 0 xy plane
    mask = cell_center[:, 2] < 0
    cell_ind = mask.nonzero()[0]
    subgrid = plotmesh.extract_cells(cell_ind)

    # advanced plotting
    plotter = pv.Plotter()
    plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True, scalars='sol')
    plotter.add_mesh(surf, 'r', 'wireframe')
    plotter.add_legend([[' Input Mesh ', 'r'], [' Tessellated Mesh ', 'black']])
    plotter.show() """

    header = 'MeshVersionFormatted \n2\n\nDimension 3\n\nSolAtVertices\n' + str(length) + '\n1 1\n'

    np.savetxt(dir, density.reshape(-1,1), header=header, footer='End', comments='', fmt='%f')
    return(mesh)

def write_sol():

    return