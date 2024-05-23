import pyvista as pv
import numpy as np
import os
import os.path as osp



def get_bl_nodes(mesh, surf, dist, dlow, dhigh, dir, plot=False):
    print('start tagging')

    closest_points = surf.find_closest_cell(mesh.points, return_closest_point=True)[1]
    d_exact = np.linalg.norm(mesh.points - closest_points, axis=1)
    density = np.where(d_exact > dist, dlow, dhigh)
    mesh["sol"] = density
    length = len(density)
    
    if plot:
        clipped = mesh.clip('z', crinkle=True)
        plotter = pv.Plotter()
        plotter.add_mesh(clipped, 'lightgrey', lighting=True, show_edges=True, scalars='sol')
        plotter.add_mesh(surf, 'r', 'wireframe')
        plotter.add_legend([[' Input Mesh ', 'r'], [' Tessellated Mesh ', 'black']])
        plotter.add_text('Initial 3D mesh with local tags')
        plotter.show()

    header = 'MeshVersionFormatted \n2\n\nDimension 3\n\nSolAtVertices\n' + str(length) + '\n1 1\n'

    np.savetxt(dir, density.reshape(-1,1), header=header, footer='End', comments='', fmt='%f')
    return(mesh)

def write_sol():

    return