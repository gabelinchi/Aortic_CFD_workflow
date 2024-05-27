import pyvista as pv

def meshreport(mesh, text):
    # Calculate mesh stats and quality
    jac = mesh.compute_cell_quality(quality_measure='scaled_jacobian')['CellQuality']
    aspect = mesh.compute_cell_quality(quality_measure='aspect_ratio')['CellQuality']
    num_points = mesh.number_of_points
    num_cells = mesh.number_of_cells

    # Print report
    print('--------------------------------------------------------------')
    print(text)
    print('num points:', num_points)
    print('num_cells: ', num_cells)
    print('Mesh quality (mean scaled jacobian)')
    print('min:', jac.min())
    print('max:', jac.max())
    print('avg:', jac.mean())
    print('Mesh quality (aspect ratio)')
    print('min:', aspect.min())
    print('max:', aspect.max())
    print('avg:', aspect.mean())
    print('--------------------------------------------------------------')
    return{'jac':jac,'aspect':aspect, 'points':num_points, 'cells':num_cells}

def clip_plot(mesh, text):
    # Plots a geometry clipped along the z axis to show the inside mesh
    surf = mesh.extract_surface()
    clipped = mesh.clip('z', crinkle=True)

    plotter = pv.Plotter()
    plotter.add_mesh(clipped, 'lightgrey', lighting=True, show_edges=True)
    plotter.add_mesh(surf, 'r', 'wireframe')
    plotter.add_legend([['Mesh surface', 'r'], ['Clipped mesh', 'black']])
    plotter.add_text(text)
    plotter.show()