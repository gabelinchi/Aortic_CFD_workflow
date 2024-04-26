import pyvista as pv
inlet_mesh = pv.read("inlet.stl")

inlet = pv.plot (inlet_mesh)

