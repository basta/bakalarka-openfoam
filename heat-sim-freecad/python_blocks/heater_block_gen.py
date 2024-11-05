import os

import classy_blocks as cb
from classy_blocks.construct.flat.sketches.disk import DiskBase

DiskBase.core_ratio = 0.7  # Default is 0.8

START = 0
END = 5
RADIUS = 10

n_cells_radial = 5
n_cells_tangential = 5
boundary_layer_thickness = 2
axial_cell_size = 1

mesh = cb.Mesh()
inlet = cb.Cylinder([START, 0, 0], [END, 0, 0], [0, 0, RADIUS])
inlet.chop_radial(count=n_cells_radial, end_size=boundary_layer_thickness)
inlet.chop_axial(start_size=axial_cell_size, end_size=2*axial_cell_size)
inlet.chop_tangential(count=n_cells_tangential)

inlet.set_start_patch('inlet')
inlet.set_outer_patch('wall')
inlet.set_end_patch('outlet')


mesh.add(inlet)

mesh.write("../system/blockMeshDict", debug_path="debug.vtk")