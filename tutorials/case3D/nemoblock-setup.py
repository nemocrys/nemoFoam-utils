"""Mesh for UDV-test-setup simulation"""
from nemoblock import *
import shutil
import numpy as np

mesh = Mesh()
####################
# For mesh optimization

one_mesh_only = False
sim_dir = "."
dimension = 3

r_melt = 0.1
h_melt = 0.1
r_crys = 0.05
####################
# Mes structure  (left: symmetry axis)
# fmt: off
#
# coordinate system:
#
# ^ z-direction
# |
# |
# |------->
#   r-direction
#
#       
# .__ ___ ___  
# |  |   |   |
# .  |   |   |
# |c1| r1| r3|
# .__|___|___|
# |  |   |   |
# .c2| r2| r4|  
# |__|___|___|
#
# fmt: on

####################
# Geometry

# cylinder c1
c1_r_top = r_crys
c1_z_top = h_melt
c1_r_bt = r_crys
c1_z_bt = 0

# ring r1
r1_r_top = r_melt
r1_z_top = h_melt
r1_r_bt = r_melt
r1_z_bt = 0
####################
# Mesh sizes
res_phi = 40

smallest_element_crucible = h_melt/40
layer_thickness_crucible = h_melt/80
growth_rate_melt = 1.4

smallest_element_side = smallest_element_crucible
layer_thickness_side = r_melt*0.05
growth_rate_side = 1.5
####################
res_z_c1, grading_bottom = boundary_layer(
    h_melt, "xmin", smallest_element_crucible, layer_thickness_crucible, growth_rate_melt
)

res_r_c1, grading_crys_rad = boundary_layer(
    r_melt, "xmax", smallest_element_side, layer_thickness_side, growth_rate_side
)

res_r_r1, grading_crucible = boundary_layer(
    h_melt, "xmin", smallest_element_crucible, layer_thickness_crucible, growth_rate_melt
)
####################
# Blocks (defined as cylinders & rings)
c1 = create_cylinder(
    mesh, [c1_r_top, c1_z_top], [c1_r_bt, c1_z_bt], res_r_c1, res_phi, res_z_c1, radius_ratio=0.25
)
r1 = create_ring(
    mesh,
    [r1_r_top, r1_z_top],
    [r1_r_bt, r1_z_bt],
    c1.surf_rad,
    res_r_r1,
    res_phi,
    res_z_c1,
)
####################
# Grading
c1.set_grading_axial(grading_bottom)
c1.set_grading_radial(grading_crys_rad)
r1.set_grading_axial(grading_bottom)
####################
# Patches
free_surf = Patch(mesh, "wall freeSurf")
free_surf.faces += r1.surf_top

crys_surf = Patch(mesh, "wall crysSurf")
crys_surf.faces += c1.surf_top

side_surf = Patch(mesh, "wall side")
side_surf.faces += r1.surf_rad

bottom_surf = Patch(mesh, "wall bottom")
bottom_surf.faces += c1.surf_bt
bottom_surf.faces += r1.surf_bt

os.remove(f"{sim_dir}/system/blockMeshDict")
mesh.write()
shutil.move("./system/blockMeshDict", f"{sim_dir}/system/blockMeshDict")

