import shutil
import matplotlib.pyplot as plt
import numpy as np
from nemoblock import *
import yaml

########################################################################
# input 1

sim_dir = "."
dimension = 2

####################
# Mesh structure  (left: symmetry axis)
# fmt: off
#
# ^ y-direction
# |
# |
# |------->
#   x-direction
# __|__|__
# |  |  |  |
#  
# __|__|__
# |  |  |  |
# .c1|r1|r2|
# |__|__|__|
# |  |  |  |
# .c2|r3|r4|
# |__|__|__|
#
# fmt: on

####################


########################################################################
# Geometry
r_melt = 0.1
r_crys = 0.05
h_melt = 0.1

# cylinder midCyl1
midCyl_r_top = r_crys
midCyl_z_top = h_melt
midCyl_r_bt = r_crys
midCyl_z_bt = 0
# ring midRing1

midRing1_r_top = r_melt
midRing1_z_top = h_melt
midRing1_r_bt = r_melt
midRing1_z_bt = 0

########################################################################
# mesh resolution
res_r_mid = 50
res_r_r1 = 50
res_h_mid = 100
    
# wedge angle
phi1 = -2 
phi2 = 2
########################################################################

# nemoblock meshing
mesh = Mesh()

def evaluate_spline(r_start, r_end, phi, spline, res=100):
    return [cartesian(r, phi, spline(r)) for r in np.linspace(r_start, r_end, res, endpoint=False)]

def permuteCartesian(point):
    # [x, y, z] -> [x, z, y]
    return [point[0], point[2], point[1]]

midCyl = Block(
    mesh,
    permuteCartesian(cartesian(0,            phi1, midCyl_z_bt)),
    permuteCartesian(cartesian(midCyl_r_bt,  phi1, midCyl_z_bt)),
    permuteCartesian(cartesian(midCyl_r_top, phi1, midCyl_z_top)),
    permuteCartesian(cartesian(0,            phi1, midCyl_z_top)),
    permuteCartesian(cartesian(0,            phi2, midCyl_z_bt)),
    permuteCartesian(cartesian(midCyl_r_bt,  phi2, midCyl_z_bt)),
    permuteCartesian(cartesian(midCyl_r_top, phi2, midCyl_z_top)),
    permuteCartesian(cartesian(0,            phi2, midCyl_z_top))
)
midCyl.p3 = "p7" # wedge
midCyl.p0 = "p4"
midCyl.set_number_of_cells(res_r_mid, res_h_mid, 1)
midCyl.create()

# midRing1 (right of midCyl)
midRing1 = Block(mesh)
midRing1.set_connection(midCyl, "left") # needs points 1,2,5,6
midRing1.p1 = permuteCartesian(cartesian(midRing1_r_bt, phi1, midRing1_z_bt))
midRing1.p2 = permuteCartesian(cartesian(midRing1_r_top, phi1, midRing1_z_top))
midRing1.p5 = permuteCartesian(cartesian(midRing1_r_bt, phi2, midRing1_z_bt))
midRing1.p6 = permuteCartesian(cartesian(midRing1_r_top, phi2, midRing1_z_top))
midRing1.cells_x1 = res_r_r1
midRing1.create()

# ####################
# # Grading
midCyl.grading = f"simpleGrading (4 1 0.25)"
midRing1.grading = f"simpleGrading (0.25 1 0.25)"

# ####################
# Patches
melt_surf = Patch(mesh, "wall freeSurf")
melt_surf.add_face(midRing1.face_back)

crys_surf = Patch(mesh, "wall crysSurf")
crys_surf.add_face(midCyl.face_back)

side_surf = Patch(mesh, "wall side")
side_surf.add_face(midRing1.face_right)

bottom_surf = Patch(mesh, "wall bottom")
bottom_surf.add_face(midCyl.face_front)
bottom_surf.add_face(midRing1.face_front)

mesh_front = Patch(mesh, "wedge front")
mesh_front.add_face(midRing1.face_top)
mesh_front.add_face(midCyl.face_top)

mesh_back = Patch(mesh, "wedge back")
mesh_back.add_face(midRing1.face_bottom)
mesh_back.add_face(midCyl.face_bottom)

os.remove(f"{sim_dir}/system/blockMeshDict")
mesh.write()
shutil.move("./system/blockMeshDict", f"{sim_dir}/system/blockMeshDict")
