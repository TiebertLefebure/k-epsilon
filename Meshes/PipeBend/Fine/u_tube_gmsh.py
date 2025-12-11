import gmsh
import math
import sys

gmsh.initialize()
gmsh.model.add("U_bend_2D")

# -----------------------------
# Geometry parameters (meters)
# -----------------------------
mm = 1e-3

R_pipe = 14.0 * mm      # pipe radius
R_curv = 125.0 * mm     # centerline radius of bend

# Inner/outer radii of the bend walls
R_inner = R_curv - R_pipe
R_outer = R_curv + R_pipe

D = 2.0 * R_pipe        # pipe diameter
H_leg = 10.0 * D        # straight-leg length (change if you like)

lc = D / 20.0           # nominal mesh size

### WALL-UNIT SCALES ###

# Reynolds number: Re= 4.5 x 10^4
# Blasius for smooth pipe (turbulent): C_f = 0.3164 * Re^(-1/4) = 0.0218
# Friction velocity: u_τ = U * sqrt(C_f/2) = 0.148 m/s

# Wall units: y+ = u_τ * y / ν => y = y+ * ν / u_τ
# First cell center at y+ = 2-3
# Definitely y+ < 5

# y+ = 1 => y ≈ 0.006 mm
# y+ = 3 => y ≈ 0.018 mm 
# y+ = 5 => y ≈ 0.03 mm 

# First cell center at y = 0.02 mm => y+ ≈ 3.3
# First cell center at y = 0.03 mm => y+ ≈ 5


# For low-Re Lam–Bremhorst integrated to the wall, aiming for y+ < 3 is nice
# y+ < 5 is still acceptable if DOFs explode
# h_wall ≈ 0.02 mm is a good target



# Center of curvature at origin
cx, cy = 0.0, 0.0

# Convenience functions
def add_point(x, y):
    return gmsh.model.geo.addPoint(x, y, 0.0, lc)

# -----------------------------
# Points
# -----------------------------
# Bottom points on outer/inner walls (theta = 0 for right, pi for left)
p_out_R_bot = add_point(cx + R_outer, cy)      # (Ro, 0)
p_out_L_bot = add_point(cx - R_outer, cy)      # (-Ro, 0)
p_in_R_bot  = add_point(cx + R_inner, cy)      # (Ri, 0)
p_in_L_bot  = add_point(cx - R_inner, cy)      # (-Ri, 0)

# Mid points to force the *bottom* half of the circle
p_out_mid   = add_point(cx, cy - R_outer)
p_in_mid    = add_point(cx, cy - R_inner)

# Top of legs (shift up by H_leg)
p_out_R_top = add_point(cx + R_outer, cy + H_leg)
p_in_R_top  = add_point(cx + R_inner, cy + H_leg)

p_out_L_top = add_point(cx - R_outer, cy + H_leg)
p_in_L_top  = add_point(cx - R_inner, cy + H_leg)

# Center point for arcs
p_center = add_point(cx, cy)

# -----------------------------
# Lines and arcs
# -----------------------------
# Outer walls (legs)
L_out_left  = gmsh.model.geo.addLine(p_out_L_top, p_out_L_bot)
L_out_right = gmsh.model.geo.addLine(p_out_R_bot, p_out_R_top)

# Inner walls (legs)
L_in_right  = gmsh.model.geo.addLine(p_in_R_top, p_in_R_bot)
L_in_left   = gmsh.model.geo.addLine(p_in_L_bot, p_in_L_top)

# Inlet and outlet cuts (for 3D extrusion later if you want)
L_inlet  = gmsh.model.geo.addLine(p_in_L_top,  p_out_L_top)   # left cross-section
L_outlet = gmsh.model.geo.addLine(p_out_R_top, p_in_R_top)    # right cross-section

# Outer bend: left bottom -> mid bottom -> right bottom (bottom semicircle)
A_out_1 = gmsh.model.geo.addCircleArc(p_out_L_bot, p_center, p_out_mid)
A_out_2 = gmsh.model.geo.addCircleArc(p_out_mid,   p_center, p_out_R_bot)

# Inner bend: right bottom -> mid bottom -> left bottom (bottom semicircle)
A_in_1  = gmsh.model.geo.addCircleArc(p_in_R_bot,  p_center, p_in_mid)
A_in_2  = gmsh.model.geo.addCircleArc(p_in_mid,    p_center, p_in_L_bot)

# -----------------------------
# Curve loop (single fluid region)
# -----------------------------
loop = gmsh.model.geo.addCurveLoop([
    L_inlet,        # inlet (top left)
    L_out_left,     # outer left wall down
    A_out_1,        # outer bend left->mid
    A_out_2,        # outer bend mid->right
    L_out_right,    # outer right wall up
    L_outlet,       # outlet (top right)
    L_in_right,     # inner right wall down
    A_in_1,         # inner bend right->mid
    A_in_2,         # inner bend mid->left
    L_in_left       # inner left wall up
])

surf = gmsh.model.geo.addPlaneSurface([loop])

gmsh.model.geo.synchronize()

# -----------------------------
# Physical groups (for FEniCS, etc.)
# -----------------------------
gmsh.model.addPhysicalGroup(2, [surf], name="Fluid")

gmsh.model.addPhysicalGroup(1, [L_inlet],  name="Inlet")
gmsh.model.addPhysicalGroup(1, [L_outlet], name="Outlet")
gmsh.model.addPhysicalGroup(
    1,
    [L_out_left, L_out_right, L_in_left, L_in_right,
     A_out_1, A_out_2, A_in_1, A_in_2],
    name="Walls"
)

# -----------------------------
# Mesh size control for low-Re k-epsilon (integrated to the wall)
# -----------------------------

# Global mesh size bounds (in meters)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.005 * mm)   # 0.005 mm
gmsh.option.setNumber("Mesh.MeshSizeMax", 5.0   * mm)   # 5 mm (core)

# Curves that form the walls (same as "Walls" physical group)
wall_curves = [
    L_out_left, L_out_right, L_in_left, L_in_right,
    A_out_1, A_out_2, A_in_1, A_in_2
]

# Distance field to walls
dist = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(dist, "EdgesList", wall_curves)
gmsh.model.mesh.field.setNumber(dist, "Sampling", 150) # ()"Sampling", 400)

# Threshold field: distance -> element size
# h_wall ~ 0.02 mm => y+ ~ 3.3 for your case (U=1.42, nu=8.9e-7, D=0.028)
h_wall = 0.08 * mm  # 0.02 * mm    # fine resolution near the wall
h_bulk = 2.0  * mm      # coarser elements in the core

# Refined strip thickness
d_min = 0.5 * mm  # 1.0 mm     # up to 1 mm from wall ~ h_wall
d_max = 3.0 * mm  # 5.0 mm     # by 5 mm from wall ~ h_bulk

th = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(th, "InField",  dist)
gmsh.model.mesh.field.setNumber(th, "SizeMin",  h_wall)
gmsh.model.mesh.field.setNumber(th, "SizeMax",  h_bulk)
gmsh.model.mesh.field.setNumber(th, "DistMin",  d_min)
gmsh.model.mesh.field.setNumber(th, "DistMax",  d_max)

gmsh.model.mesh.field.setAsBackgroundMesh(th)

# -----------------------------
# Mesh + export
# -----------------------------
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # nice for FEniCS/meshio
gmsh.model.mesh.generate(2)
gmsh.write("u_bend_2d.msh")

gmsh.finalize()
