import gmsh
import math
import sys

# Initialize gmsh
gmsh.initialize()
gmsh.model.add("bumpRect")

# --- Parameters ---
cl    = 0.075         # mesh size
A     = 0.1           # bump amplitude
x0    = 1.0           # bump center in x
ADVX  = 0.0           # ADVECTION_SPEED_X (set to 0 for static geometry)
t     = 0.0           # time parameter
sigma = 0.3           # bump width
numPoints = 20        # number of segments for the bump (including endpoints)

# --- Define Corner Points ---
# Bottom left
p1 = gmsh.model.geo.addPoint(0, 0, 0, cl)
# Bottom right
p2 = gmsh.model.geo.addPoint(5, 0, 0, cl)
# Top right
p3 = gmsh.model.geo.addPoint(5, 2, 0, cl)
# Top left (will be used as the end of the spline)
p4 = gmsh.model.geo.addPoint(0, 2, 0, cl)

# --- Define Edges for the Rectangle ---
# Bottom edge from p1 to p2
l1 = gmsh.model.geo.addLine(p1, p2)
# Right vertical edge from p2 to p3
l2 = gmsh.model.geo.addLine(p2, p3)

# --- Create the Bumpy Top Edge as a Spline ---
# We want the top edge to start at p3 (x=5) and end at p4 (x=0) with a bump.
# We will generate intermediate points between p3 and p4.
points_top = []
points_top.append(p3)  # starting point (at x=5, y=2)

for i in range(1, numPoints):
    # x value decreases linearly from 5 to 0
    x_val = 5 - i * (5.0 / numPoints)
    # Compute the bumped y-value using the provided function:
    y_val = 2 + A * math.exp(-((x_val - (x0 + ADVX * t))**2) / (2 * sigma**2))
    pt = gmsh.model.geo.addPoint(x_val, y_val, 0, cl)
    points_top.append(pt)

points_top.append(p4)  # ending point (at x=0, y=2)

# Create a spline that passes through the points on the top edge.
l3 = gmsh.model.geo.addSpline(points_top)

# --- Left vertical edge from p4 to p1 ---
l4 = gmsh.model.geo.addLine(p4, p1)

# --- Create a Closed Curve Loop and Plane Surface ---
cloop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
surface = gmsh.model.geo.addPlaneSurface([cloop])

# Optionally, you can set a transfinite meshing algorithm if desired:
# gmsh.model.geo.mesh.setTransfiniteSurface(surface)

# --- Synchronize the CAD model with gmsh ---
gmsh.model.geo.synchronize()

# --- Set Options for ASCII Mesh, Version 2 ---
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   # Use MSH version 2.2

# --- Define Physical Groups ---
# "Walls" for the bottom (l1) and left (l4) edges.
gmsh.model.addPhysicalGroup(1, [l1, l4], 1)
gmsh.model.setPhysicalName(1, 1, "Walls")
# "Outflow" for the right vertical edge (l2).
gmsh.model.addPhysicalGroup(1, [l2], 2)
gmsh.model.setPhysicalName(1, 2, "Outflow")
# "Wallm" for the bumped top edge (l3).
gmsh.model.addPhysicalGroup(1, [l3], 3)
gmsh.model.setPhysicalName(1, 3, "Wallm")
# "Domain" for the surface.
gmsh.model.addPhysicalGroup(2, [surface], 9)
gmsh.model.setPhysicalName(2, 9, "Domain")

# --- Generate Mesh and Save ---
gmsh.model.mesh.generate(2)
gmsh.write("bumpRect.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
