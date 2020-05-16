#file-type-1
reset-program
#
# This command file illustrates how the hyperbolic grid generator
# works. We will present two examples: A grid around a ship propeller
# and a grid around the stern of a ship.
#
# We begin with the propeller. First, we read in surface grids on
# PLOT3D format.
#
pause
make-surface-mapping edge
	discrete-point
		plot3d-ascii edge-surface.3dsurf
	exit
make-surface-mapping upper
	discrete-point
		plot3d-ascii upper-surface.3dsurf
	exit
#
# We put a cylindrical grid around the hub (propeller axis) 
#
pause
make-volume-mapping hub
	cylindrical-mapping
		min-radius .2
		max-radius .45
		r1-points 30   #12
		r2-points 60   #50
		r3-points 32   #24
		starting-angle -90
		ending-angle    270
		x-thickness 1.2
		x-center    -.6
	exit
#
# To be able to project some grid surfaces onto the hub, we make a
# surface grid by extracting a grid face from the hub grid.
#
pause
make-surface-from-volume hub
        hub
                direction 1
        exit
#
# We proceed by making the volume grid outside the edge. For this
# purpose, we use the hyperbolic grid generator.
#
pause
make-volume-mapping edge
	hyperbolic-mapping edge
		velocity-threshold .25
		thickness .2
		r3-points 15
		curvature-coefficient 1.e-2
		compute-mapping
		volume-plot-mode
			specify-r2-surface .5
		exit
#
# The grid faces close to the hub (cylinder) needs to be projected onto the cylinder
# to follow that surface exactly. 
#
pause
		project-face hub
			project-low-r2
		project-face hub
			project-high-r2
	exit
#
# We repeat the procedure for the grid on the side of the propeller.
#
pause
make-volume-mapping upper
	hyperbolic-mapping upper
		velocity-threshold .25
		thickness .1
		r3-points 14
#
# It is possible to add a stretching in the normal direction to
# better resolve a boundary layer near the surface.
#
pause
		stretching
 			hyperbolic-tangent-stretching
			start-grid-size .005
 			end-grid-size .13        #0.10
			show-parameters
		exit
		curvature-coefficient 0
		compute-mapping
#
# This grid has only one face that need to be projected onto the
# cylinder.
#
pause
		project-face hub
			project-low-r2
	exit
#
# This concludes the propeller part of this demo.
#
pause
reset-program
#
# We will now generate a grid on a ship hull. This case is tricky
# since the surface contains both concave and convex parts.
#
# We start by reading in the surface from a PLOT3D file.
# 
pause
make-surface-mapping aft-bulb
	discrete-point
		plot3d-ascii slim-stern.3dsurf
	exit
#
# Lets make a volume grid with the hyperbolic grid generator.
#
pause
make-volume-mapping aft-bulb
	hyperbolic-mapping aft-bulb
		thickness 1.7
		curvature-coefficient .17
		averaging-coefficient .5
		velocity-threshold 0.5
		r3-points 12
		stretching
 			hyperbolic-tangent-stretching
			start-grid-size .005
 			end-grid-size .13
		exit
		boundary-condition
			low-r2 y-constant
		exit
		compute-mapping
	exit
#
# This is the end of the hyperbolic demo. Thanks for your attention!
#
pause

