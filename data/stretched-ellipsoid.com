#file-type-1
reset-program
#
# This command file generates a grid around an ellipsoid in a channel. 
# To make the body-fitted grids, we will use the normal-surface
# grid generator.
#
# We begin by making surface grids on the ellipsoid.
#
pause
make-surface-mapping tail
	ellipsoid
		polar-axis 3
		starting-r-angle 135
		ending-r-angle 225
		starting-s-angle 45
		ending-s-angle 135
		r-points 14
		s-points 14
		x-semi 0.5
		y-semi 0.08333333333
		z-semi 0.08333333333
		center
			0.0
			0.0
			0.0
		orientation
	exit
make-surface-mapping nose
	ellipsoid
		polar-axis 3
		starting-r-angle -45
		ending-r-angle 45
		starting-s-angle 45
		ending-s-angle 135
		r-points 14
		s-points 14
		x-semi 0.5
		y-semi 0.08333333333
		z-semi 0.08333333333
		center
			0.0
			0.0
			0.0
		orientation
	exit
make-surface-mapping hull
	ellipsoid
		polar-axis 1
		starting-r-angle -180
		ending-r-angle 180
		starting-s-angle 25
		ending-s-angle 155
		r-points 30
		s-points 25
		x-semi 0.5
		y-semi 0.08333333333
		z-semi 0.08333333333
		center
			0.0
			0.0
			0.0
		orientation
	exit
#
# Next, we will make the volume grids. We begin by constructing the background grid 
# that fills the channel.
#
pause
make-volume-mapping background
	cartesian-mapping
		x-max 0.7
		x-min -0.7
		y-max 0.3
		y-min -0.3
		z-max 0.3
		z-min -0.3
		r1-points 45
		r3-points 25
		r2-points 25
	exit
#
# We proceed by making body-fitted volume grids by growing out normals
# from the surface grids on the ellipsoid.
#
pause
make-volume-mapping hull
	normal-surface-mapping hull
		width 0.12
		r1-points 31
		r2-points 25
		r3-points 17
#
# Lets introduce a stretching in the normal direction to resolve the 
# boundary layer that is likely to occur if the grid is used to simulate
# high reynolds number flow.
#
pause
		stretching
		 	hyperbolic-tangent-stretching
			start-grid-size 1e-2
			end-grid-size 0.17
		exit
	exit
make-volume-mapping nose
	normal-surface-mapping nose
		width 0.12
		r1-points 15
		r2-points 15
		r3-points 17
		stretching
		 	hyperbolic-tangent-stretching
			start-grid-size 1e-2
			end-grid-size 0.17
		exit
	exit
make-volume-mapping tail
	normal-surface-mapping tail
		width 0.12
		r1-points 15
		r2-points 15
		r3-points 17
		stretching
		 	hyperbolic-tangent-stretching
			start-grid-size 1e-2
			end-grid-size 0.17
		exit
	exit
#
# All grid faces that are aligned with the physical boundary must be assigned
# a positive surface label. In this case, the domain is simply connected, so we 
# use the same number for all those grid faces.
#
pause
surface-label nose
	set-low-r3 1
exit
surface-label tail
	set-low-r3 1
exit
surface-label hull
	set-low-r3 1
exit
surface-label background
	set-low-r1 2
	set-high-r1 2
	set-low-r2 2
	set-high-r2 2
	set-low-r3 2
	set-high-r3 2
exit
#
# To check the surface label, we plot all faces with surface label one.
#
pause
volume-plot-mode
	color-surface-label 1
exit
#
# To make the grid useful for a PDE solver, we assign boundary conditions to
# the faces that are aligned with the boundary of the computational domain.
# These values have no bearing on the overlapping grid.
#
pause
boundary-condition nose
	set-low-r3 1
exit
boundary-condition tail
	set-low-r3 1
exit
boundary-condition hull
	set-low-r3 1
exit
boundary-condition background
	set-low-r1 2
	set-high-r1 21
	set-low-r2 5
	set-high-r2 5
	set-low-r3 5
	set-high-r3 5
exit
#
# We are now ready to make the overlapping grid.
#
pause
make-overlapping-grid agard-ellipsoid
#
# Lets modify the default parameters to get implicit interpolation
#
pause
	overlap-parameters
		implicit-interpolation
	exit
#
# Lets check the value of the parameters governing the overlap algorithm.
#
pause
	show-parameters
#
# Compute the overlapping grid
#
pause
	compute-overlap
exit
#
# End of the stretched-ellipsoid demo. Thanks for your attention.
#
pause
