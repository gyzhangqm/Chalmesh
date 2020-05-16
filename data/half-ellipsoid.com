#file-type-1
reset-program
#
# This command file generates a grid around half an ellipsoid in a channel. 
# To make the body-fitted grids, we will use the normal-surface
# grid generator.
#
# We begin by making surface grids on the ellipsoid.
#
pause
make-surface-mapping main-ellipsoid
	ellipsoid
		polar-axis 2
		ending-r-angle -5
		starting-r-angle -175
		starting-s-angle 90
		ending-s-angle 165
		y-semi 2.0
		center
			0.0
			0.0
			0.0
		orientation
		s-points 30
		r-points 15
		x-semi .5
	exit
make-surface-mapping pole-patch
	ellipsoid
		polar-axis 1
		starting-r-angle 90
		ending-r-angle 270
		starting-s-angle 70
		ending-s-angle 110
		center
			0.0
			0.0
			0.0
		orientation
		surface-plot-mode
			standard-view
		exit
		s-points 14
		y-semi 2.0
		x-semi .5
		s-points 30
	exit
make-surface-mapping ellipse-2
	ellipsoid
		polar-axis 2
		starting-r-angle 5
		ending-r-angle 175
		starting-s-angle 90
		ending-s-angle 165
		y-semi 2
		orientation
		r-points 15
		s-points 15
		x-semi .5
	exit
#
# Next, we will make the volume grids. We begin by constructing the background grid 
# that fills the channel.
#
pause
make-volume-mapping background
	cartesian-mapping
		x-max 2
		x-min -2
		y-max 0
		y-min -3
		z-max 2
		z-min -2
		volume-plot-mode
			standard-view
		exit
		r1-points 25
		r3-points 25
		r2-points 35
	exit
#
# We proceed by making body-fitted volume grids by growing out normals
# from the surface grids on the ellipsoid.
#
pause
make-volume-mapping ellipsoid
	normal-surface-mapping main-ellipsoid
		width 0.65
		r1-points 31
		r2-points 35
		r3-points 10
	exit
make-volume-mapping northpole
	normal-surface-mapping pole-patch
		r1-points 15
		r3-points 10
		width .6
		r2-points 30
	exit
make-volume-mapping southpole
	normal-surface-mapping ellipse-2
		width .5
	exit
#
# All grid faces that are aligned with the physical boundary must be assigned
# a positive surface label. In this case, the domain is simply connected, so we 
# use the same number for all those grid faces.
#
pause
surface-label northpole
	set-low-r3 1
	set-low-r2 1
	set-high-r2 1
exit
surface-label ellipsoid
	set-low-r3 1
	set-low-r1 1
exit
surface-label background
	set-low-r1 1
	set-high-r1 1
	set-low-r2 1
	set-high-r2 1
	set-low-r3 1
	set-high-r3 1
exit
surface-label southpole
	set-low-r3 1
	set-low-r1 1
exit
# 
# To remove the surface grid on the background grid, inside the ellipsiod, we
# assign an edge label to the grid edges that are on the intersection between 
# the ellipsoid and the background grid.
#
pause
edge-label ellipsoid
	set-r1=0-r3=0 1
exit
edge-label ellipsoid
	set-r1=0-r3=0 1
exit
edge-label southpole
	set-r1=0-r3=0 1
exit
edge-label northpole
	set-r2=0-r3=0 1
	set-r2=1-r3=0 1
exit
#
# We are now ready to make the overlapping grid.
#
pause
make-overlapping-grid over3d-1
#
# Lets modify the default parameters to get implicit interpolation
#
pause
	overlap-parameters
		implicit-interpolation
		interpolation-width 2
	exit
#
# We change the default priority of the component grids to move the northpole grid first
#
pause
	change-member-list
		move-first northpole
	exit
#
# We can select to view the hole-cutting surface as well as the 
# result of the hole-cutting algorithm.
#
pause
	show-physical-boundaries yes
	show-holes yes
#
# Lets check the value of the parameters governing the overlap algorithm.
#
pause
	show-parameters
#
# Make the overlapping grid. This might take a couple of minutes.
#
pause
	compute-overlap
pause
	proceed
pause
	proceed
exit
#
# End of the half-ellipsoid demo. Thanks for your attention.
#
pause
