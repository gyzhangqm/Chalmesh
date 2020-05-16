#file-type-1
reset-program
#
# This command file generates a grid around a quarter of an ellipsoid in a channel. 
# To make the body-fitted grids, we will use the normal-surface
# grid generator.
#
# We begin by making surface grids on the ellipsoid.
#
pause
make-surface-mapping main-ellipsoid
	ellipsoid
		x-semi 2
		y-semi 1
		z-semi 1
		polar-axis 1
		center
			1.000000e+00
			0.000000e+00
			1.000000e+00
		ending-r-angle 0
		starting-r-angle -90
		starting-s-angle 90
		ending-s-angle 155
		orientation
		s-points 30
		r-points 15
	exit
make-surface-mapping pole-patch
	ellipsoid
		x-semi 2
		y-semi 1
		polar-axis 2
		starting-r-angle -135
		ending-r-angle -90
		ending-s-angle 90
		starting-s-angle 45
		orientation
		s-points 14
		center
			1.000000e+00
			0.000000e+00
			1.000000e+00
	exit
#
# Next, we will make the volume grids. We begin by constructing the background grid 
# that fills the channel.
#
pause
make-volume-mapping background
	cartesian-mapping
		x-min -2
		x-max 1
		y-min 0
		y-max 2
		z-min -1
		z-max 1
		r1-points 25
		r2-points 35
		r3-points 25
	exit
#
# We proceed by making body-fitted volume grids by growing out normals
# from the surface grids on the ellipsoid.
#
pause
make-volume-mapping ellipsoid
	normal-surface-mapping main-ellipsoid
		width 0.7
		volume-plot-mode
			standard-view
		exit
		r1-points 31
		r2-points 35
		r3-points 10
	exit
make-volume-mapping northpole
	normal-surface-mapping pole-patch
		width 0.8
		volume-plot-mode
			standard-view
		exit
		r1-points 15
		r2-points 15
		r3-points 10
	exit
#
# All grid faces that are aligned with the physical boundary must be assigned
# a positive surface label. In this case, the domain is simply connected, so we 
# use the same number for all those grid faces.
#
pause
surface-label northpole
	set-high-r1 1
	set-high-r2 1
	set-low-r3 1
exit
surface-label ellipsoid
	set-low-r3 1
	set-low-r1 1
	set-low-r2 1
	set-high-r2 1
exit
surface-label background
	set-low-r1 1
	set-high-r1 1
	set-low-r2 1
	set-high-r2 1
	set-low-r3 1
	set-high-r3 1
exit
# 
# To remove the surface grid on the background grid, inside the ellipsiod, we
# assign an edge label to the grid edges that are on the intersection between 
# the ellipsoid and the background grid. Observe that there are three pieces
# corresponding to the intersection with the three planes of the background grid.
# Each piece has its own edge-label.
#
pause
edge-label ellipsoid
	set-r1=0-r3=0 1
	set-r2=0-r3=0 2
	set-r2=1-r3=0 3
exit
edge-label northpole
	set-r1=1-r3=0 2
	set-r2=1-r3=0 3
exit
#
# We are now ready to make the overlapping grid. This might take a minute...
#
pause
make-overlapping-grid test
	compute-overlap
exit
#
# End of the quarter-ellipsoid demo. Thanks for your attention.
#
pause
