#file-type-1
reset-program
#
# This command file generates a grid around half an ellipsoid that is sitting on a
# cylinder. To do this, we use the hyperbolic grid generator together with
# the projection feature.
#
# We begin by making surface grids on the ellipsoid.
#
pause
make-surface-mapping upper
	ellipsoid
		polar-axis 2
		ending-r-angle 85
		starting-r-angle -85
		starting-s-angle 90
		ending-s-angle 165
		y-semi 2.0
		z-semi .5
		center
			0.0
			0.0
			0.0
		orientation
		s-points 30
		r-points 15
	exit
make-surface-mapping lower
	ellipsoid
		polar-axis 2
		starting-r-angle 95
		ending-r-angle 265
		starting-s-angle 90
		ending-s-angle 165
		y-semi 2
		z-semi .5
		orientation
		r-points 15
		s-points 30
	exit
make-surface-mapping edge
	ellipsoid
		polar-axis 3
		ending-r-angle 360
		starting-r-angle 180
		starting-s-angle 60
		ending-s-angle 120
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
		z-semi .5
		s-points 30
	exit
#
# Next we make a volume grid around the cylinder
#
pause
make-volume-mapping cylinder
	cylindrical-mapping
		max-radius 4
		min-radius 1
		starting-angle 90
		ending-angle 270
		x-thickness 4
		y-center 1
		x-center -2
		r1-points 20
		r2-points 25
		r3-points 30
	exit
#
# To be able to project grid faces onto the surface of the cylinder, we create 
# a surface grid from the face of the volume grid that is on the cylinder.
#
pause
make-surface-from-volume plane
	cylinder
		side 1
		direction 1
	exit
#
# We are now ready to grow volume grids out from the ellipsoid.
#
pause
make-volume-mapping upper
	hyperbolic-mapping upper
		thickness 0.5
		r3-points 25
#
# Lets introduce a stretching in the normal direction to resolve the 
# boundary layer that is likely to occur if the grid is used to simulate
# high reynolds number flow.
#
pause
		stretching
		 	hyperbolic-tangent-stretching
			start-grid-size 1e-4
			end-grid-size 0.17
		exit
		compute-mapping
#
# The grid face close to the cylinder needs to be projected onto the cylinder
# surface. Note that the projection also changes the ellipsoidal surface slightly
# to make sure that a grid line follows the intersection line between the two
# surfaces.
#
pause
		project-face plane
			project-low-r1
	exit
#
# We now repeat the procedure for the other two patches on the ellipsiod.
#
pause
make-volume-mapping edge
	hyperbolic-mapping edge
		r1-points 15
		r2-points 30
		r3-points 25
		thickness .6
		stretching
		 	hyperbolic-tangent-stretching
			start-grid-size 1e-4
			end-grid-size 0.17
		exit
		compute-mapping
		project-face plane
			project-low-r2
		project-face plane
			project-high-r2
	exit
make-volume-mapping lower
	hyperbolic-mapping lower
		thickness .5
		r3-points 25
		stretching
		 	hyperbolic-tangent-stretching
			start-grid-size 1e-4
			end-grid-size 0.17
		exit
		compute-mapping
		project-face plane
			project-low-r1
	exit
#
# All grid faces that are aligned with the physical boundary must be assigned
# a positive surface label. In this case, the domain is simply connected, so we 
# use the same number for all those grid faces.
#
pause
surface-label edge
	set-low-r3 1
	set-low-r2 1
	set-high-r2 1
exit
surface-label upper
	set-low-r3 1
	set-low-r1 1
exit
surface-label cylinder
	set-low-r1 1
	set-high-r1 1
	set-low-r2 1
	set-high-r2 1
	set-low-r3 1
	set-high-r3 1
exit
surface-label lower
	set-low-r3 1
	set-low-r1 1
exit
# 
# To remove the surface grid on the cylinder, inside of the ellipsiod, we
# assign an edge label to the grid edges that are on the intersection between 
# the ellipsoid and the cylinder.
#
pause
edge-label upper
	set-r1=0-r3=0 1
exit
edge-label upper
	set-r1=0-r3=0 1
exit
edge-label lower
	set-r1=0-r3=0 1
exit
edge-label edge
	set-r2=0-r3=0 1
	set-r2=1-r3=0 1
exit
#
# We speed up the algorithm by not attempting to cut any holes in the volume grids
# around the ellipsoid.
#
pause
cut-holes edge
	no
cut-holes lower
	no
cut-holes upper
	no
#
# We are now ready to make the overlapping grid.
#
pause
make-overlapping-grid over3d-1
#
# Lets set the parameters for implicit interpolation and interpolation width = 2
#
pause
	overlap-parameters
		implicit-interpolation
		interpolation-width 2
	exit
#
# We change the default priority of the component grids to move the cylinder grid
# last and the edge grid first.
#
pause
	change-member-list
		move-first edge
		move-last cylinder
	exit
#
# We do want the interpolation points to be corrected for boundary mismatch,
# but we are not interested in seeing the hybrid grid on the physical boundary
# or the result of the hole-cutting algorithm.
#
pause
	mismatch-correction yes	
	show-physical-boundaries no
	show-holes no
	show-parameters
#
# Lets make the overlapping grid. This will take a couple of minutes.
#
pause
	compute-overlap
#
# You are now looking at the physical and interpolating boundaries of
# the overlapping grid, intersected by a cutting plane through the center
# of the ellipsoid.
#
pause
	inspect-surface-grid surface-label-1
		parameter-arrows
		remove-one-patch cylinder-side=2-dir=1
#
# This is the surface grid on all physical boundaries, without
# the parameter arrows and the outermost grid surface to
# easier see the surface grid on the ellipsoid.
#
pause
	exit
exit
#
# End of the projection demo. Thanks for your attention.
#
pause



