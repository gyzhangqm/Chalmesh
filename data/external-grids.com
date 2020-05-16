#file-type-1
reset-program
#
# This command file generates a grid around a sphere in a box. 
# Some volume grids are read in from a PLOT3D file, while other grids
# are generated by the normal-surface grid generator starting from a
# surface grid that is read from a PLOT3D file.
#
# We start by reading in external volume grids that cover the two poles on the sphere. 
#
pause
make-volume-mapping south-pole
	discrete-point-mapping
		plot3d-ascii south.3dvol
	exit
make-volume-mapping north-pole
	discrete-point-mapping
		plot3d-ascii north.3dvol
	exit
#
# For the main part of the sphere, we read in an external surface grid
#
pause
make-surface-mapping main-sphere
	discrete-point
		plot3d-ascii main-sphere.3dsurf
	exit
#
# We make a volume grid in the vicinity of the surface grid by growing 
# normals out from the surface
#
pause
make-volume-mapping main-globe
	normal-surface-mapping main-sphere
		width .45
		r3-points 7
	exit
#
# Lets make a Cartesian background grid 
#
pause
make-volume-mapping box
	cartesian-mapping
		x-min -1.5
		x-max 1.5
		y-min -1.5
		y-max 1.5
		z-min -1.5
		z-max 1.5
		r1-points 34
		r2-points 34
		r3-points 34
	exit
#
# To make the grid useful for a PDE solver, we assign boundary conditions to
# the faces that are aligned with the boundary of the computational domain.
# These values have no bearing on the overlapping grid.
#
pause
boundary-condition main-globe
	set-low-r3 1
exit
boundary-condition north-pole
	set-low-r3 1
exit
boundary-condition south-pole
	set-low-r3 1
exit
#
# All grid faces that are aligned with the physical boundary must be assigned
# a positive surface label. In this case, the domain is simply connected, so we 
# use the same number for all those grid faces.
#
pause
surface-label main-globe
	set-low-r3 1
exit
surface-label north-pole
	set-low-r3 1
exit
surface-label south-pole
	set-low-r3 1
exit
surface-label box
	set-low-r1 2
	set-high-r1 2
	set-low-r2 2
	set-high-r2 2
	set-low-r3 2
	set-high-r3 2
exit
#
# We are now ready to make the overlapping grid.
#
pause
make-overlapping-grid external
#
# We change the default priority of the component grids to move the box grid last
#
pause
	change-member-list
		move-last box
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
# Compute the overlapping grid. This might take a minute...
#
pause
	compute-overlap
#
# You are looking at the hole-cutting surface
#
pause
	proceed
#
# This is the grid after the holes have been cut.
#
pause
	proceed
exit
#
# This concludes the external-grids demo. Thanks for your attention!
#
pause