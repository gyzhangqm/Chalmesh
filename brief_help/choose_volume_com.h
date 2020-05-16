/*     Chalmesh, 3-D overlapping grid generator */
/*     Copyright (C) 1997 Anders Petersson */

/*     This program is free software; you can redistribute it and/or modify */
/*     it under the terms of the GNU General Public License as published by */
/*     the Free Software Foundation; either version 2 of the License, or */
/*     (at your option) any later version. */

/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/*     You can contact the author of Chalmesh by email at andersp@na.chalmers.se,  */
/*     or by paper-mail to  */

/*     Anders Petersson */
/*     Hydromechanics division */
/*     Department of Naval Architecture and Ocean Engineering */
/*     Chalmers University of Technology */
/*     412 96 Gothenburg */
/*     SWEDEN */
enum {CARTESIAN, CYLINDRICAL, NORMAL_SPHEROID, NORMAL_SURFACE, 
      DISCRETE_POINT, HYPERBOLIC, HELP, CANCEL};
const int LAST_COM=CANCEL, LEVEL=1, NO_INDENT=0;

char *COMMAND[10], *BRIEF_HELP[10];
int ARGUMENT[10], *SAVE_ON_COPY=NULL; 

COMMAND[CARTESIAN] = "cartesian-mapping";  
/*file[CARTESIAN] = "set_cartesian_com.h";*/
ARGUMENT[CARTESIAN] = 0;
BRIEF_HELP[CARTESIAN] = "Define the mapping to be Cartesian. This mapping has "
"an analytical inverse.";

COMMAND[CYLINDRICAL] = "cylindrical-mapping";  
/*file[CYLINDRICAL] = "set_cylindrical_com.h";*/
ARGUMENT[CYLINDRICAL] = 0;
BRIEF_HELP[CYLINDRICAL] = "Define the mapping to be cylindrical. This mapping has "
"an analytical inverse.";

COMMAND[NORMAL_SPHEROID] = "normal-spheroid-mapping";  
/*file[NORMAL_SPHEROID] = "set_normal_spheroid_com.h";*/
ARGUMENT[NORMAL_SPHEROID] = 0;
BRIEF_HELP[NORMAL_SPHEROID] = "Let the mapping be defined by growing normals out "
"from a spheroid. This mapping has "
"an analytical inverse.";

COMMAND[NORMAL_SURFACE] = "normal-surface-mapping";  
/*file[NORMAL_SURFACE] = "set_normal_surface_com.h";*/
ARGUMENT[NORMAL_SURFACE] = 1;
BRIEF_HELP[NORMAL_SURFACE] = "Let the mapping be defined by growing normals out "
"from a surface. This mapping does NOT have a known analytical inverse.";

COMMAND[DISCRETE_POINT] = "discrete-point-mapping";  
/*file[DISCRETE_POINT] = "init_volume_point_com.h";*/
ARGUMENT[DISCRETE_POINT] = 0;
BRIEF_HELP[DISCRETE_POINT] = "Define the mapping to interpolate tri-linearly "
"between a discrete structured set of grid points, which are read from a file that "
"for instance could be the output from an external grid generator. This mapping does "
"NOT have a known analytical inverse.";

COMMAND[HYPERBOLIC] = "hyperbolic-mapping";  
/*file[HYPERBOLIC] = "set_hyp_grid_com.h";*/
ARGUMENT[HYPERBOLIC] = 1;
BRIEF_HELP[HYPERBOLIC] = "Define the mapping by solving a Hamilton-Jacobi equation "
"in the direction normal to a surface. (The mapping is called hyperbolic since it "
"belongs to the class of grid generation techniques that commonly is known as "
"hyperbolic). This mapping does NOT have a known analytical inverse. The mapping "
"function is\n"
"\n"
"d/dt (x, y, z) = ds/dt * v(kappa) * unit_normal\n"
"\n"
"where t and s are the parameter and arclength, respectively, along a grid line "
"normal to the surface. kappa is the mean curvature of the surface and\n"
"\n"
"v(kappa) = max( v_min, 1 - eps * kappa )\n"
"\n"
"is the speed function. Here v_min is the velocity threshold that dictates the "
"propagation speed close to highly curved convex surfaces. "
"The quantity v(kappa) * unit_normal, can be averaged to make the grid smoother. "
"In this case the algebraic average of the nearest eight neighbors is computed and "
"weighted together with the point value.";


COMMAND[HELP]   = "help";  
ARGUMENT[HELP]   = 0;
BRIEF_HELP[HELP]   = NULL;

COMMAND[CANCEL] = "cancel";       
ARGUMENT[CANCEL] = 0;
BRIEF_HELP[CANCEL] = "Go directly back to the previous command level, without "
"defining a new volume mapping.";


