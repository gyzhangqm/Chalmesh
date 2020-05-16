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
enum {COMPUTE, CURV_FACTOR, SMOOTH_FACTOR, V_MIN, WIDTH, BOUNDARY_CONDITION,
      R1_LINES, R2_LINES, R3_LINES, PROJECT_FACE, STRETCHING, NO_STRETCHING, 
      SHOW, PLOT_MODE, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[16], *BRIEF_HELP[16];
int ARGUMENT[16], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[PROJECT_FACE]  = "project-face";
/*file[PROJECT_FACE]  = "compute_projection_com.h"; */
ARGUMENT[PROJECT_FACE]  = 1;
BRIEF_HELP[PROJECT_FACE] = "Project one face of the mapping onto a surface.";

COMMAND[BOUNDARY_CONDITION]  = "boundary-condition";
/*file[BOUNDARY_CONDITION]  = "set_hyp_bc_com.h"; */
ARGUMENT[BOUNDARY_CONDITION]  = 0;
BRIEF_HELP[BOUNDARY_CONDITION] = "Set the boundary condition for the hyperbolic "
"grid generation.";

COMMAND[COMPUTE]    = "compute-mapping";  
ARGUMENT[COMPUTE]   = 0;
BRIEF_HELP[COMPUTE] = "Solve the hyperbolic PDE with the current parameter values.";

COMMAND[CURV_FACTOR]    = "curvature-coefficient";  
ARGUMENT[CURV_FACTOR]   = 1;
BRIEF_HELP[CURV_FACTOR] = "Set the coefficient e in (1 - e*k), which is "
"the propagation speed of the grid surface. In this expression, k=(k1+k2)/2 "
"is the mean curvature of the grid surface at the grid point. Also see the command "
"`velocity-threshold'.";

COMMAND[SMOOTH_FACTOR]    = "averaging-coefficient";  
ARGUMENT[SMOOTH_FACTOR]   = 1;
BRIEF_HELP[SMOOTH_FACTOR] = "Set the coefficient for the spacial "
"averaging. The averaging is done on the velocity * normal. A factor of 0.0 "
"means no averaging, so the velocity * normal is determined pointwise. A value "
"of 1.0 implies that the point value is replaced by the algebraic average of the "
"nearest neighbors. A value inbetween 0 and 1 implies a linear interpolation "
"between the two extrema.";

COMMAND[V_MIN]    = "velocity-threshold";  
ARGUMENT[V_MIN]   = 1;
BRIEF_HELP[V_MIN] = "Set the minimum propagation speed of the grid surface. This "
"value is used near convex corners, where the propagation speed otherwise could "
"go negative.";

COMMAND[WIDTH]    = "thickness";  
ARGUMENT[WIDTH]   = 1;
BRIEF_HELP[WIDTH] = "Set the approximate thickness of the grid in the r3-direction. "
"The actual thickness will vary and be smaller where the surface is convex and be "
"larger where it is concave.";

COMMAND[R1_LINES]    = "r1-points";
ARGUMENT[R1_LINES]   = 1;
BRIEF_HELP[R1_LINES] = "Change the number of grid lines in the r1-direction.";

COMMAND[R2_LINES]    = "r2-points";
ARGUMENT[R2_LINES]   = 1;
BRIEF_HELP[R2_LINES] = "Change the number of grid lines in the r2-direction.";

COMMAND[R3_LINES]    = "r3-points";
ARGUMENT[R3_LINES]   = 1;
BRIEF_HELP[R3_LINES] = "Change the number of grid lines in the r3-direction.";

COMMAND[STRETCHING]    = "stretching";
/*file[STRETCHING] = "choose_stretching_com.h"; */
ARGUMENT[STRETCHING]   = 1;
BRIEF_HELP[STRETCHING] = "Introduce a stretching function in the normal direction, "
"or, if a stretching function already is in use, modify its parameters.";

COMMAND[NO_STRETCHING]    = "no-stretching";
ARGUMENT[NO_STRETCHING]   = 0;
BRIEF_HELP[NO_STRETCHING] = "Disable the stretching function in the normal direction "
"and use a uniform distribution of grid points instead.";

COMMAND[SHOW]    = "show-parameters";
ARGUMENT[SHOW]   = 0;
BRIEF_HELP[SHOW] = "Display the parameters that control the mapping.";

COMMAND[PLOT_MODE]    = "volume-plot-mode";
/*file[PLOT_MODE] = "volume_plot_mode_com.h";*/
ARGUMENT[PLOT_MODE]   = 0;
BRIEF_HELP[PLOT_MODE] = "Modify the graphical contents of the "
"`Volume mappings' window.";

COMMAND[HELP]    = "help";
ARGUMENT[HELP]   = 0;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT]    = "exit";
ARGUMENT[EXIT]   = 0;
BRIEF_HELP[EXIT] = "Stop modifying the parameters and exit to the previous "
"command level.";


