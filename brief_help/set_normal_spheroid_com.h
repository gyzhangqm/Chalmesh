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
enum {THETA_MIN, THETA_MAX, PHI_MIN, PHI_MAX, X_SEMI, YZ_SEMI, CENTER, P_AXIS, 
      WIDTH, R1_LINES, R2_LINES, R3_LINES, STRETCHING, NO_STRETCHING, 
      SHOW, PLOT_MODE, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[18], *BRIEF_HELP[18];
int ARGUMENT[18], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[THETA_MIN] = "starting-r-angle";  
ARGUMENT[THETA_MIN] = 1;
BRIEF_HELP[THETA_MIN] = "Change the starting value of the theta angle. The theta angle "
"corresponds to the r-direction and measures the angle in the plane orthogonal to the "
"polar axis. In particular if the x-axis is polar axis, theta is the angle from the "
"y-axis. If the y-axis is polar axis, the angle is measured from the z-axis";

COMMAND[THETA_MAX] = "ending-r-angle";  
ARGUMENT[THETA_MAX] = 1;
BRIEF_HELP[THETA_MAX] = "Change the ending value of the theta angle. The theta angle "
"corresponds to the r-direction and measures the angle in the plane orthogonal to the "
"polar axis. In particular if the x-axis is polar axis, theta is the angle from the "
"y-axis. If the y-axis is polar axis, the angle is measured from the z-axis";

COMMAND[PHI_MIN]  = "starting-s-angle";  
ARGUMENT[PHI_MIN]  = 1;
BRIEF_HELP[PHI_MIN]  = "Change the starting value of the phi angle. The phi angle "
"corresponds to the s-direction and measures the angle from the polar-axis.";

COMMAND[PHI_MAX]  = "ending-s-angle";  
ARGUMENT[PHI_MAX]  = 1;
BRIEF_HELP[PHI_MAX]  = "Change the ending value of the phi angle. The phi angle "
"corresponds to the s-direction and measures the angle from the polar-axis.";

COMMAND[X_SEMI]    = "x-semi";  
ARGUMENT[X_SEMI]    = 1;
BRIEF_HELP[X_SEMI]    = "Change the semi-axis in the x-direction.";

COMMAND[YZ_SEMI]    = "yz-semi";  
ARGUMENT[YZ_SEMI]    = 1;
BRIEF_HELP[YZ_SEMI]    = "Change the semi-axis in the y and z-directions.";

COMMAND[CENTER]    = "center";  
ARGUMENT[CENTER]    = 0;
BRIEF_HELP[CENTER]    = "Move the center of the spheroid.";

COMMAND[P_AXIS]    = "polar-axis";  
ARGUMENT[P_AXIS]    = 1;
BRIEF_HELP[P_AXIS]    = "Change the polar axis. Give one integer in the range 1-2 "
"where 1 corresponds to the x-axis and 2 to the y-axis.";

COMMAND[WIDTH]    = "width";  
ARGUMENT[WIDTH]   = 1;
BRIEF_HELP[WIDTH] = "Set the width of the grid in the r3-direction.";

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
ARGUMENT[HELP]   = 1;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT]    = "exit";
ARGUMENT[EXIT]   = 0;
BRIEF_HELP[EXIT] = "Stop modifying the parameters and exit to the previous "
"command level.";


