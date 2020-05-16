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
enum {ORIENTATION, R_LINES, S_LINES, ALPHA_MIN, ALPHA_MAX, BETA_MIN, BETA_MAX, X_SEMI, 
      Y_SEMI, Z_SEMI, CENTER, P_AXIS, SHOW, PLOT_MODE, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[16], *BRIEF_HELP[16];
int ARGUMENT[16], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[ORIENTATION] = "orientation";
ARGUMENT[ORIENTATION] = 0;
BRIEF_HELP[ORIENTATION] = "Wrap the surface inside out, i.e., swap the properties of "
"the front and back sides of the surface and make the surface normal point in the "
"opposite direction.";

COMMAND[R_LINES]   = "r-points";       
ARGUMENT[R_LINES]   = 1;
BRIEF_HELP[R_LINES]   = "Change the number of grid lines in the r-direction.";

COMMAND[S_LINES]   = "s-points";       
ARGUMENT[S_LINES]   = 1;
BRIEF_HELP[S_LINES]   = "Change the number of grid lines in the s-direction.";

COMMAND[ALPHA_MIN] = "starting-r-angle";  
ARGUMENT[ALPHA_MIN] = 1;
BRIEF_HELP[ALPHA_MIN] = "Change the starting value of the theta angle. The theta angle "
"corresponds to the r-direction and measures the angle in the plane orthogonal to the "
"polar axis. In particular if the z-axis is polar axis, theta is the angle from the "
"x-axis. If the x-axis is polar axis, the angle is measured from the y-axis and if "
"the y-axis is polar axis, the angle is taken from the z-axis.";

COMMAND[ALPHA_MAX] = "ending-r-angle";  
ARGUMENT[ALPHA_MAX] = 1;
BRIEF_HELP[ALPHA_MAX] = "Change the ending value of the theta angle. The theta angle "
"corresponds to the r-direction and measures the angle in the plane orthogonal to the "
"polar axis. In particular if the z-axis is polar axis, theta is the angle from the "
"x-axis. If the x-axis is polar axis, the angle is measured from the y-axis and if "
"the y-axis is polar axis, the angle is taken from the z-axis.";

COMMAND[BETA_MIN]  = "starting-s-angle";  
ARGUMENT[BETA_MIN]  = 1;
BRIEF_HELP[BETA_MIN]  = "Change the starting value of the phi angle. The phi angle "
"corresponds to the s-direction and measures the angle from the polar axis.";

COMMAND[BETA_MAX]  = "ending-s-angle";  
ARGUMENT[BETA_MAX]  = 1;
BRIEF_HELP[BETA_MAX]  = "Change the ending value of the phi angle. The phi angle "
"corresponds to the s-direction and measures the angle from the polar axis.";

COMMAND[X_SEMI]    = "x-semi";  
ARGUMENT[X_SEMI]    = 1;
BRIEF_HELP[X_SEMI]    = "Change the length of the semi-axis in the x-direction.";

COMMAND[Y_SEMI]    = "y-semi";  
ARGUMENT[Y_SEMI]    = 1;
BRIEF_HELP[Y_SEMI]    = "Change the length of the semi-axis in the y-direction.";

COMMAND[Z_SEMI]    = "z-semi";  
ARGUMENT[Z_SEMI]    = 1;
BRIEF_HELP[Z_SEMI]    = "Change the length of the semi-axis in the z-direction.";

COMMAND[CENTER]    = "center";  
ARGUMENT[CENTER]    = 0;
BRIEF_HELP[CENTER]    = "Move the center of the ellipsoid.";

COMMAND[P_AXIS]    = "polar-axis";  
ARGUMENT[P_AXIS]    = 1;
BRIEF_HELP[P_AXIS]    = "Change the polar axis. Give one integer in the range 1-3 "
"where 1 corresponds to the x-axis, 2 to the y-axis, and 3 to the z-axis.";

COMMAND[SHOW]      = "show-parameters";
ARGUMENT[SHOW]      = 0;
BRIEF_HELP[SHOW]      = "Show the value of the parameters that control this surface "
"mapping.";

COMMAND[PLOT_MODE] = "surface-plot-mode";
/*file[PLOT_MODE] = "surface_plot_mode_com.h";*/
ARGUMENT[PLOT_MODE] = 0;
BRIEF_HELP[PLOT_MODE] = "Modify the graphical contents of the `Surface mappings' "
"window.";

COMMAND[HELP]      = "help";
ARGUMENT[HELP]      = 1;
BRIEF_HELP[HELP]      = NULL;

COMMAND[EXIT]      = "exit";
ARGUMENT[EXIT]      = 0;
BRIEF_HELP[EXIT]      = "Stop modifying the parameters and exit to the previous "
"command level.";


