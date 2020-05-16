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
enum {THICKNESS, R_MIN, R_MAX, THETA_MIN, THETA_MAX, X_CENTER, Y_CENTER, Z_CENTER, 
      R1_LINES, R2_LINES, R3_LINES, SHOW, PLOT_MODE, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[16], *BRIEF_HELP[16];
int ARGUMENT[16], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[THICKNESS]  = "x-thickness";
ARGUMENT[THICKNESS]  = 1;
BRIEF_HELP[THICKNESS]  = "Change the thickness of the cylinder.";

COMMAND[R_MIN]  = "min-radius";
ARGUMENT[R_MIN]  = 1;
BRIEF_HELP[R_MIN]  = "Change the minimum radius of the cylinder.";

COMMAND[R_MAX]  = "max-radius";
ARGUMENT[R_MAX]  = 1;
BRIEF_HELP[R_MAX]  = "Change the maximum radius of the cylinder.";

COMMAND[THETA_MIN]  = "starting-angle";
ARGUMENT[THETA_MIN]  = 1;
BRIEF_HELP[THETA_MIN]  = "Change the starting angle.";

COMMAND[THETA_MAX]  = "ending-angle";
ARGUMENT[THETA_MAX]  = 1;
BRIEF_HELP[THETA_MAX]  = "Change the ending angle.";

COMMAND[X_CENTER]  = "x-center";
ARGUMENT[X_CENTER]  = 1;
BRIEF_HELP[X_CENTER]  = "Change the x-coordinate of the origin.";

COMMAND[Y_CENTER]  = "y-center";
ARGUMENT[Y_CENTER]  = 1;
BRIEF_HELP[Y_CENTER]  = "Change the y-coordinate of the origin.";

COMMAND[Z_CENTER]  = "z-center";
ARGUMENT[Z_CENTER]  = 1;
BRIEF_HELP[Z_CENTER]  = "Change the z-coordinate of the origin.";

COMMAND[R1_LINES]   = "r1-points";       
ARGUMENT[R1_LINES]   = 1;
BRIEF_HELP[R1_LINES]   = "Change the number of grid lines in the r1-direction.";

COMMAND[R2_LINES]   = "r2-points";       
ARGUMENT[R2_LINES]   = 1;
BRIEF_HELP[R2_LINES]   = "Change the number of grid lines in the r2-direction.";

COMMAND[R3_LINES]   = "r3-points";       
ARGUMENT[R3_LINES]   = 1;
BRIEF_HELP[R3_LINES]   = "Change the number of grid lines in the r3-direction.";

COMMAND[SHOW] = "show-parameters";
ARGUMENT[SHOW] = 0;
BRIEF_HELP[SHOW] = "Show the parameters that control the mapping.";

COMMAND[PLOT_MODE] = "volume-plot-mode";
/*file[PLOT_MODE] = "volume_plot_mode_com.h";*/
ARGUMENT[PLOT_MODE] = 0;
BRIEF_HELP[PLOT_MODE] = "Modify the graphical contents of the "
"`Volume mappings' window.";

COMMAND[HELP]      = "help";
ARGUMENT[HELP]      = 1;
BRIEF_HELP[HELP]      = NULL;

COMMAND[EXIT]      = "exit";
ARGUMENT[EXIT]      = 0;
BRIEF_HELP[EXIT]      = "Stop modifying the parameters and exit to the previous "
"command level.";



