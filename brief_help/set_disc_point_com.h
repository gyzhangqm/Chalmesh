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
enum {ORIENTATION, R_LINES, S_LINES, SAVE_SUBGRID, PLOT_MODE, SHOW, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[8], *BRIEF_HELP[8];
int ARGUMENT[8], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[ORIENTATION] = "orientation";
ARGUMENT[ORIENTATION] = 0;
BRIEF_HELP[ORIENTATION] = "Wrap the surface inside out, i.e., swap the properties of "
"the front and back sides of the surface and make the surface normal point in the "
"opposite direction.";

COMMAND[R_LINES]   = "r-lines";       
ARGUMENT[R_LINES]   = 1;
BRIEF_HELP[R_LINES]   = "Change the number of grid lines in the r-direction.";

COMMAND[S_LINES]   = "s-lines";       
ARGUMENT[S_LINES]   = 1;
BRIEF_HELP[S_LINES]   = "Change the number of grid lines in the s-direction.";

COMMAND[SAVE_SUBGRID]   = "save-subgrid";       
ARGUMENT[SAVE_SUBGRID]   = 1;
BRIEF_HELP[SAVE_SUBGRID]   = "Save a sub-grid of the present grid on a PLOT3D "
"ASCII file.";

COMMAND[PLOT_MODE] = "surface-plot-mode";
/*file[PLOT_MODE] = "surface_plot_mode_com.h";*/
ARGUMENT[PLOT_MODE] = 0;
BRIEF_HELP[PLOT_MODE] = "Modify the graphical contents of the `Surface mappings' "
"window.";

COMMAND[SHOW] = "show";
ARGUMENT[SHOW] = 0;
BRIEF_HELP[SHOW] = "Show some data for the mapping.";

COMMAND[HELP]      = "help";
ARGUMENT[HELP]      = 1;
BRIEF_HELP[HELP]      = NULL;

COMMAND[EXIT]      = "exit";
ARGUMENT[EXIT]      = 0;
BRIEF_HELP[EXIT]      = "Stop modifying the parameters and exit to the previous "
"command level.";


