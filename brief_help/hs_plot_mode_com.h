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
enum {BOUNDARIES, ACTIVE_POINTS, ALL_LINES, CURVE_LABEL, EDGE_LABEL, 
	ARROWS, FLAG, COORD_CAGE, STD_VIEW, 
	SHOW_ONE_COMPONENT, ADD_ONE_COMPONENT, REMOVE_ONE_COMPONENT, 
	SHOW_ALL_COMPONENTS, REMOVE_ALL_COMPONENTS, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[16], *BRIEF_HELP[16];
int ARGUMENT[16], *SAVE_ON_COPY=NULL;
const int LEVEL=2;

COMMAND[BOUNDARIES]  = "highlight-boundaries";
ARGUMENT[BOUNDARIES] = 0;
BRIEF_HELP[BOUNDARIES] = "Toggle the boundaries of the component grids.";

COMMAND[ACTIVE_POINTS] = "active-points";
ARGUMENT[ACTIVE_POINTS] = 0;
BRIEF_HELP[ACTIVE_POINTS] = "Toggle the active (non-dead) points in the "
"overlapping surface grid";

COMMAND[ALL_LINES] = "all-grid-lines";
ARGUMENT[ALL_LINES] = 0;
BRIEF_HELP[ALL_LINES] = "Toggle the plotting of all grid lines.";

COMMAND[CURVE_LABEL] = "curve-label";
ARGUMENT[CURVE_LABEL] = 0;
BRIEF_HELP[CURVE_LABEL] = "Toggle the plotting of the curve label.";

COMMAND[EDGE_LABEL] = "edge-label";
ARGUMENT[EDGE_LABEL] = 0;
BRIEF_HELP[EDGE_LABEL] = "Toggle the plotting of the edge labels.";

COMMAND[FLAG] = "color-flag";
ARGUMENT[FLAG] = 0;
BRIEF_HELP[FLAG] = "Toggle the plotting of the color coded flag array.";

COMMAND[ARROWS]      = "parameter-arrows";
ARGUMENT[ARROWS] = 0;
BRIEF_HELP[ARROWS] = "Toggle arrows showing the parameter directions.";

COMMAND[COORD_CAGE]  = "coord-cage";
ARGUMENT[COORD_CAGE] = 0;
BRIEF_HELP[COORD_CAGE]  = "Toggle the coordinate axes.";

COMMAND[STD_VIEW]    = "standard-view";
ARGUMENT[STD_VIEW] = 0;
BRIEF_HELP[STD_VIEW]  = "Reset the view point and focal point to show everything in "
"the bounding box.";

COMMAND[SHOW_ONE_COMPONENT] = "show-one-patch";
ARGUMENT[SHOW_ONE_COMPONENT] = 1;
BRIEF_HELP[SHOW_ONE_COMPONENT] = "Specify which single patch to plot.";

COMMAND[ADD_ONE_COMPONENT] = "add-one-patch";
ARGUMENT[ADD_ONE_COMPONENT] = 1;
BRIEF_HELP[ADD_ONE_COMPONENT] = "Add one patch to the plot.";

COMMAND[REMOVE_ONE_COMPONENT] = "remove-one-patch";
ARGUMENT[REMOVE_ONE_COMPONENT] = 1;
BRIEF_HELP[REMOVE_ONE_COMPONENT] = "Remove on patch from the plot.";

COMMAND[SHOW_ALL_COMPONENTS] = "show-all-patches";
ARGUMENT[SHOW_ALL_COMPONENTS] = 0;
BRIEF_HELP[SHOW_ALL_COMPONENTS] = "Show all patches.";

COMMAND[REMOVE_ALL_COMPONENTS] = "remove-all-patches";
ARGUMENT[REMOVE_ALL_COMPONENTS] = 0;
BRIEF_HELP[REMOVE_ALL_COMPONENTS] = "Don't show any patches.";

COMMAND[HELP]        = "help";
ARGUMENT[HELP] = 0;
BRIEF_HELP[HELP]      = NULL;

COMMAND[EXIT]        = "exit";
ARGUMENT[EXIT] = 0;
BRIEF_HELP[EXIT]      = "Stop modifying the plotting parameters and exit to the "
"previous command level.";


