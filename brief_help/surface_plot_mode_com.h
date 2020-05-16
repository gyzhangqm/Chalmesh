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
enum {GRID_LINES, SURFACE, BOUNDARIES, ARROWS, COORD_CAGE, STD_VIEW, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[8], *BRIEF_HELP[8];
int *ARGUMENT=NULL, *SAVE_ON_COPY=NULL;
const int LEVEL=3;

COMMAND[GRID_LINES]  = "grid-lines";
BRIEF_HELP[GRID_LINES] = "Toggle the grid lines.";

COMMAND[SURFACE]     = "surface";
BRIEF_HELP[SURFACE] = "Toggle the opaque surface.";

COMMAND[BOUNDARIES]  = "highlight-boundaries";
BRIEF_HELP[BOUNDARIES] = "Toggle highlighting of the boundaries.";

COMMAND[ARROWS]      = "parameter-arrows";
BRIEF_HELP[ARROWS] = "Toggle arrows showing the parameter directions.";

COMMAND[COORD_CAGE]  = "coord-cage";
BRIEF_HELP[COORD_CAGE]  = "Toggle the coordinate axes.";

COMMAND[STD_VIEW]    = "standard-view";
BRIEF_HELP[STD_VIEW]  = "Reset the view point and focal point to show everything in "
"the bounding box.";

COMMAND[HELP]        = "help";
BRIEF_HELP[HELP]      = NULL;

BRIEF_HELP[EXIT]      = "Stop modifying the plotting parameters and exit to the "
"previous command level.";
COMMAND[EXIT]        = "exit";


