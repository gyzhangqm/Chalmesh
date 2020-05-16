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
enum {SIDE, DIRECTION, PLOT_MODE, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[5], *BRIEF_HELP[5];
int ARGUMENT[5], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[SIDE]      = "side";  
COMMAND[DIRECTION] = "direction";  

COMMAND[PLOT_MODE] = "surface-plot-mode";
COMMAND[HELP]      = "help";
COMMAND[EXIT]      = "exit";

ARGUMENT[SIDE]    = 1;
ARGUMENT[DIRECTION] = 1;

ARGUMENT[PLOT_MODE] = 0;
ARGUMENT[HELP]      = 1;
ARGUMENT[EXIT]      = 0;

BRIEF_HELP[SIDE] = "Choose the side of the volume mapping that should be used "
"in the definition of the surface mapping.";
BRIEF_HELP[DIRECTION] = "Choose the direction of the volume mapping that "
"should be used in the definition of the surface mapping.";

BRIEF_HELP[PLOT_MODE] = "Look at different features of the grid.";
BRIEF_HELP[HELP]      = NULL;
BRIEF_HELP[EXIT]      = "Stop modifying the parameters and exit to the previous "
"command level.";


