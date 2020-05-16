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
enum {LEFT, RIGHT, LOWER, UPPER, MY_FAR, NEAR, SHOW, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[9], *BRIEF_HELP[9];
int ARGUMENT[9], *SAVE_ON_COPY=NULL;
const int LEVEL=1, NO_INDENT=0;

COMMAND[LEFT]  ="set-low-r1";  
ARGUMENT[LEFT]  = 1;
BRIEF_HELP[LEFT] = 
"Change the value of the boundary condition on the parameter face r1=0. "
"The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is helpful information "
"for a PDE solver that uses the grid.";

COMMAND[RIGHT] ="set-high-r1"; 
ARGUMENT[RIGHT] = 1;
BRIEF_HELP[RIGHT] = 
"Change the value of the boundary condition on the parameter face r1=1. "
"The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is helpful information "
"for a PDE solver that uses the grid.";

COMMAND[LOWER] ="set-low-r2"; 
ARGUMENT[LOWER] = 1;
BRIEF_HELP[LOWER] = 
"Change the value of the boundary condition on the parameter face r2=0. "
"The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is helpful information "
"for a PDE solver that uses the grid.";

COMMAND[UPPER] ="set-high-r2"; 
ARGUMENT[UPPER] = 1;
BRIEF_HELP[UPPER] = 
"Change the value of the boundary condition on the parameter face r2=1. "
"The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is helpful information "
"for a PDE solver that uses the grid.";

COMMAND[NEAR]  ="set-low-r3";  
ARGUMENT[NEAR]  = 1;
BRIEF_HELP[NEAR] = 
"Change the value of the boundary condition on the parameter face r3=0. "
"The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is helpful information "
"for a PDE solver that uses the grid.";

COMMAND[MY_FAR]="set-high-r3";
ARGUMENT[MY_FAR] = 1;
BRIEF_HELP[MY_FAR] = 
"Change the value of the boundary condition on the parameter face r3=1. "
"The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is helpful information "
"for a PDE solver that uses the grid.";

COMMAND[SHOW]  ="show";       
ARGUMENT[SHOW]  = 0;
BRIEF_HELP[SHOW] = "Print the boundary condition on all faces of the grid.";

COMMAND[HELP]  ="help";
ARGUMENT[HELP]  = 0;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT]  ="exit";
ARGUMENT[EXIT]  = 0;
BRIEF_HELP[EXIT] = "Stop changing the boundary condition and exit from this "
"command level.";
