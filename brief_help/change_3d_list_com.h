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
enum {INSERT_FIRST, MOVE_FIRST, MOVE_LAST, REMOVE, INSERT_ALL, SHOW, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[8], *BRIEF_HELP[8];
int ARGUMENT[8], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[INSERT_FIRST] = "insert-first";  
COMMAND[MOVE_FIRST] = "move-first";  
COMMAND[MOVE_LAST] = "move-last";  
COMMAND[REMOVE]    = "remove-component";  
COMMAND[INSERT_ALL] = "insert-all-components";  
COMMAND[SHOW]      = "show-list";  
COMMAND[HELP]      = "help";
COMMAND[EXIT]      = "exit";

ARGUMENT[INSERT_FIRST] = 1;
ARGUMENT[MOVE_FIRST]   = 1;
ARGUMENT[MOVE_LAST]    = 1;
ARGUMENT[REMOVE]       = 1;
ARGUMENT[INSERT_ALL] = 0;
ARGUMENT[SHOW]       = 0;
ARGUMENT[HELP]       = 0;
ARGUMENT[EXIT]       = 0;

BRIEF_HELP[INSERT_FIRST] = "Insert one component grid at the beginning of the "
"list of components to be used in the overlapping surface grid. Observe that a "
"component can only occur once in the list.";  
BRIEF_HELP[MOVE_FIRST] = "Move one component grid to the beginning of the "
"list of components to be used in the overlapping surface grid. Observe that a "
"component can only occur once in the list.";  
BRIEF_HELP[MOVE_LAST] = "Move one component grid to the end of the "
"list of components to be used in the overlapping surface grid. Observe that a "
"component can only occur once in the list.";  
BRIEF_HELP[REMOVE]    = "Remove one component grid from the list of components "
"which will be included in the overlapping surface grid. Observe that this command "
"only removes the link in the list and NOT the entire component grid structure.";  
BRIEF_HELP[INSERT_ALL] = "Insert all defined surface grids into the overlapping "
"grid list";
BRIEF_HELP[SHOW]    = "Show the list of component grids in the overlapping grid.";
BRIEF_HELP[HELP]      = NULL;
BRIEF_HELP[EXIT]      = "Stop modifying the parameters and exit to the previous "
"command level.";



