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
enum {EDGE_1, EDGE_2, EDGE_3, EDGE_4, EDGE_5, EDGE_6, EDGE_7, EDGE_8, EDGE_9, 
      EDGE_10, EDGE_11, EDGE_12, SHOW, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[15], *BRIEF_HELP[15];
int ARGUMENT[15], *SAVE_ON_COPY=NULL;
const int LEVEL=1, NO_INDENT=0;

COMMAND[EDGE_1]  ="set-r1=0-r2=0";  
ARGUMENT[EDGE_1]  = 1;
BRIEF_HELP[EDGE_1] = "Change the edge label on the edge where the parameter face "
"r1=0 intersects the parameter face r2=0.";

COMMAND[EDGE_2]  ="set-r1=1-r2=0";  
ARGUMENT[EDGE_2]  = 1;
BRIEF_HELP[EDGE_2] = "Change the edge label on the edge where the parameter face "
"r1=1 intersects the parameter face r2=0.";

COMMAND[EDGE_3]  ="set-r1=0-r2=1";  
ARGUMENT[EDGE_3]  = 1;
BRIEF_HELP[EDGE_3] = "Change the edge label on the edge where the parameter face "
"r1=0 intersects the parameter face r2=1.";

COMMAND[EDGE_4]  ="set-r1=1-r2=1";  
ARGUMENT[EDGE_4]  = 1;
BRIEF_HELP[EDGE_4] = "Change the edge label on the edge where the parameter face "
"r1=1 intersects the parameter face r2=1.";

COMMAND[EDGE_5]  ="set-r1=0-r3=0";  
ARGUMENT[EDGE_5]  = 1;
BRIEF_HELP[EDGE_5] = "Change the edge label on the edge where the parameter face "
"r1=0 intersects the parameter face r3=0.";

COMMAND[EDGE_6]  ="set-r1=1-r3=0";  
ARGUMENT[EDGE_6]  = 1;
BRIEF_HELP[EDGE_6] = "Change the edge label on the edge where the parameter face "
"r1=1 intersects the parameter face r3=0.";

COMMAND[EDGE_7]  ="set-r1=0-r3=1";  
ARGUMENT[EDGE_7]  = 1;
BRIEF_HELP[EDGE_7] = "Change the edge label on the edge where the parameter face "
"r1=0 intersects the parameter face r3=1.";

COMMAND[EDGE_8]  ="set-r1=1-r3=1";  
ARGUMENT[EDGE_8]  = 1;
BRIEF_HELP[EDGE_8] = "Change the edge label on the edge where the parameter face "
"r1=1 intersects the parameter face r3=1.";

COMMAND[EDGE_9]  ="set-r2=0-r3=0";  
ARGUMENT[EDGE_9]  = 1;
BRIEF_HELP[EDGE_9] = "Change the edge label on the edge where the parameter face "
"r2=0 intersects the parameter face r3=0.";

COMMAND[EDGE_10]  ="set-r2=1-r3=0";  
ARGUMENT[EDGE_10]  = 1;
BRIEF_HELP[EDGE_10] = "Change the edge label on the edge where the parameter face "
"r2=1 intersects the parameter face r3=0.";

COMMAND[EDGE_11]  ="set-r2=0-r3=1";  
ARGUMENT[EDGE_11]  = 1;
BRIEF_HELP[EDGE_11] = "Change the edge label on the edge where the parameter face "
"r2=0 intersects the parameter face r3=1.";

COMMAND[EDGE_12]  ="set-r2=1-r3=1";  
ARGUMENT[EDGE_12]  = 1;
BRIEF_HELP[EDGE_12] = "Change the edge label on the edge where the parameter face "
"r2=1 intersects the parameter face r3=1.";

COMMAND[SHOW]  ="show";       
ARGUMENT[SHOW]  = 0;
BRIEF_HELP[SHOW] = "Print the surface label on all faces of the grid.";

COMMAND[HELP]  ="help";       
ARGUMENT[HELP]  = 0;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT]  ="exit";       
ARGUMENT[EXIT]  = 0;
BRIEF_HELP[EXIT] = "Stop changing the surface label and exit from this "
"command level.";
