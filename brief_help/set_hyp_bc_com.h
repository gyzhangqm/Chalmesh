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
enum {LOW_R1, HIGH_R1, LOW_R2, HIGH_R2, SHOW, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[7], *BRIEF_HELP[7];
int ARGUMENT[7], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;

COMMAND[LOW_R1]  ="low-r1";  
ARGUMENT[LOW_R1]  = 1;
BRIEF_HELP[LOW_R1] = 
"Change the boundary condition on the side of the surface grid where the parameter "
"r1=0.\n"
"The following boundary conditions are available:\n"
"Free\n"
"Constant-x\n"
"Constant-y\n"
"Constant-z\n";

COMMAND[HIGH_R1] ="high-r1"; 
ARGUMENT[HIGH_R1] = 1;
BRIEF_HELP[HIGH_R1] = 
"Change the boundary condition on the side of the surface grid where the parameter "
"r1=1.\n"
"The following boundary conditions are available:\n"
"Free\n"
"Constant-x\n"
"Constant-y\n"
"Constant-z\n";

COMMAND[LOW_R2] ="low-r2"; 
ARGUMENT[LOW_R2] = 1;
BRIEF_HELP[LOW_R2] = 
"Change the boundary condition on the side of the surface grid where the parameter "
"r2=0.\n"
"The following boundary conditions are available:\n"
"Free\n"
"Constant-x\n"
"Constant-y\n"
"Constant-z\n";

COMMAND[HIGH_R2] ="high-r2"; 
ARGUMENT[HIGH_R2] = 1;
BRIEF_HELP[HIGH_R2] = 
"Change the boundary condition on the side of the surface grid where the parameter "
"r2=1.\n"
"The following boundary conditions are available:\n"
"Free\n"
"Constant-x\n"
"Constant-y\n"
"Constant-z\n";

COMMAND[SHOW]  ="show";      
ARGUMENT[SHOW]  = 0;
BRIEF_HELP[SHOW] = "Print the boundary condition on all sides of the surface grid.";

COMMAND[HELP]  ="help";      
ARGUMENT[HELP]  = 0;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT]  ="exit";      
ARGUMENT[EXIT]  = 0;
BRIEF_HELP[EXIT] = "Stop changing the boundary condition and exit from this "
"command level.";
