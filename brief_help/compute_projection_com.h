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
enum {LEFT, RIGHT, LOWER, UPPER, HELP, CANCEL};
const int LAST_COM = CANCEL;

char *COMMAND[6], *BRIEF_HELP[6];
int ARGUMENT[6], *SAVE_ON_COPY=NULL;
const int LEVEL=3;

COMMAND[LEFT]   = "project-low-r1";       
ARGUMENT[LEFT]   = 0;
BRIEF_HELP[LEFT]  = "Project the parameter face r1=0 onto the surface.";

COMMAND[RIGHT]  = "project-high-r1";       
ARGUMENT[RIGHT]  = 0;
BRIEF_HELP[RIGHT] = "Project the parameter face r1=1 onto the surface.";

COMMAND[LOWER]  = "project-low-r2";       
ARGUMENT[LOWER]  = 0;
BRIEF_HELP[LOWER] = "Project the parameter face r2=0 onto the surface.";

COMMAND[UPPER]  = "project-high-r2";       
ARGUMENT[UPPER]  = 0;
BRIEF_HELP[UPPER] = "Project the parameter face r2=1 onto the surface.";

COMMAND[HELP]   = "help";
ARGUMENT[HELP]   = 1;
BRIEF_HELP[HELP]  = NULL;

COMMAND[CANCEL] = "cancel";
ARGUMENT[CANCEL] = 0;
BRIEF_HELP[CANCEL] = "Do not project any face onto the surface. Instead exit to "
"the previous command level.";


