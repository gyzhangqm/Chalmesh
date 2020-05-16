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

char *COMMAND[7], *BRIEF_HELP[7];
int ARGUMENT[7], *SAVE_ON_COPY=NULL; 

enum {ROT_AXIS, PRE_TRANS, POST_TRANS, RESET_TRANS, SHOW, 
      HELP, EXIT};

const int LAST_COM=EXIT, LEVEL=1, NO_INDENT=0;

COMMAND[ROT_AXIS] = "rotation";  
ARGUMENT[ROT_AXIS] = 0;
BRIEF_HELP[ROT_AXIS] = "Specify the axis around which the mapping should be rotated "
"as well as the rotation angle around that axis. Note that multiple rotations are "
"accumulated and that the rotation axis is specified with respect to the current "
"coordinate system (which is shown on the screen). Furthermore, the mapping is always "
"rotatated around the origin in the current (x,y,z) system.";

COMMAND[PRE_TRANS] = "pre-translation";  
ARGUMENT[PRE_TRANS] = 0;
BRIEF_HELP[PRE_TRANS] = "Specify the translation vector for the mapping (before the "
"mapping has been rotated).";

COMMAND[POST_TRANS] = "post-translation";  
ARGUMENT[POST_TRANS] = 0;
BRIEF_HELP[POST_TRANS] = "Specify the translation vector for the mapping (after the "
"mapping has been rotated).";

COMMAND[RESET_TRANS] = "reset-transformation";  
ARGUMENT[RESET_TRANS] = 0;
BRIEF_HELP[RESET_TRANS] = "Reset the rotation matrix to the identity matrix and the "
"translation vectors to zero.";

COMMAND[SHOW] = "show-transformation";  
ARGUMENT[SHOW] = 0;
BRIEF_HELP[SHOW] = "Show the rotation matrix and the translation "
"vectors.";

COMMAND[HELP]   = "help";  
ARGUMENT[HELP]   = 0;
BRIEF_HELP[HELP]   = NULL;

COMMAND[EXIT] = "exit";       
ARGUMENT[EXIT] = 0;
BRIEF_HELP[EXIT] = "Exit to the previous command level.";


