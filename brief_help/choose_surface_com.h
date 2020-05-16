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
enum {SPHERE, DISCRETE, HELP, CANCEL};
const int LAST_COM=CANCEL, LEVEL=1;

char *COMMAND[4], *BRIEF_HELP[4];
int *ARGUMENT=NULL, *SAVE_ON_COPY=NULL; 

COMMAND[SPHERE] = "ellipsoid";
/*file[SPHERE] = "sphere_patch_com.h";*/ 
BRIEF_HELP[SPHERE] = "Define the surface to be a patch on an ellipsoid. Proceed "
"to modify the default parameters.";

COMMAND[DISCRETE] = "discrete-point";  
/*file[DISCRETE] = "init_disc_point_com.h";*/ 
BRIEF_HELP[DISCRETE] = "Define the mapping to interpolate inbetween a discrete set of "
"grid points, which are read from a file that for instance could be the output from "
"an external grid generator.";

COMMAND[HELP]   = "help";  
BRIEF_HELP[HELP]   = NULL;

BRIEF_HELP[CANCEL] = "Don't define a surface grid. Instead go back to the previous"
" command level.";
COMMAND[CANCEL] = "cancel";
