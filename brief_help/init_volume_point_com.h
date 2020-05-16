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
enum {PLOT3D, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[3], *BRIEF_HELP[3];
int ARGUMENT[3], *SAVE_ON_COPY=NULL;
const int LEVEL=2;

COMMAND[PLOT3D] = "plot3d-ascii"; ARGUMENT[PLOT3D] = 1;
/*file[PLOT3D] = "set_volume_point_com.h";*/
COMMAND[HELP]   = "help";         ARGUMENT[HELP]   = 0;
COMMAND[EXIT]   = "cancel";       ARGUMENT[EXIT]   = 0;

BRIEF_HELP[PLOT3D] = "Read the volume grid points from a PLOT3D ascii file with the "
"following format:\n" 
" n1 n2 n3\n"
" x(1,1,1)   x(2,1,1)   ... x(n1,1,1)\n"
" x(1,2,1)   x(2,2,1)   ... x(n1,2,1)\n"
" ...\n"
" x(1,n2,1)  x(2,n2,1)  ... x(n1,n2,1)\n"
" x(1,1,2)   x(2,1,2)   ... x(n1,1,2)\n"
" ...\n"
" x(1,n2,n3) x(2,n2,n3) ... x(n1,n2,n3)\n"
" y(1,1,1)   y(2,1,1)   ... y(n1,1,1)\n"
" y(1,2,1)   y(2,2,1)   ... y(n1,2,1)\n"
" ...\n"
" y(1,n2,1)  y(2,n2,1)  ... y(n1,n2,1)\n"
" y(1,1,2)   y(2,1,2)   ... y(n1,1,2)\n"
" ...\n"
" y(1,n2,n3) y(2,n2,n3) ... y(n1,n2,n3)\n"
" z(1,1,1)   z(2,1,1)   ... z(n1,1,1)\n"
" z(1,2,1)   z(2,2,1)   ... z(n1,2,1)\n"
" ...\n"
" z(1,n2,1)  z(2,n2,1)  ... z(n1,n2,1)\n"
" z(1,1,2)   z(2,1,2)   ... z(n1,1,2)\n"
" ...\n"
" z(1,n2,n3) z(2,n2,n3) ... z(n1,n2,n3)";

BRIEF_HELP[HELP] = NULL;

BRIEF_HELP[EXIT] = "Exit from this command level without reading a grid point file.";
