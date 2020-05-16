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
/*file[PLOT3D] = "set_disc_point_com.h";*/ 
COMMAND[HELP]   = "help";         ARGUMENT[HELP]   = 0;
COMMAND[EXIT]   = "cancel";       ARGUMENT[EXIT]   = 0;

BRIEF_HELP[PLOT3D] = "Read the surface grid points from a PLOT3D ascii file with the "
"following format:\n" 
" nii njj nkk\n"
" x(1,1)  x(2,1)  ... x(nr,1)\n"
" x(1,2)  x(2,2)  ... x(nr,2)\n"
" ...\n"
" x(1,ns) x(2,ns) ... x(nr,ns)\n"
" y(1,1)  y(2,1)  ... y(nr,1)\n"
" y(1,2)  y(2,2)  ... y(nr,2)\n"
" ...\n"
" y(1,ns) y(2,ns) ... y(nr,ns)\n"
" z(1,1)  z(2,1)  ... z(nr,1)\n"
" z(1,2)  z(2,2)  ... z(nr,2)\n"
" ...\n"
" z(1,ns) z(2,ns) ... z(nr,ns)\n"
"Because chalmesh expects a surface grid, it is required that at least one of nii,"
" njj and nkk equal ONE. The first non-1 of nii and njj is set to nr and the remaining"
" non-1 is set to ns.\n"
"Example: nii = 3, njj = 1, nkk = 4 yields nr = 3 and ns = 4.";

BRIEF_HELP[HELP] = NULL;

BRIEF_HELP[EXIT] = "Exit from this command level without reading a grid point file.";
