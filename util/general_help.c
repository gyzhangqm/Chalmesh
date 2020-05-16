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
#include <stdio.h>
extern int printf(const char *format, ...);

void general_help( void ){
  printf(
"The user interface employs a tcsh-style command interpreter with the following features:\n"
" o  TAB completion of commands.\n"
" o  Upwards arrow or ^P traverses up in the command history list.\n"
" o  Downwards arrow or ^N goes down in the command history list.\n"
" o  Only the unique part between each `-' of a command needs to be given.\n"
" o  ? lists all possible commands.\n"
" o  !cmnd executes the shell command `cmnd'\n\n"
"The graphics windows have the following built in functionality:\n"
" o  Left mouse button zooms in.\n"
" o  Middle mouse button shows all.\n"
" o  Right mouse button zooms out.\n\n"
);
}
