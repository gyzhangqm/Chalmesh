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
#define ncom 6
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] = "change-strength";        argument[0] = 1;
command[1] = "change-width";           argument[1] = 1;
command[2] = "default-strength-width"; argument[2] = 0;

command[3] = "grid-lines";             argument[3] = 1;

command[4] = "show-parameters";        argument[4] = 0;
command[5] = "help";                   argument[5] = 0;
command[6] = "exit";                   argument[6] = 0;

brief_help[0] = "Change the strength of the layer stretching at one node point.";

brief_help[1] = "Change the width of the layer stretching at one node point.";

brief_help[2] = "Use the default parameters for the layer stretching.";

brief_help[3] = "Change the number of grid lines.";

brief_help[4] = "Show the strength and width of the layer stretching at each node "
"point.";

brief_help[5] = NULL;

brief_help[6] = "Stop changing the parameters of the layer stretching and exit "
"from this command level.";

