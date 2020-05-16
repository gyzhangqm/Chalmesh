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
#include "volume_internal.h"

void 
set_volume_bc(input_output *io_ptr, volume_mapping *volume){
  char prompt[80];
  int icom, quit=0, replot=1, plot_mode;

#include "volume_bc_com.h"

  sprintf(prompt, "%s: set boundary condition>", volume->name);

  plot_mode = 2 | 4 | 128 | 256;

  do{

    if (replot){
/* plot the mapping */
      if (ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
			 ogl_length_scale(volume->bb))){
	draw_volume( volume, plot_mode, NULL );
	ogl_end_plot();
	replot = 0;
      }
    }

    switch (get_command( io_ptr, prompt, COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT)){

    case LEFT: 
      bc_vol(1,1) = get_int(io_ptr, "Enter boundary condition: ", 
			bc_vol(1,1), NO_INDENT);
      replot = 1;
      break;

    case RIGHT: 
      bc_vol(2,1) = get_int(io_ptr, "Enter boundary condition: ", 
			bc_vol(2,1), NO_INDENT);
      replot = 1;
      break;

    case LOWER: 
      bc_vol(1,2) = get_int(io_ptr, "Enter boundary condition: ", 
			bc_vol(1,2), NO_INDENT);
      replot = 1;
      break;

    case UPPER:
      bc_vol(2,2) = get_int(io_ptr, "Enter boundary condition: ", 
			bc_vol(2,2), NO_INDENT);
      replot = 1;
      break;

    case NEAR: 
      bc_vol(1,3) = get_int(io_ptr, "Enter boundary condition: ", 
			bc_vol(1,3), NO_INDENT);
      replot = 1;
      break;

    case MY_FAR: 
      bc_vol(2,3) = get_int(io_ptr, "Enter boundary condition: ", 
			bc_vol(2,3), NO_INDENT);
      replot = 1;
      break;

    case SHOW:
      printf("Boundary condition at r1=0: %i, r1=1: %i, r2=0: %i, r2=1: %i,\n"
	     "r3=0: %i, r3=1: %i.\n",
	     bc_vol(1,1), bc_vol(2,1), 
	     bc_vol(1,2), bc_vol(2,2),
	     bc_vol(1,3), bc_vol(2,3));
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
      quit = 1;
      break;

    default:
      break;
    }
  }
  while(!quit);
}

