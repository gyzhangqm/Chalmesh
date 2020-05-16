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
set_surf_label(input_output *io_ptr, volume_mapping *volume){
  char prompt[80];
  int icom, quit=0, replot=1, plot_mode;

#include "surf_label_com.h"

  sprintf(prompt, "%s: set surface label>", volume->name);

  plot_mode = 2 | 4 | 64 | 256;

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
      surf_vol(1,1) = get_int(io_ptr, "Enter surface label: ", 
			surf_vol(1,1), NO_INDENT);
      replot = 1;
      break;

    case RIGHT: 
      surf_vol(2,1) = get_int(io_ptr, "Enter surface label: ", 
			surf_vol(2,1), NO_INDENT);
      replot = 1;
      break;

    case LOWER: 
      surf_vol(1,2) = get_int(io_ptr, "Enter surface label: ", 
			surf_vol(1,2), NO_INDENT);
      replot = 1;
      break;

    case UPPER:
      surf_vol(2,2) = get_int(io_ptr, "Enter surface label: ", 
			surf_vol(2,2), NO_INDENT);
      replot = 1;
      break;

    case NEAR: 
      surf_vol(1,3) = get_int(io_ptr, "Enter surface label: ", 
			surf_vol(1,3), NO_INDENT);
      replot = 1;
      break;

    case MY_FAR: 
      surf_vol(2,3) = get_int(io_ptr, "Enter surface label: ", 
			surf_vol(2,3), NO_INDENT);
      replot = 1;
      break;

    case SHOW:
      printf("Surface label at r1=0: %i, r1=1: %i, r2=0: %i, r2=1: %i,\n"
	     "r3=0: %i, r3=1: %i.\n",
	     surf_vol(1,1), surf_vol(2,1), 
	     surf_vol(1,2), surf_vol(2,2), 
	     surf_vol(1,3), surf_vol(2,3));
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

void 
set_edge_label(input_output *io_ptr, volume_mapping *volume){
  char prompt[80];
  int icom, quit=0, replot=1, plot_mode;

#include "edge_label_com.h"

  sprintf(prompt, "%s: set edge label>", volume->name);

  plot_mode = 2 | 4 | 2048;

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

    case EDGE_1: 
      edge_vol(1) = get_int(io_ptr, "Enter edge label: ", edge_vol(1), NO_INDENT);
      replot = TRUE;
      break;

    case EDGE_2: 
      printf("Sorry, this edge has not yet been implemented!\n");
/*       edge_vol(2) = get_int(io_ptr, "Enter edge label: ", edge_vol(2), NO_INDENT); */
/*       replot = TRUE; */
      break;

    case EDGE_3: 
      edge_vol(3) = get_int(io_ptr, "Enter edge label: ", edge_vol(3), NO_INDENT);
      replot = TRUE;
      break;

    case EDGE_4: 
      printf("Sorry, this edge has not yet been implemented!\n");
/*       edge_vol(4) = get_int(io_ptr, "Enter edge label: ", edge_vol(4), NO_INDENT); */
/*       replot = TRUE; */
      break;

    case EDGE_5: 
      edge_vol(5) = get_int(io_ptr, "Enter edge label: ", edge_vol(5), NO_INDENT);
      replot = TRUE;
      break;

    case EDGE_6: 
      edge_vol(6) = get_int(io_ptr, "Enter edge label: ", edge_vol(6), NO_INDENT);
      replot = TRUE;
      break;

    case EDGE_7: 
      edge_vol(7) = get_int(io_ptr, "Enter edge label: ", edge_vol(7), NO_INDENT);
      replot = TRUE;
      break;

    case EDGE_8: 
      printf("Sorry, this edge has not yet been implemented!\n");
/*       edge_vol(8) = get_int(io_ptr, "Enter edge label: ", edge_vol(8), NO_INDENT); */
/*       replot = TRUE; */
      break;

    case EDGE_9: 
      edge_vol(9) = get_int(io_ptr, "Enter edge label: ", edge_vol(9), NO_INDENT);
      replot = TRUE;
      break;

    case EDGE_10: 
      edge_vol(10) = get_int(io_ptr, "Enter edge label: ", edge_vol(10), NO_INDENT);
      replot = TRUE;
      break;

    case EDGE_11: 
      printf("Sorry, this edge has not yet been implemented!\n");
/*       edge_vol(11) = get_int(io_ptr, "Enter edge label: ", edge_vol(11), NO_INDENT); */
/*       replot = TRUE; */
      break;

    case EDGE_12: 
      printf("Sorry, this edge has not yet been implemented!\n");
/*       edge_vol(12) = get_int(io_ptr, "Enter edge label: ", edge_vol(12), NO_INDENT); */
/*       replot = TRUE; */
      break;

    case SHOW:
      printf("Edge label at\n"
	     "r1=0, r2=0: %i, r1=1, r2=0: %i, r1=0, r2=1: %i, r1=1, r2=1: %i,\n"
	     "r1=0, r3=0: %i, r1=1, r3=0: %i, r1=0, r3=1: %i, r1=1, r3=1: %i,\n"
	     "r2=0, r3=0: %i, r2=1, r3=0: %i, r2=0, r3=1: %i, r2=1, r3=1: %i,\n",
	     edge_vol(1), edge_vol(2), edge_vol(3), edge_vol(4), 
	     edge_vol(5), edge_vol(6), edge_vol(7), edge_vol(8), 
	     edge_vol(9), edge_vol(10), edge_vol(11), edge_vol(12));

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

