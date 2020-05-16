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
#include "overlap_internal.h"

void 
overlap_3d_parameters(input_output *io_ptr, overlapping_3d_grid *over3d){
  char *prompt;
  int icom, quit, all_widths;
  real new_max_b_dist;

#include "overlap_3d_parameters_com.h"

  prompt = "overlap parameters>";
  quit = 0;
  do{

    switch (get_command(io_ptr, prompt, COMMAND, LAST_COM, LEVEL,
		       SAVE_ON_COPY, ARGUMENT)) {
/* show the current parameters and list */
    case SHOW:
      show_overlap_parameters( over3d );
      break;

    case CORNER_W:
      over3d->corner_width = 
	int_max(1, get_int(io_ptr, "corner width > 0: ", 
			   over3d->corner_width, NO_INDENT));
      break;

    case DISC_W:
      over3d->disc_width = 
	int_max(1, get_int(io_ptr, 
			   "discretization width > 0: ", 
			   over3d->disc_width, NO_INDENT));
      break;

    case EXTRA:
      over3d->extra = 
	int_max(0, get_int(io_ptr, 
			   "number of ghostpoints >= 0: ", 
			   over3d->extra, NO_INDENT));
      break;

    case EXPLICIT: 
      over3d->interp_type = 'e';
      break;

    case IMPLICIT:
      over3d->interp_type = 'i';
      break;

    case INTERP_W:
      over3d->interp_width = 
	int_max(1, get_int(io_ptr, 
			   "interpolation width > 0: ", 
			   over3d->interp_width, NO_INDENT));
      break;

    case NORMAL_W:
      over3d->normal_width = 
	int_max(1, get_int( io_ptr, "normal width > 0: ", 
			   over3d->normal_width, NO_INDENT));
      break;

    case TANG_W:
      over3d->tangent_width = 
	int_max(1, get_int(io_ptr, "tangential width > 0: ", 
			   over3d->tangent_width, NO_INDENT));
      break;

    case ALL_W:
      all_widths = (over3d->disc_width + over3d->interp_width + 
		    over3d->normal_width + over3d->tangent_width +
		    over3d->corner_width)/5.0;
      all_widths = 
	int_max(1, get_int(io_ptr, "all widths > 0: ", all_widths, NO_INDENT));
/* assign */
      over3d->disc_width    =  all_widths;
      over3d->tangent_width =  all_widths;
      over3d->normal_width  =  all_widths;
      over3d->corner_width  =  all_widths;
      over3d->interp_width  =  all_widths;

      break;

    case NORMAL_DIST:
      new_max_b_dist = get_real( io_ptr, "extra boundary missfit tolerance >= 0: ", 
				over3d->max_b_dist, NO_INDENT);
      while( new_max_b_dist < 0 ){
	printf("Input error: the boundary missfit tolerance must be positive!\n");
	new_max_b_dist = get_real( io_ptr, "extra boundary missfit tolerance >= 0: ", 
				  over3d->max_b_dist, NO_INDENT);
      }
      over3d->max_b_dist = new_max_b_dist;
      break;

    case HELP:
      while ( (icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				  COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1 );
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

  } while(!quit);
}

void 
show_overlap_parameters( overlapping_3d_grid *over3d ){
  printf("\n");
  printf("Corner width: %i\n",over3d->corner_width);
  printf("Discretization width: %i\n",over3d->disc_width);
  printf("Number of ghostpoints: %i\n",over3d->extra);
  printf("Interpolation width: %i\n",over3d->interp_width);
  printf("Normal width: %i\n",over3d->normal_width);
  printf("Tangential width: %i\n",over3d->tangent_width);
  
  printf("\n");
  printf("The interpolation type is ");
  if (over3d->interp_type == 'e')
    printf("explicit\n");
  else if (over3d->interp_type == 'i')
    printf("implicit\n");
  else
    printf("unknown\n");

  printf("\n");
  printf("Boundary missfit tolerance: %e\n", over3d->max_b_dist);

}

