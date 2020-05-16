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

int
overlapping_3d_plot_mode(input_output *io_ptr, overlapping_3d_grid *over3d, 
			 int plot_mode){
  int quit=0, replot=0, icom, i;
  linked_list_member *this_volume, *this_vgrid;
  component_vgrid *vgrid;

#include "overlapping_3d_plot_mode_com.h"

  do{

    if (replot){
/* plot the overlapping grid */
      if (ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, 
			 ogl_length_scale(over3d->bb))){
 	if (plot_mode & 2){ 
 	  for (this_volume = over3d->mapping_list->first; this_volume != NULL; 
 	       this_volume = this_volume->next){ 
 	    draw_volume( (volume_mapping *) this_volume->data, 2, NULL ); 
 	  } 
 	} 
	draw_overlapping_3d_grid( over3d, plot_mode );
	if (plot_mode & 256) ogl_bb_cage( over3d->bb, OGL_WHITE );
	ogl_end_plot();
	replot = 0;
      }
    }

    switch( get_command( io_ptr, "overlapping plot mode>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL, SAVE_ON_COPY, ARGUMENT)) 
	     == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
      quit = 1;
      break;

    case SHOW_PLOT_MODE:
      printf("\nCurrent plot_mode:\n");
      printf("%s is %s\n", COMMAND[USED_BOUNDARY_FACES], (plot_mode & 1)? "ON": "off");
      printf("%s is %s\n", COMMAND[BOUNDARY_LINES], (plot_mode & 2)? "ON": "off");
      printf("%s is %s\n", COMMAND[CHANGE_SURFACE], (plot_mode & 4)? "ON": "off");
      printf("%s is %s\n", COMMAND[USED_CELLS], (plot_mode & 8)? "ON": "off");
      printf("%s is %s\n", COMMAND[USED_FACES], (plot_mode & 16)? "ON": "off");
      printf("%s is %s\n", COMMAND[USED_GRIDLINES], (plot_mode & 32)? "ON": "off");
      printf("%s is %s\n", COMMAND[INTERP_FACES], (plot_mode & 64)? "ON": "off");
      printf("%s is %s\n", COMMAND[BAD_POINTS], (plot_mode & 128)? "ON": "off");
      printf("%s is %s\n", COMMAND[COORD_CAGE], (plot_mode & 256)? "ON": "off");
      printf("%s is %s\n", COMMAND[INTERP_POINTS], (plot_mode & 512)? "ON": "off");
/*      printf("%s is %s\n", COMMAND[DUMMY_POINTS], (plot_mode & 1024)? "ON": "off");*/
      printf("%s is %s\n", COMMAND[DEAD_POINTS], (plot_mode & 2048)? "ON": "off");
      printf("%s is %s\n", COMMAND[POSITIVE_BC], (plot_mode & 4096)? "ON": "off");
      printf("%s is %s\n", COMMAND[WORLD_POINTS], (plot_mode & 8192)? "ON": "off");
      printf("%s is %s\n", COMMAND[EDGE_LABEL], (plot_mode & 16384)? "ON": "off");
      break;

    case STD_VIEW:
      ogl_standard_view( OVERLAP_W, over3d->bb );
      replot = 1;
      break;

    case SPHERE_SIZE:
      over3d->sphere_size = real_max(1.e-7, 
				     get_real(io_ptr, "relative size of spheres:", 
				     over3d->sphere_size, LEVEL+1));
      replot = 1;
      break;

    case USED_BOUNDARY_FACES:
      if (plot_mode & 1)
	plot_mode &= ~1;
      else
	plot_mode |= 1;
      printf("Toggled %s.\n", (plot_mode & 1)? "on" : "off");
      replot = 1;
      break;

    case BOUNDARY_LINES:
      if (plot_mode & 2)
	plot_mode &= ~2;
      else
	plot_mode |= 2;
      printf("Toggled %s.\n", (plot_mode & 2)? "on" : "off");
      replot = 1;
      break;

    case CHANGE_SURFACE:
      if (plot_mode & 4)
	plot_mode &= ~4;
      else
	plot_mode |= 4;
      printf("Toggled %s.\n", (plot_mode & 4)? "on" : "off");
      replot = 1;
      break;
      
    case USED_CELLS:
      if (plot_mode & 8)
	plot_mode &= ~8;
      else
	plot_mode |= 8;
      printf("Toggled %s.\n", (plot_mode & 8)? "on" : "off");
      replot = 1;
      break;
      
    case USED_FACES:
      if (plot_mode & 16)
	plot_mode &= ~16;
      else
	plot_mode |= 16;
      printf("Toggled %s.\n", (plot_mode & 16)? "on" : "off");
      replot = 1;
      break;
      
    case USED_GRIDLINES:
      if (plot_mode & 32)
	plot_mode &= ~32;
      else
	plot_mode |= 32;
      printf("Toggled %s.\n", (plot_mode & 32)? "on" : "off");
      replot = 1;
      break;
      
    case INTERP_FACES:
      if (plot_mode & 64)
	plot_mode &= ~64;
      else
	plot_mode |= 64;
      printf("Toggled %s.\n", (plot_mode & 64)? "on" : "off");
      replot = 1;
      break;

    case BAD_POINTS:
      if (plot_mode & 128)
	plot_mode &= ~128;
      else{
	plot_mode |= 128;
	printf("The color coding is as follows:\n"
	       "%s: Bad discretization points,\n"
	       "%s: Dead interpolation locations,\n"
	       "%s: Non-explicit interpolation locations,\n"
	       "%s: Bad interpolation point becaue of dead interpolation locations,\n"
	       "%s: Bad interpolation point because of non-explicit interpolation "
	       "locations.\n",
	       ogl_color_name( BAD_DISC ), 
	       ogl_color_name( DEAD_INTERP_LOC ),
	       ogl_color_name( NON_EXPLICIT_LOC ), 
	       ogl_color_name( BAD_INTERP_DEAD_LOC ),
	       ogl_color_name( BAD_INTERP_NON_EXPLICIT ));

      }
      printf("Toggled %s.\n", (plot_mode & 128)? "on" : "off");
      replot = 1;
      break;

    case INTERP_POINTS:
      if (plot_mode & 512)
	plot_mode &= ~512;
      else
	plot_mode |= 512;
      printf("Toggled %s.\n", (plot_mode & 512)? "on" : "off");
      replot = 1;
      break;

/*     case DUMMY_POINTS: */
/*       if (plot_mode & 1024) */
/* 	plot_mode &= ~1024; */
/*       else */
/* 	plot_mode |= 1024; */
/*       printf("Toggled %s.\n", (plot_mode & 1024)? "on" : "off"); */
/*       replot = 1; */
/*       break; */

    case DEAD_POINTS:
      if (plot_mode & 2048)
	plot_mode &= ~2048;
      else
	plot_mode |= 2048;
      printf("Toggled %s.\n", (plot_mode & 2048)? "on" : "off");
      replot = 1;
      break;

    case POSITIVE_BC:
      if (plot_mode & 4096)
	plot_mode &= ~4096;
      else
	plot_mode |= 4096;
      printf("Toggled %s.\n", (plot_mode & 4096)? "on" : "off");
      replot = 1;
      break;

    case WORLD_POINTS:
      if (plot_mode & 8192)
	plot_mode &= ~8192;
      else
	plot_mode |= 8192;
      printf("Toggled %s.\n", (plot_mode & 8192)? "on" : "off");
      replot = 1;
      break;

    case EDGE_LABEL:
      printf("Sorry, not implemented. Use the edge-label under `inspect-surface-grid' "
	     "instead.\n");
/*       if (plot_mode & 16384) */
/* 	plot_mode &= ~16384; */
/*       else */
/* 	plot_mode |= 16384; */
/*       printf("Toggled %s.\n", (plot_mode & 16384)? "on" : "off"); */
/*       replot = 1; */
      break;

    case CLIP:
      over3d->n_clip_planes = 
	int_min(4, int_max(0, get_int(io_ptr, "Number of clip "
				      "planes (>=0, <=4):", 
				      over3d->n_clip_planes, NO_INDENT)));
      for (i=0; i<over3d->n_clip_planes; i++){
	printf("Specify plane %i:\n", i);
	over3d->clip_plane[i][0] = 
	  get_real(io_ptr, "x-plane:", over3d->clip_plane[i][0], LEVEL+1);
	over3d->clip_plane[i][1] = 
	  get_real(io_ptr, "y-plane:", over3d->clip_plane[i][1], LEVEL+1);
	over3d->clip_plane[i][2] = 
	  get_real(io_ptr, "z-plane:", over3d->clip_plane[i][2], LEVEL+1);
      }
      replot = 1;
      break;

    case NO_CLIP:
      over3d->n_clip_planes = 0;
      replot = 1;
      break;

    case SHOW_ONE_COMPONENT:
/* first remove all components */
      for (this_vgrid = over3d->grid_list->first; this_vgrid != NULL;
	   this_vgrid = this_vgrid->next){
	vgrid = this_vgrid->data;
	vgrid->plot_it = 0;
      }
/* then add one */
      if ((this_vgrid = get_vgrid_ptr( io_ptr, over3d->grid_list )) != NULL){
	vgrid = this_vgrid->data;
	vgrid->plot_it = 1;
	replot = 1;
      }
      break;      

    case ADD_ONE_COMPONENT:
      if ((this_vgrid = get_vgrid_ptr( io_ptr, over3d->grid_list )) != NULL){
	vgrid = this_vgrid->data;
	vgrid->plot_it = 1;
	replot = 1;
      }
      break;      

    case REMOVE_ONE_COMPONENT:
      if ((this_vgrid = get_vgrid_ptr( io_ptr, over3d->grid_list )) != NULL){
	vgrid = this_vgrid->data;
	vgrid->plot_it = 0;
	replot = 1;
      }
      break;      

    case SHOW_ALL_COMPONENTS:
      for (this_vgrid = over3d->grid_list->first; this_vgrid != NULL;
	   this_vgrid = this_vgrid->next){
	vgrid = this_vgrid->data;
	vgrid->plot_it = 1;
      }
      replot = 1;
      break;      

    case REMOVE_ALL_COMPONENTS:
      for (this_vgrid = over3d->grid_list->first; this_vgrid != NULL;
	   this_vgrid = this_vgrid->next){
	vgrid = this_vgrid->data;
	vgrid->plot_it = 0;
      }
      replot = 1;
      break;      

    case COORD_CAGE:
      if (plot_mode & 256)
	plot_mode &= ~256;
      else
	plot_mode |= 256;
      printf("Toggled %s.\n", (plot_mode & 256)? "on" : "off");
      replot = 1;
      break;

    default:
      break;
    }
  } while (!quit);

  return plot_mode;

/* done */
}
