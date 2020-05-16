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

int
volume_plot_mode(input_output *io_ptr, linked_list *volume_list,
		 volume_mapping *volume, int plot_mode, clip_planes *clip){
  int quit=0, replot=1, icom, plot_surf_label=0, plot_bc=0, i;
  linked_list_member *this_link;
  bounding_box global_bb;

#include "volume_plot_mode_com.h"

/* compute the bounding box if we have a list */
  if (volume_list != NULL)
    volume_list_bb( &global_bb, volume_list );

  do{

    if (replot){
/* redraw the volume grids */
      if (volume_list != NULL){
	if (ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
			 ogl_length_scale(&global_bb))){
	  for (this_link = volume_list->first; this_link != NULL; 
	       this_link = this_link->next)
	    draw_volume( this_link->data, plot_mode, clip);
	  if (plot_mode & 256) ogl_bb_cage( &global_bb, OGL_WHITE );
	  ogl_end_plot();
	  replot = 0;
	}
      }
      else{
	if (ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
			   ogl_length_scale(volume->bb))){
	  draw_volume( volume, plot_mode, clip );
	  if (plot_mode & 256) ogl_bb_cage( volume->bb, OGL_WHITE );
	  ogl_end_plot();
	  replot = 0;
	}
      }
    }

    switch( get_command( io_ptr, "set volume plot mode>", COMMAND, LAST_COM, LEVEL, 
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

    case GRID_POINTS:
      if (plot_mode & 1)
	plot_mode &= ~1;
      else
	plot_mode |= 1;
      replot = 1;
      break;

    case BOUNDARIES:
      if (plot_mode & 2)
	plot_mode &= ~2;
      else
	plot_mode |= 2;
      replot = 1;
      break;

    case ARROWS:
      if (plot_mode & 4)
	plot_mode &= ~4;
      else
	plot_mode |= 4;
      replot = 1;
      break;

    case STD_VIEW:
      ogl_standard_view( VOLUME_W, (volume_list == NULL)? volume->bb : 
			&global_bb );
      replot = 1;
      break;

    case CLIP:
      if (clip)
	{
	  clip->n_clip_planes = 
	    int_min(4, int_max(0, get_int(io_ptr, "Number of clip "
					  "planes (>=0, <=4):", 
					  clip->n_clip_planes, NO_INDENT)));
	  for (i=0; i<clip->n_clip_planes; i++){
	    printf("Specify plane %i:\n", i);
	    clip->clip_plane[i][0] = 
	      get_real(io_ptr, "x-plane:", clip->clip_plane[i][0], LEVEL+1);
	    clip->clip_plane[i][1] = 
	      get_real(io_ptr, "y-plane:", clip->clip_plane[i][1], LEVEL+1);
	    clip->clip_plane[i][2] = 
	      get_real(io_ptr, "z-plane:", clip->clip_plane[i][2], LEVEL+1);
	  }
	  replot = 1;
	}
      break;

    case NO_CLIP:
      if (clip)
	{
	  clip->n_clip_planes = 0;
	  replot = 1;
	}
      break;

    case R1_SURFACE:
      if (plot_mode & 8)
	plot_mode &= ~8;
      else
	plot_mode |= 8;
      replot = 1;
      break;

    case SPECIFY_R1:
      if (volume_list == NULL){
	volume->r1_surface = get_real( io_ptr, "r1-surface:", volume->r1_surface, 
				     NO_INDENT);
	volume->r1_surface = real_max(0.0, volume->r1_surface);
	volume->r1_surface = real_min(1.0, volume->r1_surface);
/* turn the surface on */
	plot_mode |= 8;
	replot = 1;
      }
      else{
	printf("Sorry, not implemented for lists of volume grids\n");
      }
      break;

    case INCREASE_R1:
      if (volume_list == NULL){
	volume->r1_surface += 1.0 / ((real) volume->r1_points - 1);
	volume->r1_surface = real_min(1.0, volume->r1_surface);
      }
      else{
	for (this_link = volume_list->first; this_link != NULL;
	     this_link = this_link->next){
	  volume = this_link->data;
	  volume->r1_surface += 1.0 / ((real) volume->r1_points - 1);
	  volume->r1_surface = real_min(1.0, volume->r1_surface);
	}
      }
/* turn the surface on */
      plot_mode |= 8;
      replot = 1;
      break;

    case DECREASE_R1:
      if (volume_list == NULL){
	volume->r1_surface -= 1.0 / ((real) volume->r1_points - 1);
	volume->r1_surface = real_max(0.0, volume->r1_surface);
      }
      else{
	for (this_link = volume_list->first; this_link != NULL;
	     this_link = this_link->next){
	  volume = this_link->data;
	  volume->r1_surface -= 1.0 / ((real) volume->r1_points - 1);
	  volume->r1_surface = real_max(0.0, volume->r1_surface);
	}
      }
/* turn the surface on */
      plot_mode |= 8;
      replot = 1;
      break;

    case R2_SURFACE:
      if (plot_mode & 16)
	plot_mode &= ~16;
      else
	plot_mode |= 16;
      replot = 1;
      break;

    case SPECIFY_R2:
      if (volume_list == NULL){
	volume->r2_surface = get_real( io_ptr, "r2-surface:", volume->r2_surface, 
				     NO_INDENT);
	volume->r2_surface = real_max(0.0, volume->r2_surface);
	volume->r2_surface = real_min(1.0, volume->r2_surface);
/* turn the surface on */
	plot_mode |= 16;
	replot = 1;
      }
      else{
	printf("Sorry, not implemented for lists of volume grids\n");
      }
      break;

    case INCREASE_R2:
      if (volume_list == NULL){
	volume->r2_surface += 1.0 / ((real) volume->r2_points - 1);
	volume->r2_surface = real_min(1.0, volume->r2_surface);
      }
      else{
	for (this_link = volume_list->first; this_link != NULL;
	     this_link = this_link->next){
	  volume = this_link->data;
	  volume->r2_surface += 1.0 / ((real) volume->r2_points - 1);
	  volume->r2_surface = real_min(1.0, volume->r2_surface);
	}
      }
/* turn the surface on */
      plot_mode |= 16;
      replot = 1;
      break;

    case DECREASE_R2:
      if (volume_list == NULL){
	volume->r2_surface -= 1.0 / ((real) volume->r2_points - 1);
	volume->r2_surface = real_max(0.0, volume->r2_surface);
      }
      else{
	for (this_link = volume_list->first; this_link != NULL;
	     this_link = this_link->next){
	  volume = this_link->data;
	  volume->r2_surface -= 1.0 / ((real) volume->r2_points - 1);
	  volume->r2_surface = real_max(0.0, volume->r2_surface);
	}
      }
/* turn the surface on */
      plot_mode |= 16;
      replot = 1;
      break;

    case R3_SURFACE:
      if (plot_mode & 32)
	plot_mode &= ~32;
      else
	plot_mode |= 32;
      replot = 1;
      break;

    case SPECIFY_R3:
      if (volume_list == NULL){
	volume->r3_surface = get_real( io_ptr, "r3-surface:", volume->r3_surface, 
				     NO_INDENT);
	volume->r3_surface = real_max(0.0, volume->r3_surface);
	volume->r3_surface = real_min(1.0, volume->r3_surface);
/* turn the surface on */
	plot_mode |= 32;
	replot = 1;
      }
      else{
	printf("Sorry, not implemented for lists of volume grids\n");
      }
      break;

    case INCREASE_R3:
      if (volume_list == NULL){
	volume->r3_surface += 1.0 / ((real) volume->r3_points - 1);
	volume->r3_surface = real_min(1.0, volume->r3_surface);
      }
      else{
	for (this_link = volume_list->first; this_link != NULL;
	     this_link = this_link->next){
	  volume = this_link->data;
	  volume->r3_surface += 1.0 / ((real) volume->r3_points - 1);
	  volume->r3_surface = real_min(1.0, volume->r3_surface);
	}
      }
/* turn the surface on */
      plot_mode |= 32;
      replot = 1;
      break;

    case DECREASE_R3:
      if (volume_list == NULL){
	volume->r3_surface -= 1.0 / ((real) volume->r3_points - 1);
	volume->r3_surface = real_max(0.0, volume->r3_surface);
      }
      else{
	for (this_link = volume_list->first; this_link != NULL;
	     this_link = this_link->next){
	  volume = this_link->data;
	  volume->r3_surface -= 1.0 / ((real) volume->r3_points - 1);
	  volume->r3_surface = real_max(0.0, volume->r3_surface);
	}
      }
/* turn the surface on */
      plot_mode |= 32;
      replot = 1;
      break;

     case SURF_LABEL: 
       if (plot_mode & 64) 
 	plot_mode &= ~64; 
       else 
 	plot_mode |= 64; 
       replot = 1; 
       break; 

     case BOUNDARY_CONDITION: 
       if (plot_mode & 128) 
 	plot_mode &= ~128; 
       else 
 	plot_mode |= 128; 
       replot = 1; 
       break; 

    case COORD_CAGE:
      if (plot_mode & 256)
	plot_mode &= ~256;
      else
	plot_mode |= 256;
      replot = 1;
      break;
      
    case CONST_SURF_LABEL: 
      plot_mode |= 512; 
      if (volume_list != NULL){
	plot_surf_label = get_int( io_ptr, "surface label to plot:", 
				  plot_surf_label, 
				  NO_INDENT);
	for (this_link = volume_list->first; this_link != NULL; 
	     this_link = this_link->next){
	  volume = this_link->data;
	  volume->plot_surf_label = plot_surf_label;
	}
      }
      else{
	volume->plot_surf_label = get_int( io_ptr, "surface label to plot:", 
					  volume->plot_surf_label, 
					  NO_INDENT);
      }
      replot = 1; 
      break; 
      
    case NO_CONST_SURF_LABEL:
      plot_mode &= ~512; 
      replot = 1;
      break;

    case CONST_BC: 
      plot_mode |= 1024; 
      if (volume_list != NULL){
	plot_bc = get_int( io_ptr, "boundary condition to plot:", 
			  plot_bc, NO_INDENT);
	for (this_link = volume_list->first; this_link != NULL; 
	     this_link = this_link->next){
	  volume = this_link->data;
	  volume->plot_bc = plot_bc;
	}
      }
      else{
	volume->plot_bc = get_int( io_ptr, "boundary condition to plot:", 
				  volume->plot_bc, 
				  NO_INDENT);
      }
      replot = 1; 
      break; 

    case NO_CONST_BC:
      plot_mode &= ~1024; 
      replot = 1;
      break;

    case EDGE_LABEL: 
      if (plot_mode & 2048) 
 	plot_mode &= ~2048; 
      else 
 	plot_mode |= 2048; 
      replot = 1; 
      break; 

    default:
      break;
    }
  } while (!quit);

/* pass back the new plot_mode */
  return plot_mode;

/* done */
}



