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

/* private member functions */
static void 
set_normal_surface(input_output *io_ptr, volume_mapping *volume,
		   linked_list *surface_list);
static void
delete_normal_surface( void *volume_data_ptr );
/* mappings for smoothed polygon grids */

static void 
set_normal_surface(input_output *io_ptr, volume_mapping *volume,
		   linked_list *surface_list){
  normal_surface_data *info;
  int quit=0, replot=1, icom;
  const real eps=1.e-4;
  clip_planes clip_info = {0, {{5.0, 0.0, 0.0}, {5.0, 0.0, 0.0}, {5.0, 0.0, 0.0}, 
				 {5.0, 0.0, 0.0}}};
#include "set_normal_surface_com.h"

/* cast the data pointer to the right type */
  info = (normal_surface_data *) volume->volume_data_ptr;

  do{

    if (replot){
/* plot the surface grid */
      if (ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
			 ogl_length_scale(volume->bb))){
	draw_volume( volume, volume->plot_mode, &clip_info );
	if (volume->plot_mode & 256) ogl_bb_cage( volume->bb, OGL_WHITE );
	ogl_end_plot();
	replot = 0;
      }
    }

    switch( get_command( io_ptr, "normal surface>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case SHOW:
      printf("This mapping is based on the surface mapping `%s'\n", 
	     info->surface->name);
      printf("width: %e\n", info->width);
      printf("r1-points: %i\nr2-points: %i\nr3-points: %i\n",
	     volume->r1_points, volume->r2_points, volume->r3_points);
      printf("\n");
      break;

    case WIDTH:
      info->width = real_max( eps, get_real( io_ptr, "width:", 
					    info->width, NO_INDENT));
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case R1_LINES:
      volume->r1_points = int_max(2, get_int( io_ptr, "Number of r1-lines:", 
				 volume->r1_points, NO_INDENT));
      replot = 1;
      break;

    case R2_LINES:
      volume->r2_points = int_max(2, get_int( io_ptr, "Number of r2-lines:", 
				 volume->r2_points, NO_INDENT));
      replot = 1;
      break;

    case R3_LINES:
      volume->r3_points = int_max(2, get_int( io_ptr, "Number of r3-lines:", 
				 volume->r3_points, NO_INDENT));
      replot = 1;
      break;

    case STRETCHING:
/* normal stretching */
      if (info->normal_stretch == NULL){
	info->normal_stretch = choose_stretching(io_ptr, &volume->r3_points );
	volume->plot_mode |= 8;
	replot = TRUE;
      }
      else{
	if (set_stretching( io_ptr, info->normal_stretch, &volume->r3_points ) ){
	  replot = TRUE;
	}
	else{
	  printf("The %s cannot be changed. If you want another type of stretching,\n"
		 "delete the present stretching and try again.\n", 
		 stretching_name( info->normal_stretch ) );
	}
      }
      break;

    case NO_STRETCHING:
/* no normal stretching */
      if (info->normal_stretch != NULL){
	info->normal_stretch = delete_generic_stretching( info->normal_stretch );
	volume->plot_mode |= 8;
	replot = TRUE;
      }
      else{
	printf("There is no stretching in the r3-direction.\n");
      }
      break;

/* change the plot mode */
    case PLOT_MODE:
       volume->plot_mode = volume_plot_mode(io_ptr, NULL, volume, 
					   volume->plot_mode, &clip_info); 
/* don't need to redraw the surface grids! */
      break;

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

    default:
      break;

    }
  } while (!quit);
}


static int 
normal_surface_mapping( grid_point *gp, void *data_ptr){
  normal_surface_data *info;
  grid_point surface_point, sp;
  real n_vec[3], n_vec_r[3], n_vec_rm[3], n_vec_s[3], n_vec_sm[3], 
    t_stretch, t_ratio, dummy;
  int i;
  const real eps=1.e-4;
/* cast the data to the right type */
  info = (normal_surface_data *) data_ptr;

/* assign the parameters in the surface structure */
  surface_point.r = gp->r;
  surface_point.s = gp->s;

/* Evaluate the surface */
  if (!forward_surface_mapping( &surface_point, info->surface ))
    return ERROR;

/* compute the normal */
  if (!compute_surface_normal( n_vec, &surface_point ))
    return ERROR;

/* stretching in the normal direction? */
  uniform_to_stretch(gp->t, info->normal_stretch, &t_stretch, &t_ratio, &dummy );

/* compute the derivative of the surface normal */
  if (!derivative_surface_normal(n_vec_r, n_vec_s, gp, info->surface))
    return ERROR;

/* Cartesian coordinates */
  gp->x = surface_point.x + n_vec[0] * info->width * t_stretch;
  gp->y = surface_point.y + n_vec[1] * info->width * t_stretch;
  gp->z = surface_point.z + n_vec[2] * info->width * t_stretch;

/* Jacobian */
  gp->xr = surface_point.xr + n_vec_r[0] * info->width * t_stretch;
  gp->xs = surface_point.xs + n_vec_s[0] * info->width * t_stretch;
  gp->xt = n_vec[0] * info->width * t_ratio;

  gp->yr = surface_point.yr + n_vec_r[1] * info->width * t_stretch;
  gp->ys = surface_point.ys + n_vec_s[1] * info->width * t_stretch;
  gp->yt = n_vec[1] * info->width * t_ratio;

  gp->zr = surface_point.zr + n_vec_r[2] * info->width * t_stretch;
  gp->zs = surface_point.zs + n_vec_s[2] * info->width * t_stretch;
  gp->zt = n_vec[2] * info->width * t_ratio;

  return OK;
}


normal_surface_data *
init_normal_surface(input_output *io_ptr, volume_mapping *volume, surface_mapping *surface){
  normal_surface_data *info;
  
/* allocate memory */
  info = (normal_surface_data *) malloc( sizeof(normal_surface_data) );

/* general data */
  volume->type = "Normal surface"; 
  volume->volume_data_ptr = (void *) info; 
  volume->forward_mapping = normal_surface_mapping; 
  volume->inverse_mapping = NULL;
  volume->set_volume = set_normal_surface; 
  volume->delete_volume_data = delete_normal_surface;
  
/* copy dimensions from the surface */
  volume->r1_points = surface->r_points;
  volume->r2_points = surface->s_points;
  volume->r3_points = 5;

/* and periodicity */
  volume->r1_period = surface->r_period;
  volume->r2_period = surface->s_period;
  volume->r3_period = 0;

/* specific data */
  surface->used_by_volume++;
  info->surface = surface;
  info->width = 1.0; 
  info->normal_stretch = NULL;

  volume_bb( volume );

  return info;
}

static void
delete_normal_surface( void *volume_data_ptr ){
  normal_surface_data *info;
  
/* cast the pointer */
  info = (normal_surface_data *) volume_data_ptr;

/* mark that the surface is no longer used by the volume mapping */
  info->surface->used_by_volume--;

/* free the stretching */
  info->normal_stretch = delete_generic_stretching( info->normal_stretch );

/* free the structure */
  free(info);
}

