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

/* prototypes */
static void 
set_cartesian_mapping(input_output *io_ptr, volume_mapping *volume, 
		      linked_list *surface_list );
static int
cartesian_mapping( grid_point *gp_ptr, void *tuut );
static int
inverse_cartesian_mapping( grid_point *gp_ptr, void *tuut );
static void
delete_cartesian_data( void *volume_data_ptr );
/* end prototypes */

cartesian_grid_info *
init_cartesian_grid(volume_mapping *volume){
  cartesian_grid_info *info;
  
  info = (cartesian_grid_info *) malloc( sizeof(cartesian_grid_info ) );

/* general data */
  volume->type = "Cartesian grid"; 
  volume->inverse_known = TRUE;
  volume->volume_data_ptr = (void *) info; 
  volume->forward_mapping = cartesian_mapping; 
  volume->inverse_mapping = inverse_cartesian_mapping;
  volume->set_volume = set_cartesian_mapping; 
  volume->delete_volume_data = delete_cartesian_data;

  volume->r1_points = 10;
  volume->r2_points = 10;
  volume->r3_points = 10;

/* specific data */

  info->x_min = 0.0;
  info->x_max = 1.0;
  info->y_min = 0.0;
  info->y_max = 1.0;
  info->z_min = 0.0;
  info->z_max = 1.0;

  volume_bb( volume );

  return info;
}

static void
delete_cartesian_data( void *volume_data_ptr ){
  cartesian_grid_info *info;
/* cast the data pointer to the right type */
  info = (cartesian_grid_info *) volume_data_ptr;

/* free everything */
  free(info);
}

static void 
set_cartesian_mapping(input_output *io_ptr, volume_mapping *volume, 
		      linked_list *surface_list ){
  cartesian_grid_info *info;
  int quit=0, replot=1, icom;
  real dx1, dx2, dx3;
#include "set_cartesian_com.h"

/* cast the data pointer to the right type */
  info = (cartesian_grid_info *) volume->volume_data_ptr;

  do{

    if (replot){
/* plot the surface grid */
      if (ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
			 ogl_length_scale(volume->bb))){
	draw_volume( volume, volume->plot_mode, NULL);
	if (volume->plot_mode & 256) ogl_bb_cage( volume->bb, OGL_WHITE );
	ogl_end_plot();
	replot = 0;
      }
    }

    switch( get_command( io_ptr, "cartesian>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case SHOW:
      printf("x-min: %e, x-max: %e\n", info->x_min, info->x_max);
      printf("y-min: %e, y-max: %e\n", info->y_min, info->y_max);
      printf("z-min: %e, z-max: %e\n", info->z_min, info->z_max);
      printf("r1-points: %i\nr2-points: %i\nr3-points: %i\n",
	     volume->r1_points, volume->r2_points, volume->r3_points);
      dx1 = (info->x_max - info->x_min)/(volume->r1_points-1);
      dx2 = (info->y_max - info->y_min)/(volume->r2_points-1);
      dx3 = (info->z_max - info->z_min)/(volume->r3_points-1);
      printf("Size of space diagonal: %e\n", sqrt( dx1*dx1 + dx2*dx2 + dx3*dx3 ));
      printf("\n");
      break;

    case X_MIN:
      info->x_min = get_real( io_ptr, "x-min:", info->x_min, NO_INDENT);
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case X_MAX:
      info->x_max = get_real( io_ptr, "x-max:", info->x_max, NO_INDENT);
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case Y_MIN:
      info->y_min = get_real( io_ptr, "y-min:", info->y_min, NO_INDENT);
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case Y_MAX:
      info->y_max = get_real( io_ptr, "y-max:", info->y_max, NO_INDENT);
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case Z_MIN:
      info->z_min = get_real( io_ptr, "z-min:", info->z_min, NO_INDENT);
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case Z_MAX:
      info->z_max = get_real( io_ptr, "z-max:", info->z_max, NO_INDENT);
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

/* change the plot mode */
    case PLOT_MODE:
       volume->plot_mode = volume_plot_mode(io_ptr, NULL, volume, 
					   volume->plot_mode, NULL); 
/* don't need to redraw the surface grids! */
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
cartesian_mapping( grid_point *gp_ptr, void *tuut ){
  cartesian_grid_info *info;

/* cast the data pointer to the right type */
  info = (cartesian_grid_info *) tuut;

  gp_ptr->x = info->x_min + gp_ptr->r * (info->x_max - info->x_min);
  gp_ptr->y = info->y_min + gp_ptr->s * (info->y_max - info->y_min);
  gp_ptr->z = info->z_min + gp_ptr->t * (info->z_max - info->z_min);

  gp_ptr->xr = info->x_max - info->x_min;
  gp_ptr->xs = 0.0;
  gp_ptr->xt = 0.0;

  gp_ptr->yr = 0.0;
  gp_ptr->ys = info->y_max - info->y_min;
  gp_ptr->yt = 0.0;

  gp_ptr->zr = 0.0;
  gp_ptr->zs = 0.0;
  gp_ptr->zt = info->z_max - info->z_min;

  return OK;
}

static int
inverse_cartesian_mapping( grid_point *gp_ptr, void *tuut ){
  cartesian_grid_info *info;

/* cast the data pointer to the right type */
  info = (cartesian_grid_info *) tuut;

  gp_ptr->r = (gp_ptr->x - info->x_min) / (info->x_max - info->x_min);
  gp_ptr->s = (gp_ptr->y - info->y_min) / (info->y_max - info->y_min);
  gp_ptr->t = (gp_ptr->z - info->z_min) / (info->z_max - info->z_min);

  gp_ptr->xr = info->x_max - info->x_min;
  gp_ptr->xs = 0.0;
  gp_ptr->xt = 0.0;

  gp_ptr->yr = 0.0;
  gp_ptr->ys = info->y_max - info->y_min;
  gp_ptr->yt = 0.0;

  gp_ptr->zr = 0.0;
  gp_ptr->zs = 0.0;
  gp_ptr->zt = info->z_max - info->z_min;

  return OK;
}
