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
#include "volume.h"
#include "volume_face.h" /* surface from face of a volume mapping */

/* private member functions */
static void
volume_face_bb( surface_mapping *surface );
static int 
volume_face_mapping(grid_point *gp, surface_mapping *surface);
static void 
set_volume_face(input_output *io_ptr, surface_mapping *surface);
static void
delete_volume_face( void *surface_data_ptr );
/* mappings for a face of a volume grid mapping */

static void 
set_volume_face(input_output *io_ptr, surface_mapping *surface){
  volume_face_data *info;
  int quit=0, replot=1, icom, side, dir;
#include "set_volume_face_com.h"

/* cast the data pointer to the right type */
  info = (volume_face_data *) surface->surface_data_ptr;

  do{

    if (replot){
/* plot the surface grid */
      if (ogl_start_plot(OGL_NEW_PLOT, SURFACE_W, 
			 ogl_length_scale(surface->bb))){
	draw_surface( surface, surface->plot_mode );
	if (surface->plot_mode & 256) ogl_bb_cage( surface->bb, OGL_WHITE );
	ogl_end_plot();
	replot = 0;
      }
    }

    switch( get_command( io_ptr, "Face of volume mapping>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case SIDE:
      side = get_int( io_ptr, "Side (1 is lower, 2 is upper):", info->side, NO_INDENT);
      info->side = int_max(1, int_min(2, side));
      volume_face_bb( surface );
      replot = 1;
      break;

    case DIRECTION:
      dir = get_int( io_ptr, "Parameter direction (1, 2 or 3):", info->dir, NO_INDENT);
      info->dir = int_max(1, int_min(3, dir));

/* the correspondence between the directions in the mother grid and the component */
/* surface grid are as follows */
/* dir  r_dim   s_dim  */
/*  1   r2_dim  r3_dim */
/*  2   r3_dim  r1_dim */
/*  3   r1_dim  r2_dim */

/* copy number of grid points and periodicity */
      if (info->dir == 1){
	surface->r_points = info->volume->r2_points;
	surface->s_points = info->volume->r3_points;
	surface->r_period = info->volume->r2_period;
	surface->s_period = info->volume->r3_period;
      }
      else if (info->dir == 2){
	surface->r_points = info->volume->r3_points;
	surface->s_points = info->volume->r1_points;
	surface->r_period = info->volume->r3_period;
	surface->s_period = info->volume->r1_period;
      }
      else{
	surface->r_points = info->volume->r1_points;
	surface->s_points = info->volume->r2_points;
	surface->r_period = info->volume->r1_period;
	surface->s_period = info->volume->r2_period;
      }
/* recompute the bounding box */
      volume_face_bb( surface );
      replot = 1;
      break;

/* change the plot mode */
    case PLOT_MODE:
       surface->plot_mode = surface_plot_mode(io_ptr, NULL, surface, 
					      surface->plot_mode); 
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
volume_face_mapping(grid_point *gp, surface_mapping *surface){
  volume_face_data *info;
  grid_point volume_point;
  int status;
/* cast the data to the right type */
  info = (volume_face_data *) surface->surface_data_ptr;

/* the correspondence between the directions in the mother grid and the component */
/* surface grid are as follows */
/* dir  r_dim   s_dim  */
/*  1   r2_dim  r3_dim */
/*  2   r3_dim  r1_dim */
/*  3   r1_dim  r2_dim */

  if (info->dir == 1){
    volume_point.r = (info->side == 1)? 0.0 : 1.0;
    volume_point.s = gp->r;
    volume_point.t = gp->s;
  }
  else if (info->dir == 2){
    volume_point.r = gp->s;
    volume_point.s = (info->side == 1)? 0.0 : 1.0;
    volume_point.t = gp->r;
  }
  else{ /* info->dir == 3 */
    volume_point.r = gp->r;
    volume_point.s = gp->s;
    volume_point.t = (info->side ==1)? 0.0 : 1.0;
  }
/* Evaluate the surface */
  status = forward_volume_mapping( &volume_point, info->volume );

/* copy the result */
  gp->x = volume_point.x;
  gp->y = volume_point.y;
  gp->z = volume_point.z;

  if (info->dir == 1){
    gp->xr = volume_point.xs;
    gp->xs = volume_point.xt;
    gp->yr = volume_point.ys;
    gp->ys = volume_point.yt;
    gp->zr = volume_point.zs;
    gp->zs = volume_point.zt;
  }
  else if (info->dir == 2){
    gp->xr = volume_point.xt;
    gp->xs = volume_point.xr;
    gp->yr = volume_point.yt;
    gp->ys = volume_point.yr;
    gp->zr = volume_point.zt;
    gp->zs = volume_point.zr;
  }
  else{ /* info->dir == 3 */
    gp->xr = volume_point.xr;
    gp->xs = volume_point.xs;
    gp->yr = volume_point.yr;
    gp->ys = volume_point.ys;
    gp->zr = volume_point.zr;
    gp->zs = volume_point.zs;
  }

  return status;
}


volume_face_data *
init_volume_face(surface_mapping *surface, volume_mapping *volume){
  volume_face_data *info;
  
/* the correspondence between the directions in the mother grid and the component */
/* surface grid are as follows */
/* dir  r_dim   s_dim  */
/*  1   r2_dim  r3_dim */
/*  2   r3_dim  r1_dim */
/*  3   r1_dim  r2_dim */

/* allocate memory */
  info = (volume_face_data *) malloc( sizeof(volume_face_data) );

/* general data */
  surface->type = "Face of volume grid"; 
  surface->surface_data_ptr = (void *) info; 
  surface->forward_mapping = volume_face_mapping; 
  surface->inverse_mapping = NULL;
  surface->set_surface = set_volume_face; 
  surface->delete_surface_data = delete_volume_face;
  
/* mark that this surface relies on the volume for its definition */
  volume->used_by_surface++;

/* the default face is r3=0 */
  info->volume = volume;
  info->dir = 3; info->side = 1;

/* copy dimensions from the volume */
  surface->r_points = volume->r1_points;
  surface->s_points = volume->r2_points;
/* copy periodicity flag */
  surface->r_period = volume->r1_period;
  surface->s_period = volume->r2_period;

/* copy color from volume */
  surface->color = volume->color;

  volume_face_bb( surface );

  return info;
}

static void
delete_volume_face( void *surface_data_ptr ){
  volume_face_data *info;
  
/* cast the pointer */
  info = (volume_face_data *) surface_data_ptr;

/* mark that the volume is no longer used by the surface mapping */
  info->volume->used_by_surface--;

/* free the structure */
  free(info);
}

static void
volume_face_bb( surface_mapping *surface ){
  grid_point gp;
  const real step=0.05;

/* initialize bounding box */
  surface->bb->x_min =  1.e10;
  surface->bb->x_max = -1.e10;
  surface->bb->y_min =  1.e10;
  surface->bb->y_max = -1.e10;
  surface->bb->z_min =  1.e10;
  surface->bb->z_max = -1.e10;

  for (gp.s = 0.0; gp.s <= 1.0+0.5*step; gp.s += step)
    for (gp.r=0.0; gp.r<=1.0+0.5*step; gp.r+=step){
      volume_face_mapping( &gp, surface );
      surface->bb->x_min = real_min( surface->bb->x_min, gp.x );
      surface->bb->x_max = real_max( surface->bb->x_max, gp.x );
      surface->bb->y_min = real_min( surface->bb->y_min, gp.y );
      surface->bb->y_max = real_max( surface->bb->y_max, gp.y );
      surface->bb->z_min = real_min( surface->bb->z_min, gp.z );
      surface->bb->z_max = real_max( surface->bb->z_max, gp.z );
    }
}
