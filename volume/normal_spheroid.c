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
set_normal_spheroid(input_output *io_ptr, volume_mapping *volume,
		   linked_list *surface_list);
static void
delete_normal_spheroid( void *volume_data_ptr );
static int 
normal_spheroid_mapping( grid_point *gp, void *data_ptr);
static int 
inverse_normal_spheroid_mapping( grid_point *gp, void *data_ptr);
/* mappings for smoothed polygon grids */

static void 
set_normal_spheroid(input_output *io_ptr, volume_mapping *volume,
		    linked_list *surface_list){
  normal_spheroid_data *info;
  int quit=0, replot=1, icom;
  const real eps=1.e-4;
#include "set_normal_spheroid_com.h"

/* cast the data pointer to the right type */
  info = (normal_spheroid_data *) volume->volume_data_ptr;

  do{

    if (replot){
/* plot the surface grid */
      if (ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
			 ogl_length_scale(volume->bb))){
	draw_volume( volume, volume->plot_mode, NULL );
	if (volume->plot_mode & 256) ogl_bb_cage( volume->bb, OGL_WHITE );
	ogl_end_plot();
	replot = 0;
      }
    }

    switch( get_command( io_ptr, "normal spheroid>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case SHOW:
      printf("width: %e\n", info->width);
      printf("r1-points: %i\nr2-points: %i\nr3-points: %i\n",
	     volume->r1_points, volume->r2_points, volume->r3_points);
      printf("\n");
      break;

    case CENTER:
      info->x_center = get_real( io_ptr, "X-center:", info->x_center, LEVEL+1);
      info->y_center = get_real( io_ptr, "Y-center:", info->y_center, LEVEL+1);
      info->z_center = get_real( io_ptr, "Z-center:", info->z_center, LEVEL+1);
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case THETA_MIN:
      replot = 1;
      info->theta_min = get_real( io_ptr, "Starting theta angle:", 
				 info->theta_min*180.0/((double)M_PI), NO_INDENT) * ((double)M_PI)/180.0;
/* recompute bounding box */
      volume_bb( volume );
/* check periodicity */
      volume->r1_period = (fabs(info->theta_max - info->theta_min) + 1.e-7 >= 
			 2.0*((double)M_PI))? 1 : 0;
      break;

    case THETA_MAX:
      replot = 1;
      info->theta_max = get_real( io_ptr, "Ending theta angle:", 
				 info->theta_max*180.0/((double)M_PI), NO_INDENT) * ((double)M_PI)/180.0;
/* recompute bounding box */
      volume_bb( volume );
/* check periodicity */
      volume->r1_period = (fabs(info->theta_max - info->theta_min) + 1.e-7 >= 
			 2.0*((double)M_PI))? 1 : 0;
      break;

    case PHI_MIN:
      replot = 1;
      info->phi_min = real_max( eps, get_real( io_ptr, "Starting phi angle > 0:", 
				 info->phi_min*180.0/((double)M_PI), NO_INDENT) * ((double)M_PI)/180.0 );
      volume_bb( volume );
      break;

    case PHI_MAX:
      replot = 1;
      info->phi_max = real_min( ((double)M_PI)-eps, get_real( io_ptr, "Ending phi angle < 180:", 
				 info->phi_max*180.0/((double)M_PI), NO_INDENT) * ((double)M_PI)/180.0 );
      volume_bb( volume );
      break;

    case X_SEMI:
      replot = 1;
      info->x_semi = real_max( eps, get_real( io_ptr, "Semi-axis in x-direction:", 
				 info->x_semi, NO_INDENT) );
      volume_bb( volume );
      break;

    case YZ_SEMI:
      replot = 1;
      info->yz_semi = real_max( eps, get_real( io_ptr, "Semi-axis in y-direction:", 
				 info->yz_semi, NO_INDENT) );
      volume_bb( volume );
      break;

    case P_AXIS:
      replot = 1;
      info->polar_axis = get_int( io_ptr, "Polar axis 1 => x, 2 => y:", 
				 info->polar_axis, NO_INDENT);
      info->polar_axis = int_min( 2, int_max(1,info->polar_axis));
      volume_bb( volume );
      break;

    case WIDTH:
      info->width = real_max( eps, get_real( io_ptr, "width > 0:", 
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
					   volume->plot_mode, NULL);
/* don't need to redraw the volume grids! */
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
normal_spheroid_mapping( grid_point *gp, void *data_ptr){
  normal_spheroid_data *info;
  real t_stretch, t_ratio, dummy, scale, theta, phi, theta_r, phi_s;
#define sq(x) ((x)*(x))
/* cast the data to the right type */
  info = (normal_spheroid_data *) data_ptr;

/* angles */
  theta = info->theta_min + (info->theta_max - info->theta_min) * gp->r;
  phi = info->phi_min + (info->phi_max - info->phi_min) * gp->s;

  theta_r = (info->theta_max - info->theta_min);
  phi_s = (info->phi_max - info->phi_min);

/* evaluate the stretching in the normal direction */
  uniform_to_stretch(gp->t, info->normal_stretch, &t_stretch, &t_ratio, &dummy );

  if (info->polar_axis == 1){
    scale = 1.0/sqrt(sq(info->x_semi*sin(phi)) + sq(info->yz_semi*cos(phi)));

/* Cartesian coordinates */
    gp->x = info->x_center + 
      (info->x_semi + info->yz_semi*info->width*t_stretch*scale) * cos(phi);
    gp->y = info->y_center + 
      (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * sin(phi) * cos(theta);
    gp->z = info->z_center + 
      (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * sin(phi) * sin(theta);


/* Jacobian (neglect the variation of scale with phi) */
    gp->xr = 0.0;
    gp->xs =-(info->x_semi + info->yz_semi*info->width*t_stretch*scale) * 
      sin(phi) * phi_s;
    gp->xt = info->yz_semi*info->width*t_ratio*scale * cos(phi);

    gp->yr =-(info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      sin(phi) * sin(theta) * theta_r;
    gp->ys = (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      cos(phi) * cos(theta) * phi_s;
    gp->yt = info->x_semi*info->width*t_ratio*scale * sin(phi) * cos(theta);

    gp->zr = (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      sin(phi) * cos(theta) * theta_r;
    gp->zs = (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      cos(phi) * sin(theta) * phi_s;
    gp->zt = info->x_semi*info->width*t_ratio*scale * sin(phi) * sin(theta);
  }
  else if (info->polar_axis == 2){
    scale = 1.0/sqrt( (sq(info->x_semi*cos(theta)) + sq(info->yz_semi*sin(theta))) * 
		     sq(sin(phi)) + sq(info->x_semi*cos(phi)) );

/* Cartesian coordinates */
    gp->x = info->x_center + 
      (info->x_semi + info->yz_semi*info->width*t_stretch*scale) * sin(phi) * sin(theta);
    gp->y = info->y_center + 
      (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * cos(phi);
    gp->z = info->z_center + 
      (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * sin(phi) * cos(theta);

/* Jacobian (neglect the variation of scale with theta and phi) */
    gp->xr = (info->x_semi + info->yz_semi*info->width*t_stretch*scale) * 
      sin(phi) * cos(theta) * theta_r;
    gp->xs = (info->x_semi + info->yz_semi*info->width*t_stretch*scale) * 
      cos(phi) * sin(theta) * phi_s;
    gp->xt = info->yz_semi*info->width*t_ratio*scale * sin(phi) * sin(theta);

    gp->yr = 0.0;
    gp->ys =-(info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      sin(phi) * phi_s;
    gp->yt = info->x_semi*info->width*t_ratio*scale * cos(phi);

    gp->zr =-(info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      sin(phi) * sin(theta) * theta_r;
    gp->zs = (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      cos(phi) * cos(theta) * phi_s;
    gp->zt = info->x_semi*info->width*t_ratio*scale * sin(phi) * cos(theta);
  }
  else{
    printf("ERROR: normal_spheroid_mapping: unknown polar_axis = %i\n", 
	   info->polar_axis);
    return ERROR;
  }

  return OK;
#undef sq
}

static int 
inverse_normal_spheroid_mapping( grid_point *gp, void *data_ptr){
  normal_spheroid_data *info;
  real t_stretch, t_ratio, u_ratio, dummy, scale, theta, phi, theta_r, phi_s
    , x, y, z, r, res, xs, ys, zs;
  int itr;
  const int max_itr = 10;
  const real eps = 1.e-10;
#define sq(x) ((x)*(x))
#define func(phi) ((sq(info->x_semi) - sq(info->yz_semi))*0.5*sin(2.0*phi) + \
		   r*info->yz_semi*cos(phi) - x*info->x_semi*sin(phi))
#define d_func(phi) ((sq(info->x_semi) - sq(info->yz_semi))*cos(2.0*phi) - \
		   r*info->yz_semi*sin(phi) - x*info->x_semi*cos(phi))
/* cast the data to the right type */
  info = (normal_spheroid_data *) data_ptr;

/* center the coordinates */
  x = gp->x - info->x_center;
  y = gp->y - info->y_center;
  z = gp->z - info->z_center;

  r = sqrt(sq(y) + sq(z));
/* initial guess for phi */
  phi = 0.5 * ((double)M_PI);
/* Newton iteration to compute phi */
  itr = 0;
  do{
    res = func( phi );
    if (fabs(res) <= eps) break;
    phi -= res / d_func( phi );
  } while(++itr < max_itr);
  if (fabs(res=func(phi)) > eps)
    printf("Warning: inverse_normal_speroidal_mapping: Poor convergence in Newton.\n"
	   "x = %e, r = %e, res = %e\n", x, y, res);
  scale = sqrt( sq(info->x_semi*sin(phi)) + sq(info->yz_semi*cos(phi)) );
  if (fabs(sin(phi)) > .1){
    t_stretch = (r/sin(phi) - info->yz_semi) * scale / (info->x_semi * info->width);
  }
  else{
    t_stretch = (x/cos(phi) - info->x_semi) * scale / (info->yz_semi * info->width);
  }
  theta = atan2(z, y);
  if (theta < info->theta_min - 1.e-4) theta += 2.0*((double)M_PI);
  if (theta > info->theta_max + 1.e-4) theta -= 2.0*((double)M_PI);

  if (info->polar_axis == 2){
/* get the point on the surface */
    xs = info->x_semi * cos(phi);
    ys = info->yz_semi * sin(phi) * cos(theta);
    zs = info->yz_semi * sin(phi) * sin(theta);

/* compute the new surface parameters */
    phi = acos(ys/info->yz_semi);
    theta = atan2(xs/info->x_semi, zs/info->yz_semi);
    if (theta < info->theta_min - 1.e-4) theta += 2.0*((double)M_PI);
    if (theta > info->theta_max + 1.e-4) theta -= 2.0*((double)M_PI);

/* compute the new normal distance */
    scale = sqrt( (sq(info->x_semi*cos(theta)) + sq(info->yz_semi*sin(theta))) * 
		 sq(sin(phi)) + sq(info->x_semi*cos(phi)) );
    if (fabs(cos(phi)) > 0.1){
      t_stretch = (y/cos(phi) - info->yz_semi) * scale / (info->x_semi * info->width);
    }
    else if (fabs(sin(theta)) > 0.1){
      t_stretch = (x/sin(phi)/sin(theta) - info->x_semi) * scale / 
	(info->yz_semi * info->width);
    }
    else{
      t_stretch = (z/sin(phi)/cos(theta) - info->yz_semi) * scale /
	(info->x_semi * info->width);
    }
  } /* end if polar_axis == 2 */

/* compute r, s, t from theta, phi, t_stretch */
  gp->r = (theta - info->theta_min) / (info->theta_max - info->theta_min);
  gp->s = (phi - info->phi_min) / (info->phi_max - info->phi_min);

/* evaluate the stretching in the normal direction */
  stretch_to_uniform(t_stretch, info->normal_stretch, &(gp->t), &u_ratio, &dummy );


/* Jacobian */
  theta_r = (info->theta_max - info->theta_min);
  phi_s = (info->phi_max - info->phi_min);
  t_ratio = 1.0/u_ratio;
  scale = 1.0/scale;

  if (info->polar_axis == 1){
    gp->xr = 0.0;
    gp->xs =-(info->x_semi + info->yz_semi*info->width*t_stretch*scale) * 
      sin(phi) * phi_s;
    gp->xt = info->yz_semi*info->width*t_ratio*scale * cos(phi);

    gp->yr =-(info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      sin(phi) * sin(theta) * theta_r;
    gp->ys = (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      cos(phi) * cos(theta) * phi_s;
    gp->yt = info->x_semi*info->width*t_ratio*scale * sin(phi) * cos(theta);

    gp->zr = (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      sin(phi) * cos(theta) * theta_r;
    gp->zs = (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      cos(phi) * sin(theta) * phi_s;
    gp->zt = info->x_semi*info->width*t_ratio*scale * sin(phi) * sin(theta);
  }
  else if (info->polar_axis == 2){
    gp->xr = (info->x_semi + info->yz_semi*info->width*t_stretch*scale) * 
      sin(phi) * cos(theta) * theta_r;
    gp->xs = (info->x_semi + info->yz_semi*info->width*t_stretch*scale) * 
      cos(phi) * sin(theta) * phi_s;
    gp->xt = info->yz_semi*info->width*t_ratio*scale * sin(phi) * sin(theta);

    gp->yr = 0.0;
    gp->ys =-(info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      sin(phi) * phi_s;
    gp->yt = info->x_semi*info->width*t_ratio*scale * cos(phi);

    gp->zr =-(info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      sin(phi) * sin(theta) * theta_r;
    gp->zs = (info->yz_semi + info->x_semi*info->width*t_stretch*scale) * 
      cos(phi) * cos(theta) * phi_s;
    gp->zt = info->x_semi*info->width*t_ratio*scale * sin(phi) * cos(theta);
  }
  else{
    printf("ERROR: normal_spheroid_mapping: unknown polar_axis = %i\n", 
	   info->polar_axis);
    return ERROR;
  }



  return OK;
#undef sq
}


normal_spheroid_data *
init_normal_spheroid(volume_mapping *volume){
  normal_spheroid_data *info;
  
/* allocate memory */
  info = (normal_spheroid_data *) malloc( sizeof(normal_spheroid_data) );

/* general data */
  volume->type = "Normal spheroid"; 
  volume->inverse_known = TRUE; 
  volume->volume_data_ptr = (void *) info; 
  volume->forward_mapping = normal_spheroid_mapping; 
  volume->inverse_mapping = inverse_normal_spheroid_mapping;
  volume->set_volume = set_normal_spheroid; 
  volume->delete_volume_data = delete_normal_spheroid;
  
/* grid dimensions */
  volume->r1_points = 10;
  volume->r2_points = 10;
  volume->r3_points = 5;

/* periodicity */
  volume->r1_period = 0;
  volume->r2_period = 0;
  volume->r3_period = 0;

/* specific data */
  info->x_center = 0.0;
  info->y_center = 0.0;
  info->z_center = 0.0;
  info->x_semi = 2.0;
  info->yz_semi = 1.0;
  info->theta_min = 0.0;
  info->theta_max = ((double)M_PI);
  info->phi_min = 0.25 * ((double)M_PI);
  info->phi_max = 0.75 * ((double)M_PI);
  info->width = 1.0; 

  info->polar_axis = 1;

  info->normal_stretch = NULL;

/* initialize the bounding box */
  volume_bb( volume );

  return info;
}

static void
delete_normal_spheroid( void *volume_data_ptr ){
  normal_spheroid_data *info;
  
/* cast the pointer */
  info = (normal_spheroid_data *) volume_data_ptr;

/* free the stretching */
  info->normal_stretch = delete_generic_stretching( info->normal_stretch );

/* free the structure */
  free(info);
}

