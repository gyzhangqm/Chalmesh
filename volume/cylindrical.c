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
set_cylindrical(input_output *io_ptr, volume_mapping *volume,
		linked_list *surface_list);
static void
delete_cylindrical_data( void *volume_data_ptr );
static int 
cylindrical_mapping( grid_point *gp, void *data_ptr);
static int 
inverse_cylindrical_mapping( grid_point *gp, void *data_ptr);
/* mappings for smoothed polygon grids */

static void 
set_cylindrical(input_output *io_ptr, volume_mapping *volume,
		linked_list *surface_list){
  cylindrical_data *info;
  int quit=0, replot=1, icom;
  const real eps=1.e-4;
    
#include "set_cylindrical_com.h"

/* cast the data pointer to the right type */
  info = (cylindrical_data *) volume->volume_data_ptr;

  do{

    if (replot){
/* plot the grid */
      if (ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
			 ogl_length_scale(volume->bb))){
	draw_volume( volume, volume->plot_mode, NULL );
	if (volume->plot_mode & 256) ogl_bb_cage( volume->bb, OGL_WHITE );
	ogl_end_plot();
	replot = 0;
      }
    }

    switch( get_command( io_ptr, "cylindrical mapping>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

/*     case TEST_INVERSE: */
/*       dr = 1.0/(volume->r1_points - 1); */
/*       ds = 1.0/(volume->r2_points - 1); */
/*       dt = 1.0/(volume->r3_points - 1); */
/*       error = FALSE; */
/*       for (gp.r = 0.0; gp.r <= 1.0 + 0.5*dr; gp.r+=dr) */
/* 	for (gp.s = 0.0; gp.s <= 1.0 + 0.5*ds; gp.s+=ds) */
/* 	  for (gp.t = 0.0; gp.t <= 1.0 + 0.5*dt; gp.t+=dt){ */
/* call the forward mapping */
/* 	    cylindrical_mapping( &gp, info ); */
/* 	    gp2.x = gp.x; gp2.y = gp.y; gp2.z = gp.z; */
/* call the inverse mapping */
/*	    inverse_cylindrical_mapping( &gp2, info );*/
/* compare */
/* 	    if (fabs(gp2.r-gp.r) + fabs(gp2.s-gp.s) + fabs(gp2.t-gp.t) >  */
/* 		3*NEWTON_EPS){ */
/* 	      error = TRUE; */
/* 	      printf("Warning, inaccurate inverse at (r,s,t)= (%e,%e,%e).\n" */
/* 		     "Inverse mapping gave (%e,%e,%e).\n", gp.r, gp.s, gp.t,  */
/* 		     gp2.r, gp2.s, gp2.t); */
/* 	    } */
/* 	  } */
/*       if (!error) printf("The inverse mapping seems to be accuracte!\n"); */
	    
/*       break; */

    case SHOW:
      printf("Minimum radius: %e, Maximum radius: %e\n", info->r_min, info->r_max);
      printf("Thickness: %e\n", info->thickness);
      printf("Starting angle: %e, Ending angle: %e\n", info->theta_min*180.0/((double)M_PI), 
	     info->theta_max*180.0/((double)M_PI));
      printf("x-center: %e, y-center: %e, z-center: %e\n", 
	     info->x_center, info->y_center, info->z_center);
      printf("r1-points: %i\nr2-points: %i\nr3-points: %i\n",
	     volume->r1_points, volume->r2_points, volume->r3_points);
      printf("\n");
      break;

    case THICKNESS:
      info->thickness = real_max( eps, get_real( io_ptr, "thickness:", 
						info->thickness, NO_INDENT));
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case R_MIN:
      info->r_min = real_max( eps, get_real( io_ptr, "Minimum radius:", 
						info->r_min, NO_INDENT));
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case R_MAX:
      info->r_max = real_max( 2*eps, get_real( io_ptr, "Maximum radius:", 
						info->r_max, NO_INDENT));
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case THETA_MIN:
      info->theta_min = get_real(io_ptr, "Starting angle (degrees):", 
				 info->theta_min*180.0/((double)M_PI), NO_INDENT)*((double)M_PI)/180.0;
/* recompute bounding box */
      volume_bb( volume );
/* check periodicity */
      volume->r2_period = (fabs(info->theta_max - info->theta_min) + 1.e-7 >= 
			   2.0*((double)M_PI))? 1 : 0;
      replot = 1;
      break;

    case THETA_MAX:
      info->theta_max = get_real(io_ptr, "Ending angle (degrees):", 
				 info->theta_max*180.0/((double)M_PI), NO_INDENT)*((double)M_PI)/180.0;
/* recompute bounding box */
      volume_bb( volume );
/* check periodicity */
      volume->r2_period = (fabs(info->theta_max - info->theta_min) + 1.e-7 >= 
			   2.0*((double)M_PI))? 1 : 0;
      replot = 1;
      break;

    case X_CENTER:
      info->x_center = get_real( io_ptr, "X-center:", info->x_center, NO_INDENT);
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case Y_CENTER:
      info->y_center = get_real( io_ptr, "Y-center:", info->y_center, NO_INDENT);
/* recompute bounding box */
      volume_bb( volume );
      replot = 1;
      break;

    case Z_CENTER:
      info->z_center = get_real( io_ptr, "Z-center:", info->z_center, NO_INDENT);
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

/* change the plot mode */
    case PLOT_MODE:
       volume->plot_mode = volume_plot_mode(io_ptr, NULL, volume, 
					   volume->plot_mode, NULL); 
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
cylindrical_mapping( grid_point *gp, void *data_ptr){
  cylindrical_data *info;
  real r, theta, r_r, theta_s;
/* cast the data to the right type */
  info = (cylindrical_data *) data_ptr;

  r = info->r_min + (info->r_max - info->r_min)*gp->r;
  theta = info->theta_min + (info->theta_max - info->theta_min)*gp->s;

  r_r = (info->r_max - info->r_min);
  theta_s = (info->theta_max - info->theta_min);

/* Cartesian coordinates */
  gp->y = info->y_center + r * cos( theta );
  gp->z = info->z_center + r * sin( theta );
  gp->x = info->x_center + info->thickness * gp->t;

/* Jacobian */
  gp->yr = r_r * cos( theta );
  gp->ys =-r * sin( theta ) * theta_s;
  gp->yt = 0.0;

  gp->zr = r_r * sin( theta );
  gp->zs = r * cos( theta ) * theta_s;
  gp->zt = 0.0;

  gp->xr = 0.0;
  gp->xs = 0.0;
  gp->xt = info->thickness;

  return OK;
}

static int 
inverse_cylindrical_mapping( grid_point *gp, void *data_ptr){
  cylindrical_data *info;
  real r, theta, r_r, theta_s;
/* cast the data to the right type */
  info = (cylindrical_data *) data_ptr;

/* Polar coordinates */
  r = sqrt((gp->y - info->y_center) * (gp->y - info->y_center) +
	   (gp->z - info->z_center) * (gp->z - info->z_center));
/* theta is 2-pi periodic, and by the following construction, 0<= theta <= 2*pi */
  if (r>NEWTON_EPS){
    if (gp->z - info->z_center >= 0.0)
      theta = acos( (gp->y - info->y_center)/r );
    else
      theta = 2.0*((double)M_PI) - acos( (gp->y - info->y_center)/r );
  }
  else{
    theta = 0.0;
  }
/* sometimes, theta_min and/or theta_max are negative, we then need to add/subtract */
/* 2*pi from theta */
  if (info->theta_min < info->theta_max){
    if (info->theta_max+NEWTON_EPS < theta)
      theta -= 2.0*((double)M_PI);
    else if (theta < info->theta_min-NEWTON_EPS)
      theta += 2.0*((double)M_PI);
  }
  else{
    if (info->theta_min+NEWTON_EPS < theta)
      theta -= 2.0*((double)M_PI);
    else if (theta < info->theta_max-NEWTON_EPS)
      theta += 2.0*((double)M_PI);
  }

/* get the parameters */
  gp->s = (theta - info->theta_min)/(info->theta_max - info->theta_min);
  gp->r = (r - info->r_min)/(info->r_max - info->r_min);
  gp->t = (gp->x - info->x_center)/info->thickness;

/* Jacobian */
  r_r = (info->r_max - info->r_min);
  theta_s = (info->theta_max - info->theta_min);

  gp->yr = r_r * cos( theta );
  gp->ys =-r * sin( theta ) * theta_s;
  gp->yt = 0.0;

  gp->zr = r_r * sin( theta );
  gp->zs = r * cos( theta ) * theta_s;
  gp->zt = 0.0;

  gp->xr = 0.0;
  gp->xs = 0.0;
  gp->xt = info->thickness;

  return OK;
}


cylindrical_data *
init_cylindrical_grid(volume_mapping *volume){
  cylindrical_data *info;
  
/* allocate memory */
  info = (cylindrical_data *) malloc( sizeof(cylindrical_data) );

/* general data */
  volume->type = "Cylindrical grid"; 
  volume->inverse_known = TRUE;
  volume->volume_data_ptr = (void *) info; 
  volume->forward_mapping = cylindrical_mapping; 
  volume->inverse_mapping = inverse_cylindrical_mapping;
  volume->set_volume = set_cylindrical; 
  volume->delete_volume_data = delete_cylindrical_data;
  
/* initial number of grid points */

  volume->r1_points = 5;
  volume->r2_points = 10;
  volume->r3_points = 5;

/* specific data */
  info->thickness = 1.0; 
  info->r_min = 0.1;
  info->r_max = 1.0;
  info->theta_min = 0.0;
  info->theta_max = ((double)M_PI);
  info->x_center = 0.0;
  info->y_center = 0.0;
  info->z_center = 0.0;

  volume_bb( volume );

  return info;
}

static void
delete_cylindrical_data( void *volume_data_ptr ){
  cylindrical_data *info;
/* cast the data pointer to the right type */
  info = (cylindrical_data *) volume_data_ptr;

/* free everything */
  free(info);
}

