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
#include "surface_internal.h"

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

/* prototypes */
static int
sphere_patch_hessian(hessian_data *h_ptr, struct surface_mapping *surface);
static void
delete_sphere_patch( void *tuut );
/* end prototypes */

sphere_patch_data *
init_sphere_patch(surface_mapping *surface){
  sphere_patch_data *info;
  
  info = (sphere_patch_data *) malloc( sizeof(sphere_patch_data ) );

/* general data */
  surface->type = "Ellipsoidal patch"; 
  surface->surface_data_ptr = (void *) info; 
  surface->forward_mapping = sphere_patch; 
  surface->hessian = sphere_patch_hessian; 
  surface->inverse_mapping = NULL;
  surface->set_surface = set_sphere_patch; 
  surface->delete_surface_data = delete_sphere_patch;

  surface->r_points = 10;
  surface->s_points = 10;

/* specific data */

  info->polar_axis = 3; /* x_3-axis has polar singularity */
  info->x_semi = info->y_semi = info->z_semi = 1.0;
  info->alpha_min = 0.0;
  info->alpha_max = ((double)M_PI);
  info->beta_min  = 0.25*((double)M_PI);
  info->beta_max  = 0.75*((double)M_PI);

  info->x_center = 0.0;
  info->y_center = 0.0;
  info->z_center = 0.0;

  surface_bb( surface );

  return info;
}

static void
delete_sphere_patch( void *tuut ){
  sphere_patch_data *info;
/* cast the data pointer to the right type */
  info = (sphere_patch_data *) tuut;

/* purge everything */
  free(info);
}

void 
set_sphere_patch( input_output *io_ptr, surface_mapping *surface ){
  sphere_patch_data *info;
  int quit=0, replot=1, icom, tmp;
  const real eps=1.e-4;
#include "sphere_patch_com.h"

/* cast the data pointer to the right type */
  info = (sphere_patch_data *) surface->surface_data_ptr;

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

    switch( get_command( io_ptr, "ellipsoid>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case SHOW:
      printf("r-lines: %i\ns-lines: %i\nx-center: %e\ny-center: %e\nz-center: %e\n",
	     surface->r_points, surface->s_points, 
	     info->x_center, info->y_center, info->z_center);
      printf("alpha-min: %e\nalpha-max: %e\nbeta-min: %e\nbeta-max: %e\n",
	     info->alpha_min*180.0/((double)M_PI), info->alpha_max*180.0/((double)M_PI), 
	     info->beta_min*180.0/((double)M_PI), info->beta_max*180.0/((double)M_PI));
      printf("polar-axis: %i\nx-semi: %e\ny-semi: %e\nz-semi: %e\n",
	     info->polar_axis, info->x_semi, info->y_semi, info->z_semi);
      break;

    case ORIENTATION:
      surface->orientation = !surface->orientation;
/* swap grid point and periodicity information */
      tmp = surface->r_points;
      surface->r_points = surface->s_points; surface->s_points = tmp;
      tmp = surface->r_period;
      surface->r_period = surface->s_period; surface->s_period = tmp;
      replot = 1;
      break;

    case CENTER:
      info->x_center = get_real( io_ptr, "X-center:", info->x_center, LEVEL+1);
      info->y_center = get_real( io_ptr, "Y-center:", info->y_center, LEVEL+1);
      info->z_center = get_real( io_ptr, "Z-center:", info->z_center, LEVEL+1);
/* recompute bounding box */
      surface_bb( surface );
      replot = 1;
      break;

/* change the plot mode */
    case PLOT_MODE:
      surface->plot_mode = surface_plot_mode(io_ptr, NULL, surface,
						 surface->plot_mode);
/* don't need to redraw the surface grids! */
      break;

    case R_LINES:
      surface->r_points = int_max(2, get_int( io_ptr, "Number of r-lines:", 
				 surface->r_points, NO_INDENT));
      replot = 1;
      break;

    case S_LINES:
      surface->s_points = int_max(2, get_int( io_ptr, "Number of s-lines:", 
				 surface->s_points, NO_INDENT));
      replot = 1;
      break;

    case EXIT:
      quit = 1;
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

    case ALPHA_MIN:
      replot = 1;
      info->alpha_min = get_real( io_ptr, "Starting alpha angle:", 
				 info->alpha_min*180.0/((double)M_PI), NO_INDENT) * 
				   ((double)M_PI)/180.0;
/* recompute bounding box */
      surface_bb( surface );
/* check periodicity */
      surface->r_period = (fabs(info->alpha_max - info->alpha_min) + 1.e-7 >= 
			 2.0*((double)M_PI))? 1 : 0;
      break;

    case ALPHA_MAX:
      replot = 1;
      info->alpha_max = get_real( io_ptr, "Ending alpha angle:", 
				 info->alpha_max*180.0/((double)M_PI), NO_INDENT) * 
				   ((double)M_PI)/180.0;
/* recompute bounding box */
      surface_bb( surface );
/* check periodicity */
      surface->r_period = (fabs(info->alpha_max - info->alpha_min) + 1.e-7 >= 
			 2.0*((double)M_PI))? 1 : 0;
      break;

    case BETA_MIN:
      replot = 1;
      info->beta_min = real_max( eps, get_real( io_ptr, "Starting beta angle > 0:", 
				info->beta_min*180.0/((double)M_PI), NO_INDENT) * 
				((double)M_PI)/180.0 );
      surface_bb( surface );
      break;

    case BETA_MAX:
      replot = 1;
      info->beta_max = real_min( ((double)M_PI)-eps, get_real( io_ptr, 
							      "Ending beta angle < 180:", 
				info->beta_max*180.0/((double)M_PI), NO_INDENT) * 
				((double)M_PI)/180.0 );
      surface_bb( surface );
      break;

    case X_SEMI:
      replot = 1;
      info->x_semi = real_max( eps, get_real( io_ptr, "Semi-axis in x-direction:", 
				 info->x_semi, NO_INDENT) );
      surface_bb( surface );
      break;

    case Y_SEMI:
      replot = 1;
      info->y_semi = real_max( eps, get_real( io_ptr, "Semi-axis in y-direction:", 
				 info->y_semi, NO_INDENT) );
      surface_bb( surface );
      break;

    case Z_SEMI:
      replot = 1;
      info->z_semi = real_max( eps, get_real( io_ptr, "Semi-axis in z-direction:", 
				 info->z_semi, NO_INDENT) );
      surface_bb( surface );
      break;

    case P_AXIS:
      replot = 1;
      info->polar_axis = get_int( io_ptr, "Polar axis:", 
				 info->polar_axis, NO_INDENT);
      info->polar_axis = int_min( 3, int_max(1,info->polar_axis));
      surface_bb( surface );
      break;

    default:
      break;

    }
  } while (!quit);
}

int
sphere_patch( grid_point *gp_ptr, surface_mapping *surface ){
  sphere_patch_data *info;
  real r, s, alpha, beta, alpha_r, beta_s, x1, x1_r, x1_s, x2, x2_r, x2_s
    , x3, x3_r, x3_s, x1_semi=1.0, x2_semi=1.0, x3_semi=1.0;

/* cast the data pointer to the right type */
  info = (sphere_patch_data *) surface->surface_data_ptr;

/* orient the semi-axes */
  if (info->polar_axis == 1){
    x1_semi = info->y_semi;
    x2_semi = info->z_semi;
    x3_semi = info->x_semi;
  }
  else if (info->polar_axis == 2){
    x1_semi = info->z_semi;
    x2_semi = info->x_semi;
    x3_semi = info->y_semi;
  }
  else if (info->polar_axis == 3){
    x1_semi = info->x_semi;
    x2_semi = info->y_semi;
    x3_semi = info->z_semi;
  }

/* swap r and s ? */
  if (surface->orientation){
    r = gp_ptr->s;
    s = gp_ptr->r;
  }
  else{
    r = gp_ptr->r;
    s = gp_ptr->s;
  }

/* angles */
  alpha = info->alpha_min + r * (info->alpha_max - info->alpha_min);
  beta = info->beta_min + s * (info->beta_max - info->beta_min);
  alpha_r = info->alpha_max - info->alpha_min;
  beta_s  = info->beta_max  - info->beta_min ;

/* function */
  x1 = x1_semi * cos( alpha ) * sin( beta );
  x2 = x2_semi * sin( alpha ) * sin( beta );
  x3 = x3_semi * cos( beta );

/* jacobian */
  if (surface->orientation){
    x1_s = -x1_semi * sin( alpha ) * sin( beta ) * alpha_r;
    x1_r =  x1_semi * cos( alpha ) * cos( beta ) * beta_s;
    x2_s =  x2_semi * cos( alpha ) * sin( beta ) * alpha_r;
    x2_r =  x2_semi * sin( alpha ) * cos( beta ) * beta_s;
    x3_s =  0.0;
    x3_r = -x3_semi * sin( beta ) * beta_s;
  }
  else{
    x1_r = -x1_semi * sin( alpha ) * sin( beta ) * alpha_r;
    x1_s =  x1_semi * cos( alpha ) * cos( beta ) * beta_s;
    x2_r =  x2_semi * cos( alpha ) * sin( beta ) * alpha_r;
    x2_s =  x2_semi * sin( alpha ) * cos( beta ) * beta_s;
    x3_r =  0.0;
    x3_s = -x3_semi * sin( beta ) * beta_s;
  }

  if (info->polar_axis == 1){
    gp_ptr->x = x3; gp_ptr->y = x1; gp_ptr->z = x2;
    gp_ptr->xr = x3_r; gp_ptr->yr = x1_r; gp_ptr->zr = x2_r;
    gp_ptr->xs = x3_s; gp_ptr->ys = x1_s; gp_ptr->zs = x2_s;
  }
  else if (info->polar_axis == 2){
    gp_ptr->x = x2; gp_ptr->y = x3; gp_ptr->z = x1;
    gp_ptr->xr = x2_r; gp_ptr->yr = x3_r; gp_ptr->zr = x1_r;
    gp_ptr->xs = x2_s; gp_ptr->ys = x3_s; gp_ptr->zs = x1_s;
  }
  else if (info->polar_axis == 3){
    gp_ptr->x = x1; gp_ptr->y = x2; gp_ptr->z = x3;
    gp_ptr->xr = x1_r; gp_ptr->yr = x2_r; gp_ptr->zr = x3_r;
    gp_ptr->xs = x1_s; gp_ptr->ys = x2_s; gp_ptr->zs = x3_s;
  }
  else{
    printf("Error in sphere_patch_mapping: Unknown value of polar_axis = %i\n",
	   info->polar_axis);
    return ERROR;
  }
/* translate */
  gp_ptr->x += info->x_center;
  gp_ptr->y += info->y_center;
  gp_ptr->z += info->z_center;

  return OK;
}

static int
sphere_patch_hessian(hessian_data *h_ptr, struct surface_mapping *surface){
  sphere_patch_data *info;
  real r, s, alpha, beta, alpha_r, beta_s, x1_rr, x1_rs, x1_ss, x2_rr, x2_rs, x2_ss
    , x3_rr, x3_rs, x3_ss, x1_semi=1.0, x2_semi=1.0, x3_semi=1.0;

/* cast the data pointer to the right type */
  info = (sphere_patch_data *) surface->surface_data_ptr;

/* orient the semi-axes */
  if (info->polar_axis == 1){
    x1_semi = info->y_semi;
    x2_semi = info->z_semi;
    x3_semi = info->x_semi;
  }
  else if (info->polar_axis == 2){
    x1_semi = info->z_semi;
    x2_semi = info->x_semi;
    x3_semi = info->y_semi;
  }
  else if (info->polar_axis == 3){
    x1_semi = info->x_semi;
    x2_semi = info->y_semi;
    x3_semi = info->z_semi;
  }

/* swap r and s ? */
  if (surface->orientation){
    r = h_ptr->s;
    s = h_ptr->r;
  }
  else{
    r = h_ptr->r;
    s = h_ptr->s;
  }

/* angles */
  alpha = info->alpha_min + r * (info->alpha_max - info->alpha_min);
  beta = info->beta_min + s * (info->beta_max - info->beta_min);
  alpha_r = info->alpha_max - info->alpha_min;
  beta_s  = info->beta_max  - info->beta_min ;

/* function */
/*   x1 = x1_semi * cos( alpha ) * sin( beta ); */
/*   x2 = x2_semi * sin( alpha ) * sin( beta ); */
/*   x3 = x3_semi * cos( beta ); */

/* jacobian */
  if (surface->orientation){
/*     x1_s = -x1_semi * sin( alpha ) * sin( beta ) * alpha_r; */
    x1_ss = -x1_semi * cos( alpha ) * sin( beta ) * sqr(alpha_r);
/*     x1_r =  x1_semi * cos( alpha ) * cos( beta ) * beta_s; */
    x1_rr = -x1_semi * cos( alpha ) * sin( beta ) * sqr(beta_s);
    x1_rs = -x1_semi * sin( alpha ) * cos( beta ) * beta_s * alpha_r;

/*     x2_s =  x2_semi * cos( alpha ) * sin( beta ) * alpha_r; */
    x2_ss =  -x2_semi * sin( alpha ) * sin( beta ) * sqr(alpha_r);
/*     x2_r =  x2_semi * sin( alpha ) * cos( beta ) * beta_s; */
    x2_rr = -x2_semi * sin( alpha ) * sin( beta ) * sqr(beta_s);
    x2_rs =  x2_semi * cos( alpha ) * cos( beta ) * beta_s*alpha_r;

    x3_ss =  0.0;
/*     x3_s =  0.0; */
/*     x3_r = -x3_semi * sin( beta ) * beta_s; */
    x3_rr = -x3_semi * cos( beta ) * sqr(beta_s);
    x3_rs = 0.0;
  }
  else{
    /*    x1_r = -x1_semi * sin( alpha ) * sin( beta ) * alpha_r; */
    x1_rr = -x1_semi * cos( alpha ) * sin( beta ) * sqr(alpha_r);
    /*    x1_s =  x1_semi * cos( alpha ) * cos( beta ) * beta_s; */
    x1_ss = -x1_semi * cos( alpha ) * sin( beta ) * sqr(beta_s);
    x1_rs = -x1_semi * sin( alpha ) * cos( beta ) * beta_s * alpha_r;

    /*    x2_r =  x2_semi * cos( alpha ) * sin( beta ) * alpha_r;*/
    x2_rr =  -x2_semi * sin( alpha ) * sin( beta ) * sqr(alpha_r);
    /*    x2_s =  x2_semi * sin( alpha ) * cos( beta ) * beta_s;*/
    x2_ss = -x2_semi * sin( alpha ) * sin( beta ) * sqr(beta_s);
    x2_rs =  x2_semi * cos( alpha ) * cos( beta ) * beta_s*alpha_r;

    x3_rr =  0.0;
    /*    x3_s = -x3_semi * sin( beta ) * beta_s;*/
    x3_ss = -x3_semi * cos( beta ) * sqr(beta_s);
    x3_rs = 0.0;
  }

  if (info->polar_axis == 1){
    h_ptr->xrr = x3_rr; 
    h_ptr->yrr = x1_rr; 
    h_ptr->zrr = x2_rr;

    h_ptr->xss = x3_ss; 
    h_ptr->yss = x1_ss; 
    h_ptr->zss = x2_ss;

    h_ptr->xrs = x3_rs; 
    h_ptr->yrs = x1_rs; 
    h_ptr->zrs = x2_rs;
  }
  else if (info->polar_axis == 2){
    h_ptr->xrr = x2_rr; 
    h_ptr->yrr = x3_rr; 
    h_ptr->zrr = x1_rr;

    h_ptr->xss = x2_ss; 
    h_ptr->yss = x3_ss; 
    h_ptr->zss = x1_ss;

    h_ptr->xrs = x2_rs; 
    h_ptr->yrs = x3_rs; 
    h_ptr->zrs = x1_rs;
  }
  else if (info->polar_axis == 3){
    h_ptr->xrr = x1_rr; 
    h_ptr->yrr = x2_rr; 
    h_ptr->zrr = x3_rr;

    h_ptr->xss = x1_ss; 
    h_ptr->yss = x2_ss; 
    h_ptr->zss = x3_ss;

    h_ptr->xrs = x1_rs; 
    h_ptr->yrs = x2_rs; 
    h_ptr->zrs = x3_rs;
  }
  else{
    printf("Error in sphere_patch_hessian: Unknown value of polar_axis = %i\n",
	   info->polar_axis);
    return ERROR;
  }
  return OK;
}
