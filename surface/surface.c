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

static int
numerical_hessian( hessian_data *h_ptr, surface_mapping *surface );

surface_mapping *
new_surface_mapping( char *name){
  surface_mapping *surface;
  static int n_calls=0;

  n_calls++;

  surface = (surface_mapping *) malloc(sizeof(surface_mapping));
  
/* initialize all fields */
  surface->used_by_volume = 0;

/* copy the name */
  surface->name = (char *) malloc( (strlen(name)+1)*sizeof(char) );
  surface->name = strcpy( surface->name, name );

  surface->inverse_known = 0;
  surface->forward_mapping = NULL;
  surface->hessian = numerical_hessian;
  surface->inverse_mapping = NULL;
  surface->set_surface = NULL;
  surface->delete_surface_data = NULL;
  surface->write_surface_data = NULL;
  
  surface->surface_data_ptr = NULL;

  surface->r_points = 0;
  surface->s_points = 0;
  surface->r_period = 0;
  surface->s_period = 0;

  surface->bb = (bounding_box *) malloc( sizeof(bounding_box) );
  surface->bb->x_min = 0.0;
  surface->bb->x_max = 1.0;
  surface->bb->y_min = 0.0;
  surface->bb->y_max = 1.0;
  surface->bb->z_min = 0.0;
  surface->bb->z_max = 1.0;

/* set the grid type and the periodicity flag */
  surface->type= NULL; /* unknown grid type */

/* set the default plot mode */
  surface->plot_mode = 1 | 4 | 8 | 256;

/* set the color */
  surface->color = n_calls % 8;

/* assume the orientation is clock-wise */
  surface->orientation = 0;

  return surface;
}

surface_mapping *
delete_surface_mapping( surface_mapping *surface ){
  if (surface == NULL) return NULL;

/* check if the mapping is used by a volume mapping*/
  if (surface->used_by_volume) return surface;

/* delete the specific data */
  delete_specific_surface( surface );

/* free the name */
  free( surface->name );

/* free the bounding box */
  free( surface->bb );

/* free the structure itself */
  free( surface );
  return NULL;
}

int
compute_surface_normal(real *n_vec, grid_point *gp_ptr){
  real length;
  int i;

  n_vec[0] = (gp_ptr->yr * gp_ptr->zs - gp_ptr->zr * gp_ptr->ys);
  n_vec[1] = (gp_ptr->zr * gp_ptr->xs - gp_ptr->xr * gp_ptr->zs);
  n_vec[2] = (gp_ptr->xr * gp_ptr->ys - gp_ptr->yr * gp_ptr->xs);

/* normalize */
  length = sqrt( n_vec[0] * n_vec[0] + n_vec[1] * n_vec[1] + n_vec[2] * n_vec[2] );
  if (length > NEWTON_EPS){
    for (i=0; i<3; i++) n_vec[i] = n_vec[i]/length;
    return OK;
  }
  else
    return ERROR;
}

static int
numerical_hessian( hessian_data *h_ptr, surface_mapping *surface ){
  grid_point gp_p, gp_m;
  const real eps = 1.e-5;

  /* form the derivatives numerically */
  gp_p.r = h_ptr->r + eps;
  gp_p.s = h_ptr->s;
  forward_surface_mapping(&gp_p, surface);

  gp_m.r = h_ptr->r - eps;
  gp_m.s = h_ptr->s;
  forward_surface_mapping(&gp_m, surface);

  h_ptr->xrr = 0.5*(gp_p.xr - gp_m.xr)/eps;
  h_ptr->yrr = 0.5*(gp_p.yr - gp_m.yr)/eps;
  h_ptr->zrr = 0.5*(gp_p.zr - gp_m.zr)/eps;

  h_ptr->xrs = 0.5*(gp_p.xs - gp_m.xs)/eps;
  h_ptr->yrs = 0.5*(gp_p.ys - gp_m.ys)/eps;
  h_ptr->zrs = 0.5*(gp_p.zs - gp_m.zs)/eps;

  gp_p.r = h_ptr->r;
  gp_p.s = h_ptr->s + eps;
  forward_surface_mapping(&gp_p, surface);

  gp_m.r = h_ptr->r;
  gp_m.s = h_ptr->s - eps;
  forward_surface_mapping(&gp_m, surface);

  h_ptr->xss = 0.5*(gp_p.xs - gp_m.xs)/eps;
  h_ptr->yss = 0.5*(gp_p.ys - gp_m.ys)/eps;
  h_ptr->zss = 0.5*(gp_p.zs - gp_m.zs)/eps;

  return OK;
}

int
derivative_surface_normal(real n_vec_r[3], real n_vec_s[3], grid_point *gp_ptr, 
			  surface_mapping *surface){
/* computes the derivative of the surface normal in the two parameter directions */
  real length, n_vec[3], dldr, dlds;
  hessian_data hs;
  int i;

/* get the first derivatives */
  forward_surface_mapping(gp_ptr, surface);

  n_vec[0] = (gp_ptr->yr * gp_ptr->zs - gp_ptr->zr * gp_ptr->ys);
  n_vec[1] = (gp_ptr->zr * gp_ptr->xs - gp_ptr->xr * gp_ptr->zs);
  n_vec[2] = (gp_ptr->xr * gp_ptr->ys - gp_ptr->yr * gp_ptr->xs);

/* normalize */
  length = sqrt( n_vec[0] * n_vec[0] + n_vec[1] * n_vec[1] + n_vec[2] * n_vec[2] );

  if (length < NEWTON_EPS){
    for (i=0; i<3; i++){
      n_vec_r[i] = n_vec_s[i] = 0.0;
    }
    return ERROR;
  }
/* get the second derivatives */
  hs.r = gp_ptr->r; hs.s = gp_ptr->s;
  surface_hessian( &hs, surface);

/* compute the derivative of the normal */
  n_vec_r[0] = (hs.yrr * gp_ptr->zs + gp_ptr->yr * hs.zrs - 
		hs.zrr * gp_ptr->ys - gp_ptr->zr * hs.yrs);
  n_vec_r[1] = (hs.zrr * gp_ptr->xs + gp_ptr->zr * hs.xrs -
		hs.xrr * gp_ptr->zs - gp_ptr->xr * hs.zrs);
  n_vec_r[2] = (hs.xrr * gp_ptr->ys + gp_ptr->xr * hs.yrs - 
		hs.yrr * gp_ptr->xs - gp_ptr->yr * hs.xrs);

  n_vec_s[0] = (hs.yrs * gp_ptr->zs + gp_ptr->yr * hs.zss - 
		hs.zrs * gp_ptr->ys - gp_ptr->zr * hs.yss);
  n_vec_s[1] = (hs.zrs * gp_ptr->xs + gp_ptr->zr * hs.xss -
		hs.xrs * gp_ptr->zs - gp_ptr->xr * hs.zss);
  n_vec_s[2] = (hs.xrs * gp_ptr->ys + gp_ptr->xr * hs.yss - 
		hs.yrs * gp_ptr->xs - gp_ptr->yr * hs.xss);

  dldr = (n_vec[0]*n_vec_r[0] + n_vec[1]*n_vec_r[1] + n_vec[2]*n_vec_r[2])/length;
  dlds = (n_vec[0]*n_vec_s[0] + n_vec[1]*n_vec_s[1] + n_vec[2]*n_vec_s[2])/length;

  for (i=0; i<3; i++)
    {
      n_vec_r[i] = n_vec_r[i]/length - n_vec[i]*dldr/sqr(length);
      n_vec_s[i] = n_vec_s[i]/length - n_vec[i]*dlds/sqr(length);
    }
  return OK;
}


int 
forward_surface_mapping( grid_point *gp_ptr, surface_mapping *surface ){
/* take care of any errors in the specific mapping routine and pass them on */
  return (*surface->forward_mapping)( gp_ptr, surface );
}

int 
inverse_surface_mapping( grid_point *gp_ptr, surface_mapping *surface ){
/* take care of any errors in the specific mapping routine and pass them on */
  return (*surface->inverse_mapping)( gp_ptr, surface );
}

int
surface_hessian( hessian_data *h_ptr, surface_mapping *surface ){
/* take care of any errors in the specific mapping routine and pass them on */
  return (*surface->hessian)( h_ptr, surface );
}

void 
set_surface( input_output *io_ptr, surface_mapping *surface ){
  if (surface == NULL || surface->set_surface == NULL){
    printf("ERROR: set_surface was called with inconsistent pointers.\n");
  }
  else
    (*surface->set_surface)( io_ptr, surface );
}

void 
delete_specific_surface( surface_mapping *surface ){
  if (surface == NULL || surface->delete_surface_data == NULL){
    printf("ERROR: delete_specific_surface was called with inconsistent pointers.\n");
  }
  else
    (*surface->delete_surface_data)( surface->surface_data_ptr );
}

void 
write_specific_surface( int fd, surface_mapping *surface ){
  if (surface == NULL || surface->write_surface_data == NULL){
    printf("ERROR: write_specific_surface was called with inconsistent pointers.\n");
  }
  else
    (*surface->write_surface_data)( fd, surface->surface_data_ptr );
}


surface_mapping *
choose_surface( input_output *io_ptr, char *name ){
  char prompt[80];
  int icom, quit;
  surface_mapping *surface=NULL;

#include "choose_surface_com.h"

  sprintf(prompt,"set surface mapping type %s>", name);

  quit = 0;
  do{

    icom = get_command(io_ptr, prompt, COMMAND, LAST_COM, LEVEL,
		       SAVE_ON_COPY, ARGUMENT);

    switch (icom) {

/* surface grids */
    case SPHERE: 
/* make surface */
      surface = new_surface_mapping(name);
/* patch on a sphere */
      init_sphere_patch( surface );
      quit = 1;
      break;

    case DISCRETE: 
/* make surface */
      surface = new_surface_mapping(name);
/* read the discrete gridpoint from a file */
      if (init_disc_point( io_ptr, surface ) != NULL){
	quit = 1;
      }
      else{
	surface = delete_surface_mapping( surface );
      }
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

    case CANCEL:
/* cancel, unsucessful surface grid generation */
      return NULL;
      break;

    default:
      break;
    }
  }
  while(!quit);

/* set the view point and focal points to make the new object visible */
  ogl_standard_view( SURFACE_W, surface->bb );

/* call the specific setup routine to enter all the curve specific data */
  set_surface( io_ptr, surface ); 

/* sucessful curve generation */
  return surface;
}

surface_mapping_link *
get_surface( input_output *io_ptr, surface_mapping_list *head, int level){
  surface_mapping *surface;
  surface_mapping_link *this_link;
  char **command;
  int i, icom, ncom, *save_on_copy=NULL, *argument=NULL, quit;

  if (head->n_members == 0){
    printf("The list of surface mappings is empty.\n");
    return NULL;
  }

  ncom = head->n_members + 1;
  command = (char **) malloc( (ncom+1)*sizeof(char*) );

  i=0;
/* add color information to the name! */
  for (this_link = head->first; this_link != NULL; 
       this_link = this_link->next){
    surface  = this_link->data;
    command[i] = (char *) malloc( (strlen(surface->name) + 
				   strlen(ogl_color_name(surface->color)) + 3) 
				 * sizeof(char) );
    sprintf( command[i], "%s(%s)", surface->name, ogl_color_name(surface->color) );
    i++;
  }
  command[ncom-1] = "help";
  command[ncom]   = "cancel";

  quit = 0;
  do{
    icom = get_command(io_ptr, "surface mapping name: ", command, ncom, level, 
		       save_on_copy, argument);

/* default is to quit */
    quit = 1;

    if (icom == -1)
/* not a valid command */
      quit = 0;
    else if (icom == ncom-1){
/* help */
      printf(
" o  Enter one of the surface names if you wish to proceed.\n"
" o  Enter cancel if you do not wish to proceed. You will then return to the\n"
"    previous command level.\n");
      quit = 0;
    }

  }
  while (!quit);

/* return the pointer to the selected curve */
  i=0;
  for (this_link = head->first; this_link != NULL; 
       this_link = this_link->next){
    if (i==icom) break;
    i++;
  }

/* free the strings */
  for (i=0; i<= ncom-2; i++)
    free(command[i]);

/* free the pointers to the strings */
  free(command);

  return this_link;

}

void
surface_coordinates(surface_mapping *surface, real_array_2d **x_, 
		    real_array_2d **y_, real_array_2d **z_)
{
  int i, j, i_min, j_min, i_max, j_max, status=ERROR;
  real r_step, s_step;
  grid_point gp;

#define x_in(i, j) compute_index_2d((*x_), i, j)
#define y_in(i, j) compute_index_2d((*y_), i, j)
#define z_in(i, j) compute_index_2d((*z_), i, j)

/*   if (surface->r_period) */
/*     i_min = 2; */
/*   else */
  i_min = 1;
  i_max = surface->r_points;

/*   if (surface->s_period) */
/*     j_min = 2; */
/*   else */
  j_min = 1;
  j_max = surface->s_points;

  r_step = 1.0/(surface->r_points - 1);
  s_step = 1.0/(surface->s_points - 1);

/* define the multi-dimensional arrays */
  *x_ = create_real_array_2d(i_max, j_max);
  *y_ = create_real_array_2d(i_max, j_max);
  *z_ = create_real_array_2d(i_max, j_max);

/* assign x, y, z. */
  for (i = 1; i <= i_max; i++)
    for (j = 1; j <= j_max; j++)
      {
	gp.r = (i-i_min) * r_step;
	gp.s = (j-j_min) * s_step;
	status = forward_surface_mapping(&gp, surface); 
	x_in(i,j) = gp.x;
	y_in(i,j) = gp.y;
	z_in(i,j) = gp.z;
      }
  /* error check */
  if (status != OK)
    printf("define_grid: error in grid `%s'.\n",surface->name);
}


void
surface_list_bb(bounding_box *bb, surface_mapping_list *surface_list)
{
  surface_mapping_link *surface_link;
  surface_mapping *surface;

  if (!surface_list  || surface_list->n_members == 0)
    {
      bb->x_min = 0.0;
      bb->x_max = 1.0;
      bb->y_min = 0.0;
      bb->y_max = 1.0;
      bb->z_min = 0.0;
      bb->z_max = 1.0;
    }
  else
    {
      bb->x_min = 1.0e10;
      bb->x_max =-1.0e10;
      bb->y_min = 1.0e10;
      bb->y_max =-1.0e10;
      bb->z_min = 1.0e10;
      bb->z_max =-1.0e10;
      for (surface_link = surface_list->first; surface_link != NULL; 
	   surface_link = surface_link->next)
	{
	  surface = surface_link->data;
	  bb->x_min = real_min(bb->x_min, surface->bb->x_min);
	  bb->x_max = real_max(bb->x_max, surface->bb->x_max);
	  bb->y_min = real_min(bb->y_min, surface->bb->y_min);
	  bb->y_max = real_max(bb->y_max, surface->bb->y_max);
	  bb->z_min = real_min(bb->z_min, surface->bb->z_min);
	  bb->z_max = real_max(bb->z_max, surface->bb->z_max);
	}
    } /* end if surface_list */
}

void
check_surface(input_output *io_ptr, surface_mapping *surface){
  grid_point gp, gp_p, gp_m;
  hessian_data hs;
  real xr[3], xs[3], xrr[3], xrs[3], xss[3], dr, ds, error, max_err, 
    r_max, s_max, h_error, max_h_err, r_h_max, s_h_max, eps, n_vec_r[3], n_vec_s[3];
  int NO_INDENT = 0, i, j;

  eps = real_max( 1.e-10, get_real(io_ptr, "Enter eps for numerical "
				   "differentiation (>0) :", 1.e-4, NO_INDENT));

  dr = 1.0/(surface->r_points-1);
  ds = 1.0/(surface->s_points-1);

  max_h_err = 0.0;
  r_h_max = -1; s_h_max = -1;
  max_err = 0.0;
  r_max = -1; s_max = -1;

  for (i=1; i<=surface->r_points; i++)
    for (j=1; j<=surface->s_points; j++)
	{
	  gp.r = (i-1)*dr; gp.s = (j-1)*ds; 

/* evaluate the analytic jacobian */
	  if (!forward_surface_mapping(&gp, surface))
	    printf("forward_surface_mapping returned an error at parameter (%e, %e)\n",
		   gp.r, gp.s);

/* estimate the jacobian numerically */
	  gp_p.r = gp.r+eps;
	  gp_p.s = gp.s;
	  forward_surface_mapping(&gp_p, surface);
	  gp_m.r = gp.r-eps;
	  gp_m.s = gp.s;
	  forward_surface_mapping(&gp_m, surface);
	  xr[0] = 0.5*(gp_p.x - gp_m.x)/eps;
	  xr[1] = 0.5*(gp_p.y - gp_m.y)/eps;
	  xr[2] = 0.5*(gp_p.z - gp_m.z)/eps;

	  xrr[0] = 0.5*(gp_p.xr - gp_m.xr)/eps;
	  xrr[1] = 0.5*(gp_p.yr - gp_m.yr)/eps;
	  xrr[2] = 0.5*(gp_p.zr - gp_m.zr)/eps;

	  xrs[0] = 0.5*(gp_p.xs - gp_m.xs)/eps;
	  xrs[1] = 0.5*(gp_p.ys - gp_m.ys)/eps;
	  xrs[2] = 0.5*(gp_p.zs - gp_m.zs)/eps;

	  gp_p.r = gp.r;
	  gp_p.s = gp.s+eps;
	  forward_surface_mapping(&gp_p, surface);
	  gp_m.r = gp.r;
	  gp_m.s = gp.s-eps;
	  forward_surface_mapping(&gp_m, surface);
	  xs[0] = 0.5*(gp_p.x - gp_m.x)/eps;
	  xs[1] = 0.5*(gp_p.y - gp_m.y)/eps;
	  xs[2] = 0.5*(gp_p.z - gp_m.z)/eps;

	  xss[0] = 0.5*(gp_p.xs - gp_m.xs)/eps;
	  xss[1] = 0.5*(gp_p.ys - gp_m.ys)/eps;
	  xss[2] = 0.5*(gp_p.zs - gp_m.zs)/eps;

	  error = sqrt( sqr(gp.xr-xr[0]) + sqr(gp.yr-xr[1]) + sqr(gp.zr-xr[2]) + 
			sqr(gp.xs-xs[0]) + sqr(gp.ys-xs[1]) + sqr(gp.zs-xs[2]));
	  if (error > max_err)
	    {
	      max_err = error;
	      r_max = gp.r;
	      s_max = gp.s;
	    }

/* evaluate the analytic hessian */
	  hs.r = gp.r; hs.s = gp.s;
	  surface_hessian( &hs, surface );
	  h_error = sqrt(sqr(hs.xrr-xrr[0]) + sqr(hs.yrr-xrr[1]) + sqr(hs.zrr-xrr[2]) + 
			 sqr(hs.xrs-xrs[0]) + sqr(hs.yrs-xrs[1]) + sqr(hs.zrs-xrs[2]) + 
			 sqr(hs.xss-xss[0]) + sqr(hs.yss-xss[1]) + sqr(hs.zss-xss[2]));
	  if (h_error > max_h_err)
	    {
	      max_h_err = h_error;
	      r_h_max = gp.r;
	      s_h_max = gp.s;
	    }
	} /* end for all grid points */
  printf("The largest error in the Jacobian was %e and occured at (%e, %e)\n", 
	 max_err, r_max, s_max);

  printf("The largest error in the Hessian was %e and occured at (%e, %e)\n", 
	 max_h_err, r_h_max, s_h_max);

  printf("\nChecking the derivative of the surface normal\n");
  for (i=1; i<=surface->r_points; i++)
    for (j=1; j<=surface->s_points; j++)
	{
	  gp.r = (i-1)*dr; gp.s = (j-1)*ds; 

/* evaluate the derivative of the surface normal */
	  if (!derivative_surface_normal(n_vec_r, n_vec_s, &gp, surface))
	    printf("derivative_surface_normal returned an error at parameter "
		   "(%e, %e)\n", gp.r, gp.s);
	}

/*   gp.r = r_max; */
/*   gp.s = s_max; */
/*   while( (gp.r = get_real( io_ptr, "r ( exit with r <= -1 ) :", gp.r, NO_INDENT )) > -1) */
/*     { */
/*       gp.s = get_real( io_ptr, "s:", gp.s, NO_INDENT ); */
/*       forward_surface_mapping(&gp, surface); */

/* estimate the jacobian numerically */
/*       gp_p.r = gp.r+eps; */
/*       gp_p.s = gp.s; */
/*       forward_surface_mapping(&gp_p, surface); */
/*       gp_m.r = gp.r-eps; */
/*       gp_m.s = gp.s; */
/*       forward_surface_mapping(&gp_m, surface); */
/*       xr[0] = 0.5*(gp_p.x - gp_m.x)/eps; */
/*       xr[1] = 0.5*(gp_p.y - gp_m.y)/eps; */
/*       xr[2] = 0.5*(gp_p.z - gp_m.z)/eps; */

/*       xrr[0] = 0.5*(gp_p.xr - gp_m.xr)/eps; */
/*       xrr[1] = 0.5*(gp_p.yr - gp_m.yr)/eps; */
/*       xrr[2] = 0.5*(gp_p.zr - gp_m.zr)/eps; */

/*       xrs[0] = 0.5*(gp_p.xs - gp_m.xs)/eps; */
/*       xrs[1] = 0.5*(gp_p.ys - gp_m.ys)/eps; */
/*       xrs[2] = 0.5*(gp_p.zs - gp_m.zs)/eps; */

/*       gp_p.r = gp.r; */
/*       gp_p.s = gp.s+eps; */
/*       forward_surface_mapping(&gp_p, surface); */
/*       gp_m.r = gp.r; */
/*       gp_m.s = gp.s-eps; */
/*       forward_surface_mapping(&gp_m, surface); */
/*       xs[0] = 0.5*(gp_p.x - gp_m.x)/eps; */
/*       xs[1] = 0.5*(gp_p.y - gp_m.y)/eps; */
/*       xs[2] = 0.5*(gp_p.z - gp_m.z)/eps; */

/*       xss[0] = 0.5*(gp_p.xs - gp_m.xs)/eps; */
/*       xss[1] = 0.5*(gp_p.ys - gp_m.ys)/eps; */
/*       xss[2] = 0.5*(gp_p.zs - gp_m.zs)/eps; */

/* evaluate the analytic hessian */
/*       hs.r = gp.r; hs.s = gp.s; */
/*       surface_hessian( &hs, surface ); */

/* print the result */
/*       printf("Analytic jacobian:\n" */
/* 	     "xr = ( %e, %e, %e )\nxs = ( %e, %e, %e )\n\n",  */
/* 	     gp.xr, gp.yr, gp.zr, gp.xs, gp.ys, gp.zs); */
/*       printf("Difference between analytic and numerical jacobian:\n" */
/* 	     "xr = ( %e, %e, %e )\nxs = ( %e, %e, %e )\n\n",  */
/* 	     gp.xr-xr[0], gp.yr-xr[1], gp.zr-xr[2],  */
/* 	     gp.xs-xs[0], gp.ys-xs[1], gp.zs-xs[2]); */
/*       printf("Analytic hessian:\n" */
/* 	     "xrr = ( %e, %e, %e )\nxrs = ( %e, %e, %e )\nxss = ( %e, %e, %e )\n\n",  */
/* 	     hs.xrr, hs.yrr, hs.zrr, hs.xrs, hs.yrs, hs.zrs, hs.xss, hs.yss, hs.zss); */
/*       printf("Discrete hessian:\n" */
/* 	     "xrr = ( %e, %e, %e )\nxrs = ( %e, %e, %e )\nxss = ( %e, %e, %e )\n\n",  */
/* 	     xrr[0], xrr[1], xrr[2], xrs[0], xrs[1], xrs[2], xss[0], xss[1], xss[2]); */
/*       printf("Difference between analytic and numerical hessian:\n" */
/* 	     "xrr = ( %e, %e, %e )\nxrs = ( %e, %e, %e )\nxss = ( %e, %e, %e )\n\n",  */
/* 	     hs.xrr-xrr[0], hs.yrr-xrr[1], hs.zrr-xrr[2],  */
/* 	     hs.xrs-xrs[0], hs.yrs-xrs[1], hs.zrs-xrs[2],  */
/* 	     hs.xss-xss[0], hs.yss-xss[1], hs.zss-xss[2]); */
	     
/*     }     */

}

static void
update_bb( bounding_box *bb, grid_point *gp ){
  bb->x_min = real_min( bb->x_min, gp->x );
  bb->x_max = real_max( bb->x_max, gp->x );
  bb->y_min = real_min( bb->y_min, gp->y );
  bb->y_max = real_max( bb->y_max, gp->y );
  bb->z_min = real_min( bb->z_min, gp->z );
  bb->z_max = real_max( bb->z_max, gp->z );
}

void 
surface_bb( surface_mapping *surface ){
  grid_point gp;
  const real step=0.05;

/* reset bounding box */
  surface->bb->x_min = 1.e10;
  surface->bb->x_max = -1.e10;
  surface->bb->y_min = 1.e10;
  surface->bb->y_max = -1.e10;
  surface->bb->z_min = 1.e10;
  surface->bb->z_max = -1.e10;

  gp.s = 0.0;
  for (gp.r=0.0; gp.r<=1.0+0.5*step; gp.r+=step){
    forward_surface_mapping( &gp, surface );
    update_bb( surface->bb, &gp );
  }

  gp.s = 0.5;
  for (gp.r=0.0; gp.r<=1.0+0.5*step; gp.r+=step){
    forward_surface_mapping( &gp, surface );
    update_bb( surface->bb, &gp );
  }

  gp.s = 1.0;
  for (gp.r=0.0; gp.r<=1.0+0.5*step; gp.r+=step){
    forward_surface_mapping( &gp, surface );
    update_bb( surface->bb, &gp );
  }

  gp.r = 0.0;
  for (gp.s=0.0; gp.s<=1.0+0.5*step; gp.s+=step){
    forward_surface_mapping( &gp, surface );
    update_bb( surface->bb, &gp );
  }

  gp.r = 0.5;
  for (gp.s=0.0; gp.s<=1.0+0.5*step; gp.s+=step){
    forward_surface_mapping( &gp, surface );
    update_bb( surface->bb, &gp );
  }

  gp.r = 1.0;
  for (gp.s=0.0; gp.s<=1.0+0.5*step; gp.s+=step){
    forward_surface_mapping( &gp, surface );
    update_bb( surface->bb, &gp );
  }

}

