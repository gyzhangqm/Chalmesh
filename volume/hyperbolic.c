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

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

typedef struct{
  int n1, n2;
  real_array_2d *x_, *y_, *z_;
} grid_function;

typedef struct{
  real lambda_1r, lambda_1i, lambda_2r, lambda_2i;
  grid_function *x1, *x2, *Ax1, *Ax2, *u_eps, *u_eps_t;
} eigen_val;

/* private member functions */
static int
grid_normal(real *n_vec, grid_function *u, int i, int j, hyp_grid_data *info);
static void
set_hyp_bc( input_output *io_ptr, hyp_grid_data *info );
static void
compute_projection( input_output *io_ptr, volume_mapping *volume, surface_mapping *project );
static void
power_method_1( grid_function *u, grid_function *u_t, real t, hyp_grid_data *info,
	     eigen_val *eval );
static int
mean_curvature( grid_function *u, real_array_2d *kappa_, real *max_curvature );
static int 
hyp_grid_mapping( grid_point *gp, void *data_ptr);
static grid_function *
new_grid_function(int n1, int n2);
static grid_function *
delete_grid_function(grid_function *u);
static int
time_derivative(grid_function *u, real t, grid_function *u_t, hyp_grid_data *info);
static void 
compute_hyp_grid(input_output *io_ptr, volume_mapping *volume);
static void 
set_hyp_grid(input_output *io_ptr, volume_mapping *volume,
		   linked_list *surface_list);
static void
delete_hyp_grid( void *volume_data_ptr );
/* mappings for smoothed polygon grids */

static void 
set_hyp_grid(input_output *io_ptr, volume_mapping *volume,
	     linked_list *surface_list){
  hyp_grid_data *info;
  int quit=0, replot=TRUE, icom;
  const real eps=1.e-4;
  surface_mapping *project;
  linked_list_member *project_link;
  char *bc_names[4];
#include "set_hyp_grid_com.h"

  bc_names[0] = "free";
  bc_names[1] = "x-constant";
  bc_names[2] = "y-constant";
  bc_names[3] = "z-constant";

/* cast the data pointer to the right type */
  info = (hyp_grid_data *) volume->volume_data_ptr;

  do{

    if (replot){
/* plot the surface grid */
      if (ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
			 ogl_length_scale(volume->bb))){
/* if info->status != 1, only draw the surface */
	if (info->status == 1){
	  draw_volume( volume, volume->plot_mode, NULL );
	  if (volume->plot_mode & 256) 
	    ogl_bb_cage( volume->bb, OGL_WHITE );
	}
	else{
	  draw_surface( info->surface, info->surface->plot_mode );
	  if (info->surface->plot_mode & 256) 
	    ogl_bb_cage( info->surface->bb, OGL_WHITE );
	}
	ogl_end_plot();
	replot = FALSE;
      }
    }

    switch( get_command( io_ptr, "hyperbolic mapping>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case COMPUTE:
      if (info->status == 1){
	printf("The mapping is up to date with the parameters, \n"
	       "so there is no need to recompute it!\n");
      }
      else{
	compute_hyp_grid( io_ptr, volume );
	if (info->status == 1)
	  volume_bb( volume ); 
	replot = TRUE;
      }
      break;

    case CURV_FACTOR:
      info->curv_factor = real_max( 0.0, get_real(io_ptr, "curvature factor:", 
						  info->curv_factor, NO_INDENT));
      info->status = 0;
      break;

/*     case TIME_FACTOR: */
/*       info->time_factor = real_max( 0.0, get_real(io_ptr, "time step factor:",  */
/* 						  info->time_factor, NO_INDENT)); */
/*       info->status = 0; */
/*       break; */

    case SMOOTH_FACTOR:
      info->smooth_factor = 
	real_min( 1.0, real_max( 0.0, get_real(io_ptr, "averaging factor (>=0, <=1):", 
					       info->smooth_factor, NO_INDENT)) );
      info->status = 0;
      break;

    case V_MIN:
      info->v_min = real_max( 0.0, get_real(io_ptr, "smallest normal velocity:", 
					    info->v_min, NO_INDENT));
      info->status = 0;
      break;

    case SHOW:
      printf("This mapping is based on the surface mapping `%s'\n", 
	     info->surface->name);
      printf("averaging coefficient: %e\n", info->smooth_factor);
      printf("curvature coefficient: %e\n", info->curv_factor);
      printf("velocity threshold: %e\n", info->v_min);
      printf("thickness: %e\n", info->width);
      printf("Boundary condition at r1=0: %s, r1=1: %s, r2=0: %s, r2=1: %s.\n",
	     bc_names[hyp_bc(1,1)], bc_names[hyp_bc(2,1)], bc_names[hyp_bc(1,2)], 
	     bc_names[hyp_bc(2,2)]);
      printf("r1-points: %i\nr2-points: %i\nr3-points: %i\n",
	     volume->r1_points, volume->r2_points, volume->r3_points);
      printf("\n");
      break;

    case BOUNDARY_CONDITION:
      set_hyp_bc( io_ptr, info );
/* recompute the mapping */
      info->status = 0;
      break;

    case WIDTH:
      info->width = real_max( eps, get_real( io_ptr, "thickness:", 
					    info->width, NO_INDENT));
/* recompute recompute the mapping */
      info->status = 0;
      break;

    case R1_LINES:
      volume->r1_points = int_max(2, get_int(io_ptr, "Number of r1-lines:", 
					     volume->r1_points, NO_INDENT));
      info->status = 0;
      break;

    case R2_LINES:
      volume->r2_points = int_max(2, get_int(io_ptr, "Number of r2-lines:", 
					     volume->r2_points, NO_INDENT));
      info->status = 0;
      break;

    case R3_LINES:
      volume->r3_points = int_max(2, get_int(io_ptr, "Number of r3-lines:", 
					     volume->r3_points, NO_INDENT));
      info->status = 0;
      break;

    case PROJECT_FACE:
      if (info->status == 1){
	printf("Identify the surface to project onto:\n");
	project_link = get_surface( io_ptr, surface_list, NO_INDENT );
/* compute the projection */  
	if (project_link){
	  project = project_link->data;
	  compute_projection(io_ptr, volume, project);
	  replot = TRUE;
	}
      }
      else{
	printf("You must update the hyperbolic grid before you can project "
	       "the grid points. This is done with the command `compute-mapping'.\n");
      }
      break;

    case STRETCHING:
/* normal stretching */
      if (info->normal_stretch == NULL){
	info->normal_stretch = choose_stretching(io_ptr, &volume->r3_points );
	volume->plot_mode |= 8;
	info->status = 0;
      }
      else{
	if (set_stretching( io_ptr, info->normal_stretch, &volume->r3_points ) ){
	  info->status = 0;
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
	info->status = 0;
      }
      else{
	printf("There is no stretching in the r3-direction.\n");
      }
      break;

/* change the plot mode */
    case PLOT_MODE:
      if (info->status == 1){
	volume->plot_mode = volume_plot_mode(io_ptr, NULL, volume, volume->plot_mode, 
					     NULL); 
/* don't need to redraw the surface grids! */
      }
      else{
	printf("Sorry, you must compute the hyperbolic grid before you can look at it!\n"
	       "This is done with the command `compute-mapping'.\n");
      }
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL, NULL, NULL)) 
	     == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
/* make sure that we exit with a valid grid */
      if (info->status == 1){
	quit = 1;
      }
      else{
	printf(
"Sorry, you can only exit after the mapping has been successfully computed.\n"
"If you have trouble making the grid, I suggest that you make the thickness\n"
"very small, increase the averaging coefficient and decrease the curvature\n"
"coefficient.\n");
      }
      break;

    default:
      break;

    }
  } while (!quit);
}

static void
set_hyp_bc( input_output *io_ptr, hyp_grid_data *info ){
  char *bc_names[4];
  int icom, quit=0;

#include "set_hyp_bc_com.h"

  bc_names[0] = "free";
  bc_names[1] = "x-constant";
  bc_names[2] = "y-constant";
  bc_names[3] = "z-constant";

  do{

    switch (get_command( io_ptr, "boundary condition>", COMMAND, LAST_COM, LEVEL+1, 
			SAVE_ON_COPY, ARGUMENT)){

    case LOW_R1:
      hyp_bc(1,1) = get_command( io_ptr, "Low r1 bc>", bc_names, 3, NO_INDENT, 
				NULL, NULL);
      break;

    case HIGH_R1: 
      hyp_bc(2,1) = get_command( io_ptr, "High r1 bc>", bc_names, 3, NO_INDENT, 
				NULL, NULL);
      break;

    case LOW_R2: 
      hyp_bc(1,2) = get_command( io_ptr, "Low r2 bc>", bc_names, 3, NO_INDENT, 
				NULL, NULL);
      break;

    case HIGH_R2:
      hyp_bc(2,2) = get_command( io_ptr, "High r2 bc>", bc_names, 3, NO_INDENT, 
				NULL, NULL);
      break;

    case SHOW:
      printf("Boundary condition at r1=0: %s, r1=1: %s, r2=0: %s, r2=1: %s.\n",
	     bc_names[hyp_bc(1,1)], bc_names[hyp_bc(2,1)], bc_names[hyp_bc(1,2)], 
	     bc_names[hyp_bc(2,2)]);
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
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
  }
  while(!quit);
}


static void
compute_projection( input_output *io_ptr, volume_mapping *volume, surface_mapping *project ){

  char *prompt;
  int icom, quit=0, side = 0, dir = 0, n1, n2, n3, i, j, k;
  discrete_surface *sgrid_ptr;
  real_array_2d *x_gap_, *y_gap_, *z_gap_, *x_project, *y_project, *z_project;
  hyp_grid_data *info;
  inverse_point *interp_ptr;
  grid_point gp;
  real s, w, r_loc, s_loc, n_loc, l_scale;
  const real alpha=18.42; /* a pretty universial number */

#define x_gap(i,j) compute_index_2d(x_gap_,i,j)
#define y_gap(i,j) compute_index_2d(y_gap_,i,j)
#define z_gap(i,j) compute_index_2d(z_gap_,i,j)

#include "compute_projection_com.h"

  prompt = "Project face>";

  do{

    switch( get_command( io_ptr, prompt, COMMAND, LAST_COM, LEVEL, 
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

    case LEFT:
      side = 1; dir = 1;
      quit = 1;
      break;

    case RIGHT:
      side = 2; dir = 1;
      quit = 1;
      break;

    case LOWER:
      side = 1; dir = 2;
      quit = 1;
      break;

    case UPPER:
      side = 2; dir = 2;
      quit = 1;
      break;

/*     case NEAR: */
/*       printf("Sorry, projection of the near (r3=0) face is not implemented yet\n"); */
/*       side = 1; dir = 3; */
/*       break; */

/*     case MY_FAR: */
/*       printf("Sorry, projection of the upper (r3=1) face is not implemented yet\n"); */
/*       side = 2; dir = 3; */
/*       break; */

    case CANCEL:
      return;
      break;

    default:
      break;

    }
  } while (!quit);

/* do the projection for the selected face */

/* construct a discrete surface to approximately invert the surface */
/* mapping `project' */
  surface_coordinates(project, &x_project, &y_project, &z_project);
  sgrid_ptr = new_discrete_surface(project->name, x_project, y_project, z_project, 0, NULL);

/* make the tolerance independent of the mesh size in the `project' patch */
  l_scale = sqrt(sqr(sgrid_ptr->bb->x_max - sgrid_ptr->bb->x_min) +
		 sqr(sgrid_ptr->bb->y_max - sgrid_ptr->bb->y_min) +
		 sqr(sgrid_ptr->bb->z_max - sgrid_ptr->bb->z_min));
  sgrid_ptr->max_n_dist = 0.1 * l_scale;
/* tmp */
/*  printf("Mismatch tolerance for projection surface: %e\n", sgrid_ptr->max_n_dist);*/
/* end tmp */

/* assumes that the mapping is of type hyperbolic */
  info = (hyp_grid_data *) volume->volume_data_ptr;

/* the correspondence between the directions in the mother grid and the component */
/* surface grid are as follows */
/* dir  r_dim   s_dim  */
/*  1   r2_dim  r3_dim */
/*  2   r3_dim  r1_dim */
/*  3   r1_dim  r2_dim */

/* left or right face */
  if (dir == 1){
    n2 = info->n2;
    n3 = info->n3;
/* pick side */
    if (side == 1)
      i = 1;
    else
      i = info->n1;

/* allocate space for the projected surface */
    x_gap_ = create_real_array_2d( n2, n3 );
    y_gap_ = create_real_array_2d( n2, n3 );
    z_gap_ = create_real_array_2d( n2, n3 );

/* do the projection */
    printf("Projecting...\n");
    for (j=1; j<=n2; j++){
      k = 1;
      if ((interp_ptr = search_quad_tree(x3_vol(i,j,k), y3_vol(i,j,k), z3_vol(i,j,k), 
					 sgrid_ptr->quad_tree, sgrid_ptr, 0.0, 0, 0))){
/* save the solution to use as initial guess for next iteration */
	r_loc = interp_ptr->r_loc;
	s_loc = interp_ptr->s_loc;
	n_loc = interp_ptr->n_loc;
/* free the space taken by the interpolation point */
	free(interp_ptr);

/* evaluate the discrete surface */
	gp.r = r_loc; gp.s = s_loc;
	bi_linear(&gp, sgrid_ptr);
/* save the gap */
	x_gap(j,k) = gp.x - x3_vol(i,j,k);
	y_gap(j,k) = gp.y - y3_vol(i,j,k);
	z_gap(j,k) = gp.z - z3_vol(i,j,k);

	for (k=2; k<=n3; k++){
	  newton_search(x3_vol(i,j,k), y3_vol(i,j,k), z3_vol(i,j,k), 
			&r_loc, &s_loc, &n_loc, sgrid_ptr);
/* evaluate the discrete surface */
	  gp.r = r_loc; gp.s = s_loc;
	  bi_linear(&gp, sgrid_ptr);
/* save the gap */
	  x_gap(j,k) = gp.x - x3_vol(i,j,k);
	  y_gap(j,k) = gp.y - y3_vol(i,j,k);
	  z_gap(j,k) = gp.z - z3_vol(i,j,k);
	}
      }
      else{
	printf("compute_projection: Warning: unable to project the point (%e, %e, %e) "
	       "onto the surface mapping `%s'\n", x3_vol(i,j,k), y3_vol(i,j,k), 
	       z3_vol(i,j,k), project->name);
	x_gap(j,k) = 0.0;
	y_gap(j,k) = 0.0;
	z_gap(j,k) = 0.0;
      }
    } /* end for i */

/* update the grid coordinates */
    for (i=1; i<=info->n1; i++){
/* scaled parameter */
      s = (i-1)/((real) info->n1 - 1);
/* use an exponential weighting function */
      if (side == 1){
/* the weighting function is 1 at s=0 and 0 at s=1 */
	w = (exp(-alpha*s) - exp(-alpha))/(1.0-exp(-alpha));
      }
      else{
/* the weighting function is 0 at s=0 and 1 at s=1 */
	w = (exp(-alpha*(1.0-s)) - exp(-alpha))/(1.0-exp(-alpha));
      }

/* Move the k=1 surface eventhough it defines the shape of the surface. In this way */
/* the grid can stay smooth all the way up to the boundary when it is stretched in  */
/* the k-direction */

      for (k=1; k<=info->n3; k++){
	for (j=1; j<=info->n2; j++){
	  x3_vol(i,j,k) = x3_vol(i,j,k) + w*x_gap(j,k);
	  y3_vol(i,j,k) = y3_vol(i,j,k) + w*y_gap(j,k);
	  z3_vol(i,j,k) = z3_vol(i,j,k) + w*z_gap(j,k);
	}
      }
    } /* end for i */

/* cleanup */
    delete_real_array_2d( x_gap_ );
    delete_real_array_2d( y_gap_ );
    delete_real_array_2d( z_gap_ );

  } /* end left or right face */
/* lower or upper face */
  else if (dir == 2){
    n1 = info->n1;
    n3 = info->n3;
/* pick side */
    if (side == 1)
      j = 1;
    else
      j = info->n2;

/* allocate space for the projected surface */
    x_gap_ = create_real_array_2d( n1, n3 );
    y_gap_ = create_real_array_2d( n1, n3 );
    z_gap_ = create_real_array_2d( n1, n3 );

/* do the projection */
    printf("Projecting...\n");
    for (i=1; i<=n1; i++){
      k = 1;
      if ((interp_ptr = search_quad_tree(x3_vol(i,j,k), y3_vol(i,j,k), z3_vol(i,j,k), 
					 sgrid_ptr->quad_tree, sgrid_ptr, 0.0, 0, 0))){
/* save the solution to use as initial guess for next iteration */
	r_loc = interp_ptr->r_loc;
	s_loc = interp_ptr->s_loc;
	n_loc = interp_ptr->n_loc;
/* free the space taken by the interpolation point */
	free(interp_ptr);

/* evaluate the discrete surface */
	gp.r = r_loc; gp.s = s_loc;
	bi_linear(&gp, sgrid_ptr);
/* save the gap */
	x_gap(i,k) = gp.x - x3_vol(i,j,k);
	y_gap(i,k) = gp.y - y3_vol(i,j,k);
	z_gap(i,k) = gp.z - z3_vol(i,j,k);

	for (k=1; k<=n3; k++){
	  newton_search(x3_vol(i,j,k), y3_vol(i,j,k), z3_vol(i,j,k), 
			&r_loc, &s_loc, &n_loc, sgrid_ptr);
/* evaluate the discrete surface */
	  gp.r = r_loc; gp.s = s_loc;
	  bi_linear(&gp, sgrid_ptr);
/* save the gap */
	  x_gap(i,k) = gp.x - x3_vol(i,j,k);
	  y_gap(i,k) = gp.y - y3_vol(i,j,k);
	  z_gap(i,k) = gp.z - z3_vol(i,j,k);
	}
      }
      else{
	printf("compute_projection: Warning: unable to project the point (%e, %e, %e) "
	       "onto the surface mapping `%s'\n", x3_vol(i,j,k), y3_vol(i,j,k), 
	       z3_vol(i,j,k), project->name);
	x_gap(i,k) = 0.0;
	y_gap(i,k) = 0.0;
	z_gap(i,k) = 0.0;
      }
    } /* end for i */

/* update the grid coordinates */
    for (j=1; j<=info->n2; j++){
/* scaled parameter */
      s = (j-1)/((real) info->n2 - 1);
/* use an exponential weighting function */
      if (side == 1){
/* the weighting function is 1 at s=0 and 0 at s=1 */
	w = (exp(-alpha*s) - exp(-alpha))/(1.0-exp(-alpha));
      }
      else{
/* the weighting function is 0 at s=0 and 1 at s=1 */
	w = (exp(-alpha*(1.0-s)) - exp(-alpha))/(1.0-exp(-alpha));
      }

/* Move the k=1 surface eventhough it defines the shape of the surface. In this way */
/* the grid can stay smooth all the way up to the boundary when it is stretched in  */
/* the k-direction */

      for (k=1; k<=info->n3; k++){
	for (i=1; i<=info->n1; i++){
	  x3_vol(i,j,k) = x3_vol(i,j,k) + w*x_gap(i,k);
	  y3_vol(i,j,k) = y3_vol(i,j,k) + w*y_gap(i,k);
	  z3_vol(i,j,k) = z3_vol(i,j,k) + w*z_gap(i,k);
	}
      }
    }

/* cleanup */
    delete_real_array_2d( x_gap_ );
    delete_real_array_2d( y_gap_ );
    delete_real_array_2d( z_gap_ );

  } /* end lower or upper face */

/* delete the component grid used to invert the project_onto mapping */
  delete_discrete_surface( sgrid_ptr );
}

static grid_function *
new_grid_function(int n1, int n2){
  grid_function *u;
  u = (grid_function *) malloc( sizeof(grid_function) );
  u->n1 = n1;
  u->n2 = n2;
  u->x_ = create_real_array_2d(n1, n2);
  u->y_ = create_real_array_2d(n1, n2);
  u->z_ = create_real_array_2d(n1, n2);
  return u;
}

static grid_function *
delete_grid_function(grid_function *u){
  if( u ){
    delete_real_array_2d(u->x_);
    delete_real_array_2d(u->y_);
    delete_real_array_2d(u->z_);
    free(u);
  }
  return NULL;
}

static void 
compute_hyp_grid(input_output *io_ptr, volume_mapping *volume){
  hyp_grid_data *info;
  grid_function *u, *uh, *u_t;
  grid_point surface_point;
  int i,j,k, itr;
  real dr, ds, dt, t, t_max, max_curvature, first_dt=-1.0; 
  const real eps=1.e-10;
  real_array_2d *kappa_;
  eigen_val *eval;
#define u1(i,j) compute_index_2d(u->x_, i, j)
#define u2(i,j) compute_index_2d(u->y_, i, j)
#define u3(i,j) compute_index_2d(u->z_, i, j)

#define uh1(i,j) compute_index_2d(uh->x_, i, j)
#define uh2(i,j) compute_index_2d(uh->y_, i, j)
#define uh3(i,j) compute_index_2d(uh->z_, i, j)

#define u_t1(i,j) compute_index_2d(u_t->x_, i, j)
#define u_t2(i,j) compute_index_2d(u_t->y_, i, j)
#define u_t3(i,j) compute_index_2d(u_t->z_, i, j)

#define c_1(u,i,j) compute_index_2d(u->x_, i, j)
#define c_2(u,i,j) compute_index_2d(u->y_, i, j)
#define c_3(u,i,j) compute_index_2d(u->z_, i, j)

/* cast the data pointer to the right type */
  info = (hyp_grid_data *) volume->volume_data_ptr;

/*   0. free old x, y, z arrays */
  info->x_ = delete_real_array_3d( info->x_ );
  info->y_ = delete_real_array_3d( info->y_ );
  info->z_ = delete_real_array_3d( info->z_ );

/*   1. allocate new x, y, z arrays */
  info->n1 = volume->r1_points;
  info->n2 = volume->r2_points;
  info->n3 = volume->r3_points;

  dr = info->d1 = 1.0/(info->n1-1);
  ds = info->d2 = 1.0/(info->n2-1);
  info->d3 = 1.0/(info->n3-1);

  info->x_ = create_real_array_3d( info->n1, info->n2, info->n3 );
  info->y_ = create_real_array_3d( info->n1, info->n2, info->n3 );
  info->z_ = create_real_array_3d( info->n1, info->n2, info->n3 );

/*   2. allocate 3 gridfunctions, 1 for the present solution, 1 for the middle step, */
/*      and 1 for the time-derivative. */
  u   = new_grid_function(info->n1, info->n2);
  uh  = new_grid_function(info->n1, info->n2);
  u_t = new_grid_function(info->n1, info->n2);

/*   2.5 assign initial data */
  for (j=1; j<=info->n2; j++)
    for (i=1; i<=info->n1; i++){
/* assign the parameters in the surface structure */
      surface_point.r = dr*(i-1);
      surface_point.s = ds*(j-1);
/* Evaluate the surface */
      forward_surface_mapping( &surface_point, info->surface );
/* copy values */
      u1(i,j) = x_hyp(i,j,1) = surface_point.x;
      u2(i,j) = y_hyp(i,j,1) = surface_point.y;
      u3(i,j) = z_hyp(i,j,1) = surface_point.z;
    }

  kappa_ = create_real_array_2d( u->n1, u->n2 );

  if (mean_curvature( u, kappa_, &max_curvature ))
    {
      printf("Max kappa from initial data: %e\n", max_curvature);
      delete_real_array_2d( kappa_ );
    }
  else
    {
      printf("ERROR: compute_hyp_grid: Initial surface singular\n");
      info->status = -1;
/* de-allocate */
      delete_real_array_2d( kappa_ );
      delete_grid_function(u);
      delete_grid_function(uh);
      delete_grid_function(u_t);
      return;
    }

/* initialize the eigenvalue structure*/

  dr = 2.0*((double)M_PI)/(u->n1-1);
  ds = 2.0*((double)M_PI)/(u->n2-1);

/* allocate space */
  eval = (eigen_val *) malloc( sizeof(eigen_val) );
  eval->lambda_1r = eval->lambda_1i = eval->lambda_2r = eval->lambda_2i = 1.e10;
  eval->u_eps   = new_grid_function( u->n1, u->n2 );
  eval->u_eps_t = new_grid_function( u->n1, u->n2 );
  eval->x1    = new_grid_function( u->n1, u->n2 );
  eval->x2    = new_grid_function( u->n1, u->n2 );
  eval->Ax1   = new_grid_function( u->n1, u->n2 );
  eval->Ax2   = new_grid_function( u->n1, u->n2 );

/* initial guesses for Ax1 and Ax2 */
  for (j=1; j<=u->n2; j++)
    for (i=1; i<=u->n1; i++){
      c_1(eval->Ax1,i,j) = (1.0 - 2.0*((i+j)%2)) * cos( (i-1)*dr + (j-1)*ds );
      c_2(eval->Ax1,i,j) = (1.0 - 2.0*((i+j)%2)) * cos( (i-1)*dr + (j-1)*ds );
      c_3(eval->Ax1,i,j) = (1.0 - 2.0*((i+j)%2)) * cos( (i-1)*dr + (j-1)*ds );
      c_1(eval->Ax2,i,j) = (1.0 - 2.0*((i+j)%2)) * sin( (i-1)*dr + (j-1)*ds );
      c_2(eval->Ax2,i,j) = (1.0 - 2.0*((i+j)%2)) * sin( (i-1)*dr + (j-1)*ds );
      c_3(eval->Ax2,i,j) = (1.0 - 2.0*((i+j)%2)) * sin( (i-1)*dr + (j-1)*ds );
    }

/* done initializing the eigenvalue structure */

/*   3. for layer = 2 to layer = info->n3, step 1: */
  for (k=2; k<=info->n3; k++){
    printf("Working on layer %i\n", k);
/* copy previous layer to the time-integration variables */
    for (j=1; j<= info->n2; j++)
      for (i=1; i<= info->n1; i++){
	u1(i,j) = x_hyp(i,j,k-1);
	u2(i,j) = y_hyp(i,j,k-1);
	u3(i,j) = z_hyp(i,j,k-1);
      }
/*   4. Integrate u from time t_{k-1} to t_{k} with RK-2. */
    t     = (k-2)/((real) info->n3-1);
    t_max = (k-1)/((real) info->n3-1);
/*    printf("Working on layer %i. t0 = %e and tmax = %e.\n", k, t, t_max);*/
    itr = 0;
    do{
/* first stage */
      if (!time_derivative( u, t, u_t, info ))
	{
	  /* unsuccessful computation */
	  info->status = -1;
	  
	  /*  delete the gridfunctions the present solution, the middle step, */
	  /*  and the time-derivative. */
	  delete_grid_function(u);
	  delete_grid_function(uh);
	  delete_grid_function(u_t);
	  
	  /* free the eigenvalue structure */
	  delete_grid_function( eval->u_eps );
	  delete_grid_function( eval->u_eps_t );
	  delete_grid_function( eval->x1 );
	  delete_grid_function( eval->x2 );
	  delete_grid_function( eval->Ax1 );
	  delete_grid_function( eval->Ax2 );
	  free(eval);
	  return;
	}

/* compute the eigenvalues. We could if necessary do this in every time-step, if */
/* it turns out to be hard to estimate them analytically. */
      power_method_1( u, u_t, t, info, eval );
      if (eval->lambda_1r > 0.0)
	printf("Warning: Positive dominating eigenvalue = %e\n", eval->lambda_1r); 

/* estimate largest stable dt for R-K-2 */
      dt = 2.0*info->time_factor / fabs(eval->lambda_1r);

      if (itr == 0){
	first_dt = dt;
      }
      else if (dt < first_dt/100.){
	printf(
"Warning: The grid generation was prematurely halted because the time-integration \n"
"seemed to go unstable. Please try again with different parameters. It usually helps \n"
"to decrease the curvature coefficient and increase the averaging coefficient.\n"
"Furthermore, it also helps to make the grid thinner.\n");
/* unsuccessful computation */
	info->status = -1;

/*  delete the gridfunctions the present solution, the middle step, */
/*  and the time-derivative. */
	delete_grid_function(u);
	delete_grid_function(uh);
	delete_grid_function(u_t);

/* free the eigenvalue structure */
	delete_grid_function( eval->u_eps );
	delete_grid_function( eval->u_eps_t );
	delete_grid_function( eval->x1 );
	delete_grid_function( eval->x2 );
	delete_grid_function( eval->Ax1 );
	delete_grid_function( eval->Ax2 );
	free(eval);
	return;
      }
/* estimate largest stable dt for R-K-2 */
      dt = 2.0*info->time_factor / fabs(eval->lambda_1r);

/* do not integrate past t_max */
      dt = real_min( t_max-t, dt );

/* update */
      for (j=1; j<=info->n2; j++)
	for (i=1; i<=info->n1; i++){
	  uh1(i,j) = u1(i,j) + 0.5 * dt * u_t1(i,j);
	  uh2(i,j) = u2(i,j) + 0.5 * dt * u_t2(i,j);
	  uh3(i,j) = u3(i,j) + 0.5 * dt * u_t3(i,j);
	}
/* apply bc */

/* second stage */
      if (!time_derivative( uh, t+0.5*dt, u_t, info ))
	{
	  /* unsuccessful computation */
	  info->status = -1;
	  
	  /*  delete the gridfunctions the present solution, the middle step, */
	  /*  and the time-derivative. */
	  delete_grid_function(u);
	  delete_grid_function(uh);
	  delete_grid_function(u_t);
	  
	  /* free the eigenvalue structure */
	  delete_grid_function( eval->u_eps );
	  delete_grid_function( eval->u_eps_t );
	  delete_grid_function( eval->x1 );
	  delete_grid_function( eval->x2 );
	  delete_grid_function( eval->Ax1 );
	  delete_grid_function( eval->Ax2 );
	  free(eval);
	  return;
	}

/* update */
      for (j=1; j<=info->n2; j++)
	for (i=1; i<=info->n1; i++){
	  u1(i,j) = u1(i,j) + dt * u_t1(i,j);
	  u2(i,j) = u2(i,j) + dt * u_t2(i,j);
	  u3(i,j) = u3(i,j) + dt * u_t3(i,j);
	}
/* apply bc */

/* update time */
      t += dt;
/* increment the iteration counter */
      itr++;
/*      printf("Present time = %e\n", t);*/
    } while (t < t_max-eps);
    

/* assign the new layer from the time-integration variables */
    for (j=1; j<= info->n2; j++)
      for (i=1; i<= info->n1; i++){
	x_hyp(i,j,k) = u1(i,j);
	y_hyp(i,j,k) = u2(i,j);
	z_hyp(i,j,k) = u3(i,j);
      }
  } /* end for all k */

/* successful computation */
  info->status = 1;

/*  delete the gridfunctions the present solution, the middle step, */
/*  and the time-derivative. */
  delete_grid_function(u);
  delete_grid_function(uh);
  delete_grid_function(u_t);

/* free the eigenvalue structure */
  delete_grid_function( eval->u_eps );
  delete_grid_function( eval->u_eps_t );
  delete_grid_function( eval->x1 );
  delete_grid_function( eval->x2 );
  delete_grid_function( eval->Ax1 );
  delete_grid_function( eval->Ax2 );
  free(eval);

} /* end compute_hyp_grid */


static real
scalar_prod( grid_function *u1, grid_function *u2 ){
  real sp;
  int i, j;

  sp = 0.0;
  for (j=1; j<= u1->n2; j++)
    for (i=1; i<= u1->n1; i++){
      sp += c_1(u1,i,j) * c_1(u2,i,j) + c_2(u1,i,j) * c_2(u2,i,j) + 
	c_3(u1,i,j) * c_3(u2,i,j);
    }

  return sp;
}

static void
power_method_1( grid_function *u, grid_function *u_t, real t, hyp_grid_data *info,
	     eigen_val *eval ){
  const real eps=1.e-4;
  real Ax1_norm, lambda_1r, lambda_1i, rel_diff;
  int i, j, n1, n2, itr=0;
  grid_function *x1, *Ax1, *u_eps, *u_eps_t;

  Ax1 = eval->Ax1;
  x1 = eval->x1;
  u_eps = eval->u_eps;
  u_eps_t = eval->u_eps_t;

  n1 = u->n1; n2 = u->n2;

/* copy the eigenvalues to get the iteration going */
  lambda_1r = eval->lambda_1r; lambda_1i = eval->lambda_1i;

/* Iterate: */
  do{
/* save the eigenvalues from the previous iteration */
    eval->lambda_1r = lambda_1r; eval->lambda_1i = lambda_1i;

/* Normalize Ax1, put the result in x1 */
    Ax1_norm = sqrt( scalar_prod( Ax1, Ax1 ) );
    for (j=1; j<=n2; j++)
      for (i=1; i<=n1; i++){
	c_1(x1,i,j) = c_1(Ax1,i,j)/Ax1_norm;
	c_2(x1,i,j) = c_2(Ax1,i,j)/Ax1_norm;
	c_3(x1,i,j) = c_3(Ax1,i,j)/Ax1_norm;
      }

/* Assign u_eps = u + eps*x1, compute the time_derivative corresponding to u_eps */
/* => u_eps_t */
    for (j=1; j<=n2; j++)
      for (i=1; i<=n1; i++){
	c_1(u_eps,i,j) = c_1(u,i,j) + eps * c_1(x1,i,j);
	c_2(u_eps,i,j) = c_2(u,i,j) + eps * c_2(x1,i,j);
	c_3(u_eps,i,j) = c_3(u,i,j) + eps * c_3(x1,i,j);
      }
    if (!time_derivative( u_eps, t, u_eps_t, info ))
      {
	/* unsuccessful computation */
	info->status = -1;
	return;
      }

/* Ax1 = (u_eps_t - u_t)/eps. */
    for (j=1; j<=n2; j++)
      for (i=1; i<=n1; i++){
	c_1(Ax1,i,j) = (c_1(u_eps_t,i,j) - c_1(u_t,i,j))/eps;
	c_2(Ax1,i,j) = (c_2(u_eps_t,i,j) - c_2(u_t,i,j))/eps;
	c_3(Ax1,i,j) = (c_3(u_eps_t,i,j) - c_3(u_t,i,j))/eps;
      }
/* Compute the eigenvalue */
    lambda_1r = scalar_prod( x1, Ax1 );
    
/*     printf("lambda_1 = %e + i %e\n", lambda_1r, lambda_1i); */

    rel_diff = fabs(lambda_1r - eval->lambda_1r) / fabs(lambda_1r);
/* iterate until the eigenvalues have converged, but not for more than 10 iterations */
  } while( ++itr <= 10 && rel_diff > 0.05 );

/* display the eigenvalues */
/*  printf("Dominating eigenvalue converged after %i iterations:\n", itr); */
/*  printf("lambda_1 = %e\n", lambda_1r); */

}

static int
grid_normal(real *n_vec, grid_function *u, int i, int j, hyp_grid_data *info){
  real length, xr, xs, yr, ys, zr, zs;
  int q;
  const real eps = 1.e-10;

/* r-differences */
  if (i == 1){
/*     xr = -3.0*u1(i,j) + 4.0*u1(i+1,j) - u1(i+2,j);  */
/*     yr = -3.0*u2(i,j) + 4.0*u2(i+1,j) - u2(i+2,j); */
/*     zr = -3.0*u3(i,j) + 4.0*u3(i+1,j) - u3(i+2,j); */
    xr = -u1(i,j) + u1(i+1,j);
    yr = -u2(i,j) + u2(i+1,j);
    zr = -u3(i,j) + u3(i+1,j);
  }
  else if (i == u->n1){
/*     xr = 3.0*u1(i,j) - 4.0*u1(i-1,j) + u1(i-2,j);  */
/*     yr = 3.0*u2(i,j) - 4.0*u2(i-1,j) + u2(i-2,j); */
/*     zr = 3.0*u3(i,j) - 4.0*u3(i-1,j) + u3(i-2,j); */
    xr = u1(i,j) - u1(i-1,j);
    yr = u2(i,j) - u2(i-1,j);
    zr = u3(i,j) - u3(i-1,j);
  }
  else{
    xr = u1(i+1,j) - u1(i-1,j);
    yr = u2(i+1,j) - u2(i-1,j);
    zr = u3(i+1,j) - u3(i-1,j);
  }

/* s-differences */
  if (j == 1){
    xs = -u1(i,j) + u1(i,j+1);
    ys = -u2(i,j) + u2(i,j+1);
    zs = -u3(i,j) + u3(i,j+1);
/*     xs = -3.0*u1(i,j) + 4.0*u1(i,j+1) - u1(i,j+2); */
/*     ys = -3.0*u2(i,j) + 4.0*u2(i,j+1) - u2(i,j+2); */
/*     zs = -3.0*u3(i,j) + 4.0*u3(i,j+1) - u3(i,j+2); */
  }
  else if (j == u->n2){
    xs = u1(i,j) - u1(i,j-1);
    ys = u2(i,j) - u2(i,j-1);
    zs = u3(i,j) - u3(i,j-1);
/*     xs = 3.0*u1(i,j) - 4.0*u1(i,j-1) + u1(i,j-2); */
/*     ys = 3.0*u2(i,j) - 4.0*u2(i,j-1) + u2(i,j-2); */
/*     zs = 3.0*u3(i,j) - 4.0*u3(i,j-1) + u3(i,j-2); */
  }
  else{
    xs = u1(i,j+1) - u1(i,j-1);
    ys = u2(i,j+1) - u2(i,j-1);
    zs = u3(i,j+1) - u3(i,j-1);
  }

/* compute normal */
  n_vec[0] = yr * zs - zr * ys;
  n_vec[1] = zr * xs - xr * zs;
  n_vec[2] = xr * ys - yr * xs;

/* apply boundary conditions */
  if (i==1 && hyp_bc(1,1) != 0){
    if (hyp_bc(1,1) == 1)
      n_vec[0] = 0.0;
    else if (hyp_bc(1,1) == 2)
      n_vec[1] = 0.0;
    else if (hyp_bc(1,1) == 3)
      n_vec[2] = 0.0;
  }
  else if (i==info->n1 && hyp_bc(2,1) != 0){
    if (hyp_bc(2,1) == 1)
      n_vec[0] = 0.0;
    else if (hyp_bc(2,1) == 2)
      n_vec[1] = 0.0;
    else if (hyp_bc(2,1) == 3)
      n_vec[2] = 0.0;
  }

  if (j==1 && hyp_bc(1,2) != 0){
    if (hyp_bc(1,2) == 1)
      n_vec[0] = 0.0;
    else if (hyp_bc(1,2) == 2)
      n_vec[1] = 0.0;
    else if (hyp_bc(1,2) == 3)
      n_vec[2] = 0.0;
  }
  else if (j==info->n2 && hyp_bc(2,2) != 0){
    if (hyp_bc(2,2) == 1)
      n_vec[0] = 0.0;
    else if (hyp_bc(2,2) == 2)
      n_vec[1] = 0.0;
    else if (hyp_bc(2,2) == 3)
      n_vec[2] = 0.0;
  }

/* normalize */
  length = sqrt( n_vec[0] * n_vec[0] + n_vec[1] * n_vec[1] + n_vec[2] * n_vec[2] );
  if (length > eps){
    for (q=0; q<3; q++) n_vec[q] = n_vec[q]/length;
    return OK;
  }
  else
    return ERROR;
}

static int
mean_curvature( grid_function *u, real_array_2d *kappa_, real *max_value ){
  real e, f, g, l, m, n, length;
  real x1[3], x2[3], x3[3], x11[3], x12[3], x22[3];
  real dr, ds, dr2, ds2;
  int q, i0, j0, i, j;

#define sp( x, y ) (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])
#define kappa(i,j) compute_index_2d(kappa_, i, j)

  dr = 1.0/(u->n1-1); dr2 = dr*dr;
  ds = 1.0/(u->n2-1); ds2 = ds*ds;

  *max_value = -1.e10;
  for (j0=1; j0<= u->n2; j0++)
    for (i0=1; i0<= u->n1; i0++){

      i = i0; j = j0;

/* r-differences */
      if (i == 1)
	{
	  x1[0] = (u1(i+1,j) - u1(i,j))/dr;
	  x1[1] = (u2(i+1,j) - u2(i,j))/dr;
	  x1[2] = (u3(i+1,j) - u3(i,j))/dr;
	}
      else if (i == u->n1)
	{
	  x1[0] = (u1(i,j) - u1(i-1,j))/dr;
	  x1[1] = (u2(i,j) - u2(i-1,j))/dr;
	  x1[2] = (u3(i,j) - u3(i-1,j))/dr;
	}
      else
	{
	  x1[0] = 0.5*(u1(i+1,j) - u1(i-1,j))/dr;
	  x1[1] = 0.5*(u2(i+1,j) - u2(i-1,j))/dr;
	  x1[2] = 0.5*(u3(i+1,j) - u3(i-1,j))/dr;
	}
      
/* s-differences */
      if (j == 1)
	{
	  x2[0] = (u1(i,j+1) - u1(i,j))/ds;
	  x2[1] = (u2(i,j+1) - u2(i,j))/ds;
	  x2[2] = (u3(i,j+1) - u3(i,j))/ds;
	}
      else if (j == u->n2)
	{
	  x2[0] = (u1(i,j) - u1(i,j-1))/ds;
	  x2[1] = (u2(i,j) - u2(i,j-1))/ds;
	  x2[2] = (u3(i,j) - u3(i,j-1))/ds;
	}
      else
	{
	  x2[0] = 0.5*(u1(i,j+1) - u1(i,j-1))/ds;
	  x2[1] = 0.5*(u2(i,j+1) - u2(i,j-1))/ds;
	  x2[2] = 0.5*(u3(i,j+1) - u3(i,j-1))/ds;
	}

/* shift i and j away from any boundaries */
      i = real_max(2, i); i = real_min(u->n1-1, i);
      j = real_max(2, j); j = real_min(u->n2-1, j);

/* rr-difference */
      x11[0] = (u1(i+1,j) - 2.0*u1(i,j) + u1(i-1,j))/dr2;
      x11[1] = (u2(i+1,j) - 2.0*u2(i,j) + u2(i-1,j))/dr2;
      x11[2] = (u3(i+1,j) - 2.0*u3(i,j) + u3(i-1,j))/dr2;

/* ss-difference */ 
      x22[0] = (u1(i,j+1) - 2.0*u1(i,j) + u1(i,j-1))/dr2;
      x22[1] = (u2(i,j+1) - 2.0*u2(i,j) + u2(i,j-1))/dr2;
      x22[2] = (u3(i,j+1) - 2.0*u3(i,j) + u3(i,j-1))/dr2;

/* rs-differences */
      x12[0] = 0.25*(u1(i+1,j+1) - u1(i-1,j+1) - u1(i+1,j-1) + u1(i-1,j-1))/(dr*ds);
      x12[1] = 0.25*(u2(i+1,j+1) - u2(i-1,j+1) - u2(i+1,j-1) + u2(i-1,j-1))/(dr*ds);
      x12[2] = 0.25*(u3(i+1,j+1) - u3(i-1,j+1) - u3(i+1,j-1) + u3(i-1,j-1))/(dr*ds);

/* compute the normal (x3) */
      x3[0] = x1[1] * x2[2] - x1[2] * x2[1];
      x3[1] = x1[2] * x2[0] - x1[0] * x2[2];
      x3[2] = x1[0] * x2[1] - x1[1] * x2[0];

/* normalize x3 */
      length = sqrt( sp( x3, x3 ) );

      if (length > NEWTON_EPS)
	{
	  for (q=0; q<3; q++) x3[q] = x3[q]/length;

	  e = sp( x1, x1 );
	  f = sp( x1, x2 );
	  g = sp( x2, x2 );

	  l = sp( x11, x3 );
	  m = sp( x12, x3 );
	  n = sp( x22, x3 );

/* the sign seems to be wrong */
	  kappa(i0,j0) = -0.5 * (e*n - 2.0*f*m + g*l) / (e*g - f*f);
	}
      else
	{
	  printf("Warning: mean_curvature: The mapping is singular at i=%i, j=%i\n", 
		 i0, j0);
	  return ERROR;
	}
      if (kappa(i0,j0) > *max_value) *max_value = kappa(i0,j0);

    }
  return OK;
#undef sp
}

static int
time_derivative(grid_function *u, real t, grid_function *u_t, hyp_grid_data *info){
  int i, j;
  real n_vec[3], v, ds_dt, t_stretch, dummy, max_value;
  grid_function *u_tmp;
  real_array_2d *kappa_;

#define u_tmp1(i,j) compute_index_2d(u_tmp->x_, i, j)
#define u_tmp2(i,j) compute_index_2d(u_tmp->y_, i, j)
#define u_tmp3(i,j) compute_index_2d(u_tmp->z_, i, j)

  u_tmp  = new_grid_function( u->n1, u->n2 );
  kappa_ = create_real_array_2d( u->n1, u->n2 );

/* 2. compute the mean curvature */
  if (!mean_curvature( u, kappa_, &max_value ))
    return ERROR;

/* stretching in the normal direction? */
  uniform_to_stretch(t, info->normal_stretch, &t_stretch, &ds_dt, &dummy );
/* The time-derivative of (x,y,z)^T is computed as follows: */
  for (j=1; j<=u->n2; j++)
    for (i=1; i<=u->n1; i++){
/* 1. Compute the normals */
      grid_normal( n_vec, u, i, j, info );
/* 3. Compute the velocity function */
      v = info->width * ds_dt *
	real_max( info->v_min, 1 - info->curv_factor * kappa(i,j));
      u_tmp1(i,j) = v * n_vec[0];
      u_tmp2(i,j) = v * n_vec[1];
      u_tmp3(i,j) = v * n_vec[2];
    }
/* 4. smooth the time-derivative for all interior points */
   for (j=2; j<=u->n2-1; j++) 
     for (i=2; i<=u->n1-1; i++){ 
       u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	 info->smooth_factor *
	   (u_tmp1(i-1,j) + u_tmp1(i+1,j) + u_tmp1(i,j) + u_tmp1(i,j-1) + 
	    u_tmp1(i,j+1))/5.0;
       u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	 info->smooth_factor *
	   (u_tmp2(i-1,j) + u_tmp2(i+1,j) + u_tmp2(i,j) + u_tmp2(i,j-1) + 
	    u_tmp2(i,j+1))/5.0;
       u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	 info->smooth_factor *
	   (u_tmp3(i-1,j) + u_tmp3(i+1,j) + u_tmp3(i,j) + u_tmp3(i,j-1) + 
	    u_tmp3(i,j+1))/5.0;
     } 

/* smooth along the boundary for the boundary values */
/* low i */
  i = 1;
  if (hyp_bc(1,1) == 0){
    for (j=2; j<=u->n2-1; j++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i,j-1) + u_tmp1(i,j+1) + 2.0*u_tmp1(i+1,j) )/4.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i,j-1) + u_tmp2(i,j+1) + 2.0*u_tmp2(i+1,j) )/4.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i,j-1) + u_tmp3(i,j+1) + 2.0*u_tmp3(i+1,j) )/4.0;
    } /* end for */
  } /* end free bc */
  else if (hyp_bc(1,1) == 1){
    for (j=2; j<=u->n2-1; j++){ 
      u_t1(i,j) = 0.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i,j-1) + u_tmp2(i,j+1) + 2.0*u_tmp2(i+1,j) )/4.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i,j-1) + u_tmp3(i,j+1) + 2.0*u_tmp3(i+1,j) )/4.0;
    } /* end for */
  } /* end x=constant bc */
  else if (hyp_bc(1,1) == 2){
    for (j=2; j<=u->n2-1; j++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i,j-1) + u_tmp1(i,j+1) + 2.0*u_tmp1(i+1,j) )/4.0;
      u_t2(i,j) = 0.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i,j-1) + u_tmp3(i,j+1) + 2.0*u_tmp3(i+1,j) )/4.0;
    } /* end for */
  } /* end y=constant bc */
  else if (hyp_bc(1,1) == 3){
    for (j=2; j<=u->n2-1; j++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i,j-1) + u_tmp1(i,j+1) + 2.0*u_tmp1(i+1,j) )/4.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i,j-1) + u_tmp2(i,j+1) + 2.0*u_tmp2(i+1,j) )/4.0;
      u_t3(i,j) = 0.0;
    } /* end for */
  } /* end z=constant bc */

/* high i */
  i = u->n1;
  if (hyp_bc(2,1) == 0){
    for (j=2; j<=u->n2-1; j++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i,j-1) + u_tmp1(i,j+1) + 2.0*u_tmp1(i-1,j) )/4.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i,j-1) + u_tmp2(i,j+1) + 2.0*u_tmp2(i-1,j) )/4.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i,j-1) + u_tmp3(i,j+1) + 2.0*u_tmp3(i-1,j) )/4.0;
    } /* end for */
  }
  else if (hyp_bc(2,1) == 1){
    for (j=2; j<=u->n2-1; j++){ 
      u_t1(i,j) = 0.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i,j-1) + u_tmp2(i,j+1) + 2.0*u_tmp2(i-1,j) )/4.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i,j-1) + u_tmp3(i,j+1) + 2.0*u_tmp3(i-1,j) )/4.0;
    } /* end for */
  }
  else if (hyp_bc(2,1) == 2){
    for (j=2; j<=u->n2-1; j++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i,j-1) + u_tmp1(i,j+1) + 2.0*u_tmp1(i-1,j) )/4.0;
      u_t2(i,j) = 0.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i,j-1) + u_tmp3(i,j+1) + 2.0*u_tmp3(i-1,j) )/4.0;
    } /* end for */
  }
  else if (hyp_bc(2,1) == 3){
    for (j=2; j<=u->n2-1; j++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i,j-1) + u_tmp1(i,j+1) + 2.0*u_tmp1(i-1,j) )/4.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i,j-1) + u_tmp2(i,j+1) + 2.0*u_tmp2(i-1,j) )/4.0;
      u_t3(i,j) = 0.0;
    } /* end for */
  }


/* low j */
  j = 1;
  if (hyp_bc(1,2) == 0){
    for (i=2; i<=u->n1-1; i++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i-1,j) + u_tmp1(i+1,j) + 2.0*u_tmp1(i,j+1) )/4.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i-1,j) + u_tmp2(i+1,j) + 2.0*u_tmp2(i,j+1) )/4.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i-1,j) + u_tmp3(i+1,j) + 2.0*u_tmp3(i,j+1) )/4.0;
    } /* end for */
  }
  else if (hyp_bc(1,2) == 1){
    for (i=2; i<=u->n1-1; i++){ 
      u_t1(i,j) = 0.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i-1,j) + u_tmp2(i+1,j) + 2.0*u_tmp2(i,j+1) )/4.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i-1,j) + u_tmp3(i+1,j) + 2.0*u_tmp3(i,j+1) )/4.0;
    } /* end for */
  }
  else if (hyp_bc(1,2) == 2){
    for (i=2; i<=u->n1-1; i++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i-1,j) + u_tmp1(i+1,j) + 2.0*u_tmp1(i,j+1) )/4.0;
      u_t2(i,j) = 0.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i-1,j) + u_tmp3(i+1,j) + 2.0*u_tmp3(i,j+1) )/4.0;
    } /* end for */
  }
  else if (hyp_bc(1,2) == 3){
    for (i=2; i<=u->n1-1; i++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i-1,j) + u_tmp1(i+1,j) + 2.0*u_tmp1(i,j+1) )/4.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i-1,j) + u_tmp2(i+1,j) + 2.0*u_tmp2(i,j+1) )/4.0;
      u_t3(i,j) = 0.0;
    } /* end for */
  }


/* high j */
  j = u->n2;
  if (hyp_bc(2,2) == 0){
    for (i=2; i<=u->n1-1; i++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i-1,j) + u_tmp1(i+1,j) + 2.0*u_tmp1(i,j-1) )/4.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i-1,j) + u_tmp2(i+1,j) + 2.0*u_tmp2(i,j-1) )/4.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i-1,j) + u_tmp3(i+1,j) + 2.0*u_tmp3(i,j-1) )/4.0;
    }
  }
  else if (hyp_bc(2,2) == 1){
    for (i=2; i<=u->n1-1; i++){ 
      u_t1(i,j) = 0.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i-1,j) + u_tmp2(i+1,j) + 2.0*u_tmp2(i,j-1) )/4.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i-1,j) + u_tmp3(i+1,j) + 2.0*u_tmp3(i,j-1) )/4.0;
    }
  }
  else if (hyp_bc(2,2) == 2){
    for (i=2; i<=u->n1-1; i++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i-1,j) + u_tmp1(i+1,j) + 2.0*u_tmp1(i,j-1) )/4.0;
      u_t2(i,j) = 0.0;
      u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
	info->smooth_factor *
	  (u_tmp3(i-1,j) + u_tmp3(i+1,j) + 2.0*u_tmp3(i,j-1) )/4.0;
    }
  }
  else if (hyp_bc(2,2) == 3){
    for (i=2; i<=u->n1-1; i++){ 
      u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
	info->smooth_factor *
	  (u_tmp1(i-1,j) + u_tmp1(i+1,j) + 2.0*u_tmp1(i,j-1) )/4.0;
      u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
	info->smooth_factor *
	  (u_tmp2(i-1,j) + u_tmp2(i+1,j) + 2.0*u_tmp2(i,j-1) )/4.0;
      u_t3(i,j) = 0.0;
    }
  }


/* smooth around the corners */
/* lower left corner */
  i=1; j=1;
  if (hyp_bc(1,1) == 1 || hyp_bc(1,2) == 1)
    u_t1(i,j) = 0.0;
  else
    u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp1(i+1,j) + 2.0*u_tmp1(i,j+1) )/4.0;

  if (hyp_bc(1,1) == 2 || hyp_bc(1,2) == 2)
    u_t2(i,j) = 0.0;
  else
    u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp2(i+1,j) + 2.0*u_tmp2(i,j+1) )/4.0;

  if (hyp_bc(1,1) == 3 || hyp_bc(1,2) == 3)
    u_t3(i,j) = 0.0;
  else      
    u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp3(i+1,j) + 2.0*u_tmp3(i,j+1) )/4.0;

/* lower right corner */
  i = u->n1; j=1;
  if (hyp_bc(2,1) == 1 || hyp_bc(1,2) == 1)
    u_t1(i,j) = 0.0;
  else
    u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp1(i,j+1) + 2.0*u_tmp1(i-1,j) )/4.0;

  if (hyp_bc(2,1) == 2 || hyp_bc(1,2) == 2)
    u_t2(i,j) = 0.0;
  else
    u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp2(i,j+1) + 2.0*u_tmp2(i-1,j) )/4.0;

  if (hyp_bc(2,1) == 3 || hyp_bc(1,2) == 3)
    u_t3(i,j) = 0.0;
  else
    u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp3(i,j+1) + 2.0*u_tmp3(i-1,j) )/4.0;

/* upper left corner */
  i = 1; j=u->n2;
  if (hyp_bc(1,1) == 1 || hyp_bc(2,2) == 1)
    u_t1(i,j) = 0.0;
  else
    u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp1(i+1,j) + 2.0*u_tmp1(i,j-1) )/4.0;

  if (hyp_bc(1,1) == 2 || hyp_bc(2,2) == 2)
    u_t2(i,j) = 0.0;
  else
    u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp2(i+1,j) + 2.0*u_tmp2(i,j-1) )/4.0;

  if (hyp_bc(1,1) == 3 || hyp_bc(2,2) == 3)
    u_t3(i,j) = 0.0;
  else
    u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp3(i+1,j) + 2.0*u_tmp3(i,j-1) )/4.0;

/* upper right corner */
  i = u->n1; j = u->n2;
  if (hyp_bc(2,1) == 1 || hyp_bc(2,2) == 1)
    u_t1(i,j) = 0.0;
  else
    u_t1(i,j) = (1.0-info->smooth_factor) * u_tmp1(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp1(i-1,j) + 2.0*u_tmp1(i,j-1) )/4.0;

  if (hyp_bc(2,1) == 2 || hyp_bc(2,2) == 2)
    u_t2(i,j) = 0.0;
  else
    u_t2(i,j) = (1.0-info->smooth_factor) * u_tmp2(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp2(i-1,j) + 2.0*u_tmp2(i,j-1) )/4.0;

  if (hyp_bc(2,1) == 3 || hyp_bc(2,2) == 3)
    u_t3(i,j) = 0.0;
  else
    u_t3(i,j) = (1.0-info->smooth_factor) * u_tmp3(i,j) +  
      info->smooth_factor *
	(2.0*u_tmp3(i-1,j) + 2.0*u_tmp3(i,j-1) )/4.0;

/* free the memory */
  delete_grid_function( u_tmp );
  delete_real_array_2d( kappa_ );

  return OK;
}

hyp_grid_data *
init_hyp_grid(input_output *io_ptr, volume_mapping *volume, surface_mapping *surface){
  hyp_grid_data *info;
  
/* allocate memory */
  info = (hyp_grid_data *) malloc( sizeof(hyp_grid_data) );

/* general data */
  volume->type = "Hyperbolic"; 
  volume->volume_data_ptr = (void *) info; 
  volume->forward_mapping = hyp_grid_mapping; 
  volume->inverse_mapping = NULL;
  volume->set_volume = set_hyp_grid; 
  volume->delete_volume_data = delete_hyp_grid;
  
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
  info->width = 0.5; 
  info->normal_stretch = NULL;
  info->v_min = 0.25;
  info->curv_factor = 0.25;
  info->smooth_factor = 0.5;
  info->time_factor = 1.0;
  info->status = 0;
  info->hyp_bc_ = create_int_array_2d(2,2);
/* initialize bc */
  hyp_bc(1,1) = hyp_bc(2,1) = hyp_bc(1,2) = hyp_bc(2,2) = 0;

/* initially, copy the bounding box from the surface */
  volume->bb->x_min = surface->bb->x_min;
  volume->bb->x_max = surface->bb->x_max; 
  volume->bb->y_min = surface->bb->y_min; 
  volume->bb->y_max = surface->bb->y_max; 
  volume->bb->z_min = surface->bb->z_min; 
  volume->bb->z_max = surface->bb->z_max; 

/* no volume arrays to begin with */
  info->x_ = NULL;
  info->y_ = NULL;
  info->z_ = NULL;

  return info;
}

static void
delete_hyp_grid( void *volume_data_ptr ){
  hyp_grid_data *info;
  
/* cast the pointer */
  info = (hyp_grid_data *) volume_data_ptr;

/* mark that the surface is no longer used by the volume mapping */
  info->surface->used_by_volume--;

/* free the boundary condition array */
  info->hyp_bc_ = delete_int_array_2d(info->hyp_bc_);

/* free the stretching */
  info->normal_stretch = delete_generic_stretching( info->normal_stretch );

/* free the arrays */
  info->x_ = delete_real_array_3d( info->x_ );
  info->y_ = delete_real_array_3d( info->y_ );
  info->z_ = delete_real_array_3d( info->z_ );

/* free the structure */
  free(info);
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


#define R0 ((i_loc-1)   * info->d1)
#define R1 (i_loc * info->d1)
#define ALPHA_0(r)   ((R1 - r)/info->d1)
#define D_ALPHA_0(r) (-1.0/info->d1)
#define ALPHA_1(r)   ((r - R0)/info->d1)
#define D_ALPHA_1(r) (1.0/info->d1)

#define S0 ((j_loc-1)  * info->d2)
#define S1 (j_loc * info->d2)
#define BETA_0(s)   ((S1 - s)/info->d2)
#define D_BETA_0(s) (-1.0/info->d2)
#define BETA_1(s)   ((s - S0)/info->d2)
#define D_BETA_1(s) (1.0/info->d2)

#define T0 ((k_loc-1)  * info->d3)
#define T1 (k_loc * info->d3)
#define GAMMA_0(t)   ((T1 - t)/info->d3)
#define D_GAMMA_0(t) (-1.0/info->d3)
#define GAMMA_1(t)   ((t - T0)/info->d3)
#define D_GAMMA_1(t) (1.0/info->d3)

#define LIN(v,r,s,t) (ALPHA_0(r) * BETA_0(s) * GAMMA_0(t) * v(i_loc  ,j_loc  ,k_loc  ) +  \
		      ALPHA_1(r) * BETA_0(s) * GAMMA_0(t) * v(i_loc+1,j_loc  ,k_loc  ) +  \
		      ALPHA_0(r) * BETA_1(s) * GAMMA_0(t) * v(i_loc  ,j_loc+1,k_loc  ) +  \
		      ALPHA_1(r) * BETA_1(s) * GAMMA_0(t) * v(i_loc+1,j_loc+1,k_loc  ) +  \
		      ALPHA_0(r) * BETA_0(s) * GAMMA_1(t) * v(i_loc  ,j_loc  ,k_loc+1) +  \
		      ALPHA_1(r) * BETA_0(s) * GAMMA_1(t) * v(i_loc+1,j_loc  ,k_loc+1) +  \
		      ALPHA_0(r) * BETA_1(s) * GAMMA_1(t) * v(i_loc  ,j_loc+1,k_loc+1) +  \
		      ALPHA_1(r) * BETA_1(s) * GAMMA_1(t) * v(i_loc+1,j_loc+1,k_loc+1))
#define LIN_R(v,r,s,t) (D_ALPHA_0(r) * BETA_0(s) * GAMMA_0(t) * v(i_loc  ,j_loc  ,k_loc  ) +  \
			D_ALPHA_1(r) * BETA_0(s) * GAMMA_0(t) * v(i_loc+1,j_loc  ,k_loc  ) +  \
			D_ALPHA_0(r) * BETA_1(s) * GAMMA_0(t) * v(i_loc  ,j_loc+1,k_loc  ) +  \
			D_ALPHA_1(r) * BETA_1(s) * GAMMA_0(t) * v(i_loc+1,j_loc+1,k_loc  ) +  \
			D_ALPHA_0(r) * BETA_0(s) * GAMMA_1(t) * v(i_loc  ,j_loc  ,k_loc+1) +  \
			D_ALPHA_1(r) * BETA_0(s) * GAMMA_1(t) * v(i_loc+1,j_loc  ,k_loc+1) +  \
			D_ALPHA_0(r) * BETA_1(s) * GAMMA_1(t) * v(i_loc  ,j_loc+1,k_loc+1) +  \
			D_ALPHA_1(r) * BETA_1(s) * GAMMA_1(t) * v(i_loc+1,j_loc+1,k_loc+1))
#define LIN_S(v,r,s,t) (ALPHA_0(r) * D_BETA_0(s) * GAMMA_0(t) * v(i_loc  ,j_loc  ,k_loc  ) +  \
			ALPHA_1(r) * D_BETA_0(s) * GAMMA_0(t) * v(i_loc+1,j_loc  ,k_loc  ) +  \
			ALPHA_0(r) * D_BETA_1(s) * GAMMA_0(t) * v(i_loc  ,j_loc+1,k_loc  ) +  \
			ALPHA_1(r) * D_BETA_1(s) * GAMMA_0(t) * v(i_loc+1,j_loc+1,k_loc  ) +  \
			ALPHA_0(r) * D_BETA_0(s) * GAMMA_1(t) * v(i_loc  ,j_loc  ,k_loc+1) +  \
			ALPHA_1(r) * D_BETA_0(s) * GAMMA_1(t) * v(i_loc+1,j_loc  ,k_loc+1) +  \
			ALPHA_0(r) * D_BETA_1(s) * GAMMA_1(t) * v(i_loc  ,j_loc+1,k_loc+1) +  \
			ALPHA_1(r) * D_BETA_1(s) * GAMMA_1(t) * v(i_loc+1,j_loc+1,k_loc+1))
#define LIN_T(v,r,s,t) (ALPHA_0(r) * BETA_0(s) * D_GAMMA_0(t) * v(i_loc  ,j_loc  ,k_loc  ) +  \
			ALPHA_1(r) * BETA_0(s) * D_GAMMA_0(t) * v(i_loc+1,j_loc  ,k_loc  ) +  \
			ALPHA_0(r) * BETA_1(s) * D_GAMMA_0(t) * v(i_loc  ,j_loc+1,k_loc  ) +  \
			ALPHA_1(r) * BETA_1(s) * D_GAMMA_0(t) * v(i_loc+1,j_loc+1,k_loc  ) +  \
			ALPHA_0(r) * BETA_0(s) * D_GAMMA_1(t) * v(i_loc  ,j_loc  ,k_loc+1) +  \
			ALPHA_1(r) * BETA_0(s) * D_GAMMA_1(t) * v(i_loc+1,j_loc  ,k_loc+1) +  \
			ALPHA_0(r) * BETA_1(s) * D_GAMMA_1(t) * v(i_loc  ,j_loc+1,k_loc+1) +  \
			ALPHA_1(r) * BETA_1(s) * D_GAMMA_1(t) * v(i_loc+1,j_loc+1,k_loc+1))

static int 
hyp_grid_mapping( grid_point *gp, void *data_ptr){
  hyp_grid_data *info;
  int i_loc, j_loc, k_loc;

/* cast the data pointer to the right type */
  info = (hyp_grid_data *) data_ptr;

/* get the nearest grid point below and to the left of (r,s) */
  i_loc = int_max( 1, 1 + (int) ( (info->n1-1)*gp->r ) );
  j_loc = int_max( 1, 1 + (int) ( (info->n2-1)*gp->s ) );
  k_loc = int_max( 1, 1 + (int) ( (info->n3-1)*gp->t ) );

/* prevent overflow */
  i_loc = int_min( i_loc, info->n1-1 );
  j_loc = int_min( j_loc, info->n2-1 );
  k_loc = int_min( k_loc, info->n3-1 );

/* evaluate the tri-linear mapping */
  gp->x  = LIN(x_hyp, gp->r, gp->s, gp->t);
  gp->y  = LIN(y_hyp, gp->r, gp->s, gp->t);
  gp->z  = LIN(z_hyp, gp->r, gp->s, gp->t);

  gp->xr = LIN_R(x_hyp, gp->r, gp->s, gp->t);
  gp->yr = LIN_R(y_hyp, gp->r, gp->s, gp->t);
  gp->zr = LIN_R(z_hyp, gp->r, gp->s, gp->t);

  gp->xs = LIN_S(x_hyp, gp->r, gp->s, gp->t);
  gp->ys = LIN_S(y_hyp, gp->r, gp->s, gp->t);
  gp->zs = LIN_S(z_hyp, gp->r, gp->s, gp->t);

  gp->xt = LIN_T(x_hyp, gp->r, gp->s, gp->t);
  gp->yt = LIN_T(y_hyp, gp->r, gp->s, gp->t);
  gp->zt = LIN_T(z_hyp, gp->r, gp->s, gp->t);

  return OK;
}
