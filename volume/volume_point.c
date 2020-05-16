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

/* system prototype */
extern int sscanf( const char *line, const char *format, ...);

/* private member functions */
static void
compute_projection(input_output *io_ptr, volume_mapping *volume, 
		   surface_mapping *project );
static void
alloc_volume_point_arrays( volume_point_info *info );
static void
free_volume_point_arrays( volume_point_info *info );
static volume_point_info *
read_plot3d_file( FILE *fp, input_output *io_ptr );
static void 
set_volume_point(input_output *io_ptr, volume_mapping *volume,
		 linked_list *surface_list);
static void
delete_volume_point( void *volume_data_ptr );
static int 
volume_point_mapping(grid_point *gp, void *data_ptr );
/* end private member functions */

static void
alloc_volume_point_arrays( volume_point_info *info ){
  info->x_ = create_real_array_3d(info->n1, info->n2, info->n3);
  info->y_ = create_real_array_3d(info->n1, info->n2, info->n3);
  info->z_ = create_real_array_3d(info->n1, info->n2, info->n3);
}

static void
free_volume_point_arrays( volume_point_info *info ){
  info->x_ = delete_real_array_3d(info->x_);
  info->y_ = delete_real_array_3d(info->y_);
  info->z_ = delete_real_array_3d(info->z_);
}

static volume_point_info *
read_plot3d_file( FILE *fp, input_output *io_ptr ){
  int n1, n2, n3, i, j, k, n_grids, which_grid=1, dum, read_past;
  real x_dum;
  char question[80], first_line[120];
  volume_point_info *info;

/* formatted plot3d file */
/* first try to determine if more than one grid is present in the file */
  read_past = 0;
  fgets( first_line, 120, fp );
  if (sscanf( first_line, "%i%i%i", &n1, &n2, &n3 ) == 1){
    n_grids = n1;
    if (n_grids == 1){
      which_grid = 1;
      fscanf( fp, "%i%i%i", &n1, &n2, &n3 );
    }
    else{
      sprintf(question, "There are %i grids in this file. Which one do "
	      "you want to read: ", n_grids);
      which_grid = get_int( io_ptr, question, 1, 0);
      which_grid = int_max(1, int_min( n_grids, which_grid) );
      for (i=1; i<which_grid; i++){
	fscanf( fp, "%i%i%i", &n1, &n2, &n3 );
	read_past += n1*n2*n3;
      }
      fscanf( fp, "%i%i%i", &n1, &n2, &n3 );
      for (i=which_grid+1; i<=n_grids; i++)
	fscanf( fp, "%i%i%i", &dum, &dum, &dum );
    }
  }
/* tmp */
  printf("The dimensions were n1 = %i, n2 = %i, n3 = %i\n", 
	 n1, n2, n3);
/* end tmp */

/* check if it is 3d */
  if( n1 <= 1 || n2 <= 1 || n3 <= 1 ){
    printf("The dimensions in the file make no sense. This is supposed to be"
	   " a 3-d grid,\n"
	   "but n1 = %i, n2 = %i, n3 = %i.\n", n1, n2, n3);
    fclose( fp );
    return ERROR;
  }
  
/* read past the first which_grid-1 surfaces */  
  for (i=1; i<=3*read_past; i++){
    if ((fscanf( fp, "%lg", &x_dum ) ) == EOF){
	  printf("Input error in read_plot3d while reading past the "
		 "first %i grids.\n", which_grid - 1);
	  fclose( fp );
	  return NULL;
	}
  }

  info = (volume_point_info *) malloc( sizeof(volume_point_info) );

/* set the dimensions */
  info->n1 = n1;
  info->n2 = n2;
  info->n3 = n3;

/* the grid step is used by local_map */
  info->d1 = 1.0/((real) n1-1);
  info->d2 = 1.0/((real) n2-1);
  info->d3 = 1.0/((real) n3-1);

/* allocate x, y, z arrays */
  alloc_volume_point_arrays( info );

/* read the x coordinates */
  for (k = 1; k <= info->n3; k++)
    for (j = 1; j <= info->n2; j++)
      for (i = 1; i <= info->n1; i++){
	if ((fscanf( fp, "%lg", &(x3_vol(i,j,k)) )) == EOF){
	  printf("Input error in read_plot3d while reading x-coordinates at "
		 "i = %i, j= %i, k= %i.\n", i, j, k);
	  free_volume_point_arrays( info );
	  fclose( fp );
	  return NULL;
	}
      }  

/* read the y coordinates */
  for (k = 1; k <= info->n3; k++)
    for (j = 1; j <= info->n2; j++)
      for (i = 1; i <= info->n1; i++){
	if ((fscanf( fp, "%lg", &(y3_vol(i,j,k)) )) == EOF){
	  printf("Input error in read_plot3d while reading y-coordinates at "
		 "i = %i, j= %i, k= %i.\n", i, j, k);
	  free_volume_point_arrays( info );
	  fclose( fp );
	  return NULL;
	}
      }

/* read the z coordinates */
  for (k = 1; k <= info->n3; k++)
    for (j = 1; j <= info->n2; j++)
      for (i = 1; i <= info->n1; i++){
	if ((fscanf( fp, "%lg", &(z3_vol(i,j,k)) )) == EOF){
	  printf("Input error in read_plot3d while reading z-coordinates at "
		 "i = %i, j= %i, k= %i.\n", i, j, k);
	  free_volume_point_arrays( info );
	  fclose( fp );
	  return NULL;
	}
      }
  
/* successful completion */
  fclose( fp );

  return info;
}


static void 
set_volume_point(input_output *io_ptr, volume_mapping *volume,
		 linked_list *surface_list){
  char prompt[80];
  int icom, quit=0, replot=1, i, j, k, i_max, j_max, k_max;
  surface_mapping *project;
  linked_list_member *project_link;
  grid_point gp0, gp1;
  real diag, max_diag, dx, dy, dz;

#include "set_volume_point_com.h"

  sprintf(prompt, "%s: volume point mapping>", volume->name);

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

    switch( get_command( io_ptr, prompt, COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

/* change the plot mode */
    case PLOT_MODE:
      volume->plot_mode = volume_plot_mode(io_ptr, NULL, volume, volume->plot_mode,
					   NULL);
/* don't need to redraw the grid! */
      break;

    case SPACE_DIAGONAL:
      max_diag = 0.0;
      i_max = j_max = k_max = 0;
      for (k=1; k<volume->r3_points; k++)
	for (j=1; j<volume->r2_points; j++)
	  for (i=1; i<volume->r1_points; i++){
	    gp0.r = (i-1)/((real) volume->r1_points - 1);
	    gp0.s = (j-1)/((real) volume->r2_points - 1);
	    gp0.t = (k-1)/((real) volume->r3_points - 1);
	    volume_point_mapping(&gp0, volume->volume_data_ptr);

	    gp1.r = i/((real) volume->r1_points - 1);
	    gp1.s = j/((real) volume->r2_points - 1);
	    gp1.t = k/((real) volume->r3_points - 1);
	    volume_point_mapping(&gp1, volume->volume_data_ptr);

	    dx = gp1.x - gp0.x;
	    dy = gp1.y - gp0.y;
	    dz = gp1.z - gp0.z;
	    diag = sqrt(dx*dx+dy*dy+dz*dz);
	    if (diag > max_diag){
	      max_diag = diag;
	      i_max = i;
	      j_max = j;
	      k_max = k;
	    } /* end if new largest */
	  } /* end for all cells */	    
      printf("Largest space diagonal = %e, occured in grid cell (%i, %i, %i)\n", 
	     max_diag, i_max, j_max, k_max);
      break;

    case R1_LINES:
      volume->r1_points = int_max(2, get_int(io_ptr, "Number of r1-lines:", 
					     volume->r1_points, NO_INDENT));
      replot = 1;
      break;

    case R2_LINES:
      volume->r2_points = int_max(2, get_int(io_ptr, "Number of r2-lines:", 
					     volume->r2_points, NO_INDENT));
      replot = 1;
      break;

    case R3_LINES:
      volume->r3_points = int_max(2, get_int(io_ptr, "Number of r3-lines:", 
					     volume->r3_points, NO_INDENT));
      replot = 1;
      break;

    case PROJECT_FACE:
      replot = 1;
      printf("Identify the surface to project onto:\n");
      project_link = get_surface( io_ptr, surface_list, NO_INDENT );
/* compute the projection */  
      if (project_link){
	project = project_link->data;
	compute_projection(io_ptr, volume, project);
      }
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

    default:
      break;

    }
  } while (!quit);
}

static void
compute_projection( input_output *io_ptr, volume_mapping *volume, surface_mapping *project ){

  char *prompt;
  int icom, quit=0, side = 0, dir = 0, n1, n2, n3, i, j, k;
  discrete_surface *sgrid_ptr;
  real_array_2d *x_gap_, *y_gap_, *z_gap_, *x_project, *y_project, *z_project;
  volume_point_info *info;
  inverse_point *interp_ptr;
  grid_point gp;
  real s, w, r_loc, s_loc, n_loc, l_scale;
  const real alpha=18.42; /* makes the weight function satisfy w(0.5) = 0.0001 */

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

/*    case  NEAR: */
/*       printf("Sorry, projection of the near (r3=0) face is not implemented yet\n"); */
/*       side = 1; dir = 3; */
/*       break; */

/*    case MY_FAR: */
/*      printf("Sorry, projection of the upper (r3=1) face is not implemented yet\n"); */
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
  sgrid_ptr->max_n_dist = 0.01 * l_scale;

/* assumes that the mapping is of type volume_point */
  info = (volume_point_info *) volume->volume_data_ptr;

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

void *
init_volume_point(input_output *io_ptr, volume_mapping *volume ){
  FILE *fp;
  char *file_name, *prompt;
  int icom, quit=0, i, j, k;
  const int save_command=1;
  volume_point_info *info=NULL;
  real box_size, dist;

#include "init_volume_point_com.h"

  prompt = "select file format for reading the grid points>";

  do{

    switch (get_command(io_ptr, prompt, COMMAND, LAST_COM, LEVEL,
		       SAVE_ON_COPY, ARGUMENT)) {

    case PLOT3D:
/* plot3d ascii format */
      if ( (fp = open_ascii_file( io_ptr, "Enter plot3d ascii file: ", 
				 "test.plot3d", &file_name, 'r', 0, 
				 save_command)) == NULL ||
	  (info = read_plot3d_file( fp, io_ptr )) == NULL )
	return NULL;
      quit = 1;
      break;

    case HELP:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
/* cancel */
      return NULL;
      break;

    default:
      ;

    }
  } while( !quit );


/* general data */
  volume->type = "Volume point"; 
  volume->volume_data_ptr = (void *) info; 
  volume->forward_mapping = volume_point_mapping; 
  volume->inverse_mapping = NULL; 
  volume->set_volume = set_volume_point; 
  volume->delete_volume_data = delete_volume_point;

  volume->r1_points = info->n1;
  volume->r2_points = info->n2;
  volume->r3_points = info->n3;

/* compute bounding box */
  volume_bb( volume );
  box_size = (volume->bb->x_max - volume->bb->x_min +
	      volume->bb->y_max - volume->bb->y_min +
	      volume->bb->z_max - volume->bb->z_min)/3.0;

/* check for periodicity in direction 1 */
  dist = 0;
  for (k=1; k<=volume->r3_points; k++)
    for (j=1; j<=volume->r2_points; j++)
      dist += fabs(x3_vol(1,j,k) - x3_vol(volume->r1_points,j,k)) +
	fabs(y3_vol(1,j,k) - y3_vol(volume->r1_points,j,k)) +
	  fabs(z3_vol(1,j,k) - z3_vol(volume->r1_points,j,k));

  volume->r1_period = (dist <= 1.e-4 * box_size)? 1 : 0;

/* check for periodicity in direction 2 */
  dist = 0;
  for (k=1; k<=volume->r3_points; k++)
    for (i=1; i<=volume->r1_points; i++)
      dist += fabs(x3_vol(i,1,k) - x3_vol(i,volume->r2_points,k)) +
	fabs(y3_vol(i,1,k) - y3_vol(i,volume->r2_points,k)) +
	  fabs(z3_vol(i,1,k) - z3_vol(i,volume->r2_points,k));

  volume->r2_period = (dist <= 1.e-4 * box_size)? 1 : 0;

/* check for periodicity in direction 3 */
  dist = 0;
  for (j=1; j<=volume->r2_points; j++)
    for (i=1; i<=volume->r1_points; i++)
      dist += fabs(x3_vol(i,j,1) - x3_vol(i,j,volume->r3_points)) +
	fabs(y3_vol(i,j,1) - y3_vol(i,j,volume->r3_points)) +
	  fabs(z3_vol(i,j,1) - z3_vol(i,j,volume->r3_points));

  volume->r3_period = (dist <= 1.e-4 * box_size)? 1 : 0;

/* tmp */
  if (volume->r1_period) printf("The grid seems to be periodic in r1\n");
  if (volume->r2_period) printf("The grid seems to be periodic in r2\n");
  if (volume->r3_period) printf("The grid seems to be periodic in r3\n");
/* end tmp */

/* done */
  return (void *) volume;
}

static void
delete_volume_point( void *volume_data_ptr ){
  volume_point_info *info;
  
/* cast the pointer */
  info = (volume_point_info *) volume_data_ptr;

/* free the large arrays */
  free_volume_point_arrays( info );

/* free the structure */
  free(info);
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
volume_point_mapping(grid_point *gp, void *data_ptr ){
  volume_point_info *info;
  int i_loc, j_loc, k_loc;

/* cast the data pointer to the right type */
  info = (volume_point_info *) data_ptr;

/* get the nearest grid point below and to the left of (r,s) */
  i_loc = int_max( 1, 1 + (int) ( (info->n1-1)*gp->r ) );
  j_loc = int_max( 1, 1 + (int) ( (info->n2-1)*gp->s ) );
  k_loc = int_max( 1, 1 + (int) ( (info->n3-1)*gp->t ) );

/* prevent overflow */
  i_loc = int_min( i_loc, info->n1-1 );
  j_loc = int_min( j_loc, info->n2-1 );
  k_loc = int_min( k_loc, info->n3-1 );

/* evaluate the tri-linear mapping */
  gp->x  = LIN(x3_vol, gp->r, gp->s, gp->t);
  gp->y  = LIN(y3_vol, gp->r, gp->s, gp->t);
  gp->z  = LIN(z3_vol, gp->r, gp->s, gp->t);

  gp->xr = LIN_R(x3_vol, gp->r, gp->s, gp->t);
  gp->yr = LIN_R(y3_vol, gp->r, gp->s, gp->t);
  gp->zr = LIN_R(z3_vol, gp->r, gp->s, gp->t);

  gp->xs = LIN_S(x3_vol, gp->r, gp->s, gp->t);
  gp->ys = LIN_S(y3_vol, gp->r, gp->s, gp->t);
  gp->zs = LIN_S(z3_vol, gp->r, gp->s, gp->t);

  gp->xt = LIN_T(x3_vol, gp->r, gp->s, gp->t);
  gp->yt = LIN_T(y3_vol, gp->r, gp->s, gp->t);
  gp->zt = LIN_T(z3_vol, gp->r, gp->s, gp->t);

  return OK;
}
