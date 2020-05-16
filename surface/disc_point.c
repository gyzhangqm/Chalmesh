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

/* system prototype */
extern int sscanf( const char *line, const char *format, ...);

/* private member functions */
static int
disc_point_hessian( hessian_data *h_ptr, surface_mapping *surface );
static void
save_sub_grid( int i_min, int i_max, int j_min, int j_max, surface_mapping *surface, 
	      FILE *out_file);
static disc_point_info *
read_plot3d_file( FILE *fp, input_output *io_ptr );
/* end private member functions */

void
alloc_disc_point_arrays( disc_point_info *info ){
  info->x_ = create_real_array_2d(info->nr, info->ns);
  info->y_ = create_real_array_2d(info->nr, info->ns);
  info->z_ = create_real_array_2d(info->nr, info->ns);
}

void
free_disc_point_arrays( disc_point_info *info ){
  info->x_ = delete_real_array_2d(info->x_);
  info->y_ = delete_real_array_2d(info->y_);
  info->z_ = delete_real_array_2d(info->z_);
}

static disc_point_info *
read_plot3d_file( FILE *fp, input_output *io_ptr ){
  int nii, njj, nkk, nr, ns, i, j, left_handed, n_grids, which_grid=1, dum
    , read_past;
  real x_dum;
  char question[80], first_line[120];
  disc_point_info *info;
  
/* formatted plot3d file */
/* first try to determine if more than one surface is present in the file */
  read_past = 0;
  fgets( first_line, 120, fp );
  if (sscanf( first_line, "%i%i%i", &nii, &njj, &nkk ) == 1){
    n_grids = nii;
    if (n_grids == 1){
      which_grid = 1;
      fscanf( fp, "%i%i%i", &nii, &njj, &nkk );
    }
    else{
      sprintf(question, "There are %i surfaces in this file. Which one do "
	      "you want to read: ", n_grids);
      which_grid = get_int( io_ptr, question, 1, 0);
      which_grid = int_max(1, int_min( n_grids, which_grid) );
      for (i=1; i<which_grid; i++){
	fscanf( fp, "%i%i%i", &nii, &njj, &nkk );
	read_past += nii*njj*nkk;
      }
      fscanf( fp, "%i%i%i", &nii, &njj, &nkk );
      for (i=which_grid+1; i<=n_grids; i++)
	fscanf( fp, "%i%i%i", &dum, &dum, &dum );
    }
  }
/* tmp */
  printf("The dimensions were nii = %i, njj = %i, nkk = %i\n", 
	 nii, njj, nkk);
/* end tmp */
/* check if it is 3d */
  if( nii > 1 && njj > 1 && nkk > 1 ){
    printf("The input file contains a truely 3-d grid, but a surface grid was "
	   "expected\n");
    fclose( fp );
    return ERROR;
  }
  else if( nkk == 1 ){          
    nr = nii;
    ns = njj;
    left_handed = 0;
  } 
  else if( njj == 1 ){
    nr = nii;
    ns = nkk;
    left_handed = 1;
  }
  else if( nkk == 1 ){
    nr = njj;
    ns = nkk;
    left_handed = 0;
  }
  else{
    printf("The dimensions in the file make no sense. nii = %i, njj = %i, "
	   "nkk = %i.\n", nii, njj, nkk);
    fclose( fp );
    return ERROR;
  }

/* read past the first which_grid-1 surfaces */  
  for (i=1; i<=3*read_past; i++){
    if ((fscanf( fp, "%lg", &x_dum ) ) == EOF){
	  printf("Input error in read_plot3d while reading past the "
		 "first %i surfaces.\n", which_grid - 1);
	  fclose( fp );
	  return NULL;
	}
  }


  info = (disc_point_info *) malloc( sizeof(disc_point_info) );

/* set the dimensions */
  info->nr = nr;
  info->ns = ns;

/* the grid step is used by local_map */
  info->dr = 1.0/((real) nr-1);
  info->ds = 1.0/((real) ns-1);

/* allocate x, y, z arrays */
  alloc_disc_point_arrays( info );

/* read the x coordinates */
  for (j = 1; j <= info->ns; j++)
    if (left_handed)
      for (i = info->nr; i >= 1; i--){
	if ((fscanf( fp, "%lg", &(x_surf(i,j)) )) == EOF){
	  printf("Input error in read_plot3d while reading x-coordinates at "
		 "i = %i, j= %i.\n", info->nr-i+1, j);
	  free_disc_point_arrays( info );
	  fclose( fp );
	  return NULL;
	}
      }
    else
      for (i = 1; i <= info->nr; i++){
	if ((fscanf( fp, "%lg", &(x_surf(i,j)) )) == EOF){
	  printf("Input error in read_plot3d while reading x-coordinates at "
		 "i = %i, j= %i.\n", i, j);
	  free_disc_point_arrays( info );
	  fclose( fp );
	  return NULL;
	}
      }  

/* read the y coordinates */
  for (j = 1; j <= info->ns; j++)
    if (left_handed)
      for (i = info->nr; i >= 1; i--){
	if ((fscanf( fp, "%lg", &(y_surf(i,j)) )) == EOF){
	  printf("Input error in read_plot3d while reading y-coordinates at "
		 "i = %i, j= %i.\n", info->nr-i+1, j);
	  free_disc_point_arrays( info );
	  fclose( fp );
	  return NULL;
	}
      }
    else
      for (i = 1; i <= info->nr; i++){
	if ((fscanf( fp, "%lg", &(y_surf(i,j)) )) == EOF){
	  printf("Input error in read_plot3d while reading y-coordinates at "
		 "i = %i, j= %i.\n", i, j);
	  free_disc_point_arrays( info );
	  fclose( fp );
	  return NULL;
	}
      }

/* read the z coordinates */
  for (j = 1; j <= info->ns; j++)
    if (left_handed)
      for (i = info->nr; i >= 1; i--){
	if ((fscanf( fp, "%lg", &(z_surf(i,j)) )) == EOF){
	  printf("Input error in read_plot3d while reading z-coordinates at "
		 "i = %i, j= %i.\n", info->nr-i+1, j);
	  free_disc_point_arrays( info );
	  fclose( fp );
	  return NULL;
	}
      }
    else
      for (i = 1; i <= info->nr; i++){
	if ((fscanf( fp, "%lg", &(z_surf(i,j)) )) == EOF){
	  printf("Input error in read_plot3d while reading z-coordinates at "
		 "i = %i, j= %i.\n", i, j);
	  free_disc_point_arrays( info );
	  fclose( fp );
	  return NULL;
	}
      }
  
/* successful completion */
  fclose( fp );

  return info;
}


void 
set_disc_point(input_output *io_ptr, surface_mapping *surface){
  char prompt[80], *file_name;
  int icom, quit=0, replot=1, tmp, i_min, i_max, j_min, j_max;
  const int save_command=1;
  FILE *out_file;
  disc_point_info *info;

#include "set_disc_point_com.h"

  info = (disc_point_info *) surface->surface_data_ptr;
  sprintf(prompt, "%s: discrete point mapping>", surface->name);

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

    switch( get_command( io_ptr, prompt, COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case ORIENTATION:
      surface->orientation = !surface->orientation;
/* swap grid point and periodicity information */
      tmp = surface->r_points;
      surface->r_points = surface->s_points; surface->s_points = tmp;
      tmp = surface->r_period;
      surface->r_period = surface->s_period; surface->s_period = tmp;
      replot = 1;
      break;

/* change the plot mode */
    case PLOT_MODE:
      surface->plot_mode = surface_plot_mode(io_ptr, NULL, surface, 
						 surface->plot_mode);
/* don't need to redraw the surface grids! */
      break;

    case SAVE_SUBGRID:
      if ( (out_file = open_ascii_file( io_ptr, "Enter plot3d output file: ", 
				       "test.plot3d", &file_name, 'w', 0, 
				       save_command)) != NULL){
	i_min = int_max(1, get_int( io_ptr, "Starting index in r:", 1, LEVEL+1));
	i_max = int_min(surface->r_points, get_int( io_ptr, "Ending index in r:", 
						   surface->r_points, LEVEL+1));
	j_min = int_max(1, get_int( io_ptr, "Starting index in s:", 1, LEVEL+1));
	j_max = int_min(surface->s_points, get_int( io_ptr, "Ending index in s:", 
						   surface->s_points, LEVEL+1));
	if (i_max-i_min > 0 && j_max-j_min > 0)
	  save_sub_grid( i_min, i_max, j_min, j_max, surface, out_file);
	else
	  printf("Sorry, can not save the sub grid because the range in i or j is "
		 "empty\n");
	fclose( out_file );
      }
      else
	printf("Sorry, could not open the file %s with write permission.\n", file_name);
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

    case SHOW:
      printf("Number of grid points in r = %i, and in s = %i.\n", surface->r_points, 
	     surface->s_points);
      printf("This surface is based on %i discrete points in r and %i discrete "
	     "points in s.\n", info->nr, info->ns);
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
save_sub_grid( int i_min, int i_max, int j_min, int j_max, surface_mapping *surface, 
	      FILE *out_file){
  int i, j;
  disc_point_info *info;

  info = (disc_point_info *) surface->surface_data_ptr;
/* header */
  fprintf( out_file, "%i %i %i\n", i_max-i_min+1, j_max-j_min+1, 1);
  if (surface->orientation == 0){
/* write the x-coordinates */
    for (j=j_min; j<=j_max; j++)
      for (i=i_min; i<=i_max; i++){
	fprintf( out_file, "%e\n", x_surf(i,j));
      }
/* write the y-coordinates */
    for (j=j_min; j<=j_max; j++)
      for (i=i_min; i<=i_max; i++){
	fprintf( out_file, "%e\n", y_surf(i,j));
      }
/* write the z-coordinates */
    for (j=j_min; j<=j_max; j++)
      for (i=i_min; i<=i_max; i++){
	fprintf( out_file, "%e\n", z_surf(i,j));
      }
  } /* end if orientation == 0 */
  else{
/* write the x-coordinates */
    for (j=j_min; j<=j_max; j++)
      for (i=i_min; i<=i_max; i++){
	fprintf( out_file, "%e\n", x_surf(j,i));
      }
/* write the y-coordinates */
    for (j=j_min; j<=j_max; j++)
      for (i=i_min; i<=i_max; i++){
	fprintf( out_file, "%e\n", y_surf(j,i));
      }
/* write the z-coordinates */
    for (j=j_min; j<=j_max; j++)
      for (i=i_min; i<=i_max; i++){
	fprintf( out_file, "%e\n", z_surf(j,i));
      }
  } /* end if orientation == 1 */
}

void *
init_disc_point(input_output *io_ptr, 
		surface_mapping *surface ){
  FILE *fp;
  char *file_name, *prompt;
  int icom, quit=0, i, j;
  const int save_command=1;
  real n_vec[3], box_size, dist;
  grid_point gp;
  disc_point_info *info=NULL;

#include "init_disc_point_com.h"

  prompt = "select file format for reading grid points>";

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
  surface->type = "Discrete point"; 
  surface->surface_data_ptr = (void *) info; 
  surface->forward_mapping = disc_point_mapping; 
  surface->hessian = disc_point_hessian;
  surface->inverse_mapping = NULL; 
  surface->set_surface = set_disc_point; 
  surface->delete_surface_data = delete_disc_point;

  surface->r_points = info->nr;
  surface->s_points = info->ns;

/* check for singularities */
  for (i=1; i<= info->nr; i++){
    for (j=1; j<= info->ns; j++){
      gp.r = (i-1)*info->dr;
      gp.s = (j-1)*info->ds;
      forward_surface_mapping( &gp, surface );
      if (!compute_surface_normal( n_vec, &gp ))
	printf("Warning: The jacobian is singular at (r, s) = (%f, %f)\n", gp.r, gp.s);
    }
  }

/* compute bounding box */
  surface_bb( surface );
  box_size = (surface->bb->x_max - surface->bb->x_min +
	      surface->bb->y_max - surface->bb->y_min +
	      surface->bb->z_max - surface->bb->z_min)/3.0;

/* check for periodicity in the r-direction */
  dist = 0;
  for (j=1; j<=surface->s_points; j++)
    dist += fabs(x_surf(1,j) - x_surf(surface->r_points,j)) +
	fabs(y_surf(1,j) - y_surf(surface->r_points,j)) +
	  fabs(z_surf(1,j) - z_surf(surface->r_points,j));

  surface->r_period = (dist <= 1.e-4 * box_size)? 1 : 0;

/* check for periodicity in the s-direction */
  dist = 0;
    for (i=1; i<=surface->r_points; i++)
      dist += fabs(x_surf(i,1) - x_surf(i,surface->s_points)) +
	fabs(y_surf(i,1) - y_surf(i,surface->s_points)) +
	  fabs(z_surf(i,1) - z_surf(i,surface->s_points));

  surface->s_period = (dist <= 1.e-4 * box_size)? 1 : 0;

/* tmp */
  if (surface->r_period) printf("The grid seems to be periodic in r\n");
  if (surface->s_period) printf("The grid seems to be periodic in s\n");
/* end tmp */

/* done */
  return (void *) surface;
}

void
delete_disc_point( void *tuut ){
  disc_point_info *info;
/* cast the data pointer to the right type */
  info = (disc_point_info *) tuut;

/* free the arrays */
  free_disc_point_arrays( info );
/* free everything else */
  free(info);
}


int 
disc_point_mapping(grid_point *gp_ptr, surface_mapping *surface){
  disc_point_info *info;
  int i_loc, j_loc, i_close, j_close, j_m, j_p;
  const int iw_low=1, iw_high=2, iw1=3;
  real r, s;
  const real eps = 1.e-3;;

/*#include "linear-interp.h"*/
#include "cubic-interp.h"

/* cast the pointer to the right type */
  info = (disc_point_info *) surface->surface_data_ptr;

/* swap directions ? */
  if (surface->orientation){
    r = gp_ptr->s;
    s = gp_ptr->r;
  }
  else{
    r = gp_ptr->r;
    s = gp_ptr->s;
  }

  if (surface->r_period)
    {
      if (r < 0.0) r += 1.0;
      else if (r > 1.0) r -= 1.0;
    }

  if (surface->s_period)
    {
      if (s < 0.0) s += 1.0;
      else if (s > 1.0) s -= 1.0;
    }

/* bi-cubic interpolation */
/* get the nearest grid point below and to the left of (r,s) */
  i_close = 1 + (int) ( (info->nr-1)*r );
  j_close = 1 + (int) ( (info->ns-1)*s );

/* find the interpolation location */
  if (i_close-iw_low < 1)
    i_loc = 1;
  else if (i_close+iw_high > info->nr)
    i_loc = info->nr-iw1;
  else
    i_loc = i_close-iw_low;

  if (j_close-iw_low < 1)
    {
      if (surface->s_period)
	{
	  j_m = info->ns-1;
	  j_loc = 0;
	  j_p = 3;
	}
      else
	{
	  j_m = j_loc = 1;
	  j_p = 4;
	}
    }
  else if (j_close+iw_high > info->ns)
    {
      if (surface->s_period)
	{
	  j_p = 2;
	  j_m = j_loc = info->ns-2;
	}
      else
	{
	  j_m = j_loc = info->ns-iw1;
	  j_p = info->ns;
	}
    }
  else
    {
      j_m = j_loc = j_close-iw_low;
      j_p = j_m + 3;
    }

  /* check alpha and beta */
/*   printf("r: %e\n", r); */
/*   printf("D2_ALPHA_0(r) = %e, (D_ALPHA_0(r+eps) - D_ALPHA_0(r-eps))/(2eps) = %e\n",  */
/* 	 D2_ALPHA_0(r), (D_ALPHA_0(r+eps) - D_ALPHA_0(r-eps))/(2.*eps)); */
/*   printf("D2_ALPHA_1(r) = %e, (D_ALPHA_1(r+eps) - D_ALPHA_1(r-eps))/(2eps) = %e\n",  */
/* 	 D2_ALPHA_1(r), (D_ALPHA_1(r+eps) - D_ALPHA_1(r-eps))/(2.*eps)); */
/*   printf("D2_ALPHA_2(r) = %e, (D_ALPHA_2(r+eps) - D_ALPHA_2(r-eps))/(2eps) = %e\n",  */
/* 	 D2_ALPHA_2(r), (D_ALPHA_2(r+eps) - D_ALPHA_2(r-eps))/(2.*eps)); */
/*   printf("D2_ALPHA_3(r) = %e, (D_ALPHA_3(r+eps) - D_ALPHA_3(r-eps))/(2eps) = %e\n",  */
/* 	 D2_ALPHA_3(r), (D_ALPHA_3(r+eps) - D_ALPHA_3(r-eps))/(2.*eps)); */

/*   printf("s: %e\n", s); */
/*   printf("D2_BETA_0(s) = %e, (D_BETA_0(s+eps) - D_BETA_0(s-eps))/(2eps) = %e\n",  */
/* 	 D2_BETA_0(s), (D_BETA_0(s+eps) - D_BETA_0(s-eps))/(2.*eps)); */
/*   printf("D2_BETA_1(s) = %e, (D_BETA_1(s+eps) - D_BETA_1(s-eps))/(2eps) = %e\n",  */
/* 	 D2_BETA_1(s), (D_BETA_1(s+eps) - D_BETA_1(s-eps))/(2.*eps)); */
/*   printf("D2_BETA_2(s) = %e, (D_BETA_2(s+eps) - D_BETA_2(s-eps))/(2eps) = %e\n",  */
/* 	 D2_BETA_2(s), (D_BETA_2(s+eps) - D_BETA_2(s-eps))/(2.*eps)); */
/*   printf("D2_BETA_3(s) = %e, (D_BETA_3(s+eps) - D_BETA_3(s-eps))/(2eps) = %e\n",  */
/* 	 D2_BETA_3(s), (D_BETA_3(s+eps) - D_BETA_3(s-eps))/(2.*eps)); */

/* evaluate the bi-cubic transformation */
  gp_ptr->x  = CUB(r, s, x_surf);
  gp_ptr->y  = CUB(r, s, y_surf);
  gp_ptr->z  = CUB(r, s, z_surf);

  if (surface->orientation){
    gp_ptr->xr = CUB_S(r, s, x_surf);
    gp_ptr->yr = CUB_S(r, s, y_surf);
    gp_ptr->zr = CUB_S(r, s, z_surf);

    gp_ptr->xs = CUB_R(r, s, x_surf);
    gp_ptr->ys = CUB_R(r, s, y_surf);
    gp_ptr->zs = CUB_R(r, s, z_surf);
  }
  else{
    gp_ptr->xr = CUB_R(r, s, x_surf);
    gp_ptr->yr = CUB_R(r, s, y_surf);
    gp_ptr->zr = CUB_R(r, s, z_surf);

    gp_ptr->xs = CUB_S(r, s, x_surf);
    gp_ptr->ys = CUB_S(r, s, y_surf);
    gp_ptr->zs = CUB_S(r, s, z_surf);
  }

/* bi-linear interpolation */
/* get the nearest grid point below and to the left of (r,s) */
/*   i_loc = int_max( 1, 1 + (int) ( (info->nr-1)*r ) ); */
/*   j_loc = int_max( 1, 1 + (int) ( (info->ns-1)*s ) ); */

/* prevent overflow */
/*   i_loc = int_min( i_loc, info->nr-1 ); */
/*   j_loc = int_min( j_loc, info->ns-1 ); */

/* evaluate the bi-linear transformation */
/*   gp_ptr->x  = LIN(r, s, x_surf); */
/*   gp_ptr->y  = LIN(r, s, y_surf); */
/*   gp_ptr->z  = LIN(r, s, z_surf); */

/*   if (surface->orientation){ */
/*     gp_ptr->xr = LIN_S(r, s, x_surf); */
/*     gp_ptr->yr = LIN_S(r, s, y_surf); */
/*     gp_ptr->zr = LIN_S(r, s, z_surf); */

/*     gp_ptr->xs = LIN_R(r, s, x_surf); */
/*     gp_ptr->ys = LIN_R(r, s, y_surf); */
/*     gp_ptr->zs = LIN_R(r, s, z_surf); */
/*   } */
/*   else{ */
/*     gp_ptr->xr = LIN_R(r, s, x_surf); */
/*     gp_ptr->yr = LIN_R(r, s, y_surf); */
/*     gp_ptr->zr = LIN_R(r, s, z_surf); */

/*     gp_ptr->xs = LIN_S(r, s, x_surf); */
/*     gp_ptr->ys = LIN_S(r, s, y_surf); */
/*     gp_ptr->zs = LIN_S(r, s, z_surf); */
/*   } */

  return OK;
}

static int
disc_point_hessian( hessian_data *h_ptr, surface_mapping *surface ){
  disc_point_info *info;
  const int iw_low=1, iw_high=2, iw1=3;
  int i_loc, j_loc, i_close, j_close, j_m, j_p;
  real r, s;

/* cast the pointer to the right type */
  info = (disc_point_info *) surface->surface_data_ptr;

/* swap directions ? */
  if (surface->orientation){
    r = h_ptr->s;
    s = h_ptr->r;
  }
  else{
    r = h_ptr->r;
    s = h_ptr->s;
  }

  if (surface->r_period)
    {
      if (r < 0.0) r += 1.0;
      else if (r > 1.0) r -= 1.0;
    }

  if (surface->s_period)
    {
      if (s < 0.0) s += 1.0;
      else if (s > 1.0) s -= 1.0;
    }

/* bi-cubic interpolation */
/* get the nearest grid point below and to the left of (r,s) */
  i_close = 1 + (int) ( (info->nr-1)*r );
  j_close = 1 + (int) ( (info->ns-1)*s );

/* find the interpolation location */
  if (i_close-iw_low < 1)
    i_loc = 1;
  else if (i_close+iw_high > info->nr)
    i_loc = info->nr-iw1;
  else
    i_loc = i_close-iw_low;

  if (j_close-iw_low < 1)
    {
      if (surface->s_period)
	{
	  j_m = info->ns-1;
	  j_loc = 0;
	  j_p = 3;
	}
      else
	{
	  j_m = j_loc = 1;
	  j_p = 4;
	}
    }
  else if (j_close+iw_high > info->ns)
    {
      if (surface->s_period)
	{
	  j_p = 2;
	  j_m = j_loc = info->ns-2;
	}
      else
	{
	  j_m = j_loc = info->ns-iw1;
	  j_p = info->ns;
	}
    }
  else
    {
      j_m = j_loc = j_close-iw_low;
      j_p = j_m + 3;
    }

  if (surface->orientation){
    h_ptr->xss = CUB_RR(r, s, x_surf);
    h_ptr->yss = CUB_RR(r, s, y_surf);
    h_ptr->zss = CUB_RR(r, s, z_surf);

    h_ptr->xrs = CUB_RS(r, s, x_surf);
    h_ptr->yrs = CUB_RS(r, s, y_surf);
    h_ptr->zrs = CUB_RS(r, s, z_surf);

    h_ptr->xrr = CUB_SS(r, s, x_surf);
    h_ptr->yrr = CUB_SS(r, s, y_surf);
    h_ptr->zrr = CUB_SS(r, s, z_surf);
  }
  else{
    h_ptr->xrr = CUB_RR(r, s, x_surf);
    h_ptr->yrr = CUB_RR(r, s, y_surf);
    h_ptr->zrr = CUB_RR(r, s, z_surf);

    h_ptr->xrs = CUB_RS(r, s, x_surf);
    h_ptr->yrs = CUB_RS(r, s, y_surf);
    h_ptr->zrs = CUB_RS(r, s, z_surf);

    h_ptr->xss = CUB_SS(r, s, x_surf);
    h_ptr->yss = CUB_SS(r, s, y_surf);
    h_ptr->zss = CUB_SS(r, s, z_surf);
  }

  return OK;
}
