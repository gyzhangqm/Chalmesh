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
#ifndef surface_h
#define surface_h

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>

#include <version.h>
#include <stupid_compiler.h>
#include <real.h>
#include <c_array.h>
#include <linked_list.h>
#include <nrutil.h>
#include <min_max.h>
#include <ogl_plot.h>
#include <input_output.h>
#include <linalg.h>

#include <window_defs.h>
#include <grid_point.h>

/* typedefs */

/* define new types for lists of surface_mappings */
typedef linked_list surface_mapping_list;
typedef linked_list_member surface_mapping_link;

/* data structure for reporting a hessian */
typedef struct hessian_data{
  real r, s, xrr, xrs, xss, yrr, yrs, yss, zrr, zrs, zss;
} hessian_data;

/* data structure for a surface mapping */
typedef struct surface_mapping{

  int used_by_volume;       /* number of volume grids that base their mapping */
                            /* on this surface mapping */
  char *name;               /* name of this component grid */
  int inverse_known;        /* inverse_known = 1 if the inverse of the mapping 
			       is defined in the grid function */

/* pointer to the function that defines the forward mapping */
  int (*forward_mapping)(grid_point *gp_ptr, struct surface_mapping *surface); 

/* pointer to the function that defines the inverse mapping 
   set to NULL if the inverse mapping is not known */
  int (*inverse_mapping)(grid_point *gp_ptr, struct surface_mapping *surface);

/* pointer to the function that interactively changes the default parameters */
  void (*set_surface)(input_output *io_ptr, struct surface_mapping *surface); 

/* pointer to the function that deletes the specific data */
  void (*delete_surface_data)(void *surface_data_ptr);

/* pointers to the functions for saving and recovering a surface from a binary file */
  void (*write_surface_data)( int fd, void *surface_data_ptr );

/* pointer to the function that defines the forward mapping */
  int (*hessian)(hessian_data *h_ptr, struct surface_mapping *surface); 

/* pointer to the data necessary to compute the mapping */
  void *surface_data_ptr;

/* number of non-fictitious points in the r and s directions, respectively */
  int r_points,s_points;    
  
  int r_period,s_period;    /* r_period and s_period refer to the r and 
			       s-directions, respectively. They mean:
			       period=1: periodic
			       period=0: not periodic */

  bounding_box *bb; /* bounding box */

  char *type;            
/* specific type of grid:
   1. Spherical patch.
   2. Interpolate bi-linearly between discrete grid points.
   3. Interpolate bi-cubically between discrete grid points.
*/

  int plot_mode;
  int color;

/* to get inside and outside right */
  int orientation;
} surface_mapping;


/* prototypes for the functions in this library */
surface_mapping *
new_surface_mapping( char *name);
surface_mapping *
delete_surface_mapping( surface_mapping *surface );
void
save_surface_mapping( surface_mapping *surface, FILE *fp );
int 
forward_surface_mapping( grid_point *gp_ptr, surface_mapping *surface );
int 
inverse_surface_mapping( grid_point *gp_ptr, surface_mapping *surface );
void 
set_surface( input_output *io_ptr, surface_mapping *surface );
void 
delete_specific_surface( surface_mapping *surface );
void 
write_specific_surface( int fd, surface_mapping *surface );
surface_mapping *
choose_surface( input_output *io_ptr, char *name );
void
draw_surface( surface_mapping *surface, int plot_mode );
int
surface_plot_mode(input_output *io_ptr, surface_mapping_list *surface_list,
		  surface_mapping *surface, int plot_mode);
void
surface_coordinates(surface_mapping *surface, real_array_2d **x_, 
		    real_array_2d **y_, real_array_2d **z_ );
surface_mapping_link *
get_surface( input_output *io_ptr, surface_mapping_list *head, int level);
int
compute_surface_normal(real *n_vec, grid_point *gp_ptr);
int
derivative_surface_normal(real *n_vec_r, real *n_vec_s, grid_point *gp_ptr, 
			  surface_mapping *surface);
int
surface_hessian( hessian_data *h_ptr, surface_mapping *surface );
void
surface_list_bb(bounding_box *bb, surface_mapping_list *surface_list);
void 
surface_bb( surface_mapping *surface );
void
check_surface(input_output *io_ptr, surface_mapping *surface);

/* end of surfgrid.h */
#endif
