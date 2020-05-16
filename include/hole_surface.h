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
#ifndef hole_surface_h
#define hole_surface_h

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

/* data structure for returning an inverse */
typedef struct inverse_point{

/* the indices of the interpolation point in this grid */
  int i_point, j_point;
/* the pointer to the other surface grid (for a hole_surface) */
  struct discrete_surface *sgrid_loc;
/* r and s coordinates for the interpolation point in the other grids
   coordinate system. n_loc is the (small) normal distance to the discrete surface. */
  real r_loc, s_loc, n_loc;
/* the indices of the lower left gridpoint in the enclosing grid cell in the */
/* other grid. */
  int i_loc, j_loc;
} inverse_point;

/* data structure for describing a boundary of a component grid */
typedef struct subdomain_info{
  struct subdomain_info *low_left, *low_right, *upp_left, *upp_right;
  int i_min, i_max, j_min, j_max;
  real x_min, x_max, y_min, y_max, z_min, z_max;
} subdomain_info;

typedef struct triangle{
  real x[3], y[3], z[3];
  int i[3], j[3];
  struct discrete_surface *sgrid[3];
} triangle;

/* data structure for a hole-cutting curve */
typedef struct hole_cutting_curve
{
  int n_points;
  real *x_edge, *y_edge, *z_edge;
  real *inside_point[2][3];
  real *inside_normal[2][3];
  int periodic;
  int edge;
  real max_dist;
} hole_cutting_curve;

/* define new types for lists of hole_surfaces and discrete_surfaces */
typedef linked_list terminator_list;
typedef linked_list_member terminator_link;
typedef linked_list hole_cutting_curve_list;
typedef linked_list_member hole_cutting_curve_link;
typedef linked_list discrete_surface_list;
typedef linked_list_member discrete_surface_link;
typedef linked_list hole_surface_list;
typedef linked_list_member hole_surface_link;
typedef linked_list tri_array_node_list;
typedef linked_list_member tri_array_node_link;

/* data structure for a patched surface */
typedef struct hole_surface{

/* this surface grid corresponds to faces of 3d-grids with the following */
/* surface label */
  int surface_label;

/* is the computational domain inside or outside of this surface? */
  int inside; 

/* pointer to the linked list of discrete_surfaces and corresponding surface mappings */
  discrete_surface_list *grid_list;

/* list of hole_cutting_curves */
  hole_cutting_curve_list *edge_curves;

/* bounding box */
  bounding_box *bb;

/* boundary curve missmatch tolerance */
  real max_b_dist;

  char *name;

  int plot_mode;

} hole_surface;
  

/* data structure for a component sgrid */
typedef struct discrete_surface{

  int priority;             /* priority of this grid in the hole_surface */
  char *name;               /* name of the discrete_surface */

/* dimensions of arrays in the r and s directions, respectively */
  int r_dim,s_dim;
  
  int surface_label;
  struct component_vgrid *vgrid_mother;
  int side, dir;

/* grid sizes in the r and s-directions, */
/* respectively. r_step = 1/(r_points-1) */
  real r_step,s_step; 

  int r_period,s_period;    /* r_period and s_period refer to the r and 
			       s-directions, respectively. They mean:
			       period=1: periodic
			       period=0: not periodic */


  int_array_2d *curve_ptr;  /* curve(is,id): side is=1,2, direction 
			       id=1,2. curve(is,id)= a number identifying 
			       each side of the grid with a boundary curve. 
			       Used to identify sides that correspond to edges
			       that join faces with surface_label != 0 */

  int_array_2d *hole_curve_; /* hole_curve(is,id): side is=1,2, direction 
				id=1,2. If hole_curve(is,id) != 0, the side 
				will be used to cut holes in other discrete 
				surfaces in the same hole_surface */

  int_array_2d *flag_ptr; /* only used by the patching algorithm */
/* flag(i,j) is the flag of gridpoint i,j. */
/* flag(i,j)=grid_number: interior point,  */
/* flag(i,j)=-k: interpolate from grid k   */
/* flag(i,j)=0: unused point.              */

  real_array_2d *x_ptr, *y_ptr, *z_ptr;/* only used by the patching algorithm */
/* x(i,j), y(i,j), z(i,j) are the x, y, z-coordinates, */
/* respectively of grid point i,j.                     */

/* bounding box */
  bounding_box *bb; 

/* domain subdivision information */
  subdomain_info *quad_tree;

/* surface smoothness estimate */
  real max_n_dist, max_angle, max_t_dist;

  int plot_mode, plot_it;

} discrete_surface;

/* define some macros for computing indices of multi-dimensional arrays */
#define curve(is,id) compute_index_2d(sgrid_ptr->curve_ptr,is,id)
#define hole_curve(is,id) compute_index_2d(sgrid_ptr->hole_curve_,is,id)

#define flag(i,j)  compute_index_2d(sgrid_ptr->flag_ptr,i,j)

#define x(i,j)     compute_index_2d(sgrid_ptr->x_ptr,i,j)
#define y(i,j)     compute_index_2d(sgrid_ptr->y_ptr,i,j)
#define z(i,j)     compute_index_2d(sgrid_ptr->z_ptr,i,j)

/* prototypes for the functions in this library */
void
merge_hole_curves(hole_surface *overg);
void
cut_surface_holes(input_output *io_ptr, hole_surface *overg);
hole_cutting_curve *
delete_hole_curve(hole_cutting_curve *hole_curve);
hole_cutting_curve *
new_hole_curve(real *x_edge, real *y_edge, real *z_edge, int n1, 
	       real *inside_point[2][3], real *inside_normal[2][3], int edge);
void
insert_hole_curve(hole_surface *patch, hole_cutting_curve *hole_curve );
hole_surface *
new_hole_surface( char *name, int surface_label );
hole_surface *
delete_hole_surface( hole_surface *patch );
discrete_surface *
new_discrete_surface(char *name, real_array_2d *x_, 
		     real_array_2d *y_, real_array_2d *z_, int present_surface, 
		     struct component_vgrid *mother);
discrete_surface *
delete_discrete_surface( discrete_surface *sgrid_ptr );
void
insert_discrete_surface(hole_surface *patch, discrete_surface *sgrid_ptr);

int
hole_plot_mode(input_output *io_ptr, hole_surface *patch, int plot_mode);
void
draw_hole_surface( hole_surface *patch, int plot_mode );

void 
bi_linear(grid_point *gp_ptr, discrete_surface *sgrid_ptr);
int
inside_surface(hole_surface *patch, real xp[], int present_surface, int present_curve, int debug);
int
triangulate(input_output *io_ptr, hole_surface *patch);
hole_surface_link *
get_hole_surface( input_output *io_ptr, hole_surface_list *head);
discrete_surface_link *
get_discrete_surface( input_output *io_ptr, discrete_surface_list *head);
int
newton_search(real xp, real yp, real zp, real *r_loc, real *s_loc, real *n_loc, 
	      discrete_surface *project_onto);
inverse_point *
search_quad_tree(real xp, real yp, real zp, subdomain_info *info, 
		 discrete_surface *sgrid_ptr, real max_b_dist,
		 int query_curve, int always_use_max_b_dist );




/* end of overlap.h */
#endif
