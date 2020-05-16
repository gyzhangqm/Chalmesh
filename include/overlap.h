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
#ifndef overlap_h
#define overlap_h

#include <window_defs.h>
#include <grid_point.h>
#include <surface.h>
#include <hole_surface.h>
#include <volume.h>

/* define new types for lists of surface_grids */
typedef linked_list component_vgrid_list;
typedef linked_list_member component_vgrid_link;
typedef linked_list overlapping_3d_list;
typedef linked_list_member overlapping_3d_link;

enum {BAD_DISC, DEAD_INTERP_LOC, NON_EXPLICIT_LOC, BAD_INTERP_DEAD_LOC
	, BAD_INTERP_NON_EXPLICIT};

/* data structure for an 3-d interpolation point */
typedef struct interp_3d_point{

/* pointers to the previous interpolation point 
   for the first point in the list, prev == NULL */
  struct interp_3d_point *prev;

/* a flag to indicate if the interpolation point is still in use */
  int active;

/* this flag indicates wether the weights have been corrected for boundary mismatch */
  int corrected[3];

/* the indices of the interpolation point in this grid */
  int i_point, j_point, k_point;
/* the pointer to the other surface grid (for overlapping sgrids) */
  struct component_vgrid *vgrid_loc;
/* r and s coordinates for the interpolation point in the other grids
   coordinate system. n_loc is the (small) normal distance to the discrete surface. */
  real r_loc, s_loc, t_loc;
/* the indices of the lower left gridpoint in the other grid. The */
/*   interpolation stencil will use the points */
/*   i_loc <= i <= i_loc + interp_width - 1, */
/*   j_loc <= j <= j_loc + interp_width - 1, */
/*   k_loc <= k <= k_loc + interp_width - 1, */
  int i_loc, j_loc, k_loc;
} interp_3d_point;

/* data structure for the oct-search-tree */
typedef struct oct_tree_info{
  int i_min, i_max, j_min, j_max, k_min, k_max;
  real x_min, x_max, y_min, y_max, z_min, z_max;
  struct oct_tree_info **children;
} oct_tree_info;


/* data structure for a bad grid point */
typedef struct bad_3d_point{
  struct bad_3d_point *prev;
  int i, j, k;
  int type; 
/*  the value of type can be one of the following:  */
/* (one day, replace the color coding with bitmap symbols) */
/*    BAD_DISC: Bad discretization point,  */
/*    DEAD_INTERP_LOC: Dead interpolation location, */
/*    NON_EXPLICIT_LOC: Non-explicit interpolation location, */
/*    BAD_INTERP_DEAD_LOC: Bad interpolation point because of dead */
/*                         interpolation location, */
/*    BAD_INTERP_NON_EXPLICIT: Bad interpolation point because of */
/*                             non-explicit interpolation location. */
} bad_3d_point;


/* data structure for a component in the 3d overlapping grid */
typedef struct component_vgrid{

  char *name;               /* name of this component grid */

  int grid_type; /* 1 for Cartesian, otherwise curvilinear */

  int priority;

  int no_hole; /* 1 means that the hole-cutting algorithm will not be performed */

  int orientation; /* =1 for a right-handed system, =-1 for a left-handed system */

/* the correspondence between the directions in the mother grid and the component */
/* surface grid are as follows */
/* dir  r_dim   s_dim  */
/*  1   r2_dim  r3_dim */
/*  2   r3_dim  r1_dim */
/*  3   r1_dim  r2_dim */

/* dimensions of arrays in the r1, r2, and r3 directions */
  int r1_dim, r2_dim, r3_dim;
  
/*  Hence, if r_period == 0: r_dim = r_points + 2*extra 
    and if    r_period == 1: r_dim = r_points + 2*extra_period - 1
    See the description of range below for more details.
*/
  int r1_period, r2_period, r3_period; /* Periodicity flag meaning:
			       period=1: periodic
			       period=0: not periodic */

  int_array_2d *range_ptr;  /* range of non-fictitious points in the arrays
			       range(is,id) refers to side is=1,2 and direction
			       id=1,3. */
/* If r_period == 0, s_period == 0:
   range(1,1)=extra+1, 
   range(2,1)=r_points+extra, 
   range(1,2)=extra+1,
   range(2,2)=s_points+extra. 

   If r_period == 1:
   range(1,1) = extra_period + 1,
   range(2,1) = r_points + extra_period,

   If s_period == 1:
   range(1,2) = extra_period + 1,
   range(2,2) = s_points + extra_period,
*/
  int_array_2d *bc_ptr;     /* bc3(is,id): boundary condition on side is=1,2
			       in direction id=1,2,3 */
  int_array_2d *surf_ptr;   /* surf3(is,id): side is=1,2, direction 
			       id=1,2,3. surface(is,id)= a number identifying 
			       each side of the volume grid with a boundary surface. 
			       Used when sides of two or more grids share the
			       same physical boundary surface */
  int_array_1d *edge_curve_; /* edge_curve(i) corresponds to grid edges which join */
  /* grid faces with surface_label != 0. The index `i' is related to the */
  /* intersection line between the faces as follows: */
  /* i   edge_curve(i) */
  /* 1   r1=0, r2=0    */
  /* 2   r1=1, r2=0    */
  /* 3   r1=0, r2=1    */
  /* 4   r1=1, r2=1    */
  /* 5   r1=0, r3=0    */
  /* 6   r1=1, r3=0    */
  /* 7   r1=0, r3=1    */
  /* 8   r1=1, r3=1    */
  /* 9   r2=0, r3=0    */
  /* 10  r2=1, r3=0    */
  /* 11  r2=0, r3=1    */
  /* 12  r2=1, r3=1    */

/* grid sizes in the ri-directions, i=1,2,3.*/
/* r1_step = 1/(r1_points-1) etc. */
  real r1_step, r2_step, r3_step; /* this information is redundant and should 
				     only be used by the overlapping algorithm */

  int_array_3d *flag_ptr; /* only used by the overlapping algorithm */
/* flag(i,j,k) is the flag of gridpoint i,j,k. */
/* flag(i,j,k)=grid_number: interior point,  */
/* flag(i,j,k)=-q: interpolate from the grid with priority q, */
/* flag(i,j,k)=0: unused point.              */

  real_array_3d *x_ptr, *y_ptr, *z_ptr;/* only used by the overlapping algorithm */
/* x(i,j,k), y(i,j,k), z(i,j,k) are the x, y, z-coordinates, */
/* respectively of grid point i,j,k.                     */

  int n_interp;             /* number of interpolation points in this 
			       component grid */
  interp_3d_point *last_interp; /* pointer to the last interpolation point. 
				We save the interpolation points on a stack. */

  bad_3d_point *last_bad_point; /* pointer to the stack of bad grid points */

/* domain subdivision information */
  oct_tree_info *oct_tree;

/* hole_surface for searching for intersections with the ray */
  hole_surface *bounding_surface;

/* talkativity */
  int verbose;

/* surface missfit tolerance */
  real max_s_dist;

/* edge missfit tolerance */
  real max_e_dist;

  int analytical_inverse; /* TRUE if an analytical inverse function is known, */
                          /* otherwise FALSE */

/* pointer to the forward mapping function */
  int (*forward_mapping)(grid_point *gp_ptr, void *vgrid_data_ptr);
/* pointer to the inverse function. NULL if the inverse mapping is not known */
  int (*inverse_mapping)(grid_point *gp_ptr, void *vgrid_data_ptr);
/* data to evaluate the mapping */
  void *vgrid_data_ptr;

/* bounding box */
  bounding_box *bb; 

/* length scale (used for scaling tolerances in the Newton iteration) */
  real length_scale;

  int plot_mode;
  int color;
  int plot_it;

/* rotation matrix and translation vectors */
  real_array_2d *rotation_;
  real translation[3], pre_trans[3];

/* flag that indicates wether the mapping is translated and/or rotated */
  int transformed;

} component_vgrid;

#define rot_vg(i,j) compute_index_2d(vgrid->rotation_,i,j)

/* the overlapping 3-d grid structure */
typedef struct overlapping_3d_grid{

/* is the overlap valid? */
  int valid;

/* pointer to the linked list of discrete_surfaces and corresponding surface mappings */
  component_vgrid_list *grid_list;
  volume_mapping_list *mapping_list;

/* list of hole-surfaces corresponding to all component grid faces */
/* with matching surface labels > 0 */
  hole_surface_list *surface_grid_list;

/* the correspondence between the directions in the mother grid and the component */
/* surface grid are as follows */
/* dir  r_dim   s_dim  */
/*  1   r2_dim  r3_dim */
/*  2   r3_dim  r1_dim */
/*  3   r1_dim  r2_dim */

  int extra;          /* number of fictitious grid points outside each boundary.*/
  int extra_period;   /* number of extra grid points at a periodic boundary. */
/* 
   Observe that extra_period grid points are added at the lower bound but
   extra_period-1 grid points are added at the upper bound. See the above 
   descriptions of range and dim for details.
*/
  int disc_width;     /* width of interior discretisation stencil. Must be odd. */
  int normal_width;   /* width of boundary discretization stencil in the 
			 direction normal to the boundary */
  int tangent_width;  /* width of boundary discretization stencil in the 
			 direction tangential to the boundary */
  int corner_width;   /* width of the stencil close to a corner */
  int interp_width;   /* width of the interpolation stencil. Must be odd. */
  char interp_type;   /* type of interpolation, 
			 ='e': explicit (uncoupled) interpolation 
			 ='i': implicit (coupled) interpolation   */

/* trimming style: minimize overlap = 0, maximize overlap = 1 */
  int trim_style;

/* the amount of warnings during the overlap algorithm is controlled by: */
  int verbose;

/* boundary surface missmatch tolerance */
  real max_b_dist;

/* bounding box */
  bounding_box *bb;

  char *name;

  int plot_mode;

  int show_boundary, show_volume_holes, correct_mismatch;

/* clipping */
  int n_clip_planes;
  GLdouble clip_plane[4][3];

/* relative size of spheres (used for plotting) */
  real sphere_size;

} overlapping_3d_grid;

/* prototypes for the functions in this library */
void 
overlap_3d_parameters(input_output *io_ptr, overlapping_3d_grid *over3d);
linked_list_member *
get_over3d_ptr( input_output *io_ptr, linked_list *head);
overlapping_3d_grid *
new_overlapping_3d_grid( char *name, linked_list *volume_mappings );
overlapping_3d_grid *
delete_over3d_grid( overlapping_3d_grid *over3d );
void 
set_3d_overlap(input_output *io_ptr, overlapping_3d_grid *over3d, 
	       linked_list *all_sgrids );
void
draw_overlapping_3d_grid( overlapping_3d_grid *over3d, int plot_mode );
int
overlapping_3d_plot_mode(input_output *io_ptr, overlapping_3d_grid *over3d, 
			 int plot_mode);
void
overlap_list_bb(bounding_box *bb, overlapping_3d_list *over3d_list);
int 
forward_vgrid_mapping( grid_point *gp_ptr, component_vgrid *vgrid );
int 
inverse_vgrid_mapping( grid_point *gp_ptr, component_vgrid *vgrid );

#endif
