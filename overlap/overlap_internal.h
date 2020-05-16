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
#ifndef overlap_internal_h
#define overlap_internal_h

#include <overlap.h>

/* define some macros for computing indices of multi-dimensional arrays */
#define bc3(i,j)    compute_index_2d(vgrid->bc_ptr,i,j)
#define range3(i,j) compute_index_2d(vgrid->range_ptr,i,j)
#define surf3(is,id) compute_index_2d(vgrid->surf_ptr,is,id)
#define edge_curve(i) compute_index_1d(vgrid->edge_curve_,i)

#define flag3(i,j,k)  compute_index_3d(vgrid->flag_ptr,i,j,k)
#define other_flag3(i,j,k)  compute_index_3d(other_vgrid->flag_ptr,i,j,k)

#define x3(i,j,k)     compute_index_3d(vgrid->x_ptr,i,j,k)
#define y3(i,j,k)     compute_index_3d(vgrid->y_ptr,i,j,k)
#define z3(i,j,k)     compute_index_3d(vgrid->z_ptr,i,j,k)

/* prototypes for the functions in this library */
void 
show_overlap_parameters( overlapping_3d_grid *over3d );
discrete_surface *
new_sgrid_from_vgrid( component_vgrid *vgrid, int side, int dir, real max_b_dist );
interp_3d_point *
invert_grid(real xp, real yp, real zp, component_vgrid *vgrid, 
	    int query_surface, int query_curve, interp_3d_point *initial_guess);
void
cleanup_and_finish(component_vgrid *vgrid, overlapping_3d_grid *over3d );
void
trim_interpolation_points(component_vgrid *vgrid, overlapping_3d_grid *over3d);
int
consistency_check(component_vgrid *vgrid, overlapping_3d_grid *over3d);
void
explicit(component_vgrid *vgrid, overlapping_3d_grid *over3d);
void
initial_classification(linked_list_member *vgrid_link, overlapping_3d_grid *over3d);
int
grid_face_normal(real *n_vec, component_vgrid *vgrid, int i, int j, int k, 
		 int direction, int orientation);
int 
boundary_cell( int i0, int j0, int k0, component_vgrid *vgrid );
int 
good_interp_loc(interp_3d_point *new_interp3,
		component_vgrid *vgrid,
		overlapping_3d_grid *over3d);
void 
init_3d_arrays(volume_mapping *volume, component_vgrid *vgrid, overlapping_3d_grid *over3d);
void 
new_killed_point(int i, int j, int k, component_vgrid *vgrid);
bad_3d_point *
delete_bad_3d_list( bad_3d_point *last_bad_point );
void 
new_3d_bad_point(int i, int j, int k, int type, component_vgrid *vgrid);
void 
remove_deactive_3d_interp(component_vgrid *vgrid);
oct_tree_info *
delete_oct_tree( oct_tree_info *info );
interp_3d_point *
new_3d_interp(int i_loc, int j_loc, int k_loc, component_vgrid *vgrid_loc, 
	      real r_loc, real s_loc, real t_loc);
interp_3d_point *
delete_3d_interpolation_list( interp_3d_point *last_interp );
void
insert_3d_interp_point( interp_3d_point *new_interp, component_vgrid *vgrid );
component_vgrid *
delete_vgrid( component_vgrid *vgrid );
void
copy_from_volume_mapping( component_vgrid *vgrid, volume_mapping *volume );
void
prepare_3d_overlap( overlapping_3d_grid *over3d );
linked_list_member *
get_vgrid_ptr( input_output *io_ptr, linked_list *head);
component_vgrid *
new_vgrid( volume_mapping *volume );

#endif
