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
#ifndef volume_h
#define volume_h

#include <window_defs.h>
#include <grid_point.h>
#include <surface.h>
#include <hole_surface.h>

/* typedefs */

/* define new types for lists of surface_mappings */
typedef linked_list volume_mapping_list;
typedef linked_list_member volume_mapping_link;

/* data structure for a 3-D component mapping */
typedef struct volume_mapping{

  int used_by_surface;      /* number of times this volume mapping is used by */
                            /* surface mappings that rely on it */

  char *name;               /* name of this component grid */

  int inverse_known;        /* inverse_known = 1 if the inverse of the mapping 
			       is defined in the grid function */
  int no_hole;

/* pointer to the function that defines the forward mapping */
  int (*forward_mapping)(grid_point *gp_ptr, void *volume_data_ptr); 

/* pointer to the function that defines the inverse mapping 
   set to NULL if the inverse mapping is not known */
  int (*inverse_mapping)(grid_point *gp_ptr, void *volume_data_ptr);

/* pointer to the function that interactively changes the default parameters */
  void (*set_volume)(input_output *io_ptr, struct volume_mapping *volume_ptr, 
		    surface_mapping_list *surface_list); 

/* rotation matrix and translation vectors */
  real_array_2d *rotation_;
  real axis[3], angle;
  real translation[3], pre_trans[3];

/* flag that indicates wether the mapping is translated and/or rotated */
  int transformed;

/* pointer to the function that deletes the specific data */
  void (*delete_volume_data)(void *volume_data_ptr);

/* pointers to the functions for saving and recovering a volume from a binary file */
  void (*write_volume_data)( int fd, void *volume_data_ptr );

/* pointer to the data necessary to compute the mapping */
  void *volume_data_ptr;

/* number of non-fictitious points in the ri directions, i=1,2,3 */
  int r1_points, r2_points, r3_points;    
  
  int r1_period, r2_period, r3_period; /* Periodicity flag meaning:
			       period=1: periodic
			       period=0: not periodic */


  int_array_2d *bc_ptr;      /* bc_vol(is,id): boundary condition on side is=1,2
			       in direction id=1,2,3 */
  int_array_2d *surf_ptr;  /* surf_vol(is,id): side is=1,2, direction 
			       id=1,2,3. surface(is,id)= a number identifying 
			       each side of the volume grid with a boundary surface. 
			       Used when sides of two or more grids share the
			       same physical boundary surface */
  int_array_1d *edge_curve_; /* edge_vol(i) corresponds to grid edges which join */
  /* grid faces with surface_label != 0. The index `i' is related to the */
  /* intersection line between the faces as follows: */
  /* i   edge_vol(i) */
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


/* bounding box */
  bounding_box *bb; 

  char *type;            
/* specific type of grid:
   1. Cartesian.
   2. Normal surface.
   ...
*/

  int plot_mode;
  int color;

/* to keep track of which grid surfaces to plot */
  real r1_surface, r2_surface, r3_surface;

/* to keep track of which surfaces to plot */
  int plot_surf_label, plot_bc;

} volume_mapping;

/* define some macros for computing indices of multi-dimensional arrays */
#define bc_vol(i,j)     compute_index_2d(volume->bc_ptr,i,j)
#define surf_vol(is,id) compute_index_2d(volume->surf_ptr,is,id)
#define edge_vol(i)     compute_index_1d(volume->edge_curve_,i)
#define rot_vol(i,j)    compute_index_2d(volume->rotation_,i,j)

/* prototypes for the functions in this library */
volume_mapping *
new_volume_mapping( char *name);
volume_mapping *
delete_volume_mapping( volume_mapping *volume );
int
volume_plot_mode(input_output *io_ptr, volume_mapping_list *volume_list,
		    volume_mapping *volume, int plot_mode, clip_planes *clip);
void
draw_volume( volume_mapping *volume, int plot_mode, clip_planes *clip );
volume_mapping *
choose_volume(input_output *io_ptr, char *name, 
	      surface_mapping_list *surface_mappings );
void 
set_volume(input_output *io_ptr, volume_mapping *volume, 
	   surface_mapping_list *surface_mappings );
void 
delete_specific_volume( volume_mapping *volume );
void 
write_specific_volume( int fd, volume_mapping *volume );
int 
inverse_volume_mapping( grid_point *gp_ptr, volume_mapping *volume );
int 
forward_volume_mapping( grid_point *gp_ptr, volume_mapping *volume );
void
volume_list_bb(bounding_box *bb, volume_mapping_list *volume_list);
volume_mapping_link *
get_volume_mapping( input_output *io_ptr, volume_mapping_list *head, int level );
void 
set_surf_label(input_output *io_ptr, volume_mapping *volume);
void 
set_volume_bc(input_output *io_ptr, volume_mapping *volume);
void 
set_edge_label(input_output *io_ptr, volume_mapping *volume);
void
check_jacobian(input_output *io_ptr, volume_mapping *volume);
void
transform_volume(input_output *io_ptr, volume_mapping *volume);
void
volume_bb( volume_mapping *volume );
void
file_volume( input_output *io_ptr, volume_mapping *volume );

/* end of volume.h */
#endif
