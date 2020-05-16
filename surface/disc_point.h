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
#ifndef disc_point_h
#define disc_point_h

typedef struct {
  int nr, ns;
  real dr, ds;
  real_array_2d *x_, *y_, *z_;
} disc_point_info;

/* prototypes */
void *
init_disc_point(input_output *io_ptr, surface_mapping *surface );
void
alloc_disc_point_arrays( disc_point_info *info );
void
free_disc_point_arrays( disc_point_info *info );
int 
disc_point_mapping(grid_point *gp_ptr, surface_mapping *surface);
void 
set_disc_point(input_output *io_ptr, surface_mapping *surface);
void
delete_disc_point( void *tuut );
/* end prototypes */

/* macros for evaluating the grid point coordinates */
#define x_surf(i,j)     compute_index_2d(info->x_,i,j)
#define y_surf(i,j)     compute_index_2d(info->y_,i,j)
#define z_surf(i,j)     compute_index_2d(info->z_,i,j)

#endif


