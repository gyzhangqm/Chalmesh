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
#ifndef volume_point_h
#define volume_point_h

typedef struct {
  int n1, n2, n3;
  real d1, d2, d3;
  real_array_3d *x_, *y_, *z_;
} volume_point_info;

/* prototypes */
void *
init_volume_point(input_output *io_ptr, volume_mapping *volume );
/* end prototypes */

/* macros for evaluating the grid point coordinates */
#define x3_vol(i,j,k)     compute_index_3d(info->x_,i,j,k)
#define y3_vol(i,j,k)     compute_index_3d(info->y_,i,j,k)
#define z3_vol(i,j,k)     compute_index_3d(info->z_,i,j,k)

#endif


