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
#ifndef hyp_grid_h
#define hyp_grid_h

#include "stretchings.h"

typedef struct{
  surface_mapping *surface;
  real width, curv_factor, v_min, smooth_factor, time_factor;
  int n1, n2, n3;
  real d1, d2, d3;
  real_array_3d *x_, *y_, *z_;

/* boundary conditions */
  int_array_2d *hyp_bc_;

/* flag to indicate if the mapping is up to date (=1), needs to be (re-)computed (=0), */
/* or if the computation was unsuccessful (=-1). */
  int status;

/* pointer to the normal stretching structures */
  generic_stretching *normal_stretch;
} hyp_grid_data;

/* prototypes */
hyp_grid_data *
init_hyp_grid(input_output *io_ptr, volume_mapping *volume, surface_mapping *surface);

/* macros for evaluating the grid point coordinates */
#define x_hyp(i,j,k)     compute_index_3d(info->x_,i,j,k)
#define y_hyp(i,j,k)     compute_index_3d(info->y_,i,j,k)
#define z_hyp(i,j,k)     compute_index_3d(info->z_,i,j,k)

/* macro for the boundary condition */
#define hyp_bc(i,j) compute_index_2d(info->hyp_bc_,i,j)
#endif
