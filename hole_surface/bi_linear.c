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
#include "hole_surface_internal.h"

#define R0 ((i_loc-1)   * sgrid_ptr->r_step)
#define R1 ((i_loc+1-1) * sgrid_ptr->r_step)
#define ALPHA_0(r)   ((R1 - r)/sgrid_ptr->r_step)
#define D_ALPHA_0(r) (-1.0/sgrid_ptr->r_step)
#define ALPHA_1(r)   ((r - R0)/sgrid_ptr->r_step)
#define D_ALPHA_1(r) (1.0/sgrid_ptr->r_step)

#define S0 ((j_loc-1)  * sgrid_ptr->s_step)
#define S1 ((j_loc+1-1)* sgrid_ptr->s_step)
#define BETA_0(s)   ((S1 - s)/sgrid_ptr->s_step)
#define D_BETA_0(s) (-1.0/sgrid_ptr->s_step)
#define BETA_1(s)   ((s - S0)/sgrid_ptr->s_step)
#define D_BETA_1(s) (1.0/sgrid_ptr->s_step)

#define LIN_X(r,s) (ALPHA_0(r) * BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		    ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		    ALPHA_0(r) * BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		    ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) )
#define LIN_X_R(r,s) (D_ALPHA_0(r) * BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		      D_ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		      D_ALPHA_0(r) * BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		      D_ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) )
#define LIN_X_S(r,s) (ALPHA_0(r) * D_BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		      ALPHA_1(r) * D_BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		      ALPHA_0(r) * D_BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		      ALPHA_1(r) * D_BETA_1(s) * x(i_loc+1,j_loc+1) )

#define LIN_Y(r,s) (ALPHA_0(r) * BETA_0(s) * y(i_loc  ,j_loc  ) +  \
		    ALPHA_1(r) * BETA_0(s) * y(i_loc+1,j_loc  ) +  \
		    ALPHA_0(r) * BETA_1(s) * y(i_loc  ,j_loc+1) +  \
		    ALPHA_1(r) * BETA_1(s) * y(i_loc+1,j_loc+1) )
#define LIN_Y_R(r,s) (D_ALPHA_0(r) * BETA_0(s) * y(i_loc  ,j_loc  ) +  \
		      D_ALPHA_1(r) * BETA_0(s) * y(i_loc+1,j_loc  ) +  \
		      D_ALPHA_0(r) * BETA_1(s) * y(i_loc  ,j_loc+1) +  \
		      D_ALPHA_1(r) * BETA_1(s) * y(i_loc+1,j_loc+1) )
#define LIN_Y_S(r,s) (ALPHA_0(r) * D_BETA_0(s) * y(i_loc  ,j_loc  ) +  \
		      ALPHA_1(r) * D_BETA_0(s) * y(i_loc+1,j_loc  ) +  \
		      ALPHA_0(r) * D_BETA_1(s) * y(i_loc  ,j_loc+1) +  \
		      ALPHA_1(r) * D_BETA_1(s) * y(i_loc+1,j_loc+1) )

#define LIN_Z(r,s) (ALPHA_0(r) * BETA_0(s) * z(i_loc  ,j_loc  ) +  \
		    ALPHA_1(r) * BETA_0(s) * z(i_loc+1,j_loc  ) +  \
		    ALPHA_0(r) * BETA_1(s) * z(i_loc  ,j_loc+1) +  \
		    ALPHA_1(r) * BETA_1(s) * z(i_loc+1,j_loc+1) )
#define LIN_Z_R(r,s) (D_ALPHA_0(r) * BETA_0(s) * z(i_loc  ,j_loc  ) +  \
		      D_ALPHA_1(r) * BETA_0(s) * z(i_loc+1,j_loc  ) +  \
		      D_ALPHA_0(r) * BETA_1(s) * z(i_loc  ,j_loc+1) +  \
		      D_ALPHA_1(r) * BETA_1(s) * z(i_loc+1,j_loc+1) )
#define LIN_Z_S(r,s) (ALPHA_0(r) * D_BETA_0(s) * z(i_loc  ,j_loc  ) +  \
		      ALPHA_1(r) * D_BETA_0(s) * z(i_loc+1,j_loc  ) +  \
		      ALPHA_0(r) * D_BETA_1(s) * z(i_loc  ,j_loc+1) +  \
		      ALPHA_1(r) * D_BETA_1(s) * z(i_loc+1,j_loc+1) )

void 
bi_linear(grid_point *gp_ptr, discrete_surface *sgrid_ptr){
  int i_loc, j_loc;
  real r, s;

  r = gp_ptr->r;
  s = gp_ptr->s;

/* periodicity */
  if (sgrid_ptr->r_period && r < 0.0) r += 1.0;
  else if (sgrid_ptr->r_period && r > 1.0) r -= 1.0;

  if (sgrid_ptr->s_period && s < 0.0) s += 1.0;
  else if (sgrid_ptr->s_period && s > 1.0) s -= 1.0;

/* compute the index of the cell corresponding to the interpolation location */
  i_loc = 1 + r / sgrid_ptr->r_step;
  j_loc = 1 + s / sgrid_ptr->s_step;
/* avoid out of bounds */
  i_loc = int_max( 1, int_min( sgrid_ptr->r_dim-1, i_loc ) );
  j_loc = int_max( 1, int_min( sgrid_ptr->s_dim-1, j_loc ) );

  gp_ptr->x  = LIN_X(r, s);
  gp_ptr->y  = LIN_Y(r, s);
  gp_ptr->z  = LIN_Z(r, s);

  gp_ptr->xr = LIN_X_R(r, s);
  gp_ptr->yr = LIN_Y_R(r, s);
  gp_ptr->zr = LIN_Z_R(r, s);

  gp_ptr->xs = LIN_X_S(r, s);
  gp_ptr->ys = LIN_Y_S(r, s);
  gp_ptr->zs = LIN_Z_S(r, s);
}
