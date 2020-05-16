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
#include "overlap_internal.h"

#define R0 ((i_cell-range3(1,1))   * vgrid->r1_step)
#define R1 ((i_cell+1-range3(1,1)) * vgrid->r1_step)
#define ALPHA_0(r)   ((R1 - r)/vgrid->r1_step)
#define D_ALPHA_0(r) (-1.0/vgrid->r1_step)
#define ALPHA_1(r)   ((r - R0)/vgrid->r1_step)
#define D_ALPHA_1(r) (1.0/vgrid->r1_step)

#define S0 ((j_cell-range3(1,2))  * vgrid->r2_step)
#define S1 ((j_cell+1-range3(1,2))* vgrid->r2_step)
#define BETA_0(s)   ((S1 - s)/vgrid->r2_step)
#define D_BETA_0(s) (-1.0/vgrid->r2_step)
#define BETA_1(s)   ((s - S0)/vgrid->r2_step)
#define D_BETA_1(s) (1.0/vgrid->r2_step)

#define T0 ((k_cell-range3(1,3))  * vgrid->r3_step)
#define T1 ((k_cell+1-range3(1,3))* vgrid->r3_step)
#define GAMMA_0(t)   ((T1 - t)/vgrid->r3_step)
#define D_GAMMA_0(t) (-1.0/vgrid->r3_step)
#define GAMMA_1(t)   ((t - T0)/vgrid->r3_step)
#define D_GAMMA_1(t) (1.0/vgrid->r3_step)

#define LIN(v,r,s,t) (ALPHA_0(r) * BETA_0(s) * GAMMA_0(t) * v(i_cell  ,j_cell  ,k_cell  ) +  \
		      ALPHA_1(r) * BETA_0(s) * GAMMA_0(t) * v(i_cell+1,j_cell  ,k_cell  ) +  \
		      ALPHA_0(r) * BETA_1(s) * GAMMA_0(t) * v(i_cell  ,j_cell+1,k_cell  ) +  \
		      ALPHA_1(r) * BETA_1(s) * GAMMA_0(t) * v(i_cell+1,j_cell+1,k_cell  ) +  \
		      ALPHA_0(r) * BETA_0(s) * GAMMA_1(t) * v(i_cell  ,j_cell  ,k_cell+1) +  \
		      ALPHA_1(r) * BETA_0(s) * GAMMA_1(t) * v(i_cell+1,j_cell  ,k_cell+1) +  \
		      ALPHA_0(r) * BETA_1(s) * GAMMA_1(t) * v(i_cell  ,j_cell+1,k_cell+1) +  \
		      ALPHA_1(r) * BETA_1(s) * GAMMA_1(t) * v(i_cell+1,j_cell+1,k_cell+1))
#define LIN_R(v,r,s,t) (D_ALPHA_0(r) * BETA_0(s) * GAMMA_0(t) * v(i_cell  ,j_cell  ,k_cell  ) +  \
			D_ALPHA_1(r) * BETA_0(s) * GAMMA_0(t) * v(i_cell+1,j_cell  ,k_cell  ) +  \
			D_ALPHA_0(r) * BETA_1(s) * GAMMA_0(t) * v(i_cell  ,j_cell+1,k_cell  ) +  \
			D_ALPHA_1(r) * BETA_1(s) * GAMMA_0(t) * v(i_cell+1,j_cell+1,k_cell  ) +  \
			D_ALPHA_0(r) * BETA_0(s) * GAMMA_1(t) * v(i_cell  ,j_cell  ,k_cell+1) +  \
			D_ALPHA_1(r) * BETA_0(s) * GAMMA_1(t) * v(i_cell+1,j_cell  ,k_cell+1) +  \
			D_ALPHA_0(r) * BETA_1(s) * GAMMA_1(t) * v(i_cell  ,j_cell+1,k_cell+1) +  \
			D_ALPHA_1(r) * BETA_1(s) * GAMMA_1(t) * v(i_cell+1,j_cell+1,k_cell+1))
#define LIN_S(v,r,s,t) (ALPHA_0(r) * D_BETA_0(s) * GAMMA_0(t) * v(i_cell  ,j_cell  ,k_cell  ) +  \
			ALPHA_1(r) * D_BETA_0(s) * GAMMA_0(t) * v(i_cell+1,j_cell  ,k_cell  ) +  \
			ALPHA_0(r) * D_BETA_1(s) * GAMMA_0(t) * v(i_cell  ,j_cell+1,k_cell  ) +  \
			ALPHA_1(r) * D_BETA_1(s) * GAMMA_0(t) * v(i_cell+1,j_cell+1,k_cell  ) +  \
			ALPHA_0(r) * D_BETA_0(s) * GAMMA_1(t) * v(i_cell  ,j_cell  ,k_cell+1) +  \
			ALPHA_1(r) * D_BETA_0(s) * GAMMA_1(t) * v(i_cell+1,j_cell  ,k_cell+1) +  \
			ALPHA_0(r) * D_BETA_1(s) * GAMMA_1(t) * v(i_cell  ,j_cell+1,k_cell+1) +  \
			ALPHA_1(r) * D_BETA_1(s) * GAMMA_1(t) * v(i_cell+1,j_cell+1,k_cell+1))
#define LIN_T(v,r,s,t) (ALPHA_0(r) * BETA_0(s) * D_GAMMA_0(t) * v(i_cell  ,j_cell  ,k_cell  ) +  \
			ALPHA_1(r) * BETA_0(s) * D_GAMMA_0(t) * v(i_cell+1,j_cell  ,k_cell  ) +  \
			ALPHA_0(r) * BETA_1(s) * D_GAMMA_0(t) * v(i_cell  ,j_cell+1,k_cell  ) +  \
			ALPHA_1(r) * BETA_1(s) * D_GAMMA_0(t) * v(i_cell+1,j_cell+1,k_cell  ) +  \
			ALPHA_0(r) * BETA_0(s) * D_GAMMA_1(t) * v(i_cell  ,j_cell  ,k_cell+1) +  \
			ALPHA_1(r) * BETA_0(s) * D_GAMMA_1(t) * v(i_cell+1,j_cell  ,k_cell+1) +  \
			ALPHA_0(r) * BETA_1(s) * D_GAMMA_1(t) * v(i_cell  ,j_cell+1,k_cell+1) +  \
			ALPHA_1(r) * BETA_1(s) * D_GAMMA_1(t) * v(i_cell+1,j_cell+1,k_cell+1))

int 
tri_linear(grid_point *gp, component_vgrid *vgrid){
  int i_cell, j_cell, k_cell;
  real r, s, t;

  r = gp->r;
  s = gp->s;
  t = gp->t;

/* periodicity */
  if (vgrid->r1_period && r < 0.0) r += 1.0;
  else if (vgrid->r1_period && r > 1.0) r -= 1.0;

  if (vgrid->r2_period && s < 0.0) s += 1.0;
  else if (vgrid->r2_period && s > 1.0) s -= 1.0;

  if (vgrid->r3_period && t < 0.0) t += 1.0;
  else if (vgrid->r3_period && t > 1.0) t -= 1.0;

/* compute the index of the cell corresponding to the interpolation location */
  i_cell = range3(1,1) + r / vgrid->r1_step;
  j_cell = range3(1,2) + s / vgrid->r2_step;
  k_cell = range3(1,3) + t / vgrid->r3_step;

/* avoid out of bounds */
  i_cell = int_max( range3(1,1), int_min( range3(2,1)-1, i_cell ) );
  j_cell = int_max( range3(1,2), int_min( range3(2,2)-1, j_cell ) );
  k_cell = int_max( range3(1,3), int_min( range3(2,3)-1, k_cell ) );
  
/* evaluate the tri-linear mapping */
  gp->x  = LIN(x3, r, s, t);
  gp->y  = LIN(y3, r, s, t);
  gp->z  = LIN(z3, r, s, t);

  gp->xr = LIN_R(x3, r, s, t);
  gp->yr = LIN_R(y3, r, s, t);
  gp->zr = LIN_R(z3, r, s, t);

  gp->xs = LIN_S(x3, r, s, t);
  gp->ys = LIN_S(y3, r, s, t);
  gp->zs = LIN_S(z3, r, s, t);

  gp->xt = LIN_T(x3, r, s, t);
  gp->yt = LIN_T(y3, r, s, t);
  gp->zt = LIN_T(z3, r, s, t);

  return OK;
}
