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
#ifndef hole_surface_internal_h
#define hole_surface_internal_h

#include "hole_surface.h"

typedef struct intersection_point
{
  real ic[3], sin_theta;
  discrete_surface *sgrid;
  int unique;
  struct intersection_point *prev_ip;
} intersection_point;

typedef struct stack_of_ip
{
  int n_ip;
  intersection_point *last_ip;
} stack_of_ip;

/* tri_array_node * */
/* delete_tri_array_node( tri_array_node *node ); */
stack_of_ip *
new_stack_of_ip(void);
stack_of_ip *
delete_stack_of_ip(stack_of_ip *stack);
intersection_point *
new_intersection_point(real xp, real yp, real zp, real sin_theta,
		       discrete_surface *sgrid_ptr);
void
push_ip(intersection_point *ip_ptr, stack_of_ip *stack);
int
compute_grid_normal(real n_vec[3], int i, int j, discrete_surface *sgrid_ptr);
int
compute_cell_normal(real n_vec[3], int ic, int jc, discrete_surface *sgrid_ptr);
void
compute_continuous_normal(real normal[3], real normal_r[3], real normal_s[3],
			  real r, real s, discrete_surface *sgrid_ptr);

#endif
