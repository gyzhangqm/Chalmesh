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

int
newton_search(real xp, real yp, real zp, real *r_loc, real *s_loc, real *n_loc, 
	      discrete_surface *project_onto){
  int iter, info;
  const int max_iter=15;
  grid_point gp;
/* normal and temporary stuff*/
  real normal[3], normal_r[3], normal_s[3], res;
/* linear system stuff */
  static real_array_2d *b_ = NULL, *A_ = NULL;
  long int pivot[3];

#define A(i,j) compute_index_2d(A_,i,j)
#define b(i) compute_index_2d(b_,i,1)

/* allocate memory on the first call (and never deallocate = small memory leak) */
  if (A_==NULL){
    A_ = create_real_array_2d( 3, 3 );
    b_ = create_real_array_2d( 3, 1 );
  }

  iter = 0;
  do{
/* evaluate the surface */
    gp.r = *r_loc; gp.s = *s_loc;
    bi_linear(&gp, project_onto);
/* compute the normal */
    compute_continuous_normal(normal, normal_r, normal_s, gp.r, gp.s, project_onto);

/* compute the residual */
    b(1) = xp - *n_loc*normal[0] - gp.x;
    b(2) = yp - *n_loc*normal[1] - gp.y;
    b(3) = zp - *n_loc*normal[2] - gp.z;

    if (( res = sqrt( b(1)*b(1) + b(2)*b(2) + b(3)*b(3)) ) <= NEWTON_EPS) break;

/* assign the Jacobian matrix */
    A(1,1) = gp.xr + *n_loc*normal_r[0]; 
    A(2,1) = gp.yr + *n_loc*normal_r[1]; 
    A(3,1) = gp.zr + *n_loc*normal_r[2]; 
    
    A(1,2) = gp.xs + *n_loc*normal_s[0]; 
    A(2,2) = gp.ys + *n_loc*normal_s[1]; 
    A(3,2) = gp.zs + *n_loc*normal_s[2]; 

    A(1,3) = normal[0];
    A(2,3) = normal[1];
    A(3,3) = normal[2];

/* Solve the system */
    if ((info=LUsolve(A_, pivot, b_)))
      {
	if (info > 1)
	  printf("Newton: Singular jacobian!\n");
	else
	  printf("Newton: Invalid matrix!\n");
	return ERROR;
      }

/* update */
    *r_loc += b(1);
    *s_loc += b(2);
    *n_loc += b(3);

  } while (++iter < max_iter);

  if (res > NEWTON_EPS){
/*#ifdef D_PROJ*/
    printf("Warning: Poor convergence in Newton iteration while projecting onto "
	   "grid `%s' at x=%f, y=%f, z=%f.\n"
	   "Last residual = %e\n", project_onto->name, xp, yp, zp, res);
/*#endif*/
    return ERROR;
  }

  return OK;
}

