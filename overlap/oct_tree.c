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

/* static prototypes */
static oct_tree_info *
new_oct_tree(int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, 
	     component_vgrid *vgrid);
static interp_3d_point *
inside_cube(real xp, real yp, real zp, oct_tree_info *leaf,
	    component_vgrid *vgrid, real max_b_dist, 
	    int query_surface, int query_curve );
static interp_3d_point *
search_oct_tree(real xp, real yp, real zp, oct_tree_info *info, component_vgrid *vgrid, 
		real max_b_dist, int query_surface, int query_curve );

/* end static prototypes */

static oct_tree_info *
new_oct_tree(int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, 
	     component_vgrid *vgrid){
  int i, j, k;
  oct_tree_info *info;

/* don't save any empty leafs */
  if (i_min == i_max || j_min == j_max || k_min == k_max) 
    return NULL;

  info = (oct_tree_info *) malloc( sizeof(oct_tree_info) );

  info->i_min = i_min;
  info->i_max = i_max;
  info->j_min = j_min;
  info->j_max = j_max;
  info->k_min = k_min;
  info->k_max = k_max;

/* initialize bounding box */
  info->x_min =  1.e10;
  info->x_max = -1.e10;
  info->y_min =  1.e10;
  info->y_max = -1.e10;
  info->z_min =  1.e10;
  info->z_max = -1.e10;

  k=k_min;
  for (i=i_min; i<=i_max; i++){
    for (j=j_min; j<=j_max; j++){
      info->x_min = real_min( info->x_min, x3(i,j,k) );
      info->x_max = real_max( info->x_max, x3(i,j,k) );
      info->y_min = real_min( info->y_min, y3(i,j,k) );
      info->y_max = real_max( info->y_max, y3(i,j,k) );
      info->z_min = real_min( info->z_min, z3(i,j,k) );
      info->z_max = real_max( info->z_max, z3(i,j,k) );
    }
  }

  k=k_max;
  for (i=i_min; i<=i_max; i++){
    for (j=j_min; j<=j_max; j++){
      info->x_min = real_min( info->x_min, x3(i,j,k) );
      info->x_max = real_max( info->x_max, x3(i,j,k) );
      info->y_min = real_min( info->y_min, y3(i,j,k) );
      info->y_max = real_max( info->y_max, y3(i,j,k) );
      info->z_min = real_min( info->z_min, z3(i,j,k) );
      info->z_max = real_max( info->z_max, z3(i,j,k) );
    }
  }

  i=i_min;
  for (k=k_min; k<=k_max; k++){
    for (j=j_min; j<=j_max; j++){
      info->x_min = real_min( info->x_min, x3(i,j,k) );
      info->x_max = real_max( info->x_max, x3(i,j,k) );
      info->y_min = real_min( info->y_min, y3(i,j,k) );
      info->y_max = real_max( info->y_max, y3(i,j,k) );
      info->z_min = real_min( info->z_min, z3(i,j,k) );
      info->z_max = real_max( info->z_max, z3(i,j,k) );
    }
  }

  i=i_max;
  for (k=k_min; k<=k_max; k++){
    for (j=j_min; j<=j_max; j++){
      info->x_min = real_min( info->x_min, x3(i,j,k) );
      info->x_max = real_max( info->x_max, x3(i,j,k) );
      info->y_min = real_min( info->y_min, y3(i,j,k) );
      info->y_max = real_max( info->y_max, y3(i,j,k) );
      info->z_min = real_min( info->z_min, z3(i,j,k) );
      info->z_max = real_max( info->z_max, z3(i,j,k) );
    }
  }

  j=j_min;
  for (i=i_min; i<=i_max; i++){
    for (k=k_min; k<=k_max; k++){
      info->x_min = real_min( info->x_min, x3(i,j,k) );
      info->x_max = real_max( info->x_max, x3(i,j,k) );
      info->y_min = real_min( info->y_min, y3(i,j,k) );
      info->y_max = real_max( info->y_max, y3(i,j,k) );
      info->z_min = real_min( info->z_min, z3(i,j,k) );
      info->z_max = real_max( info->z_max, z3(i,j,k) );
    }
  }

  j=j_max;
  for (i=i_min; i<=i_max; i++){
    for (k=k_min; k<=k_max; k++){
      info->x_min = real_min( info->x_min, x3(i,j,k) );
      info->x_max = real_max( info->x_max, x3(i,j,k) );
      info->y_min = real_min( info->y_min, y3(i,j,k) );
      info->y_max = real_max( info->y_max, y3(i,j,k) );
      info->z_min = real_min( info->z_min, z3(i,j,k) );
      info->z_max = real_max( info->z_max, z3(i,j,k) );
    }
  }

/* are we on a leaf ? */
  if (i_max-i_min <= 1 && j_max-j_min <= 1 && 
      k_max-k_min <= 1){
    info->children = NULL;
  }
  else{
/* recurse */
    info->children = (oct_tree_info **) malloc( 8*sizeof(oct_tree_info *) );
    info->children[0] = new_oct_tree(i_min, (i_max+i_min)/2, j_min, (j_max+j_min)/2,
				     k_min, (k_max+k_min)/2, vgrid );
    info->children[1] = new_oct_tree((i_max+i_min)/2, i_max, j_min, (j_max+j_min)/2,
				     k_min, (k_max+k_min)/2, vgrid );
    info->children[2] = new_oct_tree(i_min, (i_max+i_min)/2, (j_max+j_min)/2, j_max,
				     k_min, (k_max+k_min)/2, vgrid );
    info->children[3] = new_oct_tree((i_max+i_min)/2, i_max, (j_max+j_min)/2, j_max,
				     k_min, (k_max+k_min)/2, vgrid );
    info->children[4] = new_oct_tree(i_min, (i_max+i_min)/2, j_min, (j_max+j_min)/2,
				     (k_max+k_min)/2, k_max, vgrid );
    info->children[5] = new_oct_tree((i_max+i_min)/2, i_max, j_min, (j_max+j_min)/2,
				     (k_max+k_min)/2, k_max, vgrid );
    info->children[6] = new_oct_tree(i_min, (i_max+i_min)/2, (j_max+j_min)/2, j_max,
				     (k_max+k_min)/2, k_max, vgrid );
    info->children[7] = new_oct_tree((i_max+i_min)/2, i_max, (j_max+j_min)/2, j_max,
				     (k_max+k_min)/2, k_max, vgrid );
  }
/* done */
  return info;
}

oct_tree_info *
delete_oct_tree( oct_tree_info *info ){
  int i;

  if (info == NULL)
    return NULL;

/* delete subtrees */
  if (info->children){
    for (i=0; i<8; i++)
      delete_oct_tree( info->children[i] );
/* delete the pointers to the subtrees */
    free(info->children);
  }

/*   delete_oct_tree( info->low_left_near ); */
/*   delete_oct_tree( info->low_right_near); */
/*   delete_oct_tree( info->upp_left_near ); */
/*   delete_oct_tree( info->upp_right_near); */

/*   delete_oct_tree( info->low_left_far ); */
/*   delete_oct_tree( info->low_right_far); */
/*   delete_oct_tree( info->upp_left_far ); */
/*   delete_oct_tree( info->upp_right_far); */

/* delete the structure itself */
  free(info);

  return NULL;
}

static interp_3d_point *
newton(real xp, real yp, real zp, component_vgrid *vgrid, real max_b_dist,
       interp_3d_point *initial_guess, int query_surface, int query_curve)
{
  interp_3d_point *solution=NULL;
  int iter, i_cell, j_cell, k_cell, info;
  const int max_iter=10;
  grid_point gp;
  real r_loc, s_loc, t_loc, res, prev_diff, diff
    , r_dist_low = NEWTON_EPS*vgrid->length_scale
      , r_dist_high = NEWTON_EPS*vgrid->length_scale
	, s_dist_low = NEWTON_EPS*vgrid->length_scale
	  , s_dist_high = NEWTON_EPS*vgrid->length_scale
	    , t_dist_low = NEWTON_EPS*vgrid->length_scale
	      , t_dist_high = NEWTON_EPS*vgrid->length_scale;
/* linear system stuff */
  static real_array_2d *b_ = NULL, *A_ = NULL;
  long int pivot[3];

#define A(i,j) compute_index_2d(A_,i,j)
#define b(i) compute_index_2d(b_,i,1)

/* allocate memory on the first call (and never deallocate = small memory leak) */
  if (A_==NULL)
    {
      A_ = create_real_array_2d( 3, 3 );
      b_ = create_real_array_2d( 3, 1 );
    }

/* initial guess */
  r_loc = initial_guess->r_loc;
  s_loc = initial_guess->s_loc;
  t_loc = initial_guess->t_loc;

/* get the iteration going */
  prev_diff = 1.e10;

  iter = 0;
  do
    {
/* evaluate the volume mapping */
    gp.r = r_loc; gp.s = s_loc; gp.t = t_loc;
    if (!forward_vgrid_mapping(&gp, vgrid))
      return NULL;

    b(1) = xp - gp.x;
    b(2) = yp - gp.y;
    b(3) = zp - gp.z;

    if (( res = sqrt( b(1)*b(1) + b(2)*b(2) + b(3)*b(3)) ) <= 
	NEWTON_EPS*vgrid->length_scale) break;

    A(1,1) = gp.xr; A(1,2) = gp.xs; A(1,3) = gp.xt;
    A(2,1) = gp.yr; A(2,2) = gp.ys; A(2,3) = gp.yt;
    A(3,1) = gp.zr; A(3,2) = gp.zs; A(3,3) = gp.zt;
    
/* Solve the system */
    if ((info=LUsolve(A_, pivot, b_)))
      {
	if (info > 1)
	  printf("Newton: Singular jacobian!\n");
	else
	  printf("Newton: Invalid matrix!\n");
	return NULL;
      }

/* update */
    r_loc = r_loc + b(1);
    s_loc = s_loc + b(2);
    t_loc = t_loc + b(3);


/* make sure we don't get overflow while computing the difference and
   that the increment is decreasing at least a little */
    if (fabs(b(1)) > 1.e10 || fabs(b(2)) > 1.e10 || fabs(b(3)) > 1.e10 || 
	(diff = b(1)*b(1) + b(2)*b(2) + b(3)*b(3)) > 0.95*prev_diff)
      {
/* tmp */
/*	printf("Newton diverged starting from the initial guess!\n");*/
	return NULL;
      }
    prev_diff = diff;

  } 
  while (++iter < max_iter);

/* tmp */
/*  printf("Newton converged starting from the initial guess!\n");*/


/* periodicity */
  if (vgrid->r1_period)
    {
      if (r_loc < 0.0) r_loc += 1.0;
      else if (r_loc > 1.0) r_loc -= 1.0;
    }

  if (vgrid->r2_period)
    {
      if (s_loc < 0.0) s_loc += 1.0;
      else if (s_loc > 1.0) s_loc -= 1.0;
    }

  if (vgrid->r3_period)
    {
      if (t_loc < 0.0) t_loc += 1.0;
      else if (t_loc > 1.0) t_loc -= 1.0;
    }

/* get the cell number */
  i_cell = range3(1,1) + r_loc / vgrid->r1_step;
  j_cell = range3(1,2) + s_loc / vgrid->r2_step;
  k_cell = range3(1,3) + t_loc / vgrid->r3_step;

/* avoid out of bounds */
  i_cell = int_max( range3(1,1), int_min( range3(2,1)-1, i_cell ) );
  j_cell = int_max( range3(1,2), int_min( range3(2,2)-1, j_cell ) );
  k_cell = int_max( range3(1,3), int_min( range3(2,3)-1, k_cell ) );

/* compute the tolerance in the (r,s,t) directions if the cell is an edge cell */
  if (query_curve > 0){
/* a solution at a singular point is considered to be invalid! */
    if (fabs(gp.xr) + fabs(gp.yr) + fabs(gp.zr) < 1.e-10 ||
	fabs(gp.xs) + fabs(gp.ys) + fabs(gp.zs) < 1.e-10 ||
	fabs(gp.xt) + fabs(gp.yt) + fabs(gp.zt) < 1.e-10)
      return NULL;

    if (query_curve == edge_curve(1) && 
	i_cell == range3(1,1) && j_cell == range3(1,2)){
      r_dist_low  += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      s_dist_low  += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_curve == edge_curve(2) && 
	     i_cell == range3(2,1)-1 && j_cell == range3(1,2)){
      r_dist_high += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      s_dist_low  += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_curve == edge_curve(3) && 
	     i_cell == range3(1,1) && j_cell == range3(2,2)-1){
      r_dist_low  += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      s_dist_high += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_curve == edge_curve(4) && 
	     i_cell == range3(2,1)-1 && j_cell == range3(2,2)-1){
      r_dist_high += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      s_dist_high += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_curve == edge_curve(5) && 
	     i_cell == range3(1,1) && k_cell == range3(1,3)){
      r_dist_low  += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      t_dist_low  += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(6) && 
	     i_cell == range3(2,1)-1 && k_cell == range3(1,3)){
      r_dist_high += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      t_dist_low  += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(7) && 
	     i_cell == range3(1,1) && k_cell == range3(2,3)-1){
      r_dist_low  += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      t_dist_high += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(8) && 
	     i_cell == range3(2,1)-1 && k_cell == range3(2,3)-1){
      r_dist_high += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      t_dist_high += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(9) && 
	     j_cell == range3(1,2) && k_cell == range3(1,3)){
      s_dist_low  += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      t_dist_low  += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(10) && 
	     j_cell == range3(2,2)-1 && k_cell == range3(1,3)){
      s_dist_high += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      t_dist_low  += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(11) && 
	     j_cell == range3(1,2) && k_cell == range3(2,3)-1){
      s_dist_low  += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      t_dist_high += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(12) && 
	     j_cell == range3(2,2)-1 && k_cell == range3(2,3)-1){
      s_dist_high += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      t_dist_high += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
  }   

/* compute the tolerance in the (r,s,t) directions if the cell is on a face with */
/* non-zero surface-label */
/* r1(r)-direction */
  if (query_surface > 0){
/* a solution at a singular point is considered to be invalid! */
    if (fabs(gp.xr) + fabs(gp.yr) + fabs(gp.zr) < 1.e-10 ||
	fabs(gp.xs) + fabs(gp.ys) + fabs(gp.zs) < 1.e-10 ||
	fabs(gp.xt) + fabs(gp.yt) + fabs(gp.zt) < 1.e-10)
      return NULL;

    if (query_surface == surf3(1,1) && i_cell == range3(1,1) ){
      r_dist_low += max_b_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
    }
    else if (query_surface == surf3(2,1) && i_cell == range3(2,1)-1 ){
      r_dist_high += max_b_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
    }
/* r2(s)-direction */
    if (query_surface == surf3(1,2) && j_cell == range3(1,2)){
      s_dist_low += max_b_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_surface == surf3(2,2) && j_cell == range3(2,2)-1){
      s_dist_high += max_b_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
/* r3(t)-direction */
    if (query_surface == surf3(1,3) && k_cell == range3(1,3)){
      t_dist_low += max_b_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_surface == surf3(2,3) && k_cell == range3(2,3)-1){
      t_dist_high += max_b_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
  }
  
/* is (r_loc, s_loc, t_loc) inside the slightly enlarged unit cube? */
  if (-r_dist_low < r_loc && r_loc < 1.0 + r_dist_high &&
      -s_dist_low < s_loc && s_loc < 1.0 + s_dist_high &&
      -t_dist_low < t_loc && t_loc < 1.0 + t_dist_high){
/* make an interp_3d_point to report the result */
    solution = new_3d_interp( i_cell, j_cell, k_cell, vgrid, r_loc, s_loc, t_loc );
    if (res > NEWTON_EPS*vgrid->length_scale)
      printf("Warning: Poor convergence in newton in grid `%s' for "
	     "the Cartesian \ncoordinate (%e, %e, %e).\n"
	     "After %i iterations, the residual was %e\n"
	     "and the solution (%e, %e, %e)\n", vgrid->name,
	     xp, yp, zp, iter, res, r_loc, s_loc, t_loc);
  }
  return solution;
}

static interp_3d_point *
inside_cube(real xp, real yp, real zp, oct_tree_info *leaf,
	    component_vgrid *vgrid, real max_b_dist, 
	    int query_surface, int query_curve ){
  interp_3d_point *solution=NULL;
  int iter, i_cell, j_cell, k_cell;
  const int max_iter=10;
  grid_point gp;
  real r_loc, s_loc, t_loc, res, prev_res, tolerance
    , r_0, r_1, s_0, s_1, t_0, t_1, r_prev, s_prev, t_prev, prev_diff, diff
      , r_dist_low = NEWTON_EPS*vgrid->length_scale
	, r_dist_high = NEWTON_EPS*vgrid->length_scale
	  , s_dist_low = NEWTON_EPS*vgrid->length_scale
	    , s_dist_high = NEWTON_EPS*vgrid->length_scale
	      , t_dist_low = NEWTON_EPS*vgrid->length_scale
		, t_dist_high = NEWTON_EPS*vgrid->length_scale;
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

/* for now, assume the leaf-size to be one */
  i_cell = leaf->i_min; j_cell = leaf->j_min; k_cell = leaf->k_min;

/* initial guess */
  r_prev = r_loc = vgrid->r1_step * (i_cell - range3(1,1) + 0.5); 
  s_prev = s_loc = vgrid->r2_step * (j_cell - range3(1,2) + 0.5);
  t_prev = t_loc = vgrid->r3_step * (k_cell - range3(1,3) + 0.5);

/* get the iteration going */
  prev_diff = 1.e10;

  iter = 0;
  do{
/* evaluate the volume mapping */
    gp.r = r_loc; gp.s = s_loc; gp.t = t_loc;
    if (!forward_vgrid_mapping(&gp, vgrid))
      return NULL;

    b(1) = xp - gp.x;
    b(2) = yp - gp.y;
    b(3) = zp - gp.z;

    if (( res = sqrt( b(1)*b(1) + b(2)*b(2) + b(3)*b(3)) ) <= 
	NEWTON_EPS*vgrid->length_scale) break;

    A(1,1) = gp.xr; A(1,2) = gp.xs; A(1,3) = gp.xt;
    A(2,1) = gp.yr; A(2,2) = gp.ys; A(2,3) = gp.yt;
    A(3,1) = gp.zr; A(3,2) = gp.zs; A(3,3) = gp.zt;
    
/* Solve the system */
    LUsolve(A_, pivot, b_); 

/* save previous guess and residual in case the residual increases */
    r_prev = r_loc; s_prev = s_loc; t_prev = t_loc;
    prev_res = res;

/* update */
    r_loc = r_loc + b(1);
    s_loc = s_loc + b(2);
    t_loc = t_loc + b(3);

/* make sure we don't get overflow while computing the difference and
   that the increment is decreasing at least a little */
    if (fabs(b(1)) > 1.e10 || fabs(b(2)) > 1.e10 || fabs(b(3)) > 1.e10 || 
	(diff = b(1)*b(1) + b(2)*b(2) + b(3)*b(3)) > 0.95*prev_diff)
      {
/* Hypothesis: Newton diverges because (xp, yp, zp) is outside the cell */
      return NULL;
    }
    prev_diff = diff;

  } while (++iter < max_iter);

/* periodicity */
  if (vgrid->r1_period)
    {
      if (r_loc < 0.0) r_loc += 1.0;
      else if (r_loc > 1.0) r_loc -= 1.0;
    }

  if (vgrid->r2_period)
    {
      if (s_loc < 0.0) s_loc += 1.0;
      else if (s_loc > 1.0) s_loc -= 1.0;
    }

  if (vgrid->r3_period)
    {
      if (t_loc < 0.0) t_loc += 1.0;
      else if (t_loc > 1.0) t_loc -= 1.0;
    }

/* get the cell number */
  i_cell = range3(1,1) + r_loc / vgrid->r1_step;
  j_cell = range3(1,2) + s_loc / vgrid->r2_step;
  k_cell = range3(1,3) + t_loc / vgrid->r3_step;

/* avoid out of bounds */
  i_cell = int_max( range3(1,1), int_min( range3(2,1)-1, i_cell ) );
  j_cell = int_max( range3(1,2), int_min( range3(2,2)-1, j_cell ) );
  k_cell = int_max( range3(1,3), int_min( range3(2,3)-1, k_cell ) );
  
/* compute the tolerance in the (r,s,t) directions if the cell is an edge cell */
  if (query_curve > 0){
    if (query_curve == edge_curve(1) && 
	i_cell == range3(1,1) && j_cell == range3(1,2)){
      r_dist_low  += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      s_dist_low  += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_curve == edge_curve(2) && 
	     i_cell == range3(2,1)-1 && j_cell == range3(1,2)){
      r_dist_high += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      s_dist_low  += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_curve == edge_curve(3) && 
	     i_cell == range3(1,1) && j_cell == range3(2,2)-1){
      r_dist_low  += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      s_dist_high += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_curve == edge_curve(4) && 
	     i_cell == range3(2,1)-1 && j_cell == range3(2,2)-1){
      r_dist_high += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      s_dist_high += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_curve == edge_curve(5) && 
	     i_cell == range3(1,1) && k_cell == range3(1,3)){
      r_dist_low  += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      t_dist_low  += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(6) && 
	     i_cell == range3(2,1)-1 && k_cell == range3(1,3)){
      r_dist_high += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      t_dist_low  += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(7) && 
	     i_cell == range3(1,1) && k_cell == range3(2,3)-1){
      r_dist_low  += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      t_dist_high += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(8) && 
	     i_cell == range3(2,1)-1 && k_cell == range3(2,3)-1){
      r_dist_high += vgrid->max_e_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      t_dist_high += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(9) && 
	     j_cell == range3(1,2) && k_cell == range3(1,3)){
      s_dist_low  += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      t_dist_low  += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(10) && 
	     j_cell == range3(2,2)-1 && k_cell == range3(1,3)){
      s_dist_high += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      t_dist_low  += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(11) && 
	     j_cell == range3(1,2) && k_cell == range3(2,3)-1){
      s_dist_low  += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      t_dist_high += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_curve == edge_curve(12) && 
	     j_cell == range3(2,2)-1 && k_cell == range3(2,3)-1){
      s_dist_high += vgrid->max_e_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      t_dist_high += vgrid->max_e_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
  }   

/* compute the tolerance in the (r,s,t) directions if the cell is on a face with */
/* non-zero surface-label */
/* r1(r)-direction */
  if (query_surface > 0){
    if (query_surface == surf3(1,1) && i_cell == range3(1,1) ){
      r_dist_low += max_b_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
    }
    else if (query_surface == surf3(2,1) && i_cell == range3(2,1)-1 ){
      r_dist_high += max_b_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
    }
/* r2(s)-direction */
    if (query_surface == surf3(1,2) && j_cell == range3(1,2)){
      s_dist_low += max_b_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if (query_surface == surf3(2,2) && j_cell == range3(2,2)-1){
      s_dist_high += max_b_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
/* r3(t)-direction */
    if (query_surface == surf3(1,3) && k_cell == range3(1,3)){
      t_dist_low += max_b_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
    else if (query_surface == surf3(2,3) && k_cell == range3(2,3)-1){
      t_dist_high += max_b_dist / sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    }
  }
  
/* is (r_loc, s_loc, t_loc) inside the slightly enlarged unit cube? */
  if (-r_dist_low < r_loc && r_loc < 1.0 + r_dist_high &&
      -s_dist_low < s_loc && s_loc < 1.0 + s_dist_high &&
      -t_dist_low < t_loc && t_loc < 1.0 + t_dist_high){
/* make an interp_3d_point to report the result */
    solution = new_3d_interp( i_cell, j_cell, k_cell, vgrid, r_loc, s_loc, t_loc );
    if (res > NEWTON_EPS*vgrid->length_scale)
      printf("Warning: Poor convergence in Newton iteration in inside_cube for "
	     "cell (%i, %i, %i) in grid `%s'.\n"
	     "After %i iterations, the residual was %e\n", i_cell, j_cell, k_cell, 
	     vgrid->name, iter, res);
  }


  return solution;

/* compute the cell boundaries in the parameter plane */
/*   r_0 = vgrid->r1_step * (i_cell - range3(1,1)) - r_dist_low;  */
/*   r_1 = vgrid->r1_step * (i_cell + 1 - range3(1,1)) + r_dist_high;  */
/*   s_0 = vgrid->r2_step * (j_cell - range3(1,2)) - s_dist_low; */
/*   s_1 = vgrid->r2_step * (j_cell + 1 - range3(1,2)) + s_dist_high; */
/*   t_0 = vgrid->r3_step * (k_cell - range3(1,3)) - t_dist_low; */
/*   t_1 = vgrid->r3_step * (k_cell + 1 - range3(1,3)) + t_dist_high; */

/* check if the point is inside the cell */
/*   if (r_0 <= r_loc && r_loc <= r_1 &&  */
/*       s_0 <= s_loc && s_loc <= s_1 && */
/*       t_0 <= t_loc && t_loc <= t_1){ */
/*     solution = new_3d_interp( i_cell, j_cell, k_cell, vgrid, r_loc, s_loc, t_loc ); */
/*   } */
/*   return solution; */
}


static interp_3d_point *
search_oct_tree(real xp, real yp, real zp, oct_tree_info *info, component_vgrid *vgrid, 
		real max_b_dist, int query_surface, int query_curve ){
  interp_3d_point *interp3=NULL;
/*  real max_b_dist;*/
  int i;

  if (info == NULL) return NULL;

/* check if we are inside the (slightly enlarged) bounding box */
  if (info->x_min-max_b_dist <= xp && xp < info->x_max+max_b_dist &&
      info->y_min-max_b_dist <= yp && yp < info->y_max+max_b_dist &&
      info->z_min-max_b_dist <= zp && zp < info->z_max+max_b_dist){
/* check if we are on a leaf */
    if (info->children == NULL){
/* perform a closer investigation to see if the point really is in the cell */
      interp3 = inside_cube(xp, yp, zp, info, vgrid, max_b_dist, 
			    query_surface, query_curve);
    }
/* Recurse. Observe that the bounding boxes can be overlapping. The different */
/* branches are therefore not mutually exclusive */
    else{
      for (i=0; i<8; i++){
	if ((interp3 = search_oct_tree(xp, yp, zp, info->children[i], vgrid, 
				       max_b_dist, query_surface, query_curve )) != NULL)
	  return interp3;
      }
    }
  }

/* return either NULL or the newly created interpolation point */
  return interp3;
}


interp_3d_point *
invert_grid(real xp, real yp, real zp, component_vgrid *vgrid, 
	    int query_surface, int query_curve, interp_3d_point *initial_guess){
  grid_point gp;
  interp_3d_point *solution=NULL, surface_solution;
  real r_dist_low = NEWTON_EPS*vgrid->length_scale
    , r_dist_high = NEWTON_EPS*vgrid->length_scale
      , s_dist_low = NEWTON_EPS*vgrid->length_scale
	, s_dist_high = NEWTON_EPS*vgrid->length_scale
	  , t_dist_low = NEWTON_EPS*vgrid->length_scale
	    , t_dist_high = NEWTON_EPS*vgrid->length_scale;
  real point[3];
  int i_cell, j_cell, k_cell;
  inverse_point *result;
  discrete_surface_link *sgrid_link;
  discrete_surface *sgrid_ptr;

  if (!vgrid){
    printf("FATAL ERROR: invert_grid was called with vgrid == NULL.\n");
    exit(-1);
  }

  if (!vgrid->analytical_inverse){
/* make the oct tree on the first call */
    if (!vgrid->oct_tree) 
      vgrid->oct_tree = new_oct_tree(range3(1,1), range3(2,1), 
				     range3(1,2), range3(2,2), 
				     range3(1,3), range3(2,3), 
				     vgrid);
    point[0] = xp; point[1] = yp; point[2] = zp;
/* first try the the initial guess (if it applies to the current vgrid) */
    if (initial_guess && initial_guess->vgrid_loc == vgrid &&
	(solution = newton(xp, yp, zp, vgrid, vgrid->max_s_dist, initial_guess, 
			   query_surface, query_curve)))
      return solution;

/* try the get a better initial guess from the bounding surface */
    if (query_surface != 0)
      {
	surface_solution.vgrid_loc = NULL;
/* check each patch */
	for (sgrid_link = vgrid->bounding_surface->grid_list->first; sgrid_link != NULL; 
	     sgrid_link = sgrid_link->next)
	  {
	    sgrid_ptr = sgrid_link->data;
	    if (sgrid_ptr->surface_label == query_surface &&
		(result = search_quad_tree(xp, yp, zp, sgrid_ptr->quad_tree, 
					   sgrid_ptr, sgrid_ptr->max_n_dist, 
					   query_curve, 0)))
	      {
/* translate the result to the volume grid */
		surface_solution.vgrid_loc = vgrid;
		if (sgrid_ptr->dir == 1)
		  {
		    surface_solution.s_loc = result->r_loc;
		    surface_solution.t_loc = result->s_loc;
		    surface_solution.r_loc = (sgrid_ptr->side==1)? 0.0:1.0;
		  }
		else if (sgrid_ptr->dir == 2)
		  {
		    surface_solution.t_loc = result->r_loc;
		    surface_solution.r_loc = result->s_loc;
		    surface_solution.s_loc = (sgrid_ptr->side==1)? 0.0:1.0;
		  }
		else if (sgrid_ptr->dir == 3)
		  {
		    surface_solution.r_loc = result->r_loc;
		    surface_solution.s_loc = result->s_loc;
		    surface_solution.t_loc = (sgrid_ptr->side==1)? 0.0:1.0;
		  }
		else
		  {
		    printf("Warning: newton: found a patch in the bounding surface "
			   "with dir = %i\n", sgrid_ptr->dir);
		    surface_solution.vgrid_loc = NULL;
		  }
		free(result);
		break; /* one initial guess is enough! */
	      }
	  } /* end for */
	if (surface_solution.vgrid_loc &&
	    (solution = newton(xp, yp, zp, vgrid, vgrid->max_s_dist, &surface_solution, 
			       query_surface, query_curve)))
	  {
/*	    printf("The initial guess from the surface grid worked!\n");*/
	    return solution;
	  }

      } /* end if query_surface != 0 */

/* The newton iteration around the initial guess did not converge. Use the robust */
/* ray + octree methods to invert the mapping */
/* first check if (xp, yp, zp) is inside the grid */
    if (inside_surface(vgrid->bounding_surface, point, query_surface, query_curve, 0))
/* search the oct-tree */
      return search_oct_tree(xp, yp, zp, vgrid->oct_tree, vgrid, vgrid->max_s_dist, 
			     query_surface, query_curve);
    else
      return NULL;
  } /* end if no analytical inversion routine */
  else{
/* call the analytical inversion routine */
    gp.x = xp; gp.y = yp; gp.z = zp;
    inverse_vgrid_mapping( &gp, vgrid );

/* compute the tolerance in the (r,s,t) directions if the cell is on a face with */
/* non-zero surface-label */
/* r1(r)-direction */
    if (query_surface > 0){
/* a solution at a singular point is considered to be invalid! */
      if (fabs(gp.xr) + fabs(gp.yr) + fabs(gp.zr) < 1.e-10 ||
	  fabs(gp.xs) + fabs(gp.ys) + fabs(gp.zs) < 1.e-10 ||
	  fabs(gp.xt) + fabs(gp.yt) + fabs(gp.zt) < 1.e-10)
	return NULL;

      if (query_surface == surf3(1,1) && fabs(gp.r) < vgrid->r1_step ){
	r_dist_low += vgrid->max_s_dist / 
	  sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      }
      else if (query_surface == surf3(2,1) && fabs(1.0-gp.r) < vgrid->r1_step ){
	r_dist_high += vgrid->max_s_dist / 
	  sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
      }
/* r2(s)-direction */
      if (query_surface == surf3(1,2) && fabs(gp.s) < vgrid->r2_step){
	s_dist_low += vgrid->max_s_dist / 
	  sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      }
      else if (query_surface == surf3(2,2) && fabs(1.0-gp.s) < vgrid->r2_step){
	s_dist_high += vgrid->max_s_dist / 
	  sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
      }
/* r3(t)-direction */
      if (query_surface == surf3(1,3) && fabs(gp.t) < vgrid->r3_step){
	t_dist_low += vgrid->max_s_dist / 
	  sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
      }
      else if (query_surface == surf3(2,3) && fabs(1.0-gp.t) < vgrid->r3_step){
	t_dist_high += vgrid->max_s_dist / 
	  sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
      }
    }
  
/* check if the point is inside the (slightly enlarged unit cube in the */
/* parameter plane */
    if (-r_dist_low <= gp.r && gp.r <= 1.0 + r_dist_high && 
	-s_dist_low <= gp.s && gp.s <= 1.0 + s_dist_high &&
	-t_dist_low <= gp.t && gp.t <= 1.0 + t_dist_high){

/* get the cell number */
      i_cell = range3(1,1) + gp.r / vgrid->r1_step;
      j_cell = range3(1,2) + gp.s / vgrid->r2_step;
      k_cell = range3(1,3) + gp.t / vgrid->r3_step;

/* avoid out of bounds */
      i_cell = int_max( range3(1,1), int_min( range3(2,1)-1, i_cell ) );
      j_cell = int_max( range3(1,2), int_min( range3(2,2)-1, j_cell ) );
      k_cell = int_max( range3(1,3), int_min( range3(2,3)-1, k_cell ) );
  
      solution = new_3d_interp( i_cell, j_cell, k_cell, vgrid, gp.r , gp.s , gp.t );
    }
    return solution;
  } /* end if analytical inversion routine */

}
