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
#include <sys/types.h>
#include <sys/times.h>

#include "hole_surface_internal.h"

#define sqr(x) (x)*(x)

/* static prototypes */
static int 
get_edge_label( int i, int j, discrete_surface *sgrid_ptr, real max_dist );
static int
inside_edge_curve(real x, real y, real z, int edge_label, real n_vec[3], real sign, 
		  hole_cutting_curve *hole_curve, int i_seg);
static int
periodic_contiguous(inverse_point *a, inverse_point *b, discrete_surface *sgrid_ptr);
static void
mark_around_segment(real segment[2][3], hole_cutting_curve *hole_curve, int i_seg, 
		    discrete_surface *sgrid_ptr, int sign, int level);
static void
ray_intersect_patch(discrete_surface *sgrid_ptr, subdomain_info *info, real xp[], 
		    int dir, stack_of_ip *stack);
static intersection_point *
ray_intersect_triangle(triangle t, real xp[], real ray[], discrete_surface *sgrid_ptr);
static int 
get_curve( int i, int j, discrete_surface *sgrid_ptr, real max_dist );
static inverse_point *
new_inverse_point(int i_loc, int j_loc, discrete_surface *sgrid_loc, 
		  real r_loc, real s_loc, real n_loc);
static inverse_point *
inside_cell(real xp, real yp, real zp, int i_cell, int j_cell, 
	    discrete_surface *sgrid_ptr, real max_b_dist, int query_curve, 
	    int always_use_max_b_dist );
/* end prototypes */


static inverse_point *
new_inverse_point(int i_loc, int j_loc, discrete_surface *sgrid_loc, 
		  real r_loc, real s_loc, real n_loc)
{
  inverse_point *interp_ptr;
  
  interp_ptr = (inverse_point *) malloc( sizeof(inverse_point) );
  
  /* fill in the interpolation point info later */
  interp_ptr->i_point = -999;
  interp_ptr->j_point = -999;

  /* copy the location info from the arguments */
  interp_ptr->sgrid_loc = sgrid_loc;

  interp_ptr->r_loc = r_loc;
  interp_ptr->s_loc = s_loc;
  interp_ptr->n_loc = n_loc;

  interp_ptr->i_loc = i_loc;
  interp_ptr->j_loc = j_loc;
  
  return interp_ptr;
}

static int 
get_curve( int i, int j, discrete_surface *sgrid_ptr, real max_dist ){
  int this_curve=0;
  real dist_2, left_dist, right_dist, lower_dist, upper_dist;

  dist_2 = max_dist * max_dist;

  left_dist = right_dist = lower_dist = upper_dist = 10.0 * dist_2;

  if (curve(1,1) != 0)
    {
      if (i == 1)
	left_dist = 0.0;
      else
	left_dist = sqr(x(i,j)-x(1,j)) + 
	  sqr(y(i,j)-y(1,j)) + sqr(z(i,j)-z(1,j));
    } /* end if curve(1,1) != 0 */

  if (curve(2,1) != 0)
    {
      if (i == sgrid_ptr->r_dim)
	right_dist = 0.0;
      else
	right_dist = sqr(x(i,j)-x(sgrid_ptr->r_dim,j)) +
	  sqr(y(i,j)-y(sgrid_ptr->r_dim,j)) + sqr(z(i,j)-z(sgrid_ptr->r_dim,j));
    } /* end if curve(2,1) != 0 */

  if (curve(1,2) != 0)
    {
      if (j == 1)
	lower_dist = 0.0;
      else
	lower_dist = sqr(x(i,j)-x(i,1)) +
	  sqr(y(i,j)-y(i,1)) + sqr(z(i,j)-z(i,1));
    } /* end if curve(1,2) != 0 */

  if (curve(2,2) != 0)
    {
      if (j == sgrid_ptr->s_dim)
	upper_dist = 0.0;
      else
	upper_dist = sqr(x(i,j)-x(i,sgrid_ptr->s_dim)) +
	  sqr(y(i,j)-y(i,sgrid_ptr->s_dim)) + sqr(z(i,j)-z(i,sgrid_ptr->s_dim));
    } /* end if curve(2,2) != 0 */
    
  if (left_dist <= dist_2 && curve(1,1) != 0)
    return curve(1,1);
  if (right_dist <= dist_2 && curve(2,1) != 0) 
    return curve(2,1);

  if (lower_dist <= dist_2 && curve(1,2) != 0)
    return curve(1,2);
  if (upper_dist <= dist_2 && curve(2,2) != 0)
    return curve(2,2);

/* if we got this far, the present point is not close to any surface with */
/* non-zero surface label */
  return 0;
}


static int 
get_edge_label( int i, int j, discrete_surface *sgrid_ptr, real max_dist ){
  int this_curve=0;
  real dist_2, left_dist, right_dist, lower_dist, upper_dist;

  dist_2 = max_dist * max_dist;

  left_dist = right_dist = lower_dist = upper_dist = 10.0 * dist_2;

  if (hole_curve(1,1) != 0)
    {
      if (i == 1)
	left_dist = 0.0;
      else
	left_dist = sqr(x(i,j)-x(1,j)) + 
	  sqr(y(i,j)-y(1,j)) + sqr(z(i,j)-z(1,j));
    } /* end if hole_curve(1,1) != 0 */

  if (hole_curve(2,1) != 0)
    {
      if (i == sgrid_ptr->r_dim)
	right_dist = 0.0;
      else
	right_dist = sqr(x(i,j)-x(sgrid_ptr->r_dim,j)) +
	  sqr(y(i,j)-y(sgrid_ptr->r_dim,j)) + sqr(z(i,j)-z(sgrid_ptr->r_dim,j));
    } /* end if hole_curve(2,1) != 0 */

  if (hole_curve(1,2) != 0)
    {
      if (j == 1)
	lower_dist = 0.0;
      else
	lower_dist = sqr(x(i,j)-x(i,1)) +
	  sqr(y(i,j)-y(i,1)) + sqr(z(i,j)-z(i,1));
    } /* end if hole_curve(1,2) != 0 */

  if (hole_curve(2,2) != 0)
    {
      if (j == sgrid_ptr->s_dim)
	upper_dist = 0.0;
      else
	upper_dist = sqr(x(i,j)-x(i,sgrid_ptr->s_dim)) +
	  sqr(y(i,j)-y(i,sgrid_ptr->s_dim)) + sqr(z(i,j)-z(i,sgrid_ptr->s_dim));
    } /* end if hole_curve(2,2) != 0 */
    
  if (left_dist <= dist_2 && hole_curve(1,1) != 0)
    return hole_curve(1,1);
  if (right_dist <= dist_2 && hole_curve(2,1) != 0) 
    return hole_curve(2,1);

  if (lower_dist <= dist_2 && hole_curve(1,2) != 0)
    return hole_curve(1,2);
  if (upper_dist <= dist_2 && hole_curve(2,2) != 0)
    return hole_curve(2,2);

/* if we got this far, the present point is not close to any surface with */
/* non-zero surface label */
  return 0;
}


static int
periodic_contiguous(inverse_point *a, inverse_point *b, discrete_surface *sgrid_ptr)
{
  int i1_dist, i2_dist;
  if (sgrid_ptr->r_period) 
    i1_dist = int_min(abs(a->i_loc - b->i_loc), 
		      sgrid_ptr->r_dim - 1 - abs(a->i_loc - b->i_loc));
  else
    i1_dist = abs(a->i_loc - b->i_loc);

  if (sgrid_ptr->s_period) 
    i2_dist = int_min(abs(a->j_loc - b->j_loc), 
		      sgrid_ptr->s_dim - 1 - abs(a->j_loc - b->j_loc));
  else
    i2_dist = abs(a->j_loc - b->j_loc);

  return (i1_dist <= 1 && i2_dist <= 1);
}

void
cut_surface_holes(input_output *io_ptr, hole_surface *overg)
{
  int sign, ip, ip0, i, j, go_right, go_up, ih, jh, q, q0, present_edge, dir;
  real n_vec[3], segment[2][3], save_n_dist, sp, surf_normal[3];
  inverse_point *point_1, *tmp_point;
  discrete_surface_link *patch_link;
  discrete_surface *sgrid_ptr;
  hole_cutting_curve_link *hole_curve_link;
  hole_cutting_curve *hole_curve;
/* STEP 1 */
/* Mark the cells that are intersected by the hole-cutting curve */
  printf("\nMaking surface holes in `%s'\n", overg->name);
   for (patch_link = overg->grid_list->first; patch_link != NULL;    
        patch_link = patch_link->next)   
    {
      sgrid_ptr = patch_link->data;
      printf("\tInvestigating patch `%s'\n", sgrid_ptr->name);

/* save the normal tolerance, so we can reset it at the end of this loop */
      save_n_dist = sgrid_ptr->max_n_dist;

/* initialize the flag array */
      for (i = 1; i <= sgrid_ptr->r_dim; i++)
	for (j = 1; j <= sgrid_ptr->s_dim; j++)
	  flag(i,j) = sgrid_ptr->priority;

/* check each edge curve */
      for (hole_curve_link = overg->edge_curves->first; hole_curve_link != NULL;
	   hole_curve_link = hole_curve_link->next)
	{
	  hole_curve = hole_curve_link->data;

/* change the normal tolerance temporarily to reflect the roughness of the edge curve */
	  sgrid_ptr->max_n_dist = real_max(sgrid_ptr->max_n_dist, hole_curve->max_dist);

/* check if one or more test points are on the discrete surface sgrid_ptr */
/* if more than one point is found, pick the closest one */
	  point_1 = NULL;
	  for (ip0=0; ip0<2; ip0++) 
	    {
	    for (q0=0; q0<hole_curve->n_points; q0++)
	      if ((tmp_point = search_quad_tree(hole_curve->inside_point[ip0][0][q0],
					      hole_curve->inside_point[ip0][1][q0],
					      hole_curve->inside_point[ip0][2][q0],
					      sgrid_ptr->quad_tree, sgrid_ptr,
					      sgrid_ptr->max_n_dist, 0, 0)))
		{
/* compute the surface normal at tmp_point */
		  compute_cell_normal(surf_normal, tmp_point->i_loc, tmp_point->j_loc, 
				      sgrid_ptr);
/* form the scalar product between the normal at the test point and at tmp_point */
		  for (dir=0, sp=0.0; dir<3; dir++)
		    sp += surf_normal[dir]*hole_curve->inside_normal[ip0][dir][q0];
		  		      
/* NEW: Also require that the surface normal is close to the normal at the test point */
		  if (fabs(sp) > 0.85 && (!point_1 || 
				      fabs(tmp_point->n_loc) < fabs(point_1->n_loc)))
		    {
/* save point_1, ip, and q */		  
		      if (point_1) free(point_1); /* cleanup the previous point_1 */
		      point_1 = tmp_point;
		      ip = ip0; q = q0;
		    } 
		  else
		    {
		      free(tmp_point); /* cleanup */
		    } /* end if */
		} /* end if, for q0... */
	  } /* end for ip0 */
	  if (point_1) 
	    {
/* compute cell normal, edge tangent, and the cross product of the two */
		compute_cell_normal(n_vec, point_1->i_loc, point_1->j_loc, sgrid_ptr);
/* ansatz: sign==1 is inside */
		present_edge = 0;
		if (inside_edge_curve(hole_curve->inside_point[ip][0][q], 
				      hole_curve->inside_point[ip][1][q], 
				      hole_curve->inside_point[ip][2][q], 
				      present_edge, n_vec, 1, hole_curve, q))
		  sign = 1;
		else
		  sign = -1;
		/* cleanup */
		free(point_1);
/* tmp */
/*    		printf("Inside point %i in segment %i of hole-cutting curve %i\n" */
/*  		       "lies on the patch. sign = %i\n",  */
/*    		       ip, q, hole_curve->edge, sign);    */
/*  		printf("Coordinates (%e, %e, %e)\n",   */
/*  		       hole_curve->inside_point[ip][0][q],   */
/*  		       hole_curve->inside_point[ip][1][q],   */
/*  		       hole_curve->inside_point[ip][2][q]);  */
/* end tmp */
	
/* check each segment */
	      printf("\t\tMarking the patch around hole-cutting curve %i\n", 
		     hole_curve->edge);
#ifdef SURFACE_HOLE_CUTTING
      	      if (ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, ogl_length_scale(overg->bb))){ 
    		ogl_marker(hole_curve->inside_point[ip][0][q], 
    			   hole_curve->inside_point[ip][1][q], 
    			   hole_curve->inside_point[ip][2][q], 0.005, /* 0.1 */  
   			   OGL_GOLD); 
      		ogl_end_plot(); 
      	      } 
#endif

/* tmp */
/*    	      printf("There are %i segments\n", hole_curve->n_points-2);    */
/*    	      while( (i = get_int(io_ptr,"Segment number:", 0, 0)) >= 0)    */
/*    		{    */
/*    		  segment[0][0] = hole_curve->x_edge[i];    */
/*    		  segment[0][1] = hole_curve->y_edge[i];    */
/*    		  segment[0][2] = hole_curve->z_edge[i];    */
/*    		  segment[1][0] = hole_curve->x_edge[i+1];  */
/*    		  segment[1][1] = hole_curve->y_edge[i+1];  */
/*    		  segment[1][2] = hole_curve->z_edge[i+1];  */

/*   		  ogl_start_plot(OGL_OLD_PLOT, OVERLAP_W, ogl_length_scale(overg->bb));  */
/*    		  mark_around_segment(segment, hole_curve, i, sgrid_ptr, sign, 0);  */
/*   		  ogl_end_plot();     */
/*    		}    */
/* end tmp */

#ifdef SURFACE_HOLE_CUTTING
 	      ogl_start_plot(OGL_OLD_PLOT, OVERLAP_W, ogl_length_scale(overg->bb)); 
#endif
	      for (i=0; i<=hole_curve->n_points-2; i++) 
/* the last array element is [n_points-1] */
 		{ 
 		  segment[0][0] = hole_curve->x_edge[i]; 
 		  segment[0][1] = hole_curve->y_edge[i]; 
 		  segment[0][2] = hole_curve->z_edge[i]; 
 		  segment[1][0] = hole_curve->x_edge[i+1]; 
 		  segment[1][1] = hole_curve->y_edge[i+1]; 
 		  segment[1][2] = hole_curve->z_edge[i+1]; 
 		  mark_around_segment(segment, hole_curve, i, sgrid_ptr, sign, 0); 
		} /* end for each segment */

#ifdef SURFACE_HOLE_CUTTING
	      ogl_end_plot();
/* get a chance to change the view point! */
	      ogl_wait_for_key("Press <return> to continue:");
#endif
	    } /* end if point_1 */
	} /* end for hole_curve_link... */

/* restore the normal tolerance to its original value */
      sgrid_ptr->max_n_dist = save_n_dist;

      
/* tmp */
/*       printf("r_period = %s, s_period = %s\n", (sgrid_ptr->r_period? "TRUE": "FALSE"),  */
/* 	     (sgrid_ptr->s_period? "TRUE": "FALSE")); */
/*       printf("Flag array after the first stage of the surface hole-cutting\n"); */
/*       print_int_array_2d(sgrid_ptr->flag_ptr); */

    } /* end STEP 1 */
  

  /* STEP 2 */
  /* Fill in the holes */
  printf("\n");
  for (patch_link = overg->grid_list->first; patch_link != NULL; 
       patch_link = patch_link->next)
    {
      sgrid_ptr = patch_link->data;
      printf("Filling the surface holes in patch `%s'\n", sgrid_ptr->name);
      for (j=1; j<=sgrid_ptr->s_dim; j++)
	for (i=1; i<=sgrid_ptr->r_dim; i++)
	  {
	    /* check if it is a hole point */
	    if (flag(i,j) == 'h')
	      {
		/* check the neighbors */
		go_right = (i == 1 || flag(i-1,j) == 'w' || flag(i-1,j) == 'h');
		go_up = (j == 1 || flag(i,j-1) == 'w' || flag(i,j-1) == 'h');
		/* horizontal sweep */
		if (go_right)
		  {
		    for (ih = i+1; ih <= sgrid_ptr->r_dim && flag(ih,j) != 'h' && 
			 flag(ih,j) != 'w'; ih++)
		      flag(ih,j) = 0;
		  }
		else
		  {
		    for (ih = i-1; ih >= 1 && flag(ih,j) != 'h' && flag(ih,j) != 'w';
			 ih--)
		      flag(ih,j) = 0;
		  }
		/* vertical sweep */
		if (go_up)
		  {
		    for (jh = j+1; jh <= sgrid_ptr->s_dim && flag(i,jh) != 'h' && 
			 flag(i,jh) != 'w'; jh++)
		      flag(i,jh) = 0;
		  }
		else
		  {
		    for (jh = j-1; jh >= 1 && flag(i,jh) != 'h' && flag(i,jh) != 'w';
			 jh--)
		      flag(i,jh) = 0;
		  }
	      } /* end if 'h' point */
	  } /* end for all grid points */

      /* change 'w' to the  grid priority and 'h' to 0 */
      for (j=1; j<=sgrid_ptr->s_dim; j++)
	for (i=1; i<=sgrid_ptr->r_dim; i++)
	  {
	    if (flag(i,j) == 'w')
	      flag(i,j) = sgrid_ptr->priority;
	    else if (flag(i,j) == 'h')
	      flag(i,j) = 0;
	  } /* end for all grid points */

   } /* end for all grids */

}

static void
mark_around_segment(real segment[2][3], hole_cutting_curve *hole_curve, int i_seg, 
		    discrete_surface *sgrid_ptr, int sign, int level)
{
  int marked[2] = {0, 0};
  inverse_point *point[2]={NULL, NULL};
  int ip, i, j, i_min, j_min, i_max, j_max, present_edge;
  real n_vec[3], sub_segment[2][3], r_0, r_1, s_0, s_1, radius = 0.005 /* 0.1 */;
  char num_str[10];
/* don't repeat the recursion ab inito */
  if (level > 10) return;

/* check both endpoints of the segment */
  for (ip=0; ip<2; ip++)
    if ((point[ip] = search_quad_tree(segment[ip][0], segment[ip][1], segment[ip][2], 
				      sgrid_ptr->quad_tree, sgrid_ptr,
				      sgrid_ptr->max_n_dist, 0, 0)) )
      {
	marked[ip] = 1;
	/* compute the cell normal */
	compute_cell_normal(n_vec, point[ip]->i_loc, point[ip]->j_loc, sgrid_ptr);

#ifdef SURFACE_HOLE_CUTTING
/* plot the point in black */
   	ogl_marker(segment[ip][0], segment[ip][1], segment[ip][2], radius, OGL_BLACK);
	sprintf(num_str,"%i", i_seg+ip);
	ogl_drawstr(segment[ip], n_vec, num_str);
#endif

/* it might be necessary mark mark more points than just the ones in the cell,
   when the point is very close to (or slightly outside of) the cell boundary */
/* compute the cell boundaries in the parameter plane */
	r_0 = sgrid_ptr->r_step * (point[ip]->i_loc - 1);
	r_1 = sgrid_ptr->r_step * (point[ip]->i_loc);
	s_0 = sgrid_ptr->s_step * (point[ip]->j_loc - 1);
	s_1 = sgrid_ptr->s_step * (point[ip]->j_loc);

/* 	printf("r_0 = %e, r_1 = %e, r_loc = %e\n", r_0, r_1, point[ip]->r_loc); */
/* 	printf("s_0 = %e, s_1 = %e, s_loc = %e\n", s_0, s_1, point[ip]->s_loc); */
	if (point[ip]->r_loc <= r_0 + NEWTON_EPS) 
	  i_min = int_max(1, point[ip]->i_loc - 1);
	else
	  i_min = point[ip]->i_loc;
	if (point[ip]->r_loc >= r_1 - NEWTON_EPS) 
	  i_max = int_min(sgrid_ptr->r_dim, point[ip]->i_loc + 2);
	else
	  i_max = point[ip]->i_loc+1;

	if (point[ip]->s_loc <= s_0 + NEWTON_EPS) 
	  j_min = int_max(1, point[ip]->j_loc - 1);
	else
	  j_min = point[ip]->j_loc;
	if (point[ip]->s_loc >= s_1 - NEWTON_EPS) 
	  j_max = int_min(sgrid_ptr->s_dim, point[ip]->j_loc + 2);
	else
	  j_max = point[ip]->j_loc+1;

/* check each vertex in the enclosing cell */
	for (i=i_min; i<=i_max; i++)
	  for (j=j_min; j<=j_max; j++)
	      /* only necessary to mark each vertex once */
	      if (flag(i, j) == sgrid_ptr->priority)
		{
		  present_edge = get_edge_label(i, j, sgrid_ptr, hole_curve->max_dist);
		  if (inside_edge_curve(x(i,j), y(i,j), z(i,j), present_edge, n_vec, 
					sign, hole_curve, i_seg))
		    {
		      flag(i,j) = 'w'; /* inside of the curve */

#ifdef SURFACE_HOLE_CUTTING
		      ogl_marker(x(i,j), y(i,j), z(i,j), radius, OGL_GREEN); 
#endif
		      if (sgrid_ptr->r_period && i == 1){
/* 			printf("Assigning r-periodic mirror point i=r_dim, j=%i\n", j); */
			flag(sgrid_ptr->r_dim,j) = 'w';
		      }
		      else if (sgrid_ptr->r_period && i == sgrid_ptr->r_dim){
/* 			printf("Assigning r-periodic mirror point i=1, j=%i\n", j); */
			flag(1,j) = 'w';
		      }
		      else if (sgrid_ptr->s_period && j == 1){
/* 			printf("Assigning s-periodic mirror point i=%i, j=s_dim\n", i); */
			flag(i,sgrid_ptr->s_dim) = 'w';
		      }
		      else if (sgrid_ptr->s_period && j == sgrid_ptr->s_dim){
/* 			printf("Assigning s-periodic mirror point i=%i, j=1\n", i); */
			flag(i,1) = 'w';
		      }
		    }
		  else
		    {
		      flag(i,j) = 'h'; /* outside of the curve */
#ifdef SURFACE_HOLE_CUTTING
		      ogl_marker(x(i,j), y(i,j), z(i,j), radius, OGL_RED); 
#endif
		      if (sgrid_ptr->r_period && i == 1){
/* 			printf("Assigning r-periodic mirror point i=r_dim, j=%i\n", j); */
			flag(sgrid_ptr->r_dim,j) = 'h';
		      }
		      else if (sgrid_ptr->r_period && i == sgrid_ptr->r_dim){
/* 			printf("Assigning r-periodic mirror point i=1, j=%i\n", j); */
			flag(1,j) = 'h';
		      }
		      else if (sgrid_ptr->s_period && j == 1){
/* 			printf("Assigning s-periodic mirror point i=%i, j=s_dim\n", i); */
			flag(i,sgrid_ptr->s_dim) = 'h';
		      }
		      else if (sgrid_ptr->s_period && j == sgrid_ptr->s_dim){
/* 			printf("Assigning s-periodic mirror point i=%i, j=1\n", i); */
			flag(i,1) = 'h';
		      }
		    }
		} /* end for each unmarked vertex */
      } /* end if the starting point is inside and separated */

/* recurse if point[0] and point[1] are non-contiguous */
  if (marked[0] && marked[1])
    {
      if (!periodic_contiguous(point[0], point[1], sgrid_ptr))
	{
/* setup sub-segment 1 */
	  for (i=0; i<3; i++)
	    {
	      sub_segment[0][i] = segment[0][i];
	      sub_segment[1][i] = 0.5*(segment[0][i] + segment[1][i]);
	    }
/* recursive call for sub-segment 1 */
	  mark_around_segment(sub_segment, hole_curve, i_seg, sgrid_ptr, sign, level+1);
/* setup sub-segment 2 */
	  for (i=0; i<3; i++)
	    {
	      sub_segment[0][i] = 0.5*(segment[0][i] + segment[1][i]);
	      sub_segment[1][i] = segment[1][i];
	    }
/* recursive call for sub-segment 2 */
	  mark_around_segment(sub_segment, hole_curve, i_seg, sgrid_ptr, sign, level+1);
	}
    } /* end if both point[0] and point[1] are alive */

/* cleanup */
  for (ip=0; ip<2; ip++)
    if (point[ip]) free(point[ip]);

}

static void
edge_tangent(hole_cutting_curve *hole_curve, int i_seg, real t_vec[3])
{
  real length;
  int ip, i;
  t_vec[0] = hole_curve->x_edge[i_seg+1] - hole_curve->x_edge[i_seg];
  t_vec[1] = hole_curve->y_edge[i_seg+1] - hole_curve->y_edge[i_seg];
  t_vec[2] = hole_curve->z_edge[i_seg+1] - hole_curve->z_edge[i_seg];
  length = sqr(t_vec[0]) + sqr(t_vec[1]) + sqr(t_vec[2]); 
  if (length > 1.e-10)
    for (i=0; i<3; i++) t_vec[i] /= length;
  else
    {
      ip = i_seg+2;
      if (ip >= hole_curve->n_points) ip = 1;
      t_vec[0] = hole_curve->x_edge[ip] - hole_curve->x_edge[i_seg];
      t_vec[1] = hole_curve->y_edge[ip] - hole_curve->y_edge[i_seg];
      t_vec[2] = hole_curve->z_edge[ip] - hole_curve->z_edge[i_seg];
      length = sqr(t_vec[0]) + sqr(t_vec[1]) + sqr(t_vec[2]); 
      if (length > 1.e-10)
	for (i=0; i<3; i++) t_vec[i] /= length;
      else
	printf("edge_tangent: Error: undefined tangent for segment %i "
	       "in edge curve %i\n", i_seg, hole_curve->edge);
    }
}

static int
inside_edge_curve(real x, real y, real z, int edge_label, real n_vec[3], real sign, 
		  hole_cutting_curve *hole_curve, int i_seg)
{
/* the edge_label corresponds to the point (x,y,z) */
  real t_vec[3], n_x_t[3], dr[3], r, current_dist, min_dist;
  int direction=0, alternate=0, i_m, i_min, i;
#define dist2(i) sqr(x - hole_curve->x_edge[i]) + \
  sqr(y - hole_curve->y_edge[i]) + sqr(z - hole_curve->z_edge[i])

/* A point with matching edge_label is by definition on the inside */
/* of the cutting curve */
  if (hole_curve->edge == edge_label) 
    return TRUE;

/* find the closest grid point on the edge curve */
  min_dist = dist2(0);
  i_min = 0;
  for (i=1; i<=hole_curve->n_points-2; i++)
    if ((current_dist = dist2(i)) < min_dist)
      {
	i_min = i;
	min_dist = current_dist;
      }
/* find the appropriate segment */
  i_seg = i_min;
  do
    {
/* compute the tangent (length = 1/|x_{i+1} - x_i| to make the computation of r easy)*/
      edge_tangent( hole_curve, i_seg, t_vec );
      r = t_vec[0]*(x - hole_curve->x_edge[i_seg]) +
	t_vec[1]*(y - hole_curve->y_edge[i_seg]) + 
	  t_vec[2]*(z - hole_curve->z_edge[i_seg]);
      if (r < 0.0)
	{
	  if (direction <= 0)
	    {
	      i_seg--;
	      direction = -1;
	      if (i_seg < 0)
		{
		  if (hole_curve->periodic)
		    i_seg = hole_curve->n_points - 2;
		  else
		    {
		      i_seg = 0;
		      break;
		    }
		}
	    }
	  else
	    alternate = TRUE;
	}
      else if (r > 1.0)
	{
	  if (direction >= 0)
	    {
	      i_seg++;
	      direction = 1;
	      if (i_seg > hole_curve->n_points-2)
		{
		  if (hole_curve->periodic)
		    i_seg = 0;
		  else
		    {
		      i_seg = hole_curve->n_points-2; /* could go one step further 
							 in the non-periodic case */
		      break;
		    }
		}
	    }
	  else
	    alternate = TRUE;
	}
    }
  while (!alternate && (r < 0.0 || r > 1.0));

/* tmp */
/*  printf("inside_edge_curve: i_seg = %i\n", i_seg);*/

/* if i_seg alternates, take the average of the two tangents and let i_seg 
   be the vertex in the middle */
  if (alternate)
    {
      i_m = i_seg-1;
      if (i_m < 0)
	{
	  if (hole_curve->periodic)
	    i_m = hole_curve->n_points - 2;
	  else
	    i_m = 0;
	}
      t_vec[0] = hole_curve->x_edge[i_seg+1] - hole_curve->x_edge[i_m];
      t_vec[1] = hole_curve->y_edge[i_seg+1] - hole_curve->y_edge[i_m];
      t_vec[2] = hole_curve->z_edge[i_seg+1] - hole_curve->z_edge[i_m];
    }

  n_x_t[0] = n_vec[1]*t_vec[2] - n_vec[2]*t_vec[1];
  n_x_t[1] = n_vec[2]*t_vec[0] - n_vec[0]*t_vec[2];
  n_x_t[2] = n_vec[0]*t_vec[1] - n_vec[1]*t_vec[0];

  dr[0] = x - hole_curve->x_edge[i_seg];
  dr[1] = y - hole_curve->y_edge[i_seg];
  dr[2] = z - hole_curve->z_edge[i_seg];

  /* check the scalar product of dr and n_x_t */
  if (sign*(dr[0]*n_x_t[0] + dr[1]*n_x_t[1] + dr[2]*n_x_t[2]) > 0.0)
    return TRUE;
  else
    return FALSE;
#undef dist2
}



static intersection_point *
ray_intersect_triangle(triangle t, real xp[], real ray[], 
		       discrete_surface *sgrid_ptr){
/* **************************************************** */
/* The ray is assumed to be normalized (ray*ray = 1)    */
/* **************************************************** */
/* Test if the ray that starts at xp and tends to infinity */
/* intersects the triangle t. If it does, return the pointer to the, */
/* intersetion point, otherwise, NULL. */

  real beta3, n_vec[3], cos_theta, n_len, bn, det;
  real_array_2d *A_, *b_;
  long int pivot[3];

  intersection_point *ip_ptr = NULL;

#define A(i,j) compute_index_2d(A_,i,j)
#define b(i) compute_index_2d(b_,i,1)

/* allocate memory for the matrix and the right hand side */
  A_ = create_real_array_2d( 3, 3 );
  b_ = create_real_array_2d( 3, 1 );

  A(1,1) = t.x[0] - t.x[2];
  A(2,1) = t.y[0] - t.y[2];
  A(3,1) = t.z[0] - t.z[2];

  A(1,2) = t.x[1] - t.x[2];
  A(2,2) = t.y[1] - t.y[2];
  A(3,2) = t.z[1] - t.z[2];

  A(1,3) = -ray[0];
  A(2,3) = -ray[1];
  A(3,3) = -ray[2];

  /* setup the right hand side */
  b(1) = xp[0] - t.x[2];
  b(2) = xp[1] - t.y[2];
  b(3) = xp[2] - t.z[2];
  
  /* Compute the normal */
  n_vec[0] = A(2,1)*A(3,2) - A(3,1)*A(2,2);
  n_vec[1] = A(3,1)*A(1,2) - A(1,1)*A(3,2);
  n_vec[2] = A(1,1)*A(2,2) - A(2,1)*A(1,2);
  n_len = sqrt(sqr(n_vec[0])+sqr(n_vec[1])+sqr(n_vec[2]));

/* check the determinant, i.e., the scalar product of ray and the normal */
/* of the triangle */
  det = n_vec[0]*ray[0] + n_vec[1]*ray[1] + n_vec[2]*ray[2];
  if (fabs(det) > NEWTON_EPS)
    {

      /* compute angle of intersection */
      cos_theta = det/n_len;

      /* solve the system */
      LUsolve(A_, pivot, b_);

      /* intersection? */
      beta3 = 1.0 - b(1) - b(2);
      if ((b(1) >= -NEWTON_EPS && b(1) <= 1.0+NEWTON_EPS) && 
	  (b(2) >= -NEWTON_EPS && b(2) <= 1.0+NEWTON_EPS) && 
	  (beta3 >= -NEWTON_EPS && beta3 <= 1.0+NEWTON_EPS) && 
	  b(3) > 0.0) 
	{
	  ip_ptr = new_intersection_point(b(1)*t.x[0] + b(2)*t.x[1] + beta3*t.x[2],
					  b(1)*t.y[0] + b(2)*t.y[1] + beta3*t.y[2],
					  b(1)*t.z[0] + b(2)*t.z[1] + beta3*t.z[2],
					  fabs(cos_theta), sgrid_ptr);
/* 	  printf("ray_intersect_triangle: xp = (%e, %e, %e) intersects\n" */
/* 		 "the triangle at (%e, %e, %e)\n", xp[0], xp[1], xp[2],  */
/* 		 ip_ptr->ic[0], ip_ptr->ic[1], ip_ptr->ic[2]); */
	}
    }
  else
    {
/* Check if xp lies in the same plane as the triangle */
      bn = b(1)*n_vec[0] + b(2)*n_vec[1] + b(3)*n_vec[2];
      if (fabs(bn) <= NEWTON_EPS)
	{
/* Define the intersection to occur at corner 0 of the triangle. */
/* The intersection angle is zero. */
	  ip_ptr = new_intersection_point(t.x[0], t.y[0], t.z[0], 0.0, sgrid_ptr);
	}
    }

/* cleanup */
  delete_real_array_2d(A_);
  delete_real_array_2d(b_);

  return ip_ptr;
}

static void
set_ray( real ray[], int dir ){
  ray[0] = ray[1] = ray[2] = 0.0;
  if (dir == 1)
    ray[0] = -1.0;
  else if (dir == 2)
    ray[0] = 1.0;
  else if (dir == 3)
    ray[1] = -1.0;
  else if (dir == 4)
    ray[1] = 1.0;
  else if (dir == 5)
    ray[2] = -1.0;
  else if (dir == 6)
    ray[2] = 1.0;
  else
    printf("ERROR: set_ray: Impossible direction: %i\n", dir);
}

static void
ray_intersect_patch(discrete_surface *sgrid_ptr, subdomain_info *info, real xp[], 
		    int dir, stack_of_ip *stack)
{

  triangle test;
  int i, j, inside_bb=FALSE;
  intersection_point *ip_ptr;
  real ray[3];

  if (info == NULL) return;

  /* check if the ray intersects the bounding box */

  if (dir == 1 && info->x_min < xp[0] &&
      info->y_min <= xp[1] && xp[1] <= info->y_max &&
      info->z_min <= xp[2] && xp[2] <= info->z_max)
    inside_bb = TRUE;
  else if (dir == 2 && info->x_max > xp[0] &&
      info->y_min <= xp[1] && xp[1] <= info->y_max &&
      info->z_min <= xp[2] && xp[2] <= info->z_max)
    inside_bb = TRUE;
  else if (dir == 3 && info->y_min < xp[1] &&
      info->x_min <= xp[0] && xp[0] <= info->x_max &&
      info->z_min <= xp[2] && xp[2] <= info->z_max)
    inside_bb = TRUE;
  else if (dir == 4 && info->y_max > xp[1] &&
      info->x_min <= xp[0] && xp[0] <= info->x_max &&
      info->z_min <= xp[2] && xp[2] <= info->z_max)
    inside_bb = TRUE;
  else if (dir == 5 && info->z_min < xp[2] &&
      info->x_min <= xp[0] && xp[0] <= info->x_max &&
      info->y_min <= xp[1] && xp[1] <= info->y_max)
    inside_bb = TRUE;
  else if (dir == 6 && info->z_max > xp[2] &&
      info->x_min <= xp[0] && xp[0] <= info->x_max &&
      info->y_min <= xp[1] && xp[1] <= info->y_max)
    inside_bb = TRUE;



  if (inside_bb)
    {
      /* check if we are on a leaf */
      if (info->i_max - info->i_min == 1 && info->j_max - info->j_min == 1)
	{
	  i = info->i_min;
	  j = info->j_min;
	  set_ray( ray, dir );
	  /* CHECK cell (i,j) if all vertices are alive */
	  if (flag(i,j) != 0 && flag(i+1,j) != 0 && flag(i,j+1) != 0 && 
	      flag(i+1,j+1) != 0)
	    {
	      /* first triangle */
	      test.x[0] = x(i,j);   test.y[0] = y(i,j);   test.z[0] = z(i,j);
	      test.x[1] = x(i+1,j); test.y[1] = y(i+1,j); test.z[1] = z(i+1,j);
	      test.x[2] = x(i,j+1); test.y[2] = y(i,j+1); test.z[2] = z(i,j+1);
	      if ((ip_ptr = ray_intersect_triangle(test, xp, ray, sgrid_ptr)))
		push_ip(ip_ptr, stack);
	      /* second triangle */
	      test.x[0] = x(i+1,j+1); test.y[0] = y(i+1,j+1); test.z[0] = z(i+1,j+1);
	      test.x[1] = x(i+1,j);   test.y[1] = y(i+1,j);   test.z[1] = z(i+1,j);
	      test.x[2] = x(i,j+1);   test.y[2] = y(i,j+1);   test.z[2] = z(i,j+1);
	      if ((ip_ptr = ray_intersect_triangle(test, xp, ray, sgrid_ptr)))
		push_ip(ip_ptr, stack);
	    } /* end if the cell is alive */
	} /* end if leaf */
      else
	/* recurse */
	{
	  ray_intersect_patch(sgrid_ptr, info->low_left,  xp, dir, stack);
	  ray_intersect_patch(sgrid_ptr, info->low_right, xp, dir, stack);
	  ray_intersect_patch(sgrid_ptr, info->upp_left,  xp, dir, stack);
	  ray_intersect_patch(sgrid_ptr, info->upp_right, xp, dir, stack);
	}
      
    } /* end if inside the bounding box */
}

int
inside_surface(hole_surface *overg, real xp[], int present_surface, int present_curve, 
	       int debug){
  discrete_surface_link *sgrid_link;
  discrete_surface *sgrid_ptr;
  int intersect=0, strictly_inside, within_tolerance = FALSE, i
    , tangential_intersection, dir;
  real x_tolerance, ray[3];

  stack_of_ip *stack = NULL;
  intersection_point *ip_ptr, *ip2_ptr;
  inverse_point *result;

/* loop through different directions of the ray until all tangential intersections */
/* are gone */
  dir = 1;
  do{
    if (stack) delete_stack_of_ip(stack);
    stack = new_stack_of_ip();
/* check each patch */
    for (sgrid_link = overg->grid_list->first; sgrid_link != NULL; 
	 sgrid_link = sgrid_link->next)
      {
	sgrid_ptr = sgrid_link->data;
	ray_intersect_patch(sgrid_ptr, sgrid_ptr->quad_tree, xp, dir, stack);
      } /* end for all patches */

/* check for small intersection angles */
    tangential_intersection = FALSE;
    for (ip_ptr = stack->last_ip; ip_ptr != NULL; ip_ptr = ip_ptr->prev_ip)
      {
	if (fabs(ip_ptr->sin_theta) < 0.1) tangential_intersection = TRUE;
      }

#ifdef TEST_INSIDE_GLOBAL_DOMAIN
    if (tangential_intersection)
      {
	set_ray( ray, dir );
	printf("A tangential intersection occured for point (%f, %f, %f) in\n"
	       "the direction (%f, %f, %f)\n", xp[0], xp[1], xp[2], 
	       ray[0], ray[1], ray[2]);
      }
#endif

  } while(tangential_intersection && ++dir <= 6);

  if (tangential_intersection)
    printf("WARNING: The ray from (%f, %f, %f) has an (almost) tangential\n"
	   "intersection with the hole-cutting surface `%s'.\n", 
	   xp[0], xp[1], xp[2], overg->name);
  
#ifdef TEST_INSIDE_GLOBAL_DOMAIN
  if (debug){
    printf("Raw intersection count = %i\n", stack->n_ip); 
    for (ip_ptr = stack->last_ip, i=0; ip_ptr != NULL; ip_ptr = ip_ptr->prev_ip, i++)
      {
	printf("#%i: (%f, %f, %f), sin(alpha) = %f\n", i, ip_ptr->ic[0], ip_ptr->ic[1], 
	       ip_ptr->ic[2], ip_ptr->sin_theta);
      }

  }
#endif

/* remove redundant intersections */ 
  for (ip_ptr = stack->last_ip; ip_ptr != NULL; ip_ptr = ip_ptr->prev_ip)
    {
      if (fabs(ip_ptr->sin_theta) < 0.1) tangential_intersection = TRUE;

      if (!ip_ptr->unique) continue;
      for (ip2_ptr = stack->last_ip; ip2_ptr != NULL; ip2_ptr = ip2_ptr->prev_ip)
	{
	  if (ip2_ptr == ip_ptr || !ip2_ptr->unique) continue;


/* if the intersections occured on the same patch (or patches that originates from the */
/* same vgrid), we can use a small tolerance */
	  if (ip_ptr->sgrid->vgrid_mother == ip2_ptr->sgrid->vgrid_mother)
	    x_tolerance = 10.0*NEWTON_EPS;
	  else
/* Otherwise,  we must use a more generous tolerance */
	    x_tolerance = 
/* The tolerance ought to be max(max_n_dist_1/sin_theta_1, max_n_dist_2/sin_theta_2) */
	      real_max(ip_ptr->sgrid->max_n_dist, ip2_ptr->sgrid->max_n_dist) / 
	      real_min(ip_ptr->sin_theta, ip2_ptr->sin_theta);

/* compare the intersection coordinates */
	  if (fabs(ip2_ptr->ic[0] - ip_ptr->ic[0]) < x_tolerance &&
	      (fabs(ip2_ptr->ic[1] - ip_ptr->ic[1]) + 
	       fabs(ip2_ptr->ic[2] - ip_ptr->ic[2])) < 10.0*NEWTON_EPS)
	    ip2_ptr->unique = FALSE;
	}
    }
  
/* count the number of unique intersection points */
  intersect = 0;
  for (ip_ptr = stack->last_ip; ip_ptr != NULL; ip_ptr = ip_ptr->prev_ip)
    {
      intersect += ip_ptr->unique;
    }

#ifdef TEST_INSIDE_GLOBAL_DOMAIN
  if (debug){
    printf("#uniqe intersections = %i\n", intersect);
    for (ip_ptr = stack->last_ip, i=0; ip_ptr != NULL; ip_ptr = ip_ptr->prev_ip, i++){
      if (!ip_ptr->unique) continue;
      printf("#%i: (%f, %f, %f), sin(alpha) = %f\n", i, ip_ptr->ic[0], ip_ptr->ic[1], 
	     ip_ptr->ic[2], ip_ptr->sin_theta);
    }
  }
#endif

/* cleanup */
  delete_stack_of_ip( stack );
  strictly_inside = (intersect % 2);
  
/* account for points that are less than the tolerance outside of the bounding surface */
  if (!strictly_inside && present_surface != 0)
    {
/* check each patch */
      for (sgrid_link = overg->grid_list->first; sgrid_link != NULL; 
	   sgrid_link = sgrid_link->next)
	{
	  sgrid_ptr = sgrid_link->data;
	  if (sgrid_ptr->surface_label == present_surface &&
	      (result = search_quad_tree(xp[0], xp[1], xp[2], sgrid_ptr->quad_tree, 
					 sgrid_ptr, sgrid_ptr->max_n_dist, 
					 present_curve, 0)))
	    {
	      within_tolerance = TRUE;
	      free(result);
	      break;
	    }
	} /* end for */
    } /* end if */
  
  return (strictly_inside || within_tolerance);
}


static inverse_point *
inside_cell(real xp, real yp, real zp, int i_cell, int j_cell, 
	    discrete_surface *sgrid_ptr, real max_b_dist,
	    int query_curve, int always_use_max_b_dist ){
  inverse_point *interp_ptr=NULL;
  int iter, info;
  const int max_iter=20;
  grid_point gp;
/* normal and temporary stuff*/
  real normal[3], normal_r[3], normal_s[3], r_loc, s_loc, n_loc, r_0, r_1, s_0, s_1
    , r_dist_low, r_dist_high, s_dist_low, s_dist_high, res, r_prev, s_prev
      , n_prev, prev_res, diff, prev_diff;
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

/* initial guess */
  r_prev = r_loc = sgrid_ptr->r_step * (i_cell - 1 + 0.5); 
  s_prev = s_loc = sgrid_ptr->s_step * (j_cell - 1 + 0.5);
  n_prev = n_loc = 0.0;

/* get the iteration going */
  prev_diff = 1.e10;

  iter = 0;
  do{
/* evaluate the surface */
    gp.r = r_loc; gp.s = s_loc;
    bi_linear(&gp, sgrid_ptr);
/* compute the normalized normal */
    compute_continuous_normal(normal, normal_r, normal_s, r_loc, s_loc, sgrid_ptr);

/* compute the residual */
    b(1) = xp + n_loc*normal[0] - gp.x;
    b(2) = yp + n_loc*normal[1] - gp.y;
    b(3) = zp + n_loc*normal[2] - gp.z;

    if (( res = sqrt( b(1)*b(1) + b(2)*b(2) + b(3)*b(3) )) <= NEWTON_EPS) break;

/* assign the matrix */
    A(1,1) = gp.xr - n_loc*normal_r[0]; A(1,2) = gp.xs - n_loc*normal_s[0]; 
    A(2,1) = gp.yr - n_loc*normal_r[1]; A(2,2) = gp.ys - n_loc*normal_s[1]; 
    A(3,1) = gp.zr - n_loc*normal_r[2]; A(3,2) = gp.zs - n_loc*normal_s[2]; 

    A(1,3) = -normal[0];
    A(2,3) = -normal[1];
    A(3,3) = -normal[2];    

/* Solve the system */
    if ((info=LUsolve(A_, pivot, b_)))
      {
	if (info > 1)
	  printf("inside_cell: Singular jacobian!\n");
	else
	  printf("inside_cell: Invalid matrix!\n");
	return NULL;
      }

/* save previous guess and residual in case the size of the increment increases */
    r_prev = r_loc; s_prev = s_loc; n_prev = n_loc;
    prev_res = res;

/* update */
    r_loc = r_loc + b(1);
    s_loc = s_loc + b(2);
    n_loc = n_loc + b(3);

/* norm of increment */
    diff = b(1)*b(1) + b(2)*b(2) + b(3)*b(3);

/* make sure the increment is decreasing at least a little */
    if (diff > 0.95*prev_diff){
/* hypothesis: The iteration diverges because the point is outside of the cell */
      return NULL;
    }
    prev_diff = diff;

  } while (++iter < max_iter);

  if (res > NEWTON_EPS)
    printf("Warning: Poor convergence in Newton iteration in `inside_cell' for "
	   "cell (%i, %i) in grid `%s'.\n"
	   "After %i iterations, the residual was = %e\n", i_cell, j_cell, 
	   sgrid_ptr->name, iter, res);

/* standard tolerance */
  r_dist_low = r_dist_high = sgrid_ptr->max_t_dist /  
    sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
  s_dist_low = s_dist_high = sgrid_ptr->max_t_dist / 
    sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );

/* tolerance in the (r,s) directions for the tangent plane if */
/* the cell is a boundary cell */
  if (query_curve > 0 || always_use_max_b_dist){
    if ((i_cell == 1 && query_curve == curve(1,1)) || always_use_max_b_dist){
      r_dist_low = max_b_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
    }
    else if ((i_cell == sgrid_ptr->r_dim-1 && query_curve == curve(2,1))
	     || always_use_max_b_dist){
      r_dist_high = max_b_dist / sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
    }
    if ((j_cell == 1 && query_curve == curve(1,2))
	|| always_use_max_b_dist){
      s_dist_low = max_b_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
    else if ((j_cell == sgrid_ptr->s_dim-1 && query_curve == curve(2,2))
	     || always_use_max_b_dist){
      s_dist_high = max_b_dist / sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    }
  }

  
/* compute the cell boundaries in the parameter plane */
  r_0 = sgrid_ptr->r_step * (i_cell - 1) - r_dist_low; 
  r_1 = sgrid_ptr->r_step * (i_cell + 1 - 1) + r_dist_high; 
  s_0 = sgrid_ptr->s_step * (j_cell - 1) - s_dist_low;
  s_1 = sgrid_ptr->s_step * (j_cell + 1 - 1) + s_dist_high;

/* check if the point is inside the cell in the tangent plane */
  if (r_0 <= r_loc && r_loc <= r_1 && 
      s_0 <= s_loc && s_loc <= s_1){
/* check if the point is sufficiently close in the normal direction */
    if (fabs(n_loc) <= sgrid_ptr->max_n_dist)
/* create an interpolation point and store the location */
      interp_ptr = new_inverse_point(i_cell, j_cell, sgrid_ptr,
				    r_loc, s_loc, n_loc);
/*     else */
/*       printf("The projection of the point (%f, %f, %f) onto the\n" */
/* 	     "tangentplane is inside cell (%i, %i) in component `%s',\n" */
/* 	     "but the normal distance %e exceeds the tolerance %e\n",  */
/*        xp, yp, zp, i_cell, j_cell, sgrid_ptr->name, n_loc, sgrid_ptr->max_n_dist); */
  }
/*   else */
/*     printf("The projection of the point (%f, %f, %f) onto the\n" */
/* 	   "tangentplane is outside cell (%i, %i) in component `%s'\n", */
/* 	   xp, yp, zp, i_cell, j_cell, sgrid_ptr->name); */

  return interp_ptr;
}

inverse_point *
search_quad_tree(real xp, real yp, real zp, subdomain_info *info, 
		 discrete_surface *sgrid_ptr, real max_b_dist,
		 int query_curve, int always_use_max_b_dist ){
  inverse_point *interp_ptr=NULL;
  real max_n_dist;

  if (info == NULL) return NULL;

/* check if we are inside the (slightly enlarged) bounding box */
  max_n_dist = sgrid_ptr->max_n_dist;
  if (info->x_min - max_n_dist <= xp && xp < info->x_max + max_n_dist &&
      info->y_min - max_n_dist <= yp && yp < info->y_max + max_n_dist &&
      info->z_min - max_n_dist <= zp && zp < info->z_max + max_n_dist){
/* check if we are on a leaf */
    if (info->i_max - info->i_min == 1 && info->j_max - info->j_min == 1){
/* perform a closer investigation to see if the point really is in the cell */
      interp_ptr = inside_cell( xp, yp, zp, info->i_min, info->j_min, 
			       sgrid_ptr, max_b_dist, query_curve, 
			       always_use_max_b_dist );
    }
/* recurse */
    else{
      if ((interp_ptr = search_quad_tree( xp, yp, zp, info->low_left, sgrid_ptr, 
					 max_b_dist, query_curve, 
					 always_use_max_b_dist )) != NULL)
	return interp_ptr;
      if ((interp_ptr = search_quad_tree( xp, yp, zp, info->low_right, sgrid_ptr, 
					 max_b_dist, query_curve, 
					 always_use_max_b_dist)) != NULL)
	return interp_ptr;
      if ((interp_ptr = search_quad_tree( xp, yp, zp, info->upp_left, sgrid_ptr, 
					 max_b_dist, query_curve, 
					 always_use_max_b_dist )) != NULL)
	return interp_ptr;
      if ((interp_ptr = search_quad_tree( xp, yp, zp, info->upp_right, sgrid_ptr, 
					 max_b_dist, query_curve, 
					 always_use_max_b_dist)) != NULL)
	return interp_ptr;
    }
  }

/* return either NULL or the newly created interpolation point */
  return interp_ptr;
}
