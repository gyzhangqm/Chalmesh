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

/* private member functions */
static void
swap_points(int p, int q, hole_cutting_curve *hole_curve);
static int
sub_sort(hole_cutting_curve *hole_curve, real length_scale);
static void
sort_coordinates(hole_cutting_curve *hole_curve);
static subdomain_info *
new_subdomain_info(int i_min, int i_max, int j_min, int j_max, 
		   discrete_surface *sgrid_ptr); 
static subdomain_info *
delete_subdomain_info( subdomain_info *info );
static void
update_hs_bb( hole_surface *patch );
/* end private member functions */


hole_surface *
new_hole_surface( char *name, int surface_label )
{
  hole_surface *overg;

  overg = (hole_surface *) malloc(sizeof(hole_surface));
  
  /* surface label */
  overg->surface_label = surface_label;

  overg->inside = 0;

  /* pointers to the linked list of component surface grids */
  overg->grid_list = new_linked_list();

  /* tolerance for two curves or two surfaces to be regarded as the same curve */
  /* or surface */
  overg->max_b_dist    = 1.e-7;
    
  /* empty list of hole cutting curves */
  overg->edge_curves = new_linked_list();

  /* allocate space for the bounding box */
  overg->bb = (bounding_box *) malloc( sizeof(bounding_box) );

  /* initialize bounding box */
  overg->bb->x_min = 0.0;
  overg->bb->x_max = 1.0;
  overg->bb->y_min = 0.0;
  overg->bb->y_max = 1.0;
  overg->bb->z_min = 0.0;
  overg->bb->z_max = 1.0;

  /* copy the name */
  overg->name = (char *) malloc( (strlen(name)+1)*sizeof(char) );
  strcpy( overg->name, name );

  /* plot grid boundaries and the coordinate cage */
  overg->plot_mode = 4 | 256;

  return overg;
}

hole_surface *
delete_hole_surface( hole_surface *overg )
{
  linked_list_member *this_link;
  discrete_surface *sgrid_ptr;

  if (overg == NULL) return NULL;

  /* free the bounding box */
  free( overg->bb );

  /* delete the discrete surfaces */
  for (this_link = overg->grid_list->first; this_link != NULL; 
       this_link = this_link->next)
    {
      sgrid_ptr = this_link->data;
      delete_discrete_surface( sgrid_ptr );
    }

  /* delete the linked list with surface grids */
  delete_linked_list( overg->grid_list );

  /* delete the hole cutting curves here */
  for (this_link = overg->edge_curves->first; this_link != NULL; 
       this_link = this_link->next)
    delete_hole_curve( this_link->data );

  /* delete the linked list of hole_cutting_curves */
  delete_linked_list( overg->edge_curves );

/* free the name string */
  free(overg->name);

  free( overg );
  return NULL;
}

hole_cutting_curve *
new_hole_curve(real *x_edge, real *y_edge, real *z_edge, int n1, 
	       real *inside_point[2][3], real *inside_normal[2][3], int edge)
{
  hole_cutting_curve *hole_curve;
  int i, ip;

  hole_curve = (hole_cutting_curve *) malloc( sizeof(hole_cutting_curve) );

  hole_curve->n_points = n1;
  hole_curve->x_edge = x_edge;
  hole_curve->y_edge = y_edge;
  hole_curve->z_edge = z_edge;
  for (ip=0; ip<2; ip++)
    for (i=0; i<3; i++){
      hole_curve->inside_point[ip][i]  = inside_point[ip][i];
      hole_curve->inside_normal[ip][i] = inside_normal[ip][i];
    }
  hole_curve->edge = edge;
  hole_curve->periodic = FALSE;
  hole_curve->max_dist = 0.0;
  return hole_curve;
}

hole_cutting_curve *
delete_hole_curve(hole_cutting_curve *hole_curve)
{
  int ip, q;
  free(hole_curve->x_edge);
  free(hole_curve->y_edge);
  free(hole_curve->z_edge);
  for (ip=0; ip<2; ip++)
    for (q=0; q<3; q++){
      free(hole_curve->inside_point[ip][q]);
      free(hole_curve->inside_normal[ip][q]);
    }

  free(hole_curve);
  return NULL;
}

void
insert_hole_curve(hole_surface *overg, hole_cutting_curve *hole_curve )
{
  new_link(overg->edge_curves)->data = hole_curve;
}

#define sqr(x) ((x)*(x))
static void
hole_curve_eps( hole_cutting_curve *hc )
{
  real dx[3], t[3], length, t_dist, n_dist;
  int i, i_max;

  hc->max_dist = 0.0;
  i_max = 0;
  for (i=1; i<hc->n_points-1; i++)
    {
      t[0] = hc->x_edge[i+1] - hc->x_edge[i-1];
      t[1] = hc->y_edge[i+1] - hc->y_edge[i-1];
      t[2] = hc->z_edge[i+1] - hc->z_edge[i-1];
      length = sqrt(sqr(t[0]) + sqr(t[1]) + sqr(t[2]));
      if (length < 1.e-10) continue;
      t[0] /= length;
      t[1] /= length;
      t[2] /= length;

      dx[0] = hc->x_edge[i+1] - hc->x_edge[i];
      dx[1] = hc->y_edge[i+1] - hc->y_edge[i];
      dx[2] = hc->z_edge[i+1] - hc->z_edge[i];

      t_dist = sqrt(sqr(t[0]*dx[0])+sqr(t[1]*dx[1])+sqr(t[2]*dx[2]));
      n_dist = sqrt(sqr(dx[0]-t[0]*t_dist) + sqr(dx[1]-t[1]*t_dist) +
		    sqr(dx[2]-t[2]*t_dist));

      if (hc->max_dist < n_dist)
	{
	  hc->max_dist = n_dist;
	  i_max = i;
	}
    }
  printf("Mismatch tolerance for edge curve %i estimated to %e at\n"
	 "point %i, n_points = %i.\n", 
	 hc->edge, hc->max_dist, i_max, hc->n_points);
}

void
merge_hole_curves(hole_surface *overg)
{
  hole_cutting_curve_link *edge_curve_link, *next_victim;
  hole_cutting_curve *edge_curve, *one_curve=NULL, *new_curve;
  int i, j, n_points, ie, q, ip;
  real *x_edge, *y_edge, *z_edge, *inside_point[2][3], *inside_normal[2][3];
/* assume that there are at most 10 different edge curves in each hole-surface */
  int edge_value[10], n_edge_values=0;

  printf("\nMerging edge curves in hole surface `%s'. There are %i edge(s).\n", 
	 overg->name, overg->edge_curves->n_members);

  if (overg->edge_curves->n_members <= 1)
    {
/* just estimate the mismatch tolerance and quit */
      if (overg->edge_curves->first) 
	hole_curve_eps(overg->edge_curves->first->data);
      return;
    }

/* record all edge_labels in the hole_surface */
  for (edge_curve_link = overg->edge_curves->first; edge_curve_link != NULL;
       edge_curve_link = edge_curve_link->next)
    {
      edge_curve = edge_curve_link->data;
/* check if the edge value is present in the array */
      for (ie=0; ie<n_edge_values; ie++)
	if (edge_value[ie] == edge_curve->edge) break;
/* record the new value */
      if (ie == n_edge_values)
	{
	  if (n_edge_values < 10)
	    {
	      printf("Recorded a new edge value = %i\n", edge_curve->edge);
	      edge_value[n_edge_values++] = edge_curve->edge;
	    }
	  else
	    {
	      printf("Warning: found more than 10 edge values in the hole surface %s.\n"
		     "Only taking the first 10 into account.\n", overg->name);
	      break;
	    }
	} /* end if new edge value */
    }

/* loop over all edge values */
  for (ie = 0; ie < n_edge_values; ie++)
    {
      printf("Doing edge value %i\n", edge_value[ie]);
      n_points = 0;
      for (edge_curve_link = overg->edge_curves->first; edge_curve_link != NULL;
	   edge_curve_link = edge_curve_link->next)
	{
	  edge_curve = edge_curve_link->data;
	  if (edge_curve->edge == edge_value[ie])
	    n_points += edge_curve->n_points;
	}
      /* allocate memory for all points (one extra to put a copy of the first point) */
      x_edge = (real *) malloc( (n_points+1)*sizeof(real) );
      y_edge = (real *) malloc( (n_points+1)*sizeof(real) );
      z_edge = (real *) malloc( (n_points+1)*sizeof(real) );
      for (ip=0; ip<2; ip++)
	for (q=0; q<3; q++){
	  inside_point[ip][q]  = (real *) malloc( (n_points+1)*sizeof(real) );
	  inside_normal[ip][q] = (real *) malloc( (n_points+1)*sizeof(real) );
	}

      /* copy the values */
      for (edge_curve_link = overg->edge_curves->first, i=0; edge_curve_link != NULL;
	   edge_curve_link = edge_curve_link->next)
	{
	  edge_curve = edge_curve_link->data;
	  if (edge_curve->edge == edge_value[ie])
	    {
	      one_curve = edge_curve;
	      for (j=0; j<edge_curve->n_points; j++)
		{
		  x_edge[i] = edge_curve->x_edge[j];
		  y_edge[i] = edge_curve->y_edge[j];
		  z_edge[i] = edge_curve->z_edge[j];
		  for (ip=0; ip<2; ip++)
		    for (q=0; q<3; q++){
		      inside_point[ip][q][i]  = edge_curve->inside_point[ip][q][j];
		      inside_normal[ip][q][i] = edge_curve->inside_normal[ip][q][j];
		    }
		  i++;
		}
	    }
	} /* end for all edge curves */

      /* make a new hole_curve, copy the edge information from the one_curve */
      new_curve = new_hole_curve(x_edge, y_edge, z_edge, n_points, 
				 inside_point, inside_normal, one_curve->edge);

      /* delete the old hole curves with the present edge value and */
      /* remove them from the list */
      edge_curve_link = overg->edge_curves->first;
      while ( edge_curve_link != NULL )
	{
	  next_victim = edge_curve_link->next;
	  edge_curve = edge_curve_link->data;
	  if (edge_curve->edge == edge_value[ie])
	    {
	      /* free the edge_curve */
	      delete_hole_curve(edge_curve_link->data);
	      /* free the link and fix up the pointers */
	      delete_link(edge_curve_link, overg->edge_curves );
	    }
	  edge_curve_link = next_victim;
	} /* end while */

      /* insert the new hole curve into the list */
      insert_hole_curve( overg, new_curve );

      /* sort the nodes in the new hole curve */
      sort_coordinates( new_curve );

      hole_curve_eps( new_curve );
      /* tmp */
      /*   printf("The sorted hole-cutting curve for edge = %i\n", new_curve->edge); */
      /*   for (i=0; i<new_curve->n_points; i++) */
      /*     { */
      /*       printf("Point %i: (%e, %e, %e)\n", i, x_edge[i], y_edge[i], z_edge[i]); */
      /*     } */
      /*   ogl_wait_for_key("Press RETURN to continue:"); */

    } /* end for all edge values */
}

static real
dist2(int i, int j, real t[3], hole_cutting_curve *hole_curve)
{
  real t_dist, n_dist, d, d2;
  t_dist = t[0] * (hole_curve->x_edge[j] - hole_curve->x_edge[i]) +
    t[1] * (hole_curve->y_edge[j] - hole_curve->y_edge[i]) +
      t[2] * (hole_curve->z_edge[j] - hole_curve->z_edge[i]);
  d2 = sqr(hole_curve->x_edge[j] - hole_curve->x_edge[i]) +
    sqr(hole_curve->y_edge[j] - hole_curve->y_edge[i]) +
      sqr(hole_curve->z_edge[j] - hole_curve->z_edge[i]);
  n_dist = sqrt( d2 - sqr(t_dist) );

/* Elliptic weighted distance. */
  if (t_dist >= 0.0)
    {
/* Normal semi-axis 1, tangential semi-axis 5. */
      d = sqr( t_dist/5.0 ) + sqr( n_dist ) ;
    }
  else
    {
/* Normal semi-axis 1, tangential semi-axis 1/5.0. */
      d = sqr( t_dist*5.0 ) + sqr( n_dist ) ;
    }
  return d;
}

static void
sort_coordinates( hole_cutting_curve *hole_curve )
{
  int i_p, q, ic;
  real x_min, x_max
    , y_min, y_max, z_min, z_max, length_scale;

/* Get the length scale of the curve to determine if the edge curve is periodic or not */
  x_max = x_min = hole_curve->x_edge[0];
  y_max = y_min = hole_curve->y_edge[0];
  z_max = z_min = hole_curve->z_edge[0];
  for (i_p=1; i_p<hole_curve->n_points; i_p++)
    {
      if (hole_curve->x_edge[i_p] < x_min) x_min = hole_curve->x_edge[i_p];
      if (hole_curve->x_edge[i_p] > x_max) x_max = hole_curve->x_edge[i_p];
      if (hole_curve->y_edge[i_p] < y_min) y_min = hole_curve->y_edge[i_p];
      if (hole_curve->y_edge[i_p] > y_max) y_max = hole_curve->y_edge[i_p];
      if (hole_curve->z_edge[i_p] < z_min) z_min = hole_curve->z_edge[i_p];
      if (hole_curve->z_edge[i_p] > z_max) z_max = hole_curve->z_edge[i_p];
    }
  length_scale = sqrt(sqr(x_max-x_min) + sqr(y_max-y_min) + sqr(z_max-z_min));

  if ((i_p = sub_sort(hole_curve, length_scale)) > 0)
    {
/* swap i_p and 0 */
      swap_points(i_p, 0, hole_curve);

/* swap i_p-1 and 1 */
      swap_points(i_p-1, 1, hole_curve);

/* redo the sort */
      printf("Redoing the sort\n");
      sub_sort(hole_curve, length_scale);
    }

/* check the distance between the first and last points */
  i_p = hole_curve->n_points - 1;
  hole_curve->periodic = (sqrt(sqr(hole_curve->x_edge[i_p] - hole_curve->x_edge[0]) +
			       sqr(hole_curve->y_edge[i_p] - hole_curve->y_edge[0]) +
			       sqr(hole_curve->z_edge[i_p] - hole_curve->z_edge[0]))
			  <= 10*length_scale/hole_curve->n_points);

  printf("The edge curve %i seems to be %s\n", hole_curve->edge,
	 (hole_curve->periodic? "periodic": "NON-periodic"));

/* put a copy of the first point last (there is space for it) */
  if (hole_curve->periodic)
    {
      i_p = hole_curve->n_points++;
      hole_curve->x_edge[i_p] = hole_curve->x_edge[0];
      hole_curve->y_edge[i_p] = hole_curve->y_edge[0];
      hole_curve->z_edge[i_p] = hole_curve->z_edge[0];

      for (q=0; q<2; q++)
	for (ic=0; ic<3; ic++){
	  hole_curve->inside_point[q][ic][i_p]  = hole_curve->inside_point[q][ic][0];
	  hole_curve->inside_normal[q][ic][i_p] = hole_curve->inside_normal[q][ic][0];
	}
    }
}

static int
sub_sort(hole_cutting_curve *hole_curve, real length_scale)
{
  int i_p, i, j, i_swap;
  real min_dist, t[3], t_length, new_dist;
  /* compute the first local tangent and normal */
  t[0] = hole_curve->x_edge[1] - hole_curve->x_edge[0];
  t[1] = hole_curve->y_edge[1] - hole_curve->y_edge[0];
  t[2] = hole_curve->z_edge[1] - hole_curve->z_edge[0];
  t_length = sqrt(sqr(t[0]) + sqr(t[1]) + sqr(t[2]));
  /* normalize */
  for (i=0; i<3; i++) t[i] /= t_length;
  
  for (i_p=0; i_p<hole_curve->n_points-1; i_p++)
    {
      /* look for the closest point in the rest of the array */
      min_dist = dist2(i_p, i_p+1, t, hole_curve);
      i_swap = i_p+1;
      for (j=i_p+2; j<hole_curve->n_points; j++)
	{
	  if ((new_dist = dist2(i_p, j, t, hole_curve)) < min_dist)
	    {
	      min_dist = new_dist;
	      i_swap = j;
	    }
	} /* end for j */
      /* i_p is an end point if min_dist > 10*length_scale/n_points */
      if (sqrt(min_dist) > 10*length_scale/hole_curve->n_points)
	{
	  printf("Warning: found an internal end point in edge curve %i at index %i\n", 
		 hole_curve->edge, i_p);
	  return i_p;
	}
      /* swap i_p+1 and i_swap */
      swap_points(i_p+1, i_swap, hole_curve);

      /* only compute a new tangent if the distance was finite */
      if (min_dist > 0.0)
	{
	  /* compute next local tangent */
	  t[0] = hole_curve->x_edge[i_p+1] - hole_curve->x_edge[i_p];
	  t[1] = hole_curve->y_edge[i_p+1] - hole_curve->y_edge[i_p];
	  t[2] = hole_curve->z_edge[i_p+1] - hole_curve->z_edge[i_p];
	  t_length = sqrt(sqr(t[0]) + sqr(t[1]) + sqr(t[2]));
	  /* normalize */
	  for (i=0; i<3; i++) t[i] /= t_length;
	}
    } /* end for i_p */
  return 0;
}

#undef sqr

static void
swap_points(int p, int q, hole_cutting_curve *hole_curve)
{
  real x_tmp, y_tmp, z_tmp;
  int ip, ic;

  if (p != q)
    {
      x_tmp = hole_curve->x_edge[q];
      y_tmp = hole_curve->y_edge[q];
      z_tmp = hole_curve->z_edge[q];

      hole_curve->x_edge[q] = hole_curve->x_edge[p];
      hole_curve->y_edge[q] = hole_curve->y_edge[p];
      hole_curve->z_edge[q] = hole_curve->z_edge[p];
      hole_curve->x_edge[p] = x_tmp;
      hole_curve->y_edge[p] = y_tmp;
      hole_curve->z_edge[p] = z_tmp;

/* also swap the inside points & normals */
      for (ip=0; ip<2; ip++)
	for (ic=0; ic<3; ic++)
	  {
	    x_tmp = hole_curve->inside_point[ip][ic][q];
	    hole_curve->inside_point[ip][ic][q] = hole_curve->inside_point[ip][ic][p];
	    hole_curve->inside_point[ip][ic][p] = x_tmp;

	    y_tmp = hole_curve->inside_normal[ip][ic][q];
	    hole_curve->inside_normal[ip][ic][q] = hole_curve->inside_normal[ip][ic][p];
	    hole_curve->inside_normal[ip][ic][p] = y_tmp;
	  }
    }
}

void
insert_discrete_surface(hole_surface *overg, discrete_surface *sgrid_ptr)
{
/* insert the discrete surface into the list of surfaces */
  new_link(overg->grid_list)->data = sgrid_ptr;
/* update the bounding box for the overg */
  update_hs_bb( overg ); 
/* assign the priority */
  sgrid_ptr->priority = overg->grid_list->n_members;
}

discrete_surface *
new_discrete_surface(char *name, real_array_2d *x_, real_array_2d *y_, 
		     real_array_2d *z_, int present_surface, struct component_vgrid *mother)
{
  real n_vec[3], i_dist, l_scale, r_gap, s_gap, middle[3];
  discrete_surface *sgrid_ptr;
  int i, j, n1, n2;
#define x_in(i,j) compute_index_2d(x_, i, j)
#define y_in(i,j) compute_index_2d(y_, i, j)
#define z_in(i,j) compute_index_2d(z_, i, j)

/* check that the dimensions of x, y and z make sense */
  if (! (x_->n1 > 1 && x_->n2 > 1 && 
	 x_->n1 == y_->n1 && x_->n2 == y_->n2 && 
	 z_->n1 == y_->n1 && z_->n2 == y_->n2))
    {
      printf("ERROR: new_discrete_surface: dimension mismatch\n");
      return NULL;
    }

  sgrid_ptr = (discrete_surface *) malloc(sizeof(discrete_surface));
  
  /* initialize all fields */

  /* copy the name */
  sgrid_ptr->name = (char *) malloc( (strlen(name)+1)*sizeof(char) );
  sgrid_ptr->name = strcpy( sgrid_ptr->name, name );

  sgrid_ptr->priority = 0;

  n1 = sgrid_ptr->r_dim = x_->n1;
  n2 = sgrid_ptr->s_dim = x_->n2;

  sgrid_ptr->bb = (bounding_box *) malloc( sizeof(bounding_box) );
  sgrid_ptr->bb->x_min = 1.e4;
  sgrid_ptr->bb->x_max = -1.e4;
  sgrid_ptr->bb->y_min = 1.e4;
  sgrid_ptr->bb->y_max = -1.e4;
  sgrid_ptr->bb->z_min = 1.e4;
  sgrid_ptr->bb->z_max = -1.e4;
  for (i=1; i<=x_->n1; i++)
    for (j=1; j<=x_->n2; j++)
      {
	sgrid_ptr->bb->x_min = real_min(sgrid_ptr->bb->x_min, x_in(i,j));
	sgrid_ptr->bb->x_max = real_max(sgrid_ptr->bb->x_max, x_in(i,j));
	sgrid_ptr->bb->y_min = real_min(sgrid_ptr->bb->y_min, y_in(i,j));
	sgrid_ptr->bb->y_max = real_max(sgrid_ptr->bb->y_max, y_in(i,j));
	sgrid_ptr->bb->z_min = real_min(sgrid_ptr->bb->z_min, z_in(i,j));
	sgrid_ptr->bb->z_max = real_max(sgrid_ptr->bb->z_max, z_in(i,j));
      }

#define sqr(x) (x)*(x)
  l_scale = sqrt(sqr(sgrid_ptr->bb->x_max-sgrid_ptr->bb->x_min) +
		 sqr(sgrid_ptr->bb->y_max-sgrid_ptr->bb->y_min) +
		 sqr(sgrid_ptr->bb->z_max-sgrid_ptr->bb->z_min));

/* determine periodicity by computing the gap between the first and last lines */
  r_gap = 0.0;
  for (j=1; j<=n2; j++)
    r_gap += fabs(x_in(1,j)-x_in(n1,j)) + 
      fabs(y_in(1,j)-y_in(n1,j)) + 
	fabs(z_in(1,j)-z_in(n1,j));
/*     r_gap += fabs(x_in(2,j)-x_in(n1,j)) +  */
/*       fabs(y_in(2,j)-y_in(n1,j)) +  */
/* 	fabs(z_in(2,j)-z_in(n1,j)); */
  r_gap /= n2;

  s_gap = 0.0;
  for (i=1; i<=n1; i++)
    s_gap += fabs(x_in(i,1)-x_in(i,n2)) + 
      fabs(y_in(i,1)-y_in(i,n2)) + 
	fabs(z_in(i,1)-z_in(i,n2));
/*     s_gap += fabs(x_in(i,2)-x_in(i,n2)) +  */
/*       fabs(y_in(i,2)-y_in(i,n2)) +  */
/* 	fabs(z_in(i,2)-z_in(i,n2)); */
  s_gap /= n1;

  sgrid_ptr->r_period = (r_gap <= 1.e-7*l_scale);
  sgrid_ptr->s_period = (s_gap <= 1.e-7*l_scale);

/*   printf("The discrete surface `%s' seems to be %s periodic in r and %s " */
/* 	 "periodic in s.\n", sgrid_ptr->name,  */
/* 	 (sgrid_ptr->r_period? "" : "NON"), (sgrid_ptr->s_period? "" : "NON")); */
  
  sgrid_ptr->curve_ptr   = create_int_array_2d(2,2);
  sgrid_ptr->hole_curve_ = create_int_array_2d(2,2);

/* initialize curve and hole_curve */
  curve(1,1) = 0; hole_curve(1,1) = 0;
  curve(2,1) = 0; hole_curve(2,1) = 0;
  curve(1,2) = 0; hole_curve(1,2) = 0;
  curve(2,2) = 0; hole_curve(2,2) = 0;

/* initialize the surface label */
  sgrid_ptr->surface_label = present_surface;
  sgrid_ptr->vgrid_mother = mother;
  sgrid_ptr->side = 0;
  sgrid_ptr->dir = 0;

  sgrid_ptr->r_step = 0.0;
  sgrid_ptr->s_step = 0.0;

  sgrid_ptr->flag_ptr = create_int_array_2d(n1, n2);
  for (i=1; i<=n1; i++)
    for (j=1; j<=n2; j++)
      flag(i,j) = 1;      /* makes the inside_surface routine happy! */

  sgrid_ptr->x_ptr = x_;
  sgrid_ptr->y_ptr = y_;
  sgrid_ptr->z_ptr = z_;

  sgrid_ptr->quad_tree = NULL;

  /* set the default plot mode */
  sgrid_ptr->plot_mode = 1 | 4 | 8 | 256;
  sgrid_ptr->plot_it = TRUE;

  sgrid_ptr->r_step = 1.0/((real) (sgrid_ptr->r_dim - 1));
  sgrid_ptr->s_step = 1.0/((real) (sgrid_ptr->s_dim - 1));

/* ensure perfect periodicity if r_period = 1 or s_period = 1 */
  if (sgrid_ptr->r_period)
    for (j=1; j<= sgrid_ptr->s_dim; j++){
      x(sgrid_ptr->r_dim,j) = x(1,j);
      y(sgrid_ptr->r_dim,j) = y(1,j);
      z(sgrid_ptr->r_dim,j) = z(1,j);
    }

  if (sgrid_ptr->s_period)
    for (i=1; i<= sgrid_ptr->r_dim; i++){
      x(i,sgrid_ptr->s_dim) = x(i,1);
      y(i,sgrid_ptr->s_dim) = y(i,1);
      z(i,sgrid_ptr->s_dim) = z(i,1);
    }

  /* compute quad search tree */
  sgrid_ptr->quad_tree = new_subdomain_info(1, sgrid_ptr->r_dim, 1, 
					    sgrid_ptr->s_dim, sgrid_ptr );

  /* set default missmatch tolerance */
  sgrid_ptr->max_angle  = 10.0 * ((double)M_PI)/180.0;

  /* estimate how large max_n_dist needs to be. */
  sgrid_ptr->max_n_dist = 1.e-7; /* even a flat surface needs some tolerance
				    to make the search_quad_tree algorithm work */
  for (i = 2; i < sgrid_ptr->r_dim; i++)
    for (j = 2; j < sgrid_ptr->s_dim; j++)
      {
	compute_grid_normal( n_vec, i, j, sgrid_ptr );
	middle[0] = 0.25*(x(i+1,j) + x(i,j+1) + x(i-1,j) + x(i,j-1));
	middle[1] = 0.25*(y(i+1,j) + y(i,j+1) + y(i-1,j) + y(i,j-1));
	middle[2] = 0.25*(z(i+1,j) + z(i,j+1) + z(i-1,j) + z(i,j-1));
	i_dist = fabs(n_vec[0]*(x(i,j) - middle[0]) + n_vec[1]*(y(i,j) - middle[1]) + 
		      n_vec[2]*(z(i,j) - middle[2]) );
	sgrid_ptr->max_n_dist = real_max(sgrid_ptr->max_n_dist, i_dist);
      }
  sgrid_ptr->max_t_dist = real_max(1.e-8, sgrid_ptr->max_n_dist * 
				   tan( sgrid_ptr->max_angle ));

  return sgrid_ptr;
#undef x_in
#undef y_in
#undef z_in
}

discrete_surface *
delete_discrete_surface( discrete_surface *sgrid_ptr )
{
  if (sgrid_ptr == NULL) return NULL;

  delete_int_array_2d( sgrid_ptr->curve_ptr );
  delete_int_array_2d( sgrid_ptr->hole_curve_ );
  delete_int_array_2d( sgrid_ptr->flag_ptr );
  delete_real_array_2d( sgrid_ptr->x_ptr );
  delete_real_array_2d( sgrid_ptr->y_ptr );
  delete_real_array_2d( sgrid_ptr->z_ptr );

/* free the search tree */
  delete_subdomain_info( sgrid_ptr->quad_tree );

/* free the bounding box */
  free( sgrid_ptr->bb );

/* free the name */
  free( sgrid_ptr->name );

/* free the structure itself */
  free( sgrid_ptr );
  return NULL;
}


static subdomain_info *
new_subdomain_info(int i_min, int i_max, int j_min, int j_max, 
		   discrete_surface *sgrid_ptr )
{
  int i, j;
  subdomain_info *info;

  /* don't save any empty leafs */
  if (i_min == i_max || j_min == j_max) 
    return NULL;

  info = (subdomain_info *) malloc( sizeof(subdomain_info) );

  info->i_min = i_min;
  info->i_max = i_max;
  info->j_min = j_min;
  info->j_max = j_max;

  /* initialize bounding box */
  info->x_min =  1.e10;
  info->x_max = -1.e10;
  info->y_min =  1.e10;
  info->y_max = -1.e10;
  info->z_min =  1.e10;
  info->z_max = -1.e10;

  for (i=i_min; i<=i_max; i++)
    for (j=j_min; j<=j_max; j++)
      {
	info->x_min = real_min( info->x_min, x(i,j) );
	info->x_max = real_max( info->x_max, x(i,j) );
	info->y_min = real_min( info->y_min, y(i,j) );
	info->y_max = real_max( info->y_max, y(i,j) );
	info->z_min = real_min( info->z_min, z(i,j) );
	info->z_max = real_max( info->z_max, z(i,j) );
      }

  /* are we at a leaf ? */
  if (i_max-i_min == 1 && j_max-j_min == 1)
    {
      info->low_left  = NULL;
      info->low_right = NULL;
      info->upp_left  = NULL;
      info->upp_right = NULL;
    }
  else
    {
      /* recurse */
      info->low_left  = new_subdomain_info(i_min, (i_max+i_min)/2, 
					   j_min, (j_max+j_min)/2, sgrid_ptr );
      info->low_right = new_subdomain_info((i_max+i_min)/2, i_max, 
					   j_min, (j_max+j_min)/2, sgrid_ptr );
      info->upp_left  = new_subdomain_info(i_min, (i_max+i_min)/2, 
					   (j_max+j_min)/2, j_max, sgrid_ptr );
      info->upp_right = new_subdomain_info((i_max+i_min)/2, i_max, 
					   (j_max+j_min)/2, j_max, sgrid_ptr );
    }
  /* done */
  return info;
}

intersection_point *
new_intersection_point(real xp, real yp, real zp, real sin_theta,
		       discrete_surface *sgrid_ptr)
{
  intersection_point *ip_ptr;
  ip_ptr = (intersection_point *) malloc( sizeof(intersection_point) );
  ip_ptr->ic[0] = xp;
  ip_ptr->ic[1] = yp;
  ip_ptr->ic[2] = zp;
  ip_ptr->sin_theta = sin_theta;
  ip_ptr->sgrid = sgrid_ptr;
  ip_ptr->unique = TRUE;
  ip_ptr->prev_ip = NULL;

  return ip_ptr;
}

stack_of_ip *
new_stack_of_ip(void)
{
  stack_of_ip *stack;

  stack = (stack_of_ip *) malloc( sizeof(stack_of_ip) );
  stack->n_ip = 0;
  stack->last_ip = NULL;

  return stack;
}

stack_of_ip *
delete_stack_of_ip(stack_of_ip *stack)
{
  intersection_point *ip_ptr, *next_victim;

  ip_ptr = stack->last_ip;
  while (ip_ptr != NULL)
    {      
      next_victim = ip_ptr->prev_ip;
      free(ip_ptr);
      ip_ptr= next_victim;
    }

  free(stack);
  return NULL;
}

void
push_ip(intersection_point *ip_ptr, stack_of_ip *stack)
{
  ip_ptr->prev_ip = stack->last_ip;
  stack->last_ip = ip_ptr;
  stack->n_ip++;
}

int
compute_grid_normal(real n_vec[3], int i, int j, discrete_surface *sgrid_ptr)
{
  double xr, yr, zr, xs, ys, zs, length, sign;
  const real eps = 1.e-10;
  int i1, i0, j1, j0;

  /*  sign = (sgrid_ptr->orientation)? -1.0 : 1.0;*/
  sign = 1.0;
  /* r-differences */
  if (sgrid_ptr->r_period)
    {
      if (i == 1 || i == sgrid_ptr->r_dim)
	{
	  i1 = 1+1; i0 = sgrid_ptr->r_dim-1;
	}
      else
	{
	  i1 = i+1; i0 = i-1;
	}
    }
  else
    {
      if (i == 1)
	{
	  i1 = i+1; i0 = i;
	}
      else if (i == sgrid_ptr->r_dim)
	{
	  i1 = i; i0 = i-1;
	}
      else
	{
	  i1 = i+1; i0 = i-1;
	}
    }

  xr = x(i1,j) - x(i0,j);
  yr = y(i1,j) - y(i0,j);
  zr = z(i1,j) - z(i0,j);

  /* s-differences */
  if (sgrid_ptr->s_period)
    {
      if (j == 1 || j == sgrid_ptr->s_dim)
	{
	  j1 = 1+1; j0 = sgrid_ptr->s_dim-1;
	}
      else
	{
	  j1 = j+1; j0 = j-1;
	}
    }
  else
    {
      if (j == 1)
	{
	  j1 = j+1; j0 = j;
	}
      else if (j == sgrid_ptr->s_dim)
	{
	  j1 = j; j0 = j-1;
	}
      else
	{
	  j1 = j+1; j0 = j-1;
	}
    }

  xs = x(i,j1) - x(i,j0);
  ys = y(i,j1) - y(i,j0);
  zs = z(i,j1) - z(i,j0);

  /* normals */
  n_vec[0] =  sign * (yr * zs - zr * ys);
  n_vec[1] =  sign * (zr * xs - xr * zs);
  n_vec[2] =  sign * (xr * ys - yr * xs);

  /* normalize */
  length = sqrt( n_vec[0] * n_vec[0] + n_vec[1] * n_vec[1] + n_vec[2] * n_vec[2] );
  if (length > eps)
    {
      for (i=0; i<3; i++) n_vec[i] = n_vec[i]/length;
      return OK;
    }
  else
    return ERROR;
}


int
compute_cell_normal(real n_vec[3], int ic, int jc, discrete_surface *sgrid_ptr)
{
  int i;
  double xr, yr, zr, xs, ys, zs, length, sign;
  const real eps = 1.e-10;

  sign = 1.0;

  if (ic < 1)
    ic = 1;
  else if (ic > sgrid_ptr->r_dim-1)
    ic = sgrid_ptr->r_dim-1;

  if (jc < 1)
    jc = 1;
  else if (jc > sgrid_ptr->s_dim-1)
    jc = sgrid_ptr->s_dim-1;

  /* r-differences */
  xr = x(ic+1,jc) - x(ic,jc);
  yr = y(ic+1,jc) - y(ic,jc);
  zr = z(ic+1,jc) - z(ic,jc);

  /* s-differences */
  xs = x(ic,jc+1) - x(ic,jc);
  ys = y(ic,jc+1) - y(ic,jc);
  zs = z(ic,jc+1) - z(ic,jc);

  /* normals */
  n_vec[0] =  sign * (yr * zs - zr * ys);
  n_vec[1] =  sign * (zr * xs - xr * zs);
  n_vec[2] =  sign * (xr * ys - yr * xs);

  /* normalize */
  length = sqrt( n_vec[0] * n_vec[0] + n_vec[1] * n_vec[1] + n_vec[2] * n_vec[2] );
  if (length > eps)
    {
      for (i=0; i<3; i++) n_vec[i] = n_vec[i]/length;
      return OK;
    }
  else
    return ERROR;
}

void
compute_continuous_normal(real normal[3], real normal_r[3], real normal_s[3],
			  real r, real s, discrete_surface *sgrid_ptr)
{
  real n_vec0[3], n_vec1[3], n_vec2[3], n_vec3[3];
  int i_cell, j_cell, q;
  double length;
  const real eps = 1.e-10;

  /* periodicity */
  if (sgrid_ptr->r_period && r < 0.0) r += 1.0;
  else if (sgrid_ptr->r_period && r > 1.0) r -= 1.0;

  if (sgrid_ptr->s_period && s < 0.0) s += 1.0;
  else if (sgrid_ptr->s_period && s > 1.0) s -= 1.0;

  /* find cell */
  i_cell = 1 + r/sgrid_ptr->r_step;
  j_cell = 1 + s/sgrid_ptr->s_step;
  /* avoid out of bounds */
  i_cell = int_max( 1, int_min( sgrid_ptr->r_dim-1, i_cell ) );
  j_cell = int_max( 1, int_min( sgrid_ptr->s_dim-1, j_cell ) );

  /* normals at the four corners of the cell */
  compute_grid_normal( n_vec0, i_cell  , j_cell  , sgrid_ptr );
  compute_grid_normal( n_vec1, i_cell+1, j_cell  , sgrid_ptr );
  compute_grid_normal( n_vec2, i_cell  , j_cell+1, sgrid_ptr );
  compute_grid_normal( n_vec3, i_cell+1, j_cell+1, sgrid_ptr );

  /* define linear base functions */
#define R0 ((i_cell-1)   * sgrid_ptr->r_step)
#define R1 ((i_cell+1-1) * sgrid_ptr->r_step)
#define ALPHA_0(r)   ((R1 - r)/sgrid_ptr->r_step)
#define ALPHA_1(r)   ((r - R0)/sgrid_ptr->r_step)
#define ALPHA_0r(r)   (-1.0/sgrid_ptr->r_step)
#define ALPHA_1r(r)   (1.0/sgrid_ptr->r_step)

#define S0 ((j_cell-1)  * sgrid_ptr->s_step)
#define S1 ((j_cell+1-1)* sgrid_ptr->s_step)
#define BETA_0(s)   ((S1 - s)/sgrid_ptr->s_step)
#define BETA_1(s)   ((s - S0)/sgrid_ptr->s_step)
#define BETA_0s(s)   (-1.0/sgrid_ptr->s_step)
#define BETA_1s(s)   (1.0/sgrid_ptr->s_step)

  /* weight the normals together */
  for (q=0; q<3; q++)
    {
      normal[q] = 
	ALPHA_0(r) * BETA_0(s) * n_vec0[q] +
	  ALPHA_1(r) * BETA_0(s) * n_vec1[q] + 
	    ALPHA_0(r) * BETA_1(s) * n_vec2[q] + 
	      ALPHA_1(r) * BETA_1(s) * n_vec3[q];
      normal_r[q] = 
	ALPHA_0r(r) * BETA_0(s) * n_vec0[q] +
	  ALPHA_1r(r) * BETA_0(s) * n_vec1[q] + 
	    ALPHA_0r(r) * BETA_1(s) * n_vec2[q] + 
	      ALPHA_1r(r) * BETA_1(s) * n_vec3[q];
      normal_s[q] = 
	ALPHA_0(r) * BETA_0s(s) * n_vec0[q] +
	  ALPHA_1(r) * BETA_0s(s) * n_vec1[q] + 
	    ALPHA_0(r) * BETA_1s(s) * n_vec2[q] + 
	      ALPHA_1(r) * BETA_1s(s) * n_vec3[q];
    }
#undef R0
#undef R1
#undef ALPHA_0
#undef ALPHA_1
#undef ALPHA_0r
#undef ALPHA_1r

#undef S0
#undef S1
#undef BETA_0
#undef BETA_1
#undef BETA_0s
#undef BETA_1s

  /* normalize */
  length = sqrt( normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2] );
  if (length > eps)
    {
      for (q=0; q<3; q++)
	{
	  normal[q] = normal[q]/length;
	  normal_r[q] = normal_r[q]/length;
	  normal_s[q] = normal_s[q]/length;
	}
    }
}


static subdomain_info *
delete_subdomain_info( subdomain_info *info )
{
  if (info == NULL)
    return NULL;

/* delete subtrees */
  delete_subdomain_info( info->low_left );
  delete_subdomain_info( info->low_right);
  delete_subdomain_info( info->upp_left );
  delete_subdomain_info( info->upp_right);

/* delete the structure itself */
  free(info);

  return NULL;
}

hole_surface_link *
get_hole_surface( input_output *io_ptr, hole_surface_list *head)
{
  hole_surface *overg;
  linked_list_member *this_link;
  char **command;
  int i, icom, last_com, *save_on_copy=NULL, *argument=NULL, level=0, quit;

  if (head->n_members == 0)
    {
      printf("The list of hole-surfaces is empty.\n");
      return NULL;
    }

  last_com = head->n_members + 1;
  command = (char **) malloc( (last_com+1)*sizeof(char*) );

  i=0;
  for (this_link = head->first; this_link != NULL; 
       this_link = this_link->next)
    {
      overg = this_link->data;
      command[i] = overg->name;
      i++;
    }
  command[last_com-1] = "help";
  command[last_com]   = "cancel";

  quit = 0;
  do{
      icom = get_command(io_ptr, "hole-surface name: ", command, 
			 last_com, level, save_on_copy, argument);

      /* default is to quit */
      quit = 1;

      if (icom == -1)
	/* not a valid command */
	quit = 0;
      else if (icom == last_com-1)
	{
	  /* help */
	  printf(
		 " o  Enter one of the hole surface names if you wish to proceed.\n"
		 " o  Enter cancel if you do not wish to proceed. You will then return\n"
		 "    to the previous command level.\n");
	  quit = 0;
	}

    }
  while (!quit);

  /* return the pointer to the selected curve */
  i=0;
  for (this_link = head->first; this_link != NULL; 
       this_link = this_link->next)
    {
      if (i==icom) break;
      i++;
    }

  /* free the pointers to the strings */
  free(command);

  return this_link;
}

discrete_surface_link *
get_discrete_surface( input_output *io_ptr, discrete_surface_list *head)
{
  discrete_surface *sgrid_ptr;
  linked_list_member *this_link;
  char **command;
  int i, icom, last_com, *save_on_copy=NULL, *argument=NULL, level=0, quit;

  if (head->n_members == 0)
    {
      printf("The list of discrete surfaces is empty.\n");
      return NULL;
    }

  last_com = head->n_members + 1;
  command = (char **) malloc( (last_com+1)*sizeof(char*) );

  i=0;
/* add color information to the name! */
  for (this_link = head->first; this_link != NULL; 
       this_link = this_link->next)
    {
      sgrid_ptr = this_link->data;
      command[i] = (char *) malloc( (strlen(sgrid_ptr->name) + 
				     strlen(ogl_color_name(sgrid_ptr->priority)) + 3) 
				   * sizeof(char) );
      sprintf( command[i], "%s(%s)", sgrid_ptr->name, ogl_color_name(sgrid_ptr->priority) );

      i++;
    }
  command[last_com-1] = "help";
  command[last_com]   = "cancel";

  quit = 0;
  do{
      icom = get_command(io_ptr, "discrete surface name: ", command, 
			 last_com, level, save_on_copy, argument);

      /* default is to quit */
      quit = 1;

      if (icom == -1)
	/* not a valid command */
	quit = 0;
      else if (icom == last_com-1)
	{
	  /* help */
	  printf(
		 " o  Enter one of the discrete surface names if you wish to proceed.\n"
		 " o  Enter cancel if you do not wish to proceed. You will then return\n"
		 "    to the previous command level.\n");
	  quit = 0;
	}

    }
  while (!quit);

  /* return the pointer to the selected curve */
  i=0;
  for (this_link = head->first; this_link != NULL; 
       this_link = this_link->next)
    {
      if (i==icom) break;
      i++;
    }

/* free the strings */
  for (i=0; i<= last_com-2; i++)
    free(command[i]);

  /* free the pointers to the strings */
  free(command);

  return this_link;
}

static void
update_hs_bb( hole_surface *overg )
{
  linked_list_member *sgrid_link;
  discrete_surface *sgrid_ptr;
  if (overg->grid_list == NULL || overg->grid_list->n_members == 0)
    {
      overg->bb->x_min = 0.0;
      overg->bb->x_max = 1.0;
      overg->bb->y_min = 0.0;
      overg->bb->y_max = 1.0;
      overg->bb->z_min = 0.0;
      overg->bb->z_max = 1.0;
    }
  else
    {
      overg->bb->x_min = 1.0e10;
      overg->bb->x_max =-1.0e10;
      overg->bb->y_min = 1.0e10;
      overg->bb->y_max =-1.0e10;
      overg->bb->z_min = 1.0e10;
      overg->bb->z_max =-1.0e10;
      for (sgrid_link = overg->grid_list->first; sgrid_link != NULL; 
	   sgrid_link = sgrid_link->next)
	{
	  sgrid_ptr = sgrid_link->data;
	  overg->bb->x_min = real_min(overg->bb->x_min, sgrid_ptr->bb->x_min);
	  overg->bb->x_max = real_max(overg->bb->x_max, sgrid_ptr->bb->x_max);
	  overg->bb->y_min = real_min(overg->bb->y_min, sgrid_ptr->bb->y_min);
	  overg->bb->y_max = real_max(overg->bb->y_max, sgrid_ptr->bb->y_max);
	  overg->bb->z_min = real_min(overg->bb->z_min, sgrid_ptr->bb->z_min);
	  overg->bb->z_max = real_max(overg->bb->z_max, sgrid_ptr->bb->z_max);
	}
    }
}

