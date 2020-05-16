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
#define sqr(x) ((x)*(x))
/* static prototypes */
linked_list_member *
get_vgrid_ptr( input_output *io_ptr, linked_list *head);
/* end static prototypes */


overlapping_3d_grid *
new_overlapping_3d_grid(char *name, linked_list *volume_mappings){
  overlapping_3d_grid *over3d;
  component_vgrid *vgrid;
  linked_list_member *this_link;
  int i, j;

/* don't call this rutine with an empty `volume_mappings' list */

  over3d = (overlapping_3d_grid *) malloc(sizeof(overlapping_3d_grid));
  
/* initially invalid overlap */
  over3d->valid = 0;

/* pointers to the linked list of component grids */
  over3d->grid_list = new_linked_list();
  over3d->mapping_list = new_linked_list();

/* pointer to the linked list of overlapping surface grids */
  over3d->surface_grid_list = new_linked_list();

/* make new links and copy the pointer to the data from the links in the */
/* volume_mappings list */
  for (this_link = volume_mappings->last; this_link != NULL;
       this_link = this_link->prev){
    vgrid = this_link->data;
    new_link( over3d->mapping_list )->data = vgrid;
  }

/* number of ghost points */
  over3d->extra = 0;

/* discretization width quantities */
  over3d->disc_width    = 3;
  over3d->normal_width  = 2;
  over3d->tangent_width = 3;
  over3d->corner_width  = 2;
    
/* interpolation width */
  over3d->interp_width  = 3;
    
/* interpolation type */
  over3d->interp_type   = 'i';
    
/* boundary surface missmatch tolerance */
  over3d->max_b_dist    = 1.e-7;

/* number of extra points at periodic boundaries */
  over3d->extra_period = 0;

/* default trimming style is to minimize the overlap */
  over3d->trim_style = 0;

/* default is to tell the user some things during the overlap algorithm */
  over3d->verbose = 2;

/* allocate space for the bounding box */
  over3d->bb = (bounding_box *) malloc( sizeof(bounding_box) );

/* bounding box */
  volume_list_bb( over3d->bb, over3d->mapping_list );

/* copy the name */
  over3d->name = (char *) malloc( (strlen(name)+1)*sizeof(char) );
  strcpy( over3d->name, name );

/* plot grid boundaries and the coordinate cage */
  over3d->plot_mode = 2 | 256;

/* The default is to not show the physical boundary representation or the */
/* volume holes during the overlap algorithm */
  over3d->show_boundary = FALSE;
  over3d->show_volume_holes = FALSE;

/* The default is to NOT correct mismatch at physical boundaries */
  over3d->correct_mismatch = FALSE;

/* one clip plane */
  over3d->n_clip_planes = 1;

/* define default clip planes */
  for (i=0; i<3; i++){
    for (j=0; j<3; j++)
      over3d->clip_plane[i][j] = 0.0;
/* the homogeneous coordinate is always 1, so that the first 3 components are the */
/* cartesian coordinates */
    over3d->clip_plane[i][3] = 1.0;
  }
/* take the normal of the clip plane to be x=1, y=0, z=0 */
  over3d->clip_plane[0][0] = 1.0;

/* default sphere size */
  over3d->sphere_size = 1.e-2;
  return over3d;
}

overlapping_3d_grid *
delete_over3d_grid( overlapping_3d_grid *over3d ){
  linked_list_member *vgrid_link, *overg_link;

  if (over3d == NULL) return NULL;

/* delete all component_vgrid's */
  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL;
       vgrid_link = vgrid_link->next){
    delete_vgrid( (component_vgrid *) vgrid_link->data );
  }
  delete_linked_list( over3d->grid_list );

/* delete the linked list with overlapping surface grids */
  for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;
       overg_link = overg_link->next){
    delete_hole_surface( (hole_surface *) overg_link->data );
  }
  delete_linked_list( over3d->surface_grid_list );

/* free the linked list of volume mappings */
  delete_linked_list( over3d->mapping_list );

/* free the bounding box and the name */
  free( over3d->bb );
  free( over3d->name );

  free( over3d );
  return NULL;
}

interp_3d_point *
new_3d_interp(int i_loc, int j_loc, int k_loc, component_vgrid *vgrid_loc, 
	      real r_loc, real s_loc, real t_loc){
  interp_3d_point *interp3;
  
  interp3 = (interp_3d_point *) malloc( sizeof(interp_3d_point) );
  
  interp3->prev = NULL;
  interp3->active = 1;
  interp3->corrected[0] = FALSE;
  interp3->corrected[1] = FALSE;
  interp3->corrected[2] = FALSE;

/* fill in the interpolation point info later */
  interp3->i_point = -999;
  interp3->j_point = -999;
  interp3->k_point = -999;

/* copy the location info from the arguments */
  interp3->vgrid_loc = vgrid_loc; 

  interp3->r_loc = r_loc;
  interp3->s_loc = s_loc;
  interp3->t_loc = t_loc;

  interp3->i_loc = i_loc;
  interp3->j_loc = j_loc;
  interp3->k_loc = k_loc;
  
  return interp3;
}

void 
remove_deactive_3d_interp(component_vgrid *vgrid){
  interp_3d_point *interp3, *prev_victim, *next_victim;

/* are there any interpolation points? */
  if (vgrid->last_interp != NULL){
    interp3 = vgrid->last_interp;
    prev_victim = NULL;
    do{
      next_victim = interp3->prev;
      if (!interp3->active){
/* remove this interpolation point */
	if (prev_victim == NULL)
	  vgrid->last_interp = interp3->prev;
	else
	  prev_victim->prev = interp3->prev;
	free(interp3);
      }
      else
	prev_victim = interp3;

      interp3 = next_victim;
    }
    while (interp3 != NULL);
  } 
}

interp_3d_point *
delete_3d_interpolation_list( interp_3d_point *last_interp ){
  interp_3d_point *interp3, *next_victim;

  interp3 = last_interp;
  while (interp3 != NULL){
    next_victim = interp3->prev;
    free(interp3);
    interp3= next_victim;
  }

  return NULL;
}

void
insert_3d_interp_point( interp_3d_point *new_interp, component_vgrid *vgrid ){
  int i;
  if (new_interp == NULL) return;

/* update the pointers */
  new_interp->prev = vgrid->last_interp;
  vgrid->last_interp = new_interp;
/* increase the count */
  vgrid->n_interp++;
}

void 
new_3d_bad_point(int i, int j, int k, int type, component_vgrid *vgrid){
  bad_3d_point *new_bad_point;

  new_bad_point = (bad_3d_point *) malloc( sizeof(bad_3d_point) );
  new_bad_point->i = i;
  new_bad_point->j = j;
  new_bad_point->k = k;
  new_bad_point->type = type;

  if (vgrid->last_bad_point == NULL)
    new_bad_point->prev = NULL;
  else
    new_bad_point->prev = vgrid->last_bad_point;
    
  vgrid->last_bad_point = new_bad_point;
}

bad_3d_point *
delete_bad_3d_list( bad_3d_point *last_bad_point ){
  bad_3d_point *bad_point_ptr, *next_victim;

  bad_point_ptr = last_bad_point;
  while (bad_point_ptr != NULL){
    next_victim = bad_point_ptr->prev;
    free(bad_point_ptr);
    bad_point_ptr = next_victim;
  }
  
  return NULL;
}

component_vgrid *
new_vgrid( volume_mapping *volume ){
  component_vgrid *vgrid;
  int i, j;

  if (!volume)
    {
      printf("ERROR: new_vgrid called with volume==NULL\n");
      return NULL;
    }

  vgrid = (component_vgrid *) malloc(sizeof(component_vgrid));
  
  vgrid->no_hole = volume->no_hole;

  vgrid->name = (char *) malloc( (strlen(volume->name)+1)*sizeof(char) );
  vgrid->name = strcpy( vgrid->name, volume->name );

/* initialize all fields */
  if (volume->type)
/* vgrid->grid_type is set to 1 for Cartesian grids. */
    vgrid->grid_type = !strcmp(volume->type, "Cartesian grid"); 
  else
    vgrid->grid_type = 0;

  vgrid->priority = 0;
  vgrid->orientation = 0;

  vgrid->r1_dim = 0;
  vgrid->r2_dim = 0;
  vgrid->r3_dim = 0;

  vgrid->r1_period = volume->r1_period;
  vgrid->r2_period = volume->r2_period;
  vgrid->r3_period = volume->r3_period;

  vgrid->range_ptr = create_int_array_2d(2,3);

  vgrid->bc_ptr      = create_int_array_2d(2,3);
  vgrid->surf_ptr    = create_int_array_2d(2,3);
  vgrid->edge_curve_ = create_int_array_1d(12);

  vgrid->r1_step = 0.0;
  vgrid->r2_step = 0.0;
  vgrid->r3_step = 0.0;

  vgrid->flag_ptr = NULL;
  vgrid->x_ptr = NULL;
  vgrid->y_ptr = NULL;
  vgrid->z_ptr = NULL;

  vgrid->n_interp = 0;
  vgrid->last_interp = NULL;

  vgrid->last_bad_point = NULL;

  vgrid->oct_tree = NULL;

  vgrid->verbose = 0;

/* missmatch tolerances */
  vgrid->max_s_dist = 0.0;
  vgrid->max_e_dist = 0.0;

  vgrid->bb = (bounding_box *) malloc( sizeof(bounding_box) );

  vgrid->bb->x_min = volume->bb->x_min;
  vgrid->bb->x_max = volume->bb->x_max;
  vgrid->bb->y_min = volume->bb->y_min;
  vgrid->bb->y_max = volume->bb->y_max;
  vgrid->bb->z_min = volume->bb->z_min;
  vgrid->bb->z_max = volume->bb->z_max;

/* compute lengthscale */
  vgrid->length_scale = sqrt( sqr(volume->bb->x_max-volume->bb->x_min) +
			      sqr(volume->bb->y_max-volume->bb->y_min) +
			      sqr(volume->bb->z_max-volume->bb->z_min) );

/* set the default plot mode */
  vgrid->plot_mode = 1 | 4 | 8 | 256;

/* set the color */
  vgrid->color = volume->color;

/* by default, plot this component */
  vgrid->plot_it = 1;

  vgrid->bounding_surface = NULL;

  vgrid->analytical_inverse    = volume->inverse_known;
  vgrid->forward_mapping = volume->forward_mapping;
  vgrid->inverse_mapping = volume->inverse_mapping;
  vgrid->vgrid_data_ptr        = volume->volume_data_ptr;

/* copy translation / rotation info */
  vgrid->rotation_ = create_real_array_2d(3,3);
  vgrid->transformed = volume->transformed;

/* copy the rotation matrix */
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++){
      rot_vg(i,j) = rot_vol(i,j);
    }

/* copy the translation vectors */
  for (i=0; i<3; i++)
    {
      vgrid->translation[i] = volume->translation[i];
      vgrid->pre_trans[i] = volume->pre_trans[i];
    }

/* initialize the 2-d arrays */
  for (i=1; i<=2; i++)
    for (j=1; j<=3; j++){
      bc3(i,j)    = bc_vol(i,j);
      range3(i,j) = 0;
      surf3(i,j)  = surf_vol(i,j);
    }

/* initialize the edge curve info */
  for (i=1; i<=12; i++)
    edge_curve(i) = edge_vol(i);

  return vgrid;
}


component_vgrid *
delete_vgrid( component_vgrid *vgrid ){
  if (vgrid == NULL) return NULL;

  delete_int_array_2d( vgrid->range_ptr );
  delete_int_array_2d( vgrid->bc_ptr );
  delete_int_array_2d( vgrid->surf_ptr );
  delete_int_array_1d( vgrid->edge_curve_ );
  delete_int_array_3d( vgrid->flag_ptr );
  delete_real_array_3d( vgrid->x_ptr );
  delete_real_array_3d( vgrid->y_ptr );
  delete_real_array_3d( vgrid->z_ptr );
  delete_real_array_2d( vgrid->rotation_ );

/* free the bounding surface */
  vgrid->bounding_surface = delete_hole_surface( vgrid->bounding_surface );

/* free the list of interpolation points */
  vgrid->last_interp = delete_3d_interpolation_list( vgrid->last_interp ); 

/* free the list of bad points */
  vgrid->last_bad_point = delete_bad_3d_list( vgrid->last_bad_point ); 

/* free the search tree */
  delete_oct_tree( vgrid->oct_tree ); 

/* free the bounding box */
  free( vgrid->bb );

/* free the name */
  free( vgrid->name );

/* free the structure itself */
  free( vgrid );
  return NULL;
}

int 
forward_vgrid_mapping( grid_point *gp_ptr, component_vgrid *vgrid ){
  grid_point gp;
  int msg;
  real x,y,z;

  if (vgrid->transformed){
    gp.r = gp_ptr->r;
    gp.s = gp_ptr->s;
    gp.t = gp_ptr->t;

/* take care of any errors in the specific mapping routine and pass them on */
    msg = (*vgrid->forward_mapping)( &gp, vgrid->vgrid_data_ptr );

/* pre-translation */
    x = gp.x + vgrid->pre_trans[0];
    y = gp.y + vgrid->pre_trans[1];
    z = gp.z + vgrid->pre_trans[2];

/* rotate and translate coordinates */
    gp_ptr->x = vgrid->translation[0] + x * rot_vg(1,1) + y * rot_vg(1,2) + 
      z * rot_vg(1,3);
    gp_ptr->y = vgrid->translation[1] + x * rot_vg(2,1) + y * rot_vg(2,2) + 
      z * rot_vg(2,3);
    gp_ptr->z = vgrid->translation[2] + x * rot_vg(3,1) + y * rot_vg(3,2) + 
      z * rot_vg(3,3);

/* rotate the Jacobian */
    gp_ptr->xr = gp.xr * rot_vg(1,1) + gp.yr * rot_vg(1,2) + gp.zr * rot_vg(1,3);
    gp_ptr->yr = gp.xr * rot_vg(2,1) + gp.yr * rot_vg(2,2) + gp.zr * rot_vg(2,3);
    gp_ptr->zr = gp.xr * rot_vg(3,1) + gp.yr * rot_vg(3,2) + gp.zr * rot_vg(3,3);

    gp_ptr->xs = gp.xs * rot_vg(1,1) + gp.ys * rot_vg(1,2) + gp.zs * rot_vg(1,3);
    gp_ptr->ys = gp.xs * rot_vg(2,1) + gp.ys * rot_vg(2,2) + gp.zs * rot_vg(2,3);
    gp_ptr->zs = gp.xs * rot_vg(3,1) + gp.ys * rot_vg(3,2) + gp.zs * rot_vg(3,3);

    gp_ptr->xt = gp.xt * rot_vg(1,1) + gp.yt * rot_vg(1,2) + gp.zt * rot_vg(1,3);
    gp_ptr->yt = gp.xt * rot_vg(2,1) + gp.yt * rot_vg(2,2) + gp.zt * rot_vg(2,3);
    gp_ptr->zt = gp.xt * rot_vg(3,1) + gp.yt * rot_vg(3,2) + gp.zt * rot_vg(3,3);
  }
  else{
/* take care of any errors in the specific mapping routine and pass them on */
    msg = (*vgrid->forward_mapping)( gp_ptr, vgrid->vgrid_data_ptr );
  }

  return msg;
}

int 
inverse_vgrid_mapping( grid_point *gp_ptr, component_vgrid *vgrid ){
  real x, y, z;

/* take care of any errors in the specific mapping routine and pass them on */
  if (vgrid == NULL || vgrid->inverse_mapping == NULL){
    printf("ERROR: inverse_vgrid_mapping was called with inconsistent pointers.\n");
    return -1;
  }

  if (vgrid->transformed){
/* invert the rotation and / or translation */
    x = gp_ptr->x - vgrid->translation[0];
    y = gp_ptr->y - vgrid->translation[1];
    z = gp_ptr->z - vgrid->translation[2];

    gp_ptr->x = x*rot_vg(1,1) + y*rot_vg(2,1) + z*rot_vg(3,1) - vgrid->pre_trans[0];
    gp_ptr->y = x*rot_vg(1,2) + y*rot_vg(2,2) + z*rot_vg(3,2) - vgrid->pre_trans[1];
    gp_ptr->z = x*rot_vg(1,3) + y*rot_vg(2,3) + z*rot_vg(3,3) - vgrid->pre_trans[2];
  }

/* pass any errors from the specific mapping routine to the calling function */
  return (*vgrid->inverse_mapping)( gp_ptr, vgrid->vgrid_data_ptr );
}

linked_list_member *
get_vgrid_ptr( input_output *io_ptr, linked_list *head){
  component_vgrid *vgrid;
  linked_list_member *this_link;
  char **command;
  int i, icom, ncom, *save_on_copy=NULL, *argument=NULL, level=0, quit;

  if (head->n_members == 0){
    printf("The volume grid list is empty.\n");
    return NULL;
  }

  ncom = head->n_members + 1;
  command = (char **) malloc( (ncom+1)*sizeof(char*) );

  i=0;
/* add color information to the name! */
  for (this_link = head->first; this_link != NULL; 
       this_link = this_link->next){
    vgrid = this_link->data;
    command[i] = (char *) malloc( (strlen(vgrid->name) + 
				   strlen(ogl_color_name(vgrid->color)) + 3) 
				 * sizeof(char) );
    sprintf( command[i], "%s(%s)", vgrid->name, ogl_color_name(vgrid->color) );
    i++;
  }
  command[ncom-1] = "help";
  command[ncom]   = "cancel";

  quit = 0;
  do{
    icom = get_command(io_ptr, "volume grid name: ", command, ncom, level, 
		       save_on_copy, argument);

/* default is to quit */
    quit = 1;

    if (icom == -1)
/* not a valid command */
      quit = 0;
    else if (icom == ncom-1){
/* help */
      printf(
" o  Enter one of the surface names if you wish to proceed.\n"
" o  Enter cancel if you do not wish to proceed. You will then return to the\n"
"    previous command level.\n");
      quit = 0;
    }

  }
  while (!quit);

/* return the pointer to the selected curve */
  i=0;
  for (this_link = head->first; this_link != NULL; 
       this_link = this_link->next){
    if (i==icom) break;
    i++;
  }

/* free the strings */
  for (i=0; i<= ncom-2; i++)
    free(command[i]);

/* free the pointers to the strings */
  free(command);

  return this_link;

}

linked_list_member *
get_over3d_ptr( input_output *io_ptr, linked_list *head){
  overlapping_3d_grid *over3d;
  linked_list_member *over3d_link;
  char **command;
  int i, icom, last_com, *save_on_copy=NULL, *argument=NULL, level=0, quit;

  if (head->n_members == 0){
    printf("The list of overlapping grids is empty.\n");
    return NULL;
  }

  last_com = head->n_members + 1;
  command = (char **) malloc( (last_com+1)*sizeof(char*) );

  i=0;
  for (over3d_link = head->first; over3d_link != NULL; 
       over3d_link = over3d_link->next){
    over3d = over3d_link->data;
    command[i] = over3d->name;
    i++;
  }
  command[last_com-1] = "help";
  command[last_com]   = "cancel";

  quit = 0;
  do{
    icom = get_command(io_ptr, "overlapping grid name: ", command, 
		       last_com, level, save_on_copy, argument);

/* default is to quit */
    quit = 1;

    if (icom == -1)
/* not a valid command */
      quit = 0;
    else if (icom == last_com-1){
/* help */
      printf(
" o  Enter one of the overlapping grid names if you wish to proceed.\n"
" o  Enter cancel if you do not wish to proceed. You will then return to the\n"
"    previous command level.\n");
      quit = 0;
    }

  }
  while (!quit);

/* return the pointer to the selected overlapping volume grid */
  i=0;
  for (over3d_link = head->first; over3d_link != NULL; 
       over3d_link = over3d_link->next){
    if (i==icom) break;
    i++;
  }

/* free the pointers to the strings */
  free(command);

  return over3d_link;
}

void
copy_from_volume_mapping( component_vgrid *vgrid, volume_mapping *volume ){
  int i, j, k, status=OK;
  grid_point gp;

/* delete any existing arrays */
  vgrid->x_ptr = delete_real_array_3d(vgrid->x_ptr);
  vgrid->y_ptr = delete_real_array_3d(vgrid->y_ptr);
  vgrid->z_ptr = delete_real_array_3d(vgrid->z_ptr);
  vgrid->flag_ptr = delete_int_array_3d(vgrid->flag_ptr);

/* define the multi-dimensional arrays */
  vgrid->x_ptr = create_real_array_3d(vgrid->r1_dim,vgrid->r2_dim,vgrid->r3_dim);
  vgrid->y_ptr = create_real_array_3d(vgrid->r1_dim,vgrid->r2_dim,vgrid->r3_dim);
  vgrid->z_ptr = create_real_array_3d(vgrid->r1_dim,vgrid->r2_dim,vgrid->r3_dim);
  vgrid->flag_ptr = create_int_array_3d(vgrid->r1_dim,vgrid->r2_dim,vgrid->r3_dim);

/* assign x, y, xr, yr, xs and ys. */
  for (i = 1; i<= vgrid->r1_dim; i++){
    for (j = 1; j<= vgrid->r2_dim; j++){
      for (k = 1; k<= vgrid->r3_dim; k++){
	gp.r = (i-range3(1,1)) * vgrid->r1_step;
	gp.s = (j-range3(1,2)) * vgrid->r2_step;
	gp.t = (k-range3(1,3)) * vgrid->r3_step;
	status = forward_volume_mapping(&gp, volume); 
	x3(i,j,k) = gp.x;
	y3(i,j,k) = gp.y;
	z3(i,j,k) = gp.z;
      }
    }
  }
/* error check */
  if (status != OK)
    printf("Warning: copy_from_volume_mapping: error in grid `%s'\n",vgrid->name);

}

void
overlap_list_bb(bounding_box *bb, overlapping_3d_list *over3d_list)
{
  overlapping_3d_link *over3d_link;
  overlapping_3d_grid *over3d;

  if (!over3d_list  || over3d_list->n_members == 0)
    {
      bb->x_min = 0.0;
      bb->x_max = 1.0;
      bb->y_min = 0.0;
      bb->y_max = 1.0;
      bb->z_min = 0.0;
      bb->z_max = 1.0;
    }
  else
    {
      bb->x_min = 1.0e10;
      bb->x_max =-1.0e10;
      bb->y_min = 1.0e10;
      bb->y_max =-1.0e10;
      bb->z_min = 1.0e10;
      bb->z_max =-1.0e10;
      for (over3d_link = over3d_list->first; over3d_link != NULL; 
	   over3d_link = over3d_link->next)
	{
	  over3d = over3d_link->data;
	  bb->x_min = real_min(bb->x_min, over3d->bb->x_min);
	  bb->x_max = real_max(bb->x_max, over3d->bb->x_max);
	  bb->y_min = real_min(bb->y_min, over3d->bb->y_min);
	  bb->y_max = real_max(bb->y_max, over3d->bb->y_max);
	  bb->z_min = real_min(bb->z_min, over3d->bb->z_min);
	  bb->z_max = real_max(bb->z_max, over3d->bb->z_max);
	}
    } /* end if over3d_list */
}

