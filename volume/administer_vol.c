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
#include "volume_internal.h"

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

volume_mapping *
new_volume_mapping( char *name){
  volume_mapping *volume;
  int is, id, i, j;
  static int n_calls=0;

  n_calls++;

  volume = (volume_mapping *) malloc(sizeof(volume_mapping));
  
/* initialize all fields */
  volume->used_by_surface = 0;
  volume->no_hole = FALSE;

/* copy the name */
  volume->name = (char *) malloc( (strlen(name)+1)*sizeof(char) );
  volume->name = strcpy( volume->name, name );

/* save the pointers to the data and the functions for
   evaluating the mapping */
  volume->inverse_known = 0;
  volume->forward_mapping = NULL;
  volume->inverse_mapping = NULL;
  volume->set_volume = NULL;
  volume->delete_volume_data = NULL;
  volume->write_volume_data = NULL;
  volume->volume_data_ptr = NULL;

/* allocate space for the rotation */
  volume->rotation_ = create_real_array_2d(3,3);

/* initially, the mapping is not transformed */
  volume->transformed = FALSE;

/* initialize the rotation */
  for (i=1; i<=3; i++){
    for (j=1; j<=3; j++){
      rot_vol(i,j) = 0.0;
    }
    rot_vol(i,i) = 1.0;
  }
  volume->axis[0] = 1.0;
  volume->axis[1] = 0.0;
  volume->axis[2] = 0.0;
  volume->angle   = 0.0;

/* initialize the translation */
  volume->translation[0] = 0.0;
  volume->translation[1] = 0.0;
  volume->translation[2] = 0.0;

  volume->pre_trans[0] = 0.0;
  volume->pre_trans[1] = 0.0;
  volume->pre_trans[2] = 0.0;

/* set the number of grid lines */
  volume->r1_points = 0;
  volume->r2_points = 0;
  volume->r3_points = 0;

  volume->r1_period = 0;
  volume->r2_period = 0;
  volume->r3_period = 0;

  volume->bc_ptr   = create_int_array_2d(2,3);
  volume->surf_ptr = create_int_array_2d(2,3);
  volume->edge_curve_ = create_int_array_1d(12);

  volume->bb = (bounding_box *) malloc( sizeof(bounding_box) );
  volume->bb->x_min = 0.0;
  volume->bb->x_max = 1.0;
  volume->bb->y_min = 0.0;
  volume->bb->y_max = 1.0;
  volume->bb->z_min = 0.0;
  volume->bb->z_max = 1.0;

/* set the grid type and the periodicity flag */
  volume->type= NULL; /* unknown grid type */

/* set the default plot mode: thick edges, parameter arrows and a coordinate cage */
  volume->plot_mode = 2 | 4 | 256;

/* set the color */
  volume->color = n_calls % 8;

/* initialize the parameter surfaces to be plotted */
  volume->r1_surface = 0.0;
  volume->r2_surface = 0.0;
  volume->r3_surface = 0.0;

  volume->plot_surf_label = 0;
  volume->plot_bc         = 0;

/* initialize the 2-d arrays */
  for (is=1; is<=2; is++)
    for (id=1; id<=3; id++){
      bc_vol(is,id) = 0;
      surf_vol(is,id) = 0;
    }

/* initialize the edge-array */
  for (i=1; i<=12; i++)
    edge_vol(i) = 0;

  return volume;
}

volume_mapping *
delete_volume_mapping( volume_mapping *volume ){
  if (volume == NULL) return NULL;

/* check if there is a surface mapping that relies on this */
/* volume mapping */
  if (volume->used_by_surface) return volume;

/* delete the specific data */
  delete_specific_volume( volume );

/* free name and bounding box */
  free( volume->name );
  free( volume->bb );

/* free the arrays */
  volume->rotation_ = delete_real_array_2d(volume->rotation_);

  volume->bc_ptr    = delete_int_array_2d(volume->bc_ptr);
  volume->surf_ptr  = delete_int_array_2d(volume->surf_ptr);
  volume->edge_curve_ = delete_int_array_1d(volume->edge_curve_);

/* free the structure itself */
  free( volume );
  return NULL;
}

int 
forward_volume_mapping( grid_point *gp_ptr, volume_mapping *volume ){
  grid_point gp;
  int msg;
  real x,y,z;

  if (volume->transformed){
    gp.r = gp_ptr->r;
    gp.s = gp_ptr->s;
    gp.t = gp_ptr->t;

/* take care of any errors in the specific mapping routine and pass them on */
    msg = (*volume->forward_mapping)( &gp, volume->volume_data_ptr );

/* pre-translation */
    x = gp.x + volume->pre_trans[0];
    y = gp.y + volume->pre_trans[1];
    z = gp.z + volume->pre_trans[2];

/* rotate and translate coordinates */
    gp_ptr->x = volume->translation[0] + x * rot_vol(1,1) + y * rot_vol(1,2) + 
      z * rot_vol(1,3);
    gp_ptr->y = volume->translation[1] + x * rot_vol(2,1) + y * rot_vol(2,2) + 
      z * rot_vol(2,3);
    gp_ptr->z = volume->translation[2] + x * rot_vol(3,1) + y * rot_vol(3,2) + 
      z * rot_vol(3,3);

/* rotate the Jacobian */
    gp_ptr->xr = gp.xr * rot_vol(1,1) + gp.yr * rot_vol(1,2) + gp.zr * rot_vol(1,3);
    gp_ptr->yr = gp.xr * rot_vol(2,1) + gp.yr * rot_vol(2,2) + gp.zr * rot_vol(2,3);
    gp_ptr->zr = gp.xr * rot_vol(3,1) + gp.yr * rot_vol(3,2) + gp.zr * rot_vol(3,3);

    gp_ptr->xs = gp.xs * rot_vol(1,1) + gp.ys * rot_vol(1,2) + gp.zs * rot_vol(1,3);
    gp_ptr->ys = gp.xs * rot_vol(2,1) + gp.ys * rot_vol(2,2) + gp.zs * rot_vol(2,3);
    gp_ptr->zs = gp.xs * rot_vol(3,1) + gp.ys * rot_vol(3,2) + gp.zs * rot_vol(3,3);

    gp_ptr->xt = gp.xt * rot_vol(1,1) + gp.yt * rot_vol(1,2) + gp.zt * rot_vol(1,3);
    gp_ptr->yt = gp.xt * rot_vol(2,1) + gp.yt * rot_vol(2,2) + gp.zt * rot_vol(2,3);
    gp_ptr->zt = gp.xt * rot_vol(3,1) + gp.yt * rot_vol(3,2) + gp.zt * rot_vol(3,3);
  }
  else{
    msg = (*volume->forward_mapping)( gp_ptr, volume->volume_data_ptr );
  }

  return msg;
}


int 
inverse_volume_mapping( grid_point *gp_ptr, volume_mapping *volume ){
  real x, y, z;

  if (volume == NULL || volume->inverse_mapping == NULL){
    printf("ERROR: inverse_volume_mapping was called with inconsistent pointers.\n");
    return -1;
  }

  if (volume->transformed){
/* invert the rotation and / or translation */
    x = gp_ptr->x - volume->translation[0];
    y = gp_ptr->y - volume->translation[1];
    z = gp_ptr->z - volume->translation[2];

    gp_ptr->x = x*rot_vol(1,1) + y*rot_vol(2,1) + z*rot_vol(3,1) - volume->pre_trans[0];
    gp_ptr->y = x*rot_vol(1,2) + y*rot_vol(2,2) + z*rot_vol(3,2) - volume->pre_trans[1];
    gp_ptr->z = x*rot_vol(1,3) + y*rot_vol(2,3) + z*rot_vol(3,3) - volume->pre_trans[2];
  }

/* pass any errors from the specific mapping routine to the calling function */
  return (*volume->inverse_mapping)( gp_ptr, volume->volume_data_ptr );
}


void 
set_volume(input_output *io_ptr, volume_mapping *volume, 
	   surface_mapping_list *surface_list ){
  if (volume == NULL || volume->set_volume == NULL){
    printf("ERROR: set_volume was called with inconsistent pointers.\n");
  }
  else
    (*volume->set_volume)( io_ptr, volume, surface_list );
}

void 
delete_specific_volume( volume_mapping *volume ){
/* don't try to delete before the pointers are assigned */
  if (volume == NULL || volume->delete_volume_data == NULL) return;

  (*volume->delete_volume_data)( volume->volume_data_ptr );
}

void 
write_specific_volume( int fd, volume_mapping *volume ){
  if (volume == NULL || volume->write_volume_data == NULL){
    printf("ERROR: write_specific_volume was called with inconsistent pointers.\n");
  }
  else
    (*volume->write_volume_data)( fd, volume->volume_data_ptr );
}

volume_mapping *
choose_volume(input_output *io_ptr, char *name, surface_mapping_list *surface_mappings){
  char prompt[80];
  int icom, quit;
  volume_mapping *volume=NULL;
  surface_mapping *surface;
  surface_mapping_link *surface_link;

#include "choose_volume_com.h"

  sprintf(prompt,"set volume mapping type %s>", name);

  quit = 0;
  do{

    icom = get_command(io_ptr, prompt, COMMAND, LAST_COM, LEVEL,
		       SAVE_ON_COPY, ARGUMENT);

    switch (icom) {

/* volume grids */
    case CARTESIAN: 
/* make volume */
      volume = new_volume_mapping(name);
/* initialize the cartesian mapping */
      init_cartesian_grid( volume );
      quit = 1;
      break;

    case CYLINDRICAL: 
/* make volume */
      volume = new_volume_mapping(name);
/* initialize the cylindrical mapping */
      init_cylindrical_grid( volume );
      quit = 1;
      break;

    case NORMAL_SPHEROID: 
/* make volume */
      volume = new_volume_mapping(name);
/* initialize the cylindrical mapping */
      init_normal_spheroid( volume );
      quit = 1;
      break;

    case NORMAL_SURFACE:
      if (surface_mappings->n_members <= 0){
	printf("This type of mapping needs a surface, but the surface list is "
	       "empty.\n");
      }
      else if ((surface_link = get_surface( io_ptr, surface_mappings, NO_INDENT )) != NULL){
	surface = surface_link->data;
/* initialize the grid */  
	volume = new_volume_mapping( name );
	init_normal_surface( io_ptr, volume, surface );
	quit = 1;
      }

      break;


    case DISCRETE_POINT: 
/* make volume grid */
      volume = new_volume_mapping(name);
/* read the discrete gridpoints from a file */
      if (init_volume_point( io_ptr, volume ) != NULL){
	quit = 1;
      }
      else{
	volume = delete_volume_mapping( volume );
      }
      break;

    case HYPERBOLIC:
      if (surface_mappings->n_members <= 0){
	printf("Sorry, this type of mapping needs a surface, but the surface list is "
	       "empty.\n");
      }
      else if ((surface_link = get_surface( io_ptr, surface_mappings, NO_INDENT )) != NULL){
	surface = surface_link->data;
/* initialize the grid */  
	volume = new_volume_mapping( name );
	init_hyp_grid( io_ptr, volume, surface );
	quit = 1;
      }

      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL, SAVE_ON_COPY, ARGUMENT)) 
	     == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case CANCEL:
/* cancel, unsucessful surface grid generation */
      return NULL;
      break;

    default:
      break;
    }
  }
  while(!quit);

/* set the view point and focal points to make the new object visible */
  ogl_standard_view( VOLUME_W, volume->bb );

/* call the specific setup routine to enter all the curve specific data */
  set_volume( io_ptr, volume, surface_mappings ); 

/* sucessful curve generation */
  return volume;
}

void
volume_list_bb(bounding_box *bb, volume_mapping_list *volume_list)
{
  volume_mapping_link *volume_link;
  volume_mapping *volume;

  if (!volume_list  || volume_list->n_members == 0)
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
      for (volume_link = volume_list->first; volume_link != NULL; 
	   volume_link = volume_link->next)
	{
	  volume = volume_link->data;
	  bb->x_min = real_min(bb->x_min, volume->bb->x_min);
	  bb->x_max = real_max(bb->x_max, volume->bb->x_max);
	  bb->y_min = real_min(bb->y_min, volume->bb->y_min);
	  bb->y_max = real_max(bb->y_max, volume->bb->y_max);
	  bb->z_min = real_min(bb->z_min, volume->bb->z_min);
	  bb->z_max = real_max(bb->z_max, volume->bb->z_max);
	}
    } /* end if volume_list */
}

volume_mapping_link *
get_volume_mapping( input_output *io_ptr, volume_mapping_list *head, int level){
  volume_mapping *volume;
  volume_mapping_link *this_link;
  char **command;
  int i, icom, ncom, *save_on_copy=NULL, *argument=NULL, quit;

  if (head->n_members == 0){
    printf("The list of volume mappings is empty.\n");
    return NULL;
  }

  ncom = head->n_members + 1;
  command = (char **) malloc( (ncom+1)*sizeof(char*) );

  i=0;
/* add color information to the name! */
  for (this_link = head->first; this_link != NULL; 
       this_link = this_link->next){
    volume = this_link->data;
    command[i] = (char *) malloc( (strlen(volume->name) + 
				   strlen(ogl_color_name(volume->color)) + 3) 
				 * sizeof(char) );
    sprintf( command[i], "%s(%s)", volume->name, ogl_color_name(volume->color) );
    i++;
  }
  command[ncom-1] = "help";
  command[ncom]   = "cancel";

  quit = 0;
  do{
    icom = get_command(io_ptr, "volume mapping name: ", command, ncom, level, 
		       save_on_copy, argument);

/* default is to quit */
    quit = 1;

    if (icom == -1)
/* not a valid command */
      quit = 0;
    else if (icom == ncom-1){
/* help */
      printf(
" o  Enter one of the volume names if you wish to proceed.\n"
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

void
check_jacobian(input_output *io_ptr, volume_mapping *volume){
  grid_point gp, gp_p, gp_m;
  real xr[3], xs[3], xt[3], dr, ds, dt, error, max_err, r_max, s_max, t_max, eps;
  int NO_INDENT = 0, i, j, k;

  eps = real_max( 1.e-10, get_real(io_ptr, "Enter eps for numerical "
				   "differentiation (>0) :", 1.e-4, NO_INDENT));

  dr = 1.0/(volume->r1_points-1);
  ds = 1.0/(volume->r2_points-1);
  dt = 1.0/(volume->r3_points-1);

  max_err = 0.0;
  r_max = -1; s_max = -1; t_max = -1;
  for (i=1; i<=volume->r1_points; i++)
    for (j=1; j<=volume->r2_points; j++)
      for (k=1; k<=volume->r3_points; k++)
	{
	  gp.r = (i-1)*dr; gp.s = (j-1)*ds; gp.t = (k-1)*dt;

/* evaluate the analytic jacobian */
	  forward_volume_mapping(&gp, volume);

/* estimate the jacobian numerically */
	  gp_p.r = gp.r+eps;
	  gp_p.s = gp.s;
	  gp_p.t = gp.t;
	  forward_volume_mapping(&gp_p, volume);
	  gp_m.r = gp.r-eps;
	  gp_m.s = gp.s;
	  gp_m.t = gp.t;
	  forward_volume_mapping(&gp_m, volume);
	  xr[0] = 0.5*(gp_p.x - gp_m.x)/eps;
	  xr[1] = 0.5*(gp_p.y - gp_m.y)/eps;
	  xr[2] = 0.5*(gp_p.z - gp_m.z)/eps;

	  gp_p.r = gp.r;
	  gp_p.s = gp.s+eps;
	  gp_p.t = gp.t;
	  forward_volume_mapping(&gp_p, volume);
	  gp_m.r = gp.r;
	  gp_m.s = gp.s-eps;
	  gp_m.t = gp.t;
	  forward_volume_mapping(&gp_m, volume);
	  xs[0] = 0.5*(gp_p.x - gp_m.x)/eps;
	  xs[1] = 0.5*(gp_p.y - gp_m.y)/eps;
	  xs[2] = 0.5*(gp_p.z - gp_m.z)/eps;

	  gp_p.r = gp.r;
	  gp_p.s = gp.s;
	  gp_p.t = gp.t+eps;
	  forward_volume_mapping(&gp_p, volume);
	  gp_m.r = gp.r;
	  gp_m.s = gp.s;
	  gp_m.t = gp.t-eps;
	  forward_volume_mapping(&gp_m, volume);
	  xt[0] = 0.5*(gp_p.x - gp_m.x)/eps;
	  xt[1] = 0.5*(gp_p.y - gp_m.y)/eps;
	  xt[2] = 0.5*(gp_p.z - gp_m.z)/eps;

	  error = sqrt( sqr(gp.xr-xr[0]) + sqr(gp.yr-xr[1]) + sqr(gp.zr-xr[2]) + 
			sqr(gp.xs-xs[0]) + sqr(gp.ys-xs[1]) + sqr(gp.zs-xs[2]) + 
			sqr(gp.xt-xt[0]) + sqr(gp.yt-xt[1]) + sqr(gp.zt-xt[2]));
	  if (error > max_err)
	    {
	      max_err = error;
	      r_max = gp.r;
	      s_max = gp.s;
	      t_max = gp.t;
	    }
	} /* end for all grid points */
  printf("The largest error in the Jacobian was %e and occured at (%e, %e, %e)\n", 
	 max_err, r_max, s_max, t_max);

  gp.r = gp.s = gp.t = 0.0;
  while( (gp.r = get_real( io_ptr, "r ( exit with r <= -1 ) :", gp.r, NO_INDENT )) > -1)
    {
      gp.s = get_real( io_ptr, "s:", gp.s, NO_INDENT );
      gp.t = get_real( io_ptr, "t:", gp.t, NO_INDENT );
      forward_volume_mapping(&gp, volume);

/* estimate the jacobian numerically */
      gp_p.r = gp.r+eps;
      gp_p.s = gp.s;
      gp_p.t = gp.t;
      forward_volume_mapping(&gp_p, volume);
      gp_m.r = gp.r-eps;
      gp_m.s = gp.s;
      gp_m.t = gp.t;
      forward_volume_mapping(&gp_m, volume);
      xr[0] = 0.5*(gp_p.x - gp_m.x)/eps;
      xr[1] = 0.5*(gp_p.y - gp_m.y)/eps;
      xr[2] = 0.5*(gp_p.z - gp_m.z)/eps;

      gp_p.r = gp.r;
      gp_p.s = gp.s+eps;
      gp_p.t = gp.t;
      forward_volume_mapping(&gp_p, volume);
      gp_m.r = gp.r;
      gp_m.s = gp.s-eps;
      gp_m.t = gp.t;
      forward_volume_mapping(&gp_m, volume);
      xs[0] = 0.5*(gp_p.x - gp_m.x)/eps;
      xs[1] = 0.5*(gp_p.y - gp_m.y)/eps;
      xs[2] = 0.5*(gp_p.z - gp_m.z)/eps;

      gp_p.r = gp.r;
      gp_p.s = gp.s;
      gp_p.t = gp.t+eps;
      forward_volume_mapping(&gp_p, volume);
      gp_m.r = gp.r;
      gp_m.s = gp.s;
      gp_m.t = gp.t-eps;
      forward_volume_mapping(&gp_m, volume);
      xt[0] = 0.5*(gp_p.x - gp_m.x)/eps;
      xt[1] = 0.5*(gp_p.y - gp_m.y)/eps;
      xt[2] = 0.5*(gp_p.z - gp_m.z)/eps;

/* print the result */
      printf("Analytic jacobian:\n"
	     "xr = ( %e, %e, %e )\nxs = ( %e, %e, %e )\nxt = ( %e, %e, %e )\n\n", 
	     gp.xr, gp.yr, gp.zr, gp.xs, gp.ys, gp.zs, gp.xt, gp.yt, gp.zt);
      printf("Difference between analytic and numerical jacobian:\n"
	     "xr = ( %e, %e, %e )\nxs = ( %e, %e, %e )\nxt = ( %e, %e, %e )\n\n", 
	     gp.xr-xr[0], gp.yr-xr[1], gp.zr-xr[2], 
	     gp.xs-xs[0], gp.ys-xs[1], gp.zs-xs[2], 
	     gp.xt-xt[0], gp.yt-xt[1], gp.zt-xt[2]);
	     
    }    
}

void
transform_volume(input_output *io_ptr, volume_mapping *volume){
  char prompt[120];
  int quit=FALSE, replot=TRUE, icom, i, j;
  real q2[3][3], q3[3][3], t[3][3], a1, a2, a3, length;
  real_array_2d *q2_, *q3_, *t_, *new_rot_;
#define q2(i,j) compute_index_2d(q2_, i, j)  
#define q3(i,j) compute_index_2d(q3_, i, j)  
#define t(i,j)  compute_index_2d(t_ , i, j)  
#define new_rot(i,j)  compute_index_2d(new_rot_ , i, j)  

#include "transform_volume_com.h"

  q2_ = create_real_array_2d(3,3);
  q3_ = create_real_array_2d(3,3);
  t_  = create_real_array_2d(3,3);
  new_rot_  = create_real_array_2d(3,3);

  sprintf(prompt,"transform volume mapping `%s'>", volume->name);

  do{

    if (replot){
/* plot the surface grid */
      if (ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
			 ogl_length_scale(volume->bb))){
	draw_volume( volume, volume->plot_mode, NULL);
	if (volume->plot_mode & 256) ogl_bb_cage( volume->bb, OGL_WHITE );
	ogl_end_plot();
	replot = 0;
      }
    }

    icom = get_command(io_ptr, prompt, COMMAND, LAST_COM, LEVEL,
		       SAVE_ON_COPY, ARGUMENT);

    switch (icom) {

/* volume grids */
    case ROT_AXIS: 
      do {
	a1 = volume->axis[0] = get_real( io_ptr, "axis(x):", volume->axis[0], 
					 NO_INDENT );
	a2 = volume->axis[1] = get_real( io_ptr, "axis(y):", volume->axis[1], 
					 NO_INDENT );
	a3 = volume->axis[2] = get_real( io_ptr, "axis(z):", volume->axis[2], 
					 NO_INDENT );
	length = sqrt(sqr(a1)+sqr(a2)+sqr(a3));
	if (length < 1.e-5)
	  printf("Sorry, the axis must have a length > 1.e-5. Please try again.\n");
      } while(length < 1.e-5);
      a1 /= length; a2 /= length; a3 /= length; 

      volume->angle   = get_real( io_ptr, "angle (degrees):", 180.0*volume->angle/
				 ((double)M_PI), NO_INDENT ) * ((double)M_PI)/180.0;

/* recompute the rotation matrix... */
      q2(1,1) = cos(volume->angle);
      q2(1,2) =-sin(volume->angle);
      q2(1,3) = 0.0;
      q2(2,1) = sin(volume->angle);
      q2(2,2) = cos(volume->angle);
      q2(2,3) = 0.0;
      q2(3,1) = 0.0;
      q2(3,2) = 0.0;
      q2(3,3) = 1.0;

      q3(1,3) = a1;
      q3(2,3) = a2;
      q3(3,3) = a3;
      if (fabs(a3) >= real_max(fabs(a2), fabs(a1))){
	  q3(1,1) = 1.0;
	  q3(2,1) = 1.0;
	  q3(3,1) = -(a1+a2)/a3;
	}
      else if (fabs(a2) >= real_max(fabs(a1), fabs(a3)))
	{
	  q3(1,1) = 1.0;
	  q3(2,1) = -(a1+a3)/a2;
	  q3(3,1) = 1.0;
	}
      else
	{
	  q3(1,1) = -(a2+a3)/a1;
	  q3(2,1) = 1.0;
	  q3(3,1) = 1.0;
	}
/* normalize */
      length = sqrt(sqr(q3(1,1)) + sqr(q3(2,1)) + sqr(q3(3,1)));
      for (i=1; i<=3; i++)
	q3(i,1) /= length;

/* second column must be orthogonal to both the first and the third */
/* lets make a cross product */
      q3(1,2) = q3(2,3)*q3(3,1) - q3(3,3)*q3(2,1);
      q3(2,2) = q3(3,3)*q3(1,1) - q3(1,3)*q3(3,1);
      q3(3,2) = q3(1,3)*q3(2,1) - q3(2,3)*q3(1,1);

      length = 0.0;
      for (i=1; i<=3; i++)
	length += sqr(q3(i,2));
      for (i=1; i<=3; i++)
	q3(i,2) /= length;

/* We first multiply q3 by q2 */
      for (i=1; i<=3; i++)
	for (j=1; j<=3; j++)
	  {
	    t(i,j) = q3(i,1)*q2(1,j) + q3(i,2)*q2(2,j) + q3(i,3)*q2(3,j);
	  }
/* we then multiply t by q3^T */
      for (i=1; i<=3; i++)
	for (j=1; j<=3; j++)
	  {
	    new_rot(i,j)  = t(i,1)*q3(j,1) + t(i,2)*q3(j,2) + t(i,3)*q3(j,3);
	  }

/* pre-multiply the new rotation matrix with the existing one. Store it temporarily */
/* in q3 */
      for (i=1; i<=3; i++)
	for (j=1; j<=3; j++)
	  {
	    q3(i,j) = new_rot(i,1)*rot_vol(1,j) + new_rot(i,2)*rot_vol(2,j) + 
	      new_rot(i,3)*rot_vol(3,j);
	  }

/* save the new rotation matrix in rot_vol */
      for (i=1; i<=3; i++)
	for (j=1; j<=3; j++)
	  {
	    rot_vol(i,j) = q3(i,j);
	  }

/* set the transformation flag */
      volume->transformed = TRUE;
/* update the bounding box */
      volume_bb( volume );
/* redraw the grid */
      replot = TRUE;
      break;

    case PRE_TRANS: 
      volume->pre_trans[0] = get_real( io_ptr, "x-translation:", 
				       volume->pre_trans[0], NO_INDENT );
      volume->pre_trans[1] = get_real( io_ptr, "y-translation:", 
				       volume->pre_trans[1], NO_INDENT );
      volume->pre_trans[2] = get_real( io_ptr, "z-translation:", 
				       volume->pre_trans[2], NO_INDENT );
/* set the transformation flag */
      volume->transformed = TRUE;
/* update the bounding box */
      volume_bb( volume );
/* redraw the grid */
      replot = TRUE;
      break;

    case POST_TRANS: 
      volume->translation[0] = get_real( io_ptr, "x-translation:", 
					 volume->translation[0], NO_INDENT );
      volume->translation[1] = get_real( io_ptr, "y-translation:", 
					 volume->translation[1], NO_INDENT );
      volume->translation[2] = get_real( io_ptr, "z-translation:", 
					 volume->translation[2], NO_INDENT );
/* set the transformation flag */
      volume->transformed = TRUE;
/* update the bounding box */
      volume_bb( volume );
/* redraw the grid */
      replot = TRUE;
      break;

    case RESET_TRANS: 
      volume->axis[0] = 1.0;
      volume->axis[1] = 0.0;
      volume->axis[2] = 0.0;
      volume->angle   = 0.0;
      for (i=1; i<=3; i++){
	for (j=1; j<=3; j++){
	  rot_vol(i,j) = 0.0;
	}
	rot_vol(i,i) = 1.0;
	volume->translation[i-1] = 0.0;
	volume->pre_trans[i-1] = 0.0;
      }
      volume->transformed = FALSE;
/* update the bounding box */
      volume_bb( volume );
      replot = TRUE;
      break;

    case SHOW: 
      printf("\nRotation matrix:\n(% e  % e  % e)\n(% e  % e  % e)\n(% e  % e  % e)\n\n",
	     rot_vol(1,1), rot_vol(1,2), rot_vol(1,3), 
	     rot_vol(2,1), rot_vol(2,2), rot_vol(2,3), 
	     rot_vol(3,1), rot_vol(3,2), rot_vol(3,3));
      printf("Translation before rotation (x,y,z) = (%e, %e, %e)\n\n", 
	     volume->pre_trans[0], volume->pre_trans[1], volume->pre_trans[2]);
      printf("Translation after rotation (x,y,z) = (%e, %e, %e)\n\n", 
	     volume->translation[0], volume->translation[1], volume->translation[2]);
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL, SAVE_ON_COPY, ARGUMENT)) 
	     == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
      quit = TRUE;
      break;

    }
  }
  while (!quit);

  delete_real_array_2d(q2_);
  delete_real_array_2d(q3_);
  delete_real_array_2d(t_);
  delete_real_array_2d(new_rot_);
}

static void
update_bb( bounding_box *bb, grid_point *gp ){
  bb->x_min = real_min( bb->x_min, gp->x );
  bb->x_max = real_max( bb->x_max, gp->x );
  bb->y_min = real_min( bb->y_min, gp->y );
  bb->y_max = real_max( bb->y_max, gp->y );
  bb->z_min = real_min( bb->z_min, gp->z );
  bb->z_max = real_max( bb->z_max, gp->z );
}

void
volume_bb( volume_mapping *volume ){
  grid_point gp;
  const real step=0.05;

/* initialize bounding box */
  volume->bb->x_min =  1.e10;
  volume->bb->x_max = -1.e10;
  volume->bb->y_min =  1.e10;
  volume->bb->y_max = -1.e10;
  volume->bb->z_min =  1.e10;
  volume->bb->z_max = -1.e10;

/* near surface */
  gp.t = 0.0;

  gp.s = 0.0;
  for (gp.r=0.0; gp.r<=1.0+0.5*step; gp.r+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

  gp.s = 1.0;
  for (gp.r=0.0; gp.r<=1.0+0.5*step; gp.r+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

  gp.r = 0.0;
  for (gp.s=0.0; gp.s<=1.0+0.5*step; gp.s+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

  gp.r = 1.0;
  for (gp.s=0.0; gp.s<=1.0+0.5*step; gp.s+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

/* far surface */
  gp.t = 1.0;

  gp.s = 0.0;
  for (gp.r=0.0; gp.r<=1.0+0.5*step; gp.r+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

  gp.s = 1.0;
  for (gp.r=0.0; gp.r<=1.0+0.5*step; gp.r+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

  gp.r = 0.0;
  for (gp.s=0.0; gp.s<=1.0+0.5*step; gp.s+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

  gp.r = 1.0;
  for (gp.s=0.0; gp.s<=1.0+0.5*step; gp.s+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

/* edges between the near and far surfaces */
  gp.r = 0.0;
  gp.s = 0.0;
  for (gp.t=0.0; gp.t<=1.0+0.5*step; gp.t+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

  gp.r = 1.0;
  gp.s = 0.0;
  for (gp.t=0.0; gp.t<=1.0+0.5*step; gp.t+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

  gp.r = 0.0;
  gp.s = 1.0;
  for (gp.t=0.0; gp.t<=1.0+0.5*step; gp.t+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

  gp.r = 1.0;
  gp.s = 1.0;
  for (gp.t=0.0; gp.t<=1.0+0.5*step; gp.t+=step){
    forward_volume_mapping( &gp, volume );
    update_bb( volume->bb, &gp );
  }

}

void
file_volume( input_output *io_ptr, volume_mapping *volume ){
  FILE *fp;
  char *fname;
  grid_point gp;
  int i, j, k;
  real dr, ds, dt;

  fp = open_ascii_file( io_ptr, "Enter file name: ", "test.plot3d", &fname, 'w', 
			0, TRUE );

  if (!fp) return;

  dr = 1.0/(volume->r1_points - 1);
  ds = 1.0/(volume->r2_points - 1);
  dt = 1.0/(volume->r3_points - 1);

  fprintf( fp, "%i %i %i\n", volume->r1_points, volume->r2_points, volume->r3_points ); 
/* read the x coordinates */
  for (k = 1; k <= volume->r3_points; k++)
    for (j = 1; j <= volume->r2_points; j++)
      for (i = 1; i <= volume->r1_points; i++){
	gp.r = (i-1)*dr;
	gp.s = (j-1)*ds;
	gp.t = (k-1)*dt;
	forward_volume_mapping( &gp, volume );
	fprintf( fp, "%e\n", gp.x);
      }  

/* read the y coordinates */
  for (k = 1; k <= volume->r3_points; k++)
    for (j = 1; j <= volume->r2_points; j++)
      for (i = 1; i <= volume->r1_points; i++){
	gp.r = (i-1)*dr;
	gp.s = (j-1)*ds;
	gp.t = (k-1)*dt;
	forward_volume_mapping( &gp, volume );
	fprintf( fp, "%e\n", gp.y);
      }

/* read the z coordinates */
  for (k = 1; k <= volume->r3_points; k++)
    for (j = 1; j <= volume->r2_points; j++)
      for (i = 1; i <= volume->r1_points; i++){
	gp.r = (i-1)*dr;
	gp.s = (j-1)*ds;
	gp.t = (k-1)*dt;
	forward_volume_mapping( &gp, volume );
	fprintf( fp, "%e\n", gp.z);
      }
  
/* successful completion */
  fclose( fp );

}
