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
#include "overlap_internal.h"

/* static prototypes */
static int
on_physical_boundary(int i, int j, int k, component_vgrid *vgrid, int *dir);
static int
get_surface_label(int i, int j, int k, component_vgrid *vgrid, int *side, int *dir);
static int
get_curve_label(int i, int j, int k, component_vgrid *vgrid);
static void
recompute_one_interp(int present_surf, component_vgrid *vgrid, interp_3d_point *interp);
static void
recompute_interpolation_location(/* FILE *err_msg,*/ interp_3d_point *interp, 
				 component_vgrid *vgrid, 
				 overlapping_3d_grid *over3d, int dir, 
				 int iw1, real x_dist, real y_dist, real z_dist);
static int
boundary_diff(int present_surf, component_vgrid *vgrid, interp_3d_point *interp,
	      real *x_dist, real *y_dist, real *z_dist);
static void
copy_surface_holes(discrete_surface *sgrid_ptr);
static void
compute_hole_curves(overlapping_3d_grid *over3d);
static hole_surface *
get_existing_hs( overlapping_3d_grid *over3d, int surface_label );
static hole_surface *
new_bounding_surface(component_vgrid *vgrid);
static int
inside_global_domain(real xp, real yp, real zp, int surface_label, 
		     overlapping_3d_grid *over3d, int debug);
static void
mark_around_tri(triangle tri, interp_3d_point *interp[3],
		int surface_label, component_vgrid *vgrid, 
		overlapping_3d_grid *over3d, interp_3d_point *initial_guess, int level);
static int 
cut_volume_holes(input_output *io_ptr, overlapping_3d_grid *over3d);
static void 
vgrid_face_coordinates(component_vgrid *vgrid, int side, int dir, 
		       real_array_2d **x_surf_, 
		       real_array_2d **y_surf_, 
		       real_array_2d **z_surf_);
static void 
show_overlap_param( overlapping_3d_grid *over3d );
static void 
change_3d_list( input_output *io_ptr, overlapping_3d_grid *over3d, 
	       linked_list *all_volumes );
static void
show_overlap_list( overlapping_3d_grid *over3d );
static int 
mark_needed_interp(overlapping_3d_grid *over3d);
static int 
needed_by_disc(int i_point, int j_point, int k_point, component_vgrid *vgrid,
	       overlapping_3d_grid *over3d);
static void 
change_sign(int i_loc, int j_loc, int k_loc, int iw1, component_vgrid *vgrid);
static int 
overlap_3d(input_output *io_ptr, overlapping_3d_grid *over3d);
static int
interior_in_i(int i, int j, int k, component_vgrid *vgrid, 
	      int dw2, int nw1, int tw2);
static int
interior_in_j(int i, int j, int k, component_vgrid *vgrid, 
	      int dw2, int nw1, int tw2);
static int
interior_in_k(int i, int j, int k, component_vgrid *vgrid, 
	      int dw2, int nw1, int tw2);
static interp_3d_point *
interp_from_lower(real xp, real yp, real zp, 
		  int i_point, int j_point, int k_point,
		  interp_3d_point *initial_guess, int query_surface, int query_curve,
		  component_vgrid *own_vgrid, 
		  int lower_than, overlapping_3d_grid *over3d);
static interp_3d_point *
interp_from_higher(real xp, real yp, real zp, 
		   int i_point, int j_point, int k_point,
		   interp_3d_point *initial_guess, int query_surface, int query_curve,
		   linked_list_member *other_link, overlapping_3d_grid *over3d);
static int
disc_point(int i, int j, int k, component_vgrid *vgrid, overlapping_3d_grid *over3d);
static hole_surface *
get_or_make_hs( overlapping_3d_grid *over3d, int surface_label );
static void
prepare_hole_surfaces(input_output *io_ptr, overlapping_3d_grid *over3d );
static int 
initialize_flag(overlapping_3d_grid *over3d);
/* end static prototypes */


void 
set_3d_overlap(input_output *io_ptr, overlapping_3d_grid *over3d, 
	       linked_list *all_volumes ){
  int quit=0, replot=0, replot_patch=0, icom;
  linked_list_member *volume_link, *patch_link;
  hole_surface *patch;

#include "set_3d_overlap_com.h"

  do{

/* plot the overlapping 3-d grid */
    if (replot && ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, 
				 ogl_length_scale(over3d->bb))){ 
/* highlight grid boundaries */
      if (over3d->plot_mode & 2){ 
	for (volume_link = over3d->mapping_list->first; volume_link != NULL; 
	     volume_link = volume_link->next){ 
	  draw_volume( (volume_mapping *) volume_link->data, 2, NULL ); 
	} 
      } 
      draw_overlapping_3d_grid( over3d, over3d->plot_mode ); 
      if (over3d->plot_mode & 256) ogl_bb_cage( over3d->bb, OGL_WHITE ); 
      ogl_end_plot(); 
      replot = FALSE; 
    } 

/* redraw the overlapping surface grids */
    if (replot_patch && ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, 
				       ogl_length_scale(over3d->bb))){
      for (patch_link = over3d->surface_grid_list->first; patch_link != NULL; 
	   patch_link = patch_link->next){
	patch = patch_link->data;
	draw_hole_surface(patch, patch->plot_mode);
      }
      if (over3d->plot_mode & 256) ogl_bb_cage( over3d->bb, OGL_WHITE );
      ogl_end_plot();
      replot_patch = FALSE;
    }

    switch( get_command( io_ptr, "overlapping grid>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

/* change the overlap parameters */
    case PARAMS:
      overlap_3d_parameters( io_ptr, over3d );
      break;

    case PHYSICAL_BNDRY:
      over3d->show_boundary = 
	get_yes_no(io_ptr, "Display the global boundary during the "
		   "overlap algorithm? ", LEVEL, TRUE);
      break;

    case MISMATCH:
      over3d->correct_mismatch = 
	get_yes_no(io_ptr, "Correct interpolation locations for boundary mismatch? ", 
		   LEVEL, TRUE);
      break;

    case SHOW_HOLES:
      over3d->show_volume_holes = 
	get_yes_no(io_ptr, "Display the result of the hole-cutting during the "
		   "overlap algorithm? ", LEVEL, TRUE);
      break;

/* toggle trimming style */
    case TRIM_STYLE: 
      over3d->trim_style = !over3d->trim_style; 
      printf("The trimming style is to %s the overlap.\n", 
 	     (over3d->trim_style)? "maximize": "minimize"); 
      break; 

/* change the list of components in the overlapping grid */
    case LIST:
      change_3d_list( io_ptr, over3d, all_volumes );
/* can't guarentee that the overlap still is valid after the list */
/* has been changed */
      over3d->valid = 0;
      break;

/* apply the overlap algorithm to compute the intergrid communication data */
    case COMPUTE:

      if (overlap_3d(io_ptr, over3d)){
/* successful */
	printf("\nThe overlap algorithm completed successfully!\n\n");
/* plot interpolation faces, the coordinate cage, and positive surface-label faces. */
	over3d->plot_mode = 64 | 256 | 4096;
	replot = 1;
      }
      else{
/* unsuccessful */
	printf("\n"
"Error: The overlap algorithm failed. Follow the following steps\n"
"to improve the situation:\n"
"\n"
" o  Make sure that all parts of the physical domain are covered by grids.\n"
"\n"
" o  Make sure that the same positive surface label has been assigned\n"
"     to all faces of the grids that are aligned with the same part of\n"
"    the physical boundary.\n"
"\n"
" o  Make sure that a positive edge label is assigned to edges that\n"
"    should cut holes in faces of other grids.\n"
"\n"
" o  Make sure there are sufficiently many grid points in the overlap\n"
"    regions.\n"
"\n"
" o  If you use explicit interpolation, you can begin by making a grid\n"
"    with implicit interpolation, which requires less overlap between\n"
"    the grids. You can also decrease the discretization and interpolation\n"
"    widths to decrease the required overlap.\n"
"\n"
	       );
      }
      break;

/* show the current parameters and list */
    case SHOW:
      show_overlap_list( over3d );
      printf("\n");
      show_overlap_param( over3d );
      printf("\n");
      printf("The trimming style is to %s the overlap.\n", 
	     (over3d->trim_style)? "maximize": "minimize");
      printf("Interpolation locations %s be compensated for boundary mismatch.\n", 
	     (over3d->correct_mismatch)? "WILL": "WILL NOT");

/*       printf("\n"); */
/*       printf("This overlapping grid will be designed to be %s.\n",  */
/* 	     (over3d->overlapping_sub_grid)?  */
/* 	     "a possibly incomplete\npart of a larger overlapping grid": "complete"); */
      break;

/* change the plot mode */
    case PLOT_MODE:
      over3d->plot_mode = overlapping_3d_plot_mode(io_ptr, over3d, over3d->plot_mode);
/* don't need to redraw the volume grids! */
      break;

    case VERBOSE:
      over3d->verbose = int_max(0, get_int(io_ptr, "How verbose >= 0: ", 
					   over3d->verbose, 2));
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL, NULL, NULL)) 
	     == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case INSPECT_OVER_SURF:
      if ((patch_link = get_hole_surface( io_ptr, over3d->surface_grid_list))){
	patch = patch_link->data;
	patch->plot_mode = hole_plot_mode(io_ptr, patch, patch->plot_mode );
	replot_patch = TRUE;
      }
/* don't need to redraw the surface grids! */
      break;

    case EXIT:
      quit = 1;
      break;

    default:
      break;

    }
  } while (!quit);
}

static void 
change_3d_list( input_output *io_ptr, overlapping_3d_grid *over3d, 
	       linked_list *all_volumes )
{
  int replot = 1, quit = 0, icom;
  linked_list_member *this_link, *old_link, *volume_link;
  volume_mapping *volume;

#include "change_3d_list_com.h"

  do{

    if (replot){
/* plot the overlapping grid */
      if (ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, 
			 ogl_length_scale(over3d->bb))){
/* draw the outline of the active component grids and the bounding box cage */
	for (volume_link = over3d->mapping_list->first; volume_link != NULL;
	     volume_link = volume_link->next){
	  draw_volume( (volume_mapping *) volume_link->data, 2, NULL);
	}
	ogl_bb_cage( over3d->bb, OGL_WHITE );
	ogl_end_plot();
	replot = 0;
      }
    }

    switch( get_command( io_ptr, "change-list>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case INSERT_FIRST:
      if ((old_link = get_volume_mapping( io_ptr, all_volumes, NO_INDENT )) != NULL){
/* check if this component is already present */
	for (this_link = over3d->mapping_list->first; this_link != NULL;
	     this_link = this_link->next){
	  if (this_link->data == old_link->data){
	    printf("Sorry, that grid is already in the list!\n");
	    break;
	  }
	}
/* make a new link, insert it fist into the list and copy the pointer to the data */
	if (this_link == NULL){
	  new_link( over3d->mapping_list )->data = volume = old_link->data;
/* update the plotting window */
	  replot = 1;
	}
/* update bounding box */
	volume_list_bb( over3d->bb, over3d->mapping_list );
      }
      break;

    case MOVE_FIRST:
      if ((old_link = get_volume_mapping( io_ptr, over3d->mapping_list, NO_INDENT )) != NULL){
/* save the pointer to the data block */
	volume = old_link->data;
/* delete the old link */
	delete_link( old_link, over3d->mapping_list );
/* make a new link, insert it first in the list and copy the pointer to the data */
	new_link( over3d->mapping_list )->data = volume;
/* update the plotting window */
	replot = 1;
      }
      break;

    case MOVE_LAST:
      if ((old_link = get_volume_mapping( io_ptr, over3d->mapping_list, NO_INDENT )) != NULL){
/* save the pointer to the data block */
	volume = old_link->data;
/* delete the old link */
	delete_link( old_link, over3d->mapping_list );
/* make a new link, insert it last in the list and copy the pointer to the data */
	new_last_link( over3d->mapping_list )->data = volume;
/* update the plotting window */
	replot = 1;
      }
      break;

    case REMOVE:
      if ((this_link = get_volume_mapping( io_ptr, over3d->mapping_list, NO_INDENT )) != NULL){
	volume = this_link->data;
	delete_link( this_link, over3d->mapping_list );
/* update the plotting window */
	replot = 1;
/* update bounding box */
	volume_list_bb( over3d->bb, over3d->mapping_list );
      }
      break;


    case INSERT_ALL:
/* make new links and copy the pointers to the data from the links in the */
/* all_volumes list */
      for (this_link = all_volumes->last; this_link != NULL;
	   this_link = this_link->prev){
	volume = this_link->data;
	new_link( over3d->mapping_list )->data = volume;
      }
/* update the plotting window */
      replot = 1;
/* update bounding box */
      volume_list_bb( over3d->bb, over3d->mapping_list );

      break;

/* show the current parameters and list */
    case SHOW:
      show_overlap_list( over3d );
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL, NULL, NULL)) 
	     == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
      quit = 1;
      break;

    default:
      break;

    }
  } while (!quit);
}

static void
show_overlap_list( overlapping_3d_grid *over3d )
{
  linked_list_member *this_link;
  volume_mapping *volume;
  int grid_points=0;

  printf("The following volume grids are used in the overlapping 3-d grid:\n"
	 "(The first grid has the highest priority)\n"
	 "Color\tName\tType\tGrid points\n");
  for (this_link = over3d->mapping_list->first; this_link != NULL; 
       this_link = this_link->next){
    volume = this_link->data;
    grid_points += volume->r1_points*volume->r2_points*volume->r3_points;
    printf("%s\t%s\t%s\t(%i x %i x %i) = %i\n", ogl_color_name(volume->color),
	   volume->name, volume->type, volume->r1_points, volume->r2_points, 
	   volume->r3_points, volume->r1_points*volume->r2_points*volume->r3_points);
  }
  printf("\nTotal number of grid points: %i\n", grid_points);
}


void
prepare_3d_overlap( overlapping_3d_grid *over3d )
{
  volume_mapping *volume;
  component_vgrid *vgrid;
  int priority;
  linked_list_member *vgrid_link, *patch_link;

  if (over3d->verbose >= 2)
    printf("Deallocating any previous component grids...\n");

/* delete any pre-existing component_vgrid's */
  if (over3d->grid_list->n_members > 0){
    for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL;
	 vgrid_link = vgrid_link->next){
      delete_vgrid( (component_vgrid *) vgrid_link->data );
    }
    delete_linked_list( over3d->grid_list );
/* make a new linked list */
    over3d->grid_list = new_linked_list();
  }

/* delete the list of hole surfaces */
  if (over3d->surface_grid_list->n_members > 0){
    for (patch_link = over3d->surface_grid_list->first; patch_link != NULL;
	 patch_link = patch_link->next){
      delete_hole_surface( (hole_surface *) patch_link->data );
    }
    delete_linked_list( over3d->surface_grid_list );
/* make a new linked list */
    over3d->surface_grid_list = new_linked_list();
  }
  if (over3d->verbose >= 2)
    printf("Done\n");

  if (over3d->verbose >= 2)
    printf("Evaluating the mapping functions at the grid points and initializing the\n"
	   "data structure...\n");
  fflush(stdout);

  priority = 0;

/* insert component grids corresponding to the mappings in the list */
  for (vgrid_link = over3d->mapping_list->last; vgrid_link != NULL; 
       vgrid_link = vgrid_link->prev){
    volume = vgrid_link->data;

/* make a new component_vgrid by using the volume mapping as a template */
    vgrid = new_vgrid( volume );

/* zero basic missfit tolerance */
    vgrid->max_s_dist = 0.0;

/* insert it first in the list */
    new_link( over3d->grid_list)->data = vgrid;

/* assign the priority and color index */
    vgrid->priority = ++priority;
    vgrid->verbose = over3d->verbose;

/* allocate space for the large arrays and initialize them */
/* (pre-existing arrays are first deleted )*/
    init_3d_arrays( volume, vgrid, over3d );

/* construct a bounding surface for the ray-method */
    printf("Constructing a bounding surface for component %s\n", vgrid->name);
    vgrid->bounding_surface = new_bounding_surface(vgrid);

  } /* end for all component vgrids */

  if (over3d->verbose >= 2)
    printf("Done\n");
}

static hole_surface *
new_bounding_surface(component_vgrid *vgrid)
{
  hole_surface *this_node;
  int side, dir;
  discrete_surface *sgrid_ptr;

  /* Make new discrete_surfaces from the grid points of all faces */
  /* of the current vgrid. */
  this_node = new_hole_surface( vgrid->name, 0 );

 
/* check all sides */
    for (side = 1; side <= 2; side++)
      for (dir = 1; dir <= 3; dir++)
	{
/* don't use periodic faces */
	  if ((vgrid->r1_period && dir == 1) ||
	      (vgrid->r2_period && dir == 2) ||
	      (vgrid->r3_period && dir == 3)) continue; 
/* Make a new discrete_surface from the grid points of the current face */
/* of the current vgrid. */
	  sgrid_ptr = new_sgrid_from_vgrid( vgrid, side, dir, 0.0 );

/* insert the new discrete_surface first in the list */
	  insert_discrete_surface( this_node, sgrid_ptr );
	}
  return this_node;
}


static void
prepare_hole_surfaces(input_output *io_ptr, overlapping_3d_grid *over3d )
{
  discrete_surface *sgrid_ptr, *other_sgrid;
  component_vgrid *vgrid=NULL;
  int side, dir, i, j;
  linked_list_member *overg_link;
  discrete_surface_link *sgrid_link, *sgrid2_link;
  component_vgrid_link *vgrid_link;
  hole_surface *overg;
  real x_test[3], max_n_dist;
  inverse_point *surf_point;

  if (over3d->verbose >= 2) printf("\n");

/* loop over all faces in all 3-d components and group faces with matching */
/* surface labels into hole surfaces  */
  for (vgrid_link = over3d->grid_list->last; vgrid_link != NULL;
       vgrid_link = vgrid_link->prev){
    vgrid = vgrid_link->data;

/* check all sides */
    for (side = 1; side <= 2; side++)
      for (dir = 1; dir <= 3; dir++)
	if (surf3(side,dir) != 0){

/* get an existing overlapping surface grid with the current surface label, or if */
/* non-existent, make a new one.. */
	  overg = get_or_make_hs( over3d, surf3(side,dir) );

/* Make a new discrete_surface from the grid points of the current face */
/* of the current vgrid. */
	  sgrid_ptr = new_sgrid_from_vgrid( vgrid, side, dir, 0.0 );

	  if (vgrid->verbose >= 2)
	    printf("Grid `%s' (Surface label = %i):\n"
		   "surface roughness estimated to %e\n", 
		   sgrid_ptr->name, surf3(side, dir), sgrid_ptr->max_n_dist);

/* use the largest sgrid_ptr->max_n_dist as edge and surface tolerances */
/* in the volume grid */
	  vgrid->max_s_dist = real_max( vgrid->max_s_dist, sgrid_ptr->max_n_dist );
	  vgrid->max_e_dist = real_max( vgrid->max_e_dist, sgrid_ptr->max_n_dist );

/* insert the new discrete_surface first in the list */
	  insert_discrete_surface( overg, sgrid_ptr );
	}
/* add extra tolerance */
    vgrid->max_s_dist += over3d->max_b_dist;
    vgrid->max_e_dist += over3d->max_b_dist;
  }

/* set the surface tolerance for the hole_surface to the largest tolerance 
   in any of the discrete surfaces */
  for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;
       overg_link = overg_link->next)
    {
      overg = overg_link->data;
      max_n_dist = 0.0;
      for (sgrid_link = overg->grid_list->first; sgrid_link != NULL;
	   sgrid_link = sgrid_link->next)
	{
	  sgrid_ptr = sgrid_link->data;
	  max_n_dist = real_max(sgrid_ptr->max_n_dist, max_n_dist);
	} /* end for all surface patches */

/* use the largest normal distance as boundary tolerance for the hole_surface */
/* 1.0 is the fudge factor */
      overg->max_b_dist = 1.0*(max_n_dist + over3d->max_b_dist); 
      printf("Tolerance for hole-cutting surface %s is: %e\n", overg->name, 
	     overg->max_b_dist);

    } /* end for all hole surfaces */


/* form the hole cutting curves cooresponding to non-zero edge_curve values */
  compute_hole_curves(over3d);

  for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;
       overg_link = overg_link->next)
    {
      overg = overg_link->data;
      overg->plot_mode = 2 | 64 | 512;

/* remove the surface points inside the hole cutting curves */
      cut_surface_holes(io_ptr, overg);

/* tmp */
/*       for (sgrid_link = overg->grid_list->first; sgrid_link != NULL; */
/* 	   sgrid_link = sgrid_link->next) */
/* 	{ */
/* 	  sgrid_ptr = sgrid_link->data; */
/* 	  printf("\nThis is the flag-array in `%s' after cut_surface_holes\n", sgrid_ptr->name); */
/* 	  print_int_array_2d(sgrid_ptr->flag_ptr); */
/* 	} */
#ifdef SURFACE_HOLE_CUTTING
      printf("After the surface holes have been made\n");
/* plot the hole-surfaces */
      if (ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, 
			 ogl_length_scale(overg->bb))){ 
	draw_hole_surface( overg, 64 | 256 | 1024 ); 
	ogl_end_plot(); 
      } 
/* get the chance to look at the grid from a different view-point */
      overg->plot_mode = hole_plot_mode(io_ptr, overg, 64 | 256 | 1024);
#endif
    }


/* copy the surface holes onto the corresponding faces of the volume grids */
/*    for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;  */
/*         overg_link = overg_link->next)  */
/*      {  */
/*        overg = overg_link->data;  */
/*        for (sgrid_link = overg->grid_list->first; sgrid_link != NULL;  */
/*  	   sgrid_link = sgrid_link->next)  */
/*  	{  */
/*  	  sgrid_ptr = sgrid_link->data;  */
/*  	  copy_surface_holes(sgrid_ptr);  */
/*  	}  */
/*      }  */

/*   for (overg_link = over3d->surface_grid_list->first; overg_link != NULL; */
/*        overg_link = overg_link->next) */
/*     { */
/*       overg = overg_link->data; */
/* triangulate the hole_surface */
/*       triangulate(io_ptr, overg); */
/*     } */
  
  /* check whether the computational domain is on the inside or outside of the */
  /* hole-cutting surface */
  printf("\n");
  for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;
       overg_link = overg_link->next)
    {
      overg = overg_link->data;
      /* find a component grid that has a side with matching surface-label */
      for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL; 
	   vgrid_link = vgrid_link->next)
	{
	  vgrid = vgrid_link->data;
	  if (abs(surf3(1,1)) == overg->surface_label ||
	      abs(surf3(1,2)) == overg->surface_label ||
	      abs(surf3(1,3)) == overg->surface_label ||
	      abs(surf3(2,1)) == overg->surface_label ||
	      abs(surf3(2,2)) == overg->surface_label ||
	      abs(surf3(2,3)) == overg->surface_label)
	    {
	      x_test[0] = x3((range3(1,1)+range3(2,1))/2, (range3(1,2)+range3(2,2))/2,
			     (range3(1,3)+range3(2,3))/2);
	      x_test[1] = y3((range3(1,1)+range3(2,1))/2, (range3(1,2)+range3(2,2))/2,
			     (range3(1,3)+range3(2,3))/2);
	      x_test[2] = z3((range3(1,1)+range3(2,1))/2, (range3(1,2)+range3(2,2))/2,
			     (range3(1,3)+range3(2,3))/2);
	      break;
	    }
	} /* end for all grids */
      if (!vgrid)
	printf("ERROR: prepare_hole_surfaces: Could not find a component "
	       "with surface_label = %i\n", 
	       overg->surface_label);
      else
	{
	  overg->inside = inside_surface(overg, x_test, 0, 0, 0);
	  printf("The computational domain appears to be %s of surface_label %i\n",
		 (overg->inside? "inside": "outside"), overg->surface_label);
	}
    } /* end for all hole-surfaces */
  
}

static void
copy_surface_holes(discrete_surface *sgrid_ptr)
{
  component_vgrid *vgrid;
  discrete_surface_link *face_link;
  discrete_surface *face;
  int i, j;
#define vgrid_face_flag(i,j) compute_index_2d(face->flag_ptr, i, j)

  vgrid = sgrid_ptr->vgrid_mother;
  if (!vgrid || !vgrid->bounding_surface || !vgrid->bounding_surface->grid_list)
    {
      printf("Error: copy_surface_holes: unable to locate the surface grid list\n"
	     "in the mother volume grid of `%s'\n", sgrid_ptr->name);
      return;
    }
  for (face_link = vgrid->bounding_surface->grid_list->first;
       face_link != NULL; face_link = face_link->next)
    {
      face = face_link->data;
      if (face->side == sgrid_ptr->side && face->dir == sgrid_ptr->dir)
	{
	  for (i=1; i <= sgrid_ptr->r_dim; i++)
	    for (j=1; j <= sgrid_ptr->s_dim; j++)
	      vgrid_face_flag(i,j) = flag(i,j);
	  /* tmp */
	  /*  	  printf("\nThis is the flag-array in `%s'\n", face->name);  */
	  /*  	  print_int_array_2d(face->flag_ptr);  */
	}
    }
#undef vgrid_face_flag
}

#define sqr(x) (x)*(x)


static void
compute_hole_curves(overlapping_3d_grid *over3d)
{
  component_vgrid_link *vgrid_link;
  component_vgrid *vgrid;
  hole_cutting_curve *hole_curve;
  hole_surface_link *overg_link;
  hole_surface *overg;
  real *inside_points[2][3], *inside_normal[2][3], n_vec[3];
  int n_points, i, j, k, ip, q, i_const, j_const, k_const;
  real *x_edge, *y_edge, *z_edge, arc, darc, length;

/* estimate the length scale of the grid */
  length = sqrt(sqr(over3d->bb->x_max - over3d->bb->x_min) + 
		sqr(over3d->bb->y_max - over3d->bb->y_min) + 
		sqr(over3d->bb->z_max - over3d->bb->z_min));

  /* loop over all faces in all 3-d components and group faces with matching */
  /* surface labels into hole surfaces  */
  printf("\n");
  for (vgrid_link = over3d->grid_list->last; vgrid_link != NULL;
       vgrid_link = vgrid_link->prev)
    {
      vgrid = vgrid_link->data;

      /* don't do this for dummies */
      /*    if (vgrid->dummy_background) continue;*/

      /* check all edges */
      if (edge_curve(1))	/* r1=0, r2=0 */
	{
	  printf("Constructing a hole-cutting curve from edge 1 of component `%s'\n", 
		 vgrid->name);
	  /* get an existing overlapping surface grid with the current surface label */
	  if ((overg = get_existing_hs( over3d, surf3(1,1) )))
	    {
	      n_points = range3(2,3) - range3(1,3) + 1;
	      x_edge = (real *) malloc( n_points*sizeof(real) );
	      y_edge = (real *) malloc( n_points*sizeof(real) );
	      z_edge = (real *) malloc( n_points*sizeof(real) );
	      for (ip=0; ip<2; ip++)
		for (q=0; q<3; q++){
		  inside_points[ip][q] = (real *) malloc( n_points*sizeof(real) );
		  inside_normal[ip][q] = (real *) malloc( n_points*sizeof(real) );
		}
/* get appropriate grid lines for test points */
	      j = range3(1,2);
	      k = (range3(1,3) + range3(2,3))/2;
	      arc = 0.0;
	      for (i=range3(1,1)+1; i<=range3(2,1); i++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i-1,j,k)) + 
			    sqr(y3(i,j,k)-y3(i-1,j,k)) + 
			    sqr(z3(i,j,k)-z3(i-1,j,k)));
	      darc = 0.0; 
	      for (i_const=range3(1,1)+1; darc < 0.005*length && i_const<=range3(2,1); 
		   i_const++)
		darc += sqrt(sqr(x3(i_const,j,k)-x3(i_const-1,j,k)) + 
			    sqr(y3(i_const,j,k)-y3(i_const-1,j,k)) + 
			    sqr(z3(i_const,j,k)-z3(i_const-1,j,k)));
		    
	      i = range3(1,1);
	      k = (range3(1,3) + range3(2,3))/2;
	      arc = 0.0;
	      for (j=range3(1,2)+1; j<=range3(2,2); j++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i,j-1,k)) + 
			    sqr(y3(i,j,k)-y3(i,j-1,k)) + 
			    sqr(z3(i,j,k)-z3(i,j-1,k)));
	      darc = 0.0; 
	      for (j_const=range3(1,2)+1; darc < 0.005*length && j_const<=range3(2,2); 
		   j_const++)
		darc += sqrt(sqr(x3(i,j_const,k)-x3(i,j_const-1,k)) + 
			    sqr(y3(i,j_const,k)-y3(i,j_const-1,k)) + 
			    sqr(z3(i,j_const,k)-z3(i,j_const-1,k)));
		    
	      printf("Gridlines for test points: i_const = %i, j_const = %i\n", 
		     i_const, j_const);
	      for (k = 0; k < n_points; k++)
		{
		  x_edge[k] = x3(range3(1,1), range3(1,2), k+range3(1,3));
		  y_edge[k] = y3(range3(1,1), range3(1,2), k+range3(1,3));
		  z_edge[k] = z3(range3(1,1), range3(1,2), k+range3(1,3));

		  inside_points[0][0][k] = x3(i_const, range3(1,2), k+range3(1,3));
		  inside_points[0][1][k] = y3(i_const, range3(1,2), k+range3(1,3));
		  inside_points[0][2][k] = z3(i_const, range3(1,2), k+range3(1,3));

		  if (grid_face_normal(n_vec, vgrid, i_const, range3(1,2), 
				       k+range3(1,3), 2, 1)){
		    inside_normal[0][0][k] = n_vec[0];
		    inside_normal[0][1][k] = n_vec[1];
		    inside_normal[0][2][k] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 1);
		    inside_normal[0][0][k] = 0;
		    inside_normal[0][1][k] = 0;
		    inside_normal[0][2][k] = 0;
		  }

		  inside_points[1][0][k] = x3(range3(1,1), j_const, k+range3(1,3));
		  inside_points[1][1][k] = y3(range3(1,1), j_const, k+range3(1,3));
		  inside_points[1][2][k] = z3(range3(1,1), j_const, k+range3(1,3));

		  if (grid_face_normal(n_vec, vgrid, range3(1,1), j_const, 
				       k+range3(1,3), 1, 1)){
		    inside_normal[1][0][k] = n_vec[0];
		    inside_normal[1][1][k] = n_vec[1];
		    inside_normal[1][2][k] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 1);
		    inside_normal[1][0][k] = 0;
		    inside_normal[1][1][k] = 0;
		    inside_normal[1][2][k] = 0;
		  }
		} /* end for */
	      /* fix the end points */
/* 	      for (ip=0; ip<2; ip++) */
/* 		for (q=0; q<3; q++) */
/* 		  { */
/* 		    inside_points[ip][q][0] = inside_points[ip][q][1]; */
/* 		    inside_points[ip][q][n_points-1] = inside_points[ip][q][n_points-2]; */
/* 		  } */

	      hole_curve = new_hole_curve(x_edge, y_edge, z_edge, n_points, 
					  inside_points, inside_normal, edge_curve(1));
	      insert_hole_curve( overg, hole_curve );
	    } /* end if get_existing_hs */
	  else
	    printf("compute_hole_curves: Error: could not find a hole_surface with "
		   "surface label = %i\n", surf3(1,1));
	} /* end if edge_curve(1) */

      if (edge_curve(3))	/* r1=0, r2=1 */
	{
	  printf("Constructing a hole-cutting curve from edge 3 of component `%s'\n", 
		 vgrid->name);
	  /* get an existing overlapping surface grid with the current surface label */
	  if ((overg = get_existing_hs( over3d, surf3(1,1) )))
	    {
	      n_points = range3(2,3) - range3(1,3) + 1;
	      x_edge = (real *) malloc( n_points*sizeof(real) );
	      y_edge = (real *) malloc( n_points*sizeof(real) );
	      z_edge = (real *) malloc( n_points*sizeof(real) );
	      for (ip=0; ip<2; ip++)
		for (q=0; q<3; q++){
		  inside_points[ip][q] = (real *) malloc( n_points*sizeof(real) );
		  inside_normal[ip][q] = (real *) malloc( n_points*sizeof(real) );
		}
/* get appropriate grid lines for test points */
	      j = range3(2,2);
	      k = (range3(1,3) + range3(2,3))/2;
	      arc = 0.0;
	      for (i=range3(1,1)+1; i<=range3(2,1); i++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i-1,j,k)) + 
			    sqr(y3(i,j,k)-y3(i-1,j,k)) + 
			    sqr(z3(i,j,k)-z3(i-1,j,k)));
	      darc = 0.0; 
	      for (i_const=range3(1,1)+1; darc < 0.005*length && i_const<=range3(2,1); 
		   i_const++)
		darc += sqrt(sqr(x3(i_const,j,k)-x3(i_const-1,j,k)) + 
			    sqr(y3(i_const,j,k)-y3(i_const-1,j,k)) + 
			    sqr(z3(i_const,j,k)-z3(i_const-1,j,k)));
		    
	      i = range3(1,1);
	      k = (range3(1,3) + range3(2,3))/2;
	      arc = 0.0;
	      for (j=range3(1,2)+1; j<=range3(2,2); j++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i,j-1,k)) + 
			    sqr(y3(i,j,k)-y3(i,j-1,k)) + 
			    sqr(z3(i,j,k)-z3(i,j-1,k)));
	      darc = 0.0; 
	      for (j_const=range3(2,2)-1; darc < 0.005*length && j_const>=range3(1,2); 
		   j_const--)
		darc += sqrt(sqr(x3(i,j_const,k)-x3(i,j_const+1,k)) + 
			    sqr(y3(i,j_const,k)-y3(i,j_const+1,k)) + 
			    sqr(z3(i,j_const,k)-z3(i,j_const+1,k)));
		    
	      printf("Gridlines for test points: i_const = %i, j_const = %i\n", 
		     i_const, j_const);
	      for (k = 0; k < n_points; k++)
		{
		  x_edge[k] = x3(range3(1,1), range3(2,2), k+range3(1,3));
		  y_edge[k] = y3(range3(1,1), range3(2,2), k+range3(1,3));
		  z_edge[k] = z3(range3(1,1), range3(2,2), k+range3(1,3));

		  inside_points[0][0][k] = x3(i_const, range3(2,2), k+range3(1,3));
		  inside_points[0][1][k] = y3(i_const, range3(2,2), k+range3(1,3));
		  inside_points[0][2][k] = z3(i_const, range3(2,2), k+range3(1,3));

		  if (grid_face_normal(n_vec, vgrid, i_const, range3(2,2), 
				       k+range3(1,3), 2, 1)){
		    inside_normal[0][0][k] = n_vec[0];
		    inside_normal[0][1][k] = n_vec[1];
		    inside_normal[0][2][k] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 3);
		    inside_normal[0][0][k] = 0;
		    inside_normal[0][1][k] = 0;
		    inside_normal[0][2][k] = 0;
		  }

		  inside_points[1][0][k] = x3(range3(1,1), j_const, k+range3(1,3));
		  inside_points[1][1][k] = y3(range3(1,1), j_const, k+range3(1,3));
		  inside_points[1][2][k] = z3(range3(1,1), j_const, k+range3(1,3));

		  if (grid_face_normal(n_vec, vgrid, range3(1,1), j_const, 
				       k+range3(1,3), 1, 1)){
		    inside_normal[1][0][k] = n_vec[0];
		    inside_normal[1][1][k] = n_vec[1];
		    inside_normal[1][2][k] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 3);
		    inside_normal[1][0][k] = 0;
		    inside_normal[1][1][k] = 0;
		    inside_normal[1][2][k] = 0;
		  }
		} /* end for */
	      /* fix the end points */
/* 	      for (ip=0; ip<2; ip++) */
/* 		for (q=0; q<3; q++) */
/* 		  { */
/* 		    inside_points[ip][q][0] = inside_points[ip][q][1]; */
/* 		    inside_points[ip][q][n_points-1] = inside_points[ip][q][n_points-2]; */
/* 		  } */

	      hole_curve = new_hole_curve(x_edge, y_edge, z_edge, n_points, 
					  inside_points, inside_normal, edge_curve(3));
	      insert_hole_curve( overg, hole_curve );
	    } /* end if get_existing_hs */
	  else
	    printf("compute_hole_curves: Error: could not find a hole_surface with "
		   "surface label = %i\n", surf3(1,1));
	} /* end if edge_curve(3) */

      if (edge_curve(5))	/* r1=0, r3=0 */
	{
	  printf("Constructing a hole-cutting curve from edge 5 of component `%s'\n", 
		 vgrid->name);
	  /* get an existing overlapping surface grid with the current surface label */
	  if ((overg = get_existing_hs( over3d, surf3(1,1) )))
	    {
	      n_points = range3(2,2) - range3(1,2) + 1;
	      x_edge = (real *) malloc( n_points*sizeof(real) );
	      y_edge = (real *) malloc( n_points*sizeof(real) );
	      z_edge = (real *) malloc( n_points*sizeof(real) );
	      for (ip=0; ip<2; ip++)
		for (q=0; q<3; q++){
		  inside_points[ip][q] = (real *) malloc( n_points*sizeof(real) );
		  inside_normal[ip][q] = (real *) malloc( n_points*sizeof(real) );
		}
/* get appropriate grid lines for test points */
	      j = (range3(1,2) + range3(2,2))/2;
	      k = range3(1,3);
	      arc = 0.0;
	      for (i=range3(1,1)+1; i<=range3(2,1); i++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i-1,j,k)) + 
			    sqr(y3(i,j,k)-y3(i-1,j,k)) + 
			    sqr(z3(i,j,k)-z3(i-1,j,k)));
	      darc = 0.0; 
	      for (i_const=range3(1,1)+1; darc < 0.005*length && i_const<=range3(2,1); 
		   i_const++)
		darc += sqrt(sqr(x3(i_const,j,k)-x3(i_const-1,j,k)) + 
			    sqr(y3(i_const,j,k)-y3(i_const-1,j,k)) + 
			    sqr(z3(i_const,j,k)-z3(i_const-1,j,k)));
		    
	      arc = 0.0;
	      j = (range3(1,2) + range3(2,2))/2;
	      i = range3(1,1);
	      for (k=range3(1,3)+1; k<=range3(2,3); k++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i,j,k-1)) + 
			    sqr(y3(i,j,k)-y3(i,j,k-1)) + 
			    sqr(z3(i,j,k)-z3(i,j,k-1)));
	      darc = 0.0; 
	      for (k_const=range3(1,3)+1; darc < 0.005*length && k_const<=range3(2,3); 
		   k_const++)
		darc += sqrt(sqr(x3(i,j,k_const)-x3(i,j,k_const-1)) + 
			    sqr(y3(i,j,k_const)-y3(i,j,k_const-1)) + 
			    sqr(z3(i,j,k_const)-z3(i,j,k_const-1)));

	      printf("Gridlines for test points: i_const = %i, k_const = %i\n", 
		     i_const, k_const);
	      for (j = 0; j < n_points; j++)
		{
		  x_edge[j] = x3(range3(1,1), j+range3(1,2), range3(1,3));
		  y_edge[j] = y3(range3(1,1), j+range3(1,2), range3(1,3));
		  z_edge[j] = z3(range3(1,1), j+range3(1,2), range3(1,3));

		  inside_points[0][0][j] = x3(i_const, j+range3(1,2), range3(1,3));
		  inside_points[0][1][j] = y3(i_const, j+range3(1,2), range3(1,3));
		  inside_points[0][2][j] = z3(i_const, j+range3(1,2), range3(1,3));

		  if (grid_face_normal(n_vec, vgrid, i_const, j+range3(1,2), 
				       range3(1,3), 3, 1)){
		    inside_normal[0][0][j] = n_vec[0];
		    inside_normal[0][1][j] = n_vec[1];
		    inside_normal[0][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 5);
		    inside_normal[0][0][j] = 0;
		    inside_normal[0][1][j] = 0;
		    inside_normal[0][2][j] = 0;
		  }

		  inside_points[1][0][j] = x3(range3(1,1), j+range3(1,2), k_const);
		  inside_points[1][1][j] = y3(range3(1,1), j+range3(1,2), k_const);
		  inside_points[1][2][j] = z3(range3(1,1), j+range3(1,2), k_const);

		  if (grid_face_normal(n_vec, vgrid, range3(1,1), j+range3(1,2), 
				       k_const, 1, 1)){
		    inside_normal[1][0][j] = n_vec[0];
		    inside_normal[1][1][j] = n_vec[1];
		    inside_normal[1][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 5);
		    inside_normal[1][0][j] = 0;
		    inside_normal[1][1][j] = 0;
		    inside_normal[1][2][j] = 0;
		  }
		} /* end for */
	      /* fix the end points */
/* 	      for (ip=0; ip<2; ip++) */
/* 		for (q=0; q<3; q++) */
/* 		  { */
/* 		    inside_points[ip][q][0] = inside_points[ip][q][1]; */
/* 		    inside_points[ip][q][n_points-1] = inside_points[ip][q][n_points-2]; */
/* 		  } */
	      hole_curve = new_hole_curve(x_edge, y_edge, z_edge, n_points, 
					  inside_points, inside_normal, edge_curve(5));
	      insert_hole_curve( overg, hole_curve );
	    } /* end if get_existing_hs */
	  else
	    printf("compute_hole_curves: Error: could not find a hole_surface with "
		   "surface label = %i\n", surf3(1,1));
	} /* end if edge_curve(5) */

      if (edge_curve(6))	/* r1=1, r3=0 */
	{
	  printf("Constructing a hole-cutting curve from edge 6 of component `%s'\n", 
		 vgrid->name);
	  /* get an existing overlapping surface grid with the current surface label */
	  if ((overg = get_existing_hs( over3d, surf3(2,1) )))
	    {
	      n_points = range3(2,2) - range3(1,2) + 1;
	      x_edge = (real *) malloc( n_points*sizeof(real) );
	      y_edge = (real *) malloc( n_points*sizeof(real) );
	      z_edge = (real *) malloc( n_points*sizeof(real) );
	      for (ip=0; ip<2; ip++)
		for (q=0; q<3; q++){
		  inside_points[ip][q] = (real *) malloc( n_points*sizeof(real) );
		  inside_normal[ip][q] = (real *) malloc( n_points*sizeof(real) );
		}
/* get appropriate grid lines for test points */
	      j = (range3(1,2) + range3(2,2))/2;
	      k = range3(1,3);
	      arc = 0.0;
	      for (i=range3(1,1)+1; i<=range3(2,1); i++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i-1,j,k)) + 
			    sqr(y3(i,j,k)-y3(i-1,j,k)) + 
			    sqr(z3(i,j,k)-z3(i-1,j,k)));
	      darc = 0.0; 
	      for (i_const=range3(2,1)-1; darc < 0.005*length && i_const>=range3(1,1); 
		   i_const--)
		darc += sqrt(sqr(x3(i_const,j,k)-x3(i_const+1,j,k)) + 
			    sqr(y3(i_const,j,k)-y3(i_const+1,j,k)) + 
			    sqr(z3(i_const,j,k)-z3(i_const+1,j,k)));
		    
	      arc = 0.0;
	      j = (range3(1,2) + range3(2,2))/2;
	      i = range3(2,1);
	      for (k=range3(1,3)+1; k<=range3(2,3); k++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i,j,k-1)) + 
			    sqr(y3(i,j,k)-y3(i,j,k-1)) + 
			    sqr(z3(i,j,k)-z3(i,j,k-1)));
	      darc = 0.0; 
	      for (k_const=range3(1,3)+1; darc < 0.005*length && k_const<=range3(2,3); 
		   k_const++)
		darc += sqrt(sqr(x3(i,j,k_const)-x3(i,j,k_const-1)) + 
			    sqr(y3(i,j,k_const)-y3(i,j,k_const-1)) + 
			    sqr(z3(i,j,k_const)-z3(i,j,k_const-1)));

	      printf("Gridlines for test points: i_const = %i, k_const = %i\n", 
		     i_const, k_const);
	      for (j = 0; j < n_points; j++)
		{
		  x_edge[j] = x3(range3(2,1), j+range3(1,2), range3(1,3));
		  y_edge[j] = y3(range3(2,1), j+range3(1,2), range3(1,3));
		  z_edge[j] = z3(range3(2,1), j+range3(1,2), range3(1,3));

		  inside_points[0][0][j] = x3(i_const, j+range3(1,2), range3(1,3));
		  inside_points[0][1][j] = y3(i_const, j+range3(1,2), range3(1,3));
		  inside_points[0][2][j] = z3(i_const, j+range3(1,2), range3(1,3));

		  if (grid_face_normal(n_vec, vgrid, i_const, j+range3(1,2), 
				       range3(1,3), 3, 1)){
		    inside_normal[0][0][j] = n_vec[0];
		    inside_normal[0][1][j] = n_vec[1];
		    inside_normal[0][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 6);
		    inside_normal[0][0][j] = 0;
		    inside_normal[0][1][j] = 0;
		    inside_normal[0][2][j] = 0;
		  }

		  inside_points[1][0][j] = x3(range3(2,1), j+range3(1,2), k_const);
		  inside_points[1][1][j] = y3(range3(2,1), j+range3(1,2), k_const);
		  inside_points[1][2][j] = z3(range3(2,1), j+range3(1,2), k_const);

		  if (grid_face_normal(n_vec, vgrid, range3(2,1), j+range3(1,2), 
				       k_const, 1, 1)){
		    inside_normal[1][0][j] = n_vec[0];
		    inside_normal[1][1][j] = n_vec[1];
		    inside_normal[1][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 6);
		    inside_normal[1][0][j] = 0;
		    inside_normal[1][1][j] = 0;
		    inside_normal[1][2][j] = 0;
		  }
		} /* end for */
	      /* fix the end points */
/* 	      for (ip=0; ip<2; ip++) */
/* 		for (q=0; q<3; q++) */
/* 		  { */
/* 		    inside_points[ip][q][0] = inside_points[ip][q][1]; */
/* 		    inside_points[ip][q][n_points-1] = inside_points[ip][q][n_points-2]; */
/* 		  } */
	      hole_curve = new_hole_curve(x_edge, y_edge, z_edge, n_points, 
					  inside_points, inside_normal, edge_curve(6));
	      insert_hole_curve( overg, hole_curve );
	    } /* end if get_existing_hs */
	  else
	    printf("compute_hole_curves: Error: could not find a hole_surface with "
		   "surface label = %i\n", surf3(2,1));
	} /* end if edge_curve(6) */

      if (edge_curve(7))	/* r1=0, r3=1 */
	{
	  printf("Constructing a hole-cutting curve from edge 7 of component `%s'\n", 
		 vgrid->name);
	  /* get an existing overlapping surface grid with the current surface label */
	  if ((overg = get_existing_hs( over3d, surf3(1,1) )))
	    {
	      n_points = range3(2,2) - range3(1,2) + 1;
	      x_edge = (real *) malloc( n_points*sizeof(real) );
	      y_edge = (real *) malloc( n_points*sizeof(real) );
	      z_edge = (real *) malloc( n_points*sizeof(real) );
	      for (ip=0; ip<2; ip++)
		for (q=0; q<3; q++){
		  inside_points[ip][q] = (real *) malloc( n_points*sizeof(real) );
		  inside_normal[ip][q] = (real *) malloc( n_points*sizeof(real) );
		}
/* get appropriate grid lines for test points */
	      j = (range3(1,2) + range3(2,2))/2;
	      k = range3(2,3);
	      arc = 0.0;
	      for (i=range3(1,1)+1; i<=range3(2,1); i++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i-1,j,k)) + 
			    sqr(y3(i,j,k)-y3(i-1,j,k)) + 
			    sqr(z3(i,j,k)-z3(i-1,j,k)));
	      darc = 0.0; 
	      for (i_const=range3(1,1)+1; darc < 0.005*length && i_const<=range3(2,1); 
		   i_const++)
		darc += sqrt(sqr(x3(i_const,j,k)-x3(i_const-1,j,k)) + 
			    sqr(y3(i_const,j,k)-y3(i_const-1,j,k)) + 
			    sqr(z3(i_const,j,k)-z3(i_const-1,j,k)));
		    
	      arc = 0.0;
	      j = (range3(1,2) + range3(2,2))/2;
	      i = range3(1,1);
	      for (k=range3(1,3)+1; k<=range3(2,3); k++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i,j,k-1)) + 
			    sqr(y3(i,j,k)-y3(i,j,k-1)) + 
			    sqr(z3(i,j,k)-z3(i,j,k-1)));
	      darc = 0.0; 
	      for (k_const=range3(2,3)-1; darc < 0.005*length && k_const>=range3(1,3); 
		   k_const--)
		darc += sqrt(sqr(x3(i,j,k_const)-x3(i,j,k_const+1)) + 
			     sqr(y3(i,j,k_const)-y3(i,j,k_const+1)) + 
			     sqr(z3(i,j,k_const)-z3(i,j,k_const+1)));

	      printf("Gridlines for test points: i_const = %i, k_const = %i\n", 
		     i_const, k_const);
	      for (j = 0; j < n_points; j++)
		{
		  x_edge[j] = x3(range3(1,1), j+range3(1,2), range3(2,3));
		  y_edge[j] = y3(range3(1,1), j+range3(1,2), range3(2,3));
		  z_edge[j] = z3(range3(1,1), j+range3(1,2), range3(2,3));

		  inside_points[0][0][j] = x3(i_const, j+range3(1,2), range3(2,3));
		  inside_points[0][1][j] = y3(i_const, j+range3(1,2), range3(2,3));
		  inside_points[0][2][j] = z3(i_const, j+range3(1,2), range3(2,3));

		  if (grid_face_normal(n_vec, vgrid, i_const, j+range3(1,2), 
				       range3(2,3), 3, 1)){
		    inside_normal[0][0][j] = n_vec[0];
		    inside_normal[0][1][j] = n_vec[1];
		    inside_normal[0][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 7);
		    inside_normal[0][0][j] = 0;
		    inside_normal[0][1][j] = 0;
		    inside_normal[0][2][j] = 0;
		  }

		  inside_points[1][0][j] = x3(range3(1,1), j+range3(1,2), k_const);
		  inside_points[1][1][j] = y3(range3(1,1), j+range3(1,2), k_const);
		  inside_points[1][2][j] = z3(range3(1,1), j+range3(1,2), k_const);

		  if (grid_face_normal(n_vec, vgrid, range3(1,1), j+range3(1,2), 
				       k_const, 1, 1)){
		    inside_normal[1][0][j] = n_vec[0];
		    inside_normal[1][1][j] = n_vec[1];
		    inside_normal[1][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 7);
		    inside_normal[1][0][j] = 0;
		    inside_normal[1][1][j] = 0;
		    inside_normal[1][2][j] = 0;
		  }
		} /* end for */
	      /* fix the end points */
/* 	      for (ip=0; ip<2; ip++) */
/* 		for (q=0; q<3; q++) */
/* 		  { */
/* 		    inside_points[ip][q][0] = inside_points[ip][q][1]; */
/* 		    inside_points[ip][q][n_points-1] = inside_points[ip][q][n_points-2]; */
/* 		  } */
	      hole_curve = new_hole_curve(x_edge, y_edge, z_edge, n_points, 
					  inside_points, inside_normal, edge_curve(7));
	      insert_hole_curve( overg, hole_curve );
	    } /* end if get_existing_hs */
	  else
	    printf("compute_hole_curves: Error: could not find a hole_surface with "
		   "surface label = %i\n", surf3(1,1));
	} /* end if edge_curve(7) */

      if (edge_curve(9))	/* r2=0, r3=0 */
	{
	  printf("Constructing a hole-cutting curve from edge 9 of component `%s'\n", 
		 vgrid->name);
	  /* get an existing overlapping surface grid with the current surface label */
	  if ((overg = get_existing_hs( over3d, surf3(1,2) )))
	    {
	      n_points = range3(2,1) - range3(1,1) + 1;
	      x_edge = (real *) malloc( n_points*sizeof(real) );
	      y_edge = (real *) malloc( n_points*sizeof(real) );
	      z_edge = (real *) malloc( n_points*sizeof(real) );
	      for (ip=0; ip<2; ip++)
		for (q=0; q<3; q++){
		  inside_points[ip][q] = (real *) malloc( n_points*sizeof(real) );
		  inside_normal[ip][q] = (real *) malloc( n_points*sizeof(real) );
		}
/* get appropriate grid lines for test points */
	      i = (range3(1,1) + range3(2,1))/2;
	      k = range3(1,3);
	      arc = 0.0;
	      for (j=range3(1,2)+1; j<=range3(2,2); j++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i,j-1,k)) + 
			    sqr(y3(i,j,k)-y3(i,j-1,k)) + 
			    sqr(z3(i,j,k)-z3(i,j-1,k)));
	      darc = 0.0; 
	      for (j_const=range3(1,2)+1; darc < 0.005*length && j_const<=range3(2,2); 
		   j_const++)
		darc += sqrt(sqr(x3(i,j_const,k)-x3(i,j_const-1,k)) + 
			    sqr(y3(i,j_const,k)-y3(i,j_const-1,k)) + 
			    sqr(z3(i,j_const,k)-z3(i,j_const-1,k)));
		    
	      arc = 0.0;
	      i = (range3(1,1) + range3(2,1))/2;
	      j = range3(1,2);
	      for (k=range3(1,3)+1; k<=range3(2,3); k++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i,j,k-1)) + 
			    sqr(y3(i,j,k)-y3(i,j,k-1)) + 
			    sqr(z3(i,j,k)-z3(i,j,k-1)));
	      darc = 0.0; 
	      for (k_const=range3(1,3)+1; darc < 0.005*length && k_const<=range3(2,3); 
		   k_const++)
		darc += sqrt(sqr(x3(i,j,k_const)-x3(i,j,k_const-1)) + 
			    sqr(y3(i,j,k_const)-y3(i,j,k_const-1)) + 
			    sqr(z3(i,j,k_const)-z3(i,j,k_const-1)));

	      printf("Gridlines for test points: j_const = %i, k_const = %i\n", 
		     j_const, k_const);

	      for (j = 0; j < n_points; j++)
		{
		  x_edge[j] = x3( j+range3(1,1), range3(1,2), range3(1,3));
		  y_edge[j] = y3( j+range3(1,1), range3(1,2), range3(1,3));
		  z_edge[j] = z3( j+range3(1,1), range3(1,2), range3(1,3));

		  inside_points[0][0][j] = x3( j+range3(1,1), j_const, range3(1,3));
		  inside_points[0][1][j] = y3( j+range3(1,1), j_const, range3(1,3));
		  inside_points[0][2][j] = z3( j+range3(1,1), j_const, range3(1,3));

		  if (grid_face_normal(n_vec, vgrid, j+range3(1,1), j_const, 
				       range3(1,3), 3, 1)){
		    inside_normal[0][0][j] = n_vec[0];
		    inside_normal[0][1][j] = n_vec[1];
		    inside_normal[0][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 9);
		    inside_normal[0][0][j] = 0;
		    inside_normal[0][1][j] = 0;
		    inside_normal[0][2][j] = 0;
		  }

		  inside_points[1][0][j] = x3( j+range3(1,1), range3(1,2), k_const);
		  inside_points[1][1][j] = y3( j+range3(1,1), range3(1,2), k_const);
		  inside_points[1][2][j] = z3( j+range3(1,1), range3(1,2), k_const);

		  if (grid_face_normal(n_vec, vgrid, j+range3(1,1), range3(1,2), 
				       k_const, 2, 1)){
		    inside_normal[1][0][j] = n_vec[0];
		    inside_normal[1][1][j] = n_vec[1];
		    inside_normal[1][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 9);
		    inside_normal[1][0][j] = 0;
		    inside_normal[1][1][j] = 0;
		    inside_normal[1][2][j] = 0;
		  }
		} /* end for */
	      /* fix the end points */
/* 	      for (ip=0; ip<2; ip++) */
/* 		for (q=0; q<3; q++) */
/* 		  { */
/* 		    inside_points[ip][q][0] = inside_points[ip][q][1]; */
/* 		    inside_points[ip][q][n_points-1] = inside_points[ip][q][n_points-2]; */
/* 		  } */
	      hole_curve = new_hole_curve(x_edge, y_edge, z_edge, n_points, 
					  inside_points, inside_normal, edge_curve(9));
	      insert_hole_curve( overg, hole_curve );
	    } /* end if get_existing_hs */
	  else
	    printf("compute_hole_curves: Error: could not find a hole_surface with "
		   "surface label = %i\n", surf3(1,2));
	} /* end if edge_curve(9) */

      if (edge_curve(10))	/* r2=1, r3=0 */
	{
	  printf("Constructing a hole-cutting curve from edge 10 of component `%s'\n", 
		 vgrid->name);
	  /* get an existing overlapping surface grid with the current surface label */
	  if ((overg = get_existing_hs( over3d, surf3(2,2) )))
	    {
	      n_points = range3(2,1) - range3(1,1) + 1;
	      x_edge = (real *) malloc( n_points*sizeof(real) );
	      y_edge = (real *) malloc( n_points*sizeof(real) );
	      z_edge = (real *) malloc( n_points*sizeof(real) );
	      for (ip=0; ip<2; ip++)
		for (q=0; q<3; q++){
		  inside_points[ip][q] = (real *) malloc( n_points*sizeof(real) );
		  inside_normal[ip][q] = (real *) malloc( n_points*sizeof(real) );
		}
/* get appropriate grid lines for test points */
	      i = (range3(1,1) + range3(2,1))/2;
	      k = range3(1,3);
	      arc = 0.0;
	      for (j=range3(1,2)+1; j<=range3(2,2); j++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i,j-1,k)) + 
			    sqr(y3(i,j,k)-y3(i,j-1,k)) + 
			    sqr(z3(i,j,k)-z3(i,j-1,k)));
	      darc = 0.0; 
	      for (j_const=range3(2,2)-1; darc < 0.005*length && j_const>=range3(1,2); 
		   j_const--)
		darc += sqrt(sqr(x3(i,j_const,k)-x3(i,j_const+1,k)) + 
			    sqr(y3(i,j_const,k)-y3(i,j_const+1,k)) + 
			    sqr(z3(i,j_const,k)-z3(i,j_const+1,k)));
		    
	      arc = 0.0;
	      i = (range3(1,1) + range3(2,1))/2;
	      j = range3(2,2);
	      for (k=range3(1,3)+1; k<=range3(2,3); k++)
		arc += sqrt(sqr(x3(i,j,k)-x3(i,j,k-1)) + 
			    sqr(y3(i,j,k)-y3(i,j,k-1)) + 
			    sqr(z3(i,j,k)-z3(i,j,k-1)));
	      darc = 0.0; 
	      for (k_const=range3(1,3)+1; darc < 0.005*length && k_const<=range3(2,3); 
		   k_const++)
		darc += sqrt(sqr(x3(i,j,k_const)-x3(i,j,k_const-1)) + 
			    sqr(y3(i,j,k_const)-y3(i,j,k_const-1)) + 
			    sqr(z3(i,j,k_const)-z3(i,j,k_const-1)));

	      printf("Gridlines for test points: j_const = %i, k_const = %i\n", 
		     j_const, k_const);


	      for (j = 0; j < n_points; j++)
		{
		  x_edge[j] = x3( j+range3(1,1), range3(2,2), range3(1,3));
		  y_edge[j] = y3( j+range3(1,1), range3(2,2), range3(1,3));
		  z_edge[j] = z3( j+range3(1,1), range3(2,2), range3(1,3));

		  inside_points[0][0][j] = x3( j+range3(1,1), j_const, range3(1,3));
		  inside_points[0][1][j] = y3( j+range3(1,1), j_const, range3(1,3));
		  inside_points[0][2][j] = z3( j+range3(1,1), j_const, range3(1,3));

		  if (grid_face_normal(n_vec, vgrid, j+range3(1,1), j_const, 
				       range3(1,3), 3, 1)){
		    inside_normal[0][0][j] = n_vec[0];
		    inside_normal[0][1][j] = n_vec[1];
		    inside_normal[0][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 10);
		    inside_normal[0][0][j] = 0;
		    inside_normal[0][1][j] = 0;
		    inside_normal[0][2][j] = 0;
		  }

		  inside_points[1][0][j] = x3( j+range3(1,1), range3(2,2), k_const);
		  inside_points[1][1][j] = y3( j+range3(1,1), range3(2,2), k_const);
		  inside_points[1][2][j] = z3( j+range3(1,1), range3(2,2), k_const);

		  if (grid_face_normal(n_vec, vgrid, j+range3(1,1), range3(2,2), 
				       k_const, 2, 1)){
		    inside_normal[1][0][j] = n_vec[0];
		    inside_normal[1][1][j] = n_vec[1];
		    inside_normal[1][2][j] = n_vec[2];
		  }
		  else{
		    printf("Warning: unable to compute normal in compute_hole_curves, "
			   "vgrid = `%s', edge_curve(%i)\n", vgrid->name, 10);
		    inside_normal[1][0][j] = 0;
		    inside_normal[1][1][j] = 0;
		    inside_normal[1][2][j] = 0;
		  }
		} /* end for */
	      /* fix the end points */
/* 	      for (ip=0; ip<2; ip++) */
/* 		for (q=0; q<3; q++) */
/* 		  { */
/* 		    inside_points[ip][q][0] = inside_points[ip][q][1]; */
/* 		    inside_points[ip][q][n_points-1] = inside_points[ip][q][n_points-2]; */
/* 		  } */
	      hole_curve = new_hole_curve(x_edge, y_edge, z_edge, n_points, 
					  inside_points, inside_normal, edge_curve(10));
	      insert_hole_curve( overg, hole_curve );
	    } /* end if get_existing_hs */
	  else
	    printf("compute_hole_curves: Error: could not find a hole_surface with "
		   "surface label = %i\n", surf3(2,2));
	} /* end if edge_curve(10) */

    } /* end for all vgrids... */
  for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;
       overg_link = overg_link->next)
    {
      overg = overg_link->data;
      merge_hole_curves(overg);
    } /* end for all hole surfaces */
  
}

static hole_surface *
get_or_make_hs( overlapping_3d_grid *over3d, int surface_label )
{
  linked_list_member *overg_link;
  hole_surface *overg;
  char overg_name[80];

  for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;
       overg_link = overg_link->next){
				       overg = overg_link->data;
				       if (overg->surface_label == surface_label)
					 return overg;
				     }

  /* if we got this far, we need to make a new overlapping surface grid */
  sprintf( overg_name, "surface-label-%i", surface_label );
  overg = new_hole_surface( overg_name, surface_label );

  /* insert the new hole_surface into the corresponding list */
  new_link(over3d->surface_grid_list)->data = overg;

  return overg;
}
     
static hole_surface *
get_existing_hs( overlapping_3d_grid *over3d, int surface_label )
{
  linked_list_member *overg_link;
  hole_surface *overg;

  for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;
       overg_link = overg_link->next)
    {
      overg = overg_link->data;
      if (overg->surface_label == surface_label)
	return overg;
    }

  /* if we got this far there is no hole_surface with that surface label */
  return NULL;
}

discrete_surface *
new_sgrid_from_vgrid( component_vgrid *vgrid, int side, int dir, real global_tolerance )
{
  discrete_surface *sgrid_ptr;
  real_array_2d *x_surf_, *y_surf_, *z_surf_;
  char sgrid_name[120];

  /* get the coordinates of the face */
  vgrid_face_coordinates(vgrid, side, dir, &x_surf_, &y_surf_, &z_surf_);

  /* make a discrete surface from the coordinates */
  sprintf(sgrid_name,"%s-side=%i-dir=%i", vgrid->name, side, dir);
  sgrid_ptr = new_discrete_surface(sgrid_name, x_surf_, y_surf_, z_surf_, 
				   surf3(side,dir), vgrid);

  sgrid_ptr->side = side;
  sgrid_ptr->dir = dir;

  /* the correspondence between the directions in the mother grid and the component */
  /* surface grid is as follows */
  /* dir  r_dim   s_dim  */
  /*  1   r2_dim  r3_dim */
  /*  2   r3_dim  r1_dim */
  /*  3   r1_dim  r2_dim */

  /* set the curve-label in the discrete surface */
  if (dir == 1)
    {
      if (surf3(1,2) > 0)
	curve(1,1) = surf3(side,dir) + surf3(1,2);
      if (surf3(2,2) > 0)
	curve(2,1) = surf3(side,dir) + surf3(2,2);
      if (surf3(1,3) > 0)
	curve(1,2) = surf3(side,dir) + surf3(1,3);
      if (surf3(2,3) > 0)
	curve(2,2) = surf3(side,dir) + surf3(2,3);
      if (side == 1)
	{
	  hole_curve(1,1) = edge_curve(1); /* r1=0, r2=0 */
	  hole_curve(2,1) = edge_curve(3); /* r1=0, r2=1 */
	  hole_curve(1,2) = edge_curve(5); /* r1=0, r3=0 */
	  hole_curve(2,2) = edge_curve(7); /* r1=0, r3=1 */
	}
      else			/* side == 2 */
	{
	  hole_curve(1,1) = edge_curve(2); /* r1=1, r2=0 */
	  hole_curve(2,1) = edge_curve(4); /* r1=1, r2=1 */
	  hole_curve(1,2) = edge_curve(6); /* r1=1, r3=0 */
	  hole_curve(2,2) = edge_curve(8); /* r1=1, r3=1 */
	}
    }
  else if (dir == 2)
    {
      if (surf3(1,1) > 0)
	curve(1,2) = surf3(side,dir) + surf3(1,1);
      if (surf3(2,1) > 0)
	curve(2,2) = surf3(side,dir) + surf3(2,1);
      if (surf3(1,3) > 0)
	curve(1,1) = surf3(side,dir) + surf3(1,3);
      if (surf3(2,3) > 0)
	curve(2,1) = surf3(side,dir) + surf3(2,3);
      if (side == 1)
	{
	  hole_curve(1,1) = edge_curve(9); /* r2=0, r3=0 */
	  hole_curve(2,1) = edge_curve(11); /* r2=0, r3=1 */
	  hole_curve(1,2) = edge_curve(1); /* r2=0, r1=0 */
	  hole_curve(2,2) = edge_curve(2); /* r2=0, r1=1 */
	}
      else			/* side == 2 */
	{
	  hole_curve(1,1) = edge_curve(10); /* r2=1, r3=0 */
	  hole_curve(2,1) = edge_curve(12); /* r2=1, r3=1 */
	  hole_curve(1,2) = edge_curve(3); /* r2=1, r1=0 */
	  hole_curve(2,2) = edge_curve(4); /* r2=1, r1=1 */
	}
    }
  else
    {				/* dir == 3 */
      if (surf3(1,1) > 0)
	curve(1,1) = surf3(side,dir) + surf3(1,1);
      if (surf3(2,1) > 0)
	curve(2,1) = surf3(side,dir) + surf3(2,1);
      if (surf3(1,2) > 0)
	curve(1,2) = surf3(side,dir) + surf3(1,2);
      if (surf3(2,2) > 0)
	curve(2,2) = surf3(side,dir) + surf3(2,2);
      if (side == 1)
	{
	  hole_curve(1,1) = edge_curve(5); /* r3=0, r1=0 */
	  hole_curve(2,1) = edge_curve(6); /* r3=0, r1=1 */
	  hole_curve(1,2) = edge_curve(9); /* r3=0, r2=0 */
	  hole_curve(2,2) = edge_curve(10); /* r3=0, r2=1 */
	}
      else			/* side == 2 */
	{
	  hole_curve(1,1) = edge_curve(7); /* r3=1, r1=0 */
	  hole_curve(2,1) = edge_curve(8); /* r3=1, r1=1 */
	  hole_curve(1,2) = edge_curve(11); /* r3=1, r2=0 */
	  hole_curve(2,2) = edge_curve(12); /* r3=1, r2=1 */
	}
    }
  
  return sgrid_ptr;
}


static void 
vgrid_face_coordinates(component_vgrid *vgrid, int side, int dir, 
		       real_array_2d **x_surf_, 
		       real_array_2d **y_surf_, 
		       real_array_2d **z_surf_)
{
  int i, j, constant, i_start, j_start, n1, n2;

#define x_surf(i, j) compute_index_2d((*x_surf_), i, j)
#define y_surf(i, j) compute_index_2d((*y_surf_), i, j)
#define z_surf(i, j) compute_index_2d((*z_surf_), i, j)

  /* the correspondence between the directions in the mother grid and the discrete */
  /* surface are as follows */
  /* dir  r_dim   s_dim  */
  /*  1   r2_dim  r3_dim */
  /*  2   r3_dim  r1_dim */
  /*  3   r1_dim  r2_dim */

  /* assign x, y, z, flag */
  if (dir == 1)
    {
      /* get the index of the grid surface in the constant direction of the 3-d vgrid */
      constant = range3(side, dir);

      /*       if (vgrid->r2_period) */
      /* 	i_start = range3(1,2)-2; */
      /*       else */
      i_start =  range3(1,2)-1;
      n1 = range3(2,2) - i_start;

      /*       if (vgrid->r3_period) */
      /* 	j_start = range3(1,3)-2; */
      /*       else */
      j_start =  range3(1,3)-1;
      n2 = range3(2,3) - j_start;
      
      /* define the coordinate arrays */
      *x_surf_ = create_real_array_2d(n1, n2);
      *y_surf_ = create_real_array_2d(n1, n2);
      *z_surf_ = create_real_array_2d(n1, n2);

      for (i = 1; i <= n1; i++)
	for (j = 1; j <= n2; j++)
	  {
	    x_surf(i,j)  = x3(constant,i+i_start,j+j_start);
	    y_surf(i,j)  = y3(constant,i+i_start,j+j_start);
	    z_surf(i,j)  = z3(constant,i+i_start,j+j_start);
	  }
    }
  else if (dir == 2)
    {
      /* get the index of the grid surface in the constant direction of the 3-d vgrid */
      constant = range3(side, dir);

      /*       if (vgrid->r3_period) */
      /* 	i_start = range3(1,3)-2; */
      /*       else */
      i_start =  range3(1,3)-1;
      n1 = range3(2,3) - i_start;

      /*       if (vgrid->r1_period) */
      /* 	j_start = range3(1,1)-2; */
      /*       else */
      j_start =  range3(1,1)-1;
      n2 = range3(2,1) - j_start;
      
      /* define the coordinate arrays */
      *x_surf_ = create_real_array_2d(n1, n2);
      *y_surf_ = create_real_array_2d(n1, n2);
      *z_surf_ = create_real_array_2d(n1, n2);

      for (i = 1; i <= n1; i++)
	for (j = 1; j <= n2; j++)
	  {
	    x_surf(i,j)  = x3(j+j_start,constant,i+i_start);
	    y_surf(i,j)  = y3(j+j_start,constant,i+i_start);
	    z_surf(i,j)  = z3(j+j_start,constant,i+i_start);
	  }

    }
  else
    {				/* dir == 3 */
      /* get the index of the grid surface in the constant direction of the 3-d vgrid */
      constant = range3(side, dir);

      /*       if (vgrid->r1_period) */
      /* 	i_start = range3(1,1)-2; */
      /*       else */
      i_start =  range3(1,1)-1;
      n1 = range3(2,1) - i_start;

      /*       if (vgrid->r2_period) */
      /* 	j_start = range3(1,2)-2; */
      /*       else */
      j_start =  range3(1,2)-1;
      n2 = range3(2,2) - j_start;
      
      /* define the coordinate arrays */
      *x_surf_ = create_real_array_2d(n1, n2);
      *y_surf_ = create_real_array_2d(n1, n2);
      *z_surf_ = create_real_array_2d(n1, n2);

      for (i = 1; i <= n1; i++)
	for (j = 1; j <= n2; j++)
	  {
	    x_surf(i,j)  = x3(i+i_start,j+j_start,constant);
	    y_surf(i,j)  = y3(i+i_start,j+j_start,constant);
	    z_surf(i,j)  = z3(i+i_start,j+j_start,constant);
	  }

    }

} /* end vgrid_face_coordinates */

void 
init_3d_arrays(volume_mapping *volume, component_vgrid *vgrid, 
	       overlapping_3d_grid *over3d)
{
  int i, j, k;
  real jac;

  /* number of extra points at periodic boundaries */
  over3d->extra_period = int_max((over3d->disc_width - 1)/2, 
				 (over3d->interp_width - 1)/2);

  if (vgrid->r1_period)
    {
      range3(1,1) = over3d->extra_period + 1; 
      range3(2,1) = volume->r1_points + over3d->extra_period;
      vgrid->r1_dim = volume->r1_points + 2*over3d->extra_period - 1;
    }
  else{
	range3(1,1) = over3d->extra + 1; 
	range3(2,1) = volume->r1_points + over3d->extra;
	vgrid->r1_dim = volume->r1_points + 2*over3d->extra;
      }

  if (vgrid->r2_period){
			 range3(1,2) = over3d->extra_period + 1; 
			 range3(2,2) = volume->r2_points + over3d->extra_period;
			 vgrid->r2_dim = volume->r2_points + 2*over3d->extra_period - 1;
		       }
  else{
	range3(1,2) = over3d->extra + 1; 
	range3(2,2) = volume->r2_points + over3d->extra;
	vgrid->r2_dim = volume->r2_points + 2*over3d->extra;
      }

  if (vgrid->r3_period){
			 range3(1,3) = over3d->extra_period + 1; 
			 range3(2,3) = volume->r3_points + over3d->extra_period;
			 vgrid->r3_dim = volume->r3_points + 2*over3d->extra_period - 1;
		       }
  else{
	range3(1,3) = over3d->extra + 1; 
	range3(2,3) = volume->r3_points + over3d->extra;
	vgrid->r3_dim = volume->r3_points + 2*over3d->extra;
      }

  vgrid->r1_step = 1.0/((real) (volume->r1_points-1));
  vgrid->r2_step = 1.0/((real) (volume->r2_points-1));
  vgrid->r3_step = 1.0/((real) (volume->r3_points-1));

/* create and assign the x, y, z arrays and create the flag array */
  copy_from_volume_mapping( vgrid, volume );

/* check the sign of the jacobian */
  i = range3(1,1); j = range3(1,2); k = range3(1,3);
/* the following statement says */
/* jac = xr*(ys*zt-yt*zs) - xs*(yr*zt-yt*zr) + xt*(yr*zs-ys*zr) */
  jac = 
    (x3(i+1,j,k)-x3(i,j,k)) * 
      ((y3(i,j+1,k)-y3(i,j,k))*(z3(i,j,k+1)-z3(i,j,k)) - 
       (y3(i,j,k+1)-y3(i,j,k))*(z3(i,j+1,k)-z3(i,j,k))) - 
	 (x3(i,j+1,k)-x3(i,j,k)) * 
	   ((y3(i+1,j,k)-y3(i,j,k))*(z3(i,j,k+1)-z3(i,j,k)) - 
	    (y3(i,j,k+1)-y3(i,j,k))*(z3(i+1,j,k)-z3(i,j,k))) +
	      (x3(i,j,k+1)-x3(i,j,k)) * 
		((y3(i+1,j,k)-y3(i,j,k))*(z3(i,j+1,k)-z3(i,j,k)) - 
		 (y3(i,j+1,k)-y3(i,j,k))*(z3(i+1,j,k)-z3(i,j,k)));
  if (jac>0.0)
    vgrid->orientation = 1;
  else
    vgrid->orientation = -1;
  /* report the orientation */
  if (over3d->verbose >= 2)
    printf("The component `%s' is %s.\n", vgrid->name, 
	   (vgrid->orientation == 1)? "right-handed":"left-handed");

  /* ensure perfect periodicity if r1_period || r2_period || r3_period */
  if (vgrid->r1_period)
    for (k=1; k<=vgrid->r3_dim; k++)
      {
	for (j=1; j<=vgrid->r2_dim; j++)
	  {
	    x3(range3(2,1),j,k) = x3(range3(1,1),j,k);
	    y3(range3(2,1),j,k) = y3(range3(1,1),j,k);
	    z3(range3(2,1),j,k) = z3(range3(1,1),j,k);
	  }
      }
  
  if (vgrid->r2_period) 
    for (k=1; k<=vgrid->r3_dim; k++)
      { 
	for (i=1; i<= vgrid->r1_dim; i++)
	  { 
	    x3(i,range3(2,2),k) = x3(i,range3(1,2),k); 
	    y3(i,range3(2,2),k) = y3(i,range3(1,2),k); 
	    z3(i,range3(2,2),k) = z3(i,range3(1,2),k); 
	  } 
      } 
  
  if (vgrid->r3_period)
    for (j=1; j<=vgrid->r2_dim; j++)
      {
	for (i=1; i<= vgrid->r1_dim; i++)
	  {
	    x3(i,j,range3(2,3)) = x3(i,j,range3(1,3));
	    y3(i,j,range3(2,3)) = y3(i,j,range3(1,3));
	    z3(i,j,range3(2,3)) = z3(i,j,range3(1,3));
	  }
      }
  
} /* end define grid */
     
static int 
overlap_3d(input_output *io_ptr, overlapping_3d_grid *over3d)
{
  int n_grid_points, ok, n_members, icom, side, dir, i, j, k, dum
    , present_curve, present_surface;
  struct tms time0, time1, time3, time4, time5, time6, time7, time8, time9;
  hole_surface *overg;
  component_vgrid *vgrid, *donor_vgrid;
  linked_list_member *vgrid_link, *overg_link;
  char * break_com[3];
  const int LEVEL=1;
  interp_3d_point *interp, *new_interp;
  real x_dist, y_dist, z_dist;
  /*  FILE *err_msg;*/

  break_com[0] = "proceed";
  break_com[1] = "break";
  break_com[2] = "plot-mode";
       
  /* get starting time */
  times(&time0);
  
  /* to begin with, the grid is regarded as invalid */
  over3d->valid = 0;
       
  /* Evaluate the mappings at the grid points and compute the octree */
  prepare_3d_overlap( over3d );
  
  /* Setup the hole-cutting surfaces corresponding to each surface label. */
  prepare_hole_surfaces(io_ptr, over3d);
       
  if (over3d->show_boundary && ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, 
					      ogl_length_scale(over3d->bb))) 
    { 
      for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;
	   overg_link = overg_link->next)
	{
	  overg = overg_link->data;
	  draw_hole_surface( overg, 64 | 512 | 2048 ); /* turn on the triangles */
	} /* end for */
      ogl_end_plot(); 

      /* break or proceed? */    
      while ((icom=get_command(io_ptr, "You are looking at the physical boundary>", 
			       break_com, 2, LEVEL, NULL, NULL)) != 0)
	{
	  switch (icom)
	    {
	    case 1: /* exit to calling routine */
	      return ERROR;
	      break;
	    case 2:  /* modify the plot mode for one surface label */
	      if (overg_link = get_hole_surface( io_ptr, over3d->surface_grid_list))
		{
		  overg = overg_link->data;
		  overg->plot_mode = hole_plot_mode(io_ptr, overg, 64 | 512 | 2048 );
		}
	      break;
	    default:
	      ;
	    } /* end switch */
	} /* end while */
    } /* end if show_boundary */

  initialize_flag(over3d);
  times(&time1);
       
  cut_volume_holes(io_ptr, over3d);
  times(&time3); 
  
  if (over3d->show_volume_holes && ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, 
						  ogl_length_scale(over3d->bb)))
    {    
      draw_overlapping_3d_grid( over3d, 4 | 4096 );    
      ogl_end_plot();    

      /* break or proceed? */    
      while ((icom=get_command(io_ptr, "You are looking at the volume holes>", 
			       break_com, 2, LEVEL, NULL, NULL)) != 0)
	{
	  switch (icom)
	    {
	    case 1: /* exit to calling routine */
	      return ERROR;
	      break;
	    case 2:  /* modify the plot mode */
	      over3d->plot_mode = overlapping_3d_plot_mode(io_ptr, over3d, 4 | 4096 );
	      break;
	    default:
	      ;
	    } /* end switch */
	} /* end while */
    } /* end if show_volume_holes */

  if (over3d->verbose >= 2)
    printf("\n");

/* test the invert_grid routine */
#ifdef TEST_INVERT_GRID
  vgrid = NULL;
  donor_vgrid = NULL;
  printf("Select the component grid to be tested\n");     
  if (vgrid_link = get_vgrid_ptr(io_ptr, over3d->grid_list))   
    {   
      vgrid = vgrid_link->data;     
    }
  printf("Select donor grid\n");     
  if (vgrid_link = get_vgrid_ptr(io_ptr, over3d->grid_list))   
    {   
      donor_vgrid = vgrid_link->data;     
    }
  if (vgrid && donor_vgrid)
    {
      printf("Grid bounds:\n%i <= i <= %i\n%i <= i <= %i\n%i <= i <= %i\n", 
	     range3(1,1), range3(2,1), range3(1,2), range3(2,2), 
	     range3(1,3), range3(2,3));
      i=j=k=1;
      do
	{
	  i = get_int(io_ptr, "i (quit with -1): ", i, 2);
	  if (i>0)
	    {
	      j = get_int(io_ptr, "j: ", j, 2);
	      k = get_int(io_ptr, "k: ", k, 2);
	      present_surface = get_surface_label(i, j, k, vgrid, &dum, &dum);
	      present_curve   = get_curve_label(i, j, k, vgrid);
	      printf("Surface label = %i\n", present_surface);
	      if (new_interp = invert_grid(x3(i,j,k), y3(i,j,k), z3(i,j,k), donor_vgrid, 
					   present_surface, present_curve, NULL))
		{
		  printf("The point is inside `%s'\n", donor_vgrid->name);
		  free(new_interp);
		} 
	      else
		{
		  printf("The point is outside `%s'\n", donor_vgrid->name);
		}/* end if */
	    }
	} while (i>0); /* end while */
    } /* end if */
#endif

  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL; 
       vgrid_link = vgrid_link->next){
    vgrid = vgrid_link->data;
    initial_classification(vgrid_link, over3d);
  } /* end for vgrid ... */
  times(&time4); 
  
#ifdef STOP_AFTER_CLASSIFICATION
/* plot the overlapping 3-d grid */
  printf("initial_classification completed.\n");
  if (ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, ogl_length_scale(over3d->bb))){    
    draw_overlapping_3d_grid( over3d, over3d->plot_mode );    
    ogl_end_plot();    
  }
/* get the chance to look at the grid differently! */
  over3d->plot_mode = overlapping_3d_plot_mode(io_ptr, over3d, over3d->plot_mode);
#endif

  if (over3d->interp_type == 'e'){
/* do the iteration */
     for (vgrid_link = over3d->grid_list->last; vgrid_link != NULL; 
  	 vgrid_link = vgrid_link->prev){ 
       explicit(vgrid_link->data, over3d); 
     }  
/* remove deactivated interpolation points */
    for (vgrid_link = over3d->grid_list->last; vgrid_link != NULL; 
	 vgrid_link = vgrid_link->prev){
      remove_deactive_3d_interp( vgrid_link->data );
    } 
  }
  times(&time5); 

/* delete any existing bad points */
  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL; 
       vgrid_link = vgrid_link->next){
    vgrid = vgrid_link->data;
    vgrid->last_bad_point = delete_bad_3d_list( vgrid->last_bad_point );
  }

/* check the consistency */
  ok = TRUE;
  if (over3d->verbose >= 2)
    printf("\nChecking the consistency of all grid points...\n");
  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL; 
       vgrid_link = vgrid_link->next){
    if (!consistency_check(vgrid_link->data, over3d)) ok = FALSE;
  }
  if (ok){
    if (over3d->verbose >= 2)
      printf("All discretization and interpolation points OK.\n");
  }
  else
    return ERROR;

  times(&time6); 

  printf("\n");  
  if (over3d->interp_type == 'i'){
    mark_needed_interp(over3d);
    for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL;
	 vgrid_link = vgrid_link->next)
      trim_interpolation_points(vgrid_link->data, over3d);
  }
  else{
    for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL;
	 vgrid_link = vgrid_link->next)
      trim_interpolation_points(vgrid_link->data, over3d);
  }

  times(&time7); 

  if (over3d->correct_mismatch){
    printf("\n");
/*    err_msg = fopen("mismatch.dat","w"); */

    for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL; 
	 vgrid_link = vgrid_link->next){
      vgrid = vgrid_link->data;
      printf("Mismatch correction for grid `%s'\n",vgrid->name);
/*      fprintf(err_msg, "\nMismatch correction for grid `%s'\n",vgrid->name);*/

/* first check all interpolation points on physical boundries */
      for (interp = vgrid->last_interp; interp != NULL;
	   interp = interp->prev){

	if (!interp->active) continue;

/* must also take the curve_label into account! */

/* direction 3: update the whole line */
	if (surf3(1,3) && interp->k_point == range3(1,3) &&
	    boundary_diff(surf3(1,3), vgrid, interp, &x_dist, &y_dist, &z_dist))
	  {
	    recompute_interpolation_location(/*err_msg,*/ interp, vgrid, over3d, 3,
					     over3d->interp_width - 1, x_dist, y_dist, z_dist);
	  }
	else if (surf3(2,3) && interp->k_point == range3(2,3) &&
	    boundary_diff(surf3(2,3), vgrid, interp, &x_dist, &y_dist, &z_dist))
	  {
	    recompute_interpolation_location(/*err_msg,*/ interp, vgrid, over3d, 3, 
					     over3d->interp_width - 1, x_dist, y_dist, z_dist);
	  }

/* direction 1: just update one point */
  	if (surf3(1,1) && interp->i_point == range3(1,1))
	  recompute_one_interp(surf3(1,1), vgrid, interp);
  	else if (surf3(2,1) && interp->i_point == range3(2,1))
	  recompute_one_interp(surf3(2,1), vgrid, interp);

/* direction 2: just update one point */
  	if (surf3(1,2) && interp->j_point == range3(1,2))
	  recompute_one_interp(surf3(1,2), vgrid, interp);
  	else if (surf3(2,2) && interp->j_point == range3(2,2))
	  recompute_one_interp(surf3(2,2), vgrid, interp);

/* direction 3: just update one point (to remove roundoff errors from */
/* recompute_interpolation_location above */
  	if (surf3(1,3) && interp->k_point == range3(1,3))
	  recompute_one_interp(surf3(1,3), vgrid, interp);
  	else if (surf3(2,3) && interp->k_point == range3(2,3))
	  recompute_one_interp(surf3(2,3), vgrid, interp);

      } /* end for all interpolation points */


/* check the mismatch */
      for (interp = vgrid->last_interp; interp != NULL;
	   interp = interp->prev)
	if (interp->active && 
	    (!interp->vgrid_loc->r1_period &&
	     (interp->r_loc < -NEWTON_EPS || interp->r_loc > 1+NEWTON_EPS)) ||
	    (!interp->vgrid_loc->r2_period &&
	     (interp->s_loc < -NEWTON_EPS || interp->s_loc > 1+NEWTON_EPS)) ||
	    (!interp->vgrid_loc->r3_period &&
	     (interp->t_loc < -NEWTON_EPS || interp->t_loc > 1+NEWTON_EPS)))
	  {
	    printf("Warning: Interpolation point (%i, %i, %i) in grid %s extrapolates "
		   "from donor grid %s\n interpolation location (%e, %e, %e)\n", 
		   interp->i_point, interp->j_point, interp->k_point, vgrid->name, 
		   interp->vgrid_loc->name, interp->r_loc, interp->s_loc, interp->t_loc);
/* 	    fprintf(err_msg,  */
/* 		    "Warning: Interpolation point (%i, %i, %i) in grid %s extrapolates " */
/* 		   "from donor grid %s\n interpolation location (%e, %e, %e)\n",  */
/* 		   interp->i_point, interp->j_point, interp->k_point, vgrid->name,  */
/* 		   interp->vgrid_loc->name, interp->r_loc, interp->s_loc, interp->t_loc); */
	  }
    } /* end for all component grids */
/*    fclose(err_msg);*/

  } /* end if correct_mismatch */

  times(&time8); 

  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL;
       vgrid_link = vgrid_link->next)
    cleanup_and_finish(vgrid_link->data, over3d);

  times(&time9); 

  if (over3d->verbose >= 1){
    printf("\nTiming information for the 3-d overlap algorithm:\n");
    printf("Initialization of data structures: %.1f seconds\n", 
	   (time1.tms_utime-time0.tms_utime)/60.0);
    printf("Hole-cutting:                      %.1f seconds\n", 
	   (time3.tms_utime-time1.tms_utime)/60.0);
    printf("Main classification:               %.1f seconds\n", 
	   (time4.tms_utime-time3.tms_utime)/60.0);
    if (over3d->interp_type == 'e')
      printf("Ensuring explicit interpolation:   %.1f seconds\n", 
	   (time5.tms_utime-time4.tms_utime)/60.0);
    printf("Consitency check:                  %.1f seconds\n", 
	   (time6.tms_utime-time5.tms_utime)/60.0);
    printf("Trimming interpolation points:     %.1f seconds\n", 
	   (time7.tms_utime-time6.tms_utime)/60.0);
    if (over3d->correct_mismatch)
      printf("Mismatch correction:               %.1f seconds\n", 
	     (time8.tms_utime-time7.tms_utime)/60.0);
    printf("Finishing touches:                 %.1f seconds\n", 
	   (time9.tms_utime-time8.tms_utime)/60.0);
    printf("\n");
    printf("Total CPU time:                    %.1f seconds\n", 
	   (time9.tms_utime-time0.tms_utime)/60.0);
  }

/* add up the number of grid points */
  n_grid_points = 0;
  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL; 
       vgrid_link = vgrid_link->next){
    vgrid = vgrid_link->data;
/*    if (vgrid->dummy_background) continue;*/
    n_grid_points += (range3(2,1) - range3(1,1) + 1) * 
      (range3(2,2) - range3(1,2) + 1) * (range3(2,3) - range3(1,3) + 1);
  }
  if (over3d->verbose >= 1){
    n_members = /*(over3d->overlapping_sub_grid)? over3d->grid_list->n_members-1 : */
      over3d->grid_list->n_members;
    printf("\n");
    printf("This overlapping grid has %i components and a total of %i grid points,\n"
	   "yielding a CPU time per grid point and component of %e seconds.\n", 
	   n_members, n_grid_points, 
	   (time9.tms_utime-time0.tms_utime)/(60.0 * n_members * n_grid_points) );
  }
  
/* mark the overlapping grid to be valid */
  over3d->valid = 1;
  
  return over3d->valid;
}

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif
#define surf_d(i,j) compute_index_2d( donor_vgrid->surf_ptr, i, j)
#define range_d(i,j) compute_index_2d( donor_vgrid->range_ptr, i, j)
#define x_d(i,j,k) compute_index_3d( donor_vgrid->x_ptr, i, j, k)
#define y_d(i,j,k) compute_index_3d( donor_vgrid->y_ptr, i, j, k)
#define z_d(i,j,k) compute_index_3d( donor_vgrid->z_ptr, i, j, k)

static int
boundary_diff(int present_surf, component_vgrid *vgrid, interp_3d_point *interp,
	      real *x_dist, real *y_dist, real *z_dist){
  real x_dist_tmp, y_dist_tmp, z_dist_tmp;
  int found_closest_bndry = FALSE;
  grid_point gp;
  component_vgrid *donor_vgrid;

  if (present_surf <= 0) return FALSE;

  donor_vgrid = interp->vgrid_loc;

  *x_dist = *y_dist = *z_dist = 1.e30;
/* low r in donor grid */
  if (surf_d(1,1) == present_surf){
    gp.r = 0.0;
    gp.s = interp->s_loc;
    gp.t = interp->t_loc;
    forward_vgrid_mapping( &gp, donor_vgrid );
    *x_dist = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    *y_dist = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    *z_dist = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    found_closest_bndry = TRUE;
  }
/* high r in donor grid */
  if (surf_d(2,1) == present_surf){
    gp.r = 1.0;
    gp.s = interp->s_loc;
    gp.t = interp->t_loc;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(*x_dist) + sqr(*y_dist) + sqr(*z_dist))){
      found_closest_bndry = TRUE;
      *x_dist = x_dist_tmp;
      *y_dist = y_dist_tmp;
      *z_dist = z_dist_tmp;
    }	    
  }
/* low s in donor grid */
  if (surf_d(1,2) == present_surf){
    gp.r = interp->r_loc;
    gp.s = 0.0;
    gp.t = interp->t_loc;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(*x_dist) + sqr(*y_dist) + sqr(*z_dist))){
      found_closest_bndry = TRUE;
      *x_dist = x_dist_tmp;
      *y_dist = y_dist_tmp;
      *z_dist = z_dist_tmp;
    }	    
  }
/* high s in donor grid */
  if (surf_d(2,2) == present_surf){
    gp.r = interp->r_loc;
    gp.s = 1.0;
    gp.t = interp->t_loc;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(*x_dist) + sqr(*y_dist) + sqr(*z_dist))){
      found_closest_bndry = TRUE;
      *x_dist = x_dist_tmp;
      *y_dist = y_dist_tmp;
      *z_dist = z_dist_tmp;
    }	    
  }
/* low t in donor grid */
  if (surf_d(1,3) == present_surf){
    gp.r = interp->r_loc;
    gp.s = interp->s_loc;
    gp.t = 0.0;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(*x_dist) + sqr(*y_dist) + sqr(*z_dist))){
      found_closest_bndry = TRUE;
      *x_dist = x_dist_tmp;
      *y_dist = y_dist_tmp;
      *z_dist = z_dist_tmp;
    }	    
  }
/* high t in donor grid */
  if (surf_d(2,3) == present_surf){
    gp.r = interp->r_loc;
    gp.s = interp->s_loc;
    gp.t = 1.0;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(*x_dist) + sqr(*y_dist) + sqr(*z_dist))){
      found_closest_bndry = TRUE;
      *x_dist = x_dist_tmp;
      *y_dist = y_dist_tmp;
      *z_dist = z_dist_tmp;
    }	    
  }
  return found_closest_bndry;
}

static void
recompute_one_interp(int present_surf, component_vgrid *vgrid, interp_3d_point *interp){
  real x_dist, y_dist, z_dist, x_dist_tmp, y_dist_tmp, z_dist_tmp;
  int found_closest_bndry = FALSE, side=0, dir=0;
  grid_point gp;
  component_vgrid *donor_vgrid;

  if (present_surf <= 0) return;

  donor_vgrid = interp->vgrid_loc;

  x_dist = y_dist = z_dist = 1.e30;
/* low r in donor grid */
  if (surf_d(1,1) == present_surf){
    gp.r = 0.0;
    gp.s = interp->s_loc;
    gp.t = interp->t_loc;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    side = 1; dir = 1;
    found_closest_bndry = TRUE;
  }
/* high r in donor grid */
  if (surf_d(2,1) == present_surf){
    gp.r = 1.0;
    gp.s = interp->s_loc;
    gp.t = interp->t_loc;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(x_dist) + sqr(y_dist) + sqr(z_dist))){
      side = 2; dir = 1;
      found_closest_bndry = TRUE;
      x_dist = x_dist_tmp;
      y_dist = y_dist_tmp;
      z_dist = z_dist_tmp;
    }	    
  }
/* low s in donor grid */
  if (surf_d(1,2) == present_surf){
    gp.r = interp->r_loc;
    gp.s = 0.0;
    gp.t = interp->t_loc;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(x_dist) + sqr(y_dist) + sqr(z_dist))){
      side = 1; dir = 2;
      found_closest_bndry = TRUE;
      x_dist = x_dist_tmp;
      y_dist = y_dist_tmp;
      z_dist = z_dist_tmp;
    }	    
  }
/* high s in donor grid */
  if (surf_d(2,2) == present_surf){
    gp.r = interp->r_loc;
    gp.s = 1.0;
    gp.t = interp->t_loc;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(x_dist) + sqr(y_dist) + sqr(z_dist))){
      side = 2; dir = 2;
      found_closest_bndry = TRUE;
      x_dist = x_dist_tmp;
      y_dist = y_dist_tmp;
      z_dist = z_dist_tmp;
    }	    
  }
/* low t in donor grid */
  if (surf_d(1,3) == present_surf){
    gp.r = interp->r_loc;
    gp.s = interp->s_loc;
    gp.t = 0.0;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(x_dist) + sqr(y_dist) + sqr(z_dist))){
      side = 1; dir = 3;
      found_closest_bndry = TRUE;
      x_dist = x_dist_tmp;
      y_dist = y_dist_tmp;
      z_dist = z_dist_tmp;
    }	    
  }
/* high t in donor grid */
  if (surf_d(2,3) == present_surf){
    gp.r = interp->r_loc;
    gp.s = interp->s_loc;
    gp.t = 1.0;
    forward_vgrid_mapping( &gp, donor_vgrid );
    x_dist_tmp = gp.x - x3(interp->i_point, interp->j_point, interp->k_point);
    y_dist_tmp = gp.y - y3(interp->i_point, interp->j_point, interp->k_point);
    z_dist_tmp = gp.z - z3(interp->i_point, interp->j_point, interp->k_point);
    if ((sqr(x_dist_tmp) + sqr(y_dist_tmp) + sqr(z_dist_tmp)) < 
	(sqr(x_dist) + sqr(y_dist) + sqr(z_dist))){
      side = 2; dir = 3;
      found_closest_bndry = TRUE;
      x_dist = x_dist_tmp;
      y_dist = y_dist_tmp;
      z_dist = z_dist_tmp;
    }	    
  }

  if (found_closest_bndry)
    {
      if (dir == 1 && side == 1)
	interp->r_loc = 0.0;
      else if (dir == 1 && side == 2)
	interp->r_loc = 1.0;
      else if (dir == 2 && side == 1)
	interp->s_loc = 0.0;
      else if (dir == 2 && side == 2)
	interp->s_loc = 1.0;
      else if (dir == 3 && side == 1)
	interp->t_loc = 0.0;
      else if (dir == 3 && side == 2)
	interp->t_loc = 1.0;
    }
}

static void
recompute_interpolation_location(/* FILE *err_msg,*/ interp_3d_point *interp, 
				 component_vgrid *vgrid, 
				 overlapping_3d_grid *over3d, int dir,
				 int iw1, real x_dist, real y_dist, real z_dist){
  grid_point gp;
  interp_3d_point *interp2, *new_interp;
  int present_surface, present_curve, dum, proj_surf;
  component_vgrid *donor_vgrid;

/* tmp */
/*   if (real_max(fabs(z_dist), real_max(fabs(x_dist), fabs(y_dist))) > vgrid->max_s_dist) */
/*       fprintf(err_msg, "Warning: Correcting large mismatch at interpolation point (%i,%i,%i)\n" */
/* 	     "The mismatch vector is (%e,%e,%e)\n",  */
/* 	 interp->i_point, interp->j_point, interp->k_point, x_dist, y_dist, z_dist); */

  proj_surf = 0;
  if (dir == 3 && interp->i_point == range3(1,1) && surf3(1,1) > 0)
    proj_surf = surf3(1,1);
  else if (dir == 3 && interp->i_point == range3(2,1) && surf3(2,1) > 0)
    proj_surf = surf3(2,1);

  if (dir == 3 && interp->j_point == range3(1,2) && surf3(1,2) > 0)
    proj_surf = surf3(1,2);
  else if (dir == 3 && interp->j_point == range3(2,2) && surf3(2,2) > 0)
    proj_surf = surf3(2,2);


/* change all points along the same grid line */
  for (interp2 = vgrid->last_interp; interp2 != NULL;
       interp2 = interp2->prev){

/* don't do the same point twice! */
    if (interp2->corrected[dir-1]) continue;

/* check if interp2 is on the same grid line as interp */
    if ((dir == 1 &&
	interp2->k_point == interp->k_point && interp2->j_point == interp->j_point) ||
	(dir == 2 &&
	interp2->i_point == interp->i_point && interp2->k_point == interp->k_point) ||
	(dir == 3 &&
	interp2->i_point == interp->i_point && interp2->j_point == interp->j_point)){
/* recompute the location*/
      gp.x = x3(interp2->i_point, interp2->j_point, interp2->k_point) + x_dist; 
      gp.y = y3(interp2->i_point, interp2->j_point, interp2->k_point) + y_dist; 
      gp.z = z3(interp2->i_point, interp2->j_point, interp2->k_point) + z_dist; 

/* set the surface label if the present point is on a real boundary */
      present_surface = get_surface_label(interp2->i_point, interp2->j_point, 
					  interp2->k_point, vgrid, &dum, &dum);
      present_curve   = get_curve_label(interp2->i_point, interp2->j_point,
					interp2->k_point, vgrid);

/* set the correction flag */
      interp2->corrected[dir-1] = TRUE;

      if (new_interp = invert_grid(gp.x, gp.y, gp.z, interp2->vgrid_loc, 
				   present_surface, present_curve, interp2)){
/* save the indices of the interpolation point */
	new_interp->i_point = interp2->i_point;
	new_interp->j_point = interp2->j_point;
	new_interp->k_point = interp2->k_point;
/* check if (i,j,k) is a valid interpolation location*/
	if (!good_interp_loc(new_interp, interp2->vgrid_loc, over3d))
	  {
	    printf("Warning, the corrected interpolation location for ip (%i,%i,%i) is close\n" 
		   "to the boundary of the donor grid `%s`. r_loc=(%f, %f, %f)\n",  
		   interp2->i_point, interp2->j_point, interp2->k_point, vgrid->name,
		   new_interp->r_loc, new_interp->s_loc, new_interp->t_loc); 
	  }
/* update the coefficients */
	    interp2->r_loc = new_interp->r_loc;
	    interp2->s_loc = new_interp->s_loc;
	    interp2->t_loc = new_interp->t_loc;

	    interp2->i_loc = new_interp->i_loc;
	    interp2->j_loc = new_interp->j_loc;
	    interp2->k_loc = new_interp->k_loc;
	    free(new_interp);

/* check if the location needs to be projected onto a grid surface */
	    if (proj_surf>0)
	      {
		donor_vgrid = interp2->vgrid_loc;

		if (interp2->i_loc == range_d(1,1) && surf_d(1,1) == proj_surf)
		  interp2->r_loc = 0.0;
		else if (interp2->i_loc+iw1 == range_d(2,1) && surf_d(2,1) == proj_surf)
		  interp2->r_loc = 1.0;

		if (interp2->j_loc == range_d(1,2) && surf_d(1,2) == proj_surf)
		  interp2->s_loc = 0.0;
		else if (interp2->j_loc+iw1 == range_d(2,2) && surf_d(2,2) == proj_surf)
		  interp2->s_loc = 1.0;

		if (interp2->k_loc == range_d(1,3) && surf_d(1,3) == proj_surf)
		  interp2->t_loc = 0.0;
		else if (interp2->k_loc+iw1 == range_d(2,3) && surf_d(2,3) == proj_surf)
		  interp2->t_loc = 1.0;

	      }

/* set the correction flag */
	    interp2->corrected[dir-1] = TRUE;

/*           } */
/*  	else  */
/*  	  printf("Warning, the interpolation location for ip (%i,%i,%i) was not "  */
/*  		 "updated because the new location was bad\n",   */
/*  		 interp2->i_point, interp2->j_point, interp2->k_point);  */
      }
      else 
 	printf("Warning, the interpolation location for ip (%i,%i,%i) was not " 
 	       "updated because the new coordinate was outside of the donor " 
 	       "grid %s.\n", interp2->i_point, interp2->j_point, interp2->k_point,  
 		interp2->vgrid_loc->name); 
    } /* end if along the same line */
  } /* end for all interpolation points */
}


static int 
initialize_flag(overlapping_3d_grid *over3d){
  component_vgrid *vgrid;
  linked_list_member *vgrid_link;
  int i, j, k, i_min, i_max, j_min, j_max, k_min, k_max;

  if (over3d->verbose >= 2)
    printf("\n");
  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL; 
       vgrid_link = vgrid_link->next){
    vgrid = vgrid_link->data;

/* make sure this is not a dummy grid */
/*    if (vgrid->dummy_background) continue;*/

    if (over3d->verbose >= 2)
      printf("Initializing component volume grid `%s'\n", vgrid->name);

/* initialize the flag */
    if (vgrid->r1_period){
      i_min = range3(1,1) - over3d->extra_period;
      i_max = range3(2,1) + over3d->extra_period - 1;
    }
    else{
      i_min = range3(1,1);
      i_max = range3(2,1);
    }

    if (vgrid->r2_period){
      j_min = range3(1,2) - over3d->extra_period;
      j_max = range3(2,2) + over3d->extra_period - 1;
    }
    else{
      j_min = range3(1,2);
      j_max = range3(2,2);
    }

    if (vgrid->r3_period){
      k_min = range3(1,3) - over3d->extra_period;
      k_max = range3(2,3) + over3d->extra_period - 1;
    }
    else{
      k_min = range3(1,3);
      k_max = range3(2,3);
    }

/* first initialize the whole array including ghostpoints to zero */
    for (i = 1; i<= vgrid->r1_dim; i++)
      for (j = 1; j<= vgrid->r2_dim; j++)
	for (k = 1; k<= vgrid->r3_dim; k++)
	  flag3(i,j,k) = 0;
  
/* now initialize the real points to discretization points */
/* = the priority of the component grid */
    for (i = i_min; i <= i_max; i++)
      for (j = j_min; j <= j_max; j++)
	for (k = k_min; k <= k_max; k++)
	  flag3(i,j,k) = vgrid->priority;

  } /* end for all components */

  /* done */
  return 1;
}

int 
boundary_cell( int i0, int j0, int k0, component_vgrid *vgrid ){

  if ( ( (i0 == range3(1,1) || i0+1 == range3(1,1)) && surf3(1,1) != 0) ||
       ( (i0 == range3(2,1) || i0+1 == range3(2,1)) && surf3(2,1) != 0) || 
       ( (j0 == range3(1,2) || j0+1 == range3(1,2)) && surf3(1,2) != 0) || 
       ( (j0 == range3(2,2) || j0+1 == range3(2,2)) && surf3(2,2) != 0) || 
       ( (k0 == range3(1,3) || k0+1 == range3(1,3)) && surf3(1,3) != 0) || 
       ( (k0 == range3(2,3) || k0+1 == range3(2,3)) && surf3(2,3) != 0) )
    return 1;

  return 0;
}


static int 
cut_volume_holes(input_output *io_ptr, overlapping_3d_grid *over3d)
{
  int i, j, k, i0, ih, jh, kh, go_right, go_up, go_far;
  hole_surface_link *overg_link;
  hole_surface *overg;
  discrete_surface_link *sgrid_link;
  discrete_surface *sgrid_ptr;
  triangle test;
  component_vgrid_link *vgrid_link;
  component_vgrid *vgrid;
  interp_3d_point *interp[3], prev_sol;
#ifdef TEST_INSIDE_GLOBAL_DOMAIN
  int ip, jp, kp, present_surface, LEVEL=2, dum;
  real radius, n_vec[3];
#endif

/* STEP 1 */
/* Mark the cells that are intersected by the hole-cutting surfaces */

  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL;
       vgrid_link = vgrid_link->next)
    {
      vgrid = vgrid_link->data;
/* do not perform the hole cutting if the no_hole flag is set */
      if (vgrid->no_hole) continue;

      printf("\nMarking component grid `%s'\n", vgrid->name);
      /* check each surface label */
      for (overg_link = over3d->surface_grid_list->first; overg_link!=NULL;
	   overg_link = overg_link->next)
	{
	  overg = overg_link->data;
	  printf("\tMarking cells intersected by `%s'\n", 
		 overg->name);
	  /* check each patch */
	  for (sgrid_link = overg->grid_list->first; sgrid_link != NULL; 
	       sgrid_link = sgrid_link->next)
	    {
	      sgrid_ptr = sgrid_link->data;
/* Since we assume that all mappings are one to one, it is not possible for a */
/* boundary face in a grid to be inside of another part of the same grid */
/* We therefore check where each vertex of the triangle originates from. */
	      if (sgrid_ptr->vgrid_mother != vgrid)
		{
		  printf("\t\tMarking cells intersected by patch `%s'\n", 
			 sgrid_ptr->name);
		  /* check each cell that is alive */
		  for (i=1; i<=sgrid_ptr->r_dim-1; i++)
		    {
		      /* the first point in each column has no initial guess */
		      prev_sol.vgrid_loc = NULL;
		      for (j=1; j<=sgrid_ptr->s_dim-1; j++)
			if (flag(i,j) != 0 && flag(i+1,j) != 0 && 
			    flag(i,j+1) != 0 && flag(i+1,j+1) != 0)
			  {
			    /* printf("\t\tChecking cell (%i,%i)\n", i, j);*/
			    /* first triangle */
			    test.x[0] = x(i,j);   test.y[0] = y(i,j);   
			    test.z[0] = z(i,j);
			    test.x[1] = x(i+1,j); test.y[1] = y(i+1,j); 
			    test.z[1] = z(i+1,j);
			    test.x[2] = x(i,j+1); test.y[2] = y(i,j+1); 
			    test.z[2] = z(i,j+1);
			    /* initialize */
			    for (i0=0; i0<3; i0++) 
			      {
				interp[i0] = NULL; 
				test.sgrid[i0] = sgrid_ptr;
			      }
			    mark_around_tri(test, interp, overg->surface_label, 
					    vgrid, over3d, &prev_sol, 0);
/* get a good initial guess for next triangle */
			    for (i0=0; i0<3; i0++)
			      if (interp[i0])
				{
				  prev_sol.r_loc = interp[i0]->r_loc;
				  prev_sol.s_loc = interp[i0]->s_loc;
				  prev_sol.t_loc = interp[i0]->t_loc;
				  prev_sol.vgrid_loc = interp[i0]->vgrid_loc;
				  break; /* one initial guess is enough */
				}
			    if (interp[0]) free(interp[0]); /* cleanup */
/* second triangle (vertices [1] and [2] are the same as above) */
			    test.x[0] = x(i+1,j+1); 
			    test.y[0] = y(i+1,j+1); 
			    test.z[0] = z(i+1,j+1);
			    interp[0] = NULL; /* initialize */
			    mark_around_tri(test, interp, overg->surface_label, 
					    vgrid, over3d, &prev_sol, 0);
			    for (i0=0; i0<3; i0++) 
			      if (interp[i0]) free(interp[i0]); /* cleanup */
			  } /* end for all cells that are alive */
		    } /* end for i */
		} /* end if sgrid_ptr->vgrid_mother != vgrid */
	    } /* end for all patches */
	} /* end for all surface labels */
    } /* end for all vgrids */

#ifdef VOLUME_HOLES_1
/* plot the coordinate cage and the face with flag = 'w'. */
    if (ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, ogl_length_scale(over3d->bb)))    
       {    
         draw_overlapping_3d_grid(over3d, 256 | 8192);    
         ogl_bb_cage( over3d->bb, OGL_WHITE );    
       }    
    ogl_end_plot();
  printf("After the first stage in the volume hole-cutting algorithm\n");
/* give the opportunity to modify the parameters */
  over3d->plot_mode = overlapping_3d_plot_mode(io_ptr, over3d, 256 | 8192); 
#endif

#ifdef TEST_INSIDE_GLOBAL_DOMAIN
/* test the inside_global_domain routine */
/* plot the hole-surfaces */
/* set the viewpoint to show the grid */
     ogl_standard_view( OVERLAP_W, over3d->bb );     
     if (ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, ogl_length_scale(overg->bb)))     
         {     
           for (overg_link = over3d->surface_grid_list->first; overg_link != NULL;
     	   overg_link = overg_link->next)     
     	{    
     	  overg = overg_link->data;     
     	  draw_hole_surface( overg, 64 | 512 | 2048 );    
     	}    
           ogl_bb_cage( over3d->bb, OGL_WHITE );     
           ogl_end_plot();     
         }     
     printf("Select the component grid to be tested\n");     
     if (vgrid_link = get_vgrid_ptr(io_ptr, over3d->grid_list))   
       {   
         vgrid = vgrid_link->data;     
         radius = 0.01 * real_min(real_min(vgrid->bb->x_max - vgrid->bb->x_min ,      
   					vgrid->bb->y_max - vgrid->bb->y_min),      
   			       vgrid->bb->z_max - vgrid->bb->z_min);     
         printf("Index ranges: %i <= i <= %i, %i <= j <= %i, %i <= k <= %i\n",  
 	       range3(1,1), range3(2,1), range3(1,2), range3(2,2),  
 	       range3(1,3), range3(2,3));     
         while( (ip=get_int(io_ptr, "Enter i-index:", range3(2,1), LEVEL))    
   	    >= range3(1,1) )     
   	{     
   	  jp=get_int(io_ptr, "Enter j-index:", range3(2,2), LEVEL);     
   	  kp=get_int(io_ptr, "Enter k-index:", range3(2,3), LEVEL);     
   	  present_surface = get_surface_label(ip, jp, kp, vgrid, &dum, &dum);
   	  printf("present_surface = %i\n", present_surface);     
   	  ogl_start_plot(OGL_OLD_PLOT, OVERLAP_W, 1.0);     
   	  if (inside_global_domain(x3(ip,jp,kp), y3(ip,jp,kp), z3(ip,jp,kp),      
   				   present_surface, over3d, 1))     
   	    {     
   	      printf("The point (%e, %e, %e) seems to be inside\n",     
   		     x3(ip,jp,kp), y3(ip,jp,kp), z3(ip,jp,kp));    
   	      ogl_marker(x3(ip,jp,kp), y3(ip,jp,kp), z3(ip,jp,kp), radius, OGL_GREEN); 
   	    }    
   	  else   
   	    {    
   	      printf("The point (%e, %e, %e) seems to be outside\n",  
   		     x3(ip,jp,kp), y3(ip,jp,kp), z3(ip,jp,kp));     
   	      ogl_marker(x3(ip,jp,kp), y3(ip,jp,kp), z3(ip,jp,kp), radius, OGL_RED);
   	    }    
   	  glBegin(GL_LINES);   
   	  grid_face_normal(n_vec, vgrid, ip, jp, kp, 3, 1);   
   	  glNormal3dv(n_vec);     
   	  glVertex3d(x3(ip,jp,kp), y3(ip,jp,kp), z3(ip,jp,kp));     
   	  glVertex3d(over3d->bb->x_min - 0.1*(over3d->bb->x_max - over3d->bb->x_min), 
   		     y3(ip,jp,kp), z3(ip,jp,kp)); 
   	  glEnd(); 
   	  ogl_end_plot(); 
   	} 
       } 
#endif

/* STEP 2 */
/* Fill in the holes */

  printf("\n");
  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL;
       vgrid_link = vgrid_link->next)
    {
      vgrid = vgrid_link->data;
      /* Fill in the holes */
      printf("Filling the holes in grid `%s'\n", vgrid->name);
      for (k=range3(1,3); k<=range3(2,3); k++)
	for (j=range3(1,2); j<=range3(2,2); j++)
	  for (i=range3(1,1); i<=range3(2,1); i++)
	    {
	    /* check if it is a hole point */
	    if (flag3(i,j,k) == 'h')
	      {
		/* check the neighbors */
		go_right = (i == range3(1,1) || flag3(i-1,j,k) == 'w' || 
			    flag3(i-1,j,k) == 'h');
		go_up    = (j == range3(1,2) || flag3(i,j-1,k) == 'w' || 
			    flag3(i,j-1,k) == 'h');
		go_far   = (k == range3(1,3) || flag3(i,j,k-1) == 'w' || 
			    flag3(i,j,k-1) == 'h');
		/* horizontal sweep */
		if (go_right)
		  {
		    for (ih = i+1; ih <= range3(2,1) && flag3(ih,j,k) != 'h' && 
			 flag3(ih,j,k) != 'w'; ih++)
		      flag3(ih,j,k) = 0;
		  }
		else
		  {
		    for (ih = i-1; ih >= range3(1,1) && flag3(ih,j,k) != 'h' && 
			 flag3(ih,j,k) != 'w'; ih--)
		      flag3(ih,j,k) = 0;
		  }
		/* vertical sweep */
		if (go_up)
		  {
		    for (jh = j+1; jh <= range3(2,2) && flag3(i,jh,k) != 'h' && 
			 flag3(i,jh,k) != 'w'; jh++)
		      flag3(i,jh,k) = 0;
		  }
		else
		  {
		    for (jh = j-1; jh >= range3(1,2) && flag3(i,jh,k) != 'h' && 
			 flag3(i,jh,k) != 'w'; jh--)
		      flag3(i,jh,k) = 0;
		  }
		/* depth sweep */
 		if (go_far) 
 		  { 
 		    for (kh = k+1; kh <= range3(2,3) && flag3(i,j,kh) != 'h' &&  
 			 flag3(i,j,kh) != 'w'; kh++) 
 		      flag3(i,j,kh) = 0; 
 		  } 
 		else
 		  { 
 		    for (kh = k-1; kh >= range3(1,3) && flag3(i,j,kh) != 'h' &&  
 			 flag3(i,j,kh) != 'w'; kh--) 
 		      flag3(i,j,kh) = 0; 
 		  } 
	      } /* end if 'h' point */
	  } /* end for all grid points */

      /* change 'w' to the  grid priority and 'h' to 0 */
      for (k=range3(1,3); k<=range3(2,3); k++)
	for (j=range3(1,2); j<=range3(2,2); j++)
	  for (i=range3(1,1); i<=range3(2,1); i++)
	    {
	      if (flag3(i,j,k) == 'w')
		flag3(i,j,k) = vgrid->priority;
	      else if (flag3(i,j,k) == 'h')
		flag3(i,j,k) = 0;
	    } /* end for all grid points */
    } /* end for all component grids */

/* STEP 3: Copy surface holes from the faces */
#ifdef COPY_SURFACE_HOLES
  printf("\n");
  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL;
       vgrid_link = vgrid_link->next)
    {
      vgrid = vgrid_link->data;
      printf("Copying surface holes in grid `%s'\n", vgrid->name);
      for (sgrid_link = vgrid->bounding_surface->grid_list->first;
	   sgrid_link != NULL; sgrid_link = sgrid_link->next)
	{
	  sgrid_ptr = sgrid_link->data;
	  if (sgrid_ptr->dir == 1)
	    {
	      if (sgrid_ptr->side == 1)
		i = range3(1,1);
	      else
		i = range3(2,1);
	      for (j=1; j <= sgrid_ptr->r_dim; j++)
		for (k=1; k <= sgrid_ptr->s_dim; k++)
		  if (flag(j,k) == 0)
		    flag3(i, range3(1,2)-1+j, range3(1,3)-1+k) = 0;
	    }
	  else if (sgrid_ptr->dir == 2)
	    {
	      if (sgrid_ptr->side == 1)
		j = range3(1,2);
	      else
		j = range3(2,2);
	      for (i=1; i <= sgrid_ptr->s_dim; i++)
		for (k=1; k <= sgrid_ptr->r_dim; k++)
		  if (flag(k,i) == 0)
		    flag3(range3(1,1)-1+i, j, range3(1,3)-1+k) = 0;
	    }
	  else if (sgrid_ptr->dir == 3)
	    {
	      if (sgrid_ptr->side == 1)
		k = range3(1,3);
	      else
		k = range3(2,3);
	      for (i=1; i <= sgrid_ptr->r_dim; i++)
		for (j=1; j <= sgrid_ptr->s_dim; j++)
		  if (flag(i,j) == 0)
		    flag3(range3(1,1)-1+i, range3(1,2)-1+j, k) = 0;
	    }
	  else
	    printf("Warning: cut_volume_holes: impossible direction = %i\n", 
		   sgrid_ptr->dir);
	} /* end for all faces */
    } /* end for all vgrids */
#endif

/* done */
  return OK;
}


static int
inside_global_domain(real xp, real yp, real zp, int surface_label, 
		     overlapping_3d_grid *over3d, int debug )
{
  hole_surface_link *overg_link;
  hole_surface *overg;
  int inside_this_part;
  real point[3];

  for (overg_link = over3d->surface_grid_list->first;
       overg_link != NULL; overg_link = overg_link->next)
    {
      overg = overg_link->data;
      /* a point with matching surface-label is by definition on the */
      /* inside of the cutting surface */
      if (surface_label != overg->surface_label)
	{
	  point[0] = xp;
	  point[1] = yp;
	  point[2] = zp;
	  inside_this_part = inside_surface( overg, point, 0, 0, debug );
/* tmp */
/* 	  printf("inside_global_domain: %s of %s (%s)\n",  */
/* 		 (inside_this_part? "Inside": "Outside"), overg->name,  */
/* 		 (overg->inside? "outer surface" : "inner surface")); */
	  /* the point (xp, yp, zp) is outside of the computational domain if it is */
	  /* not inside this part, but supposed to be, or */
	  /* inside this part, but not supposed to be */
	  if ((!inside_this_part &&  overg->inside) ||
	      (inside_this_part  && !overg->inside) )
	    return FALSE;
	} /* if NOT matching surface_label */
    }
  /* the point is inside of the computational domain if it wasn't */
  /* outside of any part of the boundary */
  return TRUE;
}


static int
periodic_contiguous(interp_3d_point *a, interp_3d_point *b, component_vgrid *vgrid){
  int i1_dist, i2_dist, i3_dist;
  if (vgrid->r1_period) 
    i1_dist = int_min(abs(a->i_loc - b->i_loc), 
		      range3(2,1) - range3(1,1) - abs(a->i_loc - b->i_loc));
  else
    i1_dist = abs(a->i_loc - b->i_loc);

  if (vgrid->r2_period) 
    i2_dist = int_min(abs(a->j_loc - b->j_loc), 
		      range3(2,2) - range3(1,2) - abs(a->j_loc - b->j_loc));
  else
    i2_dist = abs(a->j_loc - b->j_loc);

  if (vgrid->r3_period) 
    i3_dist = int_min(abs(a->k_loc - b->k_loc), 
		      range3(2,3) - range3(1,3) - abs(a->k_loc - b->k_loc));
  else
    i3_dist = abs(a->k_loc - b->k_loc);

  return (i1_dist <= 1 && i2_dist <= 1 && i3_dist <= 1);
}

static void
mark_around_tri(triangle tri, interp_3d_point *interp[3],
		int surface_label, component_vgrid *vgrid, 
		overlapping_3d_grid *over3d, interp_3d_point *initial_guess,
		int level)
{
/* Mark the cells in grid vgrid that are intersected by the triangle tri. */
/* the surface_label corresponds to the triangle tri. */
  interp_3d_point *sub_interp[3], *tmp_interp_1, *tmp_interp_2;
  int i, ip, jp, kp, i0, dum;
  triangle sub_tri;

/* don't repeat the recursion ab inito */
  if (level > 10) return;

  /* investigate each vertex of the triangle to see if it is inside a cell in vgrid */
  for (i=0; i<3; i++)
    if (!interp[i] && 
	(interp[i] = invert_grid(tri.x[i], tri.y[i], tri.z[i], vgrid, surface_label, 0, 
				 initial_guess)))
      {
	/* check if each vertex of the enclosing cell is inside the computational domain */
	for (ip = interp[i]->i_loc; ip <= interp[i]->i_loc+1; ip++)
	  for (jp = interp[i]->j_loc; jp <= interp[i]->j_loc+1; jp++)
	    for (kp = interp[i]->k_loc; kp <= interp[i]->k_loc+1; kp++)
	      {
		/* only necessary to mark each vertex once */
		if (flag3(ip,jp,kp) == vgrid->priority)
		  {
		    if (inside_global_domain(x3(ip,jp,kp), y3(ip,jp,kp), z3(ip,jp,kp), 
					     get_surface_label(ip, jp, kp, vgrid, 
							       &dum, &dum), 
					     over3d, 0))
		      flag3(ip,jp,kp) = 'w';
		    else
		      flag3(ip,jp,kp) = 'h';
		  } /* end if unmarked */
	      } /* end for all vertices in the enclosing cell */
	initial_guess = interp[i];
      } /* end if separated_from_surface, for all corners of the triangle */
  /* all corners inside of vgrid */
  if (interp[0] && interp[1] && interp[2])
    {
/* check if the cells are contiguous (take periodicity into account! )*/
      if (!periodic_contiguous(interp[0], interp[1], vgrid) ||
	  !periodic_contiguous(interp[0], interp[2], vgrid) ||
	  !periodic_contiguous(interp[1], interp[2], vgrid))
	{
	  /* subtriangle 1 */
/*	  printf("Marking sub-triangle 1\n");*/
	  sub_tri.x[0] = tri.x[0];   sub_tri.y[0] = tri.y[0];   sub_tri.z[0] = tri.z[0];
	  sub_interp[0] = interp[0];
	  sub_tri.x[1] = 0.5*(tri.x[0] + tri.x[1]); 
	  sub_tri.y[1] = 0.5*(tri.y[0] + tri.y[1]); 
	  sub_tri.z[1] = 0.5*(tri.z[0] + tri.z[1]);
	  sub_interp[1] = NULL;
	  sub_tri.x[2] = 0.5*(tri.x[0] + tri.x[2]); 
	  sub_tri.y[2] = 0.5*(tri.y[0] + tri.y[2]); 
	  sub_tri.z[2] = 0.5*(tri.z[0] + tri.z[2]);
	  sub_interp[2] = NULL;
	  /* recursive call to mark sub-triangle 1 */
	  mark_around_tri(sub_tri, sub_interp, surface_label, vgrid, over3d, 
			  interp[0], level+1);
	  /* saving important information! */
	  tmp_interp_1 = sub_interp[1]; 
	  tmp_interp_2 = sub_interp[2]; 
	  /* subtriangle 2 */
/*	  printf("Marking sub-triangle 2\n");*/
	  sub_tri.x[0] = 0.5*(tri.x[1] + tri.x[2]); 
	  sub_tri.y[0] = 0.5*(tri.y[1] + tri.y[2]); 
	  sub_tri.z[0] = 0.5*(tri.z[1] + tri.z[2]);
	  sub_interp[0] = NULL;
	  /* corner [1] ok from sub-triangle 1 */
	  sub_tri.x[2] = tri.x[1];
	  sub_tri.y[2] = tri.y[1]; 
	  sub_tri.z[2] = tri.z[1];
	  sub_interp[2] = interp[1]; 
	  /* recursive call to mark sub-triangle 2 */
	  mark_around_tri(sub_tri, sub_interp, surface_label, vgrid, over3d, 
			  interp[1], level+1);
	  /* subtriangle 3 */
/*	  printf("Marking sub-triangle 3\n");*/
	  /* corner [0] ok from sub-triangle 2 */
	  sub_tri.x[1] = tri.x[2];
	  sub_tri.y[1] = tri.y[2]; 
	  sub_tri.z[1] = tri.z[2];
	  sub_interp[1] = interp[2]; 
	  sub_tri.x[2] = 0.5*(tri.x[0] + tri.x[2]); 
	  sub_tri.y[2] = 0.5*(tri.y[0] + tri.y[2]); 
	  sub_tri.z[2] = 0.5*(tri.z[0] + tri.z[2]);
	  sub_interp[2] = tmp_interp_2; 
	  /* recursive call to mark sub-triangle 3 */
	  mark_around_tri(sub_tri, sub_interp, surface_label, vgrid, over3d, 
			  interp[2], level+1);
	  /* subtriangle 4 */
/*	  printf("Marking sub-triangle 4\n");*/
	  /* corner [0] and [2] ok from sub-triangle 3 */
	  sub_tri.x[1] = 0.5*(tri.x[0] + tri.x[1]); 
	  sub_tri.y[1] = 0.5*(tri.y[0] + tri.y[1]); 
	  sub_tri.z[1] = 0.5*(tri.z[0] + tri.z[1]);
	  sub_interp[1] = tmp_interp_1; 
	  /* recursive call to mark sub-triangle 4 */
	  mark_around_tri(sub_tri, sub_interp, surface_label, vgrid, over3d, 
			  interp[2], level+1);
/* cleanup */
	  for (i0 = 0; i0 < 3; i0++)
	    if (sub_interp[i0]) free(sub_interp[i0]);
	}
/*       else */
/* 	printf("Found a triangle that DON'T need to be subdivided\n"); */
    }
  
} /* end mark_around_tri */


/* classify all gridpoints as either interpolation points , 
   discretization points or exterior points */
void
initial_classification(linked_list_member *vgrid_link, overlapping_3d_grid *over3d){
  int i, j, k, present_surface, present_curve, dum;
  interp_3d_point *new_interp3, prev_interp;
  component_vgrid *vgrid;

  vgrid = vgrid_link->data;
/* don't do this for dummy grids */
/*  if (vgrid->dummy_background) return;*/

  if (over3d->verbose >= 2)
    printf("Main classification of all points in grid `%s'\n",vgrid->name);
  for (i = range3(1,1); i <= range3(2,1); i++){
    for (j = range3(1,2); j <= range3(2,2); j++){
/* the first point in each column has no initial guess */
      prev_interp.r_loc = 0.0; prev_interp.s_loc = 0.0; prev_interp.t_loc = 0.0;
      prev_interp.vgrid_loc = NULL;
/* it might be advantageous to do the k loop first, since that corresponds to the */
/* stretched direction in grids generated by the normal surface or hyperbolic */
/* grid generators */
      for (k = range3(1,3); k <= range3(2,3); k++)
/* don't mess with the points that previously were flagged dead */
	if (flag3(i,j,k) != 0){
/* set the curve label if the present point is on a real boundary */
	    present_surface = get_surface_label(i, j, k, vgrid, &dum, &dum);
	    present_curve   = get_curve_label(i, j, k, vgrid);

/* test if (i,j,k) can interpolalate from a higher grid */
	    if ((new_interp3 = 
		 interp_from_higher(x3(i,j,k), y3(i,j,k), z3(i,j,k), i, j, k,
				    &prev_interp, present_surface, present_curve,
				    vgrid_link, over3d)) != NULL)
	      flag3(i,j,k) = new_interp3->vgrid_loc->priority;
/* test if it can be an interior point */
	    else if ( disc_point( i, j, k, vgrid, over3d) )
	      flag3(i,j,k) = vgrid->priority;
/* test if (i,j,k) can interpolate from a lower grid */
	    else if ((new_interp3 = 
		      interp_from_lower(x3(i,j,k), y3(i,j,k), z3(i,j,k), i, j, k,
					&prev_interp, present_surface, present_curve,
					vgrid, vgrid->priority, over3d))
		     != NULL)
	      flag3(i,j,k) = new_interp3->vgrid_loc->priority;
	    else
/* dead point */
	      flag3(i,j,k) = 0;

/* insert the new interpolation point in the list */
	    if (new_interp3 != NULL){
	      insert_3d_interp_point( new_interp3, vgrid );
/* use the location as an initial guess for next point */
	      prev_interp.r_loc = new_interp3->r_loc; 
	      prev_interp.s_loc = new_interp3->s_loc; 
	      prev_interp.t_loc = new_interp3->t_loc; 
	      prev_interp.vgrid_loc = new_interp3->vgrid_loc;
	    }
	  
	  } /* end if flag != 0 */
	} /* end for k */
      } /* end for i,j */
}


void
explicit(component_vgrid *vgrid, overlapping_3d_grid *over3d){
  component_vgrid *donor_vgrid;
  int i, j, k, iw1, implicit, points_removed, convert_dp, lower_interp, present_surface
    , present_curve, did_lower, new_flag, convert_ip, points_killed
      , i_loc_min, i_loc_max, j_loc_min, j_loc_max
	, k_loc_min, k_loc_max, unfixable, dum;
  interp_3d_point *interp3, *donor_interp3, *new_interp3;

#define donor_flag3(i,j,k)  compute_index_3d(donor_vgrid->flag_ptr,i,j,k)
#define donor_range(i,j)  compute_index_2d(donor_vgrid->range_ptr,i,j)

  iw1 = over3d->interp_width - 1;

  if (over3d->verbose >= 2)
    printf("\nExplicit> Grid `%s'.\n",vgrid->name);
  points_removed = 0;
  points_killed = 0;
  convert_dp = 0;
  convert_ip = 0;
  lower_interp = 0;
  unfixable = 0;
    for (interp3 = vgrid->last_interp; interp3 != NULL; 
	 interp3 = interp3->prev){

/* only consider active interpolation points */
      if (!interp3->active) continue;
/* make sure all the points in the interpolation formula are discretization */
/* points */
	do {
	  donor_vgrid = interp3->vgrid_loc;
	  did_lower = FALSE;
	  implicit = FALSE;

	  i_loc_min = interp3->i_loc;
	  i_loc_max = interp3->i_loc + iw1;
	  j_loc_min = interp3->j_loc;
	  j_loc_max = interp3->j_loc + iw1;
	  k_loc_min = interp3->k_loc;
	  k_loc_max = interp3->k_loc + iw1;

	  for (i=i_loc_min; i<=i_loc_max; i++)
	    for (j=j_loc_min; j<=j_loc_max; j++)
	      for (k=k_loc_min; k<=k_loc_max; k++){
		if (abs(donor_flag3(i,j,k)) != donor_vgrid->priority){
/* first try to convert (i,j,k) in the donor grid to a discretization point */
		  if (disc_point( i, j, k, donor_vgrid, over3d)){
/* find the interpolation point in the donor grids interpolation list and
		     deactivate it */
		    convert_dp++;
		    for (donor_interp3 = donor_vgrid->last_interp;
			 donor_interp3 != NULL; donor_interp3 = donor_interp3->prev){
		      if (donor_interp3->active && 
			  donor_interp3->i_point == i && 
			  donor_interp3->j_point == j && 
			  donor_interp3->k_point == k){
			donor_interp3->active = 0;
			donor_vgrid->n_interp--;
/* also change the point's flag. Just as a test, make it positive, so that it isn't */
/* a candidate for an orphan point */
			donor_flag3(i,j,k) = donor_vgrid->priority;
		      }
		    }
		  }
		  else
/* this implicit interpolation was not possible to repair by changing the */
/* interpolation location to a discretization point */
		    implicit = TRUE;
		} /* end if */
	      } /* end for all donor points for the present interpolation point */
	
/* change the interpolation point to discretization, lower interpolate */
/* or unused if implicit interpolation was encountered */
	  if (implicit){
	    i = interp3->i_point; j = interp3->j_point; k = interp3->k_point;

/* set the surface label if the present point is on a real boundary */
	    present_surface = get_surface_label( i, j, k, vgrid, &dum, &dum);
	    present_curve   = get_curve_label(i, j, k, vgrid);

/* it should first be checked if the point can interpolate from a grid which is */
/* lower than the present donor but higher than itself */

	    new_interp3 = interp_from_lower(x3(i,j,k), y3(i,j,k), z3(i,j,k), i, j, k,
					    NULL, present_surface, present_curve, vgrid, 
					    interp3->vgrid_loc->priority, over3d);
	    new_flag = (new_interp3 == NULL)? 0: new_interp3->vgrid_loc->priority;

/* if the new donor priority is lower than the present grid, try to turn it */
/* into a discretization point instead */
	    if ( new_flag < vgrid->priority &&
		disc_point(i, j, k, vgrid, over3d) ){
	      convert_ip++;
	      flag3(i, j, k) = vgrid->priority;

/* deactivate this interpolation point */
	      interp3->active = 0;
	      vgrid->n_interp--;
	      points_removed++;
#ifdef D_EXPLICIT
 	      printf("The implicit interpolation point (%i, %i, %i) "
		     "was reclassified into a discretization point.\n", 
		     i, j, k);
	      if (new_interp3 != NULL)
		printf("It could also have interpolated from grid `%s'. I'm saving"
		       "that as a backup solution.\n", new_interp3->vgrid_loc->name); 
#endif
	    }
	    else if ( new_flag != 0 ){
/* try to interpolate from a grid with lower priority, or remove the point */
	      flag3(i, j, k) = new_flag;
/* update the interpolation information*/
	      interp3->vgrid_loc = new_interp3->vgrid_loc;
	      interp3->r_loc = new_interp3->r_loc;
	      interp3->s_loc = new_interp3->s_loc;
	      interp3->t_loc = new_interp3->t_loc;
	      interp3->i_loc = new_interp3->i_loc;
	      interp3->j_loc = new_interp3->j_loc;
	      interp3->k_loc = new_interp3->k_loc;
	      lower_interp++;
	      did_lower = TRUE;
/* free the temporary interpolation point */
	      free( new_interp3 );
#ifdef D_EXPLICIT
 	      printf("The implicit interpolation point (%i, %i, %i) "
		     "now interpolates from grid `%s'\n", i, j, k, 
		     interp3->vgrid_loc->name); 
#endif
	    } 
	    else{
	      unfixable++;
	    }

	  }
	} while (did_lower);

    } /* end for all interpolation points */

/* print a report of the changes on the present grid */
  if (over3d->verbose >= 2){
    if (convert_dp > 0)
      printf("Converted %i implicit donor points to discretization points\n", 
	     convert_dp);
    if (lower_interp > 0)
      printf("Lowered the interpolation location for %i implicit interpolation "
	     "points\n", lower_interp);
    if (convert_ip > 0)
      printf("Converted %i implicit interpolation points to discretization points\n", 
	     convert_ip);
    if (points_killed > 0)
      printf("Converted %i implicit interpolation points to exterior points\n\n", 
	     points_killed);
    if (unfixable > 0)
      printf("Warning: there were %i unfixable points\n", unfixable);
  }
#undef donor_range
#undef donor_flag3
}

int
consistency_check(component_vgrid *vgrid, overlapping_3d_grid *over3d){
  int iw1, ok=1, bad_disc, bad_interp, bad_explicit, n_interp, i, j, k, ok_interp
    , ok_explicit, i_loc_min, i_loc_max, j_loc_min, j_loc_max, k_loc_min, k_loc_max;
  component_vgrid *donor_vgrid;
  interp_3d_point *interp3;

  iw1 = over3d->interp_width - 1;
/* don't do this for dummy grids */
/*  if (vgrid->dummy_background) return ok;*/

/* interior points */
  bad_disc = 0;
  for (i = range3(1,1); i <= range3(2,1); i++)
    for (j = range3(1,2); j <= range3(2,2); j++)
      for (k = range3(1,3); k <= range3(2,3); k++){
	if (flag3(i,j,k) == vgrid->priority && 
	    !disc_point( i, j, k, vgrid, over3d)){
	  new_3d_bad_point(i, j, k, BAD_DISC, vgrid);
	  bad_disc++;
	}
      }
/* interpolation points */
  bad_interp = 0;
  bad_explicit = 0;
  n_interp = 0;

#define donor_range(i,j)  compute_index_2d(donor_vgrid->range_ptr,i,j)
#define donor_flag3(i,j,k)  compute_index_3d(donor_vgrid->flag_ptr,i,j,k)

  for (interp3 = vgrid->last_interp; interp3 != NULL; 
       interp3 = interp3->prev){
      
    if (interp3->active){
      n_interp++;

      donor_vgrid = interp3->vgrid_loc;
      ok_interp = 1;
      ok_explicit = 1;

      i_loc_min = interp3->i_loc;
      i_loc_max = interp3->i_loc + iw1;
      j_loc_min = interp3->j_loc;
      j_loc_max = interp3->j_loc + iw1;
      k_loc_min = interp3->k_loc;
      k_loc_max = interp3->k_loc + iw1;

      for (i = i_loc_min; i <= i_loc_max; i++)
	for (j = j_loc_min; j <= j_loc_max; j++)
	  for (k = k_loc_min; k <= k_loc_max; k++){
	    if (donor_flag3(i,j,k) == 0){
	      new_3d_bad_point(i, j, k, DEAD_INTERP_LOC, donor_vgrid);
	      ok_interp = 0;
	    }
/* check for non-explicit interpolation */
	    else if (over3d->interp_type == 'e' && 
		     abs(donor_flag3(i,j,k)) != donor_vgrid->priority){ 
	      new_3d_bad_point(i, j, k, NON_EXPLICIT_LOC, donor_vgrid);
	      ok_explicit = 0; 
	    } 
	  }	    
      if (!ok_interp){
	new_3d_bad_point(interp3->i_point, interp3->j_point, interp3->k_point, 
			 BAD_INTERP_DEAD_LOC, vgrid); 
	bad_interp++;
      }
      else if (!ok_explicit){
	new_3d_bad_point(interp3->i_point, interp3->j_point, interp3->k_point,
			 BAD_INTERP_NON_EXPLICIT, vgrid); 
	bad_explicit++;
      }
    }
    else{
      printf("Warning: A deactivated interpolation point found while checking ");
      printf("         the constistency\n");
    }
  }
/* report if there are any errors */
    if (bad_disc > 0 || bad_interp > 0 || bad_explicit > 0){
      ok = 0;
      if (over3d->verbose >= 2){
	printf("\nInconsistencies found in grid `%s'.\n", vgrid->name);
	if (bad_disc>0)
	  printf("It has %i bad discretization points\n", bad_disc);
	if (bad_interp>0)
	  printf("It has %i bad interpolation points because of dead donors\n", 
		 bad_interp);
	if (over3d->interp_type == 'e' && bad_explicit>0)
	  printf("It has %i non-explicit interpolation points\n", bad_explicit);
      }
    }
    
    if (n_interp != vgrid->n_interp)
      printf("Warning: n_interp missmatch in grid `%s'\n",vgrid->name);

  return ok;

#undef donor_flag3
#undef donor_range
}


void
trim_interpolation_points(component_vgrid *vgrid, overlapping_3d_grid *over3d){
  interp_3d_point *interp3;
  int iw1, points_removed, i, j, k;
  component_vgrid_link *other_link;
  component_vgrid *other_vgrid;
  
/* don't do this for dummy grids */
/*  if (vgrid->dummy_background) return;*/

  iw1 = over3d->interp_width - 1;
  
  if (over3d->verbose >= 2)
    printf("Trimming grid `%s'\n", vgrid->name);
    
/****************************/
/* trimstyle == 1 starts here */
/****************************/
  if (over3d->trim_style == 1){
/* if this interpolation point can be an interior point instead, 
   set flag3(i_point,j_point) = vgrid->priority and delete this 
   interpolation point from the stack */
      points_removed = 0;

      for (interp3 = vgrid->last_interp; interp3 != NULL; 
	   interp3 = interp3->prev){
      
	if (interp3->active &&
	    disc_point(interp3->i_point, interp3->j_point, interp3->k_point, vgrid, 
		       over3d)){
	
/* change the flag3 and remove this interpolation point from the stack */
	
	  flag3(interp3->i_point,interp3->j_point,interp3->k_point) = vgrid->priority;
	  interp3->active = 0;
	  points_removed++;
	}
      }

/* update the number of interpolation points */
      vgrid->n_interp += -points_removed;
      if (over3d->verbose >= 2 && points_removed > 0)
	printf("\tChanged %i interpolation points to discretization points\n", 
	       points_removed);
    }

/****************************/
/* trimstyle == 1 ends here */
/****************************/

/* search all interpolation points in all other grids for donors in the */
/* present grid */
    if (over3d->interp_type == 'i')
      {
	for (other_link = over3d->grid_list->first; other_link != NULL;
	     other_link = other_link->next)
	  {
	    other_vgrid = other_link->data;
	    if (other_vgrid == vgrid) continue;
	    for (interp3 = other_vgrid->last_interp; interp3 != NULL;
		 interp3 = interp3->prev)
	      {
		if (interp3->vgrid_loc == vgrid)
		  {
/* mark the donor points */
		    for (i=interp3->i_loc; i<=interp3->i_loc+iw1; i++)
		      for (j=interp3->j_loc; j<=interp3->j_loc+iw1; j++)
			for (k=interp3->k_loc; k<=interp3->k_loc+iw1; k++)
			  flag3(i,j,k) = -abs(flag3(i,j,k));
		  } /* end if donor is vgrid */
	      } /* end for all interpolation points */
	  } /* end for all grids except vgrid */
      } /* end if implicit interpolation */

    points_removed = 0;
    for (interp3 = vgrid->last_interp; interp3 != NULL; 
	 interp3 = interp3->prev){
      
/* don't mess with unused or marked points, i.e. points where flag <= 0 */
      if (flag3(interp3->i_point,interp3->j_point,interp3->k_point) > 0){

/* if this interpolation point is not needed, set flag3(i_point,j_point) = 0
   and delete this interpolation point from the stack */

	if (!needed_by_disc(interp3->i_point, interp3->j_point, interp3->k_point,
			    vgrid, over3d) 
/* 	    &&  */
/* 	    !needed_by_interp(interp3->i_point, interp3->j_point, interp3->k_point, */
/* 			      vgrid->priority, over3d)  */
	    ){
/* remove this interpolation point */
	  flag3(interp3->i_point,interp3->j_point,interp3->k_point) = 0;
	  interp3->active = 0;
	  points_removed++;
	}
      }
    }

/* update the number of interpolation points */
    vgrid->n_interp += -points_removed;
    if (over3d->verbose >= 2 && points_removed > 0) 
      printf("\t%i interpolation points not needed\n", points_removed);

/* minimize the width of the overlap */
    if (over3d->trim_style == 0){
/* if this interpolation point can be an interior point instead, 
   set flag3(i_point,j_point,k_point) = vgrid->priority and delete this 
   interpolation point from the stack */
      points_removed = 0;

      for (interp3 = vgrid->last_interp; interp3 != NULL; 
	   interp3 = interp3->prev){
      
	if (interp3->active &&
	    disc_point(interp3->i_point, interp3->j_point, interp3->k_point, vgrid, 
		       over3d)){
	
/* change the flag and remove this interpolation point from the stack */
	
	  flag3(interp3->i_point,interp3->j_point,interp3->k_point) = vgrid->priority;
	  interp3->active = 0;
	  points_removed++;
	}
      }

/* update the number of interpolation points */
      vgrid->n_interp += -points_removed;
      if (over3d->verbose >= 2 && points_removed > 0)
	printf("\tChanged %i interpolation points to discretization points\n", 
	       points_removed);
    }

    /* remove the de-activated interpolation points */
    remove_deactive_3d_interp( vgrid );

    /* the marking is only necessary when the interpolation type is implicit */
    if (over3d->interp_type == 'i'){
      for (interp3 = vgrid->last_interp; interp3 != NULL;
	   interp3 = interp3->prev)
	if (interp3->vgrid_loc->priority > vgrid->priority)
/* set flag = - |flag| for all points in the interpolation formula */
	  change_sign(interp3->i_loc, interp3->j_loc, interp3->k_loc, iw1,
		      interp3->vgrid_loc);
    }
}

/* change sign of the entries in the flag array so that discretization points */
/* are positive and interpolation points are negative. Also remove deactivated */
/* interpolation points */

void
cleanup_and_finish( component_vgrid *vgrid, overlapping_3d_grid *over3d ){
  interp_3d_point *interp3;
  int i, j, k, n_interp;

/* don't do this for dummy grids */
/*  if (vgrid->dummy_background) return;*/

/* update the sign of the flag array */
    for (i = range3(1,1); i <= range3(2,1); i++)
      for (j = range3(1,2); j <= range3(2,2); j++)
	for (k = range3(1,3); k <= range3(2,3); k++)
	  if (abs(flag3(i,j,k)) == vgrid->priority)
	    flag3(i,j,k) = abs(flag3(i,j,k));
	  else if (flag3(i,j,k) != 0)
	    flag3(i,j,k) = -abs(flag3(i,j,k));

/* take care of the periodic case */
    if (vgrid->r1_period){
/* remove any spurious signs from the flag array for the periodic points */
      for (j = 1; j <= vgrid->r2_dim; j++)
	for (k = 1; k <= vgrid->r3_dim; k++){
	  for (i = 1; i <= range3(1,1)-1; i++)
	    if (abs(flag3(i,j,k)) == vgrid->priority)
	      flag3(i,j,k) = abs(flag3(i,j,k));
	    else if (flag3(i,j,k) != 0)
	      flag3(i,j,k) = -abs(flag3(i,j,k));
	  for (i = vgrid->r1_dim; i >= range3(2,1); i--)
	    if (abs(flag3(i,j,k)) == vgrid->priority)
	      flag3(i,j,k) = abs(flag3(i,j,k));
	    else if (flag3(i,j,k) != 0)
	      flag3(i,j,k) = -abs(flag3(i,j,k));
	}

/* remove any interpolation points that have i_point=range3(2,1) */
      for (interp3 = vgrid->last_interp; interp3 != NULL;
	   interp3 = interp3->prev){

	if (interp3->i_point == range3(2,1)){
/* change the flag and remove this interpolation point from the stack */
	  flag3(interp3->i_point,interp3->j_point,interp3->k_point) = vgrid->priority;
	  interp3->active = 0;  
	  vgrid->n_interp--;
	}
      }

/* make sure the flag array becomes periodic */
      for (j = 1; j <= vgrid->r2_dim; j++)
	for (k = 1; k <= vgrid->r3_dim; k++)
	  flag3(range3(2,1),j,k) = flag3(range3(1,1),j,k);
    }
    if (vgrid->r2_period){
/* remove any spurious signs from the flag array for the periodic points */
      for (i = 1; i <= vgrid->r1_dim; i++)
	for (k = 1; k <= vgrid->r3_dim; k++){
	  for (j = 1; j <= range3(1,2)-1; j++)
	    if (abs(flag3(i,j,k)) == vgrid->priority)
	      flag3(i,j,k) = abs(flag3(i,j,k));
	    else if (flag3(i,j,k) != 0)
	      flag3(i,j,k) = -abs(flag3(i,j,k));
	  for (j = vgrid->r2_dim; j >= range3(2,2); j--)
	    if (abs(flag3(i,j,k)) == vgrid->priority)
	      flag3(i,j,k) = abs(flag3(i,j,k));
	    else if (flag3(i,j,k) != 0)
	      flag3(i,j,k) = -abs(flag3(i,j,k));
	}
      
/* remove any interpolation points that have j_point=range3(2,2) */
      for (interp3 = vgrid->last_interp; interp3 != NULL;
	   interp3 = interp3->prev){

	if (interp3->j_point == range3(2,2)){
/* change the flag and remove this interpolation point from the stack */
	  flag3(interp3->i_point,interp3->j_point,interp3->k_point) = vgrid->priority;
	  interp3->active = 0;
	  vgrid->n_interp--;
	}
      }
/* make sure the flag array becomes periodic */
      for (i = 1; i <= vgrid->r1_dim; i++)
	for (k = 1; k <= vgrid->r3_dim; k++)
	  flag3(i,range3(2,2),k) = flag3(i,range3(1,2),k);
    }
    if (vgrid->r3_period){
/* remove any spurious signs from the flag array for the periodic points */
      for (i = 1; i <= vgrid->r1_dim; i++)
	for (j = 1; j <= vgrid->r2_dim; j++){
	  for (k = 1; k <= range3(1,3)-1; k++)
	    if (abs(flag3(i,j,k)) == vgrid->priority)
	      flag3(i,j,k) = abs(flag3(i,j,k));
	    else if (flag3(i,j,k) != 0)
	      flag3(i,j,k) = -abs(flag3(i,j,k));
	  for (k = vgrid->r3_dim; k >= range3(2,3); k--)
	    if (abs(flag3(i,j,k)) == vgrid->priority)
	      flag3(i,j,k) = abs(flag3(i,j,k));
	    else if (flag3(i,j,k) != 0)
	      flag3(i,j,k) = -abs(flag3(i,j,k));
	}
      
/* remove any interpolation points that have j_point=range3(2,2) */
      for (interp3 = vgrid->last_interp; interp3 != NULL;
	   interp3 = interp3->prev){

	if (interp3->k_point == range3(2,3)){
/* change the flag and remove this interpolation point from the stack */
	  flag3(interp3->i_point,interp3->j_point,interp3->k_point) = vgrid->priority;
	  interp3->active = 0;
	  vgrid->n_interp--;
	}
      }
/* make sure the flag array becomes periodic */
      for (i = 1; i <= vgrid->r1_dim; i++)
	for (j = 1; j <= vgrid->r2_dim; j++)
	  flag3(i,j,range3(2,3)) = flag3(i,j,range3(1,3));
    }

/* remove any de-activated interpolation points */
    remove_deactive_3d_interp( vgrid );

/* check the number of interpolation points */
    n_interp = 0;
    for (interp3 = vgrid->last_interp; interp3 != NULL;
	 interp3 = interp3->prev)
      n_interp++;
    if (n_interp != vgrid->n_interp){
      printf(
"Warning: Correcting the number of interpolation points in grid `%s'\n"
	     ,vgrid->name);
      printf(
"         vgrid->n_interp = %i, while the actual numer was %i\n"
	     ,vgrid->n_interp, n_interp);
      vgrid->n_interp = n_interp;
    }

/* count the number of dummy interpolation points */
/*     if (over3d->overlapping_sub_grid){ */
/*       vgrid->n_dummies = 0; */
/*       for (interp3 = vgrid->last_interp; interp3 != NULL; */
/* 	   interp3 = interp3->prev) */
/* 	if (interp3->vgrid_loc->dummy_background){ */
/* 	  over3d->n_dummy_interp++; */
/* 	  vgrid->n_dummies++; */
/* 	} */
/*    } */ /* end if overlapping_sub_grid */

} /* end cleanup_and_finish */

static void 
change_sign(int i_loc, int j_loc, int k_loc, int iw1, component_vgrid *vgrid){
  int i,j,k;
  for (i=i_loc; i<=i_loc+iw1; i++)
    for (j=j_loc; j<=j_loc+iw1; j++)
      for (k=k_loc; k<=k_loc+iw1; k++)
	flag3(i,j,k) = - abs(flag3(i,j,k));
}


static int 
needed_by_disc(int i_point, int j_point, int k_point, component_vgrid *vgrid,
	       overlapping_3d_grid *over3d){
/* check if it is needed by a discretization point in the same grid */
  int i_min, j_min, i_max, j_max, k_min, k_max, i, j, k, dw2, nw1, tw2;

  dw2 = (over3d->disc_width - 1)/2;

/* needed by an interior point? */
  i_min = int_max(i_point - dw2, range3(1,1));
  i_max = int_min(i_point + dw2, range3(2,1));
  j_min = int_max(j_point - dw2, range3(1,2));
  j_max = int_min(j_point + dw2, range3(2,2));
  k_min = int_max(k_point - dw2, range3(1,3));
  k_max = int_min(k_point + dw2, range3(2,3));

  for (i=i_min; i<= i_max; i++)
    for (j=j_min; j<= j_max; j++)
      for (k=k_min; k<= k_max; k++)
	if (abs(flag3(i,j,k)) == vgrid->priority)
	  return 1;

/* also check if the point is needed by a boundary formula */
  nw1 = over3d->normal_width-1;
  tw2 = (over3d->tangent_width-1)/2;

/* needed by a boundary point on the left side? */
  if (surf3(1,1) != 0 && i_point <= range3(1,1) + nw1){
    i_min = range3(1,1);
    i_max = range3(1,1)+nw1;
    j_min = int_max( j_point - tw2, range3(1,2) );
    j_max = int_min( j_point + tw2, range3(2,2) );
    k_min = int_max( k_point - tw2, range3(1,3) );
    k_max = int_min( k_point + tw2, range3(2,3) );
    for (k=k_min; k<=k_max; k++)
      for (j=j_min; j<= j_max; j++)
	for (i=i_min; i<= i_max; i++)
	  if (abs(flag3(i,j,k)) == vgrid->priority)
	    return 1;
  }

/* needed by a boundary point on the right side? */
  if (surf3(2,1) != 0 && i_point >= range3(2,1) - nw1){
    i_min = range3(2,1)-nw1;
    i_max = range3(2,1);
    j_min = int_max( j_point - tw2, range3(1,2) );
    j_max = int_min( j_point + tw2, range3(2,2) );
    k_min = int_max( k_point - tw2, range3(1,3) );
    k_max = int_min( k_point + tw2, range3(2,3) );
    for (k=k_min; k<=k_max; k++)
      for (j=j_min; j<= j_max; j++)
	for (i=i_min; i<= i_max; i++)
	  if (abs(flag3(i,j,k)) == vgrid->priority)
	    return 1;
  }

/* needed by a boundary point on the lower side? */
  if (surf3(1,2) != 0 && j_point <= range3(1,2) + nw1){
    i_min = int_max( i_point - tw2, range3(1,1) );
    i_max = int_min( i_point + tw2, range3(2,1) );
    j_min = range3(1,2);
    j_max = range3(1,2) + nw1;
    k_min = int_max( k_point - tw2, range3(1,3) );
    k_max = int_min( k_point + tw2, range3(2,3) );
    for (k=k_min; k<=k_max; k++)
      for (j=j_min; j<= j_max; j++)
	for (i=i_min; i<= i_max; i++)
	  if (abs(flag3(i,j,k)) == vgrid->priority)
	    return 1;
  }

/* needed by a boundary point on the upper side? */
  if (surf3(2,2) != 0 && j_point >= range3(2,2) - nw1){
    i_min = int_max( i_point - tw2, range3(1,1) );
    i_max = int_min( i_point + tw2, range3(2,1) );
    j_min = range3(2,2) - nw1;
    j_max = range3(2,2);
    k_min = int_max( k_point - tw2, range3(1,3) );
    k_max = int_min( k_point + tw2, range3(2,3) );
    for (k=k_min; k<=k_max; k++)
      for (j=j_min; j<= j_max; j++)
	for (i=i_min; i<= i_max; i++)
	  if (abs(flag3(i,j,k)) == vgrid->priority)
	    return 1;
  }

/* needed by a boundary point on the near side? */
  if (surf3(1,3) != 0 && k_point <= range3(1,3) + nw1){
    i_min = int_max( i_point - tw2, range3(1,1) );
    i_max = int_min( i_point + tw2, range3(2,1) );
    j_min = int_max( j_point - tw2, range3(1,2) );
    j_max = int_min( j_point + tw2, range3(2,2) );
    k_min = range3(1,3);
    k_max = range3(1,3) + nw1;
    for (k=k_min; k<=k_max; k++)
      for (j=j_min; j<= j_max; j++)
	for (i=i_min; i<= i_max; i++)
	  if (abs(flag3(i,j,k)) == vgrid->priority)
	    return 1;
  }

/* needed by a boundary point on the upper side? */
  if (surf3(2,3) != 0 && k_point >= range3(2,3) - nw1){
    i_min = int_max( i_point - tw2, range3(1,1) );
    i_max = int_min( i_point + tw2, range3(2,1) );
    j_min = int_max( j_point - tw2, range3(1,2) );
    j_max = int_min( j_point + tw2, range3(2,2) );
    k_min = range3(2,3) - nw1;
    k_max = range3(2,3);
    for (k=k_min; k<=k_max; k++)
      for (j=j_min; j<= j_max; j++)
	for (i=i_min; i<= i_max; i++)
	  if (abs(flag3(i,j,k)) == vgrid->priority)
	    return 1;
  }

/* I suppose one should also check all 12 edges and 8 corners as well. I'll do */
/* that some other day... */

/* If you got this far, the point (i_point, j_point, k_point) */
/* is not needed by any discretization point. */

  return 0;
}

/* mark interpolation points which are needed for interpolation by higher 
   priority grids */

static int 
mark_needed_interp(overlapping_3d_grid *over3d){
  int iw1, i, j, k;
  component_vgrid *other_vgrid,*vgrid;
  interp_3d_point *this_interp3;
  linked_list_member *vgrid_link;

  /* only necessary to mark points when the interpolation type is implicit */
  if (over3d->interp_type == 'e')
    return 1;

  iw1 = over3d->interp_width - 1;
  for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL; 
       vgrid_link = vgrid_link->next){
    other_vgrid = vgrid_link->data;
    for(this_interp3 = other_vgrid->last_interp; this_interp3 != NULL;
	this_interp3 = this_interp3->prev)
/* check if this interpolation point interpolates from a grid with lower
   priority than other_vgrid->priority */
      if (this_interp3->vgrid_loc->priority < other_vgrid->priority){
/* set vgrid to point to the interpolation location, so that the
   inline macros for flag will work */
	vgrid = this_interp3->vgrid_loc;
/* change sign of all the points in that grid that are needed in that 
   interpolation formula */
	for (k=this_interp3->k_loc; k<=this_interp3->k_loc + iw1; k++)
	  for (j=this_interp3->j_loc; j<=this_interp3->j_loc + iw1; j++)
	    for (i=this_interp3->i_loc; i<=this_interp3->i_loc + iw1; i++)
	      flag3(i,j,k) = -abs(flag3(i,j,k));
      }
  }
  return 1;
}

static int
get_surface_label(int i, int j, int k, component_vgrid *vgrid, int *side, int *dir){
  real left_dist, right_dist, lower_dist, upper_dist, near_dist, far_dist, dist_2;

  dist_2 = vgrid->max_s_dist * vgrid->max_s_dist;

  left_dist = right_dist = lower_dist = upper_dist = near_dist = far_dist =
    10.0 * dist_2;

  if (surf3(1,1) != 0)
    {
      if (i == range3(1,1))
	left_dist = 0.0;
      else
	left_dist = sqr(x3(i,j,k)-x3(range3(1,1),j,k)) + 
	  sqr(y3(i,j,k)-y3(range3(1,1),j,k)) + sqr(z3(i,j,k)-z3(range3(1,1),j,k));
    } /* end if surf3(1,1) != 0 */

  if (surf3(2,1) != 0)
    {
      if (i == range3(2,1))
	right_dist = 0.0;
      else
	right_dist = sqr(x3(i,j,k)-x3(range3(2,1),j,k)) +
	  sqr(y3(i,j,k)-y3(range3(2,1),j,k)) + sqr(z3(i,j,k)-z3(range3(2,1),j,k));
    } /* end if surf3(2,1) != 0 */

  if (surf3(1,2) != 0)
    {
      if (j == range3(1,2))
	lower_dist = 0.0;
      else
	lower_dist = sqr(x3(i,j,k)-x3(i,range3(1,2),k)) +
	  sqr(y3(i,j,k)-y3(i,range3(1,2),k)) + sqr(z3(i,j,k)-z3(i,range3(1,2),k));
    } /* end if surf3(1,2) != 0 */

  if (surf3(2,2) != 0)
    {
      if (j == range3(2,2))
	upper_dist = 0.0;
      else
	upper_dist = sqr(x3(i,j,k)-x3(i,range3(2,2),k)) +
	  sqr(y3(i,j,k)-y3(i,range3(2,2),k)) + sqr(z3(i,j,k)-z3(i,range3(2,2),k));
    } /* end if surf3(2,2) != 0 */
    
  if (surf3(1,3) != 0)
    {
      if (k == range3(1,3))
	near_dist = 0.0;
      else
	near_dist = sqr(x3(i,j,k)-x3(i,j,range3(1,3))) +
	  sqr(y3(i,j,k)-y3(i,j,range3(1,3))) + sqr(z3(i,j,k)-z3(i,j,range3(1,3)));
    } /* end if surf3(1,3) != 0 */

  if (surf3(2,3) != 0)
    {
      if (k == range3(2,3))
	far_dist = 0.0;
      else
	far_dist = sqr(x3(i,j,k)-x3(i,j,range3(2,3))) +
	  sqr(y3(i,j,k)-y3(i,j,range3(2,3))) + sqr(z3(i,j,k)-z3(i,j,range3(2,3)));
    } /* end if surf3(1,3) != 0 */


/* direction 3 */
  if (near_dist <= dist_2 && surf3(1,3) != 0)
    {
      *side = 1; *dir = 3;
      return surf3(1,3);
    }
  if (far_dist <= dist_2 && surf3(2,3) != 0)
    {
      *side = 2; *dir = 3;
      return surf3(2,3);
    }

/* direction 2 */
  if (lower_dist <= dist_2 && surf3(1,2) != 0)
    {
      *side = 1; *dir = 2;
      return surf3(1,2);
    }
  if (upper_dist <= dist_2 && surf3(2,2) != 0)
    {
      *side = 2; *dir = 2;
      return surf3(2,2);
    }

/* direction 1 */
  if (left_dist <= dist_2 && surf3(1,1) != 0)
    {
      *side = 1; *dir = 1;
      return surf3(1,1);
    }
  if (right_dist <= dist_2 && surf3(2,1) != 0) 
    {
      *side = 2; *dir = 1;
      return surf3(2,1);
    }

/* if we got this far, the present point is not close to any surface with */
/* non-zero surface label */
  *side = 0; *dir = 0;
  return 0;
}

static int
on_physical_boundary(int i, int j, int k, component_vgrid *vgrid, int *dir){

  if ((surf3(1,3) && (k == range3(1,3))) ||
	   (surf3(2,3) && (k == range3(2,3))))
    {
      *dir = 3;
      return TRUE;
    }

  else if ((surf3(1,2) && (j == range3(1,2))) ||
	   (surf3(2,2) && (j == range3(2,2))))
    {
      *dir = 2;
      return TRUE;
    }

  else if ((surf3(1,1) && (i == range3(1,1))) ||
      (surf3(2,1) && (i == range3(2,1))))
    {
      *dir = 1;
      return TRUE;
    }

  return FALSE;
}

static int
get_curve_label(int i, int j, int k, component_vgrid *vgrid){
  int this_curve=0;

/* observe that this routine gives the last curve_label at an edge or corner */
/* where several physical edges meet */

  if (i == range3(1,1)){
    if (j == range3(1,2))
      this_curve = edge_curve(1);
    else if (j == range3(2,2))
      this_curve = edge_curve(3);
    else if (k == range3(1,3))
      this_curve = edge_curve(5);
    else if (k == range3(2,3))
      this_curve = edge_curve(7);
  }
  else if (i == range3(2,1)){
    if (j == range3(1,2))
      this_curve = edge_curve(2);
    else if (j == range3(2,2))
      this_curve = edge_curve(4);
    else if (k == range3(1,3))
      this_curve = edge_curve(6);
    else if (k == range3(2,3))
      this_curve = edge_curve(8);
  }
  else if (j == range3(1,2)){
    if (k == range3(1,3))
      this_curve = edge_curve(9);
    else if (k == range3(2,3))
      this_curve = edge_curve(11);
  }
  else if (j == range3(2,2)){
    if (k == range3(1,3))
      this_curve = edge_curve(10);
    else if (k == range3(2,3))
      this_curve = edge_curve(12);
  }

/* test all corners if this_curve == 0 to make sure a positive number is returned */
/* if any of the edges that join at the corner has a positive edge_curve value */
  if (this_curve == 0){
/* r3 = 0 */
    if (i == range3(1,1) && j == range3(1,2) && k == range3(1,3))
      this_curve = int_max( edge_curve(1), int_max( edge_curve(5), edge_curve(9) ) );
    else if (i == range3(2,1) && j == range3(1,2) && k == range3(1,3))
      this_curve = int_max( edge_curve(2), int_max( edge_curve(6), edge_curve(9) ) );
    else if (i == range3(1,1) && j == range3(2,2) && k == range3(1,3))
      this_curve = int_max( edge_curve(3), int_max( edge_curve(5), edge_curve(10) ) );
    else if (i == range3(2,1) && j == range3(2,2) && k == range3(1,3))
      this_curve = int_max( edge_curve(4), int_max( edge_curve(6), edge_curve(10) ) );
/* r3 = 1 */
    else if (i == range3(1,1) && j == range3(1,2) && k == range3(2,3))
      this_curve = int_max( edge_curve(1), int_max( edge_curve(7), edge_curve(11) ) );
    else if (i == range3(2,1) && j == range3(1,2) && k == range3(2,3))
      this_curve = int_max( edge_curve(2), int_max( edge_curve(8), edge_curve(11) ) );
    else if (i == range3(1,1) && j == range3(2,2) && k == range3(2,3))
      this_curve = int_max( edge_curve(3), int_max( edge_curve(7), edge_curve(12) ) );
    else if (i == range3(2,1) && j == range3(2,2) && k == range3(2,3))
      this_curve = int_max( edge_curve(4), int_max( edge_curve(8), edge_curve(12) ) );
  }

  return this_curve;
}

static int
interior_in_i(int i, int j, int k, component_vgrid *vgrid, 
	      int dw2, int nw1, int tw2){
  int ip, jp, kp, disc_pnt=1;

/* interior in the j-direction */
  if (j >= range3(1,2)+dw2 && j <= range3(2,2)-dw2){
/* interior in the k-direction */
    if (k >= range3(1,3)+dw2 && k <= range3(2,3)-dw2){
      for (ip=i - dw2; ip<= i + dw2; ip++)
	for (jp=j - dw2; jp<= j + dw2; jp++)
	  for (kp=k - dw2; kp<= k + dw2; kp++)
	    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
    }
    /* close to near side (k=range3(1,3)), away from edges and corners */
    else if (k < range3(1,3)+dw2){
      if (surf3(1,3) == 0)
	disc_pnt = 0;
      else{
	for (ip=i-tw2; ip<=i+tw2; ip++)
	  for (jp=j-tw2; jp<=j+tw2; jp++)
	    for (kp=range3(1,3); kp<=range3(1,3)+nw1; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      }
    }
    /* close to far side (k=range3(2,3)), away from edges and corners */
    else if (k > range3(2,3)-dw2){
      if (surf3(2,3) == 0)
	disc_pnt = 0;
      else{
	for (ip=i-tw2; ip<=i+tw2; ip++)
	  for (jp=j-tw2; jp<=j+tw2; jp++)
	    for (kp=range3(2,3)-nw1; kp<=range3(2,3); kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      }
    }
    /* this should never happen */
    else
      disc_pnt = 0;
  } /* end interior j */
  /* close to lower side (j=range3(1,2)) away from corners */
  else if (j < range3(1,2)+dw2){
    if ( surf3(1,2) == 0 )
      disc_pnt = 0;
    else{
      /* close to the lower side (j=range3(1,2) and interior in the k-direction, */
      /* away from edges and corners */
      if (k >= range3(1,3)+dw2 && k <= range3(2,3)-dw2){
	for (ip=i - tw2; ip<= i + tw2; ip++)
	  for (jp=range3(1,2); jp<= range3(1,2)+nw1; jp++)
	    for (kp=k - tw2; kp<= k + tw2; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      } /* end interior k */
      /* close to near side (k=range3(1,3)) and to lower side (j=range3(1,2), */
      /* away from corners */
      else if (k < range3(1,3)+dw2){
	if (surf3(1,3) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=i-tw2; ip<=i+tw2; ip++)
	    for (jp=range3(1,2); jp<= range3(1,2)+nw1; jp++)
	      for (kp=range3(1,3); kp<=range3(1,3)+nw1; kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on near side */
      } /* end k close to near side */
      /* close to far side (k=range3(2,3)) and to lower side (j=range3(1,2), */
      /* away from corners */
      else if (k > range3(2,3)-dw2){
	if (surf3(2,3) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=i-tw2; ip<=i+tw2; ip++)
	    for (jp=range3(1,2); jp<= range3(1,2)+nw1; jp++)
	      for (kp=range3(2,3)-nw1; kp<=range3(2,3); kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on far side */
      } /* end k close to far side */
    } /* end non-zero bc on lower side */
  } /* end j close to lower side */
  /* close to upper side (j=range3(2,2)) away from corners */
  else if (j > range3(2,2)-dw2){
    if ( surf3(2,2) == 0 )
      disc_pnt = 0;
    else{
      /* close to the upper side (j=range3(2,2) and interior in the k-direction, */
      /* away from edges and corners */
      if (k >= range3(1,3)+dw2 && k <= range3(2,3)-dw2){
	for (ip=i - tw2; ip<= i + tw2; ip++)
	  for (jp=range3(2,2)-nw1; jp<= range3(2,2); jp++)
	    for (kp=k - tw2; kp<= k + tw2; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      } /* end interior k */
      /* close to the near side (k=range3(1,3)) and to the upper side (j=range3(2,2), */
      /* away from corners */
      else if (k < range3(1,3)+dw2){
	if (surf3(1,3) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=i-tw2; ip<=i+tw2; ip++)
	    for (jp=range3(2,2)-nw1; jp<= range3(2,2); jp++)
	      for (kp=range3(1,3); kp<=range3(1,3)+nw1; kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on near side */
      } /* end k close to near side */
      /* close to the far side (k=range3(2,3)) and to the upper side (j=range3(2,2), */
      /* away from corners */
      else if (k > range3(2,3)-dw2){
	if (surf3(2,3) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=i-tw2; ip<=i+tw2; ip++)
	    for (jp=range3(2,2)-nw1; jp<= range3(2,2); jp++)
	      for (kp=range3(2,3)-nw1; kp<=range3(2,3); kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on far side*/
      } /* end k close to far side */
    } /* end non-zero bc on upper side */
  } /* end j close to upper side */

  return disc_pnt;
}

static int
interior_in_j(int i, int j, int k, component_vgrid *vgrid, 
	      int dw2, int nw1, int tw2){
  int ip, jp, kp, disc_pnt=1;

/* interior in the i-direction */
  if (i >= range3(1,1)+dw2 && i <= range3(2,1)-dw2){
/* interior in the k-direction */
    if (k >= range3(1,3)+dw2 && k <= range3(2,3)-dw2){
      for (ip=i - dw2; ip<= i + dw2; ip++)
	for (jp=j - dw2; jp<= j + dw2; jp++)
	  for (kp=k - dw2; kp<= k + dw2; kp++)
	    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
    }
    /* close to near side (k=range3(1,3)), away from edges and corners */
    else if (k < range3(1,3)+dw2){
      if (surf3(1,3) == 0)
	disc_pnt = 0;
      else{
	for (ip=i-tw2; ip<=i+tw2; ip++)
	  for (jp=j-tw2; jp<=j+tw2; jp++)
	    for (kp=range3(1,3); kp<=range3(1,3)+nw1; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      }
    }
    /* close to far side (k=range3(2,3)), away from edges and corners */
    else if (k > range3(2,3)-dw2){
      if (surf3(2,3) == 0)
	disc_pnt = 0;
      else{
	for (ip=i-tw2; ip<=i+tw2; ip++)
	  for (jp=j-tw2; jp<=j+tw2; jp++)
	    for (kp=range3(2,3)-nw1; kp<=range3(2,3); kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      }
    }
/* this should never happen */
    else
      disc_pnt = 0;
  } /* end interior i */
/* close to left side (i=range3(1,1)) away from corners */
  else if (i < range3(1,1)+dw2){
    if ( surf3(1,1) == 0 )
      disc_pnt = 0;
    else{
/* close to the left side (i=range3(1,1) and interior in the k-direction, */
/* away from edges and corners */
      if (k >= range3(1,3)+dw2 && k <= range3(2,3)-dw2){
	for (ip=range3(1,1); ip<= range3(1,1)+nw1; ip++)
	  for (jp= j - tw2; jp<= j + tw2; jp++)
	    for (kp= k - tw2; kp<= k + tw2; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      } /* end interior k */
/* close to near side (k=range3(1,3)) and to left side (i=range3(1,1), */
/* away from corners */
      else if (k < range3(1,3)+dw2){
	if (surf3(1,3) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=range3(1,1); ip<= range3(1,1)+nw1; ip++)
	    for (jp=j-tw2; jp<=j+tw2; jp++)
	      for (kp=range3(1,3); kp<=range3(1,3)+nw1; kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on near side */
      } /* end k close to near side */
/* close to far side (k=range3(2,3)) and to left side (i=range3(1,1), */
/* away from corners */
      else if (k > range3(2,3)-dw2){
	if (surf3(2,3) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=range3(1,1); ip<= range3(1,1)+nw1; ip++)
	    for (jp=j-tw2; jp<=j+tw2; jp++)
	      for (kp=range3(2,3)-nw1; kp<=range3(2,3); kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on far side */
      } /* end k close to far side */
    } /* end non-zero bc on left side */
  } /* end i close to left side */

/* close to right side (i=range3(2,1)) away from corners */
  else if (i > range3(2,1)-dw2){
    if ( surf3(2,1) == 0 )
      disc_pnt = 0;
    else{
/* close to the right side (i=range3(2,1) and interior in the k-direction, */
/* away from edges and corners */
      if (k >= range3(1,3)+dw2 && k <= range3(2,3)-dw2){
	for (jp=j - tw2; jp<=j + tw2; jp++)
	  for (ip=range3(2,1) - nw1; ip<=range3(2,1); ip++)
	    for (kp=k - tw2; kp<= k + tw2; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      } /* end interior k */
/* close to the near side (k=range3(1,3)) and to the right side (j=range3(2,2), */
/* away from corners */
      else if (k < range3(1,3)+dw2){
	if (surf3(1,3) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=range3(2,1) - nw1; ip<=range3(2,1); ip++)
	    for (jp=j - tw2; jp<=j + tw2; jp++)
	      for (kp=range3(1,3); kp<=range3(1,3)+nw1; kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on near side */
      } /* end k close to near side */
      /* close to the far side (k=range3(2,3)) and to the right side (j=range3(2,2), */
      /* away from corners */
      else if (k > range3(2,3)-dw2){
	if (surf3(2,3) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=range3(2,1) - nw1; ip<=range3(2,1); ip++)
	    for (jp=j - tw2; jp<=j + tw2; jp++)
	      for (kp=range3(2,3)-nw1; kp<=range3(2,3); kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on far side*/
      } /* end k close to far side */
    } /* end non-zero bc on right side */
  } /* end i close to right side */

  return disc_pnt;
}

static int
interior_in_k(int i, int j, int k, component_vgrid *vgrid, 
	      int dw2, int nw1, int tw2){
  int ip, jp, kp, disc_pnt=1;

/* interior in the j-direction */
  if (j >= range3(1,2)+dw2 && j <= range3(2,2)-dw2){
/* interior in the i-direction */
    if (i >= range3(1,1)+dw2 && i <= range3(2,1)-dw2){
      for (ip=i - dw2; ip<= i + dw2; ip++)
	for (jp=j - dw2; jp<= j + dw2; jp++)
	  for (kp=k - dw2; kp<= k + dw2; kp++)
	    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
    }
/* close to the left side (i=range3(1,1)), away from edges and corners */
    else if (i < range3(1,1)+dw2){
      if (surf3(1,1) == 0)
	disc_pnt = 0;
      else{
	for (ip=range3(1,1); ip<=range3(1,1) + nw1; ip++)
	  for (jp=j-tw2; jp<=j+tw2; jp++)
	    for (kp=k-tw2; kp<=k-tw2; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      }
    }
/* close to right side (i=range3(2,1)), away from edges and corners */
    else if (i > range3(2,1)-dw2){
      if (surf3(2,1) == 0)
	disc_pnt = 0;
      else{
	for (ip=range3(2,1)-nw1; ip<=range3(2,1); ip++)
	  for (jp=j-tw2; jp<=j+tw2; jp++)
	    for (kp=k-tw2; kp<=k-tw2; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      }
    }
/* this should never happen */
    else
      disc_pnt = 0;
  } /* end interior j */
/* close to lower side (j=range3(1,2)) away from corners */
  else if (j < range3(1,2)+dw2){
    if ( surf3(1,2) == 0 )
      disc_pnt = 0;
    else{
/* close to the lower side (j=range3(1,2) and interior in the i-direction, */
/* away from edges and corners */
      if (i >= range3(1,1)+dw2 && i <= range3(2,1)-dw2){
	for (ip=i-tw2; ip<=i+tw2; ip++)
	  for (jp=range3(1,2); jp<= range3(1,2)+nw1; jp++)
	    for (kp=k-tw2; kp<=k-tw2; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      } /* end interior i */
/* close to the left side (i=range3(1,1)) and to lower side (j=range3(1,2), */
/* away from corners */
      else if (i < range3(1,1)+dw2){
	if (surf3(1,1) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=range3(1,1); ip<=range3(1,1) + nw1; ip++)
	    for (jp=range3(1,2); jp<= range3(1,2)+nw1; jp++)
	      for (kp=k-tw2; kp<=k-tw2; kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on near side */
      } /* end i close to left side */
/* close to right side (i=range3(2,1)) and to lower side (j=range3(1,2), */
/* away from corners */
      else if (i > range3(2,1)-dw2){
	if (surf3(2,1) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=range3(2,1)-nw1; ip<=range3(2,1); ip++)
	    for (jp=range3(1,2); jp<= range3(1,2)+nw1; jp++)
	      for (kp=k-tw2; kp<=k-tw2; kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on far side */
      } /* end k close to far side */
    } /* end non-zero bc on lower side */
  } /* end j close to lower side */
  /* close to upper side (j=range3(2,2)) away from corners */
  else if (j > range3(2,2)-dw2){
    if ( surf3(2,2) == 0 )
      disc_pnt = 0;
    else{
/* close to the upper side (j=range3(2,2) and interior in the i-direction, */
/* away from edges and corners */
      if (i >= range3(1,1)+dw2 && i <= range3(2,1)-dw2){
	for (ip=i - tw2; ip<= i + tw2; ip++)
	  for (jp=range3(2,2)-nw1; jp<= range3(2,2); jp++)
	    for (kp=k - tw2; kp<= k + tw2; kp++)
	      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
      } /* end interior k */
/* close to the left side (i=range3(1,1)) and to the upper side (j=range3(2,2), */
/* away from corners */
      else if (i < range3(1,1)+dw2){
	if (surf3(1,1) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=range3(1,1); ip<=range3(1,1)+nw1; ip++)
	    for (jp=range3(2,2)-nw1; jp<= range3(2,2); jp++)
	      for (kp=k-tw2; kp<=k-tw2; kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on near side */
      } /* end k close to near side */
/* close to the right side (i=range3(2,1)) and to the upper side (j=range3(2,2), */
/* away from corners */
      else if (i > range3(2,1)-dw2){
	if (surf3(2,1) == 0)
	  disc_pnt = 0;
	else{
	  for (ip=range3(2,1)-nw1; ip<=range3(2,1); ip++)
	    for (jp=range3(2,2)-nw1; jp<= range3(2,2); jp++)
	      for (kp=k-tw2; kp<=k-tw2; kp++)
		if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	} /* end non-zero bc on far side*/
      } /* end k close to far side */
    } /* end non-zero bc on upper side */
  } /* end j close to upper side */

  return disc_pnt;
}

static int
disc_point(int i, int j, int k, component_vgrid *vgrid, overlapping_3d_grid *over3d){
  int disc_pnt, ip, jp, kp, dw2, cw1, tw2, nw1;

/* The point (i,j,k) is a valid discretization point if sufficiently many of its
   neighbors are either discretization points or interpolation points */

  dw2 = (over3d->disc_width -1) / 2;
  nw1 = over3d->normal_width - 1;
  cw1 = over3d->corner_width - 1;
  tw2 = (over3d->tangent_width -1) / 2;
  if (tw2 > dw2){
    printf("Disc_point> Warning: tangential width > discretization width "
	   "not implemented\n");
  }

/* discretization point until otherwise proven */
  disc_pnt = 1;

/* 
   PERIODIC IN THE R1-DIRECTION 
*/
  if (vgrid->r1_period){
/* shift `i' if necessary */
    if (i < range3(1,1))
      i += range3(2,1) - range3(1,1);
    else if (i > range3(2,1)-1)
      i += -(range3(2,1) - range3(1,1));
    disc_pnt = interior_in_i(i, j, k, vgrid, dw2, nw1, tw2);
  }
/* 
   PERIODIC IN THE R2-DIRECTION 
*/
  else if (vgrid->r2_period){
/* shift `j' if necessary */
    if (j < range3(1,2))
      j += range3(2,2) - range3(1,2);
    else if (j > range3(2,2)-1)
      j += -(range3(2,2) - range3(1,2));
    disc_pnt = interior_in_j(i, j, k, vgrid, dw2, nw1, tw2);
  }
/* 
   PERIODIC IN THE R3-DIRECTION 
*/
  else if (vgrid->r3_period){
/* shift `k' if necessary */
    if (k < range3(1,3))
      k += range3(2,3) - range3(1,3);
    else if (k > range3(2,3)-1)
      k += -(range3(2,3) - range3(1,3));
    disc_pnt = interior_in_k(i, j, k, vgrid, dw2, nw1, tw2);
  }
/*
   NON-PERIODIC CASE 
*/
  else{
/* interior in the i-direction */
    if (i >= range3(1,1)+dw2 && i <= range3(2,1)-dw2){
      disc_pnt = interior_in_i(i, j, k, vgrid, dw2, nw1, tw2);
    }
/* close to left side (i=range3(1,1)) */
    else if (i < range3(1,1)+dw2){
      if ( surf3(1,1) == 0 )
	disc_pnt = 0;
      else{
/* interior in j */
	if (j>= range3(1,2)+dw2 && j<= range3(2,2)-dw2){
/* interior in k */
	  if (k>= range3(1,3)+dw2 && k<= range3(2,3)-dw2){
	    for (ip=range3(1,1); ip<= range3(1,1)+nw1; ip++)
	      for (jp=j - tw2; jp<= j + tw2; jp++)
		for (kp=k - tw2; kp<= k + tw2; kp++)
		  if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	  } /* end interior in k */
/* close to near (k=range3(1,3)) side, i.e. left-near edge */
	  else if (k< range3(1,3)+dw2){
	    if (surf3(1,3) == 0)
	      disc_pnt = 0;
	    else{
	      for (ip=range3(1,1); ip<= range3(1,1)+nw1; ip++)
		for (jp=j - tw2; jp<= j + tw2; jp++)
		  for (kp=range3(1,3); kp<=range3(1,3) + nw1; kp++)
		    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	    } /* end non-zero bc on near side */
	  } /* end close to near (k=range3(1,3)) side */
/* close to far (k=range3(2,3)) side, i.e. left-far edge */
	  else if (k> range3(2,3)-dw2){
	    if (surf3(2,3) == 0)
	      disc_pnt = 0;
	    else{
	      for (ip=range3(1,1); ip<= range3(1,1)+nw1; ip++)
		for (jp=j - tw2; jp<= j + tw2; jp++)
		  for (kp=range3(2,3)-nw1; kp<=range3(2,3); kp++)
		    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	    } /* end non-zero bc on far side */
	  } /* end close to far (k=range3(2,3)) side */
	} /* end interior in j */
/* lower left edge */
	else if (j < range3(1,2)+dw2){
	  if ( surf3(1,2) == 0 )
	    disc_pnt = 0;
	  else{
/* interior in k, i.e. left-lower edge */
	    if (k>= range3(1,3)+dw2 && k<= range3(2,3)-dw2){
	      for (ip=range3(1,1); ip<= range3(1,1)+nw1; ip++)
		for (jp=range3(1,2); jp<=range3(1,2) + nw1; jp++)
		  for (kp=k - tw2; kp<= k + tw2; kp++)
		    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	    } /* end interior in k */
/* close to near (k=range3(1,3)) side, i.e. left-lower-near corner */
	    else if (k< range3(1,3)+dw2){
	      if (surf3(1,3) == 0)
		disc_pnt = 0;
	      else{
		for (ip=range3(1,1); ip<= range3(1,1)+cw1; ip++)
		  for (jp=range3(1,2); jp<= range3(1,2) + cw1; jp++)
		    for (kp=range3(1,3); kp<=range3(1,3) + cw1; kp++)
		      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	      } /* end non-zero bc on near side */
	    } /* end close to near (k=range3(1,3)) side */
/* close to far (k=range3(2,3)) side, i.e. left-lower-far corner */
	    else if (k> range3(2,3)-dw2){
	      if (surf3(2,3) == 0)
		disc_pnt = 0;
	      else{
		for (ip=range3(1,1); ip<= range3(1,1)+cw1; ip++)
		  for (jp= range3(1,2); jp<= range3(1,2) + cw1; jp++)
		    for (kp=range3(2,3)-cw1; kp<=range3(2,3); kp++)
		      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	      } /* end non-zero bc on far side */
	    } /* end close to lower-left-far corner */
	  } /* end non-zero bc on lower side */
	} /* end lower left edge */
/* upper left edge */
	else if (j > range3(2,2)-dw2){
	  if ( surf3(2,2) == 0 )
	    disc_pnt = 0;
	  else{
/* interior in k, i.e. left-upper edge */
	    if (k>= range3(1,3)+dw2 && k<= range3(2,3)-dw2){
	      for (ip=range3(1,1); ip<= range3(1,1)+nw1; ip++)
		for (jp=range3(2,2)-nw1; jp<=range3(2,2); jp++)
		  for (kp=k - tw2; kp<= k + tw2; kp++)
		    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	    } /* end interior in k */
/* close to near (k=range3(1,3)) side, i.e. left-upper-near corner */
	    else if (k< range3(1,3)+dw2){
	      if (surf3(1,3) == 0)
		disc_pnt = 0;
	      else{
		for (ip=range3(1,1); ip<= range3(1,1)+cw1; ip++)
		  for (jp=range3(2,2)-cw1; jp<= range3(2,2); jp++)
		    for (kp=range3(1,3); kp<=range3(1,3) + cw1; kp++)
		      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	      } /* end non-zero bc on near side */
	    } /* end close to near (k=range3(1,3)) side */
/* close to far (k=range3(2,3)) side, i.e. left-upper-far corner */
	    else if (k> range3(2,3)-dw2){
	      if (surf3(2,3) == 0)
		disc_pnt = 0;
	      else{
		for (ip=range3(1,1); ip<= range3(1,1)+cw1; ip++)
		  for (jp= range3(2,2)-cw1; jp<= range3(2,2); jp++)
		    for (kp=range3(2,3)-cw1; kp<=range3(2,3); kp++)
		      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	      } /* end non-zero bc on far side */
	    } /* end close to upper-left-far corner */
	  } /* end non-zero bc on upper side */
	} /* end upper left edge */
/* this should never happen */
	else
	  disc_pnt = 0;

      } /* end non-zero bc on left side */
    } /* end close to left side (i=range3(1,1)) */

/* close to right side (i=range3(2,1)) */
    else if (i > range3(2,1)-dw2){
      if ( surf3(2,1) == 0 )
	disc_pnt = 0;
      else{
/* interior in j */
	if (j>= range3(1,2)+dw2 && j<= range3(2,2)-dw2){
/* interior in k */
	  if (k>= range3(1,3)+dw2 && k<= range3(2,3)-dw2){
	    for (ip=range3(2,1)-nw1; ip<= range3(2,1); ip++)
	      for (jp=j - tw2; jp<= j + tw2; jp++)
		for (kp=k - tw2; kp<= k + tw2; kp++)
		  if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	  } /* end interior in k */
/* close to near (k=range3(1,3)) side, i.e. right-near edge */
	  else if (k< range3(1,3)+dw2){
	    if (surf3(1,3) == 0)
	      disc_pnt = 0;
	    else{
	      for (ip=range3(2,1)-nw1; ip<= range3(2,1); ip++)
		for (jp=j - tw2; jp<= j + tw2; jp++)
		  for (kp=range3(1,3); kp<=range3(1,3) + nw1; kp++)
		    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	    } /* end non-zero bc on near side */
	  } /* end close to near (k=range3(1,3)) side */
/* close to far (k=range3(2,3)) side, i.e. right-far edge */
	  else if (k> range3(2,3)-dw2){
	    if (surf3(2,3) == 0)
	      disc_pnt = 0;
	    else{
	      for (ip=range3(2,1)-nw1; ip<= range3(2,1); ip++)
		for (jp=j - tw2; jp<= j + tw2; jp++)
		  for (kp=range3(2,3)-nw1; kp<=range3(2,3); kp++)
		    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	    } /* end non-zero bc on far side */
	  } /* end close to far (k=range3(2,3)) side */
	} /* end interior in j */
/* lower right edge */
	else if (j < range3(1,2)+dw2){
	  if ( surf3(1,2) == 0 )
	    disc_pnt = 0;
	  else{
/* interior in k, i.e. right-lower edge */
	    if (k>= range3(1,3)+dw2 && k<= range3(2,3)-dw2){
	      for (ip=range3(2,1)-nw1; ip<= range3(2,1); ip++)
		for (jp=range3(1,2); jp<=range3(1,2) + nw1; jp++)
		  for (kp=k - tw2; kp<= k + tw2; kp++)
		    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	    } /* end interior in k */
/* close to near (k=range3(1,3)) side, i.e. right-lower-near corner */
	    else if (k< range3(1,3)+dw2){
	      if (surf3(1,3) == 0)
		disc_pnt = 0;
	      else{
		for (ip=range3(2,1)-cw1; ip<= range3(2,1); ip++)
		  for (jp=range3(1,2); jp<= range3(1,2) + cw1; jp++)
		    for (kp=range3(1,3); kp<=range3(1,3) + cw1; kp++)
		      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	      } /* end non-zero bc on near side */
	    } /* end close to near (k=range3(1,3)) side */
/* close to far (k=range3(2,3)) side, i.e. right-lower-far corner */
	    else if (k> range3(2,3)-dw2){
	      if (surf3(2,3) == 0)
		disc_pnt = 0;
	      else{
		for (ip=range3(2,1)-cw1; ip<= range3(2,1); ip++)
		  for (jp= range3(1,2); jp<= range3(1,2) + cw1; jp++)
		    for (kp=range3(2,3)-cw1; kp<=range3(2,3); kp++)
		      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	      } /* end non-zero bc on far side */
	    } /* end close to lower-right-far corner */
	  } /* end non-zero bc on lower side */
	} /* end lower right edge */
/* upper right edge */
	else if (j > range3(2,2)-dw2){
	  if ( surf3(2,2) == 0 )
	    disc_pnt = 0;
	  else{
/* interior in k, i.e. right-upper edge */
	    if (k>= range3(1,3)+dw2 && k<= range3(2,3)-dw2){
	      for (ip=range3(2,1)-nw1; ip<= range3(2,1); ip++)
		for (jp=range3(2,2)-nw1; jp<=range3(2,2); jp++)
		  for (kp=k - tw2; kp<= k + tw2; kp++)
		    if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	    } /* end interior in k */
/* close to near (k=range3(1,3)) side, i.e. right-upper-near corner */
	    else if (k< range3(1,3)+dw2){
	      if (surf3(1,3) == 0)
		disc_pnt = 0;
	      else{
		for (ip=range3(2,1)-cw1; ip<= range3(2,1); ip++)
		  for (jp=range3(2,2)-cw1; jp<= range3(2,2); jp++)
		    for (kp=range3(1,3); kp<=range3(1,3) + cw1; kp++)
		      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	      } /* end non-zero bc on near side */
	    } /* end close to near (k=range3(1,3)) side */
/* close to far (k=range3(2,3)) side, i.e. right-upper-far corner */
	    else if (k> range3(2,3)-dw2){
	      if (surf3(2,3) == 0)
		disc_pnt = 0;
	      else{
		for (ip=range3(2,1)-cw1; ip<= range3(2,1); ip++)
		  for (jp= range3(2,2)-cw1; jp<= range3(2,2); jp++)
		    for (kp=range3(2,3)-cw1; kp<=range3(2,3); kp++)
		      if (flag3(ip,jp,kp) == 0) disc_pnt = 0;
	      } /* end non-zero bc on far side */
	    } /* end close to upper-right-far corner */
	  } /* end non-zero bc on upper side */
	} /* end upper right edge */
/* this should never happen */
	else
	  disc_pnt = 0;

      } /* end non-zero bc on right side */
    } /* end close to right side (i=range3(1,1)) */
  
/* this should never happen */
    else
      disc_pnt = 0;

  }    
/* end non-periodic case */

  return disc_pnt;
}

static interp_3d_point *
interp_from_higher(real xp, real yp, real zp, 
		   int i_point, int j_point, int k_point,
		   interp_3d_point *initial_guess, int query_surface, int query_curve,
		   linked_list_member *other_link, overlapping_3d_grid *over3d){
  component_vgrid *vgrid;
  interp_3d_point *new_interp3=NULL;
  linked_list_member *vgrid_link;

  for (vgrid_link = over3d->grid_list->first; 
       vgrid_link != NULL && vgrid_link != other_link; 
       vgrid_link = vgrid_link->next){
    vgrid = vgrid_link->data;

    if ((new_interp3 = invert_grid(xp, yp, zp, vgrid, query_surface, query_curve, 
				   initial_guess)) 
	!= NULL){ 
/* save the indices of the interpolation point */
      new_interp3->i_point = i_point;
      new_interp3->j_point = j_point;
      new_interp3->k_point = k_point;
      if (good_interp_loc(new_interp3, vgrid, over3d))
	return new_interp3;
      else{
/* free the space */
	free(new_interp3);
	new_interp3 = NULL;
      }
    } /* end if new_interp3 ... */
  } /* end for */

/* no match, return NULL */
  return NULL;
}

static interp_3d_point *
interp_from_lower(real xp, real yp, real zp, 
		  int i_point, int j_point, int k_point,
		  interp_3d_point *initial_guess, int query_surface, int query_curve,
		  component_vgrid *own_vgrid, 
		  int lower_than, overlapping_3d_grid *over3d){
  component_vgrid *vgrid;
  interp_3d_point *new_interp3=NULL;
  linked_list_member *vgrid_link, *other_link;
  
  for (other_link = over3d->grid_list->first; other_link != NULL;
       other_link = other_link->next){
    if (lower_than == ((component_vgrid *) other_link->data)->priority) break;
  }

  for (vgrid_link = other_link->next; 
       vgrid_link != NULL;
       vgrid_link = vgrid_link->next){
    vgrid = vgrid_link->data;

    if (vgrid != own_vgrid){

/* it is always possible to interpolate from a dummy grid */
/*       if (over3d->overlapping_sub_grid && vgrid->dummy_background){ */
/* 	new_interp3 = new_3d_interp( 0, 0, 0, vgrid, -1.0, -1.0, -1.0 ); */
/* save the indices of the interpolation point */
/* 	new_interp3->i_point = i_point; */
/* 	new_interp3->j_point = j_point; */
/* 	new_interp3->k_point = k_point; */
/* 	return new_interp3; */
/*       } */

      if ((new_interp3 = invert_grid(xp, yp, zp, vgrid, query_surface, query_curve, 
				     initial_guess)) 
	  != NULL){
/* save the indices of the interpolation point */
	new_interp3->i_point = i_point;
	new_interp3->j_point = j_point;
	new_interp3->k_point = k_point;
/* check if (i,j) is a valid interpolation location*/
	if (good_interp_loc(new_interp3, vgrid, over3d))
	  return new_interp3;
	else{
/* free the space */
	  free(new_interp3);
	  new_interp3 = NULL;
	}
      } /* end if new_interp3 ... */
    } /* end if vgrid != own_grid , end for vgrid ... */
  } /* end for */

/* no match, return NULL */
  return NULL;
}

int 
good_interp_loc(interp_3d_point *new_interp3,
		component_vgrid *vgrid,
		overlapping_3d_grid *over3d){
  int iw1, iw_low, iw_high, i_min, i_max, j_min, j_max, ip, jp, i, j
    , k, k_min, k_max, kp, ok=TRUE;

/* first check that the point is at least 0.5 grid cells away from all interpolation
   boundaries. */
  if (!vgrid->r1_period){
/* explicit */
    if (over3d->interp_type == 'e'){
      if ((new_interp3->r_loc < 0.5*over3d->interp_width * vgrid->r1_step && 
	   surf3(1,1) == 0 ) ||
	  (new_interp3->r_loc > 1.0 - 0.5*over3d->interp_width * vgrid->r1_step && 
	   surf3(2,1) == 0 ))
	ok = FALSE;
/* 	return 0; */
    }
/* implicit */
    else{
      if ((new_interp3->r_loc < 0.5*vgrid->r1_step && surf3(1,1) == 0 ) ||
	  (new_interp3->r_loc > 1.0 - 0.5*vgrid->r1_step && surf3(2,1) == 0))
	ok = FALSE;
/* 	return 0; */
    }
  }
  if (!vgrid->r2_period){
/* explicit */
    if (over3d->interp_type == 'e'){
      if ((new_interp3->s_loc < 0.5*over3d->interp_width * vgrid->r2_step && 
	   surf3(1,2) == 0 ) ||
	  (new_interp3->s_loc > 1.0 - 0.5*over3d->interp_width * vgrid->r2_step && 
	   surf3(2,2) == 0 ))
	ok = FALSE;
/* 	return 0; */
    }
/* implicit */
    else{
      if ((new_interp3->s_loc < 0.5*vgrid->r2_step && surf3(1,2) == 0) ||
	  (new_interp3->s_loc > 1.0 - 0.5*vgrid->r2_step && surf3(2,2) == 0))
	ok = FALSE;
/* 	return 0; */
    }
  }
  if (!vgrid->r3_period){
/* explicit */
    if (over3d->interp_type == 'e'){
      if ((new_interp3->t_loc < 0.5*over3d->interp_width * vgrid->r3_step && 
	   surf3(1,3) == 0 ) ||
	  (new_interp3->t_loc > 1.0 - 0.5*over3d->interp_width * vgrid->r3_step && 
	   surf3(2,3) == 0 ))
	ok = FALSE;
/* 	return 0; */
    }
/* implicit */
    else{
      if ((new_interp3->t_loc < 0.5*vgrid->r3_step && surf3(1,3) == 0) ||
	  (new_interp3->t_loc > 1.0 - 0.5*vgrid->r3_step && surf3(2,3) == 0))
	ok = FALSE;
/* 	return 0; */
    }
  }

  iw1 =  over3d->interp_width - 1;
  iw_low = (over3d->interp_width - 1)/2;
  iw_high = over3d->interp_width/2;

/* even interpolation width */
  if (over3d->interp_width % 2 == 0){

/* get the enclosing cell */
    i = range3(1,1) + new_interp3->r_loc * (range3(2,1) - range3(1,1));
    j = range3(1,2) + new_interp3->s_loc * (range3(2,2) - range3(1,2));
    k = range3(1,3) + new_interp3->t_loc * (range3(2,3) - range3(1,3));

/* make sure the closest cell is inside the grid */
    i = int_min( int_max( i, range3(1,1) ), range3(2,1)-1 );
    j = int_min( int_max( j, range3(1,2) ), range3(2,2)-1 );
    k = int_min( int_max( k, range3(1,3) ), range3(2,3)-1 );
  }

/* odd interpolation width */
  else{

/* get the closest point to (r, s, t) */
    i = range3(1,1) + new_interp3->r_loc * (range3(2,1) - range3(1,1)) + 0.5;
    j = range3(1,2) + new_interp3->s_loc * (range3(2,2) - range3(1,2)) + 0.5;
    k = range3(1,3) + new_interp3->t_loc * (range3(2,3) - range3(1,3)) + 0.5;

/* make sure the closest point is inside the grid */
    i = int_min( int_max( i, range3(1,1) ), range3(2,1) );
    j = int_min( int_max( j, range3(1,2) ), range3(2,2) );
    k = int_min( int_max( k, range3(1,3) ), range3(2,3) );
  }

/* periodic in r1 */
  if (vgrid->r1_period){
    if (i < range3(1,1)){
      i += range3(2,1) - range3(1,1);
      new_interp3->r_loc += 1.0;
    }
    else if (i > range3(2,1)-1){
      i += -(range3(2,1) - range3(1,1));
      new_interp3->r_loc += -1.0;
    }
    i_min = i - iw_low;
  }

/* non-periodic in r1 */
  else {

/* check if (i,j,k) is close to a boundary */
/* i-direction */
    if (i >= range3(1,1)+iw_low && i <= range3(2,1)-iw_high){
      i_min = i - iw_low; 
    }
    else if (i < range3(1,1)+iw_low){
      i_min = range3(1,1); 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (surf3(1,1) == 0) ok = FALSE; /* return 0; */
    }
    else if (i > range3(2,1)-iw_high){
      i_min = range3(2,1)-iw1; 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (surf3(2,1) == 0) ok = FALSE; /* return 0; */
    }
/* this can't happen */
    else
      ok = FALSE;
/*       return 0; */
  }

/* periodic in r2 */
  if (vgrid->r2_period){
    if (j < range3(1,2)){
      j += range3(2,2) - range3(1,2);
      new_interp3->s_loc += 1.0;
    }
    else if (j > range3(2,2)-1){
      j += -(range3(2,2) - range3(1,2));
      new_interp3->s_loc += -1.0;
    }
    j_min = j - iw_low;
  }

/* non-periodic in r2 */
  else {
    if (j >= range3(1,2)+iw_low && j <= range3(2,2)-iw_high){
      j_min = j - iw_low; 
    }
    else if (j < range3(1,2)+iw_low){
      j_min = range3(1,2); 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (surf3(1,2) == 0) ok = FALSE; /* return 0; */
    }
    else if (j > range3(2,2)-iw_high){
      j_min = range3(2,2)-iw1; 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (surf3(2,2) == 0) ok = FALSE; /* return 0; */
    }
/* this can't happen */
    else
      ok = FALSE;
/*       return 0; */
  }

/* periodic in r3 */
  if (vgrid->r3_period){
    if (k < range3(1,3)){
      k += range3(2,3) - range3(1,3);
      new_interp3->t_loc += 1.0;
    }
    else if (k > range3(2,3)-1){
      k += -(range3(2,3) - range3(1,3));
      new_interp3->t_loc += -1.0;
    }
    k_min = k - iw_low;
  }

/* non-periodic in r3 */
  else {
    if (k >= range3(1,3)+iw_low && k <= range3(2,3)-iw_high){
      k_min = k - iw_low; 
    }
    else if (k < range3(1,3)+iw_low){
      k_min = range3(1,3); 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (surf3(1,3) == 0) ok = FALSE; /* return 0; */
    }
    else if (k > range3(2,3)-iw_high){
      k_min = range3(2,3)-iw1; 
/* can't have a skewed interpolation stencil close to an interpolation point */
      if (surf3(2,3) == 0) ok = FALSE; /* return 0; */
    }
/* this can't happen */
    else
      ok = FALSE;
/*       return 0; */
  }

/* save the interpolation location */
  new_interp3->i_loc = i_min;
  new_interp3->j_loc = j_min;
  new_interp3->k_loc = k_min;

  i_max = i_min + iw1;
  j_max = j_min + iw1;
  k_max = k_min + iw1;

/* determine if it is a good interpolation location. At this point it is
   not possible to assure explicit interpolation. */

  for (ip = i_min; ip <= i_max; ip++)
    for (jp = j_min; jp <= j_max; jp++)
      for (kp = k_min; kp <= k_max; kp++)
	if (flag3(ip,jp,kp) == 0)
	  ok = FALSE;
/* 	  return 0; */

/* the interpolation location was good if you got this far */

/* Return the status of the interpolation location */
  return ok;
}




static void 
show_overlap_param( overlapping_3d_grid *over3d ){
  printf("\n");
  printf("Corner width: %i\n",over3d->corner_width);
  printf("Discretization width: %i\n",over3d->disc_width);
  printf("Number of ghostpoints: %i\n",over3d->extra);
  printf("Interpolation width: %i\n",over3d->interp_width);
  printf("Normal width: %i\n",over3d->normal_width);
  printf("Tangential width: %i\n",over3d->tangent_width);
  
  printf("\n");
  printf("The interpolation type is ");
  if (over3d->interp_type == 'e')
    printf("explicit\n");
  else if (over3d->interp_type == 'i')
    printf("implicit\n");
  else
    printf("unknown\n");

  printf("\n");
  printf("Extra boundary missfit tolerance: %e\n", over3d->max_b_dist);

}
