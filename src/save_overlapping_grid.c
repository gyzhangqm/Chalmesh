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
#include "chalmesh.h"

#ifdef HDF_DUMMY
  #include "hdf_dummy.h"
#else
  #include "hdf_stuff.h"
#endif

static void
save_component_grid( component_vgrid *vgrid, int32 comp_dir );

int
save_overlapping_grid( input_output *io_ptr, overlapping_3d_grid *over3d ){
  char *name, comp_dir_name[80], creator[120];
  int32 root, over_dir, comp_dir;
  component_vgrid *vgrid;
  linked_list_member *vgrid_link;
  time_t today;
  real_array_2d *bbox_=NULL;
#define bbox(i,j) compute_index_2d(bbox_, i, j)

  name = get_word( io_ptr, "File name for overlapping grid: ", "test.hdf", 1);
  
  if ((root = open_hdf_file(name, 'i')) <= 0){
    printf("Unable to open the database file %s\n", name);
    return ERROR;
  }

  today = time(NULL);
  sprintf(creator, "chalmesh version "VERSION" generated this file on %s", 
	  ctime(&today));
  hput_string(creator, "creator", root);

/* make a directory for the overlapping grid */
  if ((over_dir = create_dir("overlapping grid", "overlapping grid", root)) != -1){

/* save all data for the overlapping grid */
#ifdef SINGLE
    hput_int(32, "precision", over_dir);
#else
    hput_int(64, "precision", over_dir);
#endif

    hput_int(3, "number of dimensions", over_dir);

    hput_int(over3d->grid_list->n_members, "n_components", over_dir);
    hput_string(over3d->name, "overlapping grid name", over_dir);
    hput_string((over3d->interp_type == 'e')? "explicit" : "implicit", 
		"interpolation type", over_dir);
    hput_int(over3d->interp_width, "interpolation width", over_dir);
    hput_int(over3d->disc_width, "discretization width", over_dir);
    hput_int(over3d->normal_width, "normal width", over_dir);
    hput_int(over3d->tangent_width, "tangent width", over_dir);
    hput_int(over3d->corner_width, "corner width", over_dir);
    hput_int(over3d->extra, "ghost points", over_dir);
    hput_int(over3d->extra_period, "periodic overlap", over_dir);

/* make a real array 2d to put the bounding box in */
    bbox_ = create_real_array_2d(2,3);
    bbox(1,1) = over3d->bb->x_min;
    bbox(2,1) = over3d->bb->x_max;
    bbox(1,2) = over3d->bb->y_min;
    bbox(2,2) = over3d->bb->y_max;
    bbox(1,3) = over3d->bb->z_min;
    bbox(2,3) = over3d->bb->z_max;
/* save the bounding box */
    hput_real_array_2d(bbox_, "bounding box", over_dir);

/* save all component grids */
    for (vgrid_link = over3d->grid_list->first; vgrid_link != NULL;
	 vgrid_link = vgrid_link->next){
      vgrid = vgrid_link->data;

      sprintf(comp_dir_name, "component grid %i", vgrid->priority);
      if ((comp_dir = create_dir(comp_dir_name, "component grid directory", 
				 over_dir)) != -1){
	save_component_grid(vgrid, comp_dir);
/* release the sub-directory */
	Vdetach(comp_dir);
      }
      else
	printf("Error: save_overlapping_grid: unable to create a directory for %s\n",
	       comp_dir_name);
    }
/* release the overlapping grid-directory */
    Vdetach(over_dir);
  }
  else
    printf("Error: save_overlapping_grid: unable to create a directory for "
	   "the overlapping grid\n");

/* close the database file */
  close_hdf_file(root);

/* deallocate the bounding box array */
  delete_real_array_2d(bbox_);

  return OK;
}


static void
save_component_grid( component_vgrid *vgrid, int32 comp_dir ){
  int_array_1d *dimension_, *periodicity_;
  real_array_1d *step_;
  real_array_2d *bbox_;
  real_array_2d *donor_par_;
  int_array_2d *donor_point_, *inter_point_;
  int i;
  interp_3d_point *interp3;
#define donor_par(i,j)   compute_index_2d(donor_par_, i, j)
#define donor_point(i,j) compute_index_2d(donor_point_, i, j)
#define inter_point(i,j) compute_index_2d(inter_point_, i, j)
#define dimension(i) compute_index_1d(dimension_, i)
#define periodicity(i) compute_index_1d(periodicity_, i)
#define step(i) compute_index_1d(step_, i)

/* make a real array 2d to put the bounding box in */
  bbox_ = create_real_array_2d(2,3);
  bbox(1,1) = vgrid->bb->x_min;
  bbox(2,1) = vgrid->bb->x_max;
  bbox(1,2) = vgrid->bb->y_min;
  bbox(2,2) = vgrid->bb->y_max;
  bbox(1,3) = vgrid->bb->z_min;
  bbox(2,3) = vgrid->bb->z_max;

  hput_string(vgrid->name, "component grid name", comp_dir);

  hput_string((vgrid->grid_type == 1)? "cartesian": "curvilinear", 
	      "grid type", comp_dir);

  dimension_ = create_int_array_1d(3);
  dimension(1) = vgrid->r1_dim;
  dimension(2) = vgrid->r2_dim;
  dimension(3) = vgrid->r3_dim;
  hput_int_array_1d(dimension_, "dimension", comp_dir);
  dimension_ = delete_int_array_1d(dimension_);

  periodicity_ = create_int_array_1d(3);
  periodicity(1) = vgrid->r1_period;
  periodicity(2) = vgrid->r2_period;
  periodicity(3) = vgrid->r3_period;
  hput_int_array_1d(periodicity_, "periodicity", comp_dir);
  periodicity_ = delete_int_array_1d(periodicity_);

  hput_int_array_2d(vgrid->range_ptr, "range", comp_dir);
  hput_int_array_2d(vgrid->bc_ptr, "boundary condition", comp_dir);
  hput_int_array_2d(vgrid->surf_ptr, "surface label", comp_dir);

  step_ = create_real_array_1d(3);
  step(1) = vgrid->r1_step;
  step(2) = vgrid->r2_step;
  step(3) = vgrid->r3_step;
  hput_real_array_1d(step_, "step", comp_dir);
  step_ = delete_real_array_1d(step_);

  hput_int(vgrid->priority, "priority", comp_dir);
  hput_real_array_2d(bbox_, "bounding box", comp_dir);
  hput_int_array_3d(vgrid->flag_ptr, "flag", comp_dir);
  hput_real_array_3d(vgrid->x_ptr, "x", comp_dir);
  hput_real_array_3d(vgrid->y_ptr, "y", comp_dir);
  hput_real_array_3d(vgrid->z_ptr, "z", comp_dir);
/* a discrete approximation of the jacobian is not implemented yet */
  hput_int(0, "save jacobian", comp_dir);

  hput_int(vgrid->n_interp, "n_interp", comp_dir);

  if (vgrid->n_interp > 0){
/* make 2-d arrays of the interpolation list */
    donor_par_ = create_real_array_2d(3, vgrid->n_interp);
    donor_point_ = create_int_array_2d(4, vgrid->n_interp);
    inter_point_ = create_int_array_2d(3, vgrid->n_interp);

    i = 0;
    for (interp3 = vgrid->last_interp; interp3 != NULL;
	 interp3 = interp3->prev){
      i++;
      
      donor_par(1, i)   = interp3->r_loc;
      donor_par(2, i)   = interp3->s_loc;
      donor_par(3, i)   = interp3->t_loc;

      donor_point(1, i) = interp3->i_loc;
      donor_point(2, i) = interp3->j_loc;
      donor_point(3, i) = interp3->k_loc;
      donor_point(4, i) = interp3->vgrid_loc->priority;

      inter_point(1, i) = interp3->i_point;
      inter_point(2, i) = interp3->j_point;
      inter_point(3, i) = interp3->k_point;
    }

/* save the interpolation information */
    hput_real_array_2d( donor_par_, "donor parameter", comp_dir );
    hput_int_array_2d( donor_point_, "donor point", comp_dir );
    hput_int_array_2d( inter_point_, "interpolation point", comp_dir );

/* release the temporary storage */
    delete_real_array_2d( donor_par_ );
    delete_int_array_2d( donor_point_ );
    delete_int_array_2d( inter_point_ );
  }

  delete_real_array_2d( bbox_ );
}
