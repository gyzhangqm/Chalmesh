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
#include "GL_GraphicsInterface.h"
#include "MappingInformation.h"
#include "HDF_DataBase.h"
#include "Overture.h"

/* the following function will be called from "C" */
extern "C"{
#include "volume.h"
#include "overmap.h"

/* static prototypes */
static int
overture_mapping( grid_point *gp_ptr, void *tuut );
static void 
set_overture_mapping(input_output *io_ptr, volume_mapping *volume, 
		      linked_list *surface_list );
static int
inverse_overture_mapping( grid_point *gp_ptr, void *tuut );
static void
delete_overture_data( void *volume_data_ptr );
/* end prototypes */

void
init_overture_maps(volume_mapping_list *volume_list, void *overture_list){
  MappingInformation *mappingInfo_ptr = (MappingInformation *)overture_list;
  grid_point gp;
  Bound b;
  volume_mapping *volume;
  volume_mapping_link *volume_link;
  const char *cname;
  char grid_name[120];
  int i, j;

/* loop through the list of existing mappings to see if any overture mappings */
/* have been deleted */
  for (volume_link = volume_list->first; volume_link != NULL; 
       volume_link = volume_link->next){
    volume = (volume_mapping *) volume_link->data;
/* check the type */
    if (strcmp(volume->type, "Overture mapping") == 0){
/* loop over all overture mappings in the list */
      for (i=0; i<mappingInfo_ptr->mappingList.getLength(); i++){
	MappingRC &map = mappingInfo_ptr->mappingList[i];
	if (volume->volume_data_ptr == &map) break;
      }
      if (i >= mappingInfo_ptr->mappingList.getLength()){
	printf("deleting the overture mapping `%s'\n", volume->name);
/* delete the volume mapping and remove it from the volume mapping list */
	if (delete_volume_mapping( volume ) == NULL)
	  delete_link( volume_link, volume_list );
      }
    } /* end if Overture mapping */
  } /* end for all volume mappings */

/* loop over all overture mappings in the list */
  for (i=0; i<mappingInfo_ptr->mappingList.getLength(); i++){
    MappingRC &map = mappingInfo_ptr->mappingList[i];

/* loop through the list of existing mappings to see if it is already there */
    for (volume_link = volume_list->first; volume_link != NULL; 
	 volume_link = volume_link->next){
      volume = (volume_mapping *) volume_link->data;
/* compare the pointers */
      if (volume->volume_data_ptr == &map) break;
    }
    if ( volume_link != NULL ){
      printf("Using an existing volume mapping for overture mapping #%i.\n", i);
    }
    else{
/* First try to use the same name as the overture mapping */
      cname = (const char *) map.getName(Mapping::mappingName);
      strcpy(grid_name, cname);
      j = 0;
      do{
/* check the name for uniqueness */
	for (volume_link = volume_list->first; volume_link != NULL; 
	     volume_link = volume_link->next){
	  volume = (volume_mapping *) volume_link->data;
	  if (strcmp(grid_name, volume->name) == 0) break;
	}
	if ( volume_link != NULL ){
	  printf("A volume mapping named `%s' already exits.\n", grid_name);
/* make a new name */
	  sprintf(grid_name, "%s-%i", cname, ++j);
	  printf("Trying `%s' instead.\n", grid_name);
	}
      } while(volume_link != NULL);

      printf("Making a new volume mapping for overture mapping #%i named `%s'.\n", 
	     i, grid_name);

/* make a new volume_mapping structure */
      volume = new_volume_mapping(grid_name);
/* general data */
      volume->type            = "Overture mapping"; 
      volume->inverse_known   = FALSE;
      volume->volume_data_ptr = &map; 
      volume->forward_mapping = overture_mapping; 
      volume->inverse_mapping = NULL;
      volume->set_volume      = set_overture_mapping; 
      volume->delete_volume_data = delete_overture_data;

/* insert the new volume_mapping in the global volume mapping list */
      new_link( volume_list )->data = volume;
    }
/* copy the number of points */
    volume->r1_points = map.getGridDimensions( 0 );
    volume->r2_points = map.getGridDimensions( 1 );
    volume->r3_points = map.getGridDimensions( 2 );

/* copy the bounding box */
    volume->bb->x_min = (real) map.getRangeBound(0,0);
    volume->bb->y_min = (real) map.getRangeBound(0,1);
    volume->bb->z_min = (real) map.getRangeBound(0,2);

    volume->bb->x_max = (real) map.getRangeBound(1,0);
    volume->bb->y_max = (real) map.getRangeBound(1,1);
    volume->bb->z_max = (real) map.getRangeBound(1,2);
  } /* end for */

}

static void 
set_overture_mapping(input_output *io_ptr, volume_mapping *volume, 
		      linked_list *surface_list ){
  int quit=0, replot=1, icom;
  real dx1, dx2, dx3;
#include "set_overture_com.h"

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

    switch( get_command( io_ptr, "overture>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case SHOW:
      printf("r1-points: %i\nr2-points: %i\nr3-points: %i\n",
	     volume->r1_points, volume->r2_points, volume->r3_points);
      printf("\n");
      break;

/* change the plot mode */
    case PLOT_MODE:
       volume->plot_mode = volume_plot_mode(io_ptr, NULL, volume, 
					   volume->plot_mode, NULL); 
/* don't need to redraw the surface grids! */
      break;

    case R1_LINES:
      volume->r1_points = int_max(2, get_int( io_ptr, "Number of r1-lines:", 
				 volume->r1_points, NO_INDENT));
      /*      map.setGridDimensions(0, volume->r1_points);*/
      replot = 1;
      break;

    case R2_LINES:
      volume->r2_points = int_max(2, get_int( io_ptr, "Number of r2-lines:", 
				 volume->r2_points, NO_INDENT));
      /*      map.setGridDimensions(1, volume->r2_points);*/
      replot = 1;
      break;

    case R3_LINES:
      volume->r3_points = int_max(2, get_int( io_ptr, "Number of r3-lines:", 
				 volume->r3_points, NO_INDENT));
      /*      map.setGridDimensions(2, volume->r3_points);*/
      replot = 1;
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
      quit = 1;
      break;

    default:
      break;

    }
  } while (!quit);
}

static void
delete_overture_data( void *volume_data_ptr ){
  /* don't do anything for now */
}

static int
overture_mapping( grid_point *gp_ptr, void *tuut ){
  static RealArray r(1,3), xp(1,3), xr(1,3,3);

/* cast the data pointer to the right type */
  MappingRC & map = * (MappingRC *) tuut;

  r(0,0) = gp_ptr->r;
  r(0,1) = gp_ptr->s;
  r(0,2) = gp_ptr->t;

  map.map(r, xp, xr);

  gp_ptr->x = xp(0,0);
  gp_ptr->y = xp(0,1);
  gp_ptr->z = xp(0,2);

  gp_ptr->xr = xr(0,0,0);
  gp_ptr->xs = xr(0,0,1);
  gp_ptr->xt = xr(0,0,2);

  gp_ptr->yr = xr(0,1,0);
  gp_ptr->ys = xr(0,1,1);
  gp_ptr->yt = xr(0,1,2);

  gp_ptr->zr = xr(0,2,0);
  gp_ptr->zs = xr(0,2,1);
  gp_ptr->zt = xr(0,2,2);

  return 1;
}
}
