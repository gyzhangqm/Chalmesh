#include "volume.h"
#include "overmap.h"

/* static prototypes */
static void 
set_overture_mapping(input_output *io_ptr, volume_mapping *volume, 
		      linked_list *surface_list );
static int
inverse_overture_mapping( grid_point *gp_ptr, void *tuut );
static void
delete_overture_data( void *volume_data_ptr );
/* end prototypes */

void
init_overture_map(volume_mapping *volume, void *overture_mappings, int grid_number){
  void *map_ptr;
  int r1_points, r2_points, r3_points;
  grid_point gp;
/* get data from C++ */
  map_ptr = get_map_ptr(overture_mappings, grid_number);
  r1_points = get_grid_dim(overture_mappings, grid_number, 1);
  r2_points = get_grid_dim(overture_mappings, grid_number, 2);
  r3_points = get_grid_dim(overture_mappings, grid_number, 3);

/* general data */
  volume->type            = "Overture mapping"; 
  volume->inverse_known   = FALSE;
  volume->volume_data_ptr = map_ptr; 
  volume->forward_mapping = overture_mapping; 
  volume->inverse_mapping = NULL;
  volume->set_volume      = set_overture_mapping; 
  volume->delete_volume_data = delete_overture_data;

  volume->r1_points = r1_points;
  volume->r2_points = r2_points;
  volume->r3_points = r3_points;

/* compute an approximate bounding box */
  gp.r = 0.0; gp.s = 0.0; gp.t = 0.0;
  forward_volume_mapping( &gp, volume );
  volume->bb->x_min = gp.x;
  volume->bb->y_min = gp.y;
  volume->bb->z_min = gp.z;

  gp.r = 1.0; gp.s = 1.0; gp.t = 1.0;
  forward_volume_mapping( &gp, volume );
  volume->bb->x_max = gp.x;
  volume->bb->y_max = gp.y;
  volume->bb->z_max = gp.z;
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
