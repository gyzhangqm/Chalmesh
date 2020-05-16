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
#include "ogen_interface.h"

void
main(int argc, char **argv){
  overlapping_3d_grid *over3d;
  surface_mapping *surface;
  volume_mapping *volume=NULL;

  surface_mapping_list *surface_mappings;
  volume_mapping_list *volume_mappings;
  overlapping_3d_list *over3d_list;

  surface_mapping_link *surface_link;
  volume_mapping_link *volume_link;

  overlapping_3d_link *over3d_link;

  char *token, *name, *file_name, new_name[80], new_file[80];
  FILE *fp;

  int i, normalized_plotting=0;
  char *display = NULL;

  GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lmodel_ambient[] = { 1.0, 1.0, 1.0, 1.0 };

  input_output io, *io_ptr=&io;
  int icom, quit = 0, replot_surface = TRUE, replot_volume = TRUE
    , replot_overlap = TRUE, f_number, tmp_fd, arg_len
/* plot grid surfaces, grid boundaries, parameter arrows and the coordinate cage */
      , s_plot_mode = 2 | 4 | 256, v_plot_mode = 2 | 4 | 256, graphics = TRUE;
  const int n_windows = OVERLAP_W + 1;
  char *titles[3], tmp_title[120], *comfile=NULL;
  double x_pos[3], y_pos[3], width[3], height[3];
  bounding_box bb_vol, bb_surf, bb_over;
  clip_planes clip_info = {0, {{5.0, 0.0, 0.0}, {5.0, 0.0, 0.0}, {5.0, 0.0, 0.0}, 
				 {5.0, 0.0, 0.0}}};
/* pointer to overture mappings */
#ifdef OVERTURE_INTERFACE
  void *overture_mappings;
#include "overmap.h"
#endif

/* set up the commands */
#include "chalmesh_com.h"

/* define the titles and position of the windows */
  sprintf( tmp_title, "Surface mappings (Chalmesh version %s)", VERSION );
  titles[SURFACE_W] = (char *) malloc((strlen(tmp_title)+1)*sizeof(char));
  strcpy(titles[SURFACE_W], tmp_title);
  x_pos[SURFACE_W]  = 0.52;
  y_pos[SURFACE_W]  = 0.48;
  width[SURFACE_W]  = 0.48;
  height[SURFACE_W] = 0.45;

  sprintf( tmp_title, "Volume mappings (Chalmesh version %s)", 
	  VERSION );
  titles[VOLUME_W] = (char *) malloc((strlen(tmp_title)+1)*sizeof(char));
  strcpy(titles[VOLUME_W], tmp_title);
  x_pos[VOLUME_W]  = 0.52;
  y_pos[VOLUME_W]  = 1.0;
  width[VOLUME_W]  = 0.48;
  height[VOLUME_W] = 0.45;

  sprintf( tmp_title, "Overlapping grids (Chalmesh version %s)", 
	  VERSION );
  titles[OVERLAP_W] = (char *) malloc((strlen(tmp_title)+1)*sizeof(char));
  strcpy(titles[OVERLAP_W], tmp_title);
  x_pos[OVERLAP_W]  = 0.0;
  y_pos[OVERLAP_W]  = 1.0;
  width[OVERLAP_W]  = 0.48;
  height[OVERLAP_W] = 0.45;

/* initialize the linked lists for the surface grids, overlapping surface grids, */
/* and volume grids */
  surface_mappings = new_linked_list();
  volume_mappings = new_linked_list();
  over3d_list  = new_linked_list();

/* initialize the bounding boxes */
  surface_list_bb( &bb_surf, surface_mappings ); 
  volume_list_bb( &bb_vol, volume_mappings ); 
  overlap_list_bb( &bb_over, over3d_list ); 

/* read the command line arguments */
  for (i = 1; i < argc; i++) {
    arg_len = strlen(argv[i]);
    if (!strncmp(argv[i], "-display", arg_len)) {
      if (++i >= argc){
	printf("follow -display option with display parameter\n");
	exit(1);
      }
      display = argv[i];
    } 
    else if (!strncmp(argv[i], "-nographics", arg_len)) {
      graphics = FALSE;
    } 
    else if (!strncmp(argv[i], "-command", arg_len)){
      if (++i >= argc){
	printf("follow -command option with a filename\n");
	exit(1);
      }
/* copy the filename for later use */
      comfile = (char *) malloc( (strlen(argv[i])+1) * sizeof(char) );
      strcpy( comfile, argv[i] );
    }
    else{
      printf("Usage:\n"
	     "%s [-display DISPLAY] [-command command_file] [-nographics]\n", argv[0]);
      exit(1);
    }
  }

  ogl_init_prompt( io_ptr );

/* initialize the interface to overture */
#ifdef OVERTURE_INTERFACE
  overture_mappings = initializeOvertureMappings();
#endif

/* initialize the X-windows connection and openGL */
  if (graphics && ogl_init( display, titles, x_pos, y_pos, width, height, n_windows )){

    for (i=0; i<n_windows; i++){
/* compute initial modelview transformation */
      ogl_standard_view( i, &bb_vol );

/*** configure the OpenGL context for rendering ***/
      glEnable(GL_DEPTH_TEST);	/* enable depth buffering */
      glEnable(GL_LIGHTING);	/* enable lighting */

      glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

      glEnable(GL_LIGHT0);
      glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
      glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
      glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

      glClearColor( 0.7, 0.7, 0.4, 1.0 );

      replot_surface = TRUE; 
      replot_volume  = TRUE;
      replot_overlap  = TRUE;
    }
  }
  else{
    printf("Running without graphics output.\n");
  }

  if (comfile != NULL){
/* try to interpret the argument as a command file */
    if ((fp = open_this_ascii_file( comfile, 'r', 1, 1)) != NULL){
/* close the previous input file */
      if (io_ptr->read_command != stdin)
	fclose( io_ptr->read_command );
/* assign the new file as the input stream */
      io_ptr->read_command = fp;
    }
    else
      printf(
"Warning: Attempts to interpret the argument `%s' as a command file failed \n"
"so the argument was ignored.\n", comfile);
/* free the memory */
    free( comfile );
  }

  welcome_to_chalmesh();

/* open a log file where all commands will be saved */
/* iterate on finding a unique file name */
  for (f_number = 1; ; f_number++){
    sprintf(new_file, "chalmesh-log-%-i.com", f_number);
/* make sure the file does not exist */
    if ( (tmp_fd = open(new_file, O_WRONLY, 0666)) == -1) break;
    close( tmp_fd );
  }
  if ( (fp = open_this_ascii_file( new_file, 'w', 1, 0 )) != NULL){
    printf("All commands are beeing saved in the file `%s'\n", new_file);
/* close the previous copy file */
    if (io.copy_command != NULL)
      fclose( io.copy_command );
/* assign the new file as the copy stream */
    io.copy_command = fp;
/* copy the file-name */
    io.copy_file_name = (char *) malloc( (strlen(new_file)+1) * sizeof(char) );
    io.copy_file_name = strcpy( io.copy_file_name, new_file );
  }

/* command loop */
  do{

/* redraw surface mappings */
    if (replot_surface && ogl_start_plot(OGL_NEW_PLOT, SURFACE_W, 
					 ogl_length_scale(&bb_surf)))
      {
/* redraw all surface grids */
	for (surface_link = surface_mappings->first; surface_link != NULL; 
	     surface_link = surface_link->next)
	  draw_surface( surface_link->data, s_plot_mode);
	if (s_plot_mode & 256 && surface_mappings->n_members > 0)
	  ogl_bb_cage( &bb_surf, OGL_WHITE );

	ogl_end_plot();
	replot_surface = FALSE;
      }
/* redraw volume mappings */
    if (replot_volume && ogl_start_plot(OGL_NEW_PLOT, VOLUME_W, 
					ogl_length_scale(&bb_vol)))
      {
/* redraw the volume grids */
	for (volume_link = volume_mappings->first; volume_link != NULL; 
	     volume_link = volume_link->next)
	  draw_volume( volume_link->data, v_plot_mode, &clip_info );
	if (v_plot_mode & 256 && volume_mappings->n_members > 0)
	  ogl_bb_cage( &bb_vol, OGL_WHITE );

	ogl_end_plot();
	replot_volume = FALSE;
      }
/* redraw the overlapping volume grids */
    if (replot_overlap && ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, 
					 ogl_length_scale(&bb_over)))
      {
	for (over3d_link = over3d_list->first; over3d_link != NULL;
	     over3d_link = over3d_link->next)
	  {
	    over3d = over3d_link->data;
	    draw_overlapping_3d_grid(over3d, over3d->plot_mode);
	    if (over3d->plot_mode & 256) 
	      ogl_bb_cage( &bb_over, OGL_WHITE );
	  }
	ogl_end_plot();
	replot_overlap = FALSE;
      }

    switch( get_command( io_ptr, "chalmesh>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

/* turn on the graphics */
    case GRAPHICS_ON:
/* initialize the X-windows connection and openGL */
      if ( ogl_init( display, titles, x_pos, y_pos, width, height, n_windows) ){

	for (i=0; i<n_windows; i++){
/* compute initial modelview transformation */
	  ogl_standard_view( i, &bb_vol );

/*** configure the OpenGL context for rendering ***/
	  glEnable(GL_DEPTH_TEST);	/* enable depth buffering */
	  glEnable(GL_LIGHTING);	/* enable lighting */

	  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
	  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

	  glEnable(GL_LIGHT0);
	  glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
	  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
	  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

	  glClearColor( 0.7, 0.7, 0.4, 1.0 );

	  replot_surface = TRUE; 
	  replot_volume  = TRUE;
	  replot_overlap = TRUE;
	}

	ogl_standard_view( SURFACE_W, &bb_surf );
	ogl_standard_view( VOLUME_W,  &bb_vol );
	ogl_standard_view( OVERLAP_W, &bb_over );
      }
    break;

    case GRAPHICS_OFF:
      ogl_close( );
      break;

/* make a new surface grid */
    case NEW_SURF:
      sprintf( new_name, "surface-%i", surface_mappings->n_members+1);
      token = get_word( io_ptr, "New surface grid name: ", new_name, 1 );
/* check for uniqueness */
      for (surface_link = surface_mappings->first; surface_link != NULL; 
	   surface_link = surface_link->next){
	surface = surface_link->data;
	if (strcmp(token, surface->name) == 0) break;
      }
      if ( surface_link != NULL ){
	printf("A surface grid with that name already exits.\n");
      }
      else{
/* copy the name */
	name = (char *) malloc( (strlen(token)+1)*sizeof(char) );
	name = strcpy( name, token );
/* set the specific type of surface grid */
	if ((surface=choose_surface( io_ptr, name)) != NULL ){
/* insert in the global surface grid list */
	  surface_link = new_link( surface_mappings );
	  surface_link->data = surface;

/* recompute the global bounding box */
	  surface_list_bb( &bb_surf, surface_mappings ); 
/* redraw the surface grids */
	  replot_surface = TRUE;
	}
	free( name );
      }
      break;
      
/* make a new surface grid */
    case NEW_SURF_FROM_VOLUME:
      sprintf( new_name, "surface-%i", surface_mappings->n_members+1 );
      token = get_word( io_ptr, "New surface grid name: ", new_name, 1 );
/* token is only a pointer to a string allocated in ogl_plot that will be overwritten */
/* by the following call to get_volume_mapping */
      strcpy(new_name, token);
/* check for the name for uniqueness */
      for (surface_link = surface_mappings->first; surface_link != NULL; 
	   surface_link = surface_link->next){
	surface = surface_link->data;
	if (strcmp(new_name, surface->name) == 0) break;
      }
      if ( surface_link != NULL ){
	printf("A surface grid with that name already exits.\n");
      }
      else if ((volume_link = get_volume_mapping(io_ptr, volume_mappings, LEVEL+1))){
	volume = volume_link->data;
/* make a surface mapping */
	surface = new_surface_mapping(new_name);
/* surface from a face of a volume grid */
	init_volume_face( surface, volume );
	set_surface( io_ptr, surface );
/* insert in the global surface grid list */
	surface_link = new_link( surface_mappings );
	surface_link->data = surface;
/* recompute the global bounding box */
	surface_list_bb( &bb_surf, surface_mappings ); 
/* redraw the surface grids */
	replot_surface = TRUE;
      }
      break;
      
/* delete a surface grid */
    case DELETE_SURF:
      if ((surface_link = get_surface( io_ptr, surface_mappings, NO_INDENT ))){
	surface = surface_link->data;
	if (delete_surface_mapping( surface ) == NULL){
/* remove from the surface grid list */
	  delete_link( surface_link, surface_mappings );
/* recompute the global bounding box */
	  surface_list_bb( &bb_surf, surface_mappings ); 
/* update the plotting window */
	  replot_surface = TRUE;
	}
	else{
	  printf("This surface grid can not be deleted because it is used by "
		 "a volume grid.\n");
	}
      }
      break;

/* change a surface grid */
    case CHANGE_SURF:
      if ((surface_link = get_surface( io_ptr, surface_mappings, NO_INDENT ))){
	surface = surface_link->data;
	set_surface( io_ptr, surface );
/* recompute the global bounding box */
	surface_list_bb( &bb_surf, surface_mappings ); 
/* update the plotting window */
	replot_surface = TRUE;
      }
      break;

/* check a surface grid */
    case CHECK_SURF:
      if ((surface_link = get_surface( io_ptr, surface_mappings, NO_INDENT ))){
	surface = surface_link->data;
	check_surface( io_ptr, surface );
      }
      break;

/* save a surface grid on file */
    case SAVE_SURF:
      if ((surface_link = get_surface( io_ptr, surface_mappings, NO_INDENT ))){
	surface = surface_link->data;
/* open an ascii file and save the component grid */
	if ( (fp = open_ascii_file( io_ptr, "surface grid file name: ",
				   "surf.3d", &file_name, 'w', 0, 
				   SAVE_ON_COPY[SAVE_SURF])) != NULL){
	  save_surface_mapping( surface, fp );
	  fclose( fp );
	}
      }
      break;

/* change the boundary condition */
    case BOUNDARY_CONDITION:
      if ((volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){
	volume = volume_link->data;
	set_volume_bc( io_ptr, volume );
/* update the plotting window */
	replot_volume = TRUE;
      }
      break;

/* set the surface label on a volume component */
    case SURF_LABEL:
      if ((volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){
	volume = volume_link->data;
	set_surf_label( io_ptr, volume );
/* update the volume plotting window */
	replot_volume = TRUE;
      }
      break;

/* set the edge label on a volume component */
    case EDGE_LABEL:
      if ((volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){
	volume = volume_link->data;
	set_edge_label( io_ptr, volume );
/* update the volume plotting window */
	replot_volume = TRUE;
      }
      break;

/* prevent / enforce hole-cutting in a component grid */
    case CUT_HOLES:
      if ((volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){
	volume = volume_link->data;
	volume->no_hole = !get_yes_no(io_ptr, 
				      "Apply the hole-cutting algorithm to this grid? ", 
				      LEVEL+1, 1);
      }
      break;

#ifdef OVERTURE_INTERFACE
    case OVERTURE:
      makeMappingsAndPlotThem( overture_mappings );
/* make/update the overture mappings */
      init_overture_maps(volume_mappings, overture_mappings);
/* recompute the global bounding box */
      volume_list_bb( &bb_vol, volume_mappings ); 
/* redraw the surface grids */
      replot_volume = TRUE;
      break;
#endif

    case NEW_VOLUME:
      sprintf( new_name, "volume-%i", volume_mappings->n_members+1);
      token = get_word( io_ptr, "New volume grid name: ", new_name, 1 );
/* check for uniqueness */
      for (volume_link = volume_mappings->first; volume_link != NULL; 
	   volume_link = volume_link->next){
	volume = volume_link->data;
	if (strcmp(token, volume->name) == 0) break;
      }
      if ( volume_link != NULL ){
	printf("A volume grid with that name already exits.\n");
      }
      else{
/* copy the name */
	name = (char *) malloc( (strlen(token)+1)*sizeof(char) );
	name = strcpy( name, token );
/* set the specific type of volume grid */
	if ((volume=choose_volume( io_ptr, name, surface_mappings))){
/* insert in the global volume grid list */
	  volume_link = new_link( volume_mappings );
	  volume_link->data = volume;

/* recompute the global bounding box */
	  volume_list_bb( &bb_vol, volume_mappings ); 
/* redraw the surface grids */
	  replot_volume = TRUE;
	}
	free( name );
      }
      break;

    case CHANGE_VOLUME:
      if ((volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){
	volume = volume_link->data;
	set_volume( io_ptr, volume, surface_mappings );
/* recompute the global bounding box */
	volume_list_bb( &bb_vol, volume_mappings ); 
/* don't update the plotting window */
      }
      break;

    case TRANSFORM_VOLUME:
      if ((volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){
	volume = volume_link->data;
	transform_volume(io_ptr, volume);
      }
      break;

    case FILE_VOLUME:
      if ((volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){
	volume = volume_link->data;
	file_volume(io_ptr, volume);
      }
      break;

    case CHECK_VOLUME:
      if ((volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){
	volume = volume_link->data;
	check_jacobian(io_ptr, volume);
      }
      break;

    case DELETE_VOLUME:
      if ((volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){
	volume = volume_link->data;
	if (delete_volume_mapping( volume ) == NULL){
/* remove from the volume mapping list */
	  delete_link( volume_link, volume_mappings );
/* recompute the global bounding box */
	  volume_list_bb( &bb_vol, volume_mappings ); 
/* update the plotting window */
	  replot_volume = TRUE;
	}
	else{
	  printf("This volume grid can not be deleted because it is either used by\n"
		 "a surface grid or because it shares arrays with a component grid\n"
		 "in the overlapping grid.\n");
	}
      }
      break;

    case NEW_3D_OVERLAP:
/* check if there are any mappings left */
      if (volume_mappings->n_members == 0){
	printf("There are no volume mappings to construct the overlapping "
	       "grid from.\n");
	break;
      }

/* ask for a name */
      sprintf( new_name, "over3d-%i", over3d_list->n_members+1);
      token = get_word( io_ptr, "New overlapping volume grid name: ", new_name, 1 );

/* check the name for uniqueness */
      for (over3d_link = over3d_list->first; over3d_link != NULL; 
	   over3d_link = over3d_link->next){
	over3d = over3d_link->data;
	if (strcmp(token, over3d->name) == 0) break;
      }

      if ( over3d_link != NULL ){
	printf("An overlapping grid with that name already exists.\n");
	break;
      }

/* define the overlapping volume grid */
      over3d = new_overlapping_3d_grid(token, volume_mappings);
/* insert it in the list */
      over3d_link = new_link( over3d_list );
      over3d_link->data = over3d;
/* set the standard viewpoint */
      ogl_standard_view( OVERLAP_W, over3d->bb );
/* modify its default parameters, the component list, and compute the intergrid */
/* communication data */
      set_3d_overlap( io_ptr, over3d, volume_mappings );
      overlap_list_bb( &bb_over, over3d_list ); 
      break;
    
    case CHANGE_3D_OVERLAP:
      if ((over3d_link = get_over3d_ptr( io_ptr, over3d_list))){
	over3d = over3d_link->data;
/* modify the overlap parameters, the component list, and compute the intergrid */
/* communication data */
	set_3d_overlap( io_ptr, over3d, volume_mappings );
	overlap_list_bb( &bb_over, over3d_list ); 
      }
      break;

    case DELETE_3D_OVERLAP:
      if ((over3d_link = get_over3d_ptr( io_ptr, over3d_list))){
	over3d = over3d_link->data;
	over3d = delete_over3d_grid( over3d );
/* remove the link from the list */
	delete_link( over3d_link, over3d_list );
	overlap_list_bb( &bb_over, over3d_list ); 
      }
      replot_overlap = TRUE;
      break;

/* check if a background grid is likely to be fine enogh for an overlapping subgrid */
/*    case CHECK_BACKGROUND:*/
/* get the overlapping volume grid and the mapping to add */
/*       printf("Identify the overlapping volume grid and the volume mapping to add:\n"); */
/*       if ((over3d_link = get_over3d_ptr( io_ptr, over3d_list)) && */
/* 	  (volume_link = get_volume_mapping( io_ptr, volume_mappings, NO_INDENT ))){ */
/* 	if (over3d_link == NULL || volume_link == NULL) break; */
/* 	over3d = over3d_link->data; */
/* check that the overlapping grid is incomplete and valid */
/* 	if (!over3d->valid || !over3d->overlapping_sub_grid){ */
/* 	  printf( */
/* "This function is intended for an incomplete, but valid overlapping grid, which can \n" */
/* "be made with the `make-overlapping-volume-grid' command with the \n" */
/* "`overlapping-sub-grid' flag set.\n"); */
/* 	  break; */
/* 	} */
/* 	background = volume_link->data; */
/* check that the background grid doesn't already participate in the overlapping grid */
/* 	for (volume_link = over3d->mapping_list->first; volume_link != NULL; */
/* 	     volume_link = volume_link->next){ */
/* 	  volume = volume_link->data; */
/* 	  if (!strcmp(volume->name, background->name)) break; */
/* 	} */
/* 	if (volume_link != NULL){ */
/* 	  printf( */
/* "This function is intended for components which are not yet included in an\n" */
/* "overlapping subgrid. The function can therefore not be applied since the\n" */
/* "background grid `%s' is already included in the overlapping grid `%s'.\n",  */
/* 		 background->name, over3d->name); */
/* 	  break; */
/* 	} */

/* 	check_background(over3d, background); */
/*       }	   */
/*       break; */

/* list all entities */
    case SHOW:
/* surfaces */
      if (surface_mappings->n_members == 0){
	printf("No surface grids are defined.\n");
      }
      else{
	printf("The following surface grids are defined:\n"
	       "Color\tName\tType\n");
	for (surface_link = surface_mappings->first; surface_link != NULL; 
	     surface_link = surface_link->next){
	  surface = surface_link->data;
	  printf("%s\t%s\t%s\n", ogl_color_name(surface->color),
		 surface->name, surface->type);
	}
      }
/* volume grids */
      printf("\n");
      if (volume_mappings->n_members == 0){
	printf("No volume grids are defined.\n");
      }
      else{
	printf("The following volume grids are defined:\n"
	       "Color\tName\tType\n");
	for (volume_link = volume_mappings->first; volume_link != NULL; 
	     volume_link = volume_link->next){
	  volume = volume_link->data;
	  printf("%s\t%s\t%s\n", ogl_color_name(volume->color),
		 volume->name, volume->type);
	}
      }
/* overlapping volume grids */
      printf("\n");
      if (over3d_list->n_members == 0){
	printf("No overlapping volume grids are defined.\n");
      }
      else{
	printf("The following overlapping volume grids are defined:\n"
	       /*"Name\tStatus\n"*/);
	for (over3d_link = over3d_list->first; over3d_link != NULL; 
	     over3d_link = over3d_link->next){
	  over3d = over3d_link->data;
	  printf("%s\n", over3d->name);
/* 	  printf("%s\t%s, %s\n", over3d->name, (over3d->valid)? "valid": "invalid",  */
/* 		 (over3d->overlapping_sub_grid)? "incomplete sub-grid": "complete grid"); */
	}
      }
      break;

    case SAVE_OVERLAP:
      if ((over3d_link = get_over3d_ptr( io_ptr, over3d_list))){
	over3d = over3d_link->data;
	if (over3d->valid /*&& !over3d->overlapping_sub_grid*/)
	  save_overlapping_grid( io_ptr, over3d );
	else /*if (!over3d->valid)*/
	  printf(
		 "Sorry, the overlapping grid `%s' is NOT VALID and\n"
		 "will therefore not be saved!\n", over3d->name);
/* 	else if (over3d->overlapping_sub_grid) */
/* 	  printf( */
/* 		 "Sorry, the overlapping grid `%s' is INCOMPLETE and\n" */
/* 		 "will therefore not be saved!\n", over3d->name); */
      }
      break;

    case TUTORIAL:
      run_tutorial( io_ptr );
      break;

    case RESET_ALL:
/* reset */
/* delete all overlapping grids */
      for (over3d_link = over3d_list->first; over3d_link != NULL;
	   over3d_link = over3d_link->next)
      {
	over3d = over3d_link->data;
	over3d = delete_over3d_grid( over3d );
      }

/* delete all volume grids */
      for (volume_link = volume_mappings->first; volume_link != NULL; 
	   volume_link = volume_link->next){
	volume_link->data = delete_volume_mapping( volume_link->data );
      }

/* delete all surface grids */
      for (surface_link = surface_mappings->first; surface_link != NULL; 
	   surface_link = surface_link->next){
	surface_link->data = delete_surface_mapping( surface_link->data );
      }

/* delete any volume grids that could not be deleted before the surface grids */
      for (volume_link = volume_mappings->first; volume_link != NULL; 
	   volume_link = volume_link->next){
	volume_link->data = delete_volume_mapping( volume_link->data );
      }

/* delete the linked lists */
      delete_linked_list( over3d_list );
      delete_linked_list( surface_mappings );
      delete_linked_list( volume_mappings );

/* make new empty lists */
      over3d_list = new_linked_list();
      surface_mappings = new_linked_list();
      volume_mappings = new_linked_list();

/* initialize the bounding boxes */
      surface_list_bb( &bb_surf, surface_mappings ); 
      volume_list_bb(  &bb_vol,  volume_mappings ); 
      overlap_list_bb( &bb_over, over3d_list ); 

      for (i=0; i<n_windows; i++)
	ogl_reset_modelview(i);

/* reset the plot modes */
      s_plot_mode = 2 | 4 | 256;
      v_plot_mode  = 2 | 4 | 256;

/* update all plotting windows */
      replot_surface = TRUE;
      replot_volume  = TRUE;
      replot_overlap = TRUE;
      break;

/* change the plot mode */
    case SURFACE_PLOT_MODE:
      s_plot_mode = surface_plot_mode(io_ptr, surface_mappings, NULL, 
				       s_plot_mode );
/* don't need to draw the surface grids */
      break;

    case VOLUME_PLOT_MODE:
      v_plot_mode = volume_plot_mode(io_ptr, volume_mappings, NULL,
				     v_plot_mode, &clip_info );
/* don't need to draw the volume grids */
      break;

    case OVER_3D_PLOT_MODE:
      if ((over3d_link = get_over3d_ptr( io_ptr, over3d_list))){
	over3d = over3d_link->data;
	over3d->plot_mode = overlapping_3d_plot_mode(io_ptr, over3d, over3d->plot_mode );
      }
/* don't need to redraw the overlapping 3-d grid! */
      break;

/* reset the viewpoint in all windows */
    case R_VIEW:
      for (i=0; i<n_windows; i++)
	ogl_reset_modelview(i);
      break;

/* toggle the use of normalized coordinates in all windows */
    case NORMALIZE:
      normalized_plotting = !normalized_plotting;
      if (normalized_plotting){
	printf("Normalized coordinates toggled ON\n");
	ogl_normalize(SURFACE_W, &bb_surf);
	ogl_normalize(VOLUME_W,  &bb_vol);
	ogl_normalize(OVERLAP_W, &bb_over);
      }
      else{
	printf("Normalized coordinates toggled OFF\n");
	ogl_unnormalize(SURFACE_W, &bb_surf);
	ogl_unnormalize(VOLUME_W,  &bb_vol);
	ogl_unnormalize(OVERLAP_W, &bb_over);
      }

/* reset the viewpoint in all windows */
      for (i=0; i<n_windows; i++)
	ogl_reset_modelview(i);

      break;

/* read command-file */
    case READ_COMMAND: 
      if ( (fp = open_ascii_file( io_ptr, "command file name: ",
				 "intro.com", &file_name, 'r', 1, 
				 SAVE_ON_COPY[READ_COMMAND])) != NULL){
/* close the previous input file */
	if (io.read_command != stdin)
	  fclose( io.read_command );
/* assign the new file as the input stream */
	io.read_command = fp;
      }

      break;


/* start saving-commands */
    case START_SAVE:
      
/* try to open the new copy file */
      if ( (fp = open_ascii_file( io_ptr, "Save commands in file named: ", 
				 "commands.com", &file_name, 'w', 1,
				 SAVE_ON_COPY[START_SAVE])) != NULL){
/* close the previous copy file */
	if (io.copy_command != NULL){
/* remind the user of the file name */
	  printf("Closing the log-file `%s'\n", io_ptr->copy_file_name);
/* close the file */
	  fclose( io.copy_command );
	  io.copy_command = NULL;
/* free the space occupied by the file-name */
	  free( io.copy_file_name );
	}

/* assign the new file as the copy stream */
	io.copy_command = fp;
/* copy the file-name */
	if ((io.copy_file_name = (char *) 
	     malloc( (strlen(file_name)+1) * sizeof(char) )) == NULL)
	  printf("memory error in xcog, copy filename\n");
	io.copy_file_name = strcpy( io.copy_file_name, file_name );
      }

      break;


/* stop saving-commands */
    case STOP_SAVE:
      if (io.copy_command != NULL){
/* remind the user of the file name */
	printf("Closing the log-file `%s'\n", io_ptr->copy_file_name);
/* close the file */
	fclose( io.copy_command );
	io.copy_command = NULL;
/* free the space occupied by the file-name */
	free( io.copy_file_name );
      }

      break;

    case WARRANTY:
      warranty_info();
      break;

    case DISTRIBUTION:
      distribution_info();
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

    case QUIT:
      quit = 1;
#ifdef CLEANUP
      printf("Deallocating memory...\n"); 
#endif
      break;
  
    default:
      break;

    }
  } while( !quit );
  
/* cleanup */
#ifdef CLEANUP
/* delete the overlapping 3d grid */
  for (over3d_link = over3d_list->first; over3d_link != NULL; 
       over3d_link = over3d_link->next){ 
    over3d = over3d_link->data; 
    over3d = delete_over3d_grid( over3d ); 
  } 

/* delete all volume grids */
  for (volume_link = volume_mappings->first; volume_link != NULL;  
       volume_link = volume_link->next){ 
    volume_link->data = delete_volume_mapping( volume_link->data ); 
  } 

/* delete all surface grids */
  for (surface_link = surface_mappings->first; surface_link != NULL;  
       surface_link = surface_link->next){ 
    surface_link->data = delete_surface_mapping( surface_link->data ); 
  } 

/* delete any volume grids that could not be deleted before the surface grids */
  for (volume_link = volume_mappings->first; volume_link != NULL;  
       volume_link = volume_link->next){ 
    volume_link->data = delete_volume_mapping( volume_link->data ); 
  } 

/* delete the linked lists */
  delete_linked_list( surface_mappings ); 
  delete_linked_list( volume_mappings ); 
  delete_linked_list( over3d_list ); 
  
  printf("Done\n"); 
#endif

  ogl_close( );
  ogl_close_prompt( io_ptr );

/* cleaning up the overture mapping list */
#ifdef OVERTURE_INTERFACE
  cleanupOvertureMappings(overture_mappings);
#endif

/* free the title strings */
  for (i=0; i<n_windows; i++) 
    free(titles[i]); 
  
  printf("Good bye\n");
}
