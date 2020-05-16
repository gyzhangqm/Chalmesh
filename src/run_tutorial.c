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

static void 
tutorial_intro( void ){
  printf("\n"
"To construct an overlapping grid with Chalmesh you go through the following\n"
"steps:\n\n"
" o) Generate overlapped surface mappings that describe the geometry.\n"
" o) Create volume mappings that cover the domain close to each surface mapping.\n"
" o) Fill the remaining parts of the domain with background mappings. To make it\n"
"    possible to use an efficient solution algorithm in the partial differential\n"
"    equation (PDE) solver, the background mappings should be made as simple as\n"
"    possible, preferably orthogonal and, if possible, with constant grid sizes.\n"
"    Cartesian and cylindrical mappings are good for this purpose.\n"
" o) Mark edges of volume mappings that should cut holes in faces of other\n"
"    volume mappings during the hole-cutting part of the overlap algorithm.\n"
"    For example, consider a wing-fuselage configuration with one grid around the\n"
"    fuselage and one O-grid around the wing. In this case the edge in the wing\n"
"    grid that coincides with both the fuselage and the wing surfaces should\n"
"    have a positive edge label, so that all points inside the wing on the\n"
"    fuselage can be removed. If the hole is described by several edges in\n"
"    several volume mappings, the same positive edge-label should be used for all\n"
"    those edges.\n"
" o) Identify the boundary of the computational domain by assigning a positive\n"
"    surface-label to the faces of the volume mappings that are aligned with the\n"
"    boundary. The value of the surface-label should be the same for all faces\n"
"    that are aligned with the same part of the boundary. Hence, if the domain\n"
"    is simply connected, only one surface-label value should be used. For a\n"
"    doubly, connected domain, two values should be used, and so on.\n"
"    Note that negative surface-labels are reserved for future use.\n"
" o) Assign boundary condition values to the appropriate faces of the mappings.\n"
"    The boundary condition value makes no difference for the overlapping grid\n"
"    algorithm, but is necessary information for a PDE solver that uses the grid.\n"
" o) Order the mappings in a linear hierarcy. The first grid will get the highest\n"
"    priority in the overlapping grid. Grid points from a grid with higher\n"
"    priority will be preferred where several grids overlap each other.\n"
" o) Set the discretization and interpolation widths to suit the discretization\n"
"    of the PDE you want to solve. Also select the interpolation type.\n"
" o) Compute the overlapping grid.\n"
" o) Save the overlapping grid.\n\n");
}

void 
run_tutorial( input_output *io_ptr ){
  FILE *fp;
  int quit=0, icom;

#include "run_tutorial_com.h"

/* always present some general hints and a brief outline of the principles */
  printf(
"Welcome to the tutorial of Chalmesh version "VERSION".\n"
"\n"
"A user's guide to the program can be found at\n"
"http://www.na.chalmers.se/~andersp/chalmesh/chalmesh.html\n"
"\n");
  general_help();

  do{

    switch( get_command( io_ptr, "tutorial>", COMMAND, LAST_COM, 
		       LEVEL, SAVE_ON_COPY, ARGUMENT) ){

    case INTRO:
/* overlapping grid steps */
      tutorial_intro();
      break;

    case INTERACT:
/* command interpreter */
      general_help();
      break;

    case STRETCHED_ELLIPSOID:
/* stretched-ellipsoid */
      if ((fp = open_this_ascii_file( "stretched-ellipsoid.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
/* exit to the top command level */
	quit = 1;
      }
      break;

    case PROJECTION:
/* stretched-ellipsoid */
      if ((fp = open_this_ascii_file( "projection.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
/* exit to the top command level */
	quit = 1;
      }
      break;

    case HALF_ELLIPSOID:
/* half-ellipsoid */
      if ((fp = open_this_ascii_file( "half-ellipsoid.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
/* exit to the top command level */
	quit = 1;
      }
      break;

    case QUARTER_ELLIPSOID:
/* quarter-ellipsoid */
      if ((fp = open_this_ascii_file( "quarter-ellipsoid.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
/* exit to the top command level */
	quit = 1;
      }
      break;

    case EXTERNAL_GRIDS:
/* externally defined grids */
      if ((fp = open_this_ascii_file( "external-grids.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
/* exit to the top command level */
	quit = 1;
      }
      break;

    case HYPERBOLIC:
/* externally defined grids */
      if ((fp = open_this_ascii_file( "hyperbolic.com", 'r', 1, 0)) != NULL){
/* close the previous input file */
	if (io_ptr->read_command != stdin)
	  fclose( io_ptr->read_command );
/* assign the new file as the input stream */
	io_ptr->read_command = fp;
/* exit to the top command level */
	quit = 1;
      }
      break;

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );

      break;

    case EXIT:
/* exit tutorial */
      quit = 1;
      break;

    default:
      break;

    }
  } while (!quit);

}

