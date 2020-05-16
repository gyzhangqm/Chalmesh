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
enum {NEW_VOLUME, CHANGE_VOLUME, CHECK_VOLUME, TRANSFORM_VOLUME, FILE_VOLUME, 
      DELETE_VOLUME, OVERTURE,
      NEW_SURF, NEW_SURF_FROM_VOLUME, CHANGE_SURF, CHECK_SURF, DELETE_SURF, 
      SAVE_SURF, BOUNDARY_CONDITION, SURF_LABEL, EDGE_LABEL, CUT_HOLES,
      NEW_3D_OVERLAP, CHANGE_3D_OVERLAP, DELETE_3D_OVERLAP, SAVE_OVERLAP, 
      SURFACE_PLOT_MODE, VOLUME_PLOT_MODE, OVER_3D_PLOT_MODE, 
      SHOW, TUTORIAL, R_VIEW, NORMALIZE,
      READ_COMMAND, START_SAVE, STOP_SAVE, RESET_ALL, 
      GRAPHICS_ON, GRAPHICS_OFF, WARRANTY, DISTRIBUTION, HELP, QUIT};
const int LAST_COM = QUIT, LEVEL=0, NO_INDENT=0;

char *COMMAND[38], *BRIEF_HELP[38];
int ARGUMENT[38], SAVE_ON_COPY[38];

COMMAND[WARRANTY] = "warranty-information";  
ARGUMENT[WARRANTY] = 0;
SAVE_ON_COPY[WARRANTY] = 1;
BRIEF_HELP[WARRANTY] = "Print the NO WARRANTY information.";

COMMAND[DISTRIBUTION] = "distribution-information";  
ARGUMENT[DISTRIBUTION] = 0;
SAVE_ON_COPY[DISTRIBUTION] = 1;
BRIEF_HELP[DISTRIBUTION] = "Print the conditions for distributing Chalmesh.";

COMMAND[SURFACE_PLOT_MODE] = "surface-plot-mode";  
/*file[SURFACE_PLOT_MODE] = "surface_plot_mode_com.h"; */
ARGUMENT[SURFACE_PLOT_MODE] = 0;
SAVE_ON_COPY[SURFACE_PLOT_MODE] = 1;
BRIEF_HELP[SURFACE_PLOT_MODE] = "Modify the graphical contents of the "
"`Surface mappings' window.";

COMMAND[VOLUME_PLOT_MODE] = "volume-plot-mode";  
/*file[VOLUME_PLOT_MODE] = "volume_plot_mode_com.h"; */
ARGUMENT[VOLUME_PLOT_MODE]  = 0;
SAVE_ON_COPY[VOLUME_PLOT_MODE]  = 1;
BRIEF_HELP[VOLUME_PLOT_MODE] = "Modify the graphical contents of the "
"`Volume mappings' window.";

COMMAND[OVER_3D_PLOT_MODE] = "overlapping-plot-mode";  
/*file[OVER_3D_PLOT_MODE] = "overlapping_3d_plot_mode_com.h"; */
ARGUMENT[OVER_3D_PLOT_MODE] = 0;
SAVE_ON_COPY[OVER_3D_PLOT_MODE] = 1;
BRIEF_HELP[OVER_3D_PLOT_MODE] = "Modify the graphical contents of the "
"`Overlapping grids' window.";

COMMAND[NORMALIZE]    = "normalize-plotting";                     
ARGUMENT[NORMALIZE]   = 0;
BRIEF_HELP[NORMALIZE] = "Toggle the use of normalized coordinates in the plotting. "
"If the bounding "
"box of the object is x0 <= x < x1, the normalized coordinates are x' = x/(x1-x0). "
"Normalized coordinates will make it easier to view objects with significantly "
"different length scales in the different coordinate directions.";

COMMAND[R_VIEW] = "standard-view"; 
ARGUMENT[R_VIEW]            = 0;
SAVE_ON_COPY[R_VIEW]            = 1;
BRIEF_HELP[R_VIEW] = "Reset the viewpoints in all graphics windows.";

COMMAND[NEW_SURF_FROM_VOLUME] = "make-surface-from-volume";
/*file[NEW_SURF_FROM_VOLUME] = "set_volume_face_com.h"; */
ARGUMENT[NEW_SURF_FROM_VOLUME] = 1;
SAVE_ON_COPY[NEW_SURF_FROM_VOLUME] = 1;
BRIEF_HELP[NEW_SURF_FROM_VOLUME] = "Define a new surface mapping from a face of " 
"a volume mapping. This is for instance useful for projecting "
"a face of one volume mapping onto the face of another volume mapping, since the "
"projection command only works with surfaces."; 

COMMAND[NEW_SURF] = "make-surface-mapping";       
/*file[NEW_SURF] = "choose_surface_com.h"; */
ARGUMENT[NEW_SURF]          = 1;
SAVE_ON_COPY[NEW_SURF]             = 1;
BRIEF_HELP[NEW_SURF] = "Define a new surface mapping by first choosing a name "
"for the surface, and then selecting the type of the surface mapping. "
"Thereafter, you proceed "
"directly to modify the parameters for that type of mapping, in the same way as "
"you would do with the command `change-surface-mapping'.";

COMMAND[CHANGE_SURF] = "change-surface-mapping";       
ARGUMENT[CHANGE_SURF]          = 1;
SAVE_ON_COPY[CHANGE_SURF]          = 1;
BRIEF_HELP[CHANGE_SURF] = "Modify an existing surface mapping. This command gives "
"you the same options for modifying a surface mapping as provided by the command "
"`make-surface-mapping', with the obvious exception that you do not specify the "
"type of mapping again.";

COMMAND[CHECK_SURF] = "check-surface-mapping";       
ARGUMENT[CHECK_SURF]          = 1;
SAVE_ON_COPY[CHECK_SURF]          = 1;
BRIEF_HELP[CHECK_SURF] = "Check the Jacobian and Hessian of a surface mapping.";

COMMAND[DELETE_SURF] = "delete-surface-mapping";       
ARGUMENT[DELETE_SURF]          = 1;
SAVE_ON_COPY[DELETE_SURF]          = 1;
BRIEF_HELP[DELETE_SURF] = "Remove an existing surface mapping.";

COMMAND[SAVE_SURF] = "save-surface-mapping";       
ARGUMENT[SAVE_SURF]            = 1;
SAVE_ON_COPY[SAVE_SURF]            = 1;
BRIEF_HELP[SAVE_SURF] = "Save the grid points of a surface mapping in a file in "
"PLOT3D ASCII format.";

COMMAND[BOUNDARY_CONDITION] = "boundary-condition";       
/*file[BOUNDARY_CONDITION] = "volume_bc_com.h"; */
ARGUMENT[BOUNDARY_CONDITION] = 1;
SAVE_ON_COPY[BOUNDARY_CONDITION] = 1;       
BRIEF_HELP[BOUNDARY_CONDITION] = "Assign boundary conditions to the faces of "
"a volume mapping. The value of the boundary condition has no "
"significance for the overlapping grid algorithm, but is helpful information "
"for a PDE solver that uses the grid.";       

COMMAND[SURF_LABEL] = "surface-label";       
/*file[SURF_LABEL] = "surf_label_com.h"; */
ARGUMENT[SURF_LABEL]         = 1;
SAVE_ON_COPY[SURF_LABEL]         = 1;       
BRIEF_HELP[SURF_LABEL] = "Label the faces of a volume mapping that coincide with the "
"physical surface. "      
"The surface label is used to identify the different parts of the boundary of the "
"computational domain. All faces of all component grids that are inside "
"or outside of the computational domain must have a zero surface-label "
"(which is the default value). Faces that are aligned with "
"the physical boundary must be assigned a positive surface label. The "
"value of the surface-label should be the same for all faces of all "
"grids that are aligned with the same part of the boundary. Hence, if the "
"domain is simply connected, only one surface-label value should be used. For "
"a doubly connected domain, two values should be used, and so on.";

COMMAND[EDGE_LABEL] = "edge-label";       
/*file[EDGE_LABEL] = "edge_label_com.h"; */
ARGUMENT[EDGE_LABEL]         = 1;
SAVE_ON_COPY[EDGE_LABEL]         = 1;       
BRIEF_HELP[EDGE_LABEL] = "Mark the edges of a volume mapping that should cut holes "
"in the faces of other volume mappings during the hole-cutting part of the overlap "
"algorithm. For example, "
"consider a wing-fuselage configuration with one grid around the fuselage and "
"one O-grid around the wing. In this case the edge in the wing grid that coincides "
"with both the fuselage and the wing surfaces should have a positive edge label, "
"so that all points inside the wing on the fuselage can be removed. If the hole is "
"described by several edges in several volume mappings, "
"the same positive edge-label should be used for all those edges. ";

COMMAND[CUT_HOLES] = "cut-holes";       
ARGUMENT[CUT_HOLES]         = 1;
SAVE_ON_COPY[CUT_HOLES]         = 1;       
BRIEF_HELP[CUT_HOLES] = "Toggle the flag that determines wether holes should be cut "
"out from this grid during the overlapping grid algorithm. The default is to cut "
"holes, but the algorithm will become faster if the hole-cutting is disabled for "
"grids that are known to be fully inside the computational domain.";

COMMAND[NEW_VOLUME] = "make-volume-mapping";       
/*file[NEW_VOLUME] = "choose_volume_com.h"; */
ARGUMENT[NEW_VOLUME]    = 1;
SAVE_ON_COPY[NEW_VOLUME]    = 1;
BRIEF_HELP[NEW_VOLUME] = "Define a new volume mapping. First choose a name "
"for the mapping, and then select the specific type of mapping. "
"Thereafter, you proceed "
"directly to modify the parameters for that type of mapping, in the same way as "
"you would do with the command `change-volume-mapping'.";

COMMAND[CHECK_VOLUME] = "check-volume-mapping";       
ARGUMENT[CHECK_VOLUME]    = 1;
SAVE_ON_COPY[CHECK_VOLUME]    = 1;
BRIEF_HELP[CHECK_VOLUME] = "Compare the analytic jacobian of a volume mapping to "
"a centered difference approximation.";

COMMAND[CHANGE_VOLUME] = "change-volume-mapping";       
ARGUMENT[CHANGE_VOLUME] = 1;
SAVE_ON_COPY[CHANGE_VOLUME] = 1;
BRIEF_HELP[CHANGE_VOLUME] = "Modify an existing volume mapping. This command gives "
"you the same options for modifying a volume mapping as are provided by the command "
"`make-volume-mapping', with the obvious exception that you do not specify the "
"type of mapping again.";

COMMAND[TRANSFORM_VOLUME] = "transform-volume-mapping";       
/*file[TRANSFORM_VOLUME] = "transform_volume_com.h"; */
ARGUMENT[TRANSFORM_VOLUME] = 1;
SAVE_ON_COPY[TRANSFORM_VOLUME] = 1;
BRIEF_HELP[TRANSFORM_VOLUME] = "Add a translation and/or a rotation to the mapping";

COMMAND[FILE_VOLUME] = "save-volume-mapping";       
ARGUMENT[FILE_VOLUME] = 1;
SAVE_ON_COPY[FILE_VOLUME] = 1;
BRIEF_HELP[FILE_VOLUME] = "Save the coordinates of the grid points in a volume "
"mapping on a PLOT3D ASCII file.";

COMMAND[DELETE_VOLUME] = "delete-volume-mapping";       
ARGUMENT[DELETE_VOLUME] = 1;
SAVE_ON_COPY[DELETE_VOLUME] = 1;
BRIEF_HELP[DELETE_VOLUME] = "Remove an existing volume mapping.";

COMMAND[OVERTURE] = "overture-mappings";       
ARGUMENT[OVERTURE]    = 0;
SAVE_ON_COPY[OVERTURE]    = 1;
BRIEF_HELP[OVERTURE] = "Open an overture window and manipulate the overture mappings.";

COMMAND[NEW_3D_OVERLAP] = "make-overlapping-grid";
/*file[NEW_3D_OVERLAP] = "set_3d_overlap_com.h"; */
ARGUMENT[NEW_3D_OVERLAP] = 1; 
SAVE_ON_COPY[NEW_3D_OVERLAP] = 1;
BRIEF_HELP[NEW_3D_OVERLAP]= "Create and initialize a new overlapping grid structure. "
"Before the overlapping "
"grid is computed, you will be given the opportunity to modify the "
"controlling parameters, i.e., the list of member grids, "
"the widths of the discretization and interpolation stencils, the interpolation "
"type, etc.";

COMMAND[CHANGE_3D_OVERLAP] = "change-overlapping-grid";
ARGUMENT[CHANGE_3D_OVERLAP] = 1;
SAVE_ON_COPY[CHANGE_3D_OVERLAP] = 1;
BRIEF_HELP[CHANGE_3D_OVERLAP]= "Modify an existing overlapping grid structure. "
"You will be given the same possibilities to change the parameters as with the "
"command `make-overlapping-grid'. In particular, you will be able to modify the "
"list of member grids, the widths of the discretization and interpolation stencils, "
"and the interpolation type.";

COMMAND[DELETE_3D_OVERLAP] = "delete-overlapping-grid";
ARGUMENT[DELETE_3D_OVERLAP] = 1;
SAVE_ON_COPY[DELETE_3D_OVERLAP] = 1;
BRIEF_HELP[DELETE_3D_OVERLAP]= "Delete an overlapping grid structure.";

COMMAND[SHOW] = "show-entities";       
ARGUMENT[SHOW]         = 0;
SAVE_ON_COPY[SHOW]              = 1;
BRIEF_HELP[SHOW] = "List all presently defined surface mappings, volume mappings and "
"overlapping grids.";

COMMAND[SAVE_OVERLAP] = "save-overlapping-grid";       
ARGUMENT[SAVE_OVERLAP] = 1;
SAVE_ON_COPY[SAVE_OVERLAP]      = 1;
BRIEF_HELP[SAVE_OVERLAP] = "Save a valid overlapping grid on a HDF "
"database file. See the Latex document `doc/data_base.tex' for a description of "
"the format. You might find the routines in `util/hdf_stuff.c' (util/hdf_stuff.h) "
"and `util/c_array.c' (util/c_array.h) useful for reading the data.";

COMMAND[TUTORIAL] = "tutorial";       
/*file[TUTORIAL] = "run_tutorial_com.h"; */
ARGUMENT[TUTORIAL]     = 0;
SAVE_ON_COPY[TUTORIAL]          = 1;
BRIEF_HELP[TUTORIAL] = "Get an introduction to Chalmesh and overlapping grids "
"by running commented command files that take you through the steps of the grid "
"generation process.";

COMMAND[READ_COMMAND] = "read-commands";
ARGUMENT[READ_COMMAND] = 1;
SAVE_ON_COPY[READ_COMMAND]       = 0;
BRIEF_HELP[READ_COMMAND] = "Start reading the commands from a command file "
"instead of standard input.";

COMMAND[START_SAVE] = "start-saving-commands";
ARGUMENT[START_SAVE]   = 1;
SAVE_ON_COPY[START_SAVE]  = 1;
BRIEF_HELP[START_SAVE] = "Open a new log file and start saving a copy of "
"the commands on that file.";

COMMAND[STOP_SAVE] = "stop-saving-commands";
ARGUMENT[STOP_SAVE]    = 0; 
SAVE_ON_COPY[STOP_SAVE]   = 0; 
BRIEF_HELP[STOP_SAVE] = "Close an open log file and stop saving a copy of the "
"commands on that file.";

COMMAND[RESET_ALL]  = "reset-program";
ARGUMENT[RESET_ALL]    = 0;
SAVE_ON_COPY[RESET_ALL]    = 1;
BRIEF_HELP[RESET_ALL]  = "Delete all surface mappings, volume mappings, and overlapping "
"grids. Also reset the view-point for the graphics windows.";

COMMAND[GRAPHICS_ON] = "open-graphics";
ARGUMENT[GRAPHICS_ON]  = 0;
SAVE_ON_COPY[GRAPHICS_ON] = 1;
BRIEF_HELP[GRAPHICS_ON]  = "Open all graphics windows and replot their contents."
"The windows may later be closed with the command close-graphics.";

COMMAND[GRAPHICS_OFF] = "close-graphics";
ARGUMENT[GRAPHICS_OFF] = 0;
SAVE_ON_COPY[GRAPHICS_OFF] = 1;
BRIEF_HELP[GRAPHICS_OFF] = "Close all graphics windows. This will speed up the program "
"if it is run on a machine with slow graphics. The graphics window may later be "
"re-opened with the command `open-graphics'.";

COMMAND[HELP]   = "help";
ARGUMENT[HELP]         = 0;
SAVE_ON_COPY[HELP]         = 1;
BRIEF_HELP[HELP]   = NULL;

COMMAND[QUIT]   = "quit";
ARGUMENT[QUIT]         = 0;
SAVE_ON_COPY[QUIT]         = 0;
BRIEF_HELP[QUIT]   = "Terminate the Chalmesh session. Note that if you computed "
"an overlapping grid during this session, you must do `save-overlapping-grid' to "
"use the grid outside of the program.";

