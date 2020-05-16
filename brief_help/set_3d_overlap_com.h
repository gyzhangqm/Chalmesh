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
enum {LIST, PARAMS, COMPUTE, SHOW, PLOT_MODE, INSPECT_OVER_SURF, TRIM_STYLE, 
      VERBOSE, SHOW_HOLES, PHYSICAL_BNDRY, MISMATCH, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[13], *BRIEF_HELP[13];
int ARGUMENT[13], *SAVE_ON_COPY=NULL;
const int LEVEL=1;

COMMAND[SHOW_HOLES] ="show-holes";
ARGUMENT[SHOW_HOLES] = 1;
BRIEF_HELP[SHOW_HOLES] = "Toggle the flag that determines whether the "
"used grid points in the overlapping grid should be shown after the "
"hole-cutting has been applied. If the overlap algorithm fails, it is often "
"possible to determine if it failed because of an inconsistently defined boundary, "
"or because of too few grid points in the overlap region, "
"by looking at the grid after the holes have been cut. If holes have been cut in "
"the wrong places, it is probably because some sides of some grids have been "
"given an inconsistent surface-label. But if the overlap algorithm fails "
"eventhough the hole are cut correctly, the problem is most likely "
"caused by too few grid points in the overlap region.";

COMMAND[PHYSICAL_BNDRY] ="show-physical-boundaries";
ARGUMENT[PHYSICAL_BNDRY] = 1;
BRIEF_HELP[PHYSICAL_BNDRY] = "Toggle the flag that determines if the global "
"boundary of the computational domain will be shown during the overlap algorithm.";

COMMAND[MISMATCH] ="mismatch-correction";
ARGUMENT[MISMATCH] = 1;
BRIEF_HELP[MISMATCH] = "Toggle the flag that determines whether the interpolation "
"locations should be compensated for boundary mismatch.";

COMMAND[PARAMS] = "overlap-parameters";
/*file[PARAMS] = "overlap_3d_parameters_com.h"; */
ARGUMENT[PARAMS] = 0;
BRIEF_HELP[PARAMS]  = "Set the overlap parameters, i.e., the widths of the "
"discretization and interpolation stencils.";

COMMAND[LIST]   = "change-member-list";       
/*file[LIST] = "change_3d_list_com.h"; */
ARGUMENT[LIST]   = 0;
BRIEF_HELP[LIST]    = "Change the priority of the components in the "
"overlapping grid. Also insert or remove components from the list.";

COMMAND[COMPUTE]   = "compute-overlap";       
ARGUMENT[COMPUTE]   = 0;
BRIEF_HELP[COMPUTE] = "Compute the overlapping grid with the present list of "
"component grids and the present value of the overlap parameters.";

COMMAND[SHOW]      = "show-parameters";  
ARGUMENT[SHOW]      = 0;
BRIEF_HELP[SHOW]    = "Show the list of component grids and the value "
"of the overlap parameters.";

COMMAND[PLOT_MODE] = "overlapping-plot-mode";  
/*file[PLOT_MODE] = "overlapping_3d_plot_mode_com.h"; */
ARGUMENT[PLOT_MODE] = 0;
BRIEF_HELP[PLOT_MODE] = "Modify the graphical contents of the "
"`Overlapping grids' window.";

COMMAND[INSPECT_OVER_SURF] = "inspect-surface-grid";  
/*file[INSPECT_OVER_SURF] = "hs_plot_mode_com.h"; */
ARGUMENT[INSPECT_OVER_SURF] = 1;
BRIEF_HELP[INSPECT_OVER_SURF] = "Look at the non-overlapping hybrid grid "
"corresponding to one "
"physical boundary. The hybrid grid consists of structured non-overlapping "
"components joined by unstructured triangles.";

COMMAND[TRIM_STYLE] = "trimming-style";     
ARGUMENT[TRIM_STYLE] = 0;
BRIEF_HELP[TRIM_STYLE] = "Toggle the trimming style between minimizing (default) and "
"maximizing the overlap in the overlapping grid.";

COMMAND[VERBOSE] = "verbose";     
ARGUMENT[VERBOSE] = 1;
BRIEF_HELP[VERBOSE] = "The verbose value controls the amount of information "
"that is printed during the overlap algorithm. A zero value gives no information; "
"non-zero values have the following meaning:\n"
"verbose     comment\n"
" 1          timing\n"
" 2          +stage of the algorithm\n"
" 3          +details of the overlap algorithm (debug info)\n";

COMMAND[HELP]      = "help";
ARGUMENT[HELP]      = 0;
BRIEF_HELP[HELP]      = NULL;

COMMAND[EXIT]      = "exit";
ARGUMENT[EXIT]      = 0;
BRIEF_HELP[EXIT]      = "Stop modifying the parameters and exit to the previous "
"command level.";




