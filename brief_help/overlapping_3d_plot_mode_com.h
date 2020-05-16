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
enum {BOUNDARY_LINES, CHANGE_SURFACE, USED_BOUNDARY_FACES, USED_CELLS, USED_FACES, 
	POSITIVE_BC, WORLD_POINTS, EDGE_LABEL, 
	USED_GRIDLINES, INTERP_FACES, INTERP_POINTS, BAD_POINTS,
	DEAD_POINTS, CLIP, NO_CLIP, 
	SHOW_ONE_COMPONENT, ADD_ONE_COMPONENT, REMOVE_ONE_COMPONENT, 
	SHOW_ALL_COMPONENTS, REMOVE_ALL_COMPONENTS, COORD_CAGE, STD_VIEW, SHOW_PLOT_MODE,
	SPHERE_SIZE,
	HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[26], *BRIEF_HELP[26];
int ARGUMENT[26], *SAVE_ON_COPY=NULL;
const int LEVEL=2, NO_INDENT=0;


COMMAND[SPHERE_SIZE] = "sphere-size";
ARGUMENT[SPHERE_SIZE] = 1;
BRIEF_HELP[SPHERE_SIZE] = "Set the relative size of the spheres used to represent "
"interpolation points and bad points.";

COMMAND[SHOW_PLOT_MODE] = "show-plot-mode";
ARGUMENT[SHOW_PLOT_MODE] = 0;
BRIEF_HELP[SHOW_PLOT_MODE] = "List the current plot_mode.";

COMMAND[WORLD_POINTS] = "world-points";
ARGUMENT[WORLD_POINTS] = 0;
BRIEF_HELP[WORLD_POINTS] = "Toggle the plotting of the world points (intermediate "
"result of the hole-cutting algorithm).";

COMMAND[EDGE_LABEL] = "edge-label";
ARGUMENT[EDGE_LABEL] = 0;
BRIEF_HELP[EDGE_LABEL] = "Toggle the plotting of the curves with a non-zero "
"edge label.";

COMMAND[BOUNDARY_LINES] = "highlight-boundaries";
COMMAND[CHANGE_SURFACE] = "change-surface";
COMMAND[USED_BOUNDARY_FACES]  = "used-boundary-faces";
COMMAND[USED_CELLS]  = "used-cells";
COMMAND[USED_FACES]  = "discretization-faces";
COMMAND[POSITIVE_BC]  = "physical-boundary-faces";
COMMAND[USED_GRIDLINES]  = "used-grid-lines";
COMMAND[SHOW_ONE_COMPONENT] = "show-one-component";
COMMAND[ADD_ONE_COMPONENT] = "add-one-component";
COMMAND[REMOVE_ONE_COMPONENT] = "remove-one-component";
COMMAND[SHOW_ALL_COMPONENTS] = "show-all-components";
COMMAND[REMOVE_ALL_COMPONENTS] = "remove-all-components";
COMMAND[COORD_CAGE]  = "coord-cage";
COMMAND[INTERP_FACES] = "interpolation-faces";
COMMAND[INTERP_POINTS] = "interpolation-points";
COMMAND[BAD_POINTS]  = "bad-grid-points";
/*COMMAND[DUMMY_POINTS]  = "dummy-interpolation-points";*/
COMMAND[DEAD_POINTS]  = "dead-grid-points";
COMMAND[CLIP]        = "clip-plane";
COMMAND[NO_CLIP]     = "no-clip-plane";
COMMAND[STD_VIEW]    = "standard-view";
COMMAND[HELP]        = "help";
COMMAND[EXIT]        = "exit";

ARGUMENT[BOUNDARY_LINES] = 0;
ARGUMENT[CHANGE_SURFACE] = 0;
ARGUMENT[USED_BOUNDARY_FACES]  = 0;
ARGUMENT[USED_CELLS]  = 0;
ARGUMENT[USED_FACES]  = 0;
ARGUMENT[POSITIVE_BC]  = 0;
ARGUMENT[USED_GRIDLINES]  = 0;
ARGUMENT[SHOW_ONE_COMPONENT] = 1;
ARGUMENT[ADD_ONE_COMPONENT] = 1;
ARGUMENT[REMOVE_ONE_COMPONENT] = 1;
ARGUMENT[SHOW_ALL_COMPONENTS] = 0;
ARGUMENT[REMOVE_ALL_COMPONENTS] = 0;
ARGUMENT[COORD_CAGE]  = 0;
ARGUMENT[INTERP_FACES] = 0;
ARGUMENT[INTERP_POINTS] = 0;
ARGUMENT[BAD_POINTS]  = 0;
/*ARGUMENT[DUMMY_POINTS]  = 0;*/
ARGUMENT[DEAD_POINTS]  = 0;
ARGUMENT[CLIP]        = 1;
ARGUMENT[NO_CLIP]     = 0;
ARGUMENT[STD_VIEW] = 0;
ARGUMENT[HELP]        = 1;
ARGUMENT[EXIT]        = 0;

BRIEF_HELP[CHANGE_SURFACE] = "Toggle the change-surface";
BRIEF_HELP[BOUNDARY_LINES] = "Toggle the plotting of the boundary grid lines";
BRIEF_HELP[USED_BOUNDARY_FACES] = "Toggle the boundaries of the component grids.";
BRIEF_HELP[USED_CELLS] = "Toggle the plotting of all used grid-cells";
BRIEF_HELP[USED_FACES] = "Toggle the plotting of all discretization faces";
BRIEF_HELP[POSITIVE_BC] = "Toggle the plotting of the faces with positive "
"boundary condition.";
BRIEF_HELP[USED_GRIDLINES] = "Toggle the plotting of all used grid-lines";
BRIEF_HELP[SHOW_ONE_COMPONENT] = "Specify which single component grid to plot.";
BRIEF_HELP[ADD_ONE_COMPONENT] = "Add another component grid to the plot.";
BRIEF_HELP[REMOVE_ONE_COMPONENT] = "Remove on component grid from the plot.";
BRIEF_HELP[SHOW_ALL_COMPONENTS] = "Show all component grids.";
BRIEF_HELP[REMOVE_ALL_COMPONENTS] = "Don't show any component grids.";
BRIEF_HELP[COORD_CAGE]  = "Toggle the coordinate axes.";
BRIEF_HELP[INTERP_FACES] = "Toggle the plotting of all interpolation faces.";
BRIEF_HELP[INTERP_POINTS] = "Toggle the plotting of all interpolation points.";
BRIEF_HELP[BAD_POINTS]  = "Toggle the plotting of all bad grid points.";
/*BRIEF_HELP[DUMMY_POINTS]  = "Toggle the plotting of all interpolation points "*/
/*"with dummy donors";*/
BRIEF_HELP[DEAD_POINTS]  = "Toggle the plotting of all dead (unused) grid points.";
BRIEF_HELP[CLIP]      = "Turn clipping on and specify a clip plane.";
BRIEF_HELP[NO_CLIP]   = "Turn clipping off.";
BRIEF_HELP[STD_VIEW]  = "Reset the view point and focal point to show everything in "
"the bounding box.";
BRIEF_HELP[HELP]      = NULL;
BRIEF_HELP[EXIT]      = "Stop modifying the plotting parameters and exit to the "
"previous command level.";


