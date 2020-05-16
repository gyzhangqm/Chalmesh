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
enum {GRID_POINTS, BOUNDARIES, ARROWS, SURF_LABEL, BOUNDARY_CONDITION,
      EDGE_LABEL, CLIP, NO_CLIP,
      R1_SURFACE, SPECIFY_R1, INCREASE_R1, DECREASE_R1, 
      R2_SURFACE, SPECIFY_R2, INCREASE_R2, DECREASE_R2, 
      R3_SURFACE, SPECIFY_R3, INCREASE_R3, DECREASE_R3, 
      CONST_SURF_LABEL, NO_CONST_SURF_LABEL, CONST_BC, NO_CONST_BC,
      COORD_CAGE, STD_VIEW, HELP, EXIT};
const int LAST_COM = EXIT;

char *COMMAND[28], *BRIEF_HELP[28];
int ARGUMENT[28], *SAVE_ON_COPY=NULL;
const int LEVEL=3, NO_INDENT=0;

COMMAND[GRID_POINTS]  = "grid-points";
ARGUMENT[GRID_POINTS] = 0;
BRIEF_HELP[GRID_POINTS] = "Toggle the grid points.";

COMMAND[BOUNDARIES]   = "highlight-boundaries";
ARGUMENT[BOUNDARIES] = 0;
BRIEF_HELP[BOUNDARIES] = "Toggle highlighting of the boundaries.";

COMMAND[ARROWS]       = "parameter-arrows";
ARGUMENT[ARROWS] = 0;
BRIEF_HELP[ARROWS] = "Toggle arrows showing the parameter directions.";

COMMAND[BOUNDARY_CONDITION] = "boundary-condition";
ARGUMENT[BOUNDARY_CONDITION] = 0;
BRIEF_HELP[BOUNDARY_CONDITION] = "Toggle the boundary conditions.";

COMMAND[SURF_LABEL] = "surface-label";
ARGUMENT[SURF_LABEL] = 0;
BRIEF_HELP[SURF_LABEL] = "Toggle the surface labels.";

COMMAND[EDGE_LABEL] = "edge-label";
ARGUMENT[EDGE_LABEL] = 0;
BRIEF_HELP[EDGE_LABEL] = "Toggle the edge labels.";

COMMAND[CLIP]        = "clip-plane";
ARGUMENT[CLIP]        = 1;
BRIEF_HELP[CLIP]      = "Turn clipping on and specify up to four clip planes.";

COMMAND[NO_CLIP]     = "no-clip-plane";
ARGUMENT[NO_CLIP]     = 0;
BRIEF_HELP[NO_CLIP]   = "Turn clipping off.";

COMMAND[R1_SURFACE]   = "r1-grid-lines";           
ARGUMENT[R1_SURFACE] = 0;
BRIEF_HELP[R1_SURFACE]  = "Toggle plotting of grid lines on the grid surface r1=const. The constant is set with the command `specify-r1-surface'.";

COMMAND[SPECIFY_R1]   = "specify-r1-surface";      
ARGUMENT[SPECIFY_R1] = 1;
BRIEF_HELP[SPECIFY_R1]   = "Specify which r1=const grid surface to plot.";

COMMAND[INCREASE_R1]  = "increase-r1-constant";    
ARGUMENT[INCREASE_R1] = 0;
BRIEF_HELP[INCREASE_R1]   = "Move the r1=const to the next grid surface.";

COMMAND[DECREASE_R1]  = "decrease-r1-constant";    
ARGUMENT[DECREASE_R1] = 0;
BRIEF_HELP[DECREASE_R1]   = "Move the r1=const to the previous grid surface.";

COMMAND[R2_SURFACE]   = "r2-grid-lines";           
ARGUMENT[R2_SURFACE] = 0;
BRIEF_HELP[R2_SURFACE]  = "Toggle plotting of grid lines on the grid surface r2=const. The constant is set with the command `specify-r2-surface'.";

COMMAND[SPECIFY_R2]   = "specify-r2-surface";      
ARGUMENT[SPECIFY_R2] = 1;
BRIEF_HELP[SPECIFY_R2]   = "Specify which r2=const grid surface to plot.";

COMMAND[INCREASE_R2]  = "increase-r2-constant";    
ARGUMENT[INCREASE_R2] = 0;
BRIEF_HELP[INCREASE_R2]   = "Move the r2=const to the next grid surface.";

COMMAND[DECREASE_R2]  = "decrease-r2-constant";    
ARGUMENT[DECREASE_R2] = 0;
BRIEF_HELP[DECREASE_R2]   = "Move the r2=const to the previous grid surface.";

COMMAND[R3_SURFACE]   = "r3-grid-lines";           
ARGUMENT[R3_SURFACE] = 0;
BRIEF_HELP[R3_SURFACE]  = "Toggle plotting of grid lines on the grid surface r3=const. The constant is set with the command `specify-r3-surface'.";

COMMAND[SPECIFY_R3]   = "specify-r3-surface";      
ARGUMENT[SPECIFY_R3] = 1;
BRIEF_HELP[SPECIFY_R3]   = "Specify which r3=const grid surface to plot.";

COMMAND[INCREASE_R3]  = "increase-r3-constant";    
ARGUMENT[INCREASE_R3] = 0;
BRIEF_HELP[INCREASE_R3]   = "Move the r3=const to the next grid surface.";

COMMAND[DECREASE_R3]  = "decrease-r3-constant";    
ARGUMENT[DECREASE_R3] = 0;
BRIEF_HELP[DECREASE_R3]   = "Move the r3=const to the previous grid surface.";

COMMAND[CONST_SURF_LABEL] = "color-surface-label";     
ARGUMENT[CONST_SURF_LABEL] = 1;
BRIEF_HELP[CONST_SURF_LABEL] = "Plot all grid faces corresponding to a specified "
"surface label value.";

COMMAND[NO_CONST_SURF_LABEL] = "uncolor-surface-label";
ARGUMENT[NO_CONST_SURF_LABEL] = 0;
BRIEF_HELP[NO_CONST_SURF_LABEL] = "Turn off the plotting of grid faces corresponding "
"to a specified surface label value.";

COMMAND[CONST_BC]    = "color-boundary-condition";     
ARGUMENT[CONST_BC] = 1;
BRIEF_HELP[CONST_BC]    = "Plot all grid faces corresponding to a specified "
"boundary condition value.";

COMMAND[NO_CONST_BC] = "uncolor-boundary-condition";   
ARGUMENT[NO_CONST_BC] = 0;
BRIEF_HELP[NO_CONST_BC]    = "Turn off the plotting of all grid faces corresponding "
"to a specified boundary condition value.";

COMMAND[COORD_CAGE]  = "coord-cage";                   
ARGUMENT[COORD_CAGE] = 0;
BRIEF_HELP[COORD_CAGE]  = "Toggle the coordinate axes.";

COMMAND[STD_VIEW]    = "standard-view";            
ARGUMENT[STD_VIEW] = 0;
BRIEF_HELP[STD_VIEW]  = "Reset the view point and focal point to show everything in "
"the bounding box.";

COMMAND[HELP]        = "help";                     
ARGUMENT[HELP] = 1;
BRIEF_HELP[HELP]      = NULL;

COMMAND[EXIT]        = "exit";                     
ARGUMENT[EXIT] = 0;
BRIEF_HELP[EXIT]      = "Stop modifying the plotting parameters and exit to the "
"previous command level.";


