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
enum {CORNER_W, DISC_W, EXTRA, EXPLICIT, IMPLICIT, INTERP_W, NORMAL_W, TANG_W, ALL_W, 
      NORMAL_DIST, SHOW, HELP, EXIT};
const int LAST_COM = EXIT, LEVEL=3, NO_INDENT=0;
char *COMMAND[13], *BRIEF_HELP[13];
int ARGUMENT[13], *SAVE_ON_COPY=NULL;

COMMAND[CORNER_W] ="corner-width";         
ARGUMENT[CORNER_W] = 1;
BRIEF_HELP[CORNER_W] = "Change the corner width. The corner width equals the width of "
"the discretization stencil at a corner point in the grid and at points so close "
"to the corner that the regular discretization stencil can't be used.";

COMMAND[DISC_W] ="discretization-width"; 
ARGUMENT[DISC_W] = 1;
BRIEF_HELP[DISC_W] = "Change the discretization width, which is defined as the "
"number of grid points in each grid direction in the discretization "
"stencil, away from boundaries.";

COMMAND[EXTRA] ="ghostpoints";          
ARGUMENT[EXTRA] = 1;
BRIEF_HELP[EXTRA] = "Change the number of extra grid points outside the boundary of "
"the grid. These points have no influence on the overlapping grid algorithm, but are "
"sometimes useful for discretizing a PDE close to a boundary.";

COMMAND[EXPLICIT] ="explicit-interpolation"; 
ARGUMENT[EXPLICIT] = 0;
BRIEF_HELP[EXPLICIT] = "Use explicit interpolation in the overlapping grid. When the "
"interpolation is explicit, the overlap between the grids is "
"sufficiently wide to ensure that all donor points are discretization points "
"(and not interpolation points).";

COMMAND[IMPLICIT] ="implicit-interpolation"; 
ARGUMENT[IMPLICIT] = 0;
BRIEF_HELP[IMPLICIT] = "Use implicit interpolation in the overlapping grid. When the "
"interpolation is implicit, the interpolation points are allowed to interpolate both "
"from discretization and interpolation points in the donor grid.";

COMMAND[INTERP_W] ="interpolation-width";  
ARGUMENT[INTERP_W] = 1;
BRIEF_HELP[INTERP_W] = "Change the width of the interpolation formula, which is "
"defined as the maximum number number of grid points in any coordinate direction "
"of the interpolation stencil.";

COMMAND[NORMAL_W] ="normal-width";         
ARGUMENT[NORMAL_W] = 1;
BRIEF_HELP[NORMAL_W] = "Change the width of the discretization stencil at the boundary, "
"normal to the boundary.";

COMMAND[TANG_W] ="tangential-width";     
ARGUMENT[TANG_W] = 1;
BRIEF_HELP[TANG_W] = "Change the width of the discretization stencil at the boundary, "
"tangent to the boundary.";

COMMAND[ALL_W] ="set-all-widths";       
ARGUMENT[ALL_W] = 1;
BRIEF_HELP[ALL_W] = "Set all widths at once, i.e. the discretization width, the "
"normal width, the tangential width, the corner width and the interpolation width.";

COMMAND[NORMAL_DIST] ="extra-misfit-tolerance";
ARGUMENT[NORMAL_DIST] = 1; 

BRIEF_HELP[NORMAL_DIST] = "Set the extra tolerance for surface and edge "
"missmatch that will be added to the estimated value. The surface tolerance "
"is used to allow an interpolation "
"point on a physical boundary in one grid to interpolate from a donor "
"grid, eventough the interpolation point is slightly outside of the "
"donor grid.";

COMMAND[SHOW] ="show";
ARGUMENT[SHOW] = 0;
BRIEF_HELP[SHOW] = "Show the present value of the overlap parameters.";

COMMAND[HELP] ="help";
ARGUMENT[HELP] = 0;
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT] ="exit";
ARGUMENT[EXIT] = 0;
BRIEF_HELP[EXIT] = "Stop changing the overlap parameters and exit from this command "
"level.";
