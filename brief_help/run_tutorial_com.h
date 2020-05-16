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
enum {INTRO, INTERACT, STRETCHED_ELLIPSOID, HALF_ELLIPSOID, QUARTER_ELLIPSOID,
      PROJECTION, EXTERNAL_GRIDS, HYPERBOLIC, HELP, EXIT};

const int LAST_COM = EXIT, LEVEL=1;
char *COMMAND[10], *BRIEF_HELP[10];
int *ARGUMENT = NULL, *SAVE_ON_COPY = NULL;

COMMAND[INTRO] = "overlapping-grid-steps";
BRIEF_HELP[INTRO] = "Present a brief outline of how Chalmesh is used to make "
"overlapping grids.";

COMMAND[INTERACT] = "command-interpreter-introduction";
BRIEF_HELP[INTERACT] = "Describes how 3-D objects are translated and rotated on "
"the screen, and shows the basic features of the command interpreter.";

COMMAND[STRETCHED_ELLIPSOID] = "stretched-ellipsoid";
BRIEF_HELP[STRETCHED_ELLIPSOID] = "Create an overlapping grid around a stretched "
"ellipsoid inside a box. This is the most basic type of overlapping grid.";

COMMAND[PROJECTION] = "projection";
BRIEF_HELP[PROJECTION] = "Create an overlapping grid around an ellipsoid sitting "
"on a cylinder. This command file demonstrates how to use the hyperbolic grid "
"generator and how to project grid faces onto curved surfaces.";

COMMAND[HALF_ELLIPSOID] = "half-ellipsoid";
BRIEF_HELP[HALF_ELLIPSOID] = "Create an overlapping grid around half an "
"ellipsoid inside a box. This example demonstrates the use of edge labels.";

COMMAND[QUARTER_ELLIPSOID] = "quarter-ellipsoid";
BRIEF_HELP[QUARTER_ELLIPSOID] = "Create an overlapping grid around one quarter of an "
"ellipsoid inside a box. This example shows how several edge labels can be used.";

COMMAND[EXTERNAL_GRIDS] = "external-grids";
BRIEF_HELP[EXTERNAL_GRIDS] = "Make an overlapping grid between a sphere and a box. "
"In this example, the component grids are read from file in PLOT3D ASCII format";

COMMAND[HYPERBOLIC] = "hyperbolic";
BRIEF_HELP[HYPERBOLIC] = "Demonstrate the hyperbolic grid generator.";

COMMAND[HELP] = "help";
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT] = "exit";
BRIEF_HELP[EXIT] = "Exit from the tutorial mode and return to the top command "
"level.";

