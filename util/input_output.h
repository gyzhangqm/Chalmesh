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
#ifndef input_output_h
#define input_output_h

/* data structure for i/o files */
typedef struct input_output{
  char *copy_file_name;
  FILE *read_command, *copy_command;
} input_output;

typedef struct Color {
  double r, g, b;
} color_type;

/* Marker types for X_marker */

#define DOT                 0
#define CIRCLE              1
#define BOX                 2
#define ASTERISK            3
#define CROSS               4
#define FILLED_CIRCLE       5
#define FILLED_BOX          6
#define TRIANGLE            7
#define FILLED_TRIANGLE     8
#define INV_TRIANGLE        9
#define FILLED_INV_TRIANGLE 10

#endif
