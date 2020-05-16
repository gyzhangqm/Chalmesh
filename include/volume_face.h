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
#ifndef volume_face_h
#define volume_face_h

/* the correspondence between the directions in the mother grid and the component */
/* surface grid are as follows */
/* dir  r_dim   s_dim  */
/*  1   r2_dim  r3_dim */
/*  2   r3_dim  r1_dim */
/*  3   r1_dim  r2_dim */

typedef struct{
  volume_mapping *volume;
  int side, dir;
} volume_face_data;

/* prototypes */
volume_face_data *
init_volume_face(surface_mapping *surface, volume_mapping *volume);


#endif
