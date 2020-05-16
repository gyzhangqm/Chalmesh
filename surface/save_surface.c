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
#include "surface_internal.h"

void
save_surface_mapping( surface_mapping *surface, FILE *fp ){
  int i, j;
  real dr, ds;
  grid_point gp;

  dr = 1.0/ (surface->r_points - 1);
  ds = 1.0/ (surface->s_points - 1);

/* formatted plot3d file with one grid in the file */
  fprintf( fp, "%i %i %i\n", surface->r_points, surface->s_points, 1 );

/* save the x coordinates */
  for (j = 1; j <= surface->s_points; j++){
    for (i = 1; i <= surface->r_points; i++){
      gp.r = (i-1)*dr; gp.s = (j-1)*ds;
      forward_surface_mapping( &gp, surface );
      fprintf( fp, "%.18G\n", gp.x);
    }
  }

/* save the y coordinates */
  for (j = 1; j <= surface->s_points; j++){
    for (i = 1; i <= surface->r_points; i++){
      gp.r = (i-1)*dr; gp.s = (j-1)*ds;
      forward_surface_mapping( &gp, surface );
      fprintf( fp, "%.18G\n", gp.y);
    }
  }

/* save the z coordinates */
  for (j = 1; j <= surface->s_points; j++){
    for (i = 1; i <= surface->r_points; i++){
      gp.r = (i-1)*dr; gp.s = (j-1)*ds;
      forward_surface_mapping( &gp, surface );
      fprintf( fp, "%.18G\n", gp.z);
    }
  }

} /* end save_surface_mapping */
