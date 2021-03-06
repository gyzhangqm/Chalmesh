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
#ifndef sphere_patch_h
#define sphere_patch_h

typedef struct {
  int polar_axis;
  real x_semi, y_semi, z_semi, alpha_min, alpha_max, beta_min, beta_max;
  real x_center, y_center, z_center;
} sphere_patch_data;

/* prototypes */
sphere_patch_data *
init_sphere_patch(surface_mapping *surface);
void 
set_sphere_patch( input_output *io_ptr, surface_mapping *surface );
int
sphere_patch( grid_point *gp_ptr, surface_mapping *surface );

#endif
