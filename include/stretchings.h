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
#ifndef stretchings_h
#define stretchings_h

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>

#include "version.h"
#include <stupid_compiler.h>
#include <real.h>
#include <c_array.h>
#include "linked_list.h"
#include <nrutil.h>
#include <min_max.h>
#include <ogl_plot.h>
#include <input_output.h>
#include <linalg.h>

#ifdef HDF_DUMMY
  #include "hdf_dummy.h"
#else
  #include "hdf_stuff.h"
#endif

typedef struct generic_stretching{

  int stretch_type;
/* stretch_type: */
/* 1 : arclength stretching */
/* 2 : exponential stretching */
/* 3 : layer stretching */
/* 4 : hyperbolic tangent stretching */

  void *stretch_data_ptr;
  void 
    (*uniform_to_stretch)(real uniform, struct generic_stretching *stretch_ptr, 
			  real *stretch, real *d_stretch, real *d2_stretch );
  void 
    (*stretch_to_uniform)(real stretch, struct generic_stretching *stretch_ptr, 
		   real *uniform, real *d_uniform, real *d2_uniform );
  void 
    (*write_stretch_data)( int32 dir, struct generic_stretching *stretch_ptr );
  void 
    (*set_stretching)(input_output *io_ptr, struct generic_stretching *stretch_ptr,
		      int *grid_points);
  void 
    (*cleanup_stretching)(struct generic_stretching *stretch_ptr );
  void *
    (*copy_stretching)( void *stretch_data_ptr );
} generic_stretching;

/* prototypes */
generic_stretching *
choose_stretching( input_output *io_ptr, int *grid_points );
char *
stretching_name( generic_stretching *stretch_ptr );
generic_stretching *
read_stretching( int32 dir );
void 
write_stretch( int32 dir, generic_stretching *stretch_ptr );
generic_stretching *
delete_generic_stretching( generic_stretching *stretch_ptr );
generic_stretching *
copy_stretching( generic_stretching *old_stretch_ptr );
int 
set_stretching(input_output *io_ptr, generic_stretching *stretch_ptr,
	       int *grid_points );
void 
uniform_to_stretch(real uniform, generic_stretching *stretch_ptr, 
		   real *stretch, real *d_stretch, real *d2_stretch );
void 
stretch_to_uniform(real stretch, generic_stretching *stretch_ptr, 
		   real *uniform, real *d_uniform, real *d2_uniform );

#endif
