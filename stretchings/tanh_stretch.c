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
#include "stretch_internal.h"

/* prototypes for private functions */
static void 
tanh_stretch_to_uniform(real stretch, generic_stretching *stretch_ptr, 
			real *uniform, real *d_uniform, real *d2_uniform );
static void 
tanh_uniform_to_stretch(real uniform, generic_stretching *stretch_ptr, 
			real *stretch, real *d_stretch, real *d2_stretch );
static void 
cleanup_tanh_stretch( generic_stretching *stretch_ptr );
static void 
write_tanh_stretch( int32 dir, generic_stretching *stretch_ptr );
static void 
set_tanh_stretch(input_output *io_ptr, generic_stretching *stretch_ptr, 
		 int *grid_points );
static void *
copy_tanh_stretch( void *stretch_data_ptr );
static void 
compute_tanh_stretch( tanh_stretch_info *info );
/* end private prototypes */


void 
init_tanh_stretch( generic_stretching *stretch_ptr, int grid_points ){
  tanh_stretch_info *info;

/* allocate memory for the arcl_stretch_info structure */
  info = (tanh_stretch_info *) malloc( sizeof(tanh_stretch_info) );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 4;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = tanh_uniform_to_stretch;
  stretch_ptr->stretch_to_uniform = tanh_stretch_to_uniform;
  stretch_ptr->write_stretch_data = write_tanh_stretch;
  stretch_ptr->set_stretching = set_tanh_stretch;
  stretch_ptr->cleanup_stretching = cleanup_tanh_stretch;
  stretch_ptr->copy_stretching = copy_tanh_stretch;

/* default stretching is very mild */
  info->ds_dx_0 = 0.5;
  info->ds_dx_1 = 0.5;

/* compute the coefficients */
  compute_tanh_stretch( info );
}

void 
read_tanh_stretch( int32 dir, generic_stretching *stretch_ptr ){
  tanh_stretch_info *info;

  info = (tanh_stretch_info *) malloc( sizeof(tanh_stretch_info) );

  hget_real( &(info->a),       "a",       dir );
  hget_real( &(info->delta),   "delta",   dir );
  hget_real( &(info->ds_dx_0), "ds_dx_0", dir );
  hget_real( &(info->ds_dx_1), "ds_dx_1", dir );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 4;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = tanh_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_tanh_stretch;
  stretch_ptr->set_stretching = set_tanh_stretch;
  stretch_ptr->cleanup_stretching = cleanup_tanh_stretch;
  stretch_ptr->copy_stretching = copy_tanh_stretch;

}

static void
set_tanh_stretch(input_output *io_ptr, generic_stretching *stretch_ptr, 
		 int *grid_points ){
  char *prompt;
  const real eps=1.e-7;
  int icom, quit, replot;
  const int level=3, no_indent=0;
  real sigma0, sigma1;
  tanh_stretch_info *info;

#include "set_tanh_com.h"

  info = (tanh_stretch_info *) stretch_ptr->stretch_data_ptr;

  prompt = "change tanh stretching>";
  replot = 1;

  do{

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    quit = 0;
    replot = 1;
    switch (icom) {

    case 0:
/* call the curve function to get the ratio between its parameter and the arclength */
      sigma0 = 1.0;

/* start grid size */
      info->ds_dx_0 = (*grid_points - 1) / sigma0 * 
	real_max( eps, get_real( io_ptr, "Starting grid size > 0:",
				info->ds_dx_0 * sigma0 / (*grid_points - 1), 
				no_indent ) );
      compute_tanh_stretch( info );
      break;

    case 1:
/* call the curve function to get the ratio between its parameter and the arclength */
      sigma1 = 1.0;

/* ending grid size */
      info->ds_dx_1 = (*grid_points - 1) / sigma1 * 
	real_max( eps, get_real( io_ptr, "Ending grid size > 0:",
				info->ds_dx_1 * sigma1 / (*grid_points - 1), 
				no_indent ) );
      compute_tanh_stretch( info );
      break;

    case 2:
/* default strength */
      info->ds_dx_0 = info->ds_dx_1 = 0.5;
      compute_tanh_stretch( info );
      break;

    case 3:
/* change number of grid points */
      *grid_points = int_max(2, get_int(io_ptr, 
					"Enter number of gridlines >= 2: ", 
					*grid_points,
					no_indent));
      break;

    case 4:
      sigma0 = sigma1 = 1.0;

/* show stretching parameters */
      printf("Stretching parameters:\n");
      printf("Number of grid points: %i\n", *grid_points );
      printf("Starting grid size: %f\n", info->ds_dx_0 * sigma0 / (*grid_points - 1));
      printf("Ending grid size: %f\n", info->ds_dx_1 * sigma1 / (*grid_points - 1));
      replot = 0;
      break;

    case 5:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
			 command, ncom, level+1, NULL, NULL)) == -1);
      if (brief_help[icom] == NULL)
	general_help();
      else
	print_help( brief_help[icom] );
      replot = 0;
      break;

    case 6:
/* exit */
      quit = 1;
      break;

    default:
      ;

    }
  }
  while( !quit );
}


static void 
write_tanh_stretch( int32 dir, generic_stretching *stretch_ptr ){
  tanh_stretch_info *info;

  if (stretch_ptr == NULL) return;

  info = (tanh_stretch_info *) stretch_ptr->stretch_data_ptr;

  hput_real( info->a,       "a",       dir );
  hput_real( info->delta,   "delta",   dir );
  hput_real( info->ds_dx_0, "ds_ds_0", dir );
  hput_real( info->ds_dx_1, "ds_dx_1", dir );
}


static void 
tanh_uniform_to_stretch( real uniform, generic_stretching *stretch_ptr, 
			real *stretch, real *d_stretch, real *d2_stretch ){
  tanh_stretch_info *info;
  real u, d_u, d2_u;

  info = (tanh_stretch_info *) stretch_ptr->stretch_data_ptr;

/* help function u */
  u    = 0.5 * ( 1.0 + tanh( info->delta*( uniform - 0.5 ) ) /
		tanh( 0.5*info->delta ) );
  d_u  = 0.5 * info->delta / tanh( 0.5*info->delta ) /
    cosh( info->delta*( uniform - 0.5 ) ) / 
      cosh( info->delta*( uniform - 0.5 ) );
  d2_u = - info->delta * info->delta * tanh( info->delta*( uniform - 0.5 ) )/ 
    tanh( 0.5*info->delta ) /
      cosh( info->delta*( uniform - 0.5 ) ) / 
	cosh( info->delta*( uniform - 0.5 ) );

/* hyperbolic tangent stretching function */
  *stretch    = u / ( info->a + (1.0 - info->a) * u );
  *d_stretch  = info->a * d_u / (info->a + (1.0-info->a)*u) / 
    (info->a + (1.0-info->a)*u);
  *d2_stretch = info->a * d2_u / (info->a + (1.0-info->a)*u) / 
    (info->a + (1.0-info->a)*u) - 
      2.0 * info->a * (1.0 - info->a) * d_u * d_u /
	(info->a + (1.0-info->a)*u) /
	  (info->a + (1.0-info->a)*u) /
	    (info->a + (1.0-info->a)*u);
}

static void 
tanh_stretch_to_uniform(real stretch, generic_stretching *stretch_ptr, 
			real *uniform, real *d_uniform, real *d2_uniform ){
  tanh_stretch_info *info;
  real u, d_u, d2_u, d_stretch, d2_stretch;

  info = (tanh_stretch_info *) stretch_ptr->stretch_data_ptr;

/* inverse of the tanh stretching */
  u = stretch * info->a / (1.0 - stretch*(1.0 - info->a));
  *uniform = 0.5 + atanh( (2.0*u - 1.0) * tanh( 0.5*info->delta ) ) / info->delta;

/* help function u */
  d_u  = 0.5 * info->delta / tanh( 0.5*info->delta ) /
    cosh( info->delta*( (*uniform) - 0.5 ) ) / 
      cosh( info->delta*( (*uniform) - 0.5 ) );
  d2_u = - info->delta * info->delta * tanh( info->delta*( (*uniform) - 0.5 ) )/ 
    tanh( 0.5*info->delta ) /
      cosh( info->delta*( (*uniform) - 0.5 ) ) / 
	cosh( info->delta*( (*uniform) - 0.5 ) );

/* hyperbolic tangent stretching function */
  d_stretch  = info->a * d_u / (info->a + (1.0-info->a)*u) / 
    (info->a + (1.0-info->a)*u);
  d2_stretch = info->a * d2_u / (info->a + (1.0-info->a)*u) / 
    (info->a + (1.0-info->a)*u) - 
      2.0 * info->a * (1.0 - info->a) * d_u * d_u /
	(info->a + (1.0-info->a)*u) /
	  (info->a + (1.0-info->a)*u) /
	    (info->a + (1.0-info->a)*u);

/* compute the derivatives of the inverse stretching function */
  *d_uniform = 1.0/d_stretch;
  *d2_uniform = - d2_stretch/(d_stretch*d_stretch*d_stretch);
}

static void 
compute_tanh_stretch( tanh_stretch_info *info ){
  real b, d, res;
  const real b_min = 1.0 + 1.e-4;
  int i;
  const int n_iter=10;

  info->a = sqrt( info->ds_dx_1 / info->ds_dx_0 );

  b = 1.0 / sqrt( info->ds_dx_1 * info->ds_dx_0 );

/* solve sinh( d ) = d * b */
/* this equation is only solvable if b >= 1 */

  if (b < b_min){
    printf("Error in compute_tanh_stretch: b < %f\n", b_min);
    b = b_min;
    info->ds_dx_1 = 1.0 / info->ds_dx_0 / b / b;
  }

/* use Taylor series expansion to approximate sinh d = d + d^3/6 + ... */
/* in order to get an initial guess */    
  d = sqrt( 6.0*( b - 1.0 ) );

/* newton iteration */
  for (i=1; i<=n_iter; i++){
    d = d - (sinh(d) - b*d)/(cosh(d) - b);
  }

/* check for reasonable convergence */
  if ((res=fabs( sinh(d) - b*d )) > NEWTON_EPS)
    printf("Warning: Newton iteration in compute_tanh_stretch converged "
	   "poorly.\nResidual after %i iterations: %e\n", n_iter, res);

  info->delta = d;
}

static void cleanup_tanh_stretch( generic_stretching *stretch_ptr ){
  tanh_stretch_info *info;

  if (stretch_ptr == NULL) return;

/* cast the pointer to the right type */
  info = (tanh_stretch_info *) stretch_ptr->stretch_data_ptr;

  free( info );
}

static void *copy_tanh_stretch( void *stretch_data_ptr ){
  tanh_stretch_info *info, *old_info;
  
  old_info = (tanh_stretch_info *)stretch_data_ptr;

  info = (tanh_stretch_info *) malloc( sizeof(tanh_stretch_info) );

  info->a = old_info->a;
  info->delta = old_info->delta;
  info->ds_dx_0 = old_info->ds_dx_0;
  info->ds_dx_1 = old_info->ds_dx_1;

  return (void *)info;
}
