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
exp_stretch_to_uniform(real stretch, generic_stretching *stretch_ptr, 
		       real *uniform, real *d_uniform, real *d2_uniform );
static void 
exp_uniform_to_stretch(real uniform, generic_stretching *stretch_ptr, real *stretch, 
		       real *d_stretch, real *d2_stretch );
static void 
cleanup_exp_stretch( generic_stretching *stretch_ptr );
static void 
write_exp_stretch( int32 dir, generic_stretching *stretch_ptr );
static void
set_exp_stretch(input_output *io_ptr, generic_stretching *stretch_ptr, 
		int *grid_points );
static void *
copy_exp_stretch( void *stretch_data_ptr );
static void 
compute_exp_stretch( exp_stretch_info *info );
/* end private prototypes */


void 
init_exp_stretch( generic_stretching *stretch_ptr, int grid_points ){
  exp_stretch_info *info;

/* allocate memory for the arcl_stretch_info structure */
  info = (exp_stretch_info *) malloc( sizeof(exp_stretch_info) );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 2;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = exp_uniform_to_stretch;
  stretch_ptr->stretch_to_uniform = exp_stretch_to_uniform;
  stretch_ptr->write_stretch_data = write_exp_stretch;
  stretch_ptr->set_stretching = set_exp_stretch;
  stretch_ptr->cleanup_stretching = cleanup_exp_stretch;
  stretch_ptr->copy_stretching = copy_exp_stretch;

/* focus at r=0 */
  info->focus = 0;

/* give a default starting grid step */
  info->ds_dx_0 = 0.5;

/* compute alpha */
  compute_exp_stretch( info );
}

void 
read_exp_stretch( int32 dir, generic_stretching *stretch_ptr ){
  exp_stretch_info *info;

  info = (exp_stretch_info *) malloc( sizeof(exp_stretch_info) );

  hget_int( &(info->focus),    "focus",   dir );
  hget_real( &(info->alpha),   "alpha",   dir );
  hget_real( &(info->ds_dx_0), "ds_dx_0", dir );

/* assign all fields of the generic stretching structure */
  stretch_ptr->stretch_type = 2;
  stretch_ptr->stretch_data_ptr = (void *) info;
  stretch_ptr->uniform_to_stretch = exp_uniform_to_stretch;
  stretch_ptr->write_stretch_data = write_exp_stretch;
  stretch_ptr->set_stretching = set_exp_stretch;
  stretch_ptr->cleanup_stretching = cleanup_exp_stretch;
  stretch_ptr->copy_stretching = copy_exp_stretch;

}

static void
set_exp_stretch(input_output *io_ptr, generic_stretching *stretch_ptr, 
		int *grid_points ){
  char *prompt;
  const real eps=1.e-7;
  int icom, quit, replot;
  const int level=3, no_indent=0;
  real sigma;
  exp_stretch_info *info;

#include "set_exp_com.h"

  info = (exp_stretch_info *) stretch_ptr->stretch_data_ptr;

  prompt = "change exp stretching>";

  do{

    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    quit = 0;
    replot = 1;
    switch (icom) {

    case 0:
/* call the curve function to get the ratio between its parameter and the arclength */
      sigma = 1.0;
/* change strength */
      info->ds_dx_0 = (*grid_points - 1) / sigma * 
	real_max( eps, get_real( io_ptr, "Starting grid size > 0:",
				info->ds_dx_0 * sigma/ (*grid_points - 1), 
				no_indent ) );
      compute_exp_stretch( info );
      break;

    case 1:
/* reverse parametrization */
      info->focus = !(info->focus);
      break;

    case 2:
/* change number of grid points */
      *grid_points = int_max(2, get_int(io_ptr, 
					"Enter number of gridlines >= 2: ", 
					*grid_points,
					no_indent));
      break;

    case 3:
/* default strength */
      info->ds_dx_0 = 0.5;
      compute_exp_stretch( info );
      break;

    case 4:
/* call the curve function to get the ratio between its parameter and the arclength */
      sigma = 1.0;
	
/* show stretching parameters */
      printf("Stretching parameters:\n");
      printf("The stretching is focused at %i\n", info->focus );
      printf("Number of grid points: %i\n", *grid_points );
      printf("Starting grid size: %f\n", info->ds_dx_0 * sigma / (*grid_points-1) );
      printf("Relative increase in grid size: %f.\n", 
	     info->alpha / (*grid_points - 1));
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
write_exp_stretch( int32 dir, generic_stretching *stretch_ptr ){
  exp_stretch_info *info;

  if (stretch_ptr == NULL) return;

  info = (exp_stretch_info *) stretch_ptr->stretch_data_ptr;

  hput_int( info->focus, "focus", dir );
  hput_real( info->alpha, "alpha", dir );
  hput_real( info->ds_dx_0, "ds_dx_0", dir );
}


static void 
exp_uniform_to_stretch(real uniform, generic_stretching *stretch_ptr, 
		       real *stretch, real *d_stretch, real *d2_stretch ){
  exp_stretch_info *info;

  info = (exp_stretch_info *) stretch_ptr->stretch_data_ptr;

/* flip the parameter? */
  if (info->focus == 1)
    uniform = 1.0 - uniform;

/* exponential stretching function */
  if (info->alpha > 0.0){
    *stretch = (exp( info->alpha * uniform ) - 1.0)/
      (exp( info->alpha ) - 1.0);
    *d_stretch = info->alpha * exp( info->alpha * uniform )/
      (exp( info->alpha ) - 1.0);
    *d2_stretch = info->alpha * (*d_stretch);
  }
  else{
    *stretch = uniform;
    *d_stretch = 1.0;
    *d2_stretch = 0.0;
  }

  if (info->focus == 1){
    *stretch = 1.0 - *stretch;
    *d2_stretch = - *d2_stretch;
  }
}

static void 
exp_stretch_to_uniform(real stretch, generic_stretching *stretch_ptr, 
		       real *uniform, real *d_uniform, real *d2_uniform ){
  exp_stretch_info *info;
  real u, d_stretch, d2_stretch;

  info = (exp_stretch_info *) stretch_ptr->stretch_data_ptr;

/* flip the parameter? */
  if (info->focus == 1)
    stretch = 1.0 - stretch;

/* inverse exponential stretching */
  if (info->alpha > 0.0){
    u = log( 1.0 + stretch*(exp(info->alpha) - 1.0) ) / info->alpha;

/* forward exponential stretching function */
    d_stretch  = info->alpha * exp( info->alpha * u )/ 
      (exp( info->alpha ) - 1.0);
    d2_stretch = info->alpha * d_stretch;
  }
/* no stretching function */
  else{
    u = stretch;
    d_stretch = 1.0;
    d2_stretch = 0.0;
  }

  if (info->focus == 1)
    d2_stretch = -d2_stretch;

/* compute the derivatives of the inverse stretching function */
  if (info->focus == 1)
    *uniform = 1.0 - u;
  else
    *uniform = u;
  *d_uniform = 1.0/d_stretch;
  *d2_uniform = - d2_stretch/(d_stretch*d_stretch*d_stretch);

}

static void 
cleanup_exp_stretch( generic_stretching *stretch_ptr ){
  exp_stretch_info *info;

  if (stretch_ptr == NULL) return;

/* cast the pointer to the right type */
  info = (exp_stretch_info *) stretch_ptr->stretch_data_ptr;

  free( info );
}

static void *copy_exp_stretch( void *stretch_data_ptr ){
  exp_stretch_info *info, *old_info;
  
  old_info = (exp_stretch_info *)stretch_data_ptr;

  info = (exp_stretch_info *) malloc( sizeof(exp_stretch_info) );

  info->focus = old_info->focus;
  info->alpha = old_info->alpha;

  return (void *)info;
}

static void compute_exp_stretch( exp_stretch_info *info ){
  real b, a, res;
  const real b_min = 1.0 + 1.e-4;
  int i;
  const int n_iter=10;

  b = 1.0 / info->ds_dx_0;

/* solve a: exp(a) - 1 = a * b */
/* this equation is only solvable if b > 1 */

  if (b < b_min){
    printf("Error in compute_exp_stretch: b < %f\n", b_min);
    b = b_min;
    info->ds_dx_0 = 1.0 / b;
  }

/* use Taylor series expansion to approximate exp(a)-1 = a + a^2/2 + a^3/6 + ... */
/* in order to get an initial guess */    
  a = -1.5 + sqrt( 1.5*1.5 + 6.0*( b - 1.0 ) );

/* newton iteration */
  for (i=1; i<=n_iter; i++){
    a = a - (exp(a) - 1.0 - b*a)/(exp(a) - b);
  }

/* check for reasonable convergence */
  if ((res=fabs( exp(a) - 1.0 - b*a )) > NEWTON_EPS)
    printf("Warning: Newton iteration in compute_exp_stretch converged "
	   "poorly.\nResidual after %i iterations: %e\n", n_iter, res);

  info->alpha = a;
}

