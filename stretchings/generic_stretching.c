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

generic_stretching *
choose_stretching(input_output *io_ptr, int *grid_points ){
  char prompt[80];
  int icom, quit;
  generic_stretching *stretch_ptr=NULL;

#include "choose_stretching_com.h"

  sprintf(prompt,"choose a stretching>");

  do{

    icom = get_command(io_ptr, prompt, COMMAND, LAST_COM, LEVEL,
		       SAVE_ON_COPY, ARGUMENT);

    quit = 0;
    switch (icom) {

    case EXPONENTIAL:
/* exponential stretching */
      stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );
      init_exp_stretch( stretch_ptr, *grid_points );
      quit = 1;
      break;
      
    case LAYER:
/* layer stretching */
      stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );
      init_layer_stretch( stretch_ptr );
      quit = 1;
      break;
      
    case TANH:
/* tanh stretching */ 
      stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );
      init_tanh_stretch( stretch_ptr, *grid_points );
      quit = 1;
      break;

    case HELP:
/* help */
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL+1, NULL, NULL)) == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );

      break;

    case CANCEL:
/* cancel */
      return NULL;
      break;

    default:
      ;
    }
  }
  while(!quit);

/* set the parameter value */
  set_stretching( io_ptr, stretch_ptr, grid_points );

  return stretch_ptr;
}


char *
stretching_name( generic_stretching *stretch_ptr ){
  char *name;

  if (stretch_ptr == NULL){

    name = "no stretching";

  }
  else{
    
    switch( stretch_ptr->stretch_type ){

    case 1:
      name = "arclength stretching";
      break;

    case 2:
      name = "exponential stretching";
      break;

    case 3:
      name = "layer stretching";
      break;

    default:
      name = "unknown stretching";
    }

  }

  return name;
}



generic_stretching *
read_stretching( int32 dir ){
  generic_stretching *stretch_ptr=NULL;
  int stretch_type;

  hget_int( &(stretch_type), "stretch_type", dir );

  if (stretch_type == 0)
    return NULL;
  else{
    stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );
    stretch_ptr->stretch_type = stretch_type;
  }

/* assign the appropriate curve-specific pointers */
  switch( stretch_ptr->stretch_type ){

  case 1:
/* the arclength stretching does not make any sense for surfaces */

    break;

  case 2:
/* exponential */
    read_exp_stretch( dir, stretch_ptr );
    break;

  case 3:
/* layer stretching */
    read_layer_stretch( dir, stretch_ptr );
    break;

  case 4:
/* hyperbolic tangent */
    read_tanh_stretch( dir, stretch_ptr );
    break;

  default:
    printf("Unknown stretch type: %i in read stretch\n", 
	   stretch_ptr->stretch_type);
    exit(1);
  }

  return stretch_ptr;
}

void 
write_stretch( int32 dir, generic_stretching *stretch_ptr ){

/* save the stretch type */
  hput_int( (stretch_ptr? stretch_ptr->stretch_type: 0), "stretch_type", dir );

/* only save stretching parameters if there is a stretching! */
  if (!stretch_ptr) return;

/* write the specific data (in the same directory to make it simple) */
  (*stretch_ptr->write_stretch_data)( dir, stretch_ptr );

}


generic_stretching *
delete_generic_stretching( generic_stretching *stretch_ptr ){
  if (stretch_ptr != NULL){
    (*stretch_ptr->cleanup_stretching)( stretch_ptr );
    free( stretch_ptr );
  }
  return NULL;
}

void 
uniform_to_stretch( real uniform, generic_stretching *stretch_ptr, 
			real *stretch, real *d_stretch, real *d2_stretch ){
/* *d_stretch will contain the derivataive of the stretching function with respect */
/* to the uniform parameter, and *d2_strech will contain the second derivative. */
  if ( stretch_ptr == NULL ){
    *stretch = uniform;
    *d_stretch = 1.0;
    *d2_stretch = 0.0;
  }
  else
    (*stretch_ptr->uniform_to_stretch)( uniform, stretch_ptr, stretch, 
				       d_stretch, d2_stretch );
}


void 
stretch_to_uniform(real stretch, generic_stretching *stretch_ptr, 
		   real *uniform, real *d_uniform, real *d2_uniform ){
/* *d_uniform will contain the derivataive of the inverse stretching function */
/* with respect to the stretch parameter. d2_uniform is the second derivative. */
  if ( stretch_ptr == NULL ){
    *uniform = stretch;
    *d_uniform = 1.0;
  }
  else
    (*stretch_ptr->stretch_to_uniform)(stretch, stretch_ptr,  uniform, 
				       d_uniform, d2_uniform );
}


generic_stretching *
copy_stretching( generic_stretching *old_stretch_ptr ){
  generic_stretching *stretch_ptr;

/* it is legal to make a copy of "no stretching" */
  if (old_stretch_ptr == NULL) return NULL;

  stretch_ptr = (generic_stretching *) malloc( sizeof(generic_stretching) );

  stretch_ptr->stretch_type       = old_stretch_ptr->stretch_type;
  stretch_ptr->uniform_to_stretch = old_stretch_ptr->uniform_to_stretch;
  stretch_ptr->write_stretch_data = old_stretch_ptr->write_stretch_data;
  stretch_ptr->set_stretching     = old_stretch_ptr->set_stretching;
  stretch_ptr->cleanup_stretching = old_stretch_ptr->cleanup_stretching;
  stretch_ptr->copy_stretching    = old_stretch_ptr->copy_stretching;

/* copy the specific data */
  stretch_ptr->stretch_data_ptr = 
    (*old_stretch_ptr->copy_stretching)( old_stretch_ptr->stretch_data_ptr );

  return stretch_ptr;
}


int 
set_stretching(input_output *io_ptr, generic_stretching *stretch_ptr,
	       int *grid_points ){
  if (stretch_ptr == NULL || stretch_ptr->set_stretching == NULL)
    return 0;
  else{
    (*stretch_ptr->set_stretching)( io_ptr, stretch_ptr, grid_points );
    return 1;
  }
}
