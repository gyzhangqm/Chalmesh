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

static void 
gridmaterials(void)
{
  GLfloat back_mat_diffuse[] = { 0.64, 0.36, 0.0, 1.0 };
  GLfloat back_mat_ambient[] = { 0.64, 0.36, 0.0, 1.0 };

/*   glMaterialfv(GL_FRONT, GL_DIFFUSE, front_mat_diffuse); */
/*   glMaterialfv(GL_FRONT, GL_AMBIENT, front_mat_ambient); */
   glMaterialfv(GL_BACK, GL_DIFFUSE, back_mat_diffuse); 
   glMaterialfv(GL_BACK, GL_AMBIENT, back_mat_ambient); 

}

static void 
surfacematerials(void)
{
  float front_mat_diffuse[] = { 0.2, 0.7, 0.4, 1.0 };
  float front_mat_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  float front_mat_specular[] = { 0.03, 0.03, 0.03, 1.0 };
  float back_mat_diffuse[] = { 1.0, 1.0, 0.2, 1.0 };
  float back_mat_ambient[] = { 0.1, 0.1, 0.1, 1.0 };

  glMaterialfv(GL_FRONT, GL_DIFFUSE, front_mat_diffuse);
  glMaterialfv(GL_FRONT, GL_AMBIENT, front_mat_ambient);
  glMaterialfv(GL_BACK, GL_DIFFUSE, back_mat_diffuse);
  glMaterialfv(GL_BACK, GL_AMBIENT, back_mat_ambient);

  glMaterialfv(GL_FRONT, GL_SPECULAR, front_mat_specular);
  glMateriali(GL_FRONT, GL_SHININESS, 10);
}

void
draw_surface( surface_mapping *surface, int plot_mode ){
  int i, j;
  real hu, hv, length, r_length, s_length;
  real n_vec[3], start[3], forward[3];
  grid_point gp;
  GLfloat position[] = { 0.0, 0.0, 1.0, 0.0 };
  const real eps=1.e-10;

/*  plot_mode & 1: plot the grid lines (using the forward_surface_mapping) */
/*  plot_mode & 2: plot the opaque surface (using the forward_surface_mapping) */
/*  plot_mode & 4: draw the boundaries with thick lines */
/*  plot_mode & 8: draw the parameter arrows with labels */
/*  plot_mode & 16: draw the boundary condition labels */
/*  plot_mode & 32: draw the curve labels */
/*  plot_mode & 256: plot the coordinate cage */

/* use a light which is fixed to the viewer */
  glPushMatrix();
/* position the light */
    glLoadIdentity();
    glLightfv(GL_LIGHT0, GL_POSITION, position);
  glPopMatrix();

/* Default is to use one-sided lighting */
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

/* test: cull backfacing polygons */
/*   glCullFace( GL_BACK ); */
/*   glEnable( GL_CULL_FACE ); */

/* grid sizes */
  hu = 1.0/((real) surface->r_points - 1); 
  hv = 1.0/((real) surface->s_points - 1); 

  if (plot_mode & 1 && !(plot_mode & 2)) {
/* Use two-sided lighting */
/*    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);*/
/* simple way of implementing hidden lines! Unfortunately, it doesn't work! */
   glCullFace( GL_BACK ); 
   glEnable( GL_CULL_FACE ); 
/* set surface props */
   ogl_set_color( surface->color);

/* set the orientation */
   glFrontFace(GL_CW);

    for (i=0; i<surface->r_points; i++){
/* start a new strip */
      glBegin(GL_LINE_STRIP);
      for (j=0; j<surface->s_points; j++){
	gp.r = i*hu;
	gp.s = j*hv;

	forward_surface_mapping( &gp, surface );
	compute_surface_normal(n_vec, &gp);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

      } /* end for j... */
/* end the strip */
      glEnd();
    }

    for (j=0; j<surface->s_points; j++){
/* start a new strip */
      glBegin(GL_LINE_STRIP);
      for (i=0; i<surface->r_points; i++){
	gp.r = i*hu;
	gp.s = j*hv;

	forward_surface_mapping( &gp, surface );
	compute_surface_normal(n_vec, &gp);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

      } /* end for j... */
/* end the strip */
      glEnd();
    }

    glDisable( GL_CULL_FACE ); 
  } /* end (plot_mode & 1) */


/* draw the surface with or without the grid */
  if (plot_mode & 2){
/* set the surface properties */
    if (plot_mode & 1) { /* also plot the grid lines */
/* Use one-sided lighting */
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
/* set the material props */
      ogl_set_color( surface->color );
      gridmaterials();
/* Actually, both the front and the back faces do matter */
      glPolygonMode(GL_FRONT, GL_LINE);
      glPolygonMode(GL_BACK, GL_FILL);
    }
    else{
/* Use two-sided lighting */
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
/* set surface props */
      surfacematerials();
/* both sides matter because of the two-sided lighting */
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

/* set the orientation */
    glFrontFace(GL_CW);

    for (i=0; i<surface->r_points-1; i++){
      /* start a new strip */
      glBegin(GL_QUAD_STRIP);
      for (j=0; j<surface->s_points; j++){
/* first point */
	gp.r = (i+1)*hu;
	gp.s = j*hv;

	forward_surface_mapping( &gp, surface );
	compute_surface_normal(n_vec, &gp);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	gp.r = i*hu;
	gp.s = j*hv;

	forward_surface_mapping( &gp, surface );
	compute_surface_normal(n_vec, &gp);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

      } /* end for j... */
/* end the strip */
      glEnd();
    }
  }

/* draw the boundary with thick lines */
  if (plot_mode & 4){
/* set the orientation */
    glFrontFace(GL_CW);

/* set the color */
    ogl_set_color( surface->color );
/* set the line width */
    glLineWidth(3.0);
    glEnable(GL_LINE_SMOOTH);
    glBegin(GL_LINE_STRIP);

/* grid sizes */
    hu = 1.0/((real) surface->r_points - 1); 
    hv = 1.0/((real) surface->s_points - 1); 

    j=0;
    for (i=0; i<surface->r_points; i++){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_surface_mapping( &gp, surface );
      compute_surface_normal(n_vec, &gp);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for i... */

    i=surface->r_points-1;
    for (j=0; j<surface->s_points; j++){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_surface_mapping( &gp, surface );
      compute_surface_normal(n_vec, &gp);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for j... */

    j=surface->s_points-1;
    for (i = surface->r_points-1; i>=0; i--){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_surface_mapping( &gp, surface );
      compute_surface_normal(n_vec, &gp);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for i... */

    i=0;
    for (j=surface->s_points-1; j>=0; j--){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_surface_mapping( &gp, surface );
      compute_surface_normal(n_vec, &gp);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for j... */

    glEnd();
/* reset the line width */
    glLineWidth(1.0);
    glDisable(GL_LINE_SMOOTH);
  }

/* draw the parameter directions */
  if (plot_mode & 8){
/* set orientation */
    glFrontFace( GL_CCW );

/* set the color */
    ogl_set_color( surface->color );

/* arrow length */
    length = 0.1 * real_max(surface->bb->x_max - surface->bb->x_min, 
			    real_max(surface->bb->y_max - surface->bb->y_min, 
				     surface->bb->z_max - surface->bb->z_min));
/* r-arrow */
    gp.r = 1.0;
    gp.s = 0.0;
    forward_surface_mapping( &gp, surface );
    compute_surface_normal( n_vec, &gp );
    start[0] = gp.x;
    start[1] = gp.y;
    start[2] = gp.z;
    r_length = sqrt( gp.xr*gp.xr + gp.yr*gp.yr + gp.zr*gp.zr );
    if (r_length < eps){
      r_length = 1.0;
    }
    forward[0] = gp.xr/r_length;
    forward[1] = gp.yr/r_length;
    forward[2] = gp.zr/r_length;

    ogl_draw_arrow(length, 3.0, start, forward, n_vec, "r");

/* s-arrow */
    gp.s = 1.0;
    gp.r = 0.0;
    forward_surface_mapping( &gp, surface );
    compute_surface_normal( n_vec, &gp );
    start[0] = gp.x;
    start[1] = gp.y;
    start[2] = gp.z;
    s_length = sqrt( gp.xs*gp.xs + gp.ys*gp.ys + gp.zs*gp.zs );
    if (s_length < eps){
      s_length = 1.0;
    }
    forward[0] = gp.xs/s_length;
    forward[1] = gp.ys/s_length;
    forward[2] = gp.zs/s_length;

    ogl_draw_arrow(length, 3.0, start, forward, n_vec, "s");
  } /* end plot_mode & 8 */
}

int
surface_plot_mode(input_output *io_ptr, surface_mapping_list *surface_list,
		  surface_mapping *surface, int plot_mode){
  int quit=0, replot=1, icom;
  surface_mapping_link *surf_link;
  bounding_box global_bb = {1000.0, -1000.0, 1000.0, -1000.0, 1000.0, -1000.0};
  surface_mapping *one_surf;

#include "surface_plot_mode_com.h"

/* compute the bounding box if we have a list */
  if (surface_list != NULL)
    {
      for (surf_link = surface_list->first; surf_link != NULL; 
	   surf_link = surf_link->next)
	{
	  one_surf = surf_link->data;
	  global_bb.x_min = real_min(global_bb.x_min, one_surf->bb->x_min); 
	  global_bb.x_max = real_max(global_bb.x_max, one_surf->bb->x_min); 
	  global_bb.y_min = real_min(global_bb.y_min, one_surf->bb->y_min); 
	  global_bb.y_max = real_max(global_bb.y_max, one_surf->bb->y_max); 
	  global_bb.z_min = real_min(global_bb.z_min, one_surf->bb->z_min); 
	  global_bb.z_max = real_max(global_bb.z_max, one_surf->bb->z_max); 
	}
    }

  do{

    if (replot){
/* redraw the surface grids */
      if (surface_list != NULL){
	if (ogl_start_plot(OGL_NEW_PLOT, SURFACE_W, 
			   ogl_length_scale(&global_bb))){
	  for (surf_link = surface_list->first; surf_link != NULL; 
	       surf_link = surf_link->next)
	    draw_surface( surf_link->data, plot_mode);
	  if (plot_mode & 256) ogl_bb_cage( &global_bb, OGL_WHITE );
	  ogl_end_plot();
	  replot = 0;
	}
      }
      else{
	if (ogl_start_plot(OGL_NEW_PLOT, SURFACE_W, 
			   ogl_length_scale(surface->bb))){
	  draw_surface( surface, plot_mode );
	  if (plot_mode & 256) ogl_bb_cage( surface->bb, OGL_WHITE );
	  ogl_end_plot();
	  replot = 0;
	}
      }
    }

    switch( get_command( io_ptr, "set plot mode>", COMMAND, LAST_COM, LEVEL, 
			SAVE_ON_COPY, ARGUMENT) ){

    case HELP:
      while ((icom = get_command( io_ptr, "help on subject (? lists all subjects)>", 
				 COMMAND, LAST_COM, LEVEL, SAVE_ON_COPY, ARGUMENT)) 
	     == -1);
      if (BRIEF_HELP[icom] == NULL)
	general_help();
      else
	print_help( BRIEF_HELP[icom] );
      break;

    case EXIT:
      quit = 1;
      break;

    case GRID_LINES:
      if (plot_mode & 1)
	plot_mode &= ~1;
      else
	plot_mode |= 1;
      replot = 1;
      break;

    case STD_VIEW:
      ogl_standard_view( SURFACE_W, (surface_list == NULL)? surface->bb : 
			&global_bb );
      replot = 1;
      break;

    case SURFACE:
      if (plot_mode & 2)
	plot_mode &= ~2;
      else
	plot_mode |= 2;
      replot = 1;
      break;

    case BOUNDARIES:
      if (plot_mode & 4)
	plot_mode &= ~4;
      else
	plot_mode |= 4;
      replot = 1;
      break;

    case ARROWS:
      if (plot_mode & 8)
	plot_mode &= ~8;
      else
	plot_mode |= 8;
      replot = 1;
      break;

    case COORD_CAGE:
      if (plot_mode & 256)
	plot_mode &= ~256;
      else
	plot_mode |= 256;
      replot = 1;
      break;
      
    default:
      break;
    }
  } while (!quit);

/* pass back the new plot_mode */
  return plot_mode;

/* done */
}
