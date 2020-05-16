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
#include "hole_surface_internal.h"

void
draw_hole_surface( hole_surface *overg, int plot_mode ){
  int i, j, ip, jp;
  real hu, hv, length, r_length, s_length, start[3], forward[3];
  real n_vec[3];
  const real eps=1.e-7;
  GLfloat position[] = { 0.0, 0.0, 1.0, 0.0 };
  discrete_surface *sgrid_ptr;
  discrete_surface_link *sgrid_link;
  hole_cutting_curve *hole_curve;
  hole_cutting_curve_link *hole_curve_link;

/*  plot_mode & 1: plot all grid lines */
/*  plot_mode & 2: plot parameter arrows */
/*  plot_mode & 4: plot the grid boundaries with thick lines */
/*  plot_mode & 8: color code the flag array */
/*  plot_mode & 32: draw the curve labels */
/*  plot_mode & 64: plot the grid lines with flag != 0 (using the x, y, z arrays) */
/*  plot_mode & 1024: plot all edge curves */

/* use a light which is fixed to the viewer */
  glPushMatrix();
/* position the light */
    glLoadIdentity();
    glLightfv(GL_LIGHT0, GL_POSITION, position);
  glPopMatrix();

/* Default is to use two-sided lighting */
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

/* test: cull backfacing polygons */
/*   glCullFace( GL_BACK ); */
/*   glEnable( GL_CULL_FACE ); */
  
/* loop over all component grids */
  for (sgrid_link = overg->grid_list->first; sgrid_link != NULL; 
       sgrid_link = sgrid_link->next){
    sgrid_ptr = sgrid_link->data;

    if (!sgrid_ptr->plot_it) continue;

/* grid sizes */
    hu = sgrid_ptr->r_step;
    hv = sgrid_ptr->s_step;

/* draw all grid lines */
    if (plot_mode & 1) {

/* set the color */
      ogl_set_color( sgrid_ptr->priority);

/* set the orientation */
      glFrontFace(GL_CW);

      for (i=1; i<=sgrid_ptr->r_dim; i++){
/* start a new strip */
	glBegin(GL_LINE_STRIP);
	for (j=1; j<=sgrid_ptr->s_dim; j++){
	  compute_grid_normal(n_vec, i, j, sgrid_ptr);
	  glNormal3dv(n_vec);
	  glVertex3d(x(i,j), y(i,j), z(i,j));
	} /* end for j... */
/* end the strip */
	glEnd();
      }

      for (j=1; j<=sgrid_ptr->s_dim; j++){
/* start a new strip */
	glBegin(GL_LINE_STRIP);
	for (i=1; i<=sgrid_ptr->r_dim; i++){
	  compute_grid_normal(n_vec, i, j, sgrid_ptr);
	  glNormal3dv(n_vec);
	  glVertex3d(x(i,j), y(i,j), z(i,j));
	} /* end for j... */
/* end the strip */
	glEnd();
      }
    } /* end (plot_mode & 1) */

/* draw the parameter directions */
  if (plot_mode & 2){
/* set orientation */
    glFrontFace( GL_CCW );

/* set the color */
    ogl_set_color( sgrid_ptr->priority );

/* arrow length */
    length = 0.1 * real_max(sgrid_ptr->bb->x_max - sgrid_ptr->bb->x_min, 
			    real_max(sgrid_ptr->bb->y_max - sgrid_ptr->bb->y_min, 
				     sgrid_ptr->bb->z_max - sgrid_ptr->bb->z_min));
/* r-arrow */
    compute_grid_normal(n_vec, sgrid_ptr->r_dim, 1, sgrid_ptr);
    start[0] = x(sgrid_ptr->r_dim,1);
    start[1] = y(sgrid_ptr->r_dim,1);
    start[2] = z(sgrid_ptr->r_dim,1);
    forward[0] = x(sgrid_ptr->r_dim,1) - x(sgrid_ptr->r_dim-1,1);
    forward[1] = y(sgrid_ptr->r_dim,1) - y(sgrid_ptr->r_dim-1,1);
    forward[2] = z(sgrid_ptr->r_dim,1) - z(sgrid_ptr->r_dim-1,1);
#define sqr(x) ((x)*(x))
    r_length = sqrt( sqr(forward[0]) + sqr(forward[1]) + sqr(forward[2]) );
    if (r_length < eps){
      r_length = 1.0;
    }
    forward[0] = forward[0]/r_length;
    forward[1] = forward[1]/r_length;
    forward[2] = forward[2]/r_length;

    ogl_draw_arrow(length, 3.0, start, forward, n_vec, "r");

/* s-arrow */
    compute_grid_normal(n_vec, 1, sgrid_ptr->s_dim, sgrid_ptr);
    start[0] = x(1,sgrid_ptr->s_dim);
    start[1] = y(1,sgrid_ptr->s_dim);
    start[2] = z(1,sgrid_ptr->s_dim);
    forward[0] = x(1,sgrid_ptr->s_dim) - x(1,sgrid_ptr->s_dim-1);
    forward[1] = y(1,sgrid_ptr->s_dim) - y(1,sgrid_ptr->s_dim-1);
    forward[2] = z(1,sgrid_ptr->s_dim) - z(1,sgrid_ptr->s_dim-1);
    s_length = sqrt( sqr(forward[0]) + sqr(forward[1]) + sqr(forward[2]) );
#undef sqr
    if (s_length < eps){
      s_length = 1.0;
    }
    forward[0] = forward[0]/s_length;
    forward[1] = forward[1]/s_length;
    forward[2] = forward[2]/s_length;

    ogl_draw_arrow(length, 3.0, start, forward, n_vec, "s");
  } /* end plot_mode & 2 */

/* draw the boundary with thick lines */
    if (plot_mode & 4){

/* set the color */
      ogl_set_color( sgrid_ptr->priority );

/* set the line width */
      glLineWidth(3.0);
      glEnable(GL_LINE_SMOOTH);
      glBegin(GL_LINE_STRIP);

      j=1;
      for (i=1; i<=sgrid_ptr->r_dim; i++){
	compute_grid_normal(n_vec, i, j, sgrid_ptr);
	glNormal3dv(n_vec);
	glVertex3d(x(i,j), y(i,j), z(i,j));
      } /* end for i... */

      i=sgrid_ptr->r_dim;
      for (j=1; j<=sgrid_ptr->s_dim; j++){
	compute_grid_normal(n_vec, i, j, sgrid_ptr);
	glNormal3dv(n_vec);
	glVertex3d(x(i,j), y(i,j), z(i,j));
      } /* end for j... */

      j=sgrid_ptr->s_dim;
      for (i=sgrid_ptr->r_dim; i>=1; i--){
	compute_grid_normal(n_vec, i, j, sgrid_ptr);
	glNormal3dv(n_vec);
	glVertex3d(x(i,j), y(i,j), z(i,j));
      } /* end for i... */

      i=1;
      for (j=sgrid_ptr->s_dim; j>=1; j--){
	compute_grid_normal(n_vec, i, j, sgrid_ptr);
	glNormal3dv(n_vec);
	glVertex3d(x(i,j), y(i,j), z(i,j));
      } /* end for j... */

      glEnd();
/* reset the line width */
      glLineWidth(1.0);
      glDisable(GL_LINE_SMOOTH);
    }

/* draw the grid lines where flag != 0 by using the (x, y, z) arrays */
/* check that the essential arrays are there */
  if (plot_mode & 8 && sgrid_ptr->x_ptr != NULL && sgrid_ptr->y_ptr != NULL && 
	sgrid_ptr->z_ptr != NULL && sgrid_ptr->flag_ptr != NULL){

/* set the orientation */
    glFrontFace(GL_CW);
    glPointSize(5.0);
/* start a new line segment */
    glBegin(GL_POINTS);
    for (i=1; i<=sgrid_ptr->r_dim; i++)
      for (j=1; j<sgrid_ptr->s_dim; j++)
	{
	  ogl_set_color(abs(flag(i,j)));
	  compute_grid_normal(n_vec, i, j, sgrid_ptr);
	  glNormal3dv(n_vec);
	  glVertex3d(x(i,j), y(i,j), z(i,j));
	} /* end for i,j... */
    glEnd();

    glPointSize(1.0);
  } /* end plot_mode & 8 */

/* curve label */
  if (plot_mode & 32){

/* set the line width */
    glLineWidth(3.0);
    glEnable(GL_LINE_SMOOTH);

/* left */
    if (curve(1,1) > 0){
/* set color */
      ogl_set_color( curve(1,1) );
      i = 1;
      glBegin(GL_LINE_STRIP);
      for (j=1; j<=sgrid_ptr->s_dim; j++){
	compute_grid_normal(n_vec, i, j, sgrid_ptr);
	glNormal3dv(n_vec);
	glVertex3d(x(i,j), y(i,j), z(i,j));
      } /* end for j... */
      glEnd();
    }

/* right */
    if (curve(2,1) > 0){
/* set color */
      ogl_set_color( curve(2,1) );
      i = sgrid_ptr->r_dim;
      glBegin(GL_LINE_STRIP);
      for (j=1; j<=sgrid_ptr->s_dim; j++){
	compute_grid_normal(n_vec, i, j, sgrid_ptr);
	glNormal3dv(n_vec);
	glVertex3d(x(i,j), y(i,j), z(i,j));
      } /* end for j... */
      glEnd();
    }

/* lower */
    if (curve(1,2) > 0){
/* set color */
      ogl_set_color( curve(1,2) );
      j = 1;
      glBegin(GL_LINE_STRIP);
      for (i=1; i<=sgrid_ptr->r_dim; i++){
	compute_grid_normal(n_vec, i, j, sgrid_ptr);
	glNormal3dv(n_vec);
	glVertex3d(x(i,j), y(i,j), z(i,j));
      } /* end for i... */
      glEnd();
    }

/* upper */
    if (curve(2,2) > 0){
/* set color */
      ogl_set_color( curve(2,2) );
      j = sgrid_ptr->s_dim;
      glBegin(GL_LINE_STRIP);
      for (i=1; i<=sgrid_ptr->r_dim; i++){
	compute_grid_normal(n_vec, i, j, sgrid_ptr);
	glNormal3dv(n_vec);
	glVertex3d(x(i,j), y(i,j), z(i,j));
      } /* end for i... */
      glEnd();
    }

/* reset the line width */
    glLineWidth(1.0);
    glDisable(GL_LINE_SMOOTH);
  } /* end plot_mode & 32 */

/* draw the grid lines where flag != 0 by using the (x, y, z) arrays */
/* check that the essential arrays are there */
  if (plot_mode & 64 && sgrid_ptr->x_ptr != NULL && sgrid_ptr->y_ptr != NULL && 
	sgrid_ptr->z_ptr != NULL && sgrid_ptr->flag_ptr != NULL){
/* simple way of implementing hidden lines! Unfortunately, it doesn't work!*/
    glCullFace( GL_BACK ); 
    glEnable( GL_CULL_FACE ); 

/* set surface props */
    ogl_set_color( sgrid_ptr->priority);

/* set the orientation */
    glFrontFace(GL_CW);

    for (i=1; i<=sgrid_ptr->r_dim; i++){
      for (j=1; j<sgrid_ptr->s_dim; j++){
	if (flag(i,j) != 0 && flag(i,j+1) != 0){
/* start a new line segment */
	  glBegin(GL_LINES);
/* first point */
	  compute_grid_normal(n_vec, i, j, sgrid_ptr);
	  glNormal3dv(n_vec);
	  glVertex3d(x(i,j), y(i,j), z(i,j));
/* second point */
	  compute_grid_normal(n_vec, i, j+1, sgrid_ptr);
	  glNormal3dv(n_vec);
	  glVertex3d(x(i,j+1), y(i,j+1), z(i,j+1));
/* end the line segment */
	  glEnd();
	}/* end if flag != 0 */
      } /* end for j... */
    } /* end for i... */

    for (j=1; j<=sgrid_ptr->s_dim; j++){
      for (i=1; i<sgrid_ptr->r_dim; i++){
	if (flag(i,j) != 0 && flag(i+1,j) != 0){
/* start a new line segment */
	  glBegin(GL_LINES);
/* first point */
	  compute_grid_normal(n_vec, i, j, sgrid_ptr);
	  glNormal3dv(n_vec);
	  glVertex3d(x(i,j), y(i,j), z(i,j));
/* second point */
	  compute_grid_normal(n_vec, i+1, j, sgrid_ptr);
	  glNormal3dv(n_vec);
	  glVertex3d(x(i+1,j), y(i+1,j), z(i+1,j));
/* end the line segment */
	  glEnd();
	}/* end if flag != 0 */
      } /* end for i... */
    } /* end for j... */

    glDisable( GL_CULL_FACE ); 
  } /* end plot_mode & 64 */

  } /* end for all discrete surfaces */

  /* hole curves */
  if (plot_mode & 1024)
    {

/* set the line width. Make it a bit thicker than the boundary lines */
      glLineWidth(1.0);
      glPointSize(3.0);

      for (hole_curve_link = overg->edge_curves->first; hole_curve_link != NULL;
	   hole_curve_link = hole_curve_link->next)
	{
	  hole_curve = hole_curve_link->data;
/* set color */
	  ogl_set_color( hole_curve->edge );

	  glBegin(GL_LINE_STRIP);
	  for (i=0; i<hole_curve->n_points; i++)
	    {
	      glVertex3d(hole_curve->x_edge[i], hole_curve->y_edge[i], 
			 hole_curve->z_edge[i]);
	    } /* end for i... */
	  glEnd();

	  glBegin(GL_POINTS);
	  for (i=0; i<hole_curve->n_points; i++)
	    {
	      glVertex3d(hole_curve->x_edge[i], hole_curve->y_edge[i], 
			 hole_curve->z_edge[i]);
	    } /* end for i... */
	  glEnd();

/* also mark all inside points */
	  glBegin(GL_POINTS);
	  for (ip=0; ip<2; ip++)
	    for (i=0; i<hole_curve->n_points; i++)
	      {
		glVertex3d(hole_curve->inside_point[ip][0][i], 
			   hole_curve->inside_point[ip][1][i], 
			   hole_curve->inside_point[ip][2][i]);
	      } /* end for i... */
	  glEnd();

/* draw the normals */
	  for (ip=0; ip<2; ip++)
	    for (i=0; i<hole_curve->n_points; i++)
	      {
		glBegin(GL_LINES);
		glVertex3d(hole_curve->inside_point[ip][0][i], 
			   hole_curve->inside_point[ip][1][i], 
			   hole_curve->inside_point[ip][2][i]);
		glVertex3d(hole_curve->inside_point[ip][0][i]+
			   0.1*hole_curve->inside_normal[ip][0][i], 
			   hole_curve->inside_point[ip][1][i]+ 
			   0.1*hole_curve->inside_normal[ip][1][i], 
			   hole_curve->inside_point[ip][2][i]+
			   0.1*hole_curve->inside_normal[ip][2][i]); 
		glEnd();
	      } /* end for i... */
      
	} /* end for all hole cutting curves */

      glPointSize(1.0);

    } /* end plot_mode & 1024 */

}


int
hole_plot_mode(input_output *io_ptr, hole_surface *overg, int plot_mode){
  int quit=0, replot=1, icom;
  linked_list_member *this_sgrid;
  discrete_surface *sgrid_ptr;

#include "hs_plot_mode_com.h"

  do{

    if (replot){
/* plot the patched grid */
      if (ogl_start_plot(OGL_NEW_PLOT, OVERLAP_W, 
			 ogl_length_scale(overg->bb))){
	draw_hole_surface( overg, plot_mode );
	if (plot_mode & 256) ogl_bb_cage(overg->bb, OGL_WHITE);
	ogl_end_plot();
	replot = 0;
      }
    }

    switch( get_command(io_ptr, "hole-surface plot mode>", COMMAND, 
			LAST_COM, LEVEL, SAVE_ON_COPY, ARGUMENT) ){

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

    case STD_VIEW:
      ogl_standard_view( OVERLAP_W, overg->bb );
      replot = 1;
      break;

    case ALL_LINES:
      if (plot_mode & 1)
	plot_mode &= ~1;
      else
	plot_mode |= 1;
      printf("Toggled %s.\n", (plot_mode & 1)? "on" : "off");
      replot = 1;
      break;
      
    case ARROWS:
      if (plot_mode & 2)
	plot_mode &= ~2;
      else{
	plot_mode |= 2;
      }
      printf("Toggled %s.\n", (plot_mode & 2)? "on" : "off");
      replot = 1;
      break;

    case BOUNDARIES:
      if (plot_mode & 4)
	plot_mode &= ~4;
      else{
	plot_mode |= 4;
      }
      printf("Toggled %s.\n", (plot_mode & 4)? "on" : "off");
      replot = 1;
      break;

    case FLAG:
      if (plot_mode & 8)
	plot_mode &= ~8;
      else{
	plot_mode |= 8;
      }
      printf("Toggled %s.\n", (plot_mode & 8)? "on" : "off");
      replot = 1;
      break;

    case CURVE_LABEL:
      if (plot_mode & 32)
	plot_mode &= ~32;
      else{
	plot_mode |= 32;
	plot_mode &= ~16;
      }
      printf("Toggled %s.\n", (plot_mode & 32)? "on" : "off");
      replot = 1;
      break;
      
    case ACTIVE_POINTS:
      if (plot_mode & 64)
	plot_mode &= ~64;
      else
	plot_mode |= 64;
      printf("Toggled %s.\n", (plot_mode & 64)? "on" : "off");
      replot = 1;
      break;
      
    case COORD_CAGE:
      if (plot_mode & 256)
	plot_mode &= ~256;
      else
	plot_mode |= 256;
      printf("Toggled %s.\n", (plot_mode & 256)? "on" : "off");
      replot = 1;
      break;
      
    case EDGE_LABEL:
      if (plot_mode & 1024)
	plot_mode &= ~1024;
      else{
	plot_mode |= 1024;
      }
      printf("Toggled %s.\n", (plot_mode & 1024)? "on" : "off");
      replot = 1;
      break;
      
/*     case TRIANGLES: */
/*       if (plot_mode & 2048) */
/* 	plot_mode &= ~2048; */
/*       else{ */
/* 	plot_mode |= 2048; */
/*       } */
/*       printf("Toggled %s.\n", (plot_mode & 2048)? "on" : "off"); */
/*       replot = 1; */
/*       break; */
      
    case SHOW_ONE_COMPONENT:
/* first remove all components */
      for (this_sgrid = overg->grid_list->first; this_sgrid != NULL;
	   this_sgrid = this_sgrid->next){
	sgrid_ptr = this_sgrid->data;
	sgrid_ptr->plot_it = 0;
      }
/* then add one */
      if ((this_sgrid = get_discrete_surface( io_ptr, overg->grid_list )) != NULL){
	sgrid_ptr = this_sgrid->data;
	sgrid_ptr->plot_it = 1;
	replot = 1;
      }
      break;      

    case ADD_ONE_COMPONENT:
      if ((this_sgrid = get_discrete_surface( io_ptr, overg->grid_list )) != NULL){
	sgrid_ptr = this_sgrid->data;
	sgrid_ptr->plot_it = 1;
	replot = 1;
      }
      break;      

    case REMOVE_ONE_COMPONENT:
      if ((this_sgrid = get_discrete_surface( io_ptr, overg->grid_list )) != NULL){
	sgrid_ptr = this_sgrid->data;
	sgrid_ptr->plot_it = 0;
	replot = 1;
      }
      break;      

    case SHOW_ALL_COMPONENTS:
      for (this_sgrid = overg->grid_list->first; this_sgrid != NULL;
	   this_sgrid = this_sgrid->next){
	sgrid_ptr = this_sgrid->data;
	sgrid_ptr->plot_it = 1;
      }
      replot = 1;
      break;      

    case REMOVE_ALL_COMPONENTS:
      for (this_sgrid = overg->grid_list->first; this_sgrid != NULL;
	   this_sgrid = this_sgrid->next){
	sgrid_ptr = this_sgrid->data;
	sgrid_ptr->plot_it = 0;
      }
      replot = 1;
      break;      

    default:
      break;
    }
  } while (!quit);

  return plot_mode;

/* done */
}

void
test_text_normal(void){
  double x_pos[3], n_vec[3];
  glFrontFace( GL_CCW );

/* draw a marker */
  ogl_set_color( OGL_GREEN );

  n_vec[0] = 1.0;
  n_vec[1] = 0.0;
  n_vec[2] = 0.0;
/* set the raster coordinate */
  x_pos[0] = 0.0;
  x_pos[1] = 0.0;
  x_pos[2] = 0.0;
/* draw the string */
  ogl_drawstr(x_pos, n_vec, "Test");

}
