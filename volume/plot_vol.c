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
#include "volume_internal.h"

static int
face_normal(real *n_vec, grid_point *gp_ptr, int direction, int orientation){
  real length, sign;
  int i;
  const real eps = 1.e-10;

  sign = (orientation)? -1.0 : 1.0;

  if (direction == 1){ /* r1 = const */
    n_vec[0] = sign * (gp_ptr->ys * gp_ptr->zt - gp_ptr->zs * gp_ptr->yt);
    n_vec[1] = sign * (gp_ptr->zs * gp_ptr->xt - gp_ptr->xs * gp_ptr->zt);
    n_vec[2] = sign * (gp_ptr->xs * gp_ptr->yt - gp_ptr->ys * gp_ptr->xt);
  }
  else if (direction == 2) { /* r2=const */
    n_vec[0] = sign * (gp_ptr->yt * gp_ptr->zr - gp_ptr->zt * gp_ptr->yr);
    n_vec[1] = sign * (gp_ptr->zt * gp_ptr->xr - gp_ptr->xt * gp_ptr->zr);
    n_vec[2] = sign * (gp_ptr->xt * gp_ptr->yr - gp_ptr->yt * gp_ptr->xr);
  }
  else { /* r3=const */
    n_vec[0] = sign * (gp_ptr->yr * gp_ptr->zs - gp_ptr->zr * gp_ptr->ys);
    n_vec[1] = sign * (gp_ptr->zr * gp_ptr->xs - gp_ptr->xr * gp_ptr->zs);
    n_vec[2] = sign * (gp_ptr->xr * gp_ptr->ys - gp_ptr->yr * gp_ptr->xs);
  }

/* normalize */
  length = sqrt( n_vec[0] * n_vec[0] + n_vec[1] * n_vec[1] + n_vec[2] * n_vec[2] );
  if (length > eps){
    for (i=0; i<3; i++) n_vec[i] = n_vec[i]/length;
    return OK;
  }
  else
    return ERROR;
}

/* static void  */
/* gridmaterials(void) */
/* { */
/*   GLfloat back_mat_diffuse[] = { 0.64, 0.36, 0.0, 1.0 }; */
/*   GLfloat back_mat_ambient[] = { 0.64, 0.36, 0.0, 1.0 }; */

/*    glMaterialfv(GL_BACK, GL_DIFFUSE, back_mat_diffuse);  */
/*    glMaterialfv(GL_BACK, GL_AMBIENT, back_mat_ambient);  */

/* } */

/* static void  */
/* surfacematerials(void) */
/* { */
/*   float front_mat_diffuse[] = { 0.2, 0.7, 0.4, 1.0 }; */
/*   float front_mat_ambient[] = { 0.2, 0.2, 0.2, 1.0 }; */
/*   float front_mat_specular[] = { 0.03, 0.03, 0.03, 1.0 }; */
/*   float back_mat_diffuse[] = { 1.0, 1.0, 0.2, 1.0 }; */
/*   float back_mat_ambient[] = { 0.1, 0.1, 0.1, 1.0 }; */

/*   glMaterialfv(GL_FRONT, GL_DIFFUSE, front_mat_diffuse); */
/*   glMaterialfv(GL_FRONT, GL_AMBIENT, front_mat_ambient); */
/*   glMaterialfv(GL_BACK, GL_DIFFUSE, back_mat_diffuse); */
/*   glMaterialfv(GL_BACK, GL_AMBIENT, back_mat_ambient); */

/*   glMaterialfv(GL_FRONT, GL_SPECULAR, front_mat_specular); */
/*   glMateriali(GL_FRONT, GL_SHININESS, 10); */
/* } */

void
draw_volume( volume_mapping *volume, int plot_mode, clip_planes *clip ){
  int i, j, k;
  real hu, hv, hw, length, r_length, s_length, x_pos[3];
  real n_vec[3], start[3], forward[3];
  grid_point gp;
  char label[80];
  GLfloat position[] = { 0.0, 0.0, 1.0, 0.0 };
  const real eps=1.e-10;

/*  plot_mode & 1: plot the grid points (using the forward_mapping) */
/*  plot_mode & 2: draw the boundaries with thick lines */
/*  plot_mode & 4: draw the parameter arrows with labels */
/*  plot_mode & 8: plot grid lines on an r1=const surface */
/*  plot_mode & 16: plot grid lines on an r2=const surface */
/*  plot_mode & 32: plot grid lines on an r3=const surface */
/*  plot_mode & 64: draw the surface labels */
/*  plot_mode & 128: draw the boundary conditions */
/*  plot_mode & 256: plot the coordinate cage */
/*  plot_mode & 512: plot the grid lines on all surfaces with */
/*                   surf_vol = plot_surf_label */
/*  plot_mode & 1024: plot the grid lines on all surfaces with */
/*                    bc_vol = plot_bc */
/*  plot_mode & 2048: draw the edges with edge_vol != 0 with thick lines */

/* use a light which is fixed to the viewer */
  glPushMatrix();
/* position the light */
    glLoadIdentity();
    glLightfv(GL_LIGHT0, GL_POSITION, position);
  glPopMatrix();

/* Default is to use one-sided lighting */
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

/* adjust clipping */
  if (clip)
    {
      for (i=0; i<clip->n_clip_planes; i++)
	{
	  glClipPlane(GL_CLIP_PLANE0+i, clip->clip_plane[i]); 
	  glEnable(GL_CLIP_PLANE0+i); 
	}
/* the maximum number of clip planes is 4 */
      for (i=clip->n_clip_planes; i<4; i++)
	glDisable(GL_CLIP_PLANE0+i);
    }
  else
    {
      for (i=0; i<4; i++)
	glDisable(GL_CLIP_PLANE0+i);
    }
/* grid sizes */
  hu = 1.0/((real) volume->r1_points - 1); 
  hv = 1.0/((real) volume->r2_points - 1); 
  hw = 1.0/((real) volume->r3_points - 1); 

  if (plot_mode & 1) {
/* set surface props */
    ogl_set_color( volume->color);

    glPointSize(2.0);

/* start a new strip */
    glBegin(GL_POINTS);

    for (i=0; i<volume->r1_points; i++){
      gp.r = i*hu;
      for (j=0; j<volume->r2_points; j++){
	gp.s = j*hv;
	for (k=0; k<volume->r3_points; k++){
	  gp.t = k*hw;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 3, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for k */
      } /* end for j... */
    } /* end for i */

/* end the point swarm */
    glEnd();

    glPointSize(1.0);

  } /* end (plot_mode & 1) */

/* draw the boundary with thick lines */
  if (plot_mode & 2){

/* set the color */
    ogl_set_color( volume->color );
/* set the line width */
    glLineWidth(3.0);
    glEnable(GL_LINE_SMOOTH);

    k=0;
    gp.t = 0.0;
    glBegin(GL_LINE_STRIP);
    j=0;
    for (i=0; i<volume->r1_points; i++){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 2, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for i... */

    i=volume->r1_points-1;
    for (j=0; j<volume->r2_points; j++){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 3, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for j... */

    j=volume->r2_points-1;
    for (i = volume->r1_points-1; i>=0; i--){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 2, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for i... */

    i=0;
    for (j=volume->r2_points-1; j>=0; j--){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 3, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for j... */

    glEnd();

    k=volume->r3_points-1;
    gp.t = 1.0;
    glBegin(GL_LINE_STRIP);
    j=0;
    for (i=0; i<volume->r1_points; i++){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 2, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for i... */

    i=volume->r1_points-1;
    for (j=0; j<volume->r2_points; j++){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 3, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for j... */

    j=volume->r2_points-1;
    for (i = volume->r1_points-1; i>=0; i--){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 2, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for i... */

    i=0;
    for (j=volume->r2_points-1; j>=0; j--){
      gp.r = i*hu;
      gp.s = j*hv;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 3, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for j... */

    glEnd();

    gp.r = 0.0;
    gp.s = 0.0;
    glBegin(GL_LINE_STRIP);
    for (k=0; k<volume->r3_points; k++){
      gp.t = k*hw;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 1, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for k... */
    glEnd();

    gp.r = 1.0;
    gp.s = 0.0;
    glBegin(GL_LINE_STRIP);
    for (k=0; k<volume->r3_points; k++){
      gp.t = k*hw;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 1, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for k... */
    glEnd();

    gp.r = 0.0;
    gp.s = 1.0;
    glBegin(GL_LINE_STRIP);
    for (k=0; k<volume->r3_points; k++){
      gp.t = k*hw;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 1, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for k... */
    glEnd();

    gp.r = 1.0;
    gp.s = 1.0;
    glBegin(GL_LINE_STRIP);
    for (k=0; k<volume->r3_points; k++){
      gp.t = k*hw;

      forward_volume_mapping( &gp, volume );
      face_normal(n_vec, &gp, 1, 0);
/* set current normal */
      glNormal3dv(n_vec);
/* set current vertex */
      glVertex3d(gp.x, gp.y, gp.z);

    } /* end for k... */
    glEnd();

/* reset the line width */
    glLineWidth(1.0);
    glDisable(GL_LINE_SMOOTH);
  }  /* end plot_mode & 2 */

/* draw parameter arrows */
  if (plot_mode & 4){
/* set orientation */
    glFrontFace( GL_CCW );

/* set the color */
    ogl_set_color( volume->color );

/* arrow length */
    length = 0.1 * real_max(volume->bb->x_max - volume->bb->x_min, 
			    real_max(volume->bb->y_max - volume->bb->y_min, 
				     volume->bb->z_max - volume->bb->z_min));
/* r1-arrow */
    gp.r = 1.0;
    gp.s = 0.0;
    gp.t = 0.0;
    forward_volume_mapping( &gp, volume );
    face_normal( n_vec, &gp, 1, 0 );
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

    ogl_draw_arrow(length, 3.0, start, forward, n_vec, "r1");

/* r2-arrow */
    gp.r = 0.0;
    gp.s = 1.0;
    gp.t = 0.0;
    forward_volume_mapping( &gp, volume );
    face_normal( n_vec, &gp, 2, 0 );
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

    ogl_draw_arrow(length, 3.0, start, forward, n_vec, "r2");

/* r3-arrow */
    gp.r = 0.0;
    gp.s = 0.0;
    gp.t = 1.0;
    forward_volume_mapping( &gp, volume );
    face_normal( n_vec, &gp, 3, 0 );
    start[0] = gp.x;
    start[1] = gp.y;
    start[2] = gp.z;
    s_length = sqrt( gp.xt*gp.xt + gp.yt*gp.yt + gp.zt*gp.zt );
    if (s_length < eps){
      s_length = 1.0;
    }
    forward[0] = gp.xt/s_length;
    forward[1] = gp.yt/s_length;
    forward[2] = gp.zt/s_length;

    ogl_draw_arrow(length, 3.0, start, forward, n_vec, "r3");

  }  /* end plot_mode & 4 */

/* draw an r1=const grid surface */
  if (plot_mode & 8){

/* set the material props */
    ogl_set_color( volume->color );

/* set the surface */
    gp.r = volume->r1_surface;

    for (k=0; k<volume->r3_points; k++){
/* start a new strip */
      glBegin(GL_LINE_STRIP);
      for (j=0; j<volume->r2_points; j++){
	gp.s = j*hv;
	gp.t = k*hw;

	forward_volume_mapping( &gp, volume );
	face_normal(n_vec, &gp, 1, 0);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

      } /* end for j... */
/* end the strip */
      glEnd();
    }

    for (j=0; j<volume->r2_points; j++){
/* start a new strip */
      glBegin(GL_LINE_STRIP);
      for (k=0; k<volume->r3_points; k++){
	gp.s = j*hv;
	gp.t = k*hw;

	forward_volume_mapping( &gp, volume );
	face_normal(n_vec, &gp, 1, 0);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

      } /* end for j... */
/* end the strip */
      glEnd();
    }
  } /* end plot_mode & 8 */

/* draw a r2=const grid-surface */
  if (plot_mode & 16){

/* set the material props */
    ogl_set_color( volume->color );

/* set the surface */
    gp.s = volume->r2_surface;

    for (i=0; i<volume->r1_points; i++){
/* start a new strip */
      glBegin(GL_LINE_STRIP);
      for (k=0; k<volume->r3_points; k++){
	gp.r = i*hu;
	gp.t = k*hw;

	forward_volume_mapping( &gp, volume );
	face_normal(n_vec, &gp, 2, 0);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

      } /* end for j... */
/* end the strip */
      glEnd();
    }

    for (k=0; k<volume->r3_points; k++){
/* start a new strip */
      glBegin(GL_LINE_STRIP);
      for (i=0; i<volume->r1_points; i++){
	gp.r = i*hu;
	gp.t = k*hw;

	forward_volume_mapping( &gp, volume );
	face_normal(n_vec, &gp, 2, 0);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

      } /* end for j... */
/* end the strip */
      glEnd();
    }
  } /* end plot_mode & 16 */

/* draw an r3=const grid-surface */
  if (plot_mode & 32){

/* set the material props */
    ogl_set_color( volume->color );

/* set the surface */
    gp.t = volume->r3_surface;

    for (i=0; i<volume->r1_points; i++){
/* start a new strip */
      glBegin(GL_LINE_STRIP);
      for (j=0; j<volume->r2_points; j++){
	gp.r = i*hu;
	gp.s = j*hv;

	forward_volume_mapping( &gp, volume );
	face_normal(n_vec, &gp, 3, 0);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

      } /* end for j... */
/* end the strip */
      glEnd();
    }

    for (j=0; j<volume->r2_points; j++){
/* start a new strip */
      glBegin(GL_LINE_STRIP);
      for (i=0; i<volume->r1_points; i++){
	gp.r = i*hu;
	gp.s = j*hv;

	forward_volume_mapping( &gp, volume );
	face_normal(n_vec, &gp, 3, 0);
/* set current normal */
	glNormal3dv(n_vec);
/* set current vertex */
	glVertex3d(gp.x, gp.y, gp.z);

      } /* end for j... */
/* end the strip */
      glEnd();
    }
  } /* end plot_mode & 32 */

/* draw the surface labels */
  if (plot_mode & 64){

/* set color */
    ogl_set_color( volume->color );
/* size of component grid */
    length = 0.05 * real_max(volume->bb->x_max - volume->bb->x_min, 
			     real_max(volume->bb->y_max - volume->bb->y_min, 
				      volume->bb->z_max - volume->bb->z_min));
/* left */
    gp.r = 0.0;
    gp.s = 0.5;
    gp.t = 0.5;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 1, 1);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "surface(r1=0):%i", surf_vol(1,1));
    ogl_drawstr(x_pos, n_vec, label);

/* right */
    gp.r = 1.0;
    gp.s = 0.5;
    gp.t = 0.5;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 1, 0);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "surface(r1=1):%i", surf_vol(2,1));
    ogl_drawstr(x_pos, n_vec, label);

/* lower */
    gp.r = 0.5;
    gp.s = 0.0;
    gp.t = 0.5;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 2, 1);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "surface(r2=0):%i", surf_vol(1,2));
    ogl_drawstr(x_pos, n_vec, label);

/* upper */
    gp.r = 0.5;
    gp.s = 1.0;
    gp.t = 0.5;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 2, 0);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "surface(r2=1):%i", surf_vol(2,2));
    ogl_drawstr(x_pos, n_vec, label);

/* far */
    gp.r = 0.5;
    gp.s = 0.5;
    gp.t = 0.0;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 3, 1);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "surface(r3=0):%i", surf_vol(1,3));
    ogl_drawstr(x_pos, n_vec, label);

/* near */
    gp.r = 0.5;
    gp.s = 0.5;
    gp.t = 1.0;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 3, 0);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "surface(r3=1):%i", surf_vol(2,3));
    ogl_drawstr(x_pos, n_vec, label);

  } /* end plot_mode & 64 */

/* draw the boundary conditions */
  if (plot_mode & 128){

/* set color */
    ogl_set_color( volume->color );
/* size of component grid */
    length = 0.05 * real_max(volume->bb->x_max - volume->bb->x_min, 
			     real_max(volume->bb->y_max - volume->bb->y_min, 
				      volume->bb->z_max - volume->bb->z_min));
/* left */
    gp.r = 0.0;
    gp.s = 0.5;
    gp.t = 0.5;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 1, 1);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "bc(r1=0):%i", bc_vol(1,1));
    ogl_drawstr(x_pos, n_vec, label);

/* right */
    gp.r = 1.0;
    gp.s = 0.5;
    gp.t = 0.5;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 1, 0);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "bc(r1=1):%i", bc_vol(2,1));
    ogl_drawstr(x_pos, n_vec, label);

/* lower */
    gp.r = 0.5;
    gp.s = 0.0;
    gp.t = 0.5;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 2, 1);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "bc(r2=0):%i", bc_vol(1,2));
    ogl_drawstr(x_pos, n_vec, label);

/* upper */
    gp.r = 0.5;
    gp.s = 1.0;
    gp.t = 0.5;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 2, 0);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "bc(r2=1):%i", bc_vol(2,2));
    ogl_drawstr(x_pos, n_vec, label);

/* far */
    gp.r = 0.5;
    gp.s = 0.5;
    gp.t = 0.0;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 3, 1);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "bc(r3=0):%i", bc_vol(1,3));
    ogl_drawstr(x_pos, n_vec, label);

/* near */
    gp.r = 0.5;
    gp.s = 0.5;
    gp.t = 1.0;
    forward_volume_mapping( &gp, volume );
    face_normal(n_vec, &gp, 3, 0);
/* position the text */		    
    x_pos[0] = gp.x;
    x_pos[1] = gp.y;
    x_pos[2] = gp.z;
    sprintf(label, "bc(r3=1):%i", bc_vol(2,3));
    ogl_drawstr(x_pos, n_vec, label);

  } /* end plot_mode & 128 */

/* draw all faces with surf_vol = plot_surf_label */
  if (plot_mode & 512){

/* two-sided lighting */
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
/* set surface props */
    ogl_set_color( volume->plot_surf_label );
/* both sides matter because of the two-sided lighting */
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

/* check r1=0 */
    if (surf_vol(1,1) == volume->plot_surf_label){
      glFrontFace(GL_CW);
      gp.r = 0.0;
      for (j=0; j<volume->r2_points-1; j++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (k=0; k<volume->r3_points; k++){
/* first point */
	  gp.s = (j+1)*hv;
	  gp.t = k*hw;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 1, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.s = j*hv;
	  gp.t = k*hw;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 1, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for k... */
/* end the strip */
	glEnd();
      } /* end for j... */
    }/* end r1=0 */

/* check r1=1 */
    if (surf_vol(2,1) == volume->plot_surf_label){
      glFrontFace(GL_CCW);
      gp.r = 1.0;
      for (j=0; j<volume->r2_points-1; j++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (k=0; k<volume->r3_points; k++){
/* first point */
	  gp.s = (j+1)*hv;
	  gp.t = k*hw;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 1, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.s = j*hv;
	  gp.t = k*hw;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 1, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for k... */
/* end the strip */
	glEnd();
      } /* end for j... */
    }/* end r1=1 */

/* check r2=0 */
    if (surf_vol(1,2) == volume->plot_surf_label){
      glFrontFace(GL_CW);
      gp.s = 0.0;
      for (k=0; k<volume->r3_points-1; k++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (i=0; i<volume->r1_points; i++){
/* first point */
	  gp.t = (k+1)*hw;
	  gp.r = i*hu;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 2, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.t = k*hw;
	  gp.r = i*hu;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 2, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for i... */
/* end the strip */
	glEnd();
      } /* end for k... */
    }/* end r2=0 */

/* check r2=1 */
    if (surf_vol(2,2) == volume->plot_surf_label){
      glFrontFace(GL_CCW);
      gp.s = 1.0;
      for (k=0; k<volume->r3_points-1; k++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (i=0; i<volume->r1_points; i++){
/* first point */
	  gp.t = (k+1)*hw;
	  gp.r = i*hu;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 2, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.t = k*hw;
	  gp.r = i*hu;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 2, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for i... */
/* end the strip */
	glEnd();
      } /* end for k... */
    }/* end r2=1 */

/* check r3=0 */
    if (surf_vol(1,3) == volume->plot_surf_label){
      glFrontFace(GL_CW);
      gp.t = 0.0;
      for (i=0; i<volume->r1_points-1; i++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (j=0; j<volume->r2_points; j++){
/* first point */
	  gp.r = (i+1)*hu;
	  gp.s = j*hv;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 3, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.r = i*hu;
	  gp.s = j*hv;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 3, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for j... */
/* end the strip */
	glEnd();
      } /* end for i... */
    }/* end r3=0 */

/* check r3=1 */
    if (surf_vol(2,3) == volume->plot_surf_label){
      glFrontFace(GL_CCW);
      gp.t = 1.0;
      for (i=0; i<volume->r1_points-1; i++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (j=0; j<volume->r2_points; j++){
/* first point */
	  gp.r = (i+1)*hu;
	  gp.s = j*hv;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 3, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.r = i*hu;
	  gp.s = j*hv;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 3, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for j... */
/* end the strip */
	glEnd();
      } /* end for i... */
    }/* end r3=1 */

  } /* end plot_mode & 512 */

/* draw all faces with bc_vol = plot_bc */
  if (plot_mode & 1024){

/* two-sided lighting */
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
/* set surface props */
    ogl_set_color( volume->plot_bc );
/* both sides matter because of the two-sided lighting */
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

/* check r1=0 */
    if (bc_vol(1,1) == volume->plot_bc){
      glFrontFace(GL_CW);
      gp.r = 0.0;
      for (j=0; j<volume->r2_points-1; j++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (k=0; k<volume->r3_points; k++){
/* first point */
	  gp.s = (j+1)*hv;
	  gp.t = k*hw;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 1, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.s = j*hv;
	  gp.t = k*hw;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 1, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for k... */
/* end the strip */
	glEnd();
      } /* end for j... */
    }/* end r1=0 */

/* check r1=1 */
    if (bc_vol(2,1) == volume->plot_bc){
      glFrontFace(GL_CCW);
      gp.r = 1.0;
      for (j=0; j<volume->r2_points-1; j++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (k=0; k<volume->r3_points; k++){
/* first point */
	  gp.s = (j+1)*hv;
	  gp.t = k*hw;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 1, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.s = j*hv;
	  gp.t = k*hw;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 1, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for k... */
/* end the strip */
	glEnd();
      } /* end for j... */
    }/* end r1=1 */

/* check r2=0 */
    if (bc_vol(1,2) == volume->plot_bc){
      glFrontFace(GL_CW);
      gp.s = 0.0;
      for (k=0; k<volume->r3_points-1; k++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (i=0; i<volume->r1_points; i++){
/* first point */
	  gp.t = (k+1)*hw;
	  gp.r = i*hu;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 2, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.t = k*hw;
	  gp.r = i*hu;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 2, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for i... */
/* end the strip */
	glEnd();
      } /* end for k... */
    }/* end r2=0 */

/* check r2=1 */
    if (bc_vol(2,2) == volume->plot_bc){
      glFrontFace(GL_CCW);
      gp.s = 1.0;
      for (k=0; k<volume->r3_points-1; k++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (i=0; i<volume->r1_points; i++){
/* first point */
	  gp.t = (k+1)*hw;
	  gp.r = i*hu;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 2, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.t = k*hw;
	  gp.r = i*hu;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 2, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for i... */
/* end the strip */
	glEnd();
      } /* end for k... */
    }/* end r2=1 */

/* check r3=0 */
    if (bc_vol(1,3) == volume->plot_bc){
      glFrontFace(GL_CW);
      gp.t = 0.0;
      for (i=0; i<volume->r1_points-1; i++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (j=0; j<volume->r2_points; j++){
/* first point */
	  gp.r = (i+1)*hu;
	  gp.s = j*hv;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 3, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.r = i*hu;
	  gp.s = j*hv;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 3, 0);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for j... */
/* end the strip */
	glEnd();
      } /* end for i... */
    }/* end r3=0 */

/* check r3=1 */
    if (bc_vol(2,3) == volume->plot_bc){
      glFrontFace(GL_CCW);
      gp.t = 1.0;
      for (i=0; i<volume->r1_points-1; i++){
/* start a new strip */
	glBegin(GL_QUAD_STRIP);
	for (j=0; j<volume->r2_points; j++){
/* first point */
	  gp.r = (i+1)*hu;
	  gp.s = j*hv;
	  
	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 3, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

/* second point */
	  gp.r = i*hu;
	  gp.s = j*hv;

	  forward_volume_mapping( &gp, volume );
	  face_normal(n_vec, &gp, 3, 1);
/* set current normal */
	  glNormal3dv(n_vec);
/* set current vertex */
	  glVertex3d(gp.x, gp.y, gp.z);

	} /* end for j... */
/* end the strip */
	glEnd();
      } /* end for i... */
    }/* end r3=1 */

  } /* end plot_mode & 1024 */

/* edge label */
  if (plot_mode & 2048)
    {

/* set the line width. Make it a bit thicker than the boundary lines */
      glLineWidth(4.0);
      glEnable(GL_LINE_SMOOTH);

/* edge 1: r1=0, r2=0 */
      if (edge_vol(1))
	{
/* set color */
	  ogl_set_color( edge_vol(1) );
	  gp.r = 0.0;
	  gp.s = 0.0;
	  glBegin(GL_LINE_STRIP);
	  for (k=0; k<volume->r3_points; k++)
	    {
	      gp.t = k*hw;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 1, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for k... */
	  glEnd();
	} /* end if edge_vol(1) */

/* edge 2: r1=1, r2=0 */
      if (edge_vol(2))
	{
/* set color */
	  ogl_set_color( edge_vol(2) );
	  gp.r = 1.0;
	  gp.s = 0.0;
	  glBegin(GL_LINE_STRIP);
	  for (k=0; k<volume->r3_points; k++)
	    {
	      gp.t = k*hw;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 1, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for k... */
	  glEnd();
	} /* end if edge_vol(2) */

/* edge 3: r1=0, r2=1 */
      if (edge_vol(3))
	{
/* set color */
	  ogl_set_color( edge_vol(3) );
	  gp.r = 0.0;
	  gp.s = 1.0;
	  glBegin(GL_LINE_STRIP);
	  for (k=0; k<volume->r3_points; k++)
	    {
	      gp.t = k*hw;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 1, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for k... */
	  glEnd();
	} /* end if edge_vol(3) */

/* edge 4: r1=1, r2=1 */
      if (edge_vol(4))
	{
/* set color */
	  ogl_set_color( edge_vol(4) );
	  gp.r = 1.0;
	  gp.s = 1.0;
	  glBegin(GL_LINE_STRIP);
	  for (k=0; k<volume->r3_points; k++)
	    {
	      gp.t = k*hw;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 1, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for k... */
	  glEnd();
	} /* end if edge_vol(4) */

/* edge 5: r1=0, r3=0 */
      if (edge_vol(5))
	{
/* set color */
	  ogl_set_color( edge_vol(5) );
	  gp.r = 0.0;
	  gp.t = 0.0;
	  glBegin(GL_LINE_STRIP);
	  for (j=0; j<volume->r2_points; j++)
	    {
	      gp.s = j*hv;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 3, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for j... */
	  glEnd();
	} /* end if edge_vol(5) */

/* edge 6: r1=1, r3=0 */
      if (edge_vol(6))
	{
/* set color */
	  ogl_set_color( edge_vol(6) );
	  gp.r = 1.0;
	  gp.t = 0.0;
	  glBegin(GL_LINE_STRIP);
	  for (j=0; j<volume->r2_points; j++)
	    {
	      gp.s = j*hv;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 3, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for j... */
	  glEnd();
	} /* end if edge_vol(6) */

/* edge 7: r1=1, r3=0 */
      if (edge_vol(7))
	{
/* set color */
	  ogl_set_color( edge_vol(7) );
	  gp.r = 0.0;
	  gp.t = 1.0;
	  glBegin(GL_LINE_STRIP);
	  for (j=0; j<volume->r2_points; j++)
	    {
	      gp.s = j*hv;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 3, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for j... */
	  glEnd();
	} /* end if edge_vol(7) */

/* edge 8: r1=1, r3=1 */
      if (edge_vol(8))
	{
/* set color */
	  ogl_set_color( edge_vol(8) );
	  gp.r = 1.0;
	  gp.t = 1.0;
	  glBegin(GL_LINE_STRIP);
	  for (j=0; j<volume->r2_points; j++)
	    {
	      gp.s = j*hv;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 3, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for j... */
	  glEnd();
	} /* end if edge_vol(8) */

/* edge 9: r2=0, r3=0 */
      if (edge_vol(9))
	{
/* set color */
	  ogl_set_color( edge_vol(9) );
	  gp.s = 0.0;
	  gp.t = 0.0;
	  glBegin(GL_LINE_STRIP);
	  for (i=0; i<volume->r1_points; i++)
	    {
	      gp.r = i*hu;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 2, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for i... */
	  glEnd();
	} /* end if edge_vol(9) */

/* edge 10: r2=1, r3=0 */
      if (edge_vol(10))
	{
/* set color */
	  ogl_set_color( edge_vol(10) );
	  gp.s = 1.0;
	  gp.t = 0.0;
	  glBegin(GL_LINE_STRIP);
	  for (i=0; i<volume->r1_points; i++)
	    {
	      gp.r = i*hu;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 2, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for i... */
	  glEnd();
	} /* end if edge_vol(10) */

/* edge 11: r2=0, r3=1 */
      if (edge_vol(11))
	{
/* set color */
	  ogl_set_color( edge_vol(11) );
	  gp.s = 0.0;
	  gp.t = 1.0;
	  glBegin(GL_LINE_STRIP);
	  for (i=0; i<volume->r1_points; i++)
	    {
	      gp.r = i*hu;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 2, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for i... */
	  glEnd();
	} /* end if edge_vol(11) */

/* edge 12: r2=1, r3=1 */
      if (edge_vol(12))
	{
/* set color */
	  ogl_set_color( edge_vol(12) );
	  gp.s = 1.0;
	  gp.t = 1.0;
	  glBegin(GL_LINE_STRIP);
	  for (i=0; i<volume->r1_points; i++)
	    {
	      gp.r = i*hu;

	      forward_volume_mapping( &gp, volume );
	      face_normal(n_vec, &gp, 2, 0);
	      /* set current normal */
	      glNormal3dv(n_vec);
	      /* set current vertex */
	      glVertex3d(gp.x, gp.y, gp.z);

	    } /* end for i... */
	  glEnd();
	} /* end if edge_vol(12) */

/* reset the line width */
      glLineWidth(1.0);
      glDisable(GL_LINE_SMOOTH);
    } /* end plot_mode & 2048 */

/* turn all 4 clip planes off */
  for (i=0; i<4; i++)
    glDisable(GL_CLIP_PLANE0+i);
}

