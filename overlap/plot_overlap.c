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
#include "overlap_internal.h"

int
grid_face_normal(real *n_vec, component_vgrid *vgrid, int i, int j, int k, 
		 int direction, int orientation){
  real length, sign, xr, xs, xt, yr, ys, yt, zr, zs, zt;
  int q;
  const real eps = 1.e-10;

  sign = (orientation)? -1.0 : 1.0;

/* r-differences */
  if (i == range3(1,1)){
    xr = x3(i+1,j,k) - x3(i,j,k);
    yr = y3(i+1,j,k) - y3(i,j,k);
    zr = z3(i+1,j,k) - z3(i,j,k);
  }
  else if (i== range3(2,1)){
    xr = x3(i,j,k) - x3(i-1,j,k);
    yr = y3(i,j,k) - y3(i-1,j,k);
    zr = z3(i,j,k) - z3(i-1,j,k);
  }
  else{
    xr = x3(i+1,j,k) - x3(i-1,j,k);
    yr = y3(i+1,j,k) - y3(i-1,j,k);
    zr = z3(i+1,j,k) - z3(i-1,j,k);
  }

/* s-differences */
  if (j == range3(1,2)){
    xs = x3(i,j+1,k) - x3(i,j,k);
    ys = y3(i,j+1,k) - y3(i,j,k);
    zs = z3(i,j+1,k) - z3(i,j,k);
  }
  else if (j == range3(2,2)){
    xs = x3(i,j,k) - x3(i,j-1,k);
    ys = y3(i,j,k) - y3(i,j-1,k);
    zs = z3(i,j,k) - z3(i,j-1,k);
  }
  else{
    xs = x3(i,j+1,k) - x3(i,j-1,k);
    ys = y3(i,j+1,k) - y3(i,j-1,k);
    zs = z3(i,j+1,k) - z3(i,j-1,k);
  }

/* t-differences */
  if (k == range3(1,3)){
    xt = x3(i,j,k+1) - x3(i,j,k);
    yt = y3(i,j,k+1) - y3(i,j,k);
    zt = z3(i,j,k+1) - z3(i,j,k);
  }
  else if (k == range3(2,3)){
    xt = x3(i,j,k) - x3(i,j,k-1);
    yt = y3(i,j,k) - y3(i,j,k-1);
    zt = z3(i,j,k) - z3(i,j,k-1);
  }
  else{
    xt = x3(i,j,k+1) - x3(i,j,k-1);
    yt = y3(i,j,k+1) - y3(i,j,k-1);
    zt = z3(i,j,k+1) - z3(i,j,k-1);
  }

  if (direction == 1){ /* r1 = const */
    n_vec[0] = sign * (ys * zt - zs * yt);
    n_vec[1] = sign * (zs * xt - xs * zt);
    n_vec[2] = sign * (xs * yt - ys * xt);
  }
  else if (direction == 2) { /* r2=const */
    n_vec[0] = sign * (yt * zr - zt * yr);
    n_vec[1] = sign * (zt * xr - xt * zr);
    n_vec[2] = sign * (xt * yr - yt * xr);
  }
  else { /* r3=const */
    n_vec[0] = sign * (yr * zs - zr * ys);
    n_vec[1] = sign * (zr * xs - xr * zs);
    n_vec[2] = sign * (xr * ys - yr * xs);
  }

/* normalize */
  length = sqrt( n_vec[0] * n_vec[0] + n_vec[1] * n_vec[1] + n_vec[2] * n_vec[2] );
  if (length > eps){
    for (q=0; q<3; q++) n_vec[q] = n_vec[q]/length;
    return OK;
  }
  else
    return ERROR;
}

void
draw_overlapping_3d_grid( overlapping_3d_grid *over3d, int plot_mode ){
  int i, j, k, orientation;
  real n_vec[3], bb_length;
  GLfloat position[] = { 0.0, 0.0, 1.0, 0.0 };
  bad_3d_point *bad3;
  interp_3d_point *interp3;
  component_vgrid *vgrid;
  linked_list_member *this_vgrid;

/*  plot_mode & 1: plot all used boundary faces */
/*  plot_mode & 4: plot the change surfaces (the grid surfaces separating 
                                             the used and unused grid points */
/*  plot_mode & 8: plot all used grid cells */
/*  plot_mode & 16: plot all used grid faces */
/*  plot_mode & 32: plot all used grid lines */
/*  plot_mode & 64: plot all interpolation faces */
/*  plot_mode & 128: plot all bad grid points */
/*  plot_mode & 512: plot all interpolation points */
/*  plot_mode & 1024: plot all dummy interpolation points (disabled) */
/*  plot_mode & 2048: plot all dead points */
/*  plot_mode & 4096: plot all faces with surf3>0 */
/*  plot_mode & 8192: plot the grid surface with flag = 'w' */

/* use a light which is fixed to the viewer */
  glPushMatrix();
/* position the light */
    glLoadIdentity();
    glLightfv(GL_LIGHT0, GL_POSITION, position);
  glPopMatrix();

/* Default is to use a two-sided lighting model */
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

/* orientation test: cull backfacing polygons */
/*   glCullFace( GL_BACK ); */
/*   glEnable( GL_CULL_FACE ); */
  
/* adjust clipping */
  for (i=0; i<over3d->n_clip_planes; i++){
    glClipPlane(GL_CLIP_PLANE0+i, over3d->clip_plane[i]); 
    glEnable(GL_CLIP_PLANE0+i); 
  }
/* the maximum number of clip planes is 4 */
  for (i=over3d->n_clip_planes; i<4; i++)
    glDisable(GL_CLIP_PLANE0+i);

/* compute length scale */
  bb_length = (over3d->bb->x_max - over3d->bb->x_min +
	       over3d->bb->y_max - over3d->bb->y_min +
	       over3d->bb->z_max - over3d->bb->z_min)/3.0;

/* loop over all component grids */
  for (this_vgrid = over3d->grid_list->first; this_vgrid != NULL; 
       this_vgrid = this_vgrid->next){ 
    vgrid = this_vgrid->data; 

/* don't do this for dummy grids */
/*    if (vgrid->dummy_background) continue;*/

/* check if this component should be plotted */
    if (vgrid->plot_it){

/* set the orientation */
    glFrontFace(GL_CW);

/* plot all used boundary points */
    if (plot_mode & 1){

/* set surface properties */
      ogl_set_color( vgrid->color);

/* loop over all faces */
/* face k=range3(side,3) */
      for (k=range3(1,3); k<=range3(2,3); k += (range3(2,3)-range3(1,3))){
/*	orientation = (k==range3(1,3))? 1 : 0;*/
	orientation = 1;
	for (j=range3(1,2); j<range3(2,2); j++){
	  for (i=range3(1,1); i<range3(2,1); i++){
	    if (flag3(i,j  ,k) != 0 && flag3(i+1,j  ,k) != 0 &&
		flag3(i,j+1,k) != 0 && flag3(i+1,j+1,k) != 0){
/* start a new polygon */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* end the polygon */
	      glEnd();
	    }
	  } /* end for i... */
	} /* end for j... */
      } /* end for k... */

/* face i=range3(side,1) */
      for (i=range3(1,1); i<=range3(2,1); i += (range3(2,1)-range3(1,1))){
/*	orientation = (i==range3(1,1))? 1 : 0;*/
	orientation = 1;
	for (k=range3(1,3); k<range3(2,3); k++){
	  for (j=range3(1,2); j<range3(2,2); j++){
	    if (flag3(i,j,k  ) != 0 && flag3(i,j+1,k  ) != 0 &&
		flag3(i,j,k+1) != 0 && flag3(i,j+1,k+1) != 0){
/* start a new polygon */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k+1, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* end the polygon */
	      glEnd();
	    }
	  } /* end for j... */
	} /* end for k ... */
      } /* end face i */

/* face j=range3(side,2) */
      for (j=range3(1,2); j<=range3(2,2); j += (range3(2,2)-range3(1,2))){
/*	orientation = (j==range3(1,2))? 1 : 0;*/
	orientation = 1;
	for (k=range3(1,3); k<range3(2,3); k++){
	  for (i=range3(1,1); i<range3(2,1); i++){
	    if (flag3(i  ,j,k  ) != 0 && flag3(i  ,j,k+1) != 0 &&
		flag3(i+1,j,k+1) != 0 && flag3(i+1,j,k  ) != 0){
/* start a new polygon */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* end the polygon */
	      glEnd();
	    } 
	  } /* end for i... */
	} /* end for k... */
      } /* end face j */
    }/* end plot_mode & 1 */

/* draw the surface separating the used and unused grid points */
/* check that the essential arrays are there */
    if (plot_mode & 4){

/* set surface properties */
      ogl_set_color( vgrid->color);

/* loop over all grid cells */
      for (k=range3(1,3); k<range3(2,3); k++){
	for (j=range3(1,2); j<range3(2,2); j++){
	  for (i=range3(1,1); i<range3(2,1); i++){
/* check if this cell contains both dead and alive vertices */
	    if ((flag3(i,j  ,k  ) != 0 || flag3(i+1,j  ,k  ) != 0 || 
		 flag3(i,j+1,k  ) != 0 || flag3(i+1,j+1,k  ) != 0 || 
		 flag3(i,j  ,k+1) != 0 || flag3(i+1,j  ,k+1) != 0 || 
		 flag3(i,j+1,k+1) != 0 || flag3(i+1,j+1,k+1) != 0) && 
		(flag3(i,j  ,k  ) == 0 || flag3(i+1,j  ,k  ) == 0 || 
		 flag3(i,j+1,k  ) == 0 || flag3(i+1,j+1,k  ) == 0 || 
		 flag3(i,j  ,k+1) == 0 || flag3(i+1,j  ,k+1) == 0 || 
		 flag3(i,j+1,k+1) == 0 || flag3(i+1,j+1,k+1) == 0)){
/* draw each "alive" face */
/* face k */
	      if (flag3(i,j  ,k) != 0 && flag3(i+1,j  ,k) != 0 &&
		  flag3(i,j+1,k) != 0 && flag3(i+1,j+1,k) != 0){
/* start a new polygon */
		glBegin(GL_POLYGON);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i+1, j, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* third point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* fourth point */
		grid_face_normal(n_vec, vgrid, i, j+1, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* end the polygon */
		glEnd();
	      } /* end face k */
/* face k+1 */
	      if (flag3(i,j  ,k+1) != 0 && flag3(i+1,j  ,k+1) != 0 &&
		  flag3(i,j+1,k+1) != 0 && flag3(i+1,j+1,k+1) != 0){
/* start a new polygon */
		glBegin(GL_POLYGON);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j, k+1, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* second point */
		grid_face_normal(n_vec, vgrid, i+1, j, k+1, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* third point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k+1, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k+1), y3(i+1,j+1,k+1), z3(i+1,j+1,k+1));
/* fourth point */
		grid_face_normal(n_vec, vgrid, i, j+1, k+1, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* end the polygon */
		glEnd();
	      } /* end face k+1 */
/* face i */
	      if (flag3(i,j  ,k  ) != 0 && flag3(i,j+1,k  ) != 0 &&
		  flag3(i,j+1,k+1) != 0 && flag3(i,j  ,k+1) != 0){
/* start a new polygon */
		glBegin(GL_POLYGON);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j, k, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i, j+1, k, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* third point */
		grid_face_normal(n_vec, vgrid, i, j+1, k+1, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* fourth point */
		grid_face_normal(n_vec, vgrid, i, j, k+1, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* end the polygon */
		glEnd();
	      } /* end face i */
/* face i+1 */
	      if (flag3(i+1,j  ,k  ) != 0 && flag3(i+1,j+1,k  ) != 0 &&
		  flag3(i+1,j+1,k+1) != 0 && flag3(i+1,j  ,k+1) != 0){
/* start a new polygon */
		glBegin(GL_POLYGON);
/* first point */
		grid_face_normal(n_vec, vgrid, i+1, j, k, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* third point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k+1, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k+1), y3(i+1,j+1,k+1), z3(i+1,j+1,k+1));
/* fourth point */
		grid_face_normal(n_vec, vgrid, i+1, j, k+1, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* end the polygon */
		glEnd();
	      } /* end face i+1 */
/* face j */
	      if (flag3(i  ,j,k  ) != 0 && flag3(i  ,j,k+1) != 0 &&
		  flag3(i+1,j,k+1) != 0 && flag3(i+1,j,k  ) != 0){
/* start a new polygon */
		glBegin(GL_POLYGON);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j, k, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i, j, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* third point */
		grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* fourth point */
		grid_face_normal(n_vec, vgrid, i+1, j, k, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* end the polygon */
		glEnd();
	      } /* end face j */
/* face j+1 */
	      if (flag3(i  ,j+1,k  ) != 0 && flag3(i  ,j+1,k+1) != 0 &&
		  flag3(i+1,j+1,k+1) != 0 && flag3(i+1,j+1,k  ) != 0){
/* start a new polygon */
		glBegin(GL_POLYGON);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j+1, k, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i, j+1, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* third point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k+1), y3(i+1,j+1,k+1), z3(i+1,j+1,k+1));
/* fourth point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* end the polygon */
		glEnd();
	      } /* end face j+1 */
	    }/* end if flag != 0 */
	  } /* end for i... */
	} /* end for j... */
      } /* end for k... */

    } /* end plot_mode & 4 */

/* plot all used grid cells */
    if (plot_mode & 8){

/* set surface properties */
      ogl_set_color( vgrid->color);

/* loop over all grid cells */
      for (k=range3(1,3); k<range3(2,3); k++){
	for (j=range3(1,2); j<range3(2,2); j++){
	  for (i=range3(1,1); i<range3(2,1); i++){
/* check if all vertices in this cell are alive */
	    if ((flag3(i,j  ,k  ) != 0 && flag3(i+1,j  ,k  ) != 0 && 
		 flag3(i,j+1,k  ) != 0 && flag3(i+1,j+1,k  ) != 0 && 
		 flag3(i,j  ,k+1) != 0 && flag3(i+1,j  ,k+1) != 0 && 
		 flag3(i,j+1,k+1) != 0 && flag3(i+1,j+1,k+1) != 0)){
/* draw each face */
/* face k */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* end the polygon */
	      glEnd();
/* face k+1 */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 3, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* second point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k+1, 3, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k+1, 3, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k+1), y3(i+1,j+1,k+1), z3(i+1,j+1,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k+1, 3, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* end the polygon */
	      glEnd();
/* face i */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k+1, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* end the polygon */
	      glEnd();
/* face i+1 */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 1, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k, 1, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k+1, 1, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k+1), y3(i+1,j+1,k+1), z3(i+1,j+1,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k+1, 1, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* end the polygon */
	      glEnd();
/* face j */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* end the polygon */
	      glEnd();
/* face j+1 */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 2, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k+1, 2, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k+1, 2, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k+1), y3(i+1,j+1,k+1), z3(i+1,j+1,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k, 2, 0);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* end the polygon */
	      glEnd();
	    }/* end if flag != 0 */

	  } /* end for i... */
	} /* end for j... */
      } /* end for k... */

    } /* end plot_mode & 8 */

/* plot all used faces */
    if (plot_mode & 16){

/* set surface properties */
      ogl_set_color( vgrid->color);

/* loop over all iso-k faces */
      for (k=range3(1,3); k<=range3(2,3); k++)
	for (j=range3(1,2); j<range3(2,2); j++)
	  for (i=range3(1,1); i<range3(2,1); i++){
/* check if all vertices in the face are discretization points */
	    if (abs(flag3(i  ,j  ,k  )) == vgrid->priority && 
		abs(flag3(i+1,j  ,k  )) == vgrid->priority && 
		abs(flag3(i  ,j+1,k  )) == vgrid->priority && 
		abs(flag3(i+1,j+1,k  )) == vgrid->priority ){
/* draw each face */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* end the polygon */
	      glEnd();
	    } /* end if completely alive face */

	  } /* end for i, for j, for k */

/* loop over all iso-i faces */
      for (i=range3(1,1); i<=range3(2,1); i++)
	for (j=range3(1,2); j<range3(2,2); j++)
	  for (k=range3(1,3); k<range3(2,3); k++){
/* check if all vertices in the face are discretization points */
	    if (abs(flag3(i  ,j  ,k  )) == vgrid->priority && 
		abs(flag3(i  ,j+1,k  )) == vgrid->priority && 
		abs(flag3(i  ,j+1,k+1)) == vgrid->priority && 
		abs(flag3(i  ,j  ,k+1)) == vgrid->priority ){
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k+1, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* end the polygon */
	      glEnd();
	    } /* end if completely alive face */
	  } /* end for i,j,k */

/* loop over all iso-j faces */
      for (j=range3(1,2); j<=range3(2,2); j++)
	for (i=range3(1,1); i<range3(2,1); i++)
	  for (k=range3(1,3); k<range3(2,3); k++){
/* check if all vertices in the face are discretization points */
	    if (abs(flag3(i  ,j  ,k  )) == vgrid->priority && 
		abs(flag3(i  ,j  ,k+1)) == vgrid->priority && 
		abs(flag3(i+1,j  ,k+1)) == vgrid->priority && 
		abs(flag3(i+1,j  ,k  )) == vgrid->priority ){
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* end the polygon */
	      glEnd();
	    }/* end if completely alive face */

	  } /* end for i, j, k */

    } /* end plot_mode & 16 */

/* plot all used grid lines */
    if (plot_mode & 32){

/* set surface properties */
      ogl_set_color( vgrid->color);
      glLineWidth(3.0);
      glEnable(GL_LINE_SMOOTH);

/* loop over all grid points */
      for (k=range3(1,3); k<range3(2,3); k++)
	for (j=range3(1,2); j<range3(2,2); j++)
	  for (i=range3(1,1); i<range3(2,1); i++){
	      if (flag3(i,j,k) != 0 && flag3(i+1,j,k) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i+1, j, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
		glEnd();
	      }
	      if (flag3(i+1,j,k) != 0 && flag3(i+1,j+1,k) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i+1, j, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
		glEnd();
	      }
	      if (flag3(i+1,j+1,k) != 0 && flag3(i,j+1,k) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i, j+1, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
		glEnd();
	      }
	      if (flag3(i,j+1,k) != 0 && flag3(i,j,k) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j+1, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i, j, k, 3, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
		glEnd();
	      }

	      if (flag3(i,j+1,k) != 0 && flag3(i,j+1,k+1) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j+1, k, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* second point */
		grid_face_normal(n_vec, vgrid, i, j+1, k+1, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
		glEnd();
	      }
	      if (flag3(i,j+1,k+1) != 0 && flag3(i,j,k+1) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j+1, k+1, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* second point */
		grid_face_normal(n_vec, vgrid, i, j, k+1, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
		glEnd();
	      }
	      if (flag3(i,j,k+1) != 0 && flag3(i,j,k) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j, k+1, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* second point */
		grid_face_normal(n_vec, vgrid, i, j, k, 1, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
		glEnd();
	      }
	      if (flag3(i,j,k+1) != 0 && flag3(i+1,j,k+1) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i, j, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* second point */
		grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
		glEnd();
	      }
	      if (flag3(i+1,j,k+1) != 0 && flag3(i+1,j,k) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* second point */
		grid_face_normal(n_vec, vgrid, i+1, j, k, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
		glEnd();
	      }
	      if (flag3(i+1,j,k+1) != 0 && flag3(i+1,j+1,k+1) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* second point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k+1), y3(i+1,j+1,k+1), z3(i+1,j+1,k+1));
		glEnd();
	      }
	      if (flag3(i+1,j+1,k+1) != 0 && flag3(i+1,j+1,k) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k+1), y3(i+1,j+1,k+1), z3(i+1,j+1,k+1));
/* second point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
		glEnd();
	      }
	      if (flag3(i+1,j+1,k+1) != 0 && flag3(i,j+1,k+1) != 0){
		glBegin(GL_LINES);
/* first point */
		grid_face_normal(n_vec, vgrid, i+1, j+1, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i+1,j+1,k+1), y3(i+1,j+1,k+1), z3(i+1,j+1,k+1));
/* second point */
		grid_face_normal(n_vec, vgrid, i, j+1, k+1, 2, 1);
		glNormal3dv(n_vec);
		glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
		glEnd();
	      }

	  } /* end for i, j, k */

/* reset the line width */
      glLineWidth(1.0);
      glDisable(GL_LINE_SMOOTH);

    } /* end plot_mode & 32 */

/* interpolation faces */
    if (plot_mode & 64){
/* set surface properties */
      ogl_set_color( vgrid->color);

/* loop over all iso-k faces */
      for (k=range3(1,3); k<=range3(2,3); k++)
	for (j=range3(1,2); j<range3(2,2); j++)
	  for (i=range3(1,1); i<range3(2,1); i++){
/* check if all vertices in the face are interpolating */
	    if (abs(flag3(i  ,j  ,k  )) != vgrid->priority && 
		abs(flag3(i+1,j  ,k  )) != vgrid->priority && 
		abs(flag3(i  ,j+1,k  )) != vgrid->priority && 
		abs(flag3(i+1,j+1,k  )) != vgrid->priority && 
		flag3(i  ,j  ,k  ) != 0 && 
		flag3(i+1,j  ,k  ) != 0 && 
		flag3(i  ,j+1,k  ) != 0 && 
		flag3(i+1,j+1,k  ) != 0){
/* draw each face */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 3, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* end the polygon */
	      glEnd();
	    } /* end if completely alive face */

	  } /* end for i, for j, for k */

/* loop over all iso-i faces */
      for (i=range3(1,1); i<=range3(2,1); i++)
	for (j=range3(1,2); j<range3(2,2); j++)
	  for (k=range3(1,3); k<range3(2,3); k++){
/* check if all vertices in face are alive */
	    if (flag3(i,j  ,k  ) != 0 && flag3(i,j+1,k  ) != 0 && 
		flag3(i,j+1,k+1) != 0 && flag3(i,j  ,k+1) != 0 &&
		abs(flag3(i,j  ,k  )) != vgrid->priority && 
		abs(flag3(i,j+1,k  )) != vgrid->priority && 
		abs(flag3(i,j+1,k+1)) != vgrid->priority && 
		abs(flag3(i,j  ,k+1)) != vgrid->priority){
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k+1, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 1, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* end the polygon */
	      glEnd();
	    } /* end if completely alive face */
	  } /* end for i,j,k */

/* loop over all iso-j faces */
      for (j=range3(1,2); j<=range3(2,2); j++)
	for (i=range3(1,1); i<range3(2,1); i++)
	  for (k=range3(1,3); k<range3(2,3); k++){
/* check if all vertices in face are alive */
	    if (flag3(i  ,j,k  ) != 0 && flag3(i  ,j,k+1) != 0 && 
		flag3(i+1,j,k+1) != 0 && flag3(i+1,j,k  ) != 0 &&
		abs(flag3(i  ,j,k  )) != vgrid->priority && 
		abs(flag3(i  ,j,k+1)) != vgrid->priority && 
		abs(flag3(i+1,j,k+1)) != vgrid->priority && 
		abs(flag3(i+1,j,k  )) != vgrid->priority){
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 2, 1);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* end the polygon */
	      glEnd();
	    }/* end if completely alive face */

	  } /* end for i, j, k */

    } /* end plot+mode & 64 */

/* bad points */
    if (plot_mode & 128){
      for (bad3 = vgrid->last_bad_point; bad3 != NULL; bad3 = bad3->prev){
	ogl_marker(x3(bad3->i, bad3->j, bad3->k), 
		   y3(bad3->i, bad3->j, bad3->k), 
		   z3(bad3->i, bad3->j, bad3->k),
		   over3d->sphere_size*bb_length, bad3->type);
      }
/* reset the orientation */
      glFrontFace(GL_CW);
    } /* end plot_mode & 128 */

/* interpolation points */
    if (plot_mode & 512){
      for (interp3 = vgrid->last_interp; interp3 != NULL; 
	   interp3 = interp3->prev){
	if (interp3->active)
	  ogl_marker(x3(interp3->i_point, interp3->j_point, interp3->k_point), 
		     y3(interp3->i_point, interp3->j_point, interp3->k_point), 
		     z3(interp3->i_point, interp3->j_point, interp3->k_point),
		     over3d->sphere_size*bb_length, interp3->vgrid_loc->color);
      }
/* reset the orientation */
      glFrontFace(GL_CW);
    } /* end plot_mode & 512 */

/* dummy interpolation points */
/*      if (plot_mode & 1024){  */
/*        for (interp3 = vgrid->last_interp; interp3 != NULL; */
/* 	    interp3 = interp3->prev){ */
/* 	 if (interp3->vgrid_loc->dummy_background){ */
/* 	   i = interp3->i_point; */
/* 	   j = interp3->j_point; */
/* 	   k = interp3->k_point; */
/* 	   ogl_marker(x3(i, j, k), y3(i, j, k), z3(i, j, k),  */
/*  		   0.01*bb_length, interp3->vgrid_loc->color);  */
/*	 }*/ /* end if dumyy donor */
/*       }*/ /* end for all interpolation points */
/* reset the orientation */
/*       glFrontFace(GL_CW);*/
/*     }*/  /* end plot_mode & 1024 */ 

/* all dead points */
    if (plot_mode & 2048){
      for (k=range3(1,3); k<=range3(2,3); k++)
	for (j=range3(1,2); j<=range3(2,2); j++)
	  for (i=range3(1,1); i<=range3(2,1); i++)
	    if (flag3(i,j,k) == 0) ogl_marker(x3(i, j, k), y3(i, j, k), 
					      z3(i, j, k), 0.01*bb_length, OGL_BLACK);
/* reset the orientation */
      glFrontFace(GL_CW);
    } /* end plot_mode & 2048 */

/* draw all faces with surfacel-label > 0 */
  if (plot_mode & 4096){

    ogl_set_color( vgrid->color );

/* face i=range3(1,1) */
    if (surf3(1,1)>0){
      i=range3(1,1);
      orientation = 1;
      for (k=range3(1,3); k<range3(2,3); k++){
	for (j=range3(1,2); j<range3(2,2); j++){
	  if (flag3(i,j,k  ) != 0 && flag3(i,j+1,k  ) != 0 &&
	      flag3(i,j,k+1) != 0 && flag3(i,j+1,k+1) != 0){
/* start a new polygon */
	    glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k+1, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* end the polygon */
	    glEnd();
	  }
	} /* end for j... */
      } /* end for k ... */
    } /* end if surf3(1,1)>0 */

/* face i=range3(2,1) */
    if (surf3(2,1)>0){
      i=range3(2,1);
      orientation = 1;

      for (k=range3(1,3); k<range3(2,3); k++){
	for (j=range3(1,2); j<range3(2,2); j++){
	  if (flag3(i,j,k  ) != 0 && flag3(i,j+1,k  ) != 0 &&
	      flag3(i,j,k+1) != 0 && flag3(i,j+1,k+1) != 0){
/* start a new polygon */
	    glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k+1, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k+1), y3(i,j+1,k+1), z3(i,j+1,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 1, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* end the polygon */
	    glEnd();
	  }
	} /* end for j... */
      } /* end for k ... */
    } /* end if surf3(2,1)>0 */

/* face j=range3(1,2) */
    if (surf3(1,2)>0){
      j=range3(1,2);
      orientation = 1;

      for (k=range3(1,3); k<range3(2,3); k++){
	for (i=range3(1,1); i<range3(2,1); i++){
	  if (flag3(i  ,j,k  ) != 0 && flag3(i  ,j,k+1) != 0 &&
	      flag3(i+1,j,k+1) != 0 && flag3(i+1,j,k  ) != 0){
/* start a new polygon */
	    glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* end the polygon */
	    glEnd();
	  } 
	} /* end for i... */
      } /* end for k... */
    } /* end if surf3(1,2)>0 */

/* face j=range3(2,2) */
    if (surf3(2,2)>0){
      j=range3(2,2);
      orientation = 1;

      for (k=range3(1,3); k<range3(2,3); k++){
	for (i=range3(1,1); i<range3(2,1); i++){
	  if (flag3(i  ,j,k  ) != 0 && flag3(i  ,j,k+1) != 0 &&
	      flag3(i+1,j,k+1) != 0 && flag3(i+1,j,k  ) != 0){
/* start a new polygon */
	    glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i, j, k+1, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k+1), y3(i,j,k+1), z3(i,j,k+1));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k+1, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k+1), y3(i+1,j,k+1), z3(i+1,j,k+1));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 2, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* end the polygon */
	    glEnd();
	  } 
	} /* end for i... */
      } /* end for k... */
    } /* end if surf3(2,2)>0 */

/* face k=range3(1,3) */
    if (surf3(1,3) > 0){
      orientation = 1;

      k = range3(1,3);
      for (j=range3(1,2); j<range3(2,2); j++){
	for (i=range3(1,1); i<range3(2,1); i++){
	  if (flag3(i,j  ,k) != 0 && flag3(i+1,j  ,k) != 0 &&
	      flag3(i,j+1,k) != 0 && flag3(i+1,j+1,k) != 0){
/* start a new polygon */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* end the polygon */
	      glEnd();
	    }
	} /* end for i... */
      } /* end for j... */
    } /* end if surf3(1,3)>0 */

/* face k=range3(2,3) */
    if (surf3(2,3) > 0){
      orientation = 1;

      k = range3(2,3);
      for (j=range3(1,2); j<range3(2,2); j++){
	for (i=range3(1,1); i<range3(2,1); i++){
	  if (flag3(i,j  ,k) != 0 && flag3(i+1,j  ,k) != 0 &&
	      flag3(i,j+1,k) != 0 && flag3(i+1,j+1,k) != 0){
/* start a new polygon */
	      glBegin(GL_POLYGON);
/* first point */
	      grid_face_normal(n_vec, vgrid, i, j, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j,k), y3(i,j,k), z3(i,j,k));
/* second point */
	      grid_face_normal(n_vec, vgrid, i+1, j, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j,k), y3(i+1,j,k), z3(i+1,j,k));
/* third point */
	      grid_face_normal(n_vec, vgrid, i+1, j+1, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i+1,j+1,k), y3(i+1,j+1,k), z3(i+1,j+1,k));
/* fourth point */
	      grid_face_normal(n_vec, vgrid, i, j+1, k, 3, orientation);
	      glNormal3dv(n_vec);
	      glVertex3d(x3(i,j+1,k), y3(i,j+1,k), z3(i,j+1,k));
/* end the polygon */
	      glEnd();
	    }
	} /* end for i... */
      } /* end for j... */
    } /* end if surf3(2,3)>0 */

  } /* end plot_mode & 4096 */

    if (plot_mode & 8192) /* plot the grid surface with flag = 'w' */
      {
	for (k=range3(1,3); k<=range3(2,3); k++)
	  for (j=range3(1,2); j<=range3(2,2); j++)
	    for (i=range3(1,1); i<=range3(2,1); i++)
	      {
		if (flag3(i,j,k) == 'w')
		  ogl_marker(x3(i,j,k), y3(i,j,k), z3(i,j,k), 0.01*bb_length, 
			     OGL_BLUE); 
		else if (flag3(i,j,k) == 'h')
		  ogl_marker(x3(i,j,k), y3(i,j,k), z3(i,j,k), 0.01*bb_length, 
			     OGL_YELLOW); 
	      } /* end for all points */

      } /* end if plot_mode & 8192 */
  
  } /* end if plot_it... */
  } /* end for all grids... */

/* turn all 4 clip planes off */
  for (i=0; i<4; i++)
    glDisable(GL_CLIP_PLANE0+i);

}

