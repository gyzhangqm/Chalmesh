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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* for cos(), sin(), and sqrt() */
#include <GL/glx.h>		/* this includes X and gl.h headers */
#include <GL/glu.h>		/* gluPerspective(), gluLookAt(), GLU polygon
				 * tesselator */
#include <X11/keysym.h>		/* for XK_Escape keysym */
#include <X11/Xatom.h>   	/* for XA_RGB_DEFAULT_MAP definition */

#include "stupid_compiler.h"
#include "min_max.h"
#include "ogl_plot.h"

static struct termio io_arg;
static void
set_raw(int on);

static float colorMaps[] = {
    0.000000, 1.000000, 0.000000, 1.000000, 0.000000, 1.000000, 
    0.000000, 1.000000, 0.333333, 0.776471, 0.443137, 0.556863, 
    0.443137, 0.556863, 0.219608, 0.666667, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.039216, 0.078431, 0.117647, 0.156863, 
    0.200000, 0.239216, 0.278431, 0.317647, 0.356863, 0.400000, 
    0.439216, 0.478431, 0.517647, 0.556863, 0.600000, 0.639216, 
    0.678431, 0.717647, 0.756863, 0.800000, 0.839216, 0.878431, 
    0.917647, 0.956863, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 0.000000, 0.000000, 
    1.000000, 1.000000, 0.000000, 0.000000, 1.000000, 1.000000, 
    0.333333, 0.443137, 0.776471, 0.556863, 0.443137, 0.219608, 
    0.556863, 0.666667, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.039216, 0.078431, 0.117647, 0.156863, 0.200000, 0.239216, 
    0.278431, 0.317647, 0.356863, 0.400000, 0.439216, 0.478431, 
    0.517647, 0.556863, 0.600000, 0.639216, 0.678431, 0.717647, 
    0.756863, 0.800000, 0.839216, 0.878431, 0.917647, 0.956863, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451, 
    0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176, 
    0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000, 
    0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726, 
    0.854902, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 0.333333, 0.443137, 
    0.443137, 0.219608, 0.776471, 0.556863, 0.556863, 0.666667, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333, 
    0.666667, 0.333333, 0.666667, 0.333333, 0.039216, 0.078431, 
    0.117647, 0.156863, 0.200000, 0.239216, 0.278431, 0.317647, 
    0.356863, 0.400000, 0.439216, 0.478431, 0.517647, 0.556863, 
    0.600000, 0.639216, 0.678431, 0.717647, 0.756863, 0.800000, 
    0.839216, 0.878431, 0.917647, 0.956863, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
    0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 
    0.498039, 0.498039, 0.498039, 0.498039, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 
    0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
    1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 
};

static GLboolean ogl_is_running = GL_FALSE;
static Display *xDisplay = 0;
static Screen *xScreen;
static int xScreen_number = 0; 
static Window wRoot = 0, brand_new_window;
static long OGL_known_events;
static int OGL_n_windows, OGL_current_window;
static WINDOW_REC *OGL_window_list[MAX_WINDOWS];

static int 
ErrorHandler(Display *xDisplay, XErrorEvent *event)
{
    char buf[80];

    printf("\nReceived X error!\n");
    printf("\tError code   : %d\n", event->error_code);
    printf("\tRequest code : %d\n", event->request_code);
    printf("\tMinor code   : %d\n\n", event->minor_code);
    XGetErrorText(xDisplay, event->error_code, buf, 80);
    printf("\tError text : '%s'\n\n", buf);
    return 0;
}

void
ogl_bind_rendering(int win_number){
  WINDOW_REC *window;

  if (!ogl_is_running || win_number < 0 || win_number >= OGL_n_windows) 
    return;

  window = OGL_window_list[win_number];

/* bind the redering to the window */
  glXMakeCurrent(xDisplay, window->wMain, window->cMain);
}

static int
window_exists(Window winid, int *pos){
/* checks if the window `winid' exists and in that case, returns the location */
/* of that window in the window array */
  *pos=0;
  
  if (winid == (Window) NULL)
    return FALSE;
  if (OGL_n_windows > 0) {
    while ((*pos < OGL_n_windows) && (winid != OGL_window_list[*pos]->wMain))
      (*pos)++;
    
    return(*pos != OGL_n_windows);
  } else
    return(FALSE);
}

static GLenum 
ogl_initdisplay(char* display)
{
    int erb, evb;
 
    xDisplay = XOpenDisplay(display);
    if (!xDisplay) {
       fprintf(stderr, "Can't connect to xDisplay!\n");
       return GL_FALSE;
    }

/*** make sure OpenGL's GLX extension supported ***/
    if (!glXQueryExtension(xDisplay, &erb, &evb)) {
       fprintf(stderr, "No glx extension!\n");
       return GL_FALSE;
    }
    xScreen = DefaultScreenOfDisplay( xDisplay );
    xScreen_number = DefaultScreen( xDisplay );
    wRoot = RootWindow( xDisplay, xScreen_number );
    XSetErrorHandler( ErrorHandler );

    OGL_known_events = StructureNotifyMask | ExposureMask | KeyPressMask |
      ButtonPressMask | Button1MotionMask | Button2MotionMask | Button3MotionMask;

    return(GL_TRUE);
}

static XVisualInfo *
FindMainVisual(GLenum type)
{
    int list[32], i;

    i = 0;

    list[i++] = GLX_LEVEL;
    list[i++] = 0;

    if (TK_IS_DOUBLE(type)) {
	list[i++] = GLX_DOUBLEBUFFER;
    }

    if (TK_IS_RGB(type)) {
	list[i++] = GLX_RGBA;
	list[i++] = GLX_RED_SIZE;
	list[i++] = 1;
	list[i++] = GLX_GREEN_SIZE;
	list[i++] = 1;
	list[i++] = GLX_BLUE_SIZE;
	list[i++] = 1;
	if (TK_HAS_ALPHA(type)) {
	    list[i++] = GLX_ALPHA_SIZE;
	    list[i++] = 1;
	}
	if (TK_HAS_ACCUM(type)) {
	    list[i++] = GLX_ACCUM_RED_SIZE;
	    list[i++] = 1;
	    list[i++] = GLX_ACCUM_GREEN_SIZE;
	    list[i++] = 1;
	    list[i++] = GLX_ACCUM_BLUE_SIZE;
	    list[i++] = 1;
	    if (TK_HAS_ALPHA(type)) {
		list[i++] = GLX_ACCUM_ALPHA_SIZE;
		list[i++] = 1;
	    }
	}
    } else if (TK_IS_INDEX(type)) {
	list[i++] = GLX_BUFFER_SIZE;
	list[i++] = 1;
    }

    if (TK_HAS_DEPTH(type)) {
	list[i++] = GLX_DEPTH_SIZE;
	list[i++] = 1;
    }

    if (TK_HAS_STENCIL(type)) {
	list[i++] = GLX_STENCIL_SIZE;
	list[i++] = 1;
    }

    list[i] = (int)None;

    return glXChooseVisual(xDisplay, xScreen_number, list);
}

static XVisualInfo *
FindBestMainVisual(GLenum type)
{
    int list[32], i;

    i = 0;

    list[i++] = GLX_LEVEL;
    list[i++] = 0;

    if (TK_IS_DOUBLE(type)) {
	list[i++] = GLX_DOUBLEBUFFER;
    }

    if (TK_IS_RGB(type)) {
	list[i++] = GLX_RGBA;
	list[i++] = GLX_RED_SIZE;
	list[i++] = 1;
	list[i++] = GLX_GREEN_SIZE;
	list[i++] = 1;
	list[i++] = GLX_BLUE_SIZE;
	list[i++] = 1;
	if (TK_HAS_ALPHA(type)) {
	    list[i++] = GLX_ALPHA_SIZE;
	    list[i++] = 1;
	}
	if (TK_HAS_ACCUM(type)) {
	    list[i++] = GLX_ACCUM_RED_SIZE;
	    list[i++] = 1;
	    list[i++] = GLX_ACCUM_GREEN_SIZE;
	    list[i++] = 1;
	    list[i++] = GLX_ACCUM_BLUE_SIZE;
	    list[i++] = 1;
	    if (TK_HAS_ALPHA(type)) {
		list[i++] = GLX_ACCUM_ALPHA_SIZE;
		list[i++] = 1;
	    }
	}
    } else if (TK_IS_INDEX(type)) {
	list[i++] = GLX_BUFFER_SIZE;
	list[i++] = 1;
    }

    if (TK_HAS_DEPTH(type)) {
	list[i++] = GLX_DEPTH_SIZE;
	list[i++] = 1;
    }

    if (TK_HAS_STENCIL(type)) {
	list[i++] = GLX_STENCIL_SIZE;
	list[i++] = 1;
    }

    list[i] = (int)None;

    return glXChooseVisual(xDisplay, xScreen_number, list);
}


static XVisualInfo *
FindExactMainVisual(GLenum type)
{
    int i, nvis, val, rval, gval, bval, aval;
    XVisualInfo *vis_list, *this_vis, *best_vis, sampleVis;
    int this_score, best_score;

    /* Get list of visuals for this screen */
    sampleVis.screen = xScreen_number;
    vis_list = XGetVisualInfo( xDisplay, VisualScreenMask, &sampleVis, &nvis);

    /* 
     * Loop through the visuals; find first one that matches the attr 
     * specified in type
     */
    best_score = -1; best_vis = NULL;
    for ( i = 0; i < nvis; i++ ) {
      this_vis = &vis_list[i];

      /* Visual must be supported by GLX */
      if ( glXGetConfig(xDisplay, this_vis, GLX_USE_GL, &val) ) continue;
      if ( !val ) continue;

      /* Visual must be in main planes which is level 0 */
      glXGetConfig(xDisplay, this_vis, GLX_LEVEL, &val);
      if ( val != 0 ) continue;

      /* Color Index or RGBA? It must match the requested value */
      glXGetConfig(xDisplay, this_vis, GLX_RGBA, &val);
      if ( TK_IS_RGB(type) && !val ) continue;
      if ( TK_IS_INDEX(type) && val ) continue;

      /* Double buffered or Single buffered? */
      glXGetConfig( xDisplay, this_vis, GLX_DOUBLEBUFFER, &val);
      if ( TK_IS_DOUBLE(type) && !val ) continue;
      if ( TK_IS_SINGLE(type) && val ) continue;

      /* If accum requested then accum rgb size must be > 0 */
      /* If alpha requested then alpha size must be > 0 */
      /* if accum & alpha requested then accum alpha size must be > 0 */
      if ( TK_IS_RGB(type) ) {
        glXGetConfig(xDisplay, this_vis, GLX_ACCUM_RED_SIZE, &rval);
        glXGetConfig(xDisplay, this_vis, GLX_ACCUM_GREEN_SIZE, &gval);
        glXGetConfig(xDisplay, this_vis, GLX_ACCUM_BLUE_SIZE, &bval);
        glXGetConfig(xDisplay, this_vis, GLX_ACCUM_ALPHA_SIZE, &aval);
        if ( TK_HAS_ACCUM(type) ) {
            if ( rval <= 0 || gval <= 0 || bval <= 0 ) continue;
        } else {
            if ( rval > 0 || gval > 0 || bval > 0 || aval > 0 ) continue;
        }

        glXGetConfig(xDisplay, this_vis, GLX_ALPHA_SIZE, &val);
        if ( TK_HAS_ALPHA(type) ) {
            if ( val <= 0 ) continue;
            if ( TK_HAS_ACCUM(type) && aval <= 0 ) continue;
        } else {
            if ( val > 0 ) continue;
        }

      }

      /* Check depth buffer */
      glXGetConfig(xDisplay, this_vis, GLX_DEPTH_SIZE, &val);
      if ( TK_HAS_DEPTH(type) ) {
            if ( val <= 0 ) continue;
      } else {
            if ( val > 0 ) continue;
      }

      /* Check stencil buffer */
      glXGetConfig( xDisplay, this_vis, GLX_STENCIL_SIZE, &val);
      if ( TK_HAS_STENCIL(type) ) {
            if ( val <= 0 ) continue;
      } else {
            if ( val > 0 ) continue;
      }

      glXGetConfig(xDisplay, this_vis, GLX_BUFFER_SIZE, &this_score);

      if (this_score > best_score ) {
          best_score = this_score;
          best_vis = this_vis;
      }

    }

    if ( best_vis ) {
        sampleVis.visualid = best_vis->visualid;
        sampleVis.screen = xScreen_number;
        if ( nvis > 0 ) XFree(vis_list);
        return XGetVisualInfo(xDisplay, VisualIDMask|VisualScreenMask, 
              &sampleVis, &nvis);
    } else {
        if ( nvis > 0 ) XFree(vis_list);
        return None;
    }

}

static GLenum 
GetMainWindowType(XVisualInfo *vi, GLXContext context)
{
    GLenum mask;
    int x, y, z;

    mask = 0;

    glXGetConfig(xDisplay, vi, GLX_DOUBLEBUFFER, &x);
    if (x) {
	mask |= TK_DOUBLE;
    } else {
	mask |= TK_SINGLE;
    }

    glXGetConfig(xDisplay, vi, GLX_RGBA, &x);
    if (x) {
	mask |= TK_RGB;
	glXGetConfig(xDisplay, vi, GLX_ALPHA_SIZE, &x);
	if (x > 0) {
	    mask |= TK_ALPHA;
	}
	glXGetConfig(xDisplay, vi, GLX_ACCUM_RED_SIZE, &x);
	glXGetConfig(xDisplay, vi, GLX_ACCUM_GREEN_SIZE, &y);
	glXGetConfig(xDisplay, vi, GLX_ACCUM_BLUE_SIZE, &z);
	if (x > 0 && y > 0 && z > 0) {
	    mask |= TK_ACCUM;
	}
    } else {
	mask |= TK_INDEX;
    }

    glXGetConfig(xDisplay, vi, GLX_DEPTH_SIZE, &x);
    if (x > 0) {
	mask |= TK_DEPTH;
    }

    glXGetConfig(xDisplay, vi, GLX_STENCIL_SIZE, &x);
    if (x > 0) {
	mask |= TK_STENCIL;
    }

    if (glXIsDirect(xDisplay, context)) {
	mask |= TK_DIRECT;
    } else {
	mask |= TK_INDIRECT;
    }

    return mask;
}

static int 
WaitForMainWindow(Display *d, XEvent *e, char *arg)
{
  if (e->type == MapNotify && e->xmap.window == brand_new_window) {
    return GL_TRUE;
  } else {
    return GL_FALSE;
  }
}


static void 
tkSetRGBMap(int size, float *rgb, WINDOW_REC *window)
{
  XColor c;
  int rShift, gShift, bShift, max, i;
  
  switch (window->vInfoMain->class) {

  case DirectColor:
    max = (size > window->vInfoMain->colormap_size) ? 
      window->vInfoMain->colormap_size : size;
    for (i = 0; i < max; i++) {
      rShift = ffs((unsigned int)window->vInfoMain->red_mask) - 1;
      gShift = ffs((unsigned int)window->vInfoMain->green_mask) - 1;
      bShift = ffs((unsigned int)window->vInfoMain->blue_mask) - 1;
      c.pixel = ((i << rShift) & window->vInfoMain->red_mask) |
	((i << gShift) & window->vInfoMain->green_mask) |
	  ((i << bShift) & window->vInfoMain->blue_mask);
      c.red = (unsigned short)(rgb[i] * 65535.0 + 0.5);
      c.green = (unsigned short)(rgb[size+i] * 65535.0 + 0.5);
      c.blue = (unsigned short)(rgb[size*2+i] * 65535.0 + 0.5);
      c.flags = DoRed | DoGreen | DoBlue;
      XStoreColor(xDisplay, window->cMapMain, &c);
    }
    break;

  case GrayScale:
  case PseudoColor:
    max = (size > window->vInfoMain->colormap_size) ? 
      window->vInfoMain->colormap_size : size;
    for (i = 0; i < max; i++) {
      c.pixel = i;
      c.red = (unsigned short)(rgb[i] * 65535.0 + 0.5);
      c.green = (unsigned short)(rgb[size+i] * 65535.0 + 0.5);
      c.blue = (unsigned short)(rgb[size*2+i] * 65535.0 + 0.5);
      c.flags = DoRed | DoGreen | DoBlue;
      XStoreColor(xDisplay, window->cMapMain, &c);
    }
    break;
  }

  XSync(xDisplay, 0);
}


static WINDOW_REC *
ogl_initwindow(char *title, double x, double y, double w, double h)
{
  WINDOW_REC *window;
  XSetWindowAttributes winAttrib;
  XWindowAttributes attrib;
  XSizeHints size_hints;
  Atom deleteWindowAtom;
  XEvent e;
  double aspect_ratio;
  int n_cmap, val, status;
  XVisualInfo *this_vis;

/* make sure the display is setup */
    if (!xDisplay) 
      return NULL;

  window = (WINDOW_REC *) malloc( sizeof(WINDOW_REC) );

/* set the size */
  window->x = size_hints.x = WidthOfScreen(xScreen) * x;
  window->y = size_hints.y = HeightOfScreen(xScreen) * (1.0-y);
  window->w = size_hints.width  = WidthOfScreen(xScreen) * w;
  window->h = size_hints.height = HeightOfScreen(xScreen) * h;

/* default length scale */
  window->scale = 1.0;

/* get aspect ratio of screen */
  XGetWindowAttributes(xDisplay, wRoot, &attrib);
  aspect_ratio = ((double) attrib.width) / ((double) attrib.height);

/* minimum size */
  size_hints.min_width = WidthOfScreen(xScreen) / 4;
  size_hints.min_height = size_hints.min_width / aspect_ratio;
  size_hints.flags = USPosition|USSize|PMinSize;

/* The XInternAtom function returns the atom identifier associated with the */
/* specified atom_name string. */
  deleteWindowAtom = XInternAtom(xDisplay, "WM_DELETE_WINDOW", False);

/* set the type */
  window->type = TK_RGB|TK_DOUBLE|TK_DIRECT|TK_DEPTH; 
  window->dmPolicy = (int) TK_MINIMUM_CRITERIA;

/*** (4) find an appropriate visual and a colormap for it ***/
  if (window->dmPolicy == TK_MINIMUM_CRITERIA)
/*    window->vInfoMain = FindBestMainVisual(window->type);*/
    window->vInfoMain = FindMainVisual(window->type); 
  else if (window->dmPolicy == TK_EXACT_MATCH)
    window->vInfoMain = FindExactMainVisual(window->type);
  if (!window->vInfoMain) {
    fprintf(stderr, "Window type not found!\n");
    free(window);
    return NULL;
  }

/*** (5) create an OpenGL rendering context  ***/ 
  window->cMain = glXCreateContext(xDisplay, window->vInfoMain, None,
				  (TK_IS_DIRECT(window->type))?GL_TRUE:GL_FALSE);
  if (!window->cMain) {
    fprintf(stderr, "Can't create a context!\n");
    free(window);
    return NULL;
  }

  window->type = GetMainWindowType(window->vInfoMain, window->cMain);

/* color stuff */
  if (TK_IS_INDEX(window->type)) {
    if (window->vInfoMain->class != StaticColor &&
	window->vInfoMain->class != StaticGray) {

      window->cMapMain = XCreateColormap(xDisplay, wRoot, 
					 window->vInfoMain->visual,
					 AllocAll);
/* fill in the color entries */
      tkSetRGBMap(256, colorMaps, window);
    } 
    else {
      window->cMapMain = XCreateColormap(xDisplay, wRoot, 
					window->vInfoMain->visual,
					AllocNone);
/* fill in the color entries */
      tkSetRGBMap(256, colorMaps, window);
    }
  }
  else {
/* Mesa might return a PseudoColor visual for RGB mode. */
    if ((n_cmap=MaxCmapsOfScreen(DefaultScreenOfDisplay(xDisplay))) == 1
	&& window->vInfoMain->visual == DefaultVisual(xDisplay, xScreen_number)) {
/* display only has one colormap, let's share it to
   prevent flashing */
      window->cMapMain = DefaultColormap(xDisplay, xScreen_number);
    } else {
/* Get our own PseudoColor colormap. */
      window->cMapMain = XCreateColormap(xDisplay, wRoot,
					 window->vInfoMain->visual, AllocNone);
/* fill in the color entries */
      tkSetRGBMap(256, colorMaps, window);
    }
  }

/* copy info onto the window attributes record */
  winAttrib.colormap = window->cMapMain;
  winAttrib.background_pixmap = None;
  winAttrib.border_pixel = 0;
  winAttrib.event_mask = OGL_known_events;

/*** (6) create an X window with selected visual and right properties ***/
/* We need to store the window pointer in brand_new_window in order for */
/* WaitForMainWindow to recognize it */
  brand_new_window = window->wMain = 
    XCreateWindow(xDisplay, wRoot, window->x, window->y, 
		  window->w, window->h, 0,
		  window->vInfoMain->depth, InputOutput,
		  window->vInfoMain->visual,
		  CWBackPixmap|CWBorderPixel|CWEventMask|CWColormap,
		  &winAttrib);

/* size hints */
  XSetWMNormalHints(xDisplay, brand_new_window, &size_hints);
	
/* title and icon strings */
  XStoreName(xDisplay, brand_new_window, title);

/*   XStringListToTextProperty(&title, 1, &textProp); */
/*   XSetWMProperties(xDisplay, window->wMain, &textProp, &textProp, NULL, 0,  */
/* 		   NULL, NULL, NULL); */

/* disable the "close" option in the window menu */
  XSetWMProtocols(xDisplay, window->wMain, &deleteWindowAtom, 1); 

/*** (10) request the X window to be displayed on the screen ***/
  XMapWindow(xDisplay, window->wMain);
  XIfEvent(xDisplay, &e, WaitForMainWindow, 0);


/* bind the redering to the window */
  glXMakeCurrent(xDisplay, window->wMain, window->cMain);

/* to get a font it is not enough to call xload queryfont once. It must be called */
/* for each window */
  if ((window->xFont = XLoadQueryFont(xDisplay, "9x15")) == NULL){
    printf("Error in ogl_initdisplay: Could not find a font!\n");
    return GL_FALSE;
  }
/* build a bitmap font from a X-font */
  window->bitmapBase = glGenLists(256);
  glXUseXFont(window->xFont->fid, 0, 255, window->bitmapBase);

/* set default values in the remaining fields */
  window->n_active_lists = 0;
  window->need_redraw = GL_TRUE;
  window->view_state.vdir[0] = 1.0;
  window->view_state.vdir[1] = 0.0;
  window->view_state.vdir[2] = 0.0;
  window->view_state.up[0] = 0.0;
  window->view_state.up[1] = 0.0;
  window->view_state.up[2] = 1.0;
  window->view_state.right[0] = 0.0;
  window->view_state.right[1] = 1.0;
  window->view_state.right[2] = 0.0;
  window->view_state.focus[0] = 0.0;
  window->view_state.focus[1] = 0.0;
  window->view_state.focus[2] = 0.0;
  window->view_state.distance = 10.0;

  window->initial_view_state = window->view_state;

  window->use_scaling = 0;
  window->x_scale = window->y_scale = window->z_scale = 1.0;

  return window;
}

#if defined(DEC) && defined(vms)

static int 
ffs (unsigned int mask)
{
    int			num;
    unsigned int	bit;

    if (mask == 0) return 0;
    
    for (num = 1, bit = 1; (mask & bit) != 0; num++, bit <<= 1) {}

    return num;

}

#endif /* DEC */

/******************************************************************************/

/**************** START ogl routines ***********************************/

void
ogl_init_prompt(input_output *io_ptr){
/* Initially read the commands from standard in and don't save the commands */
  io_ptr->read_command = stdin;
  io_ptr->copy_command = NULL;
}

void
ogl_close_prompt(input_output *io_ptr){
/* close open files */
  if (io_ptr->read_command != stdin){
    fclose(io_ptr->read_command);
    io_ptr->read_command = stdin;
  }

  if (io_ptr->copy_command != NULL){
/* remind the user of the file name */
    printf("Closing the log-file `%s'\n", io_ptr->copy_file_name);
    fclose(io_ptr->copy_command);
    io_ptr->copy_command = NULL;
/* free the file name */
    free( io_ptr->copy_file_name );
  }
}

int 
ogl_init(char *display, char **title, double *x, double *y, double *w, double *h, 
	 int n_windows){
  int i;

/* no need to do anything if ogl already is initialized */
  if (ogl_is_running) 
    return ERROR;

  if (n_windows<1 || n_windows>MAX_WINDOWS){
    printf("Fatal error: Too few or too many windows requested: %i\n", n_windows);
    return ERROR;
  }

/*** (2) open a connection to the X server and make sure that OpenGL is supported ***/
  if (!ogl_initdisplay(display)) 
    return ERROR;

/* get default settings for io */
  if (ioctl (fileno(stdin), TCGETA, &io_arg)==-1)
    perror("ioctl");

/* open all windows */
  for (i=0; i<n_windows; i++){
/*    OGL_window_list[i] = (WINDOW_REC *) malloc( sizeof(WINDOW_REC) );*/
    if ( (OGL_window_list[i] = ogl_initwindow( title[i], x[i], y[i], w[i], h[i] )) 
	== NULL ) 
      return ERROR;
  }

/* flush the output buffer. Maybe unnecessary, but it won't hurt */
  XFlush(xDisplay);

  OGL_n_windows = n_windows;
  ogl_is_running = GL_TRUE;

  return OK;
}

void 
ogl_close( void ){
  int i;
  WINDOW_REC *window;
  
/* check that ogl is running */
  if (!ogl_is_running) 
    return;

/* delete all display lists with fonts */
  for(i=0; i<OGL_n_windows; i++){
    window = OGL_window_list[i];
    glDeleteLists(window->bitmapBase, 256);
  }

/* release the current context */
  glXMakeCurrent(xDisplay, None, NULL);

/* destroy all windows */
  for (i=0; i<OGL_n_windows; i++){
    XDestroyWindow(xDisplay, OGL_window_list[i]->wMain);
  }

/* destroy the all contexts */
  for (i=0; i<OGL_n_windows; i++){
    glXDestroyContext( xDisplay, OGL_window_list[i]->cMain );
  }

/* flush the output buffer. Maybe unnecessary, but it won't hurt */
  XFlush(xDisplay);

/* Shut down X server */
  XCloseDisplay(xDisplay);

/* reset the value of some global variables to prevent the routines to do anything */
/* strange */
  xDisplay = NULL;
  OGL_n_windows = 0;
  ogl_is_running = GL_FALSE;

/* reset default settings for io */
  if (ioctl (fileno(stdin), TCSETA, &io_arg)==-1)
    perror("ioctl");

}


static int
ogl_set_modelview(int win_number,
		  GLdouble eye_x, GLdouble eye_y, GLdouble eye_z,
		  GLdouble focus_x, GLdouble focus_y, GLdouble focus_z,
		  GLdouble up_x, GLdouble up_y, GLdouble up_z)
{
  int i;
  GLdouble alpha;
  WINDOW_REC *window;

  if (!ogl_is_running || win_number < 0 || win_number >= OGL_n_windows) 
    return ERROR;

  window = OGL_window_list[ win_number ];
/* make the context current rendering */
  glXMakeCurrent(xDisplay, window->wMain, window->cMain);

/* fill in the elements in the ogl_view_state structure */
  window->view_state.focus[0] = focus_x;
  window->view_state.focus[1] = focus_y;
  window->view_state.focus[2] = focus_z;

  window->view_state.vdir[0] = focus_x - eye_x;
  window->view_state.vdir[1] = focus_y - eye_y;
  window->view_state.vdir[2] = focus_z - eye_z;

  window->view_state.distance = 
    sqrt(window->view_state.vdir[0] * 
	 window->view_state.vdir[0] +
	 window->view_state.vdir[1] * 
	 window->view_state.vdir[1] +
	 window->view_state.vdir[2] * 
	 window->view_state.vdir[2]);
/* normalize vdir */
  for (i=0; i<3; i++) window->view_state.vdir[i] /= 
    window->view_state.distance;

  window->view_state.up[0] = up_x; 
  window->view_state.up[1] = up_y; 
  window->view_state.up[2] = up_z;

/* make the upvector orthogonal to vdir */
  alpha = 0;
  for (i=0; i<3; i++) alpha += window->view_state.vdir[i] * 
    window->view_state.up[i];
  for (i=0; i<3; i++) window->view_state.up[i] -= 
    alpha*window->view_state.vdir[i];
    
/* normalize up */
  alpha = sqrt(window->view_state.up[0] * 
	       window->view_state.up[0] +
	       window->view_state.up[1] * 
	       window->view_state.up[1] +
	       window->view_state.up[2] * 
	       window->view_state.up[2]);
  for (i=0; i<3; i++) window->view_state.up[i] /= alpha;

/* compute the right vector */
  window->view_state.right[0] = 
    window->view_state.vdir[1] * 
      window->view_state.up[2] -
	window->view_state.vdir[2] *
	  window->view_state.up[1];
  window->view_state.right[1] = 
    window->view_state.vdir[2] * 
      window->view_state.up[0] -
	window->view_state.vdir[0] * 
	  window->view_state.up[2];
  window->view_state.right[2] = 
    window->view_state.vdir[0] * 
      window->view_state.up[1] -
	window->view_state.vdir[1] * 
	  window->view_state.up[0];

  glMatrixMode(GL_MODELVIEW);	/* now change the modelview matrix */
  glLoadIdentity();
  gluLookAt(window->view_state.focus[0] - 
	    window->view_state.distance * 
	    window->view_state.vdir[0], 
	    window->view_state.focus[1] - 
	    window->view_state.distance * 
	    window->view_state.vdir[1], 
	    window->view_state.focus[2] - 
	    window->view_state.distance * 
	    window->view_state.vdir[2], 
	    window->view_state.focus[0], 
	    window->view_state.focus[1], 
	    window->view_state.focus[2], 
	    window->view_state.up[0], 
	    window->view_state.up[1], 
	    window->view_state.up[2]);

/* scale for normalized or unnormalized coordinates */
  glScaled(window->x_scale, window->y_scale, window->z_scale);

/* save the view state */
  window->initial_view_state = window->view_state;

  return OK;
}

void
ogl_standard_view( int window_number, bounding_box *bb ){
  double length;
  WINDOW_REC *window;

/* check the window number */
  if (!ogl_is_running || window_number < 0 || window_number >= OGL_n_windows) return;

  window = OGL_window_list[ window_number ];

/* normalized coordinates? */
  if (window->use_scaling)
    {
/* scale-factors */
      window->x_scale = 1./(bb->x_max - bb->x_min);
      window->y_scale = 1./(bb->y_max - bb->y_min);
      window->z_scale = 1./(bb->z_max - bb->z_min);

      length = 3.0;
      ogl_set_modelview(window_number,
			0.5 - length, 
			0.5 - length, 
			0.5 + length, /* eye point */
			0.5, 0.5, 0.5, /* focal point */
			0.0, 0.0, 1.0);/* up vector */
    }
  else
    {
/* scale-factors */
      window->x_scale = 1.;
      window->y_scale = 1.;
      window->z_scale = 1.;

/* set the standard view point */ 
      length = real_max( 1.e-3*window->scale, fabs(bb->x_max - bb->x_min) +
			fabs(bb->y_max - bb->y_min) +
			fabs(bb->z_max - bb->z_min));
      ogl_set_modelview(window_number,
			0.5*(bb->x_max + bb->x_min) - length, 
			0.5*(bb->y_max + bb->y_min) - length, 
			0.5*(bb->z_max + bb->z_min) + length, /* eye point */
			0.5*(bb->x_max + bb->x_min), 
			0.5*(bb->y_max + bb->y_min), 
			0.5*(bb->z_max + bb->z_min), /* focal point */
			0.0, 0.0, 1.0);/* up vector */
    }

}

void
ogl_normalize( int window_number, bounding_box *bb ){
  WINDOW_REC *window;

/* check the window number */
  if (!ogl_is_running || window_number < 0 || window_number >= OGL_n_windows) return;

  window = OGL_window_list[ window_number ];

  window->use_scaling = 1;

  ogl_standard_view( window_number, bb);

  glEnable( GL_NORMALIZE );
}

void
ogl_unnormalize( int window_number, bounding_box *bb ){
  WINDOW_REC *window;

/* check the window number */
  if (!ogl_is_running || window_number < 0 || window_number >= OGL_n_windows) return;

  window = OGL_window_list[ window_number ];

  window->use_scaling = 0;

  ogl_standard_view( window_number, bb);

  glDisable( GL_NORMALIZE );
}

int
ogl_reset_modelview(int win_number)
{
  WINDOW_REC *window;

  if (!ogl_is_running || win_number < 0 || win_number >= OGL_n_windows) 
    return ERROR;

  window = OGL_window_list[ win_number ];

/* bind the redering to the window */
  glXMakeCurrent(xDisplay, window->wMain, window->cMain);

  window->view_state = window->initial_view_state;
  glMatrixMode(GL_MODELVIEW);	/* now change the modelview matrix */
  glLoadIdentity();
  gluLookAt(window->view_state.focus[0] - 
	    window->view_state.distance * 
	    window->view_state.vdir[0], 
	    window->view_state.focus[1] - 
	    window->view_state.distance * 
	    window->view_state.vdir[1], 
	    window->view_state.focus[2] - 
	    window->view_state.distance * 
	    window->view_state.vdir[2], 
	    window->view_state.focus[0], 
	    window->view_state.focus[1], 
	    window->view_state.focus[2], 
	    window->view_state.up[0], 
	    window->view_state.up[1], 
	    window->view_state.up[2]);

  glScaled(window->x_scale, window->y_scale, window->z_scale);

  window->need_redraw = GL_TRUE;

  return OK;
}

static void
ogl_change_modelview(WINDOW_REC *window,
		     GLdouble x_trans, GLdouble y_trans, GLdouble dist_trans, 
		     GLdouble approach, GLdouble x_rot, GLdouble y_rot, GLdouble z_rot){
  int i;
  GLdouble vtmp[3], rtmp[3], utmp[3];
  double min_distance;

  if (!ogl_is_running) 
    return;

/* bind the redering to the window */
  glXMakeCurrent(xDisplay, window->wMain, window->cMain);

/* change the distance between the eye and the focal point */
  window->view_state.distance += approach; 

  min_distance = window->scale*1.e-2;
/* make sure the distance stays positive */
  window->view_state.distance =  
    (window->view_state.distance > min_distance)?  
      window->view_state.distance : min_distance;    

/* move the focal point away from the eye */
  for (i=0; i<3; i++) window->view_state.focus[i] -=
    dist_trans * window->view_state.vdir[i];

/* horizontal translation */
  for (i=0; i<3; i++) window->view_state.focus[i] += 
    x_trans*window->view_state.right[i];
/* vertical translation */
  for (i=0; i<3; i++) window->view_state.focus[i] += 
    y_trans*window->view_state.up[i];

/* make sure the rotation angle x_rot is reasonable */
  if (x_rot < -0.9) 
    x_rot = -0.9;
  else if (x_rot > 0.9) 
    x_rot = 0.9;
/* rotation around the up vector */
  for (i=0; i<3; i++){
    vtmp[i] = 
      sqrt(1.0 - x_rot*x_rot) * window->view_state.vdir[i] - 
	x_rot * window->view_state.right[i];
    rtmp[i] = 
      x_rot * window->view_state.vdir[i] + 
	sqrt(1.0 - x_rot*x_rot) * window->view_state.right[i];
  }
  for (i=0; i<3; i++){
    window->view_state.vdir[i] = vtmp[i];
    window->view_state.right[i] = rtmp[i];
  }
  
/* make sure the rotation angle y_rot is reasonable */
  if (y_rot < -0.9) 
    y_rot = -0.9;
  else if (y_rot > 0.9) 
    y_rot = 0.9;
/* rotation around the right vector */
  for (i=0; i<3; i++){
    vtmp[i] = 
      sqrt(1.0 - y_rot*y_rot) * window->view_state.vdir[i] - 
	y_rot * window->view_state.up[i];
    utmp[i] = 
      y_rot * window->view_state.vdir[i] + 
	sqrt(1.0 - y_rot*y_rot) * window->view_state.up[i];
  }
  for (i=0; i<3; i++){
    window->view_state.vdir[i] = vtmp[i];
    window->view_state.up[i]   = utmp[i];
  }
  
/* make sure the rotation angle z_rot is reasonable */
  if (z_rot < -0.9) 
    z_rot = -0.9;
  else if (z_rot > 0.9) 
    z_rot = 0.9;
/* rotation around the vdir vector */
  for (i=0; i<3; i++){
    rtmp[i] = 
      sqrt(1.0 - z_rot*z_rot) * window->view_state.right[i] - 
	z_rot * window->view_state.up[i];
    utmp[i] = 
      z_rot * window->view_state.right[i] + 
	sqrt(1.0 - z_rot*z_rot) * window->view_state.up[i];
  }
  for (i=0; i<3; i++){
    window->view_state.right[i] = rtmp[i];
    window->view_state.up[i]    = utmp[i];
  }
  
/* change the modelview matrix */
  glMatrixMode(GL_MODELVIEW);	
  glLoadIdentity();
  gluLookAt(window->view_state.focus[0] - 
	    window->view_state.distance * window->view_state.vdir[0], 
	    window->view_state.focus[1] - 
	    window->view_state.distance * window->view_state.vdir[1], 
	    window->view_state.focus[2] - 
	    window->view_state.distance * window->view_state.vdir[2], 
	    window->view_state.focus[0], 
	    window->view_state.focus[1], 
	    window->view_state.focus[2], 
	    window->view_state.up[0], 
	    window->view_state.up[1], 
	    window->view_state.up[2]);

  glScaled(window->x_scale, window->y_scale, window->z_scale);

}

static void
ogl_redraw(WINDOW_REC *window)
{
  int i;

  if (!ogl_is_running) 
    return;

/* bind the redering to the window */
  glXMakeCurrent(xDisplay, window->wMain, window->cMain);

  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

/* replay all stored display lists */
  for (i=0; i<window->n_active_lists; i++)
    glCallList(window->current_display_list[i]);

  if (TK_IS_DOUBLE(window->type))
/* buffer swap does implicit glFlush */
    glXSwapBuffers(xDisplay, window->wMain);	
  else 
/* explicit flush for single buffered case */
    glFlush();
/* reset the redraw flag */
  window->need_redraw = GL_FALSE;
}

int
ogl_poll(void){
  XEvent event;
  XAnyEvent *any_event;
  KeySym ks;
  int win_number, dummy, i;
  WINDOW_REC *window=NULL;
  GLboolean recalcModelView = GL_FALSE;
  GLdouble aspect;
  GLfloat x_rot=0.0, y_rot=0.0, z_rot=0.0, x_trans=0.0, y_trans=0.0, dist_change=0.0
    , approach    = 0.0;
  unsigned int button_mask;
  static int lastX=0, lastY=0;

  if (!ogl_is_running) 
    return ERROR;

/* no-blocking event loop */
  while (XCheckMaskEvent(xDisplay, OGL_known_events, &event)){

/* get the window */
    any_event = (XAnyEvent *) &event;
    if (!window_exists( any_event->window, &win_number ))
      return ERROR;

    window = OGL_window_list[win_number];
/* bind the redering to the window */
    glXMakeCurrent(xDisplay, window->wMain, window->cMain);
    
    switch (event.type){

    case ConfigureNotify:
/* change the viewport */
      window->w = event.xconfigure.width;
      window->h = event.xconfigure.height;
      glViewport(0, 0, window->w, window->h);
      glMatrixMode(GL_PROJECTION);/* set up projection transform */
      glLoadIdentity();
      gluPerspective(40.0, /* field of view in degree */ 
		     ((double)window->w)/((double)window->h), /* aspect ratio */
		     1.e-2 * window->scale, /* Z near */  
		     100.0 * window->scale /* Z far */); 
/* the ratio 100/0.01 determines the z-buffer resolution */
      glMatrixMode(GL_MODELVIEW);/* modify the modelview matrix */
      window->need_redraw = GL_TRUE;

/* end ConfigureNotify */
      break;
  
    case Expose:
      window->need_redraw = GL_TRUE;
/* end Expose */
      break;

    case KeyPress:
      ks = XLookupKeysym((XKeyEvent *) & event, 0);
      recalcModelView = GL_FALSE;

/* only change the modelview if there is something to look at */
      if (window->n_active_lists > 0){
	if (ks == XK_Left){
	  recalcModelView = GL_TRUE;
	  x_trans     = 0.0;
	  y_trans     = 0.0;
	  dist_change = 0.0;
	  approach    = 0.0;
	  x_rot       = 0.0;
	  y_rot       = 0.0;
	  z_rot       = 0.0;
	  if (event.xkey.state & ShiftMask)
	    x_rot   = 0.1;
	  else if (event.xkey.state & ControlMask)
	    z_rot   = 0.1;
	  else
	    x_trans = 0.02 * window->view_state.distance;
	}
	else if (ks == XK_Right){
	  recalcModelView = GL_TRUE;
	  x_trans     = 0.0;
	  y_trans     = 0.0;
	  dist_change = 0.0;
	  approach    = 0.0;
	  x_rot       = 0.0;
	  y_rot       = 0.0;
	  z_rot       = 0.0;
	  if (event.xkey.state & ShiftMask)
	    x_rot   = -0.1;
	  else if (event.xkey.state & ControlMask)
	    z_rot   = -0.1;
	  else
	    x_trans = -0.02 * window->view_state.distance;
	}
	else if (ks == XK_Up){
	  recalcModelView = GL_TRUE;
	  x_trans     = 0.0;
	  y_trans     = 0.0;
	  dist_change = 0.0;
	  approach    = 0.0;
	  x_rot       = 0.0;
	  y_rot       = 0.0;
	  z_rot       = 0.0;

	  if (event.xkey.state & ShiftMask)
	    y_rot   = -0.1;
	  else if (event.xkey.state & ControlMask)
	    dist_change = 0.02 * window->view_state.distance;
	  else
	    y_trans     = -0.02 * window->view_state.distance;
	}
	else if (ks == XK_Down){
	  recalcModelView = GL_TRUE;
	  x_trans     = 0.0;
	  y_trans     = 0.0;
	  dist_change = 0.0;
	  approach    = 0.0;
	  x_rot       = 0.0;
	  y_rot       = 0.0;
	  z_rot       = 0.0;

	  if (event.xkey.state & ShiftMask)
	    y_rot   = 0.1;
	  else if (event.xkey.state & ControlMask)
	    dist_change = -0.02 * window->view_state.distance;
	  else
	    y_trans     = 0.02 * window->view_state.distance;
	}
      }

/* end KeyPress */
      break;

    case ButtonPress:
      lastX = event.xbutton.x;
      lastY = event.xbutton.y;
      recalcModelView = GL_FALSE;

/* reset the change */
      x_trans     = 0.0;
      y_trans     = 0.0;
      dist_change = 0.0;
      approach    = 0.0;
      x_rot       = 0.0;
      y_rot       = 0.0;
      z_rot       = 0.0;
/* end ButtonPress */
      break;

    case MotionNotify:

/* check which buttons are held down */
      XQueryPointer(xDisplay, window->wMain, (Window *) &dummy, (Window *) &dummy, 
		    (int *) &dummy, (int *) &dummy, &dummy, &dummy, &button_mask);
      
/* filter out the interesting buttons */
      button_mask &= (Button1Mask | Button2Mask | Button3Mask);

      switch (button_mask){

      case Button1Mask:
	recalcModelView = GL_TRUE;
	x_trans += 0.001 * window->view_state.distance * (lastX - event.xmotion.x);
	y_trans -= 0.001 * window->view_state.distance * (lastY - event.xmotion.y);
	lastX = event.xbutton.x;
	lastY = event.xbutton.y;
	break;

      case Button2Mask:
	recalcModelView = GL_TRUE;
	x_rot += 0.002 * (lastX - event.xmotion.x);
	y_rot -= 0.002 * (lastY - event.xmotion.y);
	lastX = event.xbutton.x;
	lastY = event.xbutton.y;
	break;

      case (Button1Mask | Button2Mask):
	recalcModelView = GL_TRUE;
	z_rot += 0.005 * (lastX - event.xmotion.x);
	lastX = event.xbutton.x;
	lastY = event.xbutton.y;
	break;

      case (Button2Mask | Button3Mask):
	recalcModelView = GL_TRUE;
	dist_change  += 0.005 * window->view_state.distance * (lastY - event.xmotion.y);
	lastX = event.xbutton.x;
	lastY = event.xbutton.y;
	break;

      case Button3Mask:
	recalcModelView = GL_TRUE;
	approach += 0.005 * window->view_state.distance * (lastY - event.xmotion.y);
	lastX = event.xbutton.x;
	lastY = event.xbutton.y;
	break;

      default:
 	break;

      }

/* end MotionNotify */
      break;

    default:
      break;
    
    } /* end switch( event.type ) */

  } /* end while( XCheckMaskEvent )*/

/* observe that window is assigned if recalcModelView == TRUE */
  if (window != NULL && recalcModelView && window->n_active_lists > 0){ 
/* ogl_change_modelview make `window' the current rendering */
    ogl_change_modelview(window, x_trans, y_trans, dist_change, approach,
			 x_rot, y_rot, z_rot );
    recalcModelView = GL_FALSE;  
    x_trans = 0.0;
    y_trans = 0.0;
    dist_change    = 0.0;
    approach    = 0.0;
    x_rot   = 0.0;
    y_rot   = 0.0;
    z_rot   = 0.0;
    window->need_redraw = GL_TRUE;
  }

/* loop over all windows and update the plot */
  for (i=0; i<OGL_n_windows; i++){
    window = OGL_window_list[i];

/* redraw the plot if necessary */
    if (window->need_redraw){
/* ogl_change_modelview make `window' the current rendering */
      ogl_redraw(window);
    }
  }

  return OK;
} /* end ogl_poll */


int
ogl_start_plot(int new_plot, int win_number, double scale){
  WINDOW_REC *window;

  if (!ogl_is_running || win_number < 0 || win_number >= OGL_n_windows) 
    return ERROR;

  window = OGL_window_list[win_number];

/* bind the redering to the window */
  glXMakeCurrent(xDisplay, window->wMain, window->cMain);

  if (new_plot) 
    {
      ogl_erase(win_number);
      if (window->use_scaling)
	window->scale = 1.0;
      else
	window->scale = scale;
    }
  else
    {
      if (window->use_scaling)
	window->scale = 1.0;
      else if (scale > window->scale)
	window->scale = scale;
    }

/* adjust the projection transform to suit the size of the object */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0,   /* field of view in degree */ 
		 ((double) window->w)/((double) window->h), /* aspect ratio */ 
		 1.e-2 * window->scale,  /* Z near */ 
		 100.0 * window->scale   /* Z far */ ); 
/* the ratio 100/0.01 determines the z-buffer resolution */

/* the default is to modify the modelview matrix */
  glMatrixMode(GL_MODELVIEW);

/* generate a new object */
  if (window->n_active_lists < OGL_ACTIVE_LISTS){
    window->current_display_list[window->n_active_lists] = glGenLists( 1 );
  }
  else{
    printf("OGL: Sorry, can only store %i display lists\n", OGL_ACTIVE_LISTS);
    return ERROR;
  }
/* start saving graphics commands */
  glNewList(window->current_display_list[window->n_active_lists], GL_COMPILE);
  window->n_active_lists++;

/* set the current window */
  OGL_current_window = win_number;

  return OK;
}

int
ogl_end_plot(void){
  if (!ogl_is_running) 
    return ERROR;

  glEndList();

  if (OGL_current_window < 0 || OGL_current_window >= OGL_n_windows){
    printf("Error in ogl_end_plot: OGL_current_window = %i out of bounds\n", 
	   OGL_current_window);
    return ERROR;
  }

/* tell ogl_poll to call ogl_redraw as soon as possible */
  OGL_window_list[OGL_current_window]->need_redraw = GL_TRUE;

/* process the redraw event */
  ogl_poll();

  return OK;
}
  
void 
ogl_erase(int win_number){
  int i;

  if (!ogl_is_running || win_number < 0 || win_number >= OGL_n_windows) 
    return;

/* delete all current lists */
  for (i=0; i<OGL_window_list[win_number]->n_active_lists; i++)
    glDeleteLists( OGL_window_list[win_number]->current_display_list[i], 1);
  OGL_window_list[win_number]->n_active_lists = 0;
/* tell ogl_poll to call ogl_redraw as soon as possible */
  OGL_window_list[win_number]->need_redraw = GL_TRUE;
}


void
ogl_marker(double x, double y, double z, double radius, int color){
  GLUquadricObj *my_quad;
  GLint orientation;
/* create a new quadric object */
   my_quad = gluNewQuadric(); 

/* save the present orientation */
   glGetIntegerv(GL_FRONT_FACE, &orientation); 
/* use counterclockwise orientation for the sphere */
   glFrontFace( GL_CCW ); 

/* save the modelview matrix */
   glPushMatrix(); 
   gluQuadricDrawStyle( my_quad, (int) GLU_FILL ); 
   gluQuadricOrientation( my_quad, (int) GLU_OUTSIDE ); 
   ogl_set_color( color ); 
   glTranslated(x, y, z); 
   gluSphere(my_quad, radius, 10, 5); 
/* restore the modelview matrix */
   glPopMatrix(); 


/* reset the orientation */
   glFrontFace( orientation ); 

/* zap the object */
   gluDeleteQuadric(my_quad);
}

void 
ogl_set_color( int color ){
  GLfloat *diffuse, *ambient, *specular;

/* define some colors */
/* black */
  GLfloat black_diffuse[]  = { 0.1, 0.1, 0.1, 1.0 };
  GLfloat black_ambient[]  = { 0.02, 0.02, 0.02, 1.0 };
  GLfloat black_specular[] = { 0.1, 0.1, 0.1, 1.0 };
/* red */
  GLfloat red_diffuse[]  = { 0.7, 0.2, 0.0, 1.0 };
  GLfloat red_ambient[]  = { 0.1, 0.04, 0.0, 1.0 };
  GLfloat red_specular[] = { 0.1, 0.1, 0.1, 1.0 };
/* green */
  GLfloat green_diffuse[]  = { 0.05, 0.7, 0.0, 1.0 };
  GLfloat green_ambient[]  = { 0.01, 0.2, 0.0, 1.0 };
  GLfloat green_specular[] = { 0.1, 0.1, 0.1, 1.0 };
/* yellow */
  GLfloat yellow_diffuse[]  = { 0.7, 0.5, 0.0, 1.0 };
  GLfloat yellow_ambient[]  = { 0.2, 0.15, 0.0, 1.0 };
  GLfloat yellow_specular[] = { 0.1, 0.1, 0.1, 1.0 };
/* blue */
  GLfloat blue_diffuse[]  = { 0.0, 0.2, 0.7, 1.0 };
  GLfloat blue_ambient[]  = { 0.0, 0.04, 0.2, 1.0 };
  GLfloat blue_specular[] = { 0.1, 0.1, 0.1, 1.0 };
/* magenta */
  GLfloat magenta_diffuse[]  = { 0.7, 0.0, 0.7, 1.0 };
  GLfloat magenta_ambient[]  = { 0.2, 0.0, 0.2, 1.0 };
  GLfloat magenta_specular[] = { 0.1, 0.1, 0.1, 1.0 };
/* cyan */
  GLfloat cyan_diffuse[]  = { 0.0, 0.7, 0.7, 1.0 };
  GLfloat cyan_ambient[]  = { 0.0, 0.2, 0.2, 1.0 };
  GLfloat cyan_specular[] = { 0.1, 0.1, 0.1, 1.0 };
/* gold */
  GLfloat gold_diffuse[]  = { 0.64, 0.36, 0.0, 1.0 };
  GLfloat gold_ambient[]  = { 0.012, 0.06, 0.0, 1.0 };
  GLfloat gold_specular[] = { 0.1, 0.1, 0.1, 1.0 };
/* white */
  GLfloat white_diffuse[]  = { 0.7, 0.7, 0.7, 1.0 };
  GLfloat white_ambient[]  = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat white_specular[] = { 0.1, 0.1, 0.1, 1.0 };

/* choose color */
  switch( color % 9 ){

  case OGL_BLACK:
    diffuse = black_diffuse;
    ambient = black_ambient;
    specular = black_specular;
    break;

  case OGL_RED:
    diffuse = red_diffuse;
    ambient = red_ambient;
    specular = red_specular;
    break;

  case OGL_GREEN:
    diffuse = green_diffuse;
    ambient = green_ambient;
    specular = green_specular;
    break;

  case OGL_YELLOW:
    diffuse = yellow_diffuse;
    ambient = yellow_ambient;
    specular = yellow_specular;
    break;

  case OGL_BLUE:
    diffuse = blue_diffuse;
    ambient = blue_ambient;
    specular = blue_specular;
    break;

  case OGL_MAGENTA:
    diffuse = magenta_diffuse;
    ambient = magenta_ambient;
    specular = magenta_specular;
    break;

  case OGL_CYAN:
    diffuse = cyan_diffuse;
    ambient = cyan_ambient;
    specular = cyan_specular;
    break;

  case OGL_GOLD:
    diffuse = gold_diffuse;
    ambient = gold_ambient;
    specular = gold_specular;
    break;

  case OGL_WHITE:
  default:
    diffuse = white_diffuse;
    ambient = white_ambient;
    specular = white_specular;
    break;

  }

/* set color */
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
  glMateriali( GL_FRONT_AND_BACK, GL_SHININESS, 10);

}

char *
ogl_color_name( int color ){
/* choose color */
  switch( color % 9 ){

  case OGL_BLACK:
    return "Black";
    break;

  case OGL_RED:
    return "Red";
    break;

  case OGL_GREEN:
    return "Green";
    break;

  case OGL_YELLOW:
    return "Yellow";
    break;

  case OGL_BLUE:
    return "Blue";
    break;

  case OGL_MAGENTA:
    return "Magenta";
    break;

  case OGL_CYAN:
    return "Cyan";
    break;

  case OGL_GOLD:
    return "Gold";
    break;

  case OGL_WHITE:
    return "White";
    break;

  default:
    break;
  }

  return "Unknown";
}

/******************************************************************************/

void 
ogl_drawstr(double x_pos[3], double n_vec[3], char *str)
{
  if (OGL_current_window < 0 || OGL_current_window >= OGL_n_windows)
    return;

/* draw a point twice to set the lighting (is this a bug or what?)*/
  glBegin( GL_POINTS );
  glNormal3dv( n_vec );
  glVertex3dv( x_pos );
  glNormal3dv( n_vec ); 
  glVertex3dv( x_pos );
  glEnd();

  glRasterPos3dv(x_pos);

  glPushAttrib(GL_LIST_BIT);
  glListBase(OGL_window_list[OGL_current_window]->bitmapBase);
  glCallLists(strlen(str), GL_UNSIGNED_BYTE, (unsigned char *)str);
  glPopAttrib();
}


static void
vector_cross_prod( double a[3], double b[3], double c[3] ){
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void
ogl_draw_arrow(double length, double thickness, double start[3], double forward[3],
	       double normal[3], char *label){
  double h_side[3], v_side[3], top_first[3], first_sec[3], n_vec[3], theta, phi;
  double height, base, x_pos[3];
  int i;
  GLint orientation;
  GLboolean was_normalized;
  WINDOW_REC *window;

/* save a pointer to the current window */
  window = OGL_window_list[OGL_current_window];
  if (window->use_scaling)
    length = 0.1;

/* shape of the arrow */
  height = 0.5 * length;
  base = 0.25*height;

/* save the present orientation */
  glGetIntegerv(GL_FRONT_FACE, &orientation);

/* use clockwise orientation for the surface normals */
  glFrontFace( GL_CW );
  was_normalized = glIsEnabled( GL_NORMALIZE );
  glEnable( GL_NORMALIZE );

/* it remains to set the normal */
/* compute orthogonal basis */
  theta = acos(forward[2]);
  if (fabs(theta) > 1.e-4 && fabs(forward[1]/sin(theta)) <= 1.0){
    if (forward[0]/sin(theta)>=0.0)
      phi = asin( forward[1]/sin(theta) );
    else
      phi = ((double)M_PI) - asin( forward[1]/sin(theta) );
  }
  else{
    phi = 0.0;
  }
/*   forward[0] = sin(theta) * cos(phi) */
/*   forward[1] = sin(theta) * sin(phi) */
/*   forward[2] = cos(theta)            */
  h_side[0]  =-sin(phi);
  h_side[1]  = cos(phi);
  h_side[2]  = 0.0;
  v_side[0]  = cos(theta) * cos(phi);
  v_side[1]  = cos(theta) * sin(phi);
  v_side[2]  =-sin(theta);

/* compute the first normal */
  for (i=0; i<3; i++){
    top_first[i] = -height*forward[i] - base*h_side[i];
    first_sec[i] = base*(h_side[i] - v_side[i]);
  }
  vector_cross_prod( first_sec, top_first, n_vec );

/* draw the stem */

/* set the line width */
  glLineWidth(thickness);
  glEnable(GL_LINE_SMOOTH);
  glBegin(GL_LINES);

/* set the normal */
  glNormal3dv( normal );
/* set current vertex */
  glVertex3dv(start);
  glVertex3d(start[0] + length*forward[0]/window->x_scale, 
	     start[1] + length*forward[1]/window->y_scale, 
	     start[2] + length*forward[2]/window->z_scale);
  glEnd();

/* reset the line width */
  glLineWidth(1.0);
  glDisable(GL_LINE_SMOOTH);

/* draw string */
  if (label != NULL){
/* set the raster position */
    x_pos[0] = start[0] + 1.2*length*forward[0]/window->x_scale;
    x_pos[1] = start[1] + 1.2*length*forward[1]/window->y_scale;
    x_pos[2] = start[2] + 1.2*length*forward[2]/window->z_scale;
    ogl_drawstr(x_pos, normal, label);
  }

/* draw the faces of the arrow */
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_TRIANGLES); 

/* set the normal */
  glNormal3dv( n_vec );

  glVertex3d(start[0] + length*forward[0]/window->x_scale, 
	     start[1] + length*forward[1]/window->y_scale,  
	     start[2] + length*forward[2]/window->z_scale); 

  glVertex3d(start[0] + ((length-height)*forward[0] - base*h_side[0])/window->x_scale,
	     start[1] + ((length-height)*forward[1] - base*h_side[1])/window->y_scale, 
	     start[2] + ((length-height)*forward[2] - base*h_side[2])/window->z_scale); 

  glVertex3d(start[0] + ((length-height)*forward[0] - base*v_side[0])/window->x_scale, 
	     start[1] + ((length-height)*forward[1] - base*v_side[1])/window->y_scale,  
	     start[2] + ((length-height)*forward[2] - base*v_side[2])/window->z_scale); 
  glEnd(); 

/* compute the second normal */
  for (i=0; i<3; i++){
    top_first[i] = -height*forward[i] - base*v_side[i];
    first_sec[i] = base*(h_side[i] + v_side[i]);
  }
  vector_cross_prod( first_sec, top_first, n_vec );

  glBegin(GL_TRIANGLES); 
/* set the normal */
  glNormal3dv( n_vec );

  glVertex3d(start[0] + length*forward[0]/window->x_scale, 
	     start[1] + length*forward[1]/window->y_scale,  
	     start[2] + length*forward[2]/window->z_scale); 

  glVertex3d(start[0] + ((length-height)*forward[0] - base*v_side[0])/window->x_scale, 
	     start[1] + ((length-height)*forward[1] - base*v_side[1])/window->y_scale,  
	     start[2] + ((length-height)*forward[2] - base*v_side[2])/window->z_scale); 

  glVertex3d(start[0] + ((length-height)*forward[0] + base*h_side[0])/window->x_scale, 
	     start[1] + ((length-height)*forward[1] + base*h_side[1])/window->y_scale,  
	     start[2] + ((length-height)*forward[2] + base*h_side[2])/window->z_scale); 
  glEnd(); 

/* compute the third normal */
  for (i=0; i<3; i++){
    top_first[i] = -height*forward[i] + base*h_side[i];
    first_sec[i] = base*(-h_side[i] + v_side[i]);
  }
  vector_cross_prod( first_sec, top_first, n_vec );

  glBegin(GL_TRIANGLES); 

/* set the normal */
  glNormal3dv( n_vec );
  
  glVertex3d(start[0] + length*forward[0]/window->x_scale,  
	     start[1] + length*forward[1]/window->y_scale,  
	     start[2] + length*forward[2]/window->z_scale); 

  glVertex3d(start[0] + ((length-height)*forward[0] + base*h_side[0])/window->x_scale, 
	     start[1] + ((length-height)*forward[1] + base*h_side[1])/window->y_scale,  
	     start[2] + ((length-height)*forward[2] + base*h_side[2])/window->z_scale); 

  glVertex3d(start[0] + ((length-height)*forward[0] + base*v_side[0])/window->x_scale,  
	     start[1] + ((length-height)*forward[1] + base*v_side[1])/window->y_scale,   
	     start[2] + ((length-height)*forward[2] + base*v_side[2])/window->z_scale);  

  glEnd(); 

/* compute the fourth normal */
  for (i=0; i<3; i++){
    top_first[i] = -height*forward[i] + base*v_side[i];
    first_sec[i] = base*(-h_side[i] - v_side[i]);
  }
  vector_cross_prod( first_sec, top_first, n_vec );

  glBegin(GL_TRIANGLES); 

/* set the normal */
  glNormal3dv( n_vec );

  glVertex3d(start[0] + length*forward[0]/window->x_scale,  
	     start[1] + length*forward[1]/window->y_scale,  
	     start[2] + length*forward[2]/window->z_scale); 

  glVertex3d(start[0] + ((length-height)*forward[0] + base*v_side[0])/window->x_scale,  
	     start[1] + ((length-height)*forward[1] + base*v_side[1])/window->y_scale,   
	     start[2] + ((length-height)*forward[2] + base*v_side[2])/window->z_scale);  

  glVertex3d(start[0] + ((length-height)*forward[0] - base*h_side[0])/window->x_scale,  
	     start[1] + ((length-height)*forward[1] - base*h_side[1])/window->y_scale,   
	     start[2] + ((length-height)*forward[2] - base*h_side[2])/window->z_scale);  
  glEnd();

/* reset the orientation */
  glFrontFace( orientation );

/* turn off the normalization if it was turned off before */
  if ( !was_normalized )
    glDisable( GL_NORMALIZE );
}

double
ogl_length_scale( bounding_box *bb )
{
  double length;
  length = (bb->x_max - bb->x_min > bb->y_max - bb->y_min)? 
    bb->x_max - bb->x_min : bb->y_max - bb->y_min;
  length = (length > bb->z_max - bb->z_min)? length: bb->z_max - bb->z_min;
  return length;
}

void 
ogl_bb_cage( bounding_box *bb, int color ){
  double start[3], forward[3], length, height, base;
  double n_vec[3];
  int was_normalized;
  char origin[80];
/* arrow length */
  length = 0.1 * ogl_length_scale( bb );
  height = 0.5 * length;
  base   = 0.5 * height;

/* set the color */
  ogl_set_color( color );

/* x-arrow */
  start[0] = bb->x_max;
  start[1] = bb->y_min;
  start[2] = bb->z_min;
  forward[0] = 1.0;
  forward[1] = 0.0;
  forward[2] = 0.0;
  n_vec[0] = 0.0;
  n_vec[1] =-1.0;
  n_vec[2] =-1.0;

  ogl_draw_arrow(length, 1.0, start, forward, n_vec, "x");

/* y-arrow */
  start[0] = bb->x_min;
  start[1] = bb->y_max;
  start[2] = bb->z_min;
  forward[0] = 0.0;
  forward[1] = 1.0;
  forward[2] = 0.0;
  n_vec[0] =-1.0;
  n_vec[1] = 0.0;
  n_vec[2] =-1.0;

  ogl_draw_arrow(length, 1.0, start, forward, n_vec, "y");

/* z-arrow*/
  start[0] = bb->x_min;
  start[1] = bb->y_min;
  start[2] = bb->z_max;
  forward[0] = 0.0;
  forward[1] = 0.0;
  forward[2] = 1.0;
  n_vec[0] =-1.0;
  n_vec[1] =-1.0;
  n_vec[2] = 0.0;

  ogl_draw_arrow(length, 1.0, start, forward, n_vec, "z");

/* draw the coordinates of the intersection of the x,y,z axes */
  start[0] = bb->x_min;
  start[1] = bb->y_min;
  start[2] = bb->z_min;
  n_vec[0] = 1.0;
  n_vec[1] = 1.0;
  n_vec[2] = 1.0;
  sprintf(origin,"(%.3e, %.3e, %.3e)", bb->x_min, bb->y_min, bb->z_min);
  ogl_drawstr( start, n_vec, origin );

  start[0] = bb->x_max;
  start[1] = bb->y_max;
  start[2] = bb->z_max;
  n_vec[0] = 1.0;
  n_vec[1] = 1.0;
  n_vec[2] = 1.0;
  sprintf(origin,"(%.3e, %.3e, %.3e)", bb->x_max, bb->y_max, bb->z_max);
  ogl_drawstr( start, n_vec, origin );

/* draw the cage */

/* let gl normalize */
  if ( !(was_normalized = glIsEnabled( GL_NORMALIZE )) )
    glEnable( GL_NORMALIZE );

/* x-direction */
  glBegin(GL_LINES);
  glNormal3f( 0.0, -1.0, -1.0 );
  glVertex3d(bb->x_min, bb->y_min, bb->z_min);
  glVertex3d(bb->x_max, bb->y_min, bb->z_min);
  glEnd();

  glBegin(GL_LINES);
  glNormal3f( 0.0, -1.0, 1.0 );
  glVertex3d(bb->x_min, bb->y_min, bb->z_max);
  glVertex3d(bb->x_max, bb->y_min, bb->z_max);
  glEnd();

  glBegin(GL_LINES);
  glNormal3f( 0.0, 1.0, -1.0 );
  glVertex3d(bb->x_min, bb->y_max, bb->z_min);
  glVertex3d(bb->x_max, bb->y_max, bb->z_min);
  glEnd();

  glBegin(GL_LINES);
  glNormal3f( 0.0, 1.0, 1.0 );
  glVertex3d(bb->x_min, bb->y_max, bb->z_max);
  glVertex3d(bb->x_max, bb->y_max, bb->z_max);
  glEnd();
 
/* y-direction */
  glBegin(GL_LINES);
  glNormal3f( -1.0, 0.0, -1.0 );
  glVertex3d(bb->x_min, bb->y_min, bb->z_min);
  glVertex3d(bb->x_min, bb->y_max, bb->z_min);
  glEnd();

  glBegin(GL_LINES);
  glNormal3f( 1.0, 0.0, -1.0 );
  glVertex3d(bb->x_max, bb->y_min, bb->z_min);
  glVertex3d(bb->x_max, bb->y_max, bb->z_min);
  glEnd();

  glBegin(GL_LINES);
  glNormal3f( -1.0, 0.0, 1.0 );
  glVertex3d(bb->x_min, bb->y_min, bb->z_max);
  glVertex3d(bb->x_min, bb->y_max, bb->z_max);
  glEnd();

  glBegin(GL_LINES);
  glNormal3f( 1.0, 0.0, 1.0 );
  glVertex3d(bb->x_max, bb->y_min, bb->z_max);
  glVertex3d(bb->x_max, bb->y_max, bb->z_max);
  glEnd();

/* z-direction */
  glBegin(GL_LINES);
  glNormal3f(-1.0, -1.0, 0.0 );
  glVertex3d(bb->x_min, bb->y_min, bb->z_min);
  glVertex3d(bb->x_min, bb->y_min, bb->z_max);
  glEnd();

  glBegin(GL_LINES);
  glNormal3f( 1.0, -1.0, 0.0 );
  glVertex3d(bb->x_max, bb->y_min, bb->z_min);
  glVertex3d(bb->x_max, bb->y_min, bb->z_max);
  glEnd();

  glBegin(GL_LINES);
  glNormal3f(-1.0, 1.0, 0.0 );
  glVertex3d(bb->x_min, bb->y_max, bb->z_min);
  glVertex3d(bb->x_min, bb->y_max, bb->z_max);
  glEnd();

  glBegin(GL_LINES);
  glNormal3f( 1.0, 1.0, 0.0 );
  glVertex3d(bb->x_max, bb->y_max, bb->z_min);
  glVertex3d(bb->x_max, bb->y_max, bb->z_max);
  glEnd();

/* turn off the normalization if it was turned off before */
  if ( !was_normalized )
    glDisable( GL_NORMALIZE );
}

/***************************************************************************
 *
 * Start of command interpreter 
 *
 **************************************************************************/



/* private functions. make them static so they can't be used outside this file */
static char *gettok(input_output *io_ptr, char *prompt,
                    char **cmd, int n);
static char *getline_no_com( char *s, int n, input_output *io_ptr );


int 
get_command(input_output *io_ptr, char *prompt,
	    char **command, int n_command,
	    int command_level, int *save_on_copy, int *argument) {
  char *token, *cnull, *first_token, *second_token, *second_command=NULL;
  char *third_token=NULL, *third_command=NULL, *par1;
  int match, unique, i, first_length, second_length, third_length;
  input_output pause_io;
  int quit, pause_icom;

  enum {PROCEED, BREAK, HELP, INTERACT};

  const int LAST_COM = INTERACT;
  char *pause_com[4], *pause_help[4];

  match = -1;
  unique = 1;

  token = gettok( io_ptr, prompt, command, n_command + 1);

  /* check for "?" */
  
  while (token[0] == '?' || token[0] == '\n' ) {
    if (token[0] == '?'){
      printf("Choose one of the following commands:\n");
      for (i=0; i<=n_command; i++){
	printf("%s\n",command[i]);
	if ((i+1)%3 == 0) printf("\n");
      }
    }
/* read a new answer */
    token = gettok( io_ptr, prompt, command, n_command + 1 );
  }

/* check for empty command */

  if (token[0] == '\n') 
    return -1;
  else{      /* replace any '\n' by '\0' */
    if ( (cnull = strchr(token,'\n')) != NULL )
      cnull[0] = '\0';
  }

/* check for a pause */
  if (strcmp(token, "pause") == 0){
    pause_io.read_command = stdin;
    pause_io.copy_command = NULL;

    pause_com[PROCEED] = "proceed";
    pause_com[BREAK] = "break";
    pause_com[HELP] = "help";
    pause_com[INTERACT] = "command-interpreter-introduction";

    pause_help[PROCEED] = "Continue with the tutorial, i.e. read next instruction from "
      "the command file.";
    pause_help[BREAK] = "Terminate the tutorial session, i.e. close the command file "
      "and start reading commands from standard input.";
    pause_help[HELP] = NULL;
    pause_help[INTERACT] = "Describes how 3-D objects are translated and rotated on "
      "the screen, and shows the basic features of the command interpreter.";

    quit = 0;
    do{
      pause_icom = get_command( &pause_io, "Pause>", pause_com, LAST_COM, 
			       command_level+1, NULL, NULL);

      switch (pause_icom){

      case PROCEED:
/* proceed */
	quit = 1;
	break;

      case BREAK:
/* break */
	fclose(io_ptr->read_command);
	io_ptr->read_command = stdin;
	quit = 1;
	break;

      case HELP:
/* help */
	while ((pause_icom = get_command( &pause_io, 
					 "help on subject (? lists all subjects)>", 
					 pause_com, 2, 0, NULL, NULL)) == -1);
	if (pause_help[pause_icom] == NULL)
	  general_help(); 
	else
	  print_help( pause_help[pause_icom] ); 
	break;


      case INTERACT:
	general_help();
	break;

      default:
	;
      }

    }
    while( !quit );
    return -1;
  }	 
  
  /* decompose the token into its sub parts */
    
  first_token = token;
  first_length = strcspn( token, "-" );

  if ( (second_token = strchr( token, '-' )) == NULL ){
    second_length = 0;
    third_length = 0;
  } else {
    second_length = strcspn( &(second_token[1]), "-" ) + 1;
    if ( (third_token = strchr( &(second_token[1]), '-' )) == NULL ){
      third_length = 0;
    }
    else{
      third_length = strcspn( &(third_token[1]), "-" ) + 1;
    }
  }
  
/* check for matches */

  for (i=0; i<=n_command; i++){
/* new algorithm */
/* decompose the command */
    if (second_length > 0){
      second_command = strchr( command[i], '-' );
      third_command = NULL;
      if ( third_length > 0){
	if (second_command == NULL)
	  third_command = NULL;
	else
	  third_command = strchr( &(second_command[1]), '-' );
      }
    }
/* the first part of the command matches, then check if either the second part
 is absent or if it also matches */
    if ( strncmp( first_token, command[i], first_length)==0 &&
	(second_length == 0 ||
	 (second_command != NULL && strncmp( second_token, second_command, 
					    second_length )==0)) &&
	(third_length == 0 ||
	 (third_command != NULL && strncmp( third_token, third_command, 
					   third_length )==0)) ){
	if (match == -1)  
	  unique = 1; 
	else 
	  unique = 0; 
	match = i; 
      }
  }
/* even if the command wasn't unique, it may be a perfect match. This happens
   for instance when commmand[0]="s", command[1]="show" and token="s". */
  if (match != -1 && !unique){
    match = -1;
    for (i=0; i<=n_command; i++)
      if (strcmp( token, command[i] )==0){
	match=i;
      }
  }


  if (match==-1){
    if (!unique)
      printf("Not an unique command: %s\n", token);
    else
      printf("Unknown command: %s\n", token);
  }

/* copy the command onto the copy stream */
  if (match != -1 && io_ptr->copy_command != NULL &&
      (save_on_copy==NULL || save_on_copy[match])){

/* indent accordning to command level */
    if (strcmp("exit", command[match])==0)
      command_level += -1;
    for (i=0; i< command_level; i++)
      fprintf( io_ptr->copy_command, "\t");

/* output the command excluding stuff between () */
    if ( (par1 = strchr( command[match], '(' )) != NULL)
      *par1 = '\0';
    fprintf( io_ptr->copy_command, "%s", command[match] );
/* can we expect an argument to this command which should be printed on the
   same line? */
    if (argument == NULL || !argument[match])
      fprintf( io_ptr->copy_command, "\n");
    else
      fprintf( io_ptr->copy_command, " ");
  }

/* echo a valid command on stdout */
  if (match != -1 ){
/* output the command */
    printf( "%s", command[match] );
    if ( io_ptr->read_command == stdin )
      printf( "\n" );
    else{
/* if the command was read from a file:
   can we expect an argument to this command which should be printed on the
   same line? */
      if (argument == NULL || !argument[match])
	printf( "\n" );
      else
	printf( " " );
    }
  }

  return match;

}


char *
get_word(input_output *io_ptr, char *prompt, char *deflt,
	  int save_on_copy){
  char *token, *cnull, prompt2[120];

  sprintf(prompt2, "%s (%s)", prompt, deflt );
  do{
    token = gettok( io_ptr, prompt2, NULL, 0 );
  }    
  while ( token[0] == '?' );

  if (token[0] == '\n') 
    token = deflt;
  else{
/* replace any '\n' by '\0' */
    if ( (cnull = strchr(token,'\n')) != NULL )
      cnull[0] = '\0';
  }

/* copy it onto the copy stream */
  if ( io_ptr->copy_command != NULL && save_on_copy )
    fprintf( io_ptr->copy_command, "%s\n", token );

/* echo the command on stdout if the command was read from file */
  if ( io_ptr->read_command != stdin )
    printf( "%s\n", token );

  return token;
}

int 
get_int(input_output *io_ptr, char *prompt, int deflt, int level){
  char *token, *cnull, prompt2[120];
  int i;

  sprintf(prompt2, "%s (%i)", prompt, deflt );
  do{
    token = gettok( io_ptr, prompt2, NULL, 0 );
/* check that token consists only of (+/-) and digits. Otherwise return NULL */
  }    
  while ( token[0] == '?' );

  if (token[0] == '\n') 
    sprintf(token,"%i", deflt);
  else{
/* replace any '\n' by '\0' */
    if ( (cnull = strchr(token,'\n')) != NULL )
      cnull[0] = '\0';
  }

/* copy it onto the copy stream */
  if ( io_ptr->copy_command != NULL ){
/* indent according to the command level */
    for (i=0; i< level; i++)
      fprintf( io_ptr->copy_command, "\t");
    fprintf( io_ptr->copy_command, "%s\n", token );
  }

/* echo the command on stdout if the command was read from file */
  if ( io_ptr->read_command != stdin )
    printf( "%s\n", token );

  return atoi(token);
}

real 
get_real(input_output *io_ptr, char *prompt, real deflt, int level){
  char *token, *cnull, prompt2[120];
  int i;

  sprintf(prompt2, "%s (%e)", prompt, deflt );
  do{
    token = gettok( io_ptr, prompt2, NULL, 0 );
/* check that token can be converted to a real number. Otherwise return NULL */
  }    
  while ( token[0] == '?' );

  if (token[0] == '\n') 
    sprintf(token,"%e", deflt);
  else{
/* replace any '\n' by '\0' */
    if ( (cnull = strchr(token,'\n')) != NULL )
      cnull[0] = '\0';
  }

/* copy it onto the copy stream */
  if ( io_ptr->copy_command != NULL ){
/* indent according to the command level */
    for (i=0; i< level; i++)
      fprintf( io_ptr->copy_command, "\t");
/* output the number */
    fprintf( io_ptr->copy_command, "%s\n", token );
  }

/* echo the command on stdout if the command was read from file */
  if ( io_ptr->read_command != stdin )
    printf( "%s\n", token );

  return atof(token);
}

float get_float(input_output *io_ptr, char *prompt, float deflt, int level){
  char *token, *cnull, prompt2[120];
  int i;

  sprintf(prompt2, "%s (%e)", prompt, deflt );
  do{
    token = gettok( io_ptr, prompt2, NULL, 0 );
/* check that token can be converted to a float number. Otherwise return NULL */
  }    
  while ( token[0] == '?' );

  if (token[0] == '\n') 
    sprintf(token,"%e", deflt);
  else{
/* replace any '\n' by '\0' */
    if ( (cnull = strchr(token,'\n')) != NULL )
      cnull[0] = '\0';
  }

/* copy it onto the copy stream */
  if ( io_ptr->copy_command != NULL ){
/* indent according to the command level */
    for (i=0; i< level; i++)
      fprintf( io_ptr->copy_command, "\t");
/* output the number */
    fprintf( io_ptr->copy_command, "%s\n", token );
  }

/* echo the command on stdout if the command was read from file */
  if ( io_ptr->read_command != stdin )
    printf( "%s\n", token );

  return atof(token);
}

int 
get_yes_no(input_output *io_ptr, char *prompt, int level, int save_answer ){
  const int ncom=1;
  char *command[2];
  int icom, quit, yes_no;
  int *argument=NULL, save_on_copy[2];

  command[0] ="yes"; save_on_copy[0] = save_answer;
  command[1] ="no";  save_on_copy[1] = save_answer;

  quit = 0;
  yes_no = 0;
  do{
    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    switch (icom) {
    case 0: 
      yes_no = 1;
      quit = 1;
      break;
    case 1:
      yes_no = 0;
      quit = 1;
      break;
    default:
      ;
    }
  }
  while(!quit);

  return yes_no;
}

int 
get_on_off(input_output *io_ptr, char *prompt){
  const int ncom=1;
  char *command[2];
  int icom, quit, on_off;
  int *argument=NULL, *save_on_copy=NULL, level=-1;

  command[0] ="on";
  command[1] ="off";

  quit = 0;
  on_off = 0;
  do{
    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    switch (icom) {
    case 0: 
      on_off = 1;
      quit = 1;
      break;
    case 1:
      on_off = 0;
      quit = 1;
      break;
    default:
      ;
    }
  }
  while(!quit);

  return on_off;
}

static char *gettok(input_output *io_ptr, char *prompt,
		    char **command, int n_command) {
  static char line[MAX_LINE_SIZE];
  char *token;
  static int first_call=TRUE;

/* when strtok is called with line as first parameter, it starts parsing from the
   first character of line and reads until it finds one of the non-tokens that are
   given as the second parameter. When strtok is called with NULL as the first 
   parameter, it continues to read from the previous line (which it stored internally)
   until it finds one of the non-tokens. We want the following functionality from
   this routine:
   when a line only contains a newline, it should be returned as a token.
   when a line contains other characters, newline should not be regarded as a token.
   This is because if the user by misstake gives the line `command <space> <return>',
   we want it to be interpreted only as `command'.
   To accomplish this, we give ' \t\n' as non-tokens when NULL is the first argument
   to strtok, while we only give ' \t' as non-tokens when 'line' is the first 
   argument. */

/* give the graphics a chance to update the plot for size changes */
  ogl_poll();

  if (first_call || (token = strtok(NULL," \t\n")) == NULL){
    first_call=FALSE;
    do{
/* output the prompt only if we read from stdin */
      if (io_ptr->read_command == stdin){
	getline_no_block( prompt, line, command, n_command );
      }
      else{
	getline_no_com( line, 120, io_ptr );
	ogl_poll();
      }
    }
    while( (token = strtok(line," \t"))==NULL );
  }
  return token;
}

static char *getline_no_com(char *line_no_comment, int buf_size, 
			    input_output *io_ptr ){
  char *first_comment;
  int i;

  if (fgets( line_no_comment, buf_size, io_ptr->read_command ) == NULL){
/* switch to stdin at end-of-file. */
    if (feof(io_ptr->read_command)){
/* If end-of-file occurs on standard in, try to reset it. */
      if (io_ptr->read_command == stdin){
	clearerr( stdin );
	printf("\n");
      }
/* if end-of-file occurs for a command file, close the command file and switch to stdin */
      else{
	fclose(io_ptr->read_command);
	io_ptr->read_command = stdin;
      }
/* important to reset the line */
      line_no_comment[0]='\0';
    }
  }

/* is there a comment ? */
  if ((first_comment=strchr( line_no_comment, '#' ))!=NULL){
/* echo the comment on stdout */
    printf("%s", first_comment);
/* overwrite the first occurance of the comment by '\0' */
/* and fill the rest of the string with blanks */
    for (i=0; first_comment[i] != '\0'; i++)
      first_comment[i] = ' ';
    *first_comment = '\0';
  }

/* return the pointer to the possibly modified string */
  return line_no_comment;
}


void print_help( char *brief_help ){
  int i, pos, length;
  const int min_length=70;

  if (brief_help == NULL) return;

  length = strlen(brief_help);

  pos = 0;
  printf("\n");
  do{
/* print one line that is at least 60 characters wide (if the string is
   sufficiently wide) and thereafter breaks the line at the first space.
*/
    for (i=0; i<min_length && pos < length && brief_help[pos] != '\n'; i++ )
      printf("%c", brief_help[pos++]);
/* wait with the newline until there is a space. Don't print the space */
    while (pos < length && brief_help[pos++] != ' ' && brief_help[pos-1] != '\n')
      printf("%c", brief_help[pos-1]);
    printf("\n");
  }
  while (pos<length);
}
	

void general_help( void ){
  printf("\n"
"The text window employs a tcsh-style command interpreter with the following \n"
"features:\n"
" o  TAB completion of commands.\n"
" o  Upwards arrow or ^P traverses up in the command history list.\n"
" o  Downwards arrow or ^N goes down in the command history list.\n"
" o  Only the unique part between each `-' of a command needs to be given.\n"
" o  ? lists all commands at the current level.\n"
" o  !cmnd executes the shell command `cmnd'\n\n"
"The graphics windows have the following functionality. Moving the pointer\n"
" o  while the LEFT mouse button is pressed translates the object in the window \n"
"    plane,\n"
" o  while the CENTER mouse button is pressed rotates the object around the \n"
"    vertical and horizontal axes, respectively,\n"
" o  horizontally, while both the LEFT and the CENTER mouse buttons are pressed,\n"
"    rotates the object in the window plane, and\n"
" o  vertically, while the RIGHT mouse button is pressed moves the view point\n"
"    closer or further away from the object.\n\n"
"The view point can also be moved with the arrow keys as follows:\n"
" o  Left, Right, Up or Down arrows moves the object in the corresponding \n"
"    direction in the window plane.\n"
" o  Holding down the SHIFT modifier while pressing Left, Right, Up or Down \n"
"    rotates the object around the vertical and the horizontal axes, \n"
"    respectively.\n"
" o  Holding down the CONTROL modifier while pressing Up or Down moves the view \n"
"    point further away from or closer to the object, while the Left or Right \n"
"    keys make the object rotate in the window plane.\n"
);
}
 
/********************************** END command interpreter ************************/

/********************************** START no-block *********************************/


static void
set_raw(int on){
  struct termio  arg ;
  static char    min ;
  static char    this_time;
  
  if (on != RAWOFF)
    {
      if (ioctl (fileno(stdin), TCGETA, &arg)==-1)
	perror("ioctl");
      arg.c_lflag &= ~ICANON & ~ECHO;
      min  = arg.c_cc[VMIN];
      this_time = arg.c_cc[VTIME];
      arg.c_cc[VMIN]  = (on == RAWON) ? 1 : 0;
/* block the read for 0.1 seconds instead of 0.0 seconds which made the program
   a real cpu eater when it was idle */
      arg.c_cc[VTIME] = 1;
      if (ioctl (fileno(stdin), TCSETA, &arg)==-1)
	perror("ioctl");
    }
  else
    {
      if (ioctl (fileno(stdin), TCGETA, &arg)==-1)
	perror("ioctl");
      arg.c_lflag |= ICANON | ECHO;
      arg.c_cc[VMIN]  = min ;
      arg.c_cc[VTIME] = this_time;
      if (ioctl (fileno(stdin), TCSETA, &arg)==-1)
	perror("ioctl");
    }
}

/* Returns how many (initial) charachters of strings s1 and s2 match */
/* i=strlen(s1)=strlen(s2) if s1==s2 */

static int
strolof(char *s1, char *s2) {
  int i=0;
  while (*s1 != 0 && *s2 != 0 && *s1++ == *s2++)
    i++;
  return i;
}


void
getline_no_block( char *prompt, char *buf, char **cmds, int n) {
  int c, i, k, cancel_flag=FALSE, last_esc, cb_pos, new_pos;
  unsigned char ch, bell=7;
  static char backsp[3] = {8, 32, 8};
  char compl[MAX_LINE_SIZE];
  static char cmd_buf[SCROLLBACK_SIZE][MAX_LINE_SIZE];
  static int cmd_buf_pos = 0, cmd_buf_full = FALSE;

  do {

  last_esc = -1;
  printf("%s", prompt);
  fflush(stdout);

  for (i=0; i<n; i++)
    if (!strcmp(cmds[i], "cancel")) {
      cancel_flag = TRUE;
      i=n;
    }

  set_raw(RAWNDELON);
    
  c=0;
  cb_pos = cmd_buf_pos;
  cmd_buf[cb_pos][0] = 0; /* Clear current slot in command buffer */

  do {
    while (read(fileno(stdin), &ch, 1) == 0)
      ogl_poll();
    
    /* Trap ESC-sequences generated by arrow keys */

    if (last_esc == 0 && ch == '[') {
      last_esc++;
      ch=0;
    } else {
      if (last_esc == 1 && ch == 'A') ch = 'P' & 31; /* UP   --> CTRL-P */
      if (last_esc == 1 && ch == 'B') ch = 'N' & 31; /* DOWN --> CTRL-N */
      last_esc = -1;
    }
    if (ch == 27)
      last_esc = 0;

    /* Print printable characters */

    if (ch >= 32 && ch < 256 && ch != 127) { /* Should be able to use 8 bit input */
      buf[c++]=ch;
      write(fileno(stdout), &ch, 1);
    }

    /* Backspace or Delete */

    if ((ch == 8 || ch == 127) && (c > 0)) {  
      write(fileno(stdout), backsp, 3);
      c--;
    }

    /* Completion! */

    if (ch == '\t') { 
      k=c-1;
      buf[c]=0;
      strcpy(compl, buf);
      for (i=0; i<n; i++) {
	if (strlen(cmds[i])> c && !strncmp(buf, cmds[i], c)) {
	  if (k==c-1) { /* First matching command */
	    strcpy(compl, cmds[i]);
	    k=strlen(compl);
	  } else { /* Subsequent matching commands */
	    k = strolof(compl, cmds[i]);
	    compl[k]=0;
	  }
	}
      }
      if (k <= c)
	write(fileno(stdout), &bell, 1);
      else {
	write(fileno(stdout), compl+c, k-c);
	strncpy(buf, compl, k);
	c=k;
      }
    }

    /* CTRL-P, previous command */
    
    if (ch == ('P' & 31)) { 
      new_pos = (cb_pos + SCROLLBACK_SIZE - 1) % SCROLLBACK_SIZE;
      if ((new_pos != SCROLLBACK_SIZE-1 || cmd_buf_full == TRUE) &&
	  new_pos != cmd_buf_pos) {
	cb_pos = new_pos;

	/* Erase current line */

	memset(buf, 8, c); write(fileno(stdout), buf, c);
	memset(buf, 32, c); write(fileno(stdout), buf, c);
	memset(buf, 8, c); write(fileno(stdout), buf, c);

	strcpy(buf, cmd_buf[cb_pos]);
	write(fileno(stdout), cmd_buf[cb_pos], strlen(cmd_buf[cb_pos]));
	c = strlen(cmd_buf[cb_pos]);
      } else
	write(fileno(stdout), &bell, 1);   /* Hit end of cmd_buf */
    }

    /* CTRL-N, next command */

    if (ch == ('N' & 31)) { 
      if (cb_pos != cmd_buf_pos) {
	cb_pos = (cb_pos + 1) % SCROLLBACK_SIZE;

	/* Erase current line */

	memset(buf, 8, c); write(fileno(stdout), buf, c);
	memset(buf, 32, c); write(fileno(stdout), buf, c);
	memset(buf, 8, c); write(fileno(stdout), buf, c);

	strcpy(buf, cmd_buf[cb_pos]);
	write(fileno(stdout), cmd_buf[cb_pos], strlen(cmd_buf[cb_pos]));
	c = strlen(cmd_buf[cb_pos]);
      } else
	write(fileno(stdout), &bell, 1);   /* Hit start of cmd_buf */
    }
      
    /* CTRL-D == cancel (if that is among commands in command array) */

    if (ch == ('D' & 31) && cancel_flag == TRUE) {    
      strcpy(buf, "cancel");
      c=6;
      ch='\n';
    }
    
  } while (ch != '\n');
  
  write(fileno(stdout), &ch, 1);  /* Print the last newline */
  
  buf[c++] = '\n';
  buf[c]   = '\0';

  /* Save command in command buffer */

  if (c != 1) { /* Not empty command -> save in command buffer */
    strncpy(cmd_buf[cmd_buf_pos], buf, c-1);
    cmd_buf[cmd_buf_pos][c-1] = 0;
    cmd_buf_pos = (cmd_buf_pos + 1) % SCROLLBACK_SIZE;
    if (cmd_buf_pos == 0)
      cmd_buf_full = TRUE;
  }
    
  set_raw(RAWOFF);

  /* Special! If first character is '!': escape to shell and
     return nothing to (hide it from) command interpreter */

  if (buf[0] == '!')
    system(buf+1);
} while (buf[0] == '!');
}

void
ogl_wait_for_key(char *prompt) 
{
  unsigned char ch;

  printf("%s", prompt);
  fflush(stdout);

  set_raw(RAWNDELON);
    
  do {
    while (read(fileno(stdin), &ch, 1) == 0)
      ogl_poll();
  } while (ch != '\n' && ch != ' ');

  ch = '\n';

  write(fileno(stdout), &ch, 1);  /* Print newline */

  set_raw(RAWOFF);

}

/********************************** END no-block *******************************/

