//=============================================================================================
//    Motif-OpenGL graphics Interface
//
//  /Purpose:
//    Define a graphics interface using Motif and OpenGL
//  /Author: WDH
//
//  
//============================================================================================

// ******** define this for now *****
#undef noGLwidget
#if !defined(__sgi) || defined(USE_MESA)
  #define noGLwidget
  #define VISUAL     XmNvisual, vi->visual,
#else
  #define VISUAL
#endif
// on the alpha, or 64bit sgi pointers are long int's
#if (defined(__alpha) || (__mips==4))
  #define POINTER_TO_INT long int
#else
  #define POINTER_TO_INT int
#endif
extern "C"
{
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <Xm/MainW.h>
#include <Xm/RowColumn.h>
#include <Xm/PushB.h>
#include <Xm/ToggleB.h>
#include <Xm/CascadeB.h>
#include <Xm/Frame.h>

#include <Xm/Form.h>

#include <X11/StringDefs.h>

#ifdef noGLwidget
#include <Xm/DrawingA.h>	/* Motif drawing area widget */
#else
/** NOTE: in IRIX 5.2, the OpenGL widget headers are mistakenly in   **/
/** <GL/GLwDrawA.h> and <GL/GlwMDraw.h> respectively.  Below are the **/
/** _official_ standard locations.                                   **/
#ifndef __sgi
#ifdef noMotifGLwidget
#include <X11/GLw/GLwDrawA.h> /* STANDARD: pure Xt OpenGL drawing area widget */
#else
#include <X11/GLw/GLwMDrawA.h> /* STANDARD: Motif OpenGL drawing area widget */
#endif
#else
#ifdef noMotifGLwidget
#include <GL/GLwDrawA.h> /* IRIX 5.2: pure Xt OpenGL drawing area widget */
#else
#include <GL/GLwMDrawA.h> /* IRIX 5.2: Motif OpenGL drawing area widget */
#endif
#endif
#endif
#include <X11/keysym.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

/* wdh: */
// *is this needed? * include "tk.h"
#include <Xm/Command.h>
#include <Xm/Text.h>
#include <Xm/LabelG.h>
#include <Xm/ToggleBG.h>
#include <Xm/PushBG.h>
#include <X11/cursorfont.h>
}

#include <assert.h>

#ifndef DOUBLE
typedef float real;
#else
typedef double real;
#endif

static GLenum doubleBuff, directRender, type;


#include "mogl.h"


// *NOTE* we need to add the info from the visual to popup windows --- otherwise this generates an error in XCreateWindow 

#if 0
// *changed for ultra: 
static int snglBuf[] = {GLX_RGBA, GLX_DEPTH_SIZE, 16, None};
static int dblBuf[] = {GLX_RGBA, GLX_DEPTH_SIZE, 16, GLX_DOUBLEBUFFER, None};
static String   fallbackResources[] = {
    "*title: hi", "*plotStuff*width: 350", "*plotStuff*height: 350", NULL
};


#else

// set GLX_RED_SIZE etc. to 0 and use default color map to get a smaller number of colour entries
static int dblBuf[] = {
    GLX_DOUBLEBUFFER, GLX_RGBA, 
    GLX_DEPTH_SIZE, 16,
    GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, 
    None
};
static int *snglBuf = &dblBuf[1];
static String   fallbackResources[] = {
#ifdef IRIX_5_2_or_higher
    "*sgiMode: true",           /* try to enable IRIX 5.2+ look & feel */
    "*useSchemes: all",         /* and SGI schemes */
#endif
    "*title: plotStuff",
    "*plotStuff*width: 350", "*plotStuff*height: 350", NULL
};

#endif

static Display     *dpy;
static GLboolean    doubleBuffer = GL_TRUE;
static XtAppContext app;
// static XtWorkProcId workId = 0;
static Widget       toplevel, mainw, menubar, menupane, btn, cascade, frame, plotStuff, promptText;
static GLXContext   cx;
static XVisualInfo *vi;
#ifdef noGLwidget
static Colormap     cmap;
#endif
static Arg          menuPaneArgs[1], args[1];

static GLboolean postDisplay = GL_FALSE;
static XSetWindowAttributes windowAttibutes;
static Cursor cursor;

static int getCursorPosition;
static real xRubberBandMin, xRubberBandMax, yRubberBandMin, yRubberBandMax;


typedef void MOGL_DISPLAY_FUNCTION(void);

void
clippingPlanesDialog(Widget pb, XtPointer client_data, XtPointer call_data);

static void 
display()
// default display function
{
}


// This points to the function used to display the screen
MOGL_DISPLAY_FUNCTION *displayFunction = display;

void
moglSetDisplayFunction( MOGL_DISPLAY_FUNCTION func )
// Define the display function callback -- called to redraw the screen
{
  displayFunction=func;
}



void 
draw(Widget w)
// draw the OpenGL frame on an expose event
{
  // set cursor to "watch symbol" if it is not already

  int watchCursor = windowAttibutes.cursor==cursor;
  if( !watchCursor )
  {
    windowAttibutes.cursor=cursor;
    XChangeWindowAttributes(dpy,XtWindow(toplevel),CWCursor,&windowAttibutes );
    XFlush(dpy);
  }
  // printf("draw function called, widget name = %s \n",XtName(w));
  displayFunction();
  if (doubleBuffer) glXSwapBuffers(dpy, XtWindow(w));
  if(!glXIsDirect(dpy, cx))
        glFinish(); /* avoid indirect rendering latency from queuing */


  if( !watchCursor )
  {   // set cursor back to normal
    windowAttibutes.cursor=None;
    XChangeWindowAttributes(dpy,XtWindow(toplevel),CWCursor,&windowAttibutes );
    XFlush(dpy);
  }
}

void
moglResetContext(void)
{
  glXMakeCurrent(dpy, XtWindow(plotStuff), cx);
}

void
moglDisplay()
// redraw the OpenGL window
{
  glXMakeCurrent(dpy, XtWindow(plotStuff), cx);
  draw(plotStuff);
}

static XtAppContext & xtAppContext = app;

void 
exposeOrResize(Widget w, XtPointer data, XtPointer callData)
// The screen has been resized or exposed
{
  XmDrawingAreaCallbackStruct *cbs = (XmDrawingAreaCallbackStruct *) callData;
  if( cbs->reason != XmCR_EXPOSE )
  {
    // printf("exposeOrResize called: RESIZE...cbs->reason=%i \n", cbs->reason);
    
    Dimension       width, height;
    XtVaGetValues(w, XmNwidth, &width, XmNheight, &height, NULL);
    // get old values

    GLint viewPort[4];
    glGetIntegerv( GL_VIEWPORT, viewPort );
    // printf("exposeOrResize: old viewPort = %i, %i, %i, %i \n",viewPort[0],viewPort[1],viewPort[2],viewPort[3]);
    
    glViewport(0, 0, (GLint) width, (GLint) height);
    // printf(" new: width=%i, height=%i \n",width,height);
    if( width<viewPort[2] && height<viewPort[3] )
      postDisplay=TRUE;
  }
  else
  {
    // printf("exposeOrResize called: EXPOSE, count = %i\n",cbs->event->xexpose.count);
    // multiple expose events can be queued, wait for a count of zero
    if( cbs->event->xexpose.count>0 )
    {
      int count =cbs->event->xexpose.count;
      // skip contiguous expose events
      XEvent event;
      for( int i=0; i<count; i++ )
        XtAppNextEvent(xtAppContext, &event);   // we wait here when no events are pending
    }
    else
    { 
      // printf("...postDisplay ...\n");
      postDisplay=TRUE;
    }
  }
  
//  XFlush(dpy);
}

void 
map_state_changed(Widget w, XtPointer data, XEvent * event, Boolean * cont)
// This routine called when window is iconified
{
/* ---
    switch (event->type) {
    case MapNotify:
	if (moving && workId != 0) workId = XtAppAddWorkProc(app, animate, NULL);
	break;
    case UnmapNotify:
	if (moving) XtRemoveWorkProc(workId);
	break;
    }
---- */
}


// * wdh
static int exitEventLoop = 0;
static int menuItemChosen=-999;   // number of the menu chosen
static char *menuNameChosen = NULL;     // name of the menu chosen

void
setMenuNameChosen( char* answer )
// Assign the global variable menuNameChosen to equal answer
{
  int length = strlen(answer);
  if( !menuNameChosen || length > strlen(menuNameChosen) )
  {
    delete menuNameChosen;
    menuNameChosen= new char[length+1];
  } 
  strcpy(menuNameChosen,answer);
}  

void
setMenuChosen( const int & menuItem, char* answer )
// Assign the global variable menuNameChosen to equal answer
{
  menuItemChosen = menuItem;
  setMenuNameChosen(answer);
  exitEventLoop=TRUE;
}  


void
exec_cmd( Widget widget, XtPointer client_data, XtPointer call_data )
// This widget is called when a command string is typed in
{
  char *message;
  message = XmTextGetString(widget);
  // printf("command = %s\n",message);
  XmTextSetString(widget,"");       // reset text to blank

  menuItemChosen = -1; // ****
  // printf("item %i chosen, name = %s \n",menuItemChosen,message);
  setMenuNameChosen(message);
  exitEventLoop=TRUE;
}



// *wdh : use this routine to create popup menu
void
inputCommand(Widget cmd_widget, XtPointer client_data, XtPointer call_data )
{
  Widget popup = (Widget) client_data; 
  XmDrawingAreaCallbackStruct *cbs = ( XmDrawingAreaCallbackStruct *) call_data;
  
  if( cbs->event->xany.type != ButtonPress || cbs->event->xbutton.button != 3 )
    return;

  // position the popup menu where the event occured
  XmMenuPosition(popup, (XButtonPressedEvent *) (cbs->event));
  XtManageChild(popup);
}

// *wdh : use this routine to create popup menu
void
postIt(Widget cmd_widget, XtPointer client_data, XEvent *event, char* )
{
  Widget popup = (Widget) client_data; 
  XButtonPressedEvent *bEvent = ( XButtonPressedEvent *) event;
  
  if( bEvent->button != 3 )
    return;

  // position the popup menu where the event occured
  XmMenuPosition(popup, bEvent );
  XtManageChild(popup);
}


// *wdh: Here is which menu item was chosen
void
popupCallback(Widget menuItem, XtPointer client_data, XtPointer call_data )
{
  menuItemChosen = (POINTER_TO_INT) client_data;
  // printf("popupCallback: item %i chosen, name = %s \n",menuItemChosen,XtName(menuItem));
  setMenuNameChosen(XtName(menuItem));
  exitEventLoop=TRUE;
}

/* ----
void
cascadeCallback(Widget menuItem, XtPointer client_data, XtPointer call_data )
{
  menuItemChosen = (POINTER_TO_INT) client_data;
  printf("cascadeCallback: item %i chosen, name = %s \n",menuItemChosen,XtName(menuItem));
  setMenuNameChosen(XtName(menuItem));
  exitEventLoop=TRUE;


}
--- */

// *wdh: call back for the radio buttons
static int toggle_item_set;
void
toggled(Widget widget, XtPointer client_data, XtPointer call_data )
{
  int which = (POINTER_TO_INT) client_data;
  XmToggleButtonCallbackStruct *state = (XmToggleButtonCallbackStruct *)call_data;
  // printf("%s: %s\n",XtName(widget), state->set ? "on" : "off" );
  if( state->set )
    toggle_item_set = which;
  else
    toggle_item_set = 0;
}



// Graphics context:
static GC gc;

//typedef void MOGL_VIEW_FUNCTION( float xa, float xb, float ya, float yb );
MOGL_VIEW_FUNCTION *viewFunction = NULL;

typedef void MOGL_VIEW_FUNCTION(const float & dx,   
				const float & dy , 
				const float & dz,
				const float & dThetaX=0.,
				const float & dThetaY=0.,
				const float & dThetaZ=0.,
				const float & magnify=1. );

void
moglSetViewFunction( MOGL_VIEW_FUNCTION viewFunction_ )
{
  viewFunction=viewFunction_;
}

static int max( int i1, int i2 ){ return i1>i2 ? i1 : i2; }

static int min( int i1, int i2 ){ return i1<i2 ? i1 : i2; }

static float max( float i1, float i2 ){ return i1>i2 ? i1 : i2; }

static float min( float i1, float i2 ){ return i1<i2 ? i1 : i2; }

static double max( double i1, double i2 ){ return i1>i2 ? i1 : i2; }

static double min( double i1, double i2 ){ return i1<i2 ? i1 : i2; }

static double max( double i1, float i2 ){ return i1>i2 ? i1 : i2; }

static double min( double i1, float i2 ){ return i1<i2 ? i1 : i2; }

static int 
rubberBand(Display *disp, Window *winid, XEvent *event0)
// The user can use the left mouse button to make a rubber-band box for zooming
// Return the relative positions of the next bounding box on [0,1]x[0,1]
{
  XEvent event =*event0;
  int x0=event.xbutton.x, y0=event.xbutton.y, x1, y1, x1f, y1f, dummy;
  x1=x0;
  y1=y0;

  Dimension       width, height;
  XtVaGetValues(plotStuff, XmNwidth, &width, XmNheight, &height, NULL);
  assert( width>0 && height>0 );

  xRubberBandMin =min(x0,x1)/float(width);
  xRubberBandMax =max(x0,x1)/float(width);
  yRubberBandMin =1.-max(y0,y1)/float(height);  // make lower left corner (0,0)
  yRubberBandMax =1.-min(y0,y1)/float(height);

  unsigned int mask;

  int first=TRUE;
  // *** define a bitmap for a "+"
  //  0010000   
  //  0010000
  //  1111100
  //  0010000
  //  0010000
  GLubyte cross[5] =
  {
    0x20,
    0x20,
    0xff,
    0x20,
    0x20
  };


  if (event.type == ButtonPress) 
  {
    // printf("rb: button press \n");
    if (!XCheckMaskEvent(disp, ButtonReleaseMask, &event)) 
    {

      glDisable(GL_DEPTH_TEST);  // **** do this !
      glShadeModel(GL_FLAT);

      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      gluOrtho2D(0, width, 0, height);
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();

#ifdef USE_MESA    
      glEnable(GL_BLEND);
      glBlendEquationEXT(GL_LOGIC_OP);
      glBlendEquationEXT(GL_COLOR_LOGIC_OP);
      glPixelStorei(GL_UNPACK_ALIGNMENT,1);   // ** reset?? 
#else
      glEnable(GL_LOGIC_OP);
#endif
      glLogicOp(GL_XOR);
      // glLogicOp(GL_INVERT);

      glDrawBuffer(GL_FRONT);  // draw into front buffer (default is back for double buffering)

      glColor3f(1.,1.,1.);
      int xa,xb,ya,yb;

      // printf("rb: (x0,y0)=(%i,%i)\n",x0,y0);
      x1=x0; y1=y0;
      x1f=x0; y1f=y0;

      while (!XCheckMaskEvent(disp, ButtonReleaseMask, &event)) 
      {
	while (XCheckWindowEvent(disp, *winid, PointerMotionMask, &event)) 
	{
	  if (event.type == MotionNotify) 
	  { 
	    x1 = event.xmotion.x;
	    y1 = event.xmotion.y;
	  }
	}
	if (x1 != x1f || y1 != y1f) 
	{

	  // printf("rb: (x1,y1)=(%i,%i)\n",x1,y1);
	  // Draw "+" at the corners of the rubberband box
	    
#ifdef USE_MESA
	  if( !first )
	  { // erase old marks by drawing with XOR
	    xa=x0, ya=height-y0, xb=x1f, yb=height-y1f;
	    glRasterPos2i(xa,yb);
	    glBitmap(5,5,0.,0.,0.,0.,cross);
	    glRasterPos2i(xb,ya);
	    glBitmap(5,5,0.,0.,0.,0.,cross);
	    glRasterPos2i(xb,yb);
	    glBitmap(5,5,0.,0.,0.,0.,cross);
           
	  }
	  else
	  {
            // mark point  at (x0,y0)
	    glRasterPos2i(x0,height-y0);  
	    glBitmap(5,5,0.,0.,0.,0.,cross);
	    first=FALSE;
	  }
	  xa=x0, ya=height-y0, xb=x1, yb=height-y1;
	  glRasterPos2i(xa,yb);
	  glBitmap(5,5,0.,0.,0.,0.,cross);
	  glRasterPos2i(xb,ya);
	  glBitmap(5,5,0.,0.,0.,0.,cross);
	  glRasterPos2i(xb,yb);
	  glBitmap(5,5,0.,0.,0.,0.,cross);
	  glFlush();

#else
	  xa=x0, ya=height-y0, xb=x1f, yb=height-y1f;
	  glBegin(GL_LINE_LOOP);
	  glVertex2i(xa, ya);
	  glVertex2i(xb, ya);
	  glVertex2i(xb,yb);
	  glVertex2i(xa,yb);
	  glEnd();
	  glFlush();   /* Added by Brian Paul */

	  xa=x0, ya=height-y0, xb=x1, yb=height-y1;
	  glBegin(GL_LINE_LOOP);
	  glVertex2i(xa, ya);
	  glVertex2i(xb, ya);
	  glVertex2i(xb,yb);
	  glVertex2i(xa,yb);
	  glEnd();
	  glFlush();   /* Added by Brian Paul */
#endif
	  x1f=x1; y1f=y1;  // save old values
	}
      }

      //  Erase marks
#ifdef USE_MESA
      xa=x0, ya=height-y0, xb=x1, yb=height-y1;
      glRasterPos2i(xa,ya);
      glBitmap(5,5,0.,0.,0.,0.,cross);
      glRasterPos2i(xa,yb);
      glBitmap(5,5,0.,0.,0.,0.,cross);
      glRasterPos2i(xb,ya);
      glBitmap(5,5,0.,0.,0.,0.,cross);
      glRasterPos2i(xb,yb);
      glBitmap(5,5,0.,0.,0.,0.,cross);
#else
      xa=x0, ya=height-y0, xb=x1, yb=height-y1;
      glBegin(GL_LINE_LOOP);
      glVertex2i(xa, ya);
      glVertex2i(xb, ya);
      glVertex2i(xb,yb);
      glVertex2i(xa,yb);
      glEnd();
      glFlush();   /* Added by Brian Paul */
#endif

      // Make sure the rectangle is big enough
	
      xRubberBandMin =min(x0,x1)/float(width);
      xRubberBandMax =max(x0,x1)/float(width);
      yRubberBandMin =1.-max(y0,y1)/float(height);  // make lower left corner (0,0)
      yRubberBandMax =1.-min(y0,y1)/float(height);

      //  printf("rb: rectangle = [%e,%e]X[%e,%e] \n",xRubberBandMin,xRubberBandMax,yRubberBandMin,yRubberBandMax);

      glDrawBuffer(GL_BACK);
      glDisable(GL_BLEND);
#ifdef USE_MESA
      glDisable(GL_COLOR_LOGIC_OP);
#else
      glDisable(GL_LOGIC_OP);
#endif
      glEnable(GL_DEPTH_TEST);
      glShadeModel(GL_SMOOTH);

      glPopMatrix();
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);

      if( getCursorPosition )
	exitEventLoop=TRUE;

      if( abs(x1-x0)>5 && abs(y1-y0)>5 )
      {
	// call the user routine that will zoom the window:
	// return new window on [0,1]x[0,1] with (0,0) bottom left corner
     
	// if( zoomFunction!=NULL )
	//   zoomFunction( xRubberBandMin,xRubberBandMax,
	//		yRubberBandMin,yRubberBandMax);

        float xShift=1.-(xRubberBandMin+xRubberBandMax);
        float yShift=1.-(yRubberBandMin+yRubberBandMax);
        float magnify=1./max(.005,max(xRubberBandMax-xRubberBandMin,yRubberBandMax-yRubberBandMin));
        // printf("rb: zoom: xShift=%e, yShift=%e, magnify=%e \n",xShift,yShift,magnify);
	if( viewFunction!=NULL )
          viewFunction( xShift,yShift,0.,0.,0.,0.,magnify);
      }
    }
  }

  return FALSE;
}

static int
setView(const int & type, Display *disp, Window *winid, XEvent *event0)
//
//   Mouse Driven rotate/translate/scale function
//
//   type=0 : translate in x/y
//       =1 : rotate about x/y axes
//       =2 : magnify
//       =3 : translate in z
{
  XEvent event =*event0;
  int x0=event.xbutton.x, y0=event.xbutton.y, x1, y1, x1f, y1f;
  x1=x0;
  y1=y0;
  float xShift=0.,yShift=0.,zShift=0.,dThetaX=0.,dThetaY=0.,dThetaZ=0.,magnify=1.;
  
  Dimension       width, height;
  XtVaGetValues(plotStuff, XmNwidth, &width, XmNheight, &height, NULL);
  assert( width>0 && height>0 );

  if (!XCheckMaskEvent(disp, ButtonReleaseMask, &event)) 
  {
    // printf("setView: (x0,y0)=(%i,%i)\n",x0,y0);
    x1=x0; y1=y0;
    x1f=x0; y1f=y0;

    while (!XCheckMaskEvent(disp, ButtonReleaseMask, &event)) 
    {
      while (XCheckWindowEvent(disp, *winid, PointerMotionMask, &event)) 
      {
	if (event.type == MotionNotify) 
	{ 
	  x1 = event.xmotion.x;
	  y1 = event.xmotion.y;
	}
      }
  
      if (x1 != x1f || y1 != y1f) 
      {
        if( type==0 )
	{ // translate
    	  // printf("shift: (x1,y1)=(%i,%i)\n",x1,y1);
          xShift= (x1-x1f)/float(width);
          yShift=-(y1-y1f)/float(height);
	}
	else if( type==1 )
	{ // rotate
    	  // printf("rotate: (x1,y1)=(%i,%i)\n",x1,y1);
          dThetaX= (y1-y1f)/float(height);
          dThetaY= (x1-x1f)/float(width);
	}
	else if( type==2 )
	{
    	  // printf("magnify (x1,y1)=(%i,%i)\n",x1,y1);
          magnify=max(.001,1.+.9*(x1-x1f)/float(width));
	}
	else if( type==3 )
	{
    	  // printf("zShift (x1,y1)=(%i,%i)\n",x1,y1);
          zShift=  (x1-x1f)/float(width);
	}
        if( viewFunction!=NULL )
          viewFunction( xShift,yShift,zShift,dThetaX,dThetaY,dThetaZ,magnify);
        moglDisplay();

	x1f=x1; y1f=y1;  // save old values
      }
    }
  }
  return 0;
}


// *wdh* detect motion and key presses in the drawing area
void
drawAreaInput(Widget w, XtPointer clientData, XtPointer callData)
{
  XmDrawingAreaCallbackStruct *cd = (XmDrawingAreaCallbackStruct *) callData;
  char            buffer[1];
  KeySym          keysym;
  int             rc;
  Position x,y;
  Window win;

  switch (cd->event->type) 
  {
  case KeyRelease:
    /*
     * It is necessary to convert the keycode to a keysym before it is
     * possible to check if it is an escape
     */
    rc = XLookupString((XKeyEvent *) cd->event, buffer, 1, &keysym, NULL);
    switch (keysym) 
    {
    case XK_Up:
      // printf("up-arrow\n");
      break;
    case XK_Down:
      // printf("down-arrow\n");
      break;
    case XK_Left:
      // printf("left-arrow\n");
      break;
    case XK_Right:
      // printf("right-arrow\n");
      break;
    case XK_S: case XK_s: /* the S key */
      // printf("s-key\n");
      break;
    case XK_Shift_L:
      // printf("key release: left shift\n");
      break;
    case XK_Control_L:
      // printf("key release: left control\n");
      break;
    case XK_Escape:
      exit(0);
    }
    break;

/* ---
  case KeyPress:
    rc = XLookupString((XKeyEvent *) cd->event, buffer, 1, &keysym, NULL);
    switch (keysym) 
    {
    case XK_Shift_L:
      printf("key press: left shift\n");
      break;
    case XK_Control_L:
      printf("key press: left control\n");
      break;
    }
---- */

  case ButtonPress:
    
    if( cd->event->xbutton.state & ShiftMask )
    {
      switch (cd->event->xbutton.button)
      {
      case Button1:
        // printf("<shift>button1 \n");
        win =XtWindow(w);
        setView(0,dpy,&win,cd->event); // <shift>button1 = translate
	break;
      case Button2:
        // printf("<shift>button2 \n");
        win =XtWindow(w);
        setView(1,dpy,&win,cd->event); // <shift>button1 = rotate
	break;
      case Button3:
        // printf("<shift>button3 \n");
	break;
      }
    }
    else if( cd->event->xbutton.state & ControlMask )
    {
      switch (cd->event->xbutton.button)
      {
      case Button1:
        // printf("<control>button1 \n");
        win =XtWindow(w);
        setView(2,dpy,&win,cd->event);  // magnify
	break;
      case Button2:
        // printf("<control>button1 \n");
        win =XtWindow(w);
        setView(3,dpy,&win,cd->event);  // shift z
	break;
      }
    }
    else if( cd->event->xbutton.button ==Button1 )
    {
      // printf("button1 press: rubber band zoom \n");
      win =XtWindow(w);
      rubberBand(dpy,&win,cd->event);
    }
    // printf("button press (x,y)=(%i,%i) \n",x,y);
    break;
  case ButtonRelease:
    // x=cd->event->xbutton.x;
    // y=cd->event->xbutton.y;
    // printf("button release (x,y)=(%i,%i) \n",x,y);
    break;
	
  }

}






// *wdh: call back for right column buttons
void
buttonCallback(Widget widget, XtPointer client_data, XtPointer call_data )
{
  // printf("button: %s \n",XtName(widget));

  menuItemChosen = -1; // ****
  // printf("item %i chosen, name = %s \n",menuItemChosen,XtName(widget));
  setMenuNameChosen(XtName(widget));
  exitEventLoop=TRUE;
}

// **** this is not used anymore *****
void 
processAllPendingXEvents() 
{
  if (xtAppContext !=NULL ) 
    while (XtAppPending(xtAppContext) ) 
    {
      XEvent event;
      XtAppNextEvent(xtAppContext, &event);
      XtDispatchEvent(&event);
    }
}

void 
eventLoop()
// Sit in this event loop until some callback sets the global variable exitEventLoop=TRUE
{
  // first set cursor back to normal from a "watch symbol"
  windowAttibutes.cursor=None;
  XChangeWindowAttributes(dpy,XtWindow(toplevel),CWCursor,&windowAttibutes );
  XFlush(dpy);

  exitEventLoop=FALSE;
  while( !exitEventLoop )
  {
// *** processAllPendingXEvents(); // why do this? it just chews up time while idling

    XEvent event;
    XtAppNextEvent(xtAppContext, &event);   // we wait here when no events are pending
/* ---
    printf("event.type = %i \n",event.type);
    if( event.type==Expose )
      printf(" event.xexpose.count =%i \n",event.xexpose.count);
--- */    
    XtDispatchEvent(&event);

    if( postDisplay )
    {
/* ----
      XtAppNextEvent(xtAppContext, &event);   // we wait here when no events are pending
      while( event.type==Expose || event.type==22 || event.type==14)
      {
        printf("eventLoop skipping EXPOSE event...,event.type = %i \n",event.type);
        XtAppNextEvent(xtAppContext, &event);   // we wait here when no events are pending
      }
      printf("event.type = %i \n",event.type);
---- */    
      moglDisplay();
      postDisplay=GL_FALSE;
    }
  }

  // set cursor back to "watch symbol"
  windowAttibutes.cursor=cursor;
  XChangeWindowAttributes(dpy,XtWindow(toplevel),CWCursor,&windowAttibutes );
  XFlush(dpy);
}


int 
getMenuItem( char* &answer )
// Wait for a callback to be called that sets the global variables menuItemChosen and menuNameChosen
//  /answer: The name of the menu chosen
//  /return value: the number of the menu chosen, <0 means a special menu item
{
  eventLoop();
  // printf(" exit event loop, menu chosen = %i, name=%s \n",menuItemChosen,menuNameChosen);
  answer=menuNameChosen;
  return menuItemChosen;
}

void 
getCursor( real & x, real & y )
// Wait for the user to choose a cursor position
{
  getCursorPosition=TRUE;
  eventLoop();
  getCursorPosition=FALSE;
  // printf(" exit event loop, menu chosen = %i, name=%s \n",menuItemChosen,menuNameChosen);
  x=xRubberBandMax;
  y=yRubberBandMax;
}

void 
getRubberBandBoxCorners( real & xMin, real & xMax, real & yMin, real & yMax )
// Wait for the user to choose a cursor position or a rubber band box
{
  // turn off view:
  MOGL_VIEW_FUNCTION *viewSave = viewFunction;
  viewFunction=NULL;

  getCursorPosition=TRUE;
  eventLoop();
  getCursorPosition=FALSE;
  // printf(" exit event loop, menu chosen = %i, name=%s \n",menuItemChosen,menuNameChosen);
  xMin=xRubberBandMin;
  xMax=xRubberBandMax;
  yMin=yRubberBandMin;
  yMax=yRubberBandMax;

  viewFunction=viewSave; // reset
}

XGCValues gcValues;
  

void 
moglInit(int & argc, char *argv[], 
	 const char *windowTitle, 
	 char *fileMenuItems[],
	 ClippingPlaneInfo & clippingPlaneInfo )
//===========================================================================
// /Purpose: Initialize the graphics interface
//
// /fileMenuItems (input): These items appear in the File menu. This array
//   of strings should be terminated by NULL.
//===========================================================================
{
    // title = "*title:"+windowTitle;    
    char *title = new char[strlen(windowTitle)+7+1];
    strcpy(title,"*title:");
    strcat(title,windowTitle);
    fallbackResources[0]=title;

    toplevel = XtAppInitialize(&app, "MotifOpenGLGI", NULL, 0, &argc, argv,
			       fallbackResources, NULL, 0 );

    delete title;  // **** is this ok?? *****

    dpy = XtDisplay(toplevel);

    XSynchronize(dpy,True);  // for debugging

    /* find an OpenGL-capable RGB visual with depth buffer */
    vi = glXChooseVisual(dpy, DefaultScreen(dpy), dblBuf);
    if (vi == NULL) {
	vi = glXChooseVisual(dpy, DefaultScreen(dpy), snglBuf);
	if (vi == NULL)
	    XtAppError(app, "no RGB visual with depth buffer");
	doubleBuffer = GL_FALSE;
    }
    /* create an OpenGL rendering context */
    cx = glXCreateContext(dpy, vi, /* no display list sharing */ None,
        /* favor direct */ GL_TRUE);
    if (cx == NULL)
	XtAppError(app, "could not create rendering context");
    /* create an X colormap since probably not using default visual */
#ifdef noGLwidget

#ifdef USE_MESA
  // use this with Mesa to prevent screen flashing
    // cmap = DefaultColormap(dpy, vi->screen);
     cmap = XCreateColormap(dpy, RootWindow(dpy, vi->screen),   vi->visual, AllocNone);
#else
    cmap = XCreateColormap(dpy, RootWindow(dpy, vi->screen),   vi->visual, AllocNone);
#endif
    // Establish the visual, depth, and colormap of the toplevel
    //  widget _before_ the widget is realized.
    XtVaSetValues(toplevel, XtNvisual, vi->visual, XtNdepth, vi->depth,
                  XtNcolormap, cmap, NULL);
#endif

    // 
    XtAddEventHandler(toplevel, StructureNotifyMask, False,  map_state_changed, NULL);

    mainw = XtVaCreateManagedWidget("main_w", xmMainWindowWidgetClass, toplevel,
               XmNcommandWindowLocation, XmCOMMAND_BELOW_WORKSPACE, NULL);


    // create menu bar
    menubar = XmCreateMenuBar(mainw, "menubar", NULL, 0);
    XtManageChild(menubar);

    // ************** File menu **********************
#ifdef noGLwidget
    /* Hack around Xt's unfortunate default visual inheritance. */
    XtSetArg(menuPaneArgs[0], XmNvisual, vi->visual);
    menupane = XmCreatePulldownMenu(menubar, "menupane", menuPaneArgs, 1);
#else
    menupane = XmCreatePulldownMenu(menubar, "menupane", NULL, 0);
#endif

    for( int i=0; fileMenuItems[i]!=NULL; i++ )
    {
      btn = XmCreatePushButton(menupane,fileMenuItems[i], NULL, 0);
      XtAddCallback(btn, XmNactivateCallback, popupCallback,(void*)(-(i+1)) );
      XtManageChild(btn);
    }
    
    XtSetArg(args[0], XmNsubMenuId, menupane);
    cascade = XmCreateCascadeButton(menubar, "File", args, 1);
    XtManageChild(cascade);

    // *************** view menu ***************
#ifdef noGLwidget
    menupane = XmCreatePulldownMenu(menubar, "menupane", menuPaneArgs, 1);
#else
    menupane = XmCreatePulldownMenu(menubar, "menupane", NULL, 0);
#endif
    btn = XmCreatePushButton(menupane, "clipping planes", NULL, 0);
    XtAddCallback(btn, XmNactivateCallback, clippingPlanesDialog,(void*) (&clippingPlaneInfo) );
    XtManageChild(btn);
    XtSetArg(args[0], XmNsubMenuId, menupane);
    cascade = XmCreateCascadeButton(menubar, "View", args, 1);
    XtManageChild(cascade);

    // *************** help menu ***************
#ifdef noGLwidget
    menupane = XmCreatePulldownMenu(menubar, "menupane", menuPaneArgs, 1);
#else
    menupane = XmCreatePulldownMenu(menubar, "menupane", NULL, 0);
#endif
    btn = XmCreatePushButton(menupane, "Detailed", NULL, 0);
    XtAddCallback(btn, XmNactivateCallback, popupCallback,(void*)(-3) );
    XtManageChild(btn);
    XtSetArg(args[0], XmNsubMenuId, menupane);
    cascade = XmCreateCascadeButton(menubar, "Help", args, 1);
    XtManageChild(cascade);


    // **** create a prompt and command area *****

    Widget bottom = XtVaCreateWidget( "rowcol",xmRowColumnWidgetClass, mainw,XmNorientation, XmVERTICAL,
                                     XmNborderWidth, 0,  NULL);
    // ------ create a text area: prompt: [......]

    Widget promptWindow = XtVaCreateWidget( "subform",xmFormWidgetClass, bottom, 
                NULL);
    Widget prompt =
    XtVaCreateManagedWidget("Prompt:",xmLabelGadgetClass, promptWindow, 
                            XmNtopAttachment,    XmATTACH_FORM,
                            XmNbottomAttachment, XmATTACH_FORM,
                            XmNleftAttachment,   XmATTACH_FORM,
                            XmNalignment,        XmALIGNMENT_BEGINNING,
			    NULL );

    promptText = XtVaCreateManagedWidget( "promptWindow",xmTextWidgetClass, promptWindow,
                            XmNleftAttachment, XmATTACH_WIDGET,
                            XmNrightAttachment, XmATTACH_FORM,
                            XmNleftWidget,      prompt,             
                            NULL); 
    XtManageChild(promptWindow);
    
    // ------ create a text area: command: [......]

    Widget command_w = XtVaCreateWidget( "subform",xmFormWidgetClass, bottom, 
                NULL);
    Widget commandText =
    XtVaCreateManagedWidget("Command:",xmLabelGadgetClass, command_w, 
                            XmNtopAttachment,    XmATTACH_FORM,
                            XmNbottomAttachment, XmATTACH_FORM,
                            XmNleftAttachment,   XmATTACH_FORM,
                            XmNalignment,        XmALIGNMENT_BEGINNING,
			    NULL );

    Widget command =  XtVaCreateManagedWidget( "commandText",xmTextWidgetClass, command_w,
                            XmNleftAttachment, XmATTACH_WIDGET,
                            XmNrightAttachment, XmATTACH_FORM,
                            XmNleftWidget,      commandText,             
                            NULL); 





    XtAddCallback(command,XmNactivateCallback, exec_cmd, NULL );
    XtManageChild(command_w);
    XtManageChild(bottom);


    // rightColumn holds all the buttons etc. down the right side
    Widget rightColumn= XtVaCreateWidget( "rightColumn",xmRowColumnWidgetClass, mainw, 
                XmNorientation, XmVERTICAL,
                XmNpacking, XmPACK_TIGHT,   // default
                 // ? XmNentryVerticalAlignment, XmALIGNMENT_CENTER,
                NULL);

    // make rotation buttons, put in a frame
    Widget rotateFrame = XtVaCreateWidget("rotateFrame",xmFrameWidgetClass, rightColumn,
					  XmNshadowType, XmSHADOW_ETCHED_IN, NULL );
    Widget rotateButtons= XtVaCreateWidget( "rotateButtons",xmRowColumnWidgetClass, rotateFrame, 
                XmNorientation, XmHORIZONTAL,  
                XmNpacking, XmPACK_COLUMN,
                XmNnumColumns, 3,
                NULL);
    Widget xpr = XtVaCreateManagedWidget("x+r",xmPushButtonWidgetClass, rotateButtons,NULL);    
    Widget xmr = XtVaCreateManagedWidget("x-r",xmPushButtonWidgetClass, rotateButtons,NULL);    
    Widget ypr = XtVaCreateManagedWidget("y+r",xmPushButtonWidgetClass, rotateButtons,NULL);    
    Widget ymr = XtVaCreateManagedWidget("y-r",xmPushButtonWidgetClass, rotateButtons,NULL);    
    Widget zpr = XtVaCreateManagedWidget("z+r",xmPushButtonWidgetClass, rotateButtons,NULL);    
    Widget zmr = XtVaCreateManagedWidget("z-r",xmPushButtonWidgetClass, rotateButtons,NULL);    
    XtAddCallback(xpr, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(xmr, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(ypr, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(ymr, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(zpr, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(zmr, XmNactivateCallback, buttonCallback, NULL);
    

    XtManageChild(rotateButtons);
    XtManageChild(rotateFrame);

    // make shift buttons, put in a frame
    Widget shiftFrame = XtVaCreateManagedWidget("shiftFrame",xmFrameWidgetClass, rightColumn,
					  XmNshadowType, XmSHADOW_ETCHED_IN, NULL );
    Widget shiftButtons= XtVaCreateWidget( "shiftButtons",xmRowColumnWidgetClass, shiftFrame, 
            //    XmNorientation, XmVERTICAL,  
            //    XmNpacking, XmPACK_COLUMN,
            //    XmNnumColumns, 2,
                XmNorientation, XmHORIZONTAL,  
                XmNpacking, XmPACK_COLUMN,
                XmNnumColumns, 3,
                NULL);
    Widget xp = XtVaCreateManagedWidget("x+",xmPushButtonWidgetClass, shiftButtons,NULL);    
    Widget xm = XtVaCreateManagedWidget("x-",xmPushButtonWidgetClass, shiftButtons,NULL);    
    Widget yp = XtVaCreateManagedWidget("y+",xmPushButtonWidgetClass, shiftButtons,NULL);    
    Widget ym = XtVaCreateManagedWidget("y-",xmPushButtonWidgetClass, shiftButtons,NULL);    
    Widget zp = XtVaCreateManagedWidget("z+",xmPushButtonWidgetClass, shiftButtons,NULL);    
    Widget zm = XtVaCreateManagedWidget("z-",xmPushButtonWidgetClass, shiftButtons,NULL);    
    XtAddCallback(xp, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(xm, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(yp, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(ym, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(zp, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(zm, XmNactivateCallback, buttonCallback, NULL);
    XtManageChild(shiftButtons);

    // make bigger/smaller buttons, put in a frame
    Widget bigSmallFrame = XtVaCreateManagedWidget("bigSmallFrame",xmFrameWidgetClass, rightColumn,
					  XmNshadowType, XmSHADOW_ETCHED_IN, NULL );
    Widget bigSmallButtons= XtVaCreateWidget( "bigSmallButtons",xmRowColumnWidgetClass, bigSmallFrame, 
                // XmNorientation, XmVERTICAL,  
                XmNorientation, XmHORIZONTAL,  
                XmNpacking, XmPACK_COLUMN,
                XmNnumColumns, 4,
                NULL);
    Widget bigger = XtVaCreateManagedWidget("bigger",xmPushButtonWidgetClass, bigSmallButtons,NULL);    
    Widget smaller = XtVaCreateManagedWidget("smaller",xmPushButtonWidgetClass, bigSmallButtons,NULL);    
    Widget reset = XtVaCreateManagedWidget("reset",xmPushButtonWidgetClass, bigSmallButtons,NULL);    
    Widget erase = XtVaCreateManagedWidget("erase",xmPushButtonWidgetClass, bigSmallButtons,NULL);    
    XtAddCallback(bigger, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(smaller, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(reset, XmNactivateCallback, buttonCallback, NULL);
    XtAddCallback(erase, XmNactivateCallback, buttonCallback, NULL);
    XtManageChild(bigSmallButtons);

    XtManageChild(rightColumn);

    // *********** create framed drawing area for OpenGL rendering ******************

    frame = XmCreateFrame(mainw, "frame", NULL, 0);
    XtManageChild(frame);
#ifdef noGLwidget
//    plotStuff = XtVaCreateManagedWidget("plotStuff", xmDrawingAreaWidgetClass,frame, XtNborderWidth, 0, NULL); // ***

    plotStuff = XtVaCreateManagedWidget("plotStuff", xmDrawingAreaWidgetClass,frame,  
             XtNborderWidth, 0, XmNvisual, vi->visual, NULL); // ***


//    plotStuff = XtVaCreateManagedWidget("plotStuff", xmDrawingAreaWidgetClass,frame,  
//             XmNborderWidth, 0, XmNdepth, vi->depth, XmNcolormap, cmap, NULL); // ***
#else
#ifdef noMotifGLwidget
    /* notice glwDrawingAreaWidgetClass lacks an 'M' */
    plotStuff = XtVaCreateManagedWidget("plotStuff", glwDrawingAreaWidgetClass,frame, GLwNvisualInfo, vi, NULL);
#else
    plotStuff = XtVaCreateManagedWidget("plotStuff", glwMDrawingAreaWidgetClass,frame, GLwNvisualInfo, vi, NULL);
#endif
#endif

  // *wdh create a graphics context for rubber band:
//   gcValues.foreground=BlackPixel(XtDisplay(plotStuff),0);
//   gcValues.background=WhitePixel(XtDisplay(plotStuff),0);
//   gc = XCreateGC( XtDisplay(plotStuff),RootWindowOfScreen(XtScreen(plotStuff)),
//         (GCForeground | GCBackground),&gcValues);

    gc = XCreateGC( dpy,RootWindowOfScreen(XtScreen(plotStuff)),NULL,NULL);
//  gc = XDefaultGC(dpy,0);

/* ---
    Window root_return;
    int x_return, y_return;
    unsigned int width_return, height_return, bw_return, depth_return;
    
    XGetGeometry(dpy,RootWindowOfScreen(XtScreen(plotStuff)),&root_return,&x_return,&y_return,
		 &width_return, &height_return, &bw_return, &depth_return);
    printf("root_return=%i, width_return=%i, height_return=%i, bw_return=%i, depth_return=%i \n",
	   root_return, width_return, height_return, bw_return, depth_return) ;
----- */

  // XtAddCallback(plotStuff, XmNexposeCallback, (XtCallbackProc)draw, NULL);
  XtAddCallback(plotStuff, XmNexposeCallback, exposeOrResize, NULL);
  XtAddCallback(plotStuff, XmNresizeCallback, exposeOrResize, NULL);
  XtAddCallback(plotStuff, XmNinputCallback, drawAreaInput, NULL);

  // set up application's window layout
  // XmMainWindowSetAreas(mainw, menubar, NULL, NULL, rightColumn, frame);
  XmMainWindowSetAreas(mainw, menubar, bottom, NULL, rightColumn, frame);



  XtRealizeWidget(toplevel);
    /*
     * Once widget is realized (ie, associated with a created X window), we
     * can bind the OpenGL rendering context to the window.
     */
    glXMakeCurrent(dpy, XtWindow(plotStuff), cx);
    // setup OpenGL state  **** should we do this ??
    // glClearDepth(1.0);
    // glClearColor(0.0, 0.0, 0.0, 0.0);

  cursor=XCreateFontCursor(dpy,XC_watch);
  windowAttibutes.cursor=cursor;
  XChangeWindowAttributes(dpy,XtWindow(toplevel),CWCursor,&windowAttibutes );
  XFlush(dpy);

}

void 
moglSetPrompt(char *prompt)
// ===================================================================================
// /Purpose:
//    Set the prompt.
// ===================================================================================
{
  if( prompt )
    XmTextSetString(promptText,prompt);       // set prompt
}

int 
moglGetMenuItem( char *menu[], char* &answer, char *prompt )
// ===================================================================================
// /Purpose:
//  Get a menu entry from an array of strings, terminated with a null string
//  To create a cascading menu begin the string with an '>'
//  To end the cascade begin the string with an '<' 
//  To end a cascade and start a new cascade begin the string with '<' followed by '>'
//  Here is an example:
//  char *menu1[] = {  "plot",
//                     ">component",
//                                  "u",
//                                  "v",
//                                  "w",
//                     "<erase",
//                     ">stuff",
//                              "s1",
//                              ">more stuff", 
//                                            "more1",
//                                            "more2", 
//                              "<s2", 
//                    "<>apples", 
//                              "apple1", 
//                    "<exit",NULL };  
//
//
// ===================================================================================
{
  Widget & popWidget = mainw; // The popup can be activated anywhere in this window.

//  Widget popupMenu = XmVaCreateSimplePopupMenu(popWidget,"popup",popupCallback,NULL);
  // *NOTE* we need to add the info from the visual --- otherwise this generates an error in XCreateWindow 
  Widget popupMenu = XmVaCreateSimplePopupMenu(mainw,"popup",popupCallback,VISUAL NULL);  // ok

  const int maxNumberOfCascadeLevels=10;
  Widget widget[maxNumberOfCascadeLevels], gidget;
  int level=0;
  int count[maxNumberOfCascadeLevels];
  for( int i=0; i<maxNumberOfCascadeLevels; i++ )
  {
    widget[i]=popupMenu;
    count[i]=-1;
  }
  
  int numberOfMenuEntries=0;
  for( i=0; menu[i]!=NULL && i<200; i++ )  // at most 200 entries !
  {
    // printf("menu[%i]=%s \n",i,menu[i]);
    numberOfMenuEntries++;
    if( menu[i][0]=='>' )
    {
      gidget = XtVaCreateManagedWidget(menu[i]+1,xmCascadeButtonWidgetClass,widget[level],NULL);
      count[level]++;
      level++;
      if( level>=maxNumberOfCascadeLevels )
	printf("moglGetMenuItem: Too many levels of cascading menus, max= %i \n",maxNumberOfCascadeLevels);

      // printf("add Pulldown, level=%i, count[l-1]=%i \n",level,count[level-1]);
      widget[level] = XmVaCreateSimplePulldownMenu(widget[level-1],"pullright",count[level-1],popupCallback,
                        VISUAL NULL );
    }
    else if( menu[i][0]=='<' )
    {
      count[level]=-1;
      level--;
      if( strlen(menu[i])>1 && menu[i][1]=='>' )
      { // This case occurs if the menu starts with "<>..." -- when one cascade menu follows another
	gidget = XtVaCreateManagedWidget(menu[i]+2,xmCascadeButtonWidgetClass,widget[level],NULL);
	count[level]++;
	level++;
	widget[level] = XmVaCreateSimplePulldownMenu(widget[level-1],"pullright",count[level-1],popupCallback,
                            VISUAL NULL );
      }
      else
      {
        gidget = XtVaCreateManagedWidget(menu[i]+1,xmPushButtonGadgetClass,widget[level],NULL);
        count[level]++;
      }
    }
    else
    { 
      gidget = XtVaCreateManagedWidget(menu[i],xmPushButtonGadgetClass,widget[level],NULL);
      count[level]++;
    }
    XtAddCallback(gidget,XmNactivateCallback,popupCallback,(void*)i);
  }

  if( numberOfMenuEntries >= 199 )
  {
    printf("Too many menu entries! \n");
    exit(1);  
  }
  XtAddEventHandler(popWidget, ButtonPressMask, FALSE, postIt, popupMenu );

  if( prompt )
    XmTextSetString(promptText,prompt);       // set prompt

  menuItemChosen=getMenuItem(answer);

  if( prompt )
    XmTextSetString(promptText,"");       // reset prompt to blank

  XtRemoveEventHandler(popWidget, ButtonPressMask, FALSE, postIt, popupMenu );
  return menuItemChosen;
}

int 
moglGetMenuItem2( char *menu[], char* &answer, char *prompt=NULL )
//
// /Purpose:
//  Get a menu entry from an array of strings, terminated with a null string
//
{
  Widget & popWidget = mainw; // The popup can be activated anywhere in this window.

  Widget popupMenu = XmVaCreateSimplePopupMenu(popWidget,"popup",popupCallback,NULL);

  int numberOfMenuEntries=0;
  for( int i=0; menu[i]!=NULL && i<200; i++ )  // at most 200 entries !
  {
    // printf("menu[%i]=%s \n",i,menu[i]);
    numberOfMenuEntries++;
    Widget widget = XtVaCreateManagedWidget(menu[i],xmPushButtonGadgetClass,popupMenu,NULL);
    XtAddCallback(widget,XmNactivateCallback,popupCallback,(void*)i);
  }
  if( numberOfMenuEntries >= 199 )
  {
    printf("Too many menu entries! \n");
    exit(1);  
  }
  XtAddEventHandler(popWidget, ButtonPressMask, FALSE, postIt, popupMenu );

  if( prompt )
    XmTextSetString(promptText,prompt);       // set prompt

  menuItemChosen=getMenuItem(answer);

  if( prompt )
    XmTextSetString(promptText,"");       // reset prompt to blank

  XtRemoveEventHandler(popWidget, ButtonPressMask, FALSE, postIt, popupMenu );
  return menuItemChosen;
}
    


void 
moglGetWindowSize( int & width, int & height )
{
  Dimension width0, height0;
  XtVaGetValues(plotStuff, XmNwidth, &width0, XmNheight, &height0, NULL);
  width=width0;
  height=height0;  
}

void 
moglPostDisplay()
{
  postDisplay=TRUE;
}


// ********************** clip ***************************

#include <Xm/SelectioB.h>
#include <Xm/MessageB.h>
#include <Xm/RowColumn.h>
#include <Xm/PushB.h>

#include <Xm/ToggleB.h>
#include <Xm/Text.h>
#include <Xm/TextF.h>
#include <Xm/LabelG.h>
#include <Xm/Frame.h>

#include <Xm/Scale.h>
#include <Xm/DialogS.h>
#include <stdio.h>

//void
//setMenuChosen( const int & menuItem, char* answer );

static ClippingPlaneInfo *clippingPlaneInfo;

#define clipPlane(plane,i) *(clippingPlaneInfo->clippingPlaneEquation+(i) \
                            +(plane)*4)

#define clipPlaneIsOn(plane) *(clippingPlaneInfo->clippingPlaneIsOn+(plane))

void 
destroyWidget(_WidgetRec* w, void* a, void* b)
{
  XtDestroyWidget(w);
}

static char buff[80];


void
clipPlaneOnOff( Widget widget, XtPointer clinet_data, XtPointer call_data )
{
  XmToggleButtonCallbackStruct *state = (XmToggleButtonCallbackStruct *) call_data;
//  printf("%s: %s\n", XtName(widget), state->set ? "on" : "off" );

  sprintf(buff,"%s %s",XtName(widget),state->set ? "on" : "off" );
//  printf("menu: %s \n",buff);
  setMenuChosen(-99,buff);
}


void
scaleCallback(Widget widget, XtPointer clinet_data, XtPointer call_data )
{
  XmScaleCallbackStruct *cbs = (XmScaleCallbackStruct *) call_data;
// printf("%s: %d\n",XtName(widget), cbs->value);

  sprintf(buff,"%s %e",XtName(widget),cbs->value/100.);
  // printf("menu: %s \n",buff);
  setMenuChosen(-99,buff);
}

void
changeNormal(Widget widget, XtPointer clinet_data, XtPointer call_data )
// This widget is called when a normal value is changed
{
  char *message;
  message = XmTextGetString(widget);
//XmTextSetString(widget,"");       // reset text to blank

//  printf("change %s: new value = %s\n",XtName(widget),message);
  
  float n0=0.,n1=0.,n2=-1.;
  printf("changeNormal message: %s \n",message);
  sscanf( message,"%e,%e,%e",&n0,&n1,&n2 ); 
  // printf("set normal=(%e,%e,%e) for plane %i \n",n0,n1,n2,plane);
  float norm=pow(n0*n0+n1*n1+n2*n2,.5);
  if( norm == 0. )
  {
    printf("changeNormal: invalid clipping plane normal, cannot be all zeroes! \n");
    n2=-1;
  }
  sprintf(buff,"%3.2f, %3.2f, %3.2f",n0/norm,n1/norm,n2/norm);
  XmTextSetString(widget,buff);

  sprintf(buff,"%s %s",XtName(widget),message);
  // printf("menu: %s \n",buff);
  setMenuChosen(-99,buff);


}

void
positionDialog(Widget dialog, XtPointer clinet_data, XtPointer call_data )
{
  Position x,y;
  Dimension w,h;
  XtVaGetValues(dialog, XmNwidth, &w, XmNheight, &h, NULL );
  
  x=WidthOfScreen(XtScreen(dialog))-w;
  x=500;
//  y=HeightOfScreen(XtScreen(dialog))-h;
  y=0;
  printf(" position dialog: w=%i, h=%i, x=%i, y=%i \n",w,h,x,y);
  
  XtVaSetValues(dialog,XmNx,x,XmNy,y, NULL );
}


/* ---
double clipPlane[6][4] = { 0.,0.,1.,0.,
			   0.,1.,0.,0.,
			   0.,1.,0.,0.,
			   0.,1.,0.,0.,
			   0.,1.,0.,0.,
			   0.,1.,0.,0.};

---- */

void
clippingPlanesDialog(Widget pb, XtPointer client_data, XtPointer call_data)
{
  Widget dialog;
  XmString t;
  Arg args[7];
  int n = 0;

  // retrieve the pointer to the clipping plane info
  clippingPlaneInfo=(ClippingPlaneInfo *)client_data;
  //printf(" clippingPlaneInfo: %i, (%e,%e,%e) \n",*clippingPlaneInfo->maximumNumberOfClippingPlanes,
  //        clipPlane(0,0),clipPlane(0,1),clipPlane(0,2),clipPlane(0,3));

    // Create the dialog -- the PushButton acts as the DialogShell's
    // parent (not the parent of the PromptDialog).  The "userData"
    // is used to store the value 

    t = XmStringCreateLocalized ("Cancel");
  XtSetArg (args[n], XmNcancelLabelString, t); n++;
  XtSetArg (args[n], XmNautoUnmanage, False); n++;
  XtSetArg (args[n], XtNtitle, "clipping planes"); n++;
#ifdef noGLwidget
  XtSetArg (args[n], XmNvisual, vi->visual ); n++;
#endif
//XtSetArg (args[n], XmNy, 200); n++;

  dialog = XmCreateTemplateDialog (pb, "clipping planes", 
				   args,n );

  XmStringFree (t); /* always destroy compound strings when done */

  XtAddCallback (dialog, XmNcancelCallback, destroyWidget, NULL);

// these did not work:

//  XtAddCallback (dialog, XmNmapCallback, positionDialog, NULL);
//  XtAddCallback (dialog, XmNpopupCallback, positionDialog, NULL);
	
//  XtVaSetValues(dialog,XmNx,200,XmNy,200, NULL );

    // arg list for scales
    XtVarArgsList arglist = XtVaCreateArgsList (NULL,
        XmNshowValue, True,
        XmNminimum, -175,
        XmNmaximum, 175,
        XmNscaleMultiple, 1,
	XmNorientation,XmHORIZONTAL,
	XmNvalue,0,
	XmNdecimalPoints, 2,
        NULL);

    // this rowColumn holds all the clipping planes
    Widget baseRC = XtVaCreateWidget ("rowcol", xmRowColumnWidgetClass, dialog,
				      XmNorientation, XmVERTICAL,
				      NULL);

    // **************************************************
    // *********  Make some cipping planes **************
    // **************************************************
    Widget rowcol, toggle, scale, normalRowCol, normal;
    for( int clip=0; clip<3; clip++ )
    {
      Widget frame0 = XtVaCreateManagedWidget("frame",
					      xmFrameWidgetClass, baseRC,
					      XmNshadowType, XmSHADOW_ETCHED_IN, NULL );

      rowcol = XtVaCreateWidget ("rowcol", xmRowColumnWidgetClass, frame0,
				 XmNorientation, XmVERTICAL,
				 NULL);

      // Create a toggle button 
      sprintf(buff,"clip %i",clip);
      toggle = XtVaCreateManagedWidget(
        buff,
	xmToggleButtonWidgetClass, rowcol,
        XmNindicatorType, XmN_OF_MANY,
        NULL);
      XmToggleButtonSetState(toggle,clipPlaneIsOn(clip),0);
      XtAddCallback(toggle, XmNvalueChangedCallback, clipPlaneOnOff, NULL  );


      sprintf(buff,"clip %i distance",clip);
      scale = XtVaCreateManagedWidget(buff,
				      xmScaleWidgetClass, rowcol,
				      XtVaNestedList, arglist,
				      XtVaTypedArg, XmNtitleString, XmRString, "distance", 8,
				      NULL);
      XtAddCallback(scale,XmNvalueChangedCallback, scaleCallback, NULL);
      // set scale value for the distance
      XmScaleSetValue(scale,(int)(clipPlane(clip,3)*100.+.5));

      // add ability to input the normal
      normalRowCol = XtVaCreateWidget("normalRowCol",
				       xmRowColumnWidgetClass, rowcol,
				       XmNorientation, XmHORIZONTAL, NULL);
      XtVaCreateManagedWidget("normal:", xmLabelGadgetClass, normalRowCol, NULL);

      sprintf(buff,"clip %i normal",clip);
      normal = XtVaCreateManagedWidget(buff,xmTextFieldWidgetClass, 
				       normalRowCol, XmNcolumns, 20, XmNmarginHeight, 2, NULL );
      XtAddCallback(normal, XmNactivateCallback, changeNormal, NULL );
      sprintf(buff,"%3.2f, %3.2f, %3.2f",clipPlane(clip,0),clipPlane(clip,1),clipPlane(clip,2));
      XmTextSetString(normal,buff);

/* ---
      for( int dir=0; dir<3; dir++ )
      {
        sprintf(buff,"normal%i[%i]",clip,dir);
        normal = XtVaCreateManagedWidget(buff,xmTextFieldWidgetClass, 
					       normalRowCol, XmNcolumns, 4, XmNmarginHeight, 2, NULL );
        XtAddCallback(normal, XmNactivateCallback, changeNormal, NULL );
        sprintf(buff,"%3.2f",clipPlane[clip][dir]);
        XmTextSetString(normal,buff);
      }
---- */

      XtManageChild(normalRowCol);
      XtManageChild (rowcol);
      
    }

    XtManageChild (baseRC);

    XtManageChild (dialog);
    XtPopup (XtParent (dialog), XtGrabNone);
}
