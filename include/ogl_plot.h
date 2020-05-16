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
#ifndef ogl_plot_h
#define ogl_plot_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <termio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>		/* this includes X and gl.h headers */

#include "real.h"
#include "input_output.h"

/*
** Window Types
*/

#define TK_RGB		0
#define TK_INDEX	1
#define TK_SINGLE	0
#define TK_DOUBLE	2
#define TK_DIRECT	0
#define TK_INDIRECT	4
#define TK_ACCUM	8
#define TK_ALPHA	16
#define TK_DEPTH	32
#define TK_OVERLAY	64
#define TK_UNDERLAY	128
#define TK_STENCIL	512

/*
** Display Mode Selection Criteria
*/

enum {
    TK_USE_ID = 1,
    TK_EXACT_MATCH,
    TK_MINIMUM_CRITERIA
};

/* 
** Window Masks
*/

#define TK_IS_RGB(x)		(((x) & TK_INDEX) == 0)
#define TK_IS_INDEX(x)		(((x) & TK_INDEX) != 0)
#define TK_IS_SINGLE(x)		(((x) & TK_DOUBLE) == 0)
#define TK_IS_DOUBLE(x)		(((x) & TK_DOUBLE) != 0)

/* There was no bug in TK_IS_DIRECT !*/
#define TK_IS_DIRECT(x)		(((x) & TK_INDIRECT) == 0)
#define TK_IS_INDIRECT(x)	(((x) & TK_INDIRECT) != 0)

#define TK_HAS_ACCUM(x)		(((x) & TK_ACCUM) != 0)
#define TK_HAS_ALPHA(x)		(((x) & TK_ALPHA) != 0)
#define TK_HAS_DEPTH(x)		(((x) & TK_DEPTH) != 0)
#define TK_HAS_OVERLAY(x)	(((x) & TK_OVERLAY) != 0)
#define TK_HAS_UNDERLAY(x)	(((x) & TK_UNDERLAY) != 0)
#define TK_HAS_STENCIL(x)	(((x) & TK_STENCIL) != 0)

/*
** Windowing System Specific Gets
*/

enum {
    TK_X_DISPLAY = 1,
    TK_X_WINDOW,
    TK_X_SCREEN,
    TK_CONTEXT
};

/*
** Event Status
*/

#define	TK_LEFTBUTTON		1
#define	TK_RIGHTBUTTON		2
#define	TK_MIDDLEBUTTON		4
#define	TK_SHIFT		1
#define	TK_CONTROL		2

/* 
** Key Codes
*/

#define TK_RETURN		0x0D
#define TK_ESCAPE		0x1B
#define TK_SPACE		0x20
#define TK_LEFT			0x25
#define TK_UP			0x26
#define TK_RIGHT		0x27
#define TK_DOWN			0x28
#define TK_A			'A'
#define TK_B			'B'
#define TK_C			'C'
#define TK_D			'D'
#define TK_E			'E'
#define TK_F			'F'
#define TK_G			'G'
#define TK_H			'H'
#define TK_I			'I'
#define TK_J			'J'
#define TK_K			'K'
#define TK_L			'L'
#define TK_M			'M'
#define TK_N			'N'
#define TK_O			'O'
#define TK_P			'P'
#define TK_Q			'Q'
#define TK_R			'R'
#define TK_S			'S'
#define TK_T			'T'
#define TK_U			'U'
#define TK_V			'V'
#define TK_W			'W'
#define TK_X			'X'
#define TK_Y			'Y'
#define TK_Z			'Z'
#define TK_a			'a'
#define TK_b			'b'
#define TK_c			'c'
#define TK_d			'd'
#define TK_e			'e'
#define TK_f			'f'
#define TK_g			'g'
#define TK_h			'h'
#define TK_i			'i'
#define TK_j			'j'
#define TK_k			'k'
#define TK_l			'l'
#define TK_m			'm'
#define TK_n			'n'
#define TK_o			'o'
#define TK_p			'p'
#define TK_q			'q'
#define TK_r			'r'
#define TK_s			's'
#define TK_t			't'
#define TK_u			'u'
#define TK_v			'v'
#define TK_w			'w'
#define TK_x			'x'
#define TK_y			'y'
#define TK_z			'z'
#define TK_0			'0'
#define TK_1			'1'
#define TK_2			'2'
#define TK_3			'3'
#define TK_4			'4'
#define TK_5			'5'
#define TK_6			'6'
#define TK_7			'7'
#define TK_8			'8'
#define TK_9			'9'

/*
** Color Macros
*/

enum {
    OGL_WHITE = 0,
    OGL_RED,
    OGL_GREEN,
    OGL_YELLOW,
    OGL_BLUE,
    OGL_MAGENTA,
    OGL_CYAN,
    OGL_BLACK,
    OGL_GOLD
};

extern float tkRGBMap[8][3];

#define TK_SETCOLOR(x, y) (TK_IS_RGB((x)) ? \
		           glColor3fv(tkRGBMap[(y)]) : glIndexf((y)))

/*
** RGB Image Structure
*/

typedef struct _TK_RGBImageRec {
    GLint sizeX, sizeY;
    unsigned char *data;
} TK_RGBImageRec;


#define OGL_ACTIVE_LISTS 100
#define OGL_NEW_PLOT 1
#define OGL_OLD_PLOT 0

typedef struct {
  GLdouble vdir[3], up[3], right[3], focus[3], distance;
} ogl_view_state;

/*
** Window manager structure
*/

typedef struct {
  int x, y, w, h;
  double scale;
  GLenum type;
  GLenum dmPolicy;
  Window wMain;
  XVisualInfo *vInfoMain;
  Colormap cMapMain;
  GLXContext cMain;
  GLuint bitmapBase;
  GLuint current_display_list[OGL_ACTIVE_LISTS];
  int n_active_lists;
  GLboolean need_redraw;
  ogl_view_state view_state, initial_view_state;
  XFontStruct *xFont;
  int use_scaling;
  double x_scale, y_scale, z_scale;
} WINDOW_REC;

/* structure to define up to four clip planes */
typedef struct{
  int n_clip_planes;
  GLdouble clip_plane[4][3];
} clip_planes;

/* Defines for X inter routines */

#define MAX_VERTICES               100
#define MAX_WINDOWS                20
#define MAX_MOUSE_BUTTONS          5
#define MAX_POPUP_ENTRIES          20
#define RESET                      2
#define TRUE                       1
#define FALSE                      0
#define OK                         TRUE
#define ERROR                      FALSE
#define BORDER_WIDTH               3
#define LINE                       0
#define POLYGON                    1
#define GENERIC                    2
#define POPUP_ENTRY_HEIGHT         20
#define TEXT_SAVE_BUF_SIZE         50
#define LINE_SAVE_BUF_SIZE         500
#define MARKER_SAVE_BUF_SIZE       100
#define COLOR_SAVE_BUF_SIZE        16

#define ConvertCoordX(x,xmin,xmax,pixels) (((pixels)*((x)-(xmin))/((xmax)-(xmin)))+0.5)
#define ConvertCoordY(y,ymin,ymax,pixels) (((pixels)*((ymax)-(y))/((ymax)-(ymin)))+0.5)

/* Types for X inter */
					     
typedef struct _textsave {
  int hjust, vjust;
  real x, y;
  char *s;
  unsigned int cind;
} TextSaveEntry;

typedef struct _linesave {
  real x0, y0, x1, y1;
  unsigned int cind;
} LineSaveEntry;

typedef struct _markersave {
  real x, y;
  int type;
  unsigned int cind;
} MarkerSaveEntry;

typedef Window winid_type;

typedef struct Windata {
  winid_type window;
  GC gc;
  int col_index;
  int width, height;
  int xpadl, xpadr, ypadt, ypadb;
  real xmin, xmax, ymin, ymax;
  real sxmin, sxmax, symin, symax;
  real aspect;
  real cx, cy;
  TextSaveEntry *ts;
  int ts_count;
  int ts_size;
  LineSaveEntry *ls;
  int ls_count;
  int ls_size;
  MarkerSaveEntry *ms;
  int ms_count;
  int ms_size;
} win_data_type;

typedef struct Entry {
  char *name;
  int ret_val;
  /*	  popup_type *next_popup;*/
} popup_entry;

typedef struct Popup {
  int n_entries;
  popup_entry *entry;
} popup_type;

typedef struct {
  double x_min, x_max;
  double y_min, y_max;
  double z_min, z_max;
} bounding_box;

/* Function prototypes for ogl functions */
void
ogl_init_prompt(input_output *io_ptr);
void
ogl_close_prompt(input_output *io_ptr);
int 
ogl_init(char *display, char **title, double *x, double *y, double *w, double *h, 
	 int n_windows);
void
ogl_close( void );
int
ogl_start_plot(int new_plot, int win_number, double scale);
int
ogl_end_plot(void);
void 
ogl_erase(int win_number);
void
ogl_standard_view( int window_number, bounding_box *bb );
int
ogl_reset_modelview(int win_number);
int
ogl_poll(void);
void 
ogl_drawstr(double x_pos[3], double n_vec[3], char *str);
void
ogl_marker(double x, double y, double z, double radius, int color);
void 
ogl_set_color( int color );
char *
ogl_color_name( int color );
void 
ogl_wait_for_key(char *prompt);
void
ogl_bind_rendering(int win_number);
void
ogl_draw_arrow(double length, double thickness, double start[3], double forward[3],
	       double normal[3], char *label);
void 
ogl_bb_cage( bounding_box *bb, int color );
double
ogl_length_scale( bounding_box *bb );
void
ogl_unnormalize( int window_number,  bounding_box *bb );
void
ogl_normalize( int window_number, bounding_box *bb );

void PL_set_same_scale( int window_number, int flag );
void PL_set_labels( int window_number, char *new_x_label, char *new_y_label );
void PL_scale(int, int);
void PL_title(int, char *);
void PL_postscript(int w, FILE *fp, int col_flag);

/* Defines for command interpreter */

#define RAWON     1
#define RAWOFF    0
#define RAWNDELON 2

/* Command buffer sizes */

#define MAX_LINE_SIZE 256
#define SCROLLBACK_SIZE 100

/* Function prototypes for command interpreter */

int get_command(input_output *io_ptr, char *prompt,
		char **command, int n_command,
		int command_level, int *save_on_copy, int *argument);
char *get_word(input_output *io_ptr, char *prompt, char *deflt,
	       int save_on_copy);
int get_int(input_output *io_ptr, char *prompt, int deflt, int level);
real get_real(input_output *io_ptr, char *prompt, real deflt, int level);
float get_float(input_output *io_ptr, char *prompt, float deflt, int level);
void getline_no_block(char *prompt, char *buf, char **cmds, int n);
void general_help( void );
void print_help( char *brief_help );
FILE *open_binary_file(input_output *io_ptr, char *prompt, char *deflt, 
		     char read_write, int file_type);
FILE *open_ascii_file( input_output *io_ptr, char *prompt, char *deflt, 
		      char **file_name, char read_write, int file_type, 
		      int save_on_copy);
FILE *open_this_ascii_file( char *file_name, char read_write, int file_type, 
			   int quiet);
FILE *open_this_binary_file( char *file_name, int file_type, int quiet);
int get_on_off( input_output *io_ptr, char *prompt );
int get_yes_no( input_output *io_ptr, char *prompt, int level, int save_answer );


#endif
