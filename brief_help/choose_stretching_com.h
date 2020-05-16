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
enum {EXPONENTIAL, LAYER, TANH, HELP, CANCEL};
const int LAST_COM = CANCEL;
const int LEVEL=3;

char *COMMAND[5], *BRIEF_HELP[5];
int *ARGUMENT=NULL, *SAVE_ON_COPY=NULL;

COMMAND[EXPONENTIAL] = "exponential-stretching";
/*file[EXPONENTIAL] = "set_exp_com.h"; */
BRIEF_HELP[EXPONENTIAL] = "Introduce an exponential stretching. This stretching "
"is often used to concentrate grid points close to the starting or ending point "
"of a curve. Let the uniform parameter be u and the stretched be s. The stretching "
"function is\n"
"\n"
"s(u) = (exp(a*u) - 1)/(exp(a) - 1),  0 <= u <= 1,\n"
"\n"
"where `a' is a constant that determines the strength of the stretching.";

COMMAND[LAYER] = "layer-stretching"; 
/*file[LAYER] = "set_layer_com.h"; */
BRIEF_HELP[LAYER] = "Use an inverse tanh stretching to attract grid points to "
"the end points. Let the uniform parameter be u and the stretched parameter be "
"s. The inverse of the stretching function is\n"
"\n"
"u(s) = s + constant * sum( 0.5*Ai*tanh(Bi*(s-Ci)) - Di ), 0 <= s <= 1,\n"
"\n"
"for 1 <= i <= 2, where C1=0, and C2=1. Ai is the strength and Bi is the "
"inverse of the width. Di and the constant are adjusted to make u(0) = 0 "
"and u(1) = 1\n"
"Observe that the above function u(s) is inverted to give the stretching "
"function s(u).";

COMMAND[TANH] = "hyperbolic-tangent-stretching"; 
/*file[TANH] = "set_tanh_com.h"; */
BRIEF_HELP[TANH] = "Introduce a hyperbolic tangent stretching function. "
"This function can be used to attract grid points to both ends of the curve. "
"Let the uniform parameter be u and the stretched be s. The stretching "
"function is\n"
"\n"
"s(u) = w(u) / (a + (1-a)*w(u)),\n"
"w(u) = 0.5*( 1 + (tanh(d*(u-0.5))/tanh(d*0.5)) )\n"
"\n"
"The parameters `a' and `d' adjusted to give specified grid sizes at u=0 and u=1.";

COMMAND[HELP] = "help";
BRIEF_HELP[HELP] = NULL;

COMMAND[CANCEL] = "cancel";
BRIEF_HELP[CANCEL] = "Go directly back to the previous command level without "
"picking a stretching.";

