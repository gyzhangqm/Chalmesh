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
#include "math.h"
#include "real.h"
#include "log_cosh.h"

real log_cosh(real slope, real r0, real r1, real sharpness, real r){
  real v, lc_0, lc_1;
  const real big_c = 20.0;

/* this function needs to be protected against overflow */
  if (fabs(sharpness*(r-r0)) < big_c)
    lc_0 = log(cosh(sharpness*(r-r0)));
  else
    lc_0 = fabs(sharpness*(r-r0)) - log(2.0);

  if (fabs(sharpness*(r-r1)) < big_c)
    lc_1 = log(cosh(sharpness*(r-r1)));
  else
    lc_1 = fabs(sharpness*(r-r1)) - log(2.0);

  v = 0.5*slope/sharpness * (lc_0 - lc_1);
  return v;
}
