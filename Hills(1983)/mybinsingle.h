/* -*- linux-c -*- */
/* binsingle.h

   Copyright (C) 2002-2004 John M. Fregeau
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#define FB_TIDALTOL 1.0e-4

#define FB_M00 (1.0e2 * FB_CONST_MSUN)
#define FB_M01 (1.0e2 * FB_CONST_MSUN)
#define FB_M1 (1.0 * FB_CONST_MSUN)

#define FB_R00 (0.0 * FB_CONST_RSUN)
#define FB_R01 (0.0 * FB_CONST_RSUN)
#define FB_R1 (0.0 * FB_CONST_RSUN)

#define FB_E0 0.0

#define FB_DT 1.0 /* approximate output dt */
#define FB_TCPUSTOP 1200.0 /* in seconds */

#define FB_ABSACC 1.0e-8 /* absolute accuracy of integrator */
#define FB_RELACC 1.0e-8 /* relative accuracy of integrator */
#define FB_NCOUNT 1 /* number of timesteps between calls to classify() */

#define FB_KS 0
#define FB_FEXP 3.0 /* expansion factor of merger product */
#define FB_DEBUG 0

#define FB_SEED 1
#define FB_BMIN 0.7
#define FB_BMAX 0.7 /* impact parameter in unit of pc */
#define FB_PINDEX 10 /* total times we need */

 
void print_usage(FILE *stream);
int calc_units(fb_obj_t *obj[2], fb_units_t *units);
