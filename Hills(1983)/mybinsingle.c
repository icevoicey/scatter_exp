/* -*- linux-c -*- */
/* binsingle.c

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

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include "fewbody.h"
#include "mybinsingle.h"


/* calculate the units used */
int calc_units(fb_obj_t *obj[2], fb_units_t *units)
{
	units->v = sqrt(FB_CONST_G*(obj[0]->m + obj[1]->m)/(obj[0]->m * obj[1]->m) * \
			(obj[0]->obj[0]->m * obj[0]->obj[1]->m / obj[0]->a));
	units->l = obj[0]->a;
	units->t = units->l / units->v;
	units->m = units->l * fb_sqr(units->v) / FB_CONST_G;
	units->E = units->m * fb_sqr(units->v);
	
	return(0);
}

/* the main attraction */
int main(int argc, char *argv[])
{
	int pindex = FB_PINDEX, pcount, j, seed;
	double a0, e0, test;
	double rtid, vtid, vinf, b, m0, m1, M, mu, t;
	double sigma, r_inf, q, v_omega,v_c;
	double x[3], x0[3], x1[3], x2[3], x3[3], v[3], v0[3], v1[3], v2[3], v3[3];
	double a_ini, e_ini, rperi, my_r;
	double my_i, my_psi, my_theta, my_alpha, my_phi, my_phi0, my_delta;
	fb_hier_t hier;
	fb_input_t input;
	fb_ret_t retval;
	fb_units_t units;
	char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH];
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;
	double t_hubble = 1.3e10*365*24*60*60; // hubble time of the universe set as 13 billion years

	/* initialize GSL rng */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);
	seed = FB_SEED;
	gsl_rng_set(rng, seed);

	/* make the result file with the seed used as ending number */
	FILE *fhead;
	char name[10] = "result";
	char name2[4];
	snprintf (name2, sizeof(name2), "%d", seed);
	strcat(name, name2);
	fhead = fopen(name,"w");
	fclose(fhead);

	/* start the loop */
	for (pcount=0; pcount < pindex; pcount++){

		sigma = 2.0e7 * pow((FB_M00 / (pow(10.0, 8.13) * FB_CONST_MSUN)), 1.0/4.02); // velocity dispersion of star 
		r_inf = FB_CONST_G * (FB_M00) / pow(sigma, 2.0); // influence radius of black hole

		v_c=sqrt(FB_CONST_G*(FB_M00+FB_M01)/(2.0 * r_inf));

		/* initialize a few things for integrator */
		t = 0.0;
		hier.nstarinit = 3;
		hier.nstar = 3;
		fb_malloc_hier(&hier);
		fb_init_hier(&hier);

		/* create binary */
		hier.hier[hier.hi[2]+0].obj[0] = &(hier.hier[hier.hi[1]+0]);
		hier.hier[hier.hi[2]+0].obj[1] = &(hier.hier[hier.hi[1]+1]);
		hier.hier[hier.hi[2]+0].t = t;

		/* give the objects some properties */
		for (j=0; j<hier.nstar; j++) {
			hier.hier[hier.hi[1]+j].ncoll = 1;
			hier.hier[hier.hi[1]+j].id[0] = j;
			snprintf(hier.hier[hier.hi[1]+j].idstring, FB_MAX_STRING_LENGTH, "%d", j);
			hier.hier[hier.hi[1]+j].n = 1;
			hier.hier[hier.hi[1]+j].obj[0] = NULL;
			hier.hier[hier.hi[1]+j].obj[1] = NULL;
			hier.hier[hier.hi[1]+j].Eint = 0.0;
			hier.hier[hier.hi[1]+j].Lint[0] = 0.0;
			hier.hier[hier.hi[1]+j].Lint[1] = 0.0;
			hier.hier[hier.hi[1]+j].Lint[2] = 0.0;
		}

		hier.hier[hier.hi[1]+0].R = FB_R00;
		hier.hier[hier.hi[1]+1].R = FB_R01;
		hier.hier[hier.hi[1]+2].R = FB_R1;

		hier.hier[hier.hi[1]+0].m = FB_M00;
		hier.hier[hier.hi[1]+1].m = FB_M01;
		hier.hier[hier.hi[1]+2].m = FB_M1;

		hier.hier[hier.hi[2]+0].m = FB_M00 + FB_M01;

		hier.hier[hier.hi[2]+0].a = 2.0 * r_inf;
		hier.hier[hier.hi[2]+0].e = FB_E0;

		hier.obj[0] = &(hier.hier[hier.hi[2]+0]);
		hier.obj[1] = &(hier.hier[hier.hi[1]+2]);
		hier.obj[2] = NULL;

		/* get the units and normalize */
		calc_units(hier.obj, &units);
		fb_normalize(&hier, units);

		/* set parameters to default values */
		input.ks = FB_KS;
		//input.tstop = 0.0001;
		input.tstop = t_hubble / units.t;
		input.Dflag = 0;
		input.dt = FB_DT;
		input.tcpustop = FB_TCPUSTOP;
		input.absacc = FB_ABSACC;
		input.relacc = FB_RELACC;
		input.ncount = FB_NCOUNT;
		input.tidaltol = FB_TIDALTOL;
		input.fexp = FB_FEXP;
		fb_debug = FB_DEBUG;

		b = pow((((pow(FB_BMAX, 2) - pow(FB_BMIN, 2)) * gsl_rng_uniform(rng)) + pow(FB_BMIN, 2.0)), 0.5) * FB_CONST_PARSEC / units.l;

		/* move hierarchies analytically in from infinity along hyperbolic orbit */
		m0 = hier.obj[0]->m;
		m1 = hier.obj[1]->m;
		M = m0 + m1;
		mu = m0 * m1 / M;

		a0 = hier.obj[0]->a;
		e0 = hier.obj[0]->e;

		vinf = 0.1 * v_c / units.v;

		rtid = pow(2.0 * (m0 + m1) / (m0 * input.tidaltol), 1.0/3.0) * a0 * (1.0+e0);

		rperi = (vinf==0.0?0.0:(M/fb_sqr(vinf)*(sqrt(1.0+fb_sqr(b*fb_sqr(vinf)/M))-1.0)));

		/* make sure r>=rperi, otherwise analytically moving the obj's below will give NANs */
		my_r = FB_MAX(rtid, rperi);

		a_ini = -mu / fb_sqr(vinf);
		e_ini = sqrt(1.0 + fb_sqr(b) * pow(vinf,4.0) / fb_sqr(mu)); // major semi-axis amd eccemtricity of the hyperbolic orbit moving from infinite to the start point of intergration (rtid)
		my_phi0 = acos(-1.0 / e_ini);

		my_i = acos(1.0 * gsl_rng_uniform(rng));
		my_phi = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng); // two angle position in the sphereical coordinate system setting the initial point of stellar binary on a sphereical surface at infinite 
		my_psi = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng);
		my_delta = FB_CONST_PI - my_i - my_phi0;
		my_theta = acos(-1.0 / e_ini + a_ini * (1.0-fb_sqr(e_ini)) / (e_ini * my_r));

		/* position coordinate at my_r(rtid) in the orbit framework */
		x0[0] = 0.0;
		x0[1] = my_r * sin(my_theta);
		x0[2] = -my_r * cos(my_theta);

		/* coordinate transform from the orbit framework to the black hole binary framework */
		fb_rotat(x0,x1,0,-my_delta-my_i);
		fb_rotat(x1,x2,2,-my_psi);
		fb_rotat(x2,x3,0,my_i);
		fb_rotat(x3,x,2,-my_phi);

		/* velocity coordinate at my_r(rtid) in the orbit framework */
		vtid = sqrt(fb_sqr(vinf) + 2.0 * mu / my_r);
		my_alpha = atan(x0[1] / (fb_sqr(e_ini) - 1.0) * (x0[2] - e_ini * a_ini)); 
		v0[0] = 0.0;
		v0[1] = -vtid * cos(my_alpha);
		v0[2] = -vtid * sin(my_alpha);

		/* coordinate transform from the orbit framework to the black hole binary framework */
		fb_rotat(v0,v1,0,-my_delta-my_i);
		fb_rotat(v1,v2,2,-my_psi);
		fb_rotat(v2,v3,0,my_i);
		fb_rotat(v3,v,2,-my_phi);

		/* coordinate of center of mass of black hole binary */
		hier.obj[0]->x[0] = 0.0;
		hier.obj[0]->x[1] = 0.0;
		hier.obj[0]->x[2] = 0.0;

		hier.obj[0]->v[0] = 0.0;
		hier.obj[0]->v[1] = 0.0;
		hier.obj[0]->v[2] = 0.0;

		/* coordinate of center of mass of star */
		hier.obj[1]->x[0] = x[0];
		hier.obj[1]->x[1] = x[1];
		hier.obj[1]->x[2] = x[2];

		hier.obj[1]->v[0] = v[0];
		hier.obj[1]->v[1] = v[1];
		hier.obj[1]->v[2] = v[2];

		/* set the initial condition of black hole binary */		
		q = hier.hier[hier.hi[2]+0].obj[1]->m / hier.hier[hier.hi[2]+0].obj[0]->m; // 0<q<1, m1<m0
		v_omega = sqrt((hier.hier[hier.hi[2]+0].obj[0]->m + hier.hier[hier.hi[2]+0].obj[1]->m) / pow(a0, 3.0)); //angle velocity of black hole binary
		
		hier.hier[hier.hi[2]+0].obj[0]->x[0] = 0.0;
		hier.hier[hier.hi[2]+0].obj[0]->x[1] = a0 / (1.0 + 1.0 / q);
		hier.hier[hier.hi[2]+0].obj[0]->x[2] = 0.0;

		hier.hier[hier.hi[2]+0].obj[0]->v[0] = v_omega * a0 / (1.0 + 1.0 / q);
		hier.hier[hier.hi[2]+0].obj[0]->v[1] = 0.0;
		hier.hier[hier.hi[2]+0].obj[0]->v[2] = 0.0;

		hier.hier[hier.hi[2]+0].obj[1]->x[0] = 0.0;
		hier.hier[hier.hi[2]+0].obj[1]->x[1] = -a0 / (1.0 + q);
		hier.hier[hier.hi[2]+0].obj[1]->x[2] = 0.0;

		hier.hier[hier.hi[2]+0].obj[1]->v[0] = -v_omega * a0 / (1.0 + q);
		hier.hier[hier.hi[2]+0].obj[1]->v[1] = 0.0;
		hier.hier[hier.hi[2]+0].obj[1]->v[2] = 0.0;

		
		/* trickle down the binary properties, then back up */
		//fb_downsync(&(hier.hier[hier.hi[2]+0]), t);
		//fb_upsync(&(hier.hier[hier.hi[2]+0]), t);

		/* call fewbody! */
		retval = fewbody(input, &hier, &t);
		test =my_r*units.l/FB_CONST_PARSEC;
		/* save the data into result file */
		
		FILE *fbody;
		fbody=fopen(name, "a");
		if (retval.retval == 1) {
			fprintf(fbody, "1,%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n",
				fb_sprint_hier(hier, string1), fb_sprint_hier_hr(hier, string2),\
				retval.x0[0]*units.l/FB_CONST_PARSEC,retval.x0[1]*units.l/FB_CONST_PARSEC,retval.x0[2]*units.l/FB_CONST_PARSEC,retval.x1[0]*units.l/FB_CONST_PARSEC,retval.x1[1]*units.l/FB_CONST_PARSEC,retval.x1[2]*units.l/FB_CONST_PARSEC,\
				retval.x2[0]*units.l/FB_CONST_PARSEC,retval.x2[1]*units.l/FB_CONST_PARSEC,retval.x2[2]*units.l/FB_CONST_PARSEC,\
				retval.v0[0]*units.v,retval.v0[1]*units.v,retval.v0[2]*units.v,retval.v1[0]*units.v,retval.v1[1]*units.v,retval.v1[2]*units.v,\
				retval.v2[0]*units.v,retval.v2[1]*units.v,retval.v2[2]*units.v,\
				b*units.l/FB_CONST_PARSEC,retval.DeltaLfrac,retval.DeltaEfrac,(retval.Nosc>=1?"resonance":"non-resonance"));
		} 

		if (retval.retval == 0) {
			fprintf(fbody, "0,%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n",
				fb_sprint_hier(hier, string1), fb_sprint_hier_hr(hier, string2),\
				retval.x0[0]*units.l,retval.x0[1]*units.l,retval.x0[2]*units.l,retval.x1[0]*units.l,retval.x1[1]*units.l,retval.x1[2]*units.l,\
				retval.x2[0]*units.l,retval.x2[1]*units.l,retval.x2[2]*units.l,\
				retval.v0[0]*units.v,retval.v0[1]*units.v,retval.v0[2]*units.v,retval.v1[0]*units.v,retval.v1[1]*units.v,retval.v1[2]*units.v,\
				retval.v2[0]*units.v,retval.v2[1]*units.v,retval.v2[2]*units.v,\
				b*units.l/FB_CONST_PARSEC,retval.DeltaLfrac,retval.DeltaEfrac,(retval.Nosc>=1?"resonance":"non-resonance"));
		} 	
		fclose(fbody);
		
		/* free our own stuff */
		fb_free_hier(hier);

		printf("%d\n", pcount);
		printf("%.6g\n",test);
	}

	/* save important values of paramaters at the end of result file */

	FILE *fend;
	fend = fopen(name, "a");
	fprintf(fend, "PARAMETERS:\n");
	fprintf(fend,  "m00=%.6g MSUN  m01=%.6g MSUN  m1=%.6g MSUN \n", FB_M00/FB_CONST_MSUN, FB_M01/FB_CONST_MSUN, FB_M1/FB_CONST_MSUN);
	fprintf(fend, "a0=%.6g AU  e0=%.3g\n", a0*units.l/FB_CONST_AU, e0);
	fprintf(fend, "vinf=%.6g m/s  bmin=%.6g pc  bmax=%.6g pc\n", vinf*units.v/100.0, FB_BMIN, FB_BMAX);
	fprintf(fend, "tidaltol=%.6g  abs_acc=%.6g  rel_acc=%.6g  ncount=%d  fexp=%.6g  seed=%d  num=%d\n", FB_TIDALTOL, FB_ABSACC, FB_RELACC, FB_NCOUNT, FB_FEXP, FB_SEED, FB_PINDEX);
	fclose(fend);

	/* free GSL stuff */
	gsl_rng_free(rng);

	/* done! */
	return(0);
}

	

