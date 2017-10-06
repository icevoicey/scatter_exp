/* -*- linux-c -*- */
/* binbin.c

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
#include "mybinbin.h"

/* calculate the units used */
int calc_units(fb_obj_t *obj[2], fb_units_t *units)
{
	units->v = sqrt(FB_CONST_G*(obj[0]->m + obj[1]->m)/(obj[0]->m * obj[1]->m) * \
			(obj[0]->obj[0]->m * obj[0]->obj[1]->m / obj[0]->a + \
			 obj[1]->obj[0]->m * obj[1]->obj[1]->m / obj[1]->a));
	units->l = obj[0]->a + obj[1]->a;
	units->t = units->l / units->v;
	units->m = units->l * fb_sqr(units->v) / FB_CONST_G;
	units->E = units->m * fb_sqr(units->v);
	
	return(0);
}


/* the main attraction */
int main(int argc, char *argv[])
{
	int pindex = FB_PINDEX, pcount;
	int i, j;
	int seed;
	double r00, r01, r10, r11, a0, a1, e0, e1;
	double rtid, vinf, b, m0, m1, M, mu, Ei, E, Lint[3], Li[3], l0[3], l1[3], L[3], r[3], t;
	double sigma, r_inf, bmax,fac;
	double x,y,z,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,vx0,vy0,vz0,vx1,vy1,vz1,vx2,vy2,vz2,vx3,vy3,vz3,vx,vy,vz,q,v0;
	double a_ini, e_ini,rperi,my_r;
	double my_i,my_psi,my_theta,my_alpha,my_phi,my_phi0,my_delta;
	//double myx0,myy0,myz0,myx1,myy1,myz1,myvx0,myvy0,myvz0,myvx1,myvy1,myvz1;
	fb_hier_t hier;
	fb_input_t input;
	fb_ret_t retval;
	fb_units_t units;
	char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH];
	gsl_rng *rng;
	const gsl_rng_type *rng_type=gsl_rng_mt19937;


	/* initialize GSL rng */
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(rng_type);
	seed = FB_SEED;
	gsl_rng_set(rng, seed);

	FILE *fhead;
	char name[10] = "result";
	char name2[2];
	snprintf (name2, sizeof(name2), "%d", seed);
	strcat(name, name2);
	fhead = fopen(name,"w");
	fclose(fhead);

	/* start the loop */
	for (pcount=0; pcount < pindex; pcount++){

		r00 = FB_R00;
		r01 = FB_R01;
		r10 = FB_R10;
		r11 = FB_R11;

		sigma = 2.0e7 * pow((FB_M00 / (pow(10, 8.13) * FB_CONST_MSUN)), 1.0/4.02);
		r_inf = FB_CONST_G * (FB_M00) / pow(sigma, 2);

		a0 = 2.0 * r_inf; 
		a1 = exp((log(FB_A34MAX) - log(FB_A34MIN)) * gsl_rng_uniform(rng) + log(FB_A34MIN));
		e0 = FB_E0;
		e1 = pow(((pow(FB_E34MAX, 2) - pow(FB_E34MIN, 2)) * gsl_rng_uniform(rng) + pow(FB_E34MIN, 2)), 0.5);	

		input.ks = FB_KS;
		input.tstop = FB_TSTOP;
		input.Dflag = 0;
		input.dt = FB_DT;
		input.tcpustop = FB_TCPUSTOP;
		input.absacc = FB_ABSACC;
		input.relacc = FB_RELACC;
		input.ncount = FB_NCOUNT;
		input.tidaltol = FB_TIDALTOL;
		input.fexp = FB_FEXP;
		fb_debug = FB_DEBUG;

		/* initialize a few things for integrator */
		t = 0.0;
		hier.nstarinit = 4;
		hier.nstar = 4;
		fb_malloc_hier(&hier);
		fb_init_hier(&hier);

		/* create binaries */
		hier.hier[hier.hi[2]+0].obj[0] = &(hier.hier[hier.hi[1]+0]);
		hier.hier[hier.hi[2]+0].obj[1] = &(hier.hier[hier.hi[1]+1]);
		hier.hier[hier.hi[2]+0].t = t;
		hier.hier[hier.hi[2]+1].obj[0] = &(hier.hier[hier.hi[1]+2]);
		hier.hier[hier.hi[2]+1].obj[1] = &(hier.hier[hier.hi[1]+3]);
		hier.hier[hier.hi[2]+1].t = t;

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

		hier.hier[hier.hi[1]+0].R = r00;
		hier.hier[hier.hi[1]+1].R = r01;
		hier.hier[hier.hi[1]+2].R = r10;
		hier.hier[hier.hi[1]+3].R = r11;

		hier.hier[hier.hi[1]+0].m = FB_M00;
		hier.hier[hier.hi[1]+1].m = FB_M01;
		hier.hier[hier.hi[1]+2].m = FB_M10;
		hier.hier[hier.hi[1]+3].m = FB_M11;

		hier.hier[hier.hi[2]+0].m = FB_M00 + FB_M01;
		hier.hier[hier.hi[2]+1].m = FB_M10 + FB_M11;

		hier.hier[hier.hi[2]+0].a = a0;
		hier.hier[hier.hi[2]+1].a = a1;

		hier.hier[hier.hi[2]+0].e = e0;
		hier.hier[hier.hi[2]+1].e = e1;

		hier.obj[0] = &(hier.hier[hier.hi[2]+0]);
		hier.obj[1] = &(hier.hier[hier.hi[2]+1]);
		hier.obj[2] = NULL;
		hier.obj[3] = NULL;

		/* get the units and normalize */
		calc_units(hier.obj, &units);
		fb_normalize(&hier, units);

		vinf = pow(FB_M00 / (pow(10, 8.13) * FB_CONST_MSUN), 1.0/4.02) * 2.0e7 / units.v;
		fac=1.0;		
		bmax = fac * FB_CONST_PARSEC / units.l;
		//b = pow((pow(bmax, 2) * gsl_rng_uniform(rng)), 0.5);
		b = bmax;
		/*	
		bmax = 1 * pow(2 * FB_CONST_G * (FB_M00 + FB_M01) * r_bt, 0.5) / (fabs(vinf * units.v) * units.l);
				
		*/

		/* move hierarchies analytically in from infinity along hyperbolic orbit */
		m0 = hier.obj[0]->m;
		m1 = hier.obj[1]->m;
		M = m0 + m1;
		mu = m0 * m1 / M;

		Ei = 0.5 * mu * fb_sqr(vinf);

		a0 = hier.obj[0]->a;
		a1 = hier.obj[1]->a;

		e0 = hier.obj[0]->e;
		e1 = hier.obj[1]->e;

		rtid = pow(2.0*(m0+m1)/input.tidaltol, 1.0/3.0) * FB_MAX(pow(m0, -1.0/3.0)*a0*(1.0+e0), pow(m1, -1.0/3.0) * a1 * (1.0+e1));

		rperi = (vinf==0.0?0.0:(M/fb_sqr(vinf)*(sqrt(1.0+fb_sqr(b*fb_sqr(vinf)/M))-1.0)));

		/* make sure r>=rperi, otherwise analytically moving the obj's below will give NANs */
		my_r = FB_MAX(rtid, rperi);

		a_ini = -mu / fb_sqr(vinf);
		e_ini = sqrt(1.0 + fb_sqr(b) * pow(vinf,4.0) / fb_sqr(mu));
		my_phi0 = acos(-1.0 / e_ini);

		my_i = acos(1.0 * gsl_rng_uniform(rng));
		my_phi = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng);
		my_psi = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng);
		my_delta = FB_CONST_PI - my_i - my_phi0;
		my_theta=acos(-1.0 / e_ini + a_ini * (1.0-fb_sqr(e_ini)) / (e_ini * my_r));

		x0 = 0.0;
		y0 = my_r * sin(my_theta);
		z0 = -my_r * cos(my_theta);

		x1 = x0;
		y1=cos(my_i+my_delta)*y0+sin(my_i+my_delta)*z0;
		z1=-sin(my_i+my_delta)*y0+cos(my_i+my_delta)*z0;

		x2=cos(my_psi)*x1+sin(my_psi)*y1;
		y2=-sin(my_psi)*x1+cos(my_psi)*y1;
		z2=z1;

		x3=x2;
		y3=cos(my_i)*y2-sin(my_i)*z2;
		z3=sin(my_i)*y2+cos(my_i)*z2;

		x=cos(my_phi)*x3+sin(my_phi)*y3;
		y=-sin(my_phi)*x3+cos(my_phi)*y3;
		z=z3;

		v0=sqrt(fb_sqr(vinf) + 2.0 * mu / my_r);
		my_alpha=atan(y0/(fb_sqr(e_ini)-1.0)*(z0-e_ini*a_ini));

		vx0=0.0;
		vy0=-v0*cos(my_alpha);
		vz0=-v0*sin(my_alpha);

		vx1 = vx0;
		vy1=cos(my_i+my_delta)*vy0+sin(my_i+my_delta)*vz0;
		vz1=-sin(my_i+my_delta)*vy0+cos(my_i+my_delta)*vz0;

		vx2=cos(my_psi)*vx1+sin(my_psi)*vy1;
		vy2=-sin(my_psi)*vx1+cos(my_psi)*vy1;
		vz2=vz1;

		vx3=vx2;
		vy3=cos(my_i)*vy2-sin(my_i)*vz2;
		vz3=sin(my_i)*vy2+cos(my_i)*vz2;

		vx=cos(my_phi)*vx3+sin(my_phi)*vy3;
		vy=-sin(my_phi)*vx3+cos(my_phi)*vy3;
		vz=vz3;


		hier.obj[0]->x[0] = 0.0;
		hier.obj[0]->x[1] = 0.0;
		hier.obj[0]->x[2] = 0.0;
	
		hier.obj[0]->v[0] = 0.0;
		hier.obj[0]->v[1] = 0.0;
		hier.obj[0]->v[2] = 0.0;

		hier.obj[1]->x[0] = x;
		hier.obj[1]->x[1] = y;
		hier.obj[1]->x[2] = z;
	
		hier.obj[1]->v[0] = vx;
		hier.obj[1]->v[1] = vy;
		hier.obj[1]->v[2] = vz;

	
		/* and check to see that we conserved energy and angular momentum */
		fb_cross(hier.obj[0]->x, hier.obj[0]->v, l0);
		fb_cross(hier.obj[1]->x, hier.obj[1]->v, l1);
	
		for (i=0; i<3; i++) {
			L[i] = (m0 * l0[i] + m1 * l1[i]);
			r[i] = hier.obj[1]->x[i] - hier.obj[0]->x[i];
		}

		E = - m0 * m1 / fb_mod(r) + 0.5 * (m0 * fb_dot(hier.obj[0]->v, hier.obj[0]->v) + \
						m1 * fb_dot(hier.obj[1]->v, hier.obj[1]->v));

		/* trickle down the binary properties, then back up */
		fb_randorient(&(hier.hier[hier.hi[2]+1]), rng);
		fb_downsync(&(hier.hier[hier.hi[2]+1]), t);
		fb_upsync(&(hier.hier[hier.hi[2]+1]), t);

		
		q=hier.hier[hier.hi[2]+0].obj[1]->m / hier.hier[hier.hi[2]+0].obj[0]->m;	
		hier.hier[hier.hi[2]+0].obj[0]->x[0] = 0;
		hier.hier[hier.hi[2]+0].obj[0]->x[1] = a0/(1.0+1.0/q);
		hier.hier[hier.hi[2]+0].obj[0]->x[2] = 0;

		hier.hier[hier.hi[2]+0].obj[0]->v[0] = sqrt(hier.hier[hier.hi[2]+0].obj[1]->m/(1.0+1.0/q));
		hier.hier[hier.hi[2]+0].obj[0]->v[1] = 0;
		hier.hier[hier.hi[2]+0].obj[0]->v[2] = 0;

		hier.hier[hier.hi[2]+0].obj[1]->x[0] = 0;
		hier.hier[hier.hi[2]+0].obj[1]->x[1] = -a0/(1.0+q);
		hier.hier[hier.hi[2]+0].obj[1]->x[2] = 0;

		hier.hier[hier.hi[2]+0].obj[1]->v[0] = -sqrt(hier.hier[hier.hi[2]+0].obj[0]->m/(1.0+q));
		hier.hier[hier.hi[2]+0].obj[1]->v[1] = 0;
		hier.hier[hier.hi[2]+0].obj[1]->v[2] = 0;


		

	
		/* store the initial energy and angular momentum*/
		fb_angmom(&(hier.hier[hier.hi[1]]), hier.nstar, Li);
		fb_angmomint(&(hier.hier[hier.hi[1]]), hier.nstar, Lint);
		/*
		myx0=hier.hier[hier.hi[2]+0].obj[0]->x[0]*units.l/FB_CONST_PARSEC;
		myy0=hier.hier[hier.hi[2]+0].obj[0]->x[1]*units.l/FB_CONST_PARSEC;
		myz0=hier.hier[hier.hi[2]+0].obj[0]->x[2]*units.l/FB_CONST_PARSEC;
		myx1=hier.hier[hier.hi[2]+0].obj[1]->x[0]*units.l/FB_CONST_PARSEC;
		myy1=hier.hier[hier.hi[2]+0].obj[1]->x[1]*units.l/FB_CONST_PARSEC;
		myz1=hier.hier[hier.hi[2]+0].obj[1]->x[2]*units.l/FB_CONST_PARSEC;

		myvx0=hier.hier[hier.hi[2]+0].obj[0]->v[0]*units.v;
		myvy0=hier.hier[hier.hi[2]+0].obj[0]->v[1]*units.v;
		myvz0=hier.hier[hier.hi[2]+0].obj[0]->v[2]*units.v;
		myvx1=hier.hier[hier.hi[2]+0].obj[1]->v[0]*units.v;
		myvy1=hier.hier[hier.hi[2]+0].obj[1]->v[1]*units.v;
		myvz1=hier.hier[hier.hi[2]+0].obj[1]->v[2]*units.v;


		myx0=hier.hier[hier.hi[1]+2].x[0]*units.l/FB_CONST_PARSEC;
		myy0=hier.hier[hier.hi[1]+2].x[1]*units.l/FB_CONST_PARSEC;
		myz0=hier.hier[hier.hi[1]+2].x[2]*units.l/FB_CONST_PARSEC;
		myx1=hier.hier[hier.hi[1]+3].x[0]*units.l/FB_CONST_PARSEC;
		myy1=hier.hier[hier.hi[1]+3].x[1]*units.l/FB_CONST_PARSEC;
		myz1=hier.hier[hier.hi[1]+3].x[2]*units.l/FB_CONST_PARSEC;

		myvx0=hier.hier[hier.hi[1]+2].v[0]*units.v;
		myvy0=hier.hier[hier.hi[1]+2].v[1]*units.v;
		myvz0=hier.hier[hier.hi[1]+2].v[2]*units.v;
		myvx1=hier.hier[hier.hi[1]+3].v[0]*units.v;
		myvy1=hier.hier[hier.hi[1]+3].v[1]*units.v;
		myvz1=hier.hier[hier.hi[1]+3].v[2]*units.v;
		*/

		/* call fewbody! */
		retval = fewbody(input, &hier, &t);

		/* save the data into file */
		///*
		FILE *fbody;
		fbody=fopen(name, "a");
		if (retval.retval == 1) {
			fprintf(fbody, "1,%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n",
				fb_sprint_hier(hier, string1), fb_sprint_hier_hr(hier, string2),\
retval.x0[0]*units.l,retval.x0[1]*units.l,retval.x0[2]*units.l,retval.x1[0]*units.l,retval.x1[1]*units.l,retval.x1[2]*units.l,retval.x2[0]*units.l,retval.x2[1]*units.l,retval.x2[2]*units.l,retval.x3[0]*units.l,retval.x3[1]*units.l,retval.x3[2]*units.l,\
retval.v0[0]*units.v,retval.v0[1]*units.v,retval.v0[2]*units.v,retval.v1[0]*units.v,retval.v1[1]*units.v,retval.v1[2]*units.v,retval.v2[0]*units.v,retval.v2[1]*units.v,retval.v2[2]*units.v,retval.v3[0]*units.v,retval.v3[1]*units.v,retval.v3[2]*units.v,\
a1*units.l/FB_CONST_AU, e1, b* units.l/FB_CONST_PARSEC,retval.DeltaLfrac,retval.DeltaEfrac,(retval.Nosc>=1?"resonance":"non-resonance"));
		} 
		if (retval.retval == 0 ) {
			fprintf(fbody, "0,%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n",
				fb_sprint_hier(hier, string1), fb_sprint_hier_hr(hier, string2),\
retval.x0[0]*units.l,retval.x0[1]*units.l,retval.x0[2]*units.l,retval.x1[0]*units.l,retval.x1[1]*units.l,retval.x1[2]*units.l,retval.x2[0]*units.l,retval.x2[1]*units.l,retval.x2[2]*units.l,retval.x3[0]*units.l,retval.x3[1]*units.l,retval.x3[2]*units.l,\
retval.v0[0]*units.v,retval.v0[1]*units.v,retval.v0[2]*units.v,retval.v1[0]*units.v,retval.v1[1]*units.v,retval.v1[2]*units.v,retval.v2[0]*units.v,retval.v2[1]*units.v,retval.v2[2]*units.v,retval.v3[0]*units.v,retval.v3[1]*units.v,retval.v3[2]*units.v,\
a1*units.l/FB_CONST_AU, e1, b* units.l/FB_CONST_PARSEC,retval.DeltaLfrac,retval.DeltaEfrac,(retval.Nosc>=1?"resonance":"non-resonance"));
		} 
		
		fclose(fbody);
		//*/

		/* free our own stuff */
		fb_free_hier(hier);

		printf("%d\n ", pcount);
	}

	/* print out values of paramaters */
	///*
	FILE *fend;
	fend = fopen(name, "a");
	fprintf(fend, "PARAMETERS:\n");
	fprintf(fend, "a0=%.6g AU  e0=%.3g\n", a0*units.l/FB_CONST_AU, e0);
	fprintf(fend,  "m00=%.6g MSUN  m01=%.6g MSUN  m10=%.6g MSUN  m11=%.6g MSUN \n", FB_M00/FB_CONST_MSUN, FB_M01/FB_CONST_MSUN, FB_M10/FB_CONST_MSUN, FB_M11/FB_CONST_MSUN);
	fprintf(fend, "a1min=%.6g AU  a1max=%.6g AU   e1min=%.6g  e1max=%.6g\n", FB_A34MIN/FB_CONST_AU, FB_A34MAX/FB_CONST_AU, FB_E34MIN, FB_E34MAX);
	fprintf(fend, "vinf=%.6g cm/s  bmax=%.6g pc  tstop=%.6g  tcpustop=%.6g\n", vinf*units.v, bmax*units.l/FB_CONST_PARSEC, FB_TSTOP, FB_TCPUSTOP);
	fprintf(fend, "tidaltol=%.6g  abs_acc=%.6g  rel_acc=%.6g  ncount=%d  fexp=%.6g  seed=%d  pindex=%d\n", FB_TIDALTOL, FB_ABSACC, FB_RELACC, FB_NCOUNT, FB_FEXP, FB_SEED, FB_PINDEX);
	fclose(fend);
	//*/

	/* free GSL stuff */
	gsl_rng_free(rng);

	/* done! */
	return(0);
}

