/* -*- linux-c -*- */
/* fewbody.h

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

#ifndef _FEWBODY_H
#define _FEWBODY_H 1

#include <stdio.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_rng.h>

/* version information */
#define FB_VERSION "0.26"
#define FB_NICK "Oblivion"
#define FB_DATE "Mon Sep  6 09:05:47 EDT 2010"

/* dimensionless constants */
#define FB_CONST_PI 3.141592653589793238462643

/* constants, in cgs units */
#define FB_CONST_MSUN 1.989e+33
#define FB_CONST_RSUN 6.9599e+10
#define FB_CONST_C 2.99792458e+10
#define FB_CONST_G 6.67259e-8
#define FB_CONST_AU 1.496e+13
#define FB_CONST_PARSEC 3.0857e+18
#define FB_CONST_YR 3.155693e+7

/* these usually shouldn't need to be changed */
#define FB_H 1.0e-2
#define FB_SSTOP GSL_POSINF
#define FB_AMIN GSL_POSINF
#define FB_RMIN GSL_POSINF
#define FB_ROOTSOLVER_MAX_ITER 100
#define FB_ROOTSOLVER_ABS_ACC 1.0e-11
#define FB_ROOTSOLVER_REL_ACC 1.0e-11
#define FB_MAX_STRING_LENGTH 2048
#define FB_MAX_LOGENTRY_LENGTH (32 * FB_MAX_STRING_LENGTH)


/* a struct containing the units used */
typedef struct{
	double v; /* velocity */
	double l; /* length */
	double t; /* time */
	double m; /* mass */
	double E; /* energy */
} fb_units_t;

/* the fundamental object */
typedef struct fb_obj{
	int ncoll; /* total number of stars collided together in this star */
	long *id; /* numeric id array */
	char idstring[FB_MAX_STRING_LENGTH]; /* string id */
	double m; /* mass */
	double R; /* radius */
	double Eint; /* internal energy (used to check energy conservation) */
	double Lint[3]; /* internal ang mom (used to check ang mom conservation) */
	double x[3]; /* position */
	double v[3]; /* velocity */
	int n; /* total number of stars in hierarchy */
	struct fb_obj *obj[2]; /* pointers to children */
	double a; /* semimajor axis */
	double e; /* eccentricity */
	double Lhat[3]; /* angular momentum vector */
	double Ahat[3]; /* Runge-Lenz vector */
	double t; /* time at which node was upsynced */
	double mean_anom; /* mean anomaly when node was upsynced */
} fb_obj_t;

/* parameters for the K-S integrator */
typedef struct{
	int nstar; /* number of actual stars */
	int kstar; /* nstar*(nstar-1)/2, number of separations */
	double *m; /* m[nstar] */
	double *M; /* M[kstar] */
	double **amat; /* amat[nstar][kstar] */
	double **Tmat; /* Tmat[kstar][kstar] */
	double Einit; /* initial energy used in integration scheme */
} fb_ks_params_t;


/* parameters for the non-regularized integrator */
typedef struct{
	int nstar; /* number of actual stars */
	double *m; /* m[nstar] */
} fb_nonks_params_t;

/* the hierarchy data structure */
typedef struct{
	int nstarinit; /* initial number of stars (may not equal nstar if there are collisions) */
	int nstar; /* number of stars */
	int nobj; /* number of binary trees */
	int *hi; /* hierarchical index array */
	int *narr; /* narr[i] = number of hierarchical objects with i elements */
	fb_obj_t *hier; /* memory location of hierarchy information */
	fb_obj_t **obj; /* array of pointers to top nodes of binary trees */
} fb_hier_t;

/* input parameters */
typedef struct{
	int ks; /* 0=no regularization, 1=K-S regularization */
	double tstop; /* stopping time, in units of t_dyn */
	int Dflag; /* 0=don't print to stdout, 1=print to stdout */
	double dt; /* time interval between printouts will always be greater than this value */
	double tcpustop; /* cpu stopping time, in units of seconds */
	double absacc; /* absolute accuracy of the integrator */
	double relacc; /* relative accuracy of the integrator */
	int ncount; /* number of integration steps between each call to fb_classify() */
	double tidaltol; /* tidal tolerance */
	char firstlogentry[FB_MAX_LOGENTRY_LENGTH]; /* first entry to put in printout log */
	double fexp; /* expansion factor for a merger product: R = f_exp (R_1+R_2) */
} fb_input_t;

/* return parameters */
typedef struct{
	long count; /* number of integration steps */
	int retval; /* return value; 1=success, 0=failure */
	long iclassify; /* number of times classify was called */
	double tcpu; /* cpu time taken */
	double DeltaE; /* change in energy */
	double DeltaEfrac; /* change in energy, as a fraction of initial energy */
	double DeltaL; /* change in ang. mom. */
	double DeltaLfrac; /* change in ang. mom., as a fraction of initial ang. mom. */
	double Rmin; /* minimum distance of close approach during interaction */
	int Rmin_i; /* index of star i participating in minimum close approach */
	int Rmin_j; /* index of star j participating in minimum close approach */
	int Nosc; /* number of oscillations of the quantity s^2 (McMillan & Hut 1996) (Nosc=Nmin-1, so resonance if Nosc>=1) */
	double x0[3];
	double x1[3];
	double x2[3];
	double x3[3];
	double v0[3];
	double v1[3];
	double v2[3];
	double v3[3];
} fb_ret_t;

/* fewbody.c */
fb_ret_t fewbody(fb_input_t input, fb_hier_t *hier, double *t);

/* myfewbody.c */
/*
int fb_gene_aeb(void);
*/

/* fewbody_classify.c */
int fb_classify(fb_hier_t *hier, double t, double tidaltol);
int fb_is_stable(fb_obj_t *obj);
int fb_is_stable_binary(fb_obj_t *obj);
int fb_is_stable_triple(fb_obj_t *obj);
int fb_is_stable_quad(fb_obj_t *obj);
int fb_mardling(fb_obj_t *obj, int ib, int is);

/* fewbody_coll.c */
int fb_is_collision(double r, double R1, double R2);
int fb_collide(fb_hier_t *hier, double f_exp);
void fb_merge(fb_obj_t *obj1, fb_obj_t *obj2, int nstarinit, double f_exp);

/* fewbody_hier.c */
void fb_malloc_hier(fb_hier_t *hier);
void fb_init_hier(fb_hier_t *hier);
void fb_free_hier(fb_hier_t hier);
void fb_trickle(fb_hier_t *hier, double t);
void fb_elkcirt(fb_hier_t *hier, double t);
int fb_create_indices(int *hi, int nstar);
int fb_n_hier(fb_obj_t *obj);
char *fb_sprint_hier(fb_hier_t hier, char string[FB_MAX_STRING_LENGTH]);
char *fb_sprint_hier_hr(fb_hier_t hier, char string[FB_MAX_STRING_LENGTH]);
void fb_upsync(fb_obj_t *obj, double t);
void fb_randorient(fb_obj_t *obj, gsl_rng *rng);
void fb_downsync(fb_obj_t *obj, double t);
void fb_objcpy(fb_obj_t *obj1, fb_obj_t *obj2);

/* fewbody_int.c */
void fb_malloc_ks_params(fb_ks_params_t *ks_params);
void fb_init_ks_params(fb_ks_params_t *ks_params, fb_hier_t hier);
void fb_free_ks_params(fb_ks_params_t ks_params);
void fb_malloc_nonks_params(fb_nonks_params_t *nonks_params);
void fb_init_nonks_params(fb_nonks_params_t *nonks_params, fb_hier_t hier);
void fb_free_nonks_params(fb_nonks_params_t nonks_params);

/* fewbody_io.c */
void fb_print_version(FILE *stream);
void fb_print_story(fb_obj_t *star, int nstar, double t, char *logentry);

/* fewbody_isolate.c */
int fb_collapse(fb_hier_t *hier, double t, double tidaltol);
int fb_expand(fb_hier_t *hier, double t, double tidaltol);

/* fewbody_ks.c */
double fb_ks_dot(double x[4], double y[4]);
double fb_ks_mod(double x[4]);
void fb_calc_Q(double q[4], double Q[4]);
void fb_calc_ksmat(double Q[4], double Qmat[4][4]);
void fb_calc_amat(double **a, int nstar, int kstar);
void fb_calc_Tmat(double **a, double *m, double **T, int nstar, int kstar);
int fb_ks_func(double s, const double *y, double *f, void *params);
double fb_ks_Einit(const double *y, fb_ks_params_t params);
void fb_euclidean_to_ks(fb_obj_t **star, double *y, int nstar, int kstar);
void fb_ks_to_euclidean(double *y, fb_obj_t **star, int nstar, int kstar);

/* fewbody_nonks.c */
int fb_nonks_func(double t, const double *y, double *f, void *params);
int fb_nonks_jac(double t, const double *y, double *dfdy, double *dfdt, void *params);
void fb_euclidean_to_nonks(fb_obj_t **star, double *y, int nstar);
void fb_nonks_to_euclidean(double *y, fb_obj_t **star, int nstar);

/* fewbody_scat.c */
void fb_init_scattering(fb_obj_t *obj0, fb_obj_t *obj1, double vinf, double b, double rtid);
void fb_normalize(fb_hier_t *hier, fb_units_t units);

/* fewbody_utils.c */
double *fb_malloc_vector(int n);
double **fb_malloc_matrix(int nr, int nc);
void fb_free_vector(double *v);
void fb_free_matrix(double **m);
double fb_sqr(double x);
double fb_cub(double x);
double fb_dot(double x[3], double y[3]);
double fb_mod(double x[3]);
int fb_cross(double x[3], double y[3], double z[3]);
int fb_rotat(double A[3], double B[3], int n, double theta);
int fb_angmom(fb_obj_t *star, int nstar, double L[3]);
void fb_angmomint(fb_obj_t *star, int nstar, double L[3]);
double fb_einttot(fb_obj_t *star, int nstar);
double fb_petot(fb_obj_t *star, int nstar);
double fb_ketot(fb_obj_t *star, int nstar);
double fb_outerpetot(fb_obj_t **obj, int nobj);
double fb_outerketot(fb_obj_t **obj, int nobj);
double fb_kepler(double e, double mean_anom);
double fb_keplerfunc(double mean_anom, void *params);
double fb_reltide(fb_obj_t *bin, fb_obj_t *single, double r);

/* fewbody_ui.c */
int fbui_new_hier(fb_hier_t *hier, int n);
int fbui_delete_hier(fb_hier_t *hier);
fb_obj_t *fbui_hierarchy_element(fb_hier_t *hier, int n, int m);
fb_obj_t *fbui_hierarchy_single(fb_hier_t *hier, int m);
fb_obj_t *fbui_hierarchy_binary(fb_hier_t *hier, int m);
fb_obj_t *fbui_hierarchy_triple(fb_hier_t *hier, int m);
fb_obj_t *fbui_hierarchy_quadruple(fb_hier_t *hier, int m);
fb_obj_t *fbui_hierarchy_quintuple(fb_hier_t *hier, int m);
fb_obj_t *fbui_hierarchy_sextuple(fb_hier_t *hier, int m);
fb_obj_t *fbui_hierarchy_septuple(fb_hier_t *hier, int m);
fb_obj_t *fbui_hierarchy_octuple(fb_hier_t *hier, int m);
fb_obj_t *fbui_hierarchy_nonuple(fb_hier_t *hier, int m);
fb_obj_t *fbui_hierarchy_decuple(fb_hier_t *hier, int m);
int fbui_make_pair(fb_obj_t *parentobj, fb_obj_t *obj1, fb_obj_t *obj2);
int fbui_initialize_single(fb_obj_t *single, long id, char idstring[FB_MAX_STRING_LENGTH]);
fb_obj_t *fbui_tree(fb_hier_t *hier, int n);
int fbui_obj_ncoll_set(fb_obj_t *obj, int ncoll);
int fbui_obj_ncoll_get(fb_obj_t *obj);
int fbui_obj_id_set(fb_obj_t *obj, int index, long id);
long fbui_obj_id_get(fb_obj_t *obj, int index);
int fbui_obj_idstring_set(fb_obj_t *obj, char idstring[FB_MAX_STRING_LENGTH]);
char *fbui_obj_idstring_get(fb_obj_t *obj);
int fbui_obj_mass_set(fb_obj_t *obj, double mass);
double fbui_obj_mass_get(fb_obj_t *obj);
int fbui_obj_radius_set(fb_obj_t *obj, double radius);
double fbui_obj_radius_get(fb_obj_t *obj);
int fbui_obj_Eint_set(fb_obj_t *obj, double Eint);
double fbui_obj_Eint_get(fb_obj_t *obj);
int fbui_obj_Lint_set(fb_obj_t *obj, double Lint[3]);
int fbui_obj_Linti_set(fb_obj_t *obj, double Lint0, double Lint1, double Lint2);
double *fbui_obj_Lint_get(fb_obj_t *obj);
int fbui_obj_x_set(fb_obj_t *obj, double x[3]);
int fbui_obj_xi_set(fb_obj_t *obj, double x0, double x1, double x2);
double *fbui_obj_x_get(fb_obj_t *obj);
int fbui_obj_v_set(fb_obj_t *obj, double v[3]);
int fbui_obj_vi_set(fb_obj_t *obj, double v0, double v1, double v2);
double *fbui_obj_v_get(fb_obj_t *obj);
int fbui_obj_n_set(fb_obj_t *obj, int n);
int fbui_obj_n_get(fb_obj_t *obj);
int fbui_obj_obj_set(fb_obj_t *obj, int i, fb_obj_t *objtopointto);
int fbui_obj_left_child_set(fb_obj_t *obj, fb_obj_t *objtopointto);
int fbui_obj_right_child_set(fb_obj_t *obj, fb_obj_t *objtopointto);
fb_obj_t *fbui_obj_obj_get(fb_obj_t *obj, int i);
fb_obj_t *fbui_obj_left_child_get(fb_obj_t *obj);
fb_obj_t *fbui_obj_right_child_get(fb_obj_t *obj);
int fbui_obj_a_set(fb_obj_t *obj, double a);
double fbui_obj_a_get(fb_obj_t *obj);
int fbui_obj_e_set(fb_obj_t *obj, double e);
double fbui_obj_e_get(fb_obj_t *obj);
int fbui_obj_Lhat_set(fb_obj_t *obj, double Lhat[3]);
int fbui_obj_Lhati_set(fb_obj_t *obj, double Lhat0, double Lhat1, double Lhat2);
double *fbui_obj_Lhat_get(fb_obj_t *obj);
int fbui_obj_Ahat_set(fb_obj_t *obj, double Ahat[3]);
int fbui_obj_Ahati_set(fb_obj_t *obj, double Ahat0, double Ahat1, double Ahat2);
double *fbui_obj_Ahat_get(fb_obj_t *obj);
int fbui_obj_t_set(fb_obj_t *obj, double t);
double fbui_obj_t_get(fb_obj_t *obj);
int fbui_obj_mean_anom_set(fb_obj_t *obj, double mean_anom);
double fbui_obj_mean_anom_get(fb_obj_t *obj);

/* macros */
/* The variadic macro syntax here conforms to the C99 standard, but for some
   reason won't compile on Mac OSX with gcc. */
/* #define fb_dprintf(...) if (fb_debug) fprintf(stderr, __VA_ARGS__) */
/* The variadic macro syntax here is the old gcc standard, and compiles on
   Mac OSX with gcc. */
#define fb_dprintf(args...) if (fb_debug) fprintf(stderr, args)
#define FB_MIN(a, b) ((a)<=(b)?(a):(b))
#define FB_MAX(a, b) ((a)>=(b)?(a):(b))
#define FB_DELTA(i, j) ((i)==(j)?1:0)
#define FB_KS_K(i, j, nstar) ((i)*(nstar)-((i)+1)*((i)+2)/2+(j))

/* there is just one global variable */
extern int fb_debug;

#endif /* fewbody.h */
