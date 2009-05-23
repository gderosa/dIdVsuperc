/* common.h 
 * 
 * Copyright (C) 2008, 2009 Guido De Rosa
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h> 
#include <gsl/gsl_multifit_nlin.h>

/* Constants and defaults */
#define   DATADIR                   "data"
#define   T0_default                4.200           /* Kelvin degrees */
#define   Gamma1_0_default          0.200
#define   Gamma2_0_default          0.100
#define   Delta1_0_default          1.700
#define   Delta2_0_default          3.200
#define   alpha1_0_default          0.667           /* alpha1+alpha2=1 */
#define   k_B_default               0.086           /* Boltzmann */
#define   ExpDataFile_default       DATADIR"/MgB2_01.dat" 
#define   Vi_default                -HUGE_VAL
#define   Vf_default                HUGE_VAL

#define   MaxExpPoints              4096

/* Plotting */
#define   ExtraPlotRatio            0.05 /* plot fit function slightly larger */

/* for num. integration: */
#define   SUBINTERVALS              10000 
#define   ABSTOL                    0.0000        
#define   RELTOL                    0.001  

/* for num. differentiation */
#define   DSTEPSIZE                 2e-7  

/* (chi squared) minimization algorithms                                      */
#define   MAX_SIMPLEX_ITER          400   /* Nelder-Mead "Simplex"            */
#define   MAX_FIT_ITER              100   /* Levenberg-Marquardt "downhill"   */
#define   MULTIMIN_TEST_SIZE        2e-5  /* Nelder-Mead "Simplex"            */
#define   MULTIFIT_TEST_DELTA       1e-7  /* Levenberg-Marquardt "downhill"   */

/* Constraints on parameters */
#define   CONSTRAINT_GAMMA1_MIN     0.000
#define   CONSTRAINT_GAMMA1_MAX     HUGE_VAL

#define   CONSTRAINT_GAMMA2_MIN     0.000
#define   CONSTRAINT_GAMMA2_MAX     HUGE_VAL

#define   CONSTRAINT_DELTA1_MIN     0.000
#define   CONSTRAINT_DELTA1_MAX     HUGE_VAL

#define   CONSTRAINT_DELTA2_MIN     0.000
#define   CONSTRAINT_DELTA2_MAX     HUGE_VAL

#define   CONSTRAINT_ALPHA1_MIN     0.666
#define   CONSTRAINT_ALPHA1_MAX     1.000    


/* Types */

#define BOOL unsigned short int /* ANSI C89 doesn't provide bool */

struct Gin_integrand_params 
  {
    double V, Gamma, Delta, T;
  };

struct Gin_params 
  {
    double Gamma, Delta, T;
  };
  
struct Gin_doubleDeltaGamma_params
  {
    double Gamma1, Gamma2, Delta1, Delta2, alpha1, T;
  };

struct Gin_squared_residuals_params 
  {
    size_t n;
    double *X, *Y;
  };

struct data { /* experimental points */
  size_t n;
  double * X;
  double * Y;
  double * sigmaY;
};


/* Globals */
#ifndef GLOBALS
extern double 
  Gamma1_0, Gamma2_0, 
  Delta1_0, Delta2_0, 
  alpha1_0, 
  T0, 
  k_B, 
  Vi, Vf;
extern char ExpDataFile[BUFSIZ];

#endif

/* Functions */

/* square of a complex number */
gsl_complex
cpow_2(const gsl_complex z);

/* Minimum and maximum element of an array -- TODO: use gsl_vector instead */
double
array_min(double *ary, size_t n);
double
array_max(double *ary, size_t n);

/* density of states */
double 
ni(const double u, const double Gamma, const double Delta);

/* Fermi distrib. */
double
fermi(const double u);

/* Fermi distribution derivative * (-1) */   
double 
fprime(const double u);

/* Differential conductance: integrand */
double
Gin_integrand_base(const double u,  const double V, const double Gamma, const double Delta, const double T);

/* Differential conductance: integrand GSL-ized */
double 
Gin_integrand(const double u, void * params);

double
Gin(const double V, const double Gamma, const double Delta, const double T);
 
double
Gin_doubleDeltaGamma(
  const double V, 
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2, 
  const double alpha1, 
  const double T);

double 
Gin_doubleDeltaGamma_Gamma1(const double Gamma1, void * params);
double 
Gin_doubleDeltaGamma_Gamma2(const double Gamma2, void * params);

double 
dGin_doubleDeltaGamma_dGamma1(
  const double V, 
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2, 
  const double alpha1, 
  const double T);
double 
dGin_doubleDeltaGamma_dGamma2(
  const double V, 
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2, 
  const double alpha1, 
  const double T);

double 
Gin_doubleDeltaGamma_Delta1(const double Delta1, void * params);

double 
dGin_doubleDeltaGamma_dDelta1(
  const double V, 
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2, 
  const double alpha1,  
  const double T);

double
Gin_doubleDeltaGamma_Delta2(const double Delta1, void * params);

double 
dGin_doubleDeltaGamma_dDelta2(
  const double V, 
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2, 
  const double alpha1,  
  const double T);

double
Gin_doubleDeltaGamma_alpha1(const double alpha1, void * params);

double 
dGin_doubleDeltaGamma_dalpha1(
  const double V, 
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2, 
  const double alpha1,  
  const double T);


int 
residuals_vector_f(const gsl_vector * params, void * data, gsl_vector * f);

double 
squared_residuals(const gsl_vector * params, void * data);

double 
squared_residuals_w_constraints(const gsl_vector * params, void * data);

int
simplex(
  const double Gamma1_init,
  const double Gamma2_init,    
  const double Delta1_init,
  const double Delta2_init,  
  const double alpha1_inti,
  struct data * d, 
  double * Gamma1_best, 
  double * Gamma2_best, 
  double * Delta1_best,
  double * Delta2_best,
  double * alpha1_best);

int
residuals_jacobian_df(const gsl_vector * params, void * data, gsl_matrix * J);
int
residuals_fdf(const gsl_vector * params, 
              void * data,
              gsl_vector * f,
              gsl_matrix * J);

int 
fit(struct data * d,    
    const double Gamma1_init,
    const double Gamma2_init,    
    const double Delta1_init,    
    const double Delta2_init,
    const double alpha1_init,
    double *Gamma1,
    double *Gamma2,      
    double *Delta1,
    double *Delta2,
    double *alpha1, 
    gsl_matrix * cov,
    double * chisquared);

void
print_fit_state (size_t iter, gsl_multifit_fdfsolver * s, size_t DoF);

int
init(void);

int
ui(void);

int 
plot(
  const double Gamma1, 
  const double Gamma2, 
  const double Delta1, 
  const double Delta2, 
  const double alpha1, 
  struct data * d);


int 
constraints
(
  double Gamma1,
  double Gamma2,  
  double Delta1,
  double Delta2,
  double alpha1
);
