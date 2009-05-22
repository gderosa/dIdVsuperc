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
#define   T0_default                4.2             /* Kelvin degrees */
#define   Gamma0_default            0.7
#define   Delta1_0_default          2.5
#define   Delta2_0_default          7.0
#define   alpha1_0_default          0.7             /* alpha1+alpha2=1 */
#define   k_B_default               0.086           /* Boltzmann */
#define   ExpDataFile_default       DATADIR"/sample.dat" 
#define   Vi_default                -HUGE_VAL
#define   Vf_default                HUGE_VAL

#define   MaxExpPoints              4096

/* Plotting */
#define   ExtraPlotRatio            0.05 /* plot fit function slightly larger */

/* 3D plotting of chi^2 as a function of Gamma and Delta (DISABLED)
 * #define   SPLOTSTEPS                72 
 * #define   SPLOTSIGMAS               10
 */

/* for num. integration: */
#define   SUBINTERVALS              10000 
#define   ABSTOL                    0.0000        
#define   RELTOL                    0.001  

/* for num. differentiation */
#define   DSTEPSIZE                 2e-7  

/* Types */

#define BOOL unsigned short int /* ANSI C89 / ISO C90 */

struct Gin_integrand_params 
  {
    double V, Gamma, Delta, T;
  };

struct Gin_params 
  {
    double Gamma, Delta, T;
  };
  
struct Gin_params_doubleDelta
  {
    double Gamma, Delta1, Delta2, alpha1, T;
  };

struct Gin_squared_residuals_params 
  {
    size_t n;
    double *X, *Y;
  };
#define Gin_squared_residuals_params_doubleDelta Gin_squared_residuals_params 

struct data { /* experimental points */
  size_t n;
  double * X;
  double * Y;
  double * sigmaY;
};


/* Globals */
#ifndef GLOBALS
extern double Gamma0, Delta1_0, Delta2_0, alpha1_0, T0, k_B, Vi, Vf;
extern char ExpDataFile[BUFSIZ];
extern BOOL PlotSquareResiduals;
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
Gin_doubleDelta(
  const double V, 
  const double Gamma, 
  const double Delta1, 
  const double Delta2, 
  const double alpha1, 
  const double T);


double 
Gin_doubleDelta_Gamma(const double Gamma, void * params);

double 
dGin_doubleDelta_dGamma(
  const double V, 
  const double Gamma, 
  const double Delta1, 
  const double Delta2, 
  const double alpha1, 
  const double T);

double 
Gin_doubleDelta_Delta1(const double Delta1, void * params);

double 
dGin_doubleDelta_dDelta1(
  const double V, 
  const double Gamma, 
  const double Delta1, 
  const double Delta2, 
  const double alpha1,  
  const double T);

double
Gin_doubleDelta_Delta2(const double Delta1, void * params);

double 
dGin_doubleDelta_dDelta2(
  const double V, 
  const double Gamma, 
  const double Delta1, 
  const double Delta2, 
  const double alpha1,  
  const double T);

double
Gin_doubleDelta_alpha1(const double alpha1, void * params);

double 
dGin_doubleDelta_dalpha1(
  const double V, 
  const double Gamma, 
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
  const double Gamma_init, 
  const double Delta1_init,
  const double Delta2_init,  
  const double alpha1_inti,
  struct data * d, 
  double * Gamma_best, 
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
    const double Gamma_init,
    const double Delta1_init,    
    const double Delta2_init,
    const double alpha1_init,
    double *Gamma, 
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
plot(const double Gamma, const double Delta1, const double Delta2, const double alpha1, struct data * d);

/* DISABLED 3D plotting of.... argh, actually it would require 5 dimensions!!
 *int
 * splot(const double Gamma, const double Delta, const double DGamma, const double DDelta, struct data * d); 
 */

/* the following were required only for chi-squared plotting, so....
 * const char * 
 * BOOL2yn(BOOL b);
 *
 * BOOL
 * yn2BOOL(char * str);
 */



