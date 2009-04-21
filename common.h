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
#define   PLOTDIR                   "plots"
#define   T0_default                4.2             /* Kelvin degrees */
#define   Gamma0_default            0.7
#define   Delta0_default            2.5
#define   k_B_default               0.086           /* Boltzmann */
#define   ExpDataFile_default       DATADIR"/sample.dat" 
#define   Vi_default                -HUGE_VAL
#define   Vf_default                HUGE_VAL

#define   MaxExpPoints              4096

#define   TheorExpPlotRatio         1.05
#define   SquaredResidualsPlotFile  PLOTDIR"/squared_residuals.plt"

/* for num. integration: */
#define   SUBINTERVALS              10000 
#define   ABSTOL                    0.0000        
#define   RELTOL                    0.001  

/* for num. differentiation */
#define   DSTEPSIZE                 2e-7  

/* to see (through a plot) how the sum of the squared residuals changes with 
 * parameters; i.e.: a rectangle is defined in the (Gamma, Delta) plane */
#define PLOTSTEPS 60
#define PlotSquareResiduals_default 0

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
extern double Gamma0, Delta0, T0, k_B, Vi, Vf;
extern char ExpDataFile[BUFSIZ];
extern BOOL PlotSquareResiduals;
#endif


/*Library*/
gsl_complex
cpow_2(const gsl_complex z);

/* Functions */

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
Gin_Gamma(const double Gamma, void * params);
double 
dGin_dGamma(const double V, const double Gamma, const double Delta, const double T);
double 
Gin_Delta(const double Delta, void * params);
double 
dGin_dDelta(const double V, const double Gamma, const double Delta, const double T);

int 
residuals_vector_f(const gsl_vector * params, void * data, gsl_vector * f);

double 
squared_residuals(const gsl_vector * params, void * data);
double 
squared_residuals_w_constraints(const gsl_vector * params, void * data);
int
simplex(const double Gamma_init, 
        const double Delta_init, 
        struct data * d, 
        double * Gamma_best, 
        double * Delta_best);

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
    const double Delta_init,    
    double *Gamma, 
    double *Delta, 
    gsl_matrix * cov);

void
print_fit_state (size_t iter, gsl_multifit_fdfsolver * s, size_t DoF);

int
init(void);

int
ui(void);

int 
plot(const double Gamma, const double Delta, struct data * d);

const char * 
BOOL2yn(BOOL b);

BOOL
yn2BOOL(char * str);




