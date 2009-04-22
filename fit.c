/* fit.c - finding the best fit curve parameters
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


/* Squared residuals ("chi squared") are first minimized through a Simplex 
 * algorithm, before being minimized by an actual nonlinear fit (which, 
 * most importantly, estimates statistical errors on best-fit parameters)
 */

#include "common.h"


int 
residuals_vector_f(const gsl_vector * params, void * data, gsl_vector * f)
/* strongly based on http://www.gnu.org/software/gsl/manual/html_node/
 * Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html        
 * (expb_f function)
 */
{
  size_t n        = ((const struct data *)data)->n;
  double *X       = ((const struct data *)data)->X;
  double *Y       = ((const struct data *)data)->Y;
  double *sigmaY  = ((const struct data *)data)->sigmaY;
  
  double Gamma = gsl_vector_get (params, 0);
  double Delta = gsl_vector_get (params, 1);
  
  double theory, experiment;
  
  size_t i;
  
  gsl_vector_set_all(f, 0.0);
  
  for (i = 0; i < n; i++)
    {
      /* fitting an interval only */ 
        /* might be Vi=-HUGE_VAL and/or Vf=HUGE_VAL for "unlimited" */
      if (X[i]<Vi || X[i]>Vf) continue; 
      
      theory = Gin(X[i], Gamma, Delta, T0);
      experiment = Y[i];
      gsl_vector_set (f, i, (theory - experiment)/sigmaY[i]);
    }  
  return GSL_SUCCESS;
}

/* squared residuals ("chi squared") are first minimized through a Simplex algorithm,
 * before being minimized by an actual nonlinear fit (which, most importantly,
 * estimates statistical errors on best-fit parameters)
 */
double 
squared_residuals(const gsl_vector * params, void * data)
{
  struct data * d = (struct data *) data;
  size_t n = d->n;
  gsl_vector * f = gsl_vector_alloc(n);
  /* int status = residuals_vector_f(params, data, f); */
  double result;
  
  residuals_vector_f(params, data, f);
  
  result = gsl_pow_2(gsl_blas_dnrm2(f));
  
  gsl_vector_free(f);
  
  return result; 
}

/* A naive implementation of constraints: just return HUGE_VAL when parameters
 * are outside the desired domain
 */
double 
squared_residuals_w_constraints(const gsl_vector * params, void * data)
{
  double Gamma = gsl_vector_get (params, 0);
  double Delta = gsl_vector_get (params, 1);
  
  if (Gamma < 0) return HUGE_VAL;
  if (Delta < 0) return HUGE_VAL;
  return squared_residuals(params, data);
}

int
simplex(const double Gamma_init, 
        const double Delta_init, 
        struct data * d, 
        double * Gamma_best, 
        double * Delta_best)
/* Inspired by http://www.gnu.org/software/gsl/manual/html_node/Multimin-Examples.html */
{
  const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex;
  /* const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2; */
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  size_t iter = 0;
  int status;
  double size;
  
    /* Starting point */
  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, Gamma_init);
  gsl_vector_set (x, 1, Delta_init);
  
  /* Set initial step sizes */
  ss = gsl_vector_alloc (2);
  /* In practice, these ones seem the best... */
  gsl_vector_set (ss, 0, Gamma_init/7. + 0.04);
  gsl_vector_set (ss, 1, Delta_init/7. + 0.04);
  
  /* Initialize method and iterate */
  minex_func.n = 2;
  minex_func.f = &squared_residuals_w_constraints;
  minex_func.params = (void *) d;
  
  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
  printf ("%5d %+.8f %+.8f chi^2 = %011.8f size =  ---\n", 
           (int) iter,
           gsl_vector_get (s->x, 0), 
           gsl_vector_get (s->x, 1), 
           squared_residuals (x, d)
         ); 
  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status) 
        break;
      
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-4);
      
      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }
      * Gamma_best = gsl_vector_get (s->x, 0);
      * Delta_best = gsl_vector_get (s->x, 1);
      
      printf ("%5d %+.8f %+.8f chi^2 = %011.8f size = %.8f\n", 
              (int) iter,
              * Gamma_best, 
              * Delta_best, 
              s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 100);
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
  
  return status;
  
}


int
residuals_jacobian_df(const gsl_vector * params, void * data, gsl_matrix * J)
/* strongly based on http://www.gnu.org/software/gsl/manual/html_node/
 * Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html        
 * (expb_df function)
 */
{
  size_t n        = ((struct data *)data)->n;
  double *X       = ((struct data *)data)->X;
  double *Y       = ((struct data *)data)->Y;
  double *sigmaY  = ((struct data *)data)->sigmaY;
  
  double Gamma = gsl_vector_get (params, 0);
  double Delta = gsl_vector_get (params, 1);
  
  size_t i;
  
  double Ji0, Ji1;
  
  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dparamsj,    */
      /* where fi = (theory - experiment)/sigmaY[i], */
      /* and paramsj = Gamma, Delta */
      
      Ji0 = ( dGin_dGamma (X[i], Gamma, Delta, T0) - Y[i] ) / sigmaY[i];
      Ji1 = ( dGin_dDelta (X[i], Gamma, Delta, T0) - Y[i] ) / sigmaY[i];
      
      gsl_matrix_set (J, i, 0, Ji0); 
      gsl_matrix_set (J, i, 1, Ji1);
    }

  return GSL_SUCCESS;  
}

int
residuals_fdf(const gsl_vector * params, 
              void * data,
              gsl_vector * f,
              gsl_matrix * J)
/* strongly based on http://www.gnu.org/software/gsl/manual/html_node/
 * Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html        
 * (expb_fdf function)
 */                     
{
  residuals_vector_f(params, data, f);
  residuals_jacobian_df(params, data, J);
  
  return GSL_SUCCESS;
}

     
void
print_fit_state (size_t iter, gsl_multifit_fdfsolver * s, size_t DoF)
{
  printf ("iter: %3u // Gamma = %.8f Delta = %.8f"
               " chi^2 = %.10f chi^2/DoF = %.10f\n",
               (unsigned int) iter,
               gsl_vector_get(s->x, 0), 
               gsl_vector_get(s->x, 1),
               gsl_pow_2(gsl_blas_dnrm2(s->f)),
               gsl_pow_2(gsl_blas_dnrm2(s->f))/DoF
          );
}

int 
fit(struct data * d,
    const double Gamma_init,
    const double Delta_init,
    double *Gamma_best, 
    double *Delta_best, 
    gsl_matrix *cov,
    double * reduced_chi_square)
/* strongly based on http://www.gnu.org/software/gsl/manual/html_node/
 * Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html        
 * (main function)
 */ 
{ 
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int iter = 0;
  const size_t n = d->n;
  const size_t p = 2;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  gsl_multifit_function_fdf f;
  double x_init[2]; /* initial guess for (Gamma, Delta) */
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  
  x_init[0] = Gamma_init;
  x_init[1] = Delta_init;
  
  f.f   = &residuals_vector_f;
  f.df  = &residuals_jacobian_df;
  f.fdf  = &residuals_fdf;
  f.n = n;
  f.p = p;
  f.params = d; /* here 'params' are actually the experiment. data, not Gamma and Delta */
  
  T = gsl_multifit_fdfsolver_lmsder;
  /* T = gsl_multifit_fdfsolver_lmder; */
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);  
  
  print_fit_state (iter, s, n);
  
  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
     
      /* printf ("status = %s\n", gsl_strerror (status)); */
     
      print_fit_state (iter, s, n);
     
      /* if (status) break; */
     
      status = gsl_multifit_test_delta (s->dx, s->x, 0, 1e-6); 
    }
  while (status == GSL_CONTINUE &&  iter < 500);
  
  gsl_multifit_covar (s->J, 0.0, covar);
  
  gsl_matrix_memcpy(cov, covar);
  
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))  
#define COV    gsl_matrix_get(covar,0,1)
#define COV1   gsl_matrix_get(covar,1,0)
  
  assert(COV==COV1); /* covariance matrix is simmetric! */
  
    { 
      double chi = gsl_blas_dnrm2(s->f);
      double dof = n - p;
      double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
      
      printf("\n\n  ***** FINAL RESULTS: *****   \n");
      

      *reduced_chi_square = pow(chi, 2.0) / dof;

      printf("chi^2/DoF = %g\n", *reduced_chi_square);
      
      *Gamma_best = FIT(0);
      *Delta_best = FIT(1);
      
      printf ("Gamma = %.5f +/- %.5f\n", *Gamma_best, c*ERR(0));
      printf ("Delta = %.5f +/- %.5f\n", *Delta_best, c*ERR(1));
      printf ("cov(Gamma, Delta) = %.5f\n", COV);
    }
  
  /* printf ("status = %s\n", gsl_strerror (status)); */
  
  printf("\n");
  
  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
        
  return GSL_SUCCESS;
}
