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
  
  double Gamma1 = gsl_vector_get (params, 0);
  double Gamma2 = gsl_vector_get (params, 1);  
  double Delta1 = gsl_vector_get (params, 2);
  double Delta2 = gsl_vector_get (params, 3);
  double alpha1 = gsl_vector_get (params, 4); /* alpha1 + alpha2 = 1 */
  
  double theory, experiment, residual;
  
  size_t i;
  
  gsl_vector_set_all(f, 0.0);
  
  for (i = 0; i < n; i++)
    {
      /* fitting an interval only */ 
        /* might be Vi=-HUGE_VAL and/or Vf=HUGE_VAL for "unlimited" */
      if (X[i]<Vi || X[i]>Vf) continue; 
      
      theory = Gin_doubleDeltaGamma(X[i], Gamma1, Gamma2, Delta1, Delta2, alpha1, T0);
      experiment = Y[i];
      residual = (theory - experiment)/sigmaY[i];
      if (!constraints(
        Gamma1,
        Gamma2,
        Delta1,
        Delta2,
        alpha1
      )) residual = residual * HUGE_VAL; 

      gsl_vector_set (f, i, residual);
    }  
  return GSL_SUCCESS;
}

/* squared residuals ("chi squared") are first minimized through a Simplex algorithm,
 * before being minimized by a derivatives-based nonlinear fit (which, most importantly,
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

double 
squared_residuals_w_constraints(const gsl_vector * params, void * data)
{
  double Gamma1 = gsl_vector_get (params, 0);
  double Gamma2 = gsl_vector_get (params, 1);  
  double Delta1 = gsl_vector_get (params, 2);
  double Delta2 = gsl_vector_get (params, 3);
  double alpha1 = gsl_vector_get (params, 4);

  if (!constraints(
    Gamma1,
    Gamma2,
    Delta1,
    Delta2,
    alpha1
  )) return HUGE_VAL;
  
  return squared_residuals(params, data);
}

int
simplex(
  const double Gamma1_init,
  const double Gamma2_init,  
  const double Delta1_init,
  const double Delta2_init,  
  const double alpha1_init,
  struct data * d, 
  double * Gamma1_best,
  double * Gamma2_best,    
  double * Delta1_best,
  double * Delta2_best,
  double * alpha1_best)
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
  
  size_t DoF = d->n - 5;
  
  /* Starting point */
  x = gsl_vector_alloc (5);  
  gsl_vector_set (x, 0, Gamma1_init);
  gsl_vector_set (x, 1, Gamma2_init);  
  gsl_vector_set (x, 2, Delta1_init);
  gsl_vector_set (x, 3, Delta2_init);
  gsl_vector_set (x, 4, alpha1_init);  
  
  /* Set initial step sizes */
  ss = gsl_vector_alloc (5);
  /* In practice, these ones seem the best... */
  gsl_vector_set (ss, 0, Gamma1_init/5. + 0.04);
  gsl_vector_set (ss, 1, Gamma2_init/5. + 0.04);  
  gsl_vector_set (ss, 2, Delta1_init/5. + 0.04);
  gsl_vector_set (ss, 3, Delta2_init/5. + 0.04);
  gsl_vector_set (ss, 4, 0.1);                    /* related to alpha1 */
  
  /* Initialize method and iterate */
  minex_func.n = 5;
  minex_func.f = &squared_residuals_w_constraints;
  minex_func.params = (void *) d;
  
  s = gsl_multimin_fminimizer_alloc (T, 5);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
  printf ("iter=%3d Gamma1=%+.8f Gamma2=%+.8f Delta1=%+.8f Delta2=%+.8f alpha1=%+.8f chi^2/DoF=%011.8f simplex_size=---\n", 
  
           (int) iter,
           
           gsl_vector_get (s->x, 0), 
           gsl_vector_get (s->x, 1), 
           gsl_vector_get (s->x, 2), 
           gsl_vector_get (s->x, 3),            
           gsl_vector_get (s->x, 4),            
                      
           squared_residuals (x, d) / DoF
         ); 
  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status) 
        break;
      
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, MULTIMIN_TEST_SIZE);
      
      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }
      * Gamma1_best = gsl_vector_get (s->x, 0);
      * Gamma2_best = gsl_vector_get (s->x, 1);      
      * Delta1_best = gsl_vector_get (s->x, 2);
      * Delta2_best = gsl_vector_get (s->x, 3);
      * alpha1_best = gsl_vector_get (s->x, 4);
            
      printf ("iter=%3d Gamma1=%+.8f Gamma2=%+.8f Delta1=%+.8f Delta2=%+.8f alpha1=%+.8f chi^2/DoF=%011.8f simplex_size=%.8f\n", 
              (int) iter,
              * Gamma1_best, 
              * Gamma2_best,
              * Delta1_best,
              * Delta2_best, 
              * alpha1_best,                              
              s->fval / DoF, 
              size
      );
    }
  while (status == GSL_CONTINUE && iter < MAX_SIMPLEX_ITER);
  
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
  
  double Gamma1 = gsl_vector_get (params, 0);
  double Gamma2 = gsl_vector_get (params, 1);  
  double Delta1 = gsl_vector_get (params, 2);
  double Delta2 = gsl_vector_get (params, 3);
  double alpha1 = gsl_vector_get (params, 4);
      
  size_t i;
  
  double Ji0, Ji1, Ji2, Ji3, Ji4;
  
  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dparamsj,    */
      /* where fi = (theory - experiment)/sigmaY[i], */
      /* and paramsj = Gamma, Delta */
      
      Ji0 = ( dGin_doubleDeltaGamma_dGamma1(X[i], Gamma1, Gamma2, Delta1, Delta2, alpha1, T0) - Y[i] ) / sigmaY[i];
      Ji1 = ( dGin_doubleDeltaGamma_dGamma2(X[i], Gamma1, Gamma2, Delta1, Delta2, alpha1, T0) - Y[i] ) / sigmaY[i];      
      Ji2 = ( dGin_doubleDeltaGamma_dDelta1(X[i], Gamma1, Gamma2, Delta1, Delta2, alpha1, T0) - Y[i] ) / sigmaY[i];
      Ji3 = ( dGin_doubleDeltaGamma_dDelta2(X[i], Gamma1, Gamma2, Delta1, Delta2, alpha1, T0) - Y[i] ) / sigmaY[i];
      Ji4 = ( dGin_doubleDeltaGamma_dalpha1(X[i], Gamma1, Gamma2, Delta1, Delta2, alpha1, T0) - Y[i] ) / sigmaY[i];

      gsl_matrix_set (J, i, 0, Ji0); 
      gsl_matrix_set (J, i, 1, Ji1);
      gsl_matrix_set (J, i, 2, Ji2); 
      gsl_matrix_set (J, i, 3, Ji3);
      gsl_matrix_set (J, i, 4, Ji4);            
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
  printf ("iter=%3d Gamma1=%+.8f Gamma2=%+.8f Delta1=%+.8f Delta2=%+.8f alpha1=%+.8f"
               " chi^2=%.10f chi^2/DoF=%.10f\n",
               
               (unsigned int) iter,
               
               gsl_vector_get(s->x, 0), 
               gsl_vector_get(s->x, 1),
               gsl_vector_get(s->x, 2), 
               gsl_vector_get(s->x, 3),               
               gsl_vector_get(s->x, 4),               
                              
               gsl_pow_2(gsl_blas_dnrm2(s->f)),
               gsl_pow_2(gsl_blas_dnrm2(s->f))/DoF
          );
}

int 
fit(struct data * d,
    const double Gamma1_init,
    const double Gamma2_init,
    const double Delta1_init,
    const double Delta2_init,
    const double alpha1_init,        
    double *Gamma1_best,
    double *Gamma2_best,  
    double *Delta1_best,
    double *Delta2_best, 
    double *alpha1_best,          
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
  const size_t p = 5;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  gsl_multifit_function_fdf f;
  double x_init[5]; /* initial guess for (Gamma1, Gamma2, Delta1, Delta2, alpha1) */
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  
  x_init[0] = Gamma1_init;
  x_init[1] = Gamma2_init;  
  x_init[2] = Delta1_init;
  x_init[3] = Delta2_init;
  x_init[4] = alpha1_init;    
  
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
  
  print_fit_state (iter, s, n-p); /* n-p = DoF */
  
  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);
     
      /* printf ("status = %s\n", gsl_strerror (status)); */
     
      print_fit_state (iter, s, n);
     
      /* if (status) break; */
     
      status = gsl_multifit_test_delta (s->dx, s->x, 0, MULTIFIT_TEST_DELTA); 
    }
  while (status == GSL_CONTINUE &&  iter < MAX_FIT_ITER);
  
  gsl_multifit_covar (s->J, 0.0, covar);
  
  gsl_matrix_memcpy(cov, covar);
  
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))  
/* now parameters are 5, not just 2...
 *
 * #define COV    gsl_matrix_get(covar,0,1)
 * #define COV1   gsl_matrix_get(covar,1,0)
 *
 * assert(COV==COV1);
 */
  
    { 
      double chi = gsl_blas_dnrm2(s->f);
      double dof = n - p;
      double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
      
      printf("\n\n  ***** FINAL RESULTS: *****   \n");
      

      *reduced_chi_square = pow(chi, 2.0) / dof;

      printf("chi^2/DoF = %g\n", *reduced_chi_square);
      
      *Gamma1_best = FIT(0);
      *Gamma2_best = FIT(1);      
      *Delta1_best = FIT(2);
      *Delta2_best = FIT(3);
      *alpha1_best = FIT(4);
                  
      printf ("Gamma1 = %.5f +/- %.5f\n", *Gamma1_best, c*ERR(0));
      printf ("Gamma2 = %.5f +/- %.5f\n", *Gamma2_best, c*ERR(1));      
      printf ("Delta1 = %.5f +/- %.5f\n", *Delta1_best, c*ERR(2));
      printf ("Delta2 = %.5f +/- %.5f\n", *Delta2_best, c*ERR(3));
      printf ("alpha1 = %.5f +/- %.5f\n", *alpha1_best, c*ERR(4));            
      /* printf ("cov(Gamma, Delta) = %.8f\n", COV); */ 
        /* now there are actually 5*4/2 = 10 independent covariances! */
    }
  
  /* printf ("status = %s\n", gsl_strerror (status)); */
  
  printf("\n");
  
  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
        
  return GSL_SUCCESS;
}

int
constraints
(
  double Gamma1,
  double Gamma2,  
  double Delta1,
  double Delta2,
  double alpha1
)
{
  return (
    Gamma1  > CONSTRAINT_GAMMA1_MIN &&
    Gamma1  < CONSTRAINT_GAMMA1_MAX &&

    Gamma2  > CONSTRAINT_GAMMA1_MIN &&
    Gamma2  < CONSTRAINT_GAMMA1_MAX &&

    Delta1  > CONSTRAINT_DELTA1_MIN &&
    Delta1  < CONSTRAINT_DELTA1_MAX &&

    Delta2  > CONSTRAINT_DELTA2_MIN &&
    Delta2  < CONSTRAINT_DELTA2_MAX &&

    alpha1  > CONSTRAINT_ALPHA1_MIN &&
    alpha1  < CONSTRAINT_ALPHA1_MAX 
  );     
}

