/* functions.c - Physically relevant functions w/ their derivatives
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

#include "common.h"

/* Density of states */
double 
ni(const double E, const double Gamma, const double Delta)
{
  gsl_complex tmp = 
    gsl_complex_div
      (
        gsl_complex_rect(E, -Gamma) , 
        gsl_complex_sqrt
          (  
            gsl_complex_sub_real
              (
                cpow_2(gsl_complex_rect(E, -Gamma)), 
                gsl_pow_2(Delta) 
              )   
          )
      )
  ;
  
  return fabs(GSL_REAL(tmp));
}

/* Fermi distribution -- u is E/(k_B*T) physically */   
double
fermi(const double u)
{
  return 1. / ( exp(u) + 1. ) ;
}

/* Fermi distribution derivative * (-1) */   
double 
fprime(const double u)
{
  return gsl_pow_2( 1. / (2.*cosh(u/2.)) );
}

/* Differential conductance: integrand */
/* V and eV are identified, so voltage is in Volt and Energy in eV */
double
Gin_integrand_base(const double E,  const double V, const double Gamma, const double Delta, const double T)
{
  return ni(E, Gamma, Delta) * fprime((E-V)/(k_B*T));
}

/* Differential conductance: integrand GSL-ized */
double 
Gin_integrand(const double E, void * params)
{
  struct Gin_integrand_params * my_params = (struct Gin_integrand_params *) params;
  
  return Gin_integrand_base(E, my_params->V, my_params->Gamma, my_params->Delta, my_params->T);
}

double
Gin(const double V, const double Gamma, const double Delta, const double T)
{
  double integral_result, integral_error;
  struct Gin_integrand_params GinIntegrandParams;
  gsl_function F;
  gsl_integration_workspace * w;
  
  GinIntegrandParams.V      = V;
  GinIntegrandParams.Gamma  = Gamma;
  GinIntegrandParams.Delta  = Delta;
  GinIntegrandParams.T      = T;
  
  F.function = &Gin_integrand;
  F.params = &GinIntegrandParams;
    
  w = gsl_integration_workspace_alloc (SUBINTERVALS);
  
  /* gsl_integration_qags (&F, -40.*Delta, 40.*Delta, ABSTOL, RELTOL, SUBINTERVALS,
                        w, &integral_result, &integral_error); */
  
  /* DEBUG */ /*fprintf(stderr, "\nV=%.6f Gamma=%.6f Delta=%.6f T=%.6f\n", V, Gamma, Delta, T);*/  /* DEBUG */
  gsl_integration_qagi(&F, ABSTOL, RELTOL, SUBINTERVALS,
                        w, &integral_result, &integral_error); 
  
  gsl_integration_workspace_free(w);
                                
  return integral_result / (k_B * T);
}

double
Gin_doubleDelta(
  const double V,
  const double Gamma,
  const double Delta1,
  const double Delta2,
  const double alpha1,
  const double T)
{
  return alpha1*Gin(V, Gamma, Delta1, T) + (1-alpha1)*Gin(V, Gamma, Delta2, T);
}

/* Gin as a function of Gamma; V, Delta1, Delta2, alpha1, T as parametrs */
double 
Gin_doubleDelta_Gamma(const double Gamma, void * params)
{
  double * p = (double *) params;
  double V      = p[0];
  double Delta1 = p[1];
  double Delta2 = p[2];
  double alpha1 = p[3];
  double T      = p[4];
  return Gin_doubleDelta(V, Gamma, Delta1, Delta2, alpha1, T);
}

double 
dGin_doubleDelta_dGamma(
  const double V, 
  const double Gamma, 
  const double Delta1, 
  const double Delta2,
  const double alpha1,
  const double T)
{
  gsl_function F;
  double params[5];
  double result, abserr;
  int status;
  
  params[0] = V;
  params[1] = Delta1;
  params[2] = Delta2;
  params[3] = alpha1;
  params[4] = T;
  
  F.function  = &Gin_doubleDelta_Gamma;
  F.params    = &params;
  
  status = gsl_deriv_central (&F, Gamma, DSTEPSIZE, &result, &abserr);
  
  return result;
}

/* Gin as a function of Delta1; V, Gamma, Delta2, alpha1, T as parametrs */
double 
Gin_doubleDelta_Delta1(const double Delta1, void * params)
{
  double * p = (double *) params;
  double V      = p[0];
  double Gamma  = p[1];
  double Delta2 = p[2];
  double alpha1 = p[3];
  double T      = p[4];
  return Gin_doubleDelta(V, Gamma, Delta1, Delta2, alpha1, T);
}

double 
dGin_doubleDelta_dDelta1(
  const double V, 
  const double Gamma, 
  const double Delta1, 
  const double Delta2,
  const double alpha1,
  const double T)
{
  gsl_function F;
  double params[5];
  double result, abserr;
  int status;
  
  params[0] = V;
  params[1] = Gamma;
  params[2] = Delta2;
  params[3] = alpha1;
  params[4] = T;
    
  F.function = &Gin_doubleDelta_Delta1;
  F.params = &params;
  
  status = gsl_deriv_central (&F, Delta1, DSTEPSIZE, &result, &abserr);
  
  return result;
}

double 
Gin_doubleDelta_Delta2(const double Delta2, void * params)
{
  double * p = (double *) params;
  double V      = p[0];
  double Gamma  = p[1];
  double Delta1 = p[2];
  double alpha1 = p[3];
  double T      = p[4];
  return Gin_doubleDelta(V, Gamma, Delta1, Delta2, alpha1, T);
}

double 
dGin_doubleDelta_dDelta2(
  const double V, 
  const double Gamma, 
  const double Delta1, 
  const double Delta2,
  const double alpha1,
  const double T)
{
  gsl_function F;
  double params[5];
  double result, abserr;
  int status;
  
  params[0] = V;
  params[1] = Gamma;
  params[2] = Delta1;
  params[3] = alpha1;
  params[4] = T;
    
  F.function = &Gin_doubleDelta_Delta2;
  F.params = &params;
  
  status = gsl_deriv_central (&F, Delta2, DSTEPSIZE, &result, &abserr);
  
  return result;
}

/* Gin as a function of alpha1; V, Gamma, Delta1, Delta2, T as parameters */
double 
Gin_doubleDelta_alpha1(const double alpha1, void * params)
{
  double * p = (double *) params;
  double V      = p[0];
  double Gamma  = p[1];
  double Delta1 = p[2];
  double Delta2 = p[3];
  double T      = p[4];
  return Gin_doubleDelta(V, Gamma, Delta1, Delta2, alpha1, T);
}

double 
dGin_doubleDelta_dalpha1(
  const double V, 
  const double Gamma, 
  const double Delta1, 
  const double Delta2,
  const double alpha1,
  const double T)
{
  gsl_function F;
  double params[5];
  double result, abserr;
  int status;
  
  params[0] = V;
  params[1] = Gamma;
  params[2] = Delta1;
  params[3] = Delta2;
  params[4] = T;
    
  F.function = &Gin_doubleDelta_alpha1;
  F.params = &params;
  
  status = gsl_deriv_central (&F, alpha1, DSTEPSIZE, &result, &abserr);
  
  return result;
}
