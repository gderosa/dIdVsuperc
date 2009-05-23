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
Gin_doubleDeltaGamma(
  const double V,
  const double Gamma1,
  const double Gamma2,  
  const double Delta1,
  const double Delta2,
  const double alpha1,
  const double T)
{
  return alpha1*Gin(V, Gamma1, Delta1, T) + (1-alpha1)*Gin(V, Gamma2, Delta2, T);
}

/* Gin as a function of Gamma1; Gamma2, V, Delta1, Delta2, alpha1, T as parametrs */
double 
Gin_doubleDeltaGamma_Gamma1(const double Gamma1, void * params)
{
  double * p = (double *) params;
  double V      = p[0];
  double Gamma2 = p[1];
  double Delta1 = p[2];
  double Delta2 = p[3];
  double alpha1 = p[4];
  double T      = p[5];
  return Gin_doubleDeltaGamma(V, Gamma1, Gamma2, Delta1, Delta2, alpha1, T);
}

double 
dGin_doubleDeltaGamma_dGamma1(
  const double V, 
  const double Gamma1,
  const double Gamma2,  
  const double Delta1, 
  const double Delta2,
  const double alpha1,
  const double T)
{
  gsl_function F;
  double params[6];
  double result, abserr;
  int status;
  
  params[0] = V;
  params[1] = Gamma2;
  params[2] = Delta1;
  params[3] = Delta2;
  params[4] = alpha1;
  params[5] = T;
  
  F.function  = &Gin_doubleDeltaGamma_Gamma1;
  F.params    = &params;
  
  status = gsl_deriv_central (&F, Gamma1, DSTEPSIZE, &result, &abserr);
  
  return result;
}

/* Gin as a function of Gamma2; Gamma1, V, Delta1, Delta2, alpha1, T as parameters */
double 
Gin_doubleDeltaGamma_Gamma2(const double Gamma2, void * params)
{
  double * p = (double *) params;
  double V      = p[0];
  double Gamma1 = p[1];
  double Delta1 = p[2];
  double Delta2 = p[3];
  double alpha1 = p[4];
  double T      = p[5];
  return Gin_doubleDeltaGamma(V, Gamma1, Gamma2, Delta1, Delta2, alpha1, T);
}

double 
dGin_doubleDeltaGamma_dGamma2(
  const double V, 
  const double Gamma1,
  const double Gamma2,  
  const double Delta1, 
  const double Delta2,
  const double alpha1,
  const double T)
{
  gsl_function F;
  double params[6];
  double result, abserr;
  int status;
  
  params[0] = V;
  params[1] = Gamma1;
  params[2] = Delta1;
  params[3] = Delta2;
  params[4] = alpha1;
  params[5] = T;
  
  F.function  = &Gin_doubleDeltaGamma_Gamma2;
  F.params    = &params;
  
  status = gsl_deriv_central (&F, Gamma2, DSTEPSIZE, &result, &abserr);
  
  return result;
}

/* Gin as a function of Delta1; V, Gamma1, Gamma2, Delta2, alpha1, T as parametrs */
double 
Gin_doubleDeltaGamma_Delta1(const double Delta1, void * params)
{
  double * p = (double *) params;
  double V      = p[0];
  double Gamma1 = p[1];
  double Gamma2 = p[2];  
  double Delta2 = p[3];
  double alpha1 = p[4];
  double T      = p[5];
  return Gin_doubleDeltaGamma(V, Gamma1, Gamma2, Delta1, Delta2, alpha1, T);
}

double 
dGin_doubleDeltaGamma_dDelta1(
  const double V, 
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2,
  const double alpha1,
  const double T)
{
  gsl_function F;
  double params[6];
  double result, abserr;
  int status;
  
  params[0] = V;
  params[1] = Gamma1;
  params[2] = Gamma2;  
  params[3] = Delta2;
  params[4] = alpha1;
  params[5] = T;
    
  F.function = &Gin_doubleDeltaGamma_Delta1;
  F.params = &params;
  
  status = gsl_deriv_central (&F, Delta1, DSTEPSIZE, &result, &abserr);
  
  return result;
}

double 
Gin_doubleDeltaGamma_Delta2(const double Delta2, void * params)
{
  double * p = (double *) params;
  double V      = p[0];
  double Gamma1 = p[1];
  double Gamma2 = p[2];  
  double Delta1 = p[3];
  double alpha1 = p[4];
  double T      = p[5];
  return Gin_doubleDeltaGamma(V, Gamma1, Gamma2, Delta1, Delta2, alpha1, T);
}

double 
dGin_doubleDeltaGamma_dDelta2(
  const double V, 
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2,
  const double alpha1,
  const double T)
{
  gsl_function F;
  double params[6];
  double result, abserr;
  int status;
  
  params[0] = V;
  params[1] = Gamma1;
  params[2] = Gamma2;  
  params[3] = Delta1;
  params[4] = alpha1;
  params[5] = T;
    
  F.function = &Gin_doubleDeltaGamma_Delta2;
  F.params = &params;
  
  status = gsl_deriv_central (&F, Delta2, DSTEPSIZE, &result, &abserr);
  
  return result;
}

/* Gin as a function of alpha1; V, Gamma1, Gamma2, Delta1, Delta2, T as parameters */
double 
Gin_doubleDeltaGamma_alpha1(const double alpha1, void * params)
{
  double * p = (double *) params;
  double V      = p[0];
  double Gamma1 = p[1];
  double Gamma2 = p[2];  
  double Delta1 = p[3];
  double Delta2 = p[4];
  double T      = p[5];
  return Gin_doubleDeltaGamma(V, Gamma1, Gamma2, Delta1, Delta2, alpha1, T);
}

double 
dGin_doubleDeltaGamma_dalpha1(
  const double V, 
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2,
  const double alpha1,
  const double T)
{
  gsl_function F;
  double params[6];
  double result, abserr;
  int status;
  
  params[0] = V;
  params[1] = Gamma1;
  params[2] = Gamma2;  
  params[3] = Delta1;
  params[4] = Delta2;
  params[5] = T;
    
  F.function = &Gin_doubleDeltaGamma_alpha1;
  F.params = &params;
  
  status = gsl_deriv_central (&F, alpha1, DSTEPSIZE, &result, &abserr);
  
  return result;
}
