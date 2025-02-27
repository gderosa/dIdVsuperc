/* main.c
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

/* used only in this file */
int write_results(
                  double Gamma1_best_final,                  
                  double sigma_Gamma1,
                  
                  double Gamma2_best_final,                  
                  double sigma_Gamma2,
                  
                  double Delta1_best_final,
                  double sigma_Delta1,
                                    
                  double Delta2_best_final,
                  double sigma_Delta2,
                                    
                  double alpha1_best_final,
                  double sigma_alpha1,
                  
                  double reduced_chi_squared,
                  
                  char * name
);

int 
main (void) 
{
  double ExpDataX[MaxExpPoints];
  double ExpDataY[MaxExpPoints];
  double ExpDataSigmaY[MaxExpPoints];
  unsigned int i = 0;
  
  double 
    Gamma1_best_step1,
    Gamma2_best_step1,      
    Delta1_best_step1,
    Delta2_best_step1,
    alpha1_best_step1;
  double 
     Gamma1_best_step2,
     Gamma2_best_step2,       
     Delta1_best_step2,
     Delta2_best_step2,
     alpha1_best_step2;
  double 
    sigma_Gamma1, 
    sigma_Gamma2,     
    sigma_Delta1,
    sigma_Delta2,
    sigma_alpha1;

  double reduced_chi_square;

  gsl_matrix * cov_Gamma1_Gamma2_Delta1_Delta2_alpha1 = gsl_matrix_alloc(5,5);
  
  struct data d;
  
  size_t ExpPoints;
  char line[BUFSIZ];
    
  FILE * f;

  init();
  
  printf("\nCopyright (C) 2008, 2009 Guido De Rosa <guidoderosa@gmail.com>\n"
         "Based on GNU Scientific Library http://www.gnu.org/software/gsl/ \n"
         "License: GPLv3\n\n"
         );  
  
  ui();
  
  f = fopen(ExpDataFile, "r");  
  
  if(!f) 
  	{
      fprintf(stderr,"Couldn't open %s!\n", ExpDataFile);
      return 127;
    }  
  while (fgets(line, sizeof line, f)) 
    {
      ExpDataX[i]       = strtod(strtok(line," \t\n\r"),NULL);
      ExpDataY[i]       = strtod(strtok(NULL," \t\n\r"),NULL);      
      ExpDataSigmaY[i]  = strtod(strtok(NULL," \t\n\r"),NULL);      
      i++;
    }          
  ExpPoints = i;    

  d.n       = ExpPoints;
  d.X       = ExpDataX;
  d.Y       = ExpDataY;
  d.sigmaY  = ExpDataSigmaY;
  
  printf("\n PASS 1: Minimization with Nelder-Mead SIMPLEX algorithm:\n");
  
  simplex(
    Gamma1_0,
    Gamma2_0,      
    Delta1_0,
    Delta2_0,
    alpha1_0,
    
    &d, 
    
    &Gamma1_best_step1,
    &Gamma2_best_step1,      
    &Delta1_best_step1,
    &Delta2_best_step1,
    &alpha1_best_step1
  );
  
  printf("\n PASS 2: (possibly) refining the fit; computing statistical errors\n" 
         "(Levenberg-Marquardt solver)\n\n");
  
  fit(
    &d, 

    Gamma1_best_step1,
    Gamma2_best_step1,      
    Delta1_best_step1, 
    Delta2_best_step1,
    alpha1_best_step1,
    
    &Gamma1_best_step2,
    &Gamma2_best_step2,      
    &Delta1_best_step2, 
    &Delta2_best_step2,
    &alpha1_best_step2,

    cov_Gamma1_Gamma2_Delta1_Delta2_alpha1, 
    
    &reduced_chi_square
  );

  sigma_Gamma1  = 
	  sqrt(gsl_matrix_get(cov_Gamma1_Gamma2_Delta1_Delta2_alpha1,0,0)) * 
	  GSL_MAX_DBL(1,sqrt(reduced_chi_square));
  sigma_Gamma2  = 
	  sqrt(gsl_matrix_get(cov_Gamma1_Gamma2_Delta1_Delta2_alpha1,1,1)) * 
	  GSL_MAX_DBL(1,sqrt(reduced_chi_square));	  
  sigma_Delta1 = 
	  sqrt(gsl_matrix_get(cov_Gamma1_Gamma2_Delta1_Delta2_alpha1,2,2)) * 
	  GSL_MAX_DBL(1,sqrt(reduced_chi_square));
  sigma_Delta2 = 
	  sqrt(gsl_matrix_get(cov_Gamma1_Gamma2_Delta1_Delta2_alpha1,3,3)) * 
	  GSL_MAX_DBL(1,sqrt(reduced_chi_square));
  sigma_alpha1 = 
	  sqrt(gsl_matrix_get(cov_Gamma1_Gamma2_Delta1_Delta2_alpha1,4,4)) * 
	  GSL_MAX_DBL(1,sqrt(reduced_chi_square));

#define Gamma1_best_final Gamma1_best_step2
#define Gamma2_best_final Gamma2_best_step2    
#define Delta1_best_final Delta1_best_step2 
#define Delta2_best_final Delta2_best_step2 
#define alpha1_best_final alpha1_best_step2 

/* #define COV gsl_matrix_get(cov_Gamma_Delta, 1, 0) */
  /* DISABLED: there are 5*4/2 = 10 independent covariances now! */
  /* TODO: display properly covariances matrix ;-) */

  write_results(
	  Gamma1_best_final, 
	  sigma_Gamma1,
	   
	  Gamma2_best_final, 
	  sigma_Gamma2, 
	  
	  Delta1_best_final, 
	  sigma_Delta1,
	  
	  Delta2_best_final, 
	  sigma_Delta2,
	  
 	  alpha1_best_final, 
	  sigma_alpha1,
	  
	  reduced_chi_square,
	  
	  ExpDataFile
  );
  
  plot(
    Gamma1_best_final,
    Gamma2_best_final,      
    Delta1_best_final, 
    Delta2_best_final, 
    alpha1_best_final,  
    &d
  );

/* DISABLED: it doesn't make sense in a space of 4 parameters!
 *  if (PlotSquareResiduals) 
 *    {
 *      printf("\n");
 *      splot(Gamma_best_final, Delta_best_final, SPLOTSIGMAS*sigma_Gamma, SPLOTSIGMAS*sigma_Delta, &d);
 *    }
 */
 
  gsl_matrix_free(cov_Gamma1_Gamma2_Delta1_Delta2_alpha1);
 
  return 0;
}


int write_results(
  double Gamma1,
  double sigma_Gamma1,
  double Gamma2,
  double sigma_Gamma2,  
  double Delta1,
  double sigma_Delta1,
  double Delta2,
  double sigma_Delta2,
  double alpha1,
  double sigma_alpha1,
  double reduced_chi_squared,
  char * name
)
{
  FILE * f;
  char str[BUFSIZ];
  char tmpstr[BUFSIZ];

  strcpy(str, name);
  sprintf(tmpstr, "%d", Mode); 
  strcat(str, ".Mode");
  strcat(str, tmpstr);
  strcat(str, ".out");

  f = fopen(str, "w");
  
  fprintf(f, "Gamma1    = %f +/- %f\n", Gamma1,             sigma_Gamma1);
  fprintf(f, "Gamma2    = %f +/- %f\n", Gamma2 * (Mode==2), sigma_Gamma2);  
  fprintf(f, "Delta1    = %f +/- %f\n", Delta1,             sigma_Delta1);
  fprintf(f, "Delta2    = %f +/- %f\n", Delta2 * (Mode!=0), sigma_Delta2);
  fprintf(f, "alpha1    = %f +/- %f\n", alpha1 * (Mode!=0), sigma_alpha1);
  fprintf(f, "chi^2/DoF = %f \n", reduced_chi_squared);
  
  /* fprintf(f, "cov(Gamma, Delta) = %.8f \n", cov); */

  fclose(f);

  return 0;
}


