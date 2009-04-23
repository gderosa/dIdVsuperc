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
                  double Gamma_best_final,
                  double sigma_Gamma,
                  double Delta_best_final,
                  double sigma_Delta,
                  double COV,
                  char * name);

int 
main (void) 
{
  double ExpDataX[MaxExpPoints];
  double ExpDataY[MaxExpPoints];
  double ExpDataSigmaY[MaxExpPoints];
  unsigned int i = 0;
  
  double Gamma_best_step1, Delta_best_step1;
  double Gamma_best_step2, Delta_best_step2;
  double sigma_Gamma, sigma_Delta;
  double reduced_chi_square;
  gsl_matrix * cov_Gamma_Delta = gsl_matrix_alloc(2,2);
  
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
  
  simplex(Gamma0, Delta0, &d, &Gamma_best_step1, &Delta_best_step1);
  
  printf("\n PASS 2: (possibly) refining the fit; computing statistical errors\n" 
         "(Levenberg-Marquardt solver)\n\n");
  
  fit(&d, Gamma_best_step1, Delta_best_step1, &Gamma_best_step2, &Delta_best_step2, cov_Gamma_Delta, &reduced_chi_square);

  sigma_Gamma = 
	  sqrt(gsl_matrix_get(cov_Gamma_Delta,0,0)) * 
	  GSL_MAX_DBL(1,sqrt(reduced_chi_square));
  sigma_Delta = 
	  sqrt(gsl_matrix_get(cov_Gamma_Delta,1,1)) * 
	  GSL_MAX_DBL(1,sqrt(reduced_chi_square));

#define Gamma_best_final Gamma_best_step2  
#define Delta_best_final Delta_best_step2 
#define COV gsl_matrix_get(cov_Gamma_Delta, 1, 0)

  write_results(
		  Gamma_best_final, 
		  sigma_Gamma, 
		  Delta_best_final, 
		  sigma_Delta,
		  COV,
		  ExpDataFile);
  
  plot(Gamma_best_final, Delta_best_final, &d);

  if (PlotSquareResiduals) 
    {
      printf("\n");
      splot(Gamma_best_final, Delta_best_final, SPLOTSIGMAS*sigma_Gamma, SPLOTSIGMAS*sigma_Delta, &d);
    }

  gsl_matrix_free(cov_Gamma_Delta);
  return 0;
}

int write_results(
                  double Gamma,
                  double sigma_Gamma,
                  double Delta,
                  double sigma_Delta,
                  double cov,
                  char * name)
{
  FILE * f;

  strcat(name,".out");
  
  f = fopen(name, "w");
  
  fprintf(f, "Gamma = %f +/- %f\n", Gamma, sigma_Gamma);
  fprintf(f, "Delta = %f +/- %f\n", Delta, sigma_Delta);
  fprintf(f, "cov(Gamma, Delta) = %.8f \n", cov);

  fclose(f);

  return 0;
}


