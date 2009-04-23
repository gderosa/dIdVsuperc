/* splot.c 3D - plot of chi squared (Gamma, Delta) 
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

/* 
 * 3D plot of [reduced] chi squared (Gamma, Delta), 
 */
int
splot(const double Gamma, const double Delta, const double DGamma, const double DDelta, struct data * d)
{
      double Gamma_var, Delta_var;
      char chi_filename[BUFSIZ];
      FILE * file;
      size_t i, j;
      gsl_vector * params = gsl_vector_alloc(2);

      double Gamma_start = Gamma - DGamma;
      double Gamma_step =  2*DGamma / (SPLOTSTEPS-1);

      double Delta_start = Delta - DDelta;
      double Delta_step =  2*DDelta / (SPLOTSTEPS-1);

      strcpy(chi_filename, ExpDataFile);
      strcat(chi_filename, ".chi");
      file = fopen(chi_filename, "w");

      printf("Saving data for plotting (chi^2(Delta, Gamma)/DoF)\n");
  
      for (i=0; i<SPLOTSTEPS; i++)
        {
          for (j=0; j<SPLOTSTEPS; j++)
            {
	            Gamma_var = Gamma_start + i*Gamma_step;
              Delta_var = Delta_start + j*Delta_step;
              
              gsl_vector_set(params, 0, Gamma_var);
              gsl_vector_set(params, 1, Delta_var);          
              fprintf(
                      file, 
                      "%.8f \t %.8f \t %.8f \n", 
                      Gamma_var, 
                      Delta_var, 
                      squared_residuals(params, d) / (d->n - 2) /* chi^2/DoF */
              );
	       
	            /* progress indicator (percent.) */
              printf(
                     "  in %s ... %03.1f%%\r",
		                 chi_filename,
                     ((double)(i*SPLOTSTEPS + j)/(SPLOTSTEPS*SPLOTSTEPS))*100.0
              );
              fflush(stdout);
             
            
	          }
          fprintf(file, "\n");
        }
      
      
      printf("  in %s ... done.  \n", chi_filename);
      
      fclose(file);
      gsl_vector_free(params);
      return 0;
}
