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
splot(const double Gamma, const double Delta, struct data * d)
{
      double Gamma_var, Delta_var;
      char chi_filename[BUFSIZ];
      FILE * file;
      size_t i = 0;
      gsl_vector * params = gsl_vector_alloc(2);

      strcpy(chi_filename, ExpDataFile);
      strcat(chi_filename, ".chi");
      file = fopen(chi_filename, "w");

      printf("Saving data for plotting (chi^2(Delta, Gamma))\n");
  
      for (Gamma_var=0; Gamma_var<2*Gamma; Gamma_var += 2*(Gamma/SPLOTSTEPS))
        {
          for (Delta_var=0; Delta_var<2*Delta; Delta_var += 2*(Delta/SPLOTSTEPS))
            {
              i++;
              printf(
                     "  in %s ... %03.1f%%\r",
		     chi_filename,
                     (100.0*((double)i))/(double)(SPLOTSTEPS*SPLOTSTEPS)
                     ); /* progress indicator (percent.) */
              
              gsl_vector_set(params, 0, Gamma_var);
              gsl_vector_set(params, 1, Delta_var);          
              fprintf(
                      file, 
                      "%.8f \t %.8f \t %.8f \n", 
                      Gamma_var, 
                      Delta_var, 
                      squared_residuals (params, d)                  
                      );
            }
          fprintf(file, "\n");
        }
      
      fclose(file);
      gsl_vector_free(params);
      return 0;
}
