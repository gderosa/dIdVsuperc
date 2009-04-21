/* plot.c
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

int 
plot(const double Gamma, const double Delta, struct data * d)
{
  double V, V_plot_i, V_plot_f;
  double V_exp_min = min(d->X,d->n); 
  double V_exp_max = max(d->X,d->n);
  double exp_data_xrange = V_exp_max - V_exp_min ;
  double Gamma_var, Delta_var;
  FILE * file;
  char fit_filename[BUFSIZ]; /* fit function datafile for plotting */
  gsl_vector * params = gsl_vector_alloc(2);  
  unsigned int i = 0;
  
  strcpy(fit_filename, ExpDataFile);
  strcat(fit_filename, ".fit");
  file = fopen(fit_filename, "w");

  printf("Saving data for plotting (dI/dV) in %s ... ", fit_filename);

  /* Plot theor. function "slightly larger" than experimental points */
  V_plot_i = V_exp_min - exp_data_xrange*ExtraPlotRatio;
  V_plot_f = V_exp_max + exp_data_xrange*ExtraPlotRatio; 
  for (V=V_plot_i; V<V_plot_f; V+=0.05) 
    {
      fprintf(file, "%.8f \t %.8f \n", V, Gin(V, Gamma, Delta, T0));       
    }  

  fclose(file);
  
  printf("done.\n");
  
  if (PlotSquareResiduals)    
    {
      file = fopen(SquaredResidualsPlotFile, "w");
  
      for (Gamma_var=0; Gamma_var<2*Gamma; Gamma_var += 2*(Gamma/PLOTSTEPS))
        {
          for (Delta_var=0; Delta_var<2*Delta; Delta_var += 2*(Delta/PLOTSTEPS))
            {
              i++;
              printf(
                     "Saving data for plotting (chi^2(Delta, Gamma))... %03.1f%%\r",
                     (100.0*((double)i))/(double)(PLOTSTEPS*PLOTSTEPS)
                     );
              
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
    }
  
  
  gsl_vector_free(params);
  
  return 0;
}
