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
plot(
  const double Gamma1,
  const double Gamma2,    
  const double Delta1, 
  const double Delta2, 
  const double alpha1, 
  struct data * d
)
{
  double V, V_plot_i, V_plot_f;
  double V_exp_min = array_min(d->X,d->n); 
  double V_exp_max = array_max(d->X,d->n);
  double exp_data_xrange = V_exp_max - V_exp_min ;
  FILE * file;
  char fit_filename[BUFSIZ]; /* fit function datafile for plotting */
  char tmpstr[BUFSIZ]; 
  unsigned int i;
  
  /* Why is it so difficult playing with strings in C ? */
  strcpy(fit_filename, ExpDataFile);
  strcat(fit_filename, ".Mode");
  sprintf(tmpstr, "%d", Mode); /* contains a number */
  strcat(fit_filename, tmpstr); 
  strcat(fit_filename, ".fit");

  file = fopen(fit_filename, "w");

  printf("Saving data for plotting (dI/dV)\n  in %s ... ", fit_filename);

  /* Plot theor. function "slightly larger" than experimental points */
  V_plot_i = V_exp_min - exp_data_xrange*ExtraPlotRatio;
  V_plot_f = V_exp_max + exp_data_xrange*ExtraPlotRatio;
  for (V=V_plot_i; V < d->X[0]; V+=0.05)
    {
      fprintf(
        file,
        "%.8f \t %.8f \n",
        V,
        Gin_doubleDeltaGamma(V, Gamma1, Gamma2, Delta1, Delta2, alpha1, T0)
      );
    }
  for (i = 0; i < d->n; i++) /* plot exactly AT the experimental points */
    {
      V = d->X[i];
      fprintf(
        file, 
        "%.8f \t %.8f \n", 
        V, 
        Gin_doubleDeltaGamma(V, Gamma1, Gamma2, Delta1, Delta2, alpha1, T0)
      );       
    }  
  for (V = d->X[d->n - 1]; V < V_plot_f; V+=0.05)
    {
      fprintf(
        file,
        "%.8f \t %.8f \n",
        V,
        Gin_doubleDeltaGamma(V, Gamma1, Gamma2, Delta1, Delta2, alpha1, T0)
      );
    }

  fclose(file);
  
  printf("done.\n");
  
  return 0;
}

