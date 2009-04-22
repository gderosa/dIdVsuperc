/* ui.c - interactive configuration of initial values and other parameters
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

/* just to avoid anooying warnings */
#ifndef fgets
#define fgets(line, n, stdin) line=fgets(line, n, stdin)
 
int
ui(void)
{
    char *line, *str;
    size_t n = BUFSIZ;
    
    /* Keep It Simple, Stupid ;-) */
    line = (char *) malloc(BUFSIZ*sizeof(char));
    str  = (char *) malloc(BUFSIZ*sizeof(char));
    
    printf("Initial value for Gamma                     [%g] ", Gamma0);
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Gamma0);
    
    printf("Initial value for Delta                     [%g] ", Delta0);
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Delta0);
    
    printf("Absolute temperature                        [%g] ", T0);
    fgets(line, n, stdin);
    sscanf(line, "%lf", &T0);
    
    printf("Experimental data file                      [%s] ", ExpDataFile);
    fgets(line, n, stdin);
    sscanf(line, "%s", ExpDataFile);

    printf("Fit interval: bias voltage lower limit(mV)  [unlimited] ");
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Vi);

    printf("Fit interval: bias voltage upper limit(mV)  [unlimited] ");
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Vf);

    /* 3D plotting of chi squared */
    printf(
          "Plot squared residuals (it may take a long time)?  [%s] ", 
           BOOL2yn(PlotSquareResiduals)
          );
    fgets(line, n, stdin);
    sscanf(line, "%s", str);
    PlotSquareResiduals = yn2BOOL(str);

    return 0;
}

#undef fgets
#endif



