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

/* just to avoid annoying warnings */
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
    
    printf("Initial value for Gamma1                    [%g] ", Gamma1_0);
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Gamma1_0);

    printf("Initial value for Gamma2                    [%g] ", Gamma2_0);
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Gamma2_0);
    
    printf("Initial value for Delta1                    [%g] ", Delta1_0);
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Delta1_0);

    printf("Initial value for Delta2                    [%g] ", Delta2_0);
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Delta2_0);
 
    printf("Initial value for alpha1                    [%g] ", alpha1_0);
    fgets(line, n, stdin);
    sscanf(line, "%lf", &alpha1_0);
    
    printf("Absolute temperature                        [%g] ", T0);
    fgets(line, n, stdin);
    sscanf(line, "%lf", &T0);
    
    printf("Experimental data file                      [%s] ", ExpDataFile);
    fgets(line, n, stdin);
    sscanf(line, "%s", ExpDataFile);

    /* DISABLED, waiting for a better implementation
    printf("Fit interval: bias voltage lower limit(mV)  [unlimited] ");
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Vi);
    
    printf("Fit interval: bias voltage upper limit(mV)  [unlimited] ");
    fgets(line, n, stdin);
    sscanf(line, "%lf", &Vf);
    */

    return 0;
}

#undef fgets
#endif



