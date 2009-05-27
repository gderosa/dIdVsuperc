/* init.c
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
init(void) 
  {
    T0          = T0_default;
    Gamma1_0    = Gamma1_0_default;
    Gamma2_0    = Gamma2_0_default;    
    Delta1_0    = Delta1_0_default;
    Delta2_0    = Delta2_0_default;
    alpha1_0    = alpha1_0_default;
    k_B         = k_B_default;
    Vi      		= Vi_default;
    Vf		      = Vf_default;
    strcpy(ExpDataFile, ExpDataFile_default);
    Mode        = Mode_default;
    
    return 0;
  }
