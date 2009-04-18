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
    Gamma0      = Gamma0_default;
    Delta0      = Delta0_default;
    k_B         = k_B_default;
    strcpy(ExpDataFile, ExpDataFile_default);
    PlotSquareResiduals = PlotSquareResiduals_default;
    
    return 0;
  }
