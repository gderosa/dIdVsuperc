/* constraints.c 
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
constraints
(
  double Gamma1,
  double Gamma2,  
  double Delta1,
  double Delta2,
  double alpha1
)
{
  return (
    Gamma1  > CONSTRAINT_GAMMA1_MIN &&
    Gamma1  < CONSTRAINT_GAMMA1_MAX &&

    Gamma2  > CONSTRAINT_GAMMA1_MIN &&
    Gamma2  < CONSTRAINT_GAMMA1_MAX &&

    Delta1  > CONSTRAINT_DELTA1_MIN &&
    Delta1  < CONSTRAINT_DELTA1_MAX &&

    Delta2  > CONSTRAINT_DELTA2_MIN &&
    Delta2  < CONSTRAINT_DELTA2_MAX &&

    alpha1  > CONSTRAINT_ALPHA1_MIN &&
    alpha1  < CONSTRAINT_ALPHA1_MAX 
  );     
}

