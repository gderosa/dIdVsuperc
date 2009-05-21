/* misc.c - Miscellaneous handy/useful functions
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
 
/* Complex functions */
        
gsl_complex
cpow_2(const gsl_complex z)
{
  return gsl_complex_mul(z, z);
}

/* Minimum and maximum element of an array -- TODO: use gsl_vector instead */

double array_min(double *ary, size_t n) 
{
  double min = HUGE_VAL;
  size_t i;

  for (i=0;i<n;i++)
    {
      if (ary[i] < min) 
        min = ary[i];
    }

  return min;
}

double array_max(double *ary, size_t n) 
{
  double max = -HUGE_VAL;
  size_t i;

  for (i=0;i<n;i++)
    {
      if (ary[i] > max) 
        max = ary[i];
    }

  return max;
}

/* User Interface */

const char * 
BOOL2yn(BOOL b) 
{
  if (b) 
    return "Y/n";
  else
    return "y/N";
}

BOOL
yn2BOOL(char * str)
{
  if (str[0] == 'y' || str[0] == 'y')
    return 1;
  return 0;
}
