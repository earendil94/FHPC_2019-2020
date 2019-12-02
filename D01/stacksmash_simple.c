
/* ────────────────────────────────────────────────────────────────────────── *
 │                                                                            │
 │ This file is part of the exercises for the Lectures on                     │
 │   "Foundations of High Performance Computing"                              │
 │ given at                                                                   │
 │   Master in HPC and                                                        │
 │   Master in Data Science and Scientific Computing                          │
 │ @ SISSA, ICTP and University of Trieste                                    │
 │                                                                            │
 │ contact: luca.tornatore@inaf.it                                            │
 │                                                                            │
 │     This is free software; you can redistribute it and/or modify           │
 │     it under the terms of the GNU General Public License as published by   │
 │     the Free Software Foundation; either version 3 of the License, or      │
 │     (at your option) any later version.                                    │
 │     This code is distributed in the hope that it will be useful,           │
 │     but WITHOUT ANY WARRANTY; without even the implied warranty of         │
 │     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          │
 │     GNU General Public License for more details.                           │
 │                                                                            │
 │     You should have received a copy of the GNU General Public License      │
 │     along with this program.  If not, see <http://www.gnu.org/licenses/>   │
 │                                                                            │
 * ────────────────────────────────────────────────────────────────────────── */

/* ----------------------------------------------------------------- *
 *                                                                   *
 * a very basic example of stack smashing                            *
 *                                                                   *
 * ----------------------------------------------------------------- */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>



int main( void )
{

  int array[ 3 ];

  //Trying to access to array[3] which does not exist in the array defined in the stack. I am trying to access a memory location which is not in the stack! If you define other stuff in the stack, you can actually get out of the boundaries of your array, as long as you go in a memory not protected by you. 
  for ( int ii = 0; ii <= 3; ii++ )
    array[ ii ] = ii;

  return 0;
}
