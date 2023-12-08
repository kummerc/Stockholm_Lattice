//
// Copyright (c) 2013 Forschungszentrum Juelich
//
// Author(s): Dirk Pleiter
//
// This software is available to you under a choice of one of two
// licenses.  You may choose to be licensed under the terms of the GNU
// General Public License (GPL) Version 2, available from the file
// COPYING in the main directory of this source tree, or the
// OpenIB.org BSD license below:
//
//     Redistribution and use in source and binary forms, with or
//     without modification, are permitted provided that the following
//     conditions are met:
//
//      - Redistributions of source code must retain the above
//        copyright notice, this list of conditions and the following
//        disclaimer.
//
//      - Redistributions in binary form must reproduce the above
//        copyright notice, this list of conditions and the following
//        disclaimer in the documentation and/or other materials
//        provided with the distribution.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//--------------------------------------------------------------------------------------------------

#include "defs.h"
#include "geom.h"
#include "rng.h"
#include "u.h"


//--------------------------------------------------------------------------------------------------
//! Main routine
//--------------------------------------------------------------------------------------------------


double t0_tot=0;
double t1_tot=0;
double t2_tot=0;
int main()
{
  int i;

  // Initialisation routines
  geom_init();
  rng_init();
  u_init();

  printf("# beta  = %f\n", BETA);
  printf("# uinit = %d\n", UINIT);
  printf("# nhit  = %d\n", METRO_NHIT);
  printf("#\n");
  printf("# update  plaquette  acceptance\n");

  // Update gauge fields and print plaquette
  for (i = 0; i < METRO_NSWEEP; i++)
  {
    double acc;
    acc = u_sweep_metro();
    printf("%6d     %.6e     %.2e\n", i, acc, u_plaq());
    fflush(stdout);
  }
  printf("t0_tot=%lf, t1_tot=%lf, t2_tot=%lf\n", t0_tot, t1_tot, t2_tot);

  return 0;
}
