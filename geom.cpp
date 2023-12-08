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

int nnp[V][4];
int nnm[V][4];

//--------------------------------------------------------------------------------------------------
//! Define size index (x is fastest, t is slowest running index)
//--------------------------------------------------------------------------------------------------

int site(int x, int y, int z, int t)
{
  assert((-1 <= x) && (x <= LS));
  assert((-1 <= y) && (y <= LS));
  assert((-1 <= z) && (z <= LS));
  assert((-1 <= t) && (t <= LT));

  // Enforce periodic boundary conditions
  x = (x+LS) % LS;
  y = (y+LS) % LS;
  z = (z+LS) % LS;
  t = (t+LT) % LT;

  // Return site index
  return (((t*(LS) + z) * (LS) + y) * (LS) + x);
}

//--------------------------------------------------------------------------------------------------
//! Define link index (mu is fastest running index)
//--------------------------------------------------------------------------------------------------

int link(int site, int mu)
{
  assert((0 <= site) && (site < V));
  assert((0 <= mu) && (mu < 4));

  return (4*site + mu);
}

//--------------------------------------------------------------------------------------------------
//! Initialize lattice geometry
//--------------------------------------------------------------------------------------------------

void geom_init(void)
{
  int x, y, z, t;

  // Initialize neighbour lists
  for (t = 0; t < LT; t++)
    for (z = 0; z < LS; z++)
      for (y = 0; y < LS; y++)
        for (x = 0; x < LS; x++)
        {
          int s = site(x, y, z, t);

          // Nearest neighbours in plus direction
          nnp[s][DX] = site(x+1, y,   z,   t);
          nnp[s][DY] = site(x,   y+1, z,   t);
          nnp[s][DZ] = site(x,   y,   z+1, t);
          nnp[s][DT] = site(x,   y,   z,   t+1);

          // Nearest neighbours in minus direction
          nnm[s][DX] = site(x-1, y,   z,   t);
          nnm[s][DY] = site(x,   y-1, z,   t);
          nnm[s][DZ] = site(x,   y,   z-1, t);
          nnm[s][DT] = site(x,   y,   z,   t-1);
        }
}
