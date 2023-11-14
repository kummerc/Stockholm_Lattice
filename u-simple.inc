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

//==================================================================================================
// Specialized gauge functions: simple version
//==================================================================================================

#include "defs.h"
#include <math.h>
#include "u.h"

//--------------------------------------------------------------------------------------------------
//! Initialize with zero
//--------------------------------------------------------------------------------------------------

void u_zero(SU3* u)
{
  int k, l;

  for (k = 0; k < NCOL; k++)
    for (l = 0; l < NCOL; l++)
      u->c[k][l] = 0.0;
}

//--------------------------------------------------------------------------------------------------
//! Generate unit SU(3) matrix
//--------------------------------------------------------------------------------------------------

void u_one(SU3* u)
{
  int k, l;

  for (k = 0; k < NCOL; k++)
    for (l = 0; l < NCOL; l++)
      u->c[k][l] = (k == l ? 1.0 : 0.0);
}

//--------------------------------------------------------------------------------------------------
//! Generate random SU(3) matrix
//--------------------------------------------------------------------------------------------------

void u_rng(SU3* u)
{
  int k, l;
  double phi, rho;

  for (k = 0; k < NCOL-1; k++)
    for (l = 0; l < NCOL; l++)
    {
      phi = 2. * M_PI * rng();
      rho = sqrt(-1.0 * log(rng()));
      u->c[k][l] = rho * (cos(phi) + _Complex_I * sin(phi));
    }

  u_reunitarise(u);

#if 0	// Testing
  double complex det = u_det(u);
  printf("u00 = (%e,%e)\n", creal(u->c[0][0]), cimag(u->c[0][0]));
  printf("det = (%e,%e)\n", creal(det), cimag(det));
  SU3 ud, uud;
  u_copy(&ud, u);
  u_dagger(&ud);
  u_mul(&uud, &ud, u);
  for (k = 0; k < NCOL; k++)
    for (l = 0; l < NCOL; l++)
      if ((k == l && cabs(uud.c[k][l]-1.) > 1.e-6) ||
          (k != l && cabs(uud.c[k][l]) > 1.e-6))
        die("%d %d %e %e\n", k, l, creal(uud.c[k][l]), cimag(uud.c[k][l]));
#endif
}

//--------------------------------------------------------------------------------------------------
//! Copy SU3 matrices
//--------------------------------------------------------------------------------------------------

void u_copy(SU3* o, SU3* i)
{
  int k, l;

  for (k = 0; k < NCOL; k++)
    for (l = 0; l < NCOL; l++)
      o->c[k][l] = i->c[k][l];
}

//--------------------------------------------------------------------------------------------------
//! Accumulate SU3 matrices
//--------------------------------------------------------------------------------------------------

void u_accum(SU3* io, SU3* i)
{
  int k, l;

  for (k = 0; k < NCOL; k++)
    for (l = 0; l < NCOL; l++)
      io->c[k][l] += i->c[k][l];
}

//--------------------------------------------------------------------------------------------------
//! Multiply SU3 matrices
//--------------------------------------------------------------------------------------------------

void u_mul(SU3* o, SU3* i1, SU3* i2)
{
  int k, l;

  for (k = 0; k < NCOL; k++)
    for (l = 0; l < NCOL; l++)
    {
      o->c[k][l] = i1->c[k][0] * i2->c[0][l] +
                   i1->c[k][1] * i2->c[1][l] +
                   i1->c[k][2] * i2->c[2][l];
    }
}

//--------------------------------------------------------------------------------------------------
//! Replace input by hermitian conjugate
//--------------------------------------------------------------------------------------------------

void u_dagger(SU3* u)
{
  int k, l;
  double complex t0, t1, t2;

  t0 = conj(u->c[0][1]);
  t1 = conj(u->c[0][2]);
  t2 = conj(u->c[1][2]);

  u->c[0][0] = conj(u->c[0][0]);
  u->c[0][1] = conj(u->c[1][0]);
  u->c[0][2] = conj(u->c[2][0]);
  u->c[1][0] = t0;
  u->c[1][1] = conj(u->c[1][1]);
  u->c[1][2] = conj(u->c[2][1]);
  u->c[2][0] = t1;
  u->c[2][1] = t2;
  u->c[2][2] = conj(u->c[2][2]);
}

//--------------------------------------------------------------------------------------------------
//! Normalize given row
//--------------------------------------------------------------------------------------------------

void u_norm_row(SU3* u, int row)
{
  double sum;
  int col;

  sum = 0.0;
  for (col = 0; col < NCOL; col++)
  {
    double complex z = u->c[row][col];
    sum += creal(z) * creal(z) + cimag(z) * cimag(z);
  }

  if ( sum <= 1.e-8)
    die("Normalization failed (sum=%e)\n", sum);

  sum = 1. / sqrt(sum);
  for (col = 0; col < NCOL; col++)
    u->c[row][col] *= sum;
}

//--------------------------------------------------------------------------------------------------
//! Orthogonalize two rows
//--------------------------------------------------------------------------------------------------

void u_orthog_rows(SU3* u, int row1, int row2)
{
  int col;
  double complex p1, p2, f;

  p1 = 0.0;
  p2 = 0.0;
  for (col = 0; col < NCOL; col++)
  {
    p1 += conj(u->c[row1][col]) * u->c[row2][col];
    p2 += conj(u->c[row1][col]) * u->c[row1][col];
  }

  f = p1 * (1.0 / creal(p2));

  for (col = 0; col < NCOL; col++)
    u->c[row2][col] -= f * u->c[row1][col];
}

//--------------------------------------------------------------------------------------------------
//! Compute ( row1 x row2 )*
//--------------------------------------------------------------------------------------------------

void u_cross_rows(SU3* u, int row1, int row2, int row3)
{
  int col;
  double complex t[3];

  t[0] = u->c[row1][1] * u->c[row2][2] - u->c[row2][1] * u->c[row1][2];
  t[1] = u->c[row2][0] * u->c[row1][2] - u->c[row1][0] * u->c[row2][2];
  t[2] = u->c[row1][0] * u->c[row2][1] - u->c[row2][0] * u->c[row1][1];

  for (col = 0; col < NCOL; col++)
    u->c[row3][col] = conj(t[col]);
}

//--------------------------------------------------------------------------------------------------
//! Re-unitarise SU(3) matrix
//--------------------------------------------------------------------------------------------------

void u_reunitarise(SU3* u)
{
  u_norm_row(u, 0);
  u_orthog_rows(u, 0, 1);
  u_norm_row(u, 1);
  u_cross_rows(u, 0, 1, 2);
}

//--------------------------------------------------------------------------------------------------
//! Compute determinant
//--------------------------------------------------------------------------------------------------

complex double u_det(SU3* u)
{
  return (u->c[0][0] * u->c[1][1] * u->c[2][2] +
          u->c[0][1] * u->c[1][2] * u->c[2][0] +
          u->c[0][2] * u->c[1][0] * u->c[2][1] -
          u->c[0][2] * u->c[1][1] * u->c[2][0] -
          u->c[0][1] * u->c[1][0] * u->c[2][2] -
          u->c[0][0] * u->c[1][2] * u->c[2][1]);
}
