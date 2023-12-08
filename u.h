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

#ifndef U_H
#define U_H

// Include selected implementation header
#if defined(U_SIMPLE)
#include "u-simple.h"
#else
#error "No implementation selected"
#endif

// Global variables
extern SU3 u[NLINK];

// Function declarations (general)
void u_init(void);
double u_plaq(void);
double u_sweep_metro(void);
void u_metro_offer(SU3*, SU3*);
int u_metro_accept(SU3*, SU3*, SU3*);

// Function declarations (specialized)
void u_zero(SU3*);
void u_one(SU3*);
void u_rng(SU3*);
void u_copy(SU3*, SU3*);
void u_accum(SU3*, SU3*);
void u_mul(SU3*, SU3*, SU3*);
void u_dagger(SU3*);
void u_norm_row(SU3*, int);
void u_orthog_rows(SU3*, int, int);
void u_cross_rows(SU3*, int, int, int);
void u_reunitarise(SU3*);
std::complex<double> u_det(SU3*);

#endif
