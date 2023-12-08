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

#ifndef DEFS_H
#define DEFS_H

#include <iostream>
#include <cstdlib>
#include <cstdarg>
#include <cassert>
#include <complex>

// Simulation parameters
#define LS		8			//!< Lattice size in spatial direction
#define LT		8			//!< Lattice size in temporal direction

#define BETA		6.0			//!< Gauge coupling constant
#define UINIT		cold			//!< hot or cold start
#define METRO_NSWEEP	50			//!< Number of Metropolis updates
#define METRO_NHIT	50			//!< Number of Metropolis hits
#define METRO_BIAS	3.0			//!< Bias used in Metropolis update

#define SEED		23984			//!< Random number generator seed

//--------------------------------------- NO CHANGES BELOW ---------------------------------------//

// Derived simulation parameters
#define V		((LS)*(LS)*(LS)*(LT))	//!< Volume
#define NLINK		((V)*4)			//!< Number of links

// Other constants
#define NCOL		3			//!< Number of colours (do not change)

// Enumerations
static enum {
  DX,						//!< x direction
  DY,						//!< y direction
  DZ,						//!< z direction
  DT						//!< t direction
} direction;

static enum {
  cold,						//!< Cold start
  hot						//!< Hot start
} uinit;

// Function declarations
void die_(const char *, int, const char *, ...);

// Error handling
#define die(ARGS...) die_(__FILE__, __LINE__, ##ARGS)

#endif
