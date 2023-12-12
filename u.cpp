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
#include "rng.h"
#include "geom.h"
#include "geom.cpp"
#include "u.h"
#include <time.h>
#include <sycl/sycl.hpp>

// Include specialized functions
#if defined(U_SIMPLE)
#include "u-simple.inc"
#else
#error "No implementation defined"
#endif

SU3 u[NLINK];

//--------------------------------------------------------------------------------------------------
//! Initialize gauge fields
//--------------------------------------------------------------------------------------------------

void u_init(void)
{
  int i;

  for (i = 0; i < NLINK; i++)
  {
    switch (UINIT)
    {
      case hot:
        u_rng(&u[i]);
        break;
      case cold:
        u_one(&u[i]);
        break;
      default:
        die_("Unknown initialisation (UINIT=%d)\n", UINIT);
    }
  }
}

//--------------------------------------------------------------------------------------------------
//! Compute plaquette: \f$\frac{1}{18 V} \mathrm{Re Tr}\sum_p W_p\f$ with
//! \f$\sum_p W_p = \sum_x \sum_\mu \sum_{\nu > \mu} U_{x,\mu} U_{x+\hat{\nu},\nu} U^\dagger_{x+\hat{\nu},\mu} U^\dagger_{x,\nu}\f$
//--------------------------------------------------------------------------------------------------

void u_plaq(void) {
  clock_t start_copy, end_copy, start_metro, end_metro;
    
  start_copy = clock();
  double plaq;
  
  // Create a buffer to store the result
  sycl::buffer<double, 1> plaqBuffer(sycl::range<1>(1));
  //sycl::buffer<double, 1> accBuffer(sycl::range<1>(1));
  
  // Create a SYCL queue to specify the device (e.g., GPU)
  sycl::queue queue(sycl::gpu_selector{});
   
  SU3 *ud = sycl::malloc_device<SU3>(4 * VOL, queue);
  int *nnpd = sycl::malloc_device<int>(4 * VOL, queue);
  //int *nnmd = sycl::malloc_device<int>(4 * VOL, queue);
  end_copy = clock();
  printf("Time u_copy_plaq(): %f s\n", ((double) (end_copy - start_copy)) / CLOCKS_PER_SEC);
  printf("##################################################################################\n");
  fflush(stdout);
  queue.copy<SU3>(u, ud, 4 * VOL);
  queue.copy<int>(&(nnp[0][0]), nnpd, 4 * VOL);
  //queue.copy<int>(&(nnm[0][0]), nnmd, 4 * VOL);
  
  // Submit a command group to the queue
  int i_METRO;
  double acc;
  for (i_METRO=0; i_METRO<METRO_NSWEEP; i_METRO++){
    start_metro = clock();
    
    //Metropolis update starts
    /*
    queue.submit([&](sycl::handler& cgm) {
      auto accAcc = accBuffer.get_access<sycl::access::mode::write>(cgm);
      
      cgm.parallel_for<class MetroKernel>(sycl::range<1>(LT * LS * LS * LS), sycl::reduction(accAcc, 0.0, std::plus<double>{}),[=](sycl::id<1> idx, auto& accLoc) {
        int t = idx / (LS * LS * LS);
        int z = (idx / (LS * LS)) % LS;
        int y = (idx / LS) % LS;
        int x = idx % LS;
        
        int s = site(x, y, z, t);
        double accd = 0.0;
         
        //int im, jm;
        int ihit;
        int mu, nu;
        int l;
        SU3 staple;
        int iacc = 0;
         
        for(int rem=0; rem<8; rem++){ //Do not update adjacent links
          for(mu=0; mu<4; mu++){
       	    l = link(s,mu);
      	    if(l % 8 == rem){
      	      // Compute staples
              u_zero(&staple);
               
              for (nu = 0; nu < 4; nu++)                  // Forward direction
       	        if (mu != nu)
                {
                  SU3* up[3];
                  SU3 t0, t1;
                  
                  up[0] = &ud[link(nnpd[4*s + mu], nu)];
                  up[1] = &ud[link(nnpd[4*s + nu], mu)];
                  up[2] = &ud[link(s, nu)];
                  
                  u_mul(&t0, up[2], up[1]);
                  u_dagger(&t0);
                  u_mul(&t1, up[0], &t0);
                  u_accum(&staple, &t1);
                }
              
     	      for (nu = 0; nu < 4; nu++)                  // Backward direction
                if (mu != nu)
                {
                  SU3* um[3];
                  SU3 t0, t1;
                  
                  um[0] = &ud[link(nnmd[4*nnpd[4*s + mu] + nu], nu)];
                  um[1] = &ud[link(nnmd[4*s + nu], mu)];
                  um[2] = &ud[link(nnmd[4*s + nu], nu)];
                  
                  u_mul(&t0, um[1], um[0]);
                  u_dagger(&t0);
                  u_mul(&t1, &t0, um[2]);
                  u_accum(&staple, &t1);
                }
                
       	      for (ihit = 0; ihit < METRO_NHIT; ihit++)
              {
                SU3 unew;
                
                u_metro_offer(&unew, &ud[l]);
                
                if (u_metro_accept(&staple, &ud[l], &unew))
                {
                  u_copy(&ud[l], &unew);
                  iacc++;
                }
              }
      	    }
            //queue.wait();	    
            accd = iacc / (double) (METRO_NHIT * NLINK);
      	  }
        } 
        accLoc.combine(accd);
	//acc = accLoc;
      });
    });
    
    auto accHostAcc = accBuffer.get_access<sycl::access::mode::read>();
    acc = accHostAcc[0];
    */
    //queue.wait();
    acc = u_sweep_metro();
    //acc = u_sweep_metro_gpu(); 
    end_metro = clock();
    printf("Time u_sweep_metro(): %f s\n", ((double) (end_metro - start_metro)) / CLOCKS_PER_SEC);
    
    //Plaquette action gets calculated from here on
    queue.submit([&](sycl::handler& cgh) {
      // Get an accessor for the buffer
      auto plaqAcc = plaqBuffer.get_access<sycl::access::mode::write>(cgh);
          
      // Execute the parallel_for algorithm on the GPU
      cgh.parallel_for<class PlaquetteKernel>(sycl::range<1>(LT * LS * LS * LS), sycl::reduction(plaqAcc, 0.0, std::plus<double>{}),[=](sycl::id<1> idx, auto& plaqLoc) {
        int t = idx / (LS * LS * LS);
        int z = (idx / (LS * LS)) % LS;
        int y = (idx / LS) % LS;
        int x = idx % LS;
          
        int s = site(x, y, z, t);
        double plaqd = 0.0;
         
        for (int mu = 0; mu < 4; mu++) {
          SU3* up[4];
          up[0] = &ud[link(s, mu)];
            
          for (int nu = mu + 1; nu < 4; nu++) {
            SU3 t0, t1;
            up[1] = &ud[link(nnpd[4 * s + mu], nu)];
            up[2] = &ud[link(nnpd[4 * s + nu], mu)];
            up[3] = &ud[link(s, nu)];
            
            u_mul(&t0, up[0], up[1]);
            u_mul(&t1, up[3], up[2]);
            u_dagger(&t1);
            
            double localPlaq = 0.0;
            for (int i = 0; i < NCOL; i++)
              for (int j = 0; j < NCOL; j++)
                localPlaq += real(t0.c[i][j] * t1.c[j][i]);
            plaqd += localPlaq;
          }
        }
           
        // Use the reduction extension to sum up the results across all threads
        plaqLoc.combine(plaqd);
      });
    });
    // Read the reduced result back to the host
    auto plaqHostAcc = plaqBuffer.get_access<sycl::access::mode::read>();
    plaq = plaqHostAcc[0];
    
    // Normalize by the number of lattice sites and the number of directions
    plaq /= 18. * VOL;

  //Parallelizing the Metropolis Algorithm
  
  printf("%6d     %.6e     %.2e\n", i_METRO, plaq, acc);
  fflush(stdout);

  }
  sycl::free(ud, queue);
  sycl::free(nnpd, queue);
  //sycl::free(nnmd, queue);
  //return plaq;
}


//--------------------------------------------------------------------------------------------------
//! Metropolis update: Update full lattice
//!
//! Computes stamples only once (and performs multiple hits per link):
//!
//!                p1
//!         +nu  <----
//!          |          ^
//!        p2|          | p0
//!          v          |
//!         site ----> +mu
//!          ^     1    |
//!        m2|          | m0        nu
//!          |          v           ^
//!         -nu  <---- -nu+mu       |
//!                m1                -> mu
//--------------------------------------------------------------------------------------------------

class MetroFunctor {
public:
  MetroFunctor(double *acc, SU3 *ud, int *nnpd, int *nnmd)
      : acc(acc), ud(ud), nnpd(nnpd), nnmd(nnmd) {}

  void operator()(id<1> idx) const {
    int l = idx[0];
    int nu;

    int s = l / 4;  // Assuming 4 links per site
    int mu = l % 4;

    // Only update every 8th link
    if (mu % 8 == 0) {
      SU3 staple;
      u_zero(&staple);

      for (nu = 0; nu < 4; nu++)                  // Forward direction
        if (mu != nu)
         {
           SU3* up[3];
           SU3 t0, t1;

           up[0] = &u[link(nnp[s][mu], nu)];
           up[1] = &u[link(nnp[s][nu], mu)];
           up[2] = &u[link(s, nu)];

           u_mul(&t0, up[2], up[1]);
           u_dagger(&t0);
           u_mul(&t1, up[0], &t0);
           u_accum(&staple, &t1);
          }

      for (nu = 0; nu < 4; nu++)                  // Backward direction
        if (mu != nu)
          {
            SU3* um[3];
            SU3 t0, t1;

            um[0] = &u[link(nnm[nnp[s][mu]][nu], nu)];
            um[1] = &u[link(nnm[s][nu], mu)];
            um[2] = &u[link(nnm[s][nu], nu)];

            u_mul(&t0, um[1], um[0]);
            u_dagger(&t0);
            u_mul(&t1, &t0, um[2]);
            u_accum(&staple, &t1);
          }


      for (int ihit = 0; ihit < METRO_NHIT; ihit++) {
        SU3 unew;
        u_metro_offer(&unew, &ud[l]);

        if (u_metro_accept(&staple, &ud[l], &unew)) {
          u_copy(&ud[l], &unew);
        }
      }
    }
  }

private:
  double *acc;
  SU3 *ud;
  int *nnpd;
  int *nnmd;
};

double u_sweep_metro_gpu(SU3 *d_ud, int *d_nnpd, int *d_nnmd) {
  double h_acc = 0.0;

  {
	  sycl::queue queue(sycl::gpu_selector{});

    // Allocate and copy data to the device
    buffer<double, 1> accBuffer(&h_acc, range<1>(1));
    buffer<SU3, 1> d_ud_buffer(d_ud, range<1>(NSITE * 4));
    buffer<int, 1> d_nnpd_buffer(d_nnpd, range<1>(NSITE * 4));
    buffer<int, 1> d_nnmd_buffer(d_nnmd, range<1>(NSITE * 4));

    // Launch the kernel
    queue.submit([&](handler &cgh) {
      auto accAcc = accBuffer.get_access<access::mode::write>(cgh);
      auto ud = d_ud_buffer.get_access<access::mode::read_write>(cgh);
      auto nnpd = d_nnpd_buffer.get_access<access::mode::read>(cgh);
      auto nnmd = d_nnmd_buffer.get_access<access::mode::read>(cgh);

      cgh.parallel_for<class MetroKernel>(range<1>(NSITE * 4), MetroFunctor(accAcc.get_pointer(), ud.get_pointer(), nnpd.get_pointer(), nnmd.get_pointer()));
    });

    // Copy results back to the host
    queue.wait_and_throw();
  }

  return h_acc / (METRO_NHIT * NLINK);
}





double u_sweep_metro(void)
{
  //int i, j;
  int ihit;
  int x, y, z, t;
  int mu, nu;
  int s, l;
  SU3 staple;
  int iacc;
  double acc;
  clock_t start, end;

  iacc = 0;

  for (t = 0; t < LT; t++)
    for (z = 0; z < LS; z++)
      for (y = 0; y < LS; y++)
        for (x = 0; x < LS; x++)
        {
          s = site(x, y, z, t);
          for (mu = 0; mu < 4; mu++)
          {
            l = link(s, mu);

            // Compute staples
            u_zero(&staple);

            for (nu = 0; nu < 4; nu++)			// Forward direction
              if (mu != nu)
              {
                SU3* up[3];
                SU3 t0, t1;

                up[0] = &u[link(nnp[s][mu], nu)];
                up[1] = &u[link(nnp[s][nu], mu)];
                up[2] = &u[link(s, nu)];

                u_mul(&t0, up[2], up[1]);
                u_dagger(&t0);
                u_mul(&t1, up[0], &t0);
                u_accum(&staple, &t1);
              }

            for (nu = 0; nu < 4; nu++)			// Backward direction
              if (mu != nu)
              {
                SU3* um[3];
                SU3 t0, t1;

                um[0] = &u[link(nnm[nnp[s][mu]][nu], nu)];
                um[1] = &u[link(nnm[s][nu], mu)];
                um[2] = &u[link(nnm[s][nu], nu)];

                u_mul(&t0, um[1], um[0]);
                u_dagger(&t0);
                u_mul(&t1, &t0, um[2]);
                u_accum(&staple, &t1);
              }

            start = clock();
            // Perform multiple hits
            for (ihit = 0; ihit < METRO_NHIT; ihit++)
            {
              SU3 unew;

              u_metro_offer(&unew, &u[l]);

              if (u_metro_accept(&staple, &u[l], &unew))
              {
                u_copy(&u[l], &unew);
                iacc++;
              }
            }
            end = clock();
            //printf("Time metro_hits(): %f s\n", ((double) (end - start)) / CLOCKS_PER_SEC);
          }
        }

  acc = iacc / (double) (METRO_NHIT * NLINK);

  return acc;
}

//--------------------------------------------------------------------------------------------------
//! Metropolis update: Offer new link variable
//--------------------------------------------------------------------------------------------------

void u_metro_offer(SU3* unew, SU3* uold)
{
  int k, l;
  SU3 uc;

  // Compute uc = bias * I + r where r is random matrix and I is unit matrix
  for (k = 0; k < NCOL-1; k++)
  {
    uc.c[k][k] = METRO_BIAS + std::complex<double>(0.0, rng() - 0.5);
    for (l = k+1; l < NCOL; l++)
    {
      uc.c[k][l] = (rng() - 0.5) + std::complex<double>(0.0, rng() - 0.5);
    }
  }
  uc.c[1][0] = -1.0 * conj(uc.c[0][1]);

  u_reunitarise(&uc);

  // Make sure that c and c+ are equally probable
  if (rng() < 0.5)
    u_dagger(&uc);

  // New link value obtained by multiplication with c
  u_mul(unew, &uc, uold);
}

//--------------------------------------------------------------------------------------------------
//! Metropolis update: Accept or reject new link variable
//--------------------------------------------------------------------------------------------------

int u_metro_accept(SU3* staple, SU3* uold, SU3* unew)
{
  int k, l;
  double aold;
  double anew;
  double adiff;
  int accept;

  // Compute change of action when using new link value
  aold = 0.0;
  anew = 0.0;
  for (k = 0; k < NCOL; k++)
    for (l = 0; l < NCOL; l++)
    {
      aold += real(uold->c[k][l] * staple->c[l][k]);
      anew += real(unew->c[k][l] * staple->c[l][k]);
    }
  adiff = -0.3333333 * (aold - anew);

  // Metropolis step
  if (log(rng()) > BETA * adiff)
    accept = 0;
  else
    accept = 1;

  return accept;
}
