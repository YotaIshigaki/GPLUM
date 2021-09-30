/* Standard headers */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>
#include <algorithm>
/* FDPS headers */
#include <particle_simulator.hpp>

/* Macros */
#define GET_POS_DOMAIN_ID_CHECK  (0)
//#define GET_POS_DOMAIN_DEBUG_PRINT
//#define GET_POS_DOMAIN_RESULT_PRINT
//#define GET_POS_DOMAIN_TIMING_PRINT

#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#define MAX(x,y) (((x) > (y)) ? (x) : (y))

/* Local constants */
const PS::F64 LARGE_FLOAT = std::numeric_limits<PS::F32>::max()*0.0625;

void GetPosDomain(PS::F64 ain, PS::F64 aout, PS::F64ort &pos_domain,
                  PS::F64 amax=LARGE_FLOAT) {
    PS::F64 starttime = MPI::Wtime();
    assert(amax > aout);

    //* Get rank number and the number of processes
    const PS::S32 myrank = PS::Comm::getRank();
    const PS::S32 nprocs = PS::Comm::getNumberOfProc();
    if (nprocs == 1) {
        pos_domain.low_.x = - amax;
        pos_domain.low_.y = - amax;
        pos_domain.low_.z = - amax;
        pos_domain.high_.x = amax;
        pos_domain.high_.y = amax;
        pos_domain.high_.z = amax;
        return;
    }
    assert(nprocs%4 == 0);

    //* Compute the ratio of the areas of the ring to the square region
    const PS::F64 PI = 4.0 * std::atan(1.0);
    const PS::F64 area_ratio = (0.25 * PI * (ain + aout) * (aout - ain))/(aout * aout);
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
        std::cout << "ain  = " << ain << std::endl;
        std::cout << "aout = " << aout << std::endl;
        std::cout << "area_ratio = " << area_ratio << std::endl;
    }
#endif

    //* Compute the number of sample points
    const PS::S32 num_points_per_side = (PS::S32)(4096/std::sqrt(area_ratio));
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK)
        std::cout << "num_points_per_side = " << num_points_per_side << std::endl;
#endif

    //* Compute the number of domains per side
    PS::S32 n_domain[3];
    n_domain[0] = std::sqrt((PS::F64)nprocs-0.000001)+1;
    while( nprocs % n_domain[0] != 0) n_domain[0]++;
    n_domain[1] = nprocs / n_domain[0];
    n_domain[2] = 1;
    assert(n_domain[0]*n_domain[1] == nprocs);
    assert(n_domain[0]%2 == 0);
    assert(n_domain[1]%2 == 0);
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
        std::cout << "n_domain[0] = " << n_domain[0] << std::endl;
        std::cout << "n_domain[1] = " << n_domain[1] << std::endl;
        std::cout << "n_domain[2] = " << n_domain[2] << std::endl;
    }
#endif
 
    //* Compute quadrant_num, ix_*, iy_*
    PS::S32 quadrant_num;
    PS::S32 ix_glb,iy_glb,ix_ref,iy_ref,ix_1q,iy_1q;
    ix_glb = myrank / n_domain[1];
    iy_glb = myrank % n_domain[1];
    if (myrank < nprocs/2) {
       if (iy_glb < n_domain[1]/2) {
           // x < 0, y < 0 (3rd quadrant) 
           quadrant_num = 3; 
           ix_ref = n_domain[0]/2 -1;
           iy_ref = n_domain[1]/2 -1;
       } else {
           // x < 0, y > 0 (2nd quadrant) 
           quadrant_num = 2;
           ix_ref = n_domain[0]/2 -1;
           iy_ref = n_domain[1]/2;
       }
    } else {
       if (iy_glb < n_domain[1]/2) {
           // x > 0, y < 0 (4th quadrant) 
           quadrant_num = 4;
           ix_ref = n_domain[0]/2;
           iy_ref = n_domain[1]/2-1;
       } else {
           // x > 0, y > 0 (1st quadrant) 
           quadrant_num = 1;
           ix_ref = n_domain[0]/2;
           iy_ref = n_domain[1]/2;
       }
    }
    ix_1q = std::abs(ix_glb - ix_ref);
    iy_1q = std::abs(iy_glb - iy_ref);
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
#if 0
    if (quadrant_num == 1) {
        std::cout << "(1st) " << myrank
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    } else if (quadrant_num == 2) {
        std::cout << "(2nd) " << myrank 
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    } else if (quadrant_num == 3) {
        std::cout << "(3rd) " << myrank 
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    } else {
        std::cout << "(4th) " << myrank 
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    }
#endif
#endif


    //* Compute pos_domain
    //** Initialize
    pos_domain.low_.z  = - amax; 
    pos_domain.high_.z =   amax; 
    //** Compute the interval of sample points
    PS::F64 dx = aout/(num_points_per_side-1);
    //** Compute the number of sample points in the ring
    PS::F64 x,y,r2;
    const PS::F64 ain2 = ain * ain, aout2 = aout * aout;
    PS::S32 num_points_tot = 0;
    for (PS::S32 i=0; i<num_points_per_side; i++) {
        x = dx * i;
        for (PS::S32 j=0; j<num_points_per_side; j++) {
            y = dx * j;
            r2 = x*x + y*y; 
            if ((ain2 <= r2) && (r2 < aout2)) 
                num_points_tot++;
        }
    }
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) 
        std::cout << "num_points_tot = " << num_points_tot << std::endl;
#endif
    //** Domain decomposition w.r.t. x-coordinate
    PS::S32 num_points_in_x_slab = num_points_tot/(n_domain[0]/2) + 1;
    PS::S32 num_points_in_my_x_slab = 0;
    PS::F64 * xdiv = new PS::F64[n_domain[0]/2 + 1];
    xdiv[0] = 0.0;
    xdiv[n_domain[0]/2] = amax;
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) 
        std::cout << "num_points_in_x_slab = " << num_points_in_x_slab << std::endl;
#endif
    PS::S32 num_points=0,k=1,next=num_points_in_x_slab;
    for (PS::S32 i=0; i<num_points_per_side; i++) {
        x = dx * i;
        for (PS::S32 j=0; j<num_points_per_side; j++) {
            y = dx * j;
            r2 = x*x + y*y; 
            if ((ain2 <= r2) && (r2 < aout2)) {
                if (k-1 == ix_1q) {
                    num_points_in_my_x_slab++;
                    // This is a rough number of points in my x-slab.
                }
                num_points++;
                if (num_points == next) { 
                    if (k == n_domain[0]/2) goto OUT1;
                    xdiv[k] = x;
                    k++;
                    next += num_points_in_x_slab;
                }
            }
        }
    }
    OUT1:
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
       std::cout << "num_points_in_my_x_slab = " << num_points_in_my_x_slab << std::endl;
       for (PS::S32 i=0; i<n_domain[0]/2 + 1; i++)
           std::cout << "xdiv[" << i << "] = " << xdiv[i] << std::endl;
    }
#endif

    //** Domain decomposition w.r.t. y-coordinate
    PS::S32 num_points_in_y_slab = num_points_in_my_x_slab/(n_domain[1]/2) + 1;
    PS::F64 * ydiv = new PS::F64[n_domain[1]/2 + 1];
    ydiv[0] = 0.0;
    ydiv[n_domain[1]/2] = amax;
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
   if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
       std::cout << "num_points_in_y_slab = " << num_points_in_y_slab << std::endl;
       std::cout << "xdiv[ix_1q]   = " << xdiv[ix_1q] << std::endl;
       std::cout << "xdiv[ix_1q+1] = " << xdiv[ix_1q+1] << std::endl;
   }
#endif
    num_points=0,k=1,next=num_points_in_y_slab;
    for (PS::S32 j=0; j<num_points_per_side; j++) {
        y = dx * j;
        for (PS::S32 i=0; i<num_points_per_side; i++) {
            x = dx * i;
            if ((xdiv[ix_1q] <= x) && (x < xdiv[ix_1q+1])) {
                r2 = x*x + y*y; 
                if ((ain2 <= r2) && (r2 < aout2)) {
                    num_points++;
                    if (num_points == next) { 
                        if (k == n_domain[1]/2) goto OUT2;
                        ydiv[k] = y;
                        k++;
                        next += num_points_in_y_slab;
                    }
                }
            }
        }
    } 
    OUT2:
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
       std::cout << "num_points = " << num_points << std::endl;
       for (PS::S32 j=0; j<n_domain[1]/2 + 1; j++)
           std::cout << "ydiv[" << j << "] = " << ydiv[j] << std::endl;
    }
#endif

    //* Set pos_domain
    if (quadrant_num == 1) {
       pos_domain.low_.x  = xdiv[ix_1q];
       pos_domain.low_.y  = ydiv[iy_1q];
       pos_domain.high_.x = xdiv[ix_1q+1];
       pos_domain.high_.y = ydiv[iy_1q+1];
    } else if (quadrant_num == 2) {
       pos_domain.low_.x  = - xdiv[ix_1q+1];
       pos_domain.low_.y  = ydiv[iy_1q];
       pos_domain.high_.x = - xdiv[ix_1q]; 
       pos_domain.high_.y = ydiv[iy_1q+1];
    } else if (quadrant_num == 3) {
       pos_domain.low_.x  = - xdiv[ix_1q+1];
       pos_domain.low_.y  = - ydiv[iy_1q+1];
       pos_domain.high_.x = - xdiv[ix_1q];
       pos_domain.high_.y = - ydiv[iy_1q];
    } else if (quadrant_num == 4) {
       pos_domain.low_.x  = xdiv[ix_1q]; 
       pos_domain.low_.y  = - ydiv[iy_1q+1]; 
       pos_domain.high_.x = xdiv[ix_1q+1]; 
       pos_domain.high_.y = - ydiv[iy_1q]; 
    }
#ifdef GET_POS_DOMAIN_RESULT_PRINT
    std::cout << myrank << "  "
              << pos_domain.low_.x << "  "
              << pos_domain.low_.y << "  "
              << std::endl;
    std::cout << myrank << "  " 
              << pos_domain.high_.x << "  "
              << pos_domain.high_.y << "  "
              << std::endl;
#endif

    //* Release memory
    delete [] xdiv;
    delete [] ydiv;

#ifdef GET_POS_DOMAIN_TIMING_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK)
        std::cout << "etime(GetPosDomain) = " << MPI::Wtime() - starttime << std::endl;
#endif
}

//-------------------------------------------------------------------
void GetPosDomain2(PS::F64 ain, PS::F64 aout, PS::F64ort &pos_domain,
                   PS::F64 amax=LARGE_FLOAT) {
    PS::F64 starttime = MPI::Wtime(), starttime_tmp;
    assert(amax > aout);

    //* Get rank number and the number of processes
    const PS::S32 myrank = PS::Comm::getRank();
    const PS::S32 nprocs = PS::Comm::getNumberOfProc();
    if (nprocs == 1) {
        pos_domain.low_.x = - amax;
        pos_domain.low_.y = - amax;
        pos_domain.low_.z = - amax;
        pos_domain.high_.x = amax;
        pos_domain.high_.y = amax;
        pos_domain.high_.z = amax;
        return;
    }
    assert(nprocs%4 == 0);

    //* Compute the ratio of the areas of the ring to the square region
    const PS::F64 PI = 4.0 * std::atan(1.0);
    const PS::F64 area_ratio = (0.25 * PI * (ain + aout) * (aout - ain))/(aout * aout);
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
        std::cout << "ain  = " << ain << std::endl;
        std::cout << "aout = " << aout << std::endl;
        std::cout << "area_ratio = " << area_ratio << std::endl;
    }
#endif

    //* Compute the number of sample points
    //const PS::S32 num_points_per_side = (PS::S32)(4096/std::sqrt(area_ratio));
    const PS::S32 num_points_per_side = (PS::S32)(16384/std::sqrt(area_ratio));
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK)
        std::cout << "num_points_per_side = " << num_points_per_side << std::endl;
#endif

    //* Compute the number of domains per side
    PS::S32 n_domain[3];
    n_domain[0] = std::sqrt((PS::F64)nprocs-0.000001)+1;
    while( nprocs % n_domain[0] != 0) n_domain[0]++;
    n_domain[1] = nprocs / n_domain[0];
    n_domain[2] = 1;
    assert(n_domain[0]*n_domain[1] == nprocs);
    assert(n_domain[0]%2 == 0);
    assert(n_domain[1]%2 == 0);
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
        std::cout << "n_domain[0] = " << n_domain[0] << std::endl;
        std::cout << "n_domain[1] = " << n_domain[1] << std::endl;
        std::cout << "n_domain[2] = " << n_domain[2] << std::endl;
    }
#endif
 
    //* Compute quadrant_num, ix_*, iy_*
    PS::S32 quadrant_num;
    PS::S32 ix_glb,iy_glb,ix_ref,iy_ref,ix_1q,iy_1q;
    ix_glb = myrank / n_domain[1];
    iy_glb = myrank % n_domain[1];
    if (myrank < nprocs/2) {
       if (iy_glb < n_domain[1]/2) {
           // x < 0, y < 0 (3rd quadrant) 
           quadrant_num = 3; 
           ix_ref = n_domain[0]/2 -1;
           iy_ref = n_domain[1]/2 -1;
       } else {
           // x < 0, y > 0 (2nd quadrant) 
           quadrant_num = 2;
           ix_ref = n_domain[0]/2 -1;
           iy_ref = n_domain[1]/2;
       }
    } else {
       if (iy_glb < n_domain[1]/2) {
           // x > 0, y < 0 (4th quadrant) 
           quadrant_num = 4;
           ix_ref = n_domain[0]/2;
           iy_ref = n_domain[1]/2-1;
       } else {
           // x > 0, y > 0 (1st quadrant) 
           quadrant_num = 1;
           ix_ref = n_domain[0]/2;
           iy_ref = n_domain[1]/2;
       }
    }
    ix_1q = std::abs(ix_glb - ix_ref);
    iy_1q = std::abs(iy_glb - iy_ref);
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
#if 0
    if (quadrant_num == 1) {
        std::cout << "(1st) " << myrank
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    } else if (quadrant_num == 2) {
        std::cout << "(2nd) " << myrank 
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    } else if (quadrant_num == 3) {
        std::cout << "(3rd) " << myrank 
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    } else {
        std::cout << "(4th) " << myrank 
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    }
#endif
#endif


    //* Compute pos_domain
    //** Initialize
    pos_domain.low_.z  = - amax; 
    pos_domain.high_.z =   amax; 
    //** Compute the interval of sample points
    PS::F64 dx = aout/(num_points_per_side-1);
    //** Compute the number of sample points in the ring
#ifdef GET_POS_DOMAIN_TIMING_PRINT
    starttime_tmp = MPI::Wtime();
#endif
    PS::F64 x,x2,y,y2,r2;
    PS::S32 jsta,jend;
    const PS::F64 fac = 2.0;
    const PS::F64 ain2 = ain * ain, aout2 = aout * aout;
    PS::S32 num_points_tot = 0;
    for (PS::S32 i=0; i<num_points_per_side; i++) {
        x = dx * i;
        x2 = x * x;
        PS::F64 tmp1 = (ain - fac*dx)*(ain - fac*dx)   - x2;
        if (tmp1 > 0.0) {
           jsta = (PS::S32) (std::sqrt(tmp1)/dx);
        } else {
           jsta = 0;
        }
        PS::F64 tmp2 = (aout + fac*dx)*(aout + fac*dx) - x2;
        if (tmp2 > 0.0) {
           jend = (PS::S32) (std::sqrt(tmp2)/dx);
        } else {
           jend = 0;
        }
        for (PS::S32 j=jsta; j<=jend; j++) {
            y = dx * j;
            r2 = x*x + y*y; 
            if ((ain2 <= r2) && (r2 < aout2)) 
                num_points_tot++;
        }
    }
    // For test run
    //athread_halt();
    //PS::Finalize();
    //std::exit(0);
#ifdef GET_POS_DOMAIN_TIMING_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK)
        std::cout << "etime(GetPosDomain2,cnt) = " << MPI::Wtime() - starttime_tmp << std::endl;
#endif
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) 
        std::cout << "num_points_tot = " << num_points_tot << std::endl;
#endif
    //** Domain decomposition w.r.t. x-coordinate
    PS::S32 num_points_in_x_slab = num_points_tot/(n_domain[0]/2) + 1;
    PS::S32 num_points_in_my_x_slab = 0;
    PS::F64 * xdiv = new PS::F64[n_domain[0]/2 + 1];
    xdiv[0] = 0.0;
    xdiv[n_domain[0]/2] = amax;
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) 
        std::cout << "num_points_in_x_slab = " << num_points_in_x_slab << std::endl;
#endif
#ifdef GET_POS_DOMAIN_TIMING_PRINT
    starttime_tmp = MPI::Wtime();
#endif
    PS::S32 num_points=0,k=1,next=num_points_in_x_slab;
    for (PS::S32 i=0; i<num_points_per_side; i++) {
        x = dx * i;
        x2 = x * x;
        PS::F64 tmp1 = (ain - fac*dx)*(ain - fac*dx)   - x2;
        if (tmp1 > 0.0) {
           jsta = (PS::S32) (std::sqrt(tmp1)/dx);
        } else {
           jsta = 0;
        }
        PS::F64 tmp2 = (aout + fac*dx)*(aout + fac*dx) - x2;
        if (tmp2 > 0.0) {
           jend = (PS::S32) (std::sqrt(tmp2)/dx);
        } else {
           jend = 0;
        }
        for (PS::S32 j=jsta; j<=jend; j++) {
            y = dx * j;
            r2 = x*x + y*y; 
            if ((ain2 <= r2) && (r2 < aout2)) {
                if (k-1 == ix_1q) {
                    num_points_in_my_x_slab++;
                    // This is a rough number of points in my x-slab.
                }
                num_points++;
                if (num_points == next) { 
                    if (k == n_domain[0]/2) goto OUT1;
                    xdiv[k] = x;
                    k++;
                    next += num_points_in_x_slab;
                }
            }
        }
    }
    OUT1:
#ifdef GET_POS_DOMAIN_TIMING_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK)
        std::cout << "etime(GetPosDomain2,xdiv) = " << MPI::Wtime() - starttime_tmp << std::endl;
#endif
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
       std::cout << "num_points_in_my_x_slab = " << num_points_in_my_x_slab << std::endl;
       for (PS::S32 i=0; i<n_domain[0]/2 + 1; i++)
           std::cout << "xdiv[" << i << "] = " << xdiv[i] << std::endl;
    }
#endif

    //** Domain decomposition w.r.t. y-coordinate
    PS::S32 num_points_in_y_slab = num_points_in_my_x_slab/(n_domain[1]/2) + 1;
    PS::F64 * ydiv = new PS::F64[n_domain[1]/2 + 1];
    ydiv[0] = 0.0;
    ydiv[n_domain[1]/2] = amax;
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
        std::cout << "num_points_in_y_slab = " << num_points_in_y_slab << std::endl;
        std::cout << "xdiv[ix_1q]   = " << xdiv[ix_1q] << std::endl;
        std::cout << "xdiv[ix_1q+1] = " << xdiv[ix_1q+1] << std::endl;
    }
#endif
#ifdef GET_POS_DOMAIN_TIMING_PRINT
    starttime_tmp = MPI::Wtime();
#endif
    PS::S32 ista,iend;
    PS::F64 xmin = MIN(xdiv[ix_1q],xdiv[ix_1q+1]);
    PS::F64 xmax = MAX(xdiv[ix_1q],xdiv[ix_1q+1]);
    xmax = MIN(xmax,1.05*aout);
    //*** Compute jsta
    x2 = xmax * xmax;
    PS::F64 tmp1 = (ain - fac*dx)*(ain - fac*dx) - x2; 
    if (tmp1 > 0.0) {
        jsta = (PS::S32) (std::sqrt(tmp1)/dx);
    } else {
       jsta = 0;
    }
    //*** Compute jsta
    x2 = xmin * xmin;
    PS::F64 tmp2 = (aout + fac*dx)*(aout + fac*dx) - x2;
    if (tmp2 > 0.0) {
       jend = (PS::S32) (std::sqrt(tmp2)/dx);
    } else {
       jend = 0;
    }
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
    std::cout << "(myrank,jsta,jend) = " 
              << myrank << "," << jsta << "," << jend << std::endl;
#endif
    num_points=0,k=1,next=num_points_in_y_slab;
    for (PS::S32 j=jsta; j<=jend; j++) {
        y = dx * j;
        y2 = y * y;
        PS::F64 tmp1 = (ain - fac*dx)*(ain - fac*dx)   - y2;
        if (tmp1 > 0.0) {
           ista = (PS::S32) (std::sqrt(tmp1)/dx);
        } else {
           ista = 0;
        }
        PS::F64 tmp2 = (aout + fac*dx)*(aout + fac*dx) - y2;
        if (tmp2 > 0.0) {
           iend = (PS::S32) (std::sqrt(tmp2)/dx);
        } else {
           iend = 0;
        }
        for (PS::S32 i=ista; i<=iend; i++) { 
            x = dx * i;
            if ((xdiv[ix_1q] <= x) && (x < xdiv[ix_1q+1])) {
                r2 = x*x + y*y; 
                if ((ain2 <= r2) && (r2 < aout2)) {
                    num_points++;
                    if (num_points == next) { 
                        if (k == n_domain[1]/2) goto OUT2;
                        ydiv[k] = y;
                        k++;
                        next += num_points_in_y_slab;
                    }
                }
            }
        }
    } 
    OUT2:
    // For test run
    //athread_halt();
    //PS::Finalize();
    //std::exit(0);
#ifdef GET_POS_DOMAIN_TIMING_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK)
        std::cout << "etime(GetPosDomain2,ydiv) = " << MPI::Wtime() - starttime_tmp << std::endl;
#endif
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    PS::Comm::barrier();
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
       std::cout << "num_points = " << num_points << std::endl;
       for (PS::S32 j=0; j<n_domain[1]/2 + 1; j++)
           std::cout << "ydiv[" << j << "] = " << ydiv[j] << std::endl;
    }
#endif

    //* Set pos_domain
    if (quadrant_num == 1) {
       pos_domain.low_.x  = xdiv[ix_1q];
       pos_domain.low_.y  = ydiv[iy_1q];
       pos_domain.high_.x = xdiv[ix_1q+1];
       pos_domain.high_.y = ydiv[iy_1q+1];
    } else if (quadrant_num == 2) {
       pos_domain.low_.x  = - xdiv[ix_1q+1];
       pos_domain.low_.y  = ydiv[iy_1q];
       pos_domain.high_.x = - xdiv[ix_1q]; 
       pos_domain.high_.y = ydiv[iy_1q+1];
    } else if (quadrant_num == 3) {
       pos_domain.low_.x  = - xdiv[ix_1q+1];
       pos_domain.low_.y  = - ydiv[iy_1q+1];
       pos_domain.high_.x = - xdiv[ix_1q];
       pos_domain.high_.y = - ydiv[iy_1q];
    } else if (quadrant_num == 4) {
       pos_domain.low_.x  = xdiv[ix_1q]; 
       pos_domain.low_.y  = - ydiv[iy_1q+1]; 
       pos_domain.high_.x = xdiv[ix_1q+1]; 
       pos_domain.high_.y = - ydiv[iy_1q]; 
    }
#ifdef GET_POS_DOMAIN_RESULT_PRINT
    std::cout << myrank << "  "
              << pos_domain.low_.x << "  "
              << pos_domain.low_.y << "  "
              << std::endl;
    std::cout << myrank << "  " 
              << pos_domain.high_.x << "  "
              << pos_domain.high_.y << "  "
              << std::endl;
#endif

    //* Release memory
    delete [] xdiv;
    delete [] ydiv;

#ifdef GET_POS_DOMAIN_TIMING_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK)
        std::cout << "etime(GetPosDomain2) = " << MPI::Wtime() - starttime << std::endl;
#endif
}

//-------------------------------------------------------------------

double cut_fan_area(double r,double x)
{
    if ( x >= r){
	return 0.0;
    }else{
	double theta=acos(x/r);
	return    (r*r*theta -r*x*sin(theta))/2;
    }
}

double band_area(double x0,double x1, double r0, double r1)
{
    return    cut_fan_area(r1,x0)-cut_fan_area(r1,x1) -
	cut_fan_area(r0,x0)+cut_fan_area(r0,x1) ;
}

double  patch_dx(double y,double x0,double x1, double r0, double r1)
{
    
    double yx0r1 = sqrt(r1*r1-x0*x0);
    double  yx0r0=0;
    double  yx1r1=sqrt(r1*r1-x1*x1);
    double  yx1r0=0;
    if (r0>x0){
	yx0r0 = sqrt(r0*r0-x0*x0);
    }
    if (r0>x1){
	yx1r0 = sqrt(r0*r0-x1*x1) ;
    }
    double  xmin=x0;
    double      xmax=x1;
    if (y < yx0r0){
	xmin = sqrt(r0*r0-y*y);
    }
    if (y> yx1r1){
	xmax = sqrt(r1*r1-y*y);
    }
    double dx = xmax-xmin;
    if (dx<0.0)  dx=0.0;
    return dx;
}

double  patch_area(double y0,
		   double y1,
		   double x0,
		   double x1,
		   double r0,
		   double r1)
{
    int  n=20;
    double  s= (patch_dx(y0,x0,x1,r0,r1)+patch_dx(y1,x0,x1,r0,r1))/2;
    int i;
    for(i=0;i<n-1;i++){
	s+=patch_dx(y0+(i+1)*(y1-y0)/n,x0,x1,r0,r1);
    }
    return   s/n*(y1-y0);
}

double  calc_x0(double r0, double r1, int n, double x1)
{
    double s= M_PI*(r1*r1-r0*r0)/n/4;
    double  a=0;
    double  b=x1;
    double  x=b;
    int i;
    for (i=0;i<50;i++){
	x= (a+b)/2;
	if (band_area(x,x1,r0,r1)-s > 0.0){
	    a=x;
	}else{
	    b=x;
	}
    }
    return x;
}

double  calc_y0(double x0, double x1, double r0, double r1, int n, double y1)
{
    //    fprintf(stderr,"calc y0 called with %e %e %e %e  %d %e\n",
    //	    x0, x1, r0, r1, n, y1);
    double s= M_PI*(r1*r1-r0*r0)/n/4;
    double  yx1r0=0;
    if (r0>x1){
	yx1r0 = sqrt(r0*r0-x1*x1) ;
    }
    double a=yx1r0;
    double  b=y1;
    double   y=b;
    int i;
    for (i=0;i<50;i++){
	y= (a+b)/2;
	if (patch_area(y,y1,x0,x1,r0,r1)-s > 0.0){
	    a=y;
	}else{
	    b=y;
	}
    }
    return y;
}

void box(int nx,
	 int ny,
	 int ix,
	 int iy,
	 double dr,
	 double *px0,
	 double *px1,
	 double *py0,
	 double *py1 )
{
    int i,j;
    double r0=1-dr/2;
    double r1=1+dr/2;
    double x1=r1;
    double x1old;
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    fprintf(stderr, "test routine with ix, iy=%d %d\n", ix, iy);
#endif
    double x0;
    for (i=nx;i>ix;i--){
	 x0=calc_x0(r0,r1,nx,x1);
	 x1old=x1;
	 x1=x0;
    }
    
    double y1 = sqrt(r1*r1-x0*x0);
    double y1old;
    double y0;
    for(j=ny;j>iy;j--){
	y0=calc_y0(x0,x1old,r0,r1,nx*ny,y1);
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
	fprintf(stderr,"y0, y1= %e %e\n", y0,y1);
#endif
	y1old=y1;
	y1=y0;
    }
    *px0=x0;
    *px1=x1old;
    *py0=y0;
    *py1=y1old;
}

    
    
void GetPosDomain3(PS::F64 delta_ax, PS::F64ort &pos_domain,
                   PS::S32 nx=-1, PS::S32 ny=-1, 
                   PS::F64 amax=LARGE_FLOAT) {
    PS::F64 starttime = MPI::Wtime();

    //* Get rank number and the number of processes
    const PS::S32 myrank = PS::Comm::getRank();
    const PS::S32 nprocs = PS::Comm::getNumberOfProc();
    if (nprocs == 1) {
        pos_domain.low_.x = - amax;
        pos_domain.low_.y = - amax;
        pos_domain.low_.z = - amax;
        pos_domain.high_.x = amax;
        pos_domain.high_.y = amax;
        pos_domain.high_.z = amax;
        return;
    }
    assert(nprocs%4 == 0);

    //* Compute the number of domains per side
    PS::S32 n_domain[3];
    if (nx == -1 || ny == -1) {
        n_domain[0] = std::sqrt((PS::F64)nprocs-0.000001)+1;
        while( nprocs % n_domain[0] != 0) n_domain[0]++;
        n_domain[1] = nprocs / n_domain[0];
    } else {
        n_domain[0] = nx;
        n_domain[1] = ny;
    }
    n_domain[2] = 1;
    assert(n_domain[0]*n_domain[1] == nprocs);
    assert(n_domain[0]%2 == 0);
    assert(n_domain[1]%2 == 0);
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK) {
        std::cout << "n_domain[0] = " << n_domain[0] << std::endl;
        std::cout << "n_domain[1] = " << n_domain[1] << std::endl;
        std::cout << "n_domain[2] = " << n_domain[2] << std::endl;
    }
#endif

    //* Compute quadrant_num, ix_*, iy_*
    PS::S32 quadrant_num;
    PS::S32 ix_glb,iy_glb,ix_ref,iy_ref,ix_1q,iy_1q;
    ix_glb = myrank / n_domain[1];
    iy_glb = myrank % n_domain[1];
    if (myrank < nprocs/2) {
       if (iy_glb < n_domain[1]/2) {
           // x < 0, y < 0 (3rd quadrant) 
           quadrant_num = 3; 
           ix_ref = n_domain[0]/2 -1;
           iy_ref = n_domain[1]/2 -1;
       } else {
           // x < 0, y > 0 (2nd quadrant) 
           quadrant_num = 2;
           ix_ref = n_domain[0]/2 -1;
           iy_ref = n_domain[1]/2;
       }
    } else {
       if (iy_glb < n_domain[1]/2) {
           // x > 0, y < 0 (4th quadrant) 
           quadrant_num = 4;
           ix_ref = n_domain[0]/2;
           iy_ref = n_domain[1]/2-1;
       } else {
           // x > 0, y > 0 (1st quadrant) 
           quadrant_num = 1;
           ix_ref = n_domain[0]/2;
           iy_ref = n_domain[1]/2;
       }
    }
    ix_1q = std::abs(ix_glb - ix_ref);
    iy_1q = std::abs(iy_glb - iy_ref);
#ifdef GET_POS_DOMAIN_DEBUG_PRINT
#if 0
    if (quadrant_num == 1) {
        std::cout << "(1st) " << myrank
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    } else if (quadrant_num == 2) {
        std::cout << "(2nd) " << myrank 
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    } else if (quadrant_num == 3) {
        std::cout << "(3rd) " << myrank 
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    } else {
        std::cout << "(4th) " << myrank 
                  << " (" << ix_1q << "," << iy_1q << ")" << std::endl;
    }
#endif
#endif

    //* Compute the positions of the domain by numerical integration
    pos_domain.low_.z  = - amax;
    pos_domain.high_.z =   amax;
    double x0,y0,x1,y1;
    box(n_domain[0]/2, n_domain[1]/2,
        ix_1q, iy_1q, delta_ax,
        &x0, &x1, &y0, &y1);

    //* Boundary processing
    if (ix_1q == 0) x0 = 0.0;
    if (ix_1q == n_domain[0]/2-1) x1 = amax;
    if (iy_1q == 0) y0 = 0.0;
    if (iy_1q == n_domain[1]/2-1) y1 = amax;

    //* Set pos_domain
    if (quadrant_num == 1) {
       pos_domain.low_.x  = x0;
       pos_domain.low_.y  = y0;
       pos_domain.high_.x = x1;
       pos_domain.high_.y = y1;
    } else if (quadrant_num == 2) {
       pos_domain.low_.x  = - x1;
       pos_domain.low_.y  = y0;
       pos_domain.high_.x = - x0;
       pos_domain.high_.y = y1;
    } else if (quadrant_num == 3) {
       pos_domain.low_.x  = - x1;
       pos_domain.low_.y  = - y1;
       pos_domain.high_.x = - x0;
       pos_domain.high_.y = - y0;
    } else if (quadrant_num == 4) {
       pos_domain.low_.x  = x0;
       pos_domain.low_.y  = - y1;
       pos_domain.high_.x = x1;
       pos_domain.high_.y = - y0;
    }
#ifdef GET_POS_DOMAIN_RESULT_PRINT
    std::cout << myrank << " "
              << pos_domain.low_.x << "  "
              << pos_domain.low_.y << "  "
              << pos_domain.high_.x << "  "
              << pos_domain.high_.y << std::endl;
#endif

#ifdef GET_POS_DOMAIN_TIMING_PRINT
    if (PS::Comm::getRank() == GET_POS_DOMAIN_ID_CHECK)
        std::cout << "etime(GetPosDomain3) = " << MPI::Wtime() - starttime << std::endl;
#endif
}

void GetPosDomainCyl(const PS::F64 delta_ax,
                     PS::F64ort & pos_domain,
                     const PS::S32 nx,
                     const PS::S32 ny){
    PS::F64 PI = 4.0 * atan(1.0);
    PS::F64 len_x = 2.0 * PI;
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 rank_x = my_rank / ny;
    PS::S32 rank_y = my_rank % ny;
    PS::F64 dx = len_x / nx;
    PS::F64 dy = delta_ax / ny;
    PS::F64 dy_offset = 1.0-delta_ax*0.5;
    pos_domain.low_.x = dx*rank_x;
    pos_domain.low_.y = dy*rank_y + dy_offset;
    pos_domain.low_.z = -PI;
    pos_domain.high_.x = dx*(rank_x+1);
    pos_domain.high_.y = dy*(rank_y+1) + dy_offset;
    pos_domain.high_.z = PI;
}
