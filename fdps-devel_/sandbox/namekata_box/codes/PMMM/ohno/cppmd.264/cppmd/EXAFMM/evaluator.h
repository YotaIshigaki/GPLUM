/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef evaluator_h
#define evaluator_h
#include "kernel.h"
#define splitFirst(Ci,Cj) Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R)

//! Interface between tree and kernel
class Evaluator : public Kernel {
protected:
  C_iter      CiB;                                              //!< icells begin per call
  C_iter      CiE;                                              //!< icells end per call
  C_iter      CjB;                                              //!< jcells begin per call
  C_iter      CjE;                                              //!< jcells end per call
  Lists       listM2L;                                          //!< M2L interaction list
  Lists       listM2P;                                          //!< M2P interaction list
  Lists       listP2P;                                          //!< P2P interaction list
  real        timeM2L;                                          //!< M2L execution time
  real        timeM2P;                                          //!< M2P execution time
  real        timeP2P;                                          //!< P2P execution time

  int         Iperiodic;                                        //!< Periodic image flag (using each bit for images)
  const int   Icenter;                                          //!< Periodic image flag at center
  Maps        flagM2L;                                          //!< Existance of periodic image for M2L
  Maps        flagM2P;                                          //!< Existance of periodic image for M2P
  Maps        flagP2P;                                          //!< Existance of periodic image for P2P

private:
//! Approximate interaction between two cells
  inline void approximate(C_iter Ci, C_iter Cj) {
#if HYBRID
    if( timeP2P*Cj->NDLEAF < timeM2P && timeP2P*Ci->NDLEAF*Cj->NDLEAF < timeM2L) {// If P2P is fastest
      evalP2P(Ci,Cj);                                           //  Evaluate on CPU, queue on GPU
    } else if ( timeM2P < timeP2P*Cj->NDLEAF && timeM2P*Ci->NDLEAF < timeM2L ) {// If M2P is fastest
      evalM2P(Ci,Cj);                                           //  Evaluate on CPU, queue on GPU
    } else {                                                    // If M2L is fastest
      evalM2L(Ci,Cj);                                           //  Evaluate on CPU, queue on GPU
    }                                                           // End if for kernel selection
#elif TREECODE
    evalM2P(Ci,Cj);                                             // Evaluate on CPU, queue on GPU
#else
    evalM2L(Ci,Cj);                                             // Evalaute on CPU, queue on GPU
#endif
  }

#if QUARK
  inline void interact(C_iter Ci, C_iter Cj, Quark *quark);     //!< interact() function using QUARK
#endif

//! Traverse a pair of trees using a queue
  void traverseQueue(Pair pair) {
    PairQueue pairQueue;                                        // Queue of interacting cell pairs
#if QUARK
    Quark *quark = QUARK_New(4);                                // Initialize QUARK object
#endif
    pairQueue.push(pair);                                       // Push pair to queue
    while( !pairQueue.empty() ) {                               // While dual traversal queue is not empty
      pair = pairQueue.front();                                 //  Get interaction pair from front of queue
      pairQueue.pop();                                          //  Pop dual traversal queue
      if(splitFirst(pair.first,pair.second)) {                  //  If first cell is larger
        C_iter C = pair.first;                                  //   Split the first cell
        for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {// Loop over first cell's children
          interact(Ci,pair.second,pairQueue);                   //    Calculate interaction between cells
        }                                                       //   End loop over fist cell's children
      } else {                                                  //  Else if second cell is larger
        C_iter C = pair.second;                                 //   Split the second cell
        for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {// Loop over second cell's children
          interact(pair.first,Cj,pairQueue);                    //    Calculate interaction betwen cells
        }                                                       //   End loop over second cell's children
      }                                                         //  End if for which cell to split
#if QUARK
      if( pairQueue.size() > 100 ) {                            //  When queue size reaches threshold
        while( !pairQueue.empty() ) {                           //   While dual traversal queue is not empty
          pair = pairQueue.front();                             //    Get interaction pair from front of queue
          pairQueue.pop();                                      //    Pop dual traversal queue
          interact(pair.first,pair.second,quark);               //    Schedule interact() task on QUARK
        }                                                       //   End while loop for dual traversal queue
      }                                                         //  End if for queue size
#endif
    }                                                           // End while loop for dual traversal queue
#if QUARK
    QUARK_Delete(quark);                                        // Delete QUARK object 
    writeTrace();                                               // Write event trace to file
#endif
  }

//! Get range of periodic images
  int getPeriodicRange() {
    int prange = 0;                                             //  Range of periodic images
    for( int i=0; i!=IMAGES; ++i ) {                            //  Loop over periodic image sublevels
      prange += int(pow(3,i));                                  //   Accumulate range of periodic images
    }                                                           //  End loop over perioidc image sublevels
    return prange;                                              // Return range of periodic images
  }

protected:
//! Get level from cell index
  int getLevel(bigint index) {
    int i = index;                                              // Copy to dummy index
    int level = -1;                                             // Initialize level counter
    while( i >= 0 ) {                                           // While cell index is non-negative
      level++;                                                  //  Increment level
      i -= 1 << 3*level;                                        //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

  void timeKernels();                                           //!< Time all kernels for auto-tuning

//! Upward phase for periodic cells
  void upwardPeriodic(Cells &jcells) {
    Cells pccells, pjcells;                                     // Periodic center cell and jcell
    pccells.push_back(jcells.back());                           // Root cell is first periodic cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      Cell cell;                                                //  New periodic cell at next sublevel
      C_iter C = pccells.end() - 1;                             //  Set previous periodic center cell as source
      pccells.pop_back();                                       //  Pop periodic center cell from vector
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              for( int cx=-1; cx<=1; ++cx ) {                   //      Loop over x periodic direction (child)
                for( int cy=-1; cy<=1; ++cy ) {                 //       Loop over y periodic direction (child)
                  for( int cz=-1; cz<=1; ++cz ) {               //        Loop over z periodic direction (child)
                    cell.X[0]  = C->X[0] + (ix * 6 + cx * 2) * C->R;//     Set new x coordinate for periodic image
                    cell.X[1]  = C->X[1] + (iy * 6 + cy * 2) * C->R;//     Set new y cooridnate for periodic image
                    cell.X[2]  = C->X[2] + (iz * 6 + cz * 2) * C->R;//     Set new z coordinate for periodic image
                    cell.M     = C->M;                          //         Copy multipoles to new periodic image
                    cell.NDLEAF = cell.NCHILD = 0;              //         Initialize NDLEAF & NCHILD
                    jcells.push_back(cell);                     //         Push cell into periodic jcell vector
                  }                                             //        End loop over z periodic direction (child)
                }                                               //       End loop over y periodic direction (child)
              }                                                 //      End loop over x periodic direction (child)
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            cell.X[0] = C->X[0] + ix * 2 * C->R;                //     Set new x coordinate for periodic image
            cell.X[1] = C->X[1] + iy * 2 * C->R;                //     Set new y cooridnate for periodic image
            cell.X[2] = C->X[2] + iz * 2 * C->R;                //     Set new z coordinate for periodic image
            cell.M = C->M;                                      //     Copy multipoles to new periodic image
            pjcells.push_back(cell);                            //     Push cell into periodic jcell vector
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      cell.X = C->X;                                            //  This is the center cell
      cell.R = 3 * C->R;                                        //  The cell size increases three times
      pccells.push_back(cell);                                  //  Push cell into periodic cell vector
      C_iter Ci = pccells.end() - 1;                            //  Set current cell as target for M2M
      Ci->CHILD = 0;                                            //  Set child cells for periodic M2M
      Ci->NCHILD = 27;                                          //  Set number of child cells for periodic M2M
      evalM2M(pccells,pjcells);                                 // Evaluate periodic M2M kernels for this sublevel
      pjcells.clear();                                          // Clear periodic jcell vector
    }                                                           // End loop over sublevels of tree
  }

//! Traverse tree for periodic cells
  void traversePeriodic(Cells &cells, Cells &jcells) {          
    Xperiodic = 0;                                              // Set periodic coordinate offset
    Iperiodic = Icenter;                                        // Set periodic flag to center
    C_iter Cj = jcells.end()-1;                                 // Initialize iterator for periodic source cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      for( int I=0; I!=26*27; ++I, --Cj ) {                     //  Loop over periodic images (exclude center)
        C_iter Ci = cells.end() - 1;                            //   Set root cell as target iterator
        evalM2L(Ci,Cj);                                         //   Perform M2P kernel
      }                                                         //  End loop over x periodic direction
    }                                                           // End loop over sublevels of tree
  }

public:
//! Constructor
  Evaluator() : Icenter(1 << 13) {}
//! Destructor
  ~Evaluator() {}

//! Add single list for kernel unit test
  void addM2L(C_iter Cj) {
    listM2L.resize(1);                                          // Resize vector of M2L interation lists
    flagM2L.resize(1);                                          // Resize vector of M2L periodic image flags
    listM2L[0].push_back(Cj);                                   // Push single cell into list
    flagM2L[0][Cj] |= Icenter;                                  // Flip bit of periodic image flag
  }

//! Add single list for kernel unit test
  void addM2P(C_iter Cj) {
    listM2P.resize(1);                                          // Resize vector of M2P interation lists
    flagM2P.resize(1);                                          // Resize vector of M2L periodic image flags
    listM2P[0].push_back(Cj);                                   // Push single cell into list
    flagM2P[0][Cj] |= Icenter;                                  // Flip bit of periodic image flag
  }

//! Create periodic images of bodies
  Bodies periodicBodies(Bodies &bodies) {
    Bodies jbodies;                                             // Vector for periodic images of bodies
    int prange = getPeriodicRange();                            // Get range of periodic images
    for( int ix=-prange; ix<=prange; ++ix ) {                   // Loop over x periodic direction
      for( int iy=-prange; iy<=prange; ++iy ) {                 //  Loop over y periodic direction
        for( int iz=-prange; iz<=prange; ++iz ) {               //   Loop over z periodic direction
          for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {//    Loop over bodies
            Body body = *B;                                     //     Copy current body
            body.X[0] += ix * 2 * R0;                           //     Shift x position
            body.X[1] += iy * 2 * R0;                           //     Shift y position
            body.X[2] += iz * 2 * R0;                           //     Shift z position
            jbodies.push_back(body);                            //     Push shifted body into jbodies
          }                                                     //    End loop over bodies
        }                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
    return jbodies;                                             // Return vector for periodic images of bodies
  }

//! Use multipole acceptance criteria to determine whether to approximate, do P2P, or subdivide
  void interact(C_iter Ci, C_iter Cj, PairQueue &pairQueue) {
    vect dX = Ci->X - Cj->X - Xperiodic;                        // Distance vector from source to target
    real Rq = std::sqrt(norm(dX));                              // Scalar distance
    if( Rq * THETA > Ci->R + Cj->R ) {                          // If distance if far enough
      approximate(Ci,Cj);                                       //  Use approximate kernels, e.g. M2L, M2P
    } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {             // Else if both cells are leafs
//      evalP2P(Ci,Cj);                                           //  Use P2P
    } else {                                                    // If cells are close but not leafs
      Pair pair(Ci,Cj);                                         //  Form a pair of cell iterators
      pairQueue.push(pair);                                     //  Push pair to queue
    }                                                           // End if for multipole acceptance
  }

//! Dual tree traversal
  void traverse(Cells &cells, Cells &jcells) {
    C_iter root = cells.end() - 1;                              // Iterator for root target cell
    C_iter jroot = jcells.end() - 1;                            // Iterator for root source cell
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      jroot = jcells.end() - 1 - 26 * 27 * (IMAGES - 1);        //  The root is not at the end
    }                                                           // Endif for periodic boundary condition
    Ci0 = cells.begin();                                        // Set begin iterator for target cells
    Cj0 = jcells.begin();                                       // Set begin iterator for source cells
    listM2L.resize(cells.size());                               // Resize M2L interaction list
    listM2P.resize(cells.size());                               // Resize M2P interaction list
    listP2P.resize(cells.size());                               // Resize P2P interaction list
    flagM2L.resize(cells.size());                               // Resize M2L periodic image flag
    flagM2P.resize(cells.size());                               // Resize M2P periodic image flag
    flagP2P.resize(cells.size());                               // Resize P2P periodic image flag
    if( IMAGES == 0 ) {                                         // If free boundary condition
      Iperiodic = Icenter;                                      //  Set periodic image flag to center
      Xperiodic = 0;                                            //  Set periodic coordinate offset
      Pair pair(root,jroot);                                    //  Form pair of root cells
      traverseQueue(pair);                                      //  Traverse a pair of trees
    } else {                                                    // If periodic boundary condition
      int I = 0;                                                //  Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //    Loop over z periodic direction
            Iperiodic = 1 << I;                                 //     Set periodic image flag
            Xperiodic[0] = ix * 2 * R0;                         //     Coordinate offset for x periodic direction
            Xperiodic[1] = iy * 2 * R0;                         //     Coordinate offset for y periodic direction
            Xperiodic[2] = iz * 2 * R0;                         //     Coordinate offset for z periodic direction
            Pair pair(root,jroot);                              //     Form pair of root cells
            traverseQueue(pair);                                //     Traverse a pair of trees
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {   //  Loop over target cells
        listM2L[Ci-Ci0].sort();                                 //  Sort interaction list
        listM2L[Ci-Ci0].unique();                               //  Eliminate duplicate periodic entries
        listM2P[Ci-Ci0].sort();                                 //  Sort interaction list
        listM2P[Ci-Ci0].unique();                               //  Eliminate duplicate periodic entries
        listP2P[Ci-Ci0].sort();                                 //  Sort interaction list
        listP2P[Ci-Ci0].unique();                               //  Eliminate duplicate periodic entries
      }                                                         //  End loop over target cells
    }                                                           // Endif for periodic boundary condition
  }

  void setSourceBody();                                         //!< Set source buffer for bodies (for GPU)
  void setSourceCell(bool isM);                                 //!< Set source buffer for cells (for GPU)
  void setTargetBody(Lists lists, Maps flags);                  //!< Set target buffer for bodies (for GPU)
  void setTargetCell(Lists lists, Maps flags);                  //!< Set target buffer for cells (for GPU)
  void getTargetBody(Lists &lists);                             //!< Get body values from target buffer (for GPU)
  void getTargetCell(Lists &lists, bool isM);                   //!< Get cell values from target buffer (for GPU)
  void clearBuffers();                                          //!< Clear GPU buffers

  void evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU=false);//!< Evaluate all P2P kernels (all pairs)
  void evalP2M(Cells &cells);                                   //!< Evaluate all P2M kernels
  void evalM2M(Cells &cells, Cells &jcells);                    //!< Evaluate all M2M kernels
  void evalM2L(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalM2L(Cells &cells);                                   //!< Evaluate queued M2L kernels
  void evalM2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalM2P(Cells &cells);                                   //!< Evaluate queued M2P kernels
  void evalP2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalP2P(Cells &cells);                                   //!< Evaluate queued P2P kernels (near field)
  void evalL2L(Cells &cells);                                   //!< Evaluate all L2L kernels
  void evalL2P(Cells &cells);                                   //!< Evaluate all L2P kernels
  void evalEwaldReal(C_iter Ci, C_iter Cj);                     //!< Evaluate on CPU, queue on GPU
  void evalEwaldReal(Cells &cells);                             //!< Evaluate queued Ewald real kernels
};

#include "cpuEvaluator.h"

#undef splitFirst
#endif
