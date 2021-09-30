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
#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "tree.h"

//! Handles all the partitioning of domains
class Partition : public MyMPI, public TreeStructure {
private:
  int numCells1D;                                               //!< Number of cells in one dimension (leaf level)

protected:
  int LEVEL;                                                    //!< Level of the MPI process binary tree
  std::vector<vect> XMIN;                                       //!< Minimum position vector of bodies
  std::vector<vect> XMAX;                                       //!< Maximum position vector of bodies
  int nprocs[64][2];                                            //!< Number of processes in the two split groups
  int offset[64][2];                                            //!< Offset of body in the two split groups
  int  color[64][3];                                            //!< Color of hypercube communicators
  int    key[64][3];                                            //!< Key of hypercube communicators
  MPI_Comm MPI_COMM[64][3];                                     //!< Hypercube communicators

public:
//! Constructor
  Partition() : TreeStructure() {
    LEVEL = int(log(MPISIZE) / M_LN2 - 1e-5) + 1;               // Level of the process binary tree
    if(MPISIZE == 1) LEVEL = 0;                                 // Level is 0 for a serial execution
    XMIN.resize(LEVEL+1);                                       // Minimum position vector at each level
    XMAX.resize(LEVEL+1);                                       // Maximum position vector at each level
    startTimer("Split comm   ");                                // Start timer
    nprocs[0][0] = nprocs[0][1] = MPISIZE;                      // Initialize number of processes in groups
    offset[0][0] = offset[0][1] = 0;                            // Initialize offset of body in groups
     color[0][0] =  color[0][1] =  color[0][2] = 0;             // Initialize color of communicators
       key[0][0] =    key[0][1] =    key[0][2] = 0;             // Initialize key of communicators
    stopTimer("Split comm   ",printNow);                        // Stop timer 
  }
//! Destructor
  ~Partition() {}

//! Set bounds of domain to be partitioned
  void setGlobDomain(Bodies &bodies, vect x0, real r0) {
    X0 = x0;                                                    //  Center is [0, 0, 0]
    R0 = r0;                                                    //  Radius is M_PI
    XMAX[0] = X0 + R0;                                          // Reposition global maximum
    XMIN[0] = X0 - R0;                                          // Reposition global minimum
  }

//! Partition by recursive octsection
  void octsection(Bodies &bodies) {
    startTimer("Partition    ");                                // Start timer
    for( int l=0; l!=LEVEL; ++l ) {                             // Loop over level of N-D hypercube
      int d = 2 - l % 3;                                        //  Dimension of subdivision
      XMIN[l+1] = XMIN[l];                                      //  Set XMAX for next subdivision
      XMAX[l+1] = XMAX[l];                                      //  Set XMIN for next subdivision
      if( (MPIRANK >> (LEVEL - l - 1)) % 2 ) {                  //  If on left side
        XMIN[l+1][d] = (XMAX[l][d]+XMIN[l][d]) / 2;             //   Set XMIN to midpoint
      } else {                                                  //  If on right side
        XMAX[l+1][d] = (XMAX[l][d]+XMIN[l][d]) / 2;             //   Set XMAX to midpoint
      }                                                         //  Endif for side
    }                                                           // End loop over levels
    stopTimer("Partition    ",printNow);                        // Stop timer 
  }

};

#endif
