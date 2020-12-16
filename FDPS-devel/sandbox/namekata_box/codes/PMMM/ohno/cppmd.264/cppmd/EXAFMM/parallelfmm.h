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
#ifndef parallelfmm_h
#define parallelfmm_h
#include "partition.h"

//! Handles all the communication of local essential trees
class ParallelFMM : public Partition {
private:
  std::vector<int>    sendCellCnt;                              //!< Vector of cell send counts
  std::vector<int>    sendCellDsp;                              //!< Vector of cell send displacements
  std::vector<int>    recvCellCnt;                              //!< Vector of cell recv counts
  std::vector<int>    recvCellDsp;                              //!< Vector of cell recv displacements
  std::vector<vect>   xminAll;                                  //!< Buffer for gathering XMIN
  std::vector<vect>   xmaxAll;                                  //!< Buffer for gathering XMAX
  JCells  sendCells;                                            //!< Send buffer for cells
  JCells  recvCells;                                            //!< Recv buffer for cells

private:
//! Gather bounds of other domain
  void gatherBounds() {
    xminAll.resize(MPISIZE);                                    // Resize buffer for gathering xmin
    xmaxAll.resize(MPISIZE);                                    // Resize buffer for gathering xmax
    sendCellCnt.resize(MPISIZE);                                // Resize vector of cell send counts
    sendCellDsp.resize(MPISIZE);                                // Resize vector of cell send displacements
    recvCellCnt.resize(MPISIZE);                                // Resize vector of cell recv counts
    recvCellDsp.resize(MPISIZE);                                // Resize vector of cell recv displacements
    MPI_Datatype MPI_TYPE = getType(XMIN[LEVEL][0]);            // Get MPI data type
    MPI_Allgather(&XMIN[LEVEL][0],3,MPI_TYPE,                   // Gather XMIN
                  &xminAll[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    MPI_Allgather(&XMAX[LEVEL][0],3,MPI_TYPE,                   // Gather XMAX
                  &xmaxAll[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
  }

//! Get disatnce to other domain
  real getDistance(C_iter C, vect xmin, vect xmax) {
    vect dist;                                                  // Distance vector
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      dist[d] = (C->X[d] + Xperiodic[d] > xmax[d])*             //  Calculate the distance between cell C and
                (C->X[d] + Xperiodic[d] - xmax[d])+             //  the nearest point in domain [xmin,xmax]^3
                (C->X[d] + Xperiodic[d] < xmin[d])*             //  Take the differnece from xmin or xmax
                (C->X[d] + Xperiodic[d] - xmin[d]);             //  or 0 if between xmin and xmax
    }                                                           // End loop over dimensions
    real R = std::sqrt(norm(dist));                             // Scalar distance
    return R;
  }

//! Determine which cells to send
  void getLET(C_iter C0, C_iter C, vect xmin, vect xmax) {
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;            // Level of local root cell
    if( MPISIZE == 1 ) level = 0;                               // Account for serial case
    for( int i=0; i!=C->NCHILD; i++ ) {                         // Loop over child cells
      C_iter CC = C0+C->CHILD+i;                                //  Iterator for child cell
      bool divide = false;                                      //  Initialize logical for dividing
      if( IMAGES == 0 ) {                                       //  If free boundary condition
        Xperiodic = 0;                                          //   Set periodic coordinate offset
        real R = getDistance(CC,xmin,xmax);                     //   Get distance to other domain
        divide |= CLET * CC->R > THETA * R - EPS2;              //   If the cell seems too close and not twig
      } else {                                                  //  If periodic boundary condition
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz ) {                     //     Loop over z periodic direction
              Xperiodic[0] = ix * 2 * R0;                       //      Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //      Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //      Coordinate offset for z periodic direction
              real R = getDistance(CC,xmin,xmax);               //      Get distance to other domain
              divide |= CLET * CC->R > THETA * R - EPS2;        //      If the cell seems too close and not twig
            }                                                   //     End loop over z periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over x periodic direction
      }                                                         //  Endif for periodic boundary condition
      divide |= R0 / (1 << level) + 1e-5 < CC->R;               //  If the cell is larger than the local root cell
      if( divide && CC->NCHILD != 0 ) {                         //  If the cell seems too close and not twig
        getLET(C0,CC,xmin,xmax);                                //   Traverse the tree further
      } else {                                                  //  If the cell is far or a twig
        assert( R0 / (1 << level) + 1e-5 > CC->R );             //   Can't send cells that are larger than local root
        JCell cell;                                             //   Set compact cell type for sending
        cell.ICELL = CC->ICELL;                                 //   Set index of compact cell type
        cell.M     = CC->M;                                     //   Set Multipoles of compact cell type
        sendCells.push_back(cell);                              //    Push cell into send buffer vector
      }                                                         //  Endif for interaction
    }                                                           // End loop over child cells
    if( C->ICELL == 0 && C->NCHILD == 0 ) {                     // If the root cell has no children
      JCell cell;                                               //  Set compact cell type for sending
      cell.ICELL = C->ICELL;                                    //  Set index of compact cell type
      cell.M     = C->M;                                        //  Set Multipoles of compact cell type
      sendCells.push_back(cell);                                //  Push cell into send buffer vector
    }                                                           // Endif for root cells children
  }

//! Turn cells to twigs
  void cells2twigs(Cells &cells, Cells &twigs, bool last) {
    while( !cells.empty() ) {                                   // While cell vector is not empty
      if( cells.back().NCHILD == 0 ) {                          //  If cell has no child
        if( cells.back().NDLEAF == 0 || !last ) {               //   If cell has no leaf or is not last iteration
          cells.back().NDLEAF = 0;                              //    Set number of leafs to 0
          twigs.push_back(cells.back());                        //    Push cell into twig vector
        }                                                       //   Endif for no leaf
      }                                                         //  Endif for no child
      cells.pop_back();                                         //  Pop last element from cell vector
    }                                                           // End while for cell vector
  }

//! Turn recv buffer to twigs
  void recv2twigs(Bodies &bodies, Cells &twigs) {
    for( JC_iter JC=recvCells.begin(); JC!=recvCells.end(); ++JC ) {// Loop over recv buffer
      Cell cell;                                                //  Cell structure
      cell.ICELL = JC->ICELL;                                   //  Set index of cell
      cell.M     = JC->M;                                       //  Set multipole of cell
      cell.NDLEAF = cell.NCHILD = 0;                            //  Set number of leafs and children
      cell.LEAF  = bodies.end();                                //  Set pointer to first leaf
      getCenter(cell);                                          //  Set center and radius
      twigs.push_back(cell);                                    //  Push cell into twig vector
    }                                                           // End loop over recv buffer
  }

//! Zip two groups of twigs that overlap
  void zipTwigs(Cells &twigs, Cells &cells, Cells &sticks, bool last) {
    startTimer("Sort resize  ");                                // Start timer
    Cells cbuffer = twigs;                                      // Sort buffer for cells
    stopTimer("Sort resize  ",printNow);                        // Stop timer 
    sortCells(twigs,cbuffer);                                   // Sort twigs in ascending order
    startTimer("Ziptwigs     ");                                // Start timer
    bigint index = -1;                                          // Initialize index counter
    while( !twigs.empty() ) {                                   // While twig vector is not empty
      if( twigs.back().ICELL != index ) {                       //  If twig's index is different from previous
        cells.push_back(twigs.back());                          //   Push twig into cell vector
        index = twigs.back().ICELL;                             //   Update index counter
      } else if ( twigs.back().NDLEAF == 0 || !last ) {         //  Elseif twig-twig collision
        cells.back().M += twigs.back().M;                       //   Accumulate the multipole
      } else if ( cells.back().NDLEAF == 0 ) {                  //  Elseif twig-body collision
        Mset M;                                                 //   Multipole for temporary storage
        M = cells.back().M;                                     //   Save multipoles from cells
        cells.back() = twigs.back();                            //   Copy twigs to cells
        cells.back().M = M;                                     //   Copy back multipoles to cells
        twigs.back().M = M - twigs.back().M;                    //   Take the difference of the two
        if( std::abs(twigs.back().M[0]/M[0]) > EPS ) {          //   If the difference is non-zero
          sticks.push_back(twigs.back());                       //    Save this difference in the sticks vector
        }                                                       //   Endif for non-zero difference
      } else {                                                  //  Else body-body collision (don't do anything)
      }                                                         //  Endif for collision type
      twigs.pop_back();                                         //  Pop last element from twig vector
    }                                                           // End while for twig vector
    stopTimer("Ziptwigs     ",printNow);                        // Stop timer 
    sortCells(cells,cbuffer);                                   // Sort cells in ascending order
    startTimer("Ziptwigs     ");                                // Start timer
    twigs = cells;                                              // Copy cells to twigs
    cells.clear();                                              // Clear cells
    stopTimer("Ziptwigs     ",printNow);                        // Stop timer 
  }

//! Re-index bodies
  void reindexBodies(Bodies &bodies, Cells &twigs, Cells &cells ,Cells &sticks) {
    startTimer("Reindex      ");                                // Start timer
    while( !twigs.empty() ) {                                   // While twig vector is not empty
      if( twigs.back().NDLEAF == 0 ) {                          //  If twig has no leafs
        cells.push_back(twigs.back());                          //   Push twig into cell vector
      }                                                         //  Endif for no leafs
      twigs.pop_back();                                         //  Pop last element from twig vector
    }                                                           // End while for twig vector
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Reindex      ",printNow);                        // Stop timer 
    sortBodies(bodies,buffer,false);                              // Sort bodies in descending order
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs
    startTimer("Reindex      ");                                // Start timer
    for( C_iter C=twigs.begin(); C!=twigs.end(); ++C ) {        // Loop over cells
      if( sticks.size() > 0 ) {                                 //  If stick vector is not empty
        if( C->ICELL == sticks.back().ICELL ) {                 //   If twig's index is equal to stick's index
          C->M += sticks.back().M;                              //    Accumulate multipole
          sticks.pop_back();                                    //    Pop last element from stick vector
        }                                                       //   Endif for twig's index
      }                                                         //  Endif for stick vector
    }                                                           // End loop over cells
    cells.insert(cells.begin(),twigs.begin(),twigs.end());      // Add twigs to the end of cell vector
    cells.insert(cells.begin(),sticks.begin(),sticks.end());    // Add remaining sticks to the end of cell vector
    sticks.clear();                                             // Clear sticks
    Cells cbuffer = cells;                                      // Sort buffer for cells
    stopTimer("Reindex      ",printNow);                        // Stop timer 
    sortCells(cells,cbuffer);                                   // Sort cells in ascending order
    startTimer("Reindex      ");                                // Start timer
    twigs = cells;                                              // Copy cells to twigs
    cells.clear();                                              // Clear cells
    stopTimer("Reindex      ",printNow);                        // Stop timer 
  }

public:
//! Constructor
  ParallelFMM() : Partition() {}
//! Destructor
  ~ParallelFMM() {}

  void gatherRootCell(Cells &cells) {
    C_iter C = cells.end() - 1;
    JCell cell;
    cell.ICELL = C->ICELL;
    cell.M = C->M;
    sendCells.push_back(cell);
    int bytes = sizeof(sendCells[0]);
    recvCells.resize(bytes*MPISIZE); 
    MPI_Allgather(&sendCells[0], bytes, MPI_BYTE,                      
                  &recvCells[0], bytes, MPI_BYTE,                
                  MPI_COMM_WORLD);
  }

//! Communicate cells in the local essential tree
  void commCells(Bodies &bodies, Cells &cells) {
    vect xmin = 0, xmax = 0;                                    // Initialize domain boundaries
    Cells twigs,sticks;                                         // Twigs and sticks are special types of cells
    startTimer("Gather bounds");                                // Start timer
    gatherBounds();                                             // Gather bounds of other domain
    stopTimer("Gather bounds",printNow);                        // Stop timer 
    startTimer("Get LET      ");                                // Start timer
    int ssize = 0;                                              // Initialize offset for send cells
    sendCellCnt.assign(MPISIZE,0);                              // Initialize cell send count
    sendCellDsp.assign(MPISIZE,0);                              // Initialize cell send displacement
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks to send to
      getLET(cells.begin(),cells.end()-1,xminAll[irank],xmaxAll[irank]);//  Determine which cells to send
      sendCellCnt[irank] = sendCells.size()-ssize;              //  Set cell send count of current rank
      sendCellDsp[irank] = ssize;                               //  Set cell send displacement of current rank
      ssize += sendCellCnt[irank];                              //  Increment offset for vector send cells
    }                                                           // End loop over ranks
    stopTimer("Get LET      ",printNow);                        // Stop timer 
    startTimer("Alltoall C   ");                                // Start timer
    MPI_Alltoall(&sendCellCnt[0],1,MPI_INT,&recvCellCnt[0],1,MPI_INT,MPI_COMM_WORLD);// Communicate the send counts
    int rsize = 0;                                              // Initialize total recv count
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks to recv from
      recvCellDsp[irank] = rsize;                               //  Set recv displacements
      rsize += recvCellCnt[irank];                              //  Accumulate recv counts
    }                                                           // End loop over ranks to recv from
    recvCells.resize(rsize);                                    // Resize recv buffer
    int bytes = sizeof(sendCells[0]);                           // Byte size of jbody structure
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendCellCnt[i] *= bytes;                                  //  Multiply by bytes
      sendCellDsp[i] *= bytes;                                  //  Multiply by bytes
      recvCellCnt[i] *= bytes;                                  //  Multiply by bytes
      recvCellDsp[i] *= bytes;                                  //  Multiply by bytes
    }                                                           // End loop over ranks
    /*
    if(MPIRANK==0){
      printf("send MPI_Alltoallv");
      for( int i=0; i!=MPISIZE; ++i ) printf(" %d",sendCellCnt[i]);
      printf("\n");
      printf("recv MPI_Alltoallv");
      for( int i=0; i!=MPISIZE; ++i ) printf(" %d",recvCellCnt[i]);
      printf("\n");
    }
    */
    MPI_Alltoallv(&sendCells[0],&sendCellCnt[0],&sendCellDsp[0],MPI_BYTE,
                  &recvCells[0],&recvCellCnt[0],&recvCellDsp[0],MPI_BYTE,MPI_COMM_WORLD);
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendCellCnt[i] /= bytes;                                  //  Divide by bytes
      sendCellDsp[i] /= bytes;                                  //  Divide by bytes
      recvCellCnt[i] /= bytes;                                  //  Divide by bytes
      recvCellDsp[i] /= bytes;                                  //  Divide by bytes
    }                                                           // End loop over ranks
    stopTimer("Alltoall C   ",printNow);                        // Stop timer 
    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs
    startTimer("Cells2twigs  ");                                // Start timer
    cells2twigs(cells,twigs,true);                              // Put cells into twig vector
    stopTimer("Cells2twigs  ",printNow);                        // Stop timer 
    startTimer("Recv2twigs   ");                                // Start timer
    recv2twigs(bodies,twigs);                                   // Put recv buffer into twig vector
    stopTimer("Recv2twigs   ",printNow);                        // Stop timer 
    zipTwigs(twigs,cells,sticks,true);                          // Zip two groups of twigs that overlap
    reindexBodies(bodies,twigs,cells,sticks);                   // Re-index bodies
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
    sendCells.clear();                                          // Clear send buffer
    recvCells.clear();                                          // Clear recv buffer
  }
};

#endif
