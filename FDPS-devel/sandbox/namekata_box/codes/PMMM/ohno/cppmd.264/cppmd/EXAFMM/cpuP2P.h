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

void Kernel::P2P(C_iter Ci, C_iter Cj) const {         // Laplace P2P kernel on CPU
#ifndef SPARC_SIMD
  for( B_iter Bi=Ci->LEAF; Bi!=Ci->LEAF+Ci->NDLEAF; ++Bi ) {    // Loop over target bodies
    for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {  //  Loop over source bodies
      vect dist = Bi->X - Bj->X - Xperiodic;                    //   Distance vector from source to target
      real R2 = norm(dist) + EPS2;                              //   R^2
      real invR = 1 / std::sqrt(R2);                            //   1 / R
      if( R2 == 0 ) invR = 0;                                   //   Exclude self interaction
      real invR3 = Bj->SRC * invR * invR * invR;                //   charge / R^3
      Bi->TRG[0] += Bj->SRC * invR;                             //   potential
      Bi->TRG[1] -= dist[0] * invR3;                            //   x component of force
      Bi->TRG[2] -= dist[1] * invR3;                            //   y component of force
      Bi->TRG[3] -= dist[2] * invR3;                            //   z component of force
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
#else
  real (* cbi)[10] =  (real (*)[10])(&(Ci->LEAF->IBODY));
  real (* cbj)[10] =  (real (*)[10])(&(Cj->LEAF->IBODY));
  real xp[3] = {Xperiodic[0],Xperiodic[1],Xperiodic[2]};
  int ni = Ci->NDLEAF;
  int nj = Cj->NDLEAF;
  int i,j;
#pragma loop norecurrence
  for(i=0;i<ni;i++){
    for(j=0;j<nj;j++){
      real dist_x = cbi[i][2] - cbj[j][2] - xp[0];
      real dist_y = cbi[i][3] - cbj[j][3] - xp[1];
      real dist_z = cbi[i][4] - cbj[j][4] - xp[2];
      real R2 = EPS2+dist_x*dist_x+dist_y*dist_y+dist_z*dist_z;
      real invR = 1.0/sqrt(R2);
      if( R2 == 0 ) invR = 0;
      real invR3 = cbj[j][5] * invR * invR * invR;
      cbi[i][6] += cbj[j][5] * invR;
      cbi[i][7] -= dist_x * invR3;
      cbi[i][8] -= dist_y * invR3;
      cbi[i][9] -= dist_z * invR3;
    }
  }
#endif
}
