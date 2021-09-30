#include <algorithm>

#include "HalfShell.h"

//! select target cells defined inside cubic region 
/*
  region [-reach,reach],[-reach,reach],[-reach,reach]
  if region is overflow(<0 or >=celldiv), it be shifted opposite side 
  and listed in shift_cell

  TODO 
    ghost cell
 */
size_t
halfshell_target(int celldiv, SpaceVector<int> ci, 
                 SpaceVector<double> cellsize, SpaceVector<double> margin,
                 double cutoff,
                 std::vector<int>& send_target_cell, 
                 std::vector<int>& recv_target_cell, 
                 ShiftCellArray& shift_cell)
{
  int reachx, reachy, reachz;
  int shiftx, shifty, shiftz;
  reachx = int(ceil((cutoff+2.0*margin.x)/cellsize.x));
  reachy = int(ceil((cutoff+2.0*margin.y)/cellsize.y));
  reachz = int(ceil((cutoff+2.0*margin.z)/cellsize.z));
  for(int rz=0;rz<=reachz;rz++){
    int z = ci.z + rz;
    shiftz = 0;
    if(z<0){
      z+=celldiv;
      shiftz = -1;
    }
    if(z>=celldiv){
      z-=celldiv;
      shiftz = 1;
    }
    int sz = ci.z - rz;
    if(sz<0)sz+=celldiv;
    if(sz>=celldiv)sz-=celldiv;
    int ry;
    if(rz==0){
      ry = 0;
    }else{
      ry = -reachy;
    }
    for(;ry<=reachy;ry++){
      int y = ci.y + ry;
      shifty = 0;
      if(y<0){
        y+=celldiv;
        shifty = -1;
      }
      if(y>=celldiv){
        y-=celldiv;
        shifty = 1;
      }
      int sy = ci.y - ry;
      if(sy<0)sy+=celldiv;
      if(sy>=celldiv)sy-=celldiv;
      int rx;
      if((rz==0)&&(ry==0)){
        rx = 1;
      }else{
        rx = -reachx; 
      }
      for(;rx<=reachx;rx++){
        int x = ci.x + rx;
        shiftx = 0;
        if(x<0){
          x+=celldiv;
          shiftx = -1;
        }
        if(x>=celldiv){
          x-=celldiv;
          shiftx = 1;
        }
        int sx = ci.x - rx;
        if(sx<0)sx+=celldiv;
        if(sx>=celldiv)sx-=celldiv;
        if(min_cell_distance(SpaceVector<int>(0,0,0),
                             SpaceVector<int>(rx,ry,rz),
                             cellsize, margin)<cutoff){
          int cj = x+y*celldiv+z*celldiv*celldiv;
          recv_target_cell.push_back(cj);
          {
            if(shiftx==1) shift_cell.plusx.push_back(cj);
            if(shiftx==-1) shift_cell.minusx.push_back(cj);
            if(shifty==1) shift_cell.plusy.push_back(cj);
            if(shifty==-1) shift_cell.minusy.push_back(cj);
            if(shiftz==1) shift_cell.plusz.push_back(cj);
            if(shiftz==-1) shift_cell.minusz.push_back(cj);
          }
          int scj = sx+sy*celldiv+sz*celldiv*celldiv;
          send_target_cell.push_back(scj);
        }
      }
    }
  }
  return recv_target_cell.size();
}

