#include <algorithm>

#include "CubicCell.h"

//! calculate relative cell-index inside cubic region 
/*
  region [-reach,reach],[-reach,reach],[-reach,reach]

  TODO 
    enhanced to any shape of region
    ghost cell
 */
#if 0
size_t
relative_cubic_target(int celldiv, int reach, 
                      std::vector< SpaceVector<int> >& relative_position)
{
  relative_position.resize(0);
  for(int rz=-reach;rz<=reach;rz++){
    for(int ry=-reach;ry<=reach;ry++){
      for(int rx=-reach;rx<=reach;rx++){
        if(!((rz==0)&&(ry==0)&&(rx==0))){
          SpaceVector<int> rc(rx,ry,rz);
          relative_position.push_back(rc);
        }
      }
    }
  }
  return relative_position.size();
}
#endif

//! calculate absolute index of target cell for ci
/*
  TODO 
    ghost cell
*/
#if 0
size_t
absolute_target(int celldiv, SpaceVector<int> ci, 
                std::vector< SpaceVector<int> >& relative_position,
                std::vector<int>& target_cell, 
                ShiftCellArray& shift_cell)
{
  int shiftx, shifty, shiftz;
  target_cell.resize(0);
  for(size_t cj=0;cj<relative_position.size();cj++){
    SpaceVector<int> ac = ci + relative_position[cj];
    shiftz = 0;
    if(ac.z<0){
      ac.z+=celldiv;
      shiftz = -1;
    }
    if(ac.z>=celldiv){
      ac.z-=celldiv;
      shiftz = 1;
    }
    shifty = 0;
    if(ac.y<0){
      ac.y+=celldiv;
      shifty = -1;
    }
    if(ac.y>=celldiv){
      ac.y-=celldiv;
      shifty = 1;
    }
    shiftx = 0;
    if(ac.x<0){
      ac.x+=celldiv;
      shiftx = -1;
    }
    if(ac.x>=celldiv){
      ac.x-=celldiv;
      shiftx = 1;
    }
    int csj = ac.x+ac.y*celldiv+ac.z*celldiv*celldiv;
    target_cell.push_back(csj);
    {
      if(shiftx==1) shift_cell.plusx.push_back(csj);
      if(shiftx==-1) shift_cell.minusx.push_back(csj);
      if(shifty==1) shift_cell.plusy.push_back(csj);
      if(shifty==-1) shift_cell.minusy.push_back(csj);
      if(shiftz==1) shift_cell.plusz.push_back(csj);
      if(shiftz==-1) shift_cell.minusz.push_back(csj);
    }
  }
  return target_cell.size();
}
#endif


//! select target cells defined inside cubic region 
/*
  region [-reach,reach],[-reach,reach],[-reach,reach]
  if region is overflow(<0 or >=celldiv), it be shifted opposite side 
  and listed in shift_cell

  TODO 
    ghost cell
 */
size_t
cubic_target(int celldiv, SpaceVector<int> ci, int reach,
             std::vector<int>& send_target_cell, 
             std::vector<int>& recv_target_cell, 
             ShiftCellArray& shift_cell)
{
  int shiftx, shifty, shiftz;
  for(int rz=-reach;rz<=reach;rz++){
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
    for(int ry=-reach;ry<=reach;ry++){
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
      for(int rx=-reach;rx<=reach;rx++){
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
        if(!((rz==0)&&(ry==0)&&(rx==0))){
          int cj = x+y*celldiv+z*celldiv*celldiv;
          std::vector<int>::iterator p;
          p = std::find(send_target_cell.begin(),send_target_cell.end(),cj);
          if(p==send_target_cell.end()){
            send_target_cell.push_back(cj);
          }
          p = std::find(recv_target_cell.begin(),recv_target_cell.end(),cj);
          if(p==recv_target_cell.end()){
            recv_target_cell.push_back(cj);
            {
              if(shiftx==1) shift_cell.plusx.push_back(cj);
              if(shiftx==-1) shift_cell.minusx.push_back(cj);
              if(shifty==1) shift_cell.plusy.push_back(cj);
              if(shifty==-1) shift_cell.minusy.push_back(cj);
              if(shiftz==1) shift_cell.plusz.push_back(cj);
              if(shiftz==-1) shift_cell.minusz.push_back(cj);
            }
          }
        }
      }
    }
  }
  return recv_target_cell.size();
}

size_t
cubic_target(SpaceVector<int> celldiv, SpaceVector<int> ci, int reach,
             std::vector<int>& send_target_cell, 
             std::vector<int>& recv_target_cell, 
             ShiftCellArray& shift_cell)
{
  int shiftx, shifty, shiftz;
  for(int rz=-reach;rz<=reach;rz++){
    int z = ci.z + rz;
    shiftz = 0;
    if(z<0){
      z+=celldiv.z;
      shiftz = -1;
    }
    if(z>=celldiv.z){
      z-=celldiv.z;
      shiftz = 1;
    }
    for(int ry=-reach;ry<=reach;ry++){
      int y = ci.y + ry;
      shifty = 0;
      if(y<0){
        y+=celldiv.y;
        shifty = -1;
      }
      if(y>=celldiv.y){
        y-=celldiv.y;
        shifty = 1;
      }
      for(int rx=-reach;rx<=reach;rx++){
        int x = ci.x + rx;
        shiftx = 0;
        if(x<0){
          x+=celldiv.x;
          shiftx = -1;
        }
        if(x>=celldiv.x){
          x-=celldiv.x;
          shiftx = 1;
        }
        if(!((rz==0)&&(ry==0)&&(rx==0))){
          int cj = x+y*celldiv.x+z*celldiv.x*celldiv.y;
          std::vector<int>::iterator p;
          p = std::find(send_target_cell.begin(),send_target_cell.end(),cj);
          if(p==send_target_cell.end()){
            send_target_cell.push_back(cj);
          }
          p = std::find(recv_target_cell.begin(),recv_target_cell.end(),cj);
          if(p==recv_target_cell.end()){
            recv_target_cell.push_back(cj);
            {
              if(shiftx==1) shift_cell.plusx.push_back(cj);
              if(shiftx==-1) shift_cell.minusx.push_back(cj);
              if(shifty==1) shift_cell.plusy.push_back(cj);
              if(shifty==-1) shift_cell.minusy.push_back(cj);
              if(shiftz==1) shift_cell.plusz.push_back(cj);
              if(shiftz==-1) shift_cell.minusz.push_back(cj);
            }
          }
        }
      }
    }
  }
  return recv_target_cell.size();
}

