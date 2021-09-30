/*!
  Provide : conversion method node ID(rank) and geometry

  In this version
  support only 3D Cartesian coordinate
  XYZ order : node ID of (X,Y,Z) is  X + size_of_x*(Y + size_of_y*Z)
  X, Y, Z and node ID start 0

  node ID 0 1 2 ... sx-1 sx sx+1 ... sx*sy-1 sx*sy ...sx*sy*sz-1
  X       0 1 2 ... sx-1  0    1 ...    sx-1     0 ...      sx-1
  Y       0 0 0 ...    0  1    1 ...    sy-1     0 ...      sy-1
  Z       0 0 0 ...    0  0    0 ...       0     1 ...      sz-1

*/

#include "Geometry.h"
#include <cstdlib>

void Geometry::calc_have_cell()
{
  if(size.x*size.y*size.z!=0){
    have_cell.x = full_cell.x/size.x;
    have_cell.y = full_cell.y/size.y;
    have_cell.z = full_cell.z/size.z;
  }else{
    have_cell.x = 0;
    have_cell.y = 0;
    have_cell.z = 0;
  }
}

void Geometry::set_size(SpaceVector<int> _size)
{
  size = _size;
  calc_have_cell();
}

void Geometry::set_full_cell(SpaceVector<int> _full_cell)
{
  full_cell = _full_cell;
  calc_have_cell();
}

SpaceVector<int> Geometry::getAbsolutePosition(int node_id)
{
  SpaceVector<int> pos;

  pos.z = node_id/(size.x*size.y);
  pos.y = (node_id-size.x*size.y*pos.z)/size.x;
  pos.x = node_id-size.x*(pos.y+size.y*pos.z);
  
  return pos;
}

SpaceVector<int> Geometry::getRelativePosition(int target_node, 
                                               int center_node)
{
  SpaceVector<int> pos_t, pos_c;

  pos_t = getAbsolutePosition(target_node);
  pos_c = getAbsolutePosition(center_node);
  return pos_t-pos_c;
}

SpaceVector<int> Geometry::getPeriodicRelativePosition(int target_node, 
                                                          int center_node)
{
  SpaceVector<int> pos_t, pos_c;

  pos_t = getAbsolutePosition(target_node);
  pos_c = getAbsolutePosition(center_node);
  pos_t = pos_t - pos_c;
  
  periodic_relative_shift(pos_t.x,size.x);
  periodic_relative_shift(pos_t.y,size.y);
  periodic_relative_shift(pos_t.z,size.z);
  
  return pos_t;
}

int Geometry::getNodeID(SpaceVector<int> position)
{
  return position.x + size.x*(position.y + size.y*position.z);
}

int Geometry::getNodeIDofRelativePosition(SpaceVector<int> position,
                                            int center_node)
{
  SpaceVector<int> pc = getAbsolutePosition(center_node);
  SpaceVector<int> pa = pc + position;
  periodic_shift(pa.x,size.x);
  periodic_shift(pa.y,size.y);
  periodic_shift(pa.z,size.z);
  
  return getNodeID(pa);
}



int GeometryXYZ::getSendPath(int target_node, int center_node)
{
  SpaceVector<int> rpos;
  SpaceVector<int> path_pos = getAbsolutePosition(center_node);

  if(periodic){
    rpos = getPeriodicRelativePosition(target_node,center_node);
  }else{
    rpos = getRelativePosition(target_node,center_node);
  }

  if(rpos.x>0){
    path_pos.x += 1;
  }else if(rpos.x<0){
    path_pos.x -= 1;
  }else if(rpos.y>0){
    path_pos.y += 1;
  }else if(rpos.y<0){
    path_pos.y -= 1;
  }else if(rpos.z>0){
    path_pos.z += 1;
  }else if(rpos.z<0){
    path_pos.z -= 1;
  }else{
    // self
  }

  periodic_shift(path_pos.x,size.x);
  periodic_shift(path_pos.y,size.y);
  periodic_shift(path_pos.z,size.z);

  return getNodeID(path_pos);
}

int GeometryXYZ::getReceivePath(int target_node, int center_node)
{
  SpaceVector<int> rpos;
  SpaceVector<int> path_pos = getAbsolutePosition(center_node);

  if(periodic){
    rpos = getPeriodicRelativePosition(target_node,center_node);
  }else{
    rpos = getRelativePosition(target_node,center_node);
  }

  if(rpos.z<0){
    path_pos.z -= 1;
  }else if(rpos.z>0){
    path_pos.z += 1;
  }else if(rpos.y<0){
    path_pos.y -= 1;
  }else if(rpos.y>0){
    path_pos.y += 1;
  }else if(rpos.x<0){
    path_pos.x -= 1;
  }else if(rpos.z>0){
    path_pos.x += 1;
  }else{
    // self
  }

  periodic_shift(path_pos.x,size.x);
  periodic_shift(path_pos.y,size.y);
  periodic_shift(path_pos.z,size.z);

  return getNodeID(path_pos);
}

void GeometryXYZ::path_depth(SpaceVector<int> &plus_depth,
                             SpaceVector<int> &minus_depth,
                             std::vector<int> target, int center_node)
{
  for(std::vector<int>::size_type i=0;i<target.size();i++){
    SpaceVector<int> rp = getPeriodicRelativePosition(target[i],center_node);
    if(rp.x<minus_depth.x)minus_depth.x=rp.x;
    if(rp.x>plus_depth.x)plus_depth.x=rp.x;
    if(rp.y<minus_depth.y)minus_depth.y=rp.y;
    if(rp.y>plus_depth.y)plus_depth.y=rp.y;
    if(rp.z<minus_depth.z)minus_depth.z=rp.z;
    if(rp.z>plus_depth.z)plus_depth.z=rp.z;
  }
}

void GeometryXYZ::pushtarget(std::vector< std::vector<int> > SendTarget, 
                             std::vector< std::vector<int> > ReceiveTarget,
                             int offset,
                             int splus_depth, int sminus_depth,
                             int rplus_depth, int rminus_depth,
                             int plus, int minus)
{
  for(int s=0;s<splus_depth;s++){
    SendTarget[offset+s].push_back(plus);
  }
  for(int s=0;s<-rminus_depth;s++){
    ReceiveTarget[offset+s].push_back(minus);
  }
  for(int s=0;s<-sminus_depth;s++){
    SendTarget[offset+s].push_back(minus);
  }
  for(int s=0;s<rplus_depth;s++){
    ReceiveTarget[offset+s].push_back(plus);
  }
}

int GeometryXYZ::getCommStage(std::vector< std::vector<int> > SendTarget, 
                              std::vector< std::vector<int> > ReceiveTarget,
                              std::vector<int> allSendTarget,
                              std::vector<int> allReceiveTarget,
                              int center_node)
{
  SpaceVector<int> sminus_depth(0,0,0);
  SpaceVector<int> splus_depth(0,0,0);
  SpaceVector<int> rminus_depth(0,0,0);
  SpaceVector<int> rplus_depth(0,0,0);
  SpaceVector<int> minus_depth(0,0,0);
  SpaceVector<int> plus_depth(0,0,0);

  path_depth(splus_depth,sminus_depth,allSendTarget,center_node);
  path_depth(rplus_depth,rminus_depth,allReceiveTarget,center_node);

  if((splus_depth.x!=-rminus_depth.x)||(sminus_depth.x!=-rplus_depth.x)
     ||(splus_depth.y!=-rminus_depth.y)||(sminus_depth.y!=-rplus_depth.y)
     ||(splus_depth.z!=-rminus_depth.z)||(sminus_depth.z!=-rplus_depth.z)){
    std::cout << "reverse receieve depth != send depth" << std::endl;
  }
  
  minus_depth.x = std::max(std::abs(sminus_depth.x),rplus_depth.x);
  minus_depth.y = std::max(std::abs(sminus_depth.y),rplus_depth.y);
  minus_depth.z = std::max(std::abs(sminus_depth.z),rplus_depth.z);
  plus_depth.x = std::max(splus_depth.x,std::abs(rminus_depth.x));
  plus_depth.y = std::max(splus_depth.y,std::abs(rminus_depth.y));
  plus_depth.z = std::max(splus_depth.z,std::abs(rminus_depth.z));

  SpaceVector<int> plusx(1,0,0);
  SpaceVector<int> minusx(-1,0,0);
  SpaceVector<int> plusy(0,1,0);
  SpaceVector<int> minusy(0,-1,0);
  SpaceVector<int> plusz(0,0,1);
  SpaceVector<int> minusz(0,0,-1);
  int plusx_node = getNodeIDofRelativePosition(plusx,center_node);
  int minusx_node = getNodeIDofRelativePosition(minusx,center_node);
  int plusy_node = getNodeIDofRelativePosition(plusy,center_node);
  int minusy_node = getNodeIDofRelativePosition(minusy,center_node);
  int plusz_node = getNodeIDofRelativePosition(plusz,center_node);
  int minusz_node = getNodeIDofRelativePosition(minusz,center_node);
  SendTarget.clear();
  ReceiveTarget.clear();

  int depth_x = std::max(plus_depth.x,minus_depth.x);
  int depth_y = std::max(plus_depth.y,minus_depth.y);
  int depth_z = std::max(plus_depth.z,minus_depth.z);
  SendTarget.resize(depth_x+depth_y+depth_z);
  ReceiveTarget.resize(depth_x+depth_y+depth_z);
  pushtarget(SendTarget,ReceiveTarget,0,
             splus_depth.x,sminus_depth.x,rplus_depth.x,rminus_depth.x,
             plusx_node,minusx_node);
  pushtarget(SendTarget,ReceiveTarget,0,
             splus_depth.y,sminus_depth.y,rplus_depth.y,rminus_depth.y,
             plusy_node,minusy_node);
  pushtarget(SendTarget,ReceiveTarget,0,
             splus_depth.z,sminus_depth.z,rplus_depth.z,rminus_depth.z,
             plusz_node,minusz_node);

  return 1;   // ???
}
