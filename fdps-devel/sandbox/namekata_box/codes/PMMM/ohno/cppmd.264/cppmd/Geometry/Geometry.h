/*!
  Provide : conversion method node ID(rank) and geometry

  In this version, support only 3D Cartesian coordinate
*/

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Common.h"

class Geometry {
public:
  SpaceVector<int> size;  //!< size rectangular parallelepiped
  bool periodic;                  //!< is periodic
  SpaceVector<int> full_cell;  //!< full cell number
  SpaceVector<int> have_cell;  //!< cell number in one node

  Geometry(){}

  Geometry(SpaceVector<int> _size, SpaceVector<int> _full_cell, bool _periodic=true)
    : size(_size),
      periodic(_periodic),
      full_cell(_full_cell)
  {
    calc_have_cell();
  }

  virtual ~Geometry(){}

  void calc_have_cell();
  void set_size(SpaceVector<int> _size);
  void set_full_cell(SpaceVector<int> _full_cell);

  virtual SpaceVector<int> getAbsolutePosition(int node_id); //! return absolute position of the node

  virtual SpaceVector<int> getRelativePosition(int target_node, int center_node); //! return relative position of target node from center_node
  
  virtual SpaceVector<int> getPeriodicRelativePosition(int target_node, int center_node); //! return relative position of target node from center_node

  virtual int getNodeIDofRelativePosition(SpaceVector<int> position, int center_node); //! return node ID at relative postition from center node (if out of bounds, correct as periodic)

  virtual int getNodeID(SpaceVector<int> position);

  virtual int getSendPath(int target_node, int center_node)
  {
    return target_node;
  }
  
  virtual int getReceivePath(int target_node, int center_node)
  {
    return target_node;
  }

  virtual int getCommStage(std::vector< std::vector<int> > SendTarget, 
                           std::vector< std::vector<int> > ReceiveTarget,
                           std::vector<int> allSendTarget,
                           std::vector<int> allReceiveTarget,
                           int center_node)
  {
    SendTarget.push_back(allSendTarget);
    ReceiveTarget.push_back(allReceiveTarget);
    return 1;
  }

  /*
    cell position to node position in x,y,z
   */
  inline int cell_x_to_node_x(int cell_x){
    //    std::cout << cell_x << "/" << have_cell.x << "->" << cell_x/have_cell.x << std::endl;
    return cell_x/have_cell.x;
  }
  inline int cell_y_to_node_y(int cell_y){
    return cell_y/have_cell.y;
  }
  inline int cell_z_to_node_z(int cell_z){
    return cell_z/have_cell.z;
  }
  inline int relative(int cell_pos, int node_pos, 
                      int have_cell, int full_cell){
    /* 
       cell position of this node = [node_pos*have_cell,(node_pos+1)*have_cell)
       relative position of cell_pos to first cell of this node = cell_pos - node_pos*have_cell
       rcell = periodic corrected relative position of cell_pos to first cell of this node
     */
    int rcell = cell_pos-node_pos*have_cell;
    periodic_relative_shift(rcell,full_cell);
    if(rcell<0){    // shift negative rcell by (have_cell-1) to give negative node position by /have_cell.
      rcell -= have_cell;
      rcell++;
    }
    return rcell/have_cell;
  }
  inline int relative_x(int cell_x, int node_x){
    return relative(cell_x,node_x,have_cell.x,full_cell.x);
  }
  inline int relative_y(int cell_y, int node_y){
    return relative(cell_y,node_y,have_cell.y,full_cell.y);
  }
  inline int relative_z(int cell_z, int node_z){
    return relative(cell_z,node_z,have_cell.z,full_cell.z);
  }
  inline int cellid_to_relative_x(int cellid, int node_x){
    return relative(cellid%full_cell.x,node_x,have_cell.x,full_cell.x);
  }
  inline int cellid_to_relative_y(int cellid, int node_y){
    return relative((cellid/full_cell.x)%full_cell.y,node_y,have_cell.y,full_cell.y);
  }
  inline int cellid_to_relative_z(int cellid, int node_z){
    return relative(cellid/(full_cell.x*full_cell.y),node_z,have_cell.z,full_cell.z);
  }
};

class GeometryXYZ : public Geometry {
public:

  GeometryXYZ(){}
  
  GeometryXYZ(SpaceVector<int> _size, SpaceVector<int> _full_cell, 
              bool _periodic=true)
    : Geometry(_size,_full_cell,_periodic)
  {
  }

  void path_depth(SpaceVector<int> &plus_depth,
                  SpaceVector<int> &minus_depth,
                  std::vector<int> target, int center_node);
  int getSendPath(int target_node, int center_node);
  int getReceivePath(int target_node, int center_node);
  void pushtarget(std::vector< std::vector<int> > SendTarget, 
                  std::vector< std::vector<int> > ReceiveTarget,
                  int offset,
                  int splus_depth, int sminus_depth,
                  int rplus_depth, int rminus_depth,
                  int plus, int minus);
  int getCommStage(std::vector< std::vector<int> > SendTarget, 
                   std::vector< std::vector<int> > ReceiveTarget,
                   std::vector<int> allSendTarget,
                   std::vector<int> allReceiveTarget,
                   int center_node);
};

#endif
