#include <iostream>
#include <algorithm>

#include "ParticleInfo.h"
#include "CovalentBondInfo.h"
#include "CellIndex.h"
#include "CubicCell.h"
#include "HalfShell.h"
#include "CellMethodImpl.h"

void
dump_cellindextype(CellIndexType cit)
{
  std::cout << "Cell Index Type ";
  if(cit==FullCube){
    std::cout << "Full Cube";
  }else if(cit==HalfShell){
    std::cout << "Half Shell";
  }else if(cit==FullShell){
    std::cout << "Full Shell";
  }else if(cit==SmallBall){
    std::cout << "Small Ball";
  }else{
    std::cout << "Unknown";
  }
  std::cout << std::endl;
}

//! check n is 2^p
bool
power_of_two(int n, int &p)
{
  int tmp=n;
  p = 0;
  while(tmp>1){
    tmp = tmp>>1;
    p++;
  }
  int p2 = 1<<p;

  if(p2!=n){
    return false;
  }
  return true;
}

//! divide p to sp[0] sp[1] sp[2], p = sp[0]+sp[1]+sp[2]
/* 
   p is calculated by power_of_two
   2^p = 2^sp[0] * 2^sp[1] * 2^sp[2]
   n = 2^p : number of unit, eg. node of MPI
   nx,ny,nz = 2^sp[0],2^sp[1],2^sp[2] : 3D unit grid 
   sp[0] <= sp[1] <= sp[2]
 */
void
split_power_of_two(int p, int sp[3])
{
  sp[0] = p/3;
  sp[1] = p/3;
  sp[2] = p/3;
  if((p%3)==2){
    sp[1]++;
    sp[2]++;
  }else if((p%3)==1){
    sp[2]++;
  }
}

bool
is_square(int np, int m[2])
{
  double p = sqrt(double(np));
  int ip = int(rint(p));
  
  if(ip*ip==np){
    m[0] = ip;
    m[1] = ip;
    return true;
  }
  return false;
}

bool
is_double_square(int np, int m[2])
{
  if(np%2==0){
    if(is_square(np/2,m)){
      m[1]*=2;
      return true;
    }
  }
  return false;
}

bool
power_of_two_of_square(int np, int m[2])
{
  if(!is_cube(np,m)){
    if(!is_double_cube(np,m)){
      return false;
    }
  }
  return true;
}

bool
is_cube(int np, int m[3])
{
  double p = pow(double(np),double(1.0/3.0));
  int ip = int(rint(p));
  
  if(ip*ip*ip==np){
    m[0] = ip;
    m[1] = ip;
    m[2] = ip;
    return true;
  }
  return false;
}

bool
is_double_cube(int np, int m[3])
{
  if(np%2==0){
    if(is_cube(np/2,m)){
      m[2]*=2;
      return true;
    }
  }
  return false;
}

bool
is_quad_cube(int np, int m[3])
{
  if(np%4==0){
    if(is_cube(np/4,m)){
      m[1]*=2;
      m[2]*=2;
      return true;
    }
  }
  return false;
}

bool
power_of_two_of_cube(int np, int m[3])
{
  if(!is_cube(np,m)){
    if(!is_double_cube(np,m)){
      if(!is_quad_cube(np,m)){
        return false;
      }
    }
  }
  return true;
}

void
CellIndex::set_cell_parameter(double bs, int cd, double mg, double ct,
                              double bond_len)
{
  cellsize.x = bs/double(cd);
  cellsize.y = cellsize.x;
  cellsize.z = cellsize.x;
  invcellsize.x = 1.0/cellsize.x;
  invcellsize.y = 1.0/cellsize.y;
  invcellsize.z = 1.0/cellsize.z;
  boxsize(bs,bs,bs);
  celldiv(cd,cd,cd);
  cutoff = ct;
  margin(mg,mg,mg);
  bond_length = bond_len;
  int tr = int(ceil((ct+mg*2.0)/cellsize.x));
  target_range(tr,tr,tr);
  int mr = int(mg/cellsize.x)+1;
  move_range(mr,mr,mr);
  int br = int(ceil((bond_len*2.0+mg*2.0)/cellsize.x));
  bond_range(br,br,br);
  calc_volume();
}
void
CellIndex::set_cell_parameter(SpaceVector<double> bs, 
                              SpaceVector<int> cd, SpaceVector<double> mg,
                              double ct,
                              double bond_len)
{
  cellsize.x = bs.x/double(cd.x);
  cellsize.y = bs.y/double(cd.y);
  cellsize.z = bs.z/double(cd.z);
  invcellsize.x = 1.0/cellsize.x;
  invcellsize.y = 1.0/cellsize.y;
  invcellsize.z = 1.0/cellsize.z;
  boxsize = bs;
  celldiv = cd;
  cutoff = ct;
  margin  = mg;
  bond_length = bond_len;
  target_range(int(ceil((ct+mg.x*2.0)/cellsize.x)),
               int(ceil((ct+mg.y*2.0)/cellsize.y)),
               int(ceil((ct+mg.z*2.0)/cellsize.z)) );
  move_range(int(mg.x/cellsize.x)+1,
             int(mg.y/cellsize.y)+1,
             int(mg.z/cellsize.z)+1 );
  bond_range(int(ceil((bond_len*2.0+mg.x*2.0)/cellsize.x)),
             int(ceil((bond_len*2.0+mg.y*2.0)/cellsize.y)),
             int(ceil((bond_len*2.0+mg.z*2.0)/cellsize.z)));
  calc_volume();
}

void
CellIndex::change_boxsize(SpaceVector<double> bs)
{
  cellsize.x = bs.x/double(celldiv.x);
  cellsize.y = bs.y/double(celldiv.y);
  cellsize.z = bs.z/double(celldiv.z);
  invcellsize.x = 1.0/cellsize.x;
  invcellsize.y = 1.0/cellsize.y;
  invcellsize.z = 1.0/cellsize.z;
  boxsize = bs;
  target_range(int(ceil((cutoff+margin.x*2.0)/cellsize.x)),
               int(ceil((cutoff+margin.y*2.0)/cellsize.y)),
               int(ceil((cutoff+margin.z*2.0)/cellsize.z)) );
  move_range(int(margin.x/cellsize.x)+1,
             int(margin.y/cellsize.y)+1,
             int(margin.z/cellsize.z)+1 );
  bond_range(int(ceil((bond_length*2.0+margin.x*2.0)/cellsize.x)),
             int(ceil((bond_length*2.0+margin.y*2.0)/cellsize.y)),
             int(ceil((bond_length*2.0+margin.z*2.0)/cellsize.z)));
  calc_volume();
}

double
CellIndex::calc_volume(){
  volume = boxsize.x*boxsize.y*boxsize.z;
  return volume;
}

//! set cubic cell position(cell index)
// order od cell_position : x y z
/* (0,0,0), (1,0,0), ... ,(celldiv-1,0,0), (0,1,0), ... ,(celldiv-1,2,0),
   (0,3,0), ... , (celldiv-1,celldiv-1,celldiv-2),
   (0,celldive-1,celldiv-1), ... ,(celldiv-1,celldiv-1,celldiv-1)
*/
size_t
CellIndex::generate_cell_index()
{
  cell_position.resize(celldiv.x*celldiv.y*celldiv.z);
  num_cell=0;
  for(int cz=0;cz<celldiv.z;cz++){
    for(int cy=0;cy<celldiv.y;cy++){
      for(int cx=0;cx<celldiv.x;cx++){
        cell_position[num_cell] = SpaceVector<int>(cx,cy,cz);
        num_cell++;
      }
    }
  }
  return num_cell;
}


//! divide cells to node
/* 
   @param[in] nodediv number of node each axis
   celldiv.{x|y|z} = nodediv.{x|y|z}*{l,m,n}
   cell_position : generated by generate_cell_index
   
   cell_id_in_node : index of cell_position for each node
   cell_id_to_node : reverse map index of cell_position to node

   index means array subscripts

   TODO
     error check celldiv is not sp[2]*integer
 */
bool
CellIndex::distribute_cell(const SpaceVector<int> nodediv)
{
  bool success = true;

  if((celldiv.x%nodediv.x>0)
     ||(celldiv.y%nodediv.y>0)
     ||(celldiv.z%nodediv.z>0)){
    std::cout << "celldiv != nodediv*(l,m,n) " << celldiv << nodediv << std::endl;
    success = false;
  }
  celldiv_node.x  = celldiv.x/nodediv.x;
  celldiv_node.y  = celldiv.y/nodediv.y;
  celldiv_node.z  = celldiv.z/nodediv.z;
  
  num_node = nodediv.x*nodediv.y*nodediv.z;
  cell_id_in_node.resize(num_node);
  cell_id_in_node.clear();
  for(int n=0;n<num_node;n++){
    std::vector<int> tmp(0);
    cell_id_in_node.push_back(tmp);
  }
  for(size_t i=0;i<cell_position.size();i++){
    int distx = cell_position[i].x/celldiv_node.x;
    int disty = cell_position[i].y/celldiv_node.y;
    int distz = cell_position[i].z/celldiv_node.z;
    int node = distx + (disty + distz*nodediv.y)*nodediv.x;
    cell_id_in_node[node].push_back(i);
    cell_id_to_node.insert(std::pair<int,int>(i,node));
  }
  node_size.x = boxsize.x/nodediv.x;
  node_size.y = boxsize.x/nodediv.y;
  node_size.z = boxsize.x/nodediv.z;
  node_volume = node_size.x*node_size.y*node_size.z;
  node_move_surface_volume = (margin.x*node_size.y*node_size.z
                              + margin.y*node_size.z*node_size.x
                              + margin.z*node_size.x*node_size.y)*2.0;
  node_move_surface_ratio = node_move_surface_volume/node_volume;
  return success;
}

void
CellIndex::dump_cell(int id)
{
  std::cout << " " << id << cell_position[id];
}

void
CellIndex::dump_cell_all()
{
  int num = cell_position.size();
  std::cout << num << " cell";
  for(int i=0;i<num;i++){
    dump_cell(i);
  }
  std::cout << std::endl;
}

void
CellIndex::dump_cell_subset(std::vector<int>& id)
{
  std::cout << id.size() << " cell";
  for(std::vector<int>::iterator i = id.begin(); i!=id.end();++i){
    dump_cell(*i);
  }
  std::cout << std::endl;
}






//! use in merge_shift_cell
inline
void addshift(std::vector<int>& set, SpaceVector<int> shift,
              std::map< int,SpaceVector<int> >& set_to_shift)
{
  for(std::vector<int>::iterator s=set.begin(); s!=set.end(); ++s){
    std::map< int,SpaceVector<int> >::iterator se = set_to_shift.find((*s));
    if(se!=set_to_shift.end()){
      se->second += shift;
    }else{
      set_to_shift.insert(std::pair<int,SpaceVector<int> >((*s),shift));
    }
  }
}

CellSubset::CellSubset(const CellIndex& ci, const int& mn)
  : cellindex(ci), mynode(mn)
{
}

CellSubset::CellSubset(const CellSubset* cs)
  : cellindex(cs->cellindex), mynode(cs->mynode)
{
} 

//! merge all shift_cell calculated by cubic_target()
/* 
   When number of node is small, 
      different ghosts of a cell will appear in one node.
   This version only warn such different ghosts and ignore.
*/
void
CellSubset::merge_shift_cell()
{
  std::vector< std::map< int,SpaceVector<int> > > set_to_shift;

  
  set_to_shift.resize(shift_cell.size());
  for(size_t ci=0;ci<shift_cell.size();ci++){
    addshift(shift_cell[ci].plusx, SpaceVector<int>(1,0,0),set_to_shift[ci]);
    addshift(shift_cell[ci].minusx, SpaceVector<int>(-1,0,0),set_to_shift[ci]);
    addshift(shift_cell[ci].plusy, SpaceVector<int>(0,1,0),set_to_shift[ci]);
    addshift(shift_cell[ci].minusy, SpaceVector<int>(0,-1,0),set_to_shift[ci]);
    addshift(shift_cell[ci].plusz, SpaceVector<int>(0,0,1),set_to_shift[ci]);
    addshift(shift_cell[ci].minusz, SpaceVector<int>(0,0,-1),set_to_shift[ci]);
  }
  // this version not support different shift for a cell, but only print out
  {
    int differ = 0;
    for(size_t ci=0;ci<shift_cell.size()-1;ci++){
      for(size_t ci2=ci+1;ci2<shift_cell.size();ci2++){
        for(std::map< int,SpaceVector<int> >::iterator s=set_to_shift[ci].begin();
            s!=set_to_shift[ci].end();++s){
          std::map< int,SpaceVector<int> >::iterator se = set_to_shift[ci2].find((s->first));
          if(se!=set_to_shift[ci2].end()){
    if(se->second!=s->second){
              std::cout << " " << s->first << ":" << ci << s->second << ":" << ci2 << se->second;
              differ++;
            }
          }
        }
      }
    }
    if(differ>0){
      std::cout << " found different shift cell" ;
      std::cout << std::endl;
    }
  }
  
  std::set<int> tmpset;
  for(size_t ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].plusx.begin(),shift_cell[ci].plusx.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.plusx.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.plusx.begin());
  
  tmpset.clear();
  for(size_t ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].minusx.begin(),shift_cell[ci].minusx.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.minusx.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.minusx.begin());

  tmpset.clear();
  for(size_t ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].plusy.begin(),shift_cell[ci].plusy.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.plusy.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.plusy.begin());

  tmpset.clear();
  for(size_t ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].minusy.begin(),shift_cell[ci].minusy.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.minusy.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.minusy.begin());

  tmpset.clear();
  for(size_t ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].plusz.begin(),shift_cell[ci].plusz.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.plusz.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.plusz.begin());

  tmpset.clear();
  for(size_t ci=0;ci<shift_cell.size();ci++){
    std::copy(shift_cell[ci].minusz.begin(),shift_cell[ci].minusz.end(),std::inserter(tmpset,tmpset.end()));
  }
  shift_cell_all.minusz.resize(tmpset.size());
  std::copy(tmpset.begin(),tmpset.end(),shift_cell_all.minusz.begin());
}

//! selecet self target
int
CellSubset::calc_self_target()
{
  int num_total_self_target=0;
  self_target_cell.resize(mycell.size());
  for(size_t ci=0;ci<mycell.size();ci++){
    self_target_cell[ci].clear();
    for(size_t cj=0;cj<recv_target_cell[ci].size();cj++){
      int tcell = recv_target_cell[ci][cj];
      std::vector<int>::iterator t = std::find(mycell.begin(),mycell.end(),tcell);
      if(t!=mycell.end()){
        self_target_cell[ci].push_back(tcell);
      }
    }
    num_total_self_target+=self_target_cell[ci].size();
  }
  return num_total_self_target;
}

//! select target nodes which has target cell, and merge target nodes for every cell in this node
/*
  target_node : merged target nodes, this node must receive from these nodes
 */
int
CellSubset::calc_target_gather(std::vector<int>& recv_target_node,
                               std::vector< std::vector<int> >& recv_target_set)
{
  std::map<int,int> target_to_index;   // map target node id to index for td, because td will be not sorted by target node id.
  std::vector< std::set<int> > td;     // target cell in targe node, node is not sorted
  std::set<int> tmptarget;             // target nodes, uniq and sorted
  int num_target=0;
  for(size_t i=0;i<recv_target_cell.size();i++){
    for(size_t j=0;j<recv_target_cell[i].size();j++){
      int tcell = recv_target_cell[i][j];
      int n = cellindex.cell_id_to_node[tcell];
      if(n!=mynode){
        std::map<int,int>::iterator ti = target_to_index.find(n);
        if(ti==target_to_index.end()){
          tmptarget.insert(n);
          target_to_index.insert(std::pair<int,int>(n,num_target));
          num_target++;
          std::set<int> newdist;
          newdist.insert(tcell);
          td.push_back(newdist);
        }else{
          td[(*ti).second].insert(tcell);
        }
      }
    }
  }
  /* 
     convert td, target cell in targe node which unsorted by node id,
     to recv_target_dist, target cell in targe node which sorted by node id,
   */ 
  for(std::set<int>::iterator t=tmptarget.begin();t!=tmptarget.end();++t){
    int target=(*t);
    int ti = target_to_index[target];
    recv_target_node.push_back(target);
    std::vector<int> dist;
    for(std::set<int>::iterator cell=td[ti].begin();cell!=td[ti].end();++cell){
      dist.push_back((*cell));
    }
    recv_target_set.push_back(dist);
  }
  return recv_target_node.size();
}

//! select target nodes which require cells in this node, and merge target nodes for every cell in this node
/*
  target_node : merged target node, this node must send to these nodes
 */
int
CellSubset::calc_target_scatter(std::vector<int>& send_target_node,
                                std::vector< std::vector<int> >& send_target_set)
{
  std::map<int,int> target_to_index;   // map target node id to index for td, because td will be not sorted by target node id.
  std::vector< std::set<int> > td;     // target cell in targe node, node is not sorted
  std::set<int> tmptarget;             // target nodes, uniq and sorted
  int num_target=0;
  for(size_t i=0;i<send_target_cell.size();i++){
    int icell = mycell[i];
    for(size_t j=0;j<send_target_cell[i].size();j++){
      int tcell = send_target_cell[i][j];
      int n = cellindex.cell_id_to_node[tcell];
      if(n!=mynode){
        std::map<int,int>::iterator ti = target_to_index.find(n);
        if(ti==target_to_index.end()){
          tmptarget.insert(n);
          target_to_index.insert(std::pair<int,int>(n,num_target));
          num_target++;
          std::set<int> newdist;
          newdist.insert(icell);
          td.push_back(newdist);
        }else{
          td[(*ti).second].insert(icell);
        }
      }
    }
  }
  /* 
     convert td, target cell in targe node which unsorted by node id,
     to send_target_dist, target cell in targe node which sorted by node id,
   */ 
  for(std::set<int>::iterator t=tmptarget.begin();t!=tmptarget.end();++t){
    int target=(*t);
    int ti = target_to_index[target];
    send_target_node.push_back(target);
    std::vector<int> dist;
    for(std::set<int>::iterator cell=td[ti].begin();cell!=td[ti].end();++cell){
      dist.push_back((*cell));
    }
    send_target_set.push_back(dist);
  }
  return send_target_node.size();
}

int
calc_gather_unsorted(const std::vector<int>& mycell,
		     const std::vector< std::vector<int> >& target_cell,
		     const std::map<int,int>& cell_id_to_node,
		     const int mynode,
		     std::vector<int>& target_node,
		     std::vector< std::vector<int> >& target_set)
{
  std::vector<int> target_to_index;   //  target node id to index for td, because td will be not sorted by target node id.
  std::vector< std::set<int> > td;     // target cell in targe node, node is not sorted
  int num_target=0;
  for(size_t i=0;i<target_cell.size();i++){
    for(size_t j=0;j<target_cell[i].size();j++){
      int tcell = target_cell[i][j];
      int n = cell_id_to_node.find(tcell)->second;
      if(n!=mynode){
        std::vector<int>::iterator ti = std::find(target_to_index.begin(),target_to_index.end(),n);
        if(ti==target_to_index.end()){
          target_to_index.push_back(n);
          num_target++;
          std::set<int> newdist;
          newdist.insert(tcell);
          td.push_back(newdist);
        }else{
          td[ti-target_to_index.begin()].insert(tcell);
        }
      }
    }
  }
  /* 
     convert td, target cell in targe node which unsorted by node id,
     to send_target_dist, target cell in targe node which sorted by node id,
   */ 
  for(int ttoi=0; ttoi<target_to_index.size();++ttoi){
    int target=target_to_index[ttoi];
    target_node.push_back(target);
    std::vector<int> dist;
    for(std::set<int>::iterator cell=td[ttoi].begin();cell!=td[ttoi].end();++cell){
      dist.push_back((*cell));
    }
    target_set.push_back(dist);
  }
  return target_node.size();
}

int
calc_scatter_unsorted(const std::vector<int>& mycell,
		      const std::vector< std::vector<int> >& target_cell,
		      const std::map<int,int>& cell_id_to_node,
		      const int mynode,
		      std::vector<int>& target_node,
		      std::vector< std::vector<int> >& target_set)
{
  std::vector<int> target_to_index;   //  target node id to index for td, because td will be not sorted by target node id.
  std::vector< std::set<int> > td;     // target cell in targe node, node is not sorted
  int num_target=0;
  for(size_t i=0;i<target_cell.size();i++){
    int icell = mycell[i];
    for(size_t j=0;j<target_cell[i].size();j++){
      int tcell = target_cell[i][j];
      int n = cell_id_to_node.find(tcell)->second;
      if(n!=mynode){
        std::vector<int>::iterator ti = std::find(target_to_index.begin(),target_to_index.end(),n);
        if(ti==target_to_index.end()){
          target_to_index.push_back(n);
          num_target++;
          std::set<int> newdist;
          newdist.insert(icell);
          td.push_back(newdist);
        }else{
          td[ti-target_to_index.begin()].insert(icell);
        }
      }
    }
  }
  /* 
     convert td, target cell in targe node which unsorted by node id,
     to send_target_dist, target cell in targe node which sorted by node id,
   */ 
  for(int ttoi=0; ttoi<target_to_index.size();++ttoi){
    int target=target_to_index[ttoi];
    target_node.push_back(target);
    std::vector<int> dist;
    for(std::set<int>::iterator cell=td[ttoi].begin();cell!=td[ttoi].end();++cell){
      dist.push_back((*cell));
    }
    target_set.push_back(dist);
  }
  return target_node.size();
}

//! calculate cell that is move target and target node that have target cell
int
CellSubset::calc_target_move(std::vector<int>& move_target_node,
                             std::vector< std::vector<int> >& move_target_set)
{
  if(DebugLog::verbose>1){
    std::cout << "calc_target_move ";
    std::cout << "celldiv " << cellindex.celldiv << " move range " << cellindex.move_range << std::endl;
  }

  std::set<int> tids;  //! ID of cells that is move target

  // find move target cells
  for(size_t s=0;s<mycell.size();s++){
    SpaceVector<int> pos = cellid_to_cell_position(mycell[s],cellindex.celldiv);
    //    std::cout << "search " << pos << std::endl;
    for(int z=pos.z-cellindex.move_range.z;z<=pos.z+cellindex.move_range.z;z++){
      int sz = z;
      periodic_shift(sz,cellindex.celldiv.z);
      for(int y=pos.y-cellindex.move_range.y;y<=pos.y+cellindex.move_range.y;y++){
        int sy = y;
        periodic_shift(sy,cellindex.celldiv.y);
        for(int x=pos.x-cellindex.move_range.x;x<=pos.x+cellindex.move_range.x;x++){
          int sx = x;
          periodic_shift(sx,cellindex.celldiv.x);
          int tid = cell_position_to_cellid(SpaceVector<int>(sx,sy,sz),cellindex.celldiv);
          std::vector<int>::iterator t = std::find(mycell.begin(),mycell.end(),tid);
          if(t==mycell.end()){
            //      std::cout << " out cell " << tid << SpaceVector<int>(sx,sy,sz) << " for " << mycell[s] << std::endl;
            tids.insert(tid);
          }
        }
      }
    }
  }

  // find target node that have target cell, and distribute target cell
  for(std::set<int>::iterator tid = tids.begin(); tid!=tids.end(); ++tid){
    int node = cellindex.cell_id_to_node[*tid];
    //    std::cout << " cell " << *tid << " to node " << node << std::endl;
    size_t tn;
    for(tn=0;tn<move_target_node.size();tn++){
      if(move_target_node[tn]==node)break;
    }
    if(tn==move_target_node.size()){  // new target
      move_target_node.push_back(node);
      std::vector<int> tsetid(1,(*tid));
      move_target_set.push_back(tsetid);
    }else{                      // exist target
      move_target_set[tn].push_back((*tid));
    }
  }
  return move_target_node.size();
}


inline
void print_shift_cell(std::vector<int>& cells)
{
  std::cout << cells.size() << " cells";
  for(std::vector<int>::size_type i=0;i<cells.size();i++){
    std::cout << " " << cells[i];
  }
  std::cout << std::endl;
}

void
CellSubset::makefulltarget(std::vector<int> subset, 
                           CellIndexType cellindextype,
                           OperationSelector operations)
{
  mycell = subset;
  num_mycell = mycell.size();
  send_target_cell.resize(num_mycell);
  recv_target_cell.resize(num_mycell);
  short_target_cell.resize(num_mycell);
  shift_cell.resize(num_mycell);
  ghost_cell_pairs.clear();
  int ghost_len = std::max(std::max(cellindex.target_range.x,cellindex.target_range.y),cellindex.target_range.z);
  cellmethod.setCellIndexType(cellindextype);
  for(int i=0;i<num_mycell;i++){
    send_target_cell[i].clear();
    recv_target_cell[i].clear();
    short_target_cell[i].clear();
    shift_cell[i].clear();
    if(cellindextype==FullCube){
      /*
      cubic_target(cellindex.celldiv.x,cellindex.cell_position[mycell[i]],
                   ghost_len,
                   send_target_cell[i],recv_target_cell[i], shift_cell[i]);
      */
      cubic_target(cellindex.celldiv,cellindex.cell_position[mycell[i]],
                   ghost_len,
                   send_target_cell[i],recv_target_cell[i], shift_cell[i]);
      short_target_cell[i] = recv_target_cell[i];
    }else{
#if 0
      halfshell_target(cellindex.celldiv.x, 
                       cellindex.cell_position[mycell[i]], 
                       cellindex.cellsize, cellindex.margin, cellindex.cutoff,
                       send_target_cell[i],
                       recv_target_cell[i], shift_cell[i]);
      short_target_cell[i] = recv_target_cell[i];
#else
      cellmethod.target(cellindex.celldiv.x, 
                        cellindex.cell_position[mycell[i]], 
                        cellindex.cellsize, cellindex.margin, cellindex.cutoff,
                        send_target_cell[i], recv_target_cell[i],
                        short_target_cell[i], shift_cell[i],
                        ghost_cell_pairs);
#endif
    }
    
    if(operations.doCovalentBondcaculation){
      int bond_reach = std::max(std::max(cellindex.bond_range.x,cellindex.bond_range.y),cellindex.bond_range.z);
      if(DebugLog::verbose>1){
        std::cout << "bond reach " << bond_reach << " cell " << std::endl;
      }
//       cubic_target(cellindex.celldiv.x,cellindex.cell_position[mycell[i]],
//                    bond_reach,
//                    send_target_cell[i],recv_target_cell[i], shift_cell[i]);
      cubic_target(cellindex.celldiv,cellindex.cell_position[mycell[i]],
                   bond_reach,
                   send_target_cell[i],recv_target_cell[i], shift_cell[i]);
    }
  }
  merge_shift_cell();
  if(DebugLog::verbose>1){
    int send_num=0, recv_num=0, short_num=0, pair_num=0;
    for(int i=0;i<num_mycell;i++){
      send_num += send_target_cell[i].size();
      recv_num += recv_target_cell[i].size();
      short_num += short_target_cell[i].size();
    }
    pair_num += ghost_cell_pairs.size();
    std::cout << "num target cell send,recv,short,pair ";
    std::cout << " " << send_num;
    std::cout << " " << recv_num;
    std::cout << " " << short_num;
    std::cout << " " << pair_num;
    std::cout << std::endl;
  }
  if(DebugLog::verbose>1)
  {
    std::cout << "shift_cell_all.plusx ";
    print_shift_cell(shift_cell_all.plusx);
    std::cout << "shift_cell_all.minusx  ";
    print_shift_cell(shift_cell_all.minusx );
    std::cout << "shift_cell_all.plusy ";
    print_shift_cell(shift_cell_all.plusy);
    std::cout << "shift_cell_all.minusy  ";
    print_shift_cell(shift_cell_all.minusy );
    std::cout << "shift_cell_all.plusz ";
    print_shift_cell(shift_cell_all.plusz);
    std::cout << "shift_cell_all.minusz  ";
    print_shift_cell(shift_cell_all.minusz );
  }
}

bool
CellSubset::makecommtarget(std::vector<int>& send_target_node,
                           std::vector< std::vector<int> >& send_target_set,
                           std::vector<int>& recv_target_node,
                           std::vector< std::vector<int> >& recv_target_set,
                           std::vector<int>& move_target_node,
                           std::vector< std::vector<int> >& move_target_set)
{

  size_t num_send_target_node
   //    = calc_target_scatter(send_target_node, send_target_set);
    = calc_scatter_unsorted(mycell, send_target_cell, cellindex.cell_id_to_node, mynode,
			   send_target_node, send_target_set);
  size_t num_recv_target_node 
  //    = calc_target_gather(recv_target_node, recv_target_set);
    = calc_gather_unsorted(mycell, recv_target_cell, cellindex.cell_id_to_node, mynode,
			   recv_target_node, recv_target_set);
  //DEBUG code
  /*
  if(mynode==0){
    size_t nt;
    std::vector<int> rtn;
    std::vector< std::vector<int> > rts;
    nt = calc_target_scatter(rtn, rts);
    for(int t=0;t<nt;t++){
      printf("send target sort %d : ",rtn[t]);
      for(int tc=0;tc<rts[t].size();tc++){
	printf(" %d",rts[t][tc]);
      }
      printf("\n");
    }
    for(int t=0;t<num_send_target_node;t++){
      printf("send target %d : ",send_target_node[t]);
      for(int tc=0;tc<send_target_set[t].size();tc++){
	printf(" %d",send_target_set[t][tc]);
      }
      printf("\n");
    }
    rtn.clear();
    rts.clear();
    nt = calc_target_gather(rtn, rts);
    for(int t=0;t<nt;t++){
      printf("recv target sort %d : ",rtn[t]);
      for(int tc=0;tc<rts[t].size();tc++){
	printf(" %d",rts[t][tc]);
      }
      printf("\n");
    }
    for(int t=0;t<num_recv_target_node;t++){
      printf("recv target %d : ",recv_target_node[t]);
      for(int tc=0;tc<recv_target_set[t].size();tc++){
	printf(" %d",recv_target_set[t][tc]);
      }
      printf("\n");
    }
  }
  */
  if(DebugLog::verbose>1){
    std::cout << "num_send_target_node " << num_send_target_node << std::endl;
    std::cout << "num_recv_target_node " << num_recv_target_node << std::endl;
  }

  calc_target_move(move_target_node, move_target_set);
   

  return true;
}

void
CellSubset::dump_cell_target(int cell_index, std::vector<int>& target_cell)
{
  std::cout << "target for cell " << cell_index;
  cellindex.dump_cell(mycell[cell_index]);
  std::cout << " : ";
  cellindex.dump_cell_subset(target_cell);
  std::cout << std::endl;
}

void
CellSubset::dump_cell_target_send(int cell_index)
{
  std::cout << "send ";
  dump_cell_target(cell_index,send_target_cell[cell_index]);
}

void
CellSubset::dump_cell_target_recv(int cell_index)
{
  std::cout << "recv ";
  dump_cell_target(cell_index,recv_target_cell[cell_index]);
}

void
CellSubset::dump_cell_target_self(int cell_index)
{
  std::cout << "self ";
  dump_cell_target(cell_index,self_target_cell[cell_index]);
}


template<class PA, class SPA>
int
CellSubset::distribute_particle_cell(PA& particlearray,
                                     std::vector<int>& particle_setid,
                                     std::vector<TypeRange>& typerangearray,
                                     WaterList& waterlist,
                                     ShakeList& shakelist,
                                     const SPA& source,
                                     const std::vector<PotentialModel>& pmodel,
                                     const WaterList& wl,
                                     const ShakeList& sl)
{
  int ncell = typerangearray.size();
  int n=0;
  SpaceVector<double> invcs(1.0/cellindex.cellsize.x,1.0/cellindex.cellsize.y,1.0/cellindex.cellsize.z);
  for(size_t i=0;i<source.size();i++){
#ifdef USE_SHAKE
    if(getatomtype(source,i) <= ATOMTYPE_WH)  // omit WH (of water) & H1 (of shake bond)
#else
    if(getatomtype(source,i) == ATOMTYPE_WH)  // omit WH (of water)
#endif
      continue;
    else{
      Particle pa = getparticle(source,i);
      Position& pos = pa.position;
      int shiftx = periodic_shift(pos.x,cellindex.boxsize.x);
      int shifty = periodic_shift(pos.y,cellindex.boxsize.y);
      int shiftz = periodic_shift(pos.z,cellindex.boxsize.z);
      SpaceVector<int> cell(int(pos.x*invcs.x),
                            int(pos.y*invcs.y),
                            int(pos.z*invcs.z));
      int ci = cell.x+ (cell.y+cell.z*cellindex.celldiv.y)*cellindex.celldiv.x;
      size_t si;
      for(si=0;si<particle_setid.size();si++){
        if(particle_setid[si]==ci)break;
      }
      if(getatomtype(source,i) == ATOMTYPE_WO){
        Particle wh1 = getparticle(source, wl.find(i)->second.h1 );
        Particle wh2 = getparticle(source, wl.find(i)->second.h2 );
        wh1.position.x += shiftx * cellindex.boxsize.x;
        wh1.position.y += shifty * cellindex.boxsize.y;
        wh1.position.z += shiftz * cellindex.boxsize.z;
        wh2.position.x += shiftx * cellindex.boxsize.x;
        wh2.position.y += shifty * cellindex.boxsize.y;
        wh2.position.z += shiftz * cellindex.boxsize.z;
        if(si<particle_setid.size()){
          add_water(particlearray, typerangearray[si], waterlist, pa, wh1, wh2);
          n+=3;
        }
      }else{ /// not Oxygen of Water 
#ifdef USE_SHAKE
        if(sl.find(i) != sl.end()){ // shake bond, not water
          int nh1 = sl.find(i)->second.nh1;
          for(int n1=0; n1<nh1; n1++){
            // H atoms of shake bond
            int h1_atomid = sl.find(i)->second.h1[n1];
            Particle h1 = getparticle(source, h1_atomid);
            h1.position.x += shiftx * cellindex.boxsize.x;
            h1.position.y += shifty * cellindex.boxsize.y;
            h1.position.z += shiftz * cellindex.boxsize.z;
            if(si<particle_setid.size()){
              // add two atoms(LJCoulombPotential) to shakelist & particlearray by every (HA,H1) unit
              add_shake(particlearray, typerangearray[si], waterlist, shakelist, n1, pa, h1, sl.find(i)->second.bondtype[n1]);
              if(n1==0){
                n+=2;
              }else{
                n++;
              }
            }
          }
        } else
#endif
	  { // not shake bond, not water
	    if(si<particle_setid.size()){
	      add_particle(particlearray, typerangearray[si], waterlist, shakelist, pa, pmodel[i]);
	      n++;
	    }
	  }
      }
      if(si<ncell-1){
	if(typerangearray[si].end==typerangearray[si+1].begin){
	  printf("cell %d reach max %d\n",(int)si,typerangearray[si].end);
	}
	if(typerangearray[si].end>typerangearray[si+1].begin){
	  printf("cell %d overflow %d\n",(int)si,typerangearray[si].end);
	}
      }
    }// atomtype != WH
  }
  return n;
}

template
int
CellSubset::distribute_particle_cell(ParticleArray& particlearray,
                                     std::vector<int>& particle_setid,
                                     std::vector<TypeRange>& typerangearray,
                                     WaterList& waterlist,
                                     ShakeList& shakelist,
                                     const ParticleArray& source,
                                     const std::vector<PotentialModel>& pmodel,
                                     const WaterList& wl,
                                     const ShakeList& sl);
template
int
CellSubset::distribute_particle_cell(CombinedParticleArray& particlearray,
                                     std::vector<int>& particle_setid,
                                     std::vector<TypeRange>& typerangearray,
                                     WaterList& waterlist,
                                     ShakeList& shakelist,
                                     const ParticleArray& source,
                                     const std::vector<PotentialModel>& pmodel,
                                     const WaterList& wl,
                                     const ShakeList& sl);
template
int
CellSubset::distribute_particle_cell(CombinedParticleArray& particlearray,
                                     std::vector<int>& particle_setid,
                                     std::vector<TypeRange>& typerangearray,
                                     WaterList& waterlist,
                                     ShakeList& shakelist,
                                     const CombinedParticleArray& source,
                                     const std::vector<PotentialModel>& pmodel,
                                     const WaterList& wl,
                                     const ShakeList& sl);


PostProcess::PostProcess()
  : delete_index(), delete_from_cellindex(), delete_rt_cache(),
    move_rt_cache(), move_to_cellid(),
    boxsize(), shift_cell()
{
  reserve(1);
  max_particle_in_cell = 0;
}

PostProcess::PostProcess(size_t psize)
  : delete_index(), delete_from_cellindex(), delete_rt_cache(),
    move_rt_cache(), move_to_cellid(),
    boxsize(), shift_cell() 
{
  reserve(psize);
}

PostProcess::PostProcess(size_t psize, 
                         SpaceVector<double> bs, ShiftCellArray& sc) 
  : delete_index(), delete_from_cellindex(), delete_rt_cache(),
    move_rt_cache(), move_to_cellid(),
    boxsize(bs), shift_cell(sc)
{
  reserve(psize);
}

void
PostProcess::reserve(size_t psize)
{
  max_particle_in_cell = psize;
  //  delete_index.reserve(psize);
  //  delete_from_cellindex.reserve(psize);
  //  delete_rt_cache.reserve(psize);

  move_cache.reserve(psize);
  move_to_cellid.reserve(psize);
  move_rt_cache.reserve(psize);
}

double
PostProcess::calc_volume_margin()
{
  SpaceVector<double> margincell = cellsize + margin;
  volume_margin = margincell.x*margincell.y*margincell.z/(cellsize.x*cellsize.y*cellsize.z);
  return volume_margin;
}

void
PostProcess::change_boxsize(SpaceVector<double> bs)
{
  boxsize = bs;
  cellsize.x = boxsize.x/celldiv.x;
  cellsize.y = boxsize.y/celldiv.y;
  cellsize.z = boxsize.z/celldiv.z;
  invcellsize.x = 1.0/cellsize.x;
  invcellsize.y = 1.0/cellsize.y;
  invcellsize.z = 1.0/cellsize.z;
}

//! listing particles move out other cell
/*!
  member store results
  number_move : number of move particles
  move_index : index of move particles
  move_cache : move particles
  rt_cache : type of move particles
  move_from_cellindex : index of set(cell) where particles move from
  move_to_cellid : id of set(cell) where particles move to

  scan each set(cell) typerange[*].begin to typerange[*].end
  move_index is grouped by set(cell) and sorted in ascending order inside each set(cell)
 */
template <typename PA>
bool
PostProcess::select_move_particle(PA& particle,
                                  std::vector<TypeRange>& typerange,
                                  std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                  WaterList& waterlist,
                                  ShakeList& shakelist,
                                  const std::vector<int>& cellid,
                                  int& over_move)
{
  bool move = false;
  over_move = 0;

  number_move = 0;
  number_move_water = 0;
  number_move_shake = 0;

  delete_index.clear();
  delete_from_cellindex.clear();
  delete_rt_cache.clear();
  move_cache.clear();
  move_rt_cache.clear();
  move_to_cellid.clear();
  move_bond_cache.clear();

  delete_wh_list.clear();
  move_cache_water.clear();
  move_to_cellid_water.clear();
  move_bond_cache_water.clear();

  move_cache_shake.clear();
  move_to_cellid_shake.clear();
  move_bond_cache_shake.clear();

//   move_index.clear();
//   rt_cache.clear();
//   move_cache.clear();
//   move_from_cellindex.clear();
//   move_to_cellid.clear();

  for(size_t ci=0;ci<cellid.size();ci++){
    int cid = cellid[ci];
    SpaceVector<double> cellboundmin, cellboundmax;
    /* for check interval */
    SpaceVector<double> cellboundmin_full, cellboundmax_full;
    cellboundmin.x = cellsize.x*cell_position[cid].x;
    cellboundmin.y = cellsize.y*cell_position[cid].y;
    cellboundmin.z = cellsize.z*cell_position[cid].z;
    /* for check interval */
    cellboundmax_full = cellboundmin + cellsize + margin;
    cellboundmin_full = cellboundmin - margin;
    cellboundmax = cellboundmin + cellsize + fine_margin;
    cellboundmin -= fine_margin;
    for(int i=typerange[ci].begin;i<typerange[ci].end;i++){
#ifdef USE_SHAKE
      if(getatomtype(particle,i) <= ATOMTYPE_WH)   // omit WH (of water) & H1 (of shake bond){
#else
	if(getatomtype(particle,i) == ATOMTYPE_WH)  // omit WH (of water)
#endif
	{
	continue;
	/*
        if( std::find(delete_wh_list.begin(),delete_wh_list.end(),i) != delete_wh_list.end() ){
//          printf("select particle (WH)= %d[%d]\n",i,particle[i].atomid);
          delete_index.insert(i);
          delete_from_cellindex.insert(std::pair<int,int>(i,ci));
          delete_rt_cache.insert(std::pair<int,PotentialModel>(i,typerange[ci].particle_index_to_type(i)));
        }
	*/
      }else{
        Position& pos = getpos(particle,i);
        if((pos.x<cellboundmin.x)||(pos.x>cellboundmax.x)
           || (pos.y<cellboundmin.y)||(pos.y>cellboundmax.y)
           || (pos.z<cellboundmin.z)||(pos.z>cellboundmax.z)){
          /* for check interval */
          if((pos.x<cellboundmin_full.x)||(pos.x>cellboundmax_full.x)
             || (pos.y<cellboundmin_full.y)||(pos.y>cellboundmax_full.y)
             || (pos.z<cellboundmin_full.z)||(pos.z>cellboundmax_full.z)){
            std::cout << "out of bouns original margin, must make small fine_margin" << std::endl;
            std::cout << "x:" << cellboundmin_full.x << "<" << pos.x << "<" << cellboundmax_full.x << std::endl;
            std::cout << "y:" << cellboundmin_full.y << "<" << pos.y << "<" << cellboundmax_full.y << std::endl;
            std::cout << "z:" << cellboundmin_full.z << "<" << pos.z << "<" << cellboundmax_full.z << std::endl;
            over_move = 1;
          }
          int shiftx = periodic_shift(pos.x,boxsize.x);
          int shifty = periodic_shift(pos.y,boxsize.y);
          int shiftz = periodic_shift(pos.z,boxsize.z);
          if(out_of_bounds(pos,boxsize)){
            std::cout << "out of bouns after periodic_shift " << i << pos << boxsize << std::endl;
            std::cout << "position" << getpos(particle,i);
            std::cout << "velocity" << getvelocity(particle,i);
            std::cout << "force" << getforce(particle,i);
            std::cout << std::endl;
          }
          SpaceVector<int> ncell(int(pos.x*invcellsize.x),
                                 int(pos.y*invcellsize.y),
                                 int(pos.z*invcellsize.z));
          int ncid = cell_position_to_cellid(ncell,celldiv);
          move = true;
          delete_index.insert(i);
          delete_from_cellindex.insert(std::pair<int,int>(i,ci));
          delete_rt_cache.insert(std::pair<int,PotentialModel>(i,typerange[ci].particle_index_to_type(i)));

          if(getatomtype(particle,i) != ATOMTYPE_WO){
#ifdef USE_SHAKE
            if(shakelist.find(i) != shakelist.end()){ // heavy atom of shake
              int nh1 = 0;
              nh1 = shakelist.find(i)->second.nh1;

              for(int inh1=0; inh1<nh1; inh1++){
                int sh1 = shakelist.find(i)->second.h1[inh1];

                delete_index.insert(sh1);                                                     
                delete_from_cellindex.insert(std::pair<int,int>(sh1,ci));
                delete_rt_cache.insert(std::pair<int,PotentialModel>(sh1,typerange[ci].particle_index_to_type(sh1)));

		//                particle[sh1].position.x += shiftx * boxsize.x;
		//                particle[sh1].position.y += shifty * boxsize.y;
		//                particle[sh1].position.z += shiftz * boxsize.z;
		getpos(particle,sh1).x += shiftx * boxsize.x;
		getpos(particle,sh1).y += shifty * boxsize.y;
		getpos(particle,sh1).z += shiftz * boxsize.z;

                move_cache_shake.push_back(getparticle(particle,i));
                move_cache_shake.push_back(getparticle(particle,sh1));
                move_to_cellid_shake.push_back(ncid);
                move_to_cellid_shake.push_back(ncid);

                number_move_shake+=2;                                                 // every (HA,H1) unit
                CovalentBondInfo::BondList bond_cache2;
                bond_cache2.clear();
                bondlistarray[ci].pick_up_bondlist(bond_cache2,getatomid(particle,i));   // pick up all bonds
                bondlistarray[ci].pick_up_bondlist(bond_cache2,getatomid(particle,sh1)); //  ceoncerned about HA,H1
                bondlistarray[ci].remove_atomid(getatomid(particle,i));                  // delete all bonds
                bondlistarray[ci].remove_atomid(getatomid(particle,sh1));                //  concerned about HA,H1
                move_bond_cache_shake.push_back(bond_cache2);
              }
              shakelist.delete_shake(i);
            }else
#endif
	      { // not shake atom
	      move_rt_cache.push_back(typerange[ci].particle_index_to_type(i));
              move_cache.push_back(getparticle(particle,i));
              move_to_cellid.push_back(ncid);
              CovalentBondInfo::BondList bond_cache;
              bond_cache.clear();
              //            dump_bondlistarray(bondlistarray);
              bondlistarray[ci].pick_up_bondlist(bond_cache,getatomid(particle,i));
              move_bond_cache.push_back(bond_cache);
              /*
              std::cout << " move_bond_cache ";
              dump_bondlistarray(move_bond_cache);
              dump_bondlistarray(bondlistarray);
              std::cout << "remve " << particle[i].atomid << " from " << ci << std::endl;
              */
              bondlistarray[ci].remove_atomid(getatomid(particle,i));
              //            dump_bondlistarray(bondlistarray);
              number_move++;
	      }
          }else{ // waterO
            //            std::cout << "found Water O " <<particle[i].atomid << std::endl;
            int wh1 = waterlist.find(i)->second.h1;
            int wh2 = waterlist.find(i)->second.h2;

            delete_index.insert(wh1);
            delete_from_cellindex.insert(std::pair<int,int>(wh1,ci));
            delete_rt_cache.insert(std::pair<int,PotentialModel>(wh1,typerange[ci].particle_index_to_type(wh1)));
            delete_index.insert(wh2);
            delete_from_cellindex.insert(std::pair<int,int>(wh2,ci));
            delete_rt_cache.insert(std::pair<int,PotentialModel>(wh2,typerange[ci].particle_index_to_type(wh2)));

            waterlist.delete_water(i);

            getpos(particle,wh1).x += shiftx * boxsize.x;
            getpos(particle,wh1).y += shifty * boxsize.y;
            getpos(particle,wh1).z += shiftz * boxsize.z;
            
            getpos(particle,wh2).x += shiftx * boxsize.x;
            getpos(particle,wh2).y += shifty * boxsize.y;
            getpos(particle,wh2).z += shiftz * boxsize.z;
            
            move_cache_water.push_back(getparticle(particle,i));
            move_cache_water.push_back(getparticle(particle,wh1));
            move_cache_water.push_back(getparticle(particle,wh2));
            
            move_to_cellid_water.push_back(ncid);
            move_to_cellid_water.push_back(ncid);
            move_to_cellid_water.push_back(ncid);
            
            number_move_water+=3;
            
            CovalentBondInfo::BondList bond_cache;
            bond_cache.clear();
            bondlistarray[ci].pick_up_bondlist(bond_cache,getatomid(particle,i));
            bondlistarray[ci].pick_up_bondlist(bond_cache,getatomid(particle,wh1));
            bondlistarray[ci].pick_up_bondlist(bond_cache,getatomid(particle,wh2));
            move_bond_cache_water.push_back(bond_cache);
            bondlistarray[ci].remove_atomid(getatomid(particle,i));
            bondlistarray[ci].remove_atomid(getatomid(particle,wh1));
            bondlistarray[ci].remove_atomid(getatomid(particle,wh2));
          }
        }
      }
    }
    if(DebugLog::verbose>1){  // debug_shake
      if(move) std::cout << ">>>>>>>>>>>>>>> select_move_particle = " << move << std::endl;
      if(move) std::cout << " ci= " << ci << " cid=" << cid << std::endl;
      if(move) std::cout << " number_move = " << number_move << std::endl;
      if(move) std::cout << " number_move_water = " << number_move_water << std::endl;
      if(move) std::cout << " number_move_shake = " << number_move_shake << std::endl;
    }
  }
  return move;
}
void 
PostProcess::dumpcell(ParticleArray& particle,
                      std::vector<TypeRange>& typerange)
{
  for(size_t s=0;s<typerange.size();s++){
    std::cout << " cell " << s;
    for(int i=typerange[s].begin;i<typerange[s].end;i++){
      std::cout << " " << i << ":" << particle[i].atomid;
    }
  }
  std::cout << std::endl;
}


///// TODO : fix two or more particles move out from one cell
//! move particle to own cell or set of move out to other node
/*! 
  particle : original array of particles
  typerange : original range of sets(cells)
  move_out_particle : particles move out from this node
  move_out_cellid : id of cells where move_out_particle move to
  move_out_type : type of move_out_particle
  setidtoindex : reverse map set(cell) id to index

  member prepared by select_move_particle
  number_move : number of move particles
  move_index : index of move particles
  move_cache : move particles
  rt_cache : type of move particles
  move_from_cellindex : index of set(cell) where particles move from
  move_to_cellid : id of set(cell) where particles move to

  !!!! move_index must be 
        grouped by cell
        sorted in ascending order at least in one cell
 */
template <typename PA, typename MPA>
bool
PostProcess::move_inside_node(PA& particle,
                              std::vector<TypeRange>& typerange,
                              std::vector<CovalentBondInfo::BondList>& bondlistarray,
                              WaterList& waterlist,
                              ShakeList& shakelist,
                              MPA& move_out_particle,
                              std::vector<int>& move_out_cellid,
                              std::vector<PotentialModel>& move_out_type,
                              std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                              std::map<int,int>& setidtoindex)
{
  // delete move out particle
  /* move particle scan reverse order
     last particle move to index of deleted intermediate particle.
     if last particle is move-particle, it must deleted before delete other particle
   */
  //for(int mi=delete_index.size()-1;mi>=0;mi--){
    //int i = delete_index[mi];
    //int fromindex = delete_from_cellindex[mi];
    //PotentialModel pm = delete_rt_cache[mi];
  for(std::set<int>::reverse_iterator mi=delete_index.rbegin();mi!=delete_index.rend();mi++){
    int i = *mi;
    int fromindex = delete_from_cellindex[i];
    PotentialModel pm = delete_rt_cache[i];
    TypeRange tr = typerange[fromindex];

    if(pm==OnlyLJPotential){
      if(DebugLog::verbose>2){  // debug_shake
        std::cout << "move_inside_node delete OnlyLJPotential i=" << i << " fromindex=" << fromindex << " tr.lj.end-1=" << tr.lj.end-1 << " atomid=" << getatomid(particle,i) << " atomtype=" << getatomtype(particle,i) << std::endl;
      }
      if(i<tr.lj.end-1){
        //        particle[i] = particle[tr.lj.end-1];
        setparticle(particle,getparticle(particle,tr.lj.end-1),i);
      }
      if(tr.ljcoulomb.begin<tr.ljcoulomb.end){
        //        particle[tr.ljcoulomb.begin-1] = particle[tr.ljcoulomb.end-1];
        setparticle(particle,getparticle(particle,tr.ljcoulomb.end-1),tr.ljcoulomb.begin-1);
        if(getatomtype(particle,tr.ljcoulomb.begin-1) == ATOMTYPE_WO){
          waterlist.move_o(tr.ljcoulomb.end-1, tr.ljcoulomb.begin-1);
	}
#ifdef USE_SHAKE
	else if(getatomtype(particle,tr.ljcoulomb.begin-1) < ATOMTYPE_WH){ // H1 atom of shake bond
          shakelist.move_h1(tr.ljcoulomb.end-1, tr.ljcoulomb.begin-1);
        }else{                                                           // Heavy atom of shake bond
          shakelist.move_ha(tr.ljcoulomb.end-1, tr.ljcoulomb.begin-1);
	}
#endif
      }
      if(tr.coulomb.begin<tr.coulomb.end){
        //        particle[tr.coulomb.begin-1] = particle[tr.coulomb.end-1];
        setparticle(particle,getparticle(particle,tr.coulomb.end-1),tr.coulomb.begin-1);
        if(getatomtype(particle,tr.coulomb.begin-1) == ATOMTYPE_WH)
          waterlist.move_h(tr.coulomb.end-1, tr.coulomb.begin-1);
      }
      tr.end--;
      tr.lj.end--;
      tr.ljcoulomb.shift(-1);
      tr.coulomb.shift(-1);
    }else if(pm==LJCoulombPotential){
      if(DebugLog::verbose>1){  // debug_shake
        std::cout << "move_inside_node delete LJCoulombPotential i=" << i << " fromindex=" << fromindex << " tr.ljcoulomb.end-1=" << tr.ljcoulomb.end-1 << " atomid=" << getatomid(particle,i) << " atomtype=" << getatomtype(particle,i) << std::endl;
      }
      if(i<tr.ljcoulomb.end-1){
        //        particle[i] = particle[tr.ljcoulomb.end-1];
        setparticle(particle,getparticle(particle,tr.ljcoulomb.end-1),i);
        if(getatomtype(particle,i) == ATOMTYPE_WO){
          waterlist.move_o(tr.ljcoulomb.end-1, i);
        }
#ifdef USE_SHAKE
	else if(getatomtype(particle,i) < ATOMTYPE_WH){ // H1 atom of shake bond
          shakelist.move_h1(tr.ljcoulomb.end-1, i);
        }else{                                        // Heavy atom of shake bond
          shakelist.move_ha(tr.ljcoulomb.end-1, i);
        }
#endif
      }
      if(tr.coulomb.begin<tr.coulomb.end){
        //        particle[tr.coulomb.begin-1] = particle[tr.coulomb.end-1];
        setparticle(particle,getparticle(particle,tr.coulomb.end-1),tr.coulomb.begin-1);
        if(getatomtype(particle,tr.coulomb.begin-1) == ATOMTYPE_WH)
          waterlist.move_h(tr.coulomb.end-1, tr.coulomb.begin-1);
      }
      tr.end--;
      tr.ljcoulomb.end--;
      tr.coulomb.shift(-1);
    }else if(pm==OnlyCoulombPotential){
      if(DebugLog::verbose>1){  // debug_shake
        std::cout << "move_inside_node delete OnlyCoulombPotential i=" << i << " fromindex=" << fromindex << " tr.coulomb.end-1=" << tr.coulomb.end-1 << " atomid=" << getatomid(particle,i) << " atomtype=" << getatomtype(particle,i) << std::endl;
      }
      if(i<tr.coulomb.end-1){
        //        particle[i] = particle[tr.coulomb.end-1];
        setparticle(particle,getparticle(particle,tr.coulomb.end-1),i);
        if(getatomtype(particle,i) == ATOMTYPE_WH)
          waterlist.move_h(tr.coulomb.end-1, i);
      }
      tr.end--;
      tr.coulomb.end--;
    }else{
      std::cout << " Unsupported Potential Model " << i << " move from " << fromindex << std::endl;
      tr.end--;
    }
    typerange[fromindex] = tr;
  }


  bool cell_over_flow = false;

  // append move particle
  for(size_t mi=0;mi<number_move;mi++){
    int tocellid = move_to_cellid[mi];
    PotentialModel pm = move_rt_cache[mi];
    std::map<int,int>::iterator toi = setidtoindex.find(tocellid);

    if(toi!=setidtoindex.end()){  // to cell in this node
      int tocellindex = toi->second;
      if(typerange[tocellindex].end-typerange[tocellindex].begin>=max_particle_in_cell){
	printf("cell(index %d) overflow %d add %lu/%d move particle\n",tocellindex,max_particle_in_cell,mi,number_move);
	cell_over_flow = true;
      }
      //      int fromindex = move_from_cellindex[mi];
      //      if(tocellindex==fromindex)continue; 
      /* may be tocellindex!=fromindex
         and particle was deleted from set(cell)
       */
      //bool aps=
      add_particle(particle,typerange[tocellindex], waterlist, shakelist, move_cache[mi],pm);
      //      if(aps==false)
      //        {
          //    std::cout << "move index " << i << " atomid " << move_cache[i].atomid << " from " << fromindex << "(" <<tr.begin << "-" << tr.end << ")" << " to " << tocellindex << std::endl;
          //     std::cout << "move index " << i << " atomid " << move_cache[mi].atomid << " from " << fromindex << " to " << tocellindex << std::endl;
      //  }
      bondlistarray[tocellindex].push_back(move_bond_cache[mi]);
        
    }else{   // to cell is not in this node
      move_out_particle.push_back(move_cache[mi]);
      move_out_cellid.push_back(tocellid);
      move_out_type.push_back(pm);
      move_out_bondlistarray.push_back(move_bond_cache[mi]);
    }
  }

#ifdef USE_SHAKE
  // append shake
  int mos=0;
  int nh1 = 0;
  int atomid = -999;
  for(size_t mi=0;mi<number_move_shake;mi+=2){  // every (HA,H1) unit
    int tocellid = move_to_cellid_shake[mi];
    std::map<int,int>::iterator toi = setidtoindex.find(tocellid);

    if(toi!=setidtoindex.end()){  // tocell in this node
      if(move_cache_shake[mi].atomtype > ATOMTYPE_WH){  // heavy atom
        int tocellindex = toi->second;
	AtomID id = move_cache_shake[mi+1].atomid;
	int bondtype=-1;
	int bcs = mi/2;
	int bi;
	for(bi=0;bi<move_bond_cache_shake[bcs].BondArray.size();bi++){
	  if((move_bond_cache_shake[bcs].BondArray[bi].id_of_atom[0]==id)
	     ||(move_bond_cache_shake[bcs].BondArray[bi].id_of_atom[1]==id)){
	    bondtype = move_bond_cache_shake[bcs].BondArray[bi].typeofbond;
	    break;
	  }
	}
        if(atomid != move_cache_shake[mi].atomid){  // 1st H
	  nh1 = 0;
	}else{   // not 1st H
	  if(bondtype==-1){   // bond not found move_bond_cache_shake[mi/2]
	    bcs=(mi-nh1*2)/2;   // index of move_bond_cache_shake 1st bond of same heavy atom
	    for(bi=0;bi<move_bond_cache_shake[bcs].BondArray.size();bi++){
	      if((move_bond_cache_shake[bcs].BondArray[bi].id_of_atom[0]==id)
		 ||(move_bond_cache_shake[bcs].BondArray[bi].id_of_atom[1]==id)){
		bondtype = move_bond_cache_shake[bcs].BondArray[bi].typeofbond;
		break;
	      }
	    }
	  }
	}
	if(bondtype==-1){
	  printf("not found bond in move_bond_cache_shake for %d %d\n",move_cache_shake[mi].atomid,id);
	}
        bool aps=
            add_shake(particle, typerange[tocellindex], waterlist, shakelist, // add atoms to particle & typerange & shakelist
                      nh1++, move_cache_shake[mi], move_cache_shake[mi+1], bondtype);   // modify waterlist, if needs
        bondlistarray[tocellindex].push_back(move_bond_cache_shake[mi/2]);
        atomid = move_cache_shake[mi].atomid;
      }
    }else{   // to cell not in this node
      move_out_particle.push_back(move_cache_shake[mi]);
      move_out_cellid.push_back(tocellid);
      move_out_type.push_back(LJCoulombPotential);

      move_out_particle.push_back(move_cache_shake[mi+1]);
      move_out_cellid.push_back(tocellid);
      move_out_type.push_back(LJCoulombPotential);

      move_out_bondlistarray.push_back(move_bond_cache_shake[mi/2]);
      move_out_bondlistarray.resize(move_out_bondlistarray.size()+1);
      mos++;
    }
  }
  if((DebugLog::verbose>2)
     ||((DebugLog::verbose==2)&&(mos>0))){
    std::cout << " move out shake " << mos << std::endl;
  }
#endif

  // append water
  int mow=0;
  for(size_t mi=0;mi<number_move_water;mi+=3){
    int tocellid = move_to_cellid_water[mi];
    std::map<int,int>::iterator toi = setidtoindex.find(tocellid);

    if(toi!=setidtoindex.end()){  // tocell in this node
      if(move_cache_water[mi].atomtype == ATOMTYPE_WO){
        int tocellindex = toi->second;
	if(typerange[tocellindex].end-typerange[tocellindex].begin>=max_particle_in_cell-2){
	  printf("cell(index %d) overflow\n",tocellindex);
	  cell_over_flow = true;
	}
        //bool aps=
            add_water(particle, typerange[tocellindex], waterlist,
                      move_cache_water[mi], move_cache_water[mi+1], move_cache_water[mi+2]);
        bondlistarray[tocellindex].push_back(move_bond_cache_water[mi/3]);
      }
    }else{   // to cell not in this node
      move_out_particle.push_back(move_cache_water[mi]);
      move_out_cellid.push_back(tocellid);
      move_out_type.push_back(LJCoulombPotential);

      move_out_particle.push_back(move_cache_water[mi+1]);
      move_out_cellid.push_back(tocellid);
      move_out_type.push_back(OnlyCoulombPotential);

      move_out_particle.push_back(move_cache_water[mi+2]);
      move_out_cellid.push_back(tocellid);
      move_out_type.push_back(OnlyCoulombPotential);

      move_out_bondlistarray.push_back(move_bond_cache_water[mi/3]);
      move_out_bondlistarray.resize(move_out_bondlistarray.size()+2);
      mow++;
    }
  }
  if((DebugLog::verbose>2)
     ||((DebugLog::verbose==2)&&(mow>0))){
    std::cout << " move out water " << mow << std::endl;
  }
  return !cell_over_flow;
}
/*
void
PostProcess::select_and_move_inside_node(ParticleArray& particlearray,
                                         std::vector<TypeRange>& typerangearray,
                                         WaterList& waterlist,
                                         ShakeList& shakelist,
                                         const std::vector<int>& cellid,
                                         ParticleArray& move_out_particles,
                                         std::vector<int>& move_out_cellid,
                                         std::vector<PotentialModel>& move_out_type,
                                         std::map<int,int>& setid_to_index)
{
  select_move_particle(particlearray, typerangearray, waterlist, shakelist, cellid);
  move_inside_node(particlearray, typerangearray, waterlist, shakelist,
                   move_out_particles, move_out_cellid, move_out_type,
                   setid_to_index);
}
*/
template <typename PA, typename MPA>
bool
PostProcess::select_and_move_inside_node(PA& particlearray,
                                         std::vector<TypeRange>& typerangearray,
                                         std::vector<CovalentBondInfo::BondList>& bondlistarray,

                                         WaterList& waterlist,
                                         ShakeList& shakelist,
                                         const std::vector<int>& cellid,
                                         MPA& move_out_particles,
                                         std::vector<int>& move_out_cellid,
                                         std::vector<PotentialModel>& move_out_type,
                                         std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                                         std::map<int,int>& setid_to_index,
                                         int& over_move)
{
  select_move_particle(particlearray, typerangearray, 
                       bondlistarray, waterlist, shakelist, cellid, over_move);
  bool not_overflow;
  not_overflow = move_inside_node(particlearray, typerangearray, 
                                  bondlistarray, waterlist, shakelist,
                                  move_out_particles, move_out_cellid, move_out_type,
                                  move_out_bondlistarray,
                                  setid_to_index);
  return not_overflow;
}
template 
bool
PostProcess::select_and_move_inside_node(ParticleArray& particlearray,
                                         std::vector<TypeRange>& typerangearray,
                                         std::vector<CovalentBondInfo::BondList>& bondlistarray,

                                         WaterList& waterlist,
                                         ShakeList& shakelist,
                                         const std::vector<int>& cellid,
                                         ParticleArray& move_out_particles,
                                         std::vector<int>& move_out_cellid,
                                         std::vector<PotentialModel>& move_out_type,
                                         std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                                         std::map<int,int>& setid_to_index,
                                         int& over_move);
template 
bool
PostProcess::select_and_move_inside_node(CombinedParticleArray& particlearray,
                                         std::vector<TypeRange>& typerangearray,
                                         std::vector<CovalentBondInfo::BondList>& bondlistarray,

                                         WaterList& waterlist,
                                         ShakeList& shakelist,
                                         const std::vector<int>& cellid,
                                         CombinedParticleArray& move_out_particles,
                                         std::vector<int>& move_out_cellid,
                                         std::vector<PotentialModel>& move_out_type,
                                         std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                                         std::map<int,int>& setid_to_index,
                                         int& over_move);

template 
bool
PostProcess::select_and_move_inside_node(CombinedParticleArray& particlearray,
                                         std::vector<TypeRange>& typerangearray,
                                         std::vector<CovalentBondInfo::BondList>& bondlistarray,

                                         WaterList& waterlist,
                                         ShakeList& shakelist,
                                         const std::vector<int>& cellid,
                                         ParticleArray& move_out_particles,
                                         std::vector<int>& move_out_cellid,
                                         std::vector<PotentialModel>& move_out_type,
                                         std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                                         std::map<int,int>& setid_to_index,
                                         int& over_move);

template<typename PA, typename MPA>
bool
PostProcess::merge(PA& particle,
                   std::vector<int>& cellid,
                   std::map<int,int>& cellidtoindex,
                   std::vector<TypeRange>& typerange,
                   std::vector<CovalentBondInfo::BondList>& bondlistarray,
                   WaterList& waterlist,
                   ShakeList& shakelist,
                   MPA& move_in_particle,
                   std::vector<int>& move_in_cellid,
                   std::vector<PotentialModel>& move_in_type,
                   std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                   const ShakeList& sl)
{
  bool cell_over_flow = false;
  int nh1 = 0;
  int atomid = -999;
  for(int i=0;i<int(move_in_particle.size());i++){
    Position pos = getpos(move_in_particle,i);
    SpaceVector<int> ncell(int(pos.x*invcellsize.x),
                           int(pos.y*invcellsize.y),
                           int(pos.z*invcellsize.z));
    int ncid = cell_position_to_cellid(ncell,celldiv);

    if(getatomtype(move_in_particle,i) == ATOMTYPE_WH){
      std::cout << "ERROR in merge == ATOMTYPE_WH: atomid=" << getatomid(move_in_particle,i) << std::endl;
      continue;
    }else{
      int ci = cellidtoindex[ncid];
      if(getatomtype(move_in_particle,i) == ATOMTYPE_WO){
	if(typerange[ci].end-typerange[ci].begin>=max_particle_in_cell-2){
	  printf("cell(index %d) overflow move in water\n",ci);
	  cell_over_flow = true;
	}
        add_water(particle, typerange[ci], waterlist,
                  getparticle(move_in_particle,i), getparticle(move_in_particle,i+1), getparticle(move_in_particle,i+2));
	i+=2;
      }
#ifdef USE_SHAKE
      /* In copy_mode, sl is not complete.
      else if(sl.find(getatomid(move_in_particle,i)) != sl.end()){  // heavy atom of shake bond
        if(DebugLog::verbose>1){  // debug_shake
          std::cout << "in merge: HIT HA: atomid=" << getatomid(move_in_particle,i) << std::endl;
          std::cout << "in merge: HIT HA: atomtype=" << getatomtype(move_in_particle,i) << std::endl;
          std::cout << "in merge:     H1: atomid=" << getatomid(move_in_particle,i+1) << std::endl;
          std::cout << "in merge:     H1: atomtype=" << getatomtype(move_in_particle,i+1) << std::endl;
        }

        if(atomid != getatomid(move_in_particle,i)) nh1 = 0;

        add_shake(particle, typerange[ci], waterlist, shakelist,
                  nh1++, getparticle(move_in_particle,i), getparticle(move_in_particle,i+1));
        bondlistarray[ci].push_back(move_out_bondlistarray[i]);

        atomid = getatomid(move_in_particle,i);
        i++;   // every (HA,H1) unit
      }else if(getatomtype(move_in_particle,i) < ATOMTYPE_WH){
        std::cout << "ERROR in merge < ATOMTYPE_WH: atomid=" << getatomid(move_in_particle,i) << std::endl;
      }
      */
      else if((i<move_in_particle.size()-1)&&(getatomtype(move_in_particle,i+1)<ATOMTYPE_WH)) {  // next atom is H, this atom is heavy atom of shake bond
        if(DebugLog::verbose>1){  // debug_shake
          std::cout << "in merge: HIT HA: atomid=" << getatomid(move_in_particle,i) << std::endl;
          std::cout << "in merge: HIT HA: atomtype=" << getatomtype(move_in_particle,i) << std::endl;
          std::cout << "in merge:     H1: atomid=" << getatomid(move_in_particle,i+1) << std::endl;
          std::cout << "in merge:     H1: atomtype=" << getatomtype(move_in_particle,i+1) << std::endl;
        }
	AtomID id = getatomid(move_in_particle,i+1);
	int bondtype=-1;
	int bcs = i;
	int bi;
	for(bi=0;bi<move_out_bondlistarray[bcs].BondArray.size();bi++){
	  if((move_out_bondlistarray[bcs].BondArray[bi].id_of_atom[0]==id)
	     ||(move_out_bondlistarray[bcs].BondArray[bi].id_of_atom[1]==id)){
	    bondtype = move_out_bondlistarray[bcs].BondArray[bi].typeofbond;
	    break;
	  }
	}
        if(atomid != getatomid(move_in_particle,i)){
	  nh1 = 0;
	}else{
	  if(bondtype==-1){
	    bcs = i-nh1*2;
	    for(bi=0;bi<move_out_bondlistarray[bcs].BondArray.size();bi++){
	      if((move_out_bondlistarray[bcs].BondArray[bi].id_of_atom[0]==id)
		 ||(move_out_bondlistarray[bcs].BondArray[bi].id_of_atom[1]==id)){
		bondtype = move_out_bondlistarray[bcs].BondArray[bi].typeofbond;
		break;
	      }
	    }
	  }
	}
	if(bondtype==-1){
	  printf("not found bond in move_out_bondlistarray for %d %d %d\n",getatomid(move_in_particle,i),id,nh1);
	}
        add_shake(particle, typerange[ci], waterlist, shakelist,
                  nh1++, getparticle(move_in_particle,i), getparticle(move_in_particle,i+1), bondtype);
	bondlistarray[ci].push_back(move_out_bondlistarray[i]);

        atomid = getatomid(move_in_particle,i);
        i++;   // every (HA,H1) unit
      }
#endif
      else {
	if(typerange[ci].end-typerange[ci].begin>=max_particle_in_cell){
	  printf("cell(index %d) overflow by move in\n",ci);
	  cell_over_flow = true;
	}
        add_particle(particle, typerange[ci], waterlist, shakelist,
                     getparticle(move_in_particle,i), move_in_type[i]);
        bondlistarray[ci].push_back(move_out_bondlistarray[i]);
      }
    }
  }
  return !cell_over_flow;
}
template
bool
PostProcess::merge(ParticleArray& particle,
                   std::vector<int>& cellid,
                   std::map<int,int>& cellidtoindex,
                   std::vector<TypeRange>& typerange,
                   std::vector<CovalentBondInfo::BondList>& bondlistarray,
                   WaterList& waterlist,
                   ShakeList& shakelist,
                   ParticleArray& move_in_particle,
                   std::vector<int>& move_in_cellid,
                   std::vector<PotentialModel>& move_in_type,
                   std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                   const ShakeList& sl
                   );
template
bool
PostProcess::merge(CombinedParticleArray& particle,
                   std::vector<int>& cellid,
                   std::map<int,int>& cellidtoindex,
                   std::vector<TypeRange>& typerange,
                   std::vector<CovalentBondInfo::BondList>& bondlistarray,
                   WaterList& waterlist,
                   ShakeList& shakelist,
                   ParticleArray& move_in_particle,
                   std::vector<int>& move_in_cellid,
                   std::vector<PotentialModel>& move_in_type,
                   std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                   const ShakeList& sl
                   );
template
bool
PostProcess::merge(CombinedParticleArray& particle,
                   std::vector<int>& cellid,
                   std::map<int,int>& cellidtoindex,
                   std::vector<TypeRange>& typerange,
                   std::vector<CovalentBondInfo::BondList>& bondlistarray,
                   WaterList& waterlist,
                   ShakeList& shakelist,
                   CombinedParticleArray& move_in_particle,
                   std::vector<int>& move_in_cellid,
                   std::vector<PotentialModel>& move_in_type,
                   std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                   const ShakeList& sl);

// shift position of particles listed in shift_cell
//  cell
void
PostProcess::postprocess_receive(GhostParticleArray& particle, 
                                 const std::vector<int>& cellid,
                                 const std::map<int,int>& cellidtoindex,
                                 const std::vector<TypeRange>& typerange)
{
  for(std::vector<int>::iterator it = shift_cell.plusx.begin();
      it!=shift_cell.plusx.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.x += boxsize.x;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusx.begin();
      it!=shift_cell.minusx.end(); ++it){
    //int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.x -= boxsize.x;
    }
  }

  for(std::vector<int>::iterator it = shift_cell.plusy.begin();
      it!=shift_cell.plusy.end(); ++it){
    //int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.y += boxsize.y;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusy.begin();
      it!=shift_cell.minusy.end(); ++it){
    //int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.y -= boxsize.y;
    }
  }

  for(std::vector<int>::iterator it = shift_cell.plusz.begin();
      it!=shift_cell.plusz.end(); ++it){
    //int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.z += boxsize.z;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusz.begin();
      it!=shift_cell.minusz.end(); ++it){
    //int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.z -= boxsize.z;
    }
  }
}

void
PostProcess::postprocess_receive_with_shiftidcheck(GhostParticleArray& particle, 
                                 const std::vector<int>& cellid,
                                 const std::map<int,int>& cellidtoindex,
                                 const std::vector<TypeRange>& typerange)
{
  for(std::vector<int>::iterator it = shift_cell.plusx.begin();
      it!=shift_cell.plusx.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.x += boxsize.x;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusx.begin();
      it!=shift_cell.minusx.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.x -= boxsize.x;
    }
  }

  for(std::vector<int>::iterator it = shift_cell.plusy.begin();
      it!=shift_cell.plusy.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.y += boxsize.y;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusy.begin();
      it!=shift_cell.minusy.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.y -= boxsize.y;
    }
  }

  for(std::vector<int>::iterator it = shift_cell.plusz.begin();
      it!=shift_cell.plusz.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.z += boxsize.z;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusz.begin();
      it!=shift_cell.minusz.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle.poscharge[i].position.z -= boxsize.z;
    }
  }
}

void
PostProcess::postprocess_receive(ParticleArray& particle, 
                                 const std::vector<int>& cellid,
                                 const std::map<int,int>& cellidtoindex,
                                 const std::vector<TypeRange>& typerange)
{
  for(std::vector<int>::iterator it = shift_cell.plusx.begin();
      it!=shift_cell.plusx.end(); ++it){
    //int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.x += boxsize.x;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusx.begin();
      it!=shift_cell.minusx.end(); ++it){
    //int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.x -= boxsize.x;
    }
  }

  for(std::vector<int>::iterator it = shift_cell.plusy.begin();
      it!=shift_cell.plusy.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.y += boxsize.y;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusy.begin();
      it!=shift_cell.minusy.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.y -= boxsize.y;
    }
  }

  for(std::vector<int>::iterator it = shift_cell.plusz.begin();
      it!=shift_cell.plusz.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.z += boxsize.z;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusz.begin();
      it!=shift_cell.minusz.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set = cellidtoindex.find(*it)->second;
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.z -= boxsize.z;
    }
  }
}

void
PostProcess::postprocess_receive_with_shiftidcheck(ParticleArray& particle, 
                                              const std::vector<int>& cellid,
                                              const std::map<int,int>& cellidtoindex,
                                              const std::vector<TypeRange>& typerange)
{
  for(std::vector<int>::iterator it = shift_cell.plusx.begin();
      it!=shift_cell.plusx.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.x += boxsize.x;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusx.begin();
      it!=shift_cell.minusx.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.x -= boxsize.x;
    }
  }

  for(std::vector<int>::iterator it = shift_cell.plusy.begin();
      it!=shift_cell.plusy.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.y += boxsize.y;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusy.begin();
      it!=shift_cell.minusy.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.y -= boxsize.y;
    }
  }

  for(std::vector<int>::iterator it = shift_cell.plusz.begin();
      it!=shift_cell.plusz.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.z += boxsize.z;
    }
  }
  for(std::vector<int>::iterator it = shift_cell.minusz.begin();
      it!=shift_cell.minusz.end(); ++it){
    //    int set = cellidtoindex[*it];
    int set;
    std::map<int,int>::const_iterator cit = cellidtoindex.find(*it);
    if(cit!=cellidtoindex.end()){
      set = cit->second;
    }else{
      std::cout << "shift_cell not found in ghost" << std::endl;
      continue;
    }
    for(int i=typerange[set].begin;i<typerange[set].end;i++){
      particle[i].position.z -= boxsize.z;
    }
  }
}

void 
PostProcess::change_fine_margin(int interval)
{
  if(interval<10){
    fine_margin.x = margin.x/interval;
    fine_margin.y = margin.y/interval;
    fine_margin.z = margin.z/interval;
  }else{
    fine_margin.x = 0.0;
    fine_margin.y = 0.0;
    fine_margin.z = 0.0;
  }
}
