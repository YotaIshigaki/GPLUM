#ifndef CELLINDEX_H
#define CELLINDEX_H

#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <iostream>

#include <cassert>

#include "Common.h"
#include "ParticleInfo.h"
#include "OperationSelector.h"

#include "Log.h"

enum CellIndexType {
  FullCube,
  HalfShell,
  FullShell,
  SmallBall
};


void
dump_cellindextype(CellIndexType cit);

//! check n is 2^m
bool
power_of_two(int n, int &p);

//! divide 2^m unit to 2^mx 2^my 2^mz
// mx <= my <= mz
void
split_power_of_two(int p, int sp[3]);

bool
is_square(int np, int m[2]);

bool
is_double_square(int np, int m[2]);

bool
power_of_two_of_square(int np, int m[2]);

bool
is_cube(int np, int m[3]);

bool
is_double_cube(int np, int m[3]);

bool
is_quad_cube(int np, int m[3]);

bool
power_of_two_of_cube(int np, int m[3]);

/*
  cell_id : x, y, z order
 */
inline int
cell_position_to_cellid(const SpaceVector<int>& cell_position, 
                        const SpaceVector<int>celldiv)
{
  return cell_position.x + (cell_position.y
                            + cell_position.z*celldiv.y)*celldiv.x;
}
    
inline
SpaceVector<int> cellid_to_cell_position(const int cellid, 
                                         const SpaceVector<int>celldiv)
{
  SpaceVector<int> pos;
  pos.z = cellid/(celldiv.x*celldiv.y);
  pos.y = (cellid-pos.z*(celldiv.x*celldiv.y))/celldiv.x;
  pos.x = cellid-(pos.y+pos.z*celldiv.y)*celldiv.x;
  return pos;
}

inline int
cellid_to_cell_position_x(const int cellid, 
                          const SpaceVector<int>celldiv)
{
  return cellid%celldiv.x;
}

inline int
cellid_to_cell_position_y(const int cellid, 
                          const SpaceVector<int>celldiv)
{
  return (cellid/celldiv.x)%celldiv.y;
}

inline int
cellid_to_cell_position_z(const int cellid, 
                          const SpaceVector<int>celldiv)
{
  return cellid/(celldiv.x*celldiv.y);
}

inline double
min_cell_distance_nomargin(SpaceVector<int> ci, SpaceVector<int> cj,
                  SpaceVector<double> cellsize)
{
  SpaceVector<int> dc = cj - ci;
  int dx_min = std::max(abs(dc.x)-1,0);
  int dy_min = std::max(abs(dc.y)-1,0);
  int dz_min = std::max(abs(dc.z)-1,0);
  double dx = cellsize.x*dx_min;
  double dy = cellsize.y*dy_min;
  double dz = cellsize.z*dz_min;
  return sqrt(dx*dx+dy*dy+dz*dz);
}

inline double
min_cell_distance(SpaceVector<int> ci, SpaceVector<int> cj,
                  SpaceVector<double> cellsize, SpaceVector<double> margin){
  SpaceVector<int> dc = cj - ci;
  int dx_min = std::max(abs(dc.x)-1,0);
  int dy_min = std::max(abs(dc.y)-1,0);
  int dz_min = std::max(abs(dc.z)-1,0);
  double dx = std::max(cellsize.x*dx_min-2.0*margin.x,0.0);
  double dy = std::max(cellsize.y*dy_min-2.0*margin.y,0.0);
  double dz = std::max(cellsize.z*dz_min-2.0*margin.z,0.0);
  return sqrt(dx*dx+dy*dy+dz*dz);
}


class CellIndex {
public:
  // common infomation for each node
  //// set by set_cell_parameter
  SpaceVector<double> boxsize;    //!< size of system
  double volume;                  //!< volume of system
  SpaceVector<int> celldiv;       //!< number of cell in each axis
  SpaceVector<double> cellsize;   //!< size of one cell
  SpaceVector<double> invcellsize;//!< 1.0/size of one cell
  double cutoff;                  //!< cutoff length
  double bond_length;             //!< bond length
  SpaceVector<double> margin;     //!< margin of move out
  SpaceVector<int> target_range;  //!< target range in one axis
  SpaceVector<int> move_range;    //!< 1 except very small cell
  SpaceVector<int> bond_range;    //!< bond target range in one axis
  //// set by generate_cell_index
  std::vector< SpaceVector<int> > cell_position;    //!< postion(indexies)
  int num_cell;                                     //!< number of cell
  //// set by distribute_cell
  int num_node;                   //!< number of node
  SpaceVector<int> celldiv_node;  //!< number of cell in each axis for one node
  std::vector< std::vector<int> > cell_id_in_node;  //!< cell id in each node
  std::map<int,int> cell_id_to_node;                //!< reverse map cell id to node

  SpaceVector<double> node_size;  //!< size of node
  double node_volume;             //! volume of one node
  double node_move_surface_volume; //!< surface volume contain move out candidate
  double node_move_surface_ratio; //!< node_move_surface_volume/node_volume

  void
  set_cell_parameter(double bs, int cd, double mg, double ct, 
                     double bond_len=1.7);
  void
  set_cell_parameter(SpaceVector<double> bs, 
                     SpaceVector<int> cd, SpaceVector<double> mg,
                     double ct, 
                     double bond_len=1.7);
  void
  change_boxsize(SpaceVector<double> bs);
  double
  calc_volume();

  //! set cubic cell position(cell index)
  // order od cell_position : x y z
  /* (0,0,0), (1,0,0), ... ,(celldiv-1,0,0), (0,1,0), ... ,(celldiv-1,2,0),
     (0,3,0), ... , (celldiv-1,celldiv-1,celldiv-2),
     (0,celldive-1,celldiv-1), ... ,(celldiv-1,celldiv-1,celldiv-1)
  */
  size_t
  generate_cell_index();

  //! distribute cells to node
  bool
  distribute_cell(const SpaceVector<int> nodediv);

  void
  dump_cell(int id);

  void
  dump_cell_all();

  void
  dump_cell_subset(std::vector<int>& id);


};

#include "CellMethod.h"

class CellSubset {
public:
  CellIndex cellindex;
  CellMethodModule::CellMethod cellmethod;

  // node depend information
  int mynode;
  // subunit depend 
  std::vector<int> mycell;  //!< cell id in this unit
  int num_mycell;           //!< number of cell in this unit
  std::vector< std::vector<int> > send_target_cell;  //!< cell id that requires mycell 
  std::vector< std::vector<int> > recv_target_cell;  //!< cell id required mycell
  std::vector< std::vector<int> > short_target_cell; //!< cell id required short range calculation
  std::vector<ShiftCellArray> shift_cell;
  ShiftCellArray shift_cell_all;                     //!< cell required periodic shift in recv_target_cell
  std::vector<std::pair<int, int> > ghost_cell_pairs;//!< cell id ghost-ghost
  std::vector< std::vector<int> > self_target_cell;  //!< part of recv_target_cell that owned this node

  CellSubset(){}

  CellSubset(const CellIndex& ci, const int& mn);

  CellSubset(const CellSubset* cs);

  void
  merge_shift_cell();

  int
  calc_self_target();

  int
  calc_target_gather(std::vector<int>& recv_target_node,
                     std::vector< std::vector<int> >& recv_target_dist);

  int
  calc_target_scatter(std::vector<int>& send_target_node,
                      std::vector< std::vector<int> >& send_target_dist);

  int
  calc_target_move(std::vector<int>& move_target_node,
                   std::vector< std::vector<int> >& move_target_set);

  void
  makefulltarget(std::vector<int> subset, 
                 CellIndexType cellindextype,
                 OperationSelector operations);

  bool
  makecommtarget(std::vector<int>& send_target_node,
                 std::vector< std::vector<int> >& send_target_set,
                 std::vector<int>& recv_target_node,
                 std::vector< std::vector<int> >& recv_target_set,
                 std::vector<int>& move_target_node,
                 std::vector< std::vector<int> >& move_target_set);

  void
  dump_cell_target(int cell_index, std::vector<int>& target_cell);

  void
  dump_cell_target_send(int cell_index);

  void
  dump_cell_target_recv(int cell_index);

  void
  dump_cell_target_self(int cell_index);

  template<class PA, class SPA>
    int
  distribute_particle_cell(PA& particlearray,
                           std::vector<int>& particle_setid,
                           std::vector<TypeRange>& typerangearray,
                           WaterList& waterlist,
                           ShakeList& shakelist,
                           const SPA& source,
                           const std::vector<PotentialModel>& pmodel,
                           const WaterList& wl,
                           const ShakeList& sl
                           );

};


class PostProcess {
public:

  int number_move;                        //!< number of move particles

  //std::vector<int> delete_index;            //!< store index of move particles
  //std::vector<PotentialModel> delete_rt_cache;        //!< store type of move particles
  //std::vector<int> delete_from_cellindex;   //!< store where particles move from

  int max_particle_in_cell;

  std::set<int> delete_index;
  std::map<int,PotentialModel> delete_rt_cache;
  std::map<int,int> delete_from_cellindex;

  ParticleArray move_cache;               //!< store move particles
  std::vector<PotentialModel> move_rt_cache;        //!< store type of move particles
  std::vector<int> move_to_cellid;        //!< store where particles move to

  // water
  std::vector<int> delete_wh_list;
  int number_move_water;                     
  ParticleArray move_cache_water;            
  std::vector<int> move_to_cellid_water;     

  // shake bond
  int number_move_shake;
  ParticleArray move_cache_shake;
  std::vector<int> move_to_cellid_shake;

  // bond
  std::vector<CovalentBondInfo::BondList> move_bond_cache;
  std::vector<CovalentBondInfo::BondList> move_bond_cache_water;
  std::vector<CovalentBondInfo::BondList> move_bond_cache_shake;
  

//   std::vector<int> move_index;            //!< store index of move particles
//   ParticleArray move_cache;               //!< store move particles
//   std::vector<PotentialModel> rt_cache;   //!< store potential model of each move particle
//   std::vector<int> move_from_cellindex;   //!< store where particles move from
//   std::vector<int> move_to_cellid;        //!< store where particles move to

  std::vector< SpaceVector<int> > cell_position;
  SpaceVector<int> celldiv;
  SpaceVector<double> cellsize;
  SpaceVector<double> invcellsize;
  SpaceVector<double> margin;
  SpaceVector<double> fine_margin; //!< small margin of move out when check interval > 1;
  SpaceVector<double> boxsize;

  double volume_margin;
  
  ShiftCellArray shift_cell;

  PostProcess();

  PostProcess(size_t psize);

  PostProcess(size_t psize, SpaceVector<double> bs, ShiftCellArray& sc);

  void reserve(size_t psize);

  double calc_volume_margin();

  void change_boxsize(SpaceVector<double> bs);

  //! listing particles which move out other cell
  /*!
    @param[in,out] particle list of particles
    @param[in] typerange list of typerage
    @param[in,out] waterlist list of rigid water
    @param[in,out] shakelist list of shake bond
    @param[in] cellid list of cell IDs
    @return true: move particle found, false: no move particle found 
    @note position of Hydorogen in rigid water is shifted followed Oxygen.
   */
  template <typename PA>
  bool select_move_particle(PA& particle,
                            std::vector<TypeRange>& typerange,
                            std::vector<CovalentBondInfo::BondList>& bondlistarray,
                            WaterList& waterlist,
                            ShakeList& shakelist,
                            const std::vector<int>& cellid,
                            int& over_move);

  void dumpcell(ParticleArray& particle,std::vector<TypeRange>& typerange);

  template <typename PA, typename MPA>
  bool move_inside_node(PA& particle,
                        std::vector<TypeRange>& typerange,
                        std::vector<CovalentBondInfo::BondList>& bondlistarray,
                        WaterList& waterlist,
                        ShakeList& shakelist,
                        MPA& move_out_particle,
                        std::vector<int>& move_out_cellid,
                        std::vector<PotentialModel>& move_out_type,
                        std::vector<CovalentBondInfo::BondList>&  move_out_bondlistarray,
                        std::map<int,int>& setidtoindex);

  /*
  void select_and_move_inside_node(ParticleArray& particlearray,
                                   std::vector<TypeRange>& typerangearray,
                                   WaterList& waterlist,
                                   ShakeList& shakelist,
                                   const std::vector<int>& cellid,
                                   ParticleArray& move_out_particles,
                                   std::vector<int>& move_out_cellid,
                                   std::vector<PotentialModel>& move_out_type,
                                   std::map<int,int>& setid_to_index);
  */

  template <typename PA, typename MPA>
  bool select_and_move_inside_node(PA& particlearray,
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
                                   int& over_move);

  template<typename PA, typename MPA>
  bool merge(PA& particle,
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
             const ShakeList& sl);

  // shift position of particles listed in shift_cell
  void postprocess_receive(GhostParticleArray& particle, 
                           const std::vector<int>& cellid,
                           const std::map<int,int>& cellidtoindex,
                           const std::vector<TypeRange>& typerange);
  void postprocess_receive_with_shiftidcheck(GhostParticleArray& particle, 
                           const std::vector<int>& cellid,
                           const std::map<int,int>& cellidtoindex,
                           const std::vector<TypeRange>& typerange);
  void postprocess_receive(ParticleArray& particle, 
                           const std::vector<int>& cellid,
                           const std::map<int,int>& cellidtoindex,
                           const std::vector<TypeRange>& typerange);
  void postprocess_receive_with_shiftidcheck(ParticleArray& particle, 
                           const std::vector<int>& cellid,
                           const std::map<int,int>& cellidtoindex,
                           const std::vector<TypeRange>& typerange);
  
  /*  for NearestXYZ
   */
  /* 
     particle position to cell position in x,y,z
   */
  inline int position_x_to_cell_x(double x){
    //    std::cout << x << "*" << invcellsize.x << "->" << int(x*invcellsize.x) << std::endl;
    return int(x*invcellsize.x);
  }
  inline int position_y_to_cell_y(double y){
    return int(y*invcellsize.y);
  }
  inline int position_z_to_cell_z(double z){
    return int(z*invcellsize.z);
  }
 
  void change_fine_margin(int interval);
};


#endif
