#ifndef SAMPLECELL_H
#define SAMPLECELL_H

#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <iostream>

#include <cassert>

#include "Common.h"
#include "ParticleInfo.h"

bool
power_of_two(int n, int &p);

void
split_power_of_two(int p, int sp[3]);

size_t
generate_cell_index(int celldiv, 
                    std::vector< SpaceVector<int> >& cell_position);

void
distribute_cell(int div[3], int celldiv, 
                std::vector< SpaceVector<int> >& cell_position,
                std::vector< std::vector<int> >& cell_index_node,
                std::map<int,int>& cell_to_node);

double
min_cell_distance(SpaceVector<int> ci, SpaceVector<int> cj, 
                  SpaceVector<double> scale);

enum ReqShift{
  PlusX,
  MinusX,
  PlusY,
  MinusY,
  PlusZ,
  MinusZ
};

struct ShiftCellArray{
  std::vector<int> plusx;
  std::vector<int> minusx;
  std::vector<int> plusy;
  std::vector<int> minusy;
  std::vector<int> plusz;
  std::vector<int> minusz;
};

size_t
cubic_target(int celldiv, SpaceVector<int> ci, int reach,
             std::vector<int>& target_cell, 
             ShiftCellArray& shift_cell);

void
merge_shift_cell(std::vector<ShiftCellArray>& shift_cell,
                 ShiftCellArray& shift_cell_all);

int
calc_target_gather(std::map<int,int>& cell_to_node,
                   int mynode, 
                   std::vector< std::vector<int> >& target_cell,
                   std::vector<int>& target_node,
                   std::vector< std::vector<int> >& target_dist
                   );

int
calc_target_scatter(std::map<int,int>& cell_to_node,
                    int mynode, std::vector<int> mycell,
                    std::vector< std::vector<int> >& target_cell,
                    std::vector<int>& target_node,
                    std::vector< std::vector<int> >& target_dist
                    );


class PostProcess {
public:
  SpaceVector<double> shift;
  ShiftCellArray shift_cell;

  PostProcess() : shift(), shift_cell() {}
  PostProcess(SpaceVector<double> sft, ShiftCellArray& sc) 
    : shift(sft), shift_cell(sc)
  {
  }
// shift position of particles listed in shift_cell
  inline void
  postprocess_receive(ParticleArray& particle, 
                      std::vector<int>& setid,
                      std::map<int,int>& setidtoindex,
                      std::vector<TypeRange>& typerange)
  {
    std::cout << " postprocess_receive" << std::endl;
    std::cout << " setid.size() " << setid.size() << std::endl;
    std::cout << " typerange.size() " << typerange.size()  << std::endl;


    std::cout << " shift_cell.*.size()" 
              << " " << shift_cell.plusx.size()
              << " " << shift_cell.minusx.size()
              << " " << shift_cell.plusy.size()
              << " " << shift_cell.minusy.size()
              << " " << shift_cell.plusz.size()
              << " " << shift_cell.minusz.size()
              << std::endl;
    for(std::vector<int>::iterator it = shift_cell.plusx.begin();
        it!=shift_cell.plusx.end(); ++it){
      int set = setidtoindex[*it];
      for(int i=typerange[set].begin;i<typerange[set].end;i++){
        particle[i].position.x += shift.x;
      }
    }
    for(std::vector<int>::iterator it = shift_cell.minusx.begin();
        it!=shift_cell.minusx.end(); ++it){
      int set = setidtoindex[*it];
      for(int i=typerange[set].begin;i<typerange[set].end;i++){
        particle[i].position.x -= shift.x;
      }
    }

    for(std::vector<int>::iterator it = shift_cell.plusy.begin();
        it!=shift_cell.plusy.end(); ++it){
      int set = setidtoindex[*it];
      for(int i=typerange[set].begin;i<typerange[set].end;i++){
        particle[i].position.y += shift.y;
      }
    }
    for(std::vector<int>::iterator it = shift_cell.minusy.begin();
        it!=shift_cell.minusy.end(); ++it){
      int set = setidtoindex[*it];
      for(int i=typerange[set].begin;i<typerange[set].end;i++){
        particle[i].position.y -= shift.y;
      }
    }

    for(std::vector<int>::iterator it = shift_cell.plusz.begin();
        it!=shift_cell.plusz.end(); ++it){
      int set = setidtoindex[*it];
      for(int i=typerange[set].begin;i<typerange[set].end;i++){
        particle[i].position.z += shift.z;
      }
    }
    for(std::vector<int>::iterator it = shift_cell.minusz.begin();
        it!=shift_cell.minusz.end(); ++it){
      int set = setidtoindex[*it];
      for(int i=typerange[set].begin;i<typerange[set].end;i++){
        particle[i].position.z -= shift.z;
      }
    }
  }
};
#endif
