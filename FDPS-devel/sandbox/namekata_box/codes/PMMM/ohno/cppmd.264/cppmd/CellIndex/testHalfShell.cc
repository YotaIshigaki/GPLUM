#include <iostream>
#include <cstring>

#include "HalfShell.h"

int
main(int argc, char **argv)
{
  int celldiv = 16;
  SpaceVector<double> cellsize(4.0,4.0,4.0);
  SpaceVector<double> margin(1.0,1.0,1.0);
  double cutoff = 9.0;

  std::vector< SpaceVector<int> > cellindex;
  size_t num_cell;
  std::vector< std::vector<int> > send_target_cell;
  std::vector< std::vector<int> > recv_target_cell;
  std::vector< ShiftCellArray >  shift_cell;

  int np=1;

  int opt;
  while ((opt = getopt(argc,argv,"c:C:p:"))!=-1){
    switch (opt){
    case 'c':
      celldiv = atol(optarg);
      break;
    case 'C':
      cutoff = strtod(optarg,(char **)NULL);
      break;
    case 'p':
      np = atol(optarg);
      break;
    }
  }
  if(np<1)np=1;
  if(np>celldiv*celldiv*celldiv)np=celldiv*celldiv*celldiv;

  int p;
  if(!power_of_two(np, p)){
    std::cout << "np != 2^m " << std::endl;
  }
  int sp[3];
  split_power_of_two(p,sp);
  int div[3] = {(1<<sp[0]),(1<<sp[1]),(1<<sp[2])};
  if(celldiv<(1<<sp[2])){
    celldiv = 1<<sp[2];
    std::cout << "celldiv recalc" << std::endl;
  }

  num_cell = generate_cell_index(celldiv,cellindex);

  std::cout << "np " << np << "  celldiv " << celldiv << "  num_cell " << num_cell  << std::endl;
  
  std::vector< std::vector<int> > cell_index_node;
  std::map<int,int> cell_to_node; 
  
  distribute_cell(div, celldiv, cellindex, cell_index_node,cell_to_node);
  
  for(int n=0;n<np;n++){
    std::cout << "node " << n << std::endl;
    int num_cell_node = cell_index_node[n].size();
    recv_target_cell.resize(num_cell_node);
    send_target_cell.resize(num_cell_node);
    shift_cell.resize(num_cell_node);
    for(size_t index=0;index<num_cell_node;index++){
      halfshell_target(celldiv, cellindex[cell_index_node[n][index]], 
                       cellsize, margin, cutoff,
                       send_target_cell[index],
                       recv_target_cell[index], shift_cell[index]);
    }
    for(size_t index=0;index<num_cell_node;index++){
      dump_target_cell(cellindex,cell_index_node[n][index],
                       send_target_cell[index]);
    }
    for(size_t index=0;index<num_cell_node;index++){
      dump_target_cell(cellindex,cell_index_node[n][index],
                       recv_target_cell[index]);
    }
    ShiftCellArray shift_cell_all;
    merge_shift_cell(shift_cell, shift_cell_all);
    
    std::vector<int> send_target_id;
    std::vector< std::vector<int> > send_setid;
    size_t num_send_target_node
      = calc_target_scatter(cell_to_node, n, 
                            cell_index_node[n], send_target_cell, 
                            send_target_id, 
                            send_setid);
    for(size_t t=0;t<send_target_id.size();t++){
      std::cout << "to target " << send_target_id[t] << " send";
      dump_cell(cellindex,send_setid[t]);
    }

    std::vector<int> recv_target_id;
    std::vector< std::vector<int> > recv_setid;
    size_t num_recv_target_node 
      = calc_target_gather(cell_to_node, n, recv_target_cell,
                           recv_target_id, recv_setid);
    for(size_t t=0;t<recv_target_id.size();t++){
      std::cout << "from target " << recv_target_id[t] << " receive";
      dump_cell(cellindex,recv_setid[t]);
    }
  }

}
