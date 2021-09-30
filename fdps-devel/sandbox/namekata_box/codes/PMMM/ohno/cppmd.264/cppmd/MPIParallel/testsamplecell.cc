#include <iterator>
#include "samplecell.h"

int
main(int argc, char **argv)
{
  int n=1;
  int celldiv=4;    // number of cell in axis
  int reach=1;      // range of neighbor cell

  if(argc>1){
    n=std::atol(argv[1]);
    if(argc>2){
      int cp = std::atol(argv[2]);
      if(cp<1)cp=1;
      celldiv = 1<<cp;
      if(argc>3){
        reach = std::atol(argv[3]);
        if(reach<1)reach=1;
        if(reach>=(celldiv>>1))reach=(celldiv>>1)-1;
      }
    }
  }

  int np;
  if(power_of_two(n,np)){
    std::cout << n << " = 2^" << np << std::endl;
  }else{
    std::cout << n << " != 2^" << np << std::endl;
    exit(1);
  }

  int sp[3];
  split_power_of_two(np,sp);
  int div[3];         // number of node in x,y,z axis
  div[0] = 1<<sp[0];
  div[1] = 1<<sp[1];
  div[2] = 1<<sp[2];
  int num_node = div[0]*div[1]*div[2];
  std::cout << num_node << " = " << div[0] << "*" << div[1] << "*" << div[2] << std::endl;
  
  if(celldiv<div[2]){
    std::cout << "cell div < " << div[2] << std::endl;
    celldiv=div[2];
  }
  
  std::vector< SpaceVector<int> > cell_position;
  size_t num_cell = generate_cell_index(celldiv,cell_position);

  std::cout << "number of cell " << celldiv << "^3 " << num_cell << std::endl;

  std::vector< std::vector<int> > cell_index_node;
  std::map<int,int> cell_to_node;
  distribute_cell(div, celldiv, cell_position, cell_index_node,cell_to_node);
  std::cout << "cell_to_node";
  for(int i=0;i<cell_position.size();i++){
    std::cout << " " << cell_to_node[i];
  }
  std::cout << std::endl;
  
  int mynode=0;


  std::vector<int> mycell = cell_index_node[mynode];
  size_t number_of_mycell = mycell.size();
  std::vector< std::vector<int> > target_cell(number_of_mycell);
  std::vector< ShiftCellArray > shift_cell(number_of_mycell);
   for(int i=0;i<number_of_mycell;i++){
    size_t numtc = cubic_target(celldiv,cell_position[mycell[i]],reach,
                                target_cell[i],shift_cell[i]);
  }
  
  std::vector<int> target_node;
  std::vector< std::vector<int> > target_dist;
  size_t num_target_node 
    = calc_target_gather(cell_to_node, mynode, target_cell,
                         target_node, target_dist);
  std::cout << "num_target_node " << num_target_node << std::endl;
  for(int tn=0;tn<num_target_node;tn++){
    std::cout << "target node " << target_node[tn] << "  cell ";
    for(int tc=0;tc<target_dist[tn].size();tc++){
      std::cout << " " << target_dist[tn][tc];
    }
    std::cout << std::endl;
  }
  
  std::ostream_iterator< int > ostream_iter( std::cout, " " );

  ShiftCellArray shift_cell_all;
  merge_shift_cell(shift_cell, shift_cell_all);

  std::cout << "plusx" << std::endl;
  std::copy(shift_cell_all.plusx.begin(), shift_cell_all.plusx.end(),ostream_iter);
  std::cout << std::endl;
  std::cout << "minusx" << std::endl;
  std::copy(shift_cell_all.minusx.begin(), shift_cell_all.minusx.end(),ostream_iter);
  std::cout << std::endl;
  std::cout << "plusy" << std::endl;
  std::copy(shift_cell_all.plusy.begin(), shift_cell_all.plusy.end(),ostream_iter);
  std::cout << std::endl;
  std::cout << "minusy" << std::endl;
  std::copy(shift_cell_all.minusy.begin(), shift_cell_all.minusy.end(),ostream_iter);
  std::cout << std::endl;
  std::cout << "plusz" << std::endl;
  std::copy(shift_cell_all.plusz.begin(), shift_cell_all.plusz.end(),ostream_iter);
  std::cout << std::endl;
  std::cout << "minusz" << std::endl;
  std::copy(shift_cell_all.minusz.begin(), shift_cell_all.minusz.end(),ostream_iter);
  std::cout << std::endl;

}
