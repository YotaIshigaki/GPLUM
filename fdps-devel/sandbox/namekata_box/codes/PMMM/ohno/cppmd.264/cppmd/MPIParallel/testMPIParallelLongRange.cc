#include <mpi.h>
#include "Common.h"
#include "CalcPreparator.h"
#include "MPIParallelLongRange.h"

int
main(int argc, char **argv)
{
  int node_id;
  int num_node;
  int short_id;
  int num_short;
  int long_id;
  int num_long_only;
  int num_long;

  MPI_Comm mpi_comm_all = MPI_COMM_WORLD;
  MPI_Comm mpi_comm_short;
  MPI_Comm mpi_comm_long;

  LongRangeMPIPlan longplan = Combine0D;

  OperationSelector ope;


  MPI_Init(&argc, &argv);

  num_long_only = 1; 

  if(argc>1) longplan = LongRangeMPIPlan(atoi(argv[1]));
  if(argc>2) num_long_only = atoi(argv[2]);

  MPI_Comm_rank(mpi_comm_all,&node_id);
  MPI_Comm_size(mpi_comm_all,&num_node);

  num_short = num_node - num_long_only;


  SpaceVector<double> cellmargin(1.0, 1.0, 1.0);
  CellIndexType citype=HalfShell;
  bool withlong = true;
  bool withbond = false;
  CalcPreparator_::makeoperations(node_id, 
                                  withlong, 
                                  withbond, 
                                  citype,
                                  longplan,
                                  num_long_only);
  int num_particle = num_short*100;
  int num_set;
  int celldiv=4;
  int max_per_cell = 150;
  SpaceVector<double> boxsize(64.0,64.0,64.0);
  double cutoff = 9.0;
  std::vector<int> node_id_of_shorts;
  num_set = celldiv*celldiv*celldiv;
  std::cout << "call CalcPreparator_::maketarget" << std::endl;
  CalcPreparator_::maketarget(node_id, short_id, node_id_of_shorts,
                              num_node, num_particle, num_set,
                              celldiv, max_per_cell, boxsize,
                              cellmargin, citype, cutoff,
                              longplan, num_long_only);
  ope.doLongrangecalculation = true;
  ope.doShortrangecalculation = true;

  int color = 0;
  if(node_id<num_short){
    color=1;
    ope.doShortrangecalculation = true;
  }else{
    color = MPI_UNDEFINED;
    ope.doShortrangecalculation = false;
  }

  std::map<int,int> shortidtorank;
  for(int i=0;i<num_short;i++)shortidtorank.insert(std::pair<int,int>(i,i));

  std::cout << "node_id " << node_id << " short_id " << short_id << " " << ope.doShortrangecalculation << std::endl;


  GeometryXYZ longgeometry;

  CalcPreparator_::makelongrangegeometry(longplan, num_long_only, 
                                         num_long, longgeometry);

  CalcPreparator_::makelongcomm(longplan,longgeometry,num_long,num_short,
                                MPI_COMM_WORLD, node_id, short_id, long_id);

  /*
  MPIParallelLongRange::assign_long(node_id, short_id, num_short, 
                                    longplan, mpi_comm_all,
                                    num_long, long_id, mpi_comm_long, ope);
  */

  std::cout << "num_node " << num_node << " num_short " << num_short << " num_long " << num_long << std::endl;
  std::cout << "node_id " << node_id << " long_id " << long_id <<  " " << ope.doLongrangecalculation << std::endl;

  std::vector<int> long_reqcell_list(0);
  double long_fringe = 5.0;

  CalcPreparator_::makelongrangerequiredcell(longgeometry,long_fringe,
                                             long_id,
                                             long_reqcell_list);

  std::cout << " receive " << long_reqcell_list.size() << " cell for long";
  for(int i=0;i<long_reqcell_list.size();i++){
    std::cout << " " << long_reqcell_list[i];
  }
  std::cout << std::endl;

  std::vector< std::vector<int> > send_to_long_list;

  CalcPreparator_::makecellidlistsendtolong(longgeometry,long_fringe,num_long,
                                            short_id,
                                            send_to_long_list);

  std::cout << " send " << send_to_long_list.size() << " long node " << std::endl;
  for(int lid=0;lid<send_to_long_list.size();lid++){
    std::cout << " to " << lid;
    for(int i=0;i<send_to_long_list[lid].size();i++){
      std::cout << " " << send_to_long_list[lid][i];
    }
    std::cout << std::endl;
  }

  std::vector<int> long_short_id;
  std::vector<int> sender_global_id_list;
  std::vector<int> reciever_global_id_list;
  std::vector<MPI_Comm> mpi_comm_ls_list;


  CalcPreparator_::make_long_short_comm(num_node, num_long, 
                                        node_id, short_id, long_id,
                                        MPI_COMM_WORLD,
                                        send_to_long_list,
                                        long_short_id,
                                        sender_global_id_list,
                                        reciever_global_id_list,
                                        mpi_comm_ls_list);
                                             
  for(int lid=0;lid<num_long;lid++){
    if(lid==long_id){
      std::cout << "recv from";
      for(int s=0;s<sender_global_id_list.size();s++){
        std::cout << " " << sender_global_id_list[s];
      }
      std::cout << std::endl;
    }
    if((short_id!=NO_SHORT_ID)&&(mpi_comm_ls_list[lid]!=MPI_COMM_NULL)){
      std::cout << "send to " << reciever_global_id_list[lid] << std::endl;
    }
  }


  MPI_Finalize();

  exit (0);
}
