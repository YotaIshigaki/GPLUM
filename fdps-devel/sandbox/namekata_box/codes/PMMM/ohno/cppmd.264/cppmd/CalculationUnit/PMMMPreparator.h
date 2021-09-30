#include "Common.h"
#include "PMMMInteraction.h"

namespace PMMMPreparator{

  void estimate_alpha(double &alpha, const SpaceVector<double> &box, const double scaled_alpha=2.4)
  {
    double len = std::min(std::min(box.x,box.y),box.z);
    alpha = scaled_alpha/len;
  }

  void make_mm_info(const int order, const int num_pp, const int num_mm,
		    std::vector<int> &m_boundary, 
		    std::vector<int> &pp_targetid, 
		    std::vector<int> &mm_targetid)
  {
    mm_decompose(m_boundary,num_mm,order);
    mm_targetid.resize(num_mm);
    for(int m=0;m<num_mm;m++)mm_targetid[m]=num_pp+m;
    pp_targetid.resize(num_pp);
    for(int m=0;m<num_pp;m++)pp_targetid[m]=m;
  }

  void set_cell_range(const MPI_Comm &comm, const MPI_Comm &pp_comm, const int num_pp, 
		      const MPI_Comm &mm_comm,
		      const CellIndex &cellindex, const CellSubset &cellsubset, const OperationSelector &operation,
		      std::vector<Position> &cell_center,
		      std::vector<CellMethodModule::CellRange> &cell_range)
  {
    int node_id;
    MPI_Comm_rank(comm,&node_id);
    cell_range.resize(num_pp);
    if(operation.doShortrangecalculation){
      int pp_id;
      MPI_Comm_rank(pp_comm,&pp_id);

      CellMethodModule::CellRange my_range;

      my_range.min = cellindex.celldiv;
      my_range.max.x = my_range.max.y = my_range.max.z = 0;
      Position cs = cellindex.cellsize;
      cell_center.resize(cellsubset.num_mycell);
      for(int i=0;i<cellsubset.num_mycell;i++){
	SpaceVector<int> pos = cellindex.cell_position[cellsubset.mycell[i]];
	if(pos.x>my_range.max.x)my_range.max.x=pos.x;
	if(pos.y>my_range.max.y)my_range.max.y=pos.y;
	if(pos.z>my_range.max.z)my_range.max.z=pos.z;
	if(pos.x<my_range.min.x)my_range.min.x=pos.x;
	if(pos.y<my_range.min.y)my_range.min.y=pos.y;
	if(pos.z<my_range.min.z)my_range.min.z=pos.z;
	cell_center[i].x = cs.x*(0.5+pos.x);
	cell_center[i].y = cs.y*(0.5+pos.y);
	cell_center[i].z = cs.z*(0.5+pos.z);
      }
      my_range.max.x++;
      my_range.max.y++;
      my_range.max.z++;
      if(DebugLog::verbose>1)printf("Gather cell_range %d\n",node_id);
      MPI_Gather(&(my_range.min[0]), 6, MPI_INTEGER,
		 &(cell_range[0].min[0]), 6, MPI_INTEGER, 0, pp_comm);
      if(pp_id==0){
	if(DebugLog::verbose>1)printf("send cell_range to %d : by %d\n",num_pp,node_id);
	MPI_Send(&(cell_range[0].min[0]), 6*num_pp, MPI_INTEGER, num_pp, 0, comm);
      }
    }
    if(operation.doLongrangecalculation){
      int mm_id;
      MPI_Comm_rank(mm_comm,&mm_id);
      if(mm_id==0){
	MPI_Status stat;
	if(DebugLog::verbose>1)printf("recv cell_range to %d : by %d\n",0, node_id);
	MPI_Recv(&(cell_range[0].min[0]), 6*num_pp,  MPI_INTEGER, 0, 0, comm, &stat);
      }
      if(DebugLog::verbose>1)printf("Bcast cell_range %d\n",node_id);
      MPI_Bcast(&(cell_range[0].min[0]), 6*num_pp,  MPI_INTEGER, 0, mm_comm);
    }

  }

  void initialize_PMMM_Interaction(const MPI_Comm &comm, const MPI_Comm &pp_comm, const MPI_Comm &mm_comm,
				   const CellIndex &cellindex,
				   const LongRangeParameter &longparam,
				   const OperationSelector &operation,
				   const std::vector<int> &m_boundary, 
				   const std::vector<int> &pp_targetid,
				   const std::vector<int> &mm_targetid, 
				   const std::vector<Position> &cell_center,
				   const std::vector<CellMethodModule::CellRange> &cell_range,
				   const double total_charge,
				   PMMMLocalInteraction &pmmm_local,
				   PMMMLongRangeInteraction &pmmm_long)
  {
    int node_id;
    MPI_Comm_rank(comm,&node_id);
    if((node_id==0)||(DebugLog::verbose>1)){
      printf("initialize_PMMM_Interaction\n");
    }
    const int celldiv_node[3] = {cellindex.celldiv_node.x,cellindex.celldiv_node.y,cellindex.celldiv_node.z};
    const int celldiv[3] = {cellindex.celldiv.x,cellindex.celldiv.y,cellindex.celldiv.z};
    if(operation.doShortrangecalculation){
      pmmm_local.settings(node_id, longparam, pp_comm);
      if(total_charge!=0.0){
	pmmm_local.not_charge_neutral(total_charge);
	if((node_id==0)||(DebugLog::verbose>1)){
	  printf("not_charge_neutral(%e)\n",total_charge);
	}
      }
      pmmm_local.initialize(cell_center, celldiv_node, 
			    comm, mm_targetid, m_boundary);
    }
    if(operation.doLongrangecalculation){
      pmmm_long.settings(node_id, longparam, mm_comm);
      int icut = std::min(std::min(cellindex.target_range.x,cellindex.target_range.y),cellindex.target_range.z);
      int mmid;
      MPI_Comm_rank(mm_comm,&mmid);
      if((mmid==0)||(DebugLog::verbose>1)){
	printf("ICUT %d\n",icut);
      }
      pmmm_long.initialize(celldiv, m_boundary, cellindex.cellsize.x,
			   3, 5, icut,
			   comm, pp_targetid, cell_range);
      if(DebugLog::verbose>1)printf("local mm %d %d\n",pmmm_long.pmmm_mm.m_begin,pmmm_long.pmmm_mm.m_end);
    }
     
  }

}
