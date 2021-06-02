/*
  move setLJCutoffEnergyCorrectionParmeter into makeparticle
*/

// longgeometry return by  makelongrangegeometry.
// It must user for other make***

#include "CalcPreparator.h"


// added by Gleb on 2011.05.23
#ifdef CPPMD_ENABLE_LONGRANGE
#include <LongRangeInteraction/EwaldBase.h>
#endif
//

#ifdef CPPMD_ENABLE_PMMM
#include "PMMMPreparator.h"
#endif

/*! TODO 
  supress duplicate ghost short-short short-long when Combind long
 */

/*
  INPUT for all
  node_id

  INPUT for makeoperations only
  bool withlong
  bool withbond

  INPUT for makeoperations and maketarget
  CellIndexType cellindextype

  INPUT for maketarget only
  int num_node :
  num_particle
  SpaceVector<double> cellmargin :
  
  INPUT for maketarget and constructcommunicator
  max_particle_in_cell

  INPUT for maketarget and constructcalculationunit
  SpaceVector<double> boxsize :
  cutoff

  INPUT-OUTPUT for maketarget
  num_total_set
  celldiv

  INPUT for makeparticle only
  ParticleArray& pa
  std::vector<PotentialModel>& pm
  //  SpaceVector<double> boxsize unused

  INPUT for makebond only
  CovalentBondList bl
  bool excludewaterbond

  INPUT for constructcommunicator only
  std::map<int,int> uidtorank

  INPUT for constructcalculationunit only
  ShortRange::CoulombType cltype

  makeoperations to maketarget, constructcalculationunit
  ope

  maketarget to makeparticle
  cellsubset

  maketarget to constructcommunicator
  sendparticle_target_id
  sendparticle_setid
  recvparticle_target_id
  recvparticle_setid
  move_target_id
  move_setid
  cellindex
  max_id
  node_div

  maketarget to constructcalculationunit
  bondlistarray
  targetset
  postprocess

  maketarget to makeparticle, constructcommunicator and constructcalculationunit
  particle_setid 
  
  maketarget to makeparticle and constructcalculationunit
  particlearray
  typerangearray

  makeparticle to constructcalculationunit
  int num_all_particle

  constructcommunicator to constructcalculationunit
  MPICommunicator communicator

  maketarget to
  num_set_per_node

  LOCAL for maketarget
  volume
  volume_margin

  LOCAL for makeparticle
  num_lj
*/

namespace CalcPreparator_ {

//! set by makeoperations, used in maketarget and constructcalculationunit
OperationSelector ope;
bool expire_reverse;

MPI_Comm mpi_comm_short;
MPI_Comm mpi_comm_long;

//! set by maketarget, used in makeparticle
CellSubset cellsubset;   //!< cell information

//! set by maketarget, used in makeparticle, constructcommunicator and constructcalculationunit
std::vector<int> particle_setid;

//! set by maketarget, used in makeparticle and constructcalculationunit
#ifdef OLDPARTICLE
ParticleArray particlearray;
#else
CombinedParticleArray particlearray;
#endif
std::vector<TypeRange> typerangearray;

//! set by maketarget, used in makebond
int bond_assign;  //!< 0: 1st atom,  1:2nd atom, 3:middle atom

//! set by maketarget, used in constructcommunicator
std::vector<int> sendparticle_target_id;
std::vector< std::vector<int> > sendparticle_setid;
std::vector<int> recvparticle_target_id;
std::vector< std::vector<int> > recvparticle_setid;
std::vector<int> move_target_id;
std::vector< std::vector<int> > move_setid;
CellIndex cellindex;
int max_id;
SpaceVector<int> node_div;
SpaceVector<int> cell_div3d;
SpaceVector<int> cell_div3d_in_node;
double boxsize_minimum_scale;

//! set by maketarget, used in constructcalculationunit
std::vector<CovalentBondInfo::BondList> bondlistarray;
std::vector<CovalentBondInfo::BondList> bondlistarray_idx;
std::vector< std::vector<int> > targetset;
std::vector< std::pair<int,int> > ghostpairs;
PostProcess postprocess;

int num_set_per_node;

//! set by makeparticle, used in makebond
std::map<int,int> particle_index_to_setindex;
std::map<AtomID,int> atomid_to_index;

//! set by makeparticle, used in constructcalculationunit
int num_all_particle;
WaterList waterlist;
ShakeList shakelist;

int num_freedom;

std::vector<int> longid_to_worldid; //!< ranks at world of long node

//! set by makelongrangeparameter, used in constructcalculationunit
/*!
  TODO : used in construct function for longrange communicator (constructcommunicator or new function).
 */
LongRangeParameter longrangeparameter;

//! set by constructcommunicator, used in constructcalculationunit
MPICommunicator communicator;

MPIReceiverForLong longreceiver;

MPISenderForLong longsender;

//! Todo: force longplan None when node number is invalide.
void makeoperations(const int node_id,
                    const bool withlong,
                    const bool withbond,
		    const bool withcorrecttrans,
                    const CellIndexType cellindextype,
                    const LongRangeMPIPlan longplan,
                    const int num_long_only)
{
  ParallelCoordinator parallelcoordinator;
  ope = parallelcoordinator.selectOperations(0,withlong,withbond);

  if(withlong==true){
    switch (longplan) {
      case Combine0D:
        if(node_id==0){
          ope.doLongrangecalculation = true;
        }else{
          ope.doLongrangecalculation = false;
        }
        ope.doShortLongCommunication = true;
        break;
      case Combine1D:
        // depend on node position, set in/after maketarget
        //ope.doLongrangecalculation = true;
        ope.doShortLongCommunication = true;
        break;
      case Combine2D:
        // depend on node position, set in/after maketarget
        //ope.doLongrangecalculation = true;
        ope.doShortLongCommunication = true;
        break;
      case Combine3D:
        ope.doLongrangecalculation = true;
        ope.doShortLongCommunication = true;
        break;
      case Separate0D:
      case Separate1D:
      case Separate2D:
      case Separate3D:
        ope.doShortLongCommunication = true;
        if(num_long_only<=0){
          std::cout << "Separate Long Range allowd only with positive number of long range only node" << std::endl;
          ope.doLongrangecalculation = false;
          ope.doShortLongCommunication = false;
        }
        // depend on node position, set in/after maketarget
        break;
      default:
        std::cout << "unknown long rage plan. long range disabled" << std::endl;
      case NoneLong:
        ope.doLongrangecalculation = false;
        ope.doShortLongCommunication = false;
        break;
    }
  }else{
    ope.doLongrangecalculation = false;
    ope.doShortLongCommunication = false;
  }
#ifdef CPPMD_ENABLE_PMMM
  ope.doShortLongCommunication = false;
#endif

  ope.doCorrectTranslation = withcorrecttrans;

  if(cellindextype==FullCube){
    ope.doExchangeForce = ope.doCovalentBondcaculation;
    ope.doExchangeForceIndexed = true;
    ope.doReturnShortForce = false;
    ope.doShortEnergyToHalf = true;
#ifdef USE_PAIRLIST
    expire_reverse = false;
#else    
    expire_reverse = true;
#endif
  }else{
    ope.doExchangeForce = true;
    ope.doExchangeForceIndexed = false;
    ope.doShortEnergyToHalf = false;
    ope.doReturnShortForce = true;
    expire_reverse = false;
  }

  if(node_id==0){
    std::cout << "ope.doIntegration " << ope.doIntegration << std::endl;
    std::cout << "ope.doShortrangecalculation " << ope.doShortrangecalculation << std::endl;
    std::cout << "ope.doLongrangecalculation " << ope.doLongrangecalculation << std::endl;
    std::cout << "ope.doCovalentBondcaculation " << ope.doCovalentBondcaculation << std::endl;
    std::cout << "ope.doExchangeForce " << ope.doExchangeForce << std::endl;
    std::cout << "ope.doExchangeForceIndexed " << ope.doExchangeForceIndexed << std::endl;
    std::cout << "ope.doShortEnergyToHalf " << ope.doShortEnergyToHalf << std::endl;
    std::cout << "ope.doShortLongCommunication " << ope.doShortLongCommunication << std::endl;
    std::cout << "ope.doReturnShortForce " << ope.doReturnShortForce << std::endl;
    std::cout << "ope.doMakePairList " << ope.doMakePairList << std::endl;
    std::cout << "ope.doEnergycalculation " << ope.doEnergycalculation << std::endl;
    std::cout << "ope.doVirialcalculation " << ope.doVirialcalculation << std::endl;
  }
}

void makempicomm(MPI_Comm all_comm, int isShortComm, int shortrank, 
                 MPI_Comm *short_comm)
{
  MPI_Comm_split(all_comm,isShortComm,shortrank,short_comm);
}

void makegeometry(const int node_id,
                  const int opts_nodediv3d[],
                  const int opts_celldiv3d_in_node[],
                  const int num_node,
                  const int num_short,
                  int &short_id,
                  std::vector<int>& node_id_of_shorts,
                  int& num_total_set,
                  int& celldiv,
                  SpaceVector<int>& ndiv,
                  SpaceVector<int>& cdiv_in_node,
                  SpaceVector<int>& cdiv)
{
//  int num_short;
//  num_short = num_node-num_long_only;
//  int num_set_per_node = num_total_set/num_short;

  // set node_div
  SpaceVector<int> div;
  if(opts_nodediv3d[0]==0){
    int p;
    if(!power_of_two(num_short, p)){
      //    std::cout << "num_short != 2^m " << std::endl;
      int m[3];
      if(power_of_two_of_cube(num_short,m)){
        div[0] = m[0];
        div[1] = m[1];
        div[2] = m[2];
      }
    }else{
      int sp[3];
      split_power_of_two(p,sp);
      div[0] = 1<<sp[0];
      div[1] = 1<<sp[1];
      div[2] = 1<<sp[2];
    }
    node_div = div;
  }else{
    node_div[0] = opts_nodediv3d[0];
    node_div[1] = opts_nodediv3d[1];
    node_div[2] = opts_nodediv3d[2];
  }
  if(node_id==0){
    std::cout << "num_short = " << node_div[0] << " * " << node_div[1] << " * " << node_div[2] << std::endl;
  }
  ndiv = node_div;

  // set celldiv
  if(opts_celldiv3d_in_node[0]!=0){
    for(int i=0;i<3;i++){
      cell_div3d_in_node[i] = opts_celldiv3d_in_node[i];
      cell_div3d[i] = opts_celldiv3d_in_node[i] * node_div[i];
    }
  }else{
    int min_div = node_div[0];
    if(node_div[1] < min_div) min_div=node_div[1];
    if(node_div[2] < min_div) min_div=node_div[2];

    if(celldiv<min_div){
      celldiv = min_div;
    }else{
      if(celldiv%min_div!=0){
        celldiv = (celldiv/min_div)*min_div;
      }
    }
    for(int i=0;i<3;i++){
      cell_div3d[i] = celldiv;
      if(celldiv%node_div[i] !=0){
        printf("celldiv %% node_div[%d] != 0\n",i);
        exit(1);
      }else{
        cell_div3d_in_node[i] = celldiv / node_div[i];
      }
    }
  }

  num_total_set = cell_div3d[0] * cell_div3d[1] * cell_div3d[2];
  if(node_id==0){
    std::cout << "cell_div = " << cell_div3d[0] << " * " << cell_div3d[1] << " * " << cell_div3d[2] <<  std::endl;
    std::cout << "cell_div_in_node = " << cell_div3d_in_node[0] << " * " << cell_div3d_in_node[1] << " * " << cell_div3d_in_node[2] <<  std::endl;
    std::cout << "num_total_set = " << num_total_set << std::endl;
  }
  cdiv_in_node = cell_div3d_in_node;
  cdiv = cell_div3d;

  int isShortComm = 0;
  if(node_id<num_short){
    isShortComm = 1;
    if((node_id==0)||(DebugLog::verbose>1)) {
      printf("node %d is short\n",node_id);
    }
  }else{
    isShortComm = MPI_UNDEFINED;
    ope.doShortrangecalculation = false;
    ope.doCovalentBondcaculation = false;
    if((node_id==0)||(DebugLog::verbose>1)) {
      printf("node %d is not short\n",node_id);
    }
  }
  makempicomm(MPI_COMM_WORLD, isShortComm, node_id, 
              &mpi_comm_short);
  int send;
  int *recv = new int[num_short];
  if(mpi_comm_short!=MPI_COMM_NULL){
    MPI_Comm_rank(mpi_comm_short,&short_id);
    send = short_id;
    MPI_Allgather(&send,1,MPI_INT,recv,1,MPI_INT,mpi_comm_short);
  // rank 0 of MPI_COMM_WORLD must be in mpi_comm_short
    if(recv[0]!=0){
      std::cout << " short0 != node0 " << std::endl;
      if(node_id==0){
        MPI_Status stat;
        MPI_Recv(recv,num_short,MPI_INT,recv[0],0,MPI_COMM_WORLD,&stat);
      }
      if(short_id==0){
        MPI_Send(recv,num_short,MPI_INT,0,0,MPI_COMM_WORLD);
      }
    }
  }else{
    short_id = NO_SHORT_ID;
  }
  MPI_Bcast(recv,num_short,MPI_INT,0,MPI_COMM_WORLD);
  node_id_of_shorts.resize(num_short);
  for(int i=0;i<num_short;i++)node_id_of_shorts[i] = recv[i];

}

void maketarget(const int node_id,
                const int short_id,
                const int num_short,
                const int num_particle,
                const int max_particle_in_cell,
                const SpaceVector<double> boxsize, 
                const SpaceVector<double> cellmargin,
                const CellIndexType cellindextype,
                const double cutoff,
		const double bms,
                int& num_total_set)
{
  boxsize_minimum_scale = bms;
  cellindex.set_cell_parameter(boxsize,
                               cell_div3d,
                               cellmargin,
                               cutoff);
  if(node_id==0){
    std::cout << "cell length " << cellindex.cellsize << std::endl;
  }
#if 0
  volume = cellindex.calc_volume();
#else
  cellindex.calc_volume();
#endif
  int num_cell = cellindex.generate_cell_index();
  postprocess.celldiv = cellindex.celldiv;
  postprocess.cell_position = cellindex.cell_position;

  num_set_per_node = num_total_set/num_short;
  if(num_cell!=num_set_per_node*num_short){
    if(node_id==0){
      std::cout << "num_cell "<< num_cell << " != num_set_per_node*num_short " 
                << num_set_per_node << "*" << num_short << " " <<  num_set_per_node*num_short << std::endl;
    }
    num_set_per_node = num_cell/num_short;
    num_total_set = num_cell;
  }

//  if(!cellindex.distribute_cell(div)){
  if(!cellindex.distribute_cell(node_div)){
    std::cout << "distribute_cell false " << std::endl;
  }

  max_id = num_total_set;
  int max_num_particle = max_particle_in_cell*num_total_set;
  int num_spu = num_set_per_node;
  //    int max_num_pps = max_particle_in_cell; // depend on simulation
  particlearray.resize(max_particle_in_cell*num_set_per_node);
  if ( (node_id==0) || (DebugLog::verbose>1) ) {
    //    std::cout << "particle " << num_particle << std::endl;
    std::cout << "max particle " << max_num_particle << std::endl;
    std::cout << "max particle in cell " << max_particle_in_cell << std::endl;
    std::cout << "num_set_per_node " << num_set_per_node << std::endl;
    std::cout << "reserved particlearray " << max_particle_in_cell*num_set_per_node << std::endl;
    std::cout << "node id " << node_id <<  " short id " << short_id << std::endl;
  }
  
  bondlistarray.clear();
  bondlistarray_idx.clear();


  if(short_id!=NO_SHORT_ID)
  {
    particle_setid.resize(num_spu);
    typerangearray.resize(num_spu);  
  // assign cell to unit and reserve particle in unit
    for(size_t s=0;s<particle_setid.size();s++){
      particle_setid[s] = cellindex.cell_id_in_node[short_id][s];
      typerangearray[s].begin = s*max_particle_in_cell;
      typerangearray[s].end = typerangearray[s].begin;
      typerangearray[s].lj.begin = typerangearray[s].begin;
      typerangearray[s].lj.end = typerangearray[s].end;
      typerangearray[s].ljcoulomb.begin = typerangearray[s].lj.end;
      typerangearray[s].ljcoulomb.end = typerangearray[s].lj.end;
      typerangearray[s].coulomb.begin = typerangearray[s].lj.end;
      typerangearray[s].coulomb.end = typerangearray[s].lj.end;
    }
  }else{
    particle_setid.resize(0);
    typerangearray.resize(0);
  }

  // set target set
  postprocess.max_particle_in_cell = max_particle_in_cell;
  postprocess.margin = cellindex.margin;
  postprocess.fine_margin = cellindex.margin;
  postprocess.boxsize = cellindex.boxsize;
  postprocess.cellsize = cellindex.cellsize;
  postprocess.invcellsize.x = 1.0/postprocess.cellsize.x;
  postprocess.invcellsize.y = 1.0/postprocess.cellsize.y;
  postprocess.invcellsize.z = 1.0/postprocess.cellsize.z;
#if 0
  volume_margin = postprocess.calc_volume_margin();
#else
  postprocess.calc_volume_margin();
#endif

  if(short_id!=NO_SHORT_ID){
    if(boxsize_minimum_scale<1.0){
      SpaceVector<double>min_box(boxsize*boxsize_minimum_scale);
      if(short_id==0){
	std::cout << "select target cell by minimum boxsize" << min_box << std::endl;
      }
      cellindex.change_boxsize(min_box);
    }
    cellsubset = CellSubset(cellindex,short_id);
    cellsubset.makefulltarget(particle_setid,cellindextype,ope);

    postprocess.shift_cell=cellsubset.shift_cell_all;

    cellsubset.makecommtarget(sendparticle_target_id, 
                              sendparticle_setid,
                              recvparticle_target_id,
                              recvparticle_setid,
                              move_target_id, move_setid);
    if(boxsize_minimum_scale<1.0){
      cellindex.change_boxsize(boxsize);
      cellsubset.cellindex = cellindex;
    }
    targetset = cellsubset.short_target_cell;
    ghostpairs = cellsubset.ghost_cell_pairs;
    bondlistarray.resize(num_set_per_node);
    bondlistarray_idx.resize(num_set_per_node);
    bond_assign = 1;
  }else{
    sendparticle_target_id.resize(0);
    sendparticle_setid.resize(0);
    recvparticle_target_id.resize(0);
    recvparticle_setid.resize(0);
    move_target_id.resize(0);
    move_setid.resize(0);
    targetset.resize(0);
    ghostpairs.resize(0);
    bondlistarray.resize(0);
    bondlistarray_idx.resize(0);
  }

  if(DebugLog::verbose>2){
    if(node_id==0){
     int nt=0;
     for(unsigned int i=0;i<recvparticle_setid.size();i++){
       nt += recvparticle_setid[i].size();
     }  
     printf("number of ghost for short0  %d\n",nt);
    }
    if(node_id==0){
      //    for(int t=0;t<sendparticle_target_id.size();t++)
      int t=0;
      {
        printf("send to %d :",sendparticle_target_id[t]);
        for(unsigned int c=0;c<sendparticle_setid[t].size();c++){
          printf(" %d",sendparticle_setid[t][c]);
        }
        printf("\n");
      }
    }
    if(node_id==1){
      //    for(int t=0;t<sendparticle_target_id.size();t++)
      int t=0;
      {
        printf("recv from %d :",recvparticle_target_id[t]);
        for(unsigned int c=0;c<recvparticle_setid[t].size();c++){
          printf(" %d",recvparticle_setid[t][c]);
        }
        printf("\n");
      }
    }
  }
}

  template<class PA>
void makeparticle(const int short_id,
                  const int total_num_particle,
                  const int total_num_freedom,
#if 1
                  const int num_copy,
                  const int num_copy_local,
#endif
                  const PA& pa, 
                  const WaterList& wl,
                  const ShakeList& sl,
                  const std::vector<PotentialModel>& pm,
                  const double ljcec)
{
  if(short_id==NO_SHORT_ID){
    particlearray.resize(0);
    return;
  }
  int num_lj;
  if((short_id==0)||(DebugLog::verbose>1)) {
    std::cout << "distribute_particle_cell" << std::endl;
  }
  //! number of cell in each axis
  int n;
  n = cellsubset.distribute_particle_cell(particlearray,
                                          particle_setid,
                                          typerangearray,
                                          waterlist,
                                          shakelist,
                                          pa, pm, wl, sl);
  if((short_id==0)||(DebugLog::verbose>1)) {
    std::cout << "make_atomid_map" << std::endl;
  }
  make_atomid_map(particlearray,typerangearray,atomid_to_index);

  if((short_id==0)||(DebugLog::verbose>1)) {
    std::cout << " short node " << short_id << " found " << n << " particle" ;
    if(DebugLog::verbose>1){
      int max_ppc=0;
      for(int s=0;s<int(typerangearray.size());s++){
        int ns = typerangearray[s].end-typerangearray[s].begin;
        if(ns>0){
          if(DebugLog::verbose>2){
            std::cout << " " << s << ":" << ns;
          }
          if(ns>max_ppc)max_ppc=ns;
        }
      }
      std::cout << " : max particle in cell " << max_ppc;
    }
    std::cout << std::endl;
  }

  num_all_particle = total_num_particle;
//  num_all_particle = 0;
  num_lj = 0;
  for(int i=0;i<int(pm.size());i++){
//    num_all_particle++;
    if((pm[i]==OnlyLJPotential)||(pm[i]==LJCoulombPotential)){
      num_lj++;
    }
  }
#if 1
  if((short_id==0)){
    printf("num_lj %d * %d / %d\n",num_lj,num_copy,num_copy_local);
  }
  num_lj = num_lj * num_copy / num_copy_local;
#endif
  setLJCutoffEnergyCorrectionParmeter(ljcec,num_lj);
  if((short_id==0)){
    printf("LJCutoffEnergyCorrectionV = %e\n",ShortRange::LJCutoffEnergyCorrectionV);
  }

  {
     double ljcecf = 0.0;
     {
       int nlj=0;
       std::vector<int> findatype;
       std::vector<int> numatype;
       for(int i=0;i<int(pm.size()-1);i++){
	 if((pm[i]==OnlyLJPotential)||(pm[i]==LJCoulombPotential)){
	   int ati = getatomtype(pa,i);
	   if(ati+1>findatype.size()){
	     findatype.resize(ati+1);
	     numatype.resize(ati+1);
	   }
	   numatype[ati]++;
	 }
       }
       for(int i=0;i<numatype.size();i++){
	 if(numatype[i]>0){
	   nlj += numatype[i];
	   for(int j=0;j<numatype.size();j++){
	     if(numatype[j]>0){
	       ljcecf += calcLJCutoffEnergyCorrection(cellindex.cutoff,i,j)*numatype[i]*numatype[j];
	     }
	   }
	 }
       }
       if((short_id==0)){
	 printf("counted ljatom  = %d\n",nlj);
       }
     }
     ljcecf = ljcecf*(num_copy/num_copy_local)*(num_copy/num_copy_local)*0.25;
     if((short_id==0)){
       printf("LJCutoffEnergyCorrectionV fine  = %e\n",ljcecf);
       printf("override LJCutoffEnergyCorrectionV\n");
     }
     ShortRange::LJCutoffEnergyCorrectionV = ljcecf;
  }
//  std::cout<< "num_all_particle (in makeparticle)" << num_all_particle << std::endl;
//  std::cout<< "num_water (in makeparticle)" << wl.size() << std::endl;
  num_freedom = total_num_freedom;
//  std::cout << "num_freedom (in makeparticle)" << num_freedom << std::endl;

}

#ifdef OLDPARTICLE
template
void makeparticle(const int short_id,
                  const int total_num_particle,
                  const int total_num_freedom,
#if 1
                  const int num_copy,
                  const int num_copy_local,
#endif
                  const ParticleArray& pa, 
                  const WaterList& wl,
                  const ShakeList& sl,
                  const std::vector<PotentialModel>& pm,
                  const double ljcec);
#else
template
void makeparticle(const int short_id,
                  const int total_num_particle,
                  const int total_num_freedom,
#if 1
                  const int num_copy,
                  const int num_copy_local,
#endif
                  const CombinedParticleArray& pa, 
                  const WaterList& wl,
                  const ShakeList& sl,
                  const std::vector<PotentialModel>& pm,
                  const double ljcec);
#endif

inline int atom_id_to_atom_index(const AtomID aid)
{
  std::map<AtomID,int>::iterator ait;
  ait = atomid_to_index.find(aid);
  if(ait!=atomid_to_index.end()){
    return ait->second;
  }
  return -1;
}

inline int atom_id_to_setid(const AtomID aid, int& atom_index)
{
  atom_index = -1;
#if 1
  std::map<AtomID,int>::iterator ait;
  ait = atomid_to_index.find(aid);
  if(ait!=atomid_to_index.end()){
    atom_index = ait->second;
    std::vector<TypeRange>::size_type s=0;
    for(s=0;s<typerangearray.size();s++){
      if((typerangearray[s].begin<=atom_index)&&(atom_index<typerangearray[s].end))break;
    }
    if(s<typerangearray.size()){
      return static_cast<int>(s);
    }
  }
#else
  for(std::vector<TypeRange>::size_type s=0;s<typerangearray.size();s++){
    for(int i=typerangearray[s].begin;i<typerangearray[s].end;i++){
      if(getatomid(particlearray,i) == aid){
        atom_index = i;
        return static_cast<int>(s);
      }
    }
  }
#endif
  return -1;
}

inline 
bool isWater(const AtomID atom)
{
  return ((atom==ATOMTYPE_WO)||(atom==ATOMTYPE_WH));
}

void makebond(const int short_id,
              CovalentBondList bl,
              bool excludewaterbond,
              const ShakeList& sl,
              int shake_type)
{
  if(short_id==NO_SHORT_ID){
    if(DebugLog::verbose>1){
      std::cout << " This node not short, skip makebond" << std::endl;
    }
    return;
  }
  int assigned_atom;
#ifndef FIX_CB_OFFSET
  assigned_atom = CovalentBondInfo::bond_assign_offset;
#else
  assigned_atom = bond_assign;
  if(bond_assign>1){
    assigned_atom = 1;
  }
#endif
  
  //  int found=0;
  CovalentBondInfo::Bond bond_idx;
  for(std::vector<CovalentBondInfo::Bond>::size_type b=0;b<bl.bond.size();b++){
    AtomID assign = bl.bond[b].id_of_atom[assigned_atom];

    int notassigned_atom;
    if(assigned_atom==0){
      notassigned_atom=1;
    }else if(assigned_atom==1){
      notassigned_atom=0;
    }else{
      std::cout << " Error: assigned_atom is invalid " << std::endl;
      return;
    }

    AtomID notassign = bl.bond[b].id_of_atom[notassigned_atom];

    //check shake bond or not
    int shake=0;
#ifdef USE_SHAKE
    if(sl.find(assign) != sl.end()){
      int nh1 = sl.find(assign)->second.nh1;

      for(int i=0; i<nh1; i++){
        int h1 = sl.find(assign)->second.h1[i];
        if(h1==notassign){
          if(shake_type > 0){
            shake = 1;    // this is shake bond
          }
          break;
        }
      }
    }
    if(shake==0 && sl.find(notassign) != sl.end()){
      int nh1 = sl.find(notassign)->second.nh1;

      for(int i=0; i<nh1; i++){
        int h1 = sl.find(notassign)->second.h1[i];
        if(h1==assign){
          if(shake_type > 0){
            shake = 1;    // this is shake bond
          }
          break;
        }
      }
    }
#endif

    int atom_index;
    int s = atom_id_to_setid(assign, atom_index);
    if(s>-1){
      if((excludewaterbond)&&(isWater(getatomtype(particlearray,atom_index) ))) continue;

      bl.bond[b].shake = shake;
      bond_idx.shake = shake;

      bondlistarray[s].BondArray.push_back(bl.bond[b]);
      for(int bi=0;bi<2;bi++){
        if(bi==assigned_atom){
          bond_idx.index_of_atom[bi] = atom_index;
        }else{
          bond_idx.index_of_atom[bi] = atom_id_to_atom_index(bl.bond[b].id_of_atom[bi]);
        }
      }
      bond_idx.typeofbond = bl.bond[b].typeofbond;
      bondlistarray_idx[s].BondArray.push_back(bond_idx);
      //      found = 1;
    }
  }
  //  found=0;
#ifndef FIX_CB_OFFSET
  assigned_atom = CovalentBondInfo::angle_assign_offset;
#endif
  CovalentBondInfo::Angle angle_idx;
  for(std::vector<CovalentBondInfo::Angle>::size_type b=0;
      b<bl.angle.size();b++){
    AtomID assign = bl.angle[b].id_of_atom[assigned_atom];
    int atom_index;
    int s = atom_id_to_setid(assign, atom_index);
    if(s>-1){
      if((excludewaterbond)&&(isWater(getatomtype(particlearray,atom_index) ))) continue;
      bondlistarray[s].AngleArray.push_back(bl.angle[b]);
      for(int bi=0;bi<3;bi++){
        if(bi==assigned_atom){
          angle_idx.index_of_atom[bi] = atom_index;
        }else{
          angle_idx.index_of_atom[bi] = atom_id_to_atom_index(bl.angle[b].id_of_atom[bi]);
        }
      }
      angle_idx.typeofangle = bl.angle[b].typeofangle;
      bondlistarray_idx[s].AngleArray.push_back(angle_idx);
      //      found = 1;
    }
  }
  //  found=0;
#ifndef FIX_CB_OFFSET
  assigned_atom = CovalentBondInfo::torsion_assign_offset;
#endif
  CovalentBondInfo::Torsion torsion_idx;
  for(std::vector<CovalentBondInfo::Torsion>::size_type b=0;
      b<bl.torsion.size();b++){
    AtomID assign = bl.torsion[b].id_of_atom[assigned_atom];
    int atom_index;
    int s = atom_id_to_setid(assign, atom_index);
    if(s>-1){
      if((excludewaterbond)&&(isWater(getatomtype(particlearray,atom_index) ))) continue;
      bondlistarray[s].TorsionArray.push_back(bl.torsion[b]);
      for(int bi=0;bi<4;bi++){
        if(bi==assigned_atom){
          torsion_idx.index_of_atom[bi] = atom_index;
        }else{
          torsion_idx.index_of_atom[bi] = atom_id_to_atom_index(bl.torsion[b].id_of_atom[bi]);
        }
      }
      torsion_idx.typeoftorsion = bl.torsion[b].typeoftorsion;
      bondlistarray_idx[s].TorsionArray.push_back(torsion_idx);
      //      found = 1;
    }
  }
#ifndef FIX_CB_OFFSET
  assigned_atom = CovalentBondInfo::improper_assign_offset;
#else
  if(bond_assign>1){
    assigned_atom = 2;
  }
#endif
  //  found=0;
  CovalentBondInfo::Improper improper_idx;
  for(std::vector<CovalentBondInfo::Improper>::size_type b=0;
      b<bl.improper.size();b++){
    AtomID assign = bl.improper[b].id_of_atom[assigned_atom];
    int atom_index;
    int s = atom_id_to_setid(assign, atom_index);
    if(s>-1){
      if((excludewaterbond)&&(isWater(getatomtype(particlearray,atom_index) ))) continue;
      bondlistarray[s].ImproperArray.push_back(bl.improper[b]);
      for(int bi=0;bi<4;bi++){
        if(bi==assigned_atom){
          improper_idx.index_of_atom[bi] = atom_index;
        }else{
          improper_idx.index_of_atom[bi] = atom_id_to_atom_index(bl.improper[b].id_of_atom[bi]);
        }
      }
      improper_idx.typeofimproper = bl.improper[b].typeofimproper;
      bondlistarray_idx[s].ImproperArray.push_back(improper_idx);
      //      found = 1;
    }
  }

  if((short_id==0)||(DebugLog::verbose>1)){
    int num_b=0, num_a=0, num_t=0, num_i=0;
    for(int s=0;s<int(typerangearray.size());s++){
      num_b += bondlistarray[s].BondArray.size();
      num_a += bondlistarray[s].AngleArray.size();
      num_t += bondlistarray[s].TorsionArray.size();
      num_i += bondlistarray[s].ImproperArray.size();
    }
    std::cout << " Number of CB on short node " << short_id;
    if((excludewaterbond)){
      std::cout << "(exclude water)";
    }else{
      std::cout << "(include water)";
    }
    std::cout << " Bond " << num_b;
    std::cout << " Angle " << num_a;
    std::cout << " Torsion " << num_t;
    std::cout << " Improper " << num_i;
    std::cout << std::endl;
    if(DebugLog::verbose>1){
      std::cout << "set bond angle torsion improper" << std::endl;
      for(int s=0;s<int(typerangearray.size());s++){
        std::cout << s;
        std::cout << " " << bondlistarray[s].BondArray.size();
        std::cout << " " << bondlistarray[s].AngleArray.size();
        std::cout << " " << bondlistarray[s].TorsionArray.size();
        std::cout << " " << bondlistarray[s].ImproperArray.size();
        std::cout << std::endl;
      }
      if(DebugLog::verbose>2){
        for(std::vector<CovalentBondInfo::BondList>::size_type s=0;
            s<bondlistarray.size();s++){
          std::cout << "bond in set " << s << " ";
          for(std::vector<CovalentBondInfo::Bond>::size_type b=0;
              b<bondlistarray[s].BondArray.size();b++){
            std::cout << "(" << bondlistarray[s].BondArray[b].id_of_atom[0] <<
                "," << bondlistarray[s].BondArray[b].id_of_atom[1] << ")";
          }
          std::cout << std::endl;
        }
        dump_bondlistarray(bondlistarray);
      }
    }
  }
}

void makelongrangegeometry(const LongRangeMPIPlan longplan,
                           const int num_long_only,
                           int& num_long,
                           GeometryXYZ& longgeometry)
{
  SpaceVector<int> longdomain;
  switch (longplan) {
    case Combine1D:
      longdomain = SpaceVector<int>(node_div.x,1,1);
      break;
    case Separate1D:
      longdomain = SpaceVector<int>(num_long_only,1,1);
      break;
    case Combine2D:
      longdomain = SpaceVector<int>(node_div.x,node_div.y,1);
      break;
    case Separate2D:
      int p;
      if(!power_of_two(num_long_only, p)){
        int m[2];
        if(power_of_two_of_square(num_long_only,m)){
          longdomain[0] = m[0];
          longdomain[1] = m[1];
          longdomain[2] = 1;
        }
      }else{
        longdomain[0] = 1<<(p/2);
        longdomain[1] = num_long_only/longdomain[0];
        longdomain[2] = 1;
      }
      break;
    case Combine3D:
      longdomain = node_div;
      break;
    case Separate3D:
      if(!power_of_two(num_long_only, p)){
        int m[3];
        if(power_of_two_of_cube(num_long_only,m)){
          longdomain[0] = m[0];
          longdomain[1] = m[1];
          longdomain[2] = m[2];
        }
      }else{
        int sp[3];
        split_power_of_two(p,sp);
        longdomain[0] = 1<<sp[0];
        longdomain[1] = 1<<sp[1];
        longdomain[2] = 1<<sp[2];
      }
      break;
    default:
        std::cout << "unknown long rage plan. force 0D" << std::endl;
    case Combine0D:
    case Separate0D:
      longdomain = SpaceVector<int>(1,1,1);
      break;
    case NoneLong:
      longdomain = SpaceVector<int>(0,0,0);
  }
  num_long = longdomain.x*longdomain.y*longdomain.z;
  //  longgeometry = GeometryXYZ(longdomain,longdomain);   // long decomposition just longnode-geometry
  longgeometry = GeometryXYZ(longdomain,cell_div3d); // short cells  decomposed to long node
}

//! TODO Combine and num_node > num_short
void makelongcomm(const LongRangeMPIPlan longplan,
                  const GeometryXYZ longgeometry,
                  const int num_long,
                  const int num_short,
                  const MPI_Comm mpi_comm_all,
                  const int node_id,
                  const int short_id,
                  int& long_id)
{
  int key = 0;
  int color = 0;
  int long0 = 0;
  int size_long;
  if(short_id==0){
    std::cout << "longgeometry.size " << longgeometry.size << std::endl;
  }
  int num_world;
  MPI_Comm_size(MPI_COMM_WORLD,&num_world);
  switch (longplan) {
    case Combine0D:
      if(short_id==0){
        key = 0;
        color = 1;
      }
      break;
    case Combine1D:
      //      if((short_id!=NO_SHORT_ID)&&((short_id/(node_div.y*node_div.z))==0)){ // 
      if((short_id!=NO_SHORT_ID)&&(short_id<node_div.x)) {
        key = short_id;
        color = 1;
      }
      break;
    case Combine2D:
      //      if((short_id!=NO_SHORT_ID)&&((short_id/node_div.z)==0)){
      if((short_id!=NO_SHORT_ID)&&(short_id<node_div.x*node_div.y)){
        key = short_id;
        color = 1;
      }
      break;
    case Combine3D:
      if((short_id!=NO_SHORT_ID)){
        key = short_id;
        color = 1;
      }
      break;
    case Separate0D:
    case Separate1D:
    case Separate2D:
    case Separate3D:
      if(short_id==NO_SHORT_ID){
        key = node_id;
        color = 1;
      }
    default:
      if(num_short<num_world){
        long0 = num_short;
      }else{
        long0 = 0;
        size_long = 0;
      }
      break;
  }
  if((node_id==0)){
    std::cout << " num_world " << num_world << " num_short "  << num_short  << std::endl;
  }
  MPI_Comm_split(mpi_comm_all,color,key,&mpi_comm_long);
  if((mpi_comm_long!=MPI_COMM_NULL)&&(color==1)){
    MPI_Comm_rank(mpi_comm_long,&long_id);
    MPI_Comm_size(mpi_comm_long, &size_long);
    longid_to_worldid.resize(size_long);
    MPI_Gather((void *)&node_id,1,MPI_INT,&(longid_to_worldid[0]),1,MPI_INT,0,mpi_comm_long);
    if(long_id==0){
      if(long0!=node_id){
        printf("Error world rank %d of long rank 0 is not %d\n",node_id, long0);
      }
      int expected_size = longgeometry.size.x*longgeometry.size.y*longgeometry.size.z;
      if(size_long!=expected_size){
        printf("Error number of long node is %d, but  expeceted %d\n",size_long,expected_size);
      }
    }
  }else{
    if(DebugLog::verbose>1){
      std::cout << node_id << " :No Long Comm " << std::endl;
    }
    mpi_comm_long=MPI_COMM_NULL;
    long_id = NO_LONG_ID;
  }
  MPI_Bcast(&size_long,1,MPI_INT,long0,MPI_COMM_WORLD);
  if(long_id!=0){
    longid_to_worldid.resize(size_long);
  }
  MPI_Bcast(&(longid_to_worldid[0]),size_long,MPI_INT,long0,MPI_COMM_WORLD);
  if((node_id==0)){
    std::cout << "size of Long " << size_long << std::endl;
  }
  if((node_id==0)&&(DebugLog::verbose>1)){
    printf("(long_id, world_id) ");
    for(int i=0;i<size_long;i++){
      printf("(%d %d)",i,longid_to_worldid[i]);
    }
    printf("\n");
  }
}

void long_corner_cell(const SpaceVector<int>ldiv,
                      const double fringe,
                      const int long_id,
                      SpaceVector<int>& cell_min,
                      SpaceVector<int>& cell_max)
{
  if(long_id!=NO_LONG_ID){
    int z = long_id/(ldiv.x*ldiv.y);
    int y = (long_id-z*(ldiv.x*ldiv.y))/ldiv.x;
    int x = long_id-(y+z*ldiv.y)*ldiv.x;
#if defined(CPPMD_ENABLE_MR3EXAFMM)
    double min_x = cellindex.boxsize.x/ldiv.x*x-1;
    double max_x = cellindex.boxsize.x/ldiv.x*(x+1)+1;
    double min_y = cellindex.boxsize.y/ldiv.y*y-1;
    double max_y = cellindex.boxsize.y/ldiv.y*(y+1)+1;
    double min_z = cellindex.boxsize.z/ldiv.z*z-1;
    double max_z = cellindex.boxsize.z/ldiv.z*(z+1)+1;

    cell_min.x = int(floor(min_x*cellindex.invcellsize.x));
    cell_max.x = int(ceil(max_x*cellindex.invcellsize.x));
    cell_min.y = int(floor(min_y*cellindex.invcellsize.y));
    cell_max.y = int(ceil(max_y*cellindex.invcellsize.y));
    cell_min.z = int(floor(min_z*cellindex.invcellsize.z));
    cell_max.z = int(ceil(max_z*cellindex.invcellsize.z));
    /*
    SpaceVector<double> min(min_x,min_y,min_z);
    SpaceVector<double> max(max_x,max_y,max_z);
    if((long_id==0)&&(DebugLog::verbose>0)){
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      std::cout << "long corner cell for " << long_id << "(" << rank << ")"<< " : " << cell_min << cell_max << " : " << min << max << " : at ldiv " << ldiv << std::endl;
    }
    */
#elif defined(CPPMD_ENABLE_GFMM)
    // rank0 of long node has all cell as ghost and bcast to others
    // other rank has no ghost and receive from rank0
    cell_min = SpaceVector<int>(0,0,0);
    if(long_id==0){
      cell_max.x = cell_div3d.x+1;
      cell_max.y = cell_div3d.y+1;
      cell_max.z = cell_div3d.z+1;
    }else{
      cell_max = SpaceVector<int>(0,0,0);
    }
#elif defined(CPPMD_ENABLE_EWALD)
    // This is correct only long node border is just on cell border
    // Othercase, ewald class must know border and select atom
    double min_x, min_y, min_z, max_x, max_y, max_z;
    if(cell_div3d.x%ldiv.x==0){
      min_x = cell_div3d.x/ldiv.x*x;
      max_x = cell_div3d.x/ldiv.x*(x+1);
    }else{
      min_x = double(cell_div3d.x)/ldiv.x*x;
      max_x = double(cell_div3d.x)/ldiv.x*(x+1);
    }
    if(cell_div3d.y%ldiv.y==0){
      min_y = cell_div3d.y/ldiv.y*y;
      max_y = cell_div3d.y/ldiv.y*(y+1);
    }else{
      min_y = double(cell_div3d.y)/ldiv.y*y;
      max_y = double(cell_div3d.y)/ldiv.y*(y+1);
    }
    if(cell_div3d.z%ldiv.z==0){
      min_z = cell_div3d.z/ldiv.z*z;
      max_z = cell_div3d.z/ldiv.z*(z+1);
    }else{
      min_z = double(cell_div3d.z)/ldiv.z*z;
      max_z = double(cell_div3d.z)/ldiv.z*(z+1);
    }
    cell_min.x = int(floor(min_x));
    cell_max.x = int(ceil(max_x));
    cell_min.y = int(floor(min_y));
    cell_max.y = int(ceil(max_y));
    cell_min.z = int(floor(min_z));
    cell_max.z = int(ceil(max_z));
#else
    SpaceVector<double>& margin = cellindex.margin;
    double min_x = cellindex.boxsize.x/ldiv.x*x-margin.x-fringe;
    double max_x = cellindex.boxsize.x/ldiv.x*(x+1)+margin.x+fringe;
    double min_y = cellindex.boxsize.y/ldiv.y*y-margin.y-fringe;
    double max_y = cellindex.boxsize.y/ldiv.y*(y+1)+margin.y+fringe;
    double min_z = cellindex.boxsize.z/ldiv.z*z-margin.z-fringe;
    double max_z = cellindex.boxsize.z/ldiv.z*(z+1)+margin.z+fringe;

    cell_min.x = int(floor(min_x*cellindex.invcellsize.x));
    cell_max.x = int(ceil(max_x*cellindex.invcellsize.x));
    cell_min.y = int(floor(min_y*cellindex.invcellsize.y));
    cell_max.y = int(ceil(max_y*cellindex.invcellsize.y));
    cell_min.z = int(floor(min_z*cellindex.invcellsize.z));
    cell_max.z = int(ceil(max_z*cellindex.invcellsize.z));
#endif
  }else{
    cell_min = SpaceVector<int>(0,0,0);
    cell_max = SpaceVector<int>(0,0,0);
  }
}

void makelongrangerequiredcell(const GeometryXYZ longgeometry,
                               const double fringe,
                               const int long_id,
                               std::vector<int>& cellid_list)
{
  if(long_id==NO_LONG_ID){
    cellid_list.resize(0);
    return;
  }
  SpaceVector<int> ldiv(longgeometry.size);
  SpaceVector<int> cell_min;
  SpaceVector<int> cell_max;

  /*
  if((long_id==0)||(DebugLog::verbose>1)){
    std::cout << "long div " << ldiv << " fringe " << fringe << std::endl;
  }
  */
  long_corner_cell(ldiv,fringe,long_id,cell_min,cell_max);
  for(int cell_z=cell_min.z;cell_z<cell_max.z;cell_z++){
    for(int cell_y=cell_min.y;cell_y<cell_max.y;cell_y++){
      for(int cell_x=cell_min.x;cell_x<cell_max.x;cell_x++){
        SpaceVector<int> pos(cell_x,cell_y,cell_z);
        periodic_shift(pos.x,cellindex.celldiv.x);
        periodic_shift(pos.y,cellindex.celldiv.y);
        periodic_shift(pos.z,cellindex.celldiv.z);
        int cid = cell_position_to_cellid(pos,cellindex.celldiv);
        cellid_list.push_back(cid);
      }
    }
  }
  std::sort(cellid_list.begin(),cellid_list.end());
  std::vector<int>::iterator end_it = std::unique(cellid_list.begin(),cellid_list.end());
  cellid_list.erase(end_it,cellid_list.end());
  //  if((long_id==0)||(DebugLog::verbose>1)){
  if((DebugLog::verbose>1)){
    printf("required cell for long %ld\n",cellid_list.size());
    fflush(stdout);
  }
}

// send_to_long_list[long_id][] : set_ids this node must send to long node with long_id,
//   not mention long node is self when Combine*D
void makecellidlistsendtolong(const GeometryXYZ longgeometry,
                              const double fringe,
                              const int num_long,
                              const int short_id,
                              std::vector< std::vector<int> >& send_to_long_list)
{
  if(short_id==NO_SHORT_ID){
    return;
  }
  send_to_long_list.resize(num_long);
  for(int long_id=0;long_id<num_long;long_id++){
    send_to_long_list[long_id].clear();
    send_to_long_list[long_id].resize(0);
  }
  std::vector<int> reqlist;
  for(int long_id=0;long_id<num_long;long_id++){
    reqlist.clear();
    makelongrangerequiredcell(longgeometry, fringe,long_id,reqlist);
    for(std::vector<int>::size_type i=0;i<reqlist.size();i++){
      std::vector<int>::iterator it = std::find(particle_setid.begin(),particle_setid.end(),reqlist[i]);
      if(it!=particle_setid.end()){
        int ri;
        int worldid=longid_to_worldid[long_id];
        int max_ri = sendparticle_target_id.size();
        for(ri=0;ri<max_ri;ri++){
          if(sendparticle_target_id[ri]==worldid)break;
        }
        if(ri<max_ri){
          std::vector<int>::iterator cell_it = std::find(sendparticle_setid[ri].begin(),sendparticle_setid[ri].end(),*it);
          if(cell_it!=sendparticle_setid[ri].end()){
            continue;
          }else{
            if(short_id==0){
              printf("cell %d is not in ghost to %d(%d)\n",*it,worldid,ri);
            }
          }
        }
        send_to_long_list[long_id].push_back(*it);
      }
    }
  }
  if(short_id==0){
    if(num_long>1){
      printf("send_to_long_list[1] :");
      for(unsigned int i=0;i<send_to_long_list[1].size();i++){
        printf(" %d",send_to_long_list[1][i]);
      }
      printf("\n");
    }
  }
}

void calc_cell_origin(const SpaceVector<int> cell_pos, SpaceVector<double>& origin)
{
  origin.x = cellindex.cellsize.x*cell_pos.x;
  origin.y = cellindex.cellsize.y*cell_pos.y;
  origin.z = cellindex.cellsize.z*cell_pos.z;
}

void calc_self_energy_box(const SpaceVector<int>ldiv, const int long_id,
                          SpaceVector<double>& min_corner,
                          SpaceVector<double>& max_corner)
{
  if(long_id!=NO_LONG_ID){
    int z = long_id/(ldiv.x*ldiv.y);
    int y = (long_id-z*(ldiv.x*ldiv.y))/ldiv.x;
    int x = long_id-(y+z*ldiv.y)*ldiv.x;
    min_corner.x = cellindex.boxsize.x/ldiv.x*x;
    max_corner.x = cellindex.boxsize.x/ldiv.x*(x+1);
    min_corner.y = cellindex.boxsize.x/ldiv.y*y;
    max_corner.y = cellindex.boxsize.x/ldiv.y*(y+1);
    min_corner.z = cellindex.boxsize.x/ldiv.z*z;
    max_corner.z = cellindex.boxsize.x/ldiv.z*(z+1);
  }else{
    min_corner = SpaceVector<double>(0.0,0.0,0.0);
    max_corner = SpaceVector<double>(0.0,0.0,0.0);
  }
}

void makeselfenergycell_list(const GeometryXYZ longgeometry,
                             const int long_id,
                             const std::vector<int>& cellid_list,
                             std::vector<int>& selfenergycell_list)
{
  selfenergycell_list.resize(0);
  if(long_id==NO_LONG_ID){
    return;
  }
  SpaceVector<int> ldiv(longgeometry.size);
  SpaceVector<double> cell_origin;
  SpaceVector<double> min_corner;
  SpaceVector<double> max_corner;
  calc_self_energy_box(ldiv,long_id,min_corner,max_corner);
  for(std::vector<int>::size_type i=0;i<cellid_list.size();i++){
    int cid = cellid_list[i];
    SpaceVector<int> pos = cellid_to_cell_position(cid, cellindex.celldiv);
    calc_cell_origin(pos,cell_origin);
    if((cell_origin.x>=min_corner.x)&&(cell_origin.x<max_corner.x)
       &&(cell_origin.y>=min_corner.y)&&(cell_origin.y<max_corner.y)
       &&(cell_origin.z>=min_corner.z)&&(cell_origin.z<max_corner.z)){
      selfenergycell_list.push_back(cid);
    }
  }
  if((long_id==0)||(DebugLog::verbose>1)) {
    std::cout << "number of cell for selfenergy " << selfenergycell_list.size() << std::endl;
  }
}

//! TODO remove sets that transferd for short when longplan is Combine?D
void make_long_short_comm(const int num_node,
                          const int num_long,
                          const int node_id,
                          const int short_id,
                          const int long_id,
                          const MPI_Comm mpi_comm_all,
                          const std::vector< std::vector<int> >& send_to_long_list,
                          const std::vector<int>& long_reqcell_list,
                          std::vector<MPI_Comm>& mpi_comm_ls_list,
                          std::vector<int>& long_short_id,
                          std::vector<int>& sender_local_id,
                          std::map<int,int>& idinlong_to_longrank,
                          std::vector<int>& receiver_local_rank,
                          std::vector< std::vector<int> >& long_recv_set_id_lists,
                          std::vector<int>& sender_global_id_list,
                          std::vector<int>& reciever_global_id_list)
{
  mpi_comm_ls_list.clear();

  int color;
  MPI_Comm comm_ls;
  int key;

  long_short_id.clear();
  sender_global_id_list.clear();
  reciever_global_id_list.clear();
#ifdef CPPMD_ENABLE_MR3EXAFMM
///// it correct when --calc-space=30 : longplan = Combine3D
///// it is not required, disable long-short communication
  color = long_id;
  key = 0;
  MPI_Comm_split(mpi_comm_all,color,key,&comm_ls);
  for(int lid=0;lid<num_long;lid++){
    if(lid==long_id){
      mpi_comm_ls_list.push_back(comm_ls);
    }else{
      mpi_comm_ls_list.push_back(MPI_COMM_NULL);
    }
    if(lid==long_id){
      int ls_id;
      MPI_Comm_rank(comm_ls,&ls_id);
      long_short_id.push_back(ls_id);
      int num_ls;
      MPI_Comm_size(comm_ls,&num_ls);
      int *ids = new int[num_ls];
      if(num_ls>1){
	printf("num_ls %d != 1\n",num_ls);
	num_ls = 1;
      }
      ids[0] = node_id;
      {  // Long node of this sub comm
	if((short_id==0)||(DebugLog::verbose>1)) {
	  std::cout << " Rank in long MPI_Comm of Long node is " << ls_id;
	  if(ls_id!=0){
	    std::cout << " , not 0";
	  }
	  std::cout << std::endl;
	  std::cout << "size long MPI_Comm is " << num_ls << std::endl;
	  std::cout << "max set size  " << long_reqcell_list.size() << std::endl;
	}
	sender_global_id_list.resize(num_ls-1);
	long_recv_set_id_lists.resize(num_ls-1);
	int *setid = new int[num_set_per_node+1];
	int num_remote_set=0;
	if((short_id==0)||(DebugLog::verbose>1)) {
	  std::cout << "long " << lid << " remote set " << num_remote_set << std::endl;
	}
      }
      receiver_local_rank.push_back(int(0));
      reciever_global_id_list.push_back(ids[0]);
    }else{
      long_short_id.push_back(NO_LONG_ID);
      receiver_local_rank.push_back(NO_LONG_ID);
      reciever_global_id_list.push_back(NO_LONG_ID);
    }
  }
#else  // CPPMD_ENABLE_MR3EXAFMM
  for(int lid=0;lid<num_long;lid++){
    color = MPI_UNDEFINED;
    key = 0;
    if(long_id==lid){  // this node is S-L group for lid, just lid
      color = lid;
    }else{
      if(short_id!=NO_SHORT_ID){
        if(send_to_long_list[lid].size()>0){
          color = lid;
          key = short_id + 1;  // make key >= 1 > 0(key of long node)
        }
      }
    }
    MPI_Comm_split(mpi_comm_all,color,key,&comm_ls);
    mpi_comm_ls_list.push_back(comm_ls);
    if(color!=MPI_UNDEFINED){
      int ls_id;
      MPI_Comm_rank(comm_ls,&ls_id);
      long_short_id.push_back(ls_id);
      int num_ls;
      MPI_Comm_size(comm_ls,&num_ls);
      int *ids = new int[num_ls];
      MPI_Allgather((void *)&node_id,1,MPI_INT,ids,1,MPI_INT,comm_ls);
      if(long_id==lid){  // Long node of this sub comm
	if((short_id==0)||(DebugLog::verbose>1)) {
	  std::cout << " Rank in long MPI_Comm of Long node is " << ls_id;
	  if(ls_id!=0){
	    std::cout << " , not 0";
	  }
	  std::cout << std::endl;
	  std::cout << "size long MPI_Comm is " << num_ls << std::endl;
	  std::cout << "max set size  " << long_reqcell_list.size() << std::endl;
	}
        sender_global_id_list.resize(num_ls-1);
        long_recv_set_id_lists.resize(num_ls-1);
        int *setid = new int[num_set_per_node+1];
        int num_remote_set=0;
        for(int i=1;i<num_ls;i++){
	  if((short_id==0)||(DebugLog::verbose>1)) {
	    std::cout << "set long " << lid << " rank " << i << " max " << num_set_per_node << std::endl;
	  }
          sender_local_id.push_back(i);
          idinlong_to_longrank.insert(std::pair<int,int>(i,i));
          sender_global_id_list[i-1] = ids[i];
          MPI_Status stat;
          if((short_id==0)||(DebugLog::verbose>1)) {
	    std::cout << "long " << lid << " rank " << int(0) << " MPI_Recv " << std::endl;
	  }
          MPI_Recv(setid,num_set_per_node+1,MPI_INT,i,i,comm_ls,&stat);
          long_recv_set_id_lists[i-1].clear();
          for(int si=0;si<setid[0];si++){
            long_recv_set_id_lists[i-1].push_back(setid[si+1]);
          }
          num_remote_set += long_recv_set_id_lists[i-1].size();
        }
        if((short_id==0)||(DebugLog::verbose>1)) {
	  std::cout << "long " << lid << " remote set " << num_remote_set << std::endl;
	}
      }else{  // Short node of this tsub comm
        int *setid = new int[send_to_long_list[lid].size()+1];
        setid[0] = send_to_long_list[lid].size();
        for(int si=0;si<setid[0];si++){
          setid[si+1] = send_to_long_list[lid][si];
        }
	if((short_id==0)||(DebugLog::verbose>1)) {
	  std::cout << "long " << lid << " rank " << ls_id << "(short " << short_id << " world " << node_id << ") MPI_Send number of set and " << setid[0] << " set"  << std::endl;
	}
        MPI_Send(setid,setid[0]+1,MPI_INT,0,ls_id,comm_ls);
      }
      receiver_local_rank.push_back(int(0));
      reciever_global_id_list.push_back(ids[0]);
    }else{
      long_short_id.push_back(NO_LONG_ID);
      receiver_local_rank.push_back(NO_LONG_ID);
      reciever_global_id_list.push_back(NO_LONG_ID);
    }
  }
#endif  // CPPMD_ENABLE_MR3EXAFMM
}

void makelongrangeparameter(double cutoff,
                            double alpha,
                            double kCutoff,
                            int surfaceDipole,
                            SpaceVector<double> boxSize,
                            SpaceVector<double> gridLengths,
                            int order,
                            PMEType pmeType,
                            SpaceVector<int> grid_num,
                            GeometryXYZ node_geometry,
                            MGType multigridType,
                            int multigridIteration,
                            int multigridFFTlevel)
{
#ifdef CPPMD_ENABLE_LONGRANGE
# ifdef CPPMD_ENABLE_PMMM
  if(alpha==0.0){
    //        printf("Estimate alpha for PMMM\n");
    PMMMPreparator::estimate_alpha(alpha, boxSize);
  }
# else
  EwaldModule::estimate_alpha_kCutoff(alpha,kCutoff,boxSize.getComponentMax(),cutoff);
# endif
#endif  // CPPMD_ENABLE_LONGRANGE
  ShortRange::alpha = alpha;
  longrangeparameter.cutoff = cutoff;
  longrangeparameter.alpha = alpha;
  longrangeparameter.kCutoff = kCutoff;
  longrangeparameter.surfaceDipole = surfaceDipole;
  longrangeparameter.boxSize = boxSize;
  longrangeparameter.gridLengths = gridLengths;
  longrangeparameter.order = order;
  longrangeparameter.pmeType = pmeType;
  longrangeparameter.grid_num = grid_num;
  longrangeparameter.node_geometry = node_geometry;
  longrangeparameter.multigridType = multigridType;
  longrangeparameter.multigridIteration = multigridIteration;
  longrangeparameter.multigridFFTlevel = multigridFFTlevel;
}

void constructcommunicator(const int short_id,
                           const int max_particle_in_cell,
                           const std::map<int,int>& shortidtorank,
                           const int num_long,
                           const int long_id,
                           const std::vector<MPI_Comm>& mpi_comm_ls_list,
                           const std::vector<int>& long_short_id,
                           const std::vector<int>& sender_local_id,
                           const std::vector< std::vector<int> >& long_recv_set_id_lists,
                           const std::vector<int>& receiver_local_rank,
                           const std::vector< std::vector<int> >& send_to_long_list,
                           const std::map<int,int>& idinlong_to_longrank,
                           const int  short_comm_pattern,
                           const int  move_comm_pattern)
{
  CommPattern short_c_p = ((short_comm_pattern==0) ? Direct : NearestXYZ);
  CommPattern move_c_p = ((move_comm_pattern==0) ? Direct : NearestXYZ);
  if(short_id==NO_SHORT_ID){
    short_c_p=Direct;
    move_c_p=Direct;
  }
  communicator = MPICommunicator(shortidtorank,
                                 mpi_comm_short,
                                 sendparticle_target_id,
                                 sendparticle_setid,
                                 particle_setid,
                                 recvparticle_target_id,
                                 recvparticle_setid,
                                 move_target_id,
                                 move_setid,
                                 max_particle_in_cell,
                                 max_id, short_id,
                                 node_div,
                                 cellindex.celldiv,
                                 cellindex.move_range,
                                 cellindex.node_move_surface_ratio,
                                 short_c_p,
                                 move_c_p
                                 );

  if(DebugLog::verbose>1){
    std::cout << "communicators unit_identifier";
    std::cout << " " << communicator.unit_id;
    std::cout << std::endl;
  }

  if(long_id!=NO_LONG_ID){
    longreceiver = MPIReceiverForLong(idinlong_to_longrank,
                                      mpi_comm_ls_list[long_id],
                                      sender_local_id,
                                      long_recv_set_id_lists,
                                      max_particle_in_cell,
                                      max_id,
                                      long_short_id[long_id]);
  }else{
    ope.doLongrangecalculation = false;
    //    ope.doShortLongCommunication = false;
  }
  if(short_id!=NO_SHORT_ID){
    std::vector<int> receiver_local_id = receiver_local_rank;
    longsender = MPISenderForLong(receiver_local_rank,
                                  mpi_comm_ls_list,
                                  receiver_local_id,
                                  send_to_long_list,
                                  particle_setid,
                                  max_particle_in_cell,
                                  max_id,
                                  long_short_id
                                  );
  }
}

template<class PA, class GPA>
CalculationUnit<PA,GPA>* constructcalculationunit(const int node_id,
                                          const int short_id,
                                          const SpaceVector<double> boxsize,
                                          const CellIndexType cellindextype,
                                          const double cutoff,
#ifdef USE_PAIRLIST
                                          const double pmargin,
#endif
                                          const Config::TempCtrl& temp_control,
                                          const ShortRange::CoulombType cltype,
                                          const CovalentBondParameterList& covalent_bond_parameter_list)
{
  //  std::cout << "new CalculationUnit " << std::endl;
#ifdef CPPMD_ENABLE_LONGRANGE
  if(node_id==0){
    printf("alpha : %e %e\n",longrangeparameter.alpha,ShortRange::alpha);
    printf("kCutoff : %e\n",longrangeparameter.kCutoff);
  }
#endif
  if(DebugLog::verbose>1){
    std::cout << "communicator unit_identifier";
    std::cout << " " << communicator.unit_id;
    std::cout << std::endl;
  }
  PA pass_particlearray(particlearray);
   // DEBUG code
//  MPI_Barrier(MPI_COMM_WORLD);
//  printf("call CalculationUnit %d\n",node_id);
//  fflush(stdout);
  
  return new CalculationUnit<PA,GPA>(pass_particlearray,
                             particle_setid,
                             typerangearray,
                             waterlist,
                             shakelist,
                             bondlistarray,
                             bondlistarray_idx,
                             covalent_bond_parameter_list,
                             targetset,
                             ghostpairs,
                             ope,communicator,
                             longreceiver,
                             longsender,
                             postprocess,
                             cltype,
                             longrangeparameter,
                             node_id,
                             short_id,
                             boxsize,
			     boxsize_minimum_scale,
                             cutoff,
#ifdef USE_PAIRLIST
                             pmargin,
#endif
                             num_all_particle,
                             num_freedom,
                             mpi_comm_long,
                             temp_control,
                             expire_reverse);
  /*
  particlearray.clear();
  ParticleArray(particlearray).swap(particlearray);
  */
}
#ifdef OLDPARTICLE
template
CalculationUnit<ParticleArray,ParticleArray>* constructcalculationunit(const int node_id,
                                          const int short_id,
                                          const SpaceVector<double> boxsize,
                                          const CellIndexType cellindextype,
                                          const double cutoff,
#ifdef USE_PAIRLIST
                                          const double pmargin,
#endif
                                          const Config::TempeCtrl& temp_control,
                                          const ShortRange::CoulombType cltype,
                                                                       const CovalentBondParameterList& covalent_bond_parameter_list);
#else
template
CalculationUnit<CombinedParticleArray,GhostParticleArray>* constructcalculationunit(const int node_id,
                                          const int short_id,
                                          const SpaceVector<double> boxsize,
                                          const CellIndexType cellindextype,
                                          const double cutoff,
#ifdef USE_PAIRLIST
                                          const double pmargin,
#endif
                                          const Config::TempCtrl& temp_control,
                                          const ShortRange::CoulombType cltype,
                                                                       const CovalentBondParameterList& covalent_bond_parameter_list);
#endif
void constructselflongset(const int long_id, const int short_id,
                          const std::vector<int>& short_set_list,
                          const std::vector< std::vector<int> >& send_to_long_list,
                          const std::vector<int>& selfenergycell_list,
                          std::vector<int>& longset_index,
                          std::vector<int>& selfenergycell_index)
{
  longset_index.clear();
  selfenergycell_index.clear();
  if((long_id!=NO_LONG_ID)&&(short_id!=NO_SHORT_ID)){
    for(std::vector< std::vector<int> >::size_type i=0;
        i<send_to_long_list[long_id].size();i++){
      std::vector<int>::size_type index;
      for(index=0;index<short_set_list.size();index++){
        if(short_set_list[index]==send_to_long_list[long_id][i])break;
      }
      if(index<short_set_list.size()){
        longset_index.push_back(static_cast<int>(index));
        std::vector<int>::const_iterator it 
            = std::find(selfenergycell_list.begin(),
                        selfenergycell_list.end(),
                        short_set_list[index]);
        if(it!=selfenergycell_list.end()){
          selfenergycell_index.push_back(static_cast<int>(index));
        }
      }
    }
    if(DebugLog::verbose>1){
      std::cout << " found " << longset_index.size() << " sets for self long " << std::endl;
    }
  }
}

void constructghostlongset(const int long_id, const int short_id,
                           const std::vector<int>& long_reqcell_list,
                           const std::vector< std::vector<int> >& long_recv_set_id_lists,
                           const std::vector<int>& short_set_list,
                           const std::vector<int>& selfenergycell_list,
                           std::vector<int>& ghostlongset,
                           std::vector<int>& ghost_selfenergycell_list)
{
  // TEST, must implement correct sets
  int ns=0, nl=0;
  std::vector<int> inghost(long_reqcell_list.size(),0);
  ghostlongset.clear();
  ghost_selfenergycell_list.clear();
  if(long_id!=NO_LONG_ID){
    if(short_id!=NO_SHORT_ID){
      for(unsigned int st=0;st<recvparticle_setid.size();st++){
        for(unsigned int ci=0;ci<recvparticle_setid[st].size();ci++){
          int cell = recvparticle_setid[st][ci];
          unsigned int lp;
          for(lp=0;lp<long_reqcell_list.size();lp++){
            if(long_reqcell_list[lp]==cell)break;
          }
          if(lp<long_reqcell_list.size()){
            ghostlongset.push_back(cell);
            ns++;
            if(inghost[lp]>0){
              printf("cell %d is duplicate in short ghost\n",cell);
            }
            inghost[lp]++;
            std::vector<int>::const_iterator it 
                = std::find(selfenergycell_list.begin(),
                            selfenergycell_list.end(),
                            cell);
            if(it!=selfenergycell_list.end()){
              ghost_selfenergycell_list.push_back(cell);
            }
          }
        }
      }
    }
    for(std::vector< std::vector<int> >::size_type t=0;
        t<long_recv_set_id_lists.size();t++){
      for(std::vector<int>::size_type i=0;
          i<long_recv_set_id_lists[t].size();i++){
        ghostlongset.push_back(long_recv_set_id_lists[t][i]);
        nl++;
        {
          int cell = long_recv_set_id_lists[t][i];
          unsigned int lp;
          for(lp=0;lp<long_reqcell_list.size();lp++){
            if(long_reqcell_list[lp]==cell)break;
          }
          if(lp<long_reqcell_list.size()){
            if(inghost[lp]>0){
              if(inghost[lp]%10>0){
                printf("cell %d is duplicate in short and long ghost\n",cell);
              }else{
                printf("cell %d is duplicate in long ghost\n",cell);
              }
            }
            inghost[lp]+=10;
          }else{
            printf("cell %d is not required for long\n",cell);
          }
        }
        std::vector<int>::const_iterator it 
            = std::find(selfenergycell_list.begin(),
                        selfenergycell_list.end(),
                        long_recv_set_id_lists[t][i]);
        if(it!=selfenergycell_list.end()){
          ghost_selfenergycell_list.push_back(long_recv_set_id_lists[t][i]);
        }
      }
    }
    if(DebugLog::verbose>1){
      std::cout << " set " << ghostlongset.size() << "(" << ns << " " <<  nl << ") ghost set for only long" << std::endl;
    }
  }
}
#ifdef CPPMD_ENABLE_FMM
  void construct_fmm_target_cell(CalcForce &cf)
  {
    std::vector< std::vector<int> > fmm_target_cell;
    std::vector<int> scell;
    std::vector<int> rcell;
    ShiftCellArray shift_cell;
    fmm_target_cell.resize(particle_setid.size());
    for(int i=0;i<particle_setid.size();i++){
      scell.clear();
      rcell.clear();
      shift_cell.clear();
      fmm_target_cell[i].clear();
      cubic_target(cellindex.celldiv,cellindex.cell_position[particle_setid[i]],1,scell,rcell,shift_cell);
      fmm_target_cell[i] = rcell;
    }
    cf.clear_fmmtarget_set();
    cf.generate_fmmtarget_set(fmm_target_cell,expire_reverse);
  }
#endif  // CPPMD_ENABLE_FMM
#ifdef CPPMD_ENABLE_PMMM
  void construct_pmmm_target_cell(CalcForce &cf)
  {
    std::vector< std::vector<int> > pmmm_target_cell;
    std::vector<int> scell;
    std::vector<int> rcell;
    ShiftCellArray shift_cell;
    pmmm_target_cell.resize(particle_setid.size());
    for(int i=0;i<particle_setid.size();i++){
      scell.clear();
      rcell.clear();
      shift_cell.clear();
      pmmm_target_cell[i].clear();
      cubic_target(cellindex.celldiv,cellindex.cell_position[particle_setid[i]],cellindex.target_range.x,scell,rcell,shift_cell);
      pmmm_target_cell[i] = rcell;
    }
    cf.clear_pmmmtarget_set();
    cf.generate_pmmmtarget_set(pmmm_target_cell,expire_reverse);
  }

  void construct_pmmm(CalcForce &cf, int num_pp, int num_mm, double total_charge=0.0)
  {
    std::vector<CellMethodModule::CellRange> cell_range;
    std::vector<Position> cell_center;
    PMMMPreparator::set_cell_range(MPI_COMM_WORLD, mpi_comm_short, num_pp, mpi_comm_long,
				   cellindex, cellsubset, ope,
				   cell_center, cell_range);
    std::vector<int> m_boundary;
    std::vector<int> mm_targetid;
    std::vector<int> pp_targetid;
    PMMMPreparator::make_mm_info(longrangeparameter.order, num_pp, num_mm,
				 m_boundary, pp_targetid, mm_targetid);
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    if(pid==0){
      printf("m_boundary[%d]",m_boundary.size());
      for(int m=0;m<m_boundary.size();m++){
	printf(" %d",m_boundary[m]);
      }
      printf("\n");
    }
    PMMMPreparator::initialize_PMMM_Interaction(MPI_COMM_WORLD, mpi_comm_short, mpi_comm_long,
						cellindex, longrangeparameter, ope,
						m_boundary, pp_targetid, mm_targetid,
						cell_center, cell_range,
						total_charge,
						cf.pmmmlocal, cf.pmmmlongrange);
    

  }
#endif  // CPPMD_ENABLE_PMMM
}
