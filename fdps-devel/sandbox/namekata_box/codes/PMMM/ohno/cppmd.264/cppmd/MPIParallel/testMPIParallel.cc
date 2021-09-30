#include <cstdlib>
#include <time.h>
#include "MPIParallel.h"
#include <set>
#include <iostream>
#include <cstring>
#include <fstream>
#include "CubicCell.h"
#include "CovalentBondInfo.h"

/*
double MPI_Wtime(){
  return (double)clock()/(double)CLOCKS_PER_SEC;
}
*/

// communication tester
class testCommunicator {
public:
  int unit_id;
  size_t psize;

  ParticleArray& myparticle;
  std::vector<TypeRange>& typerangearray;
  std::vector<CovalentBondInfo::BondList>& bondlistarray;
  std::vector<int>& setid;
  ForceArray myforce;

  ParticleArray ghost;
  std::vector<ParticleRange> targetparticlerange;
  std::vector<TypeRange> targettyperange;
  std::vector<CovalentBondInfo::BondList> targetbond;
  std::vector<int> targetsetid;
  std::vector<int> recvsetid;
  std::map<int,int> recvsetid_to_index;
  ForceArray sendforce;

  MPICommunicator& communicator;

  SpaceVector<double> boxsize;

  testCommunicator(ParticleArray& mp, std::vector<TypeRange>& tra, 
                   std::vector<CovalentBondInfo::BondList>& ba, std::vector<int>& sid,
                   MPICommunicator& comm,
                   int uid, size_t ps,
                   SpaceVector<double> bsize)
    : unit_id(uid), psize(ps),
      myparticle(mp), typerangearray(tra), bondlistarray(ba), setid(sid),
      ghost(0), targetparticlerange(0), targettyperange(0), targetbond(0), targetsetid(0),recvsetid(0),
      communicator(comm),
      boxsize(bsize)
  {
    targetparticlerange.resize(communicator.rpsf_number_of_target);
    targettyperange.resize(communicator.number_of_target_set);
    targetbond.resize(communicator.rpsf_number_of_target);
    targetsetid = communicator.alltargetsetid;
    std::cout << " number_of_target_set " << communicator.number_of_target_set << std::endl;
    recvsetid.resize(communicator.number_of_target_set);
    /*
    if(unit_id==0){
      std::cout << "unit " << unit_id << " taget set";
      for(std::vector<int>::iterator it = targetsetid.begin();
          it != targetsetid.end(); ++it){
        std::cout << " " << *it;
      }
      std::cout << std::endl;
    }
    */
    myforce.resize(myparticle.size());
  }


  void start(int loop, double& alltime, double& rectime)
  {
    for(int lc=0;lc<loop;lc++){
      std::cout << " loop " << lc << std::endl;

      /*
      communicator.move_particle(myparticle,
                                 typerangearray);
      */
      
      ghost.clear();
      targetparticlerange.clear();
      targettyperange.clear();
      targetbond.clear();
      recvsetid.clear();
      recvsetid_to_index.clear();
      communicator.exchangeParticleArraysubset(myparticle, 
                                               typerangearray, bondlistarray,
                                               ghost, targetparticlerange, 
                                               recvsetid, recvsetid_to_index,
                                               targettyperange, targetbond);
      std::cout << "exchangeParticleArraysubset done" << std::endl;
      sendforce.clear();
      sendforce.resize(ghost.size());
      // force operation
      SpaceVector<double> f(0.0,0.0,0.0);
      for(int j=0;j<ghost.size();j++){
        sendforce[j] = ghost[j].position + f;
      }
      communicator.exchangeForceArraysubset(sendforce,targetparticlerange,
                                            myforce);
      //      myparticle[0].position.z = boxsize.z - myparticle[0].position.z;
      int mps=0;
      for(int s=0;s<typerangearray.size();s++){
        mps += (typerangearray[s].end-typerangearray[s].begin);
      }
    std::cout << "self size";
    std::cout << " " << myparticle.size() <<  " " << mps; 
    std::cout << " " << typerangearray.size();
    std::cout << " " << bondlistarray.size();
    std::cout << std::endl;
    std::cout << "ghost size";
    std::cout << " " << ghost.size() ;
    std::cout << " " << targetparticlerange.size() ;
    std::cout << " " <<  targettyperange.size();
    std::cout << " " <<  targetbond.size();
    std::cout << " " <<  recvsetid.size();
    std::cout << " " <<  recvsetid_to_index.size();
    std::cout << std::endl;
    }
  }
};


int
main(int argc, char **argv)
{
  int n;
  int myrank;
  int m=1;
  int repeat=1;
  int unit_identifier;
  int mpi_return;
  std::map<int,int> unitid_to_rank;
  std::vector< std::vector<int> > full_sprf_target_id;
  std::vector< std::vector<int> > full_rpsf_target_id;

  int set_per_unit=1;
  int num_set;
  std::map<int,int> setid_to_unitid;
  std::vector< std::vector<int> > setidarray;
  std::vector< std::vector<int> > targetsetid;

  int ghost_len=1;
  int celldiv = 4;

  mpi_return = MPI_Init(&argc, &argv);
  if(mpi_return!=MPI_SUCCESS){
    std::cout << "MPI_Init return " << mpi_return << std::endl;
    exit(1);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&n);

  if(n<2){
    std::cout << "np<2 " << std::endl;
    MPI_Finalize();
    exit(1);
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  
  int comm_pat=0;

  if(argc>1){
    int cp = std::atol(argv[1]);
    if(cp<1)cp=1;
    celldiv = 1<<cp;
    if(argc>2){
      ghost_len = std::atol(argv[2]);
      if(ghost_len<1){
        ghost_len = 1;
      }
      if(argc>3){
        repeat = std::atol(argv[3]);
        if(argc>4){
          comm_pat = std::atol(argv[4]);
        }
      }
    }
  }
  
  CommPattern commpattern = NearestXYZ;
  if(comm_pat!=0){
    commpattern = Direct;
  }
  
  // check n = 2^m
  int np;             // n is np-th power of 2,  n = 2^np
  if(power_of_two(n,np)){
    if(myrank==0)std::cout << n << " = 2^" << np << std::endl;
  }else{
    if(myrank==0)std::cout << n << " != 2^" << np << std::endl;
    MPI_Finalize();
    exit(1);
  }


  // map unit_id to rank
  for(int i=0;i<n;i++){
    int id=i;
    int rank=i;
    unitid_to_rank.insert(std::pair<int,int>(id,rank));
  }

  // check myrank exist in unitid_to_rank
  {
    unit_identifier = -1;
    std::map<int,int>::iterator it = unitid_to_rank.begin();
    while( (it != unitid_to_rank.end())
           && ((*it).second!=myrank) ) {
      ++it;
    }
    if(it != unitid_to_rank.end()){
      unit_identifier = (*it).first;
    }else{
      // bad unitid_to_rank
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }

  SpaceVector<double> boxsize(64.0,64.0,64.0);


  // calc x,y,z divide to node
  int sp[3];    // separated power np
  split_power_of_two(np,sp);
  SpaceVector<int> div;         // number of node in x,y,z axis
  div[0] = 1<<sp[0];
  div[1] = 1<<sp[1];
  div[2] = 1<<sp[2];
  int num_node = div[0]*div[1]*div[2];
  if(myrank==0)std::cout << num_node << " = " << div[0] << "*" << div[1] << "*" << div[2] << std::endl;

  // force cell in one axis >= div[2]
  if(celldiv<div[2]){
    if(myrank==0)std::cout << "cell div < " << div[2] << std::endl;
    celldiv=div[2];
  }

  CellIndex cellindex;

  SpaceVector<double> cellmargin(1.0,1.0,1.0);

  double cutoff = 9.0;
  PostProcess postprocess;

  cellindex.set_cell_parameter(boxsize,
                               SpaceVector<int>(celldiv,celldiv,celldiv),
                               cellmargin,
                               cutoff);
  double volume = cellindex.calc_volume();
  int num_cell = cellindex.generate_cell_index();
  postprocess.celldiv = cellindex.celldiv;
  postprocess.cell_position = cellindex.cell_position;

  // force ghost_len < celldiv/2
  if(ghost_len>=(celldiv>>1))ghost_len=(celldiv>>1)-1;

  if(myrank==0)std::cout << "number of cell " << celldiv << "^3 " << num_cell << std::endl;

  if(!cellindex.distribute_cell(div)){
    std::cout << "distribute_cell false " << std::endl;
  }






  // set target set

  // construct postprocess
  postprocess.margin = 1.0;

  postprocess.margin = cellindex.margin;
  postprocess.fine_margin = cellindex.margin;
  postprocess.boxsize = cellindex.boxsize;
  postprocess.cellsize = cellindex.cellsize;
  postprocess.invcellsize.x = 1.0/postprocess.cellsize.x;
  postprocess.invcellsize.y = 1.0/postprocess.cellsize.y;
  postprocess.invcellsize.z = 1.0/postprocess.cellsize.z;
  double volume_margin = postprocess.calc_volume_margin();

  num_set = num_cell/num_node;
  int num_total_set = num_set;
  int max_id = num_total_set;
  int max_particle_in_cell = 8;
  int max_num_particle = max_particle_in_cell*num_total_set;
  int num_spu = num_set;
  //    int max_num_pps = max_particle_in_cell; // depend on simulation
  
  std::vector<int> particle_setid(num_set);
  std::vector<TypeRange> typerangearray(num_set);
  int node_id = myrank;
  int size_per_set = 1;
  int max_size_per_set = 8;
  int psize = num_set*max_size_per_set;

  ParticleArray myparticle(psize*2);
  {
    for(int s=0;s<num_set;s++){
      typerangearray[s].begin = s*max_size_per_set;
      typerangearray[s].end = typerangearray[s].begin;
      int ioffset = typerangearray[s].begin;
      particle_setid[s] = cellindex.cell_id_in_node[node_id][s];
      int setid = particle_setid[s];
      for(int i=0;i<size_per_set;i++){
        myparticle[i+ioffset].position.x = boxsize.x/celldiv*cellindex.cell_position[setid].x;
        myparticle[i+ioffset].position.y = boxsize.y/celldiv*cellindex.cell_position[setid].y;
        myparticle[i+ioffset].position.z = boxsize.z/celldiv*cellindex.cell_position[setid].z;
        myparticle[i+ioffset].atomid = i+particle_setid[s]*size_per_set;
        std::cout << myparticle[i+ioffset].atomid << " " << i+ioffset << myparticle[i+ioffset].position << std::endl;
        typerangearray[s].end++;
        typerangearray[s].coulomb.end++;
      }
    }
  }
  if(myrank==0){
    for(int s=0;s<num_set;s++){
      std::cout << "set " << s << ":" << particle_setid[s] << "[" << typerangearray[s].begin << "," << typerangearray[s].end << ")" << std::endl;
    }
  }

  std::cout << "size " << psize << std::endl;

  CellIndexType cellindextype = HalfShell;
  CellSubset cellsubset(cellindex,node_id);
  cellsubset.makefulltarget(particle_setid,cellindextype);

  postprocess.shift_cell=cellsubset.shift_cell_all;
  
  std::vector< std::vector<int> > targetset;
  std::vector<int> sendparticle_target_id;
  std::vector< std::vector<int> > sendparticle_setid;
  std::vector<int> recvparticle_target_id;
  std::vector< std::vector<int> > recvparticle_setid;
  std::vector<int> move_target_id;
  std::vector< std::vector<int> > move_setid;

  cellsubset.makecommtarget(sendparticle_target_id, 
                            sendparticle_setid,
                            recvparticle_target_id,
                            recvparticle_setid,
                            move_target_id, move_setid);
  //  if(myrank==0)
    {
    std::cout << "sendparticle_target_id " << sendparticle_target_id.size();
    for(int t=0;t<sendparticle_target_id.size();t++){
      std::cout << " " << sendparticle_target_id[t];
    }
    std::cout << std::endl;
    std::cout << "recvparticle_target_id " << recvparticle_target_id.size();
    for(int t=0;t<recvparticle_target_id.size();t++){
      std::cout << " " << recvparticle_target_id[t];
    }
    std::cout << std::endl;
  }
  targetset = cellsubset.recv_target_cell;
  /*
  if(myrank==0){
    std::cout << "targetset.size() " << targetset.size() << std::endl;
    for(int t=0;t<targetset.size();t++){
      cellsubset.dump_cell_target_recv(t);
    }
  }
  */
  MPI_Barrier(MPI_COMM_WORLD);


  


  SpaceVector<int> node_div = div;

  // prepare ghost
  ParticleArray ghost(0);
  std::vector<CovalentBondInfo::BondList> bondlistarray(1);

  // construct communicator
  MPICommunicator communicator(unitid_to_rank,
                               sendparticle_target_id,
                               sendparticle_setid,
                               particle_setid,
                               recvparticle_target_id,
                               recvparticle_setid,
                               move_target_id,
                               move_setid,
                               postprocess,
                               max_particle_in_cell,
                               max_id, node_id,
                               node_div,
                               cellindex.celldiv,
                               cellindex.move_range,
                               commpattern
                               );

  std::cout << " const comm " << std::endl;

  if(myrank==0){
    std::cout << "communicator.rpsf_number_of_target " << communicator.rpsf_number_of_target << std::endl;
  }
  

  // construct tester for communicator
  testCommunicator testcomm(myparticle, typerangearray, bondlistarray,
                            particle_setid,
                            communicator, unit_identifier, psize,
                            boxsize);


  double alltime;
  double rectime;


  MPI_Barrier(MPI_COMM_WORLD);

  // test communication
  testcomm.start(repeat,alltime,rectime);
  
  double allsum=0.0;
  double recsum=0.0;
  double allmin, allmax;
  double recmin, recmax;

  MPI_Barrier(MPI_COMM_WORLD);
  /*
    all reduce alltime -> allsum
    all reduce rectime -> recsum;
   */

  MPI_Reduce(&alltime,&allsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&alltime,&allmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&alltime,&allmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&rectime,&recsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&rectime,&recmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&rectime,&recmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  if(myrank==0){
    //    std::cout << "sum " << allsum << " " << recsum << std::endl;
    std::cout << "          all    receive" << std::endl;
    std::cout << "min. " << allmin << " " << recmin << std::endl;
    std::cout << "max. " << allmax << " " << recmax << std::endl;
    std::cout << "avr. " << allsum/n << " " << recsum/n << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for(int r=0;r<n;r++){
    if(r==myrank){
      std::cout << "testcomm.ghost size " <<  testcomm.ghost.size() << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }


  char filename[256];
  snprintf(filename,255,"myparticle.%d",myrank);
  std::cout << filename << std::endl;

  std::ofstream file;
  file.open(filename);
  int setindex=0;
  int bound = testcomm.typerangearray[setindex].end;
  int setid = testcomm.setid[setindex];
  for(int s=0;s<num_set;s++){
    int setid = particle_setid[s];
    for(int i=testcomm.typerangearray[s].begin;i<testcomm.typerangearray[s].end;i++){
      file << testcomm.myparticle[i].position << " " << testcomm.myparticle[i].atomid << " "  << setid << std::endl;
    }
  }
  file.close();


  char ghostfilename[256];
  snprintf(ghostfilename,255,"ghost.%d",myrank);
  //  std::cout << ghostfilename << std::endl;
  std::ofstream ghostfile;
  ghostfile.open(ghostfilename);
  setindex=0;
  bound = testcomm.targettyperange[setindex].end;
  setid = testcomm.recvsetid[setindex];
  for(int i=0;i<testcomm.ghost.size();i++){
    if(i>=bound){
      setindex++;
      setid = testcomm.recvsetid[setindex];
      bound = testcomm.targettyperange[setindex].end;
    }
    ghostfile << testcomm.ghost[i].position << " " << testcomm.ghost[i].atomid << " " << setid << std::endl;
  }
  ghostfile.close();

  char forcefilename[256];
  snprintf(forcefilename,255,"force.%d",myrank);
  //  std::cout << forcefilename << std::endl;
  std::ofstream forcefile;
  forcefile.open(forcefilename);
  for(int s=0;s<num_set;s++){
    int setid = particle_setid[s];
    for(int i=testcomm.typerangearray[s].begin;i<testcomm.typerangearray[s].end;i++){
      forcefile << testcomm.myforce[i] << " " << setid << std::endl;
    }
  }
  setindex=0;
  bound = testcomm.typerangearray[setindex].end;
  setid = testcomm.setid[setindex];
  for(int i=0;i<testcomm.myforce.size();i++){
    forcefile << testcomm.myforce[i] << std::endl;
  }
  forcefile.close();

  MPI_Barrier(MPI_COMM_WORLD);

  /*
  if(myrank==0){
    std::cout << "receive cell";
    for(int s=0;s<testcomm.recvsetid.size();s++){
      std::cout << " " << testcomm.recvsetid[s];
    }
    std::cout << std::endl;
  }
  */

  MPI_Finalize();

};
