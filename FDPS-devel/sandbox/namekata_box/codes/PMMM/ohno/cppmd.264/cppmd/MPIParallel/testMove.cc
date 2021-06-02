#include <cstdlib>
#include <time.h>
#include "MPIParallel.h"
#include <set>
#include <iostream>
#include <cstring>
#include <fstream>
#include "CubicCell.h"


class TestMove {
public:
  int unit_id;
  size_t psize;
  std::map<int,int> unitid_to_rank;

  std::vector<MPIMoveInOut> move;
  std::vector<int> move_target_id;
  int move_number_of_target;
  std::vector< std::vector<int> > move_setid;
  std::map<int,int> setid_to_move_target_index;

  ParticleArray all_move_out_particle;
  std::vector<int> all_move_out_set;
  std::vector<PotentialModel> all_move_out_type;
  std::vector<ParticleArray> move_out_particle;
  std::vector< std::vector<int> > move_out_set;
  std::vector< std::vector<PotentialModel> > move_out_type;
  ParticleArray move_in_particle;
  std::vector<int> move_in_set;
  std::vector<PotentialModel> move_in_type;


  TestMove(std::map<int,int>& uidtor,
           std::vector<int> mtid,
           std::vector< std::vector<int> > msid,
           int ps, int id)
    : unit_id(id), psize(ps), unitid_to_rank(uidtor),
      move(), move_target_id(mtid),  move_setid(msid),
      move_out_particle(), move_out_set(), move_out_type()
  {
    construct_sendrecv();
    make_idset_to_index(move_setid,setid_to_move_target_index);
  }

  TestMove(TestMove& tm)
    : unit_id(tm.unit_id), psize(tm.psize), unitid_to_rank(tm.unitid_to_rank),
      move(tm.move), move_target_id(tm.move_target_id),  move_setid(tm.move_setid),
      move_out_particle(tm.move_out_particle), move_out_set(tm.move_out_set),
      move_out_type(tm.move_out_type)
  {
    construct_sendrecv();
    make_idset_to_index(move_setid,setid_to_move_target_index);
  }
  

  void construct_sendrecv()
  {
    move_number_of_target = move_target_id.size();
    for(int i=0;i<move_number_of_target;i++){
      move.push_back(MPIMoveInOut(unit_id,psize));
      move[i].make_buffer();
      move[i].target_id = move_target_id[i];
      move[i].target_rank = unitid_to_rank[move_target_id[i]];
    }

    std::cout << "move_target_id" ;
    for(int i=0;i<move_number_of_target;i++){
      std::cout << " " << move_target_id[i] ;
    }
    std::cout << std::endl;
    
    all_move_out_particle.reserve(psize*move_number_of_target);
    all_move_out_set.reserve(psize*move_number_of_target);
    all_move_out_type.reserve(psize*move_number_of_target);
    move_out_particle.resize(move_number_of_target);
    move_out_set.reserve(move_number_of_target);
    move_out_type.reserve(move_number_of_target);
    for(int i=0;i<move_number_of_target;i++){
      move_out_particle.reserve(psize);
      move_out_set.reserve(psize);
      move_out_type.reserve(psize);
    }
    move_in_particle.reserve(psize*move_number_of_target);
    move_in_set.reserve(psize*move_number_of_target);
    move_in_type.reserve(psize*move_number_of_target);
    std::cout << "move_in_particle reserve " << psize*move_number_of_target << std::endl;
  }

  void move_particle(ParticleArray& particlearray,
                     std::vector<TypeRange>& typerangearray,
                     int r, int repeat)
  {
    std::cout << move_number_of_target << " " << move_out_particle.size() << std::endl;
    int setid=0;
    for(int t=0;t<move_number_of_target;t++){
      move_out_particle[t].clear();
      move_out_set[t].clear();
      move_out_type[t].clear();
    }
    if(r==(repeat>>1)){
      for(int t=0;t<move_number_of_target;t++){
        for(int i=typerangearray[setid].begin;i<typerangearray[setid].end;i++){
          move_out_particle[t].push_back(particlearray[i]);
          PotentialModel pm = typerangearray[setid].particle_index_to_type(i);
          move_out_type[t].push_back(pm);
        }
        typerangearray[setid].end = typerangearray[setid].begin;
        std::cout << "moveout " << move_out_particle[t].size() << " particle to " << move_target_id[t] << std::endl;
        setid++;
        if(setid>=typerangearray.size())break;
      }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    for(int t=0;t<move_number_of_target;t++){
      move[t].setmoveoutparticle(move_out_particle[t], move_out_type[t], move_out_set[t]);
    }
    for(int t=0;t<move_number_of_target;t++){
      move[t].prepare_receive();
    }
    for(int t=0;t<move_number_of_target;t++){
      move[t].send();
    }

    move_in_set.clear();
    move_in_particle.clear();
    move_in_type.clear();
    
    for(int t=0;t<move_number_of_target;t++){
      move[t].getreceive(move_in_particle, move_in_type, move_in_set);
    }

    for(int t=0;t<move_number_of_target;t++){
      move[t].wait_complete_send();
    }

    std::cout << " move in " << move_in_particle.size() << " particle ";
    for(int i=0;i<move_in_particle.size();i++){
      std::cout << " " << move_in_particle[i].atomid << move_in_particle[i].position;
    }
    std::cout << std::endl;
  }

};

class TestUnit {
public:
  ParticleArray myparticle;
  std::vector<TypeRange> typerangearray;
  
  TestMove testmove;
  
  
  TestUnit(ParticleArray mp, std::vector<TypeRange> tr, TestMove tm)
    : myparticle(mp),
      typerangearray(tr),
      testmove(tm)
  {
  }

  void start(int repeat){
    int r;
    for(int r=0;r<repeat;r++){
      std::cout << "repeat " << r << std::endl;
      testmove.move_particle(myparticle,typerangearray,r,repeat);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
};

int
main(int argc, char **argv)
{
  int n;
  int myrank;
  int npn=1;
  int nsn=1;
  int m=1;
  int repeat=1;
  int unit_identifier;
  int mpi_return;
  std::map<int,int> unitid_to_rank;

  int set_per_unit=1;
  int num_set;
  std::map<int,int> setid_to_unitid;
  std::vector< std::vector<int> > setidarray;
  std::vector< std::vector<int> > targetsetid;


  int ghost_len=1;
  int celldiv = 4;

  mpi_return = MPI_Init(&argc, &argv);
  if(mpi_return!=MPI_SUCCESS){
    exit(1);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&n);

  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  if(argc>1){
    npn = std::atol(argv[1]);
    if(argc>2){
      repeat = std::atol(argv[2]);
    }
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

  MPI_Barrier(MPI_COMM_WORLD);

  setidarray.resize(n);
  targetsetid.resize(n);

  for(int ni=0;ni<n;ni++){
    setidarray[ni].resize(nsn);
    int sid0 = ni*nsn;
    for(int s=0;s<nsn;s++){
      int sid = sid0+s;
      setidarray[ni][s] = sid;
      setid_to_unitid.insert(std::pair<int,int>(s,ni));
    }
  }
  

  int psize = npn;
  
  // set typerangearray
  // all particle is coulomb type
  std::vector<TypeRange> typerangearray(nsn);
  {
    for(int s=0;s<nsn;s++){
      int setid = setidarray[myrank][s];
      int ioffset = s*npn*2;
      typerangearray[s].begin = ioffset;
      typerangearray[s].end = ioffset;
      typerangearray[s].lj.begin = typerangearray[s].begin;
      typerangearray[s].lj.end = typerangearray[s].end;
      typerangearray[s].ljcoulomb.begin = typerangearray[s].end;
      typerangearray[s].ljcoulomb.end = typerangearray[s].end;
      typerangearray[s].coulomb.begin = typerangearray[s].end;
      typerangearray[s].coulomb.end = typerangearray[s].end;
    }
  }
  

  ParticleArray myparticle(psize*nsn*2);
  {
    for(int s=0;s<nsn;s++){
      int setid = setidarray[myrank][s];
      int ioffset = typerangearray[s].begin;
      for(int i=0;i<npn;i++){
        myparticle[i+ioffset].position = Position(double(setid),0.0,0.0);
        myparticle[i+ioffset].atomid = i+s*npn+nsn*npn*myrank;
        std::cout << i+ioffset << " " << myparticle[i+ioffset].atomid << myparticle[i+ioffset].position << std::endl;
        typerangearray[s].end++;
        typerangearray[s].lj.end++;
        typerangearray[s].ljcoulomb.shift(1);
        typerangearray[s].coulomb.shift(1);
      }
    }
  }

  std::cout << "size " << psize << std::endl;

  int num_target = n-1;



  std::vector<int> move_target_id;
  std::vector< std::vector<int> >  move_setid;

  for(int t=0;t<n-1;t++){
    int tid = myrank + t + 1;
    if(tid>=n)tid-=n;
    move_target_id.push_back(tid);
    move_setid.push_back(std::vector<int>(0));
  }

  TestMove testmove(unitid_to_rank,
                    move_target_id,
                    move_setid,
                    psize,myrank);

  TestUnit testunit(myparticle,typerangearray,testmove);

  testunit.start(repeat);

  std::cout << " end " << std::endl;

  MPI_Finalize();

};
