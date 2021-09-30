#include <cstdlib>
#include <time.h>
#include "MPIParallel.h"
#include <set>
#include <iostream>
#include <cstring>
#include <fstream>

/*
double MPI_Wtime(){
  return (double)clock()/(double)CLOCKS_PER_SEC;
}
*/

class testCommunicator {
public:
  int unit_id;
  size_t psize;

  ParticleArray& myparticle;
  std::vector<TypeRange>& typerangearray;
  std::vector<BondList>& bondlistarray;
  std::vector<int>& setid;
  ForceArray myforce;

  ParticleArray ghost;
  std::vector<ParticleRange> targetparticlerange;
  std::vector<TypeRange> targettyperange;
  std::vector<BondList> targetbond;
  std::vector<int> targetsetid;
  std::vector<int> recvsetid;
  ForceArray sendforce;

  MPICommunicator& communicator;

  testCommunicator(ParticleArray& mp, std::vector<TypeRange>& tra, 
                   std::vector<BondList>& ba, std::vector<int>& sid,
                   MPICommunicator& comm,
                   int uid, size_t ps )
    : unit_id(uid), psize(ps),
      myparticle(mp), typerangearray(tra), bondlistarray(ba), setid(sid),
      ghost(0), targetparticlerange(0), targettyperange(0), targetbond(0), targetsetid(0),recvsetid(0),
      communicator(comm)
  {
    targetparticlerange.resize(communicator.rpsf_number_of_target);
    targettyperange.resize(communicator.number_of_target_set);
    targetbond.resize(communicator.rpsf_number_of_target);
    targetsetid = communicator.alltargetsetid;
    std::cout << " number_of_target_set " << communicator.number_of_target_set << std::endl;
    if(unit_id==0){
      std::cout << "unit " << unit_id << " taget set";
      for(std::vector<int>::iterator it = targetsetid.begin();
          it != targetsetid.end(); ++it){
        std::cout << " " << *it;
      }
      std::cout << std::endl;
    }
    myforce.resize(myparticle.size());
  }


  void start(int loop, double& alltime, double& rectime)
  {
    for(int lc=0;lc<loop;lc++){
      ghost.clear();
      targetparticlerange.clear();
      targettyperange.clear();
      targetbond.clear();
      recvsetid.clear();
      communicator.exchangeParticleArraysubset(myparticle, 
                                               typerangearray, bondlistarray,
                                               ghost, targetparticlerange, recvsetid,
                                               targettyperange, targetbond, 
                                               alltime,rectime);
      sendforce.clear();
      sendforce.resize(ghost.size());
      // force operation
      SpaceVector<double> f(0.0,0.0,0.0);
      for(int j=0;j<ghost.size();j++){
        sendforce[j] = ghost[j].position + f;
      }
      communicator.exchangeForceArraysubset(sendforce,targetparticlerange,
                                            myforce,alltime,rectime);
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

  mpi_return = MPI_Init(&argc, &argv);
  if(mpi_return!=MPI_SUCCESS){
    exit(1);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&n);

  if(n<2){
    std::cout << "np<2 " << std::endl;
    MPI_Finalize();
    exit(1);
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  for(int i=0;i<n;i++){
    int id=i;
    int rank=i;
    unitid_to_rank.insert(std::pair<int,int>(id,rank));
  }

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


  if(argc>1){
    m = std::atol(argv[1]);
    if(m<1){
      m = 1;
    }
    if(argc>2){
      ghost_len = std::atol(argv[2]);
      if(ghost_len<1){
        ghost_len = 1;
      }
      if(argc>3){
        repeat = std::atol(argv[3]);
      }
    }
  }

  int p=4*m;
  for(;p*p*p<n;p*=2);
  num_set=p*p*p;
  set_per_unit=num_set/n;
  
  int set_in_unit[n];
  int dv[3]={1,1,1};
  {
    int pn=1;
    int divdir=2;
    set_in_unit[0]=num_set;
    while(pn<n){
      int i;
      dv[divdir]*=2;
      for(i=0;i<pn;i++){
        set_in_unit[pn+i] = set_in_unit[i]>>1;
        set_in_unit[i] -= set_in_unit[pn+i];
        if(pn+i==n-1){
          i++;
          break;
        }
      }
      if(divdir==0){
        divdir=2;
      }else{
        divdir--;
      }
      pn+=i;
      if(myrank==0){
        std::cout << " set in unit";
        for(i=0;i<pn;i++){
          std::cout << " " << set_in_unit[i];
        }
        std::cout << std::endl;
      }
    }
    if(myrank==0){
      std::cout << "div " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
    }
  }
  std::set<int> cell[dv[2]][dv[1]][dv[0]];
  {
    for(int cz=0;cz<dv[2];cz++){
      for(int cy=0;cy<dv[1];cy++){
        for(int cx=0;cx<dv[0];cx++){
          for(int ix=0;ix<p/dv[0];ix++){
            int x = cx*p/dv[0]+ix;
            for(int iy=0;iy<p/dv[1];iy++){
              int y = cy*p/dv[1]+iy;
              for(int iz=0;iz<p/dv[2];iz++){
                int z = cz*p/dv[2]+iz;
                cell[cz][cy][cx].insert(x+y*p+z*p*p);
              }
            }
          }
          if(myrank==0){
            std::cout << "cell " << cx << " " << cy << " " << cz << " :";
            for(std::set<int>::iterator it = cell[cz][cy][cx].begin();
              it != cell[cz][cy][cx].end(); ++it){
              std::cout << " " << (*it);
            }
            std::cout << std::endl;
          }
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  setidarray.resize(n);
  {
    int u=0;
    for(int cz=0;cz<dv[2];cz++){
      for(int cy=0;cy<dv[1];cy++){
        for(int cx=0;cx<dv[0];cx++){
          for(std::set<int>::iterator it = cell[cz][cy][cx].begin();
              it != cell[cz][cy][cx].end(); ++it){
            setidarray[u].push_back(*it);
            setid_to_unitid.insert(std::pair<int,int>(*it,u));
          }
          u++;
        }
      }
    }
  }
  int num_ghost = (ghost_len*2+1)*(ghost_len*2+1)*(ghost_len*2+1)-1;
  SpaceVector<int> relsid[num_ghost];
  {
    int i=0;
    for(int z=-ghost_len;z<=ghost_len;z++){
      for(int y=-ghost_len;y<=ghost_len;y++){
        for(int x=-ghost_len;x<=ghost_len;x++){
          if(!((x==0)&&(y==0)&&(z==0))){
            relsid[i] = SpaceVector<int>(x,y,z);
            i++;
          }
        }
      }
    }
    if(myrank==0){
      std::cout << " relative taget id";
      for(int i=0;i<num_ghost;i++){
        std::cout << " " << relsid[i];
      }
      std::cout << std::endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  std::vector< SpaceVector<double> > setid_to_cellposition(num_set);
  {
    for(int setid=0;setid<num_set;setid++){
      setid_to_cellposition[setid].x = (setid%p);
      setid_to_cellposition[setid].y = ((setid/p)%p);
      setid_to_cellposition[setid].z = (setid/(p*p));
    }
  }

  std::vector<int> send_target_unit(0);
  std::vector<int> recv_target_unit(0);
  std::map<int,int> target_to_index;
  std::vector< std::vector<int> > send_target_set(0);
  std::vector< std::vector<int> > recv_target_set(0);
  {
    std::vector<int> t_unit(0);
    std::map<int,int> t_to_index;
    std::vector< std::set<int> > send_set(0);
    std::vector< std::set<int> > recv_set(0);
    std::vector<int>& myset = setidarray[unit_identifier];
    for(std::vector<int>::iterator ms = myset.begin();ms!=myset.end();++ms){
      int setid=(*ms);
      for(int i=0;i<num_ghost;i++){
        SpaceVector<int> rel=relsid[i];
        rel.x += (setid%p);
        rel.y += ((setid/p)%p);
        rel.z += (setid/(p*p));
        if(rel.x>=p)rel.x-=p;
        if(rel.x<0)rel.x+=p;
        if(rel.y>=p)rel.y-=p;
        if(rel.y<0)rel.y+=p;
        if(rel.z>=p)rel.z-=p;
        if(rel.z<0)rel.z+=p;
        int tset = rel.x+rel.y*p+rel.z*p*p;
        if(tset<0)tset+=num_set;
        if(tset>=num_set)tset-=num_set;
        int target_uid = setid_to_unitid[tset];
        if(target_uid!=unit_identifier){
          int unitindex=0;
          for(unitindex=0;unitindex<t_unit.size();unitindex++){
            if(target_uid==t_unit[unitindex])break;
          }
          if(unitindex==t_unit.size()){
            t_to_index.insert(std::pair<int,int>(target_uid,unitindex));
            t_unit.push_back(target_uid);
            std::set<int> newsend;
            newsend.insert(setid);
            send_set.push_back(newsend);
            std::set<int> newrecv;
            newrecv.insert(tset);
            recv_set.push_back(newrecv);
          }else{
            send_set[unitindex].insert(setid);
            recv_set[unitindex].insert(tset);
          }
        }
      }
    }
    
    int index=0;
    for(std::map<int,int>::iterator tt = t_to_index.begin(); 
        tt!=t_to_index.end(); ++tt){
      int oldindex = (*tt).second;
      int tid = (*tt).first;
      send_target_unit.push_back(tid);
      recv_target_unit.push_back(tid);
      target_to_index.insert(std::pair<int,int>(tid,index));
      std::vector<int> setvec(0);
      for(std::set<int>::iterator it = send_set[oldindex].begin();
          it!=send_set[oldindex].end();++it){
        setvec.push_back(*it);
      }
      send_target_set.push_back(setvec);
      setvec.clear();
      for(std::set<int>::iterator it = recv_set[oldindex].begin();
          it!=recv_set[oldindex].end();++it){
        setvec.push_back(*it);
      }
      recv_target_set.push_back(setvec);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for(int r=0;r<n;r++){
      //      if(r==myrank){
      if((r==myrank)&&(r==0)){
      //      if((r==myranl)&&(r==1)){
        std::cout << "number of send target " << send_set.size() << std::endl;
        for(std::map<int,int>::iterator tt = t_to_index.begin();
            tt!=t_to_index.end(); ++tt){
          int ti = (*tt).second;
          std::cout << "set for target unit " << t_unit[ti] << " :";
          std::set<int>& sset = send_set[ti];
          for(std::set<int>::iterator it = sset.begin(); it!=sset.end();++it){
            std::cout << " " << *it;
          }
          std::cout << std::endl;
        }
        std::cout << "number of recv target " << recv_set.size() << std::endl;
        for(std::map<int,int>::iterator tt = t_to_index.begin();
            tt!=t_to_index.end(); ++tt){
          int ti = (*tt).second;
          std::cout << "set from target unit " << t_unit[ti] << " :";
          std::set<int>& rset = recv_set[ti];
          for(std::set<int>::iterator it = rset.begin(); it!=rset.end();++it){
            std::cout << " " << *it;
          }
          std::cout << std::endl;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    
  }
  {
    if(myrank==0){
      for(int st=0;st<send_target_set.size();st++){
        std::vector<int>& sset = send_target_set[st];
        std::cout << "set for target unit " << send_target_unit[st] << " :";
        for(std::vector<int>::iterator it = sset.begin(); it!=sset.end();++it){
          std::cout << " " << *it;
        }
        std::cout << std::endl;
      }
      for(int st=0;st<send_target_set.size();st++){
        std::vector<int>& rset = recv_target_set[st];
        std::cout << "set from target unit " << recv_target_unit[st] << " :";
        for(std::vector<int>::iterator it = rset.begin(); it!=rset.end();++it){
          std::cout << " " << *it;
        }
        std::cout << std::endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  {
    for(int r=0;r<n;r++){
      if(r==myrank){
        std::cout << "unit " << unit_identifier << " set";
        for(std::vector<int>::iterator it = setidarray[unit_identifier].begin();
            it != setidarray[unit_identifier].end(); ++it){
          std::cout << " " << *it;
        }
        std::cout << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  {
    /*
      even node SendParticleReceiveForce
      odd node  ReceiveParticleSendForce
    */
    full_sprf_target_id.resize(n);
    full_rpsf_target_id.resize(n);
    for(int i=0;i<n;i+=2){
      for(int j=1;j<n;j+=2){
        full_sprf_target_id[i].push_back(j);
      }
    }

    for(int i=0;i<n;i++){
      for(std::vector<int>::iterator it = full_sprf_target_id[i].begin();
          it != full_sprf_target_id[i].end(); ++it){
        int j = (*it);
        full_rpsf_target_id[j].push_back(i);
      }
    }
  }

  {
    for(std::vector<int>::iterator it = full_rpsf_target_id[unit_identifier].begin();
        it != full_rpsf_target_id[unit_identifier].end(); ++it){
      int t = (*it);
      targetsetid.push_back(setidarray[t]);
    }
  }
  int psize = num_set*8/n;
  
  ParticleArray myparticle(psize);
  {
    std::vector<int>& myset = setidarray[unit_identifier];
    int size_per_set = psize/myset.size();
    for(int i=0;i<psize;i++){
      int setindex = i/size_per_set;
      int setid = myset[setindex];
      myparticle[i].position = setid_to_cellposition[setid];
    }
    
  }

  std::cout << "size " << psize << std::endl;
  
  std::vector<TypeRange> typerangearray(set_per_unit);
  {
    if(myrank==0){
      std::cout << "typerangearray.all";
    }
    for(int s=0;s<set_per_unit;s++){
      typerangearray[s].coulomb.begin = psize*(s)/set_per_unit;
      typerangearray[s].coulomb.end = psize*(s+1)/set_per_unit;
      typerangearray[s].all.begin = psize*(s)/set_per_unit;
      typerangearray[s].all.end = psize*(s+1)/set_per_unit;
      if(myrank==0){
        std::cout << " " << typerangearray[s].all.begin << "--" << typerangearray[s].all.end ;
      }
    }
    if(myrank==0){
      std::cout << std::endl;
    }
  }
  

  std::vector<BondList> bondlistarray(1);

  int num_rpsf_target = full_rpsf_target_id[unit_identifier].size();
  //  ParticleArray ghost(psize*num_rpsf_target);
  ParticleArray ghost(0);
  std::vector<ParticleRange> targetparticlerange(num_rpsf_target);
  for(size_t t=0;t<full_rpsf_target_id[unit_identifier].size();t++){
    targetparticlerange[t].begin = 0;
    targetparticlerange[t].end = 0;
  }
  std::vector<TypeRange> targettyperange(num_rpsf_target);
  std::vector<BondList> targetbond(num_rpsf_target);


  /*
  MPICommunicator communicator(unitid_to_rank, 
                               full_sprf_target_id[unit_identifier], 
                               send_target_set,
                               setidarray[unit_identifier],
                               full_rpsf_target_id[unit_identifier], 
                               targetsetid,
                               psize, n, unit_identifier);
  */
  MPICommunicator communicator(unitid_to_rank, 
                               send_target_unit,
                               send_target_set,
                               setidarray[unit_identifier],
                               recv_target_unit,
                               recv_target_set,
                               psize, n, unit_identifier);

  std::cout << " const comm " << std::endl;


  testCommunicator testcomm(myparticle, typerangearray, bondlistarray,
                            setidarray[unit_identifier],
                            communicator, unit_identifier, psize);
  double alltime;
  double rectime;


  MPI_Barrier(MPI_COMM_WORLD);

  /*
  for(int l=0;l<repeat;l++){
    ghost.clear();
    targetparticlerange.clear();
    targettyperange.clear();
    targetbond.clear();
    communicator.exchangeParticleArray(myparticle, 
                                       typerangearray, bondlistarray,
                                       ghost, targetparticlerange,
                                       targettyperange, targetbond,
                                       alltime,rectime);
  }
  */

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

  /*
  for(std::map<int,int>::iterator it = unitid_to_rank.begin();
      it != unitid_to_rank.end(); ++it) {
    MPI_Barrier(MPI_COMM_WORLD);
    if((*it).second==myrank){ 
      std::cout << "unit " << unit_identifier << std::endl;
      std::cout << "id "<< (*it).first;
      std::cout << "typerange coulomb.end" ;
      for(std::vector<TypeRange>::iterator trit = typerangearray.begin();
          trit != typerangearray.end(); ++trit){
        std::cout << " " << (*trit).coulomb.end;
      }
      std::cout << std::endl;
      std::cout << " rank "<< (*it).second << std::endl;
      for(size_t t = 0; t<full_rpsf_target_id[unit_identifier].size(); t++){
        std::cout << testcomm.ghost[t*psize].position << testcomm.ghost[(t+1)*psize-1].position;
        std::cout << testcomm.targetparticlerange[t].begin << " -- " << testcomm.targetparticlerange[t].end;
        std::cout << std::endl;
      }
      std::cout << "recvset";
      for(std::vector<int>::iterator it = testcomm.recvsetid.begin();
          it != testcomm.recvsetid.end(); ++it){
        std::cout << " " << *it ;
      }
      std::cout << std::endl;
      std::cout << "recvtyperange coulomb" ;
      for(std::vector<TypeRange>::iterator it = testcomm.targettyperange.begin();
          it != testcomm.targettyperange.end(); ++it){
        std::cout << " " << (*it).coulomb.begin << "- " << (*it).coulomb.end;
      }
      std::cout << std::endl;
    }
  }
  */
  
  char filename[256];
  snprintf(filename,255,"myparticle.%d",myrank);
  std::cout << filename << std::endl;

  std::ofstream file;
  file.open(filename);
  int setindex=0;
  int bound = testcomm.typerangearray[setindex].all.end;
  int setid = testcomm.setid[setindex];
  for(int i=0;i<testcomm.myparticle.size();i++){
    if(i>=bound){
      setindex++;
      setid = testcomm.setid[setindex];
      bound = testcomm.typerangearray[setindex].all.end;
    }
    file << testcomm.myparticle[i].position << " " << setid << std::endl;
  }
  file.close();


  char ghostfilename[256];
  snprintf(ghostfilename,255,"ghost.%d",myrank);
  //  std::cout << ghostfilename << std::endl;
  std::ofstream ghostfile;
  ghostfile.open(ghostfilename);
  setindex=0;
  bound = testcomm.targettyperange[setindex].all.end;
  setid = testcomm.recvsetid[setindex];
  for(int i=0;i<testcomm.ghost.size();i++){
    if(i>=bound){
      setindex++;
      setid = testcomm.recvsetid[setindex];
      bound = testcomm.targettyperange[setindex].all.end;
    }
    ghostfile << testcomm.ghost[i].position << " " << setid << std::endl;
  }
  ghostfile.close();

  char forcefilename[256];
  snprintf(forcefilename,255,"force.%d",myrank);
  //  std::cout << forcefilename << std::endl;
  std::ofstream forcefile;
  forcefile.open(forcefilename);
  setindex=0;
  bound = testcomm.typerangearray[setindex].all.end;
  setid = testcomm.setid[setindex];
  for(int i=0;i<testcomm.myforce.size();i++){
    if(i>=bound){
      setindex++;
      setid = testcomm.setid[setindex];
      bound = testcomm.typerangearray[setindex].all.end;
    }
    forcefile << testcomm.myforce[i] << " " << setid << std::endl;
  }
  forcefile.close();

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

};
