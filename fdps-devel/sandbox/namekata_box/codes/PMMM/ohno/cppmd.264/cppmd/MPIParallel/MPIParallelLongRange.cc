#include "MPIParallelLongRange.h"

MPIReceiverForLong::MPIReceiverForLong(const std::map<int,int>& uidtor,
                                       MPI_Comm long_short_comm,
                                       std::vector<int> rpsf_tid,
                                       std::vector< std::vector<int> > tsid,
                                       int ps, const int maxid, const int id
                                       )
  : unit_id(id), psize(ps), unitid_to_rank(uidtor),
    mpi_comm_long_short(long_short_comm),
    rpsf(), rpsf_target_id(rpsf_tid), 
    targetsetid(tsid), alltargetsetid(),
    bond()
{
  //  std::cout << " construct MPIReceiverForLong  " << unit_id << std::endl;
  construct_recv();
}

MPIReceiverForLong::MPIReceiverForLong(const MPIReceiverForLong& comm)
  : unit_id(comm.unit_id), psize(comm.psize), 
    unitid_to_rank(comm.unitid_to_rank),
    mpi_comm_long_short(comm.mpi_comm_long_short),
    rpsf(), 
    rpsf_target_id(comm.rpsf_target_id),
    targetsetid(comm.targetsetid), alltargetsetid(comm.alltargetsetid),
    bond()
{
  //  std::cout << " copy MPIReceiverForLong  " << unit_id << std::endl;
  construct_recv();
}

MPIReceiverForLong::~MPIReceiverForLong()
{
  rpsf.clear();
}

void
MPIReceiverForLong::construct_recv()
{
  rpsf_number_of_target = rpsf_target_id.size();
  rpsf.clear();
  rpsf.resize(rpsf_number_of_target,MPIReceiveParticleSendForce());
  number_of_target_set = 0;
  for(int i=0;i<rpsf_number_of_target;i++){
    rpsf[i] = MPIReceiveParticleSendForce(unit_id,
                                          psize*targetsetid[i].size(), 
                                          mpi_comm_long_short);
    rpsf[i].make_buffer(targetsetid[i].size());
    rpsf[i].target_id = rpsf_target_id[i];
    rpsf[i].target_rank = unitid_to_rank[rpsf_target_id[i]];
    rpsf[i].recv_requestp = new MPI_Request();
    number_of_target_set += targetsetid[i].size();
    for(std::vector<int>::iterator it = targetsetid[i].begin();
        it!=targetsetid[i].end(); ++it){
      alltargetsetid.push_back(*it);
    }
  }
}

void
MPIReceiverForLong::transferParticle()
{
  //  std::cout << " prepare_receive " << rpsf_number_of_target << std::endl;
  for(int t=0;t<rpsf_number_of_target;t++){
    //    std::cout << " prepare_receive " << rpsf[t].unit_identifier << " from " << rpsf[t].target_rank << ":" << rpsf[t].target_id << " in " << rpsf[t].mpi_comm_short << " " << t << "/" << rpsf_number_of_target << std::endl;
    rpsf[t].prepare_receive();
  }
}

template<class PA>
void
MPIReceiverForLong::getReceiveParticle(PA& targetparticle,
                                       std::vector<ParticleRange>& target_range,
                                       std::vector<TypeRange>& targettyperange,
                                       std::vector<int>& recvsetid,
                                       std::map<int,int>& recvsetid_to_index)
{
  //  target_range.resize(rpsf_number_of_target);
  int t;
#ifdef MIX_MPI_WAIT
  target_range_offset = target_range.size();
  for(t=0;t<rpsf_number_of_target;t++){
    ParticleRange range;
    range.begin = targetparticle.size();
    rpsf[t].getreceive(targetparticle,range,
                       targettyperange, bond, recvsetid,
                       recvsetid_to_index);
    target_range.push_back(range);
  }
#else
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].wait_complete_recv();
  }
  target_range_offset = target_range.size();
  target_range.resize(target_range_offset+rpsf_number_of_target);
  int *poffset, *setoffset, *bondoffset;
  poffset = new int[rpsf_number_of_target+1];
  setoffset = new int[rpsf_number_of_target+1];
  bondoffset = new int[rpsf_number_of_target+1];
  poffset[0] = targetparticle.size();
  setoffset[0] = targettyperange.size();
  bondoffset[0] = bond.size();
    /// serial
  for(t=0;t<rpsf_number_of_target;t++){
    int pnum=0, setnum=0, bondnum=0;
    rpsf[t].getreceive_number(pnum,setnum,bondnum);
    poffset[t+1] = poffset[t] + pnum;
    setoffset[t+1] = setoffset[t] + setnum;
    bondoffset[t+1] = bondoffset[t] + bondnum;
    target_range[target_range_offset+t].begin = poffset[t];
    target_range[target_range_offset+t].end = poffset[t+1];
  }
  targetparticle.resize(poffset[rpsf_number_of_target]);
  targettyperange.resize(setoffset[rpsf_number_of_target]);
  bond.resize(setoffset[rpsf_number_of_target]);
  recvsetid.resize(setoffset[rpsf_number_of_target]);
# ifdef MPIPARALLEL_OPENMP
#  ifdef _OPENMP
#   pragma omp parallel for private(t)
#  endif
# endif
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].getreceive_offset(targetparticle, target_range[target_range_offset+t],
                              targettyperange, bond, recvsetid,
                              setoffset[t],bondoffset[t]);
  }
  /// serial
  for(int i=setoffset[0];i<setoffset[rpsf_number_of_target];i++){
    recvsetid_to_index.insert(std::pair<int,int>(recvsetid[i],i));
  }
  delete [] bondoffset;
  delete [] setoffset;
  delete [] poffset;
#endif
//  printf("number of target_range %d\n",target_range.size());
}

template<class PA>
void
MPIReceiverForLong::getReceiveParticle_onlyPosition(PA& targetparticle)
{

  int t;

#ifdef MIX_MPI_WAIT
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].getreceivepos(targetparticle);
  }
#else
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].wait_complete_recv();
  }
# ifdef MPIPARALLEL_OPENMP
#  ifdef _OPENMP
#pragma omp parallel for private(t)
#  endif
# endif
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].getreceivepos_nowait(targetparticle);
  }
#endif

}

template<class PA>
void
MPIReceiverForLong::exchangeParticleArraysubset(PA& targetparticle,
                                                std::vector<ParticleRange>& target_range,
                                                std::vector<int>& recvsetid,
                                                std::map<int,int>& recvsetid_to_index,
                                                std::vector<TypeRange>& targettyperange)
{
  //  transferParticle();
  //  recvsetid_to_index.clear();
  //  int pre_num = targetparticle.size();
  //  std::cout << " getReceiveParticle   "  << unit_id << std::endl;
  getReceiveParticle(targetparticle,target_range,targettyperange,recvsetid,recvsetid_to_index);
  //  std::cout << " receive " << targetparticle.size()-pre_num << " particle for long" << std::endl;
  //  std::cout << " getReceiveParticle done  "  << unit_id << std::endl;
}

template<class PA>
void
MPIReceiverForLong::exchangeParticleArraysubset_onlyPosition(PA& targetparticle,
                                                             const std::vector<ParticleRange>& target_range,
                                                             const std::vector<int>& recvsetid,
                                                             const std::map<int,int>& recvsetid_to_index,
                                                             const std::vector<TypeRange>& targettyperange)
{
  //  transferParticle();
  getReceiveParticle_onlyPosition(targetparticle);
}

#ifdef OLDPARTICLE
template
void
MPIReceiverForLong::exchangeParticleArraysubset(ParticleArray& targetparticle,
                                                std::vector<ParticleRange>& target_range,
                                                std::vector<int>& recvsetid,
                                                std::map<int,int>& recvsetid_to_index,
                                                std::vector<TypeRange>& targettyperange);
template
void
MPIReceiverForLong::exchangeParticleArraysubset_onlyPosition(ParticleArray& targetparticle,
                                                             const std::vector<ParticleRange>& target_range,
                                                             const std::vector<int>& recvsetid,
                                                             const std::map<int,int>& recvsetid_to_index,
                                                             const std::vector<TypeRange>& targettyperange);
#else
template
void
MPIReceiverForLong::exchangeParticleArraysubset(GhostParticleArray& targetparticle,
                                                std::vector<ParticleRange>& target_range,
                                                std::vector<int>& recvsetid,
                                                std::map<int,int>& recvsetid_to_index,
                                                std::vector<TypeRange>& targettyperange);
template
void
MPIReceiverForLong::exchangeParticleArraysubset_onlyPosition(GhostParticleArray& targetparticle,
                                                             const std::vector<ParticleRange>& target_range,
                                                             const std::vector<int>& recvsetid,
                                                             const std::map<int,int>& recvsetid_to_index,
                                                             const std::vector<TypeRange>& targettyperange);
#endif

void
MPIReceiverForLong::setSendForce(ForceArray& sendforce,
                                 std::vector<ParticleRange>& send_range) 
{
  for(int t=0;t<rpsf_number_of_target;t++){
    if(rpsf[t].setsend(sendforce,send_range[t+target_range_offset])==false){
      printf("MPIReceiverForLong::setSendForce fail\n");
    }
  }
}

void 
MPIReceiverForLong::transferForce()
{
  for(int t=0;t<rpsf_number_of_target;t++){
    rpsf[t].send();
  }
}

void 
MPIReceiverForLong::wait_complete_Force()
{
  for(int t=0;t<rpsf_number_of_target;t++){
    rpsf[t].wait_complete_send();
  }
}

void
MPIReceiverForLong::exchangeForceArraysubset(ForceArray& sendforce,
                                             std::vector<ParticleRange>& send_range)
{
  setSendForce(sendforce,send_range);
  transferForce();
  wait_complete_Force();
}

MPISenderForLong::MPISenderForLong(const std::vector<int>& t_rank,
                                   std::vector<MPI_Comm> ls_comms,
                                   std::vector<int> sprf_tid,
                                   const std::vector< std::vector<int> >& sprf_set,
                                   std::vector<int> sid,
                                   int ps, const int maxid, 
                                   const std::vector<int>& ids
                                   )
  : unit_id_list(ids), psize(ps), target_rank(t_rank),
    mpi_comm_ls_list(ls_comms),
    sprf(), full_target_id(sprf_tid), setid(sid), full_setid(sprf_set), 
    allsendsetid(),
    bond()
{
  //  std::cout << " construct MPISenderForLong  " << std::endl;
  construct_send();
  make_id_to_index(setid,setid_to_index); 
}

MPISenderForLong::MPISenderForLong(const MPISenderForLong& comm)
  : unit_id_list(comm.unit_id_list), psize(comm.psize), 
    target_rank(comm.target_rank),
    mpi_comm_ls_list(comm.mpi_comm_ls_list),
    sprf(), 
    full_target_id(comm.full_target_id),
    sprf_target_id(),
    setid(comm.setid), setid_to_index(), 
    full_setid(comm.full_setid),
    sprf_setid(),
    allsendsetid(comm.allsendsetid),
    bond()
{
  //  std::cout << " copy MPISenderForLong  " << std::endl;
  construct_send();
  //  if(sprf.size()>0){
  //    std::cout << " sprf[0].send_requestp  " << sprf[0].send_requestp  << std::endl;
  //  }
 make_id_to_index(setid,setid_to_index); 
}

MPISenderForLong::~MPISenderForLong()
{
  sprf.clear();
}

void
MPISenderForLong::construct_send()
{
  sprf.clear();
  sprf_setid.clear();
  //size_t num_set = setid.size();
  int si = 0;
  for(std::vector<int>::size_type i=0;i<full_target_id.size();i++){
    if(mpi_comm_ls_list[i]!=MPI_COMM_NULL){
      if(unit_id_list[i]!=0){  // this node is not Long node of this long comm
        sprf.push_back(MPISendParticleReceiveForce(unit_id_list[i],
                                                   psize*full_setid[i].size(),
                                                   mpi_comm_ls_list[i]));
        sprf[si].make_buffer(full_setid[i].size());
        sprf[si].target_id = full_target_id[i];
        sprf[si].target_rank = target_rank[i];
	//	sprf[si].send_requestp = new MPI_Request();
	//	std::cout << " sprf.send_requestp " << sprf[si].unit_identifier  << " " << sprf[si].send_requestp << std::endl;
        sprf_setid.push_back(full_setid[i]);
        sprf_target_id.push_back(full_target_id[i]);
        si++;
      }
    }
  }
  sprf_number_of_target = sprf_target_id.size();
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].send_requestp = new MPI_Request();
  }
}

template<class PA>
void
MPISenderForLong::setSendParticlesubset(const PA& myparticle,
                                        const std::vector<TypeRange>& typerangearray)
{
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].setsendsubset(myparticle,typerangearray,bond,setid,sprf_setid[t]);
  }
}

template<class PA>
void
MPISenderForLong::setSendParticlesubset_onlyPosition(const PA& myparticle,
                                                     const std::vector<TypeRange>& typerangearray)
{
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].setsendsubsetpos(myparticle,typerangearray,setid,sprf_setid[t]);
  }
}

void
MPISenderForLong::transferParticle()
{
  //  std::cout << " send " << sprf_number_of_target << std::endl;
  for(int t=0;t<sprf_number_of_target;t++){
    //    std::cout << " sprf["<< t << "].send() " << sprf[t].unit_identifier << " to " << sprf[t].target_rank << ":" << sprf[t].target_id << " in "  << sprf[t].mpi_comm_short << " " << t << "/" << sprf_number_of_target << " " << sprf[t].send_requestp << std::endl;
    sprf[t].send();
  }
}

void
MPISenderForLong::wait_complete_Particle()
{
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].wait_complete_send();
  }
}

template<class PA>
void
MPISenderForLong::exchangeParticleArraysubset_top_half(const PA& myparticle,
						       const std::vector<TypeRange>& typerangearray)
{
  //  std::cout << " setSendParticlesubset  "  << std::endl;
  setSendParticlesubset(myparticle,typerangearray);
  //  std::cout << " transferParticle  "  << std::endl;
  transferParticle();
}
template<class PA>
void
MPISenderForLong::exchangeParticleArraysubset_bottom_half(const PA& myparticle,
							  const std::vector<TypeRange>& typerangearray)
{
  //  std::cout << " wait_complete_Particle  " << std::endl;
  wait_complete_Particle();
  //  std::cout << " wait_complete_Particle done  "  << std::endl;
}
template<class PA>
void
MPISenderForLong::exchangeParticleArraysubset(const PA& myparticle,
                                              const std::vector<TypeRange>& typerangearray)
{
  //  std::cout << " setSendParticlesubset  "  << std::endl;
  setSendParticlesubset(myparticle,typerangearray);
  //  std::cout << " transferParticle  "  << std::endl;
  transferParticle();
  //  std::cout << " wait_complete_Particle  " << std::endl;
  wait_complete_Particle();
  //  std::cout << " wait_complete_Particle done  "  << std::endl;
}

template<class PA>
void
MPISenderForLong::exchangeParticleArraysubset_onlyPosition_top_half(const PA& myparticle,
								    const std::vector<TypeRange>& typerangearray)
{
  setSendParticlesubset_onlyPosition(myparticle,typerangearray);
  transferParticle();
}
template<class PA>
void
MPISenderForLong::exchangeParticleArraysubset_onlyPosition_bottom_half(const PA& myparticle,
								       const std::vector<TypeRange>& typerangearray)
{
  wait_complete_Particle();
}
template<class PA>
void
MPISenderForLong::exchangeParticleArraysubset_onlyPosition(const PA& myparticle,
                                                           const std::vector<TypeRange>& typerangearray)
{
  setSendParticlesubset_onlyPosition(myparticle,typerangearray);
  transferParticle();
  wait_complete_Particle();
}

#ifdef OLDPARTICLE
template
void
MPISenderForLong::exchangeParticleArraysubset(const ParticleArray& myparticle,
                                              const std::vector<TypeRange>& typerangearray);

template
void
MPISenderForLong::exchangeParticleArraysubset_onlyPosition(const ParticleArray& myparticle,
                                                           const std::vector<TypeRange>& typerangearray);
#else
template
void
MPISenderForLong::exchangeParticleArraysubset_top_half(const CombinedParticleArray& myparticle,
						       const std::vector<TypeRange>& typerangearray);
template
void
MPISenderForLong::exchangeParticleArraysubset_bottom_half(const CombinedParticleArray& myparticle,
							  const std::vector<TypeRange>& typerangearray);
template
void
MPISenderForLong::exchangeParticleArraysubset(const CombinedParticleArray& myparticle,
                                              const std::vector<TypeRange>& typerangearray);

template
void
MPISenderForLong::exchangeParticleArraysubset_onlyPosition_top_half(const CombinedParticleArray& myparticle,
								    const std::vector<TypeRange>& typerangearray);
template
void
MPISenderForLong::exchangeParticleArraysubset_onlyPosition_bottom_half(const CombinedParticleArray& myparticle,
								       const std::vector<TypeRange>& typerangearray);
template
void
MPISenderForLong::exchangeParticleArraysubset_onlyPosition(const CombinedParticleArray& myparticle,
                                                           const std::vector<TypeRange>& typerangearray);
#endif

void 
MPISenderForLong::transferForce()
{
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].prepare_receive();
  }
}

void 
MPISenderForLong::getReceiveForcesubset(ForceArray& recvforce)
{
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].getreceivesubset(recvforce);
  }
}

void 
MPISenderForLong::exchangeForceArraysubset(ForceArray& recvforce)
{
  //  transferForce();
  getReceiveForcesubset(recvforce);
#if 0
  {
    Force sum(0.0,0.0,0.0);
    for(int i=0;i<recvforce.size();i++){
      sum += recvforce[i];
    }
    std::cout << " sum of received force " << sum << std::endl;
  }
#endif
}
