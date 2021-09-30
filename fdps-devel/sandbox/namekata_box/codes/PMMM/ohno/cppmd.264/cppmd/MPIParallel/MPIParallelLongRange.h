//! MPIParallel operator for Long-Range Interaction
#ifndef MPIPARALLELLONGRANGE_H
#define MPIPARALLELLONGRANGE_H

#include <mpi.h>
#include "CellIndex.h"
#include "LongRangeParameter.h"
#include "MPIParallel.h"

/*! MPI Receiver / Sender of Particles for Long Range calculation
  Long Range calculation node belong to one Long Range MPI_Comm as rank 0.
  Particle Holder (Short/Integrate) node belong to several Long Range MPI_Comm.
 */



class MPIReceiverForLong {
 public:
  int unit_id;                        //!< unit ID of this node
  size_t psize;                       //!< numver of particle in subset
  std::map<int,int> unitid_to_rank;   //!< reverse map unit ID to rank

  MPI_Comm mpi_comm_long_short;       //!< MPI_Comm for long-short communication
  std::vector<MPIReceiveParticleSendForce> rpsf; //!< Recv j-Particle and return reaction force
  std::vector<int> rpsf_target_id;               //!< IDs of target nodes
  int rpsf_number_of_target;                     //!< number of recv target nodes, equal to rpsf_target_id.size().
  std::vector< std::vector<int> > targetsetid;   //!< set IDs that each target node send to this node.
  std::vector<int> alltargetsetid;               //!< set IDs that some target nodes send to this node.
  int number_of_target_set;                      //!< number of particlesets from some nodes. @note for future use, used only in testcode.
  int target_range_offset;                       //!< offset of target_range, larger than 0 when this node is short-node and receive ghost for short.

  std::vector<CovalentBondInfo::BondList> bond;  //!< dummy for compatibility


  MPIReceiverForLong(){}

  //! constructor 
  MPIReceiverForLong(const std::map<int,int>& uidtor, 
                     MPI_Comm long_short_comm,
                     std::vector<int> rpsf_tid,
                     std::vector< std::vector<int> > tsid,
                     int ps, const int maxid, const int id
                     );

  //! copy constructor
  MPIReceiverForLong(const MPIReceiverForLong& comm);

  //  MPICommunicator& operator=(const MPICommunicator& comm);

  ~MPIReceiverForLong();

  void construct_recv();

  void transferParticle();

  template<class PA>
  void getReceiveParticle(PA& targetparticle,
                          std::vector<ParticleRange>& target_range,
                          std::vector<TypeRange>& targettyperange,
                          std::vector<int>& recvsetid,
                          std::map<int,int>& recvsetid_to_index);

  template<class PA>
  void getReceiveParticle_onlyPosition(PA& targetparticle);

  template<class PA>
  void exchangeParticleArraysubset(PA& targetparticle,
                                   std::vector<ParticleRange>& target_range,
                                   std::vector<int>& recvsetid,
                                   std::map<int,int>& recvsetid_to_index,
                                   std::vector<TypeRange>& targettyperange);

  template<class PA>
  void exchangeParticleArraysubset_onlyPosition(PA& targetparticle,
                                                const std::vector<ParticleRange>& target_range,
                                                const std::vector<int>& recvsetid,
                                                const std::map<int,int>& recvsetid_to_index,
                                                const std::vector<TypeRange>& targettyperange);

  void setSendForce(ForceArray& sendforce,
                    std::vector<ParticleRange>& send_range);

  void transferForce();

  void wait_complete_Force();

  void exchangeForceArraysubset(ForceArray& sendforce,
                                std::vector<ParticleRange>& send_range);

};

class MPISenderForLong {
 public:
  std::vector<int> unit_id_list;                 //!< unit ID of this node
  size_t psize;                       //!< numver of particle in subset
  std::vector<int> target_rank;                  //!< MPI_Rank of target

  std::vector<MPI_Comm> mpi_comm_ls_list;        //!< MPI_Comm for long-short communication
  std::vector<MPISendParticleReceiveForce> sprf; //!< Send Own particle and Recieve reaction force
  std::vector<int> full_target_id;               //!< IDs of target ndoes 
  std::vector<int> sprf_target_id;               //!< IDs of target ndoes only active
  int sprf_number_of_target;                     //!< number of send target nodes, equal to srpf_target_id.size()
  std::vector<int> setid;                        //!< set IDs of this node
  std::map<int,int> setid_to_index;              //!< reverse map setID to index
  std::vector< std::vector<int> > full_setid;    //!< set IDs which each target node required.
  std::vector< std::vector<int> > sprf_setid;    //!< set IDs which each active target node required.
  std::vector<int> allsendsetid;                 //!< set IDs which required some target node.  
  int number_of_send_set;                        //!< number of particlesets which required some target node. @note for future use.

  std::vector<CovalentBondInfo::BondList> bond;  //!< dummy for compatibility

  MPISenderForLong(){}

  //! constructor 
  MPISenderForLong(const std::vector<int>& t_rank, 
                   std::vector<MPI_Comm> ls_comms,
                   std::vector<int> sprf_tid,
                   const std::vector< std::vector<int> >& srpf_set,
                   std::vector<int> sid,
                   int ps, const int maxid, const std::vector<int>& id
                   );

  //! copy constructor
  MPISenderForLong(const MPISenderForLong& comm);

  //  MPICommunicator& operator=(const MPICommunicator& comm);

  ~MPISenderForLong();

  void construct_send();

  template<class PA>
  void setSendParticlesubset(const PA& myparticle,
                             const std::vector<TypeRange>& typerangearray);

  template<class PA>
  void setSendParticlesubset_onlyPosition(const PA& myparticle,
                                          const std::vector<TypeRange>& typerangearray);

  void transferParticle();

  void wait_complete_Particle();

  template<class PA>
  void exchangeParticleArraysubset_top_half(const PA& myparticle,
					    const std::vector<TypeRange>& typerangearray);
  template<class PA>
  void exchangeParticleArraysubset_bottom_half(const PA& myparticle,
					       const std::vector<TypeRange>& typerangearray);
  template<class PA>
  void exchangeParticleArraysubset(const PA& myparticle,
                                   const std::vector<TypeRange>& typerangearray);

  template<class PA>
  void exchangeParticleArraysubset_onlyPosition_top_half(const PA& myparticle,
							 const std::vector<TypeRange>& typerangearray);
  template<class PA>
  void exchangeParticleArraysubset_onlyPosition_bottom_half(const PA& myparticle,
							    const std::vector<TypeRange>& typerangearray);
  template<class PA>
  void exchangeParticleArraysubset_onlyPosition(const PA& myparticle,
                                                const std::vector<TypeRange>& typerangearray);

  void transferForce();

  void getReceiveForcesubset(ForceArray& recvforce);

  void exchangeForceArraysubset(ForceArray& recvforce);
};

#endif

