/*! @file
  @brief Communicator by MPI

  Node Type        before force-calculation  after force-calculation
  have Particles   Send Particles            Receive Forces 
                   1 set to N nodes          N sets from N nodes, reduce to 1 set

  calculate force  Receive Particles         Send Forces 
                   N sets from N nodes       N sets to N nodes

*/


#ifndef MPIPARALLEL_H
#define MPIPARALLEL_H
#include <mpi.h>

#include <vector>
#include <map>

#include "ParticleInfo.h"
#include "CovalentBondInfo.h"
#include "SetPairList.h"
#include "CellIndex.h"
#include "CubicCell.h"
#include "Geometry.h"

#include "Timer.h"


const int buffersize = 8192;

/*! 
  Required data of Particles
  Type of Force : Required Data 
  ShortRange    : position, charge, atomtype, typerange
  LongRange     : position, charge
  CovalentBond  : position, atomid, bondlist  : at least bonded to outside

  Original Particle data : AOS format
  send buffer : SOA format
  int         Sender Unit_Identifier \
  int         Number of Position     |
  int         Number of Charge       | 
  int         Number of Atomtype     | ParticleBufferHeader
  int         Number of AtomID       |
  int         Number of ParticleSet  |
  int         Number of Bondlist     /
  Position[]  Array of Position
  double[]    Array of Charge
  Atomtype[]  Array of Atomtype
  AtomID[]    Array of AtomID
  int[]       Array of ParticleSetIdentifier
  TypeRange[] Array of Typerange
  CovalentBondInfo::BondList[]  Array of Bondlist

  6*int + NP*3*double + NC*double + NA*int + NT*Typerange + NB*Bondlist

  Limmitation of this version
  structure of Bondlist is not fixed
  Number of Bondlist must be 0

*/
//! Header of Send/Receive Positions Buffer
struct PositionBufferHeader {
  int UnitID;                   //!< unit ID of sender
  int NumberOfPosition;
};

//! Header of Send/Receive Particles Buffer
/*!
  At time step when particle move, NumberOfPosition == NumberOfCharge == NumberOfAtomType == NumberOfAtomID .
  Charge, AtomType and AtomID are not sent when no particle move in/out at sender. In this case, NumberOfCharge, NumberOfAtomType and NumberOfAtomID are 0 .
*/
struct ParticleBufferHeader {
  int UnitID;                   //!< unit ID of sender
  int NumberOfPosition;
  int NumberOfCharge;
  int NumberOfAtomType;
  int NumberOfAtomID;
  int NumberOfParticleSet;
  int NumberOfBondlist;
};

//! Header of Send/Receive Force Buffer
struct ForceBufferHeader {
  int UnitID;                    //!< unit ID of sender
  int NumberOfForce;
};

//! Header of Send/Receive Move Particle Buffer
struct MoveBufferHeader {
  int UnitID;                     //!< unit ID of sender
  int NumberOfParticle;
  int NumberOfBondlist;
  int SizeOfBondlist;
};
  
enum CommPattern {
  Direct = 0,                     //!< direct to/from all target
  NearestXYZ                      //!< nearest, multistage
};

// MPITagType
#define  MPITagPosition 0
#define  MPITagParticle 1
#define  MPITagForce    2
#define  MPITagMove     3

// suport functions calculate size of buffer
//! return size of Positoins Buffer in Byte
/*!
  @param[in] number_of_position number of position
 */
size_t calc_position_buffer_size(int number_of_position);

//! return size of Particles Buffer in Byte
/*!
  @param[in] number_of_position number of position
  @param[in] number_of_charge number of charge
  @param[in] number_of_atomtype number of atomtype
  @param[in] number_of_atomid number of atomid
  @param[in] number_of_particleset number of particleset
  @param[in] number_of_bondlist number of bond
  @return size of particle buffer
 */
size_t calc_particle_buffer_size(int number_of_position,
                                 int number_of_charge,
                                 int number_of_atomtype,
                                 int number_of_atomid,
                                 int number_of_particleset,
                                 int number_of_bondlist );

//! return size of Particles Buffer in Byte
/*!
  @param[in] header header of particle buffer
  @return size of particle buffer
 */
size_t calc_particle_buffer_size(ParticleBufferHeader& header);

//! return size of Force Buffer in Byte
/*!
  @param[in] number_of_force number of force
  @return size of force buffer
 */
size_t calc_force_buffer_size(int number_of_force);

//! return size of Move Particle Buffer in Byte
/*!
  @param[in] number_of_move number of move particle
  @return size of move buffer
  @note this function for no-bond simulation
 */
size_t calc_move_buffer_size(int number_of_move);
//! return size of Move Particle Buffer in Byte
/*!
  @param[in] number_of_move number of move particle
  @param[in] bondlistarray array of bondlist
  @return size of move buffer
 */
size_t calc_move_buffer_size(int number_of_move,
                             std::vector<CovalentBondInfo::BondList>& bondlistarray);

//! return size of Move Particle Buffer in Byte
/*!
  @param[in] header header of move buffer
  @return size of move buffer
  @note this function for no-bond simulation
 */
size_t calc_move_buffer_size(MoveBufferHeader& header);
//! return size of Move Particle Buffer in Byte
/*!
  @param[in] header header of move buffer
  @param[in] bondlistarray array of bondlist
  @return size of move buffer
 */
size_t calc_move_buffer_size(MoveBufferHeader& header,
                             std::vector<CovalentBondInfo::BondList>& bondlistarray);


/*
  communication sequence
  Sender own particle               Receiver j-particle
  MPISendParticleReceiveForce       MPIReceiveParticleSendForce
  setsendsubset()
                                     prepare_receive()
  send()
                                     getreceive()
                                     PostProcess::postprocess_receive()
  wait_complete_send()

                    if use force symmetry
                                     setsend()
  prepare_receive()
                                     send()
  getreceivesubset()
                                     wait_complete_send()
    

*/

//! Send Own Particle and Receive Force to/from ONE target node
class MPISendParticleReceiveForce {
 public:
  size_t size;                         //!< default number of particle
  int target_rank;                     //!< rank of target
  int target_id;                       //!< unit ID of target
  int unit_identifier;                 //!< unit ID of this node

  MPI_Comm mpi_comm_short;             //!< MPI_Comm for short

  // TODO make array/vector for NearestXYZ(multistage)
  size_t number_of_stage;

  std::vector<size_t> number_of_particle;    //!< actual number of particle
  
  std::vector<size_t> number_of_set;   //!< number of particle-set

  /*!
    In XYZ mode, communication between next two or more is not directly, by multi stage nearest cummunication.
    MPICommunicator control this multi stage communication.
    This member is used to identify such stages.
   */
  size_t stage;

  /*!
    Pointers of Paricle data in Buffer
    Send buffer is defined as byte and actual data is mapped
    because it is defficult to make MPI data type with array which has variable size.
  */
  PositionBufferHeader *sendposbufferheader;  //!< Pointer of Header in Buffer
  ParticleBufferHeader *sendbufferheader;  //!< Pointer of Header in Buffer
  Position *sendposition;                  //!< Pointer of Postions in Buffer
  double *sendcharge;                      //!< Pointer of Charges in Buffer
  Atomtype *sendatomtype;                  //!< Pointer of AtomTypes in Buffer
  AtomID *sendatomid;                      //!< Pointer of AtomIDs in Buffer
  int *sendsetid;                          //!< Pointer of set IDs in Buffer
  TypeRange *sendtyperange;                //!< Pointer of TypeRanges in Buffer
  CovalentBondInfo::BondList *sendbond;    //!< Pointer of BondLists in Buufer

  /*!
    send_buffer is allocated by make_buffer for default number of particle.
    make_buffer is called at initialization.
    Because send_buffer is not reallcated at actual sending,
    size of actual data size is checked. 
  */
  char *send_buffer;                   //!< Pointer of Send Paricle MPI Buffer
  size_t send_buffer_size;             //!< Size of reserved Send Paricle Buffer
  size_t send_request_size;            //!< Size of Send Buffer calculate from requested number of Particle
  MPI_Request *send_requestp;              //!< used MPI_Isend and MPI_Wait
  std::vector<TypeRange> sent_typerange; //!< store TypeRange of sent Particle and used to restore receive force at same Range
  std::vector<int> set_index;            //!< index of send Particle Set, temporarily stored

  ForceBufferHeader *receivebufferheader; //!< Pointer of Header in Buffer
  Force *ReceiveForceArray;               //!< Pointer of Forces in Buffer
  int *ReceiveForceIndex;               //!< Pointer of Index of Forces in Buffer
  std::vector<int> forceindexarray;   //!< index in retrun ForceArray 
  char *receive_buffer;               //!< Pointer of Receive Force MPI Buffer
  size_t receive_buffer_size;         //!< Size of reserved Receive Force Buffer
  size_t receive_expected_size;       //!< Size of Force Buffer estimated from numver of send Particles
  MPI_Request receive_requset;        //!< used MPI_Irecv and MPI_Wait

  //! construct, default size is hard coded
  MPISendParticleReceiveForce();

  //! construct, with unit ID, default size is hard coded
  explicit MPISendParticleReceiveForce(int unitid, MPI_Comm short_comm);

  //! construct, with unit ID and specified size
  MPISendParticleReceiveForce(int unitid, size_t sz, MPI_Comm short_comm);

  //! construct, with unit ID, specified size and number of stage
  MPISendParticleReceiveForce(int unitid, size_t sz, MPI_Comm short_comm,
                              int num_stage);

  //! copy constructor
  MPISendParticleReceiveForce(const MPISendParticleReceiveForce& sprf);

  //! destructor
  ~MPISendParticleReceiveForce();

  //! assignment operator
  MPISendParticleReceiveForce& operator=(const MPISendParticleReceiveForce& sprf);

  //! set position data pointers in send MPI buffer
  /*! 
    @param[in] number_of_position
   */ 
  bool set_send_posbuffer_pointers(int number_of_position);

  //! set particle data pointers in send MPI buffer
  /*! 
    @param[in] number_of_position
    @param[in] number_of_charge
    @param[in] number_of_atomtype
    @param[in] number_of_atomid
    @param[in] number_of_particleset
    @param[in] number_of_bondlist
    @return set successfully or not
    @retval true required buffer size <= preallocated buffer size
    @retval false required buffer size > preallocated buffer size
    @todo 2nd case, realloc send_buffer and not return false
  */
  bool set_send_buffer_pointers(int number_of_position,
                                int number_of_charge,
                                int number_of_atomtype,
                                int number_of_atomid,
                                int number_of_particleset,
                                int number_of_bondlist );

  //! set force data pointers in receive MPI buffer
  /*! 
    @param[in] number_of_force number of force
    @return set successfully or not
    @retval true expected buffer size <= preallocated buffer size
    @retval false expected buffer size > preallocated buffer size
    @todo 2nd case, realloc send_buffer and not return false
  */
  bool set_receive_buffer_pointers(int number_of_force);
  bool set_receive_buffer_pointers_with_index(int number_of_force);
  //! reserve send and receive MPI buffer and set number of particle and particleset
  /*!
    @param[in] num_set number of particleset
    Set number_of_set[stage], number_of_particle[stage], send_buffer_size and receive_buffer_size.
    Construct send_buffer and receive_buffer.
    Set send/receive_buffer_pointers
    Initialize set_typerange, set_index.
  */
  void make_buffer(const size_t& num_set);

  //! set send Positions to send MPI buffer
  /*!
    @param[in] request particlesets contaion positions to be send
    @return set successfully or not. 
    @retval true success
    @retval false fail because illeagal size. See set_send_buffer_pointers.
    @note It sotore information of sent particlesets to local member. Use getreceivesubset to receive reaction force acording with such sets information.
    @todo use fast copy.
    Paritcles have no gap between subset
  */
  bool setsendpos(ParticleArray& request);

  //! set send Particles to send MPI buffer
  /*!
    @param[in] request particlesets to be send
    @param[in] rq_tr typerange of requested particlesets
    @param[in] rq_bond bond information to be send
    @param[in] setid set-ID of requeseted particlesets
    @return set successfully or not. 
    @retval true success
    @retval false fail because illeagal size. See set_send_buffer_pointers.
    @note It sotore information of sent particlesets to local member. Use getreceivesubset to receive reaction force acording with such sets information.
    @todo use fast copy.
    Paritcles have no gap between subset
  */
  bool setsend(ParticleArray& request, std::vector<TypeRange>& rq_tr,
               std::vector<CovalentBondInfo::BondList>& rq_bond, 
               std::vector<int>& setid);

  //! set send Positions of Particles Sets to send MPI buffer
  /*!
    @param[in] request particlesets to be send
    @param[in] rq_tr typerange of requested particlesets
    @param[in] setid set-ID of requeseted particlesets
    @param[in] rq_set set-ID of particlesets to be send
    @return set successfully or not. 
    @retval true success
    @retval false fail because illeagal size. See set_send_buffer_pointers.
    @note It sotore information of sent particlesets to local member. Use getreceivesubset to receive reaction force acording with such sets information.
    @todo use fast copy.
    Paritcles have gap between subset ( reserve max size of each subset )
  */
  template<class PA>
  bool setsendsubsetpos(const PA& request, 
                        const std::vector<TypeRange>& rq_tr,
                        const std::vector<int>& setid,
                        const std::vector<int>& rq_set);

  //! set send Particle Sets to send MPI buffer
  /*!
    @param[in] request particlesets to be send
    @param[in] rq_tr typerange of requested particlesets
    @param[in] rq_bond bond information to be send
    @param[in] setid set-ID of requeseted particlesets
    @param[in] rq_set set-ID of particlesets to be send
    @return set successfully or not. 
    @retval true success
    @retval false fail because illeagal size. See set_send_buffer_pointers.
    @note It sotore information of sent particlesets to local member. Use getreceivesubset to receive reaction force acording with such sets information.
    @todo use fast copy.
    Paritcles have gap between subset ( reserve max size of each subset )
  */
  template<class PA>
  bool setsendsubset(const PA& request,
                     const std::vector<TypeRange>& rq_tr,
                     const std::vector<CovalentBondInfo::BondList>& rq_bond, 
                     const std::vector<int>& setid,
                     const std::vector<int>& rq_set);

  //! get forces from MPI recv buffer
  /*!
    @param[inout] forcearray that be added force array received from other node. Its size must be larger than number of sent particles.
    @return receive successfully or not
    @retval true number of received forces <= number of sent particles.
    @retval false Fail MPI_Wait or number of received forces > number of sent particles. No received forces are added to forcearray.
    @note prepare_receive must be called before call this method.
    @note It caused no error that number of received forces is smaller than number of sent particles.
    @note It store forces cotinuously without sets information.
  */
  bool getreceive(ForceArray& forcearray);

  //! get forces from MPI recv buffer with offset
  /*!
    @param[inout] forcearray that be added force array received from other node. Its size must be larger than offset + number of sent particles .
    @param[in] offset offset where received force store from
    @return receive successfully or not
    @retval true number of received forces <= number of sent particles.
    @retval false number of received forces > number of sent particles. No received forces are added to forcearray. 
    @note prepare_receive must be called before call this method.
    @note It caused no error that number of received forces is smaller than number of sent particles.
    @note It store forces cotinuously without sets information.
  */
  bool getreceive(ForceArray& forcearray, int offset);

  //! get forces from MPI recv buffer with offset
  /*!
    @param[inout] forcearray that be added force array received from other node. It size must be larger than or equal to type ranges of sent particlesets.
    @return receive successfully or not
    @retval true number of received forces <= number of sent particles.
    @retval false number of received forces > number of sent particles. No received forces are added to forcearray. 
    @note prepare_receive must be called before call this method.
    @note It store forces acording with sets information which stored by setsendsubset.
  */
  bool getreceivesubset(ForceArray& forcearray);

  bool getreceive_with_index(ForceArray& forcearray);
  bool getreceive_indexed_nowait(ForceArray& forcearray);
  bool getreceive_indexed(ForceArray& forcearray);

  //! request MPI asyncronous receive for Force
  /*!
    @note It must be called before calling getreceive*.
   */
  void prepare_receive();

  //! request MPI send for Particle
  /*!
    If defined PARTICLE_ISEND, use asyncronous send, and must be check by wait_complete_send() .
    If nodef PARTICLE_ISEND, use syncronous send, and MPIReceiveParticleSendForce::prepare_receive() must be called before this method.
    Recommend define PARTICLE_ISEND for multiple target node because receiver may be ready not in same order of send request.
    But no warranty multiple Isends are simultaneously waiting or start transmission not in requested order.
  */
  void send();

  //! wait/check finish send
  /*!
    If defined PARTICLE_ISEND, MUST be called.
    If nodef PARTICLE_ISEND, do nothing.
  */
  void wait_complete_send();
};

//! Receive j-Particle and Send Force from/to ONE target node
class MPIReceiveParticleSendForce {
 public:
  size_t size;                         //!< default number of particle
  int target_rank;                     //!< rank of target
  int target_id;                       //!< unit ID of target
  int unit_identifier;                 //!< unit ID of this node

  MPI_Comm mpi_comm_short;             //!< MPI_Comm for short

  // TODO make array/vector for NearestXYZ(multistage)
  size_t number_of_stage;

  std::vector<size_t> number_of_particle;           //!< actual number of particle
  
  std::vector<size_t> number_of_set;

  /*!
    In XYZ mode, communication between next two or more is not directly, by multi stage nearest cummunication.
    MPICommunicator control this multi stage communication.
    This member is used to identify such stages.
   */
  size_t stage;

  /*!
    Pointers of Paricle data in Buffer
    buffer is defined as byte and actual data is mapped
    because it is defficult to make MPI data type with array which has variable size.
  */
  PositionBufferHeader *recvposbufferheader;  //!< Pointer of Header in Buffer
  Position *recvposposition;                  //!< Pointer of Postions in Buffer
  ParticleBufferHeader *recvbufferheader;  //!< Pointer of Header in Buffer
  Position *recvposition;                  //!< Pointer of Postions in Buffer
  double *recvcharge;                      //!< Pointer of Charges in Buffer
  Atomtype *recvatomtype;                  //!< Pointer of AtomTypes in Buffer
  AtomID *recvatomid;                      //!< Pointer of AtomIDs in Buffer
  int *recvsetid;                          //!< Pointer of set IDs in Buffer
  TypeRange *recvtyperange;                //!< Pointer of TypeRanges in Buffer
  CovalentBondInfo::BondList *recvbond;    //!< Pointer of BondLists in Buufer

  /*!
    recv_buffer is allocated by make_buffer for default number of particle.
    make_buffer is called at initialization.
    Because recv_buffer is not reallcated at actual receiving,
    size of actual data size is checked. 
  */
  char *recv_buffer;                  //!< Pointer of Receive Paricle MPI Buffer
  size_t recv_buffer_size;            //!< Size of reserved Recv Paricle Buffer
  size_t recv_requested_size;         //!< Size of Recv Buffer calculate from requested number of Particle
  //!< offset of received particlearray
  /*!
    getreceive set offset store particlearray to this variable
    getreceivepos use this variable as offset re particlearray
   */
  std::vector<size_t> recv_pos_offset;             
  MPI_Request *recv_requestp;           // used MPI_Irecv and MPI_Wait

  ForceBufferHeader *sendbufferheader; //!< Pointer of Header in Buffer
  Force *SendForceArray;               //!< Pointer of Forces in Buffer
  int *SendForceIndex;               //!!< Pointer of Index of Forces
  size_t indexbegin;
  size_t indexend;
  size_t indexlast;

  char *send_buffer;                 //!< Pointer of Receive Force MPI Buffer
  size_t send_buffer_size;           //!< Size of reserved Receive Force Buffer
  size_t send_expected_size;         //!< Size of Force Buffer estimated from numver of receive Particles
  MPI_Request send_request;          //!< used MPI_Isend and MPI_Wait

  //! construct, default size is hard corded
  MPIReceiveParticleSendForce();

  //! construct, with unit ID, default size is hard coded
  explicit MPIReceiveParticleSendForce(int unitid, MPI_Comm short_comm);

  //! construct, with unit ID and specified size
  MPIReceiveParticleSendForce(int unitid, size_t sz, MPI_Comm short_comm);

  //! construct, with unit ID, specified size and number of stage
  MPIReceiveParticleSendForce(int unitid, size_t sz, MPI_Comm short_comm, 
                              int num_stage);

  //! copy constructor
  MPIReceiveParticleSendForce(const MPIReceiveParticleSendForce &rpsf);

  //! assignment operator
  MPIReceiveParticleSendForce& operator=(const MPIReceiveParticleSendForce &rpsf);

  //! destructor
  ~MPIReceiveParticleSendForce();

  //! set positinos of particle data pointers in recv MPI buffer
  /*!
    @return set successfully or not
    @retval true size of received data <= size of preallocated buffer
    @retval false size of received data > size of preallocated buffer
    @todo 2nd case, realloc recv_buffer and not return false
   */
  bool set_receive_posbuffer_pointers();

  //! set particle data pointers in recv MPI buffer
  /*!
    @return set successfully or not
    @retval true size of received data <= size of preallocated buffer
    @retval false size of received data > size of preallocated buffer
    @todo 2nd case, realloc recv_buffer and not return false
   */
  bool set_receive_buffer_pointers();

  //! set force data pointers in send MPI buffer
  /*!
    @param[in] number_of_force number of force
    @return set successfully or not
    @retval true expected buffer size <= preallocated buffer size
    @retval false expected buffer size > preallocated buffer size
    @todo 2nd case, realloc send_buffer and not return false
   */
  bool set_send_buffer_pointers(int number_of_force);
  bool set_send_buffer_pointers_with_index(int number_of_force);

  //! reserve receive and send MPI buffer
  /*!
    @param[in] num_set number of particleset
    Set number_of_set[stage], number_of_particle[stage], send_buffer_size and receive_buffer_size.
    Construct send_buffer and receive_buffer.
    Set send/receive_buffer_pointers
    Initialize set_typerange, set_index.
   */
  void make_buffer(int num_set);

  bool wait_complete_recv();
  //! get receive Particles from recv MPI buffer
  /*!
    @param[inout] particles which append/overwrite to received particles
    @return set successfully or not. 
    @retval true success
    @retval false Fail MPI_Wait or illeagal size. See set_receive_buffer_pointers.
    @note Received particlesets are stored continuously, without gap between particlesets. 
    @note recv_pos_offset have be set by getreceive
    @todo use fast copy.
    range.begin : index of request where store receive particle from 
  */
  template <class GPA>
  bool getreceivepos_nowait(GPA& request);
  template <class GPA>
  bool getreceivepos(GPA& request);

  //! get receive Particles from recv MPI buffer
  /*!
    @param[inout] particles which append/overwrite to received particles
    @param[in] range Where received particles store from. Received particles are append after ragne.begin. Only range.begin is used.
    @param[inout] rq_tr typerange which append typerange of received particlesets
    @param[inout] rq_bond bond information which append received bond information
    @param[inout] rq_set set-ID which append set-ID of received particlesets
    @param[inout] recvsetid_to_index reverse map set-ID to set-index which append reverse map of received particlesets to.
    @return set successfully or not. 
    @retval true success
    @retval false Fail MPI_Wait or illeagal size. See set_receive_buffer_pointers.
    @note Received particlesets are stored continuously, without gap between particlesets. 
    @todo use fast copy.
    range.begin : index of request where store receive particle from 
  */
  bool getreceive_number(int &pnum, int &setnum, int &bondnum);
  template<class PA>
  bool getreceive_offset(PA& request, 
                         ParticleRange& range,
                         std::vector<TypeRange>& rq_tr,
                         std::vector<CovalentBondInfo::BondList>& rq_bond,
                         std::vector<int>& rq_set,
                         int setoffset, int bondoffset);
  bool getreceive_nowait(ParticleArray& request, ParticleRange& range, 
                         std::vector<TypeRange>& rq_tr,
                         std::vector<CovalentBondInfo::BondList> &rq_bond,
                         std::vector<int>& rq_set,
                         std::map<int,int>& recvsetid_to_index);
  template<class PA>
  bool getreceive(PA& request, ParticleRange& range, 
                  std::vector<TypeRange>& rq_tr,
                  std::vector<CovalentBondInfo::BondList> &rq_bond,
                  std::vector<int>& rq_set,
                  std::map<int,int>& recvsetid_to_index);

  //! set forces to MPI send buffer
  /*!
    @param[in] forcearray to be send
    @param[in] range range of forcearray to  be send
    @return set successfully or not. 
    @retval true success
    @retval false fail because illeagal size. See set_send_buffer_pointers.
    @note It not send range/sets information. These information must be treated at force reciever (that is just particle sender and shuld be hold such information).  
    @todo use fast copy.
    Paritcles have no gap between subset
  */
  bool setsend(const ForceArray& forcearray, const ParticleRange range);

  bool setsend_with_index(const ForceArray& forcearray,
                          const ParticleRange range, 
                          const std::vector<int>& forceindexset);
  bool setsend_indexed(const ForceArray& forcearray,
                       const ParticleRange range, 
                       const std::vector<int>& forceindexset);

  //! request MPI asyncronous receive for Particle
  /*!
    @note It must be called before calling getreceive*.
   */
  void prepare_receive();

  //! request MPI send for Force
  void send();

  //! wait/check finish send
  void wait_complete_send();

};

//! Send Move out Particle to ONE target node
class MPIMoveOut {
 public:
  size_t size;                         //!< default number of particle
  int target_rank;                     //!< rank of target
  int target_id;                       //!< unit ID of target
  int unit_identifier;                 //!< unit ID of this node

  MPI_Comm mpi_comm_short;             //!< MPI_Comm for short

  size_t number_of_particle;           //!< actual number of particle

  MoveBufferHeader *sendbufferheader;  //!< Pointer of Header in Buffer
  Particle *sendparticle;              //!< Pointer of Particle in Buffer
  PotentialModel *sendtype;                 //!< Pointer of ParticleType in Buffer
  int *sendcellid;                 //!< Pointer of Cellid in Buffer
  std::vector<CovalentBondInfo::BondList *> sendbondlistarray;
  std::vector<int *> sendbondpackarray;
  char *send_buffer;             //!< Pointer of Send Paricle MPI Buffer
  size_t send_buffer_size;       //!< Size of reserved Send Paricle Buffer
  size_t send_request_size;      //!< Size of Send Buffer calculate from requested number of Particle
  MPI_Request send_request;      //!< used MPI_Isend and MPI_Wait

  //! construct, default size is hard corded
  MPIMoveOut();

  //! construct, with unit ID, default size is hard coded
  explicit MPIMoveOut(int unitid, MPI_Comm short_comm);

  //! construct, with unit ID and specified size
  MPIMoveOut(int unitid, size_t sz, MPI_Comm short_comm);

  //! copy constructor
  MPIMoveOut(const MPIMoveOut& mo);
  
  //! assignment operator
  MPIMoveOut& operator=(const MPIMoveOut& mo);

  //! destructor
  ~MPIMoveOut();

  //! set particle data pointers in send MPI buffer
  /*! 
    @param[in] number_of_particle
    @return set successfully or not
    @retval true required buffer size <= preallocated buffer size
    @retval false required buffer size > preallocated buffer size
    @todo 2nd case, realloc send_buffer and not return false
  */
  bool set_send_buffer_pointers(int num_particle);
  bool set_send_buffer_pointers(int num_particle,
                                std::vector<CovalentBondInfo::BondList>& bondlistarray);

  //! reserve send MPI buffer
  /*!
    Construct send_buffer.
    Set send_buffer_pointers.
  */
  void make_buffer();

  //! set Move out Particles to MPI send buffer
  /*!
    @param[in] request particlesets to be send
    @param[in] moveouttype potentialmodell of moveout particles
    @param[in] outcellid cell-IDs where paticle move to. 
    @return set successfully or not. 
    @retval true success
    @retval false fail because illeagal size. See set_send_buffer_pointers.
    @note It sotore information of sent particlesets to local member. Use getreceivesubset to receive reaction force acording with such sets information.
    @todo support bond information: covalent-bond or rigid water.
    @todo use fast copy.
    Paritcles have no gap between subset
  */
  bool setmoveoutparticle(ParticleArray& moveout,
                          std::vector<PotentialModel>& moveouttype,
                          std::vector<int>& outcellid);
  bool setmoveoutparticle(ParticleArray& moveout,
                          std::vector<PotentialModel>& moveouttype,
                          std::vector<int>& outcellid,
                          std::vector<CovalentBondInfo::BondList>& outbondlistarray);

  //! request MPI send for Move out Particle
  void send();

  //! wait/check finish send
  /*!
    @todo return whether MPI_Wait is done successfully or not.
   */
  void wait_complete_send();
};

//! Receive Move in Particle from ONE target node
class MPIMoveIn {
 public:
  size_t size;                         //!< default number of particle
  int target_rank;                     //!< rank of target
  int target_id;                       //!< unit ID of target
  int unit_identifier;                 //!< unit ID of this node

  MPI_Comm mpi_comm_short;             //!< MPI_Comm for short

  size_t number_of_particle;           //!< actual number of particle
  
  MoveBufferHeader *recvbufferheader;  //!< Pointer of Header in Buffer
  Particle *recvparticle;              //!< Pointer of Particle in Buffer
  PotentialModel *recvtype;                 //!< Pointer of ParticleType in Buffer
  int *recvcellid;                 //!< Pointer of Cellid in Buffer
  std::vector<CovalentBondInfo::BondList *> recvbondlistarray;
  std::vector<int *> recvbondpackarray;
  char *recv_buffer;             //!< Pointer of Recv Paricle MPI Buffer
  size_t recv_buffer_size;       //!< Size of reserved Recv Paricle Buffer
  size_t recv_requested_size;    //!< Size of Recv Buffer calculate from requested number of Particle
  MPI_Request recv_request;      //!< used MPI_Irecv and MPI_Wait

  //! construct, default size is hard corded
  MPIMoveIn();

  //! construct, with unit ID, default size is hard coded
  explicit MPIMoveIn(int unitid, MPI_Comm short_comm);

  //! construct, with specified size
  MPIMoveIn(int unitid, size_t sz, MPI_Comm short_comm);

  //! construct, with unit ID and specified size
  MPIMoveIn(const MPIMoveIn& mi);

  //! assignment operator
  MPIMoveIn& operator=(const MPIMoveIn& mi);

  //! destructor
  ~MPIMoveIn();

  //! set particle data pointers in recv MPI buffer
  /*!
    @return set successfully or not
    @retval true size of received data <= size of preallocated buffer
    @retval false size of received data > size of preallocated buffer
    @todo 2nd case, realloc recv_buffer and not return false
   */
  bool set_receive_buffer_pointers();

  //! reserve receive MPI buffer
  /*!
    Construct receive_buffer.
    Set receive_buffer_pointers
   */
  void make_buffer();

  //! get receive Move in Particles from recv MPI buffer
  /*!
    @param[inout] movein append/overwrite to received particles
    @param[inout] moveintype potentialmodel of movein particle
    @param[inout] incellid cell-IDs where movein particles are store to.
    @return set successfully or not. 
    @retval true success
    @retval false Fail MPI_Wait or illeagal size. See set_receive_buffer_pointers.
    @note Received particlesets are stored continuously, without gap between particlesets. 
    @todo use fast copy.
    range.begin : index of request where store receive particle from 
  */
  bool getreceive(ParticleArray& movein, 
                  std::vector<PotentialModel>& moveintype,
                  std::vector<int>& incellid);
  bool getreceive(ParticleArray& movein,
                  std::vector<PotentialModel>& moveintype,
                  std::vector<int>& incellid,
                  std::vector<CovalentBondInfo::BondList>& inbondlistarray);
  //! request MPI asyncronous receive for Move in Particle
  /*!
    @note It must be called before calling getreceive*.
   */
  void prepare_receive();
};


//! set communication target
/*!
  @note defined for future use.
 */
class MPITargetNearestXYZ {
 public:
  void makeMoveOutTarget();
  void makeMoveInTarget();

};



//! Send/Receive any data to other MPI node
/*!
  Target node, target subset, target data are given.
  Not concern how to divide data or how to selected node
*/
class MPICommunicator {
 public:
  int unit_id;                        //!< unit ID of this node
  size_t psize;                       //!< numver of particle in subset
  std::map<int,int> unitid_to_rank;   //!< reverse map unit ID to rank

  MPI_Comm mpi_comm_short;             //!< MPI_Comm for short

  std::vector<MPISendParticleReceiveForce> sprf; //!< Send Own particle and Recieve reaction force
  std::vector<int> sprf_target_id;               //!< IDs of target ndoes
  int sprf_number_of_target;                     //!< number of send target nodes, equal to srpf_target_id.size()
  std::vector<int> setid;                        //!< set IDs of this node
  std::map<int,int> setid_to_index;              //!< reverse map setID to index
  std::vector< std::vector<int> > srpf_setid;    //!< set IDs which each target node required.
  std::vector<int> allsendsetid;                 //!< set IDs which required some target node.  
  int number_of_send_set;                        //!< number of particlesets which required some target node. @note for future use.

  std::vector<MPIReceiveParticleSendForce> rpsf; //!< Recv j-Particle and return reaction force
  std::vector<int> rpsf_target_id;               //!< IDs of target nodes
  int rpsf_number_of_target;                     //!< number of recv target nodes, equal to rpsf_target_id.size().
  std::vector< std::vector<int> > targetsetid;   //!< set IDs that each target node send to this node.
  std::vector<int> alltargetsetid;               //!< set IDs that some target nodes send to this node.
  int number_of_target_set;                      //!< number of particlesets from some nodes. @note for future use, used only in testcode.

  SpaceVector<int> recv_depth;                   //!< recv depth (longest distance to receive target node) in x,y,z direction
  SpaceVector<int> recv_direction;               //!< recv plus and/or minus x,yz direction, 0:bi-direction, +:from +, -:from -, if half shell, (0,0,+) (recv z from only +)
  std::vector<TypeRange> send_typerange;         //!< for NearestXYZ

  CommPattern send_recv_communication_pattern;   //!< Communication pattern for send/recv,  Direct or NearestXYZ

  //! target node ID for each direction when send_recv_communication_pattern=NearestXYZ

  int xplus_target_id;                     //!< target node ID for +X direction
  int xminus_target_id;                    //!< target node ID for -X direction
  int yplus_target_id;                     //!< target node ID for +Y direction
  int yminus_target_id;                    //!< target node ID for -Y direction
  int zplus_target_id;                     //!< target node ID for +Z direction
  int zminus_target_id;                    //!< target node ID for -Z direction

  std::vector<MPIMoveOut> move_out;              //!< Send Move Out
  std::vector<MPIMoveIn> move_in;                //!< Recv Move In
  std::vector<int> move_direct_target_id;        //!< IDs of move target nodes for Direct
  std::vector<int> move_out_target_id;           //!< IDs of move out target node for NearestXYZ
  std::vector<int> move_in_target_id;            //!< IDs of move in target node, for NearestXYZ
  int move_number_of_target;                     //!< number of Move target, = (depth*2+1)^3-1 for Direct, = 6 for NearestXYZ
  std::vector< std::vector<int> > move_setid;    //!< move particlesets IDs for each target node
  std::map<int,int> setid_to_move_target_index;  //!< reverse map of move_setid
  // buffers for move in/out
  ParticleArray all_move_out_particle;           //!< all move out particle
  std::vector<int> all_move_out_cellid;          //!< cell(subset)ID move out from
  std::vector<PotentialModel> all_move_out_type;      //!< Type of move out particle
  std::vector<CovalentBondInfo::BondList> all_move_out_bondlistarray; //!< all move out bondlistarray bind to each move out particle
  int number_of_move_out_particle;          //!< number of move_out_particle array, = move_number_of_target for Direct, = 6*depth for NearestXYZ
  int number_of_move_in_particle;          //!< number of move_out_particle array, = move_number_of_target for Direct, = 1 for NearestXYZ
  std::vector<ParticleArray> move_out_particle;  //!< move out partices for each node
  std::vector< std::vector<int> > move_out_cellid;  //!< cell ID of move out particles for each node
  std::vector< std::vector<PotentialModel> > move_out_type; //!< Type of move out particles for each node
  std::vector< std::vector<CovalentBondInfo::BondList> > move_out_bondlistarray; //!< BondList bind to move out particle for each node
  ParticleArray move_in_particle;                //!< move in particle
  std::vector<int> move_in_cellid;               //!< cell ID of move in particle
  std::vector<PotentialModel> move_in_type;           //!< Type of move in particle
  std::vector<CovalentBondInfo::BondList> move_in_bondlistarray; //!< bondlistarray bind to each move in particle

  /*!
    for move_communication_pattern = NearestXYZ
  */
  SpaceVector<int> node_position;     //!< position of this node
  SpaceVector<int> move_depth;  //!< move_in/out depth of x,y,z direction
  double move_surface_ratio;    //!< surface volume which contain move out candidate / volume
  /*
    move_imd_*[move_depth.x*2+move_depth.y*2+move_depth.z*2]
    for x+1,x+2,...x+move_depths.x,x-1,...x-move_depth.x,y+1,...
  */
  std::vector< ParticleArray > move_imd_particle;  //!< intemediate move particle
  std::vector< std::vector<int> > move_imd_cellid;   //!< cell ID of intemediate
  std::vector< std::vector<PotentialModel> > move_imd_type;    //!< Type of intemediate
  std::vector< std::vector<CovalentBondInfo::BondList> > move_imd_bondlistarray; //!< BondList bind to move out particle intemediate


  CommPattern move_communication_pattern;    //!< Communication pattern for move


  GeometryXYZ node_geometry;

  int lc;                                        // obsolete

  int timer_move;
  int timer_move_comm;

  MPICommunicator(){}

  //! constructor 
  MPICommunicator(const std::map<int,int>& uidtor, 
                  MPI_Comm short_comm,
                  std::vector<int> sprf_tid,
                  std::vector< std::vector<int> > srpf_set,
                  std::vector<int> sid,
                  std::vector<int> rpsf_tid,
                  std::vector< std::vector<int> > tsid,
                  std::vector<int> mtid,
                  std::vector< std::vector<int> > msid,
                  int ps, const int maxid, const int id,
                  SpaceVector<int> nodediv,
                  SpaceVector<int> celldiv,
                  SpaceVector<int> move_range,
                  double movesurfaceratio,
                  //              CommPattern srcp=Direct,
                  CommPattern srcp=NearestXYZ,
                  CommPattern mcp=NearestXYZ);

  //! copy constructor
  MPICommunicator(const MPICommunicator& comm);

  //  MPICommunicator& operator=(const MPICommunicator& comm);

  ~MPICommunicator();

  //! construct Send Particle, Recv Particle, Move Particle
  /**
     USE : rpsf_target_id, unit_id, send_recv_communication_pattern, 
     psize, unitid_to_rank, srpf_setid, targetsetid,
     move_communication_pattern, move_depth, move_direct_target_id,
     MODIFY : recv_depth, recv_direction, {x|y|z}{plus|minus}_target_id,
     sprf, rpsf, allsendsetid, 
     sprf_number_of_target, number_of_send_set, 
     rpsf_number_of_target, number_of_target_set, alltargetsetid,
     move_number_of_target, number_of_move_out_particle,
     number_of_move_in_particle, 
     move_imd_particle, move_imd_cellid, move_imd_type,
     move_out_target_id, move_in_target_id, 
     all_move_out_particle, all_move_out_cellid, all_move_out_type,
     move_out_particle, move_out_cellid, move_out_type,
     move_in_particle, move_in_cellid, move_in_type,
     move_out, move_in, 

     LOCAL : plus_depth, minus_depth, num_set, y_num_set, z_num_set, att, res,
     relative_target, 

     calculate depth
     warn depth imbalance
     
  */
  void construct_sendrecv();

  //! resize ParticleRange of j-Partice
  void resize_receive_data(std::vector<ParticleRange>& targetparticlerange);

  //! set own particle to each send buffer
  /*!
    for continuous particle
  */
  void setSendParticle(ParticleArray& myparticle,
                       std::vector<TypeRange>& typerangearray,
                       std::vector<CovalentBondInfo::BondList>& bondlistarray);

  //! set own particle to each send buffer
  /*!
    for particle with gap
  */
  template<class PA>
  void setSendParticlesubset_onlyPosition(const PA& myparticle,
                             const std::vector<TypeRange>& typerangearray);

  //! get received particle from each send buffer
  template<class GPA>
  void getReceiveParticle_onlyPosition(GPA& targetparticle);

  //! set own particle to each send buffer
  /*!
    for particle with gap
  */
  template<class PA>
  void setSendParticlesubset(const PA& myparticle,
                             const std::vector<TypeRange>& typerangearray,
                             const std::vector<CovalentBondInfo::BondList>& bondlistarray);

  //! get received particle from each send buffer
  template<class PA>
  void getReceiveParticle(PA& targetparticle,
                          std::vector<ParticleRange>& target_range,
                          std::vector<TypeRange>& targettyperange,
                          std::vector<CovalentBondInfo::BondList>& targetbond,
                          std::vector<int>& recvsetid,
                          std::map<int,int>& recvsetid_to_index);

  //! request transfer Particle
  void transferParticle();

  //! wait/check Particle Transfer
  void wait_complete_Particle();

  //! print some information
  void print();
  
  //! send own particle (no gap) and receive j-particle
  void exchangeParticleArray(ParticleArray& myparticle,
                             std::vector<TypeRange>& typerangearray,
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             ParticleArray& targetparticle,
                             std::vector<ParticleRange>& target_range,
                             std::vector<int>& recvsetid,
                             std::map<int,int>& recvsetid_to_index,
                             std::vector<TypeRange>& targettyperange,
                             std::vector<CovalentBondInfo::BondList>& targetbond);

  //! send own particle and receive j-particle with measurement of MPI Time
  void exchangeParticleArray(ParticleArray& myparticle,
                             std::vector<TypeRange>& typerangearray,
                             std::vector<CovalentBondInfo::BondList>& bondlistarray,
                             ParticleArray& targetparticle,
                             std::vector<ParticleRange>& target_range,
                             std::vector<int>& recvsetid,
                             std::map<int,int>& recvsetid_to_index,
                             std::vector<TypeRange>& targettyperange,
                             std::vector<CovalentBondInfo::BondList>& targetbond,
                             double& alltime, double& rectime);

  //! set j-force to each send
  void setSendForce(const ForceArray& sendforce,
                    const std::vector<ParticleRange>& send_range);
  void setSendForce_with_index(const ForceArray& sendforce,
                               const std::vector<ParticleRange>& send_range,
                               const std::vector<int>& forceindexset);
  void setSendForce_indexed(const ForceArray& sendforce,
                            const std::vector<ParticleRange>& send_range,
                            const std::vector<int>& forceindexset);

  //! get received force
  void getReceiveForcesubset(ForceArray& recvforce);
  void getReceiveForcesubset_with_index(ForceArray& recvforce);
  void getReceiveForcesubset_indexed(ForceArray& recvforce);

  //! request transfer Force
  void transferForce();

  //! wait/check Force Transfer
  void wait_complete_Force();

  //! send own particle (with gap) and receive j-particle
  template<class PA, class GPA>
  void exchangeParticleArraysubset_onlyPosition_top_half(const PA& myparticle,
                                   const std::vector<TypeRange>& typerangearray,
                                   const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                   GPA& targetparticle,
                                   const std::vector<ParticleRange>& target_range,
                                   const std::vector<int>& recvsetid,
                                   const std::map<int,int>& recvsetid_to_index,
                                   const std::vector<TypeRange>& targettyperange,
                                   const std::vector<CovalentBondInfo::BondList>& targetbond);
  template<class PA, class GPA>
  void exchangeParticleArraysubset_onlyPosition_bottom_half(const PA& myparticle,
                                   const std::vector<TypeRange>& typerangearray,
                                   const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                   GPA& targetparticle,
                                   const std::vector<ParticleRange>& target_range,
                                   const std::vector<int>& recvsetid,
                                   const std::map<int,int>& recvsetid_to_index,
                                   const std::vector<TypeRange>& targettyperange,
                                   const std::vector<CovalentBondInfo::BondList>& targetbond);
  template<class PA, class GPA>
  void exchangeParticleArraysubset_onlyPosition(const PA& myparticle,
                                   const std::vector<TypeRange>& typerangearray,
                                   const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                   GPA& targetparticle,
                                   const std::vector<ParticleRange>& target_range,
                                   const std::vector<int>& recvsetid,
                                   const std::map<int,int>& recvsetid_to_index,
                                   const std::vector<TypeRange>& targettyperange,
                                   const std::vector<CovalentBondInfo::BondList>& targetbond);

  template<class PA, class GPA>
  void exchangeParticleArraysubset_top_half(const PA& myparticle,
                                   const std::vector<TypeRange>& typerangearray,
                                   const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                   GPA& targetparticle,
                                   std::vector<ParticleRange>& target_range,
                                   std::vector<int>& recvsetid,
                                   std::map<int,int>& recvsetid_to_index,
                                   std::vector<TypeRange>& targettyperange,
                                   std::vector<CovalentBondInfo::BondList>& targetbond);
  template<class PA, class GPA>
  void exchangeParticleArraysubset_bottom_half(const PA& myparticle,
                                   const std::vector<TypeRange>& typerangearray,
                                   const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                   GPA& targetparticle,
                                   std::vector<ParticleRange>& target_range,
                                   std::vector<int>& recvsetid,
                                   std::map<int,int>& recvsetid_to_index,
                                   std::vector<TypeRange>& targettyperange,
                                   std::vector<CovalentBondInfo::BondList>& targetbond);
  //! send own particle (with gap) and receive j-particle
  template<class PA, class GPA>
  void exchangeParticleArraysubset(const PA& myparticle,
                                   const std::vector<TypeRange>& typerangearray,
                                   const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                   GPA& targetparticle,
                                   std::vector<ParticleRange>& target_range,
                                   std::vector<int>& recvsetid,
                                   std::map<int,int>& recvsetid_to_index,
                                   std::vector<TypeRange>& targettyperange,
                                   std::vector<CovalentBondInfo::BondList>& targetbond);

  //! send own particle (with gap) and receive j-particle with measurement MPI time
  void exchangeParticleArraysubset(const ParticleArray& myparticle,
                                   const std::vector<TypeRange>& typerangearray,
                                   const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                   ParticleArray& targetparticle,
                                   std::vector<ParticleRange>& target_range,
                                   std::vector<int>& recvsetid,
                                   std::map<int,int>& recvsetid_to_index,
                                   std::vector<TypeRange>& targettyperange,
                                   std::vector<CovalentBondInfo::BondList>& targetbond,
                                   double& alltime, double& rectime);

  //! send j-force (with gap) and receive own force
  void exchangeForceArraysubset(ForceArray& sendforce,
                                std::vector<ParticleRange>& send_range,
                                ForceArray& recvforce);
  void exchangeForceArraysubset_with_index(ForceArray& sendforce,
                                           std::vector<ParticleRange>& send_range,
                                           const std::vector<int>& forceindexset,
                                           ForceArray& recvforce);
  void exchangeForceArraysubset_indexed(ForceArray& sendforce,
                                        std::vector<ParticleRange>& send_range,
                                        const std::vector<int>& forceindexset,
                                        ForceArray& recvforce);
  
  //! send j-force (with gap) and receive own force with measurement MPI time
  void exchangeForceArraysubset(ForceArray& sendforce,
                                std::vector<ParticleRange>& send_range,
                                ForceArray& recvforce,
                                double& alltime, double& rectime);

  template<class PA, class GPA>
  void particle_send_recv_onlyPosition_Direct_top_half(const PA& myparticle,
                          const std::vector<TypeRange>& typerangearray,
                          const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                          GPA& targetparticle,
                          const std::vector<ParticleRange>& target_range,
                          const std::vector<int>& recvsetid,
                          const std::map<int,int>& recvsetid_to_index,
                          const std::vector<TypeRange>& targettyperange,
                          const std::vector<CovalentBondInfo::BondList>& targetbond);
  template<class PA, class GPA>
  void particle_send_recv_onlyPosition_Direct_bottom_half(const PA& myparticle,
                          const std::vector<TypeRange>& typerangearray,
                          const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                          GPA& targetparticle,
                          const std::vector<ParticleRange>& target_range,
                          const std::vector<int>& recvsetid,
                          const std::map<int,int>& recvsetid_to_index,
                          const std::vector<TypeRange>& targettyperange,
                          const std::vector<CovalentBondInfo::BondList>& targetbond);
  template<class PA, class GPA>
  void particle_send_recv_onlyPosition_Direct(const PA& myparticle,
                          const std::vector<TypeRange>& typerangearray,
                          const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                          GPA& targetparticle,
                          const std::vector<ParticleRange>& target_range,
                          const std::vector<int>& recvsetid,
                          const std::map<int,int>& recvsetid_to_index,
                          const std::vector<TypeRange>& targettyperange,
                          const std::vector<CovalentBondInfo::BondList>& targetbond);
  template<class GPA>
  void particle_send_recv_one_axis_onlyPosition(int depth, int recv_direction,
                                   int pf_plus, int pf_minus,
                                   int& target,
                                   GPA& targetparticle,
                                   const std::vector<ParticleRange>& target_range,
                                   const std::vector<int>& recvsetid,
                                   const std::map<int,int>& recvsetid_to_index,
                                   const std::vector<TypeRange>& targettyperange,
                                   const std::vector<CovalentBondInfo::BondList>& targetbond);
  template<class PA, class GPA>
  void particle_send_recv_onlyPosition_NearestXYZ(const PA& myparticle,
                          const std::vector<TypeRange>& typerangearray,
                          const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                          GPA& targetparticle,
                          const std::vector<ParticleRange>& target_range,
                          const std::vector<int>& recvsetid,
                          const std::map<int,int>& recvsetid_to_index,
                          const std::vector<TypeRange>& targettyperange,
                          const std::vector<CovalentBondInfo::BondList>& targetbond);

  template<class PA, class GPA>
  void particle_send_recv_Direct_top_half(const PA& myparticle,
                          const std::vector<TypeRange>& typerangearray,
                          const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                          GPA& targetparticle,
                          std::vector<ParticleRange>& target_range,
                          std::vector<int>& recvsetid,
                          std::map<int,int>& recvsetid_to_index,
                          std::vector<TypeRange>& targettyperange,
                          std::vector<CovalentBondInfo::BondList>& targetbond);
  template<class PA, class GPA>
  void particle_send_recv_Direct_bottom_half(const PA& myparticle,
                          const std::vector<TypeRange>& typerangearray,
                          const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                          GPA& targetparticle,
                          std::vector<ParticleRange>& target_range,
                          std::vector<int>& recvsetid,
                          std::map<int,int>& recvsetid_to_index,
                          std::vector<TypeRange>& targettyperange,
                          std::vector<CovalentBondInfo::BondList>& targetbond);
  template<class PA, class GPA>
  void particle_send_recv_Direct(const PA& myparticle,
                          const std::vector<TypeRange>& typerangearray,
                          const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                          GPA& targetparticle,
                          std::vector<ParticleRange>& target_range,
                          std::vector<int>& recvsetid,
                          std::map<int,int>& recvsetid_to_index,
                          std::vector<TypeRange>& targettyperange,
                          std::vector<CovalentBondInfo::BondList>& targetbond);

  template<class GPA>
  void particle_send_recv_one_axis(int depth, int recv_direction,
                                   int pf_plus, int pf_minus,
                                   int& target,
                                   GPA& targetparticle,
                                   std::vector<ParticleRange>& target_range,
                                   std::vector<int>& recvsetid,
                                   std::map<int,int>& recvsetid_to_index,
                                   std::vector<TypeRange>& targettyperange,
                                   std::vector<CovalentBondInfo::BondList>& targetbond);

template<class PA, class GPA>
void particle_send_recv_NearestXYZ(const PA& myparticle,
                          const std::vector<TypeRange>& typerangearray,
                          const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                          GPA& targetparticle,
                          std::vector<ParticleRange>& target_range,
                          std::vector<int>& recvsetid,
                          std::map<int,int>& recvsetid_to_index,
                          std::vector<TypeRange>& targettyperange,
                          std::vector<CovalentBondInfo::BondList>& targetbond);

  void force_send_recv_one_axis(int depth, int recv_direction,
                                int pf_plus, int pf_minus,
                                int& target,
                                ForceArray& sendforce,
                                const std::vector<ParticleRange>& send_range,
                                ForceArray& recvforce);

  void force_send_recv_Direct(const ForceArray& sendforce,
                              const std::vector<ParticleRange>& send_range,
                              ForceArray& recvforce);
  void force_send_recv_Direct_with_index(const ForceArray& sendforce,
                                         const std::vector<ParticleRange>& send_range,
                                         const std::vector<int>& forceindexset,
                                         ForceArray& recvforce);
  void force_send_recv_Direct_indexed(const ForceArray& sendforce,
                                      const std::vector<ParticleRange>& send_range,
                                      const std::vector<int>& forceindexset,
                                      ForceArray& recvforce);
  void force_send_recv_NearestXYZ(ForceArray& sendforce,
                       const std::vector<ParticleRange>& send_range,
                       ForceArray& recvforce);

  //! set move out particle
  void setSendMove();
   
  //! request transfer Move particle
  void transferMove();

  //! get received Move in particle
  void getReceiveMove();
  
  //! wait/check move out transfer
  void wait_complete_Move();

  template<CommPattern>
  void move_send_recv();

  template<CommPattern>
  void distribute_move_out();

  void marge_move(ParticleArray& particle, 
                  std::vector<PotentialModel> &rangetype,
                  std::vector<int> &cellid,
                  std::vector<CovalentBondInfo::BondList>& bondlist,
                  ParticleArray& add_particle, 
                  std::vector<PotentialModel> &add_rangetype,
                  std::vector<int> &add_cellid,
                  std::vector<CovalentBondInfo::BondList>& add_bondlist);
  void marge_move(ParticleArray& particle, 
                  std::vector<PotentialModel> &rangetype,
                  ParticleArray& add_particle, 
                  std::vector<PotentialModel> &add_rangetype);

  //! move inside node, send move out particle and receive move in particle
  template<typename PA>
  void move_particle(PA& particlearray,
                     std::vector<TypeRange>& typerangearray);

};



#endif // MPIPARALLEL_H

