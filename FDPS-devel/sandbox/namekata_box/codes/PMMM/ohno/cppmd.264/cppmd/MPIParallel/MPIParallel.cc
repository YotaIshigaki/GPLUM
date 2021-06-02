#include <cstring>
#include "MPIParallel.h"

#ifdef OVERLAP
int MPICHECKER::count;
MPI_Request* MPICHECKER::reqs;
int MPICHECKER::index;
int MPICHECKER::flag;
MPI_Status* MPICHECKER::status;
bool MPICHECKER::done;
void MPICHECKER::mpi_checker(){
  //  printf("MPICHECKER::count=%d\n",count);
  if(!done){
    MPI_Testall(count,reqs,&flag,status);
    //    printf("MPI_Testall flag %d\n",flag);
    if(flag){
      done = true;
      //      printf("MPI_Testall done\n");
    }
  }
  /*
    {
    printf("mpi_testall");
    for(int i=0;i<count;i++){
      printf(" : %d %d",i,status[i]);
    }
    printf("\n");
    }
  */
}
#endif
//#define TRANSFER_RAW_BONDLIST

// suport functions calculate size of buffer
//! return size of Positoins Buffer in Byte
size_t calc_position_buffer_size(int number_of_position)
{
  return ( sizeof(PositionBufferHeader)
           + sizeof(Position)*number_of_position);
}
size_t calc_position_buffer_size(PositionBufferHeader& header)
{
  return ( sizeof(PositionBufferHeader)
           + sizeof(Position)*header.NumberOfPosition);
}

//! return size of Particles Buffer in Byte
/*!
  @todo : fix incorrect size of bondlist
 */
size_t calc_particle_buffer_size(int number_of_position,
                                 int number_of_charge,
                                 int number_of_atomtype,
                                 int number_of_atomid,
                                 int number_of_particleset,
                                 int number_of_bondlist )
{
  return ( sizeof(ParticleBufferHeader)
           + sizeof(Position)*number_of_position
           + sizeof(double)*number_of_charge
           + sizeof(Atomtype)*number_of_atomtype
           + sizeof(AtomID)*number_of_atomid
           + sizeof(int)*number_of_particleset
           + sizeof(TypeRange)*number_of_particleset
           + sizeof(CovalentBondInfo::BondList)*number_of_bondlist
           );
}

//! return size of Particles Buffer in Byte
size_t calc_particle_buffer_size(ParticleBufferHeader& header)
{
  return ( sizeof(ParticleBufferHeader)
           + sizeof(Position)*header.NumberOfPosition
           + sizeof(double)*header.NumberOfCharge
           + sizeof(Atomtype)*header.NumberOfAtomType
           + sizeof(AtomID)*header.NumberOfAtomID
           + sizeof(int)*header.NumberOfParticleSet
           + sizeof(TypeRange)*header.NumberOfParticleSet
           + sizeof(CovalentBondInfo::BondList)*header.NumberOfBondlist
           );
}

//! return size of Force Buffer in Byte
size_t calc_force_buffer_size(int number_of_force)
{
  return ( sizeof(ForceBufferHeader)
           + sizeof(Force)*number_of_force
           );
}
size_t calc_force_buffer_size_with_index(int number_of_force)
{
  return ( sizeof(ForceBufferHeader)
           + sizeof(Force)*number_of_force
           + sizeof(int)*number_of_force
           );
}


//! return size of Move Particle Buffer in Byte
/*!
  @todo : fix incorrect size of bondlist
 */
size_t calc_move_buffer_size(int number_of_move)
{
  return ( sizeof(MoveBufferHeader)
           + sizeof(Particle)*number_of_move
           + sizeof(PotentialModel)*number_of_move
           + sizeof(int)*number_of_move
           );
}

size_t calc_move_buffer_size(int number_of_move,
                             std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  size_t bla_size = 0;

  for(std::vector<CovalentBondInfo::BondList>::size_type i=0;
      i<bondlistarray.size();i++){
#ifdef TRANSFER_RAW_BONDLIST
    bla_size += bondlistarray[i].size();
#else
    bla_size += bondlistarray[i].size_of_packed_int();
#endif
  }
  return ( sizeof(MoveBufferHeader)
           + sizeof(Particle)*number_of_move
           + sizeof(PotentialModel)*number_of_move
           + sizeof(int)*number_of_move
           + bla_size
           );
}

size_t calc_move_buffer_size(MoveBufferHeader& header)
{
  return ( sizeof(MoveBufferHeader)
           + sizeof(Particle)*header.NumberOfParticle
           + sizeof(PotentialModel)*header.NumberOfParticle
           + sizeof(int)*header.NumberOfParticle
           + header.SizeOfBondlist
           );
}
size_t calc_move_buffer_size(MoveBufferHeader& header,
                             std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  return calc_move_buffer_size(header.NumberOfParticle,bondlistarray);
}

#if 1
void set_buffer_pointers(char *origin,
                         ParticleBufferHeader* &header,
                         Position* &position,
                         double* &charge,
                         Atomtype* &atomtype,
                         AtomID* &atomid,
                         int* &setid,
                         TypeRange* &typerange,
                         CovalentBondInfo::BondList* &bond)
{
  void *currentpoint = (void *)origin;
  
  header = static_cast<ParticleBufferHeader *>(currentpoint);
  currentpoint = &(header[1]);
  position = static_cast<Position *>(currentpoint);
  currentpoint = &(position[header->NumberOfPosition]);
  charge = static_cast<double *>(currentpoint);
  currentpoint = &(charge[header->NumberOfCharge]);
  atomtype = static_cast<Atomtype *>(currentpoint);
  currentpoint = &(atomtype[header->NumberOfAtomType]);
  atomid = static_cast<AtomID *>(currentpoint);
  currentpoint = &(atomid[header->NumberOfAtomID]);
  setid = static_cast<int *>(currentpoint);
  currentpoint = &(setid[header->NumberOfParticleSet]);
  typerange = static_cast<TypeRange *>(currentpoint);
  currentpoint = &(typerange[header->NumberOfParticleSet]);
  bond = static_cast<CovalentBondInfo::BondList *>(currentpoint);
}
#endif

size_t copy_position_buffer(char *source, char *dest)
{
  PositionBufferHeader *header = static_cast<PositionBufferHeader *>((void *)(source));
  size_t buffer_size = calc_position_buffer_size(header[0]);
  memcpy((void*)(dest),(void*)(source),buffer_size);
  return buffer_size;
}

static size_t copy_particle_buffer(char *source, char *dest)
{
#if 1
  ParticleBufferHeader *header = static_cast<ParticleBufferHeader *>((void *)(source));
  size_t buffer_size = calc_particle_buffer_size(header[0]);
  memcpy((void*)(dest),(void*)(source),buffer_size);
#else
  ParticleBufferHeader *s_header;
  Position *s_position;
  double *s_charge;
  Atomtype *s_atomtype;
  AtomID *s_atomid;
  int *s_setid;
  TypeRange *s_typerange;
  CovalentBondInfo::BondList *s_bond;
  ParticleBufferHeader *d_header;
  Position *d_position;
  double *d_charge;
  Atomtype *d_atomtype;
  AtomID *d_atomid;
  int *d_setid;
  TypeRange *d_typerange;
  CovalentBondInfo::BondList *d_bond;
  set_buffer_pointers(source, 
                      s_header, s_position, s_charge, s_atomtype, s_atomid,
                      s_setid, s_typerange, s_bond);
  d_header =  static_cast<ParticleBufferHeader *>((void *)(dest));
  memcpy((void*)(d_header),(void*)(s_header),sizeof(ParticleBufferHeader));
  set_buffer_pointers(dest, 
                      d_header, d_position, d_charge, d_atomtype, d_atomid,
                      d_setid, d_typerange, d_bond);
  for(size_t i=0;i<s_header->NumberOfPosition;i++){
    d_position[i] = s_position[i];
  }
  for(size_t i=0;i<s_header->NumberOfCharge;i++){
    d_charge[i] = s_charge[i];
  }
  for(size_t i=0;i<s_header->NumberOfAtomType;i++){
    d_atomtype[i] = s_atomtype[i];
  }
  for(size_t i=0;i<s_header->NumberOfAtomID;i++){
    d_atomid[i] = s_atomid[i];
  }
  for(size_t i=0;i<s_header->NumberOfParticleSet;i++){
    d_setid[i] = s_setid[i];
  }
  for(size_t i=0;i<s_header->NumberOfParticleSet;i++){
    d_typerange[i] = s_typerange[i];
  }
  for(size_t i=0;i<s_header->NumberOfBondlist;i++){
    d_bond[i] = s_bond[i];
  }
  size_t buffer_size = calc_particle_buffer_size(d_header[0]);
  //  std::cout << "pos " << s_header->NumberOfPosition << " set " << s_header->NumberOfParticleSet << std::endl;
#endif
  return buffer_size;
}

/*!
  Methods of MPISendParticleReceiveForce
 */
//! construct, default size is hard corded
MPISendParticleReceiveForce::MPISendParticleReceiveForce()
  : size(buffersize),
    mpi_comm_short(MPI_COMM_WORLD),
    number_of_stage(1),
    stage(0)
{
  number_of_particle.resize(number_of_stage,0);
  number_of_set.resize(number_of_stage,0);
  receive_buffer = (char *)NULL;
  send_buffer = (char *)NULL;
}

MPISendParticleReceiveForce::MPISendParticleReceiveForce(int unitid, 
                                                         MPI_Comm short_comm)
  : size(buffersize),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm),
    number_of_stage(1),
    stage(0)
{
  number_of_particle.resize(number_of_stage,0);
  number_of_set.resize(number_of_stage,0);
  receive_buffer = (char *)NULL;
  send_buffer = (char *)NULL;
}

//! construct, with specified size
MPISendParticleReceiveForce::MPISendParticleReceiveForce(int unitid, 
                                                         size_t sz, 
                                                         MPI_Comm short_comm)
  : size(sz),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm),
    number_of_stage(1),
    stage(0)
{
  number_of_particle.resize(number_of_stage,0);
  number_of_set.resize(number_of_stage,0);
  receive_buffer = (char *)NULL;
  send_buffer = (char *)NULL;
}

MPISendParticleReceiveForce::MPISendParticleReceiveForce(int unitid, 
                                                         size_t sz, 
                                                         MPI_Comm short_comm,
                                                         int num_stage)
  : size(sz),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm),
    number_of_stage(num_stage),
    stage(0)
{
  number_of_particle.resize(number_of_stage,0);
  number_of_set.resize(number_of_stage,0);
  receive_buffer = (char *)NULL;
  send_buffer = (char *)NULL;
}

MPISendParticleReceiveForce::MPISendParticleReceiveForce(const MPISendParticleReceiveForce& sprf)
  : size(sprf.size),
    target_rank(sprf.target_rank),
    target_id(sprf.target_id),
    unit_identifier(sprf.unit_identifier),
    mpi_comm_short(sprf.mpi_comm_short),
    number_of_stage(sprf.number_of_stage),
    number_of_particle(sprf.number_of_particle),
    number_of_set(sprf.number_of_set),
    stage(sprf.stage)
{
  if(sprf.send_buffer!=(char *)NULL){
    make_buffer(sprf.number_of_set[stage]);
  }else{
    receive_buffer = (char *)NULL;
    send_buffer = (char *)NULL;
  }
  //  make_buffer(sprf.number_of_set);
}

MPISendParticleReceiveForce& MPISendParticleReceiveForce::operator=(const MPISendParticleReceiveForce& sprf)
{
  if(this != &sprf){
    size = sprf.size;
    target_rank = sprf.target_rank;
    target_id = sprf.target_id;
    unit_identifier = sprf.unit_identifier;
    mpi_comm_short = sprf.mpi_comm_short;
    number_of_stage = sprf.number_of_stage;
    number_of_particle = sprf.number_of_particle;
    number_of_set = sprf.number_of_set;
    stage = sprf.stage;
    if(sprf.send_buffer!=(char *)NULL){
      make_buffer(sprf.number_of_set[stage]);
    }else{
      receive_buffer = (char *)NULL;
      send_buffer = (char *)NULL;
    }
  }
  return *this;
}

MPISendParticleReceiveForce::~MPISendParticleReceiveForce()
{
  delete [] receive_buffer;
  delete [] send_buffer;
}


//! set position data pointers in send MPI buffer
/*! return true  : required buffer size <= preallocated buffer size
  return false : required buffer size > preallocated buffer size
  TODO :
  2nd case, realloc send_buffer and not return false
*/
bool
MPISendParticleReceiveForce::set_send_posbuffer_pointers(int number_of_position)
{
  send_request_size = calc_position_buffer_size(number_of_position);
  if( send_request_size > send_buffer_size ){
    std::cout << "send_request_size " << send_request_size << " > send_buffer_size " << send_buffer_size << " " << number_of_position << " position" << std::endl;
    return false;
  }
  void *currentpoint = send_buffer;
  sendposbufferheader = static_cast<PositionBufferHeader *>(currentpoint);
  sendposbufferheader->UnitID = unit_identifier;
  sendposbufferheader->NumberOfPosition = number_of_position;
  currentpoint = &(sendposbufferheader[1]);
  sendposition = static_cast<Position *>(currentpoint);

  return true;
}

//! set particle data pointers in send MPI buffer
/*! return true  : required buffer size <= preallocated buffer size
  return false : required buffer size > preallocated buffer size
  @todo 2nd case, realloc send_buffer and not return false
  @todo size of bondlist is not correct.
*/
bool
MPISendParticleReceiveForce::set_send_buffer_pointers(int number_of_position,
                                                      int number_of_charge,
                                                      int number_of_atomtype,
                                                      int number_of_atomid,
                                                      int number_of_particleset,
                                                      int number_of_bondlist )
{
  send_request_size = calc_particle_buffer_size(number_of_position,
                                                number_of_charge,
                                                number_of_atomtype,
                                                number_of_atomid,
                                                number_of_particleset,
                                                number_of_bondlist);
  if( send_request_size > send_buffer_size ){
    std::cout << "send_request_size " << send_request_size << " > send_buffer_size " << send_buffer_size << " " << number_of_position << " particle" << std::endl;
    return false;
  }
  void *currentpoint = send_buffer;
  sendbufferheader = static_cast<ParticleBufferHeader *>(currentpoint);
  sendbufferheader->UnitID = unit_identifier;
  sendbufferheader->NumberOfPosition = number_of_position;
  sendbufferheader->NumberOfCharge = number_of_charge;
  sendbufferheader->NumberOfAtomType = number_of_atomtype;
  sendbufferheader->NumberOfAtomID = number_of_atomid;
  sendbufferheader->NumberOfParticleSet = number_of_particleset; 
  sendbufferheader->NumberOfBondlist = number_of_bondlist;
  currentpoint = &(sendbufferheader[1]);
  sendposition = static_cast<Position *>(currentpoint);
  currentpoint = &(sendposition[number_of_position]);
  sendcharge = static_cast<double *>(currentpoint);
  currentpoint = &(sendcharge[number_of_charge]);
  sendatomtype = static_cast<Atomtype *>(currentpoint);
  currentpoint = &(sendatomtype[number_of_atomtype]);
  sendatomid = static_cast<AtomID *>(currentpoint);
  currentpoint = &(sendatomid[number_of_atomid]);
  sendsetid = static_cast<int *>(currentpoint);
  currentpoint = &(sendsetid[number_of_particleset]);
  sendtyperange = static_cast<TypeRange *>(currentpoint);
  currentpoint = &(sendtyperange[number_of_particleset]);
  sendbond = static_cast<CovalentBondInfo::BondList *>(currentpoint);

  return true;
}

//! set force data pointers in receive MPI buffer
/*! return true  : expected buffer size <= preallocated buffer size
  return false : expected buffer size > preallocated buffer size
  TODO :
  2nd case, realloc send_buffer and not return false
*/
bool
MPISendParticleReceiveForce::set_receive_buffer_pointers(int number_of_force)
{
  receive_expected_size = calc_force_buffer_size(number_of_force);
  if( receive_expected_size > receive_buffer_size ){
    std::cout << "receive_expected_size > receive_buffer_size" << std::endl;
    return false;
  }
  void *currentpoint = receive_buffer;
  receivebufferheader = static_cast<ForceBufferHeader *>(currentpoint);
  currentpoint = &(receivebufferheader[1]);
  ReceiveForceArray = static_cast<Force *>(currentpoint);
  return true;
}
// with index, bonded only
bool
MPISendParticleReceiveForce::set_receive_buffer_pointers_with_index(int number_of_force)
{
  receive_expected_size = calc_force_buffer_size(number_of_force);
  if( receive_expected_size > receive_buffer_size ){
    std::cout << "receive_expected_size > receive_buffer_size" << std::endl;
    return false;
  }
  void *currentpoint = receive_buffer;
  receivebufferheader = static_cast<ForceBufferHeader *>(currentpoint);
  currentpoint = &(receivebufferheader[1]);
  ReceiveForceArray = static_cast<Force *>(currentpoint);
  currentpoint = &(ReceiveForceArray[number_of_force]);
  ReceiveForceIndex = static_cast<int *>(currentpoint);
  return true;
}

//! reserve send and receive MPI buffer and tempral set information
void
MPISendParticleReceiveForce::make_buffer(const size_t& num_set)
{
  number_of_set[stage] = num_set;
  number_of_particle[stage] = size;
  send_buffer_size = calc_particle_buffer_size(size, size, size, size, num_set, 1);
  //  delete [] send_buffer;
  send_buffer = new char[send_buffer_size];
  set_send_posbuffer_pointers(size);
  set_send_buffer_pointers(size,size,size,size,num_set,0);
  sent_typerange.resize(num_set);
  set_index.resize(num_set);
  
  receive_buffer_size = calc_force_buffer_size(size);
  //  delete [] receive_buffer;
  receive_buffer = new char[receive_buffer_size];
  set_receive_buffer_pointers(size);

  if((unit_identifier==0)&&(DebugLog::verbose>1)){
    std::cout << "MPISendParticleReceiveForce::make_buffer particle " << number_of_particle[stage] << " num_set " << num_set << " buffer size send " << send_buffer_size << " receive " << receive_buffer_size << std::endl;
  }

}

//! set send Positions to send MPI buffer
/*!
  Paritcles have no gap between subset
*/
// it slow copy!
bool
MPISendParticleReceiveForce::setsendpos(ParticleArray& request)
{
  number_of_particle[stage] = request.size();
  // in this version send all data : number of position, charge, atomtype, atomid = number of particle
  // in this version ignore bond !!!!!
  size_t num_position  = number_of_particle[stage];
  if( set_send_posbuffer_pointers(num_position) == true ) {
    for(size_t i=0;i<num_position;i++){
      sendposition[i] = request[i].position;
      // sendposition[i].x = request[i].position.x;
      // sendposition[i].y = request[i].position.y;
      // sendposition[i].z = request[i].position.z;
    }
    return true;
  }
  return false;
}

//! set send Particles to send MPI buffer
/*!
  Paritcles have no gap between subset
*/
// it slow copy!
bool
MPISendParticleReceiveForce::setsend(ParticleArray& request,
                                     std::vector<TypeRange>& rq_tr,
                                     std::vector<CovalentBondInfo::BondList> &rq_bond, 
                                     std::vector<int>& setid)
{
  number_of_particle[stage] = request.size();
  // in this version send all data : number of position, charge, atomtype, atomid = number of particle
  // in this version ignore bond !!!!!
  size_t num_position  = number_of_particle[stage];
  size_t num_charge    = number_of_particle[stage];
  size_t num_atomtype  = number_of_particle[stage];
  size_t num_atomid  = number_of_particle[stage];
  size_t num_particleset = setid.size();
  size_t num_bond      = 0;
  if( set_send_buffer_pointers(num_position, num_charge, 
                               num_atomtype, num_atomid,
                               num_particleset, num_bond) == true ) {
    for(size_t i=0;i<num_position;i++){
      sendposition[i] = request[i].position;
    }
    for(size_t i=0;i<num_charge;i++){
      sendcharge[i] = request[i].charge;
    }
    for(size_t i=0;i<num_atomtype;i++){
      sendatomtype[i] = request[i].atomtype;
    }
    for(size_t i=0;i<num_atomid;i++){
      sendatomid[i] = request[i].atomid;
    }
    for(size_t i=0;i<num_particleset;i++){
      sendsetid[i] = setid[i];
    }
    for(size_t i=0;i<num_particleset;i++){
      sendtyperange[i] = rq_tr[i];
    }
    for(size_t i=0;i<num_bond;i++){
      sendbond[i] = rq_bond[i];
    }
    return true;
  }
  return false;
}

//! set send Positinons to send MPI buffer
/*!
  Paritcles have gap between subset ( reserve max size of each subset )
  @note : Is it correct for multi stage?
*/
// it slow copy!
template<class PA>
bool
MPISendParticleReceiveForce::setsendsubsetpos(const PA& request,
                                              const std::vector<TypeRange>& rq_tr,
                                              const std::vector<int>& setid,
                                              const std::vector<int>& rq_set)
{
  bool stat;
//  double  (* const __restrict sendposition)[3] = (double  (* const __restrict)[3]) (this->sendposition);
  size_t num_position  = number_of_particle[stage];
  size_t num_particleset = number_of_set[stage];
  if( set_send_posbuffer_pointers(num_position) == true ) {
    //      std::cout << "send particle " << num_position << std::endl;
    int i=0;
    for(size_t s=0;s<num_particleset;s++){
      TypeRange tr = rq_tr[set_index[s]];
      int offset = i-tr.begin;
      for(int si=tr.begin;si<tr.end;si++,i++){
        // sendposition[i] = request[si].position;
        sendposition[i] = getpos(request,si);
        // sendposition[i].x = getpos(request,si).x;
        // sendposition[i].y = getpos(request,si).y;
        // sendposition[i].z = getpos(request,si).z;
        //sendposition[i][0] = getpos(request,si).x;
        //sendposition[i][1] = getpos(request,si).y;
        //sendposition[i][2] = getpos(request,si).z;
      }
      tr.shift(offset);
    }
    stat = true;
  }else{
    stat =  false;
  }
  return stat;
}
template
bool
MPISendParticleReceiveForce::setsendsubsetpos(const ParticleArray& request,
                                              const std::vector<TypeRange>& rq_tr,
                                              const std::vector<int>& setid,
                                              const std::vector<int>& rq_set);
template
bool
MPISendParticleReceiveForce::setsendsubsetpos(const CombinedParticleArray& request,
                                              const std::vector<TypeRange>& rq_tr,
                                              const std::vector<int>& setid,
                                              const std::vector<int>& rq_set);

//! set send Particles to send MPI buffer
/*!
  Paritcles have gap between subset ( reserve max size of each subset )
*/
// it slow copy!
template<class PA>
bool
MPISendParticleReceiveForce::setsendsubset(const PA& request,
                                           const std::vector<TypeRange>& rq_tr,
                                           const std::vector<CovalentBondInfo::BondList> &rq_bond,
                                           const std::vector<int>& setid,
                                           const std::vector<int>& rq_set)
{
  bool stat;

  set_index.clear();
  sent_typerange.clear();
  //    number_of_particle = request.size();
  number_of_particle[stage] = 0;
  for(std::vector<int>::const_iterator it = rq_set.begin();
      it != rq_set.end(); ++it){
    int sindex=0;
    while(setid[sindex]!=(*it)){sindex++;}
    set_index.push_back(sindex);
    sent_typerange.push_back(rq_tr[sindex]);
    number_of_particle[stage] += rq_tr[sindex].end - rq_tr[sindex].begin;
  }
  number_of_set[stage] = rq_set.size();
  //  std::cout << "set " << number_of_particle[stage] << " send particles to " << target_id << std::endl;
  // in this version send all data : number of position, charge, atomtype, atomid = number of particle
  // in this version ignore bond !!!!!
  //  std::cout << "setsend number_of_particle[" << stage << "] " << number_of_particle[stage] << std::endl;
  size_t num_position  = number_of_particle[stage];
  size_t num_charge    = number_of_particle[stage];
  size_t num_atomtype  = number_of_particle[stage];
  size_t num_atomid    = number_of_particle[stage];
  size_t num_particleset = number_of_set[stage];
  size_t num_bond      = 0;
  if( set_send_buffer_pointers(num_position, num_charge, num_atomtype, num_atomid,
                               num_particleset, num_bond) == true ) {
    //      std::cout << "send particle " << num_position << std::endl;
    int i=0;
    for(size_t s=0;s<num_particleset;s++){
      TypeRange tr = rq_tr[set_index[s]];
      int offset = i-tr.begin;
      for(int si=tr.begin;si<tr.end;si++,i++){
        sendposition[i] = getpos(request,si);
        sendcharge[i] = getcharge(request,si);
        sendatomtype[i] = getatomtype(request,si);
        sendatomid[i] = getatomid(request,si);
      }
      tr.shift(offset);
      sendtyperange[s] = tr;
    }
    for(size_t s=0;s<num_particleset;s++){
      sendsetid[s] = rq_set[s];
    }
    for(size_t b=0;b<num_bond;b++){
      sendbond[b] = rq_bond[b];
    }
    stat = true;
  }else{
    stat =  false;
  }
  //  std::cout << " send_request_size " << send_request_size << " to " << target_rank << std::endl;
  return stat;
}

//! get forces from MPI recv buffer
/*!
  return continuous ForceArray
*/
bool
MPISendParticleReceiveForce::getreceive(ForceArray& forcearray)
{
  int result;
  result = MPI_Wait(&receive_requset,MPI_STATUS_IGNORE);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Wait at MPISendParticleReceiveForce::getreceive(ForceArray& forcearray) %d\n",result);
  }
  size_t receivesize = receivebufferheader->NumberOfForce;
  if(receivesize<=number_of_particle[stage]){
    forcearray.resize(number_of_particle[stage]);
    for(size_t i=0;i<size_t(receivesize);i++){
      forcearray[i] += ReceiveForceArray[i];
    }
    return true;
  }
  return false;
}
bool
MPISendParticleReceiveForce::getreceive(ForceArray& forcearray, int offset)
{
  int result;
  result = MPI_Wait(&receive_requset,MPI_STATUS_IGNORE);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Wait at MPISendParticleReceiveForce::getreceive(ForceArray& forcearray, int offset) %d\n",result);
  }
  size_t receivesize = receivebufferheader->NumberOfForce;
  if(receivesize<=number_of_particle[stage]){
    //    std::cout << " add force " << "[" << offset << "," << offset+number_of_particle[stage] << ")" << std::endl;
    for(size_t i=0;i<size_t(receivesize);i++){
      forcearray[i+offset] += ReceiveForceArray[i];
      // forcearray[i+offset].x += ReceiveForceArray[i].x;
      // forcearray[i+offset].y += ReceiveForceArray[i].y;
      // forcearray[i+offset].z += ReceiveForceArray[i].z;
    }
    return true;
  }
  std::cout << " number of received force " << receivesize << " != number of send particle[" << stage << "] " << number_of_particle[stage] << std::endl;
  return false;
}

//! get forces from MPI receive buffer
/*!
  return ForceArray with gap
*/
bool 
MPISendParticleReceiveForce::getreceivesubset(ForceArray& forcearray)
{
  int result;
  result = MPI_Wait(&receive_requset,MPI_STATUS_IGNORE);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Wait at MPISendParticleReceiveForce::getreceivesubset(ForceArray& forcearray) %d\n",result);
  }
  size_t receivesize = receivebufferheader->NumberOfForce;
  if(receivesize==number_of_particle[stage]){
    //      forcearray.resize(number_of_particle[stage]);
    size_t i=0;
    for(std::vector<TypeRange>::iterator it = sent_typerange.begin();
        it != sent_typerange.end(); ++it){
      for(int si=(*it).begin;si<(*it).end;si++,i++){
        forcearray[si] += ReceiveForceArray[i];
        // forcearray[si].x += ReceiveForceArray[i].x;
        // forcearray[si].y += ReceiveForceArray[i].y;
        // forcearray[si].z += ReceiveForceArray[i].z;
      }
    }
    //      if(unit_identifier==0)
    /*
      {
      std::cout << "receive force from " << receivebufferheader->UnitID << " size " << receivesize << std::endl;
      }
    */
    return true;
  }
  std::cout << " number of received force " << receivesize << " != number of send particle[" << stage << "] " << number_of_particle[stage] << std::endl;
  return false;
}

//! receive force with index, for sparse e.g bonded only 
// set forceindexarray, getreceive_indexed reuse it
bool 
MPISendParticleReceiveForce::getreceive_with_index(ForceArray& forcearray)
{
  int result;
  result = MPI_Wait(&receive_requset,MPI_STATUS_IGNORE);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Wait at MPISendParticleReceiveForce::getreceive_with_index(ForceArray& forcearray) %d\n",result);
  }
  size_t receivesize = receivebufferheader->NumberOfForce;
  //  printf("receivebufferheader NumberOfForce %ld\n",receivesize);
  if(receivesize<=number_of_particle[stage]){
    //      forcearray.resize(number_of_particle[stage]);
    set_receive_buffer_pointers_with_index(receivesize);
    forceindexarray.resize(receivesize);
    size_t i=0;
    size_t si=0;
    int index_offset=sent_typerange[si].begin;
    for(i=0;i<receivesize;i++){
      int index = ReceiveForceIndex[i];
      while(sent_typerange[si].end<=index+index_offset){
        index_offset += (sent_typerange[si+1].begin - sent_typerange[si].end);
        si++;
        if(si==sent_typerange.size()){
          if(sent_typerange[si].end<=index+index_offset){
            printf("too large force index %d >= %d\n",index+index_offset,
                   sent_typerange[si].end);
            break;
          }
        }
      }
      forceindexarray[i] = index+index_offset;
    }
    for(i=0;i<receivesize;i++){
      forcearray[forceindexarray[i]] += ReceiveForceArray[i];
      // forcearray[forceindexarray[i]].x += ReceiveForceArray[i].x;
      // forcearray[forceindexarray[i]].y += ReceiveForceArray[i].y;
      // forcearray[forceindexarray[i]].z += ReceiveForceArray[i].z;
    }
    //      if(unit_identifier==0)
    /*
      {
      std::cout << "receive force from " << receivebufferheader->UnitID << " size " << receivesize << std::endl;
      }
    */
    return true;
  }
  std::cout << " number of received force " << receivesize << " > number of send particle[" << stage << "] " << number_of_particle[stage] << std::endl;
  return false;
}

/*!
  unpack and add to forcearray
  @param [in,out] forcearray forces added received force
  @note indexes of received force had been listed forceindexarray
  @note any index listed in forceindexarray must not be duplicated for thread parallelization 
  @note call MPI_Wait before calling it
 */
bool 
MPISendParticleReceiveForce::getreceive_indexed_nowait(ForceArray& forcearray)
{
  size_t receivesize = receivebufferheader->NumberOfForce;
#ifdef K_SIMD
  double (*cforcearray)[3] = (double (*)[3])(&(forcearray[0].x));
  double (*cReceiveForceArray)[3] = (double (*)[3])(&(ReceiveForceArray[0].x));
  int *cforceindexarray = (int *)(&(forceindexarray[0]));
#endif
  if(receivesize==forceindexarray.size()){
    //      forcearray.resize(number_of_particle[stage]);
    size_t n = forceindexarray.size();
    size_t i=0;
#ifdef MPIPARALLEL_OPENMP
    /*
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    */
    //! TODO make slice size over 8?
#pragma loop parallel
#endif
#pragma loop norecurrence cforcearray,cReceiveForceArray,cforceindexarray
#pragma loop noalias cforcearray,cReceiveForceArray,cforceindexarray
    for(i=0;i<n;i++){
#ifdef K_SIMD
      int fi=cforceindexarray[i];
      cforcearray[fi][0] += cReceiveForceArray[i][0];
      cforcearray[fi][1] += cReceiveForceArray[i][1];
      cforcearray[fi][2] += cReceiveForceArray[i][2];
#else
      forcearray[forceindexarray[i]] += ReceiveForceArray[i];
      // forcearray[forceindexarray[i]].x += ReceiveForceArray[i].x;
      // forcearray[forceindexarray[i]].y += ReceiveForceArray[i].y;
      // forcearray[forceindexarray[i]].z += ReceiveForceArray[i].z;
#endif
    }
    //      if(unit_identifier==0)
    /*
      {
      std::cout << "receive force from " << receivebufferheader->UnitID << " size " << receivesize << std::endl;
      }
    */
    return true;
  }
  std::cout << " number of received force " << receivesize << " != forceindexarray.size() " << forceindexarray.size() << std::endl;
  return false;
}

/*!
  wait receiving force, unpack and add to forcearray
  @param [in,out] forcearray forces added received force
  @note indexes of received force had been listed forceindexarray
  @note any index listed in forceindexarray must not be duplicated for thread parallelization 
 */
bool 
MPISendParticleReceiveForce::getreceive_indexed(ForceArray& forcearray)
{
  int result;
  result = MPI_Wait(&receive_requset,MPI_STATUS_IGNORE);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Wait at MPISendParticleReceiveForce::getreceive_indexed(ForceArray& forcearray) %d\n",result);
  }
  size_t receivesize = receivebufferheader->NumberOfForce;
#ifdef K_SIMD
  double (*cforcearray)[3] = (double (*)[3])(&(forcearray[0].x));
  double (*cReceiveForceArray)[3] = (double (*)[3])(&(ReceiveForceArray[0].x));
  int *cforceindexarray = (int *)(&(forceindexarray[0]));
#endif
  if(receivesize==forceindexarray.size()){
    //      forcearray.resize(number_of_particle[stage]);
    size_t n = forceindexarray.size();
    size_t i=0;
#ifdef MPIPARALLEL_OPENMP
    /*
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    */
    //! TODO make slice size over 8?
#pragma loop parallel
#endif
#pragma loop norecurrence cforcearray,cReceiveForceArray,cforceindexarray
#pragma loop noalias cforcearray,cReceiveForceArray,cforceindexarray
    for(i=0;i<n;i++){
#ifdef K_SIMD
      int fi=cforceindexarray[i];
      cforcearray[fi][0] += cReceiveForceArray[i][0];
      cforcearray[fi][1] += cReceiveForceArray[i][1];
      cforcearray[fi][2] += cReceiveForceArray[i][2];
#else
      forcearray[forceindexarray[i]] += ReceiveForceArray[i];
      // forcearray[forceindexarray[i]].x += ReceiveForceArray[i].x;
      // forcearray[forceindexarray[i]].y += ReceiveForceArray[i].y;
      // forcearray[forceindexarray[i]].z += ReceiveForceArray[i].z;
#endif
    }
    //      if(unit_identifier==0)
    /*
      {
      std::cout << "receive force from " << receivebufferheader->UnitID << " size " << receivesize << std::endl;
      }
    */
    return true;
  }
  std::cout << " number of received force " << receivesize << " != forceindexarray.size() " << forceindexarray.size() << std::endl;
  return false;
}

//! request MPI asyncronous receive for Force
void 
MPISendParticleReceiveForce::prepare_receive()
{
  int result;
  result = MPI_Irecv(receive_buffer, receive_buffer_size, MPI_BYTE, 
                     target_rank, MPITagForce, mpi_comm_short, &receive_requset);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Irecv at MPISendParticleReceiveForce::prepare_receive() %d\n",result);
  }

}

//! request MPI send for Particle
/*!
  If defined PARTICLE_ISEND, use asyncronous send, and must be check by wait_complete_send() .
  If nodef PARTICLE_ISEND, use syncronous send, and MPIReceiveParticleSendForce::prepare_receive() must be called before this method.
  Recommend define PARTICLE_ISEND for multiple target node because receiver may be ready not in same order of send request.
  But no warranty multiple Isends are simultaneously waiting or start transmission not in requested order.
*/
void
MPISendParticleReceiveForce::send()
{
  int result;
#ifdef PARTICLE_ISEND
  result = MPI_Isend(send_buffer, send_request_size, MPI_BYTE, 
                     target_rank, MPITagParticle, 
                     mpi_comm_short, send_requestp);
#else
  result = MPI_Send(send_buffer, send_request_size, MPI_BYTE, 
                    target_rank, MPITagParticle, mpi_comm_short);
#endif
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Isend at MPISendParticleReceiveForce::send() %d\n",result);
  }
}

//! wait/check finish send
/*!
  If defined PARTICLE_ISEND, MUST be called.
  If nodef PARTICLE_ISEND, do nothing.
*/
void 
MPISendParticleReceiveForce::wait_complete_send(){
#ifdef PARTICLE_ISEND
  int result;
  result = MPI_Wait(send_requestp,MPI_STATUS_IGNORE);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Wait at MPISendParticleReceiveForce::wait_complete_send() %d\n",result);
  }
#endif
}


/*!
  Methods of MPIReceiveParticleSendForce
 */
//! construct, default size is hard corded
MPIReceiveParticleSendForce::MPIReceiveParticleSendForce()
  : size(buffersize),
    mpi_comm_short(MPI_COMM_WORLD),
    number_of_stage(1),
    stage(0)
{
  number_of_particle.resize(number_of_stage,0);
  //  std::cout << "MPIReceiveParticleSendForce::MPIReceiveParticleSendForce()" << std::endl;
  number_of_set.resize(number_of_stage,0);
  recv_pos_offset.resize(number_of_stage,0);
  send_buffer = (char *)NULL;
  recv_buffer = (char *)NULL;
}

MPIReceiveParticleSendForce::MPIReceiveParticleSendForce(int unitid, 
                                                         MPI_Comm short_comm)
  : size(buffersize),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm),
    number_of_stage(1),
    stage(0)
{
  //  std::cout << "MPIReceiveParticleSendForce::MPIReceiveParticleSendForce(int unitid)" << std::endl;
  number_of_particle.resize(number_of_stage,0);
  number_of_set.resize(number_of_stage,0);
  recv_pos_offset.resize(number_of_stage,0);
  send_buffer = (char *)NULL;
  recv_buffer = (char *)NULL;
}

//! construct, with specified size
MPIReceiveParticleSendForce::MPIReceiveParticleSendForce(int unitid, 
                                                         size_t sz,
                                                         MPI_Comm short_comm)
  : size(sz),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm),
    number_of_stage(1),
    stage(0)
{
  //  std::cout << "MPIReceiveParticleSendForce::MPIReceiveParticleSendForce(int unitid, size_t sz)" << std::endl;
  number_of_particle.resize(number_of_stage,0);
  number_of_set.resize(number_of_stage,0);
  recv_pos_offset.resize(number_of_stage,0);
  send_buffer = (char *)NULL;
  recv_buffer = (char *)NULL;
}

MPIReceiveParticleSendForce::MPIReceiveParticleSendForce(int unitid, 
                                                         size_t sz,
                                                         MPI_Comm short_comm,
                                                         int num_stage)
  : size(sz),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm),
    number_of_stage(num_stage),
    stage(0)
{
  //  std::cout << "MPIReceiveParticleSendForce::MPIReceiveParticleSendForce(int unitid, size_t sz)" << std::endl;
  number_of_particle.resize(number_of_stage,0);
  number_of_set.resize(number_of_stage,0);
  recv_pos_offset.resize(number_of_stage,0);
  send_buffer = (char *)NULL;
  recv_buffer = (char *)NULL;
}

MPIReceiveParticleSendForce::MPIReceiveParticleSendForce(const MPIReceiveParticleSendForce& rpsf)
  : size(rpsf.size),
    target_rank(rpsf.target_rank),
    target_id(rpsf.target_id),
    unit_identifier(rpsf.unit_identifier),
    mpi_comm_short(rpsf.mpi_comm_short),
    number_of_stage(rpsf.number_of_stage),
    number_of_particle(rpsf.number_of_particle),
    number_of_set(rpsf.number_of_set),
    stage(rpsf.stage),
    recv_pos_offset(rpsf.recv_pos_offset)
{
  //  std::cout << "MPIReceiveParticleSendForce::MPIReceiveParticleSendForce(const MPIReceiveParticleSendForce& rpsf)" << std::endl;
  if(rpsf.send_buffer!=(char *)NULL){
    make_buffer(rpsf.number_of_set[stage]);
  }else{
    send_buffer = (char *)NULL;
    recv_buffer = (char *)NULL;
  }
}

MPIReceiveParticleSendForce& MPIReceiveParticleSendForce::operator=(const MPIReceiveParticleSendForce &rpsf)
{
  //  std::cout << "MPIReceiveParticleSendForce::operator=(const MPIReceiveParticleSendForce &rpsf)" << std::endl;
  if(this != &rpsf){
    size = rpsf.size;
    target_rank = rpsf.target_rank;
    target_id = rpsf.target_id;
    unit_identifier = rpsf.unit_identifier;
    mpi_comm_short = rpsf.mpi_comm_short;
    number_of_stage = rpsf.number_of_stage;
    number_of_particle = rpsf.number_of_particle;
    number_of_set = rpsf.number_of_set;
    recv_pos_offset = rpsf.recv_pos_offset;
    stage = rpsf.stage;
    if(rpsf.send_buffer!=(char *)NULL){
      make_buffer(rpsf.number_of_set[stage]);
    }else{
      send_buffer = (char *)NULL;
      recv_buffer = (char *)NULL;
    }
  }
  return *this;
}

MPIReceiveParticleSendForce::~MPIReceiveParticleSendForce()
{
  //  std::cout << "MPIReceiveParticleSendForce::~MPIReceiveParticleSendForce() target_id " << target_id << std::endl;
  //  std::cout << "MPIReceiveParticleSendForce::~MPIReceiveParticleSendForce() target_id " << target_id ;
  /*
  if(send_buffer!=(char *)NULL){
    std::cout << " delete send_buffer " << size_t(send_buffer);
  }
  if(recv_buffer!=(char *)NULL){
    std::cout << " delete recv_buffer " << size_t(recv_buffer);
  }
  std::cout << std::endl;
  */
  //  std::cout << "MPIReceiveParticleSendForce::~MPIReceiveParticleSendForce() target_id " << target_id << "delete" << send_buffer << ","  << recv_buffer << std::endl;
  delete [] send_buffer;
  delete [] recv_buffer;
}

//! set positions of particle data pointers in recv MPI buffer
bool 
MPIReceiveParticleSendForce::set_receive_posbuffer_pointers()
{
  recv_requested_size = calc_position_buffer_size(recvposbufferheader[0]);
  if(recv_requested_size>recv_buffer_size){
    std::cout << "recv_requested_size>recv_buffer_size" << std::endl;
    return false;
  }
  void *currentpoint = &(recvposbufferheader[1]);
  recvposposition = static_cast<Position *>(currentpoint);
  number_of_particle[stage] = recvposbufferheader->NumberOfPosition;
  return true;
}

//! set particle data pointers in recv MPI buffer
bool 
MPIReceiveParticleSendForce::set_receive_buffer_pointers()
{
  recv_requested_size = calc_particle_buffer_size(recvbufferheader[0]);
  if(recv_requested_size>recv_buffer_size){
    std::cout << "recv_requested_size>recv_buffer_size" << std::endl;
    return false;
  }
  void *currentpoint = &(recvbufferheader[1]);
  recvposition = static_cast<Position *>(currentpoint);
  currentpoint = &(recvposition[recvbufferheader->NumberOfPosition]);
  recvcharge = static_cast<double *>(currentpoint);
  currentpoint = &(recvcharge[recvbufferheader->NumberOfCharge]);
  recvatomtype = static_cast<Atomtype *>(currentpoint);
  currentpoint = &(recvatomtype[recvbufferheader->NumberOfAtomType]);
  recvatomid = static_cast<AtomID *>(currentpoint);
  currentpoint = &(recvatomid[recvbufferheader->NumberOfAtomID]);
  recvsetid = static_cast<int *>(currentpoint);
  currentpoint = &(recvsetid[recvbufferheader->NumberOfParticleSet]);
  recvtyperange = static_cast<TypeRange *>(currentpoint);
  currentpoint = &(recvtyperange[recvbufferheader->NumberOfParticleSet]);
  recvbond = static_cast<CovalentBondInfo::BondList *>(currentpoint);
  number_of_particle[stage] = recvbufferheader->NumberOfPosition;
  return true;
}

//! set force data pointers in send MPI buffer
bool 
MPIReceiveParticleSendForce::set_send_buffer_pointers(int number_of_force)
{
  send_expected_size = calc_force_buffer_size(number_of_force);
  if(send_expected_size>send_buffer_size){
    std::cout << "set_send_buffer_pointers from " << unit_identifier << " to " << target_id;
    std::cout << " send_expected_size " << send_expected_size << " ( " << number_of_force << " force) > send_buffer_size " << send_buffer_size << " (resereved " << size << " forece)" << std::endl;
    return false;
  }
  void *currentpoint = send_buffer;
  sendbufferheader = static_cast<ForceBufferHeader *>(currentpoint);
  currentpoint = &(sendbufferheader[1]);
  SendForceArray = static_cast<Force *>(currentpoint);
  return true;
}
bool 
MPIReceiveParticleSendForce::set_send_buffer_pointers_with_index(int number_of_force)
{
  send_expected_size = calc_force_buffer_size_with_index(number_of_force);
  if(send_expected_size>send_buffer_size){
    std::cout << "set_send_buffer_pointers_with_index from " << unit_identifier << " to " << target_id;
    std::cout << " send_expected_size " << send_expected_size << ">send_buffer_size " << send_buffer_size << std::endl;
    return false;
  }
  void *currentpoint = send_buffer;
  sendbufferheader = static_cast<ForceBufferHeader *>(currentpoint);
  currentpoint = &(sendbufferheader[1]);
  SendForceArray = static_cast<Force *>(currentpoint);
  currentpoint = &(SendForceArray[number_of_force]);
  SendForceIndex = static_cast<int *>(currentpoint);
  return true;
}

//! reserve receive and send MPI buffer
void 
MPIReceiveParticleSendForce::make_buffer(int num_set)
{
  number_of_set[stage] = num_set;
  number_of_particle[stage] = size;
  recv_buffer_size = calc_particle_buffer_size(size, size, size, size, num_set, 1);
  recv_buffer = new char[recv_buffer_size];
  void *tp = recv_buffer;
  recvposbufferheader = static_cast<PositionBufferHeader *>(tp);
  recvposbufferheader->UnitID = 0;
  recvposbufferheader->NumberOfPosition = size;
  set_receive_posbuffer_pointers();
  recvbufferheader = static_cast<ParticleBufferHeader *>(tp);
  recvbufferheader->UnitID = 0;
  recvbufferheader->NumberOfPosition = size;
  recvbufferheader->NumberOfCharge = size;
  recvbufferheader->NumberOfAtomType = size;
  recvbufferheader->NumberOfAtomID = size;
  recvbufferheader->NumberOfParticleSet = num_set;
  recvbufferheader->NumberOfBondlist = 1;
  set_receive_buffer_pointers();
    
  send_buffer_size = calc_force_buffer_size(size);
  send_buffer = new char[send_buffer_size];
  set_send_buffer_pointers(size);
  //    std::cout << "MPIReceiveParticleSendForce::make_buffer particle " << number_of_particle[stage] << " buffer size send " << send_buffer_size << " receive " << recv_buffer_size << std::endl; 
  /*
  std::cout << "MPIReceiveParticleSendForce::make_buffer ";
  std::cout << " send_buffer " << size_t(send_buffer);
  std::cout << " recv_buffer " << size_t(recv_buffer);
  std::cout << std::endl;
  */
}
bool 
MPIReceiveParticleSendForce::wait_complete_recv()
{
  int ret;
  //    std::cout << "getreceive" << std::endl;
  ret = MPI_Wait(recv_requestp,MPI_STATUS_IGNORE);
  //    std::cout << "MPI_Wait" << std::endl;
  if(ret!=MPI_SUCCESS){
    std::cout << "MPIReceiveParticleSendForce::wait_complete_recv MPI_Wait return not MPI_SUCCESS " << ret << std::endl;
    return false;
  }
  return true;
}

//! get receive Particles from recv MPI buffer
template <class GPA>
bool 
MPIReceiveParticleSendForce::getreceivepos_nowait(GPA& request)
{
  if( set_receive_posbuffer_pointers() == true ) {
    //    std::cout << " recv_requested_size " << recv_requested_size << std::endl;
    //     std::cout << "recv particle " << recvbufferheader->NumberOfPosition  << std::endl;
    size_t recvnum = recvposbufferheader->NumberOfPosition;
    size_t remain = request.size()-recv_pos_offset[stage];
    if(recvnum>remain){
      std::cout << "recieved Number of Position " << recvnum << " is larger than remaining particle array " << remain;
      std::cout << " discard last received " << recvnum-remain << " position";
      std::cout << " at stage " << stage << std::endl;
      recvnum = remain;
    }
    for(size_t i=0;i<recvnum;i++){
      //      request[recv_pos_offset[stage]+i].position = recvposposition[i];
      getpos(request,recv_pos_offset[stage]+i) = recvposposition[i];
      // getpos(request,recv_pos_offset[stage]+i).x = recvposposition[i].x;
      // getpos(request,recv_pos_offset[stage]+i).y = recvposposition[i].y;
      // getpos(request,recv_pos_offset[stage]+i).z = recvposposition[i].z;
      /*
      request[recv_pos_offset[stage]+i].force.x = 0.0;
      request[recv_pos_offset[stage]+i].force.y = 0.0;
      request[recv_pos_offset[stage]+i].force.z = 0.0;
      */
    }
    return true;
  }
  return false;
}

template
bool
MPIReceiveParticleSendForce::getreceivepos_nowait(ParticleArray& request);
template
bool
MPIReceiveParticleSendForce::getreceivepos_nowait(GhostParticleArray& request);

template <class GPA>
bool 
MPIReceiveParticleSendForce::getreceivepos(GPA& request)
{
  int ret;
  //    std::cout << "getreceive" << std::endl;
  ret = MPI_Wait(recv_requestp,MPI_STATUS_IGNORE);
  //    std::cout << "MPI_Wait" << std::endl;
  if(ret!=MPI_SUCCESS){
    std::cout << "MPIReceiveParticleSendForce::getreceivepos MPI_Wait return not MPI_SUCCESS " << ret << std::endl;
    return false;
  }
  if( set_receive_posbuffer_pointers() == true ) {
    //    std::cout << " recv_requested_size " << recv_requested_size << std::endl;
    //     std::cout << "recv particle " << recvbufferheader->NumberOfPosition  << std::endl;
    size_t recvnum = recvposbufferheader->NumberOfPosition;
    size_t remain = request.size()-recv_pos_offset[stage];
    if(recvnum>remain){
      std::cout << "recieved Number of Position " << recvnum << " is larger than remaining particle array " << remain;
      std::cout << " discard last received " << recvnum-remain << " position";
      std::cout << " at stage " << stage << std::endl;
      recvnum = remain;
    }
    for(size_t i=0;i<recvnum;i++){
      //      request[recv_pos_offset[stage]+i].position = recvposposition[i];
      getpos(request,recv_pos_offset[stage]+i) = recvposposition[i];
      // getpos(request,recv_pos_offset[stage]+i).x = recvposposition[i].x;
      // getpos(request,recv_pos_offset[stage]+i).y = recvposposition[i].y;
      // getpos(request,recv_pos_offset[stage]+i).z = recvposposition[i].z;
      /*
     request[recv_pos_offset[stage]+i].force.x = 0.0;
      request[recv_pos_offset[stage]+i].force.y = 0.0;
      request[recv_pos_offset[stage]+i].force.z = 0.0;
      */
    }
    return true;
  }
  return false;
}

template
bool 
MPIReceiveParticleSendForce::getreceivepos(ParticleArray& request);
template
bool 
MPIReceiveParticleSendForce::getreceivepos(GhostParticleArray& request);

bool 
MPIReceiveParticleSendForce::getreceive_number(int &pnum,
                                               int &setnum,
                                               int &bondnum)
{
  if( set_receive_buffer_pointers() == true ) {
    pnum = recvbufferheader->NumberOfPosition;
    if(recvbufferheader->NumberOfCharge!=pnum){
      printf("NumberOfCharge != NumberOfPosition\n");
    }
    if(recvbufferheader->NumberOfAtomType!=pnum){
      printf("NumberOfAtomType != NumberOfPosition\n");
    }
    if(recvbufferheader->NumberOfAtomID!=pnum){
      printf("NumberOfAtomID != NumberOfPosition\n");
    }
    setnum = recvbufferheader->NumberOfParticleSet;
    bondnum = recvbufferheader->NumberOfBondlist;
    return true;
  }
  return false;
}
template<class PA>
bool 
MPIReceiveParticleSendForce::getreceive_offset(PA& request, 
                                               ParticleRange& range,
                                               std::vector<TypeRange>& rq_tr,
                                               std::vector<CovalentBondInfo::BondList>& rq_bond,
                                               std::vector<int>& rq_set,
                                               int setoffset, int bondoffset)
{
  {
    recv_pos_offset[stage] = range.begin;
#ifdef SEPARATE_POSC
    for(int i=0;i<recvbufferheader->NumberOfPosition;i++){
      getpos(request,range.begin+i) = recvposition[i];
      // getpos(request,range.begin+i).x = recvposition[i].x;
      // getpos(request,range.begin+i).y = recvposition[i].y;
      // getpos(request,range.begin+i).z = recvposition[i].z;
    }
    for(int i=0;i<recvbufferheader->NumberOfCharge;i++){
      getcharge(request,range.begin+i) = recvcharge[i];
    }
#else  // new particle hold (x,y,z,charge)*n
    if(recvbufferheader->NumberOfPosition!=recvbufferheader->NumberOfCharge){
      printf("NumberOfPosition != NumberOfCharge\n");
      for(int i=0;i<recvbufferheader->NumberOfPosition;i++){
	getpos(request,range.begin+i) = recvposition[i];
      }
      for(int i=0;i<recvbufferheader->NumberOfCharge;i++){
	getcharge(request,range.begin+i) = recvcharge[i];
      }
    }else{
      for(int i=0;i<recvbufferheader->NumberOfPosition;i++){
	getpos(request,range.begin+i) = recvposition[i];
	getcharge(request,range.begin+i) = recvcharge[i];
      }
    }
#endif
    for(int i=0;i<recvbufferheader->NumberOfAtomType;i++){
      getatomtype(request,range.begin+i) = recvatomtype[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfAtomID;i++){
      getatomid(request,range.begin+i) = recvatomid[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfParticleSet;i++){
      rq_set[setoffset+i] = recvsetid[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfParticleSet;i++){
      TypeRange tr = recvtyperange[i];
      tr.shift(range.begin);
      rq_tr[setoffset+i] = tr;
    }
    for(int i=0;i<recvbufferheader->NumberOfBondlist;i++){
      rq_bond[bondoffset+i] = recvbond[i];
    }
    return true;
  }
}
//! get receive Particles from recv MPI buffer
bool 
MPIReceiveParticleSendForce::getreceive_nowait(ParticleArray& request, 
                                               ParticleRange& range,
                                               std::vector<TypeRange>& rq_tr,
                                               std::vector<CovalentBondInfo::BondList>& rq_bond,
                                               std::vector<int>& rq_set,
                                               std::map<int,int>& recvsetid_to_index)
{
  if( set_receive_buffer_pointers() == true ) {
    //    std::cout << " recv_requested_size " << recv_requested_size << std::endl;
    //     std::cout << "recv particle " << recvbufferheader->NumberOfPosition  << std::endl;
    recv_pos_offset[stage] = range.begin;
    request.resize(request.size()+recvbufferheader->NumberOfPosition);
    for(int i=0;i<recvbufferheader->NumberOfPosition;i++){
      request[range.begin+i].position = recvposition[i];
      // request[range.begin+i].position.x = recvposition[i].x;
      // request[range.begin+i].position.y = recvposition[i].y;
      // request[range.begin+i].position.z = recvposition[i].z;
    }
    range.end = range.begin+recvbufferheader->NumberOfPosition;
    for(int i=0;i<recvbufferheader->NumberOfCharge;i++){
      request[range.begin+i].charge = recvcharge[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfAtomType;i++){
      request[range.begin+i].atomtype = recvatomtype[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfAtomID;i++){
      request[range.begin+i].atomid = recvatomid[i];
    }
    int setindex=rq_set.size();
    for(int i=0;i<recvbufferheader->NumberOfParticleSet;i++){
      rq_set.push_back(recvsetid[i]);
      recvsetid_to_index.insert(std::pair<int,int>(recvsetid[i],setindex));
      setindex++;
    }
    for(int i=0;i<recvbufferheader->NumberOfParticleSet;i++){
      TypeRange tr = recvtyperange[i];
      tr.shift(range.begin);
      rq_tr.push_back(tr);
    }
    for(int i=0;i<recvbufferheader->NumberOfBondlist;i++){
      rq_bond.push_back(recvbond[i]);
    }
    return true;
  }
  return false;
}
template<class PA>
bool 
MPIReceiveParticleSendForce::getreceive(PA& request, 
                                        ParticleRange& range,
                                        std::vector<TypeRange>& rq_tr,
                                        std::vector<CovalentBondInfo::BondList>& rq_bond,
                                        std::vector<int>& rq_set,
                                        std::map<int,int>& recvsetid_to_index)
{
  int ret;
  //    std::cout << "getreceive" << std::endl;
  ret = MPI_Wait(recv_requestp,MPI_STATUS_IGNORE);
  //    std::cout << "MPI_Wait" << std::endl;
  if(ret!=MPI_SUCCESS){
    std::cout << "MPIReceiveParticleSendForce::getreceive MPI_Wait return not MPI_SUCCESS " << ret << std::endl;
    return false;
  }
  if( set_receive_buffer_pointers() == true ) {
    //    std::cout << " recv_requested_size " << recv_requested_size << std::endl;
    //     std::cout << "recv particle " << recvbufferheader->NumberOfPosition  << std::endl;
    recv_pos_offset[stage] = range.begin;
    request.resize(request.size()+recvbufferheader->NumberOfPosition);
    for(int i=0;i<recvbufferheader->NumberOfPosition;i++){
      getpos(request,range.begin+i) = recvposition[i];
      // getpos(request,range.begin+i).x = recvposition[i].x;
      // getpos(request,range.begin+i).y = recvposition[i].y;
      // getpos(request,range.begin+i).z = recvposition[i].z;
    }
    range.end = range.begin+recvbufferheader->NumberOfPosition;
    for(int i=0;i<recvbufferheader->NumberOfCharge;i++){
      getcharge(request,range.begin+i) = recvcharge[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfAtomType;i++){
      getatomtype(request,range.begin+i) = recvatomtype[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfAtomID;i++){
      getatomid(request,range.begin+i) = recvatomid[i];
    }
    int setindex=rq_set.size();
    for(int i=0;i<recvbufferheader->NumberOfParticleSet;i++){
      rq_set.push_back(recvsetid[i]);
      recvsetid_to_index.insert(std::pair<int,int>(recvsetid[i],setindex));
      setindex++;
    }
    for(int i=0;i<recvbufferheader->NumberOfParticleSet;i++){
      TypeRange tr = recvtyperange[i];
      tr.shift(range.begin);
      rq_tr.push_back(tr);
    }
    for(int i=0;i<recvbufferheader->NumberOfBondlist;i++){
      rq_bond.push_back(recvbond[i]);
    }
    return true;
  }
  return false;
}
template
bool 
MPIReceiveParticleSendForce::getreceive(GhostParticleArray& request, 
                                        ParticleRange& range,
                                        std::vector<TypeRange>& rq_tr,
                                        std::vector<CovalentBondInfo::BondList>& rq_bond,
                                        std::vector<int>& rq_set,
                                        std::map<int,int>& recvsetid_to_index);

//! set forces to MPI send buffer
bool 
MPIReceiveParticleSendForce::setsend(const ForceArray& forcearray, 
                                     const ParticleRange range)
{
  size_t sendsize = range.end-range.begin;
  if( set_send_buffer_pointers(sendsize) == true ){
    sendbufferheader->UnitID = unit_identifier;
    sendbufferheader->NumberOfForce = sendsize;
    for(size_t i=0;i<sendsize;i++){
      SendForceArray[i] = forcearray[i+range.begin];
    }
    return true;
  }
  printf("MPIReceiveParticleSendForce::setsend fail\n");
  return false;
}

//! set forces to MPI send buffer
bool 
MPIReceiveParticleSendForce::setsend_with_index(const ForceArray& forcearray, 
                                                const ParticleRange range,
                                                const std::vector<int>& forceindexset)
{
  indexlast = forceindexset.size();
  if(indexlast>0){
    //    indexbegin=0;
    //    while(forceindexset[indexbegin]<range.begin)indexbegin++;
    for(indexbegin=0;indexbegin<indexlast;indexbegin++){
      if(forceindexset[indexbegin]>=range.begin) break;
    }
    if(indexbegin<indexlast){
      indexend = indexbegin;
      while((forceindexset[indexend]<range.end)&&(indexend<indexlast))indexend++;
    }else{
      indexbegin=0;
      indexend=0;
    }
  }else{
    indexbegin=0;
    indexend=0;
  }
  size_t sendsize = indexend-indexbegin;
//  printf("sendbufferheader NumberOfForce %d,%d  %ld\n",indexbegin, indexend, sendsize);
  if( set_send_buffer_pointers_with_index(sendsize) == true ){
    sendbufferheader->UnitID = unit_identifier;
    sendbufferheader->NumberOfForce = sendsize;
    for(size_t i=0;i<sendsize;i++){
      SendForceArray[i] = forcearray[forceindexset[i+indexbegin]];
    }
    for(size_t i=0;i<sendsize;i++){
      SendForceIndex[i]=forceindexset[i+indexbegin]-range.begin;
    }
    return true;
  }
  return false;
}
bool 
MPIReceiveParticleSendForce::setsend_indexed(const ForceArray& forcearray,
                                             const ParticleRange range, 
                                             const std::vector<int>& forceindexset)
{
  size_t sendsize = indexend-indexbegin;
  if( set_send_buffer_pointers(sendsize) == true ){
    sendbufferheader->UnitID = unit_identifier;
    sendbufferheader->NumberOfForce = sendsize;
    for(size_t i=0;i<sendsize;i++){
      SendForceArray[i] = forcearray[forceindexset[i+indexbegin]];
    }
    return true;
  }
  return false;
}

//! request MPI asyncronous receive for Particle
void 
MPIReceiveParticleSendForce::prepare_receive()
{
  int result;
  //  std::cout << " recv_buffer_size " << recv_buffer_size << " from " << target_rank << std::endl;
  result = MPI_Irecv(recv_buffer, recv_buffer_size, MPI_BYTE, 
            target_rank, MPITagParticle, mpi_comm_short, recv_requestp);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Irecv at MPIReceiveParticleSendForce::prepare_receive() %d at %d from %d:%d\n",result,unit_identifier,target_rank,target_id);
    fflush(stdout);
  }

  //      if(unit_identifier==0)
  /*
    {
    std::cout << "recv particle from " << target_rank << " " << unit_identifier << " recv_buffer_size " << recv_buffer_size << std::endl;
    }
  */
}

//! request MPI send for Force
void 
MPIReceiveParticleSendForce::send()
{
  int result;
  result = MPI_Isend(send_buffer, send_expected_size, MPI_BYTE, 
                     target_rank, MPITagForce, mpi_comm_short, &send_request);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Isend at MPIReceiveParticleSendForce::send() %d\n",result);
  }
  //      if(unit_identifier==0)
  /*
    {
    std::cout << "send force to " << target_rank << " " << target_id << " size " << send_expected_size << std::endl;
    }
  */
}

//! wait/check finish send
void 
MPIReceiveParticleSendForce::wait_complete_send()
{
  int result;
  result = MPI_Wait(&send_request,MPI_STATUS_IGNORE);
  if(result!=MPI_SUCCESS){
    printf("Fail MPI_Wait at MPIReceiveParticleSendForce::wait_complete_send() %d\n",result);
  }

}




/*!
  Methods of MPIMoveOut
*/
//! construct, default size is hard corded
MPIMoveOut::MPIMoveOut()
  : size(buffersize),
    target_rank(),
    target_id(),
    mpi_comm_short(MPI_COMM_WORLD)
{
  //  std::cout << "MPIMoveOut::MPIMoveOut() " << std::endl;
  send_buffer = (char *)NULL;
}
MPIMoveOut::MPIMoveOut(int unitid, MPI_Comm short_comm)
  : size(buffersize),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm)
{
  //  std::cout << "MPIMoveOut::MPIMoveOut(int unitid) " << std::endl;
  send_buffer = (char *)NULL;
}
//! construct, with specified size
MPIMoveOut::MPIMoveOut(int unitid, size_t sz, MPI_Comm short_comm)
  : size(sz),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm)
{
  //  std::cout << "MPIMoveOut::MPIMoveOut(int unitid, size_t sz) " << std::endl;
  send_buffer = (char *)NULL;
}
MPIMoveOut::MPIMoveOut(const MPIMoveOut& mo)
   : size(mo.size),
     target_rank(mo.target_rank),
     target_id(mo.target_id),
     unit_identifier(mo.unit_identifier),
     mpi_comm_short(mo.mpi_comm_short)
{
  //  std::cout << "MPIMoveOut::MPIMoveOut(const MPIMoveOut& mo) " << std::endl;
  if(mo.send_buffer!=(char *)NULL){
    make_buffer();
  }else{
    send_buffer = (char *)NULL;
  }
}
MPIMoveOut& MPIMoveOut::operator=(const MPIMoveOut& mo)
{
 if(this != &mo){
   //   std::cout << "MPIMoveOut operator= " << std::endl;
    size = mo.size;
    target_rank = mo.target_rank;
    target_id = mo.target_id;
    unit_identifier = mo.unit_identifier;
    mpi_comm_short = mo.mpi_comm_short;
    delete [] send_buffer;
    if(mo.send_buffer!=(char *)NULL){
      make_buffer();
      //      std::cout << "MPIMoveOut operator= make_buffer" << std::endl;
    }else{
      send_buffer = (char *)NULL;
    }
  }
  return *this;
}
MPIMoveOut::~MPIMoveOut()
{
  delete [] send_buffer;
}
//! set particle data (without bond) pointers in send MPI buffer
/*!
  @param[in]  num_particle number of particle
 */
bool 
MPIMoveOut::set_send_buffer_pointers(int num_particle)
{
  send_request_size = calc_move_buffer_size(num_particle);
  //    std::cout << " make move out buffer for particle " << num_particle << " " << send_request_size << " byte "  << std::endl;

  if( send_request_size > send_buffer_size ){
    std::cout << "particle " << num_particle << "  send_request_size > send_buffer_size for " << size << " particle" <<  std::endl;
    return false;
  }
  void *currentpoint = send_buffer;
  sendbufferheader = static_cast<MoveBufferHeader *>(currentpoint);
  sendbufferheader->UnitID = unit_identifier;
  sendbufferheader->NumberOfParticle = num_particle;
  sendbufferheader->NumberOfBondlist = 0;
  sendbufferheader->SizeOfBondlist = 0;
  currentpoint = &(sendbufferheader[1]);
  sendparticle = static_cast<Particle *>(currentpoint);
  currentpoint = &(sendparticle[number_of_particle]);
  sendtype = static_cast<PotentialModel *>(currentpoint);
  currentpoint = &(sendtype[number_of_particle]);
  sendcellid = static_cast<int *>(currentpoint);
  if(DebugLog::verbose>2) {
    std::cout << "UID " << unit_identifier << " move out particle for " << target_id << " " << num_particle << " particles" << std::endl;
  }
  return true;
}
//! set particle data (with bond) pointers in send MPI buffer
/*!
  @param[in] num_particle number of particle
  @param[in] bondlistarray array of bondlist
  @todo fix size calculation of bondlist
*/
bool 
MPIMoveOut::set_send_buffer_pointers(int num_particle,
                                     std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  send_request_size = calc_move_buffer_size(num_particle,bondlistarray);
  //    std::cout << " make move out buffer for particle " << number_of_particle << " " << send_request_size << " byte "  << std::endl;

  if( send_request_size > send_buffer_size ){
    std::cout << "particle " << num_particle << "  send_request_size > send_buffer_size for " << size << " particle and " << bondlistarray.size() << " bondlist" <<  std::endl;
    return false;
  }
  void *currentpoint = send_buffer;
  sendbufferheader = static_cast<MoveBufferHeader *>(currentpoint);
  sendbufferheader->UnitID = unit_identifier;
  sendbufferheader->NumberOfParticle = num_particle;
  sendbufferheader->NumberOfBondlist = bondlistarray.size();
  currentpoint = &(sendbufferheader[1]);
  sendparticle = static_cast<Particle *>(currentpoint);
  currentpoint = &(sendparticle[number_of_particle]);
  sendtype = static_cast<PotentialModel *>(currentpoint);
  currentpoint = &(sendtype[number_of_particle]);
  sendcellid = static_cast<int *>(currentpoint);
  currentpoint = &(sendcellid[number_of_particle]);
  int sob=0;
#ifdef TRANSFER_RAW_BONDLIST
  sendbondlistarray.resize(bondlistarray.size());
  for(size_t i=0;i<bondlistarray.size();i++){
    sendbondlistarray[i] = static_cast<CovalentBondInfo::BondList *>(currentpoint);
    currentpoint = &(sendbondlistarray[i][1]);
    sob += sendbondlistarray[i]->size();
  }
#else
  sendbondpackarray.resize(bondlistarray.size());
  for(size_t i=0;i<bondlistarray.size();i++){
    sendbondpackarray[i] = static_cast<int *>(currentpoint);
    currentpoint = &(sendbondpackarray[i][bondlistarray[i].size_of_packed_int()/sizeof(int)]);
    sob += bondlistarray[i].size_of_packed_int();
  }
#endif
  sendbufferheader->SizeOfBondlist = sob;

  if(DebugLog::verbose>2) {
    std::cout << "UID " << unit_identifier << " move out particle for " << target_id << " " << num_particle << " particles" << std::endl;
    std::cout << "UID " << unit_identifier << " move out bond for " << target_id << " " << bondlistarray.size() << " bonds" << std::endl;
  }

  return true;
}
//! reserve send MPI buffer
void 
MPIMoveOut::make_buffer()
{
  number_of_particle = size;
  CovalentBondInfo::BondList tmpbondlist(2,3,3,1);
  //  CovalentBondInfo::BondList tmpbondlist(8,12,12,8);
  std::vector<CovalentBondInfo::BondList> tmpbondlistarray(size,
                                                           tmpbondlist);
  
  //   std::cout << "reserve for move out " << number_of_particle << std::endl;
  send_buffer_size = calc_move_buffer_size(size*2,tmpbondlistarray);
  if(unit_identifier==0) {
    if(DebugLog::verbose>2){
      std::cout << "reserve for " << size*2 << " particle " << std::endl;
      std::cout << "sizeof(tmpbondlist) " << sizeof(tmpbondlist) << std::endl;
      std::cout << "tmpbondlist.size() " << tmpbondlist.size() << std::endl;
      std::cout << "sizeof(tmpbondlistarray) " << sizeof(tmpbondlistarray) << std::endl;
      std::cout << "tmpbondlistarray.size() " << tmpbondlistarray.size() << std::endl;
      std::cout << "tmpbondlistarray " << tmpbondlist.size()*tmpbondlistarray.size() << std::endl;
      std::cout << "send_buffer_size " << send_buffer_size << std::endl;
    }
  }
  send_buffer = new char[send_buffer_size]();
  tmpbondlistarray.resize(size);
  set_send_buffer_pointers(size,tmpbondlistarray);
}
//! set Move out Particles to MPI send buffer
bool 
MPIMoveOut::setmoveoutparticle(ParticleArray& moveout,
                               std::vector<PotentialModel>& moveouttype,
                               std::vector<int>& outcellid)
{
  number_of_particle = moveout.size();
  if( set_send_buffer_pointers(number_of_particle) == true ) {
    for(size_t i=0;i<moveout.size();i++){
      sendparticle[i] = moveout[i];
    }
    for(size_t i=0;i<moveouttype.size();i++){
      sendtype[i] = moveouttype[i];
    }
    for(size_t i=0;i<outcellid.size();i++){
      sendcellid[i] = outcellid[i];
    }    
    return true;
  }
  std::cout << "setmoveoutparticle fail" << std::endl;
  return false;
}
bool 
MPIMoveOut::setmoveoutparticle(ParticleArray& moveout,
                               std::vector<PotentialModel>& moveouttype,
                               std::vector<int>& outcellid,
                               std::vector<CovalentBondInfo::BondList>& outbondlistarray)
{
  number_of_particle = moveout.size();
  if( set_send_buffer_pointers(number_of_particle,outbondlistarray) == true ) {
    for(size_t i=0;i<moveout.size();i++){
      sendparticle[i] = moveout[i];
    }
    for(size_t i=0;i<moveouttype.size();i++){
      sendtype[i] = moveouttype[i];
    }
    for(size_t i=0;i<outcellid.size();i++){
      sendcellid[i] = outcellid[i];
    }
    if((DebugLog::verbose>2)
       ||((DebugLog::verbose==2)&&(outbondlistarray.size()>0))){
      std::cout << "outbondlistarray " << outbondlistarray.size() << std::endl;
    }
    for(size_t i=0;i<outbondlistarray.size();i++){
#ifdef TRANSFER_RAW_BONDLIST
      sendbondlistarray[i][0] = outbondlistarray[i];
      outbondlistarray[i].print();
#else
      outbondlistarray[i].pack_int_array(sendbondpackarray[i]);
#endif
    }
    return true;
  }
  std::cout << "setmoveoutparticle fail" << std::endl;
  return false;
}
//! request MPI send for Move out Particle
void 
MPIMoveOut::send()
{
  int ret = MPI_Isend(send_buffer, send_request_size, MPI_BYTE, 
                      target_rank, MPITagMove, mpi_comm_short, &send_request);
  if(ret!=MPI_SUCCESS){
    std::cout << " MPI_Isend error at MPIMoveOut::send() " << ret << std::endl;
  }
}
//! wait/check finish send
void 
MPIMoveOut::wait_complete_send()
{
  MPI_Status mpistatus;
  int ret = MPI_Wait(&send_request,&mpistatus);
  if(ret!=MPI_SUCCESS){
    std::cout << " MPI_Wait error at MPIMoveOut::wait_complete_send() " << ret << std::endl;
  }
}


/*!
  Methods of MPIMoveIn
*/
MPIMoveIn::MPIMoveIn()
  : size(buffersize),
    target_rank(),
    target_id(),
    mpi_comm_short(MPI_COMM_WORLD)
{
  recv_buffer = (char *)NULL;
}
MPIMoveIn::MPIMoveIn(int unitid, MPI_Comm short_comm)
  : size(buffersize),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm)
{
  recv_buffer = (char *)NULL;
}
//! construct, with specified size
MPIMoveIn::MPIMoveIn(int unitid, size_t sz, MPI_Comm short_comm)
  : size(sz),
    target_rank(),
    target_id(),
    unit_identifier(unitid),
    mpi_comm_short(short_comm)
{
  recv_buffer = (char *)NULL;
}
MPIMoveIn::MPIMoveIn(const MPIMoveIn& mi)
  : size(mi.size),
    target_rank(mi.target_rank),
    target_id(mi.target_id),
    unit_identifier(mi.unit_identifier),
    mpi_comm_short(mi.mpi_comm_short)
{
  if(mi.recv_buffer!=(char *)NULL){
    make_buffer();
  }else{
    recv_buffer = (char *)NULL;
  }
}
MPIMoveIn& MPIMoveIn::operator=(const MPIMoveIn& mi)
{
  if(this != &mi){
    size = mi.size;
    target_rank = mi.target_rank;
    target_id = mi.target_id;
    unit_identifier = mi.unit_identifier;
    mpi_comm_short = mi.mpi_comm_short;
    delete [] recv_buffer;
    if(mi.recv_buffer!=(char *)NULL){
      make_buffer();
    }else{
      recv_buffer = (char *)NULL;
    }
  }
  return *this;
}
MPIMoveIn::~MPIMoveIn()
{
  delete []recv_buffer;
}
//! set particle data pointers in recv MPI buffer
/*!
  @note buffer size check is unsafe, because data size of bondlistarray is not
considerd at calc_move_buffer_size(recvbufferheader[0]).
 */
bool 
MPIMoveIn::set_receive_buffer_pointers()
{
  recv_requested_size = calc_move_buffer_size(recvbufferheader[0]);
  if(recv_requested_size>recv_buffer_size){
    std::cout << "recv_requested_size>recv_buffer_size" << std::endl;
    return false;
  }
  void *currentpoint = &(recvbufferheader[1]);
  recvparticle = static_cast<Particle *>(currentpoint);
  currentpoint = &(recvparticle[recvbufferheader->NumberOfParticle]);
  recvtype = static_cast<PotentialModel *>(currentpoint);
  currentpoint = &(recvtype[recvbufferheader->NumberOfParticle]);
  recvcellid = static_cast<int *>(currentpoint);
  currentpoint = &(recvcellid[recvbufferheader->NumberOfParticle]);
#ifdef TRANSFER_RAW_BONDLIST
  recvbondlistarray.resize(recvbufferheader->NumberOfBondlist);
  for(size_t i=0;i<recvbufferheader->NumberOfBondlist;i++){
    recvbondlistarray[i] = static_cast<CovalentBondInfo::BondList *>(currentpoint);
    currentpoint = &(recvbondlistarray[i][1]);
  }
#else
  recvbondpackarray.resize(recvbufferheader->NumberOfBondlist);
  for(int i=0;i<recvbufferheader->NumberOfBondlist;i++){
    recvbondpackarray[i] = static_cast<int *>(currentpoint);
    currentpoint = &(recvbondpackarray[i][recvbondpackarray[i][0]]);
  }
#endif
  return true;
}
//! reserve receive MPI buffer
/*!
  @todo enlarge buffer size for bondlistarray
*/
void 
MPIMoveIn::make_buffer()
{
  number_of_particle = size;

  int bondsize;
  {
    CovalentBondInfo::BondList tmpbondlist(2,3,3,1);
#ifdef TRANSFER_RAW_BONDLIST
    bondsize = tmpbondlist.size()*size;
#else
    bondsize = tmpbondlist.size_of_packed_int()*size;
#endif
  }
  recv_buffer_size = calc_move_buffer_size(size*2)+bondsize;
  if((unit_identifier==0)&&(DebugLog::verbose>1)) {
    std::cout << "MPIMoveIn.recv_buffer_size " << recv_buffer_size << std::endl;
  }
  recv_buffer = new char[recv_buffer_size]();
  // recv_buffer is clean?
  void *tp = recv_buffer;
  recvbufferheader = static_cast<MoveBufferHeader *>(tp);
  recvbufferheader->UnitID = 0;
  recvbufferheader->NumberOfParticle = size;
  recvbufferheader->NumberOfBondlist = size;
  recvbufferheader->SizeOfBondlist = bondsize;
  set_receive_buffer_pointers();
}
//! get receive Move in Particles from recv MPI buffer
bool 
MPIMoveIn::getreceive(ParticleArray& movein,
                      std::vector<PotentialModel>& moveintype,
                      std::vector<int>& incellid)
{
  MPI_Status mpistatus;
  int ret = MPI_Wait(&recv_request,&mpistatus);
  if(ret!=MPI_SUCCESS){
    std::cout << " MPI_Wait error at MPIMoveIn::getreceive(...)" << ret << std::endl;
  }
  //    std::cout << " make move in buffer " << recvbufferheader->NumberOfParticle << std::endl; 
  if( set_receive_buffer_pointers() == true ) {
    int offset = movein.size();
    int toffset = moveintype.size();
    int coffset = incellid.size();
    movein.resize(offset+(recvbufferheader->NumberOfParticle));
    moveintype.resize(toffset+(recvbufferheader->NumberOfParticle));
    incellid.resize(coffset+(recvbufferheader->NumberOfParticle));
    for(int i=0;i<recvbufferheader->NumberOfParticle;i++){
      movein[offset+i] = recvparticle[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfParticle;i++){
      moveintype[toffset+i] = recvtype[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfParticle;i++){
      incellid[coffset+i] = recvcellid[i];
    }
    return true;
  }
  std::cout << "getreceive fail" << std::endl;
  return false;
}
bool 
MPIMoveIn::getreceive(ParticleArray& movein,
                      std::vector<PotentialModel>& moveintype,
                      std::vector<int>& incellid,
                      std::vector<CovalentBondInfo::BondList>& inbondlistarray)
{
  MPI_Status mpistatus;
  int ret = MPI_Wait(&recv_request,&mpistatus);
  if(ret!=MPI_SUCCESS){
    std::cout << " MPI_Wait error at MPIMoveIn::getreceive(...)" << ret << std::endl;
  }
  //    std::cout << " make move in buffer " << recvbufferheader->NumberOfParticle << std::endl; 
  if( set_receive_buffer_pointers() == true ) {
    int offset = movein.size();
    int toffset = moveintype.size();
    int coffset = incellid.size();
    int boffset = inbondlistarray.size();
    movein.resize(offset+(recvbufferheader->NumberOfParticle));
    moveintype.resize(toffset+(recvbufferheader->NumberOfParticle));
    incellid.resize(coffset+(recvbufferheader->NumberOfParticle));
    inbondlistarray.resize(boffset+(recvbufferheader->NumberOfBondlist));
    for(int i=0;i<recvbufferheader->NumberOfParticle;i++){
      movein[offset+i] = recvparticle[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfParticle;i++){
      moveintype[toffset+i] = recvtype[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfParticle;i++){
      incellid[coffset+i] = recvcellid[i];
    }
    for(int i=0;i<recvbufferheader->NumberOfBondlist;i++){
#ifdef TRANSFER_RAW_BONDLIST
      inbondlistarray[boffset+i] = *(recvbondlistarray[i]);
#else
      inbondlistarray[boffset+i].unpack_int_array(recvbondpackarray[i]);
#endif
    }
    return true;
  }
  std::cout << "getreceive fail" << std::endl;
  return false;
}
//! request MPI asyncronous receive for Move in Particle
void 
MPIMoveIn::prepare_receive()
{
  int ret = MPI_Irecv(recv_buffer, recv_buffer_size, MPI_BYTE, 
                      target_rank, MPITagMove, mpi_comm_short, &recv_request);
  if(ret!=MPI_SUCCESS){
    std::cout << " MPI_Irecv error at MPIMoveIn::prepare_receive()" << ret << std::endl;
  }
}


/*!
  given : move_out_target_id, move_setid
 */
void
MPITargetNearestXYZ::makeMoveOutTarget()
{
  
}

/*!
  given : move_in_target_id
 */
void
MPITargetNearestXYZ::makeMoveInTarget()
{
  
}

/*!
  Methods of MPICommunicator
 */
  //! constructor 
MPICommunicator::MPICommunicator(const std::map<int,int>& uidtor, 
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
                                 CommPattern srcp,
                                 CommPattern mcp) 
  : unit_id(id), psize(ps), unitid_to_rank(uidtor), 
    mpi_comm_short(short_comm),
    sprf(), sprf_target_id(sprf_tid), setid(sid), srpf_setid(srpf_set), 
    allsendsetid(),
    rpsf(), rpsf_target_id(rpsf_tid), 
    targetsetid(tsid), alltargetsetid(),
    send_recv_communication_pattern(srcp),
    move_out(), move_in(),
    move_direct_target_id(mtid),
    move_setid(msid),
    move_out_particle(), move_out_cellid(), move_out_type(),
    move_depth(move_range),
    move_surface_ratio(movesurfaceratio),
    move_communication_pattern(mcp),
    node_geometry(nodediv,celldiv)
{
  /*
  if(unit_id==NO_SHORT_ID){
    std::cout << " This node not short " << std::endl;
    return;
  }
  */
  //  std::cout << "MPICommunicator(...) " << std::endl;
  if(send_recv_communication_pattern==Direct){
    if(unit_id==0){
      std::cout << "send_recv_communication_pattern Direct " << send_recv_communication_pattern << std::endl;
    }
  }else if(send_recv_communication_pattern==NearestXYZ){
    if(unit_id==0){
      std::cout << "send_recv_communication_pattern NearestXYZ " << send_recv_communication_pattern << std::endl;
    }
  }else{
    if(unit_id==0){
      std::cout << " not implement send_recv_communication_pattern " << send_recv_communication_pattern << std::endl;
    }
    send_recv_communication_pattern = Direct;
    if(unit_id==0){
      std::cout << "force send_recv_communication_pattern Direct " << send_recv_communication_pattern << std::endl;
    }
  }
  //  std::cout << "MPICommunicator::MPICommunicator(.....)" << std::endl;
  if(move_communication_pattern==Direct){
    if(unit_id==0){
      std::cout << "move_communication_pattern Direct " << move_communication_pattern << std::endl;
    }
  }else if(move_communication_pattern==NearestXYZ){
    if(unit_id==0){
      std::cout << "move_communication_pattern NearestXYZ " << move_communication_pattern << std::endl;
    }
  }else{
    if(unit_id==0){
      std::cout << " not implement move_communication_pattern " << move_communication_pattern << std::endl;
    }
    move_communication_pattern = Direct;
    if(unit_id==0){
      std::cout << "force move_communication_pattern Direct " << move_communication_pattern << std::endl;
    }
  }
  node_position = node_geometry.getAbsolutePosition(unit_id);
#ifdef OVERLAP
  MPICHECKER::reqs = NULL;
  MPICHECKER::status = NULL;
#endif
  construct_sendrecv();
  make_id_to_index(setid,setid_to_index); 
  make_idset_to_index(move_setid,setid_to_move_target_index);
  //  std::cout << "MPICommunicator(...) done " << std::endl;
}

//! copy constructor
MPICommunicator::MPICommunicator(const MPICommunicator& comm)
  : unit_id(comm.unit_id), psize(comm.psize), 
    unitid_to_rank(comm.unitid_to_rank),
    mpi_comm_short(comm.mpi_comm_short),
    //    sprf(comm.sprf), 
    sprf_target_id(comm.sprf_target_id),
    setid(comm.setid), setid_to_index(comm.setid_to_index), 
    srpf_setid(comm.srpf_setid),
    allsendsetid(comm.allsendsetid),
    //    rpsf(comm.rpsf), 
    rpsf_target_id(comm.rpsf_target_id),
    targetsetid(comm.targetsetid), alltargetsetid(comm.alltargetsetid),
    send_recv_communication_pattern(comm.send_recv_communication_pattern),
    //    move_out(comm.move_out), move_in(comm.move_in),
    move_direct_target_id(comm.move_direct_target_id),
    move_out_target_id(comm.move_out_target_id),  
    move_in_target_id(comm.move_in_target_id),  
    move_setid(comm.move_setid),
    setid_to_move_target_index(comm.setid_to_move_target_index),
    move_out_particle(comm.move_out_particle), 
    move_out_cellid(comm.move_out_cellid), 
    move_out_type(comm.move_out_type), 
    move_depth(comm.move_depth),
    move_surface_ratio(comm.move_surface_ratio),
    move_communication_pattern(comm.move_communication_pattern),
    node_geometry(comm.node_geometry)
{
  //  std::cout << "MPICommunicator::MPICommunicator(const MPICommunicator& comm)" << std::endl;
  node_position = node_geometry.getAbsolutePosition(unit_id);
  construct_sendrecv();
  make_id_to_index(setid,setid_to_index); 
  make_idset_to_index(move_setid,setid_to_move_target_index);
}

/*
MPICommunicator& MPICommunicator::operator=(const MPICommunicator& comm)
{
  std::cout << "MPICommunicator::operator=(const MPICommunicator& comm)" << std::endl;
  if(this != &comm){
    unit_id = comm.unit_id;
    psize = comm.psize;
    unitid_to_rank = comm.unitid_to_rank;
    //    sprf = comm.sprf;
    sprf_target_id = comm.sprf_target_id;
    setid = comm.setid;
    setid_to_index = comm.setid_to_index;
    srpf_setid = comm.srpf_setid;
    allsendsetid = comm.allsendsetid;
    //    rpsf = comm.rpsf;
    rpsf_target_id = comm.rpsf_target_id;
    targetsetid = comm.targetsetid;
    alltargetsetid = comm.alltargetsetid;
    send_recv_communication_pattern = comm.send_recv_communication_pattern;
    //    move_out = comm.move_out;
    //    move_in = comm.move_in;
    move_direct_target_id = comm.move_direct_target_id;
    move_out_target_id = comm.move_out_target_id;
    move_in_target_id = comm.move_in_target_id;
    move_setid = comm.move_setid;
    setid_to_move_target_index = comm.setid_to_move_target_index;
    move_out_particle = comm.move_out_particle;
    move_out_cellid = comm.move_out_cellid;
    move_out_type = comm.move_out_type;
    move_depth = comm.move_depth;
    move_surface_ratio = comm.move_surface_ratio;
    move_communication_pattern = comm.move_communication_pattern;
    node_geometry = comm.node_geometry;

    node_position = node_geometry.getAbsolutePosition(unit_id);
    construct_sendrecv();
    make_id_to_index(setid,setid_to_index); 
    make_idset_to_index(move_setid,setid_to_move_target_index);
  }

  return *this;
}
*/

MPICommunicator::~MPICommunicator()
{
  //  std::cout << "MPICommunicator::~MPICommunicator()" << std::endl;
  rpsf.clear();
  sprf.clear();
}

//! construct Send Particle, Recv Particle, Move Particle
void 
MPICommunicator::construct_sendrecv()
{
  //  std::cout << "MPICommunicator::construct_sendrecv() " << std::endl;
  {
    SpaceVector<int> plus_depth(0,0,0);
    SpaceVector<int> minus_depth(0,0,0);
    node_geometry.path_depth(plus_depth,minus_depth,rpsf_target_id,unit_id);
    if(unit_id==0){
      if((plus_depth.x!=-minus_depth.x)
         ||(plus_depth.y!=-minus_depth.y)
         ||(plus_depth.z!=-minus_depth.z)){
        std::cout << " depth imbalance " <<plus_depth << minus_depth << std::endl;
      }
    }
    recv_depth.x = std::max(plus_depth.x,-minus_depth.x);
    recv_depth.y = std::max(plus_depth.y,-minus_depth.y);
    recv_depth.z = std::max(plus_depth.z,-minus_depth.z);
    recv_direction.x = 0;
    if(plus_depth.x==0)recv_direction.x=-1;
    if(minus_depth.x==0)recv_direction.x=1;
    recv_direction.y = 0;
    if(plus_depth.y==0)recv_direction.y=-1;
    if(minus_depth.y==0)recv_direction.y=1;
    recv_direction.z = 0;
    if(plus_depth.z==0)recv_direction.z=-1;
    if(minus_depth.z==0)recv_direction.z=1;
    if(unit_id==0){
      std::cout << " depth and direction " <<recv_depth << recv_direction << std::endl;
    }
  }
  if(send_recv_communication_pattern==NearestXYZ){

    xplus_target_id = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>(+1,0,0),unit_id);
    xminus_target_id = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>(-1,0,0),unit_id);
    yplus_target_id = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>(0,+1,0),unit_id);
    yminus_target_id = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>(0,-1,0),unit_id);
    zplus_target_id = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>(0,0,+1),unit_id);
    zminus_target_id = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>(0,0,-1),unit_id);
    size_t num_set = setid.size();
    size_t y_num_set = num_set*(recv_depth.x*2+1);
    size_t z_num_set = y_num_set*(recv_depth.y*2+1);
    if((unit_id==0)&&(DebugLog::verbose>1)){
      std::cout << "psize " << psize;
      std::cout << " num_set " << num_set;
      size_t ssets=0, rsets=0;
      for(size_t s=0;s<srpf_setid.size();s++){
        ssets += srpf_setid[s].size();
      }
      for(size_t s=0;s<targetsetid.size();s++){
        rsets += targetsetid[s].size();
      }
      std::cout << " num_send_setid " << ssets;
      std::cout << " num_receive_id " << rsets << std::endl;
    }
    send_typerange.resize(psize*num_set);
    sprf.resize(6,MPISendParticleReceiveForce());
    sprf[0] = MPISendParticleReceiveForce(unit_id,psize*num_set, 
                                          mpi_comm_short,
                                          recv_depth.x+1);
    sprf[0].make_buffer(num_set);
    sprf[0].target_id = xplus_target_id;
    sprf[0].target_rank = unitid_to_rank[xplus_target_id];
    sprf[1] = MPISendParticleReceiveForce(unit_id,psize*num_set, 
                                          mpi_comm_short,
                                          recv_depth.x+1);
    sprf[1].make_buffer(num_set);
    sprf[1].target_id = xminus_target_id;
    sprf[1].target_rank = unitid_to_rank[xminus_target_id];
    sprf[2] = MPISendParticleReceiveForce(unit_id, psize*y_num_set, 
                                          mpi_comm_short,
                                          recv_depth.y+1);
    sprf[2].make_buffer(y_num_set);
    sprf[2].target_id = yplus_target_id;
    sprf[2].target_rank = unitid_to_rank[yplus_target_id];
    sprf[3] = MPISendParticleReceiveForce(unit_id, psize*y_num_set, 
                                          mpi_comm_short,
                                          recv_depth.y+1);
    sprf[3].make_buffer(y_num_set);
    sprf[3].target_id = yminus_target_id;
    sprf[3].target_rank = unitid_to_rank[yminus_target_id];
    sprf[4] = MPISendParticleReceiveForce(unit_id,psize*z_num_set, 
                                          mpi_comm_short,
                                          recv_depth.z+1);
    sprf[4].make_buffer(z_num_set);
    sprf[4].target_id = zplus_target_id;
    sprf[4].target_rank = unitid_to_rank[zplus_target_id];
    sprf[5] = MPISendParticleReceiveForce(unit_id,psize*z_num_set, 
                                          mpi_comm_short,
                                          recv_depth.z+1);
    sprf[5].make_buffer(z_num_set);
    sprf[5].target_id = zminus_target_id;
    sprf[5].target_rank = unitid_to_rank[zminus_target_id];
    rpsf.resize(6,MPIReceiveParticleSendForce());
    rpsf[0] = MPIReceiveParticleSendForce(unit_id,psize*num_set, 
                                          mpi_comm_short,
                                          recv_depth.x+1);
    rpsf[0].make_buffer(num_set);
    rpsf[0].target_id = xplus_target_id;
    rpsf[0].target_rank = unitid_to_rank[xplus_target_id];
    rpsf[1] = MPIReceiveParticleSendForce(unit_id,psize*num_set, 
                                          mpi_comm_short,
                                          recv_depth.x+1);
    rpsf[1].make_buffer(num_set);
    rpsf[1].target_id = xminus_target_id;
    rpsf[1].target_rank = unitid_to_rank[xminus_target_id];
    rpsf[2] = MPIReceiveParticleSendForce(unit_id, psize*y_num_set, 
                                          mpi_comm_short,
                                          recv_depth.y+1);
    rpsf[2].make_buffer(y_num_set);
    rpsf[2].target_id = yplus_target_id;
    rpsf[2].target_rank = unitid_to_rank[yplus_target_id];
    rpsf[3] = MPIReceiveParticleSendForce(unit_id, psize*y_num_set, 
                                          mpi_comm_short,
                                          recv_depth.y+1);
    rpsf[3].make_buffer(y_num_set);
    rpsf[3].target_id = yminus_target_id;
    rpsf[3].target_rank = unitid_to_rank[yminus_target_id];
    rpsf[4] = MPIReceiveParticleSendForce(unit_id,psize*z_num_set, 
                                          mpi_comm_short,
                                          recv_depth.z+1);
    rpsf[4].make_buffer(z_num_set);
    rpsf[4].target_id = zplus_target_id;
    rpsf[4].target_rank = unitid_to_rank[zplus_target_id];
    rpsf[5] = MPIReceiveParticleSendForce(unit_id,psize*z_num_set, 
                                          mpi_comm_short,
                                          recv_depth.z+1);
    rpsf[5].make_buffer(z_num_set);
    rpsf[5].target_id = zminus_target_id;
    rpsf[5].target_rank = unitid_to_rank[zminus_target_id];

    /*
    for(int i=0;i<6;i++){
      std::cout << "target rank " <<sprf[i].target_rank << ":" << rpsf[i].target_rank<< "  number_of_stage " << sprf[i].number_of_particle.size() << " " << rpsf[i].number_of_particle.size() << std::endl;
    }
    */

    allsendsetid.clear();
    sprf_number_of_target = sprf_target_id.size();
    number_of_send_set = 0;
    std::set<int> att;
    for(int i=0;i<sprf_number_of_target;i++){
      for(std::vector<int>::iterator it = srpf_setid[i].begin();
          it!=srpf_setid[i].end(); ++it){
        std::pair<std::set<int>::iterator,bool> res(att.insert((*it)));
        if(res.second==true){
          number_of_send_set++;
          allsendsetid.push_back(*it);
        }
      }
    }

    rpsf_number_of_target = rpsf_target_id.size();
    number_of_target_set = 0;
    for(int i=0;i<rpsf_number_of_target;i++){
      number_of_target_set += rpsf_number_of_target;
      for(std::vector<int>::iterator it = targetsetid[i].begin();
          it!=targetsetid[i].end(); ++it){
        alltargetsetid.push_back(*it);
      }
    }
  }else{ // ! send_recv_communication_pattern==NearestXYZ
    // construct Send Particle
    sprf_number_of_target = sprf_target_id.size();
    sprf.clear();
    size_t num_set = setid.size();
    sprf.resize(sprf_number_of_target,MPISendParticleReceiveForce());
    for(int i=0;i<sprf_number_of_target;i++){
      sprf[i] = MPISendParticleReceiveForce(unit_id,psize*num_set,
                                            mpi_comm_short);
      sprf[i].make_buffer(num_set);
      sprf[i].target_id = sprf_target_id[i];
      sprf[i].target_rank = unitid_to_rank[sprf_target_id[i]];
    }
    if(unit_id==0)
    {
      std::cout << "number of send target node " << sprf_number_of_target << std::endl;
    }
    // construct Recv Particle
    rpsf_number_of_target = rpsf_target_id.size();
    rpsf.clear();
    rpsf.resize(rpsf_number_of_target,MPIReceiveParticleSendForce());
    /*
      targetparticle.resize(psize*rpsf_number_of_target);
      targetparticle.clear();
      target_range.resize(rpsf_number_of_target);
      target_range.clear();
    */
    number_of_target_set = 0;
    for(int i=0;i<rpsf_number_of_target;i++){
      rpsf[i] = MPIReceiveParticleSendForce(unit_id,
                                            psize*targetsetid[i].size(), 
                                            mpi_comm_short);
      rpsf[i].make_buffer(setid.size());
      rpsf[i].target_id = rpsf_target_id[i];
      rpsf[i].target_rank = unitid_to_rank[rpsf_target_id[i]];
      number_of_target_set += targetsetid[i].size();
      for(std::vector<int>::iterator it = targetsetid[i].begin();
          it!=targetsetid[i].end(); ++it){
        alltargetsetid.push_back(*it);
      }
    }
  }

#ifdef OVERLAP
  if(MPICHECKER::reqs==NULL){
    MPICHECKER::reqs = new MPI_Request[rpsf_number_of_target+sprf_number_of_target];
    //    printf("make MPICHECKER::reqs %d\n",rpsf_number_of_target+sprf_number_of_target);
  }
  if(MPICHECKER::status==NULL){
    MPICHECKER::status = new MPI_Status[rpsf_number_of_target+sprf_number_of_target];
    //    printf("make MPICHECKER::status %d\n",rpsf_number_of_target+sprf_number_of_target);
  }
  MPICHECKER::count=0;
  for(int t=0;t<rpsf_number_of_target;t++){
    rpsf[t].recv_requestp = &(MPICHECKER::reqs[MPICHECKER::count]);
    MPICHECKER::count++;
  }
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].send_requestp = &(MPICHECKER::reqs[MPICHECKER::count]);
    MPICHECKER::count++;
  }
#else
  for(int t=0;t<rpsf_number_of_target;t++){
    rpsf[t].recv_requestp = new MPI_Request();
  }
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].send_requestp = new MPI_Request();
  }
#endif
    
  if(move_communication_pattern==NearestXYZ){
    move_number_of_target = 6;
    number_of_move_out_particle = 2*move_depth.x+2*move_depth.y+2*move_depth.z;
    number_of_move_in_particle = 1;
    move_imd_particle.resize(number_of_move_out_particle);
    move_imd_cellid.resize(number_of_move_out_particle);
    move_imd_type.resize(number_of_move_out_particle);
    move_imd_bondlistarray.resize(number_of_move_out_particle);
    for(int i=0;i<number_of_move_out_particle;i++){
      move_imd_particle[i].reserve(psize);
      move_imd_cellid[i].reserve(psize);
      move_imd_type[i].reserve(psize);
      move_imd_bondlistarray[i].reserve(psize);
    }
    move_out_target_id.resize(move_number_of_target);
    move_in_target_id.resize(move_number_of_target);
    SpaceVector<int> relative_target(0,0,0);
    relative_target.x=+1;
    relative_target.y=0;
    relative_target.z=0;
    move_out_target_id[0] = node_geometry.getNodeIDofRelativePosition(relative_target,unit_id);
    relative_target.x=-1;
    move_out_target_id[1] = node_geometry.getNodeIDofRelativePosition(relative_target,unit_id);
    relative_target.x=0;
    relative_target.y=+1;
    move_out_target_id[2] = node_geometry.getNodeIDofRelativePosition(relative_target,unit_id);
    relative_target.y=-1;
    move_out_target_id[3] = node_geometry.getNodeIDofRelativePosition(relative_target,unit_id);
    relative_target.y=0;
    relative_target.z=+1;
    move_out_target_id[4] = node_geometry.getNodeIDofRelativePosition(relative_target,unit_id);
    relative_target.z=-1;
    move_out_target_id[5] = node_geometry.getNodeIDofRelativePosition(relative_target,unit_id);
    for(int t=0;t<move_number_of_target;t++){
      move_in_target_id[t] = move_out_target_id[t];
    }
  }else if(move_communication_pattern==Direct){
    move_number_of_target = move_direct_target_id.size();
    number_of_move_out_particle = move_number_of_target;
    number_of_move_in_particle = move_number_of_target;
  }
  // reserve move in/out buffers
  all_move_out_particle.reserve(psize*number_of_move_out_particle);
  all_move_out_cellid.reserve(psize*number_of_move_out_particle);
  all_move_out_type.reserve(psize*number_of_move_out_particle);
  all_move_out_bondlistarray.reserve(psize*number_of_move_out_particle);
  move_out_particle.resize(number_of_move_out_particle);
  move_out_cellid.resize(number_of_move_out_particle);
  move_out_type.resize(number_of_move_out_particle);
  move_out_bondlistarray.resize(number_of_move_out_particle);        
  for(int i=0;i<number_of_move_out_particle;i++){
    move_out_particle[i].reserve(psize);
    move_out_cellid[i].reserve(psize);
    move_out_type[i].reserve(psize);
    move_out_bondlistarray[i].reserve(psize);
  }
  move_in_particle.reserve(psize*number_of_move_in_particle);
  move_in_cellid.reserve(psize*number_of_move_in_particle);
  move_in_type.reserve(psize*number_of_move_in_particle);
  move_in_bondlistarray.reserve(psize*number_of_move_in_particle);
  //    std::cout << "move_in_particle reserve " << psize*move_number_of_target << std::endl;
  
// construct Move Particle
  move_out.clear();
  move_out.resize(move_number_of_target,MPIMoveOut());
  move_in.clear();
  move_in.resize(move_number_of_target,MPIMoveIn());
  if(move_communication_pattern==NearestXYZ){
    int number_of_full_move_setid=0;
    for(std::vector<int>::size_type i=0;i<move_direct_target_id.size();i++){
      number_of_full_move_setid+=move_setid[i].size();
    }
    if(unit_id==0){
      std::cout << "ratio of move candidate " << move_surface_ratio  << std::endl;
    }
    size_t move_candidate = size_t(psize*setid.size()*move_surface_ratio);
    for(int i=0;i<move_number_of_target;i++){
      move_out[i] = MPIMoveOut(unit_id,move_candidate, mpi_comm_short);
      move_out[i].make_buffer();
      move_out[i].target_id = move_out_target_id[i];
      move_out[i].target_rank = unitid_to_rank[move_out_target_id[i]];
      move_in[i] = MPIMoveIn(unit_id,move_candidate, mpi_comm_short);
      move_in[i].make_buffer();
      move_in[i].target_id = move_out_target_id[i];
      move_in[i].target_rank = unitid_to_rank[move_out_target_id[i]];
      if(i&0x1) move_candidate*=3;
    }
    /*
    if(unit_id==0){
      std::cout << " geom " << node_geometry.size << std::endl;
      std::cout << " out in " << std::endl;
      for(int i=0;i<move_number_of_target;i++){
        std::cout << move_out[i].target_id << " " << move_in[i].target_id << std::endl;
      }
    }
    */
  }else if(move_communication_pattern==Direct){    
    for(int i=0;i<move_number_of_target;i++){
      move_out[i] = MPIMoveOut(unit_id,psize*move_setid[i].size(),
                               mpi_comm_short);
      move_out[i].make_buffer();
      move_out[i].target_id = move_direct_target_id[i];
      move_out[i].target_rank = unitid_to_rank[move_direct_target_id[i]];
      move_in[i] = MPIMoveIn(unit_id,psize*move_setid[i].size(),
                             mpi_comm_short);
      move_in[i].make_buffer();
      move_in[i].target_id = move_direct_target_id[i];
      move_in[i].target_rank = unitid_to_rank[move_direct_target_id[i]];
    }
  }


#ifdef TIMER_DETAIL
  timer_move = PerfCounter::add_target(std::string("Move"));
  timer_move_comm = PerfCounter::add_target(std::string("Move_comm"));
#endif
}


//! resize ParticleRange of j-Partice
void 
MPICommunicator::resize_receive_data(std::vector<ParticleRange>& targetparticlerange)
{
  targetparticlerange.resize(rpsf_number_of_target);
}

//! set own particle to each send buffer
/*!
  for continuous particle
*/
void 
MPICommunicator::setSendParticle(ParticleArray& myparticle,
                                 std::vector<TypeRange>& typerangearray,
                                 std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].setsend(myparticle,typerangearray,bondlistarray,setid);
  }
}

//! set own particle to each send buffer
/*!
  for particle with gap
*/
template<class PA>
void 
MPICommunicator::setSendParticlesubset_onlyPosition(const PA& myparticle,
                                       const std::vector<TypeRange>& typerangearray)
{
  int t;
#ifdef MPIPARALLEL_OPENMP
#ifdef _OPENMP
#pragma omp parallel for private(t)
#endif
#endif
  for(t=0;t<sprf_number_of_target;t++){
    sprf[t].setsendsubsetpos(myparticle,typerangearray,setid,srpf_setid[t]);
  }
}

//! get received particle from each send buffer
template<class GPA>
void 
MPICommunicator::getReceiveParticle_onlyPosition(GPA& targetparticle)
{
  //    std::cout << "recieve from " << rpsf_number_of_target << " target" << std::endl;
  //    std::cout << "recieve";
  int t;

#ifdef MIX_MPI_WAIT
  //#pragma omp parallel for private(t)
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].getreceivepos(targetparticle);
    //      std::cout << " " << target_range[t].begin << "-" << target_range[t].end;
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
  //    std::cout << std::endl;
}

//! set own particle to each send buffer
/*!
  for particle with gap
*/
template<class PA>
void 
MPICommunicator::setSendParticlesubset(const PA& myparticle,
                                       const std::vector<TypeRange>& typerangearray,
                                       const std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  int t;
#ifdef MPIPARALLEL_OPENMP
#ifdef _OPENMP
#pragma omp parallel for private(t)
#endif
#endif
  for(t=0;t<sprf_number_of_target;t++){
    sprf[t].setsendsubset(myparticle,typerangearray,bondlistarray,setid,srpf_setid[t]);
  }
}

//! get received particle from each send buffer
template<class PA>
void 
MPICommunicator::getReceiveParticle(PA& targetparticle,
                                    std::vector<ParticleRange>& target_range,
                                    std::vector<TypeRange>& targettyperange,
                                    std::vector<CovalentBondInfo::BondList>& targetbond,
                                    std::vector<int>& recvsetid,
                                    std::map<int,int>& recvsetid_to_index)
{
  //    std::cout << "recieve from " << rpsf_number_of_target << " target" << std::endl;
  //    std::cout << "recieve";
  int t;
#ifdef MIX_MPI_WAIT
  target_range.resize(rpsf_number_of_target);
  for(t=0;t<rpsf_number_of_target;t++){
    target_range[t].begin = targetparticle.size();
    rpsf[t].getreceive(targetparticle,target_range[t],
                       targettyperange, targetbond, recvsetid,
                       recvsetid_to_index);
    //      std::cout << " " << target_range[t].begin << "-" << target_range[t].end;
  }
  //    std::cout << std::endl;
#else
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].wait_complete_recv();
  }
  target_range.resize(rpsf_number_of_target);
#if 0
  for(t=0;t<rpsf_number_of_target;t++){
    target_range[t].begin = targetparticle.size();
    rpsf[t].getreceive_nowait(targetparticle,target_range[t],
                              targettyperange, targetbond, recvsetid,
                              recvsetid_to_index);
  }
#else
  int *poffset, *setoffset, *bondoffset;
  poffset = new int[rpsf_number_of_target+1];
  setoffset = new int[rpsf_number_of_target+1];
  bondoffset = new int[rpsf_number_of_target+1];
  poffset[0] = targetparticle.size();
  setoffset[0] = targettyperange.size();
  bondoffset[0] = targetbond.size();
  /// serial
  for(t=0;t<rpsf_number_of_target;t++){
    int pnum=0, setnum=0, bondnum=0;
    rpsf[t].getreceive_number(pnum,setnum,bondnum);
    poffset[t+1] = poffset[t] + pnum;
    setoffset[t+1] = setoffset[t] + setnum;
    bondoffset[t+1] = bondoffset[t] + bondnum;
    target_range[t].begin = poffset[t];
    target_range[t].end = poffset[t+1];
  }
  targetparticle.resize(poffset[rpsf_number_of_target]);
  targettyperange.resize(setoffset[rpsf_number_of_target]);
  targetbond.resize(setoffset[rpsf_number_of_target]);
  recvsetid.resize(setoffset[rpsf_number_of_target]);
#ifdef MPIPARALLEL_OPENMP
#ifdef _OPENMP
#pragma omp parallel for private(t)
#endif
#endif
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].getreceive_offset(targetparticle, target_range[t],
                              targettyperange, targetbond, recvsetid,
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
#endif
}

//! request transfer Particle
void 
MPICommunicator::transferParticle()
{
  for(int t=0;t<rpsf_number_of_target;t++){
    rpsf[t].prepare_receive();
  }
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].send();
  }

}

//! wait/check Particle Transfer
void 
MPICommunicator::wait_complete_Particle()
{
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].wait_complete_send();
  }
}

//! print some information
void 
MPICommunicator::print()
{
  std::cout << "psize " << psize << std::endl;
  std::cout << "unit_id " << unit_id << std::endl;
  std::cout << sprf_target_id.size() << " sprf_target_id" ;
  for(std::vector<int>::iterator it = sprf_target_id.begin();it!=sprf_target_id.end();++it){
    std::cout << " " << (*it);
  }
  std::cout << std::endl;
    
  std::cout << rpsf_target_id.size()  << " rpsf_target_id" ;
  for(std::vector<int>::iterator it = rpsf_target_id.begin();it!=rpsf_target_id.end();++it){
    std::cout << " " << (*it);
  }
  std::cout << std::endl;
}

//! send own particle (no gap) and receive j-particle
void 
MPICommunicator::exchangeParticleArray(ParticleArray& myparticle,
                                       std::vector<TypeRange>& typerangearray,
                                       std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                       ParticleArray& targetparticle,
                                       std::vector<ParticleRange>& target_range,
                                       std::vector<int>& recvsetid,
                                       std::map<int,int>& recvsetid_to_index,
                                       std::vector<TypeRange>& targettyperange,
                                       std::vector<CovalentBondInfo::BondList>& targetbond)
{
  setSendParticle(myparticle,typerangearray,bondlistarray);
  transferParticle();
  recvsetid_to_index.clear();
  getReceiveParticle(targetparticle,target_range,targettyperange,targetbond,recvsetid,recvsetid_to_index);
  wait_complete_Particle();
}
//! send own particle and receive j-particle with measurement of MPI Time
void 
MPICommunicator::exchangeParticleArray(ParticleArray& myparticle,
                                       std::vector<TypeRange>& typerangearray,
                                       std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                       ParticleArray& targetparticle,
                                       std::vector<ParticleRange>& target_range,
                                       std::vector<int>& recvsetid,
                                       std::map<int,int>& recvsetid_to_index,
                                       std::vector<TypeRange>& targettyperange,
                                       std::vector<CovalentBondInfo::BondList>& targetbond,
                                       double& alltime, double& rectime)
{
  alltime -= MPI_Wtime();
  setSendParticle(myparticle,typerangearray,bondlistarray);
  transferParticle();
  rectime -= MPI_Wtime();
  recvsetid_to_index.clear();
  targettyperange.clear();
  getReceiveParticle(targetparticle,target_range,targettyperange,targetbond,recvsetid,recvsetid_to_index);
  rectime += MPI_Wtime();
  wait_complete_Particle();
  alltime += MPI_Wtime();
}

//! set j-force to each send
void 
MPICommunicator::setSendForce(const ForceArray& sendforce,
                              const std::vector<ParticleRange>& send_range) 
{
  for(int t=0;t<rpsf_number_of_target;t++){
    if(rpsf[t].setsend(sendforce,send_range[t])==false){
      printf("MPICommunicator::setSendForce fail\n");
    }
  }
}

void 
MPICommunicator::setSendForce_with_index(const ForceArray& sendforce,
                                         const std::vector<ParticleRange>& send_range,
                                         const std::vector<int>& forceindexset) 
{
  int t;
#ifdef MPIPARALLEL_OPENMP
#ifdef _OPENMP
#pragma omp parallel for private(t)
#endif
#endif
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].setsend_with_index(sendforce,send_range[t],forceindexset);
  }
}
void 
MPICommunicator::setSendForce_indexed(const ForceArray& sendforce,
                                      const std::vector<ParticleRange>& send_range,
                                      const std::vector<int>& forceindexset) 
{
  int t;
#ifdef MPIPARALLEL_OPENMP
#ifdef _OPENMP
#pragma omp parallel for private(t)
#endif
#endif
  for(t=0;t<rpsf_number_of_target;t++){
    rpsf[t].setsend_indexed(sendforce,send_range[t],forceindexset);
  }
}

//! get received force
void 
MPICommunicator::getReceiveForcesubset(ForceArray& recvforce)
{
  int t;
  for(t=0;t<sprf_number_of_target;t++){
    sprf[t].getreceivesubset(recvforce);
  }
}
void 
MPICommunicator::getReceiveForcesubset_with_index(ForceArray& recvforce)
{
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].getreceive_with_index(recvforce);
  }
}
void 
MPICommunicator::getReceiveForcesubset_indexed(ForceArray& recvforce)
{
#ifdef ORDERED_RECV
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].getreceive_indexed(recvforce);
  }
#else
  MPI_Request *receive_requsets = new MPI_Request[sprf_number_of_target];
  int ti;
  for(ti=0;ti<sprf_number_of_target;ti++){
    receive_requsets[ti] = sprf[ti].receive_requset;
  }
  int t;
  int result;
  for(ti=0;ti<sprf_number_of_target;ti++){
    result = MPI_Waitany(sprf_number_of_target, receive_requsets, &t, MPI_STATUS_IGNORE);
    if(result!=MPI_SUCCESS){
      printf("Fail MPI_Waitany at MPICommunicator::getReceiveForcesubset_indexed(ForceArray& forcearray) %d\n",result);
    }
    sprf[t].getreceive_indexed_nowait(recvforce);
  }
#endif
}

//! request transfer Force
void 
MPICommunicator::transferForce()
{
  for(int t=0;t<sprf_number_of_target;t++){
    sprf[t].prepare_receive();
  }
  for(int t=0;t<rpsf_number_of_target;t++){
    rpsf[t].send();
  }
}

//! wait/check Force Transfer
void 
MPICommunicator::wait_complete_Force()
{
  for(int t=0;t<rpsf_number_of_target;t++){
    rpsf[t].wait_complete_send();
  }
}

template<class PA, class GPA>
void MPICommunicator::particle_send_recv_onlyPosition_Direct_top_half(const PA& myparticle,
                                            const std::vector<TypeRange>& typerangearray,
                                            const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                            GPA& targetparticle,
                                            const std::vector<ParticleRange>& target_range,
                                            const std::vector<int>& recvsetid,
                                            const std::map<int,int>& recvsetid_to_index,
                                            const std::vector<TypeRange>& targettyperange,
                                            const std::vector<CovalentBondInfo::BondList>& targetbond)
{
  setSendParticlesubset_onlyPosition(myparticle,typerangearray);
  transferParticle();
#ifdef OVERLAP
  MPICHECKER::done = false;
#endif
}
template<class PA, class GPA>
void MPICommunicator::particle_send_recv_onlyPosition_Direct_bottom_half(const PA& myparticle,
                                            const std::vector<TypeRange>& typerangearray,
                                            const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                            GPA& targetparticle,
                                            const std::vector<ParticleRange>& target_range,
                                            const std::vector<int>& recvsetid,
                                            const std::map<int,int>& recvsetid_to_index,
                                            const std::vector<TypeRange>& targettyperange,
                                            const std::vector<CovalentBondInfo::BondList>& targetbond)
{
  getReceiveParticle_onlyPosition(targetparticle);
}
template<class PA, class GPA>
void MPICommunicator::particle_send_recv_onlyPosition_Direct(const PA& myparticle,
                                            const std::vector<TypeRange>& typerangearray,
                                            const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                            GPA& targetparticle,
                                            const std::vector<ParticleRange>& target_range,
                                            const std::vector<int>& recvsetid,
                                            const std::map<int,int>& recvsetid_to_index,
                                            const std::vector<TypeRange>& targettyperange,
                                            const std::vector<CovalentBondInfo::BondList>& targetbond)
{
  setSendParticlesubset_onlyPosition(myparticle,typerangearray);
  transferParticle();
  getReceiveParticle_onlyPosition(targetparticle);
}

template<class PA, class GPA>
void MPICommunicator::particle_send_recv_Direct_top_half(const PA& myparticle,
                                            const std::vector<TypeRange>& typerangearray,
                                            const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                            GPA& targetparticle,
                                            std::vector<ParticleRange>& target_range,
                                            std::vector<int>& recvsetid,
                                            std::map<int,int>& recvsetid_to_index,
                                            std::vector<TypeRange>& targettyperange,
                                            std::vector<CovalentBondInfo::BondList>& targetbond)
{
  setSendParticlesubset(myparticle,typerangearray,bondlistarray);
  transferParticle();
#ifdef OVERLAP
  MPICHECKER::done = false;
#endif
}
template<class PA, class GPA>
void MPICommunicator::particle_send_recv_Direct_bottom_half(const PA& myparticle,
                                            const std::vector<TypeRange>& typerangearray,
                                            const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                            GPA& targetparticle,
                                            std::vector<ParticleRange>& target_range,
                                            std::vector<int>& recvsetid,
                                            std::map<int,int>& recvsetid_to_index,
                                            std::vector<TypeRange>& targettyperange,
                                            std::vector<CovalentBondInfo::BondList>& targetbond)
{
  getReceiveParticle(targetparticle,target_range,targettyperange,targetbond,recvsetid,recvsetid_to_index);
}
template<class PA, class GPA>
void MPICommunicator::particle_send_recv_Direct(const PA& myparticle,
                                            const std::vector<TypeRange>& typerangearray,
                                            const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                            GPA& targetparticle,
                                            std::vector<ParticleRange>& target_range,
                                            std::vector<int>& recvsetid,
                                            std::map<int,int>& recvsetid_to_index,
                                            std::vector<TypeRange>& targettyperange,
                                            std::vector<CovalentBondInfo::BondList>& targetbond)
{
  setSendParticlesubset(myparticle,typerangearray,bondlistarray);
  transferParticle();
  getReceiveParticle(targetparticle,target_range,targettyperange,targetbond,recvsetid,recvsetid_to_index);
}

template<class GPA>
void 
MPICommunicator::particle_send_recv_one_axis_onlyPosition(int depth, int recv_direction,
                                             int pf_plus, int pf_minus,
                                             int& target,
                                             GPA& targetparticle,
                                             const std::vector<ParticleRange>& target_range,
                                             const std::vector<int>& recvsetid,
                                             const std::map<int,int>& recvsetid_to_index,
                                             const std::vector<TypeRange>& targettyperange,
                                             const std::vector<CovalentBondInfo::BondList>& targetbond)
{
  if(recv_direction>-1){
    sprf[pf_minus].stage=0;
    rpsf[pf_plus].stage=0;
  }
  if(recv_direction<+1){
    sprf[pf_plus].stage=0;
    rpsf[pf_minus].stage=0;
  }
  for(int i=1;i<=depth;i++){
    if(i==1){
      // set all_receive to send to -
      if(recv_direction>-1){
        sprf[pf_minus].setsendsubsetpos(targetparticle,targettyperange,
                                        recvsetid,recvsetid);
      }
      // set all_receive to send to {+1,0,0}
      if(recv_direction<+1){
        sprf[pf_plus].setsendsubsetpos(targetparticle,targettyperange,
                                       recvsetid,recvsetid);
      }
    }else{
      // set buff+ to send to {-1,0,0}
      if(recv_direction>-1){
        size_t rssize = copy_position_buffer(rpsf[pf_plus].recv_buffer, 
                                             sprf[pf_minus].send_buffer);
        PositionBufferHeader *passheader=rpsf[pf_plus].recvposbufferheader;
        sprf[pf_minus].set_send_posbuffer_pointers(passheader->NumberOfPosition);
        if(rssize>sprf[pf_minus].send_buffer_size){
          std::cout << "copy recv > send " << std::endl;
        }
        sprf[pf_minus].number_of_particle[sprf[pf_minus].stage] = sprf[pf_minus].sendposbufferheader->NumberOfPosition;
      }
      // set buff- to send to {+1,0,0}
      if(recv_direction<+1){
        size_t rssize = copy_position_buffer(rpsf[pf_minus].recv_buffer,
                                             sprf[pf_plus].send_buffer);
        PositionBufferHeader *passheader=rpsf[pf_minus].recvposbufferheader;
        sprf[pf_plus].set_send_posbuffer_pointers(passheader->NumberOfPosition);
        if(rssize>sprf[pf_plus].send_buffer_size){
          std::cout << "copy recv > send " << std::endl;
        }
        sprf[pf_plus].number_of_particle[sprf[pf_plus].stage] = sprf[pf_plus].sendposbufferheader->NumberOfPosition;
      }
    }

    // prepare receive from node {+1,0,0}
    if(recv_direction>-1){
      //      std::cout << " rpsf pf_plus " << pf_plus << " recv from " << rpsf[pf_plus].target_id << std::endl;
      rpsf[pf_plus].prepare_receive();
    }
    // prepare receive from node {-1,0,0}
    if(recv_direction<+1){
      //      std::cout << " rpsf pf_minus " << pf_minus << " recv from " << rpsf[pf_minus].target_id  << std::endl;
      rpsf[pf_minus].prepare_receive();
    }
    // send to {-1,0,0}
    if(recv_direction>-1){
      //      std::cout << " srpf pf_minus " << pf_minus << " send  number_of_particle[" << sprf[pf_minus].stage << "] " << sprf[pf_minus].number_of_particle[sprf[pf_minus].stage] << " to " << sprf[pf_minus].target_id << std::endl;
      sprf[pf_minus].send();
    }
    // send to {+1,0,0}
    if(recv_direction<+1){
      //      std::cout << " srpf pf_plus " << pf_plus << " send  number_of_particle[" << sprf[pf_plus].stage << "] " << sprf[pf_plus].number_of_particle[sprf[pf_plus].stage] << " to " << sprf[pf_plus].target_id << std::endl;
      sprf[pf_plus].send();
    }
    // receive buff+ from {+1,0,0}
    // receive buff- from {-1,0,0}
    // append buff+ to all_receive
    // append buff- to all_receive
    if(recv_direction>-1){
      rpsf[pf_plus].getreceivepos(targetparticle);
      //      std::cout << " particle_send_recv_one_axis pf_plus " << pf_plus << " target " << target << " number_of_particle[" << rpsf[pf_plus].stage << "] " << rpsf[pf_plus].number_of_particle[rpsf[pf_plus].stage] << std::endl;
      sprf[pf_minus].stage++;
      rpsf[pf_plus].stage++;
      target++;
    }
    if(recv_direction<+1){
      rpsf[pf_minus].getreceivepos(targetparticle);
      //      std::cout << " particle_send_recv_one_axis pf_minus " << pf_minus << " target "<< target << " number_of_particle[" << rpsf[pf_minus].stage << "] " << rpsf[pf_minus].number_of_particle[rpsf[pf_minus].stage] << std::endl;
      sprf[pf_plus].stage++;
      rpsf[pf_minus].stage++;
      target++;
    }
    if(recv_direction>-1){
      sprf[pf_minus].wait_complete_send();
    }
    if(recv_direction<+1){
      sprf[pf_plus].wait_complete_send();
    }
  }
}

template<class PA, class GPA>
void MPICommunicator::particle_send_recv_onlyPosition_NearestXYZ(const PA& myparticle,
                                                const std::vector<TypeRange>& typerangearray,
                                                const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                                GPA& targetparticle,
                                                const std::vector<ParticleRange>& target_range,
                                                const std::vector<int>& recvsetid,
                                                const std::map<int,int>& recvsetid_to_index,
                                                const std::vector<TypeRange>& targettyperange,
                                                const std::vector<CovalentBondInfo::BondList>& targetbond)
{
  //    clear all_receive
  //  append my node cell to all_receive  // send particle set packed dense same as ghost
  int target=0;

  ParticleRange pr;
  pr.begin=0;
  pr.end = 0;
  for(std::vector<int>::size_type si=0;si<allsendsetid.size();si++){
    int sindex=0;
    while(setid[sindex]!=allsendsetid[si]){sindex++;}
    TypeRange tr = typerangearray[sindex];
    int size_of_set = tr.end-tr.begin;
    for(int i=0;i<size_of_set;i++){
      //      targetparticle[pr.end+i] = myparticle[tr.begin+i];
      getpos(targetparticle,pr.end+i) = getpos(myparticle,tr.begin+i);
      // getpos(targetparticle,pr.end+i).x = getpos(myparticle,tr.begin+i).x;
      // getpos(targetparticle,pr.end+i).y = getpos(myparticle,tr.begin+i).y;
      // getpos(targetparticle,pr.end+i).z = getpos(myparticle,tr.begin+i).z;
    }
    pr.end += size_of_set;
  }
  target++;
  
  //  std::cout << " particle_send_recv_one_axis x done " << std::endl;
  particle_send_recv_one_axis_onlyPosition(recv_depth.x,recv_direction.x,0,1,target,
                              targetparticle,target_range,
                              recvsetid,recvsetid_to_index,
                              targettyperange,targetbond);
  //  std::cout << " particle_send_recv_one_axis x done " << std::endl;
  particle_send_recv_one_axis_onlyPosition(recv_depth.y,recv_direction.y,2,3,target,
                              targetparticle,target_range,
                              recvsetid,recvsetid_to_index,
                              targettyperange,targetbond);

  particle_send_recv_one_axis_onlyPosition(recv_depth.z,recv_direction.z,4,5,target,
                              targetparticle,target_range,
                              recvsetid,recvsetid_to_index,
                              targettyperange,targetbond);

}

template<class GPA>
void 
MPICommunicator::particle_send_recv_one_axis(int depth, int recv_direction,
                                             int pf_plus, int pf_minus,
                                             int& target,
                                             GPA& targetparticle,
                                             std::vector<ParticleRange>& target_range,
                                             std::vector<int>& recvsetid,
                                             std::map<int,int>& recvsetid_to_index,
                                             std::vector<TypeRange>& targettyperange,
                                             std::vector<CovalentBondInfo::BondList>& targetbond)
{
  if(recv_direction>-1){
    sprf[pf_minus].stage=0;
    rpsf[pf_plus].stage=0;
  }
  if(recv_direction<+1){
    sprf[pf_plus].stage=0;
    rpsf[pf_minus].stage=0;
  }
  for(int i=1;i<=depth;i++){
    if(i==1){
      // set all_receive to send to -
      if(recv_direction>-1){
        sprf[pf_minus].setsendsubset(targetparticle,targettyperange,
                                     targetbond,recvsetid,recvsetid);
      }
      // set all_receive to send to {+1,0,0}
      if(recv_direction<+1){
        sprf[pf_plus].setsendsubset(targetparticle,targettyperange,
                                     targetbond,recvsetid,recvsetid);
      }
    }else{
      // set buff+ to send to {-1,0,0}
      if(recv_direction>-1){
        size_t rssize = copy_particle_buffer(rpsf[pf_plus].recv_buffer, 
                                             sprf[pf_minus].send_buffer);
        ParticleBufferHeader *passheader=rpsf[pf_plus].recvbufferheader;
        sprf[pf_minus].set_send_buffer_pointers(passheader->NumberOfPosition,
                                                passheader->NumberOfCharge,
                                                passheader->NumberOfAtomType,
                                                passheader->NumberOfAtomID,
                                                passheader->NumberOfParticleSet,
                                                passheader->NumberOfBondlist);
        if(rssize>sprf[pf_minus].send_buffer_size){
          std::cout << "copy recv > send " << std::endl;
        }
        sprf[pf_minus].number_of_particle[sprf[pf_minus].stage] = sprf[pf_minus].sendbufferheader->NumberOfPosition;
      }
      // set buff- to send to {+1,0,0}
      if(recv_direction<+1){
        size_t rssize = copy_particle_buffer(rpsf[pf_minus].recv_buffer,
                                             sprf[pf_plus].send_buffer);
        ParticleBufferHeader *passheader=rpsf[pf_minus].recvbufferheader;
        sprf[pf_plus].set_send_buffer_pointers(passheader->NumberOfPosition,
                                                passheader->NumberOfCharge,
                                                passheader->NumberOfAtomType,
                                                passheader->NumberOfAtomID,
                                                passheader->NumberOfParticleSet,
                                                passheader->NumberOfBondlist);
	if(rssize>sprf[pf_plus].send_buffer_size){
          std::cout << "copy recv > send " << std::endl;
        }
        sprf[pf_plus].number_of_particle[sprf[pf_plus].stage] = sprf[pf_plus].sendbufferheader->NumberOfPosition;
      }
    }
    // prepare receive from node {+1,0,0}
    if(recv_direction>-1){
      //      std::cout << " rpsf pf_plus " << pf_plus << " recv from " << rpsf[pf_plus].target_id << std::endl;
      rpsf[pf_plus].prepare_receive();
    }
    // prepare receive from node {-1,0,0}
    if(recv_direction<+1){
      //      std::cout << " rpsf pf_minus " << pf_minus << " recv from " << rpsf[pf_minus].target_id  << std::endl;
      rpsf[pf_minus].prepare_receive();
    }
    // send to {-1,0,0}
    if(recv_direction>-1){
      //      std::cout << " srpf pf_minus " << pf_minus << " send  number_of_particle[" << sprf[pf_minus].stage << "] " << sprf[pf_minus].number_of_particle[sprf[pf_minus].stage] << " to " << sprf[pf_minus].target_id << std::endl;
      sprf[pf_minus].send();
    }
    // send to {+1,0,0}
    if(recv_direction<+1){
      //      std::cout << " srpf pf_plus " << pf_plus << " send  number_of_particle[" << sprf[pf_plus].stage << "] " << sprf[pf_plus].number_of_particle[sprf[pf_plus].stage] << " to " << sprf[pf_plus].target_id << std::endl;
      sprf[pf_plus].send();
    }
    // receive buff+ from {+1,0,0}
    // receive buff- from {-1,0,0}
    // append buff+ to all_receive
    // append buff- to all_receive
    if(recv_direction>-1){
      ParticleRange pr;
      pr.begin = targetparticle.size();
      rpsf[pf_plus].getreceive(targetparticle,
                               pr,
                               targettyperange, targetbond,
                               recvsetid, recvsetid_to_index);
      //      std::cout << " particle_send_recv_one_axis pf_plus " << pf_plus << " target " << target << " number_of_particle[" << rpsf[pf_plus].stage << "] " << rpsf[pf_plus].number_of_particle[rpsf[pf_plus].stage] << std::endl;
      sprf[pf_minus].stage++;
      rpsf[pf_plus].stage++;
      target_range.push_back(pr);
      target++;
    }
    if(recv_direction<+1){
      ParticleRange pr;
      pr.begin = targetparticle.size();
      rpsf[pf_minus].getreceive(targetparticle,
                                pr,
                                targettyperange, targetbond,
                                recvsetid, recvsetid_to_index);
      //      std::cout << " particle_send_recv_one_axis pf_minus " << pf_minus << " target "<< target << " number_of_particle[" << rpsf[pf_minus].stage << "] " << rpsf[pf_minus].number_of_particle[rpsf[pf_minus].stage] << std::endl;
      sprf[pf_plus].stage++;
      rpsf[pf_minus].stage++;
      target_range.push_back(pr);
      target++;
    }
    if(recv_direction>-1){
      sprf[pf_minus].wait_complete_send();
    }
    if(recv_direction<+1){
      sprf[pf_plus].wait_complete_send();
    }
    //    MPI_Barrier(mpi_comm_short);
  }
  /*
  if(unit_id==0){
    std::cout << " send to " << sprf[pf_minus].target_id;
    for(int s=0;s<sprf[pf_minus].stage;s++){
      std::cout << " " << sprf[pf_minus].number_of_particle[s] << " stage " << s;
    }
    std::cout << std::endl;
    std::cout << " send to " << sprf[pf_plus].target_id;
    for(int s=0;s<sprf[pf_plus].stage;s++){
      std::cout << " " << sprf[pf_plus].number_of_particle[s] << " stage " << s;
    }
    std::cout << std::endl;
  }
  */
}
inline void copy_particle_ghost(ParticleArray& ghost, const int j,
                           const ParticleArray& pa, const int i)
{
  ghost[j] = pa[i];
}
inline void copy_particle_ghost(GhostParticleArray& ghost, const int j,
                           const CombinedParticleArray& pa, const int i)
{
  ghost.poscharge[j] = pa.poscharge[i];
  ghost.atomtype[j] = pa.atomtype[i];
  ghost.force[j] = pa.force[i];
  // ghost.force[j].x = pa.force[j].x;
  // ghost.force[j].y = pa.force[j].y;
  // ghost.force[j].z = pa.force[j].z;
  ghost.atomid[j] = pa.atomid[i];
}
template<class PA, class GPA>
void MPICommunicator::particle_send_recv_NearestXYZ(const PA& myparticle,
                                                const std::vector<TypeRange>& typerangearray,
                                                const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                                GPA& targetparticle,
                                                std::vector<ParticleRange>& target_range,
                                                std::vector<int>& recvsetid,
                                                std::map<int,int>& recvsetid_to_index,
                                                std::vector<TypeRange>& targettyperange,
                                                std::vector<CovalentBondInfo::BondList>& targetbond)
{
  //    clear all_receive
  //  append my node cell to all_receive  // send particle set packed dense same as ghost
  int target=0;

  targetparticle.clear();  ///// different Direct (append, not clear)
  target_range.reserve((recv_depth.x*2+1)*(recv_depth.y*2+1)*(recv_depth.z*2+1));
  target_range.clear();
  ParticleRange pr;
  pr.begin=0;
  pr.end = 0;
  /*
  send_typerange.resize(setid.size());
  targettyperange.resize(setid.size());
  recvsetid.resize(setid.size());
  recvsetid_to_index.clear();
  for(int sindex=0;sindex<setid.size();sindex++){
    TypeRange tr = typerangearray[sindex];
    send_typerange[sindex] = tr;
    int size_of_set = tr.end-tr.begin;
    targetparticle.resize(targetparticle.size()+size_of_set);
    for(int i=0;i<size_of_set;i++){
      targetparticle[pr.end+i] = myparticle[tr.begin+i];
    }
    tr.shift(pr.end-tr.begin);
    pr.end += size_of_set;
    targettyperange[sindex] = tr;
    recvsetid[sindex] = setid[sindex];
    recvsetid_to_index.insert(std::pair<int,int>(recvsetid[sindex],sindex));
  }
  */
  send_typerange.resize(allsendsetid.size());
  targettyperange.resize(allsendsetid.size());
  recvsetid.resize(allsendsetid.size());
  recvsetid_to_index.clear();
  for(std::vector<int>::size_type si=0;si<allsendsetid.size();si++){
    int sindex=0;
    while(setid[sindex]!=allsendsetid[si]){sindex++;}
    TypeRange tr = typerangearray[sindex];
    send_typerange[si] = tr;
    int size_of_set = tr.end-tr.begin;
    targetparticle.resize(targetparticle.size()+size_of_set);
    for(int i=0;i<size_of_set;i++){
      //      targetparticle[pr.end+i] = myparticle[tr.begin+i];
      copy_particle_ghost(targetparticle,pr.end+i,myparticle,tr.begin+i);
    }
    tr.shift(pr.end-tr.begin);
    pr.end += size_of_set;
    targettyperange[si] = tr;
    recvsetid[si] = allsendsetid[si];
    recvsetid_to_index.insert(std::pair<int,int>(recvsetid[si],si));
  }
  target_range.push_back(pr);
  target++;
  targetbond = bondlistarray;
  
  //  std::cout << " particle_send_recv_one_axis x done " << std::endl;
  particle_send_recv_one_axis(recv_depth.x,recv_direction.x,0,1,target,
                              targetparticle,target_range,
                              recvsetid,recvsetid_to_index,
                              targettyperange,targetbond);
  //  std::cout << " particle_send_recv_one_axis x done " << std::endl;
  particle_send_recv_one_axis(recv_depth.y,recv_direction.y,2,3,target,
                              targetparticle,target_range,
                              recvsetid,recvsetid_to_index,
                              targettyperange,targetbond);

  particle_send_recv_one_axis(recv_depth.z,recv_direction.z,4,5,target,
                              targetparticle,target_range,
                              recvsetid,recvsetid_to_index,
                              targettyperange,targetbond);

  /*
  if(unit_id==0){
    std::cout << target << " target_range ";
    for(int t=0;t<target;t++){
      std::cout << "[" << target_range[t].begin << "," << target_range[t].end << ")";
    }
    std::cout << std::endl;
  }
  */
}

template<class PA, class GPA>
void 
MPICommunicator::exchangeParticleArraysubset_onlyPosition_top_half(const PA& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GPA& targetparticle,
                                             const std::vector<ParticleRange>& target_range,
                                             const std::vector<int>& recvsetid,
                                             const std::map<int,int>& recvsetid_to_index,
                                             const std::vector<TypeRange>& targettyperange,
                                             const std::vector<CovalentBondInfo::BondList>& targetbond)
{
  if(send_recv_communication_pattern==NearestXYZ){
    std::cout << "NearestXYZ not support overlap " << std::endl;
    particle_send_recv_onlyPosition_NearestXYZ(myparticle,
                                               typerangearray, bondlistarray,
                                   targetparticle, target_range, 
                                   recvsetid, recvsetid_to_index,
                                   targettyperange, targetbond);
  }else{
    particle_send_recv_onlyPosition_Direct_top_half(myparticle, typerangearray, bondlistarray,
                               targetparticle, target_range, 
                               recvsetid, recvsetid_to_index,
                               targettyperange, targetbond);
  }
}
template<class PA, class GPA>
void 
MPICommunicator::exchangeParticleArraysubset_onlyPosition_bottom_half(const PA& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GPA& targetparticle,
                                             const std::vector<ParticleRange>& target_range,
                                             const std::vector<int>& recvsetid,
                                             const std::map<int,int>& recvsetid_to_index,
                                             const std::vector<TypeRange>& targettyperange,
                                             const std::vector<CovalentBondInfo::BondList>& targetbond)
{
  if(send_recv_communication_pattern!=NearestXYZ){
    particle_send_recv_onlyPosition_Direct_bottom_half(myparticle, typerangearray, bondlistarray,
                               targetparticle, target_range, 
                               recvsetid, recvsetid_to_index,
                               targettyperange, targetbond);
    wait_complete_Particle();
  }
}
//! send own particle (with gap) and receive j-particle
/*!
  target_range
  Direct : particle ranges of all received particle, particle range of ghost
  NearestXYZ : particle range of self send particle and all received particle, particle range of ghost include a part of self particle which is ghost for others 
 */
template<class PA, class GPA>
void 
MPICommunicator::exchangeParticleArraysubset_onlyPosition(const PA& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GPA& targetparticle,
                                             const std::vector<ParticleRange>& target_range,
                                             const std::vector<int>& recvsetid,
                                             const std::map<int,int>& recvsetid_to_index,
                                             const std::vector<TypeRange>& targettyperange,
                                             const std::vector<CovalentBondInfo::BondList>& targetbond)
{
  if(send_recv_communication_pattern==NearestXYZ){
    //    std::cout << "particle_send_recv " << std::endl;
    particle_send_recv_onlyPosition_NearestXYZ(myparticle,
                                               typerangearray, bondlistarray,
                                   targetparticle, target_range, 
                                   recvsetid, recvsetid_to_index,
                                   targettyperange, targetbond);
    //    MPI_Barrier(mpi_comm_short);
    //    std::cout << "particle_send_recv done" << std::endl;
    
  }else{
    particle_send_recv_onlyPosition_Direct(myparticle, typerangearray, bondlistarray,
                               targetparticle, target_range, 
                               recvsetid, recvsetid_to_index,
                               targettyperange, targetbond);
  }
  if(send_recv_communication_pattern!=NearestXYZ){
    wait_complete_Particle();
  }
}

template
void 
MPICommunicator::exchangeParticleArraysubset_onlyPosition(const ParticleArray& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             ParticleArray& targetparticle,
                                             const std::vector<ParticleRange>& target_range,
                                             const std::vector<int>& recvsetid,
                                             const std::map<int,int>& recvsetid_to_index,
                                             const std::vector<TypeRange>& targettyperange,
                                                          const std::vector<CovalentBondInfo::BondList>& targetbond);
template
void 
MPICommunicator::exchangeParticleArraysubset_onlyPosition(const CombinedParticleArray& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GhostParticleArray& targetparticle,
                                             const std::vector<ParticleRange>& target_range,
                                             const std::vector<int>& recvsetid,
                                             const std::map<int,int>& recvsetid_to_index,
                                             const std::vector<TypeRange>& targettyperange,
                                                          const std::vector<CovalentBondInfo::BondList>& targetbond);
template
void 
MPICommunicator::exchangeParticleArraysubset_onlyPosition_top_half(const CombinedParticleArray& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GhostParticleArray& targetparticle,
                                             const std::vector<ParticleRange>& target_range,
                                             const std::vector<int>& recvsetid,
                                             const std::map<int,int>& recvsetid_to_index,
                                             const std::vector<TypeRange>& targettyperange,
                                                          const std::vector<CovalentBondInfo::BondList>& targetbond);
template
void 
MPICommunicator::exchangeParticleArraysubset_onlyPosition_bottom_half(const CombinedParticleArray& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GhostParticleArray& targetparticle,
                                             const std::vector<ParticleRange>& target_range,
                                             const std::vector<int>& recvsetid,
                                             const std::map<int,int>& recvsetid_to_index,
                                             const std::vector<TypeRange>& targettyperange,
                                                          const std::vector<CovalentBondInfo::BondList>& targetbond);


template<class PA, class GPA>
void 
MPICommunicator::exchangeParticleArraysubset_top_half(const PA& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GPA& targetparticle,
                                             std::vector<ParticleRange>& target_range,
                                             std::vector<int>& recvsetid,
                                             std::map<int,int>& recvsetid_to_index,
                                             std::vector<TypeRange>& targettyperange,
                                             std::vector<CovalentBondInfo::BondList>& targetbond)
{
  if(send_recv_communication_pattern==NearestXYZ){
    std::cout << "NearestXYZ not support overlap " << std::endl;
    particle_send_recv_NearestXYZ( myparticle, typerangearray, bondlistarray,
                                   targetparticle, target_range, 
                                   recvsetid, recvsetid_to_index,
                                   targettyperange, targetbond);
  }else{
    particle_send_recv_Direct_top_half( myparticle, typerangearray, bondlistarray,
                               targetparticle, target_range, 
                               recvsetid, recvsetid_to_index,
                               targettyperange, targetbond);
  }
}
template<class PA, class GPA>
void 
MPICommunicator::exchangeParticleArraysubset_bottom_half(const PA& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GPA& targetparticle,
                                             std::vector<ParticleRange>& target_range,
                                             std::vector<int>& recvsetid,
                                             std::map<int,int>& recvsetid_to_index,
                                             std::vector<TypeRange>& targettyperange,
                                             std::vector<CovalentBondInfo::BondList>& targetbond)
{
  if(send_recv_communication_pattern!=NearestXYZ){
    particle_send_recv_Direct_bottom_half( myparticle, typerangearray, bondlistarray,
                               targetparticle, target_range, 
                               recvsetid, recvsetid_to_index,
                               targettyperange, targetbond);
    wait_complete_Particle();
  }
}
template
void 
MPICommunicator::exchangeParticleArraysubset_top_half(const CombinedParticleArray& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GhostParticleArray& targetparticle,
                                             std::vector<ParticleRange>& target_range,
                                             std::vector<int>& recvsetid,
                                             std::map<int,int>& recvsetid_to_index,
                                             std::vector<TypeRange>& targettyperange,
						      std::vector<CovalentBondInfo::BondList>& targetbond);
template
void 
MPICommunicator::exchangeParticleArraysubset_bottom_half(const CombinedParticleArray& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GhostParticleArray& targetparticle,
                                             std::vector<ParticleRange>& target_range,
                                             std::vector<int>& recvsetid,
                                             std::map<int,int>& recvsetid_to_index,
                                             std::vector<TypeRange>& targettyperange,
							 std::vector<CovalentBondInfo::BondList>& targetbond);

//! send own particle (with gap) and receive j-particle
/*!
  target_range
  Direct : particle ranges of all received particle, particle range of ghost
  NearestXYZ : particle range of self send particle and all received particle, particle range of ghost include a part of self particle which is ghost for others 
 */
template<class PA, class GPA>
void 
MPICommunicator::exchangeParticleArraysubset(const PA& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GPA& targetparticle,
                                             std::vector<ParticleRange>& target_range,
                                             std::vector<int>& recvsetid,
                                             std::map<int,int>& recvsetid_to_index,
                                             std::vector<TypeRange>& targettyperange,
                                             std::vector<CovalentBondInfo::BondList>& targetbond)
{
  if(send_recv_communication_pattern==NearestXYZ){
    //    std::cout << "particle_send_recv " << std::endl;
    particle_send_recv_NearestXYZ( myparticle, typerangearray, bondlistarray,
                                   targetparticle, target_range, 
                                   recvsetid, recvsetid_to_index,
                                   targettyperange, targetbond);
    //    MPI_Barrier(mpi_comm_short);
    //    std::cout << "particle_send_recv done" << std::endl;
    
  }else{
    //    std::cout << " particle_send_recv_Direct  "  << std::endl;
    particle_send_recv_Direct( myparticle, typerangearray, bondlistarray,
                               targetparticle, target_range, 
                               recvsetid, recvsetid_to_index,
                               targettyperange, targetbond);
  }
  if(send_recv_communication_pattern!=NearestXYZ){
    //    std::cout << " wait_complete_Particle  "  << std::endl;
    wait_complete_Particle();
    //    std::cout << " wait_complete_Particle  done "  << std::endl;
  }
}
template
void 
MPICommunicator::exchangeParticleArraysubset(const ParticleArray& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             ParticleArray& targetparticle,
                                             std::vector<ParticleRange>& target_range,
                                             std::vector<int>& recvsetid,
                                             std::map<int,int>& recvsetid_to_index,
                                             std::vector<TypeRange>& targettyperange,
                                             std::vector<CovalentBondInfo::BondList>& targetbond);
template
void 
MPICommunicator::exchangeParticleArraysubset(const CombinedParticleArray& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             GhostParticleArray& targetparticle,
                                             std::vector<ParticleRange>& target_range,
                                             std::vector<int>& recvsetid,
                                             std::map<int,int>& recvsetid_to_index,
                                             std::vector<TypeRange>& targettyperange,
                                             std::vector<CovalentBondInfo::BondList>& targetbond);

//! send own particle (with gap) and receive j-particle with measurement MPI time
void 
MPICommunicator::exchangeParticleArraysubset(const ParticleArray& myparticle,
                                             const std::vector<TypeRange>& typerangearray,
                                             const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                             ParticleArray& targetparticle,
                                             std::vector<ParticleRange>& target_range,
                                             std::vector<int>& recvsetid,
                                             std::map<int,int>& recvsetid_to_index,
                                             std::vector<TypeRange>& targettyperange,
                                             std::vector<CovalentBondInfo::BondList>& targetbond,
                                             double& alltime, double& rectime)
{
  alltime -= MPI_Wtime();
  setSendParticlesubset(myparticle,typerangearray,bondlistarray);
  transferParticle();
  rectime -= MPI_Wtime();
  getReceiveParticle(targetparticle,target_range,targettyperange,targetbond,recvsetid,recvsetid_to_index);
  rectime += MPI_Wtime();
  wait_complete_Particle();
  alltime += MPI_Wtime();
}

void 
MPICommunicator::force_send_recv_Direct(const ForceArray& sendforce,
                                        const std::vector<ParticleRange>& send_range,
                                        ForceArray& recvforce)
{
  setSendForce(sendforce,send_range);
  transferForce();
  getReceiveForcesubset(recvforce);
  wait_complete_Force();
}

void 
MPICommunicator::force_send_recv_Direct_with_index(const ForceArray& sendforce,
                                                   const std::vector<ParticleRange>& send_range,
                                                   const std::vector<int>& forceindexset,
                                                   ForceArray& recvforce)
{
  setSendForce_with_index(sendforce,send_range,forceindexset);
  transferForce();
  getReceiveForcesubset_with_index(recvforce);
  wait_complete_Force();
}

void 
MPICommunicator::force_send_recv_Direct_indexed(const ForceArray& sendforce,
                                                const std::vector<ParticleRange>& send_range,
                                                const std::vector<int>& forceindexset,
                                                ForceArray& recvforce)
{
  setSendForce_indexed(sendforce,send_range,forceindexset);
  transferForce();
  getReceiveForcesubset_indexed(recvforce);
  wait_complete_Force();
}

void 
MPICommunicator::force_send_recv_one_axis(int depth, int recv_direction,
                                          int pf_plus, int pf_minus,
                                          int& target,
                                          ForceArray& sendforce,
                                          const std::vector<ParticleRange>& send_range,
                                          ForceArray& recvforce)
{
  for(int i=depth;i>=1;i--){
    if(recv_direction<+1){
      rpsf[pf_minus].stage--;
      sprf[pf_plus].stage--;
      rpsf[pf_minus].setsend(sendforce,send_range[target]);
      //      std::cout << " force_send_recv_one_axis pf_minus target " << pf_minus << " " << target << std::endl;
      target--;
    }
    if(recv_direction>-1){
      rpsf[pf_plus].stage--;
      sprf[pf_minus].stage--;
      rpsf[pf_plus].setsend(sendforce,send_range[target]);
      //      std::cout << " force_send_recv_one_axis pf_plus target " << pf_plus << " " << target << std::endl;
      target--;
    }
    // prepare receive from node {0,0,-1}
    if(recv_direction<+1){
      //      std::cout << " sprf " << pf_plus << " recv " << std::endl;
      sprf[pf_plus].prepare_receive();
    }
    if(recv_direction>-1){
      //      std::cout << " sprf " << pf_minus << " recv " << std::endl;
      sprf[pf_minus].prepare_receive();
    }
    // send to {0,0,+1}
    if(recv_direction<+1){
      //      std::cout << " rpsf " << pf_minus << " send " << std::endl;
      rpsf[pf_minus].send();
    }
    if(recv_direction>-1){
      //      std::cout << " rpsf " << pf_plus << " send " << std::endl;
      rpsf[pf_plus].send();
    }
    //! receive buff- from {0,0,-1}
    //! reduce buff+ to force to {[-X,+X],[-Y,+Y],+(k-1)}
    if(i>1){
      int recv_target=target;
      if(recv_direction<+1){
        sprf[pf_plus].getreceive(sendforce,send_range[recv_target].begin); 
        recv_target--;
      }
      if(recv_direction>-1){
        sprf[pf_minus].getreceive(sendforce,send_range[recv_target].begin); 
      }
    }else{
      if(recv_direction<+1){
        sprf[pf_plus].getreceivesubset(recvforce);
      }
      if(recv_direction>-1){
        sprf[pf_minus].getreceivesubset(recvforce);
      }
    }
    if(recv_direction<+1){
      rpsf[pf_minus].wait_complete_send();
    }
    if(recv_direction>-1){
      rpsf[pf_plus].wait_complete_send();
    }
  }
}

void 
MPICommunicator::force_send_recv_NearestXYZ(ForceArray& sendforce,
                                            const std::vector<ParticleRange>& send_range,
                                            ForceArray& recvforce)
{
  int target = 0;
  if(recv_direction.x==0){
    target += recv_depth.x*2;
  }else{
    target += recv_depth.x;
  }
  if(recv_direction.y==0){
    target += recv_depth.y*2;
  }else{
    target += recv_depth.y;
  }
  if(recv_direction.z==0){
    target += recv_depth.z*2;
  }else{
    target += recv_depth.z;
  }
  //  MPI_Barrier(mpi_comm_short);
    
  /*
  for(int i=send_range[0].begin;i<send_range[0].end;i++){
    sendforce[i].x=0.0;
    sendforce[i].y=0.0;
    sendforce[i].z=0.0;
  }
  */
  force_send_recv_one_axis(recv_depth.z,recv_direction.z,4,5,target,
                           sendforce,send_range,sendforce);
  force_send_recv_one_axis(recv_depth.y,recv_direction.y,2,3,target,
                           sendforce,send_range,sendforce);
  force_send_recv_one_axis(recv_depth.x,recv_direction.x,0,1,target,
                           sendforce,send_range,sendforce);
  int si=0;
  for(std::vector<TypeRange>::size_type sindex=0;
      sindex<send_typerange.size();sindex++){
    for(int i=send_typerange[sindex].begin;i<send_typerange[sindex].end;i++){
      recvforce[i] += sendforce[si];
      si++;
    }
  }
  
}

//! send j-force (with gap) and receive own force
/*!
  called by Integrator::time_integrations
  communicator.exchangeForceArraysubset(ghostforce, targetparticlerange, force);
  targetparticlerange is 5th argument for exchangeParticleArraysubset : target_range in exchangeParticleArraysubset
 */
void 
MPICommunicator::exchangeForceArraysubset(ForceArray& sendforce,
                                          std::vector<ParticleRange>& send_range,
                                          ForceArray& recvforce)
{
  if(send_recv_communication_pattern==NearestXYZ){
    //    std::cout << "force_send_recv " << std::endl;
    force_send_recv_NearestXYZ(sendforce, send_range, recvforce);
    //    std::cout << "force_send_recv done " << std::endl;
  }else{
    force_send_recv_Direct(sendforce, send_range, recvforce);
  }

  //   std::cout << "wait_complete_Force" << std::endl;
}
void 
MPICommunicator::exchangeForceArraysubset_with_index(ForceArray& sendforce,
                                                     std::vector<ParticleRange>& send_range,
                                                     const std::vector<int>& forceindexset,
                                                     ForceArray& recvforce)
{
  if(send_recv_communication_pattern==NearestXYZ){
    //    std::cout << "force_send_recv " << std::endl;
    force_send_recv_NearestXYZ(sendforce, send_range, recvforce);
    //    std::cout << "force_send_recv done " << std::endl;
  }else{
    force_send_recv_Direct_with_index(sendforce, send_range, forceindexset, recvforce);
  }

  //   std::cout << "wait_complete_Force" << std::endl;
}
void 
MPICommunicator::exchangeForceArraysubset_indexed(ForceArray& sendforce,
                                                  std::vector<ParticleRange>& send_range,
                                                  const std::vector<int>& forceindexset,
                                                  ForceArray& recvforce)
{
  if(send_recv_communication_pattern==NearestXYZ){
    //    std::cout << "force_send_recv " << std::endl;
    force_send_recv_NearestXYZ(sendforce, send_range, recvforce);
    //    std::cout << "force_send_recv done " << std::endl;
  }else{
    force_send_recv_Direct_indexed(sendforce, send_range, forceindexset, recvforce);
  }

  //   std::cout << "wait_complete_Force" << std::endl;
}
 
 
//! send j-force (with gap) and receive own force with measurement MPI time
void 
MPICommunicator::exchangeForceArraysubset(ForceArray& sendforce,
                                          std::vector<ParticleRange>& send_range,
                                          ForceArray& recvforce,
                                          double& alltime, double& rectime)
{
  alltime -= MPI_Wtime();
  setSendForce(sendforce,send_range);
  transferForce();
  rectime -= MPI_Wtime();
  getReceiveForcesubset(recvforce);
  rectime += MPI_Wtime();
  wait_complete_Force();
  alltime += MPI_Wtime();
}

//! set move out particle
void 
MPICommunicator::setSendMove()
{
  for(int t=0;t<move_number_of_target;t++){
    move_out[t].setmoveoutparticle(move_out_particle[t],move_out_type[t],move_out_cellid[t],move_out_bondlistarray[t]);
  }
}
   
//! request transfer Move particle
void 
MPICommunicator::transferMove()
{
  for(int t=0;t<move_number_of_target;t++){
    move_in[t].prepare_receive();
  }
  for(int t=0;t<move_number_of_target;t++){
    move_out[t].send();
  }
}

//! get received Move in particle
void 
MPICommunicator::getReceiveMove()
{
  for(int t=0;t<move_number_of_target;t++){
    move_in[t].getreceive(move_in_particle,move_in_type,move_in_cellid,move_in_bondlistarray);
  }
}
  
//! wait/check move out transfer
void 
MPICommunicator::wait_complete_Move()
{
  for(int t=0;t<move_number_of_target;t++){
    move_out[t].wait_complete_send();
  }
}

template<> void
MPICommunicator::move_send_recv<Direct>()
{
  setSendMove();
  transferMove();
  move_in_particle.clear();
  move_in_cellid.clear();
  move_in_type.clear();
  move_in_bondlistarray.clear();
  getReceiveMove();
}

template<> void
MPICommunicator::distribute_move_out<Direct>()
{
  for(size_t i=0;i<all_move_out_particle.size();i++){
    int tset = all_move_out_cellid[i];
    int target = setid_to_move_target_index[tset];
    move_out_particle[target].push_back(all_move_out_particle[i]);
    move_out_cellid[target].push_back(all_move_out_cellid[i]);
    move_out_type[target].push_back(all_move_out_type[i]);
    move_out_bondlistarray[target].push_back(all_move_out_bondlistarray[i]);
  }
}

template<> void
MPICommunicator::distribute_move_out<NearestXYZ>()
{
  /*
  if(all_move_out_particle.size()>0){
        std::cout << unit_id << ":"<< node_position << " found " << all_move_out_particle.size() << " move_out " << std::endl;
  }
  */
  //  std::cout << unit_id << ": all_move_out_bondlistarray.size()" << all_move_out_bondlistarray.size() << std::endl;
  for(size_t i=0;i<all_move_out_particle.size();i++){
    Particle out_particle = all_move_out_particle[i];
    int outset;
    int diffx = node_geometry.cellid_to_relative_x(all_move_out_cellid[i],node_position.x);
    if(diffx>0){         // for +x 
      if(diffx>move_depth.x){ // to far ?
        diffx = move_depth.x;
      }
      outset = diffx;
    }else if(diffx<0){      // for -x
      if(diffx<-move_depth.x){ // to far ?
        diffx = -move_depth.x;
      }
      outset = move_depth.x - diffx;
    }else{
      int diffy = node_geometry.cellid_to_relative_y(all_move_out_cellid[i],node_position.y);
      if(diffy>0){   // for +y
        if(diffy>move_depth.y){ // to far ?
          diffy = move_depth.y;
        }
        outset = move_depth.x*2 + diffy;
      }else if(diffy<0){ // for -y
        if(diffy<-move_depth.y){ // to far ?
          diffy = -move_depth.y;
        }
        outset = move_depth.x*2 + move_depth.y - diffy;
      }else{          // for this or +-z
        int diffz = node_geometry.cellid_to_relative_z(all_move_out_cellid[i],node_position.z);
        if(diffz>0){   // for +z
          if(diffz>move_depth.z){ // to far ?
            diffz = move_depth.z;
          }
          outset = move_depth.x*2 + move_depth.y*2 + diffz;
        }else if(diffz<0){ // for -z
          if(diffz<-move_depth.z){ // to far ?
            diffz = -move_depth.z;
          }
          outset = move_depth.x*2 + move_depth.y*2 + move_depth.z - diffz;
        }else{          // for this ??
          //      std::cout << unit_id << ": found " << out_particle.position << " move_out to self? " << diffx << " " << diffy << " " << diffz << std::endl;
          move_in_particle.push_back(out_particle);
          move_in_type.push_back(all_move_out_type[i]);
          move_in_cellid.push_back(all_move_out_cellid[i]);
          move_in_bondlistarray.push_back(all_move_out_bondlistarray[i]);
          continue;
        }
      }
    }
    outset--;
    //    std::cout << unit_id << ": found " << out_particle.position << " move_out to " <<  outset << ":" << move_out_target_id[outset]  << std::endl;
    move_out_particle[outset].push_back(out_particle);
    move_out_type[outset].push_back(all_move_out_type[i]);
    move_out_cellid[outset].push_back(all_move_out_cellid[i]);
    move_out_bondlistarray[outset].push_back(all_move_out_bondlistarray[i]);
  }
}

/*
  This is for move_out/in with cellid
  not for current move_out/in
*/
void
MPICommunicator::marge_move(ParticleArray& particle, 
                            std::vector<PotentialModel> &rangetype,
                            std::vector<int> &cellid,
                            std::vector<CovalentBondInfo::BondList>& bondlist, 
                            ParticleArray& add_particle, 
                            std::vector<PotentialModel> &add_rangetype,
                            std::vector<int> &add_cellid,
                            std::vector<CovalentBondInfo::BondList>& add_bondlist)
{
  size_t org_size = particle.size();
  size_t marge_size = particle.size() + add_particle.size();

  particle.resize(marge_size);
  cellid.resize(marge_size);
  rangetype.resize(marge_size);
  bondlist.resize(marge_size);

  for(ParticleArray::size_type i=0;i<add_particle.size();i++){
    particle[org_size+i] = add_particle[i];
    rangetype[org_size+i] = add_rangetype[i];
    cellid[org_size+i] = add_cellid[i];
    bondlist[org_size+i] = add_bondlist[i];
  }
}

/*
  current move_out/in not send/receive cellid
*/
void
MPICommunicator::marge_move(ParticleArray& particle, 
                            std::vector<PotentialModel> &rangetype,
                            ParticleArray& add_particle, 
                            std::vector<PotentialModel> &add_rangetype)
{
  size_t org_size = particle.size();
  size_t marge_size = particle.size() + add_particle.size();

  particle.resize(marge_size);
  rangetype.resize(marge_size);

  for(ParticleArray::size_type i=0;i<add_particle.size();i++){
    particle[org_size+i] = add_particle[i];
    rangetype[org_size+i] = add_rangetype[i];
  }
}

/*
  sequence of move NearestXYZ
  N_x, N_y, N_z : move range

  Stage   Send              Recv
  X 0   {0}->+-N to +-X  {+-1}->-+(N-1) from +-X
          marge {0}->+-(N-1) {-+1}->+-(N-1)
  X 1   {-+1,0}->+-(N-1) to +-X  {+-2,+-1}->-+(N-2) from +-X
  ...
          marge {0}->+-(N-i) {-+i,-+1}->+-(N-i)
  X i     {-+i,0}->+-(N-i) to +-X

  X N-1 {-+(N-1),0}->+-1 to +-X  {+-N,+-1}->0 from +-X 
 */

/*
  order of move_in, move_out
  +x,-x,+y,-y,+z,-z

  order move_out_particle
  2N_x,2N_y,2N_z : 0,1,...N_x-1,N_x,... for (1,*,*),(2,*,*),...(N_x,*,*),(-1,*,*),...
 */
template<> void
MPICommunicator::move_send_recv<NearestXYZ>()
{

  move_in_particle.clear();
  move_in_cellid.clear();
  move_in_type.clear();
  move_in_bondlistarray.clear();

  // x-direction stage
  // in this stage data for x contain data for (x,[-N_y,N_y],[-N_z,N_z])
  for(int x=move_depth.x-1;x>=0;x--){
    int i=x;
    int j=move_depth.x+x;
    if(x<move_depth.x-1){
      // marge received data from -x node at previous stage to moveout data for move_depth.x-i
      marge_move(move_out_particle[i], move_out_type[i], 
                 move_out_cellid[i], move_out_bondlistarray[i],
                 move_imd_particle[1], move_imd_type[1], 
                 move_imd_cellid[1], move_imd_bondlistarray[1]);
      // marge received data from +x node at previous stage to moveout data for -(move_depth.x-i)
      marge_move(move_out_particle[j], move_out_type[j], 
                 move_out_cellid[j], move_out_bondlistarray[j],
                 move_imd_particle[0], move_imd_type[0], 
                 move_imd_cellid[0], move_imd_bondlistarray[0]);
    } 
    // prepapre receive from -/+x for +/-(move_depth.x-1-i)
    move_in[1].prepare_receive();
    move_in[0].prepare_receive();
    // send to +/-x for +/-(move_depth.x-i)
    move_out[0].setmoveoutparticle(move_out_particle[i],
                                   move_out_type[i],
                                   move_out_cellid[i],
                                   move_out_bondlistarray[i]);
    move_out[1].setmoveoutparticle(move_out_particle[j],
                                   move_out_type[j],
                                   move_out_cellid[j],
                                   move_out_bondlistarray[j]);
    /*
    if(move_out_particle[i].size()>0){
      std::cout << unit_id << ": send " << move_out_particle[i].size() << " particle to " << move_out_target_id[i];
    }
    if(move_out_particle[j].size()>0){
      std::cout << unit_id << ": send " << move_out_particle[j].size() << " particle to " << move_out_target_id[j];
    }
    //    if(unit_id==0)
      {
      if(move_out_particle[i].size()>0){
        std::cout << unit_id << ": send " << move_out_particle[i].size() << " particle to " << move_out_target_id[i];
        for(int p=0;p<move_out_particle[i].size();p++){
          std::cout << " " << move_out_cellid[i][p];
        }
        std::cout << std::endl;
      }
    }
    */
    move_out[0].send();
    move_out[1].send();
    // receive  from -/+x for +/-(move_depth.x-1-i)
    move_imd_particle[1].clear();
    move_imd_type[1].clear();
    move_imd_cellid[1].clear();
    move_imd_bondlistarray[1].clear();
    move_imd_particle[0].clear();
    move_imd_type[0].clear();
    move_imd_cellid[0].clear();
    move_imd_bondlistarray[0].clear();
    move_in[1].getreceive(move_imd_particle[1],move_imd_type[1],
                          move_imd_cellid[1],move_imd_bondlistarray[1]);
    move_in[0].getreceive(move_imd_particle[0],move_imd_type[0],
                          move_imd_cellid[0],move_imd_bondlistarray[0]);
    /*
    //    if(unit_id==0)
{
      if(move_imd_particle[1].size()>0){
        std::cout << unit_id << ": recv " << move_imd_particle[1].size() << " : ";
        for(int p=0;p<move_imd_particle[1].size();p++){
          std::cout << " " << move_imd_cellid[1][p];
        }
        std::cout << std::endl;
      }
    }
    */
    // wite_send
    move_out[0].wait_complete_send();
    move_out[1].wait_complete_send();
  }

  // marge 0 to y-direction and data for y-direction in last received x-dirction
  for(int dir=0;dir<2;dir++){
    for(std::vector<ParticleArray>::size_type i=0;
        i<move_imd_particle[dir].size();i++){
      Particle imd_particle = move_imd_particle[dir][i];
      int diffy = node_geometry.cellid_to_relative_y(move_imd_cellid[dir][i],node_position.y);
      int outset;
      if(diffy>0){   // for +y
        if(diffy>move_depth.y){ // to far ?
          diffy = move_depth.y;
        }
        outset = move_depth.x*2 + diffy;
      }else if(diffy<0){ // for -y
        if(diffy<-move_depth.y){ // to far ?
          diffy = -move_depth.y;
        }
        outset = move_depth.x*2 + move_depth.y - diffy;
      }else{          // for this or +-z
        int diffz = node_geometry.cellid_to_relative_z(move_imd_cellid[dir][i],node_position.z);
        if(diffz>0){   // for +z
          if(diffz>move_depth.z){ // to far ?
            diffz = move_depth.z;
          }
          outset = move_depth.x*2 + move_depth.y*2 + diffz;
        }else if(diffz<0){ // for -z
          if(diffz<-move_depth.z){ // to far ?
            diffz = -move_depth.z;
          }
          outset = move_depth.x*2 + move_depth.y*2 + move_depth.z - diffz;
        }else{          // for this 
          move_in_particle.push_back(imd_particle);
          move_in_type.push_back(move_imd_type[dir][i]);
          move_in_cellid.push_back(move_imd_cellid[dir][i]);
          move_in_bondlistarray.push_back(move_imd_bondlistarray[dir][i]);
          continue;
        }
      }
      outset--;
      //      std::cout << unit_id << ": found " << imd_particle.position << " pass through " << outset << ":" << move_out_target_id[outset] <<  std::endl;
      move_out_particle[outset].push_back(imd_particle);
      move_out_type[outset].push_back(move_imd_type[dir][i]);
      move_out_cellid[outset].push_back(move_imd_cellid[dir][i]);
      move_out_bondlistarray[outset].push_back(move_imd_bondlistarray[dir][i]);
    }
  }

  // y-direction stage
  // in this stage data for y contain data for (0,y,[-N_z,N_z])
  for(int y=move_depth.y-1;y>=0;y--){
    int i = move_depth.x*2+y;
    int j = move_depth.x*2+move_depth.y+y;
    if(y<move_depth.y-1){
      // marge received data from -y node at previous stage to moveout data for move_depth.y-i
      marge_move(move_out_particle[i], move_out_type[i], 
                 move_out_cellid[i], move_out_bondlistarray[i],
                 move_imd_particle[3], move_imd_type[3], 
                 move_imd_cellid[3], move_imd_bondlistarray[3]);
      // marge received data from +y node at previous stage to moveout data for -(move_depth.y-i)
      marge_move(move_out_particle[j], move_out_type[j], 
                 move_out_cellid[j], move_out_bondlistarray[j],
                 move_imd_particle[2], move_imd_type[2], 
                 move_imd_cellid[2], move_imd_bondlistarray[2]);
    } 
    // prepapre receive from -/+y for +/-(move_depth.y-1-i)
    move_in[3].prepare_receive();
    move_in[2].prepare_receive();
    // send to +/-y for +/-(move_depth.y-i)
    move_out[2].setmoveoutparticle(move_out_particle[i],
                                   move_out_type[i],
                                   move_out_cellid[i],
                                   move_out_bondlistarray[i]);
    move_out[3].setmoveoutparticle(move_out_particle[j],
                                   move_out_type[j],
                                   move_out_cellid[j],
                                   move_out_bondlistarray[j]);
    move_out[2].send();
    move_out[3].send();
    // receive  from -/+y for +/-(move_depth.y-1-i)
    move_imd_particle[3].clear();
    move_imd_type[3].clear();
    move_imd_cellid[3].clear();
    move_imd_bondlistarray[3].clear();
    move_imd_particle[2].clear();
    move_imd_type[2].clear();
    move_imd_cellid[2].clear();
    move_imd_bondlistarray[2].clear();
    move_in[3].getreceive(move_imd_particle[3],move_imd_type[3],
                          move_imd_cellid[3],move_imd_bondlistarray[3]);
    move_in[2].getreceive(move_imd_particle[2],move_imd_type[2],
                          move_imd_cellid[2],move_imd_bondlistarray[2]);
    // wite_send
    move_out[2].wait_complete_send();
    move_out[3].wait_complete_send();
  }

  // marge 0 to z-direction and data for z-direction in last received z-dirction
  for(int dir=2;dir<4;dir++){
    for(ParticleArray::size_type i=0;i<move_imd_particle[dir].size();i++){
      Particle imd_particle = move_imd_particle[dir][i];
      int diffz = node_geometry.cellid_to_relative_z(move_imd_cellid[dir][i],node_position.z);
      int outset;
      if(diffz>0){   // for +z
        if(diffz>move_depth.z){ // to far ?
          diffz = move_depth.z;
        }
        outset = move_depth.x*2 + move_depth.y*2 + diffz;
      }else if(diffz<0){ // for -z
        if(diffz<-move_depth.z){ // to far ?
          diffz = -move_depth.z;
        }
        outset = move_depth.x*2 + move_depth.y*2 + move_depth.z - diffz;
      }else{          // for this 
        move_in_particle.push_back(imd_particle);
        move_in_type.push_back(move_imd_type[dir][i]);
        move_in_cellid.push_back(move_imd_cellid[dir][i]);
        move_in_bondlistarray.push_back(move_imd_bondlistarray[dir][i]);
        continue;
      }
      outset--;
      move_out_particle[outset].push_back(imd_particle);
      move_out_type[outset].push_back(move_imd_type[dir][i]);
      move_out_cellid[outset].push_back(move_imd_cellid[dir][i]);
      move_out_bondlistarray[outset].push_back(move_imd_bondlistarray[dir][i]);
    }
  }

  // z-direction stage
  // in this stage data for z contain data for (0,0,z)
  for(int z=move_depth.z-1;z>=0;z--){
    int i = move_depth.x*2+move_depth.y*2+z;
    int j = move_depth.x*2+move_depth.y*2+move_depth.z+z;
    if(z<move_depth.z-1){
      // marge received data from -z node at previous stage to moveout data for move_depth.z-i
      marge_move(move_out_particle[i], move_out_type[i],
                 move_out_cellid[i], move_out_bondlistarray[i],
                 move_imd_particle[5], move_imd_type[5],
                 move_imd_cellid[5], move_imd_bondlistarray[5]);
      // marge received data from +z node at previous stage to moveout data for -(move_depth.z-i)
      marge_move(move_out_particle[j], move_out_type[j],
                 move_out_cellid[j], move_out_bondlistarray[j],
                 move_imd_particle[4], move_imd_type[4], 
                 move_imd_cellid[4], move_imd_bondlistarray[4]);
    } 
    // prepapre receive from -/+z for +/-(move_depth.z-1-i)
    move_in[5].prepare_receive();
    move_in[4].prepare_receive();
    // send to +/-z for +/-(move_depth.z-i)
    move_out[4].setmoveoutparticle(move_out_particle[i],
                                   move_out_type[i],
                                   move_out_cellid[i], 
                                   move_out_bondlistarray[i]);
    move_out[5].setmoveoutparticle(move_out_particle[j],
                                   move_out_type[j],
                                   move_out_cellid[j], 
                                   move_out_bondlistarray[j]);
    move_out[4].send();
    move_out[5].send();
    // receive  from -/+z for +/-(move_depth.z-1-i)
    move_imd_particle[5].clear();
    move_imd_type[5].clear();
    move_imd_cellid[5].clear();
    move_imd_bondlistarray[5].clear();
    move_imd_particle[4].clear();
    move_imd_type[4].clear();
    move_imd_cellid[4].clear();
    move_imd_bondlistarray[4].clear();
    move_in[5].getreceive(move_imd_particle[5],move_imd_type[5],
                          move_imd_cellid[5],move_imd_bondlistarray[5]);
    move_in[4].getreceive(move_imd_particle[4],move_imd_type[4],
                          move_imd_cellid[4],move_imd_bondlistarray[4]);
    // wite_send
    move_out[4].wait_complete_send();
    move_out[5].wait_complete_send();
  }

  marge_move(move_in_particle, move_in_type, 
             move_in_cellid, move_in_bondlistarray,
             move_imd_particle[5], move_imd_type[5], 
             move_imd_cellid[5], move_imd_bondlistarray[5]);
  marge_move(move_in_particle, move_in_type, 
             move_in_cellid, move_in_bondlistarray,
             move_imd_particle[4], move_imd_type[4], 
             move_imd_cellid[4], move_imd_bondlistarray[4]);
}

//! send move out particle and receive move in particle
template<typename PA>
void 
MPICommunicator::move_particle(PA& particlearray,
                               std::vector<TypeRange>& typerangearray)
{
  // select move particle
#ifdef TIMER_DETAIL
  PerfCounter::start(timer_move);
#endif
  for(int t=0;t<move_number_of_target;t++){
    move_out_particle[t].clear();
    move_out_cellid[t].clear();
    move_out_type[t].clear();
    move_out_bondlistarray[t].clear();
  }
  // divide move out particles to target
  if(move_communication_pattern==Direct){
    distribute_move_out<Direct>();
  }else if(move_communication_pattern==NearestXYZ){
    distribute_move_out<NearestXYZ>();
  }else{
    std::cout << " not implement move_communication_pattern " << move_communication_pattern << std::endl;  
  }
#ifdef PROF_CHECK
# ifdef PROF_CHECK_VERB
    if(unit_id==0){
      printf("MPI_Barrier before move_send_recv\n");
    }
# endif
    if(mpi_comm_short!=MPI_COMM_NULL){
      MPI_Barrier(mpi_comm_short);
    }
#endif

#ifdef TIMER_DETAIL
  PerfCounter::start(timer_move_comm);
#endif
  if(move_communication_pattern==Direct){
    move_send_recv<Direct>();
  }else if(move_communication_pattern==NearestXYZ){
    move_send_recv<NearestXYZ>();
  }else{
    std::cout << " not implement move_communication_pattern " << move_communication_pattern << std::endl;
  }
  wait_complete_Move();
#ifdef TIMER_DETAIL
  PerfCounter::stop();
  PerfCounter::stop();
#endif
}
template
void 
MPICommunicator::move_particle(ParticleArray& particlearray,
                               std::vector<TypeRange>& typerangearray);
template
void 
MPICommunicator::move_particle(CombinedParticleArray& particlearray,
                               std::vector<TypeRange>& typerangearray);
/*
  NearestXYZ MPISendParticleReceiveForce/MPIReceiveParticleSendForce

  target {[-X,+X],[-Y,+Y],[mZ,+Z]}   : mZ = -Z for full shell, 0 for half shell

  number_of_target 6 for full, 5 for half
  target of sprf :  0:+X 1:-X 2:+Y 3:-Y 4:+Z(used only full) 5:-Z
  target of rpsf :  0:+X 1:-X 2:+Y 3:-Y 4:+Z 5:-Z(used only full)

        recv -1   recv +1   send +1  send -1
  X1 :  {-1,0,0}  {+1,0,0}  {0,0,0}  {0,0,0}
  X2 :  {-2,0,0}  {+2,0,0}  {-1,0,0} {+1,0,0}
...
  Xi :  {-i,0,0}  {+i,0,0}  {-(i-1),0,0} {+(i-1),0,0}
...
  XX :  {-X,0,0}  {+X,0,0}  {-(X-1),0,0} {+(X-1),0,0}
  
  Y1 :  {[-X,+X],-1,0} {[-X,+X],+1,0} {[-X,+X],0,0} {[-X,+X],0,0} 
  Y2 :  {[-X,+X],-2,0} {[-X,+X],+2,0} {[-X,+X],-1,0} {[-X,+X],+1,0} 
...
  Yj :  {[-X,+X],-j,0} {[-X,+X],+j,0} {[-X,+X],-(j-1),0} {[-X,+X],+(j-1),0} 
...
  YY :  {[-X,+X],-Y,0} {[-X,+X],+Y,0} {[-X,+X],-(Y-1),0} {[-X,+X],+(Y-1),0}

full shell/cube
  Z1 :  {[],[],-1} {[],[],+1} {[],[],0} {[],[],0}
  Z2 :  {[],[],-2} {[],[],+2} {[],[],-1} {[],[],+1}
...
  Zk :  {[],[],-k} {[],[],+k} {[],[],-(k-1)} {[],[],+(k-1)}
...
  ZZ :  {[],[],-Z} {[],[],+Z} {[],[],-(Z-1)} {[],[],+(Z-1)}

half shell/cube
  Z1 :  { - } {[],[],+1} { - } {[],[],0}
  Z2 :  { - } {[],[],+2} { - } {[],[],+1}
...
  Zk :  { - } {[],[],+k} { - } {[],[],+(k-1)}
...
  ZZ :  { - } {[],[],+Z} { - } {[],[],+(Z-1)}


  Xi stage receive cell {+/-i,0,0} from node {+/-1,0,0} 
           send    cell {+/-(i-1),0,0} to node {-/+1,0,0} 
                   = receieved cell at X(i-1) stage i>0
                   = my node cell                   i=1

  Yj stage receive cell {[-X,+X],+/-j,0} from node {0,+/-j,0}
           send    cell {[-X,+X],+/-(j-1),0} to node {0,-/+j,0}
                   = receieved cell at Y(j-1) stage j>0
                     receieved cell at XX stage j=1

  Zk stage receive cell {[-X,+X],[-Y,+Y],+{/-}k} from node {0,0,+{/-}j}
           send    cell {[-X,+X],[-Y,+Y],+{/-}(k-1)} to node {0,0,-{/+}j}
                   = receieved cell at Y(k-1) stage k>0
                     receieved cell at YY stage k=1

  clear all_receive
  append my node cell to all_receive  // send particle set packed dense same as ghost
  for(i=1;i<=X;i++){
    if(i==1){
      // set all_receive to send to {-1,0,0}
      sprf[1].setsend(all_receive);
      // set all_receive to send to {+1,0,0}
      sprf[0].setsend(all_receive);
    }else{
      // set buff+ to send to {-1,0,0}
      copy rpsf[0].recv_buffer to sprf[1].send_buffer
      // set buff- to send to {+1,0,0}
      copy rpsf[1].recv_buffer to sprf[0].send_buffer
    }
    // prepare receive from node {+1,0,0}
    rpsf[0].prepare_receive();
    // prepare receive from node {-1,0,0}
    rpsf[1].prepare_receive();
    // send to {-1,0,0}
    sprf[1].send();
    // send to {+1,0,0}
    sprf[0].send();
    // receive buff+ from {+1,0,0}
    // receive buff- from {-1,0,0}
    // append buff+ to all_receive
    // append buff- to all_receive
    rpsf[0].getreceive(all_receive);
    rpsf[1].getreceive(all_receive);
  }

  for(j=1;j<=Y;j++){
    if(j==1){
      // set all_receive to send to {0,-1,0}
      sprf[3].setsend(all_receive);
      // set all_receive to send to {0,+1,0}
      sprf[2].setsend(all_receive);
    }else{
      // set buff+ to send to {0,-1,0}
      copy rpsf[2].recv_buffer to sprf[3].send_buffer
      // set buff- to send to {0,+1,0}
      copy rpsf[3].recv_buffer to sprf[2].send_buffer
    }
    // prepare receive from node {0,+1,0}
    rpsf[2].prepare_receive();
    // prepare receive from node {0,-1,0}
    rpsf[3].prepare_receive();
    //send to {0,-1,0}
    sprf[3].send();
    //send to {0,+1,0}
    sprf[2].send();
    // receive buff+ from {0,+1,0}
    // receive buff- from {0,-1,0}
    // append buff+ to all_receive
    // append buff- to all_receive
    rpsf[2].getreceive(all_receive);
    rpsf[3].getreceive(all_receive);
  }

  for(k=1;k<=Z;k++){
    if(k==1){
      // set all_receive to send to {0,0,-1}
      sprf[5].setsend(all_receive);
      // if(full)set all_receive to send to {0,0,+1}
      if(full)sprf[4].setsend(all_receive);
    }else{
      // set buff+ to send to {0,0,-1}
      copy rpsf[4].recv_buffer to sprf[5].send_buffer
      //if(full)set buff- to send to {0,0,+1}
      if(full)copy rpsf[5].recv_buffer to sprf[4].send_buffer
    }
    // prepare receive from node {0,0,+1}
    rpsf[4].prepare_receive();
    //if(full)prepare receive from node {0,0,-1}
    if(full)rpsf[5].prepare_receive();
    // send to {0,0,-1}
    sprf[5].send();
    // if(full)send to {0,0,+1}
    sprf[4].send();
    //receive buff+ from {0,0,+1}
    //if(full)receive buff- from {0,0,-1}
    //append buff+ to all_receive
    //if(full)merge buff- to all_receive
    rpsf[4].getreceive(all_receive);
    if(full)rpsf[5].getreceive(all_receive);
  }

  all_receive may be packed 
  ( 0, 0, 0),(+1, 0, 0),(-1, 0, 0),(+2, 0, 0),(-2, 0, 0),...,(-X, 0, 0),  2X   (0,0,0) is not receive, copy from my cell
  ( 0,+1, 0),(+1,+1, 0),(-1,+1, 0),  ...                    ,(-X,+1, 0),  2X+1\
  ( 0,-1, 0),(+1,-1, 0),(-1,-1, 0),                         ,(-X,-1, 0),  2X+1|
  ...                                                                         |2Y
  ( 0,+Y, 0),(+1,+Y, 0),(-1,+Y, 0),                         ,(-X,+Y, 0),  2X+1|
  ( 0,-Y, 0),(+1,-Y, 0),(-1,-Y, 0),                         ,(-X,-Y, 0),  2X+1/
  ( 0, 0,+1),(+1, 0,+1),(-1, 0,+1),                         ,(-X, 0,+1),  ----\-------------\
  ( 0,+1,+1),(+1,+1,+1),(-1,+1,+1),                         ,(-X,+1,+1),      |             |
  ...                                                                         |(2X+1)*(2Y+1)|
  ( 0,-Y,+1),(+1,-Y,+1),(-1,-Y,+1),                         ,(-X,-Y,+1),  ----/             |
 [( 0, 0,-1),(+1, 0,-1),(-1, 0,-1),                         ,(-X, 0,-1),                    |
  ( 0,+1,-1),(+1,+1,-1),(-1,+1,-1),                         ,(-X,+1,-1),                    |
  ...                                                                                       |
  ( 0,-Y,-1),(+1,-Y,-1),(-1,-Y,-1),                         ,(-X,-Y,-1),]                   |
  ...                                                                                       |
  ( 0, 0,+Z),(+1, 0,+Z),(-1, 0,+Z),                         ,(-X, 0,+Z),                    |[2]Z
  ( 0,+1,+Z),(+1,+1,+Z),(-1,+1,+Z),                         ,(-X,+1,+Z),                    |
  ...                                                                                       |
  ( 0,-Y,+Z),(+1,-Y,+Z),(-1,-Y,+Z),                         ,(-X,-Y,+Z),                    |
 [( 0, 0,-Z),(+1, 0,-Z),(-1, 0,-Z),                         ,(-X, 0,-Z),                    |
  ( 0,+1,-Z),(+1,+1,-Z),(-1,+1,-Z),                         ,(-X,+1,-Z),                    |
  ...                                                                                       |
  ( 0,-Y,-Z),(+1,-Y,-Z),(-1,-Y,-Z),                         ,(-X,-Y,-Z),]     --------------/


  all_receive except (0,0,0) used as ghost


  return force phase only for half cube/shell
  reduce mean add force to same particle
  
  // sprf[5] MPI recv buffer, rpsf[4] MPI send buffer :  (2X+1)*(2Y+1) set of force
  for(k=Z;k=>1;k--){
    // set force to {[-X,+X],[-Y,+Y],+k}     
    rpsf[4].setsend(ghostforce(z=k));
    // prepare receive from node {0,0,-1}
    sprf[5].prepare_receive();
    // send to {0,0,+1}
    rpsf[4].send();
    //! receive buff- from {0,0,-1}
    //! reduce buff+ to force to {[-X,+X],[-Y,+Y],+(k-1)}  
    sprf[5].getreceive(ghostforce(z=-(k-1))); 
  }

  for(j=Y;j>=1;j--){
    // set force to {[-X,+X],+j,0}
    rpsf[2].setsend(ghostforce(y=j,z=0));
    // set force to {[-X,+X],-j,0}
    rpsf[3].setsend(ghostforce(y=-j,z=0));
    // prepare receive from node {0,-1,0}
    sprf[3].prepare_receive();
    // prepare receive from node {0,+1,0}
    sprf[2].prepare_receive();
    // send to {0,+1,0}
    rpsf[2].send();
    // send to {0,-1,0}
    rpsf[3].send();
    //receive buff- from {0,-1,0}
    //receive buff+ from {0,+1,0}
    //reduce buff- to force to {[-X,+X],+(j-1),0}
    //reduce buff+ to force to {[-X,+X],-(j-1),0}
    sprf[3].getreceive(ghostforce(y=-(j-1),z=0));
    sprf[2].getreceive(ghostforce(y=j-1,z=0));
  }

  for(i=X;i>=1;i--){
    // set force to {+i,0,0}
    rpsf[0].setsend(ghostforce(x=i,y=0,z=0));
    // set force to {-i,0,0}
    rpsf[1].setsend(ghostforce(x=-i,y=0,z=0));
    // prepare receive from node {-1,0,0}
    sprf[1].prepare_receive();
    // prepare receive from node {+1,0,0}
    sprf[0].prepare_receive();
    // send to {+1,0,0}
    rpsf[0].send();
    // send to {-1,0,0}
    rpsf[1].send();
    //receive buff- from {-1,0,0}
    //receive buff+ from {+1,0,0}
    //reduce buff- to force to {+(i-1),0,0}
    //reduce buff+ to force to {-(i-1),0,0}
    sprf[1].getreceive(ghostforce(x=-(i-1),y=0,z=0));
    sprf[0].getreceive(ghostforce(x=i-1,y=0,z=0)); 
  }
  reduce ghostforce(x=0,y=0,z=0) to force to my cell

 // in this stage receive/pass particle was filled dense, and not mention set, must be indicate offset, TODO make getreceive with offset?
 // z last (k=1) received force is for z=0 node (include this), y last (j=1) received force is for y=0 node (include this), x last (i=1) received force is for y=0, this node. It means that ghostforce include dummy area for this node (current ghost is particel/force outside of this node) to receive same operation in getreceive.
 // check number_of_particle in getreceive must be support multi stage, number_of_particle depends stage.
 // buffer_size must multiply 2X+1 for y stage(t=2,3), (2X+1)*(2Y+1) for z stage(t=4,5)
 */


/*******
NearestXYZ : typerange miss transfer ?

***********/
