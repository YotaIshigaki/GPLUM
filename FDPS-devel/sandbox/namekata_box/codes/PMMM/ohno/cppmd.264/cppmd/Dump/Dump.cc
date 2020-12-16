#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <mpi.h>
#include "Dump.h"
#include "LJAmber94.h"
#include <algorithm>
#ifdef TTHA_NCHAIN
#include <sstream>
#endif

#ifdef TTHA_VITRO
const int ttha_begin[TTHA_NCHAIN] = {0, 1010, 2020, 3030};
const int ttha_end[TTHA_NCHAIN] = {1010, 2020, 3030, 4040};
const int ttha_natoms = 418707;
#else
#ifdef TTHA_VIVO
/*
const int ttha_begin[TTHA_NCHAIN] = {47984, 97188, 146392, 195596,
                                    244800, 294004, 343208, 392412};
const int ttha_end[TTHA_NCHAIN] =   {48994, 98198, 147402, 196606,
                                    245810, 295014, 344218, 393422};
*/
const int ttha_begin[TTHA_NCHAIN] = { 47984, 216045, 384106, 552167,
                                     720228, 888289, 1056350, 1224411};
const int ttha_end[TTHA_NCHAIN] =   { 48994, 217055, 385116, 553177,
                                     721238, 889299, 1057360, 1225421};
const int ttha_natoms = 1344488;

#endif
#endif


//#define EXCLUDE_TEST

using namespace amberfile;

Dump::Dump(const SpaceVector<double>& bs,
           int np, int nn, int me, int exclude,
#ifdef DEBUG_UNITOUTPUT
           int nc,
#endif  // DEBUG_UNITOUTPUT
           int restart)
  : num_particle(np),
    num_node(nn),
    node_id(me),
#ifdef DEBUG_UNITOUTPUT
    num_copy(nc),
#endif  // DEBUG_UNITOUTPUT
    is_restart(restart),
    fp_(0),
    buffer(),
    exclude_type(exclude)
{
  setBoxSize(bs);
  setExcludeTypes();
  //  init(pa);
}

Dump::Dump(const SpaceVector<double>& bs,
           int np, int nn, int me, int exclude, const std::string& filename,
#ifdef DEBUG_UNITOUTPUT
           int nc,
#endif  // DEBUG_UNITOUTPUT
           int restart)
  : num_particle(np),
    num_node(nn),
    node_id(me),
#ifdef DEBUG_UNITOUTPUT
    num_copy(nc),
#endif  // DEBUG_UNITOUTPUT
    is_restart(restart),
    fp_(0),
    buffer(),
    exclude_type(exclude)
{
  setBoxSize(bs);
  setExcludeTypes();

  init(filename);
}

Dump::~Dump()
{
  CloseTitledFile(fp_);
  if(buffer.begin!=NULL){
    delete [] buffer.begin;
  }
}

void Dump::init(const std::string& filename)
{
#ifdef DEBUG_UNITOUTPUT
  num_particle /= num_copy;
#endif  // DEBUG_UNITOUTPUT
  amber_crd.natom = num_particle;
  amber_crd.coordinates.clear();
  amber_crd.velocities.clear();
  if(node_id==0){
    amber_crd.coordinates.reserve(3 * amber_crd.natom);
    if(is_restart!=0){
      amber_crd.velocities.reserve(3 * amber_crd.natom);
    }
#ifdef TTHA_NCHAIN
    if(is_restart==0){
      for (int chain = 0; chain < TTHA_NCHAIN; chain++) {
        std::string ttha_filename;
        std::ostringstream stream;
        stream << filename << '.' << chain;
        ttha_filename = stream.str();
        fp_ttha[chain] = PrepareTitledFile(amber_crd, ttha_filename);
      }
    }
#else
    fp_ = PrepareTitledFile(amber_crd, filename);
#endif
#if 0
    if (exclude_type > 0) {
      type_list.reserve(amber_crd.natom);
      for(int i=0;i<num_particle;i++){
        amber_crd.coordinates.push_back(pa[i].position.x);
        amber_crd.coordinates.push_back(pa[i].position.y);
        amber_crd.coordinates.push_back(pa[i].position.z);
        type_list.push_back(pa[i].atomtype);
      }
    } else {
      for(int i=0;i<num_particle;i++){
        amber_crd.coordinates.push_back(pa[i].position.x);
        amber_crd.coordinates.push_back(pa[i].position.y);
        amber_crd.coordinates.push_back(pa[i].position.z);
      }
    }
#else
    if (exclude_type > 0) type_list.resize(amber_crd.natom);
#endif
    if(is_restart==0){
      std::cout << "init amber_crd " << num_particle << std::endl;
    }else{
      std::cout << "init amber_rst " << num_particle << std::endl;
    }
  }
  buffer_size = sizeof(buffer.np[0])
    +sizeof(buffer.atomid[0])*num_particle
    +sizeof(buffer.position[0].x)*3*num_particle;
  if (exclude_type > 0) buffer_size += sizeof(buffer.atomtype[0])*num_particle;
  if(is_restart!=0){
    buffer_size += sizeof(buffer.velocity[0].x)*3*num_particle;
  }
  buffer.begin = new char[buffer_size];
  if(buffer.begin==NULL){
    std::cout << "fail new buffer buffer_size " << buffer_size << std::endl;
  }
  void *cp = buffer.begin;
  buffer.np = static_cast<int *>(cp);
  cp = &(buffer.np[1]);
  buffer.atomid = static_cast<int *>(cp);
  cp = &(buffer.atomid[num_particle]);
  buffer.position = static_cast<Position *>(cp);
  if (exclude_type > 0) {
    cp = &(buffer.position[num_particle]);
    buffer.atomtype = static_cast<Atomtype *>(cp);
  }
  if(is_restart!=0){
    if (exclude_type > 0) {
      cp = &(buffer.atomtype[num_particle]);
    }else{
      cp = &(buffer.position[num_particle]);
    }
    buffer.velocity = static_cast<Velocity *>(cp);
  }
}

template<typename PA>
void
Dump::setAmberCRD(const PA& pa,
                  const std::vector<TypeRange>& typerange)
{
  if (exclude_type > 0) {
    for(std::vector<TypeRange>::size_type tr=0;tr<typerange.size();tr++){
      for(int i=typerange[tr].begin;i<typerange[tr].end;i++){
        AtomID aid = getatomid(pa,i);
#ifdef DEBUG_UNITOUTPUT
        if(aid>=num_particle) continue;
#endif  // DEBUG_UNITOUTPUT
        amber_crd.coordinates[3 * aid    ] = getpos(pa,i).x;
        amber_crd.coordinates[3 * aid + 1] = getpos(pa,i).y;
        amber_crd.coordinates[3 * aid + 2] = getpos(pa,i).z;
        type_list[aid] = getatomtype(pa,i);
        if(is_restart!=0){
          amber_crd.velocities[3 * aid    ] = UnitParameter::denormalizeVelocity(getvelocity(pa,i).x/2045.5);
          amber_crd.velocities[3 * aid + 1] = UnitParameter::denormalizeVelocity(getvelocity(pa,i).y/2045.5);
          amber_crd.velocities[3 * aid + 2] = UnitParameter::denormalizeVelocity(getvelocity(pa,i).z/2045.5);
        }
      }
    }
  } else {
    for(std::vector<TypeRange>::size_type tr=0;tr<typerange.size();tr++){
      for(int i=typerange[tr].begin;i<typerange[tr].end;i++){
        AtomID aid = getatomid(pa,i);
#ifdef DEBUG_UNITOUTPUT
        if(aid>=num_particle) continue;
#endif  // DEBUG_UNITOUTPUT
        amber_crd.coordinates[3 * aid    ] = getpos(pa,i).x;
        amber_crd.coordinates[3 * aid + 1] = getpos(pa,i).y;
        amber_crd.coordinates[3 * aid + 2] = getpos(pa,i).z;
        if(is_restart!=0){
          amber_crd.velocities[3 * aid    ] = UnitParameter::denormalizeVelocity(getvelocity(pa,i).x/2045.5);
          amber_crd.velocities[3 * aid + 1] = UnitParameter::denormalizeVelocity(getvelocity(pa,i).y/2045.5);
          amber_crd.velocities[3 * aid + 2] = UnitParameter::denormalizeVelocity(getvelocity(pa,i).z/2045.5);
        }
      }
    }
  }
}

void
Dump::DumpAmberCRD(const ParticleArray& pa,
                   const std::vector<TypeRange>& typerange,
                   const std::string& filename)
{
  setAmberCRD(pa,typerange);
  if (fp_ == 0) fp_ = PrepareTitledFile(amber_crd, filename);
  DumpCoordinatesToFile();
}

template<typename PA>
void
Dump::DumpAmberCRD(const PA& pa, 
                   const std::vector<TypeRange>& typerange)
{
  setAmberCRD(pa,typerange);
  DumpCoordinatesToFile();
}
template
void
Dump::DumpAmberCRD(const CombinedParticleArray& pa, const std::vector<TypeRange>& typerange);

template<typename PA>
void
Dump::GatherDumpAmberCRD(const PA& pa,
                         const std::vector<TypeRange>& typerange)
{
  int np = 0;
  int transfer_size;

#ifdef DEBUG_UNITOUTPUT
  for(std::vector<TypeRange>::size_type tr=0;tr<typerange.size();tr++){
    for(int i=typerange[tr].begin;i<typerange[tr].end;i++){
      if(getatomid(pa,i)<num_particle)
        np++;
    }
  }
#else  // DEBUG_UNITOUTPUT
  for(std::vector<TypeRange>::size_type tr=0;tr<typerange.size();tr++){
    np += (typerange[tr].end-typerange[tr].begin);
  }
#endif  // DEBUG_UNITOUTPUT

  if(node_id==0){
    transfer_size = buffer_size;
    setAmberCRD(pa,typerange);
    //    std::cout << " overwrite " << np << std::endl;
  }else{
    *buffer.np = np;
    void *cp = &(buffer.atomid[np]);
    buffer.position = static_cast<Position *>(cp);
    if (exclude_type > 0) {
      cp = &(buffer.position[np]);
      buffer.atomtype = static_cast<Atomtype *>(cp);
      if(is_restart!=0){
        cp = &(buffer.atomtype[np]);
        buffer.velocity = static_cast<Velocity *>(cp);
      }
      int bp =0;
      for(std::vector<TypeRange>::size_type tr=0;tr<typerange.size();tr++){
        for(int i=typerange[tr].begin;i<typerange[tr].end;i++){
#ifdef DEBUG_UNITOUTPUT
          if(getatomid(pa,i)>=num_particle) continue;
#endif  // DEBUG_UNITOUTPUT
          buffer.atomid[bp] = getatomid(pa,i);
          buffer.position[bp] = getpos(pa,i);
          buffer.atomtype[bp] = getatomtype(pa,i);
          if(is_restart!=0){
            buffer.velocity[bp].x = UnitParameter::denormalizeVelocity(getvelocity(pa,i).x/2045.5);
            buffer.velocity[bp].y = UnitParameter::denormalizeVelocity(getvelocity(pa,i).y/2045.5);
            buffer.velocity[bp].z = UnitParameter::denormalizeVelocity(getvelocity(pa,i).z/2045.5);
          }
          bp++;
        }
      }
      transfer_size = sizeof(buffer.np[0])+sizeof(buffer.atomid[0])*np+
        sizeof(buffer.position[0].x)*3*np+sizeof(buffer.atomtype[0])*np;
      if(is_restart!=0){
        transfer_size += sizeof(buffer.velocity[0].x)*3*np;
      }
    } else {
      if(is_restart!=0){
        cp = &(buffer.position[np]);
        buffer.velocity = static_cast<Velocity *>(cp);
      }
      int bp =0;
      for(std::vector<TypeRange>::size_type tr=0;tr<typerange.size();tr++){
        for(int i=typerange[tr].begin;i<typerange[tr].end;i++){
#ifdef DEBUG_UNITOUTPUT
          if(getatomid(pa,i)>=num_particle) continue;
#endif  // DEBUG_UNITOUTPUT
          buffer.atomid[bp] = getatomid(pa,i);
          buffer.position[bp] = getpos(pa,i);
          if(is_restart!=0){
            buffer.velocity[bp].x = UnitParameter::denormalizeVelocity(getvelocity(pa,i).x/2045.5);
            buffer.velocity[bp].y = UnitParameter::denormalizeVelocity(getvelocity(pa,i).y/2045.5);
            buffer.velocity[bp].z = UnitParameter::denormalizeVelocity(getvelocity(pa,i).z/2045.5);
          }
          bp++;
        }
      }
      transfer_size = sizeof(buffer.np[0])+sizeof(buffer.atomid[0])*np+
        sizeof(buffer.position[0].x)*3*np;
      if(is_restart!=0){
        transfer_size += sizeof(buffer.velocity[0].x)*3*np;
      }
    }
    if(transfer_size>buffer_size){
      std::cout << "np " <<  *(buffer.np) << "  transfer_size " << transfer_size << " > " << buffer_size << " buffer_size " << std::endl;
    }
  }
  MPI_Request req;
  for(int n=1;n<num_node;n++){
    if(node_id==0){
#if 0
      MPI_Irecv(buffer.begin,transfer_size,MPI_BYTE,n,n,MPI_COMM_WORLD,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
#else
      MPI_Recv(buffer.begin,transfer_size,MPI_BYTE,n,n,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
#endif
      np = *(buffer.np);
      int recvsize = sizeof(buffer.np[0])+sizeof(buffer.atomid[0])*np+
          sizeof(buffer.position[0].x)*3*np;
      if (exclude_type > 0) recvsize += sizeof(buffer.atomtype[0])*np;
      if (is_restart!=0) recvsize += sizeof(buffer.velocity[0].x)*3*np;
      if (recvsize>buffer_size) {
        printf("receive number of atom %d is large\n",np);
      }
      void *cp = &(buffer.atomid[np]);
      buffer.position = static_cast<Position *>(cp);
      if (exclude_type > 0) {
        cp = &(buffer.position[np]);
        buffer.atomtype = static_cast<Atomtype *>(cp);
        if(is_restart!=0){
          cp = &(buffer.atomtype[np]);
          buffer.velocity = static_cast<Velocity *>(cp);
        }        
        for(int i=0;i<np;i++){
          //std::cout<<" overwrite "<<buffer.atomid[i]<<" "<<i<<std::endl;
          amber_crd.coordinates[3*buffer.atomid[i]  ] = buffer.position[i].x;
          amber_crd.coordinates[3*buffer.atomid[i]+1] = buffer.position[i].y;
          amber_crd.coordinates[3*buffer.atomid[i]+2] = buffer.position[i].z;
          type_list[buffer.atomid[i]] = buffer.atomtype[i];
          if(is_restart!=0){
            amber_crd.velocities[3*buffer.atomid[i]  ] = buffer.velocity[i].x;
            amber_crd.velocities[3*buffer.atomid[i]+1] = buffer.velocity[i].y;
            amber_crd.velocities[3*buffer.atomid[i]+2] = buffer.velocity[i].z;
          }
        }
        //std::cout<<" overwrite "<<np<<" from "<<n<<std::endl;
      } else {
        if(is_restart!=0){
          cp = &(buffer.position[np]);
          buffer.velocity = static_cast<Velocity *>(cp);
        }
        for(int i=0;i<np;i++){
          //std::cout<<" overwrite "<<buffer.atomid[i]<<" "<<i<<std::endl;
          amber_crd.coordinates[3*buffer.atomid[i]  ] = buffer.position[i].x;
          amber_crd.coordinates[3*buffer.atomid[i]+1] = buffer.position[i].y;
          amber_crd.coordinates[3*buffer.atomid[i]+2] = buffer.position[i].z;
          if(is_restart!=0){
            amber_crd.velocities[3*buffer.atomid[i]  ] = buffer.velocity[i].x;
            amber_crd.velocities[3*buffer.atomid[i]+1] = buffer.velocity[i].y;
            amber_crd.velocities[3*buffer.atomid[i]+2] = buffer.velocity[i].z;
          }
        }
        //std::cout<<" overwrite "<<np<<" from "<<n<<std::endl;
      }
    }else if(node_id==n){
#if 0
      MPI_Isend(buffer.begin,transfer_size,MPI_BYTE,0,n,MPI_COMM_WORLD,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
#else
      MPI_Send(buffer.begin,transfer_size,MPI_BYTE,0,n,MPI_COMM_WORLD);
#endif
    }
  }
  
  if(node_id==0){
    DumpCoordinatesToFile();
  }
}
template
void
Dump::GatherDumpAmberCRD(const ParticleArray& pa, const std::vector<TypeRange>& typerange);
template
void
Dump::GatherDumpAmberCRD(const CombinedParticleArray& pa, const std::vector<TypeRange>& typerange);

void Dump::setBoxSize(const SpaceVector<double>& bs)
{
  boxsize = bs;
  amber_crd.box.resize(6);
  amber_crd.box[0] = boxsize.x;
  amber_crd.box[1] = boxsize.y;
  amber_crd.box[2] = boxsize.z;
  amber_crd.box[3] = 90.0;      // !!!
  amber_crd.box[4] = 90.0;      // !!!
  amber_crd.box[5] = 90.0;      // !!!
}

void Dump::DumpCoordinatesToFile() {
  if (exclude_type > 0) {
#ifdef EXCLUDE_TEST
    // for debug output
    for (int s = 0; s < amber_crd.natom; ++s) {
      if (std::find(exclude_types.begin(), exclude_types.end(), type_list[s]) !=
          exclude_types.end()) {
        // if the atom[s] is to be excluded
        amber_crd.coordinates[3*s  ] = 0.0;
        amber_crd.coordinates[3*s+1] = 0.0;
        amber_crd.coordinates[3*s+2] = 0.0;
      }
    }
# if 0
    std::cout << "mdcrd : natom " << amber_crd.natom << std::endl;
# endif
    AppendCoordinatesToFile(amber_crd, fp_);
#else  // ! EXCLUDE_TEST
    if(is_restart==0){
#ifdef TTHA_NCHAIN
      int saved_natom = amber_crd.natom;
      for (int chain = 0; chain < TTHA_NCHAIN; chain++) {
        int d = 0;
        for (int s = ttha_begin[chain]; s < ttha_end[chain]; ++s) {
          amber_crd.coordinates[3*d  ] = amber_crd.coordinates[3*s  ];
          amber_crd.coordinates[3*d+1] = amber_crd.coordinates[3*s+1];
          amber_crd.coordinates[3*d+2] = amber_crd.coordinates[3*s+2];
          ++d;
        }
        amber_crd.natom = d;
        AppendCoordinatesToFile(amber_crd, fp_ttha[chain]);
      }
      amber_crd.natom = saved_natom;
# if 1
      std::cout << "mdcrd : natom " << amber_crd.natom << std::endl;
# endif
#else  // ! TTHA_NCHAIN
      int d = 0;
      for (int s = 0; s < amber_crd.natom; ++s) {
        if (std::find(exclude_types.begin(), exclude_types.end(), type_list[s]) ==
            exclude_types.end()) {
          // if the atom[s] is not to be excluded
          amber_crd.coordinates[3*d  ] = amber_crd.coordinates[3*s  ];
          amber_crd.coordinates[3*d+1] = amber_crd.coordinates[3*s+1];
          amber_crd.coordinates[3*d+2] = amber_crd.coordinates[3*s+2];
          ++d;
        }
      }
# if 1
      std::cout << "mdcrd : natom " << amber_crd.natom << " -> " << d <<
          ", (" << amber_crd.natom - d << " excluded atoms)" << std::endl;
# endif
      int saved_natom = amber_crd.natom;
      amber_crd.natom = d;
      AppendCoordinatesToFile(amber_crd, fp_);
      amber_crd.natom = saved_natom;
#endif // ! TTHA_NCHAIN
    }else{  // restart, with velocity
      WriteToFile(amber_crd, fp_, 1);
    }
#endif // ! EXCLUDE_TEST
  } else {
    if(is_restart==0){
      AppendCoordinatesToFile(amber_crd, fp_);
    }else{
      WriteToFile(amber_crd, fp_, 1);
    }
  }
}

void Dump::setExcludeTypes() {
  if (exclude_type > 0) {
    Atomtype at = 0;
    do{
      if("HW"==ljamber94parameters[at].name||
         "OW"==ljamber94parameters[at].name)
        exclude_types.push_back(at);
      at++;
    }while(ljamber94parameters[at].name!=LJAmberParameterEP.name);
  }
}

inline void bswapn(unsigned char *d, const unsigned char *s, int n)
{
  for(int c=0;c<n;c++){
    d[c] = s[n-1-c];
  }
}

inline int intswap(int i)
{
#ifdef bswap_32
  return bswap_32(i);
#else
  int is;
  bswapn((unsigned char *)(&is),(unsigned char *)(&i),4);
  return is;
#endif
}

inline long longswap(long l)
{
#if (SIZEOF_LONG==4)
  return intswap(int l);
#else
# ifdef bswap_64
  return bswap_64(l);
# else
  long ls;
  bswapn((unsigned char *)(&ls),(unsigned char *)(&l),8);
  return ls;
# endif
#endif
}

inline double doubleswap(double d)
{
#ifdef bswap_64
  return bswap_64(d);
#else
  double ds;
  bswapn((unsigned char *)(&ds),(unsigned char *)(&d),8);
  return ds;
#endif
}

template<class PA>
BinaryDump<PA>::BinaryDump(const int _mode)
  : mode(_mode), native_byteorder(true), count(0)
{
  fp = NULL;
}

template<class PA>
BinaryDump<PA>::BinaryDump(const std::string& filename, const int _mode)
  : mode(_mode), native_byteorder(true), count(0)
{
  init(filename);
}

template<class PA>
BinaryDump<PA>::~BinaryDump()
{
  if(fp!=NULL)fclose(fp);
}

template<class PA>
void BinaryDump<PA>::dumpmode()
{
  mode = BDUMP;
}

template<class PA>
void BinaryDump<PA>::restoremode()
{
  mode = BRESTORE;
}

template<class PA>
int BinaryDump<PA>::init(const std::string& filename)
{
  if(mode==BRESTORE){
    fp = fopen(filename.c_str(), "r");
  }else{
    fp = fopen(filename.c_str(), "w");
  }
  if(fp==NULL) return 1;
  return 0;
}

template<class PA>
int BinaryDump<PA>::close()
{
  if(fp==NULL) return 0;
  int ret = fclose(fp);
  fp = NULL;
  return ret;
}

template<class PA>
inline 
void BinaryDump<PA>::writeint(const int i)
{
  fwrite((const void *)(&i), sizeof(i), 1, fp);
#ifdef DUMP_COUNT
  count += sizeof(i);
#endif
}
template<class PA>
inline 
void BinaryDump<PA>::writelong(const long l)
{
  fwrite((const void *)(&l), sizeof(l), 1, fp);
#ifdef DUMP_COUNT
  count += sizeof(l);
#endif
}
template<class PA>
inline 
void BinaryDump<PA>::writedouble(const double d)
{
  fwrite((const void *)(&d), sizeof(d), 1, fp);
#ifdef DUMP_COUNT
  count += sizeof(d);
#endif
}

template<class PA>
inline
int BinaryDump<PA>::readint()
{
  int ret;
  fread((void *)(&ret), sizeof(ret), 1, fp);
#ifdef DUMP_COUNT
  count += sizeof(ret);
#endif
  if(native_byteorder){
    return ret;
  }else{
    int sret = intswap(ret);
    return sret;
  }
}
template<class PA>
inline
long BinaryDump<PA>::readlong()
{
  long ret;
  fread((void *)(&ret), sizeof(ret), 1, fp);
#ifdef DUMP_COUNT
  count += sizeof(ret);
#endif
  if(native_byteorder){
    return ret;
  }else{
    long sret = longswap(ret);
    return sret;
  }
}
template<class PA>
inline
double BinaryDump<PA>::readdouble()
{
  double ret;
  fread((void *)(&ret), sizeof(ret), 1, fp);
#ifdef DUMP_COUNT
  count += sizeof(ret);
#endif
  if(native_byteorder){
    return ret;
  }else{
    double sret = doubleswap(ret);
    return sret;
  }
}

#define BOMARK 0x01234567

template<class PA>
void BinaryDump<PA>::dump_binary_basic(const PA& particlearray,
				   const std::vector<TypeRange>& typerangearray,
				   const std::vector<CovalentBondInfo::BondList>& bondlistarray,
				   const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
				   const WaterList& waterlist,
				   const ShakeList& shakelist,
				   const long t)
{
  if(fp==NULL) return;
  writeint(BOMARK);
  int num_tr = typerangearray.size();
  writeint(num_tr);
  num_particle = 0;
  for(int tr=0;tr<num_tr;tr++){
    writeint(typerangearray[tr].begin);
    writeint(typerangearray[tr].end);
    writeint(typerangearray[tr].lj.begin);
    writeint(typerangearray[tr].lj.end);
    writeint(typerangearray[tr].ljcoulomb.begin);
    writeint(typerangearray[tr].ljcoulomb.end);
    writeint(typerangearray[tr].coulomb.begin);
    writeint(typerangearray[tr].coulomb.end);
    num_particle += (typerangearray[tr].end-typerangearray[tr].begin);
  }
  for(int tr=0;tr<num_tr;tr++){
    for(int i=typerangearray[tr].begin;i<typerangearray[tr].end;i++){
      writedouble(getpos(particlearray,i).x);
      writedouble(getpos(particlearray,i).y);
      writedouble(getpos(particlearray,i).z);
      writedouble(getcharge(particlearray,i));
      writedouble(getvelocity(particlearray,i).x);
      writedouble(getvelocity(particlearray,i).y);
      writedouble(getvelocity(particlearray,i).z);
#ifndef DUMP_WITHOUT_FORCE
      writedouble(getforce(particlearray,i).x);
      writedouble(getforce(particlearray,i).y);
      writedouble(getforce(particlearray,i).z);
#endif
      writedouble(getmass(particlearray,i));
      writedouble(getinvmass(particlearray,i));
      writeint(getatomtype(particlearray,i));
      writeint(getatomid(particlearray,i));
    }
  }
  int num_bondlist = bondlistarray.size();
  writeint(num_bondlist);
  for(int b=0;b<num_bondlist;b++){
    int num_bond =bondlistarray[b].BondArray.size();
    writeint(num_bond);
    for(int i=0;i<num_bond;i++){
      writeint(bondlistarray[b].BondArray[i].id_of_atom[0]);
      writeint(bondlistarray[b].BondArray[i].id_of_atom[1]);
      writeint(bondlistarray[b].BondArray[i].typeofbond);
      writeint(bondlistarray[b].BondArray[i].shake);
    }
    int num_angle = bondlistarray[b].AngleArray.size();
    writeint(num_angle);
    for(int i=0;i<num_angle;i++){
      writeint(bondlistarray[b].AngleArray[i].id_of_atom[0]);
      writeint(bondlistarray[b].AngleArray[i].id_of_atom[1]);
      writeint(bondlistarray[b].AngleArray[i].id_of_atom[2]);
      writeint(bondlistarray[b].AngleArray[i].typeofangle);
    }
    int num_torsion = bondlistarray[b].TorsionArray.size();
    writeint(num_torsion);
    for(int i=0;i<num_torsion;i++){
      writeint(bondlistarray[b].TorsionArray[i].id_of_atom[0]);
      writeint(bondlistarray[b].TorsionArray[i].id_of_atom[1]);
      writeint(bondlistarray[b].TorsionArray[i].id_of_atom[2]);
      writeint(bondlistarray[b].TorsionArray[i].id_of_atom[3]);
      writeint(bondlistarray[b].TorsionArray[i].typeoftorsion);
      int c14 = (bondlistarray[b].TorsionArray[i].calc14interaction ? 1 : 0);
      writeint(c14);
    }
    int num_improper = bondlistarray[b].ImproperArray.size();
    writeint(num_improper);
    for(int i=0;i<num_improper;i++){
      writeint(bondlistarray[b].ImproperArray[i].id_of_atom[0]);
      writeint(bondlistarray[b].ImproperArray[i].id_of_atom[1]);
      writeint(bondlistarray[b].ImproperArray[i].id_of_atom[2]);
      writeint(bondlistarray[b].ImproperArray[i].id_of_atom[3]);
      writeint(bondlistarray[b].ImproperArray[i].typeofimproper);
    }
  }
  int num_waterlist = waterlist.size();
  writeint(num_waterlist);
  for(WaterList::const_iterator w=waterlist.begin();w!=waterlist.end();++w){
    writeint(w->first);
    writeint(w->second.h1);
    writeint(w->second.h2);
  }
  int num_shakelist = shakelist.size();
  writeint(num_shakelist);
  for(ShakeList::const_iterator s=shakelist.begin();s!=shakelist.end();++s){
    writeint(s->first);
    int nh1=s->second.nh1;
    writeint(nh1);
    for(int h=0;h<nh1;h++){
      writeint(s->second.h1[h]);
      writeint(s->second.bondtype[h]);
    }
  }
  writelong(t);
}

template<class PA>
void BinaryDump<PA>::restore_binary_typerangearray(std::vector<TypeRange>& typerangearray)
{
  if(fp==NULL) return;
  int bo = readint();
  native_byteorder = (bo==BOMARK);
  {
    if(!native_byteorder)printf("byteorder\n");
  }
  int num_tr = readint();
  typerangearray.resize(num_tr);
  for(int t=0;t<num_tr;t++){
    typerangearray[t].begin = readint();
    typerangearray[t].end = readint();
    typerangearray[t].lj.begin = readint();
    typerangearray[t].lj.end = readint();
    typerangearray[t].ljcoulomb.begin = readint();
    typerangearray[t].ljcoulomb.end = readint();
    typerangearray[t].coulomb.begin = readint();
    typerangearray[t].coulomb.end = readint();
  }
  rewind(fp);
}

template<class PA>
void BinaryDump<PA>::restore_binary_basic(PA& particlearray,
				      std::vector<TypeRange>& typerangearray,
				      std::vector<CovalentBondInfo::BondList>& bondlistarray,
				      std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
				      WaterList& waterlist,
				      ShakeList& shakelist,
				      long& t)
{
  if(fp==NULL) return;
  int bo = readint();
  native_byteorder = (bo==BOMARK);
  {
    if(!native_byteorder)printf("byteorder\n");
  }
  int num_tr = readint();
  typerangearray.resize(num_tr);
  int max_id=0;
  num_particle = 0;
  for(int tr=0;tr<num_tr;tr++){
    typerangearray[tr].begin = readint();
    typerangearray[tr].end = readint();
    typerangearray[tr].lj.begin = readint();
    typerangearray[tr].lj.end = readint();
    typerangearray[tr].ljcoulomb.begin = readint();
    typerangearray[tr].ljcoulomb.end = readint();
    typerangearray[tr].coulomb.begin = readint();
    typerangearray[tr].coulomb.end = readint();
    if(typerangearray[tr].end>max_id)max_id=typerangearray[tr].end;
    num_particle += (typerangearray[tr].end-typerangearray[tr].begin);
  }
  if(particlearray.size()<max_id){
    printf("Max typerange %d is lager than reserved size of paritclearray %d : resize particleaaray\n",max_id,particlearray.size());
    particlearray.resize(max_id);
  }
  for(int tr=0;tr<num_tr;tr++){
    for(int i=typerangearray[tr].begin;i<typerangearray[tr].end;i++){
      getpos(particlearray,i).x = readdouble();
      getpos(particlearray,i).y = readdouble();
      getpos(particlearray,i).z = readdouble();
      getcharge(particlearray,i) = readdouble();
      getvelocity(particlearray,i).x = readdouble();
      getvelocity(particlearray,i).y = readdouble();
      getvelocity(particlearray,i).z = readdouble();
#ifndef DUMP_WITHOUT_FORCE
      getforce(particlearray,i).x = readdouble();
      getforce(particlearray,i).y = readdouble();
      getforce(particlearray,i).z = readdouble();
#endif
      getmass(particlearray,i) = readdouble();
      getinvmass(particlearray,i) = readdouble();
      getatomtype(particlearray,i) = readint();
      getatomid(particlearray,i) = readint();
    }
  }
  int num_bondlist = readint();
  bondlistarray.resize(num_bondlist);
  for(int b=0;b<num_bondlist;b++){
    int num_bond = readint();
    bondlistarray[b].BondArray.resize(num_bond);
    for(int i=0;i<num_bond;i++){
      bondlistarray[b].BondArray[i].id_of_atom[0] = readint();
      bondlistarray[b].BondArray[i].id_of_atom[1] = readint();
      bondlistarray[b].BondArray[i].typeofbond = readint();
      bondlistarray[b].BondArray[i].shake = readint();
    }
    int num_angle = readint();
    bondlistarray[b].AngleArray.resize(num_angle);
    for(int i=0;i<num_angle;i++){
      bondlistarray[b].AngleArray[i].id_of_atom[0] = readint();
      bondlistarray[b].AngleArray[i].id_of_atom[1] = readint();
      bondlistarray[b].AngleArray[i].id_of_atom[2] = readint();
      bondlistarray[b].AngleArray[i].typeofangle = readint();
    }
    int num_torsion = readint();
    bondlistarray[b].TorsionArray.resize(num_torsion);
    for(int i=0;i<num_torsion;i++){
      bondlistarray[b].TorsionArray[i].id_of_atom[0] = readint();
      bondlistarray[b].TorsionArray[i].id_of_atom[1] = readint();
      bondlistarray[b].TorsionArray[i].id_of_atom[2] = readint();
      bondlistarray[b].TorsionArray[i].id_of_atom[3] = readint();
      bondlistarray[b].TorsionArray[i].typeoftorsion = readint();
      int c14 = readint();
      bondlistarray[b].TorsionArray[i].calc14interaction = (c14==1);
    }
    int num_improper = readint();
    bondlistarray[b].ImproperArray.resize(num_improper);
    for(int i=0;i<num_improper;i++){
      bondlistarray[b].ImproperArray[i].id_of_atom[0] = readint();
      bondlistarray[b].ImproperArray[i].id_of_atom[1] = readint();
      bondlistarray[b].ImproperArray[i].id_of_atom[2] = readint();
      bondlistarray[b].ImproperArray[i].id_of_atom[3] = readint();
      bondlistarray[b].ImproperArray[i].typeofimproper = readint();
    }
  }
  int num_waterlist = readint();
  waterlist.clear();
  waterlist.reverse_list.clear();
  for(int w=0;w<num_waterlist;w++){
    int o = readint();
    int h1 = readint();
    int h2 = readint();
    waterlist.add_water(o,h1,h2);
  }
  int num_shakelist = readint();
  shakelist.clear();
  shakelist.reverse_list.clear();
  for(int s=0;s<num_shakelist;s++){
    int ha = readint();
    int nh1 = readint();
    for(int h=0;h<nh1;h++){
      int h1 = readint();
      int bondtype = readint();
      shakelist.add_shake(ha,h1,bondtype);
    }
  }
  t = readlong();
}


template<class PA>
void BinaryDump<PA>::merge_binary_basic_by_atomid(PA& particlearray,
					     CovalentBondInfo::BondList& bondlist,
					     WaterList& waterlist,
					     ShakeList& shakelist,
					     long& t)
{
  if(fp==NULL) return;

  int bo = readint();
  native_byteorder = (bo==BOMARK);
  {
    if(!native_byteorder)printf("byteorder\n");
  }

  int num_tr = readint();
  std::vector<TypeRange> typerangearray(num_tr);

  int num_merge_atom = 0;

  num_particle = particlearray.size();
  
  printf("Merge %d particle sets\n",num_tr);

  for(int tr=0;tr<num_tr;tr++){
    typerangearray[tr].begin = readint();
    typerangearray[tr].end = readint();
    typerangearray[tr].lj.begin = readint();
    typerangearray[tr].lj.end = readint();
    typerangearray[tr].ljcoulomb.begin = readint();
    typerangearray[tr].ljcoulomb.end = readint();
    typerangearray[tr].coulomb.begin = readint();
    typerangearray[tr].coulomb.end = readint();
    num_merge_atom += (typerangearray[tr].end-typerangearray[tr].begin);
  }
  for(int i=0;i<num_merge_atom;i++){
    double x = readdouble();
    double y = readdouble();
    double z = readdouble();
    double c = readdouble();
    double vx = readdouble();
    double vy = readdouble();
    double vz = readdouble();
#ifndef DUMP_WITHOUT_FORCE
    double fx = readdouble();
    double fy = readdouble();
    double fz = readdouble();
#endif
    double m = readdouble();
    double im = readdouble();
    Atomtype at = readint();
    AtomID ai = readint();
    if(ai>=num_particle){
      printf("AtomID %ld is larger than reserved number of particle %ld\n",ai,num_particle);
      printf("Resize particle array\n");
      num_particle = ai+1;
      particlearray.resize(num_particle);
    }
    getpos(particlearray,ai).x = x;
    getpos(particlearray,ai).y = y;
    getpos(particlearray,ai).z = z;
    getcharge(particlearray,ai) = c;
    getvelocity(particlearray,ai).x = vx;
    getvelocity(particlearray,ai).y = vy;
    getvelocity(particlearray,ai).z = vz;
#ifndef DUMP_WITHOUT_FORCE
    getforce(particlearray,ai).x = fx;
    getforce(particlearray,ai).y = fy;
    getforce(particlearray,ai).z = fz;
#endif
    getmass(particlearray,ai) = m;
    getinvmass(particlearray,ai) = im;
    getatomtype(particlearray,ai) = at;
    getatomid(particlearray,ai) = ai;
  }

  int num_bondlist = readint();
  printf("Merge %d bondlists\n",num_bondlist);
  for(int b=0;b<num_bondlist;b++){
    int num_original_bond = bondlist.BondArray.size();
    int num_bond = num_original_bond + readint();
    bondlist.BondArray.resize(num_bond);
    for(int i=num_original_bond;i<num_bond;i++){
      bondlist.BondArray[i].id_of_atom[0] = readint();
      bondlist.BondArray[i].id_of_atom[1] = readint();
      bondlist.BondArray[i].typeofbond = readint();
      bondlist.BondArray[i].shake = readint();
    }
    int num_original_angle = bondlist.AngleArray.size();
    int num_angle = num_original_angle + readint();
    bondlist.AngleArray.resize(num_angle);
    for(int i=num_original_angle;i<num_angle;i++){
      bondlist.AngleArray[i].id_of_atom[0] = readint();
      bondlist.AngleArray[i].id_of_atom[1] = readint();
      bondlist.AngleArray[i].id_of_atom[2] = readint();
      bondlist.AngleArray[i].typeofangle = readint();
    }
    int num_original_torsion = bondlist.TorsionArray.size();
    int num_torsion = num_original_torsion + readint();
    bondlist.TorsionArray.resize(num_torsion);
    for(int i=num_original_torsion;i<num_torsion;i++){
      bondlist.TorsionArray[i].id_of_atom[0] = readint();
      bondlist.TorsionArray[i].id_of_atom[1] = readint();
      bondlist.TorsionArray[i].id_of_atom[2] = readint();
      bondlist.TorsionArray[i].id_of_atom[3] = readint();
      bondlist.TorsionArray[i].typeoftorsion = readint();
      int c14 = readint();
      bondlist.TorsionArray[i].calc14interaction = (c14==1);
    }
    int num_original_improper = bondlist.ImproperArray.size();
    int num_improper = num_original_improper + readint();
    bondlist.ImproperArray.resize(num_improper);
    for(int i=num_original_improper;i<num_improper;i++){
      bondlist.ImproperArray[i].id_of_atom[0] = readint();
      bondlist.ImproperArray[i].id_of_atom[1] = readint();
      bondlist.ImproperArray[i].id_of_atom[2] = readint();
      bondlist.ImproperArray[i].id_of_atom[3] = readint();
      bondlist.ImproperArray[i].typeofimproper = readint();
    }
  }

  int num_waterlist = readint();
  printf("Merge %d waterlist\n",num_waterlist);
  for(int w=0;w<num_waterlist;w++){
    int o = readint();
    int h1 = readint();
    int h2 = readint();
    waterlist.add_water(o,h1,h2);
  }

  int num_shakelist = readint();
  printf("Merge %d shakelist\n",num_shakelist);
  for(int s=0;s<num_shakelist;s++){
    int ha = readint();
    int nh1 = readint();
    for(int h=0;h<nh1;h++){
      int h1 = readint();
      int bondtype = readint();
      shakelist.add_shake(ha,h1,bondtype);
    }
  }

  long tt = readlong();
  if(tt!=t){
    printf("Time step differ : original %ld, new %ld",t,tt);
    if(tt>t)t=tt;
    printf(" : use larger one %ld\n",t);
  }
}


template<class PA>
void BinaryDump<PA>::dump_binary_optional_double(const double *od, const int num)
{
  if(fp==NULL) return;
  for(int n=0;n<num;n++){
    writedouble(od[n]);
  }
}

template<class PA>
void BinaryDump<PA>::restore_binary_optional_double(double *od, const int num)
{
  if(fp==NULL) return;
  for(int n=0;n<num;n++){
    od[n] = readdouble();
  }
}

template<class PA>
void BinaryDump<PA>::dump_binary_optional_int(const int *oi, const int num)
{
  if(fp==NULL) return;
  for(int n=0;n<num;n++){
    writeint(oi[n]);
  }
#ifdef DUMP_COUNT
  printf("write count %ld\n",count);
#endif
}

template<class PA>
void BinaryDump<PA>::restore_binary_optional_int(int *oi, const int num)
{
  if(fp==NULL) return;
  for(int n=0;n<num;n++){
    oi[n] = readint();
  }
#ifdef DUMP_COUNT
  printf("read count %ld\n",count);
#endif
}

template<class PA>
int BinaryDump<PA>::exact_num_particle()
{
  return num_particle;
}

#ifdef OLDPARTICLE
template class BinaryDump<ParticleArray>;
#else
template class BinaryDump<CombinedParticleArray>;
#endif


