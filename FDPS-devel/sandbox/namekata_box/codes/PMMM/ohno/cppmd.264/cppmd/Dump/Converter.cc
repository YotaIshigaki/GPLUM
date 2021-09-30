#include <mpi.h>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include "Converter.h"


template<class PA>
Converter<PA>::Converter(const std::string& fbase, int nrsts, MPI_Comm _comm):
  filenamebase(fbase), number_of_restorefiles(nrsts)
{
  comm = _comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_rank);
  if(rank==0){
    std::cout << "Number of rank " << num_rank << std::endl;
    std::cout << "base of filename " << filenamebase << std::endl;
    std::cout << "number of file " << number_of_restorefiles << std::endl;
  }
  number_of_localrestorefiles=0;
  restorefileindexes.clear();
  if(number_of_restorefiles>0){
    /*
    for(int i=0;i<number_of_restorefiles;i++){
      if((i%num_rank)==rank){
	number_of_localrestorefiles++;
	restorefileindexes.push_back(i);
      }
    }
    */
    int m = number_of_restorefiles/num_rank;
    int l = number_of_restorefiles%num_rank;
    int first_index, last_index;
    if(rank<l){
      number_of_localrestorefiles = m+1;
      first_index = (m+1)*rank;
      last_index = (m+1)*(rank+1)-1;
    }else{
      number_of_localrestorefiles = m;
      first_index = (m+1)*l + m*(rank-l);
      last_index = (m+1)*l + m*(rank-l+1)-1;
    }
    for(int i=first_index;i<=last_index;i++){
      restorefileindexes.push_back(i);
    }
    std::cout << "Rank " << rank << " : number of file " << number_of_localrestorefiles << std::endl;
    read_rawimages();
    if(number_of_localrestorefiles>0){
      merge_particle_reindex();
    }else{
      local_number_of_particle=0;
      particlearray.resize(0);
      typerangearray.resize(0);
      local_maxid=0;
    }
    gather_size();
  }else{
    if(rank==0){
      printf("Number of Restorefile(%d) must be atleast 1\n",number_of_restorefiles);
    }
  }
}

template<class PA>
Converter<PA>::~Converter()
{
}

template<class PA>
int Converter<PA>::open_restorefiles()
{
  if(number_of_localrestorefiles<1){
    number_of_localrestorefiles = 0;
    restorefiles.resize(0);
    return 0;
  }
  int opened = 1;
  restorefiles.resize(number_of_localrestorefiles);
  char name[256];
  int i;
  int li=0;
  for(li=0;li<number_of_localrestorefiles;li++){
    i = restorefileindexes[li];
    snprintf(name,256,"%s.%05d",filenamebase.c_str(),i);
    std::string filename(name);
    std::cout << filename << std::endl;
    restorefiles[li].restoremode();
    opened = restorefiles[li].init(filename);
    if(opened!=0){
      std::cout << "Can't open Restorefile " << filename << std::endl;
      break;
    }
  }
  if(opened>0){
    li--;
    for(;li>=0;li--){
      restorefiles[li].close();
    }
  }
  return opened;
}

template<class PA>
int Converter<PA>::read_restorefiles()
{
  particlearrays.resize(number_of_localrestorefiles);
  typerangearrays.resize(number_of_localrestorefiles);
  bondlistarrays.resize(number_of_localrestorefiles);
  bondlistarray_idxs.resize(number_of_localrestorefiles);
  waterlists.resize(number_of_localrestorefiles);
  shakelists.resize(number_of_localrestorefiles);
  timesteps.resize(number_of_localrestorefiles);
  
  for(int i=0;i<number_of_localrestorefiles;i++){
    restorefiles[i].restore_binary_basic(particlearrays[i],
					 typerangearrays[i],
					 bondlistarrays[i],
					 bondlistarray_idxs[i],
					 waterlists[i],
					 shakelists[i],
					 timesteps[i]);
  }
  total_size_of_particlearray = 0;
  local_number_of_particle = 0;
  for(int i=0;i<number_of_localrestorefiles;i++){
    total_size_of_particlearray += particlearrays[i].size();
    local_number_of_particle += restorefiles[i].exact_num_particle();
  }
  printf("Rank %d : size of particlearrays %d\n",rank,total_size_of_particlearray);
  printf("Rank %d : size of typerange(exact number of particle) %d\n",rank,local_number_of_particle);

  return local_number_of_particle;
}

/* it depend Integrator.cc */
/* must be called just after read_restorefiles() */
template<class PA>
void Converter<PA>::read_optionals()
{
  kenergys.resize(number_of_localrestorefiles);
  penergys.resize(number_of_localrestorefiles);
  virials.resize(number_of_localrestorefiles);
  total_penergys.resize(number_of_localrestorefiles);
  eta_poss.resize(number_of_localrestorefiles);
  eta_vels.resize(number_of_localrestorefiles);
  eta_forces.resize(number_of_localrestorefiles);
  logv_poss.resize(number_of_localrestorefiles);
  logv_vels.resize(number_of_localrestorefiles);
  logv_forces.resize(number_of_localrestorefiles);
  boxsizes.resize(number_of_localrestorefiles);
  volumes.resize(number_of_localrestorefiles);
  tljcecs.resize(number_of_localrestorefiles);
  ris.resize(number_of_localrestorefiles);
  pis.resize(number_of_localrestorefiles);
  dcis.resize(number_of_localrestorefiles);
  dris.resize(number_of_localrestorefiles);

  for(int i=0;i<number_of_localrestorefiles;i++){
    double od[15];
    restorefiles[i].restore_binary_optional_double(od, 15);
    kenergys[i] = od[0];
    penergys[i] = od[1];
    virials[i] = od[2];
    total_penergys[i] = od[3];
    eta_poss[i] = od[4];
    eta_vels[i] = od[5];
    eta_forces[i] = od[6];
    logv_poss[i] = od[7];
    logv_vels[i] = od[8];
    logv_forces[i] = od[9];
    boxsizes[i].x = od[10];
    boxsizes[i].y = od[11];
    boxsizes[i].z = od[12];
    volumes[i] = od[13];
    tljcecs[i] = od[14];
    int oi[4];
    restorefiles[i].restore_binary_optional_int(oi, 4);
    ris[i] = oi[0];
    pis[i] = oi[1];
    dcis[i] = oi[2];
    dris[i] = oi[3];
  }
  kenergy = kenergys[0];
  penergy = penergys[0];
  virial = virials[0];
  total_penergy = total_penergys[0];
  eta_pos = eta_poss[0];
  eta_vel = eta_vels[0];
  eta_force = eta_forces[0];
  logv_pos = logv_poss[0];
  logv_vel = logv_vels[0];
  logv_force = logv_forces[0];
  boxsize = boxsizes[0];
  volume = volumes[0];
  tljcec = tljcecs[0];
  ri = ris[0];
  pi = pis[0];
  dci = dcis[0];
  dri = dris[0];
}

template<class PA>
void Converter<PA>::read_optional()
{
  double od[15];
  restorefiles[0].restore_binary_optional_double(od, 15);
  kenergy = od[0];
  penergy = od[1];
  virial = od[2];
  total_penergy = od[3];
  eta_pos = od[4];
  eta_vel = od[5];
  eta_force = od[6];
  logv_pos = od[7];
  logv_vel = od[8];
  logv_force = od[9];
  boxsize.x = od[10];
  boxsize.y = od[11];
  boxsize.z = od[12];
  volume = od[13];
  tljcec = od[14];
  int oi[4];
  restorefiles[0].restore_binary_optional_int(oi, 4);
  ri = oi[0];
  pi = oi[1];
  dci = oi[2];
  dri = oi[3];
}

template<class PA>
void Converter<PA>::read_rawimages()
{
  if(number_of_localrestorefiles>0){
    open_restorefiles();
    read_restorefiles();
    read_optionals();
  }else{
    number_of_localrestorefiles = 0;
    restorefiles.resize(number_of_localrestorefiles);
    particlearrays.resize(number_of_localrestorefiles);
    typerangearrays.resize(number_of_localrestorefiles);
    bondlistarrays.resize(number_of_localrestorefiles);
    bondlistarray_idxs.resize(number_of_localrestorefiles);
    waterlists.resize(number_of_localrestorefiles);
    shakelists.resize(number_of_localrestorefiles);
    timesteps.resize(number_of_localrestorefiles);
  }
}

template<class PA>
int Converter<PA>::read_typerangearrays()
{
  local_number_of_particle = 0;
  typerangearrays.resize(number_of_localrestorefiles);
  particlearrays.resize(number_of_localrestorefiles);
  for(int i=0;i<number_of_localrestorefiles;i++){
    restorefiles[i].restore_binary_typerangearray(typerangearrays[i]);
    AtomID num = (typerangearrays[i][typerangearrays[i].size()].end-typerangearrays[i][0].begin);
    particlearrays[i].resize(num);
    local_number_of_particle += num;
  }
}

template<class PA>
int Converter<PA>::merge_particle()
{
  int n=0;
  local_maxid=0;
  particlearray.resize(local_number_of_particle);
  for(int i=0;i<particlearrays.size();i++){
    for(int tr=0;tr<typerangearrays[i].size();tr++){
      for(int p=typerangearrays[i][tr].begin;p<typerangearrays[i][tr].end;p++){
	AtomID aid = getatomid(particlearrays[i],p);
	if(aid>local_maxid)local_maxid=aid;
	setparticle(particlearray,getparticle(particlearrays[i],p),n);
	n++;
      }
    }
  }
  printf("Rank %d : Merge %d atom, max atomID %d\n",rank,n,local_maxid);
  typerangearray.resize(1);
  typerangearray[0].begin = 0;
  typerangearray[0].end = n;
}

template<class PA>
int Converter<PA>::merge_particle_reindex()
{
  size_t local_number_of_targetrange=0;
  for(int i=0;i<particlearrays.size();i++){
    local_number_of_targetrange += typerangearrays[i].size();
  }
  typerangearray.resize(local_number_of_targetrange);
  particlearray.resize(local_number_of_particle);
  std::vector<AtomID> reindex;
  waterlist.clear();
  shakelist.clear();
  int t=0;
  int n=0;
  local_maxid=0;
  for(int i=0;i<particlearrays.size();i++){
    reindex.resize(particlearrays[i].size());
    for(int tr=0;tr<typerangearrays[i].size();tr++){
      int index_shift = n-typerangearrays[i][tr].begin;
      typerangearray[t] = typerangearrays[i][tr];
      typerangearray[t].shift(index_shift);
      t++;
      for(int p=typerangearrays[i][tr].begin;p<typerangearrays[i][tr].end;p++){
	AtomID aid = getatomid(particlearrays[i],p);
	if(aid>local_maxid)local_maxid=aid;
	setparticle(particlearray,getparticle(particlearrays[i],p),n);
	reindex[p] = n;
	n++;
      }
    }
    for(WaterList::iterator it=waterlists[i].begin();it!=waterlists[i].end();++it){
      int oindex = reindex[it->first];
      int h1index = reindex[it->second.h1];
      int h2index = reindex[it->second.h2];
      waterlist.add_water(oindex,h1index,h2index);
    }
    for(ShakeList::iterator it=shakelists[i].begin();it!=shakelists[i].end();++it){
      int haid = reindex[it->first];
      H1List h1list;
      h1list.nh1 = it->second.nh1;
      int h;
      for(h=0;h<h1list.nh1;h++){
	h1list.h1[h] = reindex[it->second.h1[h]];
	h1list.bondtype[h] = it->second.bondtype[h];
      }
      for(;h<MAX_HBOND;h++){
	h1list.h1[h] = -1;
	h1list.bondtype[h] = -1;
      }
      shakelist.add_shake(haid,h1list);
    }
  }
  
  printf("Rank %d : Merge %d atom, max atomID %d\n",rank,n,local_maxid);
}

template<class PA>
int Converter<PA>::merge_particle_by_atomid(PA& dpa, 
					    const PA& spa, const std::vector<TypeRange>& typerange)
{
  AtomID merge_maxid=dpa.size()-1;

  if(dpa.size()<spa.size())dpa.resize(spa.size());
  for(int tr=0;tr<typerange.size();tr++){
    for(int p=typerange[tr].begin;p<typerange[tr].end;p++){
      AtomID aid = getatomid(spa,p);
      if(aid>=merge_maxid){
	merge_maxid = aid;
	if(merge_maxid>=dpa.size()){
	  dpa.resize(merge_maxid+1);
	}
      }
      setparticle(dpa,getparticle(spa,p),aid);
    }
  }
}

template<class PA>
int Converter<PA>::merge_bondlistarray()
{
  bondlistarray.clear();
  int num_bla = bondlistarrays.size();
  for(int b=0;b<num_bla;b++){
    bondlistarray.insert(bondlistarray.end(),bondlistarrays[b].begin(),bondlistarrays[b].end());
    printf("bondlistarrays[%d].size() = %d\n",b,bondlistarrays[b].size());
  }
}

template<class PA>
int Converter<PA>::merge_waterlists()
{
  waterlist.clear();
  int num_wl = waterlists.size();
  for(int w=0;w<num_wl;w++){
    printf("waterlist %d : %d\n",w,waterlists[w].size());
    waterlist.append(waterlists[w]);
  }
}

template<class PA>
int Converter<PA>::merge_waterlists_atomid()
{
  waterlist.clear();
  int num_wl = waterlists.size();
  for(int w=0;w<num_wl;w++){
    printf("waterlist %d : %d\n",w,waterlists[w].size());
    for(WaterList::iterator it=waterlists[w].begin();it!=waterlists[w].end();++it){
      int oid = getatomid(particlearrays[w],it->first);
      int h1id = getatomid(particlearrays[w],it->second.h1);
      int h2id = getatomid(particlearrays[w],it->second.h2);
      waterlist.add_water(oid,h1id,h2id);
    }
  }
}

template<class PA>
int Converter<PA>::merge_shakelists_atomid()
{
  shakelist.clear();
  int num_sl = shakelists.size();
  for(int s=0;s<num_sl;s++){
    printf("shakelist %d : %d\n",s,shakelists[s].size());
    for(ShakeList::iterator it=shakelists[s].begin();it!=shakelists[s].end();++it){
      int haid = getatomid(particlearrays[s],it->first);
      H1List h1list;
      h1list.nh1 = it->second.nh1;
      int h;
      for(h=0;h<h1list.nh1;h++){
	h1list.h1[h] = getatomid(particlearrays[s],it->second.h1[h]);
	h1list.bondtype[h] = it->second.bondtype[h];
      }
      for(;h<MAX_HBOND;h++){
	h1list.h1[h] = -1;
	h1list.bondtype[h] = -1;
      }
      shakelist.add_shake(haid,h1list);
    }
  }
}


template<class PA>
int Converter<PA>::dump_amber_rst(std::string& filename)
{
  Dump dump_rst(boxsize,total_number_of_particle,num_rank,rank,0,filename,1); 
  dump_rst.GatherDumpAmberCRD(particlearray,typerangearray);
}

template<class PA>
int Converter<PA>::dump_amber_crd(std::string& filename)
{
  Dump dump_crd(boxsize,total_number_of_particle,num_rank,rank,0,filename,0); 
  dump_crd.GatherDumpAmberCRD(particlearray,typerangearray);
}

template<class PA>
int Converter<PA>::gather_size()
{
  int lnum,tnum;
  int lmid,tmid;
  lnum = local_number_of_particle;
  MPI_Allreduce(&lnum,&tnum,1,MPI_INT,MPI_SUM,comm);
  total_number_of_particle = tnum;
  lmid = local_maxid;
  MPI_Allreduce(&lmid,&tmid,1,MPI_INT,MPI_MAX,comm);
  maxid = tmid;
  if(maxid>=total_number_of_particle){
    if(rank==0){
      printf("Max AtomID %d is larger than number of atom %d\n",maxid,total_number_of_particle);
    }
    total_number_of_particle = maxid;
  }
}

template<class PA>
int Converter<PA>::gather_particle(PA& dpa, const PA& spa)
{
  int np;
  int bufsize;
  int singlesize = sizeof(Particle);
  std::vector<Particle> pbuf;
  np = spa.size();
  pbuf.resize(np);

  std::vector<TypeRange> typerange(1);
  typerange[0].begin = 0;
  typerange[0].end = spa.size();
  dpa.resize(total_number_of_particle);
  if(rank==0){
    merge_particle_by_atomid(dpa,spa,typerange);
  }else{
    for(int i=0;i<np;i++){
      pbuf[i] = spa.getparticle(i);
    }
  }
  for(int r=1;r<num_rank;r++){
    int ret;
    if(rank==0){
      MPI_Status stat;
      ret = MPI_Recv((void *)(&np), 1, MPI_INT, r, r, comm, &stat);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Recv from %d\n",ret,r);
	MPI_Abort(comm,ret);
	break;
      }
      if(np>=pbuf.size())pbuf.resize(np);
      ret = MPI_Recv((void *)(&(pbuf[0].position.x)), np*singlesize, MPI_BYTE, r, r, comm, &stat);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Recv from %d\n",ret,r);
	MPI_Abort(comm,ret);
	break;
      }
      for(int i=0;i<np;i++){
	AtomID aid = pbuf[i].atomid;
	if(aid>=dpa.size())dpa.resize(aid+1);
	dpa.setparticle(pbuf[i],aid);
      }
    }else if(rank==r){
      ret = MPI_Send((void *)(&np), 1, MPI_INT, 0, r, comm);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Send to 0\n",ret);
	MPI_Abort(comm,ret);
	break;
      }
      ret = MPI_Send((void *)(&(pbuf[0].position.x)), np*singlesize, MPI_BYTE, 0, r, comm);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Send to 0\n",ret);
	MPI_Abort(comm,ret);
	break;
      }
    }
  }
}

template<class PA>
int Converter<PA>::gather_bondlist(std::vector<CovalentBondInfo::BondList>& dbl,
				   const std::vector<CovalentBondInfo::BondList>& sbl)
{
  int nbl;
  int *bufsize;
  int *buffer = NULL;

  if(rank==0){
    dbl.insert(dbl.end(),sbl.begin(),sbl.end());
  }
  for(int r=1;r<num_rank;r++){
    int ret;
    if(rank==0){
      MPI_Status stat;

      ret = MPI_Recv((void *)(&nbl), 1, MPI_INT, r, r, comm, &stat);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Recv from %d\n",ret,r);
	MPI_Abort(comm,ret);
	break;
      }

      bufsize = new int[nbl];
      ret = MPI_Recv((void *)bufsize, nbl, MPI_INT, r, r, comm, &stat);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Recv from %d\n",ret,r);
	MPI_Abort(comm,ret);
	break;
      }
      int maxsize;
      for(int b=0;b<nbl;b++){
	if(bufsize[b]>maxsize)maxsize=bufsize[b];
      }
      buffer = (int *)realloc(buffer,sizeof(int)*maxsize);

      for(int b=0;b<nbl;b++){
	ret = MPI_Recv((void *)buffer, bufsize[b], MPI_INT, r, r, comm, &stat);
	if(ret!=MPI_SUCCESS){
	  printf("MPI Error %d : MPI_Recv from %d\n",ret,r);
	  MPI_Abort(comm,ret);
	  break;
	}
	CovalentBondInfo::BondList bl;
	bl.unpack_int_array(buffer);
	dbl.push_back(bl);
      }
    }else if(rank==r){

      nbl = sbl.size();
      ret = MPI_Send((void *)(&nbl), 1, MPI_INT, 0, r, comm);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Send to 0\n",ret);
	MPI_Abort(comm,ret);
	break;
      }

      int maxsize = 0;
      bufsize = new int[nbl];
      for(int b=0;b<nbl;b++){
	bufsize[b] = sbl[b].size_of_packed_int();
	if(bufsize[b]>maxsize)maxsize=bufsize[b];
      }
      ret = MPI_Send((void *)(bufsize), nbl, MPI_INT, 0 , r, comm);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Send to 0\n",ret);
	MPI_Abort(comm,ret);
	break;
      }
      buffer = (int *)realloc(buffer,sizeof(int)*maxsize);
      
      for(int b=0;b<nbl;b++){
	sbl[b].pack_int_array(buffer);
	ret = MPI_Send((void *)(buffer), bufsize[b], MPI_INT, 0, r, comm);
	if(ret!=MPI_SUCCESS){
	  printf("MPI Error %d : MPI_Send to 0\n",ret);
	  MPI_Abort(comm,ret);
	  break;
	}
      }
    }
  }
}

template<class PA>
int Converter<PA>::gather_waterlist(WaterList& dwl, const WaterList& swl)
{
  int nw;
  int bufsize;
  int *buffer = NULL;

  if(rank==0){
    dwl.append(swl);
  }

  for(int r=1;r<num_rank;r++){
    int ret;
    if(rank==0){
      MPI_Status stat;

      ret = MPI_Recv((void *)(&nw), 1, MPI_INT, r, r, comm, &stat);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Recv from %d\n",ret,r);
	MPI_Abort(comm,ret);
	break;
      }
      buffer = (int *)realloc(buffer,sizeof(int)*nw*3);

      ret = MPI_Recv((void *)(buffer), nw*3, MPI_INT, r, r, comm, &stat);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Recv from %d\n",ret,r);
	MPI_Abort(comm,ret);
	break;
      }
      for(int w=0;w<nw;w++){
	dwl.add_water(buffer[w*3],buffer[w*3+1],buffer[w*3+2]);
      }

    }else if(rank==r){
      nw = swl.size();
      ret = MPI_Send((void *)(&nw), 1, MPI_INT, 0, r, comm);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Send to 0\n",ret);
	MPI_Abort(comm,ret);
	break;
      }

      buffer = (int *)realloc(buffer,sizeof(int)*nw*3);
      int w=0;
      for(WaterList::const_iterator it=swl.begin();it!=swl.end();++it){
	buffer[w*3] = it->first;
	buffer[w*3+1] = it->second.h1;
	buffer[w*3+2] = it->second.h2;
	w++;
      }
      ret = MPI_Send((void *)(buffer), nw*3, MPI_INT, 0, r, comm);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Send to 0\n",ret);
	MPI_Abort(comm,ret);
	break;
      }
    }
  }
}

template<class PA>
int Converter<PA>::gather_shakelist(ShakeList& dsl, const ShakeList& ssl)
{
  int ns;
  int bufsize;
  const int singlesize = 1+sizeof(H1List);
  int *buffer = NULL;

  if(rank==0){
    dsl.append(ssl);
  }

  for(int r=1;r<num_rank;r++){
    int ret;
    if(rank==0){
      MPI_Status stat;

      ret = MPI_Recv((void *)(&ns), 1, MPI_INT, r, r, comm, &stat);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Recv from %d\n",ret,r);
	MPI_Abort(comm,ret);
	break;
      }
      buffer = (int *)realloc(buffer,sizeof(int)*ns*singlesize);

      ret = MPI_Recv((void *)(buffer), ns*singlesize, MPI_INT, r, r, comm, &stat);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Recv from %d\n",ret,r);
	MPI_Abort(comm,ret);
	break;
      }
      for(int s=0;s<ns;s++){
	H1List *h1list = (H1List *)(&buffer[s*singlesize+1]);
	dsl.add_shake(buffer[s*singlesize],*h1list);
      }

    }else if(rank==r){
      ns = ssl.size();
      ret = MPI_Send((void *)(&ns), 1, MPI_INT, 0, r, comm);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Send to 0\n",ret);
	MPI_Abort(comm,ret);
	break;
      }

      buffer = (int *)realloc(buffer,sizeof(int)*ns*singlesize);
      int s=0;
      for(ShakeList::const_iterator it=ssl.begin();it!=ssl.end();++it){
	buffer[s*singlesize] = it->first;
	H1List *h1list = (H1List *)(&buffer[s*singlesize+1]);
	*h1list = it->second;
	s++;
      }
      ret = MPI_Send((void *)(buffer), ns*singlesize, MPI_INT, 0, r, comm);
      if(ret!=MPI_SUCCESS){
	printf("MPI Error %d : MPI_Send to 0\n",ret);
	MPI_Abort(comm,ret);
	break;
      }
    }
  }
}

template class Converter<CombinedParticleArray>;
