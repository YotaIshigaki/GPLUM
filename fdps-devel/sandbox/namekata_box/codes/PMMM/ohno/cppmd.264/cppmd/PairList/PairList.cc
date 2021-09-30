#include <cstdio>
#include <algorithm>
#include "PairList.h"
#include "ShortRangeInteractionSet.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// hold j index only
bool addljpairlist(int *pid,
                    int *lj,
                    int *npair,
                    const ParticlePosCharge& particlei,
                    const Atomtype& atomi,
                    const PosChargeArray& particlej,
                    const AtomtypeArray& atomja,
                    const TypeRange& typerange,
                    const double cutoff2,
                    const double margin2,
                    const bool self,
                    const int i,
                    const int offset)
{
  bool success=true;

  int np=(*npair) + offset;

  double cut2 = cutoff2 + margin2;
  int j;
  int jmin, jmax;
  jmin = typerange.lj.begin;
  jmax = typerange.ljcoulomb.end;
  for(j=jmin;j<jmax;j++){
    if((self==false)||(i!=j)){
      Position d = particlei.position - particlej[j].position;
      double r2 = d.norm2();
      if(r2<cut2){
        pid[np] = j;
        int atomj = atomja[j];
        lj[np] = atomj;
        np++;
        if(np>=MAX_PAIR){
          printf("number of pair reach to MAX_PAIR %d\n",np);
          break;
        }
      }
    }
  }
  (*npair) = np-offset;
  if(j<typerange.ljcoulomb.end){
    printf("pair list full\n");
    success = false;
  }
  
  return success;
}
// hold j index only
bool addljcpairlist(int *pid,
                    int *lj,
                    int *npair,
                    const ParticlePosCharge& particlei,
                    const Atomtype& atomi,
                    const PosChargeArray& particlej,
                    const AtomtypeArray& atomja,
                    const TypeRange& typerange,
                    const double cutoff2,
                    const double margin2,
                    const bool self,
                    const int i,
                    const int offset)
{
  bool success=true;

  int np=(*npair) + offset;

  double cut2 = cutoff2 + margin2;
  int j;
  int jmin, jmax;
  jmin = typerange.begin;
  jmax = typerange.end;
  for(j=jmin;j<jmax;j++){
    if((self==false)||(i!=j)){
      Position d = particlei.position - particlej[j].position;
      double r2 = d.norm2();
      if(r2<cut2){
        pid[np] = j;
        int atomj = atomja[j];
        lj[np] = atomj;
        np++;
        if(np>=MAX_PAIR){
          printf("number of pair reach to MAX_PAIR %d\n",np);
          break;
        }
      }
    }
  }
  (*npair) = np-offset;
  if(j<typerange.end){
    printf("pair list full\n");
    success = false;
  }
  
  return success;
}

bool addljcpairlist(int *pid,
                    int *lj,
                    int *npair,
                    const ParticlePosCharge& particlei,
                    const Atomtype& atomi,
		    const WHPair& waterh,
                    const PosChargeArray& particlej,
                    const AtomtypeArray& atomja,
                    const TypeRange& typerange,
                    const double cutoff2,
                    const double margin2,
                    const bool self,
                    const int i,
                    const int offset)
{
  bool success=true;

  int np=(*npair) + offset;

  double cut2 = cutoff2 + margin2;
  int j;
  int jmin, jmax;
  jmin = typerange.begin;
  jmax = typerange.end;
  for(j=jmin;j<jmax;j++){
    if((self==false)||(i!=j)){
      if(self){
	if((j==waterh.h1)||(j==waterh.h2))continue;
      }
      Position d = particlei.position - particlej[j].position;
      double r2 = d.norm2();
      if(r2<cut2){
        pid[np] = j;
        int atomj = atomja[j];
        lj[np] = atomj;
        np++;
        if(np>=MAX_PAIR){
          printf("number of pair reach to MAX_PAIR %d\n",np);
          break;
        }
      }
    }
  }
  (*npair) = np-offset;
  if(j<typerange.end){
    printf("pair list full\n");
    success = false;
  }
  
  return success;
}

bool addljcpairlist(int *pid,
                    int *lj,
                    int *npair,
                    const ParticlePosCharge& particlei,
                    const Atomtype& atomi,
		    const std::vector<int>& exclude,
                    const PosChargeArray& particlej,
                    const AtomtypeArray& atomja,
		    const AtomIDArray& atomidj,
                    const TypeRange& typerange,
                    const double cutoff2,
                    const double margin2,
                    const bool self,
                    const int i,
                    const int offset)
{
  bool success=true;

  int np=(*npair) + offset;

  double cut2 = cutoff2 + margin2;
  int j;
  int jmin, jmax;
  jmin = typerange.begin;
  jmax = typerange.end;
  for(j=jmin;j<jmax;j++){
    if((self==false)||(i!=j)){
      if(self){
	std::vector<int>::const_iterator it = std::find(exclude.begin(),exclude.end(),atomidj[j]);
	if(it!=exclude.end())continue;
      }
      Position d = particlei.position - particlej[j].position;
      double r2 = d.norm2();
      if(r2<cut2){
        pid[np] = j;
        int atomj = atomja[j];
        lj[np] = atomj;
        np++;
        if(np>=MAX_PAIR){
          printf("number of pair reach to MAX_PAIR %d\n",np);
          break;
        }
      }
    }
  }
  (*npair) = np-offset;
  if(j<typerange.end){
    printf("pair list full\n");
    success = false;
  }
  
  return success;
}

template<typename GPA>
bool makeljpairlist1(int *pid,
                      int *lj,
                      int *npair,
                      const ParticlePosCharge& particlei,
                      const Atomtype& atomi,
                      const GPA& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljpairlist(pid,lj,&np,
                             particlei,
                             atomi,
                             particlej.poscharge,
                             particlej.atomtype,
                             typerange[targetcell],
                             cutoff2, margin2,self,i);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}
template<typename GPA>
bool makeljcpairlist1(int *pid,
                      int *lj,
                      int *npair,
                      const ParticlePosCharge& particlei,
                      const Atomtype& atomi,
                      const GPA& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljcpairlist(pid,lj,&np,
                             particlei,
                             atomi,
                             particlej.poscharge,
                             particlej.atomtype,
                             typerange[targetcell],
                             cutoff2, margin2,self,i);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}
template<typename GPA>
bool makeljcpairlist1(int *pid,
                      int *lj,
                      int *npair,
                      const ParticlePosCharge& particlei,
                      const Atomtype& atomi,
		      const WHPair& waterh,
                      const GPA& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljcpairlist(pid,lj,&np,
                             particlei,
                             atomi,
			     waterh,
                             particlej.poscharge,
                             particlej.atomtype,
                             typerange[targetcell],
                             cutoff2, margin2,self,i);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}
template<typename GPA>
bool makeljcpairlist1(int *pid,
                      int *lj,
                      int *npair,
                      const ParticlePosCharge& particlei,
                      const Atomtype& atomi,
		      const std::vector<int>& exclude,
                      const GPA& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljcpairlist(pid,lj,&np,
                             particlei,
                             atomi,
			     exclude,
                             particlej.poscharge,
                             particlej.atomtype,
			     particlej.atomid,
                             typerange[targetcell],
                             cutoff2, margin2,self,i);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}

template<typename PA, typename GPA>
bool makeljpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const PA& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const GPA& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].lj.begin;
    imax = typerangei[ci].ljcoulomb.end;
    inum = imax-imin;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
    for(i=0;i<inum;i++){
      success = makeljpairlist1(pid[lnum+i], lj[lnum+i],
                                 &(npair[lnum+i]),
                                 particlei.poscharge[i+imin],
                                 particlei.atomtype[i+imin],
                                 particlej,
                                 typerangej,
                                 shorttarget_index[ci],
                                 cutoff2,
                                 margin2,self,i+imin);
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
    }
    lnum += inum;
  }
  *npl = lnum;
  return success;
}
template<typename PA, typename GPA>
bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const PA& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const GPA& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    inum = imax-imin;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
    for(i=0;i<inum;i++){
      success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
                                 &(npair[lnum+i]),
                                 particlei.poscharge[i+imin],
                                 particlei.atomtype[i+imin],
                                 particlej,
                                 typerangej,
                                 shorttarget_index[ci],
                                 cutoff2,
                                 margin2,self,i+imin);
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
    }
    lnum += inum;
#ifdef OVERLAP
    if(self){
      MPICHECKER::mpi_checker();
    }
#endif
  }
  *npl = lnum;
  return success;
}
template<typename PA, typename GPA>
bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const PA& particlei,
                     const std::vector<TypeRange>& typerangei,
		     const WaterList& waterlist,
                     const GPA& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    inum = imax-imin;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
    for(i=0;i<inum;i++){
      bool inwaterlist=false;
      WHPair waterh;
      if(particlei.atomtype[i+imin]==ATOMTYPE_WO){
	WaterList::const_iterator wlit = waterlist.find(i+imin);
	if(wlit!=waterlist.end()){
	  inwaterlist=true;
	  waterh = wlit->second;
	}
      }else if(particlei.atomtype[i+imin]==ATOMTYPE_WH){
	std::map<int, int>::const_iterator wlrit = waterlist.reverse_list.find(i+imin);
	if(wlrit!=waterlist.reverse_list.end()){
	  AtomID wo = wlrit->second;
	  WaterList::const_iterator wlit = waterlist.find(wo);
	  if(wlit!=waterlist.end()){
	    inwaterlist=true;
	    waterh.h1 = wo;
	    if(wlit->second.h1==i+imin){
	      waterh.h2 = wlit->second.h2;
	    }else if(wlit->second.h2==i+imin){
	      waterh.h2 = wlit->second.h1;
	    }else{
	      inwaterlist=false;
	    }
	  }
	}
      }
      if(inwaterlist){
	success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
				   &(npair[lnum+i]),
				   particlei.poscharge[i+imin],
				   particlei.atomtype[i+imin],
				   waterh,
				   particlej,
				   typerangej,
				   shorttarget_index[ci],
				   cutoff2,
				   margin2,self,i+imin);
      }else{
	success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
				   &(npair[lnum+i]),
				   particlei.poscharge[i+imin],
				   particlei.atomtype[i+imin],
				   particlej,
				   typerangej,
				   shorttarget_index[ci],
				   cutoff2,
				   margin2,self,i+imin);
      }
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
    }
    lnum += inum;
#ifdef OVERLAP
    if(self){
      MPICHECKER::mpi_checker();
    }
#endif
  }
  *npl = lnum;
  return success;
}

template<typename PA, typename GPA>
bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const PA& particlei,
                     const std::vector<TypeRange>& typerangei,
		     const WaterList& waterlist,
		     const std::vector<int>& excludei,
		     const std::vector< std::vector<int> >& excludej,
                     const GPA& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    inum = imax-imin;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
    for(i=0;i<inum;i++){
      bool inwaterlist=false;
      WHPair waterh;
      if(particlei.atomtype[i+imin]==ATOMTYPE_WO){
	WaterList::const_iterator wlit = waterlist.find(i+imin);
	if(wlit!=waterlist.end()){
	  inwaterlist=true;
	  waterh = wlit->second;
	}
      }else if(particlei.atomtype[i+imin]==ATOMTYPE_WH){
	std::map<int, int>::const_iterator wlrit = waterlist.reverse_list.find(i+imin);
	if(wlrit!=waterlist.reverse_list.end()){
	  AtomID wo = wlrit->second;
	  WaterList::const_iterator wlit = waterlist.find(wo);
	  if(wlit!=waterlist.end()){
	    inwaterlist=true;
	    waterh.h1 = wo;
	    if(wlit->second.h1==i+imin){
	      waterh.h2 = wlit->second.h2;
	    }else if(wlit->second.h2==i+imin){
	      waterh.h2 = wlit->second.h1;
	    }else{
	      inwaterlist=false;
	    }
	  }
	}
      }
      if(inwaterlist){
	success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
				   &(npair[lnum+i]),
				   particlei.poscharge[i+imin],
				   particlei.atomtype[i+imin],
				   waterh,
				   particlej,
				   typerangej,
				   shorttarget_index[ci],
				   cutoff2,
				   margin2,self,i+imin);
      }else{
	int iex;
	for(iex=0;iex<excludei.size();iex++){
	  if(particlei.atomid[i+imin]==excludei[iex])break;
	}
	if(iex<excludei.size()){
	  success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
				     &(npair[lnum+i]),
				     particlei.poscharge[i+imin],
				     particlei.atomtype[i+imin],
				     excludej[iex],
				     particlej,
				     typerangej,
				     shorttarget_index[ci],
				     cutoff2,
				     margin2,self,i+imin);
	}else{
	  success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
				     &(npair[lnum+i]),
				     particlei.poscharge[i+imin],
				     particlei.atomtype[i+imin],
				     particlej,
				     typerangej,
				     shorttarget_index[ci],
				     cutoff2,
				     margin2,self,i+imin);
	}
      }
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
    }
    lnum += inum;
  }
  *npl = lnum;
  return success;
}

template<typename GPA>
bool makeljpairlist1(int *pid,
                      int *lj,
                      int *npair,
                      const ParticlePosCharge& particlei,
                      const Atomtype& atomi,
                      const GPA& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
                      const int selfnpair,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljpairlist(pid,lj,&np,
                             particlei,
                             atomi,
                             particlej.poscharge,
                             particlej.atomtype,
                             typerange[targetcell],
                             cutoff2, margin2,self,i,
                             selfnpair);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}
template<typename GPA>
bool makeljcpairlist1(int *pid,
                      int *lj,
                      int *npair,
                      const ParticlePosCharge& particlei,
                      const Atomtype& atomi,
                      const GPA& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
                      const int selfnpair,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljcpairlist(pid,lj,&np,
                             particlei,
                             atomi,
                             particlej.poscharge,
                             particlej.atomtype,
                             typerange[targetcell],
                             cutoff2, margin2,self,i,
                             selfnpair);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}
template<typename GPA>
bool makeljcpairlist1(int *pid,
                      int *lj,
                      int *npair,
                      const ParticlePosCharge& particlei,
                      const Atomtype& atomi,
		      const WHPair& waterh,
                      const GPA& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
                      const int selfnpair,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljcpairlist(pid,lj,&np,
                             particlei,
                             atomi,
			     waterh,
                             particlej.poscharge,
                             particlej.atomtype,
                             typerange[targetcell],
                             cutoff2, margin2,self,i,
                             selfnpair);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}

template<typename PA, typename GPA>
bool makeljpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const PA& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const GPA& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const int *selfnpair,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].lj.begin;
    imax = typerangei[ci].ljcoulomb.end;
    inum = imax-imin;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
    for(i=0;i<inum;i++){
      success = makeljpairlist1(pid[lnum+i], lj[lnum+i],
                                 &(npair[lnum+i]),
                                 particlei.poscharge[i+imin],
                                 particlei.atomtype[i+imin],
                                 particlej,
                                 typerangej,
                                 shorttarget_index[ci],
                                 selfnpair[lnum+i],
                                 cutoff2,
                                 margin2,self,i+imin);
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
    }
    lnum += inum;
  }
  *npl = lnum;
  return success;
}
template<typename PA, typename GPA>
bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const PA& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const GPA& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const int *selfnpair,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    inum = imax-imin;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
    for(i=0;i<inum;i++){
      success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
                                 &(npair[lnum+i]),
                                 particlei.poscharge[i+imin],
                                 particlei.atomtype[i+imin],
                                 particlej,
                                 typerangej,
                                 shorttarget_index[ci],
                                 selfnpair[lnum+i],
                                 cutoff2,
                                 margin2,self,i+imin);
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
    }
    lnum += inum;
  }
  *npl = lnum;
  return success;
}
template<typename PA, typename GPA>
bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const PA& particlei,
                     const std::vector<TypeRange>& typerangei,
		     const WaterList& waterlist,
                     const GPA& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const int *selfnpair,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    inum = imax-imin;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
    for(i=0;i<inum;i++){
      bool inwaterlist=false;
      WHPair waterh;
      if(particlei.atomtype[i+imin]==ATOMTYPE_WO){
	WaterList::const_iterator wlit = waterlist.find(i+imin);
	if(wlit!=waterlist.end()){
	  inwaterlist=true;
	  waterh = wlit->second;
	}
      }else if(particlei.atomtype[i+imin]==ATOMTYPE_WH){
	std::map<int, int>::const_iterator wlrit = waterlist.reverse_list.find(i+imin);
	if(wlrit!=waterlist.reverse_list.end()){
	  AtomID wo = wlrit->second;
	  WaterList::const_iterator wlit = waterlist.find(wo);
	  if(wlit!=waterlist.end()){
	    inwaterlist=true;
	    waterh.h1 = wo;
	    if(wlit->second.h1==i+imin){
	      waterh.h2 = wlit->second.h2;
	    }else if(wlit->second.h2==i+imin){
	      waterh.h2 = wlit->second.h1;
	    }else{
	      inwaterlist=false;
	    }
	  }
	}
      }
      if(inwaterlist){
	success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
				   &(npair[lnum+i]),
				   particlei.poscharge[i+imin],
				   particlei.atomtype[i+imin],
				   waterh,
				   particlej,
				   typerangej,
				   shorttarget_index[ci],
				   selfnpair[lnum+i],
				   cutoff2,
				   margin2,self,i+imin);
      }else{
	success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
				   &(npair[lnum+i]),
				   particlei.poscharge[i+imin],
				   particlei.atomtype[i+imin],
				   particlej,
				   typerangej,
				   shorttarget_index[ci],
				   selfnpair[lnum+i],
				   cutoff2,
				   margin2,self,i+imin);
      }
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
    }
    lnum += inum;
  }
  *npl = lnum;
  return success;
}

template
bool makeljpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const CombinedParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const CombinedParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self);
template
bool makeljpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const CombinedParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const GhostParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self);
template
bool makeljpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const CombinedParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const GhostParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const int *selfnpair,
                     const double cutoff2,
                     const double margin2,
                     const bool self);

template
bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const CombinedParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const CombinedParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self);
template
bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const CombinedParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
		     const WaterList& waterlist,
                     const CombinedParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self);
template
bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const CombinedParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const GhostParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self);
template
bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const CombinedParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const GhostParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const int *selfnpair,
                     const double cutoff2,
                     const double margin2,
                     const bool self);

bool addljcpairlist(int *pid,
                    int *lj,
                    int *npair,
                    const Particle& particlei,
                    const ParticleArray& particlej,
                    const TypeRange& typerange,
                    const double cutoff2,
                    const double margin2,
                    const bool self,
                    const int i,
		    const int offset)
{
  bool success=true;

  int np=(*npair)+offset;

  //  int atomi = particlei.atomtype;
  double cut2 = cutoff2 + margin2;
  int j;
  int jmin, jmax;
  jmin = typerange.begin;
  jmax = typerange.end;
  for(j=jmin;j<jmax;j++){
    if((self==false)||(i!=j)){
      Position d = particlei.position - particlej[j].position;
      double r2 = d.norm2();
      if(r2<cut2){
        pid[np] = j;
        int atomj = particlej[j].atomtype;
        lj[np] = atomj;
        np++;
        if(np==MAX_PAIR){
          printf("number of pair reach to MAX_PAIR %d\n",np);
          break;
        }
      }
    }
  }
  (*npair) = np-offset;
  if(j<typerange.end){
    printf("pair list full\n");
    success = false;
  }
  
  return success;
}

bool makeljcpairlist1(int *pid,
                      int *lj,
                      int *npair,
                      const Particle& particlei,
                      const ParticleArray& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljcpairlist(pid,lj,&np,
                             particlei,
                             particlej,
                             typerange[targetcell],
                             cutoff2, margin2,self,i);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}

bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const ParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const ParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    inum = imax-imin;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
    for(i=0;i<inum;i++){
      success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
                                 &(npair[lnum+i]),
                                 particlei[i+imin],
                                 particlej,
                                 typerangej,
                                 shorttarget_index[ci],
                                 cutoff2,
                                 margin2,self,i+imin);
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
    }
    lnum += inum;
  }
  *npl = lnum;
  return success;
}

bool makeljcpairlist1(int *pid,
                      int *lj,
                      int *npair,
                      const Particle& particlei,
                      const ParticleArray& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
		      const int selfnpair,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljcpairlist(pid,lj,&np,
                             particlei,
                             particlej,
                             typerange[targetcell],
                             cutoff2, margin2,self,i,
			     selfnpair);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}

bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     int (*lj)[MAX_PAIR],
                     int *npair,
                     int *iid,
                     int *npl,
                     const ParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const ParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
		     const int *selfnpair,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    inum = imax-imin;
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
    for(i=0;i<inum;i++){
      success = makeljcpairlist1(pid[lnum+i], lj[lnum+i],
                                 &(npair[lnum+i]),
                                 particlei[i+imin],
                                 particlej,
                                 typerangej,
                                 shorttarget_index[ci],
				 selfnpair[lnum+i],
                                 cutoff2,
                                 margin2,self,i+imin);
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
    }
    lnum += inum;
  }
  *npl = lnum;
  return success;
}

// hold value of pos, charge
bool addljcpairlist(int *pid,
                    double (*pos)[3],
                    double *charge,
                    double (*lj)[4],
                    int *npair,
                    const Particle& particlei,
                    const ParticleArray& particlej,
                    const TypeRange& typerange,
                    const double cutoff2,
                    const double margin2,
                    const bool self,
                    const int i,
		    const int offset)
{
  bool success=true;

  int np=(*npair)+offset;

  int atomi = particlei.atomtype;
  double cut2 = cutoff2 + margin2;
  int j;
  int jmin, jmax;
  jmin = typerange.begin;
  jmax = typerange.end;
  for(j=jmin;j<jmax;j++){
    if((self==false)||(i!=j)){
      Position d = particlei.position - particlej[j].position;
      double r2 = d.norm2();
      if(r2<cut2){
        pid[np] = j;
        pos[np][0] =  particlej[j].position.x;
        pos[np][1] =  particlej[j].position.y;
        pos[np][2] =  particlej[j].position.z;
        charge[np] =  particlej[j].charge;
        int atomj = particlej[j].atomtype;
        lj[np][0] = ShortRange::ljmixparameters[atomi][atomj].potFirst;
        lj[np][1] = ShortRange::ljmixparameters[atomi][atomj].potSecond;
        lj[np][2] = ShortRange::ljmixparameters[atomi][atomj].forceFirst;
        lj[np][3] = ShortRange::ljmixparameters[atomi][atomj].forceSecond;
        np++;
        if(np==MAX_PAIR){
          printf("number of pair reach to MAX_PAIR %d\n",np);
          break;
        }
      }
    }
  }
  (*npair) = np-offset;
  if(j<typerange.end){
    printf("pair list full\n");
    success = false;
  }
  
  return success;
}

bool makeljcpairlist1(int *pid,
                      double (*pos)[3],
                      double *charge,
                      double (*lj)[4],
                      int *npair,
                      const Particle& particlei,
                      const ParticleArray& particlej,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& shorttarget_index,
                      const double cutoff2,
                      const double margin2,
                      const bool self,
                      const int i)
{
  bool success = true;
  int np=0;
  int cjmax = shorttarget_index.size();
  for(int cj=0;cj<cjmax;cj++){
    int targetcell = shorttarget_index[cj];
    success = addljcpairlist(pid,pos,charge,lj,&np,
                             particlei,
                             particlej,
                             typerange[targetcell],
                             cutoff2, margin2,self,i);
    if(success==false){
      printf("pair list full at targetcell %d\n",targetcell);
      break;
    }
  }
  (*npair) = np;
  return success;
}

bool makeljcpairlist(int (*pid)[MAX_PAIR],
                     double (*pos)[MAX_PAIR][3],
                     double (*charge)[MAX_PAIR],
                     double (*lj)[MAX_PAIR][4],
                     int *npair,
                     int *iid,
                     double (*posi)[3],
                     double *chargei,
                     int *npl,
                     const ParticleArray& particlei,
                     const std::vector<TypeRange>& typerangei,
                     const ParticleArray& particlej,
                     const std::vector<TypeRange>& typerangej,
                     const std::vector<std::vector<int> >& shorttarget_index,
                     const double cutoff2,
                     const double margin2,
                     const bool self)
{
  bool success=true;
  int lnum=0;
  int cimax = typerangei.size();
  int ci;
#ifdef SMALL_CELL
  int *loffset;
  loffset = new int[cimax+1];
  loffset[0] = 0;
  for(ci=0;ci<cimax;ci++){
    loffset[ci+1] = loffset[ci] + (typerangei[ci].end - typerangei[ci].begin);
  }
#ifdef _OPENMP
#pragma omp parallel for private(ci,lnum) schedule(runtime)
#endif
#endif
  for(ci=0;ci<cimax;ci++){
    int i;
    int imin, imax, inum;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    inum = imax-imin;
#ifdef SMALL_CELL
    lnum = loffset[ci];
#else
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(runtime)
#endif    
#endif
    for(i=0;i<inum;i++){
      success = makeljcpairlist1(pid[lnum+i],
                                 pos[lnum+i],charge[lnum+i],lj[lnum+i],
                                 &(npair[lnum+i]),
                                 particlei[i+imin],
                                 particlej,
                                 typerangej,
                                 shorttarget_index[ci],
                                 cutoff2,
                                 margin2,self,i+imin);
      if(success==false){
        printf("pair list full for cell %d particle %d\n",ci,i);
      }
      iid[lnum+i] = i+imin;
      posi[lnum+i][0] = particlei[i+imin].position.x;
      posi[lnum+i][1] = particlei[i+imin].position.y;
      posi[lnum+i][2] = particlei[i+imin].position.z;
      chargei[lnum+i] = particlei[i+imin].charge;
    }
#ifndef SMALL_CELL
    lnum += inum;
#endif
  }
#ifdef SMALL_CELL
  *npl = loffset[cimax];
#else
  *npl = lnum;
#endif
  return success;
}

template<typename PA>
void updatepairlistpos1(const int *pid,
                        double (*pos)[3],
                        const int npair,
                        const PA& particlej)
{
  for(int p=0;p<npair;p++){
    int j=pid[p];
    pos[p][0] = particlej.poscharge[j].position.x;
    pos[p][1] = particlej.poscharge[j].position.y;
    pos[p][2] = particlej.poscharge[j].position.z;
  }
}
template<>
void updatepairlistpos1<ParticleArray>(const int *pid,
                                       double (*pos)[3],
                                       const int npair,
                                       const ParticleArray& particlej)
{
  for(int p=0;p<npair;p++){
    int j=pid[p];
    pos[p][0] = particlej[j].position.x;
    pos[p][1] = particlej[j].position.y;
    pos[p][2] = particlej[j].position.z;
  }
}

template<typename PA>
void updatepairlistpos(const int (*pid)[MAX_PAIR],
                       double (*pos)[MAX_PAIR][3],
                       const int *npair,
                       const int npl,
                       const PA& particlej)
{
  int pl;
#ifdef _OPENMP
#pragma omp parallel for private(pl) schedule(runtime)
#endif    
  for(pl=0;pl<npl;pl++){
    updatepairlistpos1(pid[pl],pos[pl],npair[pl],particlej);
  }
}
template
void updatepairlistpos(const int (*pid)[MAX_PAIR],
                       double (*pos)[MAX_PAIR][3],
                       const int *npair,
                       const int npl,
                       const ParticleArray& particlej);

void importpairforce(ForceArray& forcei,
                     const int *iid,
                     const double (*pforce)[3],
                     const int npl)
{
  for(int p=0;p<npl;p++){
    int i = iid[p];
    forcei[i].x += pforce[p][0];
    forcei[i].y += pforce[p][1];
    forcei[i].z += pforce[p][2];
  }
}
