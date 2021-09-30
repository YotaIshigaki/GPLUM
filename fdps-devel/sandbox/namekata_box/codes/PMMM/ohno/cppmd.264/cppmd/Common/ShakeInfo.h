#ifndef SHAKEINFO_H
#define SHAKEINFO_H

#include <cstdio>
#include "ParticleInfo.h"
#include "CovalentBondInfo.h"


/// It's correct only pa[i].atomid==i .
template<class PA>
void append_shake(ShakeList& shakelist, const PA& pa, const CovalentBondList& cb) 
{
  int id, hid;
  for(std::vector<CovalentBondInfo::Bond>::size_type b=0;b<cb.bond.size();b++){
    int is_h=-1;
    if(getatomtype(pa,cb.bond[b].id_of_atom[1])<ATOMTYPE_WH){
      hid = cb.bond[b].id_of_atom[0];
      id = cb.bond[b].id_of_atom[1];
      is_h = 1;
    }else if(getatomtype(pa,cb.bond[b].id_of_atom[0])<ATOMTYPE_WH){
      hid = cb.bond[b].id_of_atom[1];
      id = cb.bond[b].id_of_atom[0];
      is_h = 0;
    }
    if(is_h!=-1){
      shakelist.add_shake(hid, id, cb.bond[b].typeofbond);
    }
  }
}

template<class PA>
void append_shake_offset(ShakeList& shakelist, const PA& pa, const CovalentBondList& cb) 
{
  size_t np = pa.size();
  int offset = 0;
  for(std::vector<CovalentBondInfo::Bond>::size_type b=0;b<cb.bond.size();b++){
    AtomID index0, index1;
    int is_h=-1;
    AtomID tid0 = cb.bond[b].id_of_atom[0];
    AtomID tid1 = cb.bond[b].id_of_atom[1];
    AtomID i;
    AtomID idx, hidx;
    for(i=0;i<np;i++){
      index0 = (tid0+i+offset)%np;
      if(getatomid(pa,index0)==tid0)break;
    }
    if(i>=np){
      printf("atom id %d is not found\n",tid0);
    }else{
      if(offset!=index0-tid0){
	offset = index0 - tid0;
	printf("offset %d at %ld %ld for bond %ld 0\n",
         offset,static_cast<long>(i),static_cast<long>(tid0),b);
      }
      for(i=0;i<np;i++){
	index1 = (tid1+i+offset)%np;
	if(getatomid(pa,index1)==tid1)break;
      }
      if(i>=np){
	printf("atom id %d is not found\n",tid1);
      }else{
	if(offset!=index1-tid1){
	  offset = index1 - tid1;
	  printf("offset %d at %ld %ld for bond %ld 1\n",
           offset,static_cast<long>(i),static_cast<long>(tid1),b);
	}
	if(getatomtype(pa,index1)<ATOMTYPE_WH){
	  hidx = index0;
	  idx = index1;
	  is_h = 1;
	}else if(getatomtype(pa,index0)<ATOMTYPE_WH){
	  hidx = index1;
	  idx = index0;
	  is_h = 0;
	}
	if(is_h!=-1){
	  shakelist.add_shake(hidx, idx, cb.bond[b].typeofbond);
	}
      }
    }
  }
}

template<class PA>
void append_shake_offset_with_mark(ShakeList& shakelist, const PA& pa, const CovalentBondList& cb) 
{
  size_t np = pa.size();
  int offset = 0;
  for(std::vector<CovalentBondInfo::Bond>::size_type b=0;b<cb.bond.size();b++){
    if(cb.bond[b].shake>0){
      AtomID index0, index1;
      int is_h=-1;
      AtomID tid0 = cb.bond[b].id_of_atom[0];
      AtomID tid1 = cb.bond[b].id_of_atom[1];
      AtomID idx, hidx;
      index0 = tid0 + offset;
      if((index0<0)||(index0>=np)||(getatomid(pa,index0)!=tid0)){
	for(index0=0;index0<np;index0++){
	  if(getatomid(pa,index0)==tid0)break;
	}
      }
      if(index0>=np){
	printf("atom id %d is not found for bond %ld 0 (last index %d)\n",tid0,b,index0);
      }else{
	if(offset!=index0-tid0){
	  offset = index0 - tid0;
	  //	  printf("offset %d at %ld %ld for bond %ld 0\n",offset,i,tid0,b);
	}
	index1 = tid1 + offset;
	if((index1<0)||(index1>=np)||(getatomid(pa,index1)!=tid1)){
	  for(index1=0;index1<np;index1++){
	    if(getatomid(pa,index1)==tid1)break;
	  }
	}
	if(index1>=np){
	  printf("atom id %d is not found for bond %ld 1 (last index %d)\n",tid1,b,index1);
	}else{
	  if(offset!=index1-tid1){
	    int oo = offset;
	    offset = index1 - tid1;
	    printf("offset %d at %d %d for bond %ld 1, changed from %d for 0\n",offset,index1,tid1,b,oo);   
	  }
	  if(getatomtype(pa,index1)<ATOMTYPE_WH){
	    hidx = index0;
	    idx = index1;
	    is_h = 1;
	  }else if(getatomtype(pa,index0)<ATOMTYPE_WH){
	    hidx = index1;
	    idx = index0;
	    is_h = 0;
	  }
	  if(is_h!=-1){
	    shakelist.add_shake(hidx, idx, cb.bond[b].typeofbond);
	  }
	}
      }
    }
  }
}

template<class PA>
void append_shake_scanid(ShakeList& shakelist, const PA& pa, const CovalentBondList& cb) 
{
  size_t np = pa.size();
  for(std::vector<CovalentBondInfo::Bond>::size_type b=0;b<cb.bond.size();b++){
    AtomID index0, index1;
    int is_h=-1;
    AtomID tid0 = cb.bond[b].id_of_atom[0];
    AtomID tid1 = cb.bond[b].id_of_atom[1];
    AtomID idx, hidx;
    for(index0=0;index0<np;index0++){
      if(getatomid(pa,index0)==tid0)break;
    }
    if(index0>=np){
      printf("atom id %d is not found\n",tid0);
    }else{
      for(index1=0;index1<np;index1++){
	if(getatomid(pa,index1)==tid1)break;
      }
      if(index1>=np){
	printf("atom id %d is not found\n",tid1);
      }else{
	if(getatomtype(pa,index1)<ATOMTYPE_WH){
	  hidx = index0;
	  idx = index1;
	  is_h = 1;
	}else if(getatomtype(pa,index0)<ATOMTYPE_WH){
	  hidx = index1;
	  idx = index0;
	  is_h = 0;
	}
	if(is_h!=-1){
	  shakelist.add_shake(hidx, idx, cb.bond[b].typeofbond);
	}
      }
    }
  }
}

template<class PA>
void append_shake_scanid_with_mark(ShakeList& shakelist, const PA& pa, const CovalentBondList& cb) 
{
  size_t np = pa.size();
  for(std::vector<CovalentBondInfo::Bond>::size_type b=0;b<cb.bond.size();b++){
    if(cb.bond[b].shake>0){
      AtomID index0, index1;
      int is_h=-1;
      AtomID tid0 = cb.bond[b].id_of_atom[0];
      AtomID tid1 = cb.bond[b].id_of_atom[1];
      AtomID idx, hidx;
      for(index0=0;index0<np;index0++){
	if(getatomid(pa,index0)==tid0)break;
      }
      if(index0>=np){
	printf("atom id %d is not found\n",tid0);
      }else{
	for(index1=0;index1<np;index1++){
	  if(getatomid(pa,index1)==tid1)break;
	}
	if(index1>=np){
	  printf("atom id %d is not found\n",tid1);
	}else{
	  if(getatomtype(pa,index1)<ATOMTYPE_WH){
	    hidx = index0;
	    idx = index1;
	    is_h = 1;
	  }else if(getatomtype(pa,index0)<ATOMTYPE_WH){
	    hidx = index1;
	    idx = index0;
	    is_h = 0;
	  }
	  if(is_h!=-1){
	    shakelist.add_shake(hidx, idx, cb.bond[b].typeofbond);
	  }
	}
      }
    }
  }
}

template<class PA>
void append_shake1(ShakeList& shakelist, const PA& pa, const CovalentBondList& cb) {
  int i_ha = 1;  // if pa[0] = H atom, pa[1] = Heavy atom
  for(size_t i = 0; i < pa.size(); i++){
    if(getatomtype(pa,i) == ATOMTYPE_WO) { /// It's correct Water orger WO WH WH, ex AmberPrmtop
      i+=2; // WH
    }else if(getatomtype(pa,i) < ATOMTYPE_WH){
      AtomID haid=-1;
      AtomID aid = getatomid(pa,i);
      int typeofbond;
      for(std::vector<CovalentBondInfo::Bond>::size_type b=0;b<cb.bond.size();b++){
	if(cb.bond[b].id_of_atom[1]==aid){
	  haid = cb.bond[b].id_of_atom[0];
	  typeofbond = cb.bond[b].typeofbond;
	  break;
	}
	if(cb.bond[b].id_of_atom[0]==aid){
	  haid = cb.bond[b].id_of_atom[1];
	  typeofbond = cb.bond[b].typeofbond;
	  break;
	}
      }
      if(haid==-1){
	printf("atom id %d is not found in CovalentBond Bond\n",aid);
      }else{
	i_ha=-1;
	if((i<pa.size()-1)&&(pa[i+1].atomid==haid)){
	  i_ha=i+1;
	}else{
	  int j;
	  for(j=i-1;j>=0;j--){
	    if(getatomid(pa,j)==haid){
	      i_ha = j;
	      break;
	    }
	  }
	  if(i_ha==-1){
	    for(j=i+2;j<pa.size();j++){
	      if(getatomid(pa,j)==haid){
		i_ha = j;
		break;
	      }
	    }
	  }
	}
	if(i_ha==-1){
	  printf("atom id %d (heavy atom for id %d, index %d, atomtype %d) is not found\n",haid,aid,i,pa[i].atomtype);
	}else{
	  shakelist.add_shake(i_ha, i, typeofbond);
	}
      }
    }else{
      i_ha = i;
    }
  }
}

#endif
