// Copyright (c) 2008, 2009, 2010, 2011 RIKEN. All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//    1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//    2. Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//    3. The name of the author may not be used to endorse or promote
//    products derived from this software without specific prior written
//    permission.
//
// THIS SOFTWARE IS PROVIDED BY RIKEN ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H
#ifdef _OPENMP
#include <omp.h>
#endif

#include <mpi.h>

#include "cppmd/simulator.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <algorithm>
#include <limits>
#include <map>
#include <vector>
#include <sstream>

#include "cppmd/cppmd-int.h"

#include "Config.h"

#include "Common.h"
#include "LJAmber94.h"
#include "LJAr.h"
#include "CellIndex.h"
#include "CubicCell.h"
#include "HalfShell.h"
#include "NoInputSystem.h"
#include "ArFCC.h"
#include "NaClFCC.h"
#include "WaterLattice.h"
#include "CalculationUnit.h"
#include "CalcPreparator.h"
#include "UnitParameter.h"
#include "Timer.h"

#include "Dump.h"

#include "AmberFile/AmberFile.h"

#include "ShakeInfo.h"

#include "Log.h"

#ifdef CPPMD_ENABLE_PMMM
# ifdef PMMM_QPOS_OUT
#include "qpos_io.h"
# endif
#endif

#if 0
#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#endif

namespace {
  
  template<class PA>
  void mark_shake(const PA& particle, CovalentBondList& cb, std::vector<int>& heavyid)
  {
    int b;
    for(b=0;b<cb.bond.size();b++){
      int id0 = cb.bond[b].id_of_atom[0];
      int id1 = cb.bond[b].id_of_atom[1];
      if(getatomid(particle,id0)!=id0){
	printf("AtomID of particle %d is not %d\n",id0,id0);
      }
      if(getatomid(particle,id1)!=id1){
	printf("AtomID of particle %d is not %d\n",id1,id1);
      }
      if((getatomtype(particle,id0)<ATOMTYPE_WH)||
	 (getatomtype(particle,id1)<ATOMTYPE_WH)){
	cb.bond[b].shake = 1;
	if(getatomtype(particle,id0)<ATOMTYPE_WH){
	  heavyid[id0] = id1;
	}else{
	  heavyid[id1] = id0;
	}
      }
    }
  }

  // for COPY_MODE
  void ShiftAtomIdOfCovalentBondList(const std::map<int,int>& idmap,
				     CovalentBondList* cb) {
    if ((DebugLog::verbose > 1)) {
      printf("ShiftCovalentBond()\n");
    }
    typedef std::vector<CovalentBondInfo::Bond> Bonds;
    typedef std::vector<CovalentBondInfo::Angle> Angles;
    typedef std::vector<CovalentBondInfo::Torsion> Torsions;
    typedef std::vector<CovalentBondInfo::Improper> Impropers;
    typedef std::map<int,int> Idmap;

    Bonds& bonds = cb->bond;
    for (Bonds::size_type i = 0; i < bonds.size(); ++i) {
      Idmap::const_iterator idpair = idmap.find(bonds[i].get_assign_atomid());
      if (idpair != idmap.end()) {
	int ofst = idpair->second - idpair->first;
	bonds[i].id_of_atom[0] += ofst;
	bonds[i].id_of_atom[1] += ofst;
      }
    }

    Angles& angles = cb->angle;
    for (Angles::size_type i = 0; i < angles.size(); ++i) {
      Idmap::const_iterator idpair = idmap.find(angles[i].get_assign_atomid());
      if (idpair != idmap.end()) {
	int ofst = idpair->second - idpair->first;
	angles[i].id_of_atom[0] += ofst;
	angles[i].id_of_atom[1] += ofst;
	angles[i].id_of_atom[2] += ofst;
      }
    }

    Torsions& torsions = cb->torsion;
    for (Torsions::size_type i = 0; i < torsions.size(); ++i) {
      Idmap::const_iterator idpair = idmap.find(torsions[i].get_assign_atomid());
      if (idpair != idmap.end()) {
	int ofst = idpair->second - idpair->first;
	torsions[i].id_of_atom[0] += ofst;
	torsions[i].id_of_atom[1] += ofst;
	torsions[i].id_of_atom[2] += ofst;
	torsions[i].id_of_atom[3] += ofst;
      }
    }

    Impropers& impropers = cb->improper;
    for (Impropers::size_type i = 0; i < impropers.size(); ++i) {
      Idmap::const_iterator idpair = idmap.find(impropers[i].get_assign_atomid());
      if (idpair != idmap.end()) {
	int ofst = idpair->second - idpair->first;
	impropers[i].id_of_atom[0] += ofst;
	impropers[i].id_of_atom[1] += ofst;
	impropers[i].id_of_atom[2] += ofst;
	impropers[i].id_of_atom[3] += ofst;
      }
    }

    //cb->make_atomidset();
  }
  //
  // for COPY_MODE
  // ShiftParticleAndCovalentBond()
  // - shift particlearray and covalent_bond
  // - enlarge box
  // - count number of copies
  // Returns copy number
  //
  // Note: We use GeometryXYZ to calculate 3d vector to index or index to
  // 3d vector conversion.
  //
  template<class PA>
  int ShiftParticleAndCovalentBond(const int short_id,
				   const int num_particle,
				   const SpaceVector<int>& nodediv3d,
				   const SpaceVector<int>& copynum3d,
				   SpaceVector<double>* boxsize,
				   PA* pa,
				   CovalentBondList* cb,
				   std::vector<int>& heavyid,
				   int &num_copy_local) {
    if ((short_id == 0) || (DebugLog::verbose > 1)) {
      printf("[%d]: ShiftParticleAndCovalentBond()\n", short_id);
    }

    // prepare
    const SpaceVector<double>& bs = *boxsize;  // just an alias
    SpaceVector<int> unused(0, 0, 0);
    GeometryXYZ node_geometry = GeometryXYZ(nodediv3d, /*cheat on ctor*/unused);
    GeometryXYZ copy_geometry = GeometryXYZ(copynum3d, /*cheat on ctor*/unused);

    // numbers for node
    SpaceVector<int> node_pos3d = node_geometry.getAbsolutePosition(short_id);

    // numbers for copy
    // nth_... is a zero base index
    // the node computes copyies of [nth_copy3d_lb .. nth_copy3d_ub)
    const int num_copy = copynum3d.x * copynum3d.y * copynum3d.z;
    SpaceVector<int> nth_copy3d_lb(
				   copynum3d.x * node_pos3d.x / nodediv3d.x,
				   copynum3d.y * node_pos3d.y / nodediv3d.y,
				   copynum3d.z * node_pos3d.z / nodediv3d.z);
    SpaceVector<int> nth_copy3d_ub(
				   (copynum3d.x * (node_pos3d.x + 1) + (nodediv3d.x - 1)) / nodediv3d.x,
				   (copynum3d.y * (node_pos3d.y + 1) + (nodediv3d.y - 1)) / nodediv3d.y,
				   (copynum3d.z * (node_pos3d.z + 1) + (nodediv3d.z - 1)) / nodediv3d.z);

    // print information
    if ((short_id == 0) || (DebugLog::verbose > 1)) {
      std::cout <<
        "short_id: " << short_id <<
        ", node_pos3d: " << node_pos3d <<
        ", copynum3d: " << copynum3d <<
        ", nth_copy3d_lb: " << nth_copy3d_lb <<
        ", nth_copy3d_ub: " << nth_copy3d_ub <<
        ", boxsize: " << *boxsize <<
        ", num_particle: " << num_particle <<
        std::endl;
    }

    PA newpa;
    CovalentBondList newcb;

    // foreach copy in [nth_copy3d_lb .. nth_copy3d_ub)
    SpaceVector<int> nth_copy3d = nth_copy3d_lb;
    for (nth_copy3d[2] = nth_copy3d_lb[2];
	 nth_copy3d[2] < nth_copy3d_ub[2]; ++nth_copy3d[2]) {
      for (nth_copy3d[1] = nth_copy3d_lb[1];
	   nth_copy3d[1] < nth_copy3d_ub[1]; ++nth_copy3d[1]) {
	for (nth_copy3d[0] = nth_copy3d_lb[0];
	     nth_copy3d[0] < nth_copy3d_ub[0]; ++nth_copy3d[0]) {
	  int nth_copy = copy_geometry.getNodeID(nth_copy3d);
	  if ((short_id == 0) || (DebugLog::verbose > 2)) {
	    printf("make %d copy\n",nth_copy);
	  }
	  // copy boundary
	  const SpaceVector<double> shift3d(bs[0] * nth_copy3d[0],
					    bs[1] * nth_copy3d[1],
					    bs[2] * nth_copy3d[2]);
	  const SpaceVector<double> copybox_lb = shift3d;
	  const SpaceVector<double> copybox_ub = bs + shift3d;

	  SpaceVector<int> atomid_ofst(
				       num_particle * nth_copy3d[0],
				       num_particle * nth_copy3d[1] * copynum3d[0],
				       num_particle * nth_copy3d[2] * copynum3d[0] * copynum3d[1]);

	  SpaceVector<int> atomid_corr(
				       num_particle,
				       num_particle * copynum3d[0],
				       num_particle * copynum3d[0] * copynum3d[1]);

	  SpaceVector<int> div(
			       num_particle * copynum3d[0],
			       num_particle * copynum3d[0] * copynum3d[1],
			       num_particle * copynum3d[0] * copynum3d[1] * copynum3d[2]);

	  PA copypa(*pa);  // load pa into copypa
	  CovalentBondList copycb;
	  copycb.append(*cb);         // load cb into copycb
	  std::map<int,int> idmap[2];
	  std::vector<int> originalid((*pa).size());
	  for(int i=0;i<(*pa).size();i++){
	    originalid[i] = getatomid((*pa),i);
	    if(originalid[i]!=i){
	      printf("Original particle[%d] is not %d\n",i,i);
	    }
	  }
	  for (int d = 0; d < /* number of dimension */ 3; ++d) {
	    if ((short_id == 0) || (DebugLog::verbose > 1)) {
	      printf("[%d]: ShiftCovalentBond(%d,%d,%d) dim:%d\n",
		     short_id, nth_copy3d[0], nth_copy3d[1], nth_copy3d[2], d);
	    }

	    double wo_position = 0.0;
	    int num_move_wh = 0;
	    int d_nz = ((d > 0) ? 1 : 0);

	    int index=0;
	    //	    for (ParticleArray::iterator p = copypa.begin();
	    //	 p != copypa.end(); ++p)
	    for(size_t p=0;p<copypa.size();p++)
	      {
		int oldatomid = getatomid(copypa,p);
		getpos(copypa,p)[d] += shift3d[d];
		getatomid(copypa,p) += atomid_ofst[d];
		if ((getatomtype(copypa,p) == ATOMTYPE_WH) && (num_move_wh > 0)) {
		  if (getpos(copypa,p)[d] - wo_position >= bs[d] * 0.5) {
		    getpos(copypa,p)[d] -= bs[d];
		    getatomid(copypa,p) = (getatomid(copypa,p) - atomid_corr[d] + div[d]) % div[d];
		  }
		  else if (getpos(copypa,p)[d] - wo_position < -bs[d] * 0.5) {
		    getpos(copypa,p)[d] += bs[d];
		    getatomid(copypa,p) = (getatomid(copypa,p) + atomid_corr[d]) % div[d];
		  }
		  --num_move_wh;
		}
#ifdef USE_SHAKE
		else if ((getatomtype(copypa,p) < ATOMTYPE_WH)&&(heavyid[index]>=0)) {
		  int orgid=originalid[index];
		  int orghid=heavyid[index];
		  if(getatomid((*pa),orghid)!=orghid){
		    printf("Original particle[%d] is not %d\n",orghid,orghid);
		  }
		  Particle heavy = getparticle(copypa,orghid);
		  if(orghid>orgid){                 // heavy backward, not shifted this dimension
		    heavy.position[d] += shift3d[d];
		    heavy.atomid += atomid_ofst[d];
		    if (heavy.position[d] >= copybox_ub[d]) {
		      heavy.position[d] -= bs[d];
		      heavy.atomid = (heavy.atomid - atomid_corr[d] + div[d]) % div[d];
		    }
		    else if (heavy.position[d] < copybox_lb[d]) {
		      heavy.position[d] += bs[d];
		      heavy.atomid = (heavy.atomid + atomid_corr[d]) % div[d];
		    }
		  }
		  if (getpos(copypa,p)[d] - heavy.position[d] >= bs[d] * 0.5) {
		    getpos(copypa,p)[d] -= bs[d];
		    getatomid(copypa,p) = (getatomid(copypa,p) - atomid_corr[d] + div[d]) % div[d];
		  }
		  else if (getpos(copypa,p)[d] - heavy.position[d] < -bs[d] * 0.5) {
		    getpos(copypa,p)[d] += bs[d];
		    getatomid(copypa,p) = (getatomid(copypa,p) + atomid_corr[d]) % div[d];
		  }
		  {
		    int orgdiff = heavyid[index]-index;
		    int sftdiff = heavy.atomid-getatomid(copypa,p);
		    if(orgdiff!=sftdiff){
		      printf("relativ Shifted heavy AtomID %d, but original %d\n",sftdiff,orgdiff);
		    }
		  }
		}
#endif
		else {
		  if (getpos(copypa,p)[d] >= copybox_ub[d]) {
		    getpos(copypa,p)[d] -= bs[d];
		    getatomid(copypa,p) = (getatomid(copypa,p) - atomid_corr[d] + div[d]) % div[d];
		  }
		  else if (getpos(copypa,p)[d] < copybox_lb[d]) {
		    getpos(copypa,p)[d] += bs[d];
		    getatomid(copypa,p) = (getatomid(copypa,p) + atomid_corr[d]) % div[d];
		  }
		  if (getatomtype(copypa,p) == ATOMTYPE_WO) {
		    wo_position = getpos(copypa,p)[d];
		    num_move_wh = 2;  // number of hydrogen atom in a water molecule
		  }
		}
		idmap[d_nz][oldatomid] = getatomid(copypa,p);  // p->atomid is a new atomid
		index++;
	      }

	    if (d_nz) {
	      for (typename PA::size_type i = 0; i < copypa.size(); ++i)
		idmap[0][i] = idmap[1][idmap[0][i]];
	      idmap[1].clear();
	    }
	  }
	  ShiftAtomIdOfCovalentBondList(idmap[0], &copycb);

	  // append copypa to newpa
	  //	  newpa.insert(newpa.end(), copypa.begin(), copypa.end());
	  append(newpa,copypa);
	  // append copycb to newcb
	  newcb.append(copycb);
	}
      }
    }

    // set new particlearray and covalentbondlist to *pa and *cv
    pa->clear();
    cb->clear();
    pa->swap(newpa);
    cb->swap(newcb);

    // adjust box
    boxsize->x *= copynum3d.x;
    boxsize->y *= copynum3d.y;
    boxsize->z *= copynum3d.z;
  
    num_copy_local = (nth_copy3d_ub[0]-nth_copy3d_lb[0])
      *(nth_copy3d_ub[1]-nth_copy3d_lb[1])
      *(nth_copy3d_ub[2]-nth_copy3d_lb[2]);

    return num_copy;
  }

}  // namespace

//======================================================================

namespace cppmd {

Simulator::Simulator(Config& config) : config_(config) {
}

Simulator::~Simulator() {
}

//----------------------------------------------------------------------
#ifdef OLDPARTICLE
#define PACLASS ParticleArray
#define GPACLASS ParticleArray
#else
#define PACLASS CombinedParticleArray
#define GPACLASS GhostParticleArray
#endif

int Simulator::Run() {
  PerfCounter::init(20);                  // prepare PerfCounter
  int perf_target = PerfCounter::add_target(std::string("simulator"));

  int node_id = 0;
  int num_node = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_node);

  // set num_short
  int num_short = num_node;
  if (config_.md.nodediv3d[0] != 0) {   // when nodediv3d option is given
    num_short =                                     // reset num_short
        config_.md.nodediv3d[0] *
        config_.md.nodediv3d[1] *
        config_.md.nodediv3d[2];
    if (num_short == 0 || num_short > num_node) {   // error exit
      if (node_id == 0)
        printf("invalid nodediv3d given (%s). abort.\n",
               (num_short == 0) ? "num_short == 0" : "num_short > num_node");
      exit(1);
    }
  }

  // set num_long_only and long_plan
  int num_long_only = 0;
  LongRangeMPIPlan long_plan = NoneLong;
  if (config_.md.withlong) {            // when withlong option is set
    bool is_combined = ((config_.md.calc_space & 0x1) == 0);
    int dimension = (config_.md.calc_space / 10);
    switch (dimension) {
      case 3: long_plan = (is_combined ? Combine3D : Separate3D); break;
      case 2: long_plan = (is_combined ? Combine2D : Separate2D); break;
      case 1: long_plan = (is_combined ? Combine1D : Separate1D); break;
      default:
      case 0: long_plan = (is_combined ? Combine0D : Separate0D); break;
    }
    num_long_only = (is_combined ? 0 : (num_node - num_short));
    if (!is_combined && (num_long_only == 0))  // when no long-only node exists
      long_plan = Separate0D;
  }
#if 1
  if (node_id == 0) {                     // print out long plan
    switch (long_plan) {
      case Combine3D:  printf("long plan Combine3D\n"); break;
      case Combine2D:  printf("long plan Combine2D\n"); break;
      case Combine1D:  printf("long plan Combine1D\n"); break;
      case Combine0D:  printf("long plan Combine0D\n"); break;
      case Separate3D: printf("long plan Separate3D\n"); break;
      case Separate2D: printf("long plan Separate2D\n"); break;
      case Separate1D: printf("long plan Separate1D\n"); break;
      case Separate0D: printf("long plan Separate0D\n"); break;
      case NoneLong: default: break;
    }
  }
#endif

  // node and cell geometry
  CellIndexType citype;
  switch (config_.md.citype) {
    case Config::kHalfShell:
      citype = HalfShell;
      break;
    case Config::kFullShell:
      citype = FullShell;
      break;
    case Config::kSmallBall:
      citype = SmallBall;
      break;
    default:
      citype = FullCube;
      break;
  }
#ifdef USE_PAIRLIST
  if(node_id==0){
    printf("PairList support only FullCube. citype force to FullCube\n");
  }
  citype = FullCube;
#endif

  CalcPreparator_::makeoperations(node_id, 
                                  config_.md.withlong, 
                                  config_.md.withbond, 
				  config_.md.withcorrecttrans,
                                  citype,
                                  long_plan, num_long_only);
  int short_id, num_total_set;
  std::vector<int> node_id_of_shorts;
  SpaceVector<int> nodediv3d, celldiv3d_in_node, celldiv3d;
  CalcPreparator_::makegeometry(node_id, 
                                config_.md.nodediv3d,
                                config_.md.celldiv3d_in_node,
                                num_node, num_short,
                                short_id, node_id_of_shorts,
                                num_total_set, config_.md.celldiv,
                                nodediv3d, celldiv3d_in_node, celldiv3d);

  SpaceVector<int> unitnode3d(config_.md.unitnode3d[0],
                              config_.md.unitnode3d[1],
                              config_.md.unitnode3d[2]);
  SpaceVector<int> copynum3d(config_.md.copynum3d[0],
                             config_.md.copynum3d[1],
                             config_.md.copynum3d[2]);
  if ((unitnode3d[0]>0) && (unitnode3d[1]>0) && (unitnode3d[2]>0)) {
    if ((unitnode3d[0]>nodediv3d[0]) || (unitnode3d[1]>nodediv3d[1]) || (unitnode3d[2]>nodediv3d[2]) || (nodediv3d[0]%unitnode3d[0]!=0) || (nodediv3d[1]%unitnode3d[1]!=0) || (nodediv3d[2]%unitnode3d[2]!=0)) {
      printf("invalid unitnode3d\n");
    }else{
      copynum3d[0] = nodediv3d[0]/unitnode3d[0];
      copynum3d[1] = nodediv3d[1]/unitnode3d[1];
      copynum3d[2] = nodediv3d[2]/unitnode3d[2];
    }
  }
  if (node_id == 0){
    printf("Number of Node : %d\n", num_node);
    printf("nodediv3d = ");
    for(int i=0;i<3;i++) printf("%d ",nodediv3d[i]);
    printf("\n");
    printf("celldiv3d = ");
    for(int i=0;i<3;i++) printf("%d ",celldiv3d[i]);
    printf("\n");
    printf("celldiv3d_in_node = ");
    for(int i=0;i<3;i++) printf("%d ",celldiv3d_in_node[i]);
    printf("\n");
    printf("copynum3d = ");
    for(int i=0;i<3;i++) printf("%d ",copynum3d[i]);
    printf("\n");
#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
        int ompnum = omp_get_num_threads();
        printf("Number of OMP threads %d\n",ompnum);
      }
    }
#endif
  }
  if ((DebugLog::verbose = config_.md.verbose) > 0) {
    if (node_id == 0) printf("Verbose  %d\n", DebugLog::verbose);
  }
  if ((DebugLog::particle_dump = config_.md.particle_dump) > 0) {
    if (node_id == 0) printf("Particle_dump  %d\n", DebugLog::particle_dump);
  }

  int num_particle = 0, max_per_cell = 0;
  int num_set = celldiv3d[0] * celldiv3d[1] * celldiv3d[2];

  // lattice
  NoInputSystem<PACLASS>* fcc = NULL;
  switch (config_.md.input_mode) {
    case Config::kARGON:  fcc = new ArFCC<PACLASS>;        break;
    case Config::kNACL:   fcc = new NaClFCC<PACLASS>;      break;
    case Config::kWATER:  fcc = new WaterLattice<PACLASS>; break;
    case Config::kAMBER:  fcc = new WaterLattice<PACLASS>; break;
    default:              fcc = new ArFCC<PACLASS>;        break;  // for now
  }
  SpaceVector<double> boxsize = ((config_.md.bsize < fcc->side.x)
                                 ? fcc->side
                                 : SpaceVector<double>(config_.md.bsize));
  num_particle = config_.md.total_num_particle / num_node;
  max_per_cell = 8;

  CovalentBondList covalent_bond;
  CovalentBondParameterList covalent_bond_parameter_list;

  std::vector<int> heavyid; /// for SHAKE

  bool initialforce = false; 
  //DEUBG
  //  initialforce = true;

  // for COPY_MODE
  bool copy_mode = (copynum3d.x > 1 || copynum3d.y > 1 || copynum3d.z > 1);
  int total_num_particle = 0;
  double total_charge = 0.0;
  int num_copy = 1;
  int num_copy_local = 1;
  double simulation_start_time = 0.0;

  // AMBER file 
  switch (config_.md.input_mode){
    case Config::kAMBER: {  // AMBER
      SpaceVector<double> skew;
      int restart = !(config_.amberfile.restrt.empty());
#ifdef TUNE_NODEIO
# ifndef TUNE_NODEIO_MOD
#  define TUNE_NODEIO_MOD 0
# endif  // TUNE_NODEIO_MOD
      // clock_t inputdata_start_clock = clock();
      double inputdata_start_clock = getrealtime();
      {
        int read_id =
            (TUNE_NODEIO_MOD < 1) ? node_id : (node_id % TUNE_NODEIO_MOD);

        std::string prmtop_file, rescrd_file;

        int bufsiz = std::max(config_.amberfile.prmtop.size(),
                              std::max(config_.amberfile.restrt.size(),
                                       config_.amberfile.inpcrd.size()))
            + std::numeric_limits<int>::digits10 + 1/*NUL*/; // +1(sign)-2(%d)
        {
          char *buf = new char[bufsiz];

          buf[0] = '\0';
          if (std::snprintf(buf, bufsiz,
                            config_.amberfile.prmtop.c_str(),
                            read_id) < 0) {
            // error
          }
          prmtop_file = buf;

          buf[0] = '\0';
          if (std::snprintf(buf, bufsiz,
                            (restart
                             ? config_.amberfile.restrt.c_str()
                             : config_.amberfile.inpcrd.c_str()),
                            read_id) < 0) {
            // error
          }
          rescrd_file = buf;

          delete[] buf;
        }
# if 0
        printf("node %d is going to read %s and %s\n", node_id,
               prmtop_file.c_str(), rescrd_file.c_str());
# endif

        ReadAmberFile(prmtop_file,
                      rescrd_file,
                      &fcc->particle,
                      &covalent_bond, &covalent_bond_parameter_list,
                      &boxsize, &skew, restart, &simulation_start_time);
      }
      // clock_t inputdata_end_clock = clock();
      double inputdata_end_clock = getrealtime();
      if (node_id == 0) {
        // double idtime = static_cast<double>
        // (inputdata_end_clock-inputdata_start_clock) / CLOCKS_PER_SEC;
        double idtime = inputdata_end_clock-inputdata_start_clock;
        printf("input data %g sec\n",idtime);
      }
#else  // TUNE_NODEIO
# ifndef NUMIO
#  define NUMIO 16
# endif

      int num_ionode = NUMIO;
# ifdef __FCC_VERSION
      num_ionode = (num_node+31)/32;
# endif
# ifdef IO_Z0
      num_ionode = nodediv3d[0]*nodediv3d[1];
      int zsize = num_node/num_ionode;
      if((node_id==0)&&(DebugLog::verbose>0)){
	printf("ReadAmber at node Z=0\n");
	if(zsize!=nodediv3d[2]){
	  printf("z-size %d != nodediv3d[2] %d : with LongOnly node\n",zsize,nodediv3d[2]);
	}
      }
      int iog_size = zsize;
      int iog_id = node_id%num_ionode;
      int id_in_iog = node_id/num_ionode;
# else
      if(config_.md.io_node_number>0){
        num_ionode = config_.md.io_node_number;
      }
      int iog_size = (num_node-1)/num_ionode+1;
      int iog_id = node_id/iog_size;
      int id_in_iog = node_id%iog_size;
# endif
      MPI_Comm io_comm;
      MPI_Comm_split(MPI_COMM_WORLD,iog_id,id_in_iog,&io_comm);
      int num_iog_node;
      MPI_Comm_size(io_comm, &num_iog_node);
      int rank_iog;
      MPI_Comm_rank(io_comm,&rank_iog);
      int num_io_particle;
      Particle* io_particle;
      int *packed_bond;
      int pb_size;
      double *packed_bond_parm;
      int pbp_size;
      double bs[3];
      int bpnums[4];
      if((node_id==0)&&(DebugLog::verbose>0)){
        printf("Number of file load node %d\n",num_ionode);
        printf("Number of file recipient node %d %d\n",iog_size-1, num_iog_node-1);
      }
      if(rank_iog==0){
	//  clock_t inputdata_start_clock = clock();
	double inputdata_start_clock = getrealtime();
        ReadAmberFile(config_.amberfile.prmtop,
                      (restart
                       ? config_.amberfile.restrt
                       : config_.amberfile.inpcrd),
                      &fcc->particle,
                      &covalent_bond, &covalent_bond_parameter_list,
                      &boxsize, &skew, restart, &simulation_start_time);
	//  clock_t inputdata_end_clock = clock();
	double inputdata_end_clock = getrealtime();
	if (node_id == 0) {
	  //    double idtime = static_cast<double>(inputdata_end_clock-inputdata_start_clock) / CLOCKS_PER_SEC;
	  double idtime = inputdata_end_clock-inputdata_start_clock;
	  printf("input data %g sec\n",idtime);
	}
        if(num_iog_node>1){
	  double inputdatatransfer_start_clock = getrealtime();
          // bcast
          if((iog_id==0)||(DebugLog::verbose>1)){
            printf(" caster %d\n",node_id);
          }
          num_io_particle = fcc->particle.size();
          io_particle = new Particle[num_io_particle];
          pb_size = covalent_bond.size_of_packed_int();
          packed_bond = new int[pb_size];
          pbp_size = covalent_bond_parameter_list.size_of_packed_double(bpnums);
          packed_bond_parm = new double[pbp_size];
          bs[0] = boxsize.x;
          bs[1] = boxsize.y;
          bs[2] = boxsize.z;
          int sendnum[3];
          sendnum[0] = num_io_particle;
          sendnum[1] = pb_size;
          sendnum[2] = pbp_size;
          if((iog_id==0)||(DebugLog::verbose>1)){
            printf(" cast number %d %d %d\n",sendnum[0],sendnum[1],sendnum[2]);
          }
          MPI_Bcast(sendnum,3,MPI_INT,0,io_comm);
          if((iog_id==0)||(DebugLog::verbose>1)){
            printf(" send box %e %e %e\n",boxsize.x ,boxsize.y ,boxsize.z );
          }
          MPI_Bcast(bs,3,MPI_DOUBLE,0,io_comm);
          for(int i=0;i<num_io_particle;i++){
            io_particle[i] = getparticle(fcc->particle,i);
          }
          MPI_Bcast(io_particle,num_io_particle*sizeof(Particle),MPI_CHAR,0,io_comm);
          covalent_bond.pack_int_array(packed_bond);
          MPI_Bcast(packed_bond,pb_size,MPI_INT,0,io_comm);
          MPI_Bcast(bpnums,4,MPI_INT,0,io_comm);
          covalent_bond_parameter_list.pack_double_array(bpnums,packed_bond_parm);
          MPI_Bcast(packed_bond_parm,pbp_size,MPI_DOUBLE,0,io_comm);
          delete []packed_bond_parm;
          delete []packed_bond;
          delete []io_particle;
	  double inputdatatransfer_end_clock = getrealtime();
	  if (node_id == 0) {
	    //    double idtime = static_cast<double>(inputdata_end_clock-inputdata_start_clock) / CLOCKS_PER_SEC;
	    double idtime = inputdatatransfer_end_clock-inputdatatransfer_start_clock;
	    printf("transfer input data %g sec\n",idtime);
	  }
        }
      }else{
        if(num_iog_node>1){
	  double inputdatatransfer_start_clock = getrealtime();
          // recv
          if(DebugLog::verbose>1){
            printf(" receiver %d\n",node_id);
          }
          int sendnum[3];
          MPI_Bcast(sendnum,3,MPI_INT,0,io_comm);
          if(DebugLog::verbose>1){
            printf(" recv number %d %d %d\n",sendnum[0],sendnum[1],sendnum[2]);
          }
          num_io_particle = sendnum[0];
          pb_size = sendnum[1];
          pbp_size = sendnum[2];
          MPI_Bcast(bs,3,MPI_DOUBLE,0,io_comm);
          boxsize.x = bs[0];
          boxsize.y = bs[1];
          boxsize.z = bs[2];
          if(DebugLog::verbose>1){
            printf(" recv box %e %e %e\n",boxsize.x ,boxsize.y ,boxsize.z );
          }
          fcc->particle.clear();
          fcc->particle.resize(num_io_particle);
          io_particle = new Particle[num_io_particle];
          MPI_Bcast(io_particle,num_io_particle*sizeof(Particle),MPI_CHAR,0,io_comm);
          for(int i=0;i<num_io_particle;i++){
            setparticle(fcc->particle,io_particle[i],i);
          }
          packed_bond  = new int[pb_size];
          covalent_bond.clear();
          MPI_Bcast(packed_bond,pb_size,MPI_INT,0,io_comm);
          covalent_bond.unpack_int_array(packed_bond);
          packed_bond_parm = new double[pbp_size];
          covalent_bond_parameter_list.clear();
          MPI_Bcast(bpnums,4,MPI_INT,0,io_comm);
          MPI_Bcast(packed_bond_parm,pbp_size,MPI_DOUBLE,0,io_comm);
          covalent_bond_parameter_list.unpack_double_array(bpnums,packed_bond_parm);
          delete []packed_bond_parm;
          delete []packed_bond;
          delete []io_particle;
	  double inputdatatransfer_end_clock = getrealtime();
	  if (node_id == 0) {
	    //    double idtime = static_cast<double>(inputdata_end_clock-inputdata_start_clock) / CLOCKS_PER_SEC;
	    double idtime = inputdatatransfer_end_clock-inputdatatransfer_start_clock;
	    printf("transfer input data %g sec\n",idtime);
	  }
        }else{
          printf("Incorrect io group\n");
        }
      }
#endif  // TUNE_NODEIO
      if(DebugLog::verbose>1){
        printf(" size  %lu\n",fcc->particle.size());
      }

#ifdef CPPMD_ENABLE_PMMM
# ifdef PMMM_QPOS_OUT
      if(node_id==0){
	double len = std::max(std::max(boxsize.x,boxsize.y),boxsize.z);
	qpos_out(fcc->particle,fcc->particle.size(),len);
      }
# endif
#endif

      total_num_particle = fcc->particle.size();
      {
	for(int i=0;i<total_num_particle;i++){
	  total_charge += getcharge(fcc->particle,i);
	}
      }
#ifdef USE_SHAKE
      heavyid.clear();
      heavyid.resize(fcc->particle.size(),-1);
      mark_shake(fcc->particle, covalent_bond, heavyid);
#endif
      if (copy_mode) {  // for COPY_MODE
        // TODO: [AH] I'm not sure where to insert the following lines.
        // These lines must be placed after the ReadAmberFile() call.
	double inputdatacopy_start_clock = getrealtime();
        num_copy = ShiftParticleAndCovalentBond(short_id,
                                                fcc->particle.size(),
                                                nodediv3d, copynum3d, &boxsize,
                                                &fcc->particle, &covalent_bond, heavyid,
                                                num_copy_local);
	double inputdatacopy_end_clock = getrealtime();
#define DENSE_NP_10CUBE 125
	if (node_id == 0) {
	  //    double idtime = static_cast<double>(inputdata_end_clock-inputdata_start_clock) / CLOCKS_PER_SEC;
	  double idtime = inputdatacopy_end_clock-inputdatacopy_start_clock;
	  printf("copy input data %g sec\n",idtime);
	}
      } else {
        num_copy = 1;
        if(config_.md.bsize>boxsize.x*2.0){
          boxsize = SpaceVector<double>(config_.md.bsize);
        }
      }
      total_num_particle *= num_copy;
      total_charge *= num_copy;
      num_particle = total_num_particle / num_node;

      if (node_id == 0)
        printf("boxsize (%g,%g,%g)\n", boxsize.x, boxsize.y, boxsize.z);
      /*
	max_per_cell = static_cast<int>(
	0.5 * (boxsize.x * boxsize.y * boxsize.z) / num_set);
      */
      {
        int np_cell = static_cast<int>(
				       ceil(double(fcc->particle.size()) / num_set));
        int dense_np = int(boxsize.x * boxsize.y * boxsize.z*DENSE_NP_10CUBE/1000.0/num_set);
        int max_np = std::max(np_cell, dense_np);
        max_per_cell = max_np + int(ceil(sqrt(double(max_np))))*3;
	max_per_cell = (max_per_cell/4+1)*4;
      }
      break;
    }
    default: {  // NOT Config::kAMBER (for now)
      double num_density = 8;
      if (config_.md.num_lattice > 0) {
	fcc->latticeNum(config_.md.num_lattice);
	// TODO: use fcc->particle.size() or fcc->natom
	total_num_particle = fcc->setLattice();
	num_particle = total_num_particle / num_node;
	num_density = total_num_particle/(fcc->side.x*fcc->side.y*fcc->side.z);
	if ((config_.md.bsize < fcc->side.x)) {
	  boxsize = fcc->side;
	}else{
	  if (node_id == 0) {
	    printf("natom%d\n", fcc->natom);
	    printf("boxsize(%g,%g,%g) larger than lattice size(%g,%g,%g)\n",
		   boxsize.x, boxsize.y, boxsize.z,
		   fcc->side.x, fcc->side.y, fcc->side.z);
	  }
	}
      }
      int np_cell = static_cast<int>(
				     ceil(num_density * boxsize.x * boxsize.y * boxsize.z / num_set));
      /*
	max_per_cell = std::max(num_particle * num_node * 2 / num_set, 
	std::max(8, np_cell*4));
      */
      {
	int dense_np = int(boxsize.x * boxsize.y * boxsize.z*DENSE_NP_10CUBE/1000.0/num_set);
	int max_np = std::max(np_cell, dense_np);
	max_per_cell = max_np + int(ceil(sqrt(double(max_np))));
	max_per_cell = (max_per_cell/4+1)*4;
      }
      break;
    }
  }

  //  clock_t makecell_start_clock = clock();
  double makecell_start_clock = getrealtime();

  //printf(" max_per_cell %d\n", max_per_cell);
  if (node_id == 0) {
    dump_cellindextype(citype);
    printf(" number of local copy %d\n",num_copy_local);
    printf(" number of total particle %lu\n",
           static_cast<unsigned long>(total_num_particle));
    printf(" number of particle per node : average %d reserved %d\n",
           num_particle, max_per_cell);
    if (copy_mode) {  // for COPY_MODE
      printf(" number of total covalent bond "
             " Bond %lu  Angle %lu  Torsion %lu  Improper %lu\n",
             static_cast<unsigned long>(covalent_bond.bond.size())*num_copy/num_copy_local,
             static_cast<unsigned long>(covalent_bond.angle.size())*num_copy/num_copy_local,
             static_cast<unsigned long>(covalent_bond.torsion.size())*num_copy/num_copy_local,
             static_cast<unsigned long>(covalent_bond.improper.size())*num_copy/num_copy_local
            );
    } else {
      printf(" number of total covalent bond "
             " Bond %lu  Angle %lu  Torsion %lu  Improper %lu\n",
             static_cast<unsigned long>(covalent_bond.bond.size()),
             static_cast<unsigned long>(covalent_bond.angle.size()),
             static_cast<unsigned long>(covalent_bond.torsion.size()),
             static_cast<unsigned long>(covalent_bond.improper.size())
            );
    }
    printf(" total charge %e\n",total_charge);
    {  /// count periodicity of torsion and improper
      int bct[5] = {0,0,0,0,0};
      int bci[5] = {0,0,0,0,0};
      { // count 1 copy
	for (int i=0;i<covalent_bond.torsion.size();i++) {
	  int periodicity = covalent_bond_parameter_list.torsion[covalent_bond.torsion[i].typeoftorsion].periodicity;
	  if((periodicity>4)||(periodicity<1)){
	    bct[0]++;
	  }else{
	    bct[periodicity]++;
	  }
	}
	for (int i=0;i<covalent_bond.improper.size();i++) {
	  int periodicity = covalent_bond_parameter_list.improper[covalent_bond.improper[i].typeofimproper].periodicity;
	  if((periodicity>4)||(periodicity<1)){
	    bci[0]++;
	  }else{
	    bci[periodicity]++;
	  }
	}
      }
      int bc[5];
      if(copy_mode) {
	for(int i=0;i<5;i++)bc[i]=(bct[i]+bci[i])*num_copy/num_copy_local;
      }else{
	for(int i=0;i<5;i++)bc[i]=bct[i]+bci[i];
      }
      printf(" local copy : periodicity of Torsion 1 2 3 4 other : %d %d %d %d %d\n",bct[1],bct[2],bct[3],bct[4],bct[0]);
      printf(" local copy : periodicity of Improper 1 2 3 4 other : %d %d %d %d %d\n",bci[1],bci[2],bci[3],bci[4],bci[0]);
      printf(" total : periodicity of Torsion and Improper 1 2 3 4 other : %d %d %d %d %d\n",bc[1],bc[2],bc[3],bc[4],bc[0]);
    }

    printf(" unit lattice (%g,%g,%g)\n",
           fcc->latticeSpacing.x, fcc->latticeSpacing.y, fcc->latticeSpacing.z);
    printf(" num lattice (%d,%d,%d)\n",
           fcc->latticeNum.x, fcc->latticeNum.y, fcc->latticeNum.z);
    printf(" box size (%g,%g,%g)\n", boxsize.x, boxsize.y, boxsize.z);
    printf(" cutoff %g\n", config_.md.cutoff);
  }

  // shift the whole system for debug
  if (fabs(config_.md.offset) >= std::numeric_limits<double>::epsilon()) {
    Position sft(fcc->side * config_.md.offset);
    if (node_id == 0) printf(" shift (%g,%g,%g)\n", sft.x, sft.y, sft.z);
    for (size_t i = 0; i< fcc->particle.size(); ++i) {
      getpos(fcc->particle,i) += sft;
      periodic_shift(getpos(fcc->particle,i).x, boxsize.x);
      periodic_shift(getpos(fcc->particle,i).y, boxsize.y);
      periodic_shift(getpos(fcc->particle,i).z, boxsize.z);
    }
  }

  // set initial velocity for debug
  if (fabs(config_.md.velocity) >= std::numeric_limits<double>::epsilon()) {
    if (node_id == 0) printf(" initial velocity %g\n", config_.md.velocity);
    Velocity vs(config_.md.velocity);
    for(size_t i = 0; i < fcc->particle.size(); ++i) {
      getvelocity(fcc->particle,i) += vs;
    }
  }

  switch (config_.md.input_mode) {
    case Config::kARGON: {
      ShortRange::ljmixparameters.resize2d(1,1);
      ShortRange::ljmixparameters[0][0] = ljmixAr[0][0];
      break;
    }
    default: {
      LJAmber94 ljamber94;
      ljamber94.convertLJMixparameterArray(ShortRange::ljmixparameters);
      break;
    }
  }

  // LJ cutoff energy correction
  double ljcec = 0.0;
  if (config_.md.cutoff > 0.0) {
    switch (config_.md.input_mode) {
    case Config::kARGON: {    
      LJAr LJArgon;
      ljcec = LJArgon.calcLJCutoffEnergyCorrection(config_.md.cutoff);
      break;
    }
    default:{
      ljcec = calcLJCutoffEnergyCorrection(config_.md.cutoff,ATOMTYPE_WO,ATOMTYPE_WO);
      break;
    }
    //if (node_id == 0) printf(" ljcec %f\n", ljcec);
    }
  }

#ifdef USE_PAIRLIST
  double pairlistmargin = 2.0;
#endif

  /*
  //  std::cout << "Long Plan " << long_plan << std::endl;
  SpaceVector<double> cellmargin(1.0, 1.0, 1.0);

  CalcPreparator_::makeoperations(node_id, 
  config_.md.withlong, 
  config_.md.withbond, 
  citype,
  long_plan, num_long_only);
  int short_id;
  std::vector<int> node_id_of_shorts;
  CalcPreparator_::maketarget(node_id, short_id, node_id_of_shorts,
  num_node, num_particle, num_total_set,
  config_.md.celldiv, max_per_cell, boxsize,
  cellmargin, citype, config_.md.cutoff,
  long_plan, num_long_only);
  */
  double boxsize_minimun_scale=1.0;
  if(config_.tempctrl.method==Config::kANDERSEN_HOOVER){
    boxsize_minimun_scale = 0.98;
    if (node_id == 0) {
      printf("boxsize_minimun_scale %f\n",boxsize_minimun_scale);
    }
  }
  SpaceVector<double> cellmargin(1.0, 1.0, 1.0);
  CalcPreparator_::maketarget(node_id, short_id, num_short,
                              total_num_particle, max_per_cell, 
                              boxsize, cellmargin, citype, config_.md.cutoff,
			      boxsize_minimun_scale,
                              num_total_set);

  //  clock_t makecell_end_clock = clock();
  double makecell_end_clock = getrealtime();
  if (node_id == 0) {

    //    double mctime = static_cast<double>(makecell_end_clock-makecell_start_clock) / CLOCKS_PER_SEC;
    double mctime = makecell_end_clock-makecell_start_clock;
    printf("make cell %g sec\n",mctime);
  }
  //  clock_t makepa_start_clock = clock();
  double makepa_start_clock = getrealtime();

  WaterList waterlist(fcc->particle);
  int waterlist_size = waterlist.size();
  int num_freedom = (3 * fcc->particle.size()) - (3 * waterlist.size());

/////////// SHAKE
  double makeshake_start_clock = getrealtime();
//  ShakeList shakelist(fcc->particle);  // make global shake bond
  ShakeList shakelist;
#ifdef USE_SHAKE
  if(copy_mode){ /// atomid of copy particle is different from index
    //    append_shake_scanid(shakelist,fcc->particle,covalent_bond);
    //    append_shake_scanid_with_mark(shakelist,fcc->particle,covalent_bond);
    //    append_shake_offset(shakelist,fcc->particle,covalent_bond);
    if(node_id==0){
      printf("find shake in %ld bond from %ld particles\n",covalent_bond.bond.size(),fcc->particle.size());
    }
    append_shake_offset_with_mark(shakelist,fcc->particle,covalent_bond);
  }else{
    append_shake(shakelist,fcc->particle,covalent_bond);
  }
#endif
  int nw=0;
  for(int i=0;i<fcc->particle.size();i++){
    if(getatomtype(fcc->particle,i)<ATOMTYPE_WH)nw++;
  }
/////////// SHAKE
  if(node_id == 0){
    std::cout << "Number of Hydrogen (local, not water) " << nw << std::endl;
    std::cout << "Number of Shake Hydrogen(local) " << shakelist.reverse_list.size() << std::endl;
    std::cout << "Number of Shake HeavyAtom(local) " << shakelist.size() << std::endl;
    if(DebugLog::verbose>1){
      std::cout << "in simulator global shakelist=" << std::endl;
      dump_shakelist_distance(shakelist,fcc->particle);
    }
  }
  double makeshake_end_clock = getrealtime();
  if (node_id == 0) {

    double mstime = makeshake_end_clock-makeshake_start_clock;
    printf("make shakelist %g sec\n",mstime);
  }
  
  int shakelist_size = shakelist.reverse_list.size();
  if(config_.shake.type>0){
    num_freedom -= shakelist_size;
  }

  if (copy_mode) {  // for COPY_MODE
    if(short_id!=NO_SHORT_ID){
      waterlist_size = waterlist_size*num_copy/num_copy_local;
      shakelist_size = shakelist_size*num_copy/num_copy_local;
      num_freedom = num_freedom*num_copy/num_copy_local;
    }else{
      waterlist_size = 0;
      shakelist_size = 0;
      num_freedom = 0;
    }
  }

  if(node_id == 0){
    std::cout << "Number of Water Molecules " << waterlist_size << std::endl;
    if(config_.shake.type>0){
      std::cout << "Number of Shake Hydrogen " << shakelist_size << std::endl;
    }
    std::cout << "Number of Freedom " << num_freedom << std::endl;
  }
#if 1
  fcc->potentialmodel.clear();
  for (size_t i = 0; i < fcc->particle.size(); ++i) {
    /*
       switch (fcc->particle[i].atomtype) {
       case 1 :fcc->particle[i].atomtype = ATOMTYPE_WO; break;
       case 2 :fcc->particle[i].atomtype = ATOMTYPE_WH; break;
       }
       */
    switch (getatomtype(fcc->particle,i)) {
      case ATOMTYPE_WO:
        fcc->potentialmodel.push_back(LJCoulombPotential);
        break;
      case ATOMTYPE_WH:
        fcc->potentialmodel.push_back(OnlyCoulombPotential);
        break;
      default:
        fcc->potentialmodel.push_back(LJCoulombPotential);
        break;
    }
  }
#endif

/////////// SHAKE
  CalcPreparator_::makeparticle(short_id, total_num_particle, num_freedom,
#if 1
                                num_copy,
                                num_copy_local,
#endif
                                fcc->particle, waterlist, shakelist,
                                fcc->potentialmodel, ljcec);
  //  if(node_id == 0)
  //    std::cout << "Number of Water Molecules(after makeparticle)" << waterlist.size() << std::endl;
  bool excludewaterbond=true;
/////////// SHAKE
  // make local shake bond too
  CalcPreparator_::makebond(short_id, covalent_bond, excludewaterbond,
                            shakelist, config_.shake.type);
  std::map<int,int> shortidtorank;
  for (int u = 0; u < num_short; ++u) shortidtorank.insert(std::pair<int,int>(u,u));

  //  clock_t makepa_end_clock = clock();
  double makepa_end_clock = getrealtime();
  if (node_id == 0) {
    //    double mptime = static_cast<double>(makepa_end_clock-makepa_start_clock) / CLOCKS_PER_SEC;
    double mptime = makepa_end_clock-makepa_start_clock;
    printf("make subset particle and bond %g sec\n",mptime);
  }

  ShortRange::CoulombType cltype;
#ifndef FORCE_MWOLF
  if(config_.md.coulomb_type==Config::kZeroDipole){
    cltype = ShortRange::ZeroDipole;
    if(node_id == 0){
      std::cout << "Coulomb Type for ZeroDipole" << std::endl;
    }
  }else if(config_.md.coulomb_type==Config::kEwald){
    cltype = ShortRange::ForEwaldReal;
    if(node_id == 0){
      std::cout << "Coulomb Type for EwaldReal" << std::endl;
    }
  }else{
    cltype = ShortRange::OriginalCoulomb;
    if(node_id == 0){
      std::cout << "Coulomb Type OriginalCoulomb" << std::endl;
    }
  }
#else
  cltype = ShortRange::MWOLF;
  if(node_id == 0){
    std::cout << "Coulomb Type Modified Wolf" << std::endl;
  }
#endif

  //  clock_t makelong_start_clock = clock();
  double makelong_start_clock = getrealtime();

  double pme_alpha = config_.pme.alpha;
  double pme_kCutoff = 0.0;
  double pme_order = config_.pme.order;
  int pme_surfaceDipole = 0;
  SpaceVector<double> gridLengths(config_.pme.grid_length);
  PMEType pmetype = SmoothPME;
  int grid_1d(config_.pme.grid_number);
  SpaceVector<int> grid_num(64,64,64);
  GeometryXYZ long_geometry;
  MGType multigridType = FMG;
  int multigridIteration = 1;
  int multigridFFTlevel = 0;
  int num_long;
  if(config_.pme.grid_length==0.0){
    if(grid_1d>0){
      grid_num.x = grid_num.y = grid_num.z = grid_1d;
    }
    gridLengths.x = boxsize.x/grid_num.x;
    gridLengths.y = boxsize.y/grid_num.y;
    gridLengths.z = boxsize.z/grid_num.z;
  }else{
    grid_num.x = int(ceil(boxsize.x/gridLengths.x));
    grid_num.y = int(ceil(boxsize.y/gridLengths.y));
    grid_num.z = int(ceil(boxsize.z/gridLengths.z));
    if(grid_1d>0){
      if((grid_num.x<grid_1d)||(grid_num.y<grid_1d)||(grid_num.z<grid_1d)){
	double gl = (std::min(std::min(boxsize.x,boxsize.y),boxsize.z))/grid_1d;
	gridLengths.x = gridLengths.y = gridLengths.z = gl;
	grid_num.x = int(ceil(boxsize.x/gridLengths.x));
	grid_num.y = int(ceil(boxsize.y/gridLengths.y));
	grid_num.z = int(ceil(boxsize.z/gridLengths.z));
      }
    }
  }
  CalcPreparator_::makelongrangegeometry(long_plan,num_long_only,num_long,long_geometry);
  int long_id;
  CalcPreparator_::makelongcomm(long_plan,long_geometry,num_long,num_short,
                                MPI_COMM_WORLD,node_id,short_id,long_id);
  std::vector<int> long_reqcell_list(0);
  double long_fringe = 0.0;
#ifndef CPPMD_ENABLE_EWALD
  long_fringe = boxsize.x/grid_num.x*(pme_order/2);
#endif  // CPPMD_ENABLE_EWALD
  CalcPreparator_::makelongrangerequiredcell(long_geometry,long_fringe,
                                             long_id,
                                             long_reqcell_list);
  std::vector<int> selfenergycell_list(0);
  CalcPreparator_::makeselfenergycell_list(long_geometry,
                                           long_id,
                                           long_reqcell_list,
                                           selfenergycell_list);


  if(DebugLog::verbose>1){
    if(long_id!=NO_LONG_ID){
      std::cout << "node_id " << node_id << " long_id " << long_id << " required " << long_reqcell_list.size() << " cells for long";
      for(std::vector<int>::size_type i=0;i<long_reqcell_list.size();i++){
        std::cout << " " << long_reqcell_list[i];
      }
      std::cout << std::endl;
    }
  }
  std::vector< std::vector<int> > send_to_long_list;
  CalcPreparator_::makecellidlistsendtolong(long_geometry,long_fringe,num_long,
                                            short_id,
                                            send_to_long_list);

  if(DebugLog::verbose>1){
    if(short_id!=NO_SHORT_ID){
      std::cout << "node_id " << node_id << " short_id " << short_id << " send to long";
      for(std::vector< std::vector<int> >::size_type i=0;i<send_to_long_list.size();i++){
        std::cout << " " << i << ":" << send_to_long_list[i].size();
      }
      std::cout << std::endl;
    }
  }

  std::vector<MPI_Comm> mpi_comm_ls_list;
  std::vector<int> long_short_id;
  std::vector<int> sender_local_id;
  std::map<int,int> idinlong_to_longrank;
  std::vector<int> receiver_local_rank;
  std::vector< std::vector<int> > long_recv_set_id_lists;
  std::vector<int> sender_global_id_list;
  std::vector<int> reciever_global_id_list;

  CalcPreparator_::make_long_short_comm(num_node, num_long, 
                                        node_id, short_id, long_id,
                                        MPI_COMM_WORLD,
                                        send_to_long_list,
                                        long_reqcell_list,
                                        mpi_comm_ls_list,
                                        long_short_id,
                                        sender_local_id,
                                        idinlong_to_longrank,
                                        receiver_local_rank,
                                        long_recv_set_id_lists,
                                        sender_global_id_list,
                                        reciever_global_id_list);

  CalcPreparator_::makelongrangeparameter(config_.md.cutoff,
                                          pme_alpha,
                                          pme_kCutoff,
                                          pme_surfaceDipole,
                                          boxsize,
                                          gridLengths,
                                          int(pme_order),
                                          pmetype,
                                          grid_num,
                                          long_geometry,
                                          multigridType,
                                          multigridIteration,
                                          multigridFFTlevel
                                         );
  MPI_Barrier(MPI_COMM_WORLD);
  //  printf("MPI_Barrier after makelongrangeparameter %d\n", node_id);
  //  fflush(stdout);
  //  clock_t makelong_end_clock = clock();
  double makelong_end_clock = getrealtime();
  if (node_id == 0) {
    //    double mltime = static_cast<double>(makelong_end_clock-makelong_start_clock) / CLOCKS_PER_SEC;
    double mltime = makelong_end_clock-makelong_start_clock;
    printf("make long %g sec\n",mltime);
  }
  //  clock_t makecomm_start_clock = clock();
  double makecomm_start_clock = getrealtime();

  if(config_.md.coulomb_type==Config::kZeroDipole){
#ifdef ZERODIPOLE0
    ShortRange::alpha = 0.0;
    if (node_id == 0) {
      printf("ZeroDipole with zero alpha\n");
    }
#endif
    if (node_id == 0) {
      printf("ZeroDipole alpha = %f\n",ShortRange::alpha);
    }
  }

  CalcPreparator_::constructcommunicator(short_id, max_per_cell, shortidtorank,
                                         num_long, long_id,
                                         mpi_comm_ls_list, long_short_id, 
                                         sender_local_id, 
                                         long_recv_set_id_lists,
                                         receiver_local_rank,
                                         send_to_long_list,
                                         idinlong_to_longrank,
                                         config_.md.short_comm_pattern,
                                         config_.md.move_comm_pattern);
  MPI_Barrier(MPI_COMM_WORLD);
  //  printf("MPI_Barrier after constructcommunicator %d\n", node_id);
  //  fflush(stdout);
  //  clock_t makecomm_end_clock = clock();
  double makecomm_end_clock = getrealtime();
  if (node_id == 0) {
    //    double mcmtime = static_cast<double>(makecomm_end_clock-makecomm_start_clock) / CLOCKS_PER_SEC;
    double mcmtime = makecomm_end_clock-makecomm_start_clock;
    printf("make comm %g sec\n",mcmtime);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  //  printf("MPI_Barrier before constructcalculationunit %d\n", node_id);
  //  fflush(stdout);
  //clock_t makeunit_start_clock = clock();
  double makeunit_start_clock = getrealtime();

  CalculationUnit<PACLASS,GPACLASS>*calculationunit 
      =       CalcPreparator_::constructcalculationunit<PACLASS,GPACLASS>
    (node_id, short_id,
                                                boxsize, citype,
                                                config_.md.cutoff,
#ifdef USE_PAIRLIST
                                                pairlistmargin,
#endif
                                                config_.tempctrl,
                                                cltype, 
                                                covalent_bond_parameter_list);
  calculationunit->initialize(celldiv3d,nodediv3d);
  //  if(calculationunit->longsender.sprf.size()>0){
  //    std::cout << " calculationunit->longsender.sprf[0].send_requestp  " <<calculationunit->longsender.sprf[0].send_requestp  << std::endl;
  //  }
#ifdef CPPMD_ENABLE_FMM
  CalcPreparator_::construct_fmm_target_cell(calculationunit->calcforce);
#endif  // CPPMD_ENABLE_FMM
#ifdef CPPMD_ENABLE_PMMM
  CalcPreparator_::construct_pmmm_target_cell(calculationunit->calcforce);
  CalcPreparator_::construct_pmmm(calculationunit->calcforce, num_short, num_long, total_charge);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  //  printf("MPI_Barrier after constructcalculationunit %d\n", node_id);
  //  fflush(stdout);
  // TEST, must implement correct sets
  CalcPreparator_::constructselflongset(long_id,short_id,
                                        calculationunit->calcforce.shortset,
                                        send_to_long_list,
                                        selfenergycell_list,
                                        calculationunit->calcforce.self_longset_index,
                                        calculationunit->calcforce.self_selfenergycell_index);
  CalcPreparator_::constructghostlongset(long_id,short_id,
                                         long_reqcell_list,
                                         long_recv_set_id_lists,
                                         calculationunit->calcforce.shortset,
                                         selfenergycell_list,
                                         calculationunit->calcforce.ghostlongset,
                                         calculationunit->calcforce.ghost_selfenergycell_list);


  MPI_Barrier(MPI_COMM_WORLD);
  //  clock_t makeunit_end_clock = clock();
  double makeunit_end_clock = getrealtime();
  if (node_id == 0) {
    //    double mutime = static_cast<double>(makeunit_end_clock-makeunit_start_clock) / CLOCKS_PER_SEC;
    double mutime = makeunit_end_clock-makeunit_start_clock;
    printf("make u unit %g sec\n",mutime);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  Dump *dump = NULL;
  if (config_.amberfile.mdcrd_interval > 0) {
    dump = new Dump(boxsize, total_num_particle, 
                    num_node, node_id,
                    config_.md.dump_without_water,
                    config_.amberfile.mdcrd
#ifdef DEBUG_UNITOUTPUT
                    , num_copy
#endif  // DEBUG_UNITOUTPUT
                    );
    dump->GatherDumpAmberCRD(calculationunit->particlearray,
                             calculationunit->typerangearray);
  }else{
    dump = new Dump(boxsize, total_num_particle,
                    num_node, node_id,
                    config_.md.dump_without_water
#ifdef DEBUG_UNITOUTPUT
                    , num_copy
#endif  // DEBUG_UNITOUTPUT
                    );
  }
#ifdef USE_HDF
  {
    std::string hdffilename("hdfrestore.hdf5");
    int opened=0;
    if(calculationunit->operations.doShortrangecalculation){
      opened = calculationunit->integrator.open_hdfrestore(hdffilename);
    }
    int noopen = 0;
    MPI_Allreduce(&opened,&noopen,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(noopen==0){
      calculationunit->integrator.restore = true;
      if(node_id==0){
	printf("restore from hdfrestore.XXXXX.hdf5\n");
      }
    }else{
      if(node_id==0){
	printf("%d nodes cannt open hdfrestorefile\n",noopen);
      }
    }
  }
#else
#ifdef BINARY_RESTORE
  {
    char name[256];
    snprintf(name,256,"binaryrestorefile.%05d",node_id);
    std::string binfilename(name);
    int opened=0;
    if(calculationunit->operations.doShortrangecalculation){
      opened = calculationunit->integrator.open_binaryrestore(binfilename);
    }
    int noopen =0;
    MPI_Allreduce(&opened,&noopen,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(noopen==0){
      calculationunit->integrator.restore = true;
      if(node_id==0){
	printf("restore from %s\n",name);
      }
    }else{
      if(node_id==0){
	printf("%d nodes cannt open binaryrestorefile\n",noopen);
      }
    }
  }
#endif
#endif
#ifdef USE_HDF
  if(calculationunit->operations.doShortrangecalculation){
    std::string hdffilename("hdfdump.hdf5");
    calculationunit->integrator.open_hdfdump(hdffilename);
    if(node_id==0){
      printf("dump to hdfdump.XXXXX.hdf5\n");
    }
  }
#else
#ifdef BINARY_DUMP
  if(calculationunit->operations.doShortrangecalculation){
    char name[256];
    snprintf(name,256,"binarydumpfile.%05d",node_id);
    std::string binfilename(name);
    calculationunit->integrator.open_binarydump(binfilename);
    if(node_id==0){
      printf("dump to %s\n",name);
    }
  }
#endif
#endif

#ifndef LATEFREEFCC
  delete fcc, fcc = NULL;                       // is this ok here?
#endif

  Dump *dump_rst = NULL;
  int write_restart = 0;
  if (config_.amberfile.mdrst_interval > 0) {
    write_restart = 1;
    std::string mdrst_output_file = config_.amberfile.mdrst;
    mdrst_output_file.append(".save");
    dump_rst = new Dump(boxsize, total_num_particle, 
                        num_node, node_id,
                        config_.md.dump_without_water,
                        mdrst_output_file
#ifdef DEBUG_UNITOUTPUT
                        , num_copy
#endif  // DEBUG_UNITOUTPUT
                        , write_restart);
  }
  calculationunit->calcforce.covalentbond.setParameterList(covalent_bond_parameter_list);
  //====================================================================
  PerfCounter::start(perf_target);
  //  clock_t start_clock = clock();
  double start_clock = getrealtime();

  //  const double dt = 0.053763;  // 0.215107
  const double dt = (config_.md.delta_t * 1e-15) / UnitParameter::unitTime;
  if (node_id==0) {
    printf("delta_t = %g fs (%g)\n",config_.md.delta_t,dt);
  }

  if(initialforce==true){
    calculationunit->integrator.pre_calcforce = false;
  }else{
    calculationunit->integrator.pre_calcforce = true;
  }

  //! TODO move_cell_interval set by option
#ifndef MOVE_CELL_INTERVAL
#define MOVE_CELL_INTERVAL 40
#endif
  int move_cell_interval;

  if(config_.md.ci_update_interval>0){
    move_cell_interval=config_.md.ci_update_interval;
  }else{
    move_cell_interval = MOVE_CELL_INTERVAL;
  }
  if (node_id==0) {
    printf("interval inter-cell move %d step\n",move_cell_interval);
  }

  bool complete_integration;
  long last_t;
/////////// SHAKE
  complete_integration = calculationunit->startTimeIntegration(dt, config_.md.tmax,
                                       shakelist,
                                       config_.shake.type,
                                       config_.shake.max_iterate,
                                       config_.shake.tolerance,
                                       config_.md.reduce_interval,
                                       config_.md.print_interval,
                                       config_.amberfile.mdcrd_interval,
                                       config_.amberfile.mdrst_interval,
                                       *dump, *dump_rst,
                                        move_cell_interval,
							       write_restart, last_t);
  //  clock_t end_clock = clock();
  double end_clock = getrealtime();
  PerfCounter::stop();
  //====================================================================

  // print running time information
  //  double calc_time = static_cast<double>(end_clock - start_clock) / CLOCKS_PER_SEC;
  double calc_time = end_clock - start_clock;
  {
    std::vector<double> cts(num_node, 0.0);
    MPI_Gather(&calc_time,1,MPI_DOUBLE, &cts[0],1,MPI_DOUBLE, 0,MPI_COMM_WORLD);
    if (node_id == 0) {
      printf("time max %g :", *std::max_element(cts.begin(), cts.end()));
      if (DebugLog::verbose > 1)
        for (int n = 0; n < num_node; ++n) printf(" %g", cts[n]);
      putchar('\n');
    }
  }
  if (node_id == 0) { 
    puts("timer of rank 0");
    PerfCounter::print(); 
    puts("end"); 
  }

  { 
    //    PerfCounter::print(node_id); 

    int num_timer;
    double *times = PerfCounter::get_time_array(num_timer);
    double *max_times = new double[num_timer];
    MPI_Reduce(times,max_times,num_timer,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if (node_id == 0) { 
    PerfCounter::put_time_array(max_times,num_timer);
      std::cout << "max time" << std::endl;
      PerfCounter::print_time();
    }
    delete []max_times;
    delete []times;
  }
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef CPPMD_ENABLE_LONGRANGE
  if(node_id == 0){
# if defined (CPPMD_ENABLE_PME) || defined (CPPMD_ENABLE_OLDPME)
    calculationunit->calcforce.pmelongrange.lr_time_dump();
# endif  // CPPMD_ENABLE_PME || CPPMD_ENABLE_OLDPME
  }
#endif  // CPPMD_ENABLE_LONGRANGE

  if (dump_rst != NULL) delete dump_rst, dump_rst = NULL;
  if (dump != NULL) delete dump, dump = NULL;

  if(config_.tempctrl.method==Config::kANDERSEN_HOOVER){
    SpaceVector<double> oldbs = boxsize;
    boxsize = calculationunit->integrator.boxsize;
    if (node_id==0) {
      std::cout << "Boxsize change " << oldbs << " -> " << boxsize << std::endl;
    }
  }
  if(write_restart!=0){
    //    double simulation_end_time = simulation_start_time + config_.md.delta_t*(double)(config_.md.tmax)*1e-3;
    double simulation_end_time = simulation_start_time + config_.md.delta_t*(double)(last_t)*1e-3;
    MPI_Barrier(MPI_COMM_WORLD);
    Dump *restart_dump;
    std::string mdrst_output_file = config_.amberfile.mdrst;
    if(complete_integration==false){
      mdrst_output_file.append("_rollback_");
      std::ostringstream lt;
      lt << last_t;
      mdrst_output_file += lt.str();
    }
    if (node_id==0) {
      printf("dump restart %s\n",mdrst_output_file.c_str());
    }
    restart_dump = new Dump (boxsize, total_num_particle, 
                             num_node, node_id,
                             config_.md.dump_without_water,
                             mdrst_output_file,
#ifdef DEBUG_UNITOUTPUT
                             num_copy,
#endif  // DEBUG_UNITOUTPUT
                             write_restart);
    restart_dump->amber_crd.current_time = simulation_end_time;
    restart_dump->GatherDumpAmberCRD(calculationunit->particlearray,
                             calculationunit->integrator.typerangearray);
    MPI_Barrier(MPI_COMM_WORLD);

  }
  


#ifdef LATEFREEFCC
  delete fcc, fcc = NULL;                       // is this ok here?
#endif
  return 0;
}

}  // namespace cppmd
