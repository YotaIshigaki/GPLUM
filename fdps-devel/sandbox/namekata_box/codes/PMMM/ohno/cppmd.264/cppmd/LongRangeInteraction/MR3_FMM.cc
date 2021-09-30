#include "MR3_FMM.h"

#ifdef SAMPLE_FMM
# ifdef FPCOLL
#  include <fj_tool/fjsamp.h>
#  define SAMPLE_START() fpcoll_start()
#  define SAMPLE_STOP() fpcoll_stop()
# else
#  include <fj_tool/fipp.h>
#  define SAMPLE_START() fipp_start()
#  define SAMPLE_STOP() fipp_stop()
# endif
#endif
#ifndef SAMPLE_START
# define SAMPLE_START() 
# define SAMPLE_STOP() 
#endif
#ifdef FAPP
# include <fj_tool/fapp.h>
#endif
#ifdef K_PA
#include "fjcoll.h"
#endif
#include "../EXAFMM2/parallelfmm.h"


inline
void counter_start(const char *timername, const int timer_id=0)
{
#ifdef FAPP
  fapp_start(timername,timer_id,1);
#endif
#ifdef K_PA
  start_collection(timername);
#endif
}

inline
void counter_stop(const char *timername, const int timer_id=0)
{
#ifdef K_PA
  stop_collection(timername);
#endif
#ifdef FAPP
  fapp_stop(timername,timer_id,1);
#endif
}


#ifndef FMM_IMAGES
#define FMM_IMAGES 2
#endif

MR3_FMM::~MR3_FMM()
{
  parallelfmm->deallocate();
  delete parallelfmm;
}

void MR3_FMM::initialize(size_t size, SpaceVector<int> celldiv3d, SpaceVector<int> nodediv3d)
{
  parallelfmm = new ParallelFMM();
  numcell3d[0] = celldiv3d[0];
  numcell3d[1] = celldiv3d[1];
  numcell3d[2] = celldiv3d[2];
  num_cell = numcell3d[0]*numcell3d[1]*numcell3d[2];
  numlocalcell3d[0] = celldiv3d[0]/nodediv3d[0];
  numlocalcell3d[1] = celldiv3d[1]/nodediv3d[1];
  numlocalcell3d[2] = celldiv3d[2]/nodediv3d[2];
  num_local_cell = numlocalcell3d[0]*numlocalcell3d[1]*numlocalcell3d[2];
  int level = int(log(num_local_cell+1)/M_LN2/3);
  int images = FMM_IMAGES;
#if 1
  parallelfmm->allocate((int)size,level,images);
  allocated = true;
  parallelfmm->printNow = false;
#endif
  numnode3d[0] = nodediv3d[0];
  numnode3d[1] = nodediv3d[1];
  numnode3d[2] = nodediv3d[2];
  parallelfmm->partitioner(numnode3d,1); // Second input is gather level
  if(parallelfmm->MPIRANK==0){
    printf("FMM Natom %ld Node(%d,%d,%d), Cell(%d,%d,%d), LocalCell(%d,%d,%d), LocalLevel %d, Cartesian order %d\n",
	   size,
	   numnode3d[0],numnode3d[1],numnode3d[2],
	   numcell3d[0],numcell3d[1],numcell3d[2],
	   numlocalcell3d[0],numlocalcell3d[1],numlocalcell3d[2],
	   level, P);
  }
}

template<class PA, class GPA>
void MR3_FMM::calccoulomb(PA& particlearray,
			  const std::vector<TypeRange>& typerangearray,
			  std::vector<int>& self_longset_index,
			  std::vector<int>& self_selfenergycell_index,
			  GPA& ghost,
			  const std::vector<TypeRange>& ghosttyperange,
			  std::vector<int>& ghost_longset_index,
			  std::vector<int>& ghost_selfenergycell_index,
			  double& longenergy)
{
  SAMPLE_START();

  counter_start("calccoulomb");

  if(allocated==false){
    size_t size = particlearray.size();
    int level = int(log(self_longset_index.size()+1)/M_LN2/3);
    int images = FMM_IMAGES;
    parallelfmm->allocate(size,level,images);
    allocated = true;
    parallelfmm->printNow = false;
    if(parallelfmm->MPIRANK==0){
      printf("FMM Natom %ld Node(%d,%d,%d), Cell(%d,%d,%d), LocalCell(%d,%d,%d), LocalLevel %d, Cartesian order %d\n",
	     size,
	     numnode3d[0],numnode3d[1],numnode3d[2],
	     numcell3d[0],numcell3d[1],numcell3d[2],
	     numlocalcell3d[0],numlocalcell3d[1],numlocalcell3d[2],
	     level, P);
    }
  }
  
  for_3d parallelfmm->RGlob[d] = xmax[d] / 2;
  parallelfmm->R0 = xmax[0] / 2 / parallelfmm->numPartition[parallelfmm->maxGlobLevel][0];
  int ix[3];
  parallelfmm->getGlobIndex(ix,parallelfmm->MPIRANK,parallelfmm->maxGlobLevel);
  for_3d parallelfmm->X0[d] = 2 * parallelfmm->R0 * (ix[d] + .5);


  counter_start("copy");

  int rankOffset = 13 * parallelfmm->numLeafs;
  real diameter = 2 * parallelfmm->R0 / (1 << parallelfmm->maxLevel);
  for( int i=rankOffset; i<parallelfmm->numLeafs+rankOffset; i++ ) {
    parallelfmm->Leafs[i][0] = parallelfmm->Leafs[i][1] = 0;
  }
  parallelfmm->numBodies = 0;
  for( int ti=0; ti<self_longset_index.size(); ++ti ) {
    int t = self_longset_index[ti];
    real X[3] = {0, 0, 0};
    for( int i=typerangearray[t].begin; i<typerangearray[t].end; ++i ) {
      X[0] += getpos(particlearray,i).x;
      X[1] += getpos(particlearray,i).y;
      X[2] += getpos(particlearray,i).z;
    }
    for_3d X[d] /= typerangearray[t].end-typerangearray[t].begin;
    int ix[3];
    for_3d ix[d] = int((X[d] + parallelfmm->R0 - parallelfmm->X0[d]) / diameter);
    int ileaf = parallelfmm->getKey(ix,parallelfmm->maxLevel,false) + rankOffset;
    parallelfmm->Leafs[ileaf][0] = parallelfmm->numBodies;
    for( int i=typerangearray[t].begin; i<typerangearray[t].end; ++i, ++parallelfmm->numBodies ) {
      parallelfmm->Jbodies[parallelfmm->numBodies][0] = getpos(particlearray,i).x;
      parallelfmm->Jbodies[parallelfmm->numBodies][1] = getpos(particlearray,i).y;
      parallelfmm->Jbodies[parallelfmm->numBodies][2] = getpos(particlearray,i).z;
      parallelfmm->Jbodies[parallelfmm->numBodies][3] = getcharge(particlearray,i);
      for_4d parallelfmm->Ibodies[parallelfmm->numBodies][d] = 0;
    }
    parallelfmm->Leafs[ileaf][1] = parallelfmm->numBodies;
  }
  for( int i=rankOffset; i<parallelfmm->numLeafs+rankOffset; i++ ) {
    assert( parallelfmm->Leafs[i][1] != parallelfmm->Leafs[i][0] );
  }

  counter_stop("copy");


  counter_start("upward");

  parallelfmm->upwardPass();

  counter_stop("upward");


  counter_start("M2LSend");

  parallelfmm->M2LSend();

  counter_stop("M2LSend");


  counter_start("M2LRecv");

  parallelfmm->M2LRecv();

  counter_stop("M2LRecv");


  counter_start("rootGather");

  parallelfmm->rootGather();

  counter_stop("rootGather");


  counter_start("globM2M");

  parallelfmm->globM2M();

  counter_stop("globM2M");


  counter_start("globM2L");

  parallelfmm->globM2L();

  counter_stop("globM2L");


  counter_start("periodicM2L");

  parallelfmm->periodicM2L();

  counter_stop("periodicM2L");


  counter_start("globL2L");

  parallelfmm->globL2L();

  counter_stop("globL2L");



  counter_start("downward");

  parallelfmm->downwardPass();

  counter_stop("downward");

//  globDirect();



#ifdef DIPOLE_CORRECTION
  counter_start("dipole");
  
  double localdipole[3] = {0.0, 0.0, 0.0};
  for( int i=0; i<parallelfmm->numBodies; i++ ) {
    for_3d localdipole[d] += parallelfmm->Jbodies[i][3] * (parallelfmm->Jbodies[i][d] - parallelfmm->RGlob[d]);
  }
  double coef = M_PI / (6.0 * parallelfmm->RGlob[0] * parallelfmm->RGlob[1] * parallelfmm->RGlob[2]);
  double dipole[3] = {0.0,0.0,0.0};
  MPI_Allreduce(localdipole,dipole,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  int totalnumBodies =0 ;
  MPI_Allreduce(&(parallelfmm->numBodies),&totalnumBodies,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  //  double pc = .5 * coef
  //    * (dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]) / (double)(totalnumBodies);
  double pc = .5 * coef * (dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]) * (double)(parallelfmm->numBodies)/ (double)(totalnumBodies);
  for( int i=0; i<parallelfmm->numBodies; i++ ) {
    //    parallelfmm->Ibodies[i][0] += pc;
    for_3d parallelfmm->Ibodies[i][d+1] -= coef * dipole[d];
  }
  longenergy += pc;

  counter_stop("dipole");
#endif


  counter_start("addforce");

  int ibody = 0;
  for( int ti=0; ti<self_longset_index.size(); ++ti ) {
    int t = self_longset_index[ti];
    for( int i=typerangearray[t].begin; i<typerangearray[t].end; ++i, ++ibody ){
      longenergy += .5 * parallelfmm->Jbodies[ibody][3] * parallelfmm->Ibodies[ibody][0];
      Force &pif = getforce(particlearray,i);
      pif.x -= parallelfmm->Jbodies[ibody][3] * parallelfmm->Ibodies[ibody][1];
      pif.y -= parallelfmm->Jbodies[ibody][3] * parallelfmm->Ibodies[ibody][2];
      pif.z -= parallelfmm->Jbodies[ibody][3] * parallelfmm->Ibodies[ibody][3];
    }
  }

  counter_stop("addforce");

//  deallocate();


  counter_stop("calccoulomb");

  SAMPLE_STOP();

}

template
void MR3_FMM::calccoulomb(CombinedParticleArray& particlearray,
			  const std::vector<TypeRange>& typerangearray,
			  std::vector<int>& self_longset_index,
			  std::vector<int>& self_selfenergycell_index,
			  GhostParticleArray& ghost,
			  const std::vector<TypeRange>& ghosttyperange,
			  std::vector<int>& ghost_longset_index,
			  std::vector<int>& ghost_selfenergycell_index,
			  double& longenergy);
