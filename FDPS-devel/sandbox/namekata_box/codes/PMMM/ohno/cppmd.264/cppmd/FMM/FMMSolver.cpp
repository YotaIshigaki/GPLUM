/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, Inc             *
 *               Ohno Yousuke, PhD at RIKEN                  *
 *                                                           *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *         ~~~  FMM Solver Class. Source File.  ~~~          *
 *                                                           *
 *************************************************************/
#include "FMMSolver.h"
#include "string.h"                   //  for memset ()
#include <iostream>                   //  for cout, etc.



extern "C++" {

//  privided by:  upward_pass.cpp
void  upward_pass (
          Cell  ** & cell_tree,                       //  [in/out]  Tree of cells
          double position[],                          //  [in]  Array containing all particles
          double charges [],                          //  [in]  Array containing all charges
          const int & finest_level,                   //  [in]  Finest level of the tree hierarchy
          const unsigned long & total,                //  [in]  Total number of particles
          const Coefficients * coeffs                 //  [in]  Pointer to the coefficients structure (with pre-computed coefficients)
        );

//  privided by:  upward_pass.cpp
void  upward_pass_mpi (
          Cell  **  & cell_tree,                      //  [in/out]  Tree of cells
          double position[],                          //  [in]  Array containing all particles
          double charges [],                          //  [in]  Array containing all charges
          const int & myrank,                         //  [in]  rank
          const int & n_active,                       //  [in]  number of active ranks
          const int & rank_finest_level,              //  [in]  rank-finest level
          const int & finest_level,                   //  [in]  Finest level of the tree hierarchy
          const unsigned long & total,                //  [in]  Total number of particles
          const Coefficients * coeffs,                //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
          MPI_Comm fmm_world = MPI_COMM_WORLD         //  [in]  Communicator
        );

//  privided by:  downward_pass.cpp
void  downward_pass (
          double & total_potential,                   //  [out] computed total potential
#ifdef  USE_POTENTIAL_ARRAY
          double  potential [],                       //  [out] computed potential for each particle
#endif
          double  force [],                           //  [out] computed force (in triples (Fx, Fy, Fz).
          Cell  ** &  cell_tree,                      //  [in/out]  Tree of cells
          int     finest_level,                       //  [in]  Finest level of the tree hierarchy
          double  position [],                        //  [in]  Array containing all particles
          double  charges [],                         //  [in]  Array containing all charges
          const Coefficients * coeffs,                //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
          std::complex <double> * sph_harmonic_storage,    //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
          const BOUNDARY_CONDITION &
                boundary_condition =
                DEFAULT_BOUNDARY_CONDITION            //  [in]  Boundary Condition
        );

//  privided by:  downward_pass.cpp
void  downward_pass_mpi (
          double & total_potential,                   //  [out] computed total potential
#ifdef  USE_POTENTIAL_ARRAY
          double  potential [],                       //  [out] computed potential for each particle
#endif
          double  force [],                           //  [out] computed force (in triples (Fx, Fy, Fz).
          Cell  ** &  cell_tree,                      //  [in/out]  Tree of cells
          int     myrank,                             //  [in]  Rank
          int     rank_finest_level,                  //  [in]  Rank-finest level
          int     finest_level,                       //  [in]  Finest level of the tree hierarchy
          double  position [],                        //  [in]  Array containing all particles
          double  charges [],                         //  [in]  Array containing all charges
          const Coefficients * coeffs,                //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
          std::complex <double> * sph_harmonic_storage,    //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
          const BOUNDARY_CONDITION &
                boundary_condition =
                DEFAULT_BOUNDARY_CONDITION,           //  Boundary Condition
          MPI_Comm fmm_world = MPI_COMM_WORLD         //  Communicator
        );

}  //  extern "C++"



//  ====================================================================  Sequential FMM Solver  =========

FMMSolver :: FMMSolver()
{
  p = 0;
  finest_level = 0;
  cell_tree = NULL;
  coeffs = NULL;
  spherical_harmonics = NULL;
  _positions = NULL;
  _charges = NULL;
//  Added by Gleb. 2011.06.02
  set_scale (1.0);
//	Addition ends
}

FMMSolver :: FMMSolver (int L, int expansion_order) :  finest_level(L), p(expansion_order)
{

    set_scale (1.0);

    Cell :: set_expansion_order (p);
    coeffs = new Coefficients (p);
//    ::coeffs = coeffs;

    spherical_harmonics = new std::complex <double> [(2*p+1)*(2*p+1)];
//    ::sph_harmonic_storage = spherical_harmonics;

    fmm_create_tree ();

    std::cout << "Sequential FMMSolver Created :: L="<<finest_level<<", p="<<p<<std::endl;

}

void FMMSolver :: destroyParticles() {

  if (_positions) {
    delete [] _positions;
    _positions = NULL;
  }

  if (_charges) {
    delete [] _charges;
    _charges = NULL;
  }

}

FMMSolver ::  ~FMMSolver() {
    if (cell_tree)
      fmm_delete_tree ();

    if (coeffs) {
      delete coeffs;
      coeffs = NULL;
    }

    if (spherical_harmonics) {
      delete spherical_harmonics;
      spherical_harmonics = NULL;
    }


    destroyParticles();
    std::cout <<"FMMSolver Destroyed."<<std::endl;
}


//  Copies particles data to internal storage.  In case of MPI, broadcast this data
//  to other slave ranks.
void  FMMSolver :: setParticles (
            unsigned long  total,                   // [in]
            double positions     [],                // [in]
            double charges       []                 // [in]
        )
{

  _total = total;

//  Allocating memory for positions/charges

  _positions = new double [_total*3];
  _charges   = new double [_total];

//  copy it to internal storage in order to decouple from input arrays
  memcpy (_positions, positions, sizeof (double)*_total*3);
  memcpy (_charges,   charges,   sizeof (double)*_total);

}


void  FMMSolver :: setParticles (
            ParticleArray & particles               // [in]
        )
{
  _total  =  particles.size();

  _positions = new double [_total*3];
  _charges   = new double [_total];

//  Copying necessary data from "ParticleArray & particles"
    for (unsigned long  i=0; i<_total; i++) {
      Particle & particle = particles [i];

      _positions [3*i + 0] = particle.position.x*scale;
      _positions [3*i + 1] = particle.position.y*scale;
      _positions [3*i + 2] = particle.position.z*scale;
      _charges [i] = particle.charge;
    }
}



void  FMMSolver :: RunSolver (
#ifdef  USE_POTENTIAL_ARRAY
            double potential     [],                // [out]
#endif
            double forces        [],                // [out]
            double & total_potential,               // [out]
            BOUNDARY_CONDITION  boudary_condition   // [in]
        )
{

    double total_pot_test = 0;

    std::cout <<"Sequential FMM :: UPWARD PASS.       ====="<<std::endl;

    upward_pass (
                cell_tree,
                _positions,
                _charges,
                finest_level,
                _total,
                coeffs
              );


    std::cout <<"Sequential FMM :: DOWNWARD PASS.     ====="<<std::endl;

    downward_pass (
#ifdef  USE_POTENTIAL_ARRAY
                total_pot_test,               // scalar value of total potential
                potential,                    // array of potentials for each particle
#else
                total_potential,              // scalar value of total potential
#endif
                forces,
                cell_tree,
                finest_level,
                _positions,
                _charges,
                coeffs,
                spherical_harmonics,
                boudary_condition
              );



#ifdef  USE_POTENTIAL_ARRAY
    double pot = 0;
    for (unsigned long  i=0; i<_total; i++)
      pot += potential [i];

    if (pot != total_pot_test) {
      fprintf (stdout, "Sequential Solver :: ERROR.  total_pot_test = %e, pot_sum = %e. delta = %e, rel_delta = %e\n",
         total_pot_test, pot, fabs (pot - total_pot_test), fabs ((pot - total_pot_test)/pot));  fflush (stdout);
    }
    total_potential = pot;
#endif

    total_potential *= scale;
    scale_force(forces);
  fprintf (stdout, "Sequential Solver :: Total Potential = %e\n", total_potential); fflush (stdout);

}


void  FMMSolver :: fmm_create_tree ()
{
  int total_levels = finest_level + 1;
  cell_tree = new Cell* [total_levels];

  for (int l = 0; l<total_levels; l++) {
    ucellindex_t cells_per_level = 1 << (3*l);  // = 8^l
    cell_tree[l] = new Cell [cells_per_level];
  }

}


void  FMMSolver :: fmm_clear_cell_tree_particles ()
{

// Particles positions/charges are stored in the cells only on the finest level,
// so it should be only 1 loop.  Also, one must call cell::destroy_particles_data() to
// deallocate arrays for particle positions/charges

  for(ucellindex_t cl = 0; cl<1<<(3*finest_level); cl++) {
    cell_tree[finest_level][cl].particles.clear();
    cell_tree[finest_level][cl].destroy_particles_data();
  }


#if 0
  int total_levels = finest_level + 1;
  for (int l = 0; l<total_levels; l++) {
    ucellindex_t cells_per_level = 1 << (3*l);  // = 8^l
    for(ucellindex_t cl = 0;cl<cells_per_level;cl++){
      cell_tree[l][cl].particles.clear();
    }
  }
#endif

}

void  FMMSolver :: fmm_clear_cell_tree_MECs ()
{
//  for debug
  return;
//  end for debug

  int total_levels = finest_level + 1;
  for (int l = 0; l<total_levels; l++) {
    ucellindex_t cells_per_level = 1 << (3*l);  // = 8^l
    for(ucellindex_t cl = 0;cl<cells_per_level;cl++){
      int num_coeff = cell_tree[l][cl].get_num_coeff();
      for(int c=0;c<num_coeff;c++){
        cell_tree[l][cl].MEC[c]=std::complex<double>(0.0,0.0);
      }
    }
  }

}


void  FMMSolver :: fmm_delete_tree ()
{
  for (int l = 0; l < finest_level + 1 /* total levels*/; l++)
    delete [] cell_tree[l];

  delete [] cell_tree;
  cell_tree = NULL;

}

unsigned long FMMSolver :: get_total() {
  return _total;
}

void FMMSolver :: set_scale(const double  _scale) {
  scale = _scale;
  fscale = scale*scale*scale;
}

void FMMSolver :: scale_force(double *force) {
  for(unsigned long i=0; i<_total*3; i++){
    force[i] *= fscale;
  }
}

//  ====================================================================  Parallel FMM Solver  =========


FMMSolverParallel :: FMMSolverParallel (int L, int expansion_order, MPI_Comm lcomm)
{
    finest_level = L;
    p = expansion_order;
    myworld = lcomm;


//  Added by Gleb. 2011.06.02
    set_scale (1.0);
    _positions = NULL;
    _charges = NULL;
//	Addition ends



    if(myworld==MPI_COMM_NULL){
      return;
    }
    MPI_Bcast(&finest_level, 1, MPI_INT, 0 /* ROOT_RANK */, myworld);
    MPI_Bcast(&p,            1, MPI_INT, 0 /* ROOT_RANK */, myworld);

//  Determining the rank-finest level and the number of active ranks
    MPI_Comm_rank (myworld, &myrank);
    MPI_Comm_size (myworld, &numproc);

    rank_finest_level = int (log2(numproc)/3.0);
    if (rank_finest_level > finest_level)    rank_finest_level = finest_level;
    n_active = 1 << (3*rank_finest_level);
    printf("rank_finest_level %d  n_active %d  expansion_orde %d\n",rank_finest_level,n_active,p);

    Cell :: set_expansion_order (p);
    printf("fmm_create_tree\n");
    fmm_create_tree ();

    printf("new Coefficients\n");
    coeffs = new Coefficients (p);
//    ::coeffs = coeffs;

    printf("spherical_harmonics\n");
    spherical_harmonics = new std::complex <double> [(2*p+1)*(2*p+1)];
//    ::sph_harmonic_storage = spherical_harmonics;

    std::cout << "sizeof(cellindex_t) " << sizeof(cellindex_t) << std::endl;

    std::cout <<"RANK "<<myrank<<" :: Parallel FMMSolver Created :: L="<<finest_level<<", p="<<p<<std::endl;
}

FMMSolverParallel ::  ~FMMSolverParallel() {
  if(myworld!=MPI_COMM_NULL){
    std::cout <<"RANK "<<myrank<<" :: Parallel "; // the rest of the message is printed by ~FMMSolver()
  }
}


int   FMMSolverParallel :: get_n_active () { return n_active; }



//  Copies particles data to internal storage, then broadcasts this data
//  to other (slave) ranks.

void  FMMSolverParallel :: setParticles (
            unsigned long total,                 // [in]
            double * positions,                  // [in]
            double * charges                     // [in]
        )
{

    if (myrank == 0)
      _total = total;

    MPI_Bcast(&_total, 1, MPI_LONG, 0 /* ROOT_RANK */, myworld);


//  Allocating memory for positions/charges

    _positions = new double [_total*3];
    _charges   = new double [_total];

//  Copy it to internal storage in order to decouple from input arrays
    if (myrank == 0) {
      memcpy (_positions, positions, sizeof (double)*_total*3);
      memcpy (_charges,   charges,   sizeof (double)*_total);
    }

//  Broadcasting particles positions, charges from the master rank
    MPI_Bcast(_positions, _total*3, MPI_DOUBLE, 0 /*ROOT_RANK*/, myworld);
    MPI_Bcast(_charges,   _total,   MPI_DOUBLE, 0 /*ROOT_RANK*/, myworld);

//  Release positions/charges memory for the inactive ranks (ranks
//  which do not participate in computations)
    if (myrank >= n_active) {
      delete [] _positions; _positions = NULL;
      delete [] _charges;    _charges  = NULL;
    }
}


void  FMMSolverParallel :: setParticles (
            ParticleArray & particles               // [in]
        )
{
  if (myrank == 0) {
    _total  =  particles.size();
  }
    MPI_Bcast(&_total, 1, MPI_LONG, 0 /* ROOT_RANK */, myworld);


  _positions = new double [_total*3];
  _charges   = new double [_total];


//  Copying necessary data from "ParticleArray & particles"
  if (myrank == 0) {
    for (unsigned long  i=0; i<_total; i++) {
      Particle & particle = particles [i];

      _positions [3*i + 0] = particle.position.x*scale;
      _positions [3*i + 1] = particle.position.y*scale;
      _positions [3*i + 2] = particle.position.z*scale;
      _charges [i] = particle.charge;
    }
  }

//  Broadcasting particles positions, charges from the master rank
  MPI_Bcast(_positions, _total*3, MPI_DOUBLE, 0 /*ROOT_RANK*/, myworld );
  MPI_Bcast(_charges,   _total,   MPI_DOUBLE, 0 /*ROOT_RANK*/, myworld );

//  Release positions/charges memory for the inactive ranks
//  (ranks which do not participate in computations)
  if (myrank >= n_active) {
    delete [] _positions; _positions = NULL;
    delete [] _charges;   _charges   = NULL;
  }

}

template<class PA, class GPA>
unsigned long FMMSolverParallel :: allocate_particles(const PA& particlearray,
                                    const std::vector<TypeRange>& typerange,
                                    const std::vector<int>& self_longset_index,
                                    const GPA& ghost,
                                    const std::vector<TypeRange>& ghosttyperange,
                                    const std::vector<int>& ghost_longset_index)
{
  int num=0;

  if (myrank == 0) {
    for(size_t si=0;si<self_longset_index.size();si++){
      num += (typerange[si].end - typerange[si].begin);
    }
    for(size_t si=0;si<ghost_longset_index.size();si++){
      num += (ghosttyperange[si].end - ghosttyperange[si].begin);
    }
    _total = num;
  }

  MPI_Bcast(&_total, 1, MPI_LONG, 0 /* ROOT_RANK */, myworld);

  _positions = new double [_total*3];
  _charges   = new double [_total];

  update_particles(particlearray,
                   typerange,
                   self_longset_index,
                   ghost,
                   ghosttyperange,
                   ghost_longset_index);

  return _total;
}


template<class PA, class GPA>
void FMMSolverParallel :: update_particles (const PA& particlearray,
                                 const std::vector<TypeRange>& typerange,
                                 const std::vector<int>& self_longset_index,
                                 const GPA& ghost,
                                 const std::vector<TypeRange>& ghosttyperange,
                                 const std::vector<int>& ghost_longset_index)
{
  if (myrank == 0) {
    int n = 0;
    for(size_t si=0;si<self_longset_index.size();si++){
      for(int i=typerange[si].begin;i<typerange[si].end;i++){
        const Position & pos = getpos(particlearray,i);
        _positions [3*n + 0]  = pos.x*scale;
        _positions [3*n + 1]  = pos.y*scale;
        _positions [3*n + 2]  = pos.z*scale;
        if(_positions [3*n + 0]<0.0)_positions [3*n + 0] +=1.0;
        if(_positions [3*n + 0]>1.0)_positions [3*n + 0] -=1.0;
        if(_positions [3*n + 1]<0.0)_positions [3*n + 1] +=1.0;
        if(_positions [3*n + 1]>1.0)_positions [3*n + 1] -=1.0;
        if(_positions [3*n + 2]<0.0)_positions [3*n + 2] +=1.0;
        if(_positions [3*n + 2]>1.0)_positions [3*n + 2] -=1.0;
        _charges [n] = getcharge(particlearray,i);
        n++;
      }
    }
    for(size_t si=0;si<ghost_longset_index.size();si++){
      for(int i=ghosttyperange[si].begin;i<ghosttyperange[si].end;i++){
        const Position & pos = getpos(ghost,i);
        _positions [3*n + 0]  = pos.x*scale;
        _positions [3*n + 1]  = pos.y*scale;
        _positions [3*n + 2]  = pos.z*scale;
        if(_positions [3*n + 0]<0.0)_positions [3*n + 0] +=1.0;
        if(_positions [3*n + 0]>1.0)_positions [3*n + 0] -=1.0;
        if(_positions [3*n + 1]<0.0)_positions [3*n + 1] +=1.0;
        if(_positions [3*n + 1]>1.0)_positions [3*n + 1] -=1.0;
        if(_positions [3*n + 2]<0.0)_positions [3*n + 2] +=1.0;
        if(_positions [3*n + 2]>1.0)_positions [3*n + 2] -=1.0;
        _charges [n] = getcharge(ghost,i);
        n++;
      }
    }
    if(n!=_total){
      printf("update %d particle but total %d particle\n",n,_total);
    }
  }

  //  Broadcasting particles positions, charges from the master rank
  MPI_Bcast(_positions, _total*3, MPI_DOUBLE, 0 /*ROOT_RANK*/, myworld );
  MPI_Bcast(_charges,   _total,   MPI_DOUBLE, 0 /*ROOT_RANK*/, myworld );

}

#ifdef OLDPARTICLE
template
void FMMSolverParallel :: update_particles (       const ParticleArray& particlearray,
                                const std::vector<TypeRange>& typerange,
                                const std::vector<int>& self_longset_index,
                                const ParticleArray& ghost,
                                const std::vector<TypeRange>& ghosttyperange,
                                const std::vector<int>& ghost_longset_index
                              );

template
unsigned long FMMSolverParallel :: allocate_particles (     const ParticleArray& particlearray,
                                         const std::vector<TypeRange>& typerange,
                                         const std::vector<int>& self_longset_index,
                                         const ParticleArray& ghost,
                                         const std::vector<TypeRange>& ghosttyperange,
                                         const std::vector<int>& ghost_longset_index
                                         );
#else
template
void FMMSolverParallel :: update_particles (       const CombinedParticleArray& particlearray,
                                const std::vector<TypeRange>& typerange,
                                const std::vector<int>& self_longset_index,
                                const GhostParticleArray& ghost,
                                const std::vector<TypeRange>& ghosttyperange,
                                const std::vector<int>& ghost_longset_index
                              );

template
unsigned long FMMSolverParallel :: allocate_particles (     const CombinedParticleArray& particlearray,
                                         const std::vector<TypeRange>& typerange,
                                         const std::vector<int>& self_longset_index,
                                         const GhostParticleArray& ghost,
                                         const std::vector<TypeRange>& ghosttyperange,
                                         const std::vector<int>& ghost_longset_index
                                         );
#endif



void  FMMSolverParallel :: RunSolver (
#ifdef  USE_POTENTIAL_ARRAY
            double potential     [],                // [out]
#endif
            double forces        [],                // [out]
            double & total_potential,               // [out]
            BOUNDARY_CONDITION  boudary_condition   // [in]
        )
{

#ifdef  USE_POTENTIAL_ARRAY
    double * l_potential    = new double [_total];
    memset (l_potential,  0,  _total*sizeof(double));
#endif

    double   l_total_pot = 0;
    double * l_force        = new double [3*_total];

    memset (l_force,      0,  _total*3*sizeof(double));

  if (myrank < n_active) {

    //         fprintf (stdout, "RANK %6d :: UPWARD PASS.    Parallel FMM ========\n", myrank);  fflush (stdout);
      fmm_clear_cell_tree_particles();
      fmm_clear_cell_tree_MECs();
      upward_pass_mpi (
                cell_tree,
                _positions,
                _charges,
                myrank,
                n_active,
                rank_finest_level,
                finest_level,
                _total,
                coeffs,
                myworld
              );


      //            fprintf (stdout, "RANK %6d :: DOWNWARD PASS.  Parallel FMM ========\n", myrank);  fflush (stdout);
      downward_pass_mpi (
                l_total_pot,                  // scalar value of total potential
#ifdef  USE_POTENTIAL_ARRAY
                l_potential,                  // array of potentials for each particle
#endif
                l_force,
                cell_tree,
                myrank,
                rank_finest_level,
                finest_level,
                _positions,
                _charges,
                coeffs,
                spherical_harmonics,
                boudary_condition,
                myworld
              );

  }
  scale_force(l_force);

#if 0
  {
    unsigned long i;
    double minp[3]={0.0,0.0,0.0}, maxp[3]={0.0,0.0,0.0};
    double minf[3]={0.0,0.0,0.0}, maxf[3]={0.0,0.0,0.0};


    for(i=0;i<_total; i++){
      for(int j=0;j<3;j++){
        if(_positions[i*3+j]<minp[j])minp[j]=_positions[i*3+j];
        if(_positions[i*3+j]>maxp[j])maxp[j]=_positions[i*3+j];
        if(l_force[i*3+j]<minf[j])minf[j]=l_force[i*3+j];
        if(l_force[i*3+j]>maxf[j])maxf[j]=l_force[i*3+j];
      }
    }
    fprintf (stdout, "RANK %6d :: position range (%f,%f,%f)(%f,%f,%f)\n", myrank, minp[0],minp[1],minp[2],maxp[0],maxp[1],maxp[2]); fflush (stdout);
    fprintf (stdout, "RANK %6d :: force range (%f,%f,%f)(%f,%f,%f)\n", myrank, minf[0],minf[1],minf[2],maxf[0],maxf[1],maxf[2]); fflush (stdout);
  }
#endif


#ifdef  USE_POTENTIAL_ARRAY
  double pot = 0;
  for (unsigned long  i=0; i<_total; i++)
    pot += l_potential [i];

  if (pot != l_total_pot) {
    fprintf (stdout, "RANK %6d :: ERROR.  l_total_pot = %e,\tpot_sum = %e. delta = %e, rel_delta = %e\n",
      myrank, l_total_pot, pot, fabs (pot - l_total_pot), fabs ((pot - l_total_pot)/pot));  fflush (stdout);
  }
#endif


  MPI_Reduce (l_force    /*in*/,  forces    /*out*/,       _total*3,     MPI_DOUBLE,  MPI_SUM,  0 /*ROOT*/,  myworld );

  total_potential = 0.0;

#ifdef  USE_POTENTIAL_ARRAY
  MPI_Reduce (l_potential/*in*/,  potential /*out*/,         _total,     MPI_DOUBLE,  MPI_SUM,  0 /*ROOT*/,  myworld );
  pot *= scale;
  //  fprintf (stdout, "RANK %6d :: local Potential = %e\n", myrank, pot);
  MPI_Reduce (&pot/*in*/,         &total_potential /*out*/,       1,     MPI_DOUBLE,  MPI_SUM,  0 /*ROOT*/,  myworld );
#else
  l_total_pot *= scale;
  //  fprintf (stdout, "RANK %6d :: local Potential = %e\n", myrank, l_total_pot);
  MPI_Reduce (&l_total_pot/*in*/, &total_potential /*out*/,       1,     MPI_DOUBLE,  MPI_SUM,  0 /*ROOT*/,  myworld );
#endif

  total_potential *= 0.5;
  /*
  if(myrank==0){
    fprintf (stdout, "RANK %6d :: Total Potential = %e\n", myrank, total_potential); fflush (stdout);
  }
  */

#ifdef  USE_POTENTIAL_ARRAY
  delete [] l_potential;
#endif

  delete [] l_force;

}

