#ifndef FMMLONGRANGEINTERACTION_H
#define FMMLONGRANGEINTERACTION_H

#include <FMM/FMMSolver.h>

class FMMLongRangeInteraction {
 public:
  FMMLongRangeInteraction(int unitid, const LongRangeParameter& param, MPI_Comm lcomm=MPI_COMM_WORLD)
      : unit_identifier(unitid),
        myworld(lcomm),
        L(3), // TODO : get from LongRangeParameter
        expansion_order(10), // TODO get from LongRangeParameter
        fmmsolveparallel(L,expansion_order,myworld),
        //fmmsolveparallel(L,expansion_order,MPI_COMM_NULL),
        particle_allocated(false),
        boxsize(param.boxSize)
  {
    if(myworld!=MPI_COMM_NULL){
      MPI_Comm_rank (myworld, &myrank);
    }else{
      myrank = -1;
    }
    potential = NULL;
    forces    = NULL;
#ifndef FMM_PERIODIC
    boundary_condition = ZERO;
#else
    boundary_condition = PBC;
#endif
    double size = std::max(boxsize.x,boxsize.y);
    size = std::max(size,boxsize.z);
    scale = 1.0/size;
    fmmsolveparallel.set_scale(scale);
  }

  ~FMMLongRangeInteraction()
  {
    delete [] forces;
    delete [] potential;
  }

  template<class PA, class GPA>
  void calcForce(PA& particlearray,
                 const std::vector<TypeRange>& typerange,
                 const std::vector<int>& self_longset_index,
                 const std::vector<int>& self_selfenergycell_index,
                 GPA& ghost,
                 const std::vector<TypeRange>& ghosttyperange,
                 const std::vector<int>& ghost_longset_index,
                 const std::vector<int>& ghost_selfenergycell_index,
                 double& energy)
  {
#if 1
    if(particle_allocated){
      /*
      if(myrank==0){
        printf("FMM particle allcated previous step\n");
      }
      */
      fmmsolveparallel.update_particles(particlearray,
                                      typerange,
                                      self_longset_index,
                                      ghost,
                                      ghosttyperange,
                                      ghost_longset_index);
    }else{
      unsigned long np=0;
      np = fmmsolveparallel.allocate_particles(particlearray,
                                               typerange,
                                               self_longset_index,
                                               ghost,
                                               ghosttyperange,
                                               ghost_longset_index);
      particle_allocated = true;
      if(myrank==0){
        forces = new double [np*3];
        printf("number of FMM particle %ld\n",np);
      }
    }
#endif
#if 1
    fmmsolveparallel.RunSolver(
                             #ifdef  USE_POTENTIAL_ARRAY
                               potential,
                             #endif
                               forces,
                               total_potential,
                               boundary_condition);
#endif
#if 0
    if(myrank==0){
      printf("FMM force(0) (%f,%f,%f)\n",forces[0],forces[1],forces[2]);
    }
#endif
#if 1
    restore_forces(particlearray,
                   typerange,
                   self_longset_index,
                   ghost,
                   ghosttyperange,
                   ghost_longset_index);
#endif
    energy += total_potential;
  }
 private:
  int unit_identifier;
  MPI_Comm myworld;
  int myrank;
  int L;
  int expansion_order;
  FMMSolverParallel fmmsolveparallel;
/*

//  error: ISO C++ forbids initialization of member potential
  double *potential = NULL;
  double *forces = NULL;
*/
  double *potential;
  double *forces;

  double total_potential;
  BOUNDARY_CONDITION  boundary_condition;
  bool particle_allocated;

  SpaceVector<double> boxsize;

  double scale;

  template<class PA, class GPA>
  void restore_forces(PA& particle,
                      const std::vector<TypeRange>& typerange,
                      const std::vector<int>& self_longset_index,
                      GPA& ghost,
                      const std::vector<TypeRange>& ghosttyperange,
                      const std::vector<int>& ghost_longset_index)
  {
    int fi=0;
    for(size_t si=0;si<self_longset_index.size();si++){
      for(int i=typerange[si].begin;i<typerange[si].end;i++){
        /*
        particle[i].force.x += forces[fi*3];
        particle[i].force.y += forces[fi*3+1];
        particle[i].force.z += forces[fi*3+2];
        */
        Force &pif = getforce(particle,i);
        pif.x = forces[fi*3];
        pif.y = forces[fi*3+1];
        pif.z = forces[fi*3+2];
        fi++;
      }
    }
    for(size_t si=0;si<ghost_longset_index.size();si++){
      for(int i=ghosttyperange[si].begin;i<ghosttyperange[si].end;i++){
        /*
        ghost[i].force.x += forces[fi*3];
        ghost[i].force.y += forces[fi*3+1];
        ghost[i].force.z += forces[fi*3+2];
        */
        Force &gif = getforce(ghost,i);
        gif.x = forces[fi*3];
        gif.y = forces[fi*3+1];
        gif.z = forces[fi*3+2];
        fi++;
      }
    }
  }
};


#endif
