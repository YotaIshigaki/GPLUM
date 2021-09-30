/*
 *  Staggard Mesh Ewald Method for Long Range Interaction Force 
 *  (total energy is also computed) 
 *  MPI Parallelization designed
 *
 *
 *
 *    mesh configulation c: cell center mesh point, v: cell vertex mesh point
 *
 *    +---+---+---+---+
 *    | c | c | c | c |
 *    v---v---v---v---+          North
 *    | c | c | c | c |            Y   Z up
 *    v---v---v---v---+            | / 
 *    | c | c | c | c |            |/
 *    v---v---v---v---+   West ----+----> X East  
 *    | c | c | c | c |            |
 *    v---v---v---v---+          South
 *    1   2   3   4   5
 *   Make sure that the corner is starting with 1 (not 0)
 *
 *
 *   CAUTION
 *   two ghost meshes are used, i.e., overlap(1,1,1), 
 *   overlap(Ncut,Ncut,Ncut) where Ncut is a cut off meshes for 
 *   chargeassignment/backinterpolation
 *
 *
 *
 *
 *  To Do list
 *  multigrid
 *  chargeassign/backinterpolation exponential function pre-computed
 *  Poisson solver twice to once
 *  use old potential as an initial guess
 *
 *
 *
 */
#ifndef StMELONGRANGEINTERACTION_H
#define StMELONGRANGEINTERACTION_H

#include <mpi.h>
#include <cmath>
#include "LongRangeInteraction.h"

/*     0   1  2  3  4 5 6
 *   256 128 64 32 16 8 4 
 *     8   7  6  5  4 3 2
 *     
 */
typedef struct Grid_s{
  double *u, *f, *res;
  SpaceVector<int> n;
  SpaceVector<double> dmesh;
  double h1sqinv, h2sqinv, h3sqinv;
  double hhh;
  double hhxy, hhzx, hhyz;
  double hx, hy, hz;
} grid;


#define MAX_LEVEL 8
#define MAX_N  257 // MAX_N = 2^MAX_LEVEL
#define MAX_NN  (MAX_N*MAX_N)
#define MAX_NNN  (MAX_N*MAX_N*MAX_N)

#define SQR(x) ((x)*(x)) 
#define CUBE(x) ((x)*(x)*(x)) 
#define MAX(a,b) ((a > b) ? a : b) 
#define MIN(a,b) ((a < b) ? a : b) 

#define FI3D(i,j,k,n1,n2,n3) ((k) + (n3)*((j) + (n2)*(i)))
#define FI2D(i,j,n1,n2) ((j) + ((n2)*(i))

class StMELongRangeInteraction {
 public:
  double alpha; // = 0.3;//default value
  int Niter;

  StMELongRangeInteraction(int unitid, const LongRangeParameter& param, MPI_Comm lcomm=MPI_COMM_WORLD)
    :     alpha(param.alpha),
          Niter(160),
          unit_identifier(unitid),
    myworld(lcomm),
    boxsize(param.boxSize),
    node_geometry(param.node_geometry)

      {
	if(myworld!=MPI_COMM_NULL){
	  MPI_Comm_rank (myworld, &myrank);
	}else{
	  myrank = -1;
	}
	potential = NULL;
	forces    = NULL;

        /*
          node_geometry (copy of param.node_geometry)
          size : number of long node (x,y,z)
          have_cell : number of "short cell" (x,y,z) in on "long node"
          full_cell : number of "short cell" (x,y,z) in full system
          short cell information is set at 
              makelongrangegeometry in CalcPreparator.cc
         */
        local_mesh = node_geometry.have_cell*4;
        global_mesh = node_geometry.full_cell*4;
        SpaceVector<int> node_pos = node_geometry.getAbsolutePosition(myrank);
        mesh_begin.x = node_pos.x*local_mesh.x;
        mesh_begin.y = node_pos.y*local_mesh.y;
        mesh_begin.z = node_pos.z*local_mesh.z;

        std::cout << "mesh size" << local_mesh << global_mesh << std::endl;
        std::cout << "node pos" << node_pos << "  mesh_begin" << mesh_begin << std::endl;
        if(local_mesh.x>MAX_NNN){
          printf("local_mesh.x %d > MAX_NNN %d\n",local_mesh.x,MAX_NNN);
        }
        if(local_mesh.y>MAX_NNN){
          printf("local_mesh.y %d > MAX_NNN %d\n",local_mesh.y,MAX_NNN);
        }
        if(local_mesh.z>MAX_NNN){
          printf("local_mesh.z %d > MAX_NNN %d\n",local_mesh.z,MAX_NNN);
        }

        init_param();
    grd.u = new double [MAX_NNN];
    grd.f = new  double [MAX_NNN];
    grd.res = new double [MAX_NNN];
      }

    ~StMELongRangeInteraction()
      {
        delete [] grd.res;
        delete [] grd.f;
        delete [] grd.u;
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
#ifdef DEBUG_VERB
      printf("Staggard Mesh Ewald Method \n");
#endif
      total_potential = 0.0;
      chargeassignment(particlearray,
		       typerange,
		       self_longset_index,
		       ghost,
		       ghosttyperange,
		       ghost_longset_index);
    

      RBGSPoissonSolver(rhoc, potc, local_mesh, dmesh, Niter);
      RBGSPoissonSolver(rhov, potv, local_mesh, dmesh, Niter);
      //      MultigridPoissonSolver(rhoc, potc, local_mesh, dmesh, nPREsmooth, nPOSTsmooth, NCYCLE);
      //      MultigridPoissonSolver(rhov, potv, local_mesh, dmesh, nPREsmooth, nPOSTsmooth, NCYCLE);

      backinterpolation(particlearray,
			typerange,
			self_longset_index,
			ghost,
			ghosttyperange,
			ghost_longset_index);
      
      if(unit_identifier==0){
        printf("STME total_potential %e\n", total_potential);
      }

      energy += total_potential;
    }
 private:
    int unit_identifier;
    MPI_Comm myworld;
    SpaceVector<double> boxsize;
    GeometryXYZ node_geometry;
    int myrank;
    double *potential;
    double *forces;

    double total_potential;

    double rhoc[MAX_NNN];
    double potc[MAX_NNN];
    double rhov[MAX_NNN];
    double potv[MAX_NNN];
    
    SpaceVector<int> local_mesh, mesh_begin, global_mesh;
  
    double kappa, lambda;
    //  const double sigma = 3.0/sqrt(2.0);
    double sigma;
    double sigma1;
    //  const double sigma2 = sigma*sqrt(1.0 - 2.0*kappa);
    double Rcut1;
    //const double Rcut2 = lambda*sqrt(2.0)*sigma2;

    SpaceVector<double> dmesh;
    SpaceVector<int> Ncut;
    int nPREsmooth, nPOSTsmooth, NCYCLE;
  
double sendbuff[MAX_NNN], recvbuff[MAX_NNN];

grid grd;

  // used in backinterpolation
      double potv_with_cut[MAX_NNN];
      double potc_with_cut[MAX_NNN];

    void init_param()
    {
      kappa = 0.5;
      lambda = 3.0;
      //  const double sigma = 3.0/sqrt(2.0);
      sigma = 1.0/(sqrt(2.0)*alpha);
      sigma1 = sigma*sqrt(kappa);
      //  const double sigma2 = sigma*sqrt(1.0 - 2.0*kappa);
      Rcut1 = lambda*sqrt(2.0)*sigma1;
      //const double Rcut2 = lambda*sqrt(2.0)*sigma2;
      
      dmesh = 
        SpaceVector<double>(boxsize.x/(double)global_mesh.x, boxsize.y/(double)global_mesh.y, boxsize.z/(double)global_mesh.z);

      Ncut = SpaceVector<int>(int(ceil(Rcut1/dmesh.x)), int(ceil(Rcut1/dmesh.y)), int(ceil(Rcut1/dmesh.z)));

      std::cout << "dmesh" << dmesh << "  Ncut" << Ncut << " sigma1 " << sigma1 << std::endl;

      nPREsmooth = 4;
      nPOSTsmooth = 4;
      NCYCLE = 4;

    }


    void grid_pbc_parallel(double *u, SpaceVector<int> n, SpaceVector<int> overlap);
    void copy(double *aout, double *ain, SpaceVector<int> n);
    void fillzero(double *a, SpaceVector<int> n);
    void addint(double *uf, double *res, SpaceVector<int> n);
    void restr(double *uc, double *uf, grid *coarsegrid, grid *finegrid);
    void interp(double *uf, double *uc, grid *finegrid, grid *coarsegrid);
    void relax(double *u, double *f, grid *grd, int niter);
    void resid(double *res, double *u, double *f, grid *grd);
    double L2norm(double *u, SpaceVector<int> n);

    void direct_solver(double *u, double *rhs, grid *grd);

    void MultigridPoissonSolver(double *rho, double *pot,
				SpaceVector<int> n, 
				SpaceVector<double> dmesh, 
				int nPREsmooth, int nPOSTsmooth, int NCYCLE);

    void RBGSPoissonSolver(double *rho, double *pot,
			   SpaceVector<int> n, 
			   SpaceVector<double> dmesh, 
			   int niter);

    template<class PA, class GPA>
      void chargeassignment(PA& particle,
			    const std::vector<TypeRange>& typerange,
			    const std::vector<int>& self_longset_index,
			    GPA& ghost,
			    const std::vector<TypeRange>& ghosttyperange,
			    const std::vector<int>& ghost_longset_index);


    template<class PA, class GPA>
      void backinterpolation(PA& particle,
			     const std::vector<TypeRange>& typerange,
			     const std::vector<int>& self_longset_index,
			     GPA& ghost,
			     const std::vector<TypeRange>& ghosttyperange,
			     const std::vector<int>& ghost_longset_index);


};    
#endif

