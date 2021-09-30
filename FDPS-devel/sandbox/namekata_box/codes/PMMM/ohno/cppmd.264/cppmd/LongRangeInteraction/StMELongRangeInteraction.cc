#include <mpi.h>
#include <cmath>
#include "LongRangeInteraction.h"
#include "StMELongRangeInteraction.h"

#ifndef FORTH_ORDER
#ifndef SECOND_ORDER
#define SECOND_ORDER
#endif
#endif

    template<class PA, class GPA>
    void StMELongRangeInteraction::chargeassignment(PA& particle,
			    const std::vector<TypeRange>& typerange,
			    const std::vector<int>& self_longset_index,
			    GPA& ghost,
			    const std::vector<TypeRange>& ghosttyperange,
			    const std::vector<int>& ghost_longset_index)
    {
      //      printf("StMELongRangeInteraction::chargeassignment\n");
      const SpaceVector<int> overlap(1,1,1);

      SpaceVector<int> mesh, size = local_mesh + 2*overlap;
      for(mesh.x = 0; mesh.x < size.x; ++mesh.x){
	for(mesh.y = 0; mesh.y < size.y; ++mesh.y){
	  for(mesh.z = 0; mesh.z < size.z; ++mesh.z){
	    int addr = FI3D(mesh.x,mesh.y,mesh.z,size.x,size.y,size.z);
	    rhov[addr] = rhoc[addr] = 0.0;
	  }
	}    
      }


      //First charge spreading  
      //    const int Support = (Ncut.x*2+1)*(Ncut.y*2+1)*(Ncut.z*2+1);
      //printf("Ncut1=(%d,%d,%d), Support=%d\n",Ncut.x,Ncut.y,Ncut.z,Support);
    
#ifdef DEBUG_VERB
      int np=0;
      int ng=0;
      double tv=0.0, tc=0.0, oc=0.0;
#endif
      for(size_t si=0;si<self_longset_index.size();si++){
	for(int i=typerange[si].begin;i<typerange[si].end;i++){
	  Position& particle_i_pos = getpos(particle,i);
	  double& particle_i_charge = getcharge(particle,i);
#ifdef HAVE_SPACEVECTOR_RECP
	  const SpaceVector<int> src_mesh = particle_i_pos/dmesh + overlap;//world coord
#else
	  const SpaceVector<int> src_mesh(
					  int(floor(particle_i_pos.x/dmesh.x)) + overlap.x,
					  int(floor(particle_i_pos.y/dmesh.y)) + overlap.y,
					  int(floor(particle_i_pos.z/dmesh.z)) + overlap.z );
#endif
	  const SpaceVector<int> a = src_mesh - mesh_begin - Ncut;
          const SpaceVector<int> begin(MAX(a.x,0), MAX(a.y, 0), MAX(a.z, 0));
	  const SpaceVector<int> b = src_mesh - mesh_begin + Ncut + SpaceVector<int>(1,1,1);
          //          const SpaceVector<int> b = src_mesh - mesh_begin + Ncut;
          const SpaceVector<int> c = local_mesh + overlap*2;
	  const SpaceVector<int> end(MIN(b.x, c.x), MIN(b.y, c.y), MIN(b.z, c.z));
	  SpaceVector<int> mesh;//in node coord
	  for(mesh.x = begin.x; mesh.x < end.x; ++mesh.x){
	    for(mesh.y = begin.y; mesh.y < end.y; ++mesh.y){
	      for(mesh.z = begin.z; mesh.z < end.z; ++mesh.z){
		int addr =FI3D(mesh.x,mesh.y,mesh.z,size.x,size.y,size.z); 
		const SpaceVector<double> rv(
					     particle_i_pos.x - (mesh.x + mesh_begin.x) * dmesh.x,
					     particle_i_pos.y - (mesh.y + mesh_begin.y) * dmesh.y,
					     particle_i_pos.z - (mesh.z + mesh_begin.z) * dmesh.z );
		const double Gaussv = exp(-rv.norm2()/(2.0*sigma1*sigma1))/CUBE(sqrt(2.0*M_PI)*sigma1); 
		rhov[addr] += particle_i_charge * Gaussv;
		const SpaceVector<double> rc = rv - 0.5 * dmesh;

		const double Gaussc = exp(-rc.norm2()/(2.0*sigma1*sigma1))/CUBE(sqrt(2.0*M_PI)*sigma1); 
		rhoc[addr] += particle_i_charge * Gaussc;
#ifdef DEBUG_VERB
		tv += particle_i_charge * Gaussv;
                tc += particle_i_charge * Gaussc;
                ng++;
                oc += particle_i_charge;
#endif
	      }
	    }
	  }
#ifdef DEBUG_VERB
          np++;
#endif
	}
      }
      for(size_t si=0;si<ghost_longset_index.size();si++){
	for(int i=ghosttyperange[si].begin;i<ghosttyperange[si].end;i++){
	  Position& ghost_i_pos = getpos(ghost,i);
	  double& ghost_i_charge = getcharge(ghost,i);
#ifdef HAVE_SPACEVECTOR_RECP
	  const SpaceVector<int> src_mesh = ghost_i_pos/dmesh + overlap;//world coord
#else
	  const SpaceVector<int> src_mesh(
					  int(floor(ghost_i_pos.x/dmesh.x)) + overlap.x,
					  int(floor(ghost_i_pos.y/dmesh.y)) + overlap.y,
					  int(floor(ghost_i_pos.z/dmesh.z)) + overlap.z);
#endif
	  const SpaceVector<int> a = src_mesh - mesh_begin - Ncut;
          const SpaceVector<int> begin(MAX(a.x,0), MAX(a.y, 0), MAX(a.z, 0));
	  const SpaceVector<int> b = src_mesh - mesh_begin + Ncut;
          const SpaceVector<int> c = local_mesh + overlap*2;
	  const SpaceVector<int> end(MIN(b.x, c.x), MIN(b.y, c.y), MIN(b.z, c.z));
	  
	  SpaceVector<int> mesh;//in node coord
	  for(mesh.x = begin.x; mesh.x < end.x; ++mesh.x){
	    for(mesh.y = begin.y; mesh.y < end.y; ++mesh.y){
	      for(mesh.z = begin.z; mesh.z < end.z; ++mesh.z){
		int addr = FI3D(mesh.x,mesh.y,mesh.z,size.x,size.y,size.z);
		/*
		  const SpaceVector<double> rv = ghost_i_pos - (double)(mesh + mesh_begin) * dmesh;
		*/
		const SpaceVector<double> rv(
					     ghost_i_pos.x - (mesh.x + mesh_begin.x) * dmesh.x,
					     ghost_i_pos.y - (mesh.y + mesh_begin.y) * dmesh.y,
					     ghost_i_pos.z - (mesh.z + mesh_begin.z) * dmesh.z );
		const double Gaussv = exp(-rv.norm2()/(2.0*sigma1*sigma1))/CUBE(sqrt(2.0*M_PI)*sigma1); 
		rhov[addr] += ghost_i_charge * Gaussv;

		const SpaceVector<double> rc = rv - 0.5 * dmesh;

		const double Gaussc = exp(-rc.norm2()/(2.0*sigma1*sigma1))/CUBE(sqrt(2.0*M_PI)*sigma1); 
		rhoc[addr] += ghost_i_charge * Gaussc;
#ifdef DEBUG_VERB
                tv += ghost_i_charge * Gaussv;
                tc += ghost_i_charge * Gaussc;
                ng++;
                oc += ghost_i_charge;
#endif
	      }
	    }
	  }
#ifdef DEBUG_VERB
          np++;
#endif
	}
      }

#ifdef DEBUG_VERB
      printf("self+ghost assign %d particle %d grid total rhov %e rhoc %e charge %e\n",np,ng,tv,tc,oc);
#endif

    }
#ifdef OLDPARTICLE
template
void StMELongRangeInteraction::chargeassignment(ParticleArray& particle,
                                                const std::vector<TypeRange>& typerange,
                                                const std::vector<int>& self_longset_index,
                                                ParticleArray& ghost,
                                                const std::vector<TypeRange>& ghosttyperange,
                                                const std::vector<int>& ghost_longset_index);
#else
template
void StMELongRangeInteraction::chargeassignment(CombinedParticleArray& particle,
                                                const std::vector<TypeRange>& typerange,
                                                const std::vector<int>& self_longset_index,
                                                GhostParticleArray& ghost,
                                                const std::vector<TypeRange>& ghosttyperange,
                                                const std::vector<int>& ghost_longset_index);
#endif

    template<class PA, class GPA>
    void StMELongRangeInteraction::backinterpolation(PA& particle,
			     const std::vector<TypeRange>& typerange,
			     const std::vector<int>& self_longset_index,
			     GPA& ghost,
			     const std::vector<TypeRange>& ghosttyperange,
			     const std::vector<int>& ghost_longset_index)
    {
      //      printf("StMELongRangeInteraction::backinterpolation\n");

     //mesh convolution on the particles 
      const double hhh = dmesh.x*dmesh.y*dmesh.z; /// It used 137, moved

      const SpaceVector<int> overlap(1,1,1);
      const SpaceVector<int> size = local_mesh + overlap*2; //////////// node_mesh means local_mesh?   It used 192, moved out from next scope

      //total potential
      {
	//      const SpaceVector<int> start = overlap, end = node_mesh + overlap;
	const SpaceVector<int> start = overlap;
        const SpaceVector<int> end = local_mesh + overlap;
	  
	SpaceVector<int> mesh;
	for(mesh.x = start.x; mesh.x < end.x; ++mesh.x){
	  for(mesh.y = start.y; mesh.y < end.y; ++mesh.y){
	    for(mesh.z = start.z; mesh.z < end.z; ++mesh.z){
	      int addr = FI3D(mesh.x,mesh.y,mesh.z,size.x,size.y,size.z);
	      total_potential += 0.25 * hhh 
		* (rhoc[addr] * potc[addr] + rhov[addr] * potv[addr]);
	    }
	  }
	}
	double source = total_potential, dest;
	MPI_Allreduce(&source, &dest, 1, MPI_DOUBLE_PRECISION, MPI_SUM, myworld);
	total_potential = dest;
      }
      
      //force computation
      /*
      double potv_with_cut[MAX_NNN];
      double potc_with_cut[MAX_NNN];
      */
  

      SpaceVector<int> mesh;
      for(mesh.x = 0; mesh.x < local_mesh.x; mesh.x++){
	for(mesh.y = 0; mesh.y < local_mesh.y; mesh.y++){
	  for(mesh.z = 0; mesh.z < local_mesh.z; mesh.z++){
	    int src = FI3D(mesh.x+1,mesh.y+1,mesh.z+1,
			   local_mesh.x+2,local_mesh.y+2,local_mesh.z+2);
	    int dst = FI3D(mesh.x+Ncut.x,mesh.y+Ncut.y,mesh.z+Ncut.z,
			   local_mesh.x+2+Ncut.x,local_mesh.y+2*Ncut.y,local_mesh.z+2*Ncut.z);
	    potv_with_cut[dst] = potv[src]; 
	    potc_with_cut[dst] = potc[src]; 
	  }
	}    
      }
      grid_pbc_parallel(potv_with_cut, local_mesh, Ncut);
      grid_pbc_parallel(potc_with_cut, local_mesh, Ncut);
      

      for(size_t si=0;si<self_longset_index.size();si++){
	for(int i=typerange[si].begin;i<typerange[si].end;i++){
	  Position& particle_i_pos = getpos(particle,i);
	  double& particle_i_charge = getcharge(particle,i);
#ifdef HAVE_SPACEVECTOR_RECP
	  const SpaceVector<int> src_mesh = particle_i_pos/dmesh + Ncut;//world coord
#else
	  const SpaceVector<int> src_mesh(
					  int(floor(particle_i_pos.x/dmesh.x)) + Ncut.x,
					  int(floor(particle_i_pos.y/dmesh.y)) + Ncut.y,
					  int(floor(particle_i_pos.z/dmesh.z)) + Ncut.z);
#endif
	  const SpaceVector<int> a = src_mesh - mesh_begin - Ncut;
          const SpaceVector<int> begin(MAX(a.x,0), MAX(a.y, 0), MAX(a.z, 0));
	  const SpaceVector<int> b = src_mesh - mesh_begin + Ncut + SpaceVector<int>(1,1,1);
          //	  const SpaceVector<int> b = src_mesh - mesh_begin + Ncut;
          const SpaceVector<int> c = local_mesh + Ncut*2;
          const SpaceVector<int> end(MIN(b.x, c.x), MIN(b.y, c.y), MIN(b.z, c.z));


	  SpaceVector<int> mesh;//in node coord
	  for(mesh.x = begin.x; mesh.x < end.x; ++mesh.x){
	    for(mesh.y = begin.y; mesh.y < end.y; ++mesh.y){
	      for(mesh.z = begin.z; mesh.z < end.z; ++mesh.z){
		int addr = FI3D(mesh.x,mesh.y,mesh.z,size.x,size.y,size.z);
		/*
		  const SpaceVector<double> rv = particle_i_pos - (mesh + mesh_begin) * dmesh;
		*/
		const SpaceVector<double> rv(
					     particle_i_pos.x - (mesh.x + mesh_begin.x) * dmesh.x,
					     particle_i_pos.y - (mesh.y + mesh_begin.y) * dmesh.y,
					     particle_i_pos.z - (mesh.z + mesh_begin.z) * dmesh.z );
		const double Gaussv = exp(-rv.norm2()/(2.0*sigma1*sigma1)) * potv_with_cut[addr]; 
		const SpaceVector<double> rc = rv - 0.5 * dmesh;
		const double Gaussc = exp(-rc.norm2()/(2.0*sigma1*sigma1)) * potc_with_cut[addr];

		getforce(particle,i) += rv * Gaussv * particle_i_charge * hhh/(CUBE(sqrt(2.0*M_PI)*sigma1)*SQR(sigma1))*0.5;
                getforce(particle,i) += rc * Gaussc * particle_i_charge * hhh/(CUBE(sqrt(2.0*M_PI)*sigma1)*SQR(sigma1))*0.5;
	      }
	    }
	  }
	
	}
      }
    
    }

#ifdef OLDPARTICLE
template
void StMELongRangeInteraction::backinterpolation(ParticleArray& particle,
                                                 const std::vector<TypeRange>& typerange,
                                                 const std::vector<int>& self_longset_index,
                                                 ParticleArray& ghost,
                                                 const std::vector<TypeRange>& ghosttyperange,
                                                 const std::vector<int>& ghost_longset_index);
#else
template
void StMELongRangeInteraction::backinterpolation(CombinedParticleArray& particle,
                                                 const std::vector<TypeRange>& typerange,
                                                 const std::vector<int>& self_longset_index,
                                                 GhostParticleArray& ghost,
                                                 const std::vector<TypeRange>& ghosttyperange,
                                                 const std::vector<int>& ghost_longset_index);
#endif



void StMELongRangeInteraction::relax(double *u, double *f, grid *grd, int niter)  {
  //  printf("StMELongRangeInteraction::relax\n");
    const int n1 = grd->n.x, n2 = grd->n.y, n3 = grd->n.z; 
    const double hhhinv = 1.0/grd->hhh;
    const SpaceVector<int> overlap(1,1,1);

    //    printf("call grid_pbc_parallel\n");
    grid_pbc_parallel(f, grd->n, overlap);// required for 4th-order

#ifdef DEBUG_VERB
    double sumu1=0.0,sumul=0.0;
#endif
    for(int iter = 0; iter < niter; ++iter){
      for(int color = 0; color < 2; ++color){
	grid_pbc_parallel(u, grd->n, overlap);  ///// Caution or _serial without overlap?

#pragma omp parallel for default(none) shared(u,f,grd,color,hhhinv)
	for(int i=1; i < grd->n.x+1; ++i){
	  for(int j=1; j < grd->n.y+1; ++j){
	    int ksw = ((color + i + j) % 2) + 1; // ksw={1,2} for RedBlack
	    for(int k=ksw; k < grd->n.z+1; k+=2){
#ifdef SECOND_ORDER
	      u[FI3D(i,j,k,n1,n2,n3)] = 
		( ( u[FI3D(i+1,j,k,n1,n2,n3)] + u[FI3D(i-1,j,k,n1,n2,n3)]) * grd->h1sqinv 
		  +(u[FI3D(i,j+1,k,n1,n2,n3)] + u[FI3D(i,j-1,k,n1,n2,n3)]) * grd->h2sqinv
		  +(u[FI3D(i,j,k+1,n1,n2,n3)] + u[FI3D(i,j,k-1,n1,n2,n3)]) * grd->h3sqinv
		  -f[FI3D(i,j,k,n1,n2,n3)])
		/(2.0*(grd->h1sqinv + grd->h2sqinv + grd->h3sqinv));
#endif
	    // Mehrstellen discretization (forth order scheme)
#define FORTH_ORDER
#ifdef FORTH_ORDER
              u[FI3D(i,j,k,n1,n2,n3)] =  
		((u[FI3D(i+1,j,k,n1,n2,n3)] + u[FI3D(i-1,j,k,n1,n2,n3)]) * grd->hx 
		 +(u[FI3D(i,j+1,k,n1,n2,n3)] + u[FI3D(i,j-1,k,n1,n2,n3)]) * grd->hy 
		 +(u[FI3D(i,j,k+1,n1,n2,n3)] + u[FI3D(i,j,k-1,n1,n2,n3)]) * grd->hz 
		 +(u[FI3D(i+1,j+1,k,n1,n2,n3)]+u[FI3D(i-1,j-1,k,n1,n2,n3)] 
		   +u[FI3D(i+1,j-1,k,n1,n2,n3)]+u[FI3D(i-1,j+1,k,n1,n2,n3)]) * grd->hhxy 
		 +(u[FI3D(i+1,j,k+1,n1,n2,n3)]+u[FI3D(i-1,j,k-1,n1,n2,n3)] 
		   +u[FI3D(i+1,j,k-1,n1,n2,n3)]+u[FI3D(i-1,j,k+1,n1,n2,n3)]) * grd->hhzx 
		 +(u[FI3D(i,j+1,k+1,n1,n2,n3)]+u[FI3D(i,j-1,k-1,n1,n2,n3)] 
		   +u[FI3D(i,j-1,k+1,n1,n2,n3)]+u[FI3D(i,j+1,k-1,n1,n2,n3)]) * grd->hhyz 
		 -0.5*f[FI3D(i,j,k,n1,n2,n3)] 
		 -(f[FI3D(i+1,j,k,n1,n2,n3)]+f[FI3D(i-1,j,k,n1,n2,n3)]
		   +f[FI3D(i,j+1,k,n1,n2,n3)]+f[FI3D(i,j-1,k,n1,n2,n3)]
		   +f[FI3D(i,j,k+1,n1,n2,n3)]+f[FI3D(i,j,k-1,n1,n2,n3)])
		 /12.0) * hhhinv;
#endif
#ifdef DEBUG_VERB
              if(iter==0)sumu1+=u[FI3D(i,j,k,n1,n2,n3)];
              if(iter==niter-1)sumul+=u[FI3D(i,j,k,n1,n2,n3)];
#endif
	    }
	  }
	}
      }//color
    }//iter 
#ifdef DEBUG_VERB
    printf("sumu %e -> %e\n",sumu1,sumul);
#endif
  }

  void StMELongRangeInteraction::resid(double *res, double *u, double *f, grid *grd)  {
    //  printf("StMELongRangeInteraction::resid\n");
    const int n1 = grd->n.x, n2 = grd->n.y, n3 = grd->n.z; 
    const SpaceVector<int> overlap(1,1,1);
    grid_pbc_parallel(u, grd->n, overlap);

#pragma omp parallel for default(none) shared(res,u,f,grd,n1,n2,n3)
    for (int i=1;i<grd->n.x+1;i++){
      for (int j=1;j<grd->n.y+1;j++){
	for (int k=1;k<grd->n.z+1;k++){  
#ifdef SECOND_ORDER
	  res[FI3D(i,j,k,n1,n2,n3)] 
	    = -(  ( u[FI3D(i+1,j,k,n1,n2,n3)] + u[FI3D(i-1,j,k,n1,n2,n3)])*grd->h1sqinv
		  +(u[FI3D(i,j+1,k,n1,n2,n3)] + u[FI3D(i,j-1,k,n1,n2,n3)])*grd->h2sqinv
		  +(u[FI3D(i,j,k+1,n1,n2,n3)] + u[FI3D(i,j,k-1,n1,n2,n3)])*grd->h3sqinv
		  -u[FI3D(i,j,k,n1,n2,n3)]
		  *(2.0*(grd->h1sqinv + grd->h2sqinv + grd->h3sqinv)) ) 
	    + f[FI3D(i,j,k,n1,n2,n3)];
#endif

#ifdef FORTH_ORDER
	  res[FI3D(i,j,k,n1,n2,n3)] =   
	    -((u[FI3D(i+1,j,k,n1,n2,n3)] + u[FI3D(i-1,j,k,n1,n2,n3)]) * grd->hx 
	      +(u[FI3D(i,j+1,k,n1,n2,n3)] + u[FI3D(i,j-1,k,n1,n2,n3)]) * grd->hy 
	    +(u[FI3D(i,j,k+1,n1,n2,n3)] + u[FI3D(i,j,k-1,n1,n2,n3)]) * grd->hz 
	    +(u[FI3D(i+1,j+1,k,n1,n2,n3)]+u[FI3D(i-1,j-1,k,n1,n2,n3)]
	      +u[FI3D(i+1,j-1,k,n1,n2,n3)]+u[FI3D(i-1,j+1,k,n1,n2,n3)]) * grd->hhxy 
	    +(u[FI3D(i+1,j,k+1,n1,n2,n3)]+u[FI3D(i-1,j,k-1,n1,n2,n3)]
	      +u[FI3D(i+1,j,k-1,n1,n2,n3)]+u[FI3D(i-1,j,k+1,n1,n2,n3)]) * grd->hhzx 
	    +(u[FI3D(i,j+1,k+1,n1,n2,n3)]+u[FI3D(i,j-1,k-1,n1,n2,n3)]
	      +u[FI3D(i,j-1,k+1,n1,n2,n3)]+u[FI3D(i,j+1,k-1,n1,n2,n3)]) * grd->hhyz
	    -u[FI3D(i,j,k,n1,n2,n3)]*grd->hhh) 
	  +0.5*f[FI3D(i,j,k,n1,n2,n3)] 
	    +(f[FI3D(i+1,j,k,n1,n2,n3)]+f[FI3D(i-1,j,k,n1,n2,n3)]
	      +f[FI3D(i,j+1,k,n1,n2,n3)]+f[FI3D(i,j-1,k,n1,n2,n3)]
	      +f[FI3D(i,j,k+1,n1,n2,n3)]+f[FI3D(i,j,k-1,n1,n2,n3)])/12.0;
#endif
	}
      }
    }
  }

  double StMELongRangeInteraction::L2norm(double *u, SpaceVector<int> n){

    double L2 = 0.0;
    for(int x = 1; x < n.x+1; ++x){
      for(int y = 1; y < n.y+1; ++y){
	for(int z = 1; z < n.z+1; ++z){
	  int addr = FI3D(x,y,z,n.x+2,n.y+2,n.z+2);
	  L2 += SQR(u[addr]);
	}
      }
    }
    double src = L2, dest;
    MPI_Allreduce(&src, &dest, 1, MPI_DOUBLE_PRECISION, MPI_SUM, myworld);
    L2 = dest;
    return (L2);
  }

  void StMELongRangeInteraction::RBGSPoissonSolver(double *rho, double *pot,
			 SpaceVector<int> n, 
			 SpaceVector<double> dmesh, 
			 int niter)  {
    //    printf("StMELongRangeInteraction::RBGSPoissonSolver\n");

    /*
    grid grd;
    
    grd.u = new double [MAX_NNN];
    grd.f = new  double [MAX_NNN];
    grd.res = new double [MAX_NNN];
    */

    grd.n = n;

    grd.h1sqinv = 1.0/SQR(dmesh.x);
    grd.h2sqinv = 1.0/SQR(dmesh.y);
    grd.h3sqinv = 1.0/SQR(dmesh.z);

    grd.hhh = 16.0*(grd.h1sqinv + grd.h2sqinv + grd.h3sqinv)/12.0;
    grd.hhxy = (grd.h1sqinv + grd.h2sqinv)/12.0; 
    grd.hhzx = (grd.h1sqinv + grd.h3sqinv)/12.0; 
    grd.hhyz = (grd.h3sqinv + grd.h2sqinv)/12.0;
    grd.hx = (4.0 * grd.h1sqinv - grd.h2sqinv - grd.h3sqinv)/6.0;
    grd.hy = (4.0 * grd.h2sqinv - grd.h3sqinv - grd.h1sqinv)/6.0;
    grd.hz = (4.0 * grd.h3sqinv - grd.h1sqinv - grd.h2sqinv)/6.0; 

    for(int i = 1; i < grd.n.x+1; ++i){
      for(int j = 1; j < grd.n.y+1; ++j){
	for(int k = 1; k < grd.n.z+1; ++k){
	  int addr = FI3D(i,j,k,grd.n.x+2,grd.n.y+2,grd.n.z+2);
	  grd.f[addr] = -4.0 * M_PI * rho[addr];
	  grd.u[addr] = 0.0;/// modify if initial guess applied
	}
      }
    }
    relax(grd.u, grd.f, &grd, niter);
    resid(grd.res, grd.u ,grd.f, &grd);

    double L2res = L2norm(grd.res, grd.n);
    if(unit_identifier==0){
      printf("L2res=%e\n",sqrt(L2res*grd.dmesh.x*grd.dmesh.y*grd.dmesh.z));
    }

    //    double offset = 0.0; 
    //for(int i = 0; i < volume(n); ++i) offset += pot[i];
    //for(int i = 0; i < volume(n); ++i) pot[i] -= offset/volume(n);


    for(int x = 1; x < grd.n.x+1; ++x){
      for(int y = 1; y < grd.n.y+1; ++y){
	for(int z = 1; z < grd.n.z+1; ++z){
	  int addr = FI3D(x,y,z,grd.n.x+2,grd.n.y+2,grd.n.z+2);
	  pot[addr] = grd.u[addr]; 
	}
      }
    }

    /*
    delete [] grd.u;
    delete [] grd.f;
    delete [] grd.res;
    */
  }




void StMELongRangeInteraction::grid_pbc_parallel(double *u, SpaceVector<int> n, SpaceVector<int> overlap){

  //  printf("StMELongRangeInteraction::grid_pbc_parallel\n");
  //  fflush(stdout);

    const int n1 = n.x + 2*overlap.x, n2 = n.y + 2*overlap.y, n3 = n.z + 2*overlap.z;
    //    double sendbuff[MAX_NNN], recvbuff[MAX_NNN];

    int dest, source, tag = 1, count;
    MPI_Status mpi_status;

    int down = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0,0,-1),myrank);
    int up = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0,0, 1),myrank); 

    //z+
    count = 0;
    for(int i=0;i<n1;i++) {
      for(int j=0;j<n2;j++) { 
	for(int k=local_mesh.z;k<local_mesh.z + overlap.z;k++){
	  sendbuff[count]   = u[FI3D(i,j,k,n1,n2,n3)];
	  count++;
	}
      }
    }
    dest = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0,0, 1),myrank); 
    source = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0,0,-1),myrank);

    //    printf("MPI_Sendrecv to %d from %d\n",dest,source);
    MPI_Sendrecv( &sendbuff, count, MPI_DOUBLE_PRECISION, dest, tag,
		  &recvbuff, count, MPI_DOUBLE_PRECISION, source, tag,
		  myworld, &mpi_status );
    count = 0;
    for (int i=0;i<n1;i++) {
      for (int j=0;j<n2;j++) { 
	for(int k=0;k<overlap.z;k++){
	  u[FI3D(i,j,k,n1,n2,n3)] = recvbuff[count];
	  count++;
	}
      }
    }
    //z-
    count = 0;
    for(int i=0;i<n1;i++) {
      for(int j=0;j<n2;j++) { 
	for(int k=overlap.z;k<overlap.z*2;k++){
	  sendbuff[count]   = u[FI3D(i,j,k,n1,n2,n3)];
	  count++;
	}
      }
    }
    dest = down;//node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0,0,-1),myrank);
    source = up;//node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0,0, 1),myrank); 

    MPI_Sendrecv( &sendbuff, count, MPI_DOUBLE_PRECISION, dest, tag,
		  &recvbuff, count, MPI_DOUBLE_PRECISION, source, tag,
		  myworld, &mpi_status );
    count =0;
    for (int i=0;i<n1;i++) {
      for (int j=0;j<n2;j++) { 
	for(int k=local_mesh.z + overlap.z; k<local_mesh.z + overlap.z*2; k++){
	  u[FI3D(i,j,k,n1,n2,n3)] = recvbuff[count];
	  count++;
	}
      }
    }
    //y+
    count = 0;
    for(int i=0;i<n1;i++) {
      for(int j=local_mesh.y; j<local_mesh.y + overlap.y;j++){
	for(int k=0;k<n3;k++) { 
	  sendbuff[count]   = u[FI3D(i,j,k,n1,n2,n3)];
	  count++;
	}
      }
    }


    dest = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0, 1,0),myrank);
    source = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0,-1,0),myrank);

    MPI_Sendrecv( &sendbuff, count, MPI_DOUBLE_PRECISION, dest, tag,
		  &recvbuff, count, MPI_DOUBLE_PRECISION, source, tag,
		  myworld, &mpi_status );
    count =0;
    for (int i=0;i<n1;i++) {
      for(int j=0;j<overlap.y;j++){
	for (int k=0;k<n3;k++) { 
	  u[FI3D(i,j,k,n1,n2,n3)] = recvbuff[count];
	  count++;
	}
      }
    }
    //y-
    count = 0;
    for(int i=0;i<n1;i++) {
      for(int j=overlap.y;j<overlap.y*2;j++){
	for(int k=0;k<n3;k++) { 
	  sendbuff[count]   = u[FI3D(i,j,k,n1,n2,n3)];
	  count++;
	}
      }
    }

    dest = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0,-1,0),myrank);
    source = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 0, 1,0),myrank);
    MPI_Sendrecv( &sendbuff, count, MPI_DOUBLE_PRECISION, dest, tag,
		  &recvbuff, count, MPI_DOUBLE_PRECISION, source, tag,
		  myworld, &mpi_status );
    count =0;
    for (int i=0;i<n1;i++) {
      for(int j=local_mesh.y + overlap.y; j<local_mesh.y + overlap.y*2; j++){
	for (int k=0;k<n3;k++) { 
	  u[FI3D(i,j,k,n1,n2,n3)] = recvbuff[count];
	  count++;
	}
      }
    }
    //x+
    count = 0;
    for(int i=local_mesh.x;i<local_mesh.x + overlap.x;i++){
      for(int j=0;j<n2;j++) {
	for(int k=0;k<n3;k++) { 
	  sendbuff[count]   = u[FI3D(i,j,k,n1,n2,n3)];
	  count++;
	}
      }
    }

    dest  = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 1, 0,0),myrank);
    source  = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>(-1, 0,0),myrank);

    MPI_Sendrecv( &sendbuff, count, MPI_DOUBLE_PRECISION, dest, tag,
		  &recvbuff, count, MPI_DOUBLE_PRECISION, source, tag,
		  myworld, &mpi_status );
    count =0;
    for(int i=0;i<overlap.x;i++){
      for (int j=0;j<n2;j++) {
	for (int k=0;k<n3;k++) { 
	  u[FI3D(i,j,k,n1,n2,n3)] = recvbuff[count];
	  count++;
	}
      }
    }
    //x-
    count = 0;
    for(int i=overlap.x;i<overlap.x*2;i++){
      for(int j=0;j<n2;j++) { 
	  for(int k=0;k<n3;k++) {
	    sendbuff[count]   = u[FI3D(i,j,k,n1,n2,n3)];
	  count++;
	}
      }
    }

    dest  = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>(-1, 0,0),myrank);
    source  = node_geometry.getNodeIDofRelativePosition(SpaceVector<int>( 1, 0,0),myrank);

    MPI_Sendrecv( &sendbuff, count, MPI_DOUBLE_PRECISION, dest, tag,
		  &recvbuff, count, MPI_DOUBLE_PRECISION, source, tag,
		  myworld, &mpi_status );
    count =0;
    for(int i=local_mesh.x + overlap.x; i<local_mesh.x + overlap.x*2; i++){
      for (int j=0;j<n2;j++) { 
	for (int k=0;k<n3;k++) {
	  u[FI3D(i,j,k,n1,n2,n3)] = recvbuff[count];
	  count++;
	}
      }
    }



  }


void StMELongRangeInteraction::copy(double *aout, double *ain, SpaceVector<int> n){
  for(int i=0;i<n.x+2; i++){
    for(int j=0;j<n.y+2; j++){
      for(int k=0;k<n.z+2; k++){
	int addr = FI3D(i,j,k,n.x+2,n.y+2,n.z+2);
	aout[addr] = ain[addr];
      }
    }
  }
}
void StMELongRangeInteraction::fillzero(double *a, SpaceVector<int> n){
  for(int i=0;i<n.x+2; i++){
    for(int j=0;j<n.y+2; j++){
      for(int k=0;k<n.z+2; k++){
	int addr = FI3D(i,j,k,n.x+2,n.y+2,n.z+2);
	a[addr] = 0.0;
      }
    }
  }
}
void StMELongRangeInteraction::addint(double *uf, double *res, SpaceVector<int> n){
  for(int i=0;i<n.x+2; i++){
    for(int j=0;j<n.y+2; j++){
      for(int k=0;k<n.z+2; k++){
	int addr = FI3D(i,j,k,n.x+2,n.y+2,n.z+2);
	uf[addr] += res[addr];
      }
    }
  }
}

void StMELongRangeInteraction::MultigridPoissonSolver(double *rho, double *pot,
						      SpaceVector<int> n, 
						      SpaceVector<double> dmesh, 
						      int nPREsmooth, int nPOSTsmooth, int NCYCLE)
{

  grid grd[MAX_LEVEL];

  // set finest level
  int finestLevel = 0;
  grd[finestLevel].n.x = n.x;
  grd[finestLevel].n.y = n.y;
  grd[finestLevel].n.z = n.z;

  // set local number of grid for serial mode
  const int top = 4;
  int current = finestLevel, n1c, n2c, n3c;
  while(grd[current].n.x > top || grd[current].n.y > top || grd[current].n.z > top){
      if(grd[current].n.x % 2 != 0) break;
      if(grd[current].n.y % 2 != 0) break;
      if(grd[current].n.z % 2 != 0) break;
      n1c = grd[current].n.x / 2;
      n2c = grd[current].n.y / 2;
      n3c = grd[current].n.z / 2;
      current++;
      grd[current].n.x = n1c;
      grd[current].n.y = n2c;
      grd[current].n.z = n3c;
    }
  // set coarsesst level
  int coarsestLevel = current;

  const SpaceVector<double> local_side(n.x*dmesh.x, n.y*dmesh.y, n.z*dmesh.z);
  for(current = coarsestLevel; current < finestLevel; current++){

    int n = (grd[current].n.x+2)
      *(grd[current].n.z+2)
      *(grd[current].n.y+2);  

    grd[current].u = new double [n];  
    grd[current].f = new double [n];  
    grd[current].res = new double [n];  

    grd[current].h1sqinv = SQR(grd[current].n.x/local_side.x);
    grd[current].h2sqinv = SQR(grd[current].n.y/local_side.y);
    grd[current].h3sqinv = SQR(grd[current].n.z/local_side.z);

    grd[current].hhh = 16.0*(grd[current].h1sqinv + grd[current].h2sqinv + grd[current].h3sqinv)/12.0;
    grd[current].hhxy = (grd[current].h1sqinv + grd[current].h2sqinv)/12.0; 
    grd[current].hhzx = (grd[current].h1sqinv + grd[current].h3sqinv)/12.0; 
    grd[current].hhyz = (grd[current].h3sqinv + grd[current].h2sqinv)/12.0;
    grd[current].hx = (4.0 * grd[current].h1sqinv - grd[current].h2sqinv - grd[current].h3sqinv)/6.0;
    grd[current].hy = (4.0 * grd[current].h2sqinv - grd[current].h3sqinv - grd[current].h1sqinv)/6.0;
    grd[current].hz = (4.0 * grd[current].h3sqinv - grd[current].h1sqinv - grd[current].h2sqinv)/6.0; 
  }

  copy(rho, grd[finestLevel].f, grd[finestLevel].n);

  // V-cycle loop
  for (int i = 0; i < NCYCLE; i++){
      // Downward stroke of V
    int nLevel = finestLevel;
    for (int j = nLevel; j < coarsestLevel; j++){
      relax(grd[j].u, grd[j].f, &grd[j], nPREsmooth); // Pre-smoothing
      resid(grd[j].res, grd[j].u, grd[j].f, &grd[j]);
      int coarse = j + 1;
      restr(grd[coarse].f, grd[j].res, &grd[coarse], &grd[j]);
      fillzero(grd[coarse].u, grd[coarse].n); // 0 for initial guess in next relax
    }

    direct_solver(grd[coarsestLevel].u, grd[coarsestLevel].f, &grd[coarsestLevel]); 
      // Upward stroke of V
    for (int j = coarsestLevel - 1; j >= nLevel; j--){
      int coarse = j + 1;
      interp(grd[j].res, grd[coarse].u, &grd[j], &grd[coarse]);
      addint(grd[j].u, grd[j].res, grd[j].n);
      relax( grd[j].u, grd[j].f, &grd[j], nPOSTsmooth); // Post-smoothing
    }
  } // end of V-cycle loop

  copy(grd[finestLevel].u, pot, grd[finestLevel].n);

  for(current = coarsestLevel; current < finestLevel; current++){
    delete [] grd[current].u;  
    delete [] grd[current].f;  
    delete [] grd[current].res;  
  }


}


/*-----------------------------------------------------------------------
 * restriction. nc is the coase-gird dimension. The fine-grid
 * solution is input in uf[1..2*nc-1][1..2*nc-1][1..2*nc-1], the coase-grid
 * solution in uc[0..nc-1][0..nc-1][0..nc-1]
 */
void StMELongRangeInteraction::restr(double *uc, double *uf, grid *coarsegrid, grid *finegrid){
  const int n1c = coarsegrid->n.x, n2c = coarsegrid->n.y, n3c = coarsegrid->n.z;
  const int n1f = finegrid->n.x, n2f = finegrid->n.y, n3f = finegrid->n.z;  

  int iFine = 1; 
#pragma omp parallel for default(none) private(iFine) shared(uc,uf,n1c,n2c,n3c,n1f,n2f,n3f)
  for (int iCoarse=1; iCoarse <= n1c; iCoarse++)  {
    iFine = iCoarse * 2 - 1;
    for (int jFine=1,jCoarse=1; jCoarse <= n2c; jCoarse++,jFine+=2)    {
      for (int kFine=1,kCoarse=1; kCoarse <= n3c; kCoarse++,kFine+=2)      { 
        uc[FI3D(iCoarse,jCoarse,kCoarse,n1c,n2c,n3c)] = // 8-points
          ( uf[FI3D(iFine,jFine  ,kFine,n1f,n2f,n3f)]   
	       + uf[FI3D(iFine+1,jFine  ,kFine,  n1f,n2f,n3f)]
	       + uf[FI3D(iFine,  jFine  ,kFine+1,n1f,n2f,n3f)] 
	       + uf[FI3D(iFine+1,jFine  ,kFine+1,n1f,n2f,n3f)]
	       + uf[FI3D(iFine,  jFine+1,kFine,  n1f,n2f,n3f)]   
	       + uf[FI3D(iFine+1,jFine+1,kFine,  n1f,n2f,n3f)]
	       + uf[FI3D(iFine,  jFine+1,kFine+1,n1f,n2f,n3f)] 
	       + uf[FI3D(iFine+1,jFine+1,kFine+1,n1f,n2f,n3f)] 
	       ) / 8.0; 
	    }
      }
    }
    
    
  }

/*------------------------------------------------------------------------
 * Coase-to-fine prolongation by bilinear interporation. nf is the 
 * fine-gird dimension. The coase-grid solution is input as 
 * nc[0..nc-1][0..nc-1][0..nc-1], where nc=nf/2+1. the fine-gird solution is 
 * returned in uf[0..nf-1][0..nf-1][0..nf-1]
 */
void StMELongRangeInteraction::interp(double *uf, double *uc, grid *finegrid, grid *coarsegrid){
  const int n1c = coarsegrid->n.x, n2c = coarsegrid->n.y, n3c = coarsegrid->n.z;
  const int n1f = finegrid->n.x, n2f = finegrid->n.y, n3f = finegrid->n.z;  
  const double w1 = 1.0/64.0, w3 = 3.0/64.0, w9 = 9.0/64.0, w27 = 27.0/64.0;

  const SpaceVector<int> overlap(1,1,1);
  grid_pbc_parallel(uf, finegrid->n, overlap);
  grid_pbc_parallel(uc, coarsegrid->n, overlap);


  int iFine = 0;
#pragma omp parallel for default(none) private(iFine) shared(uf,uc,n1c,n2c,n3c,n1f,n2f,n3f)
  for (int iCoarse=0; iCoarse <= n1c; iCoarse++)  {
    iFine = iCoarse * 2;
    for (int jCoarse=0,jFine=0; jCoarse <= n2c; jCoarse++,jFine+=2){
      for (int kCoarse=0,kFine=0; kCoarse <= n3c; kCoarse++,kFine+=2){
        double mmm = uc[FI3D(iCoarse  ,jCoarse  ,kCoarse  ,n1c,n2c,n3c)]; 
        double pmm = uc[FI3D(iCoarse+1,jCoarse  ,kCoarse  ,n1c,n2c,n3c)];
        double mpm = uc[FI3D(iCoarse  ,jCoarse+1,kCoarse  ,n1c,n2c,n3c)]; 
        double ppm = uc[FI3D(iCoarse+1,jCoarse+1,kCoarse  ,n1c,n2c,n3c)];
        double mmp = uc[FI3D(iCoarse  ,jCoarse  ,kCoarse+1,n1c,n2c,n3c)]; 
        double pmp = uc[FI3D(iCoarse+1,jCoarse  ,kCoarse+1,n1c,n2c,n3c)];
        double mpp = uc[FI3D(iCoarse  ,jCoarse+1,kCoarse+1,n1c,n2c,n3c)]; 
        double ppp = uc[FI3D(iCoarse+1,jCoarse+1,kCoarse+1,n1c,n2c,n3c)];

        uf[FI3D(iFine,jFine,kFine,n1f,n2f,n3f)] /* mmm */
          = w27 * mmm + w9 * pmm + w9 * mmp + w3 * pmp + w9 * mpm + w3 * ppm + w3 * mpp + w1 * ppp; 

        uf[FI3D(iFine+1,jFine,kFine,n1f,n2f,n3f)] /* pmm */
          = w9 * mmm + w27 * pmm + w3 * mmp + w9 * pmp + w3 * mpm + w9 * ppm + w1 * mpp + w3 * ppp; 

        uf[FI3D(iFine,jFine+1,kFine,n1f,n2f,n3f)] /* mpm */
          = w9 * mmm + w3 * pmm + w3 * mmp + w1 * pmp + w27 * mpm + w9 * ppm + w9 * mpp + w3 * ppp; 

        uf[FI3D(iFine+1,jFine+1,kFine,n1f,n2f,n3f)] /* ppm */
          = w3 * mmm + w9 * pmm + w1 * mmp + w3 * pmp + w9 * mpm + w27 * ppm + w3 * mpp + w9 * ppp; 

        uf[FI3D(iFine,jFine,kFine+1,n1f,n2f,n3f)] /* mmp */
          = w9 * mmm + w3 * pmm + w27 * mmp + w9 * pmp + w3 * mpm + w1 * ppm + w9 * mpp + w3 * ppp; 

        uf[FI3D(iFine+1,jFine,kFine+1,n1f,n2f,n3f)] /* pmp */
          = w3 * mmm + w9 * pmm + w9 * mmp + w27 * pmp + w1 * mpm + w3 * ppm + w3 * mpp + w9 * ppp; 

        uf[FI3D(iFine,jFine+1,kFine+1,n1f,n2f,n3f)] /* mpp */
          = w3 * mmm + w1 * pmm + w9 * mmp + w3 * pmp + w9 * mpm + w3 * ppm + w27 * mpp + w9 * ppp; 

        uf[FI3D(iFine+1,jFine+1,kFine+1,n1f,n2f,n3f)] /* ppp */
          = w1 * mmm + w3 * pmm + w3 * mmp + w9 * pmp + w3 * mpm + w9 * ppm + w9 * mpp + w27 * ppp; 
      }
    }
  }

}

/*-----------------------------------------------------------------
 * Solution of the model problem on the coasest grid.
 * The rhs is input in rhs[0..2][0..2][0..2] and the solution is 
 * returned in u[0..2][0..2][0..2].
 */
void StMELongRangeInteraction::direct_solver(double *u, double *rhs, grid *grd){
  int iter = global_mesh.x + global_mesh.y + global_mesh.z;
  relax(u, rhs, grd, iter);
}
