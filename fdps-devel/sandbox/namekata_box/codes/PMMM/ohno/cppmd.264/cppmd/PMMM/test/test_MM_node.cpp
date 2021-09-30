#include <mpi.h>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <sys/time.h>
#include <fftw3.h>
#include "vector3.h"

#include "MM_Cell.h"
#include "PP.h"
#include "MPIParallelPMMM.h"


#include "ewald.h"


static void print_err(
		const std::vector<double> &err,
		const char * name,
		const int icut,
		const int p)
{
	static char fname[256];
	sprintf(fname, "%s.c%dp%d.dat", name, icut, p);
	FILE *fp = fopen(fname, "w");
	assert(fp);
	const int len = err.size();
	for(int i=0; i<len; i++){
		fprintf(fp, "%e %e\n", double(i)/len, err[i]);
	}
	fclose(fp);
}


int
main(int argc, char **argv)
{
	enum{
		NP   = 128,
		NC   = 8,
		NC3  = NC*NC*NC,
		PFMM = 5,
		ICUT = 2,
	};

	MPI::Init(argc, argv);

  int node_id = 0;
  int num_node = 1;
  node_id = MPI::COMM_WORLD.Get_rank();
  num_node = MPI::COMM_WORLD.Get_size();

  if(num_node<2){
    printf("MPI_Comm_size smaller than 2\n");
    exit(EXIT_FAILURE);
  }

  int num_pp;
  int num_mm;

  num_pp = num_node/2;
  num_mm = num_node-num_pp;

  int m_decomp[(PFMM+1)*(PFMM+1)];
  mm_decompose<PFMM>(m_decomp,num_mm,num_pp);

  const double alpha = 2.4;

  MPI::Intracomm  pmmm_comm;


  int nump = NP;

  if(node_id<num_pp){
    pmmm_comm = MPI::COMM_WORLD.Split(PMMM_PP, node_id);
    
    int pp_id;
    pp_id = node_id;
    	typedef Cell_FMM<PFMM>	Cell_t;

	static Cell_t   cell[NC][NC][NC];
	const double clen = 1.0 / NC;
	for(int k=0; k<NC; k++){
		for(int j=0; j<NC; j++){
			for(int i=0; i<NC; i++){
				cell[k][j][i].set(ivec3(i,j,k), clen);
			}
		}
	}

	CellRange cell_range;
	{
	  int min[3] = {0,0,0};
	  int max[3] = {NC,NC,NC};
	  cell_range.set_min_max(min,max);
	}
	int targets[(PFMM+1)*(PFMM+1)];
	{
	  int m;
	  for(m=0;m<(PFMM+1)*(PFMM+1);m++){
	    targets[m] = m_decomp[m];
	  }
	}
	{
	  int m;
	  for(m=0;m<(PFMM+1)*(PFMM+1);m++){
	    printf(" %d",m_decomp[m]);
	  }
	  printf("\n");
	}
	PP_MM_Comm<PFMM,NC,NC,NC> pp_mm_comm;
	pp_mm_comm.set_targets(targets);
	pp_mm_comm.set_cell_range(cell_range);
	
	std::vector<Particle> pa;
	Particle *ptcl;
#ifdef EWALD_FILE
	{
		FILE *fp = fopen("qpos.dat", "r");
		assert(fp);
		fscanf(fp, "%d", &nump);
		pa.resize(nump);
		for(int i=0; i<nump; i++){
			fscanf(fp, "%lf %lf %lf %lf",
					&pa[i].mass,
					&pa[i].pos.x,
					&pa[i].pos.y,
					&pa[i].pos.z);
			pa[i].pos.x -= floor(pa[i].pos.x);
			pa[i].pos.y -= floor(pa[i].pos.y);
			pa[i].pos.z -= floor(pa[i].pos.z);
		}
		fclose(fp);
		ptcl = &(pa[0]);
	}
#else
	pa.resize(nump);
	ptcl = &(pa[0]);
	Particle::gen_rand_dist(nump, ptcl);
#endif
	double msum = 0.0;
	for(int i=0; i<nump; i++){
		msum += ptcl[i].mass;
	}
	printf("num particle %d\n",nump);
        printf("msum %e\n",msum);
	for(int i=0; i<nump; i++){
		const ivec3 idx = cell_nearest(ptcl[i].pos, clen);
		assert(0 <= idx.x && idx.x < NC);
		assert(0 <= idx.y && idx.y < NC);
		assert(0 <= idx.z && idx.z < NC);
		cell[idx.z][idx.y][idx.x].plist.push_back(&ptcl[i]);
	}

	for(int k=0; k<NC; k++) for(int j=0; j<NC; j++) for(int i=0; i<NC; i++){
		cell[k][j][i].sanity_check();
	}

	// Gen Green at MM node

	puts("Eval PM");
	Cell_t *cell1d = cell[0][0];
	for(int i=0; i<NC3; i++){
		cell1d[i].do_P2M();
	}

	// send M to MM node
	//	send_M<PFMM,NC,NC,NC>(cell, num_node-1);
	pp_mm_comm.send_M(cell);

	puts("Eval PP");
	PP_interact_PBC<PFMM, ICUT, NC, NC, NC>  (cell);

	// Dipole correction
	dvec3 dipole(0.0);
	double quad0 = 0.0;
	for(int i=0; i<NC3; i++){
		dipole.x += cell1d[i].mm.buf[3];
		dipole.y += cell1d[i].mm.buf[1];
		dipole.z += cell1d[i].mm.buf[2];
		quad0 += cell1d[i].dispersion();
	}
	const double pi = 4.0 * atan(1.0);
	dipole *= (4./3.) * pi;
	printf("quad : %e\n", quad0);

	// reciev L from MM node
	//recv_L<PFMM,NC,NC,NC>(cell, num_node-1);
	pp_mm_comm.recv_L(cell);


	for(int i=0; i<NC3; i++){
		cell1d[i].le.buf[3] += 2.0 * dipole.x;
		cell1d[i].le.buf[1] -= 2.0 * dipole.y;
		cell1d[i].le.buf[2] += 1.0 * dipole.z;
		cell1d[i].le.buf[0] += ((2./3.) * pi) * quad0;
		// self energy correction
		cell1d[i].le.buf[0] -= 
			alpha * (2.0/sqrt(pi)) * cell1d[i].mm.buf[0];
	}

		for(int i=0; i<NC3; i++){
		cell1d[i].do_L2P();
#ifdef NON_CHARGE_NEUTRAL
		cell1d[i].do_L2P_corr(msum, alpha);
#endif
	}

	dvec3 fpp(0.0), fpm(0.0);
	for(int i=0; i<nump; i++){
		fpp += ptcl[i].mass * ptcl[i].acc_direct;
		fpm += ptcl[i].mass * ptcl[i].acc_app;
	}
	printf("PP ftot : (%e, %e, %e)\n", fpp.x, fpp.y, fpp.z);
	printf("PM ftot : (%e, %e, %e)\n", fpm.x, fpm.y, fpm.z);

        double epp=0.0, epm=0.0;
        for(int i=0; i<nump; i++){
                epp += 0.5*ptcl[i].mass*ptcl[i].phi_direct;
                epm += 0.5*ptcl[i].mass*ptcl[i].phi_app;
        }
        printf("Epp %24.16e  Epm %24.16e   sum %24.16e\n",epp,epm,epp+epm);

	for(int i=0; i<nump; i++){
		ptcl[i].move_accp();
	}

#ifndef PMMM_ONLY
	puts("Ewald sum");
	// Direct Ewald
#ifdef EWALD_FILE
	{
		FILE *fp = fopen("ewald.dat", "r");
		assert(fp);
		int n;
		fscanf(fp, "%d", &n);
		assert(nump == n);
		for(int i=0; i<nump; i++){
			fscanf(fp, "%lf %lf %lf %lf",
					&ptcl[i].phi_direct,
					&ptcl[i].acc_direct.x,
					&ptcl[i].acc_direct.y,
					&ptcl[i].acc_direct.z);
		}
		fclose(fp);
	}
#else
	const double alpha_ewald = 2.4;
	eval_k_space<5>(nump, alpha_ewald, ptcl);
	eval_r_space<3>(nump, alpha_ewald, msum, ptcl);
#endif	
	double en_app=0.0, en_dir=0.0;
	for(int i=0; i<nump; i++){
		en_app += 0.5 * ptcl[i].mass * ptcl[i].phi_app;
		en_dir += 0.5 * ptcl[i].mass * ptcl[i].phi_direct;
	}
	printf("energy : %24.16e, %24.16e\n", en_app, en_dir);

	std::vector<double> err(nump);
	for(int i=0; i<nump; i++) err[i] = ptcl[i].adiff_rel();
	std::sort(err.begin(), err.end());
#ifdef EWALD_FILE
	print_err(err, "adiffr.file", ICUT, PFMM);
#else
	print_err(err, "adiffr", ICUT, PFMM);
#endif

	for(int i=0; i<nump; i++) err[i] = ptcl[i].pdiff_rel();
	std::sort(err.begin(), err.end());
#ifdef EWALD_FILE
	print_err(err, "pdiffr.file", ICUT, PFMM);
#else
	print_err(err, "pdiffr", ICUT, PFMM);
#endif
#endif
  }else{
    pmmm_comm = MPI::COMM_WORLD.Split(PMMM_MM, node_id);

    int mm_id;
    mm_id = node_id-num_pp;
    std::vector<int> mm_list;
    {
      mm_list.clear();
      int m;
      for(m=0;m<(PFMM+1)*(PFMM+1);m++){
	if(m_decomp[m]==node_id){
	  mm_list.push_back(m);
	}
      }
    }
    {
      int m;
      for(m=0;m<mm_list.size();m++){
	printf(" %d",mm_list[m]);
      }
      printf("\n");
    }
    int mm_targets[(PFMM+1)*(PFMM+1)];
    {
      int m;
      for(m=0;m<(PFMM+1)*(PFMM+1);m++){
	mm_targets[m] = m_decomp[m]-num_pp;
      }
    }

    std::vector<int> targets;
    targets.resize(1);
    targets[0] = 0;
    std::vector<CellRange> cell_range;
    cell_range.resize(1);
    {
      int min[3] = {0,0,0};
      int max[3] = {NC,NC,NC};
      cell_range[0].set_min_max(min,max);
    }
    MM_PP_Comm<PFMM,NC,NC,NC> mm_pp_comm;
    mm_pp_comm.set_targets(targets, cell_range);

    puts("Gen Green");
	MM_Cell<PFMM,NC,NC,NC, 3, 5,ICUT> mmcell;
	const double alpha = 2.4;
	mmcell.initialize(alpha,1./NC);
	mmcell.set_mm_mm_com(pmmm_comm, mm_list, mm_targets);
	
	puts("Gather M");
	// gather M
	//	recv_M<PFMM,NC,NC,NC>(mmcell.mm,0);
	{
	  int m;
	  for(m=0;m<mm_list.size();m++){
	    mm_pp_comm.recv_M(mmcell.mm,mm_list[m]);
	  }
	}

	puts("convolutin PBC");
	mmcell.convolution_PBC();

	// broadcast L
	//send_L<PFMM,NC,NC,NC>(mmcell.le,0);
	{
	  int m;
	  for(m=0;m<mm_list.size();m++){
	    mm_pp_comm.send_L(mmcell.le,mm_list[m]);
	  }
	}
  }

  puts("Done");

  MPI::Finalize();
  return EXIT_SUCCESS;

  
}
