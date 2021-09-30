#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


typedef enum ps_boundary_condition {
				    BC_OPEN,
				    BC_PERIODIC_X,
				    BC_PERIODIC_Y,
				    BC_PERIODIC_Z,
				    BC_PERIODIC_XY,
				    BC_PERIODIC_XZ,
				    BC_PERIODIC_YZ,
				    BC_PERIODIC_XYZ,
				    BC_SHEARING_BOX,
				    BC_USER_DEFINED,
}PS_BOUNDARY_CONDITION;


typedef enum ps_interaction_list_mode{
				      MAKE_LIST,
				      MAKE_LIST_FOR_REUSE,
				      REUSE_LIST,
}PS_INTERACTION_LIST_MODE;

#include "./FDPS_vector.h"
#include "FDPS_types.h"
#include "user_defined.h"
#include "FDPS_c_if.h"

void fdps_get_psys_cptr(const int psys_num,
			Full_particle **cptr);

void  calc_gravity_c_epep(PFull_particle ep_i,
			   int n_ip,
			   PFull_particle ep_j,
			   int n_jp,
			  PFull_particle f);


void  calc_gravity_c_epsp(PFull_particle ep_i,
			   int n_ip,
			   PSPJMonopole ep_j,
			   int n_jp,
			  PFull_particle f);

void calc_force_all_and_write_back(int tree_num,
				   void (*pfunc_ep_ep)(Full_particle*, int,
						       Full_particle*, int,
						       Full_particle*),
				   void (*pfunc_ep_sp)(Full_particle*, int,
						       SPJMonopole*, int,
						       Full_particle*),
				   int psys_num,
				   int dinfo_num,
				   bool clear,
				   PS_INTERACTION_LIST_MODE list_mode)
{
    fdps_calc_force_all_and_write_back_l00000(tree_num,
					      pfunc_ep_ep,
					      pfunc_ep_sp,
					      psys_num,  
					      dinfo_num,
					      clear,
					      list_mode);
}

void dump_fullp(Full_particle p)
{
    printf("%d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e",
	   p.id,   p.mass, p.pos.x, p.pos.y, p.pos.z,
	   p.vel.x, p.vel.y, p.vel.z);
    printf("%15.7e %15.7e %15.7e %15.7e\n",
	   p.acc.x, p.acc.y, p.acc.z, p.pot);
}
void dump_fullpsys(Full_particle *p, int n)
{
    for (int i=0;i<n;i++)dump_fullp(p[i]);
}

void dump_particles(int psys_num)
{
    PFull_particle ptcl;
    fdps_get_psys_cptr(psys_num,&ptcl);
    int n =  fdps_get_nptcl_loc(psys_num);
    dump_fullpsys(ptcl, n);
}



void setup_IC(int psys_num,
	      int nptcl_glb)
{

    double  m_tot=1.0;
    double rmax=3.0;
    double r2max=rmax*rmax;
    //   Get # of MPI processes and rank number
    int nprocs = fdps_get_num_procs();
    int myrank = fdps_get_rank();
    // Make an initial condition at RANK 0
    if (myrank == 0 ){
	//Set # of local particles
	fdps_set_nptcl_loc(psys_num,nptcl_glb);
	PFull_particle ptcl;
	fdps_get_psys_cptr(psys_num,&ptcl);
	//  ** initialize Mersenne twister
	fdps_mt_init_genrand(0);
	Cvec_Float64 zerov= Cvec_Float64_zero();
	for (int i=0; i < nptcl_glb; i++){
	    PFull_particle q = ptcl+i;
	    q->id = i;
	    q->mass = m_tot/nptcl_glb;
	    double r2 = r2max*2;
	    Cvec_Float64 v;
	    while (r2 >= r2max){
		v.x= (2*fdps_mt_genrand_res53()-1.0) * rmax;
		v.y= (2*fdps_mt_genrand_res53()-1.0) * rmax;
		v.z= (2*fdps_mt_genrand_res53()-1.0) * rmax;
		r2 = v.x*v.x+v.y*v.y+v.z*v.z;
	    }
	    q->pos = v;
	    q->vel = zerov;
	    q->eps = 1.0/32;
	}
	Cvec_Float64 cm_pos = zerov;
	Cvec_Float64 cm_vel = zerov;
	double cm_mass = 0;
	for (int i=0; i < nptcl_glb; i++){
	    PFull_particle pi = ptcl+i;
	    cm_pos.x +=  pi->pos.x* pi->mass;
	    cm_pos.y +=  pi->pos.y* pi->mass;
	    cm_pos.z +=  pi->pos.z* pi->mass;
	    cm_vel.x +=  pi->vel.x* pi->mass;
	    cm_vel.y +=  pi->vel.y* pi->mass;
	    cm_vel.z +=  pi->vel.z* pi->mass;
	    cm_mass += pi->mass;
	}
	cm_pos.x /= cm_mass;
	cm_pos.y /= cm_mass;
	cm_pos.z /= cm_mass;
	cm_vel.x /= cm_mass;
	cm_vel.y /= cm_mass;
	cm_vel.z /= cm_mass;

	for (int i=0; i < nptcl_glb; i++){
	    PFull_particle q = ptcl+i;
	    q->pos.x -= cm_pos.x;
	    q->pos.y -= cm_pos.y;
	    q->pos.z -= cm_pos.z;
	    q->vel.x -= cm_vel.x;
	    q->vel.y -= cm_vel.y;
	    q->vel.z -= cm_vel.z;
	}
	//	dump_fullpsys(ptcl, nptcl_glb);
    } else{
	fdps_set_nptcl_loc(psys_num,0);
    }
}

void  calc_energy(int psys_num,
		  double *etot,
		  double *ekin,
		  double *epot)
{
    *etot = *ekin = *epot = 0;
    int nptcl_loc = fdps_get_nptcl_loc(psys_num);
    PFull_particle ptcl;
    fdps_get_psys_cptr(psys_num,&ptcl);
    
    double  ekin_loc = 0;
    double  epot_loc = 0;
    for (int i=0;i < nptcl_loc; i++){
	PFull_particle pi = ptcl+i;
	Cvec_Float64 v = pi->vel;
	ekin_loc += pi->mass * (v.x*v.x+v.y*v.y+v.z*v.z);
	epot_loc += pi->mass * (pi->pot + pi->mass/pi->eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    double etot_loc = ekin_loc + epot_loc;
    fdps_get_sum_r64(ekin_loc, ekin);
    fdps_get_sum_r64(epot_loc, epot);
    fdps_get_sum_r64(etot_loc, etot);
}

void  kick(int psys_num,
	   double dt)
{
    PFull_particle ptcl;
    fdps_get_psys_cptr(psys_num,&ptcl);
    int n =   fdps_get_nptcl_loc(psys_num);
    for (int i=0;i < n; i++){
	PFull_particle pi = ptcl+i;
	Cvec_Float64 *pv, *pa;
	pv = &(pi->vel);
	pa = &(pi->acc);
	pv->x += pa->x * dt;
	pv->y += pa->y * dt;
	pv->z += pa->z * dt;
    }
}

void  drift(int psys_num,
	   double dt)
{
    PFull_particle ptcl;
    fdps_get_psys_cptr(psys_num,&ptcl);
    int n =   fdps_get_nptcl_loc(psys_num);
    for (int i=0;i < n; i++){
	PFull_particle pi = ptcl+i;
	Cvec_Float64 *px, *pv; 
	pv = &(pi->vel);
	px = &(pi->pos);
	px->x += pv->x * dt;
	px->y += pv->y * dt;
	px->z += pv->z * dt;
    }
}


int cmain()
{
    fprintf(stderr, "FDPS on C test code\n");
    fdps_initialize();
    int dinfo_num=1;
    float coef_ema=0.3;
    fdps_create_dinfo(&dinfo_num);
    fdps_init_dinfo(dinfo_num,coef_ema);
    int psys_num = 0;
    fdps_create_psys(&psys_num,"full_particle");
    fdps_init_psys(psys_num);

    //  Create tree object
    int  tree_num = 0;
    fdps_create_tree(&tree_num, 
		     "Long,full_particle,full_particle,full_particle,Monopole");
    int  ntot=1024;
    double theta = 0.5;
    int n_leaf_limit = 8;
    int   n_group_limit = 64;
    fdps_init_tree(tree_num,ntot, theta, n_leaf_limit, n_group_limit);

    setup_IC(psys_num,ntot);
    // !* Domain decomposition and exchange particle
    fdps_decompose_domain_all(dinfo_num,psys_num,1);
    fdps_exchange_particle(psys_num,dinfo_num);
   //* Compute force at the initial time
    calc_force_all_and_write_back(tree_num, &calc_gravity_c_epep,
				  &calc_gravity_c_epsp,
				  psys_num,
				  dinfo_num,
				  1,
				  MAKE_LIST);
    //    dump_particles(psys_num);
   //* Compute energies at the initial time
    double etot0,ekin0,epot0;
    calc_energy(psys_num, &etot0, &ekin0,&epot0);
    printf( "Energies = %21.14e  %21.14e  %21.14e\n",etot0,ekin0,epot0);
   //* Time integration
    double   time_diag = 0;
    double   time_snap = 0;
    double   time_sys  = 0;
    double   time_end = 10.0;
    double   dt = 1.0/128.0;
    double   dt_diag = 1.0;
    double   dt_snap = 1.0;
    int   num_loop = 0;
    while (time_sys <= time_end){
	if (time_sys + dt/2 >= time_snap){
	    //       output(psys_num)
	    time_snap += dt_snap;
	}
	double etot1, ekin1, epot1;
	calc_energy(psys_num, &etot1,&ekin1,&epot1);
	//	printf( "Energies = %21.14e  %21.14e  %21.14e\n",etot1,ekin1,epot1);
	//	dump_particles(psys_num);
	if (fdps_get_rank() == 0){
	    if (time_sys + dt/2 >= time_diag){
		printf ("time: %10.3f, energy error: %15.7e\n",
			time_sys, (etot1-etot0)/etot0);
		time_diag = time_diag + dt_diag;
		//        fdps_get_tree_time_prof(tree_num,pointerof(tprof))
		//         p tprof
	    }
	}
	kick(psys_num,0.5*dt);
	time_sys +=  dt;
	drift(psys_num,dt);
	//   !* Domain decomposition & exchange particle
	if (num_loop%4 == 0){
	    fdps_decompose_domain_all(dinfo_num,psys_num,1);
	}
	fdps_exchange_particle(psys_num,dinfo_num);
	//#    !* Force calculation
	calc_force_all_and_write_back(tree_num, &calc_gravity_c_epep,
				      &calc_gravity_c_epsp,
				      psys_num,
				      dinfo_num,
				      1,
				      MAKE_LIST);
	kick(psys_num,0.5*dt);
	num_loop += 1;
    }
    fdps_finalize();
}
