//
// main.cpp
//
// MD simulation code for water molecules using LIPS-SW method
//
// Kentaro NOMURA 2017/12/22
//
// Known problem:
// nothing so far

#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>

#include <string>
#include <vector>

#include <unit.h>
#include <pdb_manager.h>
#include <water_params.h>

#include "user_defined_class.h"

const PS::F64 alpha[10] = {
  0.19578,
  0.0125143224110408,
  -0.603493863454666,
  11.7355819865242,
  -96.296895305654,
  216.649868508398,
  -197.409191110696,
  59.9544311773618,
  13.9564907382725,
  -8.66620089071555
};

struct CalcForceEpEp{
  const PS::F64 rc_lj,rc_cl;

  CalcForceEpEp(const PS::F64 _rc_lj,const PS::F64 _rc_cl)
    : rc_lj(_rc_lj), rc_cl(_rc_cl){}

  template<class EPJ>
  void operator () (const EP * ep_i,
		    const PS::S32 n_ip,
		    const EPJ * ep_j,
		    const PS::S32 n_jp,
		    Force * force){
    const PS::F64 rc_lj_sq = rc_lj*rc_lj;
    const PS::F64 rc_cl_sq = rc_cl*rc_cl;
    const PS::F64 rc_inv = 1.0 / rc_cl;
    const PS::F64 rc2_inv = rc_inv * rc_inv;
    const PS::F64 rc3_inv = rc2_inv * rc_inv;

    const PS::F64 bound_ele_pot = rc_inv*1.296557;

    for(PS::S32 i=0; i<n_ip; i++){
      const PS::F64vec ri = ep_i[i].pos;
      PS::F64vec flj = 0.0;
      PS::F64vec fcl = 0.0;
      PS::F64 poti = 0.0;
      const PS::S64 idi = ep_i[i].id;
      const PS::F64 sigma_i   = ep_i[i].sigma;
      const PS::F64 epsilon_i = ep_i[i].epsilon;
      const PS::F64 charge_i  = ep_i[i].charge;
      for(PS::S32 j=0; j<n_jp; j++){
	PS::F64vec rij = ri - ep_j[j].pos;

	const PS::F64 r2 = rij*rij;
	if((idi/3) == (ep_j[j].id/3) && r2 < rc_lj_sq) continue;
	const PS::F64 r2_inv = 1.0/r2;
 	if(0.0 < r2 && r2 <= rc_cl_sq){
	  const PS::F64 r = sqrt(r2);
	  const PS::F64 r_inv = r2_inv * r;
	  const PS::F64 qq = charge_i * ep_j[j].charge;
	  const PS::F64 rc = r * rc_inv;
	  const PS::F64 r2c2 = rc*rc;
	  const PS::F64 coef = (rc*rc - alpha[0]*alpha[0]);
	  const PS::F64 coef2 = coef*coef;
	  const PS::F64 utmp = r2c2 * (alpha[1] + r2c2*
				       (alpha[2] + r2c2*
					(alpha[3] + r2c2*
					 (alpha[4] + r2c2*
					  (alpha[5] + r2c2*
					   (alpha[6] + r2c2*
					    (alpha[7] + r2c2*
					     (alpha[8] + (alpha[9] * r2c2)))))))));
	  const PS::F64 ftmp = (2.0*alpha[1] + r2c2*
				(4.0*alpha[2] + r2c2*
				 (6.0*alpha[3] + r2c2*
				  (8.0*alpha[4] + r2c2*
				   (10.0*alpha[5] + r2c2*
				    (12.0*alpha[6] + r2c2*
				     (14.0*alpha[7] + r2c2*
				      (16.0*alpha[8] + (18.0*alpha[9] * r2c2)))))))));
	  poti += qq * (r_inv - 0.5*rc_inv * coef * coef2 * utmp - bound_ele_pot);
	  fcl  += qq * (r2_inv*r_inv + 0.5*rc3_inv*coef2*(6.0*utmp + coef*ftmp))*rij;
	}
	if (0.0 < r2 && r2 <= rc_lj_sq) {
	  const PS::F64 sigma    = 0.5*(sigma_i + ep_j[j].sigma);
	  const PS::F64 epsilon  = 4.0*sqrt(epsilon_i * ep_j[j].epsilon);
	  const PS::F64 rs2_inv  = sigma*sigma * r2_inv;
	  const PS::F64 rs6_inv  = rs2_inv * rs2_inv * rs2_inv;
	  const PS::F64 rs12_inv = rs6_inv * rs6_inv;
	  poti += epsilon * (rs12_inv - rs6_inv);
	  flj  += (epsilon * (12.0*rs12_inv - 6.0*rs6_inv) * r2_inv) * rij;
	}
      }
      force[i].acc += flj + fcl;
      force[i].pot += 0.5*poti;
    }
  }
};

template<class Tpsys>
void ReadPDBFile(Tpsys &psys,
		 const std::string filename){
  FileManager<PDBData> pdb;
  std::ifstream ifs_test(filename);
  assert(ifs_test.is_open());

  int natom = 0;
  PS::F64vec pos_min(1e32,1e32,1e32);
  for(std::string line;std::getline(ifs_test,line);){
    PDBData data = pdb.ReadLine(line);
    if(pos_min.x > data.coord[0]) pos_min.x = data.coord[0];
    if(pos_min.y > data.coord[1]) pos_min.y = data.coord[1];
    if(pos_min.z > data.coord[2]) pos_min.z = data.coord[2];
    if(data.record == "HETATM") natom++;
  }

  psys.setNumberOfParticleLocal(natom);
  PS::MTTS mt;
  mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());

  std::ifstream ifs(filename);
  assert(ifs.is_open());

  int count = 0;
  for(std::string line;std::getline(ifs,line);){
    PDBData data = pdb.ReadLine(line);
    if(data.record == "HETATM"){
      if(count%3==0) assert(data.resName == "O  ");
      else           assert(data.resName == "H  ");

      psys[count].pos.x = data.coord[0] - pos_min.x;
      psys[count].pos.y = data.coord[1] - pos_min.y;
      psys[count].pos.z = data.coord[2] - pos_min.z;

      const double v_max = 0.1;
      do {
	psys[count].vel[0] = (2.*mt.genrand_res53() - 1.) * v_max;
	psys[count].vel[1] = (2.*mt.genrand_res53() - 1.) * v_max;
	psys[count].vel[2] = (2.*mt.genrand_res53() - 1.) * v_max;
      }while(psys[count].vel * psys[count].vel >= v_max * v_max);

      psys[count].id      = count;
      psys[count].type    = count%3==0 ? 0 : 1;
      for(int c=0;c<3;c++){
	const int offset = (count / 3)*3;
	psys[count].cid[c] = offset + c;
      }
      psys[count].mass    = (count%3==0 ? MASS_OXY    : MASS_HYD);
      psys[count].sigma   = (count%3==0 ? SIGMA_OXY   : SIGMA_HYD);
      psys[count].epsilon = (count%3==0 ? EPSILON_OXY : EPSILON_HYD);
      psys[count].charge  = (count%3==0 ? CHARGE_OXY  : CHARGE_HYD);
      count++;
    }
  }
  assert(count == natom);
}

template<class Tpsys>
void MakeIceLattice(Tpsys &psys,
		    PS::F64vec &cell_size){

  const PS::S32 nx=8,ny=8,nz=4;
  const PS::F64vec unit(4.511,4.511,7.315);
  cell_size.x = unit.x * nx * sin(M_PI*60./180.);
  cell_size.y = unit.y * ny;
  cell_size.z = unit.z * nz;
  fprintf(stdout,"cellsize: %lf %lf %lf\n",cell_size.x,cell_size.y,cell_size.z);
  
  if(PS::Comm::getRank()==0){
    ReadPDBFile(psys,"ice_unit.pdb");

    const int nunit = psys.getNumberOfParticleLocal();
    const int natom = nunit*nx*ny*nz;
    psys.setNumberOfParticleLocal(natom);
    int count = nunit;
    for(int x=0;x<nx;x++){
      for(int y=0;y<ny;y++){
	for(int z=0;z<nz;z++){
	  if(x==0 && y==0 && z==0)continue;
	  for(int i=0;i<12;i++){
	    psys[count] = psys[i];

	    psys[count].pos.x += unit.x*x*sin(M_PI*60./180.);
	    psys[count].pos.y += unit.x*x*cos(M_PI*60./180.) + unit.y*y;
	    psys[count].pos.z += unit.z*z;
	    psys[count].id = count;
	    for(int c=0;c<3;c++){
	      const int offset = (count / 3)*3;
	      psys[count].cid[c] = offset + c;
	    }
	    count++;
	  }
	}}}
    assert(count == natom);
  }
}

template<class Tpsys,class Tdinfo>
void PeriodicBoundaryCondition(Tpsys & system,const Tdinfo &dinfo){
  PS::S32 n = system.getNumberOfParticleLocal();
  const PS::F64ort domain = dinfo.getPosRootDomain();
  const PS::F64vec length = domain.getFullLength();

  for(int i=0; i<n; i++){
    for(int k=0;k<3;k++){
      if(system[i].pos.x <  domain.low_.x){
	system[i].pos.x  += length.x;
      }
      if(system[i].pos.x >= domain.high_.x){
	system[i].pos.x  -= length.x;
      }
      if(system[i].pos.y <  domain.low_.y){
	system[i].pos.y  += length.y;
      }
      if(system[i].pos.y >= domain.high_.y){
	system[i].pos.y  -= length.y;
      }
      if(system[i].pos.z <  domain.low_.z){
	system[i].pos.z  += length.z;
      }
      if(system[i].pos.z >= domain.high_.z){
	system[i].pos.z  -= length.z;
      }
    }
  }
}

template<class Tpsys>
void RemoveTotalMomentum(Tpsys &system){
  const PS::S32 n_loc = system.getNumberOfParticleLocal();
  PS::F64vec cm_vel_loc = 0.0;
  PS::F64    cm_mass_loc = 0.0;
  for(int i=0; i<n_loc; i++){
    cm_vel_loc  += system[i].mass * system[i].vel;
    cm_mass_loc += system[i].mass;
  }
  PS::F64vec cm_vel;
  PS::F64    cm_mass;
  cm_vel  = PS::Comm::getSum(cm_vel_loc);
  cm_mass = PS::Comm::getSum(cm_mass_loc);
  cm_vel /= cm_mass;
  for(int i=0; i<n_loc; i++){
    system[i].vel -= cm_vel;
  }
}

template<class Tpsys>
void ScaleVelocity(Tpsys & system,const PS::F64 T){
  const PS::S32 natom_local = system.getNumberOfParticleLocal();
  PS::F64 ekin_loc = 0.0;
  for(PS::S32 i=0; i<natom_local; i++){
    ekin_loc += system[i].mass * system[i].vel * system[i].vel;
  }
  ekin_loc *= 0.5;
  PS::S32 natom;
  PS::F64 ekin;
  ekin  = PS::Comm::getSum(ekin_loc);
  natom = PS::Comm::getSum(natom_local);
  const PS::F64 scaler = sqrt(1.5*natom*T / ekin);
  for(PS::S32 i=0;i<natom_local;i++) system[i].vel *= scaler;

  RemoveTotalMomentum(system);
}

template<class Tpsys>
void CalcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
                const bool clear=true){
  if(clear){
    etot = ekin = epot = 0.0;
  }
  PS::F64 ekin_loc = 0.0;
  PS::F64 epot_loc = 0.0;
  const PS::S32 nbody = system.getNumberOfParticleLocal();
  for(PS::S32 i=0; i<nbody; i++){
    ekin_loc += system[i].mass * (system[i].vel * system[i].vel);
    epot_loc += system[i].pot;
  }
  ekin_loc *= 0.5;

  epot = PS::Comm::getSum(epot_loc);
  ekin = PS::Comm::getSum(ekin_loc);
  etot = epot + ekin;
}

  template <class Ttree>
  void CalcIntraForce(FP &fp,Ttree &tree,const PS::F64vec box_size){
    EP ep[3];
    PS::S32 myid;
    for(int i=0;i<3;i++){
      if(fp.cid[i]==fp.id){
	ep[i].copyFromFP(fp);
	myid = i;
      }else{
	assert(fp.cid[i] >= 0);
	EP* tmp = tree.getEpjFromId(fp.cid[i]);
	assert(tmp != NULL);
	ep[i] = *tmp;
	// shift ep if ep is image particle in image cell
	while(fp.pos.x - ep[i].pos.x >= 0.5*box_size.x) ep[i].pos.x += box_size.x;
	while(fp.pos.x - ep[i].pos.x < -0.5*box_size.x) ep[i].pos.x -= box_size.x;
	while(fp.pos.y - ep[i].pos.y >= 0.5*box_size.y) ep[i].pos.y += box_size.y;
	while(fp.pos.y - ep[i].pos.y < -0.5*box_size.y) ep[i].pos.y -= box_size.y;
	while(fp.pos.z - ep[i].pos.z >= 0.5*box_size.z) ep[i].pos.z += box_size.z;
	while(fp.pos.z - ep[i].pos.z < -0.5*box_size.z) ep[i].pos.z -= box_size.z;
      }
    }

    PS::F64vec drij[2] = {ep[1].pos-ep[0].pos, ep[2].pos-ep[0].pos};
    PS::F64     rij[2] = {sqrt(drij[0]*drij[0]),sqrt(drij[1]*drij[1])};

    // bond
    if(myid == 0){ // case that fp is Oxygen
      const PS::F64 dij = rij[0] - BOND_OH;
      const PS::F64 dkj = rij[1] - BOND_OH;
      fp.pot += 0.5*BOND_COEF*(dij*dij + dkj*dkj);
      fp.acc += (2.0*BOND_COEF * dij / rij[0]) * drij[0];
      fp.acc += (2.0*BOND_COEF * dkj / rij[1]) * drij[1];
    }else{
      int bond_id = myid - 1;
      const PS::F64 dij = rij[bond_id] - BOND_OH;
      fp.pot += 0.5*BOND_COEF*dij*dij;
      fp.acc -= (2.0*BOND_COEF*dij/rij[bond_id]) * drij[bond_id];
    }
    // angle
    const PS::F64 ip = drij[0] * drij[1]; // inner product
    const PS::F64 rr_inv = 1.0 /(rij[0]*rij[1]);
    const PS::F64 cos_ijk = ip * rr_inv;
    const PS::F64 sin_ijk = sqrt(1.0 - cos_ijk*cos_ijk);
    const PS::F64 dijk = acos(cos_ijk) - ANGLE_HOH;
    fp.pot += ANGLE_COEF*dijk*dijk / 3.0;
    if(myid == 0){
      fp.acc -= 2.0*ANGLE_COEF*dijk / (rij[0]*rij[0]*sin_ijk) * (drij[1] - cos_ijk*drij[0]);
      fp.acc -= 2.0*ANGLE_COEF*dijk / (rij[1]*rij[1]*sin_ijk) * (drij[0] - cos_ijk*drij[1]);
    }else{
      myid -= 1;
      PS::S32 yourid = (myid==0)? 1:0;
      fp.acc += 2.0*ANGLE_COEF*dijk / (rij[myid]*rij[myid]*sin_ijk) * (drij[yourid] - cos_ijk*drij[myid]);
    }
  }


int main(int argc, char *argv[]){
  PS::Initialize(argc, argv);

  PS::F64 Tbegin = PS::GetWtime();
  std::cout<<std::setprecision(15);
  std::cerr<<std::setprecision(15);

  PS::F64 dt    = 0.2 * 1e-15 / unit_time; // 2 fs
  PS::F64 temperature = 270.0 / unit_temp;

  PS::S32 n_group_limit = 64;
  PS::S32 n_leaf_limit = 8;

  PS::S32 nstep = 10;
  char dir_name[1024];

  int c;
  sprintf(dir_name,"./result");
  while((c=getopt(argc,argv,"o:s:n:h")) != -1){
    switch(c){
    case 'o':
      sprintf(dir_name,optarg);
      break;
    case 's':
      nstep = atoi(optarg);
      if(PS::Comm::getRank()==0)std::cerr<<"nstep="<<nstep<<std::endl;
      break;
    case 'n':
      n_group_limit = atoi(optarg);
      if(PS::Comm::getRank()==0)std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
      break;
    case 'h':
      if(PS::Comm::getRank()==0){
	std::cerr<<"s: number of steps (default: 1000)"<<std::endl;
	std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
	std::cerr<<"n: n_group_limit (default: 64.0)"<<std::endl;
      }
      return 0;
    }
  }

  struct stat st;
  if(stat(dir_name, &st) != 0) {
    PS::S32 rank = PS::Comm::getRank();
    PS::S32 ret_loc, ret=0;
    if(rank == 0)
      ret_loc = mkdir(dir_name, 0777);
    PS::Comm::broadcast(&ret_loc, ret);
    if(ret == 0) {
      if(rank == 0)
	fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
      fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
      PS::Abort();
      exit(0);
    }
  }

  std::ofstream fout_eng;
  std::ofstream fout_tcal;
  if(PS::Comm::getRank() == 0){
    char sout_de[1024];
    char sout_tcal[1024];
    sprintf(sout_de, "%s/t-de.dat", dir_name);
    sprintf(sout_tcal, "%s/t-tcal.dat", dir_name);
    std::cerr<<sout_de<<std::endl;
    std::cerr<<sout_tcal<<std::endl;
    fout_eng.open(sout_de);
    fout_tcal.open(sout_tcal);
  }

  PS::ParticleSystem<FP> system_water;
  system_water.initialize();

  PS::F64vec box_size;
  MakeIceLattice(system_water,box_size);

  const PS::F64 coef_ema = 0.3;
  PS::DomainInfo dinfo;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
			 PS::F64vec(box_size.x,box_size.y,box_size.z));
  PeriodicBoundaryCondition(system_water,dinfo);

  dinfo.collectSampleParticle(system_water);
  dinfo.decomposeDomain();
  system_water.exchangeParticle(dinfo);

  PS::TreeForForceShort<Force, EP, EP>::Scatter ips_tree;
  const PS::F64 rcut_lj = 4.0*SIGMA_OXY;
  const PS::F64 rcut_cl = 28.0; // 28 angstrom
  for(int i=0;i<system_water.getNumberOfParticleLocal();i++)
    system_water[i].search_radius = rcut_cl;
  ips_tree.initialize(system_water.getNumberOfParticleLocal(),0.0, n_leaf_limit, n_group_limit);
  ips_tree.calcForceAllAndWriteBack(CalcForceEpEp(rcut_lj,rcut_cl),
				    system_water,dinfo);
  for(int i=0;i<system_water.getNumberOfParticleLocal();i++){
    CalcIntraForce(system_water[i],ips_tree,box_size);
  }
  ScaleVelocity(system_water,temperature);

  PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
  CalcEnergy(system_water, Etot0, Ekin0, Epot0);
  if(PS::Comm::getRank() == 0) fprintf(stdout,"Etot = %lf, Epot = %lf, Ekin = %lf\n",Etot0,Epot0,Ekin0);

  // Main roop
  const PS::S32 nskip = 4;
  PS::F64 time_sys = 0.0;
  system_water.clearTimeProfile();
  dinfo.clearTimeProfile();
  ips_tree.clearTimeProfile();
  for(int s=0;s<nstep;s++){
    for(int i=0;i<system_water.getNumberOfParticleLocal();i++){
      system_water[i].IntegrateVel(0.5*dt);
      system_water[i].IntegratePos(dt,box_size);
    }
    PeriodicBoundaryCondition(system_water,dinfo);
    dinfo.collectSampleParticle(system_water);
    dinfo.decomposeDomain();
    system_water.exchangeParticle(dinfo);

    ips_tree.calcForceAllAndWriteBack(CalcForceEpEp(rcut_lj,rcut_cl),
				      system_water, dinfo);
    for(int i=0;i<system_water.getNumberOfParticleLocal();i++){
      CalcIntraForce(system_water[i],ips_tree,box_size);
    }

    for(int i=0;i<system_water.getNumberOfParticleLocal();i++){
      system_water[i].IntegrateVel(0.5*dt);
    }

    CalcEnergy(system_water, Etot1, Ekin1, Epot1);
    if(PS::Comm::getRank() == 0){
      fout_eng<<time_sys<<"   "<< " " << Epot1 << " " << Ekin1 << " " <<(Etot1-Etot0)/Etot0<<std::endl;
      fprintf(stderr, "%10.7f %lf %lf %+e\n",
	      time_sys, Epot1, Ekin1, (Etot1 - Etot0) / Etot0);
    }
    time_sys += dt;
  }
  PS::TimeProfile tp = ips_tree.getTimeProfile();
  fprintf(stdout,"Total time of force calculation: %e sec\n",tp.getTotalTime());
  fprintf(stdout,"\n");
  fprintf(stdout,"calc_force: %e sec\n",tp.calc_force);
  fprintf(stdout,"  calc_force__core: %e sec\n",tp.calc_force__core);
  fprintf(stdout,"    calc_force__core__walk_tree: %e sec\n",tp.calc_force__core__walk_tree);
  fprintf(stdout,"    calc_force__core__keep_list: %e sec\n",tp.calc_force__core__keep_list);
  fprintf(stdout,"    calc_force__core__copy_ep:   %e sec\n",tp.calc_force__core__copy_ep);
  fprintf(stdout,"    calc_force__core__dispatch:  %e sec\n",tp.calc_force__core__dispatch);
  fprintf(stdout,"    calc_force__core__retrieve:  %e sec\n",tp.calc_force__core__retrieve);
  fprintf(stdout,"\n");
  fprintf(stdout,"TimeProfile.dump():\n");
  tp = tp + system_water.getTimeProfile();
  tp = tp + dinfo.getTimeProfile();
  tp.dump();
  PS::Finalize();
  return 0;
}
