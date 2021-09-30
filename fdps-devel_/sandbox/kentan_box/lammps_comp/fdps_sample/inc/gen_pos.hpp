#ifndef HPP_GEN_POS
#define HPP_GEN_POS

#include <particle_simulator.hpp>
#include <cassert>

/*
Requirement:
3D Periodic boundary condition
N = 4*nx*ny*nz
*/
template<class Tpsys,class Tdinfo>
void generateFCC(Tpsys &psys,Tdinfo &dinfo,
		 const PS::S32 nx,
		 const PS::S32 ny,
		 const PS::S32 nz,
		 const PS::F64 rho
		 ){
  const PS::S32 n_in_uc = 4; // # of particle in unit cell
  const PS::S32 natoms = n_in_uc * nx * ny * nz;

  const PS::F64 l_uc = cbrt((PS::F64)natoms/(rho*nx*ny*nz));
  const PS::F64vec L(l_uc*nx,l_uc*ny,l_uc*nz);
  const PS::F64vec ul[4] = {PS::F64vec(      0.,      0.,      0.),
			    PS::F64vec(0.5*l_uc,0.5*l_uc,      0.),
			    PS::F64vec(0.5*l_uc,      0.,0.5*l_uc),
			    PS::F64vec(      0.,0.5*l_uc,0.5*l_uc)};

  const PS::F64ort root_domain(PS::F64vec(0.,0.,0.),PS::F64vec(L.x,L.y,L.z));
  dinfo.setPosRootDomain(root_domain.low_,root_domain.high_);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);

  // find lower and upper bound of my rank
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  PS::S32 ndomain[3] = {dinfo.getNDomain(0),
			dinfo.getNDomain(1),
			dinfo.getNDomain(2)};
  PS::S32 rank_3d[3];
  PS::S32 rank_glb = PS::Comm::getRank();
  for(PS::S32 d=2;d>=0;d--){
    rank_3d[d] = rank_glb % ndomain[d];
    rank_glb /= ndomain[d];
  }
  PS::F64vec domain_size(L.x/ndomain[0],L.y/ndomain[1],L.z/ndomain[2]);
  PS::F64ort my_domain(PS::F64vec(domain_size.x*rank_3d[0],
				  domain_size.y*rank_3d[1],
				  domain_size.z*rank_3d[2]),
		       PS::F64vec(domain_size.x*(rank_3d[0]+1),
				  domain_size.y*(rank_3d[1]+1),
				  domain_size.z*(rank_3d[2]+1)));
#else
  const PS::F64ort my_domain = root_domain;
#endif
  // count and set local N
  int natom_loc = 0;
  PS::F64vec pos;
  for(int x=0;x<nx;x++){
    pos.x = l_uc*x;
    for(int y=0;y<ny;y++){
      pos.y = l_uc*y;
      for(int z=0;z<nz;z++){
	pos.z = l_uc*z;
	if(my_domain.contained(pos + ul[0])) natom_loc++;
	if(my_domain.contained(pos + ul[1])) natom_loc++;
	if(my_domain.contained(pos + ul[2])) natom_loc++;
	if(my_domain.contained(pos + ul[3])) natom_loc++;
      }
    }
  }
  psys.setNumberOfParticleLocal(natom_loc);

  // generate fcc inside mydomain
  int count = 0;
  for(int x=0;x<nx;x++){
    pos.x = l_uc*x;
    for(int y=0;y<ny;y++){
      pos.y = l_uc*y;
      for(int z=0;z<nz;z++){
	pos.z = l_uc*z;
	PS::F64vec pos_tmp;
	pos_tmp = pos + ul[0];
	assert(root_domain.contained(pos_tmp));
	if(my_domain.contained(pos_tmp)) psys[count++].setPos(pos_tmp);
	pos_tmp = pos + ul[1];
	assert(root_domain.contained(pos_tmp));
	if(my_domain.contained(pos_tmp)) psys[count++].setPos(pos_tmp);
	pos_tmp = pos + ul[2];
	assert(root_domain.contained(pos_tmp));
	if(my_domain.contained(pos_tmp)) psys[count++].setPos(pos_tmp);
	pos_tmp = pos + ul[3];
	assert(root_domain.contained(pos_tmp));
	if(my_domain.contained(pos_tmp)) psys[count++].setPos(pos_tmp);
      }
    }
  }
  assert(count == natom_loc);
}
#endif
