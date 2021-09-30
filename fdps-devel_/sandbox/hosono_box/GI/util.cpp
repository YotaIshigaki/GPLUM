#include "header.h"

void DisplayInfo(void){
	if(PS::Comm::getRank() == 0){
		std::cout << "//==================================\\\\" << std::endl;
		std::cout << "||                                  ||" << std::endl;
		std::cout << "|| ::::::: ::::::. ::::::. .::::::. ||" << std::endl;
		std::cout << "|| ::      ::    : ::    : ::       ||" << std::endl;
		std::cout << "|| ::::::  ::    : ::::::'  `:::::. ||" << std::endl;
		std::cout << "|| ::      ::::::' ::      `......' ||" << std::endl;
		std::cout << "||     Framework for Developing     ||" << std::endl;
		std::cout << "||        Particle Simulator        ||" << std::endl;
		std::cout << "\\\\==================================//" << std::endl;
		std::cout << "//=====================================" << std::endl;
		std::cout << "Working on FDPS!" << std::endl;
		std::cout << "# of proc is   " << PS::Comm::getNumberOfProc() << std::endl;
		#ifdef _OPENMP
		#pragma omp parallel
		{
			#pragma omp single
			{
				std::cout << "# of thread is " << omp_get_num_threads() << std::endl;
			}
		}
		#else 
		std::cout << "No OpenMP Mode... " << std::endl;
		#endif
		std::cout << "//=====================================" << std::endl;
	}
	return ;
}

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>& sph_system, system_t sysinfo){
	PS::F64vec Mom  = 0;//total momentum
	PS::F64    Eng  = 0;//total enegry
	PS::F64    Mass = 0;//total enegry
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		Mom += sph_system[i].vel * sph_system[i].mass;
		Eng += (sph_system[i].eng + 0.5 * sph_system[i].vel * sph_system[i].vel + 0.5 * sph_system[i].pot * sysinfo.Grav) * sph_system[i].mass;
		Mass += sph_system[i].mass;
	}
	Eng  = PS::Comm::getSum(Eng);
	Mom  = PS::Comm::getSum(Mom);
	Mass = PS::Comm::getSum(Mass);
	printf("Mass = %.16e\n", Mass);
	printf("Eng  = %.16e\n", Eng);
	printf("Mom = (%.16e, %.16e, %.16e)\n", Mom.x, Mom.y, Mom.z);
}

void ShiftOrigin(PS::ParticleSystem<RealPtcl>& sph_system){
	PS::F64vec com_loc = 0;//center of mass of TARGET
	PS::F64vec mom_loc = 0;//moment of TARGET
	PS::F64 mass_loc = 1.0e-16;//mass of TARGET
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		if(sph_system[i].tag != 0) continue;
		com_loc += sph_system[i].pos * sph_system[i].mass;
		mom_loc += sph_system[i].vel * sph_system[i].mass;
		mass_loc += sph_system[i].mass;
	}
	PS::Comm::barrier();
	PS::F64vec com  = PS::Comm::getSum(com_loc);
	PS::F64vec mom  = PS::Comm::getSum(mom_loc);
	PS::F64 mass = PS::Comm::getSum(mass_loc);
	com /= mass;
	mom /= mass;
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].pos -= com;
		//OOPS!
		//sph_system[i].pos -= mom;
		sph_system[i].vel -= mom;
	}

}
