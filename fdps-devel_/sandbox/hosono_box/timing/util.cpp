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
		//omp_set_num_threads(omp_get_max_threads()/2);
		std::cout << "# of thread is " << omp_get_num_threads() << std::endl;
		#else 
		std::cout << "No OpenMP Mode... " << std::endl;
		#endif
		std::cout << "//=====================================" << std::endl;
	}
	return ;
}

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>& sph_system){
	PS::F64vec Mom = 0;//total momentum
	PS::F64    Eng = 0;//total enegry
	PS::F64    Mass = 0;//total enegry
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		Mom += sph_system[i].vel * sph_system[i].mass;
		Eng += (sph_system[i].eng + 0.5 * sph_system[i].vel * sph_system[i].vel) * sph_system[i].mass;
		Mass += sph_system[i].mass;
	}
	Eng  = PS::Comm::getSum(Eng);
	Mom  = PS::Comm::getSum(Mom);
	Mass = PS::Comm::getSum(Mass);
	printf("%.16e\n", Mass);
	printf("%.16e\n", Eng);
	printf("%.16e\n", Mom.x);
	printf("%.16e\n", Mom.y);
	printf("%.16e\n", Mom.z);
}

void ShiftOrigin(PS::ParticleSystem<RealPtcl>& sph_system){
	PS::F64vec com = 0;//center of mass OF TARGET
	PS::F64 mass = 0;
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		if(sph_system[i].tag != 0) continue;
		com += sph_system[i].pos * sph_system[i].mass;
		mass += sph_system[i].mass;
	}
	com  = PS::Comm::getSum(com);
	mass = PS::Comm::getSum(mass);
	com /= mass;
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].pos -= com;
	}

}