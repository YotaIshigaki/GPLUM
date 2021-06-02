#include "header.h"

void DisplayInfo(void){
	if(PS::Comm::getRank() == 0){
#if 0
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
#endif
		std::cout << "This is a sample program of Smoothed Particle Hydrodynamics on FDPS!" << std::endl;
		std::cout << "# of proc is   " << PS::Comm::getNumberOfProc() << std::endl;
		#ifdef _OPENMP
		std::cout << "# of thread is " << PS::Comm::getNumberOfThread() << std::endl;
		#else 
		std::cout << "No OpenMP Mode... " << std::endl;
		#endif
		std::cout << "//=====================================" << std::endl;
	}
	return ;
}

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>& sph_system){
	static PS::F64vec Mom_init = 0;//total momentum
	static PS::F64    Eng_init = 0;//total enegry
	static PS::F64    Eng_kin_init = 0;//total enegry
	static PS::F64    Eng_U_init = 0;//total enegry
	static PS::F64    Eng_pot_init = 0;//total enegry
	static PS::F64    Mass_init = 0;//total enegry
        static bool first = true;
	PS::F64vec Mom = 0;//total momentum
	PS::F64    Eng = 0;//total enegry
	PS::F64    Eng_kin = 0;//total enegry
	PS::F64    Eng_U = 0;//total enegry
	PS::F64    Eng_pot = 0;//total enegry
	PS::F64    Mass = 0;//total enegry
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		Mom += sph_system[i].vel * sph_system[i].mass;
		//Eng += (sph_system[i].eng + 0.5 * sph_system[i].vel * sph_system[i].vel) * sph_system[i].mass;
                Eng += (sph_system[i].eng + 0.5 * sph_system[i].vel * sph_system[i].vel + 0.5*sph_system[i].pot) * sph_system[i].mass;
                Eng_kin += 0.5 * sph_system[i].vel * sph_system[i].vel * sph_system[i].mass;
                Eng_U += sph_system[i].eng * sph_system[i].mass;
                Eng_pot += 0.5*sph_system[i].pot * sph_system[i].mass;
		Mass += sph_system[i].mass;
		//std::cout << sph_system[i].mass << std::endl;
	}
	Eng  = PS::Comm::getSum(Eng);
	Eng_U  = PS::Comm::getSum(Eng_U);
	Eng_kin  = PS::Comm::getSum(Eng_kin);
	Eng_pot  = PS::Comm::getSum(Eng_pot);
	Mom  = PS::Comm::getSum(Mom);
	Mass = PS::Comm::getSum(Mass);
        if(first){
            Mom_init = Mom;
            Eng_init = Eng;
            Eng_U_init = Eng_U;
            Eng_kin_init = Eng_kin;
            Eng_pot_init = Eng_pot;
            Mass_init = Mass;
            first = false;
        }
        if(PS::Comm::getRank()==0){
            std::cerr<<" (Eng-Eng_init)/Eng_init= "<<(Eng-Eng_init)/Eng_init
                     <<" Eng "<<Eng
                     <<" Eng_U= "<<Eng_U
                     <<" Eng_kin= "<<Eng_kin
                     <<" Eng_pot= "<<Eng_pot
                     <<std::endl;
            std::cout<<" (Eng-Eng_init)/Eng_init= "<<(Eng-Eng_init)/Eng_init
                     <<" Eng "<<Eng
                     <<" Eng_U= "<<Eng_U
                     <<" Eng_kin= "<<Eng_kin
                     <<" Eng_pot= "<<Eng_pot
                     <<std::endl;
        }
        /*
	printf("%.16e\n", Mass);
	printf("%.16e\n", Eng);
	printf("%.16e\n", Mom.x);
	printf("%.16e\n", Mom.y);
	printf("%.16e\n", Mom.z);
        */


}
