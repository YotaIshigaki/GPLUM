#include<particle_simulator.hpp>
#include<particle_mesh.hpp>
int main(int argc, char* argv[]){
	PS::Initialize(argc, argv);
	PS::PM::ParticleMesh pm;
	PS::Finalize();
	return 0;
}

