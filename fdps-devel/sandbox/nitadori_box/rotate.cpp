#include <cstdlib>
#include <cmath>
#include "../../src/vector3.hpp"

typedef ParticleSimulator::Vector3<double> F64vec;

F64vec rotate_fwd(const F64vec v){
	return F64vec(
			( 1./sqrt(2.)) * v.x + ( 1./sqrt(6.)) * v.y + (1./sqrt(3.)) * v.z,
			(-1./sqrt(2.)) * v.x + ( 1./sqrt(6.)) * v.y + (1./sqrt(3.)) * v.z,
			                       (-2./sqrt(6.)) * v.y + (1./sqrt(3.)) * v.z);

}
F64vec rotate_bwd(const F64vec v){
	return F64vec(
			( 1./sqrt(2.)) * v.x + (-1./sqrt(2.)) * v.y,
			( 1./sqrt(6.)) * v.x + ( 1./sqrt(6.)) * v.y + (-2./sqrt(6.)) * v.z,
			( 1./sqrt(3.)) * v.x + ( 1./sqrt(3.)) * v.y + ( 1./sqrt(3.)) * v.z);

}

int main(){
	srand48(20150415);

	F64vec v0(drand48(), drand48(), drand48());
	std::cout << "before  : " << v0 << std::endl;

	F64vec v1 = rotate_fwd(v0);
	std::cout << "after   : " << v1 << std::endl;

	F64vec v2 = rotate_bwd(v1);
	std::cout << "recover : " << v2 << std::endl;

	 return 0;
}


