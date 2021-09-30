#include <vector>
#include <particle_simulator.hpp>

struct Cmpvec{
	PS::F64 PS::F64vec::*loc;
	Cmpvec(PS::F64 PS::F64vec::*_loc) : loc(_loc) {}

	bool operator()(const PS::F64vec &lhs, const PS::F64vec &rhs){
		return (lhs.*loc < rhs.*loc);
	}
};

int main(){
	std::vector<PS::F64vec> vec(10);

	for(int i=0, n=vec.size(); i<n; i++){
		vec[i].x = drand48();
		vec[i].y = drand48();
		vec[i].z = drand48();
	}

	std::sort(vec.begin(), vec.end(), Cmpvec(&PS::F64vec::x));
	for(int i=0, n=vec.size(); i<n; i++) std::cout << vec[i] << std::endl;
	std::cout << std::endl;

	std::sort(vec.begin(), vec.end(), Cmpvec(&PS::F64vec::y));
	for(int i=0, n=vec.size(); i<n; i++) std::cout << vec[i] << std::endl;
	std::cout << std::endl;

	std::sort(vec.begin(), vec.end(), Cmpvec(&PS::F64vec::z));
	for(int i=0, n=vec.size(); i<n; i++) std::cout << vec[i] << std::endl;
	std::cout << std::endl;

	return 0;
}
