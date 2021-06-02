#include <cassert>

template <int N>
struct NaCl{
	enum{ natom = N*N*N};
	struct Atom{
		double x, y, z, q;
	}atoms[natom];

	NaCl(){
		assert(0 == N%2);
		const double dx = 1.0 / N;
		Atom *a = atoms;
		for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++){
			a->x = (0.5 + i) * dx;
			a->y = (0.5 + j) * dx;
			a->z = (0.5 + k) * dx;
			a->q = (i+j+k)%2 ? 1.0 : -1.0;
			a->q /= (N*N);
			a++;
	   	}
	}

	double madelung_energy() const{
		const double M = 1.747558;
		return M;
	}
};
