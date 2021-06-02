#include<algorithm>
#include<random>
#include<iostream>
#include"mpi.h"

struct Pos{
    double x, y, z;
};

int main(){

    int n_max = 1<<22;
    int n_min = 1<<10;

    std::mt19937 mt;
    std::uniform_real_distribution<double> gen(-10.0, 10.0);

    Pos * pos_org = new Pos[n_max];
    for(int i=0; i<n_max; i++){
        pos_org[i].x = gen(mt);
        pos_org[i].y = gen(mt);
        pos_org[i].z = gen(mt);
    }
    Pos * pos = new Pos[n_max];
    for(int n_tmp=n_min; n_tmp<=n_max; n_tmp*=2){
        const int n = n_tmp;
        for(int i=0; i<n; i++){
            pos[i] = pos_org[i];
        }
        const double t0 = MPI_Wtime();
        std::sort(pos, pos+n, [](const Pos & r, const Pos & l){return r.x < l.x;} );
        const double t = MPI_Wtime() - t0;
        std::cout<<"n= "<<n<<" t= "<<t<<" n*log2(n)= "<<n*log2(n)<<" t/(n*log2(n))= "<<t/(n*log2(n))<<std::endl;
    }
    return 0;
}
