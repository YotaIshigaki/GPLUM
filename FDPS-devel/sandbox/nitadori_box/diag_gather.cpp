#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <algorithm>

struct DiagGather{
	MPI_Comm comm_rad;
	MPI_Comm comm_phi;
	MPI_Comm comm_short;
	MPI_Comm comm_long;

	int rank_glb, size_glb;
	int rank_rad, size_rad;
	int rank_phi, size_phi;
	int rank_short, size_short;
	int rank_long,  size_long;

	int icount;
	int *counts;
	int *displs;

	DiagGather(int nrad, int nphi, MPI_Comm &&comm_glb){
		int color, key;

		MPI_Comm_rank(comm_glb, &rank_glb);
		MPI_Comm_size(comm_glb, &size_glb);
		assert(size_glb == nrad*nphi);

		rank_rad = rank_glb % nrad;
		rank_phi = rank_glb / nrad;

		MPI_Comm_split(comm_glb, color=rank_phi, key=rank_rad, &comm_rad);
		MPI_Comm_split(comm_glb, color=rank_rad, key=rank_phi, &comm_phi);

		MPI_Comm_size(comm_rad, &size_rad);
		MPI_Comm_size(comm_phi, &size_phi);

		assert(nrad == size_rad);
		assert(nphi == size_phi);

		color = rank_phi / nrad;
		key   = rank_phi % nrad;
		if(nphi % nrad){
			int ncolor = nphi / nrad;
			if(rank_phi >= nrad * ncolor){
				color--;
				key += nrad;
			}
		}
		MPI_Comm_split(comm_phi, color, key, &comm_short);
		MPI_Comm_size(comm_short, &size_short);
		MPI_Comm_rank(comm_short, &rank_short);

		color = rank_phi % nrad;
		key   = rank_phi / nrad;

		MPI_Comm_split(comm_phi, color, key, &comm_long);
		MPI_Comm_size(comm_long, &size_long);
		MPI_Comm_rank(comm_long, &rank_long);

		counts = new int[nrad];
		displs = new int[nrad+1];

		icount = size_long;
		int root = rank_rad;
		MPI_Bcast(&icount, 1,      MPI_INT,    root, comm_short);
		MPI_Allgather(
				&icount, 1, MPI_INT,
				counts,  1, MPI_INT,
				comm_rad);
		displs[0] = 0;
		for(int i=0; i<nrad; i++){
			displs[i+1] = displs[i] + counts[i];
		}
	}

	~DiagGather(){
		delete [] counts;
		delete [] displs;
	}

	void show_comm_info(){
		printf("glb:(%2d, %2d), rad:(%2d, %2d), phi:(%2d, %2d), short:(%2d, %2d), long:(%2d, %2d)\n",
				rank_glb, size_glb,
				rank_rad, size_rad,
				rank_phi, size_phi,
				rank_short, size_short,
				rank_long, size_long);
	}

	void gather_cm_org(double val, double vals[]){
		int root;
		double tmp;
		MPI_Reduce(&val, &tmp, 1, MPI_DOUBLE, MPI_SUM, root=0, comm_rad);

		if(0 == rank_rad){
			MPI_Allgather(
					&tmp, 1, MPI_DOUBLE,
					vals, 1, MPI_DOUBLE,
					comm_phi);
		}
		MPI_Bcast(vals, size_phi, MPI_DOUBLE, root=0, comm_rad);
	}

	void gather_cm_fast(double val, double vals[], double work[]){
		int root = rank_phi % size_rad;
		double sum;
		MPI_Reduce(&val, &sum, 1, MPI_DOUBLE, MPI_SUM, root, comm_rad);

		if(root == rank_rad){
			MPI_Allgather(
					&sum, 1, MPI_DOUBLE,
					work, 1, MPI_DOUBLE,
					comm_long);
		}

		root = rank_rad;
#if 0
		int icount = size_long;
		MPI_Bcast(&icount, 1,      MPI_INT,    root, comm_short);
#endif
		MPI_Bcast(work,    icount, MPI_DOUBLE, root, comm_short);

		MPI_Allgatherv(
				work, icount,        MPI_DOUBLE,
				vals, counts, displs, MPI_DOUBLE,
				comm_rad);
	}

	void reorder_result(const double *src, double *dst) const {
		const int nrad = this->size_rad;
		for(int i=0; i<nrad; i++){
			int nj   = this->counts[i];
			int joff = this->displs[i];
			for(int j=0; j<nj; j++){
				dst[i + nrad*j] = src[joff + j];
			}
		}
	}
};

int main(int ac, char **av){
	MPI_Init(&ac, &av);

	int nrad = ac>1 ? atoi(av[1]) : 3;
	int nphi = ac>2 ? atoi(av[2]) : 31;

	DiagGather dg(nrad, nphi, MPI_COMM_WORLD);

	for(int i=0; i<dg.size_glb; i++){
		if(i == dg.rank_glb){
			dg.show_comm_info();
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	srand48(dg.rank_glb);
	for(int i=0; i<100; i++) (void)drand48();
	double val = drand48();
	std::vector<double> vals1(nphi);
	dg.gather_cm_org(val, vals1.data());

	for(int i=0; i<dg.size_glb; i++){
		if(i == dg.rank_glb){
			// std::sort(vals1.begin(), vals1.end());
			for(auto v : vals1){
				printf("%e ", v);
			}
			puts("");
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if(0 == dg.rank_glb){
		puts("");
	}
	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<double> vals2(nphi), work(nphi), vals3(nphi);;
	dg.gather_cm_fast(val, vals2.data(), work.data());
	dg.reorder_result(vals2.data(), vals3.data());

	for(int i=0; i<dg.size_glb; i++){
		if(i == dg.rank_glb){
			// std::sort(vals2.begin(), vals2.end());
			for(auto v : vals3){
				printf("%e ", v);
			}
			puts("");
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	for(int i=0; i<nphi; i++){
		assert(vals1[i] == vals3[i]);
	}

	MPI_Finalize();

	return 0;
}
