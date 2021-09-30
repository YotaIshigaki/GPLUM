#include <stdio.h>
#include <unistd.h>
#include <mpi.h>

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);
	{
		int p, size, rank, namelen;
		char hostname[MPI_MAX_PROCESSOR_NAME];

		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Get_processor_name(hostname, &namelen);

		for(p=0; p<size; p++){
			if(p == rank){
				printf("Hello, my name is %s, rank %6d of %6d\n",
						hostname, rank, size);
				fflush(stdout);
			}
			usleep(1000);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	MPI_Finalize();
	return 0;
}
