#include <iostream>
#include <fftw3.h>
#include <mpi.h>

using S32 = int;
using U32 = unsigned int;

void sub(const S32 p,
         const S32 p_spj2mm,
         const U32 flag = FFTW_MEASURE,
         const bool use_mpifft_if_possible = true,
         const S32 fft_size_crit = 10000) {
    std::cout << "flag = " << flag << std::endl;
}


int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    const U32 flag1 = FFTW_MEASURE;
    const U32 flag2 = 0U;
    
    std::cout << "flag1 = " << flag1 << std::endl;
    std::cout << "flag2 = " << flag2 << std::endl;

    const S32 p = 7;
    const S32 p_spj2mm = 5;
    sub(p, p_spj2mm);

    MPI_Finalize();

    return 0;
}
