#ifndef MPIFFT3DFFTW2_H
#define MPIFFT3DFFTW2_H
#ifdef USE_FFTW2
#ifdef USE_MPI

#include "MPIFFT3D.h"
#include <rfftw_mpi.h>

namespace PMEModule {

class MPIFFT3Dfftw2: public FFT3D {
    public:
	MPIFFT3Dfftw2(int nx, int ny, int nz,
	  int local_nx, int nx_start, int local_ny, int ny_start,
	  int local_sizes);
	~MPIFFT3Dfftw2();
	void forward();
	void inverse();
	void backward();
	static MPI_Comm getCommunicator() {
            if (pComm) return pComm->communicator;
	    return MPI_COMM_WORLD;
	}
	static FFT3D::Ptr createFFT3D(int nx, int ny, int nz) {
	    rfftwnd_mpi_plan plan;
	    int local_nx, local_x_start;
	    int local_ny_after_tranpose, local_y_start_after_transpose;
	    int total_local_sizes;
	    plan = rfftw3d_mpi_create_plan(getCommunicator(),
	      nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
	    rfftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
	      &local_ny_after_tranpose, &local_y_start_after_transpose,
	      &total_local_sizes);
	    rfftwnd_mpi_destroy_plan(plan);
#ifdef CPPMD_ENABLE_LARGEMODEL
            if(DebugLog::verbose>1){
              std::cout << "nx,ny,nz " << nx << " " << ny << " " << nz << std::endl;
              std::cout << "local_nx local_x_start total_local_sizes " << local_nx << " " << local_x_start << " " << total_local_sizes << std::endl;
            }
#endif

	    return FFT3D::Ptr(new MPIFFT3Dfftw2(nx, ny, nz,
	      local_nx, local_x_start,
	      local_ny_after_tranpose, local_y_start_after_transpose,
	      total_local_sizes));
	}
	DivideType getRealDivideType() {
            return DIVIDE_X;
        }
    private:
	fftw_complex* top;
	double normalizeFactor;
	MPI_Comm fftWorld;
	rfftwnd_mpi_plan forwardPlan;
	rfftwnd_mpi_plan inversePlan;
	std::vector<double> work;
};

MPIFFT3Dfftw2::MPIFFT3Dfftw2(int _nx, int _ny, int _nz,
    int local_nx, int nx_start, int local_ny, int ny_start, int local_sizes):
    FFT3D(local_nx, _ny, _nz, 0, local_sizes),
    top((fftw_complex *)pc3d->data()), normalizeFactor(1.0/(_nx*_ny*_nz)),
    fftWorld(getCommunicator()),
    forwardPlan(rfftw3d_mpi_create_plan(fftWorld,
	    _nx, _ny, _nz, FFTW_REAL_TO_COMPLEX,
	    FFTW_ESTIMATE | FFTW_IN_PLACE)),
    inversePlan(rfftw3d_mpi_create_plan(fftWorld,
	    _nx, _ny, _nz, FFTW_COMPLEX_TO_REAL,
	    FFTW_ESTIMATE | FFTW_IN_PLACE)),
    work()
{
    size[0] = _nx;
    size[1] = _ny;
    size[2] = _nz;
#ifdef CPPMD_ENABLE_LARGEMODEL
    if(DebugLog::verbose>1){
      for(int i=0;i<3;i++){
        std::cout << " size[" <<i << "] " << size[i] ;
      }
      std::cout << std::endl;
      std::cout << "local_nx local_ny local_sizes " << local_nx << " " << local_ny << " " << local_sizes << std::endl;
    }
#endif
    complexAxis[0] = 1;
    complexAxis[1] = 0;
    complexAxis[2] = 2;
    Dims rbases = {{nx_start, 0, 0}};
    pr3d->reindex(rbases);
    Dims cshape = {{local_ny, _nx, nz/2+1}};
    pc3d->reshape(cshape);
    Dims cbases = {{ny_start, 0, 0}};
    pc3d->reindex(cbases);
    work.assign(local_sizes, 0.0);
    //    std::cout << "call setStartEnd()" << std::endl;
    setStartEnd();
}

MPIFFT3Dfftw2::~MPIFFT3Dfftw2()
{
    int finalized;
    MPI_Finalized(&finalized);
    if (finalized) return;
    rfftwnd_mpi_destroy_plan(forwardPlan);
    rfftwnd_mpi_destroy_plan(inversePlan);
}

void MPIFFT3Dfftw2::forward()
{
    rfftwnd_mpi(forwardPlan, 1, pr3d->data(), &work[0], FFTW_TRANSPOSED_ORDER);
}

void MPIFFT3Dfftw2::inverse()
{
    rfftwnd_mpi(inversePlan, 1, pr3d->data(), &work[0], FFTW_TRANSPOSED_ORDER);
    scaling(normalizeFactor);
}
void MPIFFT3Dfftw2::backward()
{
    rfftwnd_mpi(inversePlan, 1, pr3d->data(), &work[0], FFTW_TRANSPOSED_ORDER);
}
}
#endif
#endif
#endif
