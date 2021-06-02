#include "FFT.h"
#ifdef USE_SSLIIMPI
extern "C" 
void ds_v3drcf_(double *x, int *kx1, int *kx2, int *kx3p, int *n1, int *n2, int *n3,
		double *w, int *kw2p, int *kw3, int *isin, int *isn, MPI_Fint *comm, int *icon);
#else
#include <rfftw_mpi.h>
#endif
#ifdef FFT_TIMER
#include "Timer.h"
#endif

#ifdef FIPP
# include <fj_tool/fipp.h>
#endif
#ifdef FAPP
# include <fj_tool/fapp.h>
#endif
#ifdef K_PA
# include "fjcoll.h"
#endif

/*
fft method used in PME.cc
*/

void print_r_single(int realStart[3], int realEnd[3], std::vector<double> realvec, int r3d[3])
{
  int x,y,z;
  for(x=realStart[0];x<realEnd[0];x++){
    std::cout << "x = " << x;
    std::cout << "  y = " <<  realStart[1] << ":" << realEnd[1];
    std::cout << "  z = " <<  realStart[2] << ":" << realEnd[2] << std::endl;
    for(y=realStart[1];y<realEnd[1];y++){
      for(z=realStart[2];z<realEnd[2];z++){
	int index = (z-realStart[2]) + r3d[2]*((y-realStart[1]) + r3d[1]*(x-realStart[0]));
	std::cout << realvec[index] << " ";
      }
      std::cout << std::endl;
    }
  }
}

void print_r(int realStart[3], int realEnd[3], std::vector<double> &realvec, int r3dx, int r3dy, int r3dz)
{
  int r3d[3] = {r3dx, r3dy, r3dz};
  int rank;
  int nnum;
  MPI_Status stat;
  MPI_Comm fftcomm = MPI_COMM_WORLD;

  MPI_Comm_rank(fftcomm,&rank);
  MPI_Comm_size(fftcomm,&nnum);
  MPI_Barrier(fftcomm);
  
  std::vector<double> buf;
  int num = realvec.size();

  if(rank==0){
    std::cout << "real" << std::endl;
    buf.resize(num);
  }

  std::cout.setf(std::ios::fixed);
  for(int r=0;r<nnum; r++){
    MPI_Barrier(fftcomm);
    if(rank==0){
      if(r==0){
	print_r_single(realStart,realEnd,realvec,r3d);
      }else{
	MPI_Recv((void *)(&(buf[0])),num,MPI_DOUBLE,r,0,fftcomm, &stat);
	print_r_single(realStart,realEnd,buf,r3d);
      }
    }else{
      if(r==rank){
	MPI_Send((void *)(&(realvec[0])),num,MPI_DOUBLE,0,0,fftcomm);
      }
    }
  }
  if(rank==0) std::cout << std::endl;
  MPI_Barrier(fftcomm);
}


void print_c_single(int complexStart[3], int complexEnd[3], std::vector<double> &realvec, int c3d[3])
{
  int x,y,z;
  for(x=complexStart[0];x<complexEnd[0];x++){
    std::cout << "x = " << x;
    std::cout << "  y = " <<  complexStart[1] << ":" << complexEnd[1];
    std::cout << "  z = " <<  complexStart[2] << ":" << complexEnd[2] << std::endl;
    for(y=complexStart[1];y<complexEnd[1];y++){
      for(z=complexStart[2];z<complexEnd[2];z++){
	int index = ((z-complexStart[2]) + c3d[2]*((y-complexStart[1]) + c3d[1]*(x-complexStart[0])))*2;
	//	std::cout << "(" << rvec[index] << "," << rvec[index+1] << ")";
	std::cout << sqrt(realvec[index]*realvec[index]+realvec[index+1]*realvec[index+1]) << " ";
      }
      std::cout << std::endl;
    }
  }
}

void print_c(int complexStart[3], int complexEnd[3], std::vector<double> realvec, int c3dx, int c3dy, int c3dz)
{
  int c3d[3] = {c3dx, c3dy, c3dz};
  int rank;
  int nnum;
  MPI_Status stat;
  MPI_Comm fftcomm = MPI_COMM_WORLD;

  MPI_Comm_rank(fftcomm,&rank);
  MPI_Comm_size(fftcomm,&nnum);
  MPI_Barrier(fftcomm);
  
  std::vector<double> buf;
  int num = realvec.size();

  if(rank==0){
    std::cout << "complex(amp)" << std::endl;
    buf.resize(num);
  }

  std::cout.setf(std::ios::fixed);
  for(int r=0;r<nnum; r++){
    MPI_Barrier(fftcomm);
    if(rank==0){
      if(r==0){
	print_c_single(complexStart,complexEnd,realvec,c3d);
      }else{
	MPI_Recv((void *)(&(buf[0])),num,MPI_DOUBLE,r,0,fftcomm, &stat);
	print_c_single(complexStart,complexEnd,buf,c3d);
      }
    }else{
      if(r==rank){
	MPI_Send((void *)(&(realvec[0])),num,MPI_DOUBLE,0,0,fftcomm);
      }
    }
  }
  if(rank==0) std::cout << std::endl;
  MPI_Barrier(fftcomm);
}



namespace FFTModule {

  inline void timer_start(int counter_id)
  {
#ifdef FFT_TIMER
    PerfCounter::start(counter_id);
#endif
  }
  inline void timer_stop(int counter_id)
  {
#ifdef FFT_TIMER
    PerfCounter::stop();
#endif
  }


// FFT3d.h 109 - 114
void FFT3D::assign(double value) {
  Real3D& r3d = *pr3d;
  for (size_t i=0;i<r3d.num_elements();++i) {
    r3d.data()[i] = value;
  }
}

Communicator* FFT3D::pComm = 0;

// for fftw2

class MPIFFT3Dfftw2: public FFT3D {
public:
  MPIFFT3Dfftw2(int nx, int ny, int nz,
		int local_nx, int nx_start, int local_ny, int ny_start,
		int local_sizes);

  ~MPIFFT3Dfftw2();
  void forward();
  void inverse();
  void backward();
// MPIFFT3Dfftw2.h 20 - 23
  static MPI_Comm getCommunicator() {
    if (pComm) return pComm->communicator;
    return MPI_COMM_WORLD;
  }

  static FFT3D::Ptr createFFT3D(int nx, int ny, int nz);
  
  FFT3D::DivideType getRealDivideType();

#ifdef USE_SSLIIMPI
#else
  fftw_complex* top;
#endif
  double normalizeFactor;
  MPI_Comm fftWorld;

#ifdef USE_SSLIIMPI
  MPI_Fint fftWorld_f;
  int kx1;
  int kx2;
  int kx3p;
  int kw2p;
  int kw3;
  int isin_forward;
  int isn_forward;
  int isin_backward;
  int isn_backward;
  int icon;
  int np;
#else
  rfftwnd_mpi_plan forwardPlan;
  rfftwnd_mpi_plan inversePlan;
#endif
  std::vector<double> work;
};


// MPIFFT3Dfftw2.h 107 - 110
void MPIFFT3Dfftw2::forward()
{
#ifdef TIMER_BARRIER
  MPI_Barrier(fftWorld);
#endif
  timer_start(timer_forward);
#ifdef FIPP
  fipp_start();
#endif
#ifdef FAPP
  fapp_start("forward",timer_forward,1);
#endif
#ifdef K_PA
  start_collection("forward");
#endif
#ifdef USE_SSLIIMPI
//  std::cout << "pr3d " << pr3d->shape()[0] << ", "  << pr3d->shape()[1] << ", " << pr3d->shape()[2]  << std::endl;
//  std::cout << " work.size() " << work.size() << std::endl;
//  std::cout << "ds_v3drcf_ " << kx1 << ", " << kx2 << ", " << kx3p << ", " 
//	    << nx << ", " << ny << ", " << nz  << ", " 
//	    << kw2p << ", " << kw3 <<   std::endl;
  ds_v3drcf_(pr3d->data(),&kx1,&kx2,&kx3p,&(size[2]),&(size[1]),&(size[0]),
	     &work[0],&kw2p,&kw3,&isin_forward,&isn_forward,&fftWorld_f,&icon);
//  std::cout << " done" << std::endl;
#else
  rfftwnd_mpi(forwardPlan, 1, pr3d->data(), &work[0], FFTW_TRANSPOSED_ORDER);
#endif
#ifdef K_PA
  stop_collection("forward");
#endif
#ifdef FAPP
  fapp_stop("forward",timer_forward,1);
#endif
#ifdef FIPP
  fipp_stop();
#endif
  timer_stop(timer_forward);
}

// MPIFFT3Dfftw2.h  117 - 120
void MPIFFT3Dfftw2::backward()
{
#ifdef TIMER_BARRIER
  MPI_Barrier(fftWorld);
#endif
  timer_start(timer_backward);
#ifdef FIPP
  fipp_start();
#endif
#ifdef FAPP
  fapp_start("backward",timer_forward,1);
#endif
#ifdef K_PA
  start_collection("backward");
#endif
#ifdef USE_SSLIIMPI
//  std::cout << "ds_v3drcf_" << std::endl;
  ds_v3drcf_(pr3d->data(),&kx1,&kx2,&kx3p,&(size[2]),&(size[1]),&(size[0]),
	     &work[0],&kw2p,&kw3,&isin_backward,&isn_backward,&fftWorld_f,&icon);
#else
  rfftwnd_mpi(inversePlan, 1, pr3d->data(), &work[0], FFTW_TRANSPOSED_ORDER);
#endif
#ifdef K_PA
  stop_collection("backward");
#endif
#ifdef FAPP
  fapp_stop("backward",timer_forward,1);
#endif
#ifdef FIPP
  fipp_stop();
#endif
  timer_stop(timer_backward);
}

// FFT3d.h 184 - 231
double FFT3D::convolution(const FFT3D::Real3D& BC)
{
  FFT3D::Complex3D& FQ = *pc3d;  // FQ is Fourie transform of Q
  //  print_r(realStart,realEnd,rvec,r3d_nx,r3d_ny,r3d_nz);
  forward();
  //  print_c(complexStart,complexEnd,rvec,pc3d->shape()[0],pc3d->shape()[1],pc3d->shape()[2]);
  int nzh = complexEnd[2]-1;
  double convnorm = 0.0;
  switch (getComplexStorageType()) {
  case FFT3D::FULL_Z:
    for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
      for (int my=complexStart[1];my<complexEnd[1];++my) {
	for (int mz=complexStart[2];mz<complexEnd[2]; ++mz) {
	  convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);  // Ulrich Essmann 1995 (4.7)
	  FQ[mx][my][mz] *= BC[mx][my][mz];
	}
      }
    }
    convnorm *= 0.5;
    break;
  case FFT3D::HALF_Z_EVEN:
    for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
      for (int my=complexStart[1];my<complexEnd[1];++my) {
	convnorm += 0.5*BC[mx][my][0]*norm(FQ[mx][my][0]);
	FQ[mx][my][0] *= BC[mx][my][0];
	convnorm += 0.5*BC[mx][my][nzh]*norm(FQ[mx][my][nzh]);
	FQ[mx][my][nzh] *= BC[mx][my][nzh];
	for (int mz=1;mz<nzh;++mz) {
	  convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);
	  FQ[mx][my][mz] *= BC[mx][my][mz];
	}
      }
    }
    break;
  case FFT3D::HALF_Z_ODD:
    for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
      for (int my=complexStart[1];my<complexEnd[1];++my) {
	convnorm += 0.5*BC[mx][my][0]*norm(FQ[mx][my][0]);
	FQ[mx][my][0] *= BC[mx][my][0];
	for (int mz=1;mz<complexEnd[2];++mz) {
	  convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);
	  FQ[mx][my][mz] *= BC[mx][my][mz];
	}
      }
    }
    break;
  }
  backward();
  return convnorm;
}

// FFT3D.h 232-289
double FFT3D::convolution(const FFT3D::Real3D& BC, 
		   const FFT3D::Real3D& VBC, double &v)
{
  FFT3D::Complex3D& FQ = *pc3d;  // FQ is Fourie transform of Q
  forward();
  int nzh = complexEnd[2]-1;
  double convnorm = 0.0;
  double virial = 0.0;
  switch (getComplexStorageType()) {
  case FFT3D::FULL_Z:
    for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
      for (int my=complexStart[1];my<complexEnd[1];++my) {
	for (int mz=complexStart[2];mz<complexEnd[2]; ++mz) {
	  convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);  // Ulrich Essmann 1995 (4.7)
	  virial += VBC[mx][my][mz]*norm(FQ[mx][my][mz]);
	  FQ[mx][my][mz] *= BC[mx][my][mz];
	}
      }
    }
    convnorm *= 0.5;
    virial *= 0.5;
    break;
  case FFT3D::HALF_Z_EVEN:
    for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
      for (int my=complexStart[1];my<complexEnd[1];++my) {
	convnorm += 0.5*BC[mx][my][0]*norm(FQ[mx][my][0]);
	virial += 0.5*VBC[mx][my][0]*norm(FQ[mx][my][0]);
	FQ[mx][my][0] *= BC[mx][my][0];
	convnorm += 0.5*BC[mx][my][nzh]*norm(FQ[mx][my][nzh]);
	virial += 0.5*VBC[mx][my][nzh]*norm(FQ[mx][my][nzh]);
	FQ[mx][my][nzh] *= BC[mx][my][nzh];
	for (int mz=1;mz<nzh;++mz) {
	  convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);
	  virial += VBC[mx][my][mz]*norm(FQ[mx][my][mz]);
	  FQ[mx][my][mz] *= BC[mx][my][mz];
	}
      }
    }
    break;
  case FFT3D::HALF_Z_ODD:
    for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
      for (int my=complexStart[1];my<complexEnd[1];++my) {
	convnorm += 0.5*BC[mx][my][0]*norm(FQ[mx][my][0]);
	virial += 0.5*VBC[mx][my][0]*norm(FQ[mx][my][0]);
	FQ[mx][my][0] *= BC[mx][my][0];
	for (int mz=1;mz<complexEnd[2];++mz) {
	  convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);
	  virial += VBC[mx][my][mz]*norm(FQ[mx][my][mz]);
	  FQ[mx][my][mz] *= BC[mx][my][mz];
	}
      }
    }
    break;
  }
  backward();
  v = -(convnorm + 2.0*virial);
  return convnorm;
}



// MPIFFT3Dfftw2.h 47-49
FFT3D::DivideType MPIFFT3Dfftw2::getRealDivideType() {
  return FFT3D::DIVIDE_X;
}




void FFT3D::setStartEnd() 
{
#ifdef CPPMD_ENABLE_LARGEMODEL
  if(DebugLog::verbose>1){
    std::cout << "setStartEnd()" << std::endl;
  }
#endif
  for (int i=0;i<3;++i) {
    realStart[i] = pr3d->index_bases()[i];
    realEnd[i] = realStart[i]+pr3d->shape()[i];
    complexStart[i] = pc3d->index_bases()[i];
    complexEnd[i] = complexStart[i]+pc3d->shape()[i];
  }
  if (getComplexStorageType() != FULL_Z) {
    realEnd[2] = std::min(realEnd[2], nz);
  }
  std::cout << "real " << realStart[0] << ":" << realEnd[0] 
	    << " " << realStart[1] << ":" << realEnd[1] 
	    << " " << realStart[2] << ":" << realEnd[2] << std::endl;
  std::cout << "complex " << complexStart[0] << ":" << complexEnd[0] 
	    << " " << complexStart[1] << ":" << complexEnd[1] 
	    << " " << complexStart[2] << ":" << complexEnd[2]  << std::endl; 
}

FFT3D::Ptr FFT3D::createFFT3D(int nx, int ny, int nz, bool single)
{
  return MPIFFT3Dfftw2::createFFT3D(nx, ny, nz);
}

// FFT3D.h 43 - 97
FFT3D::FFT3D(int _nx, int _ny, int _nz, int memtype,
	     int rvecSize, int cvecSize):
  nx(_nx), ny(_ny), nz(_nz),
  pr3d(), pc3d(), rvec(), cvec(),
  realStart(), realEnd(),
  complexStart(), complexEnd()
{
  size[0] = nx;
  size[1] = ny;
  size[2] = nz;
  r3d_nx = nx;
  r3d_ny = ny;
  r3d_nz = nz;
  complexAxis[0] = 0;
  complexAxis[1] = 1;
  complexAxis[2] = 2;
  switch (memtype) {
  case 0: /* in-place Real 3D FFT */
    if (rvecSize == 0) {
      rvecSize = nx*ny*(2*(nz/2+1));
    }
    rvec.assign(rvecSize, 0.0);
    r3d_nz = 2*(nz/2+1);
//    std::cout << "rvec.size() " << rvec.size() << std::endl;
//    std::cout << "createReal3D(&rvec[0], " << nx << ", " << ny << ", " << 2*(nz/2+1) << " )" << std::endl;
    pr3d = createReal3D(&rvec[0], nx, ny, 2*(nz/2+1));
    pc3d = createComplex3D(pr3d->data(), nx, ny, nz/2+1);
    break;
  case 1: /* out-of-place Complex 3D FFT */
    if (rvecSize == 0) {
      rvecSize = nx*ny*nz;
    }
    if (cvecSize == 0) {
      cvecSize = nx*ny*nz;
    }
    rvec.assign(rvecSize, 0.0);
    cvec.assign(cvecSize, 0.0);
    pr3d = createReal3D(&rvec[0], nx, ny, nz);
    pc3d = createComplex3D(&cvec[0], nx, ny, nz);
    break;
  case 2: /* out-of-place Real 3D FFT */
    if (rvecSize == 0) {
      rvecSize = nx*ny*nz;
    }
    if (cvecSize == 0) {
      cvecSize = nx*ny*(nz/2+1);
    }
    rvec.assign(rvecSize, 0.0);
    cvec.assign(cvecSize, 0.0);
    pr3d = createReal3D(&rvec[0], nx, ny, nz);
    pc3d = createComplex3D(&cvec[0], nx, ny, nz/2+1);
    break;
  default: 
    throw(std::runtime_error("FFT3D::FFT3D() invalid memtype"));
  }
  setStartEnd();

#ifdef FFT_TIMER
  timer_forward = PerfCounter::add_target("FFT_forward");
  timer_backward = PerfCounter::add_target("FFT_backward");
#endif

}

// MPIFFT3Dfftw2.h 59 - 96
MPIFFT3Dfftw2::MPIFFT3Dfftw2(int _nx, int _ny, int _nz,
    int local_nx, int nx_start, int local_ny, int ny_start, int local_sizes):
    FFT3D(local_nx, _ny, _nz, 0, local_sizes),
#ifdef USE_SSLIIMPI
#else
    top((fftw_complex *)pc3d->data()), 
#endif
    normalizeFactor(1.0/(_nx*_ny*_nz)),
    fftWorld(getCommunicator()),
#ifdef USE_SSLIIMPI
#else
    forwardPlan(rfftw3d_mpi_create_plan(fftWorld,
					_nx, _ny, _nz, FFTW_REAL_TO_COMPLEX,
					FFTW_ESTIMATE | FFTW_IN_PLACE)),
    inversePlan(rfftw3d_mpi_create_plan(fftWorld,
					_nx, _ny, _nz, FFTW_COMPLEX_TO_REAL,
					FFTW_ESTIMATE | FFTW_IN_PLACE)),
#endif
    work()
{
    size[0] = _nx;
    size[1] = _ny;
    size[2] = _nz;
#ifdef USE_SSLIIMPI
    isin_forward = 1;
    isn_forward = 1;
    isin_backward = 1;
    isn_backward = -1;
    fftWorld_f = MPI_Comm_c2f(fftWorld);
    MPI_Comm_size(fftWorld,&np);
    kx1 = r3d_nz;
    kx2 = local_ny*np;
    kx3p = local_nx;
    kw2p = local_ny;
    kw3 = kx3p*np;
    r3d_ny = kx2;
    rvec.resize(kx3p*kx2*kx1);
    /*
    Dims rshape = {{kx3p,kx2,kx1}};
    pc3d->reshape(rshape);
    */
    delete pc3d;
    delete pr3d;
    pr3d = createReal3D(&rvec[0], kx3p,kx2,kx1);
    pc3d = createComplex3D(pr3d->data(), kx3p,kx2,kx1/2);
#endif
#ifdef CPPMD_ENABLE_LARGEMODEL
    if(DebugLog::verbose>1){
      for(int i=0;i<3;i++){
        std::cout << " size[" <<i << "] " << size[i] ;
      }
      std::cout << std::endl;
      std::cout << "local_nx local_ny local_sizes " << local_nx << " " << local_ny << " " << local_sizes << std::endl;
    }
#endif
    Dims rbases = {{nx_start, 0, 0}};
    pr3d->reindex(rbases);
#ifdef USE_SSLIIMPI
    std::cout << "kx1 kx2 kx3p  kw2p kw3 " << kx1 << " " << kx2  << " " << kx3p << "  " <<   kw2p  << " " << kw3 << std::endl;
    Dims cshape = {{kx3p,kx2, nz/2+1}};
    pc3d->reshape(cshape);
    Dims cbases = {{nx_start, 0, 0}};
    pc3d->reindex(cbases);
    work.assign(kx1*kw2p*kw3, 0.0);
#else
    complexAxis[0] = 1;
    complexAxis[1] = 0;
    complexAxis[2] = 2;
    Dims cshape = {{local_ny, _nx, nz/2+1}};
    pc3d->reshape(cshape);
    Dims cbases = {{ny_start, 0, 0}};
    pc3d->reindex(cbases);
    work.assign(local_sizes, 0.0);
#endif
    //    std::cout << "call setStartEnd()" << std::endl;
    setStartEnd();
}


// MPIFFT3Dfftw2.h 98 - 105
MPIFFT3Dfftw2::~MPIFFT3Dfftw2()
{
    int finalized;
    MPI_Finalized(&finalized);
    if (finalized) return;
#ifdef USE_SSLIIMPI
#else
    rfftwnd_mpi_destroy_plan(forwardPlan);
    rfftwnd_mpi_destroy_plan(inversePlan);
#endif
}

// MPIFFT3Dfftw2.h 24 - 46
FFT3D::Ptr MPIFFT3Dfftw2::createFFT3D(int nx, int ny, int nz) 
{
#ifdef USE_SSLIIMPI
#else
  rfftwnd_mpi_plan plan;
#endif
  int local_nx, local_x_start;
  int local_ny_after_tranpose, local_y_start_after_transpose;
  int total_local_sizes;
#ifdef USE_SSLIIMPI
  int pn, p;
  MPI_Comm_size(getCommunicator(),&pn);
  MPI_Comm_rank(getCommunicator(),&p);
  local_nx = (nx+pn-1)/pn;
  local_x_start = local_nx*p;
  local_ny_after_tranpose = (ny+pn-1)/pn;
  total_local_sizes = local_nx*pn*local_ny_after_tranpose*2*(nz/2+1);
#else
  plan = rfftw3d_mpi_create_plan(getCommunicator(),
				 nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  rfftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
			  &local_ny_after_tranpose, &local_y_start_after_transpose,
			  &total_local_sizes);
  rfftwnd_mpi_destroy_plan(plan);
#endif
//  std::cout << "local_nx local_x_start total_local_sizes " << local_nx << " " << local_x_start << " " << total_local_sizes << std::endl;
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

}

/*
// FFT3D.h 347 - 356
static void registerLibrary(const std::string& name, bool single,
			    FFT3D::Ptr (*func)(int nx, int ny, int nz),
			    int (*funcA)(int n)=0) {
  LibraryNameFunc lib;
  lib.libraryName = name;
  lib.func = func;
  lib.funcA = funcA;
  lib.single = single;
  libraries.push_back(lib);
}

// FFT3D.cc 32 - 61
void initialize()
{
    if (initialized) return;
#ifdef USE_MPI
#ifdef USE_FFTW2
    registerLibrary("fftw2", false, MPIFFT3Dfftw2::createFFT3D);
#endif
#ifdef USE_FFTE
    registerLibrary("ffte", false, MPIFFT3Dffte::createFFT3D,
	MPIFFT3Dffte::getAllowableLength);
#endif
#ifdef USE_MKL
    registerLibrary("mkl", false, MPIFFT3Dmkl::createFFT3D);
#endif
#endif // USE_MPI
#ifdef USE_FFTW2
    registerLibrary("fftw2", true, FFT3Dfftw2::createFFT3D);
#endif
#ifdef USE_FFTW3
    registerLibrary("fftw3", true, FFT3Dfftw3::createFFT3D);
#endif
#ifdef USE_FFTE
    registerLibrary("ffte", true, FFT3Dffte::createFFT3D,
	FFT3Dffte::getAllowableLength);
#endif
#ifdef USE_MKL
    registerLibrary("mkl", true, FFT3Dmkl::createFFT3D);
#endif
    initialized = true;
}

//FFT3D.h 156
//static Ptr createFFT3D(int nx, int ny, int nz, bool single=false);
//FFT3D.cc 75-84
FFT3D::Ptr FFT3D::createFFT3D(int nx, int ny, int nz, bool single)
{
  initialize();
  for (std::vector<LibraryNameFunc>::const_iterator it=libraries.begin();
       it != libraries.end();++it) {
    if (single && !it->single) continue;
    return it->func(nx, ny, nz);
  }
  throw std::runtime_error("no FFT libraries");
}
*/






/*
Array3D used in PME

Real3D_Ptr(FFT3D::R3D::createArray(q[t],nx,ny,nz))
Real3D_Ptr()
typedef FFT3D::Real3D* Real3D_Ptr;; // FFT3D.h 26
Real3D()
typedef R3D::ArrayRef Real3D;   // FFT3D.h 20
typedef Array3D<double> R3D; // FFT3D.h 18
struct Array3D { ... }   // Array3D.h 32-52
typedef boost::multi_array_ref<T, 3> ArrayRef; // Array3D.h 35

Real3D_Ptr reindex(base)
boost::multi_array_ref<T, 3>::reindex


R3D::createArray(q[t],nx,ny,nz)
// Array3D.h 49-51
  static ArrayRef* createArray(T* data, int nx, int ny, int nz) {
    return new ArrayRef(data, boost::extents[nx][ny][nz]);
  }

Dims()
  typedef R3D::Dims Dims; // FFT3D.h 22
// HAVE_BOOST_ARRAY_HPP
typedef boost::array<typename ArrayRef::index, 3> Dims;  // Array3D.h 43


*/
