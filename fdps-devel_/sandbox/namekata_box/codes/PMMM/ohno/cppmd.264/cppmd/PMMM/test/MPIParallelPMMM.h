#pragma once

#include <mpi.h>

#include <vector>
#include <algorithm>

#include "fmm.h"
#include "cell.h"

template <int p>
void mm_decompose(int m_decomp[(p+1)*(p+1)], int num_mm, int mm0)
{
  int size = (p+1)*(p+1);

  int nm = size/num_mm;
  int rm = size-nm*num_mm;
  int m;
  int t = mm0;
  for(m=0;m<size;){
    int n = nm;
    if(m<rm)n++;
    int nn = m+n;
    for(;m<nn;m++){
      m_decomp[m] = t; 
    }
    t++;
  }
}

enum MPIPMMM_color{
  PMMM_PP = 0,
  PMMM_MM = 1
};


template <int p, int NX, int NY, int NZ, typename real_t=double, typename cplx_t = std::complex<real_t> >
struct MM_MM_Comm{
  enum{
    LEN  = lmbuf<p>::length,
  };
  MPI::Intracomm comm;
  int mm_size;
  int mm_id;
  std::vector<int> mm_list;
  std::vector<int> mm_count;
  std::vector<int> mm_displs;
  
  void set_comm(const MPI::Intracomm &_comm)
  {
    comm = _comm;
    mm_id = comm.Get_rank();
    mm_size = comm.Get_size();
  }

  void set_mm_list(const std::vector<int> &_mm_list)
  {
    mm_list.resize(_mm_list.size());
    std::copy(_mm_list.begin(),_mm_list.end(),mm_list.begin());
  }

  void set_mm_count(const int mm_targets[LEN])
  {
    mm_count.clear();
    mm_count.resize(mm_size);
    mm_displs.clear();
    mm_displs.resize(mm_size);
    int mm;
    for(mm=0;mm<mm_size;mm++)mm_count[mm]=0;
    int size = NX*NY*(1+NX/2);
    int m;
    for(m=0;m<LEN;m++){
      mm_count[mm_targets[m]]+=size;
    }
    mm_displs[0] = 0;
    for(mm=1;mm<mm_size;mm++){
      mm_displs[mm] = mm_displs[mm-1]+mm_count[mm-1];
    }
  }

  void exchange_mm(cplx_t mm[NZ][NY][1+NX/2][LEN])
  {
    cplx_t sbuf[NZ*NY*(1+NX/2)*LEN];
    cplx_t buf[NZ*NY*(1+NX/2)*LEN];
    int size_cell = NZ*NY*(1+NX/2);
    for(int mi=0;mi<mm_list.size();mi++){
      int m = mm_list[mi];
      for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<(1+NX/2);i++){
	    sbuf[(i+((1+NX/2)*(j+NY*k)))+mi*size_cell] = mm[k][j][i][m];
	  }
    }
    {
      int m;
      for(m=0;m<mm_count.size();m++){
	printf(" %d:%d",mm_displs[m],mm_count[m]);
      }
      printf("\n");
    }
    {
      int m;
      for(m=0;m<LEN;m++){
	printf(" %e",real(mm[0][0][0][m]));
      }
      printf("\n");
    }
    {
      for(int mi=0;mi<mm_list.size();mi++){
	int m = mm_list[mi];
	printf(" %e",real(sbuf[mi*size_cell]));
      }
      printf("\n");
    }
    int m0 = mm_list[0];
    comm.Allgatherv(sbuf, size_cell*mm_list.size(), MPI::DOUBLE_COMPLEX, 
		    buf, &(mm_count[0]), &(mm_displs[0]), MPI::DOUBLE_COMPLEX);
    int mi=0;
    for(int m=0;m<LEN;m++){
      /*
      if(m==mm_list[mi]){
	mi++;
	continue;
      }
      */
      for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<(1+NX/2);i++){
	    mm[k][j][i][m] = buf[(i+((1+NX/2)*(j+NY*k)))+m*size_cell];
	  }
    }
    {
      int m;
      for(m=0;m<LEN;m++){
	printf(" %e",real(mm[0][0][0][m]));
      }
      printf("\n");
    }
  }
};

template <int p, int NX, int NY, int NZ, typename real_t=double>
struct PP_MM_Comm{
  int targets[(p+1)*(p+1)];
  CellRange cell_range;
  typedef std::vector<real_t> buf_t;
  buf_t buf[(p+1)*(p+1)];
  MPI::Request send_request[(p+1)*(p+1)];
  MPI::Request recv_request[(p+1)*(p+1)];

  void set_targets(int _targets[(p+1)*(p+1)]){
    int m;
    for(m=0;m<(p+1)*(p+1);m++){
      targets[m] = _targets[m];
    }
  }

  void set_cell_range(const CellRange _cell_range)
  {
    cell_range = _cell_range;
    int m;
    for(m=0;m<(p+1)*(p+1);m++)buf[m].resize(cell_range.size);
  }

  void send_M(const Cell_FMM<p> cell[NZ][NY][NX])
  {
    int m;
#pragma omp parallel for
    for(m=0;m<(p+1)*(p+1);m++){
      for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	    buf[m][i+(NX*(j+NY*k))] = cell[k][j][i].mm.buf[m];
	  }
    }
    for(m=0;m<(p+1)*(p+1);m++){
      send_request[m] = MPI::COMM_WORLD.Isend(&(buf[m][0]),cell_range.size,MPI::DOUBLE,targets[m],m);
    }
    //    MPI_Status status;
    for(m=0;m<(p+1)*(p+1);m++){
      send_request[m].Wait();
    }
  }

  void recv_L(Cell_FMM<p> cell[NZ][NY][NX])
  {
    int m;
    for(m=0;m<(p+1)*(p+1);m++){
      recv_request[m] = MPI::COMM_WORLD.Irecv(&(buf[m][0]),cell_range.size,MPI::DOUBLE,targets[m],m);
    }
    //    MPI_Status status;
    for(m=0;m<(p+1)*(p+1);m++){
      recv_request[m].Wait();
    }

#pragma omp parallel for
    for(m=0;m<(p+1)*(p+1);m++){
      for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	    cell[k][j][i].le.buf[m] = buf[m][i+(NX*(j+NY*k))];
	  }
    }
  }

};

template <int p, int NX, int NY, int NZ, typename real_t=double>
struct MM_PP_Comm{
  std::vector<int> targets;
  typedef std::vector<real_t> buf_t;
  std::vector< buf_t > buf;
  std::vector<MPI::Request> send_request;
  std::vector<MPI::Request> recv_request;
  std::vector<CellRange> cell_range;
  int num_target;


  void set_targets(const std::vector<int> _targets, const std::vector<CellRange> _cell_range){
    targets.resize(_targets.size());
    std::copy(_targets.begin(),_targets.end(),targets.begin());
    cell_range.resize(_cell_range.size());
    std::copy(_cell_range.begin(),_cell_range.end(),cell_range.begin());
    num_target = targets.size();
    buf.resize(num_target);
    int t;
    for(t=0;t<num_target;t++){
      buf[t].resize(cell_range[t].size);
    }
    send_request.resize(num_target);
    recv_request.resize(num_target);
  }


  void recv_M(MultipoleMoment<p> mm[NZ][NY][NX], const int m){
    int t;
    for(t=0;t<num_target;t++){
      recv_request[t] = MPI::COMM_WORLD.Irecv(&(buf[t][0]),buf[t].size(),MPI::DOUBLE,
					      targets[t], m);
    }
    //    MPI_Status status;
    for(t=0;t<num_target;t++){
      recv_request[t].Wait();
    }


    for(t=0;t<num_target;t++){
      int i,j,k, n;
      n=0;
      for(k=cell_range[t].min[2];k<cell_range[t].max[2];k++){
	for(j=cell_range[t].min[1];j<cell_range[t].max[1];j++){
	  for(i=cell_range[t].min[0];i<cell_range[t].max[0];i++){
	    mm[k][j][i].buf[m] = buf[t][n];
	    n++;
	  }
	}
      }
    }
  }
  
  void send_L(const LocalExpansion<p> le[NZ][NY][NX], const int m){
    int t;
    for(t=0;t<num_target;t++){
      int i,j,k, n;
      n=0;
      for(k=cell_range[t].min[2];k<cell_range[t].max[2];k++){
	for(j=cell_range[t].min[1];j<cell_range[t].max[1];j++){
	  for(i=cell_range[t].min[0];i<cell_range[t].max[0];i++){
	    buf[t][n] = le[k][j][i].buf[m];
	    n++;
	  }
	}
      }
    }


    for(t=0;t<num_target;t++){
      send_request[t] = MPI::COMM_WORLD.Isend(&(buf[t][0]),buf[t].size(),MPI::DOUBLE,
					      targets[t], m);
    }
    //    MPI_Status status;
    for(t=0;t<num_target;t++){
      send_request[t].Wait();
    }
  }
  
};



// for one - one

template <int p, int NX, int NY, int NZ>
  void send_M(Cell_FMM<p> cell[NZ][NY][NX], int target)
{
  MPI_Request send_request[NZ][NY][NX];
  int size = (p+1)*(p+1);

  for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	MPI_Isend(cell[k][j][i].mm.buf,size,MPI_DOUBLE,
		  target, i+(NX*(j+NY*k)), MPI_COMM_WORLD, &(send_request[k][j][i]));
      }

  MPI_Status status;
  for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	MPI_Wait(&(send_request[k][j][i]),&status);
      }
  
}
template <int p, int NX, int NY, int NZ>
  void recv_M(MultipoleMoment<p> mm[NZ][NY][NX], int target)
{
  MPI_Request recv_request[NZ][NY][NX];
  int size = (p+1)*(p+1);

  for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	MPI_Irecv(mm[k][j][i].buf,size,MPI_DOUBLE,
		  target, i+(NX*(j+NY*k)), MPI_COMM_WORLD, &(recv_request[k][j][i]));
      }

  MPI_Status status;
  for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	MPI_Wait(&(recv_request[k][j][i]),&status);
      }
  
}



template <int p, int NX, int NY, int NZ>
  void recv_L(Cell_FMM<p> cell[NZ][NY][NX], int target)
{
  MPI_Request recv_request[NZ][NY][NX];
  int size = (p+1)*(p+1);

  for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	MPI_Irecv(cell[k][j][i].le.buf,size,MPI_DOUBLE,
		  target, i+(NX*(j+NY*k)), MPI_COMM_WORLD, &(recv_request[k][j][i]));
      }

  MPI_Status status;
  for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	MPI_Wait(&(recv_request[k][j][i]),&status);
      }
  
}
template <int p, int NX, int NY, int NZ>
  void send_L(LocalExpansion<p> le[NZ][NY][NX], int target)
{
  MPI_Request send_request[NZ][NY][NX];
  int size = (p+1)*(p+1);

  for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	MPI_Isend(le[k][j][i].buf,size,MPI_DOUBLE,
		  target, i+(NX*(j+NY*k)), MPI_COMM_WORLD, &(send_request[k][j][i]));
      }

  MPI_Status status;
  for(int k=0;k<NZ;k++) for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
	MPI_Wait(&(send_request[k][j][i]),&status);
      }
  
}

