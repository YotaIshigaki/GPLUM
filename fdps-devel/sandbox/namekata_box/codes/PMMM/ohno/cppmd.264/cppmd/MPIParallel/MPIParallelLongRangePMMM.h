//! MPIParallel operator for PMMM Long-Range Interaction
#ifndef MPIPARALLELLONGRANGEPMMM_H
#define MPIPARALLELLONGRANGEPMMM_H

#include <mpi.h>
#include <vector>
#include "CellIndex.h"
#include "LongRangeParameter.h"
#include "pmmm_fmm.h"

template <typename real_t=double, typename cplx_t = std::complex<real_t> >
class MPIMexchagePMMM {
 public:
 MPI_Comm mpi_mm_comm;
 int mm_id;
 int num_long;
 int total_cell;
 std::vector<int> m_boundary;
 std::vector<int> m_count;
 std::vector<int> m_displs;

 std::vector<cplx_t> send_buf;
 std::vector<cplx_t> buf;

 void set_comm(MPI_Comm _mm_comm){
   mpi_mm_comm = _mm_comm;
   MPI_Comm_size(mpi_mm_comm, &num_long);
   MPI_Comm_rank(mpi_mm_comm, &mm_id);
 }

 void set_size(int _total_cell){
   total_cell = _total_cell;
 }

 void set_m_boundary(const std::vector<int> &_m_boundary){
   if(_m_boundary.size()!=num_long+1){
     /// error
   }
   m_boundary.resize(_m_boundary.size());
   std::copy(_m_boundary.begin(),_m_boundary.end(),m_boundary.begin());
 }
 
 void set_buf_param(){
   m_count.resize(num_long);
   m_displs.resize(num_long+1);
   int l;
   m_displs[0] = 0;
   for(l=0;l<num_long;l++){
     m_count[l] = total_cell*(m_boundary[l+1]-m_boundary[l]);
     m_displs[l+1] = m_displs[l]+m_count[l];
   }
   send_buf.resize(m_count[mm_id]);
   buf.resize(total_cell*m_boundary[num_long]);
 }
 
 void exchangeM(std::vector< MultipoleMoment<cplx_t> > &mm){
   int n = 0;
   for(int m=m_boundary[mm_id];m<m_boundary[mm_id+1];m++){
     for(int i=0;i<mm.size();i++,n++){
       send_buf[n] = mm[i].buf[m];
     }
   }

   MPI_Allgatherv(&(send_buf[0]), m_count[mm_id], MPI_DOUBLE_COMPLEX,
		  &(buf[0]), &(m_count[0]), &(m_displs[0]), MPI_DOUBLE_COMPLEX, mpi_mm_comm);
   n = 0;
   for(int m=0;m<m_boundary[num_long];m++){
     for(int i=0;i<mm.size();i++,n++){
       mm[i].buf[m] = buf[n];
     }
   }
 }
 
};

template <typename real_t=double>
class MPIReceiverForLongPMMM {
 public:
 MPI_Comm pmmm_comm;
 int mm_id;
 int pmmm_id;
 std::vector<MPI_Request> send_request;
 std::vector<MPI_Request> recv_request;
 std::vector<CellMethodModule::CellRange> cell_range;
 std::vector<int> targetid;
 int num_target;
 int num_m;
 int cell_size[3];
 
 typedef std::vector<real_t> buf_t;
 std::vector< buf_t > buf;


 void set_size(int _cell_size[3], int _num_m){
   cell_size[0] = _cell_size[0];
   cell_size[1] = _cell_size[1];
   cell_size[2] = _cell_size[2];
   num_m = _num_m;
 }

 void set_comm(const MPI_Comm &_comm, const int &_mm_id){
   pmmm_comm = _comm;
   mm_id = _mm_id;
   MPI_Comm_rank(pmmm_comm, &pmmm_id);
 }

 void set_target(const std::vector<int> &_targetid, const std::vector<CellMethodModule::CellRange> &_cell_range){
   if(_cell_range.size()!=_targetid.size()){
     // error
   }
   num_target = _targetid.size();
   targetid.resize(num_target);
   std::copy(_targetid.begin(),_targetid.end(),targetid.begin());
   cell_range.resize(num_target);
   std::copy(_cell_range.begin(),_cell_range.end(),cell_range.begin());
   buf.resize(num_target);
   for(int t=0;t<num_target;t++){
     buf[t].resize((cell_range[t].max[0]-cell_range[t].min[0])
		   *(cell_range[t].max[1]-cell_range[t].min[1])
		   *(cell_range[t].max[2]-cell_range[t].min[2])
		   *num_m);
   }
   send_request.resize(num_target);
   recv_request.resize(num_target);
 }
 
 void recvM(){
   int t;
   for(t=0;t<num_target;t++){
     //     printf("MPI_Irecv %d double from %d tag %d\n",buf[t].size(),targetid[t],pmmm_id);
     MPI_Irecv(&(buf[t][0]),buf[t].size(),MPI_DOUBLE,
	       targetid[t], pmmm_id, pmmm_comm, &(recv_request[t]));
   }
 }
 void recvM_wait(std::vector< std::vector<real_t> > &mm){
   int result;
   int ti;
   //   printf("recvM_wait\n");
   for(ti=0;ti<num_target;ti++){
     int t;
     result = MPI_Waitany(num_target, &(recv_request[0]), &t, MPI_STATUS_IGNORE);
     if(result!=MPI_SUCCESS){
       printf("Fail MPI_Waitany %d\n",result);
     }else{
       //       printf("MPI_Waitany %d\n",t);
     }
     int m;
     int i,j,k,n;
     n=0;
     for(m=0;m<num_m;m++){
       for(k=cell_range[t].min[2];k<cell_range[t].max[2];k++){
	 for(j=cell_range[t].min[1];j<cell_range[t].max[1];j++){
	   for(i=cell_range[t].min[0];i<cell_range[t].max[0];i++){
	     mm[m][i+cell_size[0]*(j+cell_size[1]*(k))] = buf[t][n];
	     n++;
	   }
	 }
       }
     }
   }
 }

 void sendL(std::vector< std::vector<real_t> > &le){
   int t;
   for(t=0;t<num_target;t++){
     int m;
     int i,j,k,n;
     n=0;
     for(m=0;m<num_m;m++){
       for(k=cell_range[t].min[2];k<cell_range[t].max[2];k++){
	 for(j=cell_range[t].min[1];j<cell_range[t].max[1];j++){
	   for(i=cell_range[t].min[0];i<cell_range[t].max[0];i++){
	     buf[t][n] = le[m][i+cell_size[0]*(j+cell_size[1]*(k))];
	     n++;
	   }
	 }
       }
     }
   }
   for(t=0;t<num_target;t++){
     MPI_Isend(&(buf[t][0]),buf[t].size(),MPI_DOUBLE,
	       targetid[t], pmmm_id, pmmm_comm, &(send_request[t]));
   }
 }
 void sendL_wait(){
   int ti;
   for(ti=0;ti<num_target;ti++){
     int t;
     MPI_Waitany(num_target, &(send_request[0]), &t, MPI_STATUS_IGNORE);
   }
 }

 
};

template <typename real_t=double>
class MPISenderForLongPMMM {
 public:
 MPI_Comm pmmm_comm;
 std::vector<MPI_Request> send_request;
 std::vector<MPI_Request> recv_request;
 std::vector<int> targetid;
 int num_target;
 std::vector<int> m_boundary;
 int cell_size[3];
 int size_m;

 typedef std::vector<real_t> buf_t;
 std::vector< buf_t > buf;

 void set_size(const int local_cell_size[3], const int _size_m){
   cell_size[0] = local_cell_size[0];
   cell_size[1] = local_cell_size[1];
   cell_size[2] = local_cell_size[2];
   size_m = _size_m;
 }

 void set_comm(const MPI_Comm &_comm){
   pmmm_comm = _comm;
 }
 
 void set_long_target(const std::vector<int> &_targetid, const std::vector<int> &_m_boundary){
   if(_m_boundary.size()!=_targetid.size()+1){
     // error
   }
   num_target = _targetid.size();
   targetid.resize(num_target);
   std::copy(_targetid.begin(),_targetid.end(),targetid.begin());
   m_boundary.resize(_m_boundary.size());
   std::copy(_m_boundary.begin(),_m_boundary.end(),m_boundary.begin());
   buf.resize(num_target);
   for(int t=0;t<num_target;t++){
     int size_tm = m_boundary[t+1]-m_boundary[t];
     buf[t].resize(cell_size[0]*cell_size[1]*cell_size[2]*size_tm);
   }
   send_request.resize(num_target);
   recv_request.resize(num_target);
 }
 
 void sendM(std::vector<MultipoleMoment<> > &mm){
   int t;
#pragma omp parallel for
   for(t=0;t<num_target;t++){
     int m,mi;
     for(m=m_boundary[t],mi=0;m<m_boundary[t+1];m++,mi++){
       for(int k=0;k<cell_size[2];k++){
	 for(int j=0;j<cell_size[1];j++){
	   for(int i=0;i<cell_size[0];i++){
	     buf[t][i+(cell_size[0]*(j+cell_size[1]*(k+cell_size[2]*mi)))] = mm[(i+(cell_size[0]*(j+cell_size[1]*k)))].buf[m];
	   }
	 }
       }
     }
   }
   
   for(t=0;t<num_target;t++){
     //     printf("MPI_Isend %d double for %d tag %d\n",buf[t].size(),targetid[t],targetid[t]);
     MPI_Isend(&(buf[t][0]),buf[t].size(),MPI_DOUBLE,targetid[t],targetid[t],pmmm_comm,&(send_request[t]));
   }
 }
 void sendM_wait(){
   int ti;
   for(ti=0;ti<num_target;ti++){ 
     int t;
    MPI_Waitany(num_target, &(send_request[0]), &t, MPI_STATUS_IGNORE);
   }
 }

 void recvL(){
   int t;
   for(t=0;t<num_target;t++){
     MPI_Irecv(&(buf[t][0]),buf[t].size(),MPI_DOUBLE,targetid[t],targetid[t],pmmm_comm,&(recv_request[t]));
   }
 }
 void recvL_wait(std::vector<LocalExpansion<> > &le){
   int ti;
   for(ti=0;ti<num_target;ti++){
     int t;
     MPI_Waitany(num_target, &(recv_request[0]), &t, MPI_STATUS_IGNORE);
     int m,mi;
     for(m=m_boundary[t],mi=0;m<m_boundary[t+1];m++,mi++){
       for(int k=0;k<cell_size[2];k++){
	 for(int j=0;j<cell_size[1];j++){
	   for(int i=0;i<cell_size[0];i++){
	     le[(i+(cell_size[0]*(j+cell_size[1]*k)))].buf[m] = buf[t][i+(cell_size[0]*(j+cell_size[1]*(k+cell_size[2]*mi)))];
	   }
	 }
       }
     }
   }
 }
};


#endif
