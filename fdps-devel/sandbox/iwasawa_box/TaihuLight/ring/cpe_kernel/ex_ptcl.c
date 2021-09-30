/* Standard C headers */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/* Headers specific to Sunway TaihuLight */
#include "slave.h"
#include "dma.h"
#include "cpe_func.h"

/* Local macros */
//#define CHECK_FIND_DEST_RANK
//#define CHECK_LVL_FIND_DEST_RANK (3)
//#define SPECIFY_CPE_IN_FIND_DEST_RANK
//#define USE_DEBUG_FLAG_FROM_MPE
//#define CHECK_MAKE_SEND_BUFS
#define RANK_CHECK (4)
#define ID_CHECK   (0)
#define CHUNK_ID_CHECK (31488)
#define PTCL_BUF_ID_CHECK (39)

extern int MY_RANK_MPI;
extern char HOST_NAME[170000][64];

extern double __ieee754_atan2 __P((double, double, double *, double *, double * ));

__thread_local int my_id, my_row_id, my_col_id;
__thread_local volatile dma_desc dma_get, dma_put;
__thread_local volatile int reply_get = 0;
__thread_local volatile int reply_put = 0;

__thread_local volatile int previous_rank;
__thread_local volatile int have_pos_domain_prev = 0;
__thread_local volatile F64ort_ pos_domain_prev;

__thread_local volatile int debug_flag = 0;


/* CPE communication functions */
static inline void cpe_bcast_int32(const int root_cpe_id, int * data) {
    int root_col_id, root_row_id;
    root_row_id = root_cpe_id/8;
    root_col_id = root_cpe_id%8;
    if (my_id == root_cpe_id) {
        REG_PUTR(*data,8);
        REG_PUTC(*data,8);
    }
    else {
        if (my_row_id == root_row_id) {
            REG_GETR(*data);
            REG_PUTC(*data,8);
        }
        else {
            REG_GETC(*data);
        }
    }
}

static inline void cpe_bcast_uint64(const int root_cpe_id, unsigned long * data) {
    int root_col_id, root_row_id;
    root_row_id = root_cpe_id/8;
    root_col_id = root_cpe_id%8;
    if (my_id == root_cpe_id) {
        REG_PUTR(*data,8);
        REG_PUTC(*data,8);
    }
    else {
        if (my_row_id == root_row_id) {
            REG_GETR(*data);
            REG_PUTC(*data,8);
        }
        else {
            REG_GETC(*data);
        }
    }
}

typedef int sumtype;
void prefix_sum(sumtype val, sumtype *beg, sumtype *end, sumtype *psum){
	const int myid = athread_get_id(-1);
	const int jcol = myid % 8;
	const int irow = myid / 8;
	sumtype sbuf, rbuf, sum;
	sumtype v0 = val;

	sync_array_();
	// __shfl_up(1) 
	{
		sbuf = val;
		REG_PUTR(sbuf, (jcol+1)%8);
		REG_GETR(rbuf);
		if(jcol < 1){
			sbuf = rbuf;
			REG_PUTC(sbuf, (irow+1)%8);
			REG_GETC(rbuf);
		}
		if(myid >= 1){
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(2) 
	{
		sbuf = val;
		REG_PUTR(sbuf, (jcol+2)%8);
		REG_GETR(rbuf);
		if(jcol < 2){
			sbuf = rbuf;
			REG_PUTC(sbuf, (irow+1)%8);
			REG_GETC(rbuf);
		}
		if(myid >= 2){
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(4) 
	{
		sbuf = val;
		REG_PUTR(sbuf, (jcol+4)%8);
		REG_GETR(rbuf);
		if(jcol < 4){
			sbuf = rbuf;
			REG_PUTC(sbuf, (irow+1)%8);
			REG_GETC(rbuf);
		}
		if(myid >= 4){
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(8) 
	{
		if(myid < 56){
			sbuf = val;
			REG_PUTC(sbuf, (irow+1)%8);
		}
		if(myid >= 8){
			REG_GETC(rbuf);
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(16) 
	{
		if(myid < 48){
			sbuf = val;
			REG_PUTC(sbuf, (irow+2)%8);
		}
		if(myid >= 16){
			REG_GETC(rbuf);
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(32) 
	{
		if(myid < 32){
			sbuf = val;
			REG_PUTC(sbuf, (irow+4)%8);
		}
		if(myid >= 32){
			REG_GETC(rbuf);
			val += rbuf;
		}
	}
	sync_array_();

	// bcast the sum
	if(7 == irow){
		if(7 == jcol){
			sum = val;
			REG_PUTR(sum, 8);
		}else{
			REG_GETR(sum);
		}
		REG_PUTC(sum, 8);
	}else{
		REG_GETC(sum);
	}
	sync_array_();

	*beg  = val-v0;
	*end  = val;
	*psum = sum;
}

static inline void IsNotInDomain(F64ort_ * pos_domain, F64vec_ * pos_ptcl, int * flag){
    F64_ x = pos_ptcl->x;
    F64_ y = pos_ptcl->y;
    F64_ z = pos_ptcl->z;
    *flag = 0;
    if(x < (pos_domain->low_).x || (pos_domain->high_).x <= x
       || y < (pos_domain->low_).y || (pos_domain->high_).y <= y
       || z < (pos_domain->low_).z || (pos_domain->high_).z <= z){
        *flag = 1;
    }
};

static inline void IsInDomain(F64ort_ * pos_domain, F64vec_ * pos_ptcl, int * flag){
    F64_ x = pos_ptcl->x;
    F64_ y = pos_ptcl->y;
    F64_ z = pos_ptcl->z;
    *flag = 0;
    if(x >= (pos_domain->low_).x && (pos_domain->high_).x > x
       || y >= (pos_domain->low_).y && (pos_domain->high_).y > y
       || z >= (pos_domain->low_).z && (pos_domain->high_).z > z){
        *flag = 1;
    }
}

void CheckIsInDomain(void * args){
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    unsigned long largs[5];
    reply_get = 0;
    dma_set_op(&dma_get,    DMA_GET);
    dma_set_mode(&dma_get,  PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get,  sizeof(largs));
    dma(dma_get, (long*)(args), (long*)(largs));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}

    int n_ptcl_ = (int)(largs[0]);
    void * adr_pos_my_domain_ = (void *)(largs[1]);
    F64ort_ pos_my_domain = *(F64ort_ *)adr_pos_my_domain_;
    /*
    if(my_id == 0){
        printf("a) n_ptcl_= %d \n", n_ptcl_);
        printf("A) pos_my_domain= %e %e %e %e %e %e \n",
               (pos_my_domain).low_.x,  (pos_my_domain).low_.y,  (pos_my_domain).low_.z,
               (pos_my_domain).high_.x, (pos_my_domain).high_.y, (pos_my_domain).high_.z);
    } 
    */
    void * adr_ptcl_ = (void *)(largs[2]);
    void * adr_adr_exchange_ = (void *)(largs[3]);
    void * adr_n_exchange_   = (void *)(largs[4]);    
    enum {
        CHUNK_SIZE = 64,
        N_ADR_EXCHANGE_LIMIT = 5000,
    };
    int n_loc = n_ptcl_/NUMBER_OF_CPE + ( (my_id < n_ptcl_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_ptcl_/NUMBER_OF_CPE)*my_id + ( (my_id < n_ptcl_ % NUMBER_OF_CPE) ? my_id : n_ptcl_ % NUMBER_OF_CPE );
    volatile fpLM ptcl[CHUNK_SIZE];
    volatile int adr_exchange[N_ADR_EXCHANGE_LIMIT];
    volatile int n_exchange[NUMBER_OF_CPE];
    //int n_diff_exchange[NUMBER_OF_CPE+1];

    double MY_PI = 3.14159265358979323846;
    double MY_PI_2 = 2.0*MY_PI;
    double * atanhi = (double *)ldm_malloc( sizeof(double)*4 );
    double * atanlo = (double *)ldm_malloc( sizeof(double)*4 );
    double * aT = (double *)ldm_malloc( sizeof(double)*11 );
    atanhi[0] = 4.63647609000806093515e-01; // atan(0.5)hi 0x3FDDAC67, 0x0561BB4F
    atanhi[1] = 7.85398163397448278999e-01; // atan(1.0)hi 0x3FE921FB, 0x54442D18
    atanhi[2] = 9.82793723247329054082e-01; // atan(1.5)hi 0x3FEF730B, 0xD281F69B
    atanhi[3] = 1.57079632679489655800e+00; // atan(inf)hi 0x3FF921FB, 0x54442D18
    atanlo[0] = 2.26987774529616870924e-17; // atan(0.5)lo 0x3C7A2B7F, 0x222F65E2
    atanlo[1] = 3.06161699786838301793e-17; // atan(1.0)lo 0x3C81A626, 0x33145C07
    atanlo[2] = 1.39033110312309984516e-17; // atan(1.5)lo 0x3C700788, 0x7AF0CBBD
    atanlo[3] = 6.12323399573676603587e-17; // atan(inf)lo 0x3C91A626, 0x33145C07
    aT[0] = 3.33333333333329318027e-01; // 0x3FD55555, 0x5555550D
    aT[1] = -1.99999999998764832476e-01; // 0xBFC99999, 0x9998EBC4
    aT[2] = 1.42857142725034663711e-01; // 0x3FC24924, 0x920083FF
    aT[3] = -1.11111104054623557880e-01; // 0xBFBC71C6, 0xFE231671
    aT[4] = 9.09088713343650656196e-02; // 0x3FB745CD, 0xC54C206E
    aT[5] = -7.69187620504482999495e-02; // 0xBFB3B0F2, 0xAF749A6D
    aT[6] = 6.66107313738753120669e-02; // 0x3FB10D66, 0xA0D03D51
    aT[7] = -5.83357013379057348645e-02; // 0xBFADDE2D, 0x52DEFD9A
    aT[8] = 4.97687799461593236017e-02; // 0x3FA97B4B, 0x24760DEB
    aT[9] = -3.65315727442169155270e-02; // 0xBFA2B444, 0x2C6A6C2F
    aT[10] = 1.62858201153657823623e-02; // 0x3F90AD3A, 0xE322DA11

    
    int i,j,k;
    int n_cnt = 0;
    for(i=0; i<n_loc; i += CHUNK_SIZE){
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get ptcl
        reply_get = 0;
        dma_set_op(&dma_get,    DMA_GET);
        dma_set_mode(&dma_get,  PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get,  sizeof(fpLM)*nn);
        dma(dma_get, (long*)(((fpLM *)adr_ptcl_) + my_offset + i), (long*)(ptcl));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for(k=0; k<nn; k++){
            assert(n_cnt < N_ADR_EXCHANGE_LIMIT-1);
            int flag = 0;
            const double x = ptcl[k].pos.x;
            const double y = ptcl[k].pos.y;
            const double z = ptcl[k].pos.z;
            double phi = __ieee754_atan2(y, x, atanhi, atanlo, aT);
            if(phi < 0.0) phi += MY_PI_2;
            double r = sqrt(x*x+y*y);
            F64vec_ pos;
            pos.x = phi;
            pos.y = r;
            pos.z = z;
            IsNotInDomain(&pos_my_domain, &pos, &flag);
            if(flag==1){
                adr_exchange[n_cnt++] = my_offset + i + k;
            }
        }
    }
    sync_array_();
    n_exchange[my_id] = n_cnt;
    sync_array_();
    /*
    for(i=0;i<my_id*1000000;i++) NOP();
    printf("my_id= %d, n_exchange[0]=%d, , n_exchange[1]=%d \n",
           my_id, n_exchange[0], n_exchange[1]);
    */
    int n_exchange_total = 0;
    int n_diff_exchange = 0;
    int id=0;
    for(id=0; id<NUMBER_OF_CPE; id++){
        sync_array_();
        cpe_bcast_int32(id, &n_exchange[id]);
        sync_array_();
        n_exchange_total += n_exchange[id];
        if(id < my_id){
            n_diff_exchange += n_exchange[id];
            //n_diff_exchange[id] += n_exchange[id];
        }
    }
    assert(n_exchange[my_id] < N_ADR_EXCHANGE_LIMIT);

    //** Put adr
    if(n_exchange[my_id] > 0){
        reply_put = 0;
        dma_set_op(&dma_put,    DMA_PUT);
        dma_set_mode(&dma_put,  PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put,  sizeof(int)*n_exchange[my_id]);
        dma(dma_put, (long*)((int *)adr_adr_exchange_ + n_diff_exchange), (long*)(adr_exchange));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }

    if(my_id == 0){
        //** Put adr
        reply_put = 0;
        dma_set_op(&dma_put,    DMA_PUT);
        dma_set_mode(&dma_put,  PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put,  sizeof(int));
        dma(dma_put, (long*)((int *)adr_n_exchange_), (long*)(&n_exchange_total));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
}


void GetCylCoord(void * args){
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    int n_ptcl_ = (int)((unsigned long*)args)[0];
    void * adr_len_peri_x_ = (void*)((unsigned long*)args)[1];
    double len_peri_x_ = *(double*)adr_len_peri_x_;
    void * adr_ptcl_ = (void *)((unsigned long*)args)[2];
    void * adr_pos_phi_ = (void *)((unsigned long*)args)[3];
    void * adr_pos_r_ = (void *)((unsigned long*)args)[4];

    enum {
        CHUNK_SIZE = 32,
    };
    int n_loc = n_ptcl_/NUMBER_OF_CPE + ( (my_id < n_ptcl_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_ptcl_/NUMBER_OF_CPE)*my_id + ( (my_id < n_ptcl_ % NUMBER_OF_CPE) ? my_id : n_ptcl_ % NUMBER_OF_CPE );
    volatile fpLM ptcl[CHUNK_SIZE];
    volatile double pos_phi[CHUNK_SIZE];
    volatile double pos_r[CHUNK_SIZE];
    double MY_PI = 3.14159265358979323846;
    double MY_PI_2 = 2.0*MY_PI;
    double * atanhi = (double *)ldm_malloc( sizeof(double)*4 );
    double * atanlo = (double *)ldm_malloc( sizeof(double)*4 );
    double * aT = (double *)ldm_malloc( sizeof(double)*11 );
    atanhi[0] = 4.63647609000806093515e-01; // atan(0.5)hi 0x3FDDAC67, 0x0561BB4F
    atanhi[1] = 7.85398163397448278999e-01; // atan(1.0)hi 0x3FE921FB, 0x54442D18
    atanhi[2] = 9.82793723247329054082e-01; // atan(1.5)hi 0x3FEF730B, 0xD281F69B
    atanhi[3] = 1.57079632679489655800e+00; // atan(inf)hi 0x3FF921FB, 0x54442D18
    atanlo[0] = 2.26987774529616870924e-17; // atan(0.5)lo 0x3C7A2B7F, 0x222F65E2
    atanlo[1] = 3.06161699786838301793e-17; // atan(1.0)lo 0x3C81A626, 0x33145C07
    atanlo[2] = 1.39033110312309984516e-17; // atan(1.5)lo 0x3C700788, 0x7AF0CBBD
    atanlo[3] = 6.12323399573676603587e-17; // atan(inf)lo 0x3C91A626, 0x33145C07
    aT[0] = 3.33333333333329318027e-01; // 0x3FD55555, 0x5555550D
    aT[1] = -1.99999999998764832476e-01; // 0xBFC99999, 0x9998EBC4
    aT[2] = 1.42857142725034663711e-01; // 0x3FC24924, 0x920083FF
    aT[3] = -1.11111104054623557880e-01; // 0xBFBC71C6, 0xFE231671
    aT[4] = 9.09088713343650656196e-02; // 0x3FB745CD, 0xC54C206E
    aT[5] = -7.69187620504482999495e-02; // 0xBFB3B0F2, 0xAF749A6D
    aT[6] = 6.66107313738753120669e-02; // 0x3FB10D66, 0xA0D03D51
    aT[7] = -5.83357013379057348645e-02; // 0xBFADDE2D, 0x52DEFD9A
    aT[8] = 4.97687799461593236017e-02; // 0x3FA97B4B, 0x24760DEB
    aT[9] = -3.65315727442169155270e-02; // 0xBFA2B444, 0x2C6A6C2F
    aT[10] = 1.62858201153657823623e-02; // 0x3F90AD3A, 0xE322DA11

    int i,j,k;
    for(i=0; i<n_loc; i += CHUNK_SIZE){
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get ptcl
        reply_get = 0;
        dma_set_op(&dma_get,    DMA_GET);
        dma_set_mode(&dma_get,  PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get,  sizeof(fpLM)*nn);
        dma(dma_get, (long*)(((fpLM *)adr_ptcl_) + my_offset + i), (long*)(ptcl));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for(j=0; j<nn; j++){
            double pos_x = ptcl[j].pos.x;
            double pos_y = ptcl[j].pos.y;
            pos_phi[j] = __ieee754_atan2(pos_y, pos_x, atanhi, atanlo, aT);
            if(pos_phi[j] < 0.0) pos_phi[j] += MY_PI_2;
            pos_r[j] = sqrt(pos_x*pos_x + pos_y*pos_y);
        }
        //** Put adr
        reply_put = 0;
        dma_set_op(&dma_put,    DMA_PUT);
        dma_set_mode(&dma_put,  PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put,  sizeof(double)*nn);
        dma(dma_put, (long*)((double *)adr_pos_phi_ + my_offset + i), (long*)(pos_phi));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
        reply_put = 0;
        dma_set_op(&dma_put,    DMA_PUT);
        dma_set_mode(&dma_put,  PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put,  sizeof(double)*nn);
        dma(dma_put, (long*)((double *)adr_pos_r_ + my_offset + i), (long*)(pos_r));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
}



void FindDomainParticleGoTo(void * args){

    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    int n_ptcl_ = (int)((unsigned long*)args)[0];
    int my_rank_ = (int)((unsigned long*)args)[1];
    int n_proc_x_ = (int)((unsigned long*)args)[2];
    int n_proc_y_ = (int)((unsigned long*)args)[3];
    void * adr_len_peri_x_ = (void*)((unsigned long*)args)[4];
    double len_peri_x_ = *(double*)adr_len_peri_x_;
    void * adr_ptcl_ = (void *)((unsigned long*)args)[5];
    void * adr_pos_domain_ = (void *)((unsigned long*)args)[6];
    void * adr_rank_send_ = (void *)((unsigned long*)args)[7]; //output: size is n_ptcl 

    enum {
        CHUNK_SIZE = 32,
    };
    int n_loc = n_ptcl_/NUMBER_OF_CPE + ( (my_id < n_ptcl_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_ptcl_/NUMBER_OF_CPE)*my_id + ( (my_id < n_ptcl_ % NUMBER_OF_CPE) ? my_id : n_ptcl_ % NUMBER_OF_CPE );
    volatile fpLM ptcl[CHUNK_SIZE];
    volatile int rank_send[CHUNK_SIZE];

    double MY_PI = 3.14159265358979323846;
    double MY_PI_2 = 2.0*MY_PI;
    double * atanhi = (double *)ldm_malloc( sizeof(double)*4 );
    double * atanlo = (double *)ldm_malloc( sizeof(double)*4 );
    double * aT = (double *)ldm_malloc( sizeof(double)*11 );
    atanhi[0] = 4.63647609000806093515e-01; // atan(0.5)hi 0x3FDDAC67, 0x0561BB4F
    atanhi[1] = 7.85398163397448278999e-01; // atan(1.0)hi 0x3FE921FB, 0x54442D18
    atanhi[2] = 9.82793723247329054082e-01; // atan(1.5)hi 0x3FEF730B, 0xD281F69B
    atanhi[3] = 1.57079632679489655800e+00; // atan(inf)hi 0x3FF921FB, 0x54442D18
    atanlo[0] = 2.26987774529616870924e-17; // atan(0.5)lo 0x3C7A2B7F, 0x222F65E2
    atanlo[1] = 3.06161699786838301793e-17; // atan(1.0)lo 0x3C81A626, 0x33145C07
    atanlo[2] = 1.39033110312309984516e-17; // atan(1.5)lo 0x3C700788, 0x7AF0CBBD
    atanlo[3] = 6.12323399573676603587e-17; // atan(inf)lo 0x3C91A626, 0x33145C07
    aT[0] = 3.33333333333329318027e-01; // 0x3FD55555, 0x5555550D
    aT[1] = -1.99999999998764832476e-01; // 0xBFC99999, 0x9998EBC4
    aT[2] = 1.42857142725034663711e-01; // 0x3FC24924, 0x920083FF
    aT[3] = -1.11111104054623557880e-01; // 0xBFBC71C6, 0xFE231671
    aT[4] = 9.09088713343650656196e-02; // 0x3FB745CD, 0xC54C206E
    aT[5] = -7.69187620504482999495e-02; // 0xBFB3B0F2, 0xAF749A6D
    aT[6] = 6.66107313738753120669e-02; // 0x3FB10D66, 0xA0D03D51
    aT[7] = -5.83357013379057348645e-02; // 0xBFADDE2D, 0x52DEFD9A
    aT[8] = 4.97687799461593236017e-02; // 0x3FA97B4B, 0x24760DEB
    aT[9] = -3.65315727442169155270e-02; // 0xBFA2B444, 0x2C6A6C2F
    aT[10] = 1.62858201153657823623e-02; // 0x3F90AD3A, 0xE322DA11

    int i,j,k;
    int n_cnt = 0;
    int rank_prev = my_rank_;
    volatile F64ort_ * pos_domain_prev = (F64ort_ *)ldm_malloc( sizeof(F64ort_) );
    reply_get = 0;
    dma_set_op(&dma_get,    DMA_GET);
    dma_set_mode(&dma_get,  PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get,  sizeof(F64ort_));
    dma(dma_get, (long*)(((F64ort_ *)adr_pos_domain_) + rank_prev), (long*)(pos_domain_prev));
    dma_wait(&reply_get, 1);
    while (reply_get != 1) {}
    for(i=0; i<n_loc; i += CHUNK_SIZE){
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get ptcl
        reply_get = 0;
        dma_set_op(&dma_get,    DMA_GET);
        dma_set_mode(&dma_get,  PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get,  sizeof(fpLM)*nn);
        dma(dma_get, (long*)(((fpLM *)adr_ptcl_) + my_offset + i), (long*)(ptcl));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for(j=0; j<nn; j++){
            double pos_x = ptcl[j].pos.x;
            double pos_y = ptcl[j].pos.y;
            double pos_z = ptcl[j].pos.z;
            double pos_phi = __ieee754_atan2(pos_y, pos_x, atanhi, atanlo, aT);
            if(pos_phi < 0.0) pos_phi += MY_PI_2;
            double pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
            rank_send[j] = SearchRank(pos_phi, pos_r, n_proc_x_, n_proc_y_,
                                      adr_pos_domain_, len_peri_x_, rank_prev,
                                      pos_domain_prev);
            rank_prev = rank_send[j];
        }
        //** Put adr
        reply_put = 0;
        dma_set_op(&dma_put,    DMA_PUT);
        dma_set_mode(&dma_put,  PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put,  sizeof(int)*nn);
        dma(dma_put, (long*)((int *)adr_rank_send_ + my_offset + i), (long*)(rank_send));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    ldm_free(pos_domain_prev, sizeof(F64ort_));
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
}

int SearchRank(double pos_phi,
               double pos_r,
               int n_proc_x,
               int n_proc_y,
               F64ort_ * adr_pos_domain,
               double len_peri_x,
               int rank_prev,
               F64ort_ * pos_domain){

    int idomain = rank_prev;
    if(pos_phi >= (pos_domain->low_).x && (pos_domain->high_).x > pos_phi
       && pos_r >= (pos_domain->low_).y && (pos_domain->high_).y > pos_r){
        // is in the same domain
        return idomain;
    }

    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int shift = n_proc_y;
    double box_cen_x = (pos_domain->high_.x + pos_domain->low_.x)*0.5;
    double box_len_x = (pos_domain->high_.x - pos_domain->low_.x)*0.5;
    double dx = pos_phi-box_cen_x;
    if(dx > box_len_x && dx > len_peri_x*0.5){
        double rank_y = (idomain/shift) % n_proc_y;
        idomain = rank_y*1 + 0*shift;
    }
    else if(-dx >= box_len_x && -dx > len_peri_x*0.5){
        double rank_y = (idomain/shift) % n_proc_y;
        idomain = rank_y*1 + (n_proc_x-1)*shift;
    }
    assert( idomain >= 0 && idomain < (n_proc_x*n_proc_y) );
    reply_get = 0;
    dma_set_op(&dma_get,    DMA_GET);
    dma_set_mode(&dma_get,  PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get,  sizeof(F64ort_));
    dma(dma_get, (long*)(((F64ort_ *)adr_pos_domain) + idomain), (long*)(pos_domain));
    dma_wait(&reply_get, 1);

    // x direction
    if(pos_domain->high_.x <= pos_phi){
        do{
            idomain += shift;
            reply_get = 0;
            dma_set_op(&dma_get,    DMA_GET);
            dma_set_mode(&dma_get,  PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get,  sizeof(F64ort_));
            dma(dma_get, (long*)(((F64ort_ *)adr_pos_domain) + idomain), (long*)(pos_domain));
            dma_wait(&reply_get, 1);
        }while(pos_domain->high_.x <= pos_phi);
    }
    else if(pos_domain->low_.x > pos_phi){
        do{
            idomain -= shift;
            reply_get = 0;
            dma_set_op(&dma_get,    DMA_GET);
            dma_set_mode(&dma_get,  PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get,  sizeof(F64ort_));
            dma(dma_get, (long*)(((F64ort_ *)adr_pos_domain) + idomain), (long*)(pos_domain));
            dma_wait(&reply_get, 1);            
        }while(pos_domain->low_.x > pos_phi);
    }

    //assert(pos_phi >= (pos_domain->low_).x && (pos_domain->high_).x > pos_phi);

    
    // y direction
    shift = 1;
    if(pos_domain->high_.y <= pos_r){
        do{
            idomain += shift;
            //** Get ptcl
            reply_get = 0;
            dma_set_op(&dma_get,    DMA_GET);
            dma_set_mode(&dma_get,  PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get,  sizeof(F64ort_));
            dma(dma_get, (long*)(((F64ort_ *)adr_pos_domain) + idomain), (long*)(pos_domain));
            dma_wait(&reply_get, 1);
        }while(pos_domain->high_.y <= pos_r);
    }
    else if(pos_domain->low_.y > pos_r){
        do{
            idomain -= shift;
            //** Get ptcl
            reply_get = 0;
            dma_set_op(&dma_get,    DMA_GET);
            dma_set_mode(&dma_get,  PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get,  sizeof(F64ort_));
            dma(dma_get, (long*)(((F64ort_ *)adr_pos_domain) + idomain), (long*)(pos_domain));
            dma_wait(&reply_get, 1);
        }while(pos_domain->low_.y <= pos_r);
    }

    return idomain; 
}


void SetParticleToSendBuffer(void * args){

    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int i,j,k,id;
    
    enum {
        CHUNK_SIZE = 32,
        N_TARGET_LIMIT = 3,
    };

    int n_ptcl_ = (int)((unsigned long*)args)[0];
    int my_rank_ = (int)((unsigned long*)args)[1];
    int n_proc_x_ = (int)((unsigned long*)args)[2];
    int n_proc_y_ = (int)((unsigned long*)args)[3];
    void * adr_len_peri_x_ = (void*)((unsigned long*)args)[4];
    double len_peri_x_ = *(double*)adr_len_peri_x_;
    void * adr_ptcl_ = (void *)((unsigned long*)args)[5];
    void * adr_pos_domain_ = (void *)((unsigned long*)args)[6];
    void * adr_rank_send_ = (void *)((unsigned long*)args)[7]; //output: size is n_ptcl 
    int n_proc_target_ = (int)((unsigned long*)args)[8];
    void * adr_adr_no_set_ = (void*)((unsigned long*)args)[9];
    void * adr_n_no_set_ = (void*)((unsigned long*)args)[10];

    void * adr_ptcl_send_[N_TARGET_LIMIT];
    for(i=0; i<N_TARGET_LIMIT; i++){
        adr_ptcl_send_[i] = (void *)((unsigned long*)args)[i+11];
    }
    void * adr_n_send_[N_TARGET_LIMIT];
    for(i=0; i<N_TARGET_LIMIT; i++){
        adr_n_send_[i] = (void *)((unsigned long*)args)[i+11+N_TARGET_LIMIT];
    }
    int rank_target_[N_TARGET_LIMIT];
    for(i=0; i<N_TARGET_LIMIT; i++){
        rank_target_[i] = (int)((unsigned long*)args)[i+11+N_TARGET_LIMIT*2];
    }
    
    volatile int n_send[N_TARGET_LIMIT][NUMBER_OF_CPE];
    volatile int n_send_total[N_TARGET_LIMIT];
    volatile int n_send_offset[N_TARGET_LIMIT];
    for(i=0; i<N_TARGET_LIMIT; i++){
        n_send_total[i] = n_send_offset[i] = 0;
    }
    volatile int n_no_set[NUMBER_OF_CPE];
    int n_no_set_total = 0;
    int n_no_set_offset = 0;
    

    int n_loc = n_ptcl_/NUMBER_OF_CPE + ( (my_id < n_ptcl_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_ptcl_/NUMBER_OF_CPE)*my_id + ( (my_id < n_ptcl_ % NUMBER_OF_CPE) ? my_id : n_ptcl_ % NUMBER_OF_CPE );
    volatile fpLM ptcl[CHUNK_SIZE]; // for input
    volatile fpLM ptcl_send[N_TARGET_LIMIT][CHUNK_SIZE]; // for output
    volatile int adr_no_set[CHUNK_SIZE];
    volatile int rank_send[CHUNK_SIZE];
        
    double MY_PI = 3.14159265358979323846;
    double MY_PI_2 = 2.0*MY_PI;
    double * atanhi = (double *)ldm_malloc( sizeof(double)*4 );
    double * atanlo = (double *)ldm_malloc( sizeof(double)*4 );
    double * aT = (double *)ldm_malloc( sizeof(double)*11 );
    atanhi[0] = 4.63647609000806093515e-01; // atan(0.5)hi 0x3FDDAC67, 0x0561BB4F
    atanhi[1] = 7.85398163397448278999e-01; // atan(1.0)hi 0x3FE921FB, 0x54442D18
    atanhi[2] = 9.82793723247329054082e-01; // atan(1.5)hi 0x3FEF730B, 0xD281F69B
    atanhi[3] = 1.57079632679489655800e+00; // atan(inf)hi 0x3FF921FB, 0x54442D18
    atanlo[0] = 2.26987774529616870924e-17; // atan(0.5)lo 0x3C7A2B7F, 0x222F65E2
    atanlo[1] = 3.06161699786838301793e-17; // atan(1.0)lo 0x3C81A626, 0x33145C07
    atanlo[2] = 1.39033110312309984516e-17; // atan(1.5)lo 0x3C700788, 0x7AF0CBBD
    atanlo[3] = 6.12323399573676603587e-17; // atan(inf)lo 0x3C91A626, 0x33145C07
    aT[0] = 3.33333333333329318027e-01; // 0x3FD55555, 0x5555550D
    aT[1] = -1.99999999998764832476e-01; // 0xBFC99999, 0x9998EBC4
    aT[2] = 1.42857142725034663711e-01; // 0x3FC24924, 0x920083FF
    aT[3] = -1.11111104054623557880e-01; // 0xBFBC71C6, 0xFE231671
    aT[4] = 9.09088713343650656196e-02; // 0x3FB745CD, 0xC54C206E
    aT[5] = -7.69187620504482999495e-02; // 0xBFB3B0F2, 0xAF749A6D
    aT[6] = 6.66107313738753120669e-02; // 0x3FB10D66, 0xA0D03D51
    aT[7] = -5.83357013379057348645e-02; // 0xBFADDE2D, 0x52DEFD9A
    aT[8] = 4.97687799461593236017e-02; // 0x3FA97B4B, 0x24760DEB
    aT[9] = -3.65315727442169155270e-02; // 0xBFA2B444, 0x2C6A6C2F
    aT[10] = 1.62858201153657823623e-02; // 0x3F90AD3A, 0xE322DA11

    int rank_prev = my_rank_;
    F64ort_ * pos_domain_prev = (F64ort_ *)ldm_malloc( sizeof(F64ort_) );
    for(i=0; i<NUMBER_OF_CPE; i++) n_no_set[i] = 0;
    reply_get = 0;
    dma_set_op(&dma_get,    DMA_GET);
    dma_set_mode(&dma_get,  PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get,  sizeof(F64ort_));
    dma(dma_get, (long*)(((F64ort_ *)adr_pos_domain_) + rank_prev), (long*)(pos_domain_prev));
    dma_wait(&reply_get, 1);
    while (reply_get != 1) {}

    int loop = 0;
    int loop_max = (((n_ptcl_/NUMBER_OF_CPE)+1)/CHUNK_SIZE)+1;
    for(i=0; ; i+=CHUNK_SIZE){
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        n_no_set[my_id] = 0;
        for(k=0; k<n_proc_target_; k++){ n_send[k][my_id] = 0; }
        if(i<n_loc){
            //** Get ptcl
            reply_get = 0;
            dma_set_op(&dma_get,    DMA_GET);
            dma_set_mode(&dma_get,  PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get,  sizeof(fpLM)*nn);
            dma(dma_get, (long*)(((fpLM *)adr_ptcl_) + my_offset + i), (long*)(ptcl));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            for(j=0; j<nn; j++){
                double pos_x = ptcl[j].pos.x;
                double pos_y = ptcl[j].pos.y;
                double pos_z = ptcl[j].pos.z;
                double pos_phi = __ieee754_atan2(pos_y, pos_x, atanhi, atanlo, aT);
                if(pos_phi < 0.0) pos_phi += MY_PI_2;
                double pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                rank_send[j] = SearchRank(pos_phi, pos_r, n_proc_x_, n_proc_y_,
                                          adr_pos_domain_, len_peri_x_, rank_prev,
                                          pos_domain_prev);
                rank_prev = rank_send[j];
                /*
                assert(pos_phi >= (pos_domain_prev->low_).x
                       && (pos_domain_prev->high_).x > pos_phi
                       && pos_r >= (pos_domain_prev->low_).y
                       && (pos_domain_prev->high_).y > pos_r);
                */
                int flag_send = -1;
                for(k=0; k<n_proc_target_; k++){
                    if(rank_send[j] == rank_target_[k]){
                        int n_tmp = n_send[k][my_id];
                        ptcl_send[k][n_tmp].pos.x = pos_x;
                        ptcl_send[k][n_tmp].pos.y = pos_y;
                        ptcl_send[k][n_tmp].pos.z = pos_z;
                        ptcl_send[k][n_tmp].vel.x = ptcl[j].vel.x;
                        ptcl_send[k][n_tmp].vel.y = ptcl[j].vel.y;
                        ptcl_send[k][n_tmp].vel.z = ptcl[j].vel.z;
                        ptcl_send[k][n_tmp].mass = ptcl[j].mass;
                        ptcl_send[k][n_tmp].id = ptcl[j].id;
                        n_send[k][my_id]++;
                        flag_send = k;
                        break;
                    }
                }
                if(flag_send == -1){
                    adr_no_set[n_no_set[my_id]] = my_offset + i + j;
                    n_no_set[my_id]++;
                }

                // PUT rank_send
                reply_put = 0;
                dma_set_op(&dma_put,    DMA_PUT);
                dma_set_mode(&dma_put,  PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put,  sizeof(int)*nn);
                dma(dma_put, (long*)(((int *)adr_rank_send_) + my_offset + i), (long*)(rank_send));
                dma_wait(&reply_put, 1);
                while (reply_put != 1) {}
            }
        }
        
        for(k=0; k<n_proc_target_; k++){
            n_send_offset[k] = n_send_total[k];
            for(id=0; id<NUMBER_OF_CPE; id++){
                sync_array_();
                cpe_bcast_int32(id, &(n_send[k][id]));
                sync_array_();
                n_send_total[k] += n_send[k][id];
                if(id < my_id){
                    n_send_offset[k] += n_send[k][id];
                }
            }
        }
        n_no_set_offset = n_no_set_total;
        for(id=0; id<NUMBER_OF_CPE; id++){
            sync_array_();
            cpe_bcast_int32(id, &(n_no_set[id]));
            sync_array_();
            n_no_set_total += n_no_set[id];
            if(id < my_id){
                n_no_set_offset += n_no_set[id];
            }
        }

        if(i<n_loc){
            for(k=0; k<n_proc_target_; k++){
                int n_tmp = n_send[k][my_id];
                if( n_tmp > 0){
                    int n_offset_tmp = n_send_offset[k];
                    reply_put = 0;
                    dma_set_op(&dma_put,    DMA_PUT);
                    dma_set_mode(&dma_put,  PE_MODE);
                    dma_set_reply(&dma_put, &reply_put);
                    dma_set_size(&dma_put,  sizeof(fpLM)*n_tmp);
                    dma(dma_put, (long*)(((fpLM*)(adr_ptcl_send_[k]))+n_offset_tmp), (long*)(ptcl_send[k]));
                    dma_wait(&reply_put, 1);
                    while (reply_put != 1) {}
                }
            }
            if( n_no_set[my_id] > 0){
                reply_put = 0;
                dma_set_op(&dma_put,    DMA_PUT);
                dma_set_mode(&dma_put,  PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put,  sizeof(int)*n_no_set[my_id]);
                dma(dma_put, (long*)(((int*)adr_adr_no_set_)+n_no_set_offset), (long*)(adr_no_set));
                dma_wait(&reply_put, 1);
                while (reply_put != 1) {}
            }
        }
        
        if(loop >= loop_max) break;
        loop++;
    } // end of loop over particle

    if(my_id == 0){
        // put n_send
        for(i=0; i<n_proc_target_; i++){
            reply_put = 0;
            dma_set_op(&dma_put,    DMA_PUT);
            dma_set_mode(&dma_put,  PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put,  sizeof(int));
            dma(dma_put, (long*)(((int *)(adr_n_send_[i]))), (long*)(((int*)n_send_total)+i));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }

        // Put n_no_set
        reply_put = 0;
        dma_set_op(&dma_put,    DMA_PUT);
        dma_set_mode(&dma_put,  PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put,  sizeof(int));
        dma(dma_put, (long*)((int *)adr_n_no_set_), (long*)((int*)(&n_no_set_total)));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
        
    }

    ldm_free(pos_domain_prev, sizeof(F64ort_));
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
}

//######################################################################################
static int searchWhichDomainParticleGoToPeriodicX(F64vec_ *pos,
                                                  int n_domain [],
                                                  void *adr_pos_domain) {
    volatile int idomain = previous_rank;
    volatile F64ort_ pos_domain;
#ifdef CHECK_FIND_DEST_RANK
    if (debug_flag) {
        printf("Enter searchWhichDomainParticleGoToPeriodicX().\n");
        printf("  pos->x = %e\n",pos->x);
        printf("  pos->y = %e\n",pos->y);
        printf("  pos->z = %e\n",pos->z);
        printf("  n_domain[0] = %d\n",n_domain[0]);
        printf("  n_domain[1] = %d\n",n_domain[1]);
        printf("  n_domain[2] = %d\n",n_domain[2]);
        printf("  previous_rank = %d\n",previous_rank);
        printf("  have_pos_domain_prev = %d\n",have_pos_domain_prev);
    }
#endif

    //* Check if this particle is also in the domain of the previous rank
    if (!have_pos_domain_prev) {
#ifdef CHECK_FIND_DEST_RANK
        if (debug_flag) {
            printf("DMA get of pos_domain_prev is started!\n");
        }
#endif
        //** Get pos_domain_prev
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get,  sizeof(F64ort_));
        dma(dma_get, (long*)((F64ort_ *)adr_pos_domain + idomain), (long*)(&pos_domain));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        //** Save pos_domain_prev
        pos_domain_prev.low_.x  = pos_domain.low_.x ;
        pos_domain_prev.low_.y  = pos_domain.low_.y ;
        pos_domain_prev.low_.z  = pos_domain.low_.z ;
        pos_domain_prev.high_.x = pos_domain.high_.x;
        pos_domain_prev.high_.y = pos_domain.high_.y;
        pos_domain_prev.high_.z = pos_domain.high_.z;
        //** Update have_pos_domain_prev
        have_pos_domain_prev = 1;
    } else {
        //** Copy from pos_domain_prev
        pos_domain.low_.x  = pos_domain_prev.low_.x ;
        pos_domain.low_.y  = pos_domain_prev.low_.y ;
        pos_domain.low_.z  = pos_domain_prev.low_.z ;
        pos_domain.high_.x = pos_domain_prev.high_.x;
        pos_domain.high_.y = pos_domain_prev.high_.y;
        pos_domain.high_.z = pos_domain_prev.high_.z;
    }
#ifdef CHECK_FIND_DEST_RANK
    if (debug_flag) {
        printf("DMA_get of pos_domain_prev is completed!\n");
        printf("  pos_domain_prev.low_.x = %e\n",pos_domain_prev.low_.x);
        printf("  pos_domain_prev.low_.y = %e\n",pos_domain_prev.low_.y);
        printf("  pos_domain_prev.low_.z = %e\n",pos_domain_prev.low_.z);
        printf("  pos_domain_prev.high_.x = %e\n",pos_domain_prev.high_.x);
        printf("  pos_domain_prev.high_.y = %e\n",pos_domain_prev.high_.y);
        printf("  pos_domain_prev.high_.z = %e\n",pos_domain_prev.high_.z);
    }
#endif
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
    if ((pos_domain.low_.x <= pos->x) && (pos->x < pos_domain.high_.x) &&
        (pos_domain.low_.y <= pos->y) && (pos->y < pos_domain.high_.y) &&
        (pos_domain.low_.z <= pos->z) && (pos->z < pos_domain.high_.z)) {
        return idomain;
    }
#else
    if ((pos_domain.low_.x <= pos->x) && (pos->x < pos_domain.high_.x) &&
        (pos_domain.low_.y <= pos->y) && (pos->y < pos_domain.high_.y)) {
        return idomain;
    }
#endif

#ifdef CHECK_FIND_DEST_RANK
    if (debug_flag) {
        printf("Start to search the destination rank using do-while loops.\n");
    }
#endif
    
    //* Search a process whose domain contains this particle
    //** Set the parameters of DMA get
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get,  sizeof(F64ort_));
    //** x direction
    int shift = n_domain[1] * n_domain[2];
    if (pos_domain.high_.x <= pos->x){
#if defined(CHECK_FIND_DEST_RANK) && (CHECK_LVL_FIND_DEST_RANK >= 2)
        if (debug_flag) {
            printf("Enter 1st if branch in x-searcing.\n");
        }
#endif
        do {
            idomain += shift;
#if defined(CHECK_FIND_DEST_RANK) && (CHECK_LVL_FIND_DEST_RANK >= 3)
            if (debug_flag) {
                printf("idomain = %d\n",idomain);
            }
#endif
            // Update pos_domain
            reply_get = 0;
            dma(dma_get, (long*)((F64ort_ *)adr_pos_domain + idomain), (long*)(&pos_domain));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
        } while (pos_domain.high_.x <= pos->x);
    }
    else if (pos_domain.low_.x > pos->x){
#if defined(CHECK_FIND_DEST_RANK) && (CHECK_LVL_FIND_DEST_RANK >= 2)
        if (debug_flag) {
            printf("Enter 2nd if branch in x-searcing.\n");
        }
#endif
        do {
            idomain -= shift;
#if defined(CHECK_FIND_DEST_RANK) && (CHECK_LVL_FIND_DEST_RANK >= 3)
            if (debug_flag) {
                printf("idomain = %d\n",idomain);
            }
#endif
            // Update pos_domain
            reply_get = 0;
            dma(dma_get, (long*)((F64ort_ *)adr_pos_domain + idomain), (long*)(&pos_domain));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
        } while (pos_domain.low_.x > pos->x);
    }
#ifdef CHECK_FIND_DEST_RANK
    if (debug_flag) {
        printf("Searching along x direction is completed!\n");
    }
#endif

    //** y direction
    shift = n_domain[2];
    if (pos_domain.high_.y <= pos->y) {
        do { 
            idomain += shift;
            // Update pos_domain
            reply_get = 0;
            dma(dma_get, (long*)((F64ort_ *)adr_pos_domain + idomain), (long*)(&pos_domain));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
        } while (pos_domain.high_.y <= pos->y);
    }
    else if (pos_domain.low_.y > pos->y) {
        do { 
            idomain -= shift;
            // Update pos_domain
            reply_get = 0;
            dma(dma_get, (long*)((F64ort_ *)adr_pos_domain + idomain), (long*)(&pos_domain));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
        } while (pos_domain.low_.y > pos->y);
    }
#ifdef CHECK_FIND_DEST_RANK
    if (debug_flag) {
        printf("Searching along y direction is completed!\n");
    }
#endif

    //** z direction
    if (pos_domain.high_.z <= pos->z) {
        do { 
            idomain++;
            // Update pos_domain
            reply_get = 0;
            dma(dma_get, (long*)((F64ort_ *)adr_pos_domain + idomain), (long*)(&pos_domain));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
        } while (pos_domain.high_.z <= pos->z);
    }
    else if (pos_domain.low_.z > pos->z) {
        do { 
            idomain--;
            // Update pos_domain
            reply_get = 0;
            dma(dma_get, (long*)((F64ort_ *)adr_pos_domain + idomain), (long*)(&pos_domain));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
        } while (pos_domain.low_.z > pos->z);
    }
#ifdef CHECK_FIND_DEST_RANK
    if (debug_flag) {
        printf("Searching along z direction is completed!\n");
    }
#endif
    
    //* Update previous_rank & pos_domain_prev
    previous_rank = idomain;
    pos_domain_prev.low_.x  = pos_domain.low_.x ;
    pos_domain_prev.low_.y  = pos_domain.low_.y ;
    pos_domain_prev.low_.z  = pos_domain.low_.z ;
    pos_domain_prev.high_.x = pos_domain.high_.x;
    pos_domain_prev.high_.y = pos_domain.high_.y;
    pos_domain_prev.high_.z = pos_domain.high_.z;

    return idomain;
}

void FindDestinationRank(void *args) {
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    //* Initialize local vars. for DMA comm.
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    //* Process the arguments
    int myrank                   = (int)(((unsigned long*)args)[0]);
    int nprocs                   = (int)(((unsigned long*)args)[1]);
    void *adr_n_domain           = (void *)((unsigned long*)args)[2];
    int nptcl_loc                = (int)(((unsigned long*)args)[3]);
    void *adr_ptcl_              = (void *)((unsigned long*)args)[4];
    void *adr_pos_domain         = (void *)((unsigned long*)args)[5];
    void *adr_dest_rank          = (void *)((unsigned long*)args)[6];
    void *adr_n_dest_rank_cpe    = (void *)((unsigned long*)args)[7];
    void *adr_dest_rank_list_cpe = (void *)((unsigned long*)args)[8];
    int debug_flag_from_mpe      = (int)(((unsigned long*)args)[9]);
#ifdef CHECK_FIND_DEST_RANK
    if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
        printf("myrank                 = %d\n",myrank);
        printf("nprocs                 = %d\n",nprocs);
        printf("adr_n_domain           = %lu\n",(unsigned long)adr_n_domain);
        printf("nptcl_loc              = %d\n",nptcl_loc);
        printf("adr_ptcl_              = %lu\n",(unsigned long)adr_ptcl_);
        printf("adr_pos_domain         = %lu\n",(unsigned long)adr_pos_domain);
        printf("adr_dest_rank          = %lu\n",(unsigned long)adr_dest_rank);
        printf("adr_n_dest_rank_cpe    = %lu\n",(unsigned long)adr_n_dest_rank_cpe);
        printf("adr_dest_rank_list_cpe = %lu\n",(unsigned long)adr_dest_rank_list_cpe);
        printf("debug_flag_from_mpe    = %d\n",debug_flag_from_mpe);
    }
#endif
    //* Get n_domain[]
    volatile int n_domain[3];
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get,  sizeof(int) * 3);
    dma(dma_get, (long*)(adr_n_domain), (long*)(n_domain));
    dma_wait(&reply_get, 1);
    while (reply_get != 1) {}
#ifdef CHECK_FIND_DEST_RANK
    if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
        printf("n_domain[0]         = %d\n",n_domain[0]);
        printf("n_domain[1]         = %d\n",n_domain[1]);
        printf("n_domain[2]         = %d\n",n_domain[2]);
    }
#endif
    //* Local variables
    //-(Loop counters)
    int i,j,jloc,k,id;
    //-(atan)
    double MY_PI = 3.14159265358979323846;
    double MY_2_PI = 2.0*MY_PI;
    double * atanhi = (double *)ldm_malloc( sizeof(double)*4 );
    double * atanlo = (double *)ldm_malloc( sizeof(double)*4 );
    double * aT = (double *)ldm_malloc( sizeof(double)*11 );
    atanhi[0] = 4.63647609000806093515e-01; // atan(0.5)hi 0x3FDDAC67, 0x0561BB4F
    atanhi[1] = 7.85398163397448278999e-01; // atan(1.0)hi 0x3FE921FB, 0x54442D18
    atanhi[2] = 9.82793723247329054082e-01; // atan(1.5)hi 0x3FEF730B, 0xD281F69B
    atanhi[3] = 1.57079632679489655800e+00; // atan(inf)hi 0x3FF921FB, 0x54442D18
    atanlo[0] = 2.26987774529616870924e-17; // atan(0.5)lo 0x3C7A2B7F, 0x222F65E2
    atanlo[1] = 3.06161699786838301793e-17; // atan(1.0)lo 0x3C81A626, 0x33145C07
    atanlo[2] = 1.39033110312309984516e-17; // atan(1.5)lo 0x3C700788, 0x7AF0CBBD
    atanlo[3] = 6.12323399573676603587e-17; // atan(inf)lo 0x3C91A626, 0x33145C07
    aT[0] = 3.33333333333329318027e-01; // 0x3FD55555, 0x5555550D
    aT[1] = -1.99999999998764832476e-01; // 0xBFC99999, 0x9998EBC4
    aT[2] = 1.42857142725034663711e-01; // 0x3FC24924, 0x920083FF
    aT[3] = -1.11111104054623557880e-01; // 0xBFBC71C6, 0xFE231671
    aT[4] = 9.09088713343650656196e-02; // 0x3FB745CD, 0xC54C206E
    aT[5] = -7.69187620504482999495e-02; // 0xBFB3B0F2, 0xAF749A6D
    aT[6] = 6.66107313738753120669e-02; // 0x3FB10D66, 0xA0D03D51
    aT[7] = -5.83357013379057348645e-02; // 0xBFADDE2D, 0x52DEFD9A
    aT[8] = 4.97687799461593236017e-02; // 0x3FA97B4B, 0x24760DEB
    aT[9] = -3.65315727442169155270e-02; // 0xBFA2B444, 0x2C6A6C2F
    aT[10] = 1.62858201153657823623e-02; // 0x3F90AD3A, 0xE322DA11
    //-(local buffers for ptcl_ & dest_rank[], etc.)
    enum {
        CHUNK_SIZE = 64,
        TBL_SIZE   = 2048,
        CYCLE_SIZE = CHUNK_SIZE * NUMBER_OF_CPE,
    };
    size_t bsize_ptcl_array = sizeof(fpLM) * CHUNK_SIZE;
    size_t bsize_rank_array = sizeof(int)  * CHUNK_SIZE;
    size_t bsize_tbl        = sizeof(int)  * TBL_SIZE;
    volatile fpLM *ptcl     = (fpLM *) ldm_malloc( bsize_ptcl_array );
    volatile int *dest_rank = (int *)  ldm_malloc( bsize_rank_array );
    volatile int n_dest_rank = 0;
    volatile int *dest_rank_list = (int *)  ldm_malloc( bsize_tbl );

    //* Examine the destination rank
#ifdef SPECIFY_CPE_IN_FIND_DEST_RANK
    if (my_id == ID_CHECK) {
#endif
    for (i=0; i<nptcl_loc; i+=CYCLE_SIZE) {
        int nrem = nptcl_loc - i;
        int nn   = nrem < CYCLE_SIZE ? nrem : CYCLE_SIZE;
        //** Compute the task of each CPE
        int n_loc, my_offset;
        n_loc = nn/NUMBER_OF_CPE + ( (my_id < nn % NUMBER_OF_CPE) ? 1 : 0 );
        my_offset = i + (nn/NUMBER_OF_CPE)*my_id + ( (my_id < nn % NUMBER_OF_CPE) ? my_id : nn % NUMBER_OF_CPE );
#ifdef CHECK_FIND_DEST_RANK
        // Output i
        if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
            printf("i         = %d\n",i);
            printf("n_loc     = %d\n",n_loc);
            printf("my_offset = %d\n",my_offset);
        }
#endif
        if (n_loc > 0) {
            //** Get ptcl_
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get,  sizeof(fpLM) * n_loc);
            dma(dma_get, (long*)((fpLM *)adr_ptcl_ + my_offset), (long*)(ptcl));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            //** Examine the destination rank
            for (k=0; k<n_loc; k++) {
#ifdef CHECK_FIND_DEST_RANK
                // Output k
                if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
                    printf("k = %d (%d)\n",k,my_offset+k);
                }
#endif
                //*** Compute \phi and r
                double pos_x = ptcl[k].pos.x;
                double pos_y = ptcl[k].pos.y;
                double pos_z = ptcl[k].pos.z;
                double pos_phi = __ieee754_atan2(pos_y, pos_x, atanhi, atanlo, aT);
                if (pos_phi < 0.0) pos_phi += MY_2_PI;
                double pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                //*** Search the destination rank
                F64vec_ pos;
                pos.x = pos_phi;
                pos.y = pos_r;
                pos.z = pos_z;
                dest_rank[k] = searchWhichDomainParticleGoToPeriodicX(&pos, 
                                                                      n_domain,
                                                                      adr_pos_domain); 
                //*** Compute n_dest_rank, dest_rank_list
                if (n_dest_rank > 0) {
                    assert(n_dest_rank < TBL_SIZE);
                    //* Check if dest_rank[k] is already registered or not.
                    jloc = -1;
                    for (j=0; j<n_dest_rank; j++) {
                        if (dest_rank_list[j] == dest_rank[k]) {
                            jloc = j;
                            break;
                        }
                    }
                    // Update tables
                    if (jloc == -1) {
                        // In this case, dest_rank[k] is a new rank.
                        dest_rank_list[n_dest_rank] = dest_rank[k];
                        n_dest_rank++;
                    }
                } else {
                    n_dest_rank++;
                    dest_rank_list[0] = dest_rank[k];
                }
            }
            //** Put dest_rank
            reply_put = 0;
            dma_set_op(&dma_put,    DMA_PUT);
            dma_set_mode(&dma_put,  PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put,  sizeof(int) * n_loc);
            dma(dma_put, (long*)((int *)adr_dest_rank + my_offset), (long*)(dest_rank));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }
    }
#ifdef SPECIFY_CPE_IN_FIND_DEST_RANK
    }
#endif
    //* Put n_dest_rank, dest_rank_list
    assert(n_dest_rank <= TBL_SIZE);
    //** n_dest_rank 
    reply_put = 0;
    dma_set_op(&dma_put,    DMA_PUT);
    dma_set_mode(&dma_put,  PE_MODE);
    dma_set_reply(&dma_put, &reply_put);
    dma_set_size(&dma_put,  sizeof(int));
    dma(dma_put, (long*)((int *)adr_n_dest_rank_cpe + my_id), (long*)(&n_dest_rank));
    dma_wait(&reply_put, 1);
    while (reply_put != 1) {}
    //** dest_rank_list
    if (n_dest_rank > 0) {
        reply_put = 0;
        dma_set_op(&dma_put,    DMA_PUT);
        dma_set_mode(&dma_put,  PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put,  sizeof(int) * n_dest_rank);
        dma(dma_put, (long*)((int *)adr_dest_rank_list_cpe + my_id * nprocs), (long*)(dest_rank_list));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }

    //* Release memory
    //-(atan2)
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
    //-(buffers)
    ldm_free(ptcl, bsize_ptcl_array);
    ldm_free(dest_rank, bsize_rank_array);
    ldm_free(dest_rank_list, bsize_tbl);

}

void MakeSendBuffers(void *args) {
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    //* Initialize local vars. for DMA comm.
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    //* Process the arguments
    int myrank                   = (int)(((unsigned long*)args)[0]);
    int nprocs                   = (int)(((unsigned long*)args)[1]);
    int nptcl_loc                = (int)(((unsigned long*)args)[2]);
    void *adr_ptcl_              = (void *)((unsigned long*)args)[3];
    void *adr_dest_rank          = (void *)((unsigned long*)args)[4];
    int n_dest_rank              = (int)(((unsigned long*)args)[5]);
    void *adr_dest_rank_list     = (void *)((unsigned long*)args)[6];
    void *adr_adr_ptcl_send_buf_ = (void *)((unsigned long*)args)[7];
    int debug_flag_from_mpe      = (int)(((unsigned long*)args)[8]);
#ifdef CHECK_MAKE_SEND_BUFS
    sync_array_();
    if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
        printf("myrank                 = %d\n",myrank);
        printf("nprocs                 = %d\n",nprocs);
        printf("nptcl_loc              = %d\n",nptcl_loc);
        printf("adr_ptcl_              = %lu\n",(unsigned long)adr_ptcl_);
        printf("adr_dest_rank          = %lu\n",(unsigned long)adr_dest_rank);
        printf("n_dest_rank            = %d\n",n_dest_rank);
        printf("adr_dest_rank_list     = %lu\n",(unsigned long)adr_dest_rank_list);
        printf("adr_adr_ptcl_send_buf_ = %lu\n",(unsigned long)adr_adr_ptcl_send_buf_);
        printf("debug_flag_from_mpe    = %d\n",debug_flag_from_mpe);
    }
    sync_array_();
#endif
    //* Local variables
    //-(Loop counters)
    int i,j,jglb,jloc,k,rank,id;
    int adr;
    //-(prefix sum calc.)
    int val,beg,end,psum;
    //-(local buffers for ptcl_ & dest_rank[], etc.)
    enum {
        CHUNK_SIZE = 64,
        TBL_SIZE = 2048,
        CYCLE_SIZE = CHUNK_SIZE * NUMBER_OF_CPE,
    };
    size_t bsize_ptcl_chunk = sizeof(fpLM) * CHUNK_SIZE;
    size_t bsize_int_chunk  = sizeof(int)  * CHUNK_SIZE;
    size_t bsize_int_tbl    = sizeof(int)  * TBL_SIZE;
    size_t bsize_long_tbl   = sizeof(unsigned long) * TBL_SIZE;
    volatile fpLM *ptcl                   = (fpLM *) ldm_malloc( bsize_ptcl_chunk );
    volatile int *dest_rank               = (int *)  ldm_malloc( bsize_int_chunk );
    //-(local tables)
    int n_dest_rank_cpe;
    volatile int *dest_rank_list_cpe      = (int *)  ldm_malloc( bsize_int_chunk );
    volatile int *n_send_cpe              = (int *)  ldm_malloc( bsize_int_chunk );
    volatile int *n_disp_send_cpe         = (int *)  ldm_malloc( bsize_int_chunk );
    //-(temporal buffers)
    volatile fpLM *ptcl_put_buf           = (fpLM *) ldm_malloc( bsize_ptcl_chunk );
    volatile int *n_disp_put_buf          = (int *)  ldm_malloc( bsize_int_chunk );
    //-(global tables)
    volatile int *dest_rank_list          = (int *) ldm_malloc( bsize_int_tbl );
    volatile int *n_disp_send             = (int *) ldm_malloc( bsize_int_tbl );
    volatile unsigned long *adr_send_bufs = (unsigned long *) ldm_malloc( bsize_long_tbl );
#ifdef CHECK_MAKE_SEND_BUFS
    if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
        printf("Total memory = %d [B]\n",
               2*bsize_ptcl_chunk + 5*bsize_int_chunk + 2*bsize_int_tbl + bsize_long_tbl);
    }
#endif

    //* Get the table data
    // n_dest_rank
    assert((0 < n_dest_rank) && (n_dest_rank <= TBL_SIZE)); 
    // Because myrank is also counted as a destination rank,
    // n_dest_rank must be larger than 0.
#ifdef CHECK_MAKE_SEND_BUFS
    sync_array_();
    if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
        printf("n_dest_rank = %d\n",n_dest_rank);
    }
    sync_array_();
#endif
    // dest_rank_list[] & adr_send_bufs[]
    if (my_id == 0) {
        // DMA get for dest_rank_list
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get,  sizeof(int) * n_dest_rank);
        dma(dma_get, (long*)((int *)adr_dest_rank_list), (long*)(dest_rank_list));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        // DMA get conf. for adr_send_bufs[]
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get,  sizeof(unsigned long));
        // DMA get for adr_send_bufs[]
        for (i=0; i<n_dest_rank; i++) {
            rank = dest_rank_list[i];
            reply_get = 0;
            dma(dma_get, (long*)((unsigned long *)adr_adr_ptcl_send_buf_ + rank), (long*)(&adr_send_bufs[i]));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
        }
    }
    sync_array_();
    // Share info.
    for (i=0; i<n_dest_rank; i++) {
        cpe_bcast_int32(0,&dest_rank_list[i]); 
        cpe_bcast_uint64(0,&adr_send_bufs[i]); 
    }
    sync_array_();
    // Initialize n_disp_send
    for (i=0; i<n_dest_rank; i++) n_disp_send[i] = 0;
#ifdef CHECK_MAKE_SEND_BUFS
    sync_array_();
    if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
        printf("DMA get of the tables is completed!\n");
    }
    sync_array_();
#endif

    //* Pack to the send buffers
    for (i=0; i<nptcl_loc; i+=CYCLE_SIZE) {
#ifdef CHECK_MAKE_SEND_BUFS
        sync_array_();
        if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
            printf("i = %d\n",i);
        }
        sync_array_();
#endif
        int nrem = nptcl_loc - i;
        int nn   = nrem < CYCLE_SIZE ? nrem : CYCLE_SIZE;
        //* Compute the task of each CPE
        int n_loc,my_offset;
        n_loc = nn/NUMBER_OF_CPE + ( (my_id < nn % NUMBER_OF_CPE) ? 1 : 0 );
        my_offset = i + (nn/NUMBER_OF_CPE)*my_id + ( (my_id < nn % NUMBER_OF_CPE) ? my_id : nn % NUMBER_OF_CPE );
#ifdef CHECK_MAKE_SEND_BUFS
        if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
            printf("i         = %d\n",i);
            printf("n_loc     = %d\n",n_loc);
            printf("my_offset = %d\n",my_offset);
        }
#endif
        if (n_loc > 0) {
            //** Get ptcl_ 
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get,  sizeof(fpLM) * n_loc);
            dma(dma_get, (long*)((fpLM *)adr_ptcl_ + my_offset), (long*)(ptcl));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            //** Get dest_rank[]
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get,  sizeof(int) * n_loc);
            dma(dma_get, (long*)((int *)adr_dest_rank + my_offset), (long*)(dest_rank));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            //** Make the list of destination ranks for this CPE and
            //   count the number of sent particles for each process
            n_dest_rank_cpe = 1;
            dest_rank_list_cpe[0] = dest_rank[0];
            n_send_cpe[0] = 1;
            for (k=1; k<n_loc; k++) {
                jloc = -1;
                for (j=0; j<n_dest_rank_cpe; j++) {
                    if (dest_rank[k] == dest_rank_list_cpe[j]) {
                        jloc = j;
                        n_send_cpe[j]++;
                        break;
                    }
                }
                if (jloc == -1) {
                    // A new rank is found.
                    dest_rank_list_cpe[n_dest_rank_cpe] = dest_rank[k];
                    n_send_cpe[n_dest_rank_cpe] = 1;
                    n_dest_rank_cpe++;
                }
            }
        } else {
            // In this case, this CPE does not perform DMA put.
            n_dest_rank_cpe = 0;
        }
        //** Compute the prefix sum (Collective CPE communication)
        for (k=0; k<n_dest_rank; k++) {
            rank = dest_rank_list[k];
            // Check if I have rank as a destination rank
            jloc = -1;
            for (j=0; j<n_dest_rank_cpe; j++) {
                if (rank == dest_rank_list_cpe[j]) {
                    jloc = j;
                    val = n_send_cpe[j];
                    break;
                }
            }
            if (jloc == -1) val = 0;
            // Prefix sum calc.
            prefix_sum(val,&beg,&end,&psum);
            if (jloc != -1) {
                n_disp_send_cpe[jloc] = n_disp_send[k] + beg;
            } 
            // Update n_disp_send[k] for the next cycle
            n_disp_send[k] += psum;
        }
        //** Set the send buffers 
        //if (my_id == 0) {// for debug
        //if (i == 0) {// for debug
        if (n_loc > 0) {
#if 0
            // Reset n_send_cpe[]
            for (j=0; j<n_dest_rank_cpe; j++) 
                n_send_cpe[j] = 0;
            // DMA put
            for (k=0; k<n_loc; k++) {
                // Find the destination rank
                jglb = -1;
                for (j=0; j<n_dest_rank; j++) {
                    if (dest_rank[k] == dest_rank_list[j]) {
                        jglb = j;
                        break;
                    }
                }
                //[DEBUG]
                //if (jglb == -1) {
                //    printf("i = %d\n",i);
                //    printf("k = %d\n",k);
                //    printf("dest_rank[k] = %d\n",dest_rank[k]);
                //    for (j=0; j<n_dest_rank; j++) 
                //        printf("j = %d, dest_rank_list[j] = %d\n",j,dest_rank_list[j]);
                //}
                assert(jglb >= 0);
                jloc = -1;
                for (j=0; j<n_dest_rank_cpe; j++) {
                    if (dest_rank[k] == dest_rank_list_cpe[j]) {
                        jloc = j;
                        break;
                    }
                }
                assert(jloc >= 0);
                // Copy ptcl_[] to ptcl_send_buf_[]
                reply_put = 0;
                dma_set_op(&dma_put,    DMA_PUT);
                dma_set_mode(&dma_put,  PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put,  sizeof(fpLM));
                dma(dma_put, (long*)((fpLM *)adr_send_bufs[jglb] + n_disp_send_cpe[jloc] + n_send_cpe[jloc]), (long*)(&ptcl[k]));
                dma_wait(&reply_put, 1);
                while (reply_put != 1) {}
                // Update n_send_cpe for the next DMA put
                n_send_cpe[jloc]++;
            }
#else
            // Compute n_disp_put_buf[]
            n_disp_put_buf[0] = 0;
            for (j=1; j<n_dest_rank_cpe; j++)
                n_disp_put_buf[j] = n_disp_put_buf[j-1] + n_send_cpe[j-1];
            // Reset n_send_cpe[]
            for (j=0; j<n_dest_rank_cpe; j++) n_send_cpe[j] = 0;
            // Pack ptcl[] to local buffers 
            for (k=0; k<n_loc; k++) {
                jloc = -1;
                for (j=0; j<n_dest_rank_cpe; j++) {
                    if (dest_rank[k] == dest_rank_list_cpe[j]) {
                        jloc = j;
                        break;
                    }
                }
                assert(jloc >= 0);
                adr = n_disp_put_buf[jloc] + n_send_cpe[jloc];
                ptcl_put_buf[adr] = ptcl[k];
                n_send_cpe[jloc]++;
            }
            // Copy ptcl_[] to ptcl_send_buf_[]
            for (k=0; k<n_dest_rank_cpe; k++) {
                jglb = -1;
                for (j=0; j<n_dest_rank; j++) {
                    if (dest_rank_list_cpe[k] == dest_rank_list[j]) {
                        jglb = j;
                        break;
                    }
                }
                if(jglb < 0){
                    printf("MY_RANK_MPI= %d HOST_NAME= %s \n", MY_RANK_MPI, HOST_NAME[MY_RANK_MPI]);
                }
                assert(jglb >= 0);
                adr = n_disp_put_buf[k];
                reply_put = 0;
                dma_set_op(&dma_put,    DMA_PUT);
                dma_set_mode(&dma_put,  PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put,  sizeof(fpLM) * n_send_cpe[k]);
                dma(dma_put, (long*)((fpLM *)adr_send_bufs[jglb] + n_disp_send_cpe[k]), (long*)(&ptcl_put_buf[adr]));
                dma_wait(&reply_put, 1);
                while (reply_put != 1) {}
            }
#endif
        }
        //}// for debug
        //}// for debug
    }

#ifdef CHECK_MAKE_SEND_BUFS
    sync_array_();
    if ((myrank == RANK_CHECK) && (my_id == ID_CHECK)) {
        printf("Making of the send buffers is completed!\n");
    }
    sync_array_();
#endif

    //* Release memory
    ldm_free(ptcl,               bsize_ptcl_chunk);
    ldm_free(dest_rank,          bsize_int_chunk);
    ldm_free(dest_rank_list_cpe, bsize_int_chunk);
    ldm_free(n_send_cpe,         bsize_int_chunk);
    ldm_free(n_disp_send_cpe,    bsize_int_chunk);

    ldm_free(ptcl_put_buf, bsize_ptcl_chunk);
    ldm_free(n_disp_put_buf, bsize_int_chunk);

    ldm_free(dest_rank_list,     bsize_int_tbl);
    ldm_free(n_disp_send,        bsize_int_tbl);
    ldm_free(adr_send_bufs,      bsize_long_tbl);
}
