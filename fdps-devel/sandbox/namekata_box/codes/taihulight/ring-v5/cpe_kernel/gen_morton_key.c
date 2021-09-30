/* Standard C headers */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/* Headers specific to Sunway TaihuLight */
#include "slave.h"
#include "dma.h"
#include "cpe_func.h"

extern double __ieee754_atan2 __P((double, double, double *, double *, double * ));

/* Macro function */
#define separateBit(s_in,s_out) \
    do { \
    (s_in) = ((s_in) | (s_in)<<32) & 0xffff00000000ffffL; \
    (s_in) = ((s_in) | (s_in)<<16) & 0x00ff0000ff0000ffL; \
    (s_in) = ((s_in) | (s_in)<<8)  & 0xf00f00f00f00f00fL; \
    (s_in) = ((s_in) | (s_in)<<4)  & 0x30c30c30c30c30c3L; \
    (s_out) = ((s_in) | (s_in)<<2) & 0x9249249249249249L; \
    } while(0)
// 11111111 11111111 00000000 00000000 00000000 00000000 11111111 11111111
// 00000000 11111111 00000000 00000000 11111111 00000000 00000000 11111111
// 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111
// 11 00 00 11 00 11 00 00 11
// 1 0 0 1 0 0 1 0 0 1 0 0 1

#define separateBit32(s_in,s_out) \
    do { \
    (s_in) = ((s_in) | (s_in)<<16) & 0xff0000ff; \
    (s_in) = ((s_in) | (s_in)<<8)  & 0x0f00f00f; \
    (s_in) = ((s_in) | (s_in)<<4)  & 0xc30c30c3; \
    (s_out) = ((s_in) | (s_in)<<2) & 0x49249249; \
    } while(0)


#ifdef KEY_3D
__attribute__((noinline))
void conv_mtn(int N, int ioffset, epjLM epj[], tpLM tp[],
              double *length_, F64vec_ *center_,
              double * ah, double * al, double * at){
    double MY_PI = 3.14159265358979323846;
    double MY_PI_2 = 2.0*MY_PI;
    int i;
#ifdef USE_96BIT_KEY
    U64_ kLevMaxUp = 21;
    U64_ kLevMaxLo = 10;
    U64_ kLevMax = 31;
#else
    int kLevMax = 21;
#endif
    double half_len_ = 0.5 * (*length_);
    double normalized_factor_ = (1.0/(half_len_*2.0))*(((U64_)1)<<kLevMax);
    double cx = center_->x - half_len_;
    double cy = center_->y - half_len_;
    double cz = center_->z - half_len_;
    
    for (i=0; i<N; i++){
        double pos_x = epj[i].pos.x;
        double pos_y = epj[i].pos.y;
        double pos_z = epj[i].pos.z;
        double pos_phi = __ieee754_atan2(pos_y, pos_x, ah, al, at);
        if(pos_phi < 0.0) pos_phi += MY_PI_2;
        double pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
        U64_ nx, ny, nz, retx, rety, retz;
        nx = (U64_)( (pos_phi - cx) * normalized_factor_);
        ny = (U64_)( (pos_r - cy) * normalized_factor_);
        nz = (U64_)( (pos_z - cz) * normalized_factor_);
#ifdef USE_96BIT_KEY
        U64_ nx_up = nx>>kLevMaxLo;
        U64_ ny_up = ny>>kLevMaxLo;
        U64_ nz_up = nz>>kLevMaxLo;
        separateBit(nx_up, retx);
        separateBit(ny_up, rety);
        separateBit(nz_up, retz);
        tp[i].key_hi_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
        U32_ nx_lo = (U32_)(nx & 0x3ff);
        U32_ ny_lo = (U32_)(ny & 0x3ff);
        U32_ nz_lo = (U32_)(nz & 0x3ff);
        U64_ retx_lo, rety_lo, retz_lo;
        separateBit32(nx_lo, retx_lo);
        separateBit32(ny_lo, rety_lo);
        separateBit32(nz_lo, retz_lo);
        tp[i].key_lo_ = ( retx_lo<<2 | rety_lo<<1 | retz_lo ); // the compile option must be -O0 or -O1.
#else
        separateBit(nx, retx);
        separateBit(ny, rety);
        separateBit(nz, retz);
        tp[i].key_hi_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
#endif
        tp[i].adr_ptcl_ = ioffset + i;
    }
}
#else
__attribute__((noinline))
void conv_mtn(int N, int ioffset, epjLM epj[], tpLM tp[], double *length_, F64vec_ *center_,
              double * ah, double * al, double * at){
    double MY_PI = 3.14159265358979323846;
    double MY_PI_2 = 2.0*MY_PI;
    int i;
    int kLevMax = 21;
    double half_len_ = 0.5 * (*length_);
    double normalized_factor_ = (1.0/(half_len_*2.0))*(1<<kLevMax);
    double cx = center_->x - half_len_;
    double cy = center_->y - half_len_;
    for (i=0; i<N; i++){
        double pos_x = epj[i].pos.x;
        double pos_y = epj[i].pos.y;
        double pos_phi = __ieee754_atan2(pos_y, pos_x, ah, al, at);
        if(pos_phi < 0.0) pos_phi += MY_PI_2;
        double pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
        U64_ nx, ny, nz, retx, rety, retz;
        nx = (U64_)( (pos_phi - cx) * normalized_factor_);
        ny = (U64_)( (pos_r   - cy) * normalized_factor_);
        separateBit(nx, retx);
        separateBit(ny, rety);
        tp[i].key_ = ( retx<<2 | rety<<1 ); // the compile option must be -O0 or -O1.
        tp[i].adr_ptcl_ = ioffset + i;
    }
}
#endif

void GenMortonKey(void *args) {
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int n_loc_tot_ = (int)(((unsigned long*)args)[0]);
    void *adr_length_ = (void *)((unsigned long*)args)[1];
    void *adr_center_ = (void *)((unsigned long*)args)[2];
    void *adr_tp_head = (void *)((unsigned long*)args)[3];
    void *adr_epj_head = (void *)((unsigned long*)args)[4];
    int myrank = (int)(((unsigned long*)args)[5]);
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_loc_tot_/NUMBER_OF_CPE + ( (my_id < n_loc_tot_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_loc_tot_/NUMBER_OF_CPE)*my_id + ( (my_id < n_loc_tot_ % NUMBER_OF_CPE) ? my_id : n_loc_tot_ % NUMBER_OF_CPE );
#if 0
    //* Data Check (1)
    if (my_id == 63) {
        printf("--------\n");
        printf("my_id         = %d\n",my_id);
        printf("n_loc_tot_    = %d\n",n_loc_tot_);
        printf("n_loc         = %d\n",n_loc);
        printf("my_offset     = %d\n",my_offset);
    }
#endif
    enum {
        CHUNK_SIZE = 64,
    };
    //* Memory allocation of local buffers
    size_t bsize;
    size_t bsize_tp  = sizeof(tpLM);
    size_t bsize_epj = sizeof(epjLM);
    size_t bsize_length_ = sizeof(double);
    size_t bsize_center_ = sizeof(F64vec_);
    volatile tpLM *tp = (tpLM *)(long)ldm_malloc( bsize_tp * CHUNK_SIZE);
    volatile epjLM *epj = (epjLM *)(long)ldm_malloc( bsize_epj * CHUNK_SIZE);
    volatile double *length_ = (double *)(long)ldm_malloc( bsize_length_ );
    volatile F64vec_ *center_ = (F64vec_ *)(long)ldm_malloc( bsize_center_ );

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
    
#if 0
    //* Data Check (2)
    if (my_id == 63) {
        printf("bsize_tp      = %d\n",bsize_tp);
        printf("bsize_epj     = %d\n",bsize_epj);
        printf("bsize_length_ = %d\n",bsize_length_);
        printf("bsize_center_ = %d\n",bsize_center_);
        printf("bsize_U64_    = %d\n",sizeof(U64_));
    }
#endif

    //if (myrank == 3) tp[CHUNK_SIZE+16384].key_ = 1; // To test
    //* Compute morton keys on each CPE
    // (loop counters)
    int i,j,k;
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (local vars. for morton key calculation)
    int kLevMax = 21;
    F64_ half_len_;
    F64_ normalized_factor_;
    U64_ nx,ny,nz,retx,rety,retz;
    //** Get length_
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, bsize_length_);
    dma(dma_get, (long)((double *)adr_length_), (long)(length_));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
    //** Get center_
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, bsize_center_);
    dma(dma_get, (long)((F64vec_ *)adr_center_), (long)(center_));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
    //* Compute half_len_ & normalized_factor_
    half_len_ = 0.5 * (*length_);
    normalized_factor_ = (1.0/(half_len_*2.0))*(1<<kLevMax);
#if 0
    //* Data Check (3)
    if (my_id == 63) {
        printf("length_       = %e\n",*length_);
        printf("center_->x    = %e\n",center_->x);
        printf("center_->y    = %e\n",center_->y);
        printf("center_->z    = %e\n",center_->z);
        printf("norm.fac.     = %e\n",normalized_factor_);
    }
#endif
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn  = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //* Get epj
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_epj * nn);
        dma(dma_get, (long)((epjLM *)adr_epj_head + my_offset + i), (long)(epj));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        conv_mtn(nn, my_offset+i, epj, tp, length_, center_, atanhi, atanlo, aT);

        //* Put tp
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tp * nn);
        dma(dma_put, (long)((tpLM *)adr_tp_head + my_offset + i), (long)(tp));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }

    //* Release memory
    ldm_free(tp, bsize_tp * CHUNK_SIZE);
    ldm_free(epj, bsize_epj * CHUNK_SIZE);
    ldm_free(length_, bsize_length_);
    ldm_free(center_, bsize_center_);

    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);    
}
