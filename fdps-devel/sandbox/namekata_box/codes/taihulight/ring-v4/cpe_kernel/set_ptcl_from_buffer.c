#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include"slave.h"
#include"dma.h"
#include"cpe_func.h"
#include"cpe_prof.h"

extern int MY_RANK_MPI;

static inline U32_ SetMsb(U32_ val){
    return val | 0x80000000;
}

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



void SetEpjTpFromBufferCpe(void * args){

    double MY_PI = 3.14159265358979323846;
    double MY_PI_2 = 2.0*MY_PI;
#ifdef USE_96BIT_KEY
    U64_ kLevMaxUp = 21;
    U64_ kLevMaxLo = 10;
    U64_ kLevMax = 31;
#else
    int kLevMax = 21;
#endif    
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
    
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n_tot = (int)(((unsigned long*)args)[0]);
    int n_disp = (int)(((unsigned long*)args)[1]);
    void * adr_tp = (void *)((unsigned long*)args)[2];
    void * adr_epj_org = (void *)((unsigned long*)args)[3];
    void * adr_epj_buf = (void *)((unsigned long*)args)[4];
    void * adr_center  = (void *)((unsigned long*)args)[5];
    void * adr_hlen    = (void *)((unsigned long*)args)[6];

    F64vec_ center = *(F64vec_*)adr_center;
    F64_    hlen   = *(F64_*)adr_hlen;
    double normalized_factor_ = (1.0/(hlen*2.0))*(((U64_)1)<<kLevMax); 
    F64_ cx = center.x - hlen;
    F64_ cy = center.y - hlen;
    F64_ cz = center.z - hlen;
    
    int bsize_tp = sizeof(tpLM);
    int bsize_ep = sizeof(epjLM);
    
    int my_n = n_tot/NUMBER_OF_CPE + ( (my_id < n_tot % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot % NUMBER_OF_CPE) ? my_id : n_tot % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    epjLM epj_buf[CHUNK_SIZE];
    tpLM  tp[CHUNK_SIZE];
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        // get epj_buf
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_ep * nn);
        dma(dma_get, (long*)((epjLM *)adr_epj_buf+my_offset+i),
            (long*)(epj_buf));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        // put epj_org
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_ep * nn);
        dma(dma_put, (long*)((epjLM *)adr_epj_org+my_offset+i+n_disp),
            (long*)(epj_buf));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
        
        // put tp
        for (j=0; j<nn; j++) {
            F64_ pos_x = epj_buf[j].pos.x;
            F64_ pos_y = epj_buf[j].pos.y;
            F64_ pos_z = epj_buf[j].pos.z;
            F64_ pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
            F64_ pos_phi = __ieee754_atan2(pos_y, pos_x, atanhi, atanlo, aT);
            if(pos_phi < 0.0) pos_phi += MY_PI_2;

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
            tp[j].key_hi_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
            U32_ nx_lo = (U32_)(nx & 0x3ff);
            U32_ ny_lo = (U32_)(ny & 0x3ff);
            U32_ nz_lo = (U32_)(nz & 0x3ff);
            U64_ retx_lo, rety_lo, retz_lo;
            separateBit32(nx_lo, retx_lo);
            separateBit32(ny_lo, rety_lo);
            separateBit32(nz_lo, retz_lo);
            tp[j].key_lo_ = ( retx_lo<<2 | rety_lo<<1 | retz_lo ); // the compile option must be -O0 or -O1.
#else
            separateBit(nx, retx);
            separateBit(ny, rety);
            separateBit(nz, retz);
            tp[j].key_hi_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
            //tp[j].key_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
#endif
            tp[j].adr_ptcl_ = n_disp + my_offset+i+j;
        }

        // put tp_glb
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tp * nn);
        dma(dma_put, (long*)((tpLM *)adr_tp+my_offset+i+n_disp),
            (long*)(tp));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
        
    }
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
}




void SetSpjTpFromBufferCpe(void * args){

    double MY_PI = 3.14159265358979323846;
    double MY_PI_2 = 2.0*MY_PI;
#ifdef USE_96BIT_KEY
    U64_ kLevMaxUp = 21;
    U64_ kLevMaxLo = 10;
    U64_ kLevMax = 31;
#else
    int kLevMax = 21;
#endif    
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n_tot = (int)(((unsigned long*)args)[0]);
    int n_disp = (int)(((unsigned long*)args)[1]);
    void * adr_tp = (void *)((unsigned long*)args)[2];
    void * adr_spj_org = (void *)((unsigned long*)args)[3];
    void * adr_spj_buf = (void *)((unsigned long*)args)[4];
    void * adr_center  = (void *)((unsigned long*)args)[5];
    void * adr_hlen    = (void *)((unsigned long*)args)[6];

    F64vec_ center = *(F64vec_*)adr_center;
    F64_    hlen   = *(F64_*)adr_hlen;
    double normalized_factor_ = (1.0/(hlen*2.0))*(((U64_)1)<<kLevMax); 
    F64_ cx = center.x - hlen;
    F64_ cy = center.y - hlen;
    F64_ cz = center.z - hlen;
    
    int bsize_tp = sizeof(tpLM);
    int bsize_sp = sizeof(spjLM);
    
    int my_n = n_tot/NUMBER_OF_CPE + ( (my_id < n_tot % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot % NUMBER_OF_CPE) ? my_id : n_tot % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    spjLM spj_buf[CHUNK_SIZE];
    tpLM  tp[CHUNK_SIZE];
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        // get spj_buf
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_sp * nn);
        dma(dma_get, (long*)((spjLM *)adr_spj_buf+my_offset+i),
            (long*)(spj_buf));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        // put spj_org
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_sp * nn);
#ifdef REDUCE_MEMORY
        dma(dma_put, (long*)((spjLM *)adr_spj_org+my_offset+i),
            (long*)(spj_buf));
#else        
        dma(dma_put, (long*)((spjLM *)adr_spj_org+my_offset+i+n_disp),
            (long*)(spj_buf));
#endif
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}

        // put tp
        for (j=0; j<nn; j++) {
            F64_ pos_phi = spj_buf[j].pos_phi;
            F64_ pos_r   = spj_buf[j].pos_r;
            F64_ pos_z   = spj_buf[j].pos.z;
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
            tp[j].key_hi_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
            U32_ nx_lo = (U32_)(nx & 0x3ff);
            U32_ ny_lo = (U32_)(ny & 0x3ff);
            U32_ nz_lo = (U32_)(nz & 0x3ff);
            U64_ retx_lo, rety_lo, retz_lo;
            separateBit32(nx_lo, retx_lo);
            separateBit32(ny_lo, rety_lo);
            separateBit32(nz_lo, retz_lo);
            tp[j].key_lo_ = ( retx_lo<<2 | rety_lo<<1 | retz_lo ); // the compile option must be -O0 or -O1.
#else
            separateBit(nx, retx);
            separateBit(ny, rety);
            separateBit(nz, retz);
            tp[j].key_hi_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
            //tp[j].key_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
#endif
#ifdef REDUCE_MEMORY
            tp[j].adr_ptcl_ = SetMsb( (U32_)(my_offset+i+j) );
#else
            tp[j].adr_ptcl_ = SetMsb( (U32_)(n_disp+my_offset+i+j) );
#endif
        }

        // put tp_glb
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tp * nn);
        dma(dma_put, (long*)((tpLM *)adr_tp+my_offset+i+n_disp),
            (long*)(tp));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}

    }
}
