#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
//#include<athread.h>
#include"slave.h"
#include"dma.h"
#include"cpe_func.h"

extern double __ieee754_atan2 __P((double, double, double *, double *, double * ));

void GetMinIpgBox(void * args){
    int my_id  = athread_get_id(-1);
    int n_ipg_ = (int)(((unsigned long*)args)[0]);
    void * adr_ipg_ = (void *)((unsigned long*)args)[1];
    void * adr_epj_ = (void *)((unsigned long*)args)[2];
    int n_ipg_loc = n_ipg_/NUMBER_OF_CPE + ( (my_id < n_ipg_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_ipg_/NUMBER_OF_CPE)*my_id + ( (my_id < n_ipg_ % NUMBER_OF_CPE) ? my_id : n_ipg_ % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 32, // 64 dose not work
    };
    size_t bsize_ipg       = sizeof(ipgLM);
    size_t bsize_epj       = sizeof(epjLM);
    size_t bsize_epj_array = bsize_epj * CHUNK_SIZE;
    ipgLM * ipg = (ipgLM *) ldm_malloc( bsize_ipg );
    epjLM * epj = (epjLM *) ldm_malloc( bsize_epj_array );

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
    
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    for (i=0; i<n_ipg_loc; i++) {
        //** Get ipg
        reply_get = 0;
        dma_set_op(&dma_get,    DMA_GET);
        dma_set_mode(&dma_get,  PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get,  bsize_ipg);
        dma(dma_get, (long*)(((ipgLM *)adr_ipg_) + my_offset + i), (long*)(ipg));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        int n_ptcl   = ipg->n_ptcl_;
        int adr_ptcl = ipg->adr_ptcl_;
        ipg->vertex_.high_.x = ipg->vertex_.high_.y = ipg->vertex_.high_.z = -99999.9;
        ipg->vertex_.low_.x  = ipg->vertex_.low_.y  = ipg->vertex_.low_.z = 99999.9;
        for(j=0; j<n_ptcl; j+=CHUNK_SIZE){
            int nrem = n_ptcl - j;
            int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;

            reply_get = 0;
            dma_set_op(&dma_get,    DMA_GET);
            dma_set_mode(&dma_get,  PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get,  bsize_epj*nn);
            dma(dma_get, (long*)(((epjLM *)adr_epj_)+adr_ptcl+j), (long*)(epj));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}

            for(k=0; k<nn; k++){
                double pos_x = epj[k].pos.x;
                double pos_y = epj[k].pos.y;
                double pos_z = epj[k].pos.z;
                double pos_phi = __ieee754_atan2(pos_y, pos_x, atanhi, atanlo, aT);
                if(pos_phi < 0.0) pos_phi += MY_PI_2;
                double pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                ipg->vertex_.high_.x = (ipg->vertex_.high_.x > pos_phi) ? ipg->vertex_.high_.x : pos_phi;
                ipg->vertex_.high_.y = (ipg->vertex_.high_.y > pos_r)   ? ipg->vertex_.high_.y : pos_r;
                ipg->vertex_.high_.z = (ipg->vertex_.high_.z > pos_z)   ? ipg->vertex_.high_.z : pos_z;
                ipg->vertex_.low_.x  = (ipg->vertex_.low_.x <= pos_phi) ? ipg->vertex_.low_.x : pos_phi;
                ipg->vertex_.low_.y  = (ipg->vertex_.low_.y <= pos_r)   ? ipg->vertex_.low_.y : pos_r;
                ipg->vertex_.low_.z  = (ipg->vertex_.low_.z <= pos_z)   ? ipg->vertex_.low_.z : pos_z;
            }
        }
        
        //** Put ipg
        reply_put = 0;
        dma_set_op(&dma_put,    DMA_PUT);
        dma_set_mode(&dma_put,  PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put,  bsize_ipg);
        dma(dma_put, (long*)(((ipgLM *)adr_ipg_) + my_offset + i), (long*)(ipg));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
        sync_array_(); // to avoid DMA bug
    }
    if(n_ipg_loc == n_ipg_/NUMBER_OF_CPE) sync_array_();
    //* Release memory
    ldm_free(ipg, bsize_ipg);
    ldm_free(epj, bsize_epj_array);
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
}
