/* Standard C headers */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/* Headers specific to Sunway TaihuLight */
#include "slave.h"
#include "dma.h"
#include "cpe_func.h"

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

__attribute__((noinline))
void conv_mtn(int N, int ioffset, epjLM epj[], tpLM tp[], double *length_, F64vec_ *center_){
    int i;
    int kLevMax = 21;
    double half_len_ = 0.5 * (*length_);
    double normalized_factor_ = (1.0/(half_len_*2.0))*(1<<kLevMax);
    double cx = center_->x - half_len_;
    double cy = center_->y - half_len_;
    double cz = center_->z - half_len_;

    for (i=0; i<N; i++){
        U64_ nx, ny, nz, retx, rety, retz;
        nx = (U64_)( (epj[i].pos.x - cx) * normalized_factor_);
        ny = (U64_)( (epj[i].pos.y - cy) * normalized_factor_);
        nz = (U64_)( (epj[i].pos.z - cz) * normalized_factor_);
        separateBit(nx, retx);
        separateBit(ny, rety);
        separateBit(nz, retz);
        tp[i].key_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
        tp[i].adr_ptcl_ = ioffset + i;
    }
}

void GenMortonKey(void *args) {
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int n_loc_tot_ = (int)(((unsigned long*)args)[0]);
    void *adr_length_ = (void *)((unsigned long*)args)[1];
    void *adr_center_ = (void *)((unsigned long*)args)[2];
    void *adr_tp_head = (void *)((unsigned long*)args)[3];
    void *adr_epj_head = (void *)((unsigned long*)args)[4];
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
    tpLM *tp = (tpLM *)(long)ldm_malloc( bsize_tp * CHUNK_SIZE);
    epjLM *epj = (epjLM *)(long)ldm_malloc( bsize_epj * CHUNK_SIZE);
    double *length_ = (double *)(long)ldm_malloc( bsize_length_ );
    F64vec_ *center_ = (F64vec_ *)(long)ldm_malloc( bsize_center_ );
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
#if 0
        //* Compute morton key
        nx = (U64_)( (epj->pos.x - center_->x + half_len_) * normalized_factor_);
        ny = (U64_)( (epj->pos.y - center_->y + half_len_) * normalized_factor_);
        nz = (U64_)( (epj->pos.z - center_->z + half_len_) * normalized_factor_);
        separateBit(nx,retx);
        separateBit(ny,rety);
        separateBit(nz,retz);
        tp->key_ = ( retx<<2 | rety<<1 | retz ); // the compile option must be -O0 or -O1.
        tp->adr_ptcl_ = my_offset + i;
#else
        conv_mtn(nn, my_offset+i, epj, tp, length_, center_);
#endif
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

}

