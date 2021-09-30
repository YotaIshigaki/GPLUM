/* Standard C headers */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
/* Headers specific to Sunway TaihuLight */
#include "slave.h"
#include "dma.h"
#include "cpe_func.h"

void CalcMoment(void *args) {
    int my_id = athread_get_id(-1);
    //* Process the arguments
    void *adr_epj_head = (void *)((unsigned long*)args)[0];
    void *adr_tc_head  = (void *)((unsigned long*)args)[1];
    int n_leaf_limit = (int)(((unsigned long*)args)[2]);
    int head = (int)(((unsigned long*)args)[3]);
    int next = (int)(((unsigned long*)args)[4]);
    void * adr_rsrch = (void *)((unsigned long*)args)[5];
    //* Compute the task of each CPE
    int num_tc,num_tc_loc,my_offset;
    num_tc = next - head;
    num_tc_loc = num_tc/NUMBER_OF_CPE + ( (my_id < num_tc % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (num_tc/NUMBER_OF_CPE)*my_id + ( (my_id < num_tc % NUMBER_OF_CPE) ? my_id : num_tc % NUMBER_OF_CPE );
#if 0
    //* Data Check (1)
    if (my_id == 63) {
        printf("--------\n");
        printf("my_id        = %d\n",my_id);
        printf("n_leaf_limit = %d\n",n_leaf_limit);
        printf("head         = %d\n",head);
        printf("next         = %d\n",next);
        printf("num_tc       = %d\n",num_tc);
        printf("num_tc_loc   = %d\n",num_tc_loc);
        printf("my_offset    = %d\n",my_offset);
    }
#endif
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tc  = sizeof(tcLM);
    size_t bsize_epj = sizeof(epjLM);
    size_t bsize_epj_array = bsize_epj * CHUNK_SIZE;
    size_t bsize_tc_array  = bsize_tc  * N_CHILDREN;
    size_t bsize_rsrch = sizeof(double);
    tcLM *tc_tgt = (tcLM *) ldm_malloc( bsize_tc );
    epjLM *epj = (epjLM *) ldm_malloc( bsize_epj_array );
    tcLM *tc = (tcLM *) ldm_malloc( bsize_tc_array );
    double *rsrch = (double *) ldm_malloc( bsize_rsrch );
#if 0
    //* Data Check (2)
    if (my_id == 63) {
        printf("bsize_tc        = %d\n",bsize_tc);
        printf("bsize_epj_array = %d\n",bsize_epj_array);
        printf("bsize_tc_array  = %d\n",bsize_tc_array);
    }
#endif
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (local vars. for moment calculation)
    F64_ mass;
    F64vec_ pos;
    double cx,cy,cz;
    //** Get the search radius
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, bsize_rsrch);
    dma(dma_get, (long*)((double *)adr_rsrch), (long*)(rsrch));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
#if 0
    if (my_id == 63) {
        printf("rsrch        = %e\n",*rsrch);
    }
#endif
    for (i=0; i<num_tc_loc; i++) {
        //* Get the target tc
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tc);
        dma(dma_get, (long*)((tcLM *)adr_tc_head + head + my_offset + i), (long*)(tc_tgt));
        dma_wait(&reply_get, 1);
        while (reply_get !=1 ) {}

        //* Moment calculation
        tc_tgt->mass = 0.0;
        tc_tgt->pos.x = 0.0;
        tc_tgt->pos.y = 0.0;
        tc_tgt->pos.z = 0.0;
#ifdef USE_QUADRUPOLE
        tc_tgt->quad.xx = 0.0;
        tc_tgt->quad.yy = 0.0;
        tc_tgt->quad.zz = 0.0;
        tc_tgt->quad.xy = 0.0;
        tc_tgt->quad.xz = 0.0;
        tc_tgt->quad.yz = 0.0;
#endif
        tc_tgt->vertex_out_.low_.x = 9999.9;
        tc_tgt->vertex_out_.low_.y = 9999.9;
        tc_tgt->vertex_out_.low_.z = 9999.9;
        tc_tgt->vertex_out_.high_.x = -9999.9;
        tc_tgt->vertex_out_.high_.y = -9999.9;
        tc_tgt->vertex_out_.high_.z = -9999.9;
        tc_tgt->vertex_in_.low_.x = 9999.9;
        tc_tgt->vertex_in_.low_.y = 9999.9;
        tc_tgt->vertex_in_.low_.z = 9999.9;
        tc_tgt->vertex_in_.high_.x = -9999.9;
        tc_tgt->vertex_in_.high_.y = -9999.9;
        tc_tgt->vertex_in_.high_.z = -9999.9;

        if (tc_tgt->n_ptcl_ == 0) {
            continue;
        } else if (tc_tgt->n_ptcl_ <= n_leaf_limit || tc_tgt->level_ == TREE_LEVEL_LIMIT){
            //* Inline expansion of tc_tmp->mom_.accumulateAtLeaf(epj[k]);
            for (j=0; j<tc_tgt->n_ptcl_; j+= CHUNK_SIZE) {
                int nrem = tc_tgt->n_ptcl_ - j;
                int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
                //** Get epj
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_epj * nn);
                dma(dma_get, (long*)((epjLM *)adr_epj_head + tc_tgt->adr_ptcl_ + j), (long*)(epj));
                dma_wait(&reply_get, 1);
                while (reply_get !=1 ) {}
                //** Moment calculation
                for (k=0; k<nn; k++){
                   tc_tgt->mass += epj[k].mass;
                   tc_tgt->pos.x += epj[k].mass * epj[k].pos.x;
                   tc_tgt->pos.y += epj[k].mass * epj[k].pos.y;
                   tc_tgt->pos.z += epj[k].mass * epj[k].pos.z;

                   tc_tgt->vertex_out_.low_.x = ( tc_tgt->vertex_out_.low_.x <= epj[k].pos.x - *rsrch) ? tc_tgt->vertex_out_.low_.x : epj[k].pos.x - *rsrch;
                   tc_tgt->vertex_out_.low_.y = ( tc_tgt->vertex_out_.low_.y <= epj[k].pos.y - *rsrch) ? tc_tgt->vertex_out_.low_.y : epj[k].pos.y - *rsrch;
                   tc_tgt->vertex_out_.low_.z = ( tc_tgt->vertex_out_.low_.z <= epj[k].pos.z - *rsrch) ? tc_tgt->vertex_out_.low_.z : epj[k].pos.z - *rsrch;
                   tc_tgt->vertex_out_.high_.x = ( tc_tgt->vertex_out_.high_.x > epj[k].pos.x + *rsrch) ? tc_tgt->vertex_out_.high_.x : epj[k].pos.x + *rsrch;
                   tc_tgt->vertex_out_.high_.y = ( tc_tgt->vertex_out_.high_.y > epj[k].pos.y + *rsrch) ? tc_tgt->vertex_out_.high_.y : epj[k].pos.y + *rsrch;
                   tc_tgt->vertex_out_.high_.z = ( tc_tgt->vertex_out_.high_.z > epj[k].pos.z + *rsrch) ? tc_tgt->vertex_out_.high_.z : epj[k].pos.z + *rsrch;

                   tc_tgt->vertex_in_.low_.x = ( tc_tgt->vertex_in_.low_.x <= epj[k].pos.x ) ? tc_tgt->vertex_in_.low_.x : epj[k].pos.x;
                   tc_tgt->vertex_in_.low_.y = ( tc_tgt->vertex_in_.low_.y <= epj[k].pos.y ) ? tc_tgt->vertex_in_.low_.y : epj[k].pos.y;
                   tc_tgt->vertex_in_.low_.z = ( tc_tgt->vertex_in_.low_.z <= epj[k].pos.z ) ? tc_tgt->vertex_in_.low_.z : epj[k].pos.z;
                   tc_tgt->vertex_in_.high_.x = ( tc_tgt->vertex_in_.high_.x > epj[k].pos.x ) ? tc_tgt->vertex_in_.high_.x : epj[k].pos.x;
                   tc_tgt->vertex_in_.high_.y = ( tc_tgt->vertex_in_.high_.y > epj[k].pos.y ) ? tc_tgt->vertex_in_.high_.y : epj[k].pos.y;
                   tc_tgt->vertex_in_.high_.z = ( tc_tgt->vertex_in_.high_.z > epj[k].pos.z ) ? tc_tgt->vertex_in_.high_.z : epj[k].pos.z;
                }
            }
            //** Inline expansion of tc_tmp->mom_.set();
            tc_tgt->pos.x /= tc_tgt->mass;
            tc_tgt->pos.y /= tc_tgt->mass;
            tc_tgt->pos.z /= tc_tgt->mass;
#ifdef USE_QUADRUPOLE
            //** Inline expansion of tc_tmp->mom_.accumulateAtLeaf2(epj[k]);
            for (j=0; j<tc_tgt->n_ptcl_; j+= CHUNK_SIZE) {
                int nrem = tc_tgt->n_ptcl_ - j;
                int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
                //* Get epj
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_epj * nn);
                dma(dma_get, (long*)((epjLM *)adr_epj_head + tc_tgt->adr_ptcl_ + j), (long*)(epj));
                dma_wait(&reply_get, 1);
                while (reply_get !=1 ) {}
                //** Moment calculation
                for (k=0; k<nn; k++) {
                    mass = epj[k].mass;
                    pos.x = epj[k].pos.x - tc_tgt->pos.x;
                    pos.y = epj[k].pos.y - tc_tgt->pos.y;
                    pos.z = epj[k].pos.z - tc_tgt->pos.z;
                    cx = mass * pos.x;
                    cy = mass * pos.y;
                    cz = mass * pos.z;
                    tc_tgt->quad.xx += cx * pos.x;
                    tc_tgt->quad.yy += cy * pos.y;
                    tc_tgt->quad.zz += cz * pos.z;
                    tc_tgt->quad.xy += cx * pos.y;
                    tc_tgt->quad.xz += cx * pos.z;
                    tc_tgt->quad.yz += cy * pos.z;
                }
            }
#endif
        } else {
            //** Get tc (childnodes)
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tc * N_CHILDREN);
            dma(dma_get, (long*)((tcLM *)adr_tc_head + tc_tgt->adr_tc_), (long*)(tc));
            dma_wait(&reply_get, 1);
            while (reply_get !=1 ) {}
            //** Inline expansion of tc_tmp->mom_.accumulate( tc_tmp_tmp->mom_ );
            for (k=0; k<N_CHILDREN; k++){
                if (tc[k].n_ptcl_ == 0) continue;
                tc_tgt->mass += tc[k].mass;
                tc_tgt->pos.x += tc[k].mass * tc[k].pos.x;
                tc_tgt->pos.y += tc[k].mass * tc[k].pos.y;
                tc_tgt->pos.z += tc[k].mass * tc[k].pos.z;

                tc_tgt->vertex_out_.low_.x = ( tc_tgt->vertex_out_.low_.x <= tc[k].vertex_out_.low_.x ) ? tc_tgt->vertex_out_.low_.x : tc[k].vertex_out_.low_.x;
                tc_tgt->vertex_out_.low_.y = ( tc_tgt->vertex_out_.low_.y <= tc[k].vertex_out_.low_.y ) ? tc_tgt->vertex_out_.low_.y : tc[k].vertex_out_.low_.y;
                tc_tgt->vertex_out_.low_.z = ( tc_tgt->vertex_out_.low_.z <= tc[k].vertex_out_.low_.z ) ? tc_tgt->vertex_out_.low_.z : tc[k].vertex_out_.low_.z;
                tc_tgt->vertex_out_.high_.x = ( tc_tgt->vertex_out_.high_.x > tc[k].vertex_out_.high_.x ) ? tc_tgt->vertex_out_.high_.x : tc[k].vertex_out_.high_.x;
                tc_tgt->vertex_out_.high_.y = ( tc_tgt->vertex_out_.high_.y > tc[k].vertex_out_.high_.y ) ? tc_tgt->vertex_out_.high_.y : tc[k].vertex_out_.high_.y;
                tc_tgt->vertex_out_.high_.z = ( tc_tgt->vertex_out_.high_.z > tc[k].vertex_out_.high_.z ) ? tc_tgt->vertex_out_.high_.z : tc[k].vertex_out_.high_.z;

                tc_tgt->vertex_in_.low_.x = ( tc_tgt->vertex_in_.low_.x <= tc[k].vertex_in_.low_.x ) ? tc_tgt->vertex_in_.low_.x : tc[k].vertex_in_.low_.x;
                tc_tgt->vertex_in_.low_.y = ( tc_tgt->vertex_in_.low_.y <= tc[k].vertex_in_.low_.y ) ? tc_tgt->vertex_in_.low_.y : tc[k].vertex_in_.low_.y;
                tc_tgt->vertex_in_.low_.z = ( tc_tgt->vertex_in_.low_.z <= tc[k].vertex_in_.low_.z ) ? tc_tgt->vertex_in_.low_.z : tc[k].vertex_in_.low_.z;
                tc_tgt->vertex_in_.high_.x = ( tc_tgt->vertex_in_.high_.x > tc[k].vertex_in_.high_.x ) ? tc_tgt->vertex_in_.high_.x : tc[k].vertex_in_.high_.x;
                tc_tgt->vertex_in_.high_.y = ( tc_tgt->vertex_in_.high_.y > tc[k].vertex_in_.high_.y ) ? tc_tgt->vertex_in_.high_.y : tc[k].vertex_in_.high_.y;
                tc_tgt->vertex_in_.high_.z = ( tc_tgt->vertex_in_.high_.z > tc[k].vertex_in_.high_.z ) ? tc_tgt->vertex_in_.high_.z : tc[k].vertex_in_.high_.z;
            }
            //** Inline expansion of tc_tmp->mom_.set();
            tc_tgt->pos.x /= tc_tgt->mass;
            tc_tgt->pos.y /= tc_tgt->mass;
            tc_tgt->pos.z /= tc_tgt->mass;
#ifdef USE_QUADRUPOLE
            //** Inline expansion of tc_tmp->mom_.accumulate2( tc_tmp_tmp->mom_ );
            for (k=0; k<N_CHILDREN; k++){
                if (tc[k].n_ptcl_ == 0) continue;
                mass = tc[k].mass;
                pos.x = tc[k].pos.x - tc_tgt->pos.x;
                pos.y = tc[k].pos.y - tc_tgt->pos.y;
                pos.z = tc[k].pos.z - tc_tgt->pos.z;
                cx = mass * pos.x;
                cy = mass * pos.y;
                cz = mass * pos.z;
                tc_tgt->quad.xx += cx * pos.x + tc[k].quad.xx;
                tc_tgt->quad.yy += cy * pos.y + tc[k].quad.yy;
                tc_tgt->quad.zz += cz * pos.z + tc[k].quad.zz;
                tc_tgt->quad.xy += cx * pos.y + tc[k].quad.xy;
                tc_tgt->quad.xz += cx * pos.z + tc[k].quad.xz;
                tc_tgt->quad.yz += cy * pos.z + tc[k].quad.yz;
            }
#endif
        }
        //* Put the target tc (LM -> MM) (overwrite the target tc on MPE)
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tc);
        dma(dma_put, (long*)((tcLM *)adr_tc_head + head + my_offset + i), (long*)(tc_tgt));
        dma_wait(&reply_put, 1);
        while (reply_put !=1 ) {}
    }

    //* Release memory
    ldm_free(tc_tgt, bsize_tc);
    ldm_free(epj, bsize_epj_array);
    ldm_free(tc, bsize_tc_array);
    ldm_free(rsrch, bsize_rsrch);
}


void CalcMomentLongGlobalTree(void *args) {
    int my_id = athread_get_id(-1);
    //* Process the arguments
    void *adr_tc_head  = (void *)((unsigned long*)args)[0];
    void *adr_tp_head  = (void *)((unsigned long*)args)[1];
    void *adr_epj_head = (void *)((unsigned long*)args)[2];
    void *adr_spj_head = (void *)((unsigned long*)args)[3];
    int n_leaf_limit = (int)(((unsigned long*)args)[4]);
    int head = (int)(((unsigned long*)args)[5]);
    int next = (int)(((unsigned long*)args)[6]);
    void * adr_rsrch = (void *)((unsigned long*)args)[7];
    //* Compute the task of each CPE
    int num_tc,num_tc_loc,my_offset;
    num_tc = next - head;
    num_tc_loc = num_tc/NUMBER_OF_CPE + ( (my_id < num_tc % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (num_tc/NUMBER_OF_CPE)*my_id + ( (my_id < num_tc % NUMBER_OF_CPE) ? my_id : num_tc % NUMBER_OF_CPE );
#if 0
    //* Data Check (1)
    if (my_id == 63) {
        printf("--------\n");
        printf("my_id        = %d\n",my_id);
        printf("n_leaf_limit = %d\n",n_leaf_limit);
        printf("head         = %d\n",head);
        printf("next         = %d\n",next);
        printf("num_tc       = %d\n",num_tc);
        printf("num_tc_loc   = %d\n",num_tc_loc);
        printf("my_offset    = %d\n",my_offset);
    }
#endif
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 32,
    };
    size_t bsize_tc  = sizeof(tcLM);
    size_t bsize_tp  = sizeof(tpLM);
    size_t bsize_epj = sizeof(epjLM);
    size_t bsize_spj = sizeof(spjLM);
    size_t bsize_tp_array  = bsize_tp  * CHUNK_SIZE;
    size_t bsize_epj_array = bsize_epj * CHUNK_SIZE;
    size_t bsize_spj_array = bsize_spj * CHUNK_SIZE;
    size_t bsize_tc_array  = bsize_tc  * N_CHILDREN;
    size_t bsize_rsrch = sizeof(double);
    tcLM *tc_tgt = (tcLM *) ldm_malloc( bsize_tc );
    tpLM *tp = (tpLM *) ldm_malloc( bsize_tp_array );
    epjLM *epj = (epjLM *) ldm_malloc( bsize_epj_array );
    spjLM *spj = (spjLM *) ldm_malloc( bsize_spj_array );
    tcLM *tc = (tcLM *) ldm_malloc( bsize_tc_array );
    double *rsrch = (double *) ldm_malloc( bsize_rsrch );
#if 0
    //* Data Check (2)
    if (my_id == 63) {
        printf("bsize_tc        = %d\n",bsize_tc);
        printf("bsize_tp        = %d\n",bsize_tp);
        printf("bsize_tp_array  = %d\n",bsize_tp_array);
        printf("bsize_epj_array = %d\n",bsize_epj_array);
        printf("bsize_spj_array = %d\n",bsize_spj_array);
        printf("bsize_tc_array  = %d\n",bsize_tc_array);
    }
#endif
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (local vars. for moment calculation)
    F64_ mass;
    F64vec_ pos;
    double cx,cy,cz;
    //** Get the search radius
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, bsize_rsrch);
    dma(dma_get, (long*)((double *)adr_rsrch), (long*)(rsrch));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
#if 0
    if (my_id == 63) {
        printf("rsrch        = %e\n",*rsrch);
    }
#endif
    for (i=0; i<num_tc_loc; i++) {
        //* Get the target tc
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tc);
        dma(dma_get, (long*)((tcLM *)adr_tc_head + head + my_offset + i), (long*)(tc_tgt));
        dma_wait(&reply_get, 1);
        while (reply_get !=1 ) {}

        //* Moment calculation
        tc_tgt->mass = 0.0;
        tc_tgt->pos.x = 0.0;
        tc_tgt->pos.y = 0.0;
        tc_tgt->pos.z = 0.0;
#ifdef USE_QUADRUPOLE
        tc_tgt->quad.xx = 0.0;
        tc_tgt->quad.yy = 0.0;
        tc_tgt->quad.zz = 0.0;
        tc_tgt->quad.xy = 0.0;
        tc_tgt->quad.xz = 0.0;
        tc_tgt->quad.yz = 0.0;
#endif
        tc_tgt->vertex_out_.low_.x = 9999.9;
        tc_tgt->vertex_out_.low_.y = 9999.9;
        tc_tgt->vertex_out_.low_.z = 9999.9;
        tc_tgt->vertex_out_.high_.x = -9999.9;
        tc_tgt->vertex_out_.high_.y = -9999.9;
        tc_tgt->vertex_out_.high_.z = -9999.9;
        tc_tgt->vertex_in_.low_.x = 9999.9;
        tc_tgt->vertex_in_.low_.y = 9999.9;
        tc_tgt->vertex_in_.low_.z = 9999.9;
        tc_tgt->vertex_in_.high_.x = -9999.9;
        tc_tgt->vertex_in_.high_.y = -9999.9;
        tc_tgt->vertex_in_.high_.z = -9999.9;

        if (tc_tgt->n_ptcl_ == 0) {
            continue;
        } else if (tc_tgt->n_ptcl_ <= n_leaf_limit || tc_tgt->level_ == TREE_LEVEL_LIMIT){
            for (j=0; j<tc_tgt->n_ptcl_; j+= CHUNK_SIZE) {
                int nrem = tc_tgt->n_ptcl_ - j;
                int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
                //* Get tp
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_tp * nn);
                dma(dma_get, (long*)((tpLM *)adr_tp_head + tc_tgt->adr_ptcl_ + j), (long*)(tp));
                dma_wait(&reply_get, 1);
                while (reply_get !=1 ) {}
                //* Get epj
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_epj * nn);
                dma(dma_get, (long*)((epjLM *)adr_epj_head + tc_tgt->adr_ptcl_ + j), (long*)(epj));
                dma_wait(&reply_get, 1);
                while (reply_get !=1 ) {}
                //* Get spj
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_spj * nn);
                dma(dma_get, (long*)((spjLM *)adr_spj_head + tc_tgt->adr_ptcl_ + j), (long*)(spj));
                dma_wait(&reply_get, 1);
                //* Moment calculation
                for (k=0; k<nn; k++){
                    if ( ((tp[k].adr_ptcl_>>31) & 0x1) == 0) {// tp[k].adr_ptcl_ is U32_ type.
                        //** Inline expansion of tc_tmp->mom_.accumulateAtLeaf(epj[k]);
                        tc_tgt->mass += epj[k].mass;
                        tc_tgt->pos.x += epj[k].mass * epj[k].pos.x;
                        tc_tgt->pos.y += epj[k].mass * epj[k].pos.y;
                        tc_tgt->pos.z += epj[k].mass * epj[k].pos.z;
    
                        tc_tgt->vertex_out_.low_.x = ( tc_tgt->vertex_out_.low_.x <= epj[k].pos.x - *rsrch) ? tc_tgt->vertex_out_.low_.x : epj[k].pos.x - *rsrch;
                        tc_tgt->vertex_out_.low_.y = ( tc_tgt->vertex_out_.low_.y <= epj[k].pos.y - *rsrch) ? tc_tgt->vertex_out_.low_.y : epj[k].pos.y - *rsrch;
                        tc_tgt->vertex_out_.low_.z = ( tc_tgt->vertex_out_.low_.z <= epj[k].pos.z - *rsrch) ? tc_tgt->vertex_out_.low_.z : epj[k].pos.z - *rsrch;
                        tc_tgt->vertex_out_.high_.x = ( tc_tgt->vertex_out_.high_.x > epj[k].pos.x + *rsrch) ? tc_tgt->vertex_out_.high_.x : epj[k].pos.x + *rsrch;
                        tc_tgt->vertex_out_.high_.y = ( tc_tgt->vertex_out_.high_.y > epj[k].pos.y + *rsrch) ? tc_tgt->vertex_out_.high_.y : epj[k].pos.y + *rsrch;
                        tc_tgt->vertex_out_.high_.z = ( tc_tgt->vertex_out_.high_.z > epj[k].pos.z + *rsrch) ? tc_tgt->vertex_out_.high_.z : epj[k].pos.z + *rsrch;
    
                        tc_tgt->vertex_in_.low_.x = ( tc_tgt->vertex_in_.low_.x <= epj[k].pos.x ) ? tc_tgt->vertex_in_.low_.x : epj[k].pos.x;
                        tc_tgt->vertex_in_.low_.y = ( tc_tgt->vertex_in_.low_.y <= epj[k].pos.y ) ? tc_tgt->vertex_in_.low_.y : epj[k].pos.y;
                        tc_tgt->vertex_in_.low_.z = ( tc_tgt->vertex_in_.low_.z <= epj[k].pos.z ) ? tc_tgt->vertex_in_.low_.z : epj[k].pos.z;
                        tc_tgt->vertex_in_.high_.x = ( tc_tgt->vertex_in_.high_.x > epj[k].pos.x ) ? tc_tgt->vertex_in_.high_.x : epj[k].pos.x;
                        tc_tgt->vertex_in_.high_.y = ( tc_tgt->vertex_in_.high_.y > epj[k].pos.y ) ? tc_tgt->vertex_in_.high_.y : epj[k].pos.y;
                        tc_tgt->vertex_in_.high_.z = ( tc_tgt->vertex_in_.high_.z > epj[k].pos.z ) ? tc_tgt->vertex_in_.high_.z : epj[k].pos.z;
                    } else {
                        //** Inline expansion of tc_tmp->mom_.accumulate(spj[k].convertToMoment());
                        tc_tgt->mass += spj[k].mass;
                        tc_tgt->pos.x += spj[k].mass * spj[k].pos.x;
                        tc_tgt->pos.y += spj[k].mass * spj[k].pos.y;
                        tc_tgt->pos.z += spj[k].mass * spj[k].pos.z;

                        //tc_tgt->vertex_out_.low_.x = ( tc_tgt->vertex_out_.low_.x <= tc[k].vertex_out_.low_.x ) ? tc_tgt->vertex_out_.low_.x : tc[k].vertex_out_.low_.x;
                        //tc_tgt->vertex_out_.low_.y = ( tc_tgt->vertex_out_.low_.y <= tc[k].vertex_out_.low_.y ) ? tc_tgt->vertex_out_.low_.y : tc[k].vertex_out_.low_.y;
                        //tc_tgt->vertex_out_.low_.z = ( tc_tgt->vertex_out_.low_.z <= tc[k].vertex_out_.low_.z ) ? tc_tgt->vertex_out_.low_.z : tc[k].vertex_out_.low_.z;
                        //tc_tgt->vertex_out_.high_.x = ( tc_tgt->vertex_out_.high_.x > tc[k].vertex_out_.high_.x ) ? tc_tgt->vertex_out_.high_.x : tc[k].vertex_out_.high_.x;
                        //tc_tgt->vertex_out_.high_.y = ( tc_tgt->vertex_out_.high_.y > tc[k].vertex_out_.high_.y ) ? tc_tgt->vertex_out_.high_.y : tc[k].vertex_out_.high_.y;
                        //tc_tgt->vertex_out_.high_.z = ( tc_tgt->vertex_out_.high_.z > tc[k].vertex_out_.high_.z ) ? tc_tgt->vertex_out_.high_.z : tc[k].vertex_out_.high_.z;

                        //tc_tgt->vertex_in_.low_.x = ( tc_tgt->vertex_in_.low_.x <= tc[k].vertex_in_.low_.x ) ? tc_tgt->vertex_in_.low_.x : tc[k].vertex_in_.low_.x;
                        //tc_tgt->vertex_in_.low_.y = ( tc_tgt->vertex_in_.low_.y <= tc[k].vertex_in_.low_.y ) ? tc_tgt->vertex_in_.low_.y : tc[k].vertex_in_.low_.y;
                        //tc_tgt->vertex_in_.low_.z = ( tc_tgt->vertex_in_.low_.z <= tc[k].vertex_in_.low_.z ) ? tc_tgt->vertex_in_.low_.z : tc[k].vertex_in_.low_.z;
                        //tc_tgt->vertex_in_.high_.x = ( tc_tgt->vertex_in_.high_.x > tc[k].vertex_in_.high_.x ) ? tc_tgt->vertex_in_.high_.x : tc[k].vertex_in_.high_.x;
                        //tc_tgt->vertex_in_.high_.y = ( tc_tgt->vertex_in_.high_.y > tc[k].vertex_in_.high_.y ) ? tc_tgt->vertex_in_.high_.y : tc[k].vertex_in_.high_.y;
                        //tc_tgt->vertex_in_.high_.z = ( tc_tgt->vertex_in_.high_.z > tc[k].vertex_in_.high_.z ) ? tc_tgt->vertex_in_.high_.z : tc[k].vertex_in_.high_.z;
                    }
                }
            }
            //** Inline expansion of tc_tmp->mom_.set();
            tc_tgt->pos.x /= tc_tgt->mass;
            tc_tgt->pos.y /= tc_tgt->mass;
            tc_tgt->pos.z /= tc_tgt->mass;
#ifdef USE_QUADRUPOLE
            for (j=0; j<tc_tgt->n_ptcl_; j+= CHUNK_SIZE) {
                int nrem = tc_tgt->n_ptcl_ - j;
                int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
                //* Get tp
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_tp * nn);
                dma(dma_get, (long*)((tpLM *)adr_tp_head + tc_tgt->adr_ptcl_ + j), (long*)(tp));
                dma_wait(&reply_get, 1);
                while (reply_get !=1 ) {}
                //* Get epj
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_epj * nn);
                dma(dma_get, (long*)((epjLM *)adr_epj_head + tc_tgt->adr_ptcl_ + j), (long*)(epj));
                dma_wait(&reply_get, 1);
                while (reply_get !=1 ) {}
                //* Get spj
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_spj * nn);
                dma(dma_get, (long*)((spjLM *)adr_spj_head + tc_tgt->adr_ptcl_ + j), (long*)(spj));
                dma_wait(&reply_get, 1);
                //* Moment calculation
                for (k=0; k<nn; k++) {
                    if ( ((tp[k].adr_ptcl_>>31) & 0x1) == 0) {// tp[k].adr_ptcl_ is U32_ type.
                        //** Inline expansion of tc_tmp->mom_.accumulateAtLeaf2(epj[k]);
                        mass = epj[k].mass;
                        pos.x = epj[k].pos.x - tc_tgt->pos.x;
                        pos.y = epj[k].pos.y - tc_tgt->pos.y;
                        pos.z = epj[k].pos.z - tc_tgt->pos.z;
                        cx = mass * pos.x;
                        cy = mass * pos.y;
                        cz = mass * pos.z;
                        tc_tgt->quad.xx += cx * pos.x;
                        tc_tgt->quad.yy += cy * pos.y;
                        tc_tgt->quad.zz += cz * pos.z;
                        tc_tgt->quad.xy += cx * pos.y;
                        tc_tgt->quad.xz += cx * pos.z;
                        tc_tgt->quad.yz += cy * pos.z;
                    } else {
                        //** Inline expansion of tc_tmp->mom_.accumulate2(spj[k].convertToMoment());
                        mass = spj[k].mass;
                        pos.x = spj[k].pos.x - tc_tgt->pos.x;
                        pos.y = spj[k].pos.y - tc_tgt->pos.y;
                        pos.z = spj[k].pos.z - tc_tgt->pos.z;
                        cx = mass * pos.x;
                        cy = mass * pos.y;
                        cz = mass * pos.z;
                        tc_tgt->quad.xx += cx * pos.x + spj[k].quad.xx;
                        tc_tgt->quad.yy += cy * pos.y + spj[k].quad.yy;
                        tc_tgt->quad.zz += cz * pos.z + spj[k].quad.zz;
                        tc_tgt->quad.xy += cx * pos.y + spj[k].quad.xy;
                        tc_tgt->quad.xz += cx * pos.z + spj[k].quad.xz;
                        tc_tgt->quad.yz += cy * pos.z + spj[k].quad.yz;
                    }
                }
            }
#endif
        } else {
            //** Get tc (childnodes)
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tc * N_CHILDREN);
            dma(dma_get, (long*)((tcLM *)adr_tc_head + tc_tgt->adr_tc_), (long*)(tc));
            dma_wait(&reply_get, 1);
            while (reply_get !=1 ) {}
            //** Inline expansion of tc_tmp->mom_.accumulate( tc_tmp_tmp->mom_ );
            for (k=0; k<N_CHILDREN; k++){
                if (tc[k].n_ptcl_ == 0) continue;
                tc_tgt->mass += tc[k].mass;
                tc_tgt->pos.x += tc[k].mass * tc[k].pos.x;
                tc_tgt->pos.y += tc[k].mass * tc[k].pos.y;
                tc_tgt->pos.z += tc[k].mass * tc[k].pos.z;

                tc_tgt->vertex_out_.low_.x = ( tc_tgt->vertex_out_.low_.x <= tc[k].vertex_out_.low_.x ) ? tc_tgt->vertex_out_.low_.x : tc[k].vertex_out_.low_.x;
                tc_tgt->vertex_out_.low_.y = ( tc_tgt->vertex_out_.low_.y <= tc[k].vertex_out_.low_.y ) ? tc_tgt->vertex_out_.low_.y : tc[k].vertex_out_.low_.y;
                tc_tgt->vertex_out_.low_.z = ( tc_tgt->vertex_out_.low_.z <= tc[k].vertex_out_.low_.z ) ? tc_tgt->vertex_out_.low_.z : tc[k].vertex_out_.low_.z;
                tc_tgt->vertex_out_.high_.x = ( tc_tgt->vertex_out_.high_.x > tc[k].vertex_out_.high_.x ) ? tc_tgt->vertex_out_.high_.x : tc[k].vertex_out_.high_.x;
                tc_tgt->vertex_out_.high_.y = ( tc_tgt->vertex_out_.high_.y > tc[k].vertex_out_.high_.y ) ? tc_tgt->vertex_out_.high_.y : tc[k].vertex_out_.high_.y;
                tc_tgt->vertex_out_.high_.z = ( tc_tgt->vertex_out_.high_.z > tc[k].vertex_out_.high_.z ) ? tc_tgt->vertex_out_.high_.z : tc[k].vertex_out_.high_.z;

                tc_tgt->vertex_in_.low_.x = ( tc_tgt->vertex_in_.low_.x <= tc[k].vertex_in_.low_.x ) ? tc_tgt->vertex_in_.low_.x : tc[k].vertex_in_.low_.x;
                tc_tgt->vertex_in_.low_.y = ( tc_tgt->vertex_in_.low_.y <= tc[k].vertex_in_.low_.y ) ? tc_tgt->vertex_in_.low_.y : tc[k].vertex_in_.low_.y;
                tc_tgt->vertex_in_.low_.z = ( tc_tgt->vertex_in_.low_.z <= tc[k].vertex_in_.low_.z ) ? tc_tgt->vertex_in_.low_.z : tc[k].vertex_in_.low_.z;
                tc_tgt->vertex_in_.high_.x = ( tc_tgt->vertex_in_.high_.x > tc[k].vertex_in_.high_.x ) ? tc_tgt->vertex_in_.high_.x : tc[k].vertex_in_.high_.x;
                tc_tgt->vertex_in_.high_.y = ( tc_tgt->vertex_in_.high_.y > tc[k].vertex_in_.high_.y ) ? tc_tgt->vertex_in_.high_.y : tc[k].vertex_in_.high_.y;
                tc_tgt->vertex_in_.high_.z = ( tc_tgt->vertex_in_.high_.z > tc[k].vertex_in_.high_.z ) ? tc_tgt->vertex_in_.high_.z : tc[k].vertex_in_.high_.z;
            }
            //** Inline expansion of tc_tmp->mom_.set();
            tc_tgt->pos.x /= tc_tgt->mass;
            tc_tgt->pos.y /= tc_tgt->mass;
            tc_tgt->pos.z /= tc_tgt->mass;
#ifdef USE_QUADRUPOLE
            //** Inline expansion of tc_tmp->mom_.accumulate2( tc_tmp_tmp->mom_ );
            for (k=0; k<N_CHILDREN; k++){
                if (tc[k].n_ptcl_ == 0) continue;
                mass = tc[k].mass;
                pos.x = tc[k].pos.x - tc_tgt->pos.x;
                pos.y = tc[k].pos.y - tc_tgt->pos.y;
                pos.z = tc[k].pos.z - tc_tgt->pos.z;
                cx = mass * pos.x;
                cy = mass * pos.y;
                cz = mass * pos.z;
                tc_tgt->quad.xx += cx * pos.x + tc[k].quad.xx;
                tc_tgt->quad.yy += cy * pos.y + tc[k].quad.yy;
                tc_tgt->quad.zz += cz * pos.z + tc[k].quad.zz;
                tc_tgt->quad.xy += cx * pos.y + tc[k].quad.xy;
                tc_tgt->quad.xz += cx * pos.z + tc[k].quad.xz;
                tc_tgt->quad.yz += cy * pos.z + tc[k].quad.yz;
            }
#endif
        }
        //* Put the target tc (LM -> MM) (overwrite the target tc on MPE)
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tc);
        dma(dma_put, (long*)((tcLM *)adr_tc_head + head + my_offset + i), (long*)(tc_tgt));
        dma_wait(&reply_put, 1);
        while (reply_put !=1 ) {}
    }

    //* Release memory
    ldm_free(tc_tgt, bsize_tc);
    ldm_free(tp, bsize_tp_array);
    ldm_free(epj, bsize_epj_array);
    ldm_free(spj, bsize_spj_array);
    ldm_free(tc, bsize_tc_array);
    ldm_free(rsrch, bsize_rsrch);
}
