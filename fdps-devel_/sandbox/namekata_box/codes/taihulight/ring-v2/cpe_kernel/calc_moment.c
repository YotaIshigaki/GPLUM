/* Standard C headers */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
/* Headers specific to Sunway TaihuLight */
#include "slave.h"
#include "dma.h"
#include "cpe_func.h"

extern double __ieee754_atan2 __P((double, double, double *, double *, double * ));

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

#ifndef REMOVE_VERTEX
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
#endif
        
#ifdef PHI_R_TREE
        tc_tgt->pos_phi = 0.0;
        tc_tgt->pos_r   = 0.0;
#endif
        
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

#ifdef PHI_R_TREE
                   double pos_phi = __ieee754_atan2(epj[k].pos.y, epj[k].pos.x, atanhi, atanlo, aT);
                   if(pos_phi < 0.0) pos_phi += MY_PI_2;
                   double pos_r = sqrt(epj[k].pos.x*epj[k].pos.x + epj[k].pos.y*epj[k].pos.y);

                   tc_tgt->pos_phi += epj[k].mass * pos_phi;
                   tc_tgt->pos_r   += epj[k].mass * pos_r;                   

    #ifndef REMOVE_VERTEX
                   tc_tgt->vertex_out_.low_.x = ( tc_tgt->vertex_out_.low_.x <= pos_phi - *rsrch) ? tc_tgt->vertex_out_.low_.x : pos_phi - *rsrch;
                   tc_tgt->vertex_out_.low_.y = ( tc_tgt->vertex_out_.low_.y <= pos_r - *rsrch) ? tc_tgt->vertex_out_.low_.y : pos_r - *rsrch;
                   tc_tgt->vertex_out_.low_.z = ( tc_tgt->vertex_out_.low_.z <= epj[k].pos.z - *rsrch) ? tc_tgt->vertex_out_.low_.z : epj[k].pos.z - *rsrch;
                   tc_tgt->vertex_out_.high_.x = ( tc_tgt->vertex_out_.high_.x > pos_phi + *rsrch) ? tc_tgt->vertex_out_.high_.x : pos_phi + *rsrch;
                   tc_tgt->vertex_out_.high_.y = ( tc_tgt->vertex_out_.high_.y > pos_r + *rsrch) ? tc_tgt->vertex_out_.high_.y : pos_r + *rsrch;
                   tc_tgt->vertex_out_.high_.z = ( tc_tgt->vertex_out_.high_.z > epj[k].pos.z + *rsrch) ? tc_tgt->vertex_out_.high_.z : epj[k].pos.z + *rsrch;

                   tc_tgt->vertex_in_.low_.x = ( tc_tgt->vertex_in_.low_.x <= pos_phi - *rsrch) ? tc_tgt->vertex_in_.low_.x : pos_phi - *rsrch;
                   tc_tgt->vertex_in_.low_.y = ( tc_tgt->vertex_in_.low_.y <= pos_r - *rsrch) ? tc_tgt->vertex_in_.low_.y : pos_r - *rsrch;
                   tc_tgt->vertex_in_.low_.z = ( tc_tgt->vertex_in_.low_.z <= epj[k].pos.z - *rsrch) ? tc_tgt->vertex_in_.low_.z : epj[k].pos.z - *rsrch;
                   tc_tgt->vertex_in_.high_.x = ( tc_tgt->vertex_in_.high_.x > pos_phi + *rsrch) ? tc_tgt->vertex_in_.high_.x : pos_phi + *rsrch;
                   tc_tgt->vertex_in_.high_.y = ( tc_tgt->vertex_in_.high_.y > pos_r + *rsrch) ? tc_tgt->vertex_in_.high_.y : pos_r + *rsrch;
                   tc_tgt->vertex_in_.high_.z = ( tc_tgt->vertex_in_.high_.z > epj[k].pos.z + *rsrch) ? tc_tgt->vertex_in_.high_.z : epj[k].pos.z + *rsrch;
    #endif
                   
#else

    #ifndef REMOVE_VERTEX
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
    #endif
#endif
                }
            }
            //** Inline expansion of tc_tmp->mom_.set();
            tc_tgt->pos.x /= tc_tgt->mass;
            tc_tgt->pos.y /= tc_tgt->mass;
            tc_tgt->pos.z /= tc_tgt->mass;
#ifdef PHI_R_TREE
            tc_tgt->pos_phi /= tc_tgt->mass;
            tc_tgt->pos_r   /= tc_tgt->mass;
#endif

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
    #ifndef REMOVE_VERTEX
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
    #endif
#ifdef PHI_R_TREE
                tc_tgt->pos_phi += tc[k].mass * tc[k].pos_phi;
                tc_tgt->pos_r   += tc[k].mass * tc[k].pos_r;
#endif
            }
            //** Inline expansion of tc_tmp->mom_.set();
            tc_tgt->pos.x /= tc_tgt->mass;
            tc_tgt->pos.y /= tc_tgt->mass;
            tc_tgt->pos.z /= tc_tgt->mass;
#ifdef PHI_R_TREE
            tc_tgt->pos_phi /= tc_tgt->mass;
            tc_tgt->pos_r /= tc_tgt->mass;
#endif

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

    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);    
}

void CalcMomentLongGlobalTree(void *args) {
    int my_id = athread_get_id(-1);
    double MY_PI = 3.14159265358979323846;
    double MY_PI_2 = 2.0*MY_PI;
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
        
#ifdef PHI_R_TREE
        tc_tgt->pos_phi = 0.0;
        tc_tgt->pos_r   = 0.0;
#endif
        
#ifdef USE_QUADRUPOLE
        tc_tgt->quad.xx = 0.0;
        tc_tgt->quad.yy = 0.0;
        tc_tgt->quad.zz = 0.0;
        tc_tgt->quad.xy = 0.0;
        tc_tgt->quad.xz = 0.0;
        tc_tgt->quad.yz = 0.0;
#endif
    #ifndef REMOVE_VERTEX
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
    #endif
        
        if (tc_tgt->n_ptcl_ == 0) {
            continue;
        }
        else if (tc_tgt->n_ptcl_ <= n_leaf_limit || tc_tgt->level_ == TREE_LEVEL_LIMIT){
            // leaf
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
#ifdef REDUCE_MEMORY
                U32_ offset_epj = 0;
                U32_ offset_spj = 0;
                int has_epj = 0;
                int has_spj = 0;
                for (k=0; k<nn; k++){
                    if ( ((tp[k].adr_ptcl_>>31) & 0x1) == 0 && has_epj == 0) {
                        offset_epj = tp[k].adr_ptcl_;
                        has_epj = 1;
                    }
                    else if ( ((tp[k].adr_ptcl_>>31) & 0x1) == 1 && has_spj == 0) {
                        offset_spj = (tp[k].adr_ptcl_ & 0x7fffffff);
                        has_spj = 1;
                    }
                }
#endif
                //* Get epj
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_epj * nn);
#ifdef REDUCE_MEMORY
                dma(dma_get, (long*)((epjLM *)adr_epj_head + offset_epj + j), (long*)(epj));
#else
                dma(dma_get, (long*)((epjLM *)adr_epj_head + tc_tgt->adr_ptcl_ + j), (long*)(epj));
#endif
                dma_wait(&reply_get, 1);
                while (reply_get !=1 ) {}
                //* Get spj
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_spj * nn);
#ifdef REDUCE_MEMORY
                dma(dma_get, (long*)((spjLM *)adr_spj_head + offset_spj + j), (long*)(spj));
#else
                dma(dma_get, (long*)((spjLM *)adr_spj_head + tc_tgt->adr_ptcl_ + j), (long*)(spj));
#endif
                dma_wait(&reply_get, 1);

#ifdef REDUCE_MEMORY
                int n_cnt_ep = 0;
                int n_cnt_sp = 0;
                int l;
                for (l=0; l<nn; l++){
                    if ( ((tp[l].adr_ptcl_>>31) & 0x1) == 0) {// tp[k].adr_ptcl_ is U32_ type.
                        k = n_cnt_ep;
                        n_cnt_ep++;
#else
                for (k=0; k<nn; k++){
                    if ( ((tp[k].adr_ptcl_>>31) & 0x1) == 0) {// tp[k].adr_ptcl_ is U32_ type.
#endif
                        //** Inline expansion of tc_tmp->mom_.accumulateAtLeaf(epj[k]);
                        tc_tgt->mass += epj[k].mass;
                        tc_tgt->pos.x += epj[k].mass * epj[k].pos.x;
                        tc_tgt->pos.y += epj[k].mass * epj[k].pos.y;
                        tc_tgt->pos.z += epj[k].mass * epj[k].pos.z;

#ifdef PHI_R_TREE
                        double pos_phi = __ieee754_atan2(epj[k].pos.y, epj[k].pos.x, atanhi, atanlo, aT);
                        if(pos_phi < 0.0) pos_phi += MY_PI_2;
                        double pos_r = sqrt(epj[k].pos.x*epj[k].pos.x + epj[k].pos.y*epj[k].pos.y);

                        tc_tgt->pos_phi += epj[k].mass * pos_phi;
                        tc_tgt->pos_r   += epj[k].mass * pos_r;

    #ifndef REMOVE_VERTEX                        
                        tc_tgt->vertex_out_.low_.x = ( tc_tgt->vertex_out_.low_.x <= pos_phi - *rsrch) ? tc_tgt->vertex_out_.low_.x : pos_phi - *rsrch;
                        tc_tgt->vertex_out_.low_.y = ( tc_tgt->vertex_out_.low_.y <= pos_r - *rsrch) ? tc_tgt->vertex_out_.low_.y : pos_r - *rsrch;
                        tc_tgt->vertex_out_.low_.z = ( tc_tgt->vertex_out_.low_.z <= epj[k].pos.z - *rsrch) ? tc_tgt->vertex_out_.low_.z : epj[k].pos.z - *rsrch;
                        tc_tgt->vertex_out_.high_.x = ( tc_tgt->vertex_out_.high_.x > pos_phi + *rsrch) ? tc_tgt->vertex_out_.high_.x : pos_phi + *rsrch;
                        tc_tgt->vertex_out_.high_.y = ( tc_tgt->vertex_out_.high_.y > pos_r + *rsrch) ? tc_tgt->vertex_out_.high_.y : pos_r + *rsrch;
                        tc_tgt->vertex_out_.high_.z = ( tc_tgt->vertex_out_.high_.z > epj[k].pos.z + *rsrch) ? tc_tgt->vertex_out_.high_.z : epj[k].pos.z + *rsrch;
    
                        tc_tgt->vertex_in_.low_.x = ( tc_tgt->vertex_in_.low_.x <= pos_phi ) ? tc_tgt->vertex_in_.low_.x : pos_phi;
                        tc_tgt->vertex_in_.low_.y = ( tc_tgt->vertex_in_.low_.y <= pos_r ) ? tc_tgt->vertex_in_.low_.y : pos_r;
                        tc_tgt->vertex_in_.low_.z = ( tc_tgt->vertex_in_.low_.z <= epj[k].pos.z ) ? tc_tgt->vertex_in_.low_.z : epj[k].pos.z;
                        tc_tgt->vertex_in_.high_.x = ( tc_tgt->vertex_in_.high_.x > pos_phi ) ? tc_tgt->vertex_in_.high_.x : pos_phi;
                        tc_tgt->vertex_in_.high_.y = ( tc_tgt->vertex_in_.high_.y > pos_r ) ? tc_tgt->vertex_in_.high_.y : pos_r;
                        tc_tgt->vertex_in_.high_.z = ( tc_tgt->vertex_in_.high_.z > epj[k].pos.z ) ? tc_tgt->vertex_in_.high_.z : epj[k].pos.z;
    #endif
                        

#else
    #ifndef REMOVE_VERTEX
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
    #endif
#endif
                    }
                    else {
#ifdef REDUCE_MEMORY
                        k = n_cnt_sp;
                        n_cnt_sp++;
#endif
                        //** Inline expansion of tc_tmp->mom_.accumulate(spj[k].convertToMoment());
                        tc_tgt->mass += spj[k].mass;
                        tc_tgt->pos.x += spj[k].mass * spj[k].pos.x;
                        tc_tgt->pos.y += spj[k].mass * spj[k].pos.y;
                        tc_tgt->pos.z += spj[k].mass * spj[k].pos.z;
#ifdef PHI_R_TREE
                        tc_tgt->pos_phi += spj[k].mass * spj[k].pos_phi;
                        tc_tgt->pos_r   += spj[k].mass * spj[k].pos_r;
#endif
                    }
                }
            }
            //** Inline expansion of tc_tmp->mom_.set();
            tc_tgt->pos.x /= tc_tgt->mass;
            tc_tgt->pos.y /= tc_tgt->mass;
            tc_tgt->pos.z /= tc_tgt->mass;
#ifdef PHI_R_TREE
            tc_tgt->pos_phi /= tc_tgt->mass;
            tc_tgt->pos_r /= tc_tgt->mass;
#endif
            
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
        }
        else {
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
                        
    #ifndef REMOVE_VERTEX
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
    #endif
                
#ifdef PHI_R_TREE
                tc_tgt->pos_phi += tc[k].mass * tc[k].pos_phi;
                tc_tgt->pos_r   += tc[k].mass * tc[k].pos_r;
#endif
            }
            //** Inline expansion of tc_tmp->mom_.set();
            tc_tgt->pos.x /= tc_tgt->mass;
            tc_tgt->pos.y /= tc_tgt->mass;
            tc_tgt->pos.z /= tc_tgt->mass;
            
#ifdef PHI_R_TREE
            tc_tgt->pos_phi /= tc_tgt->mass;
            tc_tgt->pos_r /= tc_tgt->mass;
#endif
            
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
    
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
    
}

