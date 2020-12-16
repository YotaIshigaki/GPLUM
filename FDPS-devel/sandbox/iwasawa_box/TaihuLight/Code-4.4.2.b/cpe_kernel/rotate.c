/* Standard C headers */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
/* Headers specific to Sunway TaihuLight */
#include "slave.h"
#include "dma.h"
#include "cpe_func.h"

void Rotate(void *args) {
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int n_tot = (int)(((unsigned long*)args)[0]); // total number of local particles
    void *adr_fp_head   = (void *)((unsigned long*)args)[1];
    void *adr_cos_theta = (void *)((unsigned long*)args)[2];
    void *adr_sin_theta = (void *)((unsigned long*)args)[3];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_tot/NUMBER_OF_CPE + ( (my_id < n_tot % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_tot/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot % NUMBER_OF_CPE) ? my_id : n_tot % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_fp = sizeof(fpLM);
    size_t bsize_fp_array = bsize_fp * CHUNK_SIZE;
    size_t bsize_cth = sizeof(double);
    size_t bsize_sth = sizeof(double);
    fpLM *fp = (fpLM *) ldm_malloc( bsize_fp_array );
    double *cth = (double *) ldm_malloc( bsize_cth ); // cos(\theta)
    double *sth = (double *) ldm_malloc( bsize_sth ); // sin(\theta)
    //* Compute on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (local vars. for rotation of coordinates)
    F64_ x_new,y_new,vx_new,vy_new;
    //** Get cos(\theta) 
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, bsize_cth);
    dma(dma_get, (long*)((double *)adr_cos_theta), (long*)(cth));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
    //** Get sin(\theta)
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, bsize_sth);
    dma(dma_get, (long*)((double *)adr_sin_theta), (long*)(sth));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
    //** Rotate 
    for (i=0; i<n_loc; i+= CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get fp  
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_fp * nn);
        dma(dma_get, (long)((fpLM *)adr_fp_head + my_offset + i), (long)(fp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        //* Compute
        for(k=0; k<nn; k++){
            x_new  = (*cth) * fp[k].pos.x - (*sth) * fp[k].pos.y;
            y_new  = (*sth) * fp[k].pos.x + (*cth) * fp[k].pos.y;
            vx_new = (*cth) * fp[k].vel.x - (*sth) * fp[k].vel.y;
            vy_new = (*sth) * fp[k].vel.x + (*cth) * fp[k].vel.y;
            fp[k].pos.x = x_new;
            fp[k].pos.y = y_new;
            fp[k].vel.x = vx_new;
            fp[k].vel.y = vy_new;
        }
        //** Put fp
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_fp * nn);
        dma(dma_put, (long)((fpLM *)adr_fp_head + my_offset + i), (long)(fp));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }

    //* Release memory
    ldm_free(fp, bsize_fp_array);
    ldm_free(cth, bsize_cth);
    ldm_free(sth, bsize_sth);

}
