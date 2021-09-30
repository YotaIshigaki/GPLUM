#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
//#include<athread.h>
#include"slave.h"
#include"dma.h"
#include"cpe_func.h"

extern double __ieee754_atan2 __P((double, double, double *, double *, double * ));

//#include"fdlibm.h"

void CopyTCToETC(void *args){
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int n_tc_      = (int)(((unsigned long*)args)[0]);
    void *adr_tc_  = (void *)((unsigned long*)args)[1];
    void *adr_etc_ = (void *)((unsigned long*)args)[2];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_tc_/NUMBER_OF_CPE + ( (my_id < n_tc_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_tc_/NUMBER_OF_CPE)*my_id + ( (my_id < n_tc_ % NUMBER_OF_CPE) ? my_id : n_tc_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tc        = sizeof(tcLM);
    size_t bsize_etc       = sizeof(etcLM);
    size_t bsize_tc_array  = bsize_tc  * CHUNK_SIZE;
    size_t bsize_etc_array = bsize_etc * CHUNK_SIZE;
    tcLM  *tc  = (tcLM *)  ldm_malloc( bsize_tc_array );
    etcLM *etc = (etcLM *) ldm_malloc( bsize_etc_array );
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get tc
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tc*nn);
        dma(dma_get, (long*)((tcLM *)adr_tc_ + my_offset + i), (long*)(tc));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            etc[j].adr_tc_   = tc[j].adr_tc_;
            etc[j].adr_ptcl_ = tc[j].adr_ptcl_;
            etc[j].n_ptcl_   = tc[j].n_ptcl_;
            etc[j].level_    = tc[j].level_;            
        }
        //** Put etc
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_etc * nn);
        dma(dma_put, (long*)((etcLM *)adr_etc_ + my_offset + i), (long*)(etc));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }

    //* Release memory
    ldm_free(tc,  bsize_tc_array);
    ldm_free(etc, bsize_etc_array);
}


void CopyEPISortedToEPJSortedLoc(void *args){
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int n_tot_             = (int)(((unsigned long*)args)[0]);
    void * epi_sorted_     = (void *)((unsigned long*)args)[1];
    void * epj_sorted_loc_ = (void *)((unsigned long*)args)[2];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_tot_/NUMBER_OF_CPE + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_tot_/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? my_id : n_tot_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_epi       = sizeof(epiLM);
    size_t bsize_epj       = sizeof(epjLM);
    size_t bsize_epi_array = bsize_epj * CHUNK_SIZE;
    size_t bsize_epj_array = bsize_epj * CHUNK_SIZE;
    epiLM  *epi     = (epiLM *) ldm_malloc( bsize_epi_array );
    epjLM  *epj     = (epjLM *) ldm_malloc( bsize_epj_array );
    epjLM  *epj_org = (epjLM *) ldm_malloc( bsize_epj_array );
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get epi_sorted_loc_ 
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_epi*nn);
        dma(dma_get, (long*)(((epiLM *)epi_sorted_)+my_offset+i), (long*)(epi));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_epj*nn);
        dma(dma_get, (long*)(((epjLM *)epj_sorted_loc_)+my_offset+i), (long*)(epj_org));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}        

        
        for(j=0; j<nn; j++){
            epj[j].id      = epj_org[j].id;
            epj[j].mass    = epj_org[j].mass;
            //epj[j].pos_phi = epj_org[j].pos_phi;
            //epj[j].pos_r   = epj_org[j].pos_r;
            epj[j].pos.x = epi[j].pos.x;
            epj[j].pos.y = epi[j].pos.y;
            epj[j].pos.z = epi[j].pos.z;
            epj[j].vel.x = epi[j].vel.x;
            epj[j].vel.y = epi[j].vel.y;
            epj[j].vel.z = epi[j].vel.z;
        }
        
        //** Put into epj_sorted_ 
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_epj*nn);
        dma(dma_put, (long*)(((epjLM *)epj_sorted_loc_)+my_offset+i), (long*)(epj));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    //* Release memory
    ldm_free(epi, bsize_epi_array);
    ldm_free(epj, bsize_epj_array);
    ldm_free(epj_org, bsize_epj_array);
}

#if 0
void CopyEPJLocToEPJ(void *args){ 
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int n_epj_sorted_             = (int)(((unsigned long*)args)[0]);
    void *adr_epj_sorted_loc_     = (void *)((unsigned long*)args)[1];
    void *adr_epj_sorted_         = (void *)((unsigned long*)args)[2];
    void *adr_adr_epj_sorted_loc_ = (void *)((unsigned long*)args)[3];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_epj_sorted_/NUMBER_OF_CPE + ( (my_id < n_epj_sorted_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_epj_sorted_/NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted_ % NUMBER_OF_CPE) ? my_id : n_epj_sorted_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_epj       = sizeof(epjLM);
    size_t bsize_epj_array = bsize_epj * CHUNK_SIZE;
    size_t bsize_adr_array = sizeof(int) * CHUNK_SIZE;
    epjLM  *epj = (epjLM *) ldm_malloc( bsize_epj_array );
    int    *adr = (int *)   ldm_malloc( bsize_adr_array );
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get adr_tmp on MPE
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, sizeof(int) * nn);
        dma(dma_get, (long*)((int *)adr_adr_epj_sorted_loc_ + my_offset + i), (long*)(adr));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        //** Get epj_sorted_loc_ 
        int n_get_issued=0;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_epj);
        for (j=0; j<nn; j++) {
            if (adr[j] != -1) {
                n_get_issued++;
                dma(dma_get, (long*)((epjLM *)adr_epj_sorted_loc_ + adr[j]), (long*)(&epj[j]));
            }
        }
        dma_wait(&reply_get, n_get_issued);
        while (reply_get != n_get_issued) {}
        //** Put into epj_sorted_ 
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_epj * nn);
        dma(dma_put, (long*)((epjLM *)adr_epj_sorted_ + my_offset + i), (long*)(epj));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }

    //* Release memory
    ldm_free(epj, bsize_epj_array);
    ldm_free(adr, bsize_adr_array);

}
#endif


#if 1
void CopyIndirect(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[ELEMENT_SIZE_LIMIT];
    int adr_lm[LOCAL_MEMORY_SIZE_LIMIT/4];    
    int n_element = (int)(((unsigned long*)arg)[0]);
    double * data_src = (double*)((unsigned long*)arg)[1];
    double * data_dst = (double*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int element_size = element_size_in_byte / sizeof(double);
    //assert(element_size*4 <= ELEMENT_SIZE_LIMIT); // factor 4 means 4 unroll
    //assert(element_size*8 <= ELEMENT_SIZE_LIMIT); // factor 8 means 8 unroll
    assert(element_size*16 <= ELEMENT_SIZE_LIMIT);
    int * adr_array = (int*)(((unsigned long*)arg)[4]);
    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );    
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int data_block_limit = data_size_limit_in_byte / sizeof(int);
    int data_size_rest = data_size_total;
    int data_loc = id_head;
    
    while(data_size_rest > 0 ){
        int data_size_tmp = (data_block_limit < data_size_rest) ? data_block_limit : data_size_rest;
        reply_get = reply_put = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, data_size_tmp*sizeof(int));
        dma(dma_get, (long)((int*)adr_array+data_loc), (long)(adr_lm));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        /*
        for(i=0; i+4<=data_size_tmp; i+=4){
            reply_get = reply_put = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, element_size_in_byte*4);
            dma(dma_get, (long)((double*)data_src+(data_loc*element_size)), (long)(data_lm));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            data_loc+=4;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, element_size_in_byte);
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+0]*element_size)), (long)((double*)data_lm+element_size*0));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+1]*element_size)), (long)((double*)data_lm+element_size*1));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+2]*element_size)), (long)((double*)data_lm+element_size*2));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+3]*element_size)), (long)((double*)data_lm+element_size*3));
            dma_wait(&reply_put, 4);
            while (reply_put != 4) {}
        }
        */
        /*
        for(i=0; i+8<=data_size_tmp; i+=8){
            reply_get = reply_put = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, element_size_in_byte*8);
            dma(dma_get, (long)((double*)data_src+(data_loc*element_size)), (long)(data_lm));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            data_loc+=8;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, element_size_in_byte);
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+0]*element_size)), (long)((double*)data_lm+element_size*0));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+1]*element_size)), (long)((double*)data_lm+element_size*1));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+2]*element_size)), (long)((double*)data_lm+element_size*2));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+3]*element_size)), (long)((double*)data_lm+element_size*3));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+4]*element_size)), (long)((double*)data_lm+element_size*4));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+5]*element_size)), (long)((double*)data_lm+element_size*5));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+6]*element_size)), (long)((double*)data_lm+element_size*6));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+7]*element_size)), (long)((double*)data_lm+element_size*7));            
            dma_wait(&reply_put, 8);
            while (reply_put != 8) {}
        } 
        */
        for(i=0; i+16<=data_size_tmp; i+=16){
            reply_get = reply_put = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, element_size_in_byte*16);
            dma(dma_get, (long)((double*)data_src+(data_loc*element_size)), (long)(data_lm));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            data_loc+=16;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, element_size_in_byte);
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 0]*element_size)), (long)((double*)data_lm+element_size* 0));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 1]*element_size)), (long)((double*)data_lm+element_size* 1));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 2]*element_size)), (long)((double*)data_lm+element_size* 2));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 3]*element_size)), (long)((double*)data_lm+element_size* 3));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 4]*element_size)), (long)((double*)data_lm+element_size* 4));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 5]*element_size)), (long)((double*)data_lm+element_size* 5));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 6]*element_size)), (long)((double*)data_lm+element_size* 6));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 7]*element_size)), (long)((double*)data_lm+element_size* 7));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 8]*element_size)), (long)((double*)data_lm+element_size* 8));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+ 9]*element_size)), (long)((double*)data_lm+element_size* 9));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+10]*element_size)), (long)((double*)data_lm+element_size*10));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+11]*element_size)), (long)((double*)data_lm+element_size*11));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+12]*element_size)), (long)((double*)data_lm+element_size*12));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+13]*element_size)), (long)((double*)data_lm+element_size*13));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+14]*element_size)), (long)((double*)data_lm+element_size*14));
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+15]*element_size)), (long)((double*)data_lm+element_size*15));
            dma_wait(&reply_put, 16);
            while (reply_put != 16) {}
        } 
        for(; i<data_size_tmp; i++){
            reply_get = reply_put = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, element_size_in_byte);
            dma(dma_get, (long)((double*)data_src+(data_loc*element_size)), (long)(data_lm));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            data_loc++;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, element_size_in_byte);
            dma(dma_put, (long)((double*)data_dst+(adr_lm[i+0]*element_size)), (long)((double*)data_lm+element_size*0));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }
        data_size_rest -= data_size_tmp;
    }
}
#else
// original
void CopyIndirect(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[LOCAL_MEMORY_SIZE_LIMIT/8];
    int n_element = (int)(((unsigned long*)arg)[0]);
    double * data_src = (double*)((unsigned long*)arg)[1];
    double * data_dst = (double*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int element_size = element_size_in_byte / sizeof(double);
    int * adr_array = (int*)(((unsigned long*)arg)[4]);
    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );    
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int n_loop = 0;
    for(i=id_head; i<id_head+data_size_total; i++, n_loop++){
        reply_get = reply_put = 0;
        int adr = (int)(((int*)adr_array))[i];
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, element_size_in_byte);

        dma(dma_get, (long)(((double*)data_src)+(long)((long)i*(long)element_size)), (long)(data_lm));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, element_size_in_byte);

        dma(dma_put, (long)(((double*)data_dst)+((long)adr*(long)element_size)), (long)(data_lm));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
}
#endif

#if 1
void CopyIndirectInverse(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[ELEMENT_SIZE_LIMIT];
    int adr_lm[LOCAL_MEMORY_SIZE_LIMIT/4];
    int n_element = (int)(((unsigned long*)arg)[0]);
    double * data_src = (double*)((unsigned long*)arg)[1];
    double * data_dst = (double*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int element_size = element_size_in_byte / sizeof(double);
    //assert(element_size*4 <= ELEMENT_SIZE_LIMIT); // factor 4 means 4 unroll
    //assert(element_size*8 <= ELEMENT_SIZE_LIMIT); // factor 8 means 8 unroll
    assert(element_size*16 <= ELEMENT_SIZE_LIMIT);
    unsigned int * adr_array = (unsigned int*)(((unsigned long*)arg)[4]);
    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int data_block_limit = data_size_limit_in_byte / sizeof(int);
    int data_size_rest = data_size_total;
    int data_loc = id_head;
    
    while(data_size_rest > 0 ){
        int data_size_tmp = (data_block_limit < data_size_rest) ? data_block_limit : data_size_rest;
        reply_get = reply_put = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, data_size_tmp*sizeof(int));
        dma(dma_get, (long)((int*)adr_array+data_loc), (long)(adr_lm));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        for(i=0; i+16<=data_size_tmp; i+=16, data_loc+=16){
            reply_get = reply_put = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, element_size_in_byte);
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 0]*element_size)), (long)((double*)data_lm+element_size* 0));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 1]*element_size)), (long)((double*)data_lm+element_size* 1));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 2]*element_size)), (long)((double*)data_lm+element_size* 2));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 3]*element_size)), (long)((double*)data_lm+element_size* 3));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 4]*element_size)), (long)((double*)data_lm+element_size* 4));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 5]*element_size)), (long)((double*)data_lm+element_size* 5));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 6]*element_size)), (long)((double*)data_lm+element_size* 6));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 7]*element_size)), (long)((double*)data_lm+element_size* 7));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 8]*element_size)), (long)((double*)data_lm+element_size* 8));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 9]*element_size)), (long)((double*)data_lm+element_size* 9));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+10]*element_size)), (long)((double*)data_lm+element_size*10));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+11]*element_size)), (long)((double*)data_lm+element_size*11));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+12]*element_size)), (long)((double*)data_lm+element_size*12));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+13]*element_size)), (long)((double*)data_lm+element_size*13));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+14]*element_size)), (long)((double*)data_lm+element_size*14));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+15]*element_size)), (long)((double*)data_lm+element_size*15));            
            
            dma_wait(&reply_get, 16);
            while (reply_get != 16) {}
            
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, element_size_in_byte*16);
            dma(dma_put, (long)((double*)data_dst+(data_loc*element_size)), (long)(data_lm));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        } 
        for(; i<data_size_tmp; i++, data_loc++){
            reply_get = reply_put = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, element_size_in_byte);
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+0]*element_size)), (long)((double*)data_lm+element_size*0));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, element_size_in_byte);
            dma(dma_put, (long)((double*)data_dst+(data_loc*element_size)), (long)(data_lm));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }
        data_size_rest -= data_size_tmp;
    }
}
#else
// original
void CopyIndirectInverse(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[LOCAL_MEMORY_SIZE_LIMIT/8];
    int n_element = (int)(((unsigned long*)arg)[0]);
    double * data_src = (double*)((unsigned long*)arg)[1];
    double * data_dst = (double*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int element_size = element_size_in_byte / sizeof(double);
    int * adr_array = (int*)(((unsigned long*)arg)[4]);
    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );    
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    for(i=id_head; i<id_head+data_size_total; i++){
        reply_get = reply_put = 0;
        int adr = (int)(((int*)adr_array))[i];
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, element_size_in_byte);
        dma(dma_get, (long)(((double*)data_src)+((long)adr*(long)element_size)), (long)(data_lm));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, element_size_in_byte);
        dma(dma_put, (long)(((double*)data_dst)+(long)((long)i*(long)element_size)), (long)(data_lm));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
}
#endif

#if 0
void CopyIndirectInverse2(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[ELEMENT_SIZE_LIMIT];
    int adr_lm[LOCAL_MEMORY_SIZE_LIMIT/4];
    int n_element = (int)(((unsigned long*)arg)[0]);
    double * data_src = (double*)((unsigned long*)arg)[1];
    double * data_dst = (double*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int element_size = element_size_in_byte / sizeof(double);
    //assert(element_size*4 <= ELEMENT_SIZE_LIMIT); // factor 4 means 4 unroll
    //assert(element_size*8 <= ELEMENT_SIZE_LIMIT); // factor 8 means 8 unroll
    assert(element_size*16 <= ELEMENT_SIZE_LIMIT);
    int * adr_array = (int*)(((unsigned long*)arg)[4]);
    int adr_type_in_byte = (int)(((unsigned long*)arg)[5]);
    int adr_type = adr_type_in_byte / sizeof(int);
    int adr_stride_in_byte = (int)(((unsigned long*)arg)[6]);
    int adr_stride = adr_stride_in_byte / sizeof(int);

    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );    
    dma_desc dma_get, dma_put, dma_get2;
    volatile int reply_get  = 0;
    volatile int reply_get2 = 0;
    volatile int reply_put  = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_get2, &reply_get2);
    dma_descriptor_init(&dma_put, &reply_put);
    int data_block_limit = data_size_limit_in_byte / (adr_type_in_byte+adr_stride_in_byte);
    //int data_block_limit = 10;
    int data_size_rest = data_size_total;
    int data_loc = id_head;
    
    while(data_size_rest > 0 ){
        int data_size_tmp = (data_block_limit < data_size_rest) ? data_block_limit : data_size_rest;
        reply_get2 = 0;
        dma_descriptor_init(&dma_get2, &reply_get2);
        dma_set_op(&dma_get2, DMA_GET);
        dma_set_mode(&dma_get2, PE_MODE);
        dma_set_reply(&dma_get2, &reply_get2);
        dma_set_size(&dma_get2, data_size_tmp*(adr_type_in_byte+adr_stride_in_byte));
        dma_set_bsize(&dma_get2, adr_type_in_byte);
        dma_set_stepsize(&dma_get2, adr_stride_in_byte);
        dma(dma_get2, (long)((int*)adr_array+(data_loc*(adr_type+adr_stride))), (long)(adr_lm));
        dma_wait(&reply_get2, 1);
        while (reply_get2 != 1) {}
        for(i=0; i+16<=data_size_tmp; i+=16, data_loc+=16){
            reply_get = reply_put = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, element_size_in_byte);
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 0]*element_size)), (long)((double*)data_lm+element_size* 0));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 1]*element_size)), (long)((double*)data_lm+element_size* 1));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 2]*element_size)), (long)((double*)data_lm+element_size* 2));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 3]*element_size)), (long)((double*)data_lm+element_size* 3));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 4]*element_size)), (long)((double*)data_lm+element_size* 4));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 5]*element_size)), (long)((double*)data_lm+element_size* 5));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 6]*element_size)), (long)((double*)data_lm+element_size* 6));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 7]*element_size)), (long)((double*)data_lm+element_size* 7));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 8]*element_size)), (long)((double*)data_lm+element_size* 8));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+ 9]*element_size)), (long)((double*)data_lm+element_size* 9));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+10]*element_size)), (long)((double*)data_lm+element_size*10));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+11]*element_size)), (long)((double*)data_lm+element_size*11));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+12]*element_size)), (long)((double*)data_lm+element_size*12));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+13]*element_size)), (long)((double*)data_lm+element_size*13));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+14]*element_size)), (long)((double*)data_lm+element_size*14));
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+15]*element_size)), (long)((double*)data_lm+element_size*15));            
            dma_wait(&reply_get, 16);
            while (reply_get != 16) {}
            
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, element_size_in_byte*16);
            dma(dma_put, (long)((double*)data_dst+(data_loc*element_size)), (long)(data_lm));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        } 
        for(; i<data_size_tmp; i++, data_loc++){
            reply_get = reply_put = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, element_size_in_byte);
            dma(dma_get, (long)((double*)data_src+(adr_lm[i+0]*element_size)), (long)((double*)data_lm+element_size*0));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, element_size_in_byte);
            dma(dma_put, (long)((double*)data_dst+(data_loc*element_size)), (long)(data_lm));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }
        data_size_rest -= data_size_tmp;
    }
}
#elif 0
void CopyIndirectInverse2(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[LOCAL_MEMORY_SIZE_LIMIT/8];
    int n_element = (int)(((unsigned long*)arg)[0]);
    double * data_src = (double*)((unsigned long*)arg)[1];
    double * data_dst = (double*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int element_size = element_size_in_byte / sizeof(double);
    int * adr_array = (int*)(((unsigned long*)arg)[4]);
    int adr_type_in_byte = (int)(((unsigned long*)arg)[5]);
    int adr_type = adr_type_in_byte / sizeof(int);
    int adr_stride_in_byte = (int)(((unsigned long*)arg)[6]);
    int adr_stride = adr_stride_in_byte / sizeof(int);
    
    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );    
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int * adr_lm;
    for(i=id_head; i<id_head+data_size_total; i++){
        reply_get = reply_put = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, element_size_in_byte);
        dma(dma_get, (long)((int*)adr_array+i*(adr_type+adr_stride)), (long)(adr_lm));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        dma(dma_get, (long)((double*)data_src+(adr_lm[0]*element_size)), (long)(data_lm));
        dma_wait(&reply_get, 2);
        while (reply_get != 2) {}
        
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, element_size_in_byte);
        dma(dma_put, (long)(((double*)data_dst)+(long)((long)i*(long)element_size)), (long)(data_lm));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
}
#else
void CopyIndirectInverse2(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[LOCAL_MEMORY_SIZE_LIMIT/8];
    int n_element = (int)(((unsigned long*)arg)[0]);
    double * data_src = (double*)((unsigned long*)arg)[1];
    double * data_dst = (double*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int element_size = element_size_in_byte / sizeof(double);
    int * adr_array = (int*)(((unsigned long*)arg)[4]);
    int adr_type_in_byte = (int)(((unsigned long*)arg)[5]);
    int adr_type = adr_type_in_byte / sizeof(int);
    int adr_stride_in_byte = (int)(((unsigned long*)arg)[6]);
    int adr_stride = adr_stride_in_byte / sizeof(int);
    
    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );    
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int n_loop = 0;
    for(i=id_head; i<id_head+data_size_total; i++, n_loop++){
        reply_get = reply_put = 0;
        int adr = *(((int*)adr_array)+i*(adr_stride+adr_type));
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, element_size_in_byte);
        dma(dma_get, (long)(((double*)data_src)+((long)adr*(long)element_size)), (long)(data_lm));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, element_size_in_byte);
        dma(dma_put, (long)(((double*)data_dst)+(long)((long)i*(long)element_size)), (long)(data_lm));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
}
#endif

void CopyDirect(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    //if(my_id==0) printf("direct copy \n");
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[LOCAL_MEMORY_SIZE_LIMIT/8];
    //double * data_lm = (double*)ldm_malloc(LOCAL_MEMORY_SIZE_LIMIT);
    int n_element = (int)(((unsigned long*)arg)[0]);
    double * data_src = (double*)((unsigned long*)arg)[1];
    double * data_dst = (double*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int element_size = element_size_in_byte / sizeof(double);
    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );    
    int data_size_in_byte = data_size_total * element_size_in_byte;
    int data_block_limit = data_size_limit_in_byte / element_size_in_byte;
    /*
    if(my_id == 0){
        printf("data_block_limit*element_size_in_byte=%d \n", data_block_limit*element_size_in_byte);
    }
    */
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int data_size_rest = data_size_total;
    while(data_size_rest > 0 ){
        int data_size_tmp = (data_block_limit < data_size_rest) ? data_block_limit : data_size_rest;
        /*
        for(i=0; i<my_id*1000000; i++){
            asm volatile("nop");
        }
        printf("data_size_tmp=%d \n", data_size_tmp);
        */
        reply_get = reply_put = 0;
#if 0
        athread_get(PE_MODE, (((double*)data_src)+id_head*element_size), data_lm, data_size_tmp*element_size_in_byte, (void *)&reply_get, 0, 0, 0);
        while(reply_get < 1) {}
        athread_put(PE_MODE, data_lm, (((double*)data_dst)+id_head*element_size), data_size_tmp*element_size_in_byte, (void *)&reply_put, 0, 0);
        while(reply_put < 1) {}
#else
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, data_size_tmp*element_size_in_byte);
        
        //dma_set_bsize(&dma_get, element_size_in_byte);
        //dma_set_stepsize(&dma_get, 0);
        
        dma(dma_get, (long)(((double*)data_src)+id_head*element_size), (long)(data_lm));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, data_size_tmp*element_size_in_byte);
        dma(dma_put, (long)(((double*)data_dst)+id_head*element_size), (long)(data_lm));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
#endif
        id_head += data_size_tmp;
        data_size_rest -= data_size_tmp;
    }
    /*
    for(i=0; i<my_id*1000000; i++){
        asm volatile("nop");
    }
    printf("loop end, my_id=%d \n", my_id);
    */
    //ldm_free(data_lm, LOCAL_MEMORY_SIZE_LIMIT);
}


#if 0
void CopyStride(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[LOCAL_MEMORY_SIZE_LIMIT/8];
    int n_element = (int)(((unsigned long*)arg)[0]);
    double * data_src = (double*)((unsigned long*)arg)[1];
    double * data_dst = (double*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int stride_size_in_byte = (int)(((unsigned long*)arg)[4]);
    int element_size = element_size_in_byte / sizeof(double);
    int stride_size = stride_size_in_byte / sizeof(double);
    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );
    int data_size_in_byte = data_size_total * element_size_in_byte;
    int data_block_limit = data_size_limit_in_byte / (element_size_in_byte+stride_size_in_byte);
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;    
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int data_size_rest = data_size_total;
    while(data_size_rest > 0 ){
        int data_size_tmp = (data_block_limit < data_size_rest) ? data_block_limit : data_size_rest;
        reply_get = reply_put = 0;

        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, data_size_tmp*(element_size_in_byte+stride_size_in_byte));
        dma_set_bsize(&dma_get, element_size_in_byte);
        dma_set_stepsize(&dma_get, stride_size_in_byte);
        dma(dma_get, (long)(((double*)data_src)+id_head*(element_size+stride_size)), (long)(data_lm));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, data_size_tmp*element_size_in_byte);
        dma(dma_put, (long)(((double*)data_dst)+id_head*element_size), (long)(data_lm));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}

        id_head += data_size_tmp;
        data_size_rest -= data_size_tmp;
    }
}

void CopyStrideByHand(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    int data_lm[LOCAL_MEMORY_SIZE_LIMIT/4];
    int n_element = (int)(((unsigned long*)arg)[0]);
    int * data_src = (int*)((unsigned long*)arg)[1];
    int * data_dst = (int*)((unsigned long*)arg)[2];
    int element_size_in_byte = (int)(((unsigned long*)arg)[3]);
    int stride_size_in_byte = (int)(((unsigned long*)arg)[4]);
    int element_size = element_size_in_byte / sizeof(int);
    int stride_size = stride_size_in_byte / sizeof(int);
    int data_size_total = n_element/NUMBER_OF_CPE + ( (my_id < n_element%NUMBER_OF_CPE) ? 1 : 0 );
    int id_head = (n_element/NUMBER_OF_CPE)*my_id + ( (my_id < n_element%NUMBER_OF_CPE) ? my_id : n_element%NUMBER_OF_CPE );
    int data_size_in_byte = data_size_total * element_size_in_byte;
    int data_block_limit = data_size_limit_in_byte / (element_size_in_byte+stride_size_in_byte);
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;    
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int data_size_rest = data_size_total;
    while(data_size_rest > 0 ){
        int data_size_tmp = (data_block_limit < data_size_rest) ? data_block_limit : data_size_rest;
        reply_get = reply_put = 0;

        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, data_size_tmp*(element_size_in_byte+stride_size_in_byte));
        dma_set_bsize(&dma_get, element_size_in_byte);
        dma_set_stepsize(&dma_get, stride_size_in_byte);
        dma(dma_get, (long)(((int*)data_src)+id_head*(element_size+stride_size)), (long)(data_lm));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, data_size_tmp*element_size_in_byte);
        dma(dma_put, (long)(((int*)data_dst)+id_head*element_size), (long)(data_lm));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}

        id_head += data_size_tmp;
        data_size_rest -= data_size_tmp;
    }
}
#endif

void CopyTCMomToSPJ(void * args){
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int num_tc = (int)(((unsigned long*)args)[0]);
    int common_offset = (int)(((unsigned long*)args)[1]);
    void * adr_tc_head = (void *)((unsigned long*)args)[2];
    void * adr_spj_head =(void *)((unsigned long*)args)[3];
#if 0
    //* Data check (1)
    if (my_id == 0) {
       printf("num_tc = %d\n",num_tc);
       printf("common_offset = %d\n",common_offset);
    }
#endif
    //* Compute the task of each CPE
    int num_tc_loc,my_offset;
    num_tc_loc = num_tc/NUMBER_OF_CPE + ( (my_id < num_tc % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (num_tc/NUMBER_OF_CPE)*my_id + ( (my_id < num_tc % NUMBER_OF_CPE) ? my_id : num_tc % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tc  = sizeof(tcLM);
    size_t bsize_spj = sizeof(spjLM);
    tcLM *tc = (tcLM *) ldm_malloc( bsize_tc * CHUNK_SIZE );
    spjLM *spj = (spjLM *) ldm_malloc( bsize_spj * CHUNK_SIZE );
    //* Perform copy
    // (loop counters)
    int i,j;
    // (local vars. for DMA)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    for (i=0; i<num_tc_loc; i+=CHUNK_SIZE) {
        int nrem = num_tc_loc - i;
        int nn = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get tc
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tc * nn);
        dma(dma_get, (long*)((tcLM *)adr_tc_head + my_offset + i), (long*)(tc));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
    
        //** Inline expansion of the following codes:
        //   spj_sorted_loc_[i].copyFromMoment(tc_loc_[i].mom_); @addMomentAsSpLocalTreeImpl
        //   spj_sorted_[n_spj_prev+i].copyFromMoment(tc_glb_[i].mom_); @addMomentAsSpGlobalTreeImpl
        //   (Note that common_offset corresponds to n_spj_prev)
        for (j=0; j<nn; j++) {
           spj[j].pos.x = tc[j].pos.x;
           spj[j].pos.y = tc[j].pos.y;
           spj[j].pos.z = tc[j].pos.z;
           spj[j].mass  = tc[j].mass;
#ifdef PHI_R_TREE
           spj[j].pos_phi = tc[j].pos_phi;
           spj[j].pos_r   = tc[j].pos_r;
#endif
        }
        
        //** Put spj 
        reply_put = 0; 
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_spj * nn);
        dma(dma_put, (long*)((spjLM *)adr_spj_head + common_offset + my_offset + i), (long*)(spj));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
   
    //* Release memory
    ldm_free(tc, bsize_tc * CHUNK_SIZE );
    ldm_free(spj, bsize_spj * CHUNK_SIZE );
 
}


void CopyETCMomToSPJ(void * args){
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int num_tc = (int)(((unsigned long*)args)[0]);
    int common_offset = (int)(((unsigned long*)args)[1]);
    void * adr_tc_head = (void *)((unsigned long*)args)[2];
    void * adr_spj_head =(void *)((unsigned long*)args)[3];
    //* Compute the task of each CPE
    int num_tc_loc,my_offset;
    num_tc_loc = num_tc/NUMBER_OF_CPE + ( (my_id < num_tc % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (num_tc/NUMBER_OF_CPE)*my_id + ( (my_id < num_tc % NUMBER_OF_CPE) ? my_id : num_tc % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tc  = sizeof(etcLM);
    size_t bsize_spj = sizeof(spjLM);
    etcLM *tc = (etcLM *) ldm_malloc( bsize_tc * CHUNK_SIZE );
    spjLM *spj = (spjLM *) ldm_malloc( bsize_spj * CHUNK_SIZE );
    //* Perform copy
    // (loop counters)
    int i,j;
    // (local vars. for DMA)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    for (i=0; i<num_tc_loc; i+=CHUNK_SIZE) {
        int nrem = num_tc_loc - i;
        int nn = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get tc
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tc * nn);
        dma(dma_get, (long*)((etcLM *)adr_tc_head + my_offset + i), (long*)(tc));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
    
        //** Inline expansion of the following codes:
        //   spj_sorted_loc_[i].copyFromMoment(tc_loc_[i].mom_); @addMomentAsSpLocalTreeImpl
        //   spj_sorted_[n_spj_prev+i].copyFromMoment(tc_glb_[i].mom_); @addMomentAsSpGlobalTreeImpl
        //   (Note that common_offset corresponds to n_spj_prev)
        for (j=0; j<nn; j++) {
           spj[j].pos.x = tc[j].pos.x;
           spj[j].pos.y = tc[j].pos.y;
           spj[j].pos.z = tc[j].pos.z;
           spj[j].mass  = tc[j].mass;
        }
        
        //** Put spj 
        reply_put = 0; 
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_spj * nn);
        dma(dma_put, (long*)((spjLM *)adr_spj_head + common_offset + my_offset + i), (long*)(spj));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
   
    //* Release memory
    ldm_free(tc, bsize_tc * CHUNK_SIZE );
    ldm_free(spj, bsize_spj * CHUNK_SIZE );
 
}

void CopyEpjOrgToEpiSortedEpjSorted(void *args){
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int n_tot_             = (int)(((unsigned long*)args)[0]);
    void * adr_epj_org_    = (void *)((unsigned long*)args)[1];
    void * adr_epi_sor_ = (void *)((unsigned long*)args)[2];
    void * adr_epj_sor_ = (void *)((unsigned long*)args)[3];
    void * adr_tp_loc_     = (void *)((unsigned long*)args)[4];
    //* Compute the task of each CPE
    int n_loc = n_tot_/NUMBER_OF_CPE + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot_/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? my_id : n_tot_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_epi       = sizeof(epiLM);
    size_t bsize_epj       = sizeof(epjLM);
    size_t bsize_tp        = sizeof(tpLM);
    size_t bsize_epi_array = bsize_epi * CHUNK_SIZE;
    size_t bsize_epj_array = bsize_epj * CHUNK_SIZE;
    size_t bsize_tp_array  = bsize_tp  * CHUNK_SIZE;
    epjLM * epj_org = (epjLM *) ldm_malloc( bsize_epj_array );
    epiLM * epi_sor = (epiLM *) ldm_malloc( bsize_epi_array );
    epjLM * epj_sor = (epjLM *) ldm_malloc( bsize_epj_array );
    tpLM  * tp_loc  = (tpLM *) ldm_malloc( bsize_tp_array );
    //int   * adr_arr = (int *) ldm_malloc( sizeof(int)*CHUNK_SIZE );
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp*nn);
        dma(dma_get, (long*)(((tpLM *)adr_tp_loc_)+my_offset+i), (long*)(tp_loc));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_epj);
        for(j=0; j<nn; j++){
            int adr = tp_loc[j].adr_ptcl_;
            dma(dma_get, (long*)(((epjLM*)adr_epj_org_)+adr), (long*)(((epjLM*)epj_org)+j));
        }
        dma_wait(&reply_get, nn);
        while (reply_get != nn) {}
        
        for(j=0; j<nn; j++){
            epj_sor[j].id      = epi_sor[j].id      = epj_org[j].id;
            epj_sor[j].mass    = epi_sor[j].mass    = epj_org[j].mass;
            epj_sor[j].pos.x   = epi_sor[j].pos.x   = epj_org[j].pos.x;
            epj_sor[j].pos.y   = epi_sor[j].pos.y   = epj_org[j].pos.y;
            epj_sor[j].pos.z   = epi_sor[j].pos.z   = epj_org[j].pos.z;
            epj_sor[j].vel.x   = epi_sor[j].vel.x   = epj_org[j].vel.x;
            epj_sor[j].vel.y   = epi_sor[j].vel.y   = epj_org[j].vel.y;
            epj_sor[j].vel.z   = epi_sor[j].vel.z   = epj_org[j].vel.z;
        }
        
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_epj*nn);
        dma(dma_put, (long*)(((epjLM *)adr_epj_sor_)+my_offset+i), (long*)(epj_sor));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}

        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_epi*nn);
        dma(dma_put, (long*)(((epiLM *)adr_epi_sor_)+my_offset+i), (long*)(epi_sor));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}        
    }
    //* Release memory
    ldm_free(epj_org, bsize_epj_array);
    ldm_free(epi_sor, bsize_epi_array);
    ldm_free(epj_sor, bsize_epj_array);
    ldm_free(tp_loc,  bsize_tp_array);
    //ldm_free(adr_arr, sizeof(int)*CHUNK_SIZE);
}

#if 1

void CopyEpiToFpWithCoordinateTransformCPE(void * args){
    int my_id = athread_get_id(-1);
    int n_tot_      = (int)(((unsigned long*)args)[0]);
    void * adr_epi_ = (void *)((unsigned long*)args)[1];
    void * adr_fp_  = (void *)((unsigned long*)args)[2];

    //* Compute the task of each CPE
    int n_loc = n_tot_/NUMBER_OF_CPE + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot_/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? my_id : n_tot_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_fp       = sizeof(fpLM);
    size_t bsize_epi       = sizeof(epiLM);
    size_t bsize_fp_array = bsize_fp * CHUNK_SIZE;
    size_t bsize_epi_array = bsize_epi * CHUNK_SIZE;
    epiLM * epi = (epiLM *) ldm_malloc( bsize_epi_array );
    fpLM *  fp =  (fpLM *)  ldm_malloc( bsize_fp_array );

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

    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_epi*nn);
        dma(dma_get, (long*)(((epiLM *)adr_epi_)+my_offset+i), (long*)(epi));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        for(j=0; j<nn; j++){
            const double x = epi[j].pos.x;
            const double y = epi[j].pos.y;
            const double z = epi[j].pos.z;
            double phi = __ieee754_atan2(y, x, atanhi, atanlo, aT);
            if(phi < 0.0) phi += MY_PI_2;
            const double r = sqrt(x*x+y*y);
            fp[j].pos.x = phi;
            fp[j].pos.y = r;
            fp[j].pos.z = z;
            fp[j].vel.x = epi[j].vel.x;
            fp[j].vel.y = epi[j].vel.y;
            fp[j].vel.z = epi[j].vel.z;
            fp[j].id    = epi[j].id;
            fp[j].mass  = epi[j].mass;
        }

        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_fp*nn);
        dma(dma_put, (long*)(((fpLM *)adr_fp_)+my_offset+i), (long*)(fp));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    ldm_free(epi, bsize_epi_array);
    ldm_free(fp,  bsize_fp_array);
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
}

#elif 1

void CopyEpiToFpWithCoordinateTransformCPE(void * args){
    int my_id = athread_get_id(-1);
    int n_tot_      = (int)(((unsigned long*)args)[0]);
    void * adr_epi_ = (void *)((unsigned long*)args)[1];
    void * adr_fp_  = (void *)((unsigned long*)args)[2];

    //* Compute the task of each CPE
    int n_loc = n_tot_/NUMBER_OF_CPE + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot_/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? my_id : n_tot_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_fp       = sizeof(fpLM);
    size_t bsize_epi       = sizeof(epiLM);
    size_t bsize_fp_array = bsize_fp * CHUNK_SIZE;
    size_t bsize_epi_array = bsize_epi * CHUNK_SIZE;
    epiLM * epi = (epiLM *) ldm_malloc( bsize_epi_array );
    fpLM *  fp =  (fpLM *)  ldm_malloc( bsize_fp_array );
    
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    double one   = 1.0;
    double huge   = 1.0e300;

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

    double y_x;
    double w,s1,s2,z_tmp;
    int ix,hx,id;
    
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_epi*nn);
        dma(dma_get, (long*)(((epiLM *)adr_epi_)+my_offset+i), (long*)(epi));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        for(j=0; j<nn; j++){
            const double x = epi[j].pos.x;
            const double y = epi[j].pos.y;
            const double z = epi[j].pos.z;

            double y_x = y / x;
            double phi = 0.0;

            hx = *(1+(int*)&y_x);
            ix = hx&0x7fffffff;
            if(ix>=0x44100000) {	// if |x| >= 2^66
                if(ix>0x7ff00000||
                   (ix==0x7ff00000&&(*(int*)&y_x!=0)))
                    phi = y_x+y_x;		//
                if(hx>0) phi = atanhi[3]+atanlo[3];
                else     phi = -atanhi[3]-atanlo[3];
                id = -2;
            }
            else if (ix < 0x3fdc0000) {	// |x| < 0.4375
                if (ix < 0x3e200000) {	// |x| < 2^-29 
                    if(huge+x>one) phi = y_x;	// raise inexact
                }
                id = -1;
            }
            else {
                y_x = fabs(y_x);
                if (ix < 0x3ff30000) {		// |x| < 1.1875
                    if (ix < 0x3fe60000) {	// 7/16 <=|x|<11/16
                        id = 0;
                        y_x = (2.0*y_x-one)/(2.0+y_x); 
                    }
                    else {			// 11/16<=|x|< 19/16
                        id = 1;
                        y_x  = (y_x-one)/(y_x+one); 
                    }
                }
                else {
                    if (ix < 0x40038000) {	// |x| < 2.4375
                        id = 2;
                        y_x  = (y_x-1.5)/(one+1.5*y_x);
                    }
                    else {			// 2.4375 <= |x| < 2^66
                        id = 3;
                        y_x  = -1.0/y_x;
                    }
                }
            }
            /* end of argument reduction */
            z_tmp = y_x*y_x;
            w = z_tmp*z_tmp;
            /* break sum from i=0 to 10 aT[i]z**(i+1) into odd and even poly */
            s1 = z_tmp*(aT[0]+w*(aT[2]+w*(aT[4]+w*(aT[6]+w*(aT[8]+w*aT[10])))));
            s2 = w*(aT[1]+w*(aT[3]+w*(aT[5]+w*(aT[7]+w*aT[9]))));
            //s1 = z_tmp*(aT0+w*(aT2+w*(aT4+w*(aT6+w*(aT8+w*aT10)))));
            //s2 = w*(aT1+w*(aT3+w*(aT5+w*(aT7+w*aT9))));
            if (id == -1) phi = y_x - y_x*(s1+s2);
            else if(id >= 0){
                z_tmp = atanhi[id] - ((y_x*(s1+s2) - atanlo[id]) - y_x);
                phi = (hx<0)? -z_tmp:z_tmp;
            }
            if(x < 0.0){
                if(y >= 0.0){
                    phi = PI + phi;
                }
                else{
                    phi = phi - PI;
                }
            }
            
            //double phi = __ieee754_atan2(y, x);
            //double phi = atan(y/x);
            if(phi < 0.0) phi += 2.0*PI;
            else if(phi >= 2.0*PI) phi -= 2.0*PI;
            
            const double r = sqrt(x*x+y*y);
            fp[j].pos.x = phi;
            fp[j].pos.y = r;
            fp[j].pos.z = z;
            fp[j].vel.x = epi[j].vel.x;
            fp[j].vel.y = epi[j].vel.y;
            fp[j].vel.z = epi[j].vel.z;
            fp[j].id    = epi[j].id;
            fp[j].mass  = epi[j].mass;
        }

        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_fp*nn);
        dma(dma_put, (long*)(((fpLM *)adr_fp_)+my_offset+i), (long*)(fp));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    ldm_free(epi, bsize_epi_array);
    ldm_free(fp,  bsize_fp_array);
    #if 1
    ldm_free(atanhi, sizeof(double)*4);
    ldm_free(atanlo, sizeof(double)*4);
    ldm_free(aT, sizeof(double)*11);
    #endif
    
}

#else
void CopyEpiToFpWithCoordinateTransformCPE(void * args){
    int my_id = athread_get_id(-1);
    int n_tot_      = (int)(((unsigned long*)args)[0]);
    void * adr_epi_ = (void *)((unsigned long*)args)[1];
    void * adr_fp_  = (void *)((unsigned long*)args)[2];
    double PI = 3.14159265358979323846;
    //* Compute the task of each CPE
    int n_loc = n_tot_/NUMBER_OF_CPE + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot_/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? my_id : n_tot_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_fp       = sizeof(fpLM);
    size_t bsize_epi       = sizeof(epiLM);
    size_t bsize_fp_array = bsize_fp * CHUNK_SIZE;
    size_t bsize_epi_array = bsize_epi * CHUNK_SIZE;
    epiLM * epi = (epiLM *) ldm_malloc( bsize_epi_array );
    fpLM *  fp =  (fpLM *)  ldm_malloc( bsize_fp_array );
    
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_epi*nn);
        dma(dma_get, (long*)(((epiLM *)adr_epi_)+my_offset+i), (long*)(epi));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        for(j=0; j<nn; j++){
            const double x = epi[j].pos.x;
            const double y = epi[j].pos.y;
            const double z = epi[j].pos.z;
            //double phi = __ieee754_atan2(y, x);
            double phi = atan(y/x);
            //if(phi < 0.0) phi += 2.0*PI;
            //else if(phi >= 2.0*PI) phi -= 2.0*PI;
            //double phi = 0.0;
            const double r = sqrt(x*x+y*y);
            fp[j].pos.x = phi;
            fp[j].pos.y = r;
            fp[j].pos.z = z;
            fp[j].vel.x = epi[j].vel.x;
            fp[j].vel.y = epi[j].vel.y;
            fp[j].vel.z = epi[j].vel.z;
            fp[j].id    = epi[j].id;
            fp[j].mass  = epi[j].mass;
        }

        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_fp*nn);
        dma(dma_put, (long*)(((fpLM *)adr_fp_)+my_offset+i), (long*)(fp));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    ldm_free(epi, bsize_epi_array);
    ldm_free(fp,  bsize_fp_array);    
}
#endif

void CopyFPToEPJ(void *args){
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int n_tot_      = (int)(((unsigned long*)args)[0]);
    void * adr_fp_  = (void *)((unsigned long*)args)[1];
    void * adr_epj_ = (void *)((unsigned long*)args)[2];
    //* Compute the task of each CPE
    int n_loc = n_tot_/NUMBER_OF_CPE + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot_/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? my_id : n_tot_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_fp        = sizeof(fpLM);
    size_t bsize_epj       = sizeof(epjLM);
    size_t bsize_fp_array  = bsize_fp * CHUNK_SIZE;
    size_t bsize_epj_array = bsize_epj * CHUNK_SIZE;
    fpLM  *fp  = (fpLM *) ldm_malloc( bsize_fp_array );
    epjLM *epj = (epjLM *) ldm_malloc( bsize_epj_array );
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get fp
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get,  bsize_fp*nn);
        dma(dma_get, (long*)(((fpLM *)adr_fp_)+my_offset+i), (long*)(fp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        for(j=0; j<nn; j++){
            epj[j].pos.x  = fp[j].pos.x;
            epj[j].pos.y  = fp[j].pos.y;
            epj[j].pos.z  = fp[j].pos.z;
            epj[j].mass = fp[j].mass;
            epj[j].vel.x  = fp[j].vel.x;
            epj[j].vel.y  = fp[j].vel.y;
            epj[j].vel.z  = fp[j].vel.z;
            epj[j].id   = fp[j].id;
        }
        
        //** Put into epj_sorted_ 
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_epj*nn);
        dma(dma_put, (long*)(((epjLM *)adr_epj_)+my_offset+i), (long*)(epj));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    //* Release memory
    ldm_free(fp,  bsize_fp_array);
    ldm_free(epj, bsize_epj_array);
}


void CopyEPIToFP(void *args){
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int n_tot_      = (int)((unsigned long*)args)[0];
    void * adr_epi_ = (void *)((unsigned long*)args)[1];
    void * adr_fp_  = (void *)((unsigned long*)args)[2];
    //* Compute the task of each CPE
    int n_loc = n_tot_/NUMBER_OF_CPE + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot_/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? my_id : n_tot_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_epi       = sizeof(epiLM);
    size_t bsize_fp        = sizeof(fpLM);
    size_t bsize_epi_array = bsize_epi * CHUNK_SIZE;
    size_t bsize_fp_array  = bsize_fp  * CHUNK_SIZE;
    fpLM  * fp  = (fpLM *)  ldm_malloc( bsize_fp_array );
    epiLM * epi = (epiLM *) ldm_malloc( bsize_epi_array );
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k;
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get epi
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get,  bsize_epi*nn);
        dma(dma_get, (long*)(((epiLM *)adr_epi_)+my_offset+i), (long*)(epi));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for(j=0; j<nn; j++){
            fp[j].pos.x = epi[j].pos.x;
            fp[j].pos.y = epi[j].pos.y;
            fp[j].pos.z = epi[j].pos.z;
            fp[j].mass  = epi[j].mass;
            fp[j].vel.x = epi[j].vel.x;
            fp[j].vel.y = epi[j].vel.y;
            fp[j].vel.z = epi[j].vel.z;
            fp[j].id  = epi[j].id;
        }
        //** Put fp
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_fp*nn);
        dma(dma_put, (long*)(((fpLM *)adr_fp_)+my_offset+i), (long*)(fp));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    //* Release memory
    ldm_free(fp,  bsize_fp_array);
    ldm_free(epi, bsize_epi_array);
}
