#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
//#include<athread.h>
#include"slave.h"
#include"dma.h"
#include"cpe_func.h"

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

#if 0
void CopyToGlobalTree(void * arg){
    int i;
    int my_id = athread_get_id(-1);
    int data_size_limit_in_byte = LOCAL_MEMORY_SIZE_LIMIT;
    double data_lm[LOCAL_MEMORY_SIZE_LIMIT/8];
    int n_glb_tot = (int)(((unsigned long*)arg)[0]);
    int n_loc_tot = (int)(((unsigned long*)arg)[1]);
    double * epj_org = (double*)((unsigned long*)arg)[2];
    int epj_size_in_byte = (int)(((unsigned long*)arg)[3]);
    double * spj_org = (double*)((unsigned long*)arg)[4];
    int spj_size_in_byte = (int)(((unsigned long*)arg)[5]);
    int epj_size = epj_size_in_byte / sizeof(double);
    int spj_size = spj_size_in_byte / sizeof(double);
    
    int * adr_array = (int*)(((unsigned long*)arg)[4]);
    int adr_stride_in_byte = (int)(((unsigned long*)arg)[5]);
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
        int adr = *(((int*)adr_array)+i*adr_stride);
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
