#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
//#include<athread.h>
#include"slave.h"
#include"dma.h"
#include"cpe_func.h"
#include"cpe_prof.h"

/* Local macros */
//#define CHECK_MSORT_GLB_TREE_REUSE
#define RANK_CHECK (0)
#define ID_CHECK   (62)

#define PROFILE_MSORT_GLB_TREE_REUSE

#if VERSION_MSORT_GLB_TREE_REUSE == 0
void Preproc1_CopyEPJLocToEPJGlb(void *args) {} // dummy
void Preproc2_CopyEPJLocToEPJGlb(void *args) {} // dummy
void CopyEPJLocToEPJGlb(void *args){ 
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
#elif VERSION_MSORT_GLB_TREE_REUSE == 1
void Preproc1_CopyEPJLocToEPJGlb(void *args) {
    // This function counts the number of parts of adr_epj_loc2glb_
    // where the addresses to epj_sorted_ lines up continously.
    unsigned long wt_start0 = rtc_();
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int myrank                 = (int)(((unsigned long*)args)[0]);
    int n_epj_sorted_loc_      = (int)(((unsigned long*)args)[1]);
    void *adr_adr_epj_loc2glb_ = (void *)((unsigned long*)args)[2];
    void *adr_n_groups         = (void *)((unsigned long*)args)[3];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_epj_sorted_loc_/NUMBER_OF_CPE + ( (my_id < n_epj_sorted_loc_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_epj_sorted_loc_/NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted_loc_ % NUMBER_OF_CPE) ? my_id : n_epj_sorted_loc_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_adr_array = sizeof(int) * CHUNK_SIZE;
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
    // (group information)
    int n_group = 0;
    int adr_glb_prev = -1;
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get adr_epj_loc2glb_
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, sizeof(int) * nn);
        dma(dma_get, (long*)((int *)adr_adr_epj_loc2glb_ + my_offset + i), (long*)(adr));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        //** Count the number of groups
        for (k=0; k<nn; k++) {
            if (adr[k] - adr_glb_prev != 1) n_group++;
            // In this case, adr_loc = my_offset + i + k is the head of a new group.
            adr_glb_prev = adr[k]; 
        }
    }

    //* Put the number of groups
    reply_put = 0;
    dma_set_op(&dma_put, DMA_PUT);
    dma_set_mode(&dma_put, PE_MODE);
    dma_set_reply(&dma_put, &reply_put);
    dma_set_size(&dma_put, sizeof(int));
    dma(dma_put, (long*)((int *)adr_n_groups + my_id), (long*)(&n_group));
    dma_wait(&reply_put, 1);
    while (reply_put != 1) {}

    //* Release memory
    ldm_free(adr, bsize_adr_array);

#ifdef CHECK_MSORT_GLB_TREE_REUSE
#if 0
    //* Check n_loc & my_offset
    if (myrank == RANK_CHECK) {
        sync_array_();
        for (i=0; i<NUMBER_OF_CPE; ++i) {
            if (i == my_id)
                printf("cpe_id = %d, n_loc = %d, my_offset = %d\n",
                       my_id,n_loc,my_offset);
            sync_array_();
        }
    }
#endif
#endif

#ifdef PROFILE_MSORT_GLB_TREE_REUSE
    //* Output profile
    unsigned long wt_end0 = rtc_();
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        double bw_dma_get=0.363e9, bw_dma_put=0.363e9; // bandwidth per CPE
        double cpu_clock=1.45e9;
        printf("Actual total time = %e [s]\n",(wt_end0-wt_start0)/cpu_clock);
    }
#endif

}

void Preproc2_CopyEPJLocToEPJGlb(void *args) {
    // This function constructs the group information used in the copy
    // from epj_sorted_loc_ to epj_sorted_.
    unsigned long wt_start0 = rtc_();
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int myrank                        = (int)(((unsigned long*)args)[0]);
    int n_epj_sorted_loc_             = (int)(((unsigned long*)args)[1]);
    void *adr_adr_epj_loc2glb_        = (void *)((unsigned long*)args)[2];
    void *adr_adr_epj_loc_group_head_ = (void *)((unsigned long*)args)[3];
    void *adr_adr_epj_glb_group_head_ = (void *)((unsigned long*)args)[4];
    void *adr_group_size_epj_loc_     = (void *)((unsigned long*)args)[5];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_epj_sorted_loc_/NUMBER_OF_CPE + ( (my_id < n_epj_sorted_loc_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_epj_sorted_loc_/NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted_loc_ % NUMBER_OF_CPE) ? my_id : n_epj_sorted_loc_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
        BUF_SIZE = 3584,
    };
    size_t bsize_adr_array = sizeof(int) * CHUNK_SIZE;
    size_t bsize_buf_array = sizeof(int) * BUF_SIZE;
    int *adr                    = (int *) ldm_malloc( bsize_adr_array );
    int *adr_epj_loc_group_head = (int *) ldm_malloc( bsize_buf_array );
    int *adr_epj_glb_group_head = (int *) ldm_malloc( bsize_buf_array );
    int *group_size_epj_loc     = (int *) ldm_malloc( bsize_buf_array );
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (group information)
    int group_id = -1;
    int n_group_member = 1;
    int adr_glb_prev = -1;
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get adr_epj_loc2glb_
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, sizeof(int) * nn);
        dma(dma_get, (long*)((int *)adr_adr_epj_loc2glb_ + my_offset + i), (long*)(adr));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        //** Count the number of groups
        for (k=0; k<nn; k++) {
            if (adr[k] - adr_glb_prev != 1) {
                // In this case, adr_loc = my_offset + i is the head of a new group.
                group_id++; // new group id
                assert(group_id < BUF_SIZE);
                adr_epj_loc_group_head[group_id] = my_offset + i + k;
                adr_epj_glb_group_head[group_id] = adr[k];
                if (group_id > 0) 
                    group_size_epj_loc[group_id-1] = n_group_member;
                n_group_member = 1; // reset # of group members
            } else {
                n_group_member++;
            }
            adr_glb_prev = adr[k];
        }
        group_size_epj_loc[group_id] = n_group_member;
    }
    //* for debug
    //if (myrank == 0 && my_id == 21) {
    //    k = 89;
    //    printf("----------");
    //    printf("group_id = %d\n",k);
    //    printf("adr_epj_loc_group_head = %d\n",adr_epj_loc_group_head[k]);
    //    printf("adr_epj_glb_group_head = %d\n",adr_epj_glb_group_head[k]);
    //    printf("group_size_epj_loc     = %d\n",group_size_epj_loc[k]);
    //}

    //* Get the address of DMA put
    unsigned long adr_dma_put[3];
    //** DMA setting 
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(unsigned long));
    //** (1) the address of adr_epj_loc_group_head_[my_id]
    dma(dma_get, (long*)((unsigned long *)adr_adr_epj_loc_group_head_ + my_id), (long*)(&adr_dma_put[0]));
    //** (2) the address of adr_epj_glb_group_head_[my_id]
    dma(dma_get, (long*)((unsigned long *)adr_adr_epj_glb_group_head_ + my_id), (long*)(&adr_dma_put[1]));
    //** (3) the address of group_size_epj_loc_[my_id]
    dma(dma_get, (long*)((unsigned long *)adr_group_size_epj_loc_ + my_id), (long*)(&adr_dma_put[2]));
    //** Wait all
    dma_wait(&reply_get, 3);
    while (reply_get != 3) {}
#if 0
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        printf("adr_dma_put[0] = %lx\n",adr_dma_put[0]);
        printf("adr_dma_put[1] = %lx\n",adr_dma_put[1]);
        printf("adr_dma_put[2] = %lx\n",adr_dma_put[2]);
    }
#endif


    //* Put the group information
    //** DMA setting
    reply_put = 0;
    dma_set_op(&dma_put, DMA_PUT);
    dma_set_mode(&dma_put, PE_MODE);
    dma_set_reply(&dma_put, &reply_put);
    dma_set_size(&dma_put, sizeof(int) * (group_id+1)); // n_group = group_id + 1
    //** (1) Put to adr_epj_loc_group_head_[my_id]
    dma(dma_put, (long)(adr_dma_put[0]), (long*)(adr_epj_loc_group_head));
    //** (2) Put to adr_epj_glb_group_head_[my_id]
    dma(dma_put, (long)(adr_dma_put[1]), (long*)(adr_epj_glb_group_head));
    //** (3) Put to group_size_epj_loc_[my_id]
    dma(dma_put, (long)(adr_dma_put[2]), (long*)(group_size_epj_loc));
    //** Wait all
    dma_wait(&reply_put, 3);
    while (reply_put != 3) {}

    //* Release memory
    ldm_free(adr, bsize_adr_array);
    ldm_free(adr_epj_loc_group_head, bsize_buf_array);
    ldm_free(adr_epj_glb_group_head, bsize_buf_array);
    ldm_free(group_size_epj_loc, bsize_buf_array);

#ifdef PROFILE_MSORT_GLB_TREE_REUSE
    //* Output profile
    unsigned long wt_end0 = rtc_();
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        double bw_dma_get=0.363e9, bw_dma_put=0.363e9; // bandwidth per CPE
        double cpu_clock=1.45e9;
        printf("Actual total time = %e [s]\n",(wt_end0-wt_start0)/cpu_clock);
    }
#endif

}

void CopyEPJLocToEPJGlb(void *args){ 
    unsigned long wt_start0 = rtc_();
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int myrank                        = (int)(((unsigned long*)args)[0]);
    int n_epj_sorted_loc_             = (int)(((unsigned long*)args)[1]);
    void *adr_epj_sorted_loc_         = (void *)((unsigned long*)args)[2];
    void *adr_epj_sorted_             = (void *)((unsigned long*)args)[3];
    void *adr_n_groups                = (void *)((unsigned long*)args)[4];
    void *adr_adr_epj_loc_group_head_ = (void *)((unsigned long*)args)[5];
    void *adr_adr_epj_glb_group_head_ = (void *)((unsigned long*)args)[6];
    void *adr_group_size_epj_loc_     = (void *)((unsigned long*)args)[7];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_epj_sorted_loc_/NUMBER_OF_CPE + ( (my_id < n_epj_sorted_loc_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_epj_sorted_loc_/NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted_loc_ % NUMBER_OF_CPE) ? my_id : n_epj_sorted_loc_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        EPJ_ARRAY_SIZE = 512,
        BUF_ARRAY_SIZE = 1024,
    };
    size_t bsize_epj       = sizeof(epjLM);
    size_t bsize_epj_array = bsize_epj   * EPJ_ARRAY_SIZE;
    size_t bsize_buf_array = sizeof(int) * BUF_ARRAY_SIZE;
    epjLM  *epj = (epjLM *) ldm_malloc( bsize_epj_array );
    int *adr_epj_loc_group_head = (int *) ldm_malloc( bsize_buf_array );
    int *adr_epj_glb_group_head = (int *) ldm_malloc( bsize_buf_array );
    int *group_size_epj_loc     = (int *) ldm_malloc( bsize_buf_array );
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k,ig; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (number of group)
    int n_group;
    // (addresses of group info.)
    unsigned long adr_group_info[3];
    //* Get n_group
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(int));
    dma(dma_get, (long*)((int *)adr_n_groups + my_id), (long*)(&n_group));
    dma_wait(&reply_get, 1);
    while (reply_get != 1) {}
#if 0
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        printf("n_group = %d\n",n_group);
    }
#endif
    //* Get adr_group_info[3]
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(unsigned long));
    dma(dma_get, (long*)((unsigned long *)adr_adr_epj_loc_group_head_ + my_id), (long*)(&adr_group_info[0]));
    dma(dma_get, (long*)((unsigned long *)adr_adr_epj_glb_group_head_ + my_id), (long*)(&adr_group_info[1]));
    dma(dma_get, (long*)((unsigned long *)adr_group_size_epj_loc_ + my_id), (long*)(&adr_group_info[2]));
    dma_wait(&reply_get, 3);
    while (reply_get != 3) {}
#if 0
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        printf("adr_group_info[0] = %lx\n",adr_group_info[0]);
        printf("adr_group_info[1] = %lx\n",adr_group_info[1]);
        printf("adr_group_info[2] = %lx\n",adr_group_info[2]);
    }
#endif
    //* Copy
    for (ig=0; ig<n_group; ig+=BUF_ARRAY_SIZE) {
        int nrem = n_group - ig;
        int nn   = nrem < BUF_ARRAY_SIZE ? nrem : BUF_ARRAY_SIZE;
        //** Get group information
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, sizeof(int) * nn);
        dma(dma_get, (long)(adr_group_info[0] + sizeof(int) * ig), (long*)(adr_epj_loc_group_head));
        dma(dma_get, (long)(adr_group_info[1] + sizeof(int) * ig), (long*)(adr_epj_glb_group_head));
        dma(dma_get, (long)(adr_group_info[2] + sizeof(int) * ig), (long*)(group_size_epj_loc));
        dma_wait(&reply_get, 3);
        while (reply_get != 3) {}
        //** Copy
        for (i=0; i<nn; i++) {
            //** Check the size of epj
            assert(group_size_epj_loc[i] < EPJ_ARRAY_SIZE);
            //** Get epj_sorted_loc_ 
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_epj * group_size_epj_loc[i]);
            dma(dma_get, (long*)((epjLM *)adr_epj_sorted_loc_ + adr_epj_loc_group_head[i]), (long*)(epj));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
#ifdef CHECK_MSORT_GLB_TREE_REUSE
            if (myrank == 0) {
               if (adr_epj_glb_group_head[i] == 344787) {
                   printf("----------\n");
                   printf("my_id   = %d\n",my_id);
                   printf("adr_loc = %d\n",adr_epj_loc_group_head[i]);
                   printf("adr_glb = %d\n",adr_epj_glb_group_head[i]);
                   printf("gsize   = %d\n",group_size_epj_loc[i]);
               }
            }
#endif
            //** For debug
            //** Put into epj_sorted_ 
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, bsize_epj * group_size_epj_loc[i]);
            dma(dma_put, (long*)((epjLM *)adr_epj_sorted_ + adr_epj_glb_group_head[i]), (long*)(epj));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }
    }

    //* Release memory
    ldm_free(epj, bsize_epj_array);
    ldm_free(adr_epj_loc_group_head, bsize_buf_array);
    ldm_free(adr_epj_glb_group_head, bsize_buf_array);
    ldm_free(group_size_epj_loc, bsize_buf_array);

#ifdef PROFILE_MSORT_GLB_TREE_REUSE
    //* Output profile
    sync_array_();
    unsigned long wt_end0 = rtc_();
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        double bw_dma_get=0.363e9, bw_dma_put=0.363e9; // bandwidth per CPE
        double cpu_clock=1.45e9;
        printf("Actual total time = %e [s]\n",(wt_end0-wt_start0)/cpu_clock);
    }
#endif
}
#elif VERSION_MSORT_GLB_TREE_REUSE == 2
void Preproc1_CopyEPJLocToEPJGlb(void *args) {
    // This function counts the number of parts of adr_epj_glb2loc_
    // where the addresses to epj_sorted_loc_ lines up continously.
    unsigned long wt_start0 = rtc_();
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int myrank                 = (int)(((unsigned long*)args)[0]);
    int n_epj_sorted_          = (int)(((unsigned long*)args)[1]);
    void *adr_adr_epj_glb2loc_ = (void *)((unsigned long*)args)[2];
    void *adr_n_groups         = (void *)((unsigned long*)args)[3];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_epj_sorted_/NUMBER_OF_CPE + ( (my_id < n_epj_sorted_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_epj_sorted_/NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted_ % NUMBER_OF_CPE) ? my_id : n_epj_sorted_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_adr_array = sizeof(int) * CHUNK_SIZE;
    int    *adr = (int *)   ldm_malloc( bsize_adr_array );
    //* Count the number of groups on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (group information)
    int n_group = 0;
    int adr_loc_prev = -10; // must be sufficiently larger than -1.
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get adr_epj_glb2loc_
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, sizeof(int) * nn);
        dma(dma_get, (long*)((int *)adr_adr_epj_glb2loc_ + my_offset + i), (long*)(adr));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        //** Count the number of groups
        for (k=0; k<nn; k++) {
            if (((adr_loc_prev == -1) && (adr[k] != -1)) ||
                ((adr_loc_prev != -1) && (adr[k] != adr_loc_prev+1))) n_group++;
            // In this case, adr_glb = my_offset + i + k is the head of a new group.
            adr_loc_prev = adr[k]; 
        }
    }

    //* Put the number of groups
    reply_put = 0;
    dma_set_op(&dma_put, DMA_PUT);
    dma_set_mode(&dma_put, PE_MODE);
    dma_set_reply(&dma_put, &reply_put);
    dma_set_size(&dma_put, sizeof(int));
    dma(dma_put, (long*)((int *)adr_n_groups + my_id), (long*)(&n_group));
    dma_wait(&reply_put, 1);
    while (reply_put != 1) {}

    //* Release memory
    ldm_free(adr, bsize_adr_array);

#ifdef CHECK_MSORT_GLB_TREE_REUSE
#if 0
    //* Check n_loc & my_offset
    if (myrank == RANK_CHECK) {
        sync_array_();
        for (i=0; i<NUMBER_OF_CPE; ++i) {
            if (i == my_id)
                printf("cpe_id = %d, n_loc = %d, my_offset = %d\n",
                       my_id,n_loc,my_offset);
            sync_array_();
        }
    }
#endif
#endif

#ifdef PROFILE_MSORT_GLB_TREE_REUSE
    //* Output profile
    unsigned long wt_end0 = rtc_();
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        double bw_dma_get=0.363e9, bw_dma_put=0.363e9; // bandwidth per CPE
        double cpu_clock=1.45e9;
        printf("Actual total time = %e [s]\n",(wt_end0-wt_start0)/cpu_clock);
    }
#endif

}

void Preproc2_CopyEPJLocToEPJGlb(void *args) {
    // This function constructs the group information used in the copy
    // from epj_sorted_loc_ to epj_sorted_.
    unsigned long wt_start0 = rtc_();
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int myrank                        = (int)(((unsigned long*)args)[0]);
    int n_epj_sorted_                 = (int)(((unsigned long*)args)[1]);
    void *adr_adr_epj_glb2loc_        = (void *)((unsigned long*)args)[2];
    void *adr_adr_epj_glb_group_head_ = (void *)((unsigned long*)args)[3];
    void *adr_adr_epj_loc_group_head_ = (void *)((unsigned long*)args)[4];
    void *adr_group_size_epj_glb_     = (void *)((unsigned long*)args)[5];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_epj_sorted_/NUMBER_OF_CPE + ( (my_id < n_epj_sorted_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_epj_sorted_/NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted_ % NUMBER_OF_CPE) ? my_id : n_epj_sorted_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
        //BUF_SIZE = 3584,
        BUF_SIZE = 5000,
    };
    size_t bsize_adr_array = sizeof(int) * CHUNK_SIZE;
    size_t bsize_buf_array = sizeof(int) * BUF_SIZE;
    int *adr                    = (int *) ldm_malloc( bsize_adr_array );
    int *adr_epj_glb_group_head = (int *) ldm_malloc( bsize_buf_array );
    int *adr_epj_loc_group_head = (int *) ldm_malloc( bsize_buf_array );
    int *group_size_epj_glb     = (int *) ldm_malloc( bsize_buf_array );
    //* Make group info. on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (group information)
    int group_id = -1;
    int n_group_member = 1;
    int adr_loc_prev = -10; // must be sufficiently larger than -1.
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        //** Get adr_epj_loc2glb_
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, sizeof(int) * nn);
        dma(dma_get, (long*)((int *)adr_adr_epj_glb2loc_ + my_offset + i), (long*)(adr));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        //** Make group info.
        for (k=0; k<nn; k++) {
            if (((adr_loc_prev == -1) && (adr[k] != -1)) ||
                ((adr_loc_prev != -1) && (adr[k] != adr_loc_prev+1))) {
                // In this case, adr_glb = my_offset + i is the head of a new group.
                group_id++; // new group id
                assert(group_id < BUF_SIZE);
                adr_epj_glb_group_head[group_id] = my_offset + i + k;
                adr_epj_loc_group_head[group_id] = adr[k];
                if (group_id > 0) 
                    group_size_epj_glb[group_id-1] = n_group_member;
                n_group_member = 1; // reset # of group members
            } else {
                n_group_member++;
            }
            adr_loc_prev = adr[k];
        }
        group_size_epj_glb[group_id] = n_group_member;
    }
#ifdef CHECK_MSORT_GLB_TREE_REUSE
    //* Check if group information is consistent with n_loc
    int gid,sum=0;
    for (gid=0; gid<=group_id; gid++) {
        sum += group_size_epj_glb[gid];
    }
    //printf("sum   = %d\n",sum);
    //printf("n_loc = %d\n",n_loc);
    assert(sum == n_loc);
#endif

    //* Get the address of DMA put
    unsigned long adr_dma_put[3];
    //** DMA setting 
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(unsigned long));
    //** (1) the address of adr_epj_glb_group_head_[my_id]
    dma(dma_get, (long*)((unsigned long *)adr_adr_epj_glb_group_head_ + my_id), (long*)(&adr_dma_put[0]));
    //** (2) the address of adr_epj_loc_group_head_[my_id]
    dma(dma_get, (long*)((unsigned long *)adr_adr_epj_loc_group_head_ + my_id), (long*)(&adr_dma_put[1]));
    //** (3) the address of group_size_epj_glb_[my_id]
    dma(dma_get, (long*)((unsigned long *)adr_group_size_epj_glb_ + my_id), (long*)(&adr_dma_put[2]));
    //** Wait all
    dma_wait(&reply_get, 3);
    while (reply_get != 3) {}
#if 0
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        printf("adr_dma_put[0] = %lx\n",adr_dma_put[0]);
        printf("adr_dma_put[1] = %lx\n",adr_dma_put[1]);
        printf("adr_dma_put[2] = %lx\n",adr_dma_put[2]);
    }
#endif


    //* Put the group information
    //** DMA setting
    reply_put = 0;
    dma_set_op(&dma_put, DMA_PUT);
    dma_set_mode(&dma_put, PE_MODE);
    dma_set_reply(&dma_put, &reply_put);
    dma_set_size(&dma_put, sizeof(int) * (group_id+1)); // n_group = group_id + 1
    //** (1) Put to adr_epj_glb_group_head_[my_id]
    dma(dma_put, (long)(adr_dma_put[0]), (long*)(adr_epj_glb_group_head));
    //** (2) Put to adr_epj_loc_group_head_[my_id]
    dma(dma_put, (long)(adr_dma_put[1]), (long*)(adr_epj_loc_group_head));
    //** (3) Put to group_size_epj_glb_[my_id]
    dma(dma_put, (long)(adr_dma_put[2]), (long*)(group_size_epj_glb));
    //** Wait all
    dma_wait(&reply_put, 3);
    while (reply_put != 3) {}

    //* Release memory
    ldm_free(adr, bsize_adr_array);
    ldm_free(adr_epj_glb_group_head, bsize_buf_array);
    ldm_free(adr_epj_loc_group_head, bsize_buf_array);
    ldm_free(group_size_epj_glb, bsize_buf_array);

#ifdef PROFILE_MSORT_GLB_TREE_REUSE
    //* Output profile
    unsigned long wt_end0 = rtc_();
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        double bw_dma_get=0.363e9, bw_dma_put=0.363e9; // bandwidth per CPE
        double cpu_clock=1.45e9;
        printf("Actual total time = %e [s]\n",(wt_end0-wt_start0)/cpu_clock);
    }
#endif

}

void CopyEPJLocToEPJGlb(void *args){ 
    unsigned long wt_start0 = rtc_();
    int my_id = athread_get_id(-1);
    //* Process the arguments
    int myrank                        = (int)(((unsigned long*)args)[0]);
    int n_epj_sorted_                 = (int)(((unsigned long*)args)[1]);
    void *adr_epj_sorted_loc_         = (void *)((unsigned long*)args)[2];
    void *adr_epj_sorted_             = (void *)((unsigned long*)args)[3];
    void *adr_n_groups                = (void *)((unsigned long*)args)[4];
    void *adr_adr_epj_glb_group_head_ = (void *)((unsigned long*)args)[5];
    void *adr_adr_epj_loc_group_head_ = (void *)((unsigned long*)args)[6];
    void *adr_group_size_epj_glb_     = (void *)((unsigned long*)args)[7];
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_epj_sorted_/NUMBER_OF_CPE + ( (my_id < n_epj_sorted_ % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_epj_sorted_/NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted_ % NUMBER_OF_CPE) ? my_id : n_epj_sorted_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        EPJ_ARRAY_SIZE = 256,
        BUF_ARRAY_SIZE = 1024,
    };
    size_t bsize_epj       = sizeof(epjLM);
    size_t bsize_epj_array = bsize_epj   * EPJ_ARRAY_SIZE;
    size_t bsize_buf_array = sizeof(int) * BUF_ARRAY_SIZE;
    epjLM  *epj = (epjLM *) ldm_malloc( bsize_epj_array );
    int *adr_epj_glb_group_head = (int *) ldm_malloc( bsize_buf_array );
    int *adr_epj_loc_group_head = (int *) ldm_malloc( bsize_buf_array );
    int *group_size_epj_glb     = (int *) ldm_malloc( bsize_buf_array );
#if 0
    if (myrank == 0 && my_id == 0) {
        printf("bsize_epj       = %d\n",bsize_epj);
        printf("bsize_epj_array = %d\n",bsize_epj_array);
        printf("bsize_buf_array = %d\n",bsize_buf_array);
        printf("total           = %d\n",bsize_epj_array+3*bsize_buf_array);
    }
#endif
    //* Compute moments on each CPE
    // (loop counters, etc.)
    int i,j,k;
    int ig,igl,ng;
    int chunk_id,adr_chunk,chunk_size;
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (number of group)
    int n_group;
    // (addresses of group info.)
    unsigned long adr_group_info[3];
    //* Get n_group
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(int));
    dma(dma_get, (long*)((int *)adr_n_groups + my_id), (long*)(&n_group));
    dma_wait(&reply_get, 1);
    while (reply_get != 1) {}
#if 0
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        printf("n_group = %d\n",n_group);
    }
#endif
    //* Get adr_group_info[3]
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(unsigned long));
    dma(dma_get, (long*)((unsigned long *)adr_adr_epj_glb_group_head_ + my_id), (long*)(&adr_group_info[0]));
    dma(dma_get, (long*)((unsigned long *)adr_adr_epj_loc_group_head_ + my_id), (long*)(&adr_group_info[1]));
    dma(dma_get, (long*)((unsigned long *)adr_group_size_epj_glb_ + my_id), (long*)(&adr_group_info[2]));
    dma_wait(&reply_get, 3);
    while (reply_get != 3) {}
#if 0
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        printf("adr_group_info[0] = %lx\n",adr_group_info[0]);
        printf("adr_group_info[1] = %lx\n",adr_group_info[1]);
        printf("adr_group_info[2] = %lx\n",adr_group_info[2]);
    }
#endif

    //* Copy
    //** Initialize the chunk information
    //   where a chunk is defined as a group of array elements of
    //   epj_sorted_[adr_chunk] ~ epj_sorted_[adr_chunk + chunk_size - 1].
    chunk_id   = 0;
    adr_chunk  = my_offset;
    chunk_size = n_loc < EPJ_ARRAY_SIZE ? n_loc : EPJ_ARRAY_SIZE;
    for (ig=0; ig<n_group; ig+=BUF_ARRAY_SIZE) {
        int nrem = n_group - ig;
        ng = nrem < BUF_ARRAY_SIZE ? nrem : BUF_ARRAY_SIZE;
#ifdef CHECK_MSORT_GLB_TREE_REUSE
        if (myrank == RANK_CHECK && my_id == ID_CHECK) {
           printf("#########################\n");
           printf("ig = %d\n",ig);
           printf("ng = %d\n",ng);
        }
#endif
        //** Get group information
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, sizeof(int) * ng);
        dma(dma_get, (long)(adr_group_info[0] + sizeof(int) * ig), (long*)(adr_epj_glb_group_head));
        dma(dma_get, (long)(adr_group_info[1] + sizeof(int) * ig), (long*)(adr_epj_loc_group_head));
        dma(dma_get, (long)(adr_group_info[2] + sizeof(int) * ig), (long*)(group_size_epj_glb));
        dma_wait(&reply_get, 3);
        while (reply_get != 3) {}
#ifdef CHECK_MSORT_GLB_TREE_REUSE
        int sum=0;
        for (i=0; i<ng; i++) sum += group_size_epj_glb[i];
        if (myrank == RANK_CHECK && my_id == ID_CHECK) {
            printf("the scheduled copy size = %d\n",sum);
        }
#endif
        //** Copy using the loaded group information
        igl = 0; // (local) group id; takes from 0 to nn-1.
        for (;;) {
#ifdef CHECK_MSORT_GLB_TREE_REUSE
            if (myrank == RANK_CHECK && my_id == ID_CHECK) {
                printf("=========================\n");
                printf("chunk_id   = %d\n",chunk_id);
                printf("adr_chunk  = %d\n",adr_chunk);
                printf("chunk_size = %d\n",chunk_size);
                printf("adr_plan   = %d\n",my_offset + EPJ_ARRAY_SIZE * chunk_id);
            }
#endif

            //** Get epj_sorted_loc_
            int copied_size = 0;
            for (;;) {
                int nrem = chunk_size - copied_size;
                int group_size, adr_offset;
                if (adr_chunk > adr_epj_glb_group_head[igl]) {
                   // In this case, a part of group i is already copied.
                   // Therefore, we start copy from an intermediate point.
                   adr_offset = adr_chunk - adr_epj_glb_group_head[igl];
                   group_size = group_size_epj_glb[igl] - adr_offset;
                } else {
                   // In this case, we try to copy from the head of the group.
                   group_size = group_size_epj_glb[igl];
                   adr_offset = 0;
                }
                int nn = nrem < group_size ? nrem : group_size; // nn is the number of data actually copied
#ifdef CHECK_MSORT_GLB_TREE_REUSE
                if (myrank == RANK_CHECK && my_id == ID_CHECK) {
                    printf("-------------------------\n");
                    printf("copied_size   = %d\n",copied_size);
                    printf("igl           = %d\n",igl);
                    printf("nrem          = %d\n",nrem);
                    printf("group_size    = %d\n",group_size);
                    printf("group_size[T] = %d\n",group_size_epj_glb[igl]);
                    printf("adr_offset    = %d\n",adr_offset);
                    printf("nn            = %d\n",nn);
                    printf("adr_epj_ghead = %d\n",adr_epj_glb_group_head[igl]);
                }
#endif
                if (adr_epj_loc_group_head[igl] != -1) {
                    reply_get = 0;
                    dma_set_op(&dma_get, DMA_GET);
                    dma_set_mode(&dma_get, PE_MODE);
                    dma_set_reply(&dma_get, &reply_get);
                    dma_set_size(&dma_get, bsize_epj * nn);
                    dma(dma_get,
                        (long*)((epjLM *)adr_epj_sorted_loc_ + adr_epj_loc_group_head[igl] + adr_offset),
                        (long*)((epjLM *)epj + copied_size));
                    dma_wait(&reply_get, 1);
                    while (reply_get != 1) {}
                 }
                 copied_size += nn; // the total number of copied data
                 if (nn < group_size) {
                     // In this case, the chunk is filled by a part of this group i.
                     // So, the remaining elements must be copied at the next copy.
                     // Therefore, group ID i shall NOT be updated.
                     break; 
                 }
                 igl++; // set the next group ID
                 if ((igl == ng) || (nn == nrem)) {
                     // In the former case, the number of groups is insufficient to
                     // fill the chunk. In the latter case, the chunk is filled 
                     // without splitting a group information.
                     break;
                 }
            }
           
            //** Put into epj_sorted_ 
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, bsize_epj * copied_size);
            dma(dma_put, (long*)((epjLM *)adr_epj_sorted_ + adr_chunk), (long*)(epj));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
 
            //** Setup the next chunk 
            chunk_id++;
            adr_chunk += copied_size;
            int nrem = (my_offset + n_loc) - adr_chunk;
            chunk_size = nrem < EPJ_ARRAY_SIZE ? nrem : EPJ_ARRAY_SIZE;
#ifdef CHECK_MSORT_GLB_TREE_REUSE
            if (myrank == RANK_CHECK && my_id == ID_CHECK) {
                printf("copied_data[A] = %d\n",copied_size);
                printf("chunk_size[N]  = %d\n",chunk_size);
            }
#endif

            //** Break the loop because all the group information were used.
            //   (we have to update group info.)
            if (igl == ng) break;
        }

#ifdef CHECK_MSORT_GLB_TREE_REUSE
        if (myrank == RANK_CHECK && my_id == ID_CHECK) {
            printf("the actual copy size = %d\n",adr_chunk-my_offset);
        }
#endif
    }

    //* Release memory
    ldm_free(epj, bsize_epj_array);
    ldm_free(adr_epj_glb_group_head, bsize_buf_array);
    ldm_free(adr_epj_loc_group_head, bsize_buf_array);
    ldm_free(group_size_epj_glb, bsize_buf_array);

#ifdef PROFILE_MSORT_GLB_TREE_REUSE
    //* Output profile
    sync_array_();
    unsigned long wt_end0 = rtc_();
    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
        double bw_dma_get=0.363e9, bw_dma_put=0.363e9; // bandwidth per CPE
        double cpu_clock=1.45e9;
        printf("Actual total time = %e [s]\n",(wt_end0-wt_start0)/cpu_clock);
    }
#endif
}
#else // VERSION_MSORT_GLB_TREE_REUSE
#error The value of the macro `VERSION_MSORT_GLB_TREE_REUSE` is incorrect.
#endif // VERSION_MSORT_GLB_TREE_REUSE

