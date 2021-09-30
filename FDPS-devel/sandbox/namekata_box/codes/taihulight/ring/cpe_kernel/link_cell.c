/* Standard C headers */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/* Headers specific to Sunway TaihuLight */
#include "slave.h"
#include "dma.h"
#include "cpe_func.h"
#include "cpe_prof.h"

/* Local macros */
//#define CHECK_LINK_CELL
#define RANK_CHECK (2)
#define ID_CHECK   (62)
// These macros are used to check the code.


//#define TIME_LINK_CELL
// This macro is used to time each part of the functions.


__thread_local int my_id, my_row_id, my_col_id;
__thread_local dma_desc dma_get, dma_put;
__thread_local volatile int reply_get = 0;
__thread_local volatile int reply_put = 0;

__thread_local volatile int debug_getPID = 0;

/* CPE communication functions */
static inline void cpe_bcast_int32(const int root_cpe_id, int* data) {
  int root_col_id, root_row_id;
  root_row_id = root_cpe_id/8;
  root_col_id = root_cpe_id%8;

  if (my_id == root_cpe_id) {
      REG_PUTR(*data,8);
      REG_PUTC(*data,8);
  }
  else {
    if (my_row_id == root_row_id) {
        REG_GETR(*data);
        REG_PUTC(*data,8);
    }
    else {
        REG_GETC(*data);
    }
  }
}

/* Local functions */
__attribute__((noinline))
static int getCellID(int lev, unsigned long long mkey){
    int kLevMax = 21;
    unsigned long long s = mkey >> ( (kLevMax - lev) * 3 );
    return (s & 0x7);
}

static int GetPartitionID(void *adr_tp,
                          int left,
                          int right,
                          int ref,
                          int lev) {
    //* Local variables
    tpLM tp_l,tp_r,tp_c;
    //* Set the parameters of DMA get
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(tpLM));
    //* Get tp[left] & tp[right]
    reply_get = 0;
    dma(dma_get, (long)((tpLM *)adr_tp +  left), (long)(&tp_l));
    dma(dma_get, (long)((tpLM *)adr_tp + right), (long)(&tp_r));
    dma_wait(&reply_get, 2);
    while (reply_get != 2) {}
#if defined(CHECK_LINK_CELL)
    if (debug_getPID == 1) {
        printf("tp_l.key_     = %lu\n",tp_l.key_);
        printf("tp_l.adr_ptcl = %d\n",tp_l.adr_ptcl_);
        printf("tp_r.key_     = %lu\n",tp_r.key_);
        printf("tp_r.adr_ptcl = %d\n",tp_r.adr_ptcl_);
        printf("getCID_r      = %d\n",getCellID(lev,tp_r.key_));
        printf("getCID_l      = %d\n",getCellID(lev,tp_l.key_));
    }
#endif
    if( getCellID(lev, tp_r.key_) < ref) return -1;
    if( getCellID(lev, tp_l.key_) > ref) return -1;
    //* Find with binary search
    int l = left;
    int r = right;
    while(l < r){
        //** Compute cen
        int cen = (l+r) / 2;
        //** Get tp[cen]
        reply_get = 0;
        dma(dma_get, (long)((tpLM *)adr_tp + cen), (long)(&tp_c));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
#if defined(CHECK_LINK_CELL)
        if (debug_getPID == 1) {
            printf("cen            = %d\n",cen);
            printf("tp_c.key_      = %lu\n",tp_c.key_);
            printf("tp_c.adr_ptcl_ = %d\n",tp_c.adr_ptcl_);
        }
#endif
        //** Compute val_cen
        int val_cen = getCellID(lev, tp_c.key_);
#if defined(CHECK_LINK_CELL)
        if (debug_getPID == 1) {
            printf("val_cen        = %d\n",val_cen);
        }
#endif
        if (ref > val_cen){
            l = cen+1;
        } else {
            r = cen;
        }
        //** Get tp[l]
        reply_get = 0;
        dma(dma_get, (long)((tpLM *)adr_tp + l), (long)(&tp_l));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
#if defined(CHECK_LINK_CELL)
        if (debug_getPID == 1) {
            printf("tp_l.key_      = %lu\n",tp_l.key_);
            printf("tp_l.adr_ptcl_ = %d\n",tp_l.adr_ptcl_);
            printf("l              = %d\n",l);
            printf("r              = %d\n",r);
        }
#endif
        //** Check 
        if( getCellID(lev, tp_l.key_) == ref ) return l;
    }
    return -1;
}
 
static int GetPartitionID_wo_DMA(tpLM tp[],
                                 int left,
                                 int right,
                                 int ref,
                                 int lev) {
    //* Local variables
    int l = 0,r = right - left;
    //* Get tp[left] & tp[right]
    if( getCellID(lev, tp[r].key_) < ref) return -1;
    if( getCellID(lev, tp[l].key_) > ref) return -1;
    //* Find with binary search
    while(l < r){
        //** Compute cen
        int cen = (l+r) / 2;
        //** Compute val_cen
        int val_cen = getCellID(lev, tp[cen].key_);
        if (ref > val_cen){
            l = cen+1;
        } else {
            r = cen;
        }
        //** Check 
        if( getCellID(lev, tp[l].key_) == ref ) return (l+left);
    }
    return -1;
} 


void LinkCell(void *args) {
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    //* Initialize local vars. for DMA comm.
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    //* Process the arguments
    int myrank                = (int)(((unsigned long*)args)[0]);
    void *adr_tc_array        = (void *)((unsigned long*)args)[1];
    void *adr_adr_tc_lev_part = (void *)((unsigned long*)args)[2];
    void *adr_tp              = (void *)((unsigned long*)args)[3];
    void *adr_lev_max         = (void *)((unsigned long*)args)[4];
    int n_tot                 = (int)(((unsigned long*)args)[5]);
    int n_leaf_limit          = (int)(((unsigned long*)args)[6]);
    int ca_tc_array           = (int)(((unsigned long*)args)[7]);
    void *adr_OOM_flag        = (void *)((unsigned long*)args)[8];
    void *adr_offset          = (void *)((unsigned long*)args)[9];
#if 0
    if ((myrank == 0) && (my_id == 0)) {
        printf("myrank              = %d\n",myrank);
        printf("adr_tc_array        = %lu\n",(unsigned long)adr_tc_array);
        printf("adr_adr_tc_lev_part = %lu\n",(unsigned long)adr_adr_tc_lev_part);
        printf("adr_tp              = %lu\n",(unsigned long)adr_tp);
        printf("adr_lev_max         = %lu\n",(unsigned long)adr_lev_max);
        printf("n_tot               = %d\n",n_tot);
        printf("n_leaf_limit        = %d\n",n_leaf_limit);
        printf("ca_tc_array         = %d\n",ca_tc_array);
        printf("adr_OOM_flag        = %lu\n",(unsigned long)adr_OOM_flag);
        printf("adr_offset          = %lu\n",(unsigned long)adr_offset);
    }
#endif
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
        TP_ARRAY_SIZE = 1024,
    };
    size_t bsize_tc = sizeof(tcLM);
    size_t bsize_tp = sizeof(tpLM);
    size_t bsize_tc_parent = bsize_tc * CHUNK_SIZE;
    size_t bsize_tc_child  = bsize_tc * N_CHILDREN;
    size_t bsize_tp_array  = bsize_tp * TP_ARRAY_SIZE;
    tcLM *tc_parent = (tcLM *)ldm_malloc( bsize_tc_parent );
    tcLM *tc_child  = (tcLM *)ldm_malloc( bsize_tc_child );
    tpLM *tp        = (tpLM *)ldm_malloc( bsize_tp_array );
    int *n_cell_new_tbl = (int *)ldm_malloc( sizeof(int) * NUMBER_OF_CPE);
    int *adr_tc_level_partition = (int *)ldm_malloc( sizeof(int) * (TREE_LEVEL_LIMIT+2) );
#if 0
    if ((myrank == 0) && (my_id == 0)) {
        printf("bsize_tc             = %d\n",bsize_tc);
        printf("bsize_tp             = %d\n",bsize_tp);
        printf("bsize_tc_parent      = %d\n",bsize_tc_parent);
        printf("bsize_tc_child       = %d\n",bsize_tc_child);
        printf("bsize_tc_array       = %d\n",bsize_tp_array);
        printf("bsize_tbl            = %d\n",sizeof(int)*NUMBER_OF_CPE);
        printf("bsize_adr_tc_lev_par = %d\n",sizeof(int)*(TREE_LEVEL_LIMIT+2));
        printf("Total (in bytes)     = %d\n",bsize_tc_parent + bsize_tc_child + bsize_tp_array
                                             + sizeof(int)*(NUMBER_OF_CPE+TREE_LEVEL_LIMIT+2) );
    }
#endif

    //* Local variables
    // (loop counters)
    int i,j,k,n,id;
    // (range of array index for tc_array[])
    int id_cell_left = 0, id_cell_right = N_CHILDREN-1; 
    int ista,iend;
    // (other vars.)
    int lev_max,n_cell_new,offset;
    int use_tp_buffer;
#if defined(TIME_LINK_CELL)
    // (time measurements)
    unsigned long tstart,t_blks[4];
    for (i=0; i<4; i++) t_blks[i]=0;
    unsigned long tstart_in, t_blks_in[4];
    for (i=0; i<4; i++) t_blks_in[i]=0;
#endif

    //* Get lev_max
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(int));
    dma(dma_get, (long)((int *)adr_lev_max), (long)(&lev_max));
    dma_wait(&reply_get, 1);
    while (reply_get != 1) {}
#if 0
    if ((myrank == 0) && (my_id == 0)) {
        printf("lev_max = %d\n",lev_max);
    }
#endif

    //** Compute the first task of each CPE
    int n_tc_tgt = id_cell_right - id_cell_left + 1;
    int n_loc,my_offset;
    n_loc = n_tc_tgt/NUMBER_OF_CPE + ( (my_id < n_tc_tgt % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = id_cell_left + (n_tc_tgt/NUMBER_OF_CPE)*my_id + ( (my_id < n_tc_tgt % NUMBER_OF_CPE) ? my_id : n_tc_tgt % NUMBER_OF_CPE );
#if 0
    if ((myrank == 0) && (my_id == 0)) {
        printf("lev_max       = %d\n",lev_max);
        printf("id_cell_left  = %d\n",id_cell_left);
        printf("id_cell_right = %d\n",id_cell_right);
        printf("my_id         = %d\n",my_id); 
        printf("n_loc         = %d\n",n_loc);
        printf("my_offset     = %d\n",my_offset);
    }
#endif

    //* Link tree cells
    //  assign particles to child cells and count # of particles in child cells
    //  but loop over parent cells because they have indexes of particles
    while(1){
#if defined(CHECK_LINK_CELL)
        sync_array_();
        if (myrank == RANK_CHECK && my_id == 0) {
            printf("Processing lev_max %d (limit: %d)\n",lev_max,TREE_LEVEL_LIMIT);
        }
        sync_array_();
#endif 
#if defined(TIME_LINK_CELL)
        tstart = rtc_();
#endif
        n_cell_new = 0;
        for (n=0; n<n_loc; n+= CHUNK_SIZE) {
            int nrem = n_loc - n;
            int nn = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
            ista = my_offset + n;
            iend = ista + nn - 1;
            // Check
#if defined(CHECK_LINK_CELL)
            if (myrank == RANK_CHECK && my_id == ID_CHECK) {
                printf("###### Outermost loop: n = %d (myrank,my_id,lev_max) = (%d,%d,%d)\n",n,myrank,my_id,lev_max);
                printf("ista=%d, iend=%d, nn = %d\n",ista,iend,nn);
            }
#endif
            //** Get tc_array[i ~ i+nn-1] and store them into tc_parent[]
#if defined(TIME_LINK_CELL)
            tstart_in = rtc_();
#endif
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tc * nn);
            dma(dma_get, (long)((tcLM *)adr_tc_array + ista), (long)(tc_parent));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
#if defined(TIME_LINK_CELL)
            t_blks_in[0] += rtc_() - tstart_in;
#endif
            //** Process each of tc_parent[]
            for (i=ista; i<=iend; i++) {
                int ii = i - ista;
                int n_ptcl_tmp = tc_parent[ii].n_ptcl_;
                if (n_ptcl_tmp <= n_leaf_limit) continue;
                int adr_ptcl_tmp = tc_parent[ii].adr_ptcl_;
                // Check
#if defined(CHECK_LINK_CELL)
                if (myrank == RANK_CHECK && my_id == ID_CHECK) {
                    printf("====== Middle loop: i = %d\n",i);
                    printf("n_ptcl_tmp   = %d\n",n_ptcl_tmp);
                    printf("adr_ptcl_tmp = %d\n",adr_ptcl_tmp);
                    printf("left         = %d\n",adr_ptcl_tmp);
                    printf("right        = %d\n",adr_ptcl_tmp+n_ptcl_tmp-1);
                }
#endif

                // Get tp if a search range is small
                if (n_ptcl_tmp < TP_ARRAY_SIZE) {
                    use_tp_buffer = 1;
                    reply_get = 0;
                    dma_set_op(&dma_get, DMA_GET);
                    dma_set_mode(&dma_get, PE_MODE);
                    dma_set_reply(&dma_get, &reply_get);
                    dma_set_size(&dma_get, bsize_tp * n_ptcl_tmp);
                    dma(dma_get, (long)((tpLM *)adr_tp + adr_ptcl_tmp), (long)(tp));
                    dma_wait(&reply_get, 1);
                    while (reply_get != 1) {}
                } else {
                   use_tp_buffer = 0;
                }

                // Link child tree cells
#if defined(TIME_LINK_CELL)
                tstart_in = rtc_();
#endif
                int adr[N_CHILDREN];
                int n_cnt = 0;
                for(j=0; j<N_CHILDREN; j++){
                    int adr_tc_tmp = tc_parent[ii].adr_tc_ + j;
                    // Check
#if defined(CHECK_LINK_CELL)
                    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
                        printf("---- Inner loop: j = %d\n",j);
                        printf("adr_tc_tmp    = %d\n",adr_tc_tmp);
                        printf("use_tp_buffer = %d\n",use_tp_buffer);
                    }
                    if (n == 192 && i == 21939 && j == 7) debug_getPID = 1;
                    else debug_getPID = 0;
#endif
                    if (use_tp_buffer == 0) {
                        tc_child[j].adr_ptcl_ = GetPartitionID(adr_tp, adr_ptcl_tmp, adr_ptcl_tmp+n_ptcl_tmp-1, j, lev_max+1);
                    } else {
                        tc_child[j].adr_ptcl_ = GetPartitionID_wo_DMA(tp, adr_ptcl_tmp, adr_ptcl_tmp+n_ptcl_tmp-1, j, lev_max+1);
                    }

                    // Check
#if defined(CHECK_LINK_CELL)
                    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
                        printf("adr_ptcl_    = %d\n",tc_child[j].adr_ptcl_);
                    }
#endif
                    tc_child[j].level_ = lev_max+1;
                    if((tc_child[j].adr_ptcl_>>31) & 0x1){
                        tc_child[j].n_ptcl_ = 0;
                        tc_child[j].adr_tc_ = ((U32_)0) | 0x80000000;
                    } else {
                        adr[n_cnt] = j;
                        n_cnt++;
                    }
#if defined(CHECK_LINK_CELL)
                    if (myrank == RANK_CHECK && my_id == ID_CHECK) {
                        printf("level_       = %d\n",tc_child[j].level_);
                        printf("n_cnt        = %d\n",n_cnt);
                    }
#endif
                }
                int n_ptcl_cum = 0;
                for (j=0; j<n_cnt-1; j++){
                    tc_child[adr[j]].n_ptcl_ = tc_child[adr[j+1]].adr_ptcl_ - tc_child[adr[j]].adr_ptcl_;
                    n_ptcl_cum += tc_child[adr[j]].n_ptcl_;
                }
                tc_child[adr[n_cnt-1]].n_ptcl_ = n_ptcl_tmp - n_ptcl_cum;
#if defined(TIME_LINK_CELL)
                t_blks_in[2] += rtc_() - tstart_in;
#endif
#if defined(CHECK_LINK_CELL)
                if (myrank == RANK_CHECK && my_id == ID_CHECK) {
                    printf("List of n_ptcl_:\n");
                    for (j=0; j<n_cnt; j++)
                        printf("  j = %d, n_ptcl_ = %d\n",adr[j],tc_child[adr[j]].n_ptcl_);
                }
#endif
                // Put tc_child[] to tc_array[tc_parent[ii].adr_tc_]
#if defined(TIME_LINK_CELL)
                tstart_in = rtc_();
#endif
                reply_put = 0;
                dma_set_op(&dma_put, DMA_PUT);
                dma_set_mode(&dma_put, PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put, bsize_tc * N_CHILDREN);
                dma(dma_put, (long)((tcLM *)adr_tc_array + tc_parent[ii].adr_tc_), (long)(tc_child));
                dma_wait(&reply_put, 1);
                while (reply_put != 1) {}
#if defined(TIME_LINK_CELL)
                t_blks_in[3] += rtc_() - tstart_in;
#endif
                // Update n_cell_new
                n_cell_new += N_CHILDREN;
            } // END: i-loop
        } // END: n-loop (chunk-loop)
        sync_array_();
#if defined(TIME_LINK_CELL)
        t_blks[0] += rtc_() - tstart;
#endif

        //** Allreduce n_cell_new
#if defined(TIME_LINK_CELL)
        tstart = rtc_();
#endif
        n_cell_new_tbl[my_id] = n_cell_new;
        n_cell_new = 0; // reset to 0 in order to store the sum. 
        for (id=0; id<NUMBER_OF_CPE; id++) {
            cpe_bcast_int32(id,&n_cell_new_tbl[id]);
            n_cell_new += n_cell_new_tbl[id];
        }
#if defined(TIME_LINK_CELL)
        t_blks[1] += rtc_() - tstart;
#endif
        if (n_cell_new == 0) break; // Escape from while(1) loop

        //** Update id_* to go deeper
        id_cell_left = id_cell_right + 1;
        id_cell_right += n_cell_new;
        lev_max++;
        adr_tc_level_partition[lev_max+1] = id_cell_right + 1;
        if (lev_max == TREE_LEVEL_LIMIT) break;
       
        //** Compute the next task of each CPE
        n_tc_tgt = id_cell_right - id_cell_left + 1;
        n_loc = n_tc_tgt/NUMBER_OF_CPE + ( (my_id < n_tc_tgt % NUMBER_OF_CPE) ? 1 : 0 );
        my_offset = id_cell_left + (n_tc_tgt/NUMBER_OF_CPE)*my_id + ( (my_id < n_tc_tgt % NUMBER_OF_CPE) ? my_id : n_tc_tgt % NUMBER_OF_CPE );

        //** Compute prefix sum
#if defined(TIME_LINK_CELL)
        tstart = rtc_();
#endif
        n_cell_new = 0;
        for (n=0; n<n_loc; n+= CHUNK_SIZE) {
            int nrem = n_loc - n;
            int nn = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
            ista = my_offset + n;
            iend = ista + nn - 1;
            //** Get tc_array[i ~ i+nn-1] and store them into tc_parent[]
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tc * nn);
            dma(dma_get, (long)((tcLM *)adr_tc_array + ista), (long)(tc_parent));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            //** Compute n_cell_new
            for (i=ista; i<=iend; i++){
                int ii = i - ista;
                if(!(tc_parent[ii].n_ptcl_ <= n_leaf_limit || 
                     tc_parent[ii].level_ == TREE_LEVEL_LIMIT)) {
                    n_cell_new += N_CHILDREN;
                }
            }
        }
        n_cell_new_tbl[my_id] = n_cell_new;
        offset = id_cell_right+1;
        for (id=0; id<NUMBER_OF_CPE; id++) {
            cpe_bcast_int32(id,&n_cell_new_tbl[id]);
            if (id < my_id) offset += n_cell_new_tbl[id];
        }
#if defined(TIME_LINK_CELL)
        t_blks[2] += rtc_() - tstart;
#endif
#if 0
        if ((myrank == 0) && (my_id == 1)) {
            printf("offset = %d (my_id = %d)\n",offset,my_id);
        }
#endif

        //** Set tc_array[].adr_tc_
#if defined(TIME_LINK_CELL)
        tstart = rtc_();
#endif
        for (n=0; n<n_loc; n+= CHUNK_SIZE) {
            int nrem = n_loc - n;
            int nn = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
            ista = my_offset + n;
            iend = ista + nn - 1;
            //** Get tc_array[i ~ i+nn-1] and store them into tc_parent[]
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tc * nn);
            dma(dma_get, (long)((tcLM *)adr_tc_array + ista), (long)(tc_parent));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            //** Compute tc_array[].adr_tc_
            for (i=ista; i<=iend; i++){
                int ii = i - ista;
                if(!(tc_parent[ii].n_ptcl_ <= n_leaf_limit || 
                     tc_parent[ii].level_ == TREE_LEVEL_LIMIT)) {
                    tc_parent[ii].adr_tc_ = offset;
                    offset += N_CHILDREN;
                }
            }
            // Put tc_array[] 
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, bsize_tc * nn);
            dma(dma_put, (long)((tcLM *)adr_tc_array + ista), (long)(tc_parent));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }
#if defined(TIME_LINK_CELL)
        t_blks[3] += rtc_() - tstart;
#endif

        // Check if there are sufficient memory 
        cpe_bcast_int32(63,&offset); // we must use `offset` of the last CPE.
        sync_array_();
        if (offset > ca_tc_array) {
           if (my_id == 0) {
                // Put OOM_flag
                int OOM_flag = 1;
                reply_put = 0;
                dma_set_op(&dma_put, DMA_PUT);
                dma_set_mode(&dma_put, PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put, sizeof(int));
                dma(dma_put, (long*)(adr_OOM_flag), (long*)(&OOM_flag));
                dma_wait(&reply_put, 1);
                while (reply_put !=1 ) {}
            }
            break;
        } 
    }

    //* Put lev_max, adr_tc_level_partition, and offset
    if (my_id == 0) {
        // (i) lev_max
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, sizeof(int));
        dma(dma_put, (long)(adr_lev_max), (long)(&lev_max));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
        // (ii) adr_tc_level_partition
        if (lev_max > 0) {
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, sizeof(int) * lev_max );
            dma(dma_put, (long)((int *)adr_adr_tc_lev_part + 2), (long)(adr_tc_level_partition + 2));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
            // Note that adr_tc_level_partition[0,1] are already set in MPE side.
            // Therefore, we only have to update adr_tc_level_partition[i] (i>1).
        }
        // (iii) offset
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, sizeof(int));
        dma(dma_put, (long)(adr_offset), (long)(&offset));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    
#if defined(TIME_LINK_CELL)
    if (myrank == 0 && my_id == 0) {
        printf("---- time profile for slave_LinkCell() ----\n");
        printf("1st block = %lu\n",t_blks[0]);
        printf("   get  tc_parent = %lu\n",t_blks_in[0]);
        printf("   get   tc_child = %lu\n",t_blks_in[1]);
        printf("   proc. tc_child = %lu\n",t_blks_in[2]);
        printf("   put   tc_child = %lu\n",t_blks_in[3]);
        printf("alltoall  = %lu\n",t_blks[1]);
        printf("2nd block = %lu\n",t_blks[2]);
        printf("3rd block = %lu\n",t_blks[3]);
    }
#endif

    //* Release memory 
    ldm_free(tc_parent, bsize_tc_parent );
    ldm_free(tc_child,  bsize_tc_child );
    ldm_free(tp, bsize_tp_array );
    ldm_free(n_cell_new_tbl, sizeof(int) * NUMBER_OF_CPE);
    ldm_free(adr_tc_level_partition, sizeof(int) * (TREE_LEVEL_LIMIT+2) );

}
