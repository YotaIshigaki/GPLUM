/* Standard C headers */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/* Headers specific to Sunway TaihuLight */
#include "slave.h"
#include "dma.h"
#include "cpe_func.h"

/* Local macros */
#define DEBUG_MAKELIST

__thread_local int my_id, my_row_id, my_col_id;





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

void MakeListUsingTreeWithStack(void * adr_tc_first,
                                int adr_tc,
                                void * adr_tp_first,
                                int * size_of_ep_list,
                                int * id_ep_list,
                                int * size_of_sp_list,
                                int * id_sp_list,
                                F64ort_ pos_target_box,
                                double r_crit_sq,
                                double len_peri_x,
                                int n_leaf_limit,
                                int adr_tree_sp_first // adress of first sp coming from the (global) tree.
#ifdef REMOVE_VERTEX
                                ,F64_   r_search,
                                F64_    hlen_root,
                                F64vec_ cen_root
#endif
                                ){
#ifdef REMOVE_VERTEX
    F64vec_ TREE_SHIFT_CENTER[8] = { {-0.5, -0.5, -0.5}, {-0.5, -0.5, 0.5},
                                     {-0.5,  0.5, -0.5}, {-0.5,  0.5, 0.5},
                                     { 0.5, -0.5, -0.5}, { 0.5, -0.5, 0.5},
                                     { 0.5,  0.5, -0.5}, { 0.5,  0.5, 0.5} };
#endif
    //* Prepare stacks
    enum {
       STACK_SIZE = 128,
    };
    enum { // initial value for open_bit
       OPEN_BIT_IV = 1<<(N_CHILDREN+1),
    };
    typedef struct {
        S32_ adr_tc;
        F64_ r_crit_sq;
        U32_ open_bit;
        S32_ ista;
    #ifdef REMOVE_VERTEX
        F64_    hlen_tree;
        F64vec_ cen_tree;
    #endif
    } Stack_t;
    Stack_t stack[STACK_SIZE];
    S32_ stack_num = 0;
    //* Initialize the stack
    stack_num = 1;
    stack[0].adr_tc = adr_tc;
    stack[0].r_crit_sq = r_crit_sq;
    stack[0].open_bit = OPEN_BIT_IV;
    stack[0].ista = 0;
#ifdef REMOVE_VERTEX
    stack[0].hlen_tree = hlen_root;
    stack[0].cen_tree  = cen_root;
#endif // REMOVE_VERTEX
    //* Memory allocation of local buffers
    enum {
       TP_ARRAY_SIZE = 64,
    };
    size_t bsize_tc = sizeof(tcLM);
    size_t bsize_tp = sizeof(tpLM);
    size_t bsize_tc_child = bsize_tc * N_CHILDREN;
    size_t bsize_tp_array = bsize_tp * TP_ARRAY_SIZE; 
    volatile tcLM *tc_cur   = (tcLM *) ldm_malloc( bsize_tc );
    volatile tcLM *tc_child = (tcLM *) ldm_malloc( bsize_tc_child ); 
    volatile tpLM *tp       = (tpLM *) ldm_malloc( bsize_tp_array );

    //* Make an interaction list with stack
    // (loop counters)
    S32_ i,j,k,ip;
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    // (local vars for open_bit)
    F64_ dx,dy,dz,dr2,dis0,dis1;
    U32_ f_contained;
    while(stack_num) {
        //* Pop a tc from the stack
        S32_ adr_tc_tgt = stack[stack_num-1].adr_tc;
        F64_ r_crit_sq_cur = stack[stack_num-1].r_crit_sq;
        U32_ open_bit = stack[stack_num-1].open_bit;
        S32_ ista = stack[stack_num-1].ista;
#ifdef REMOVE_VERTEX
        F64_    hlen1 = stack[stack_num-1].hlen_tree;
        F64vec_ cen1  = stack[stack_num-1].cen_tree;
#endif // REMOVE_VERTEX
        stack_num--;
        //* Get tc_cur 
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tc);
        dma(dma_get, (long*)((tcLM *)adr_tc_first + adr_tc_tgt), (long*)(tc_cur));
        dma_wait(&reply_get, 1);
        while (reply_get !=1 ) {}
        //* Initialize local aux. vars.
        F32vec_ pos_tc = tc_cur->pos; // not used actually
        S32_ n_ptcl = tc_cur->n_ptcl_;
        S32_ adr_ptcl_child = tc_cur->adr_ptcl_;
        S32_ adr_tc_child = tc_cur->adr_tc_;
        if( !(tc_cur->n_ptcl_ <= n_leaf_limit || tc_cur->level_ == TREE_LEVEL_LIMIT) ){ // not leaf
            //* Get child cells
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tc_child);
            dma(dma_get, (long*)((tcLM *)adr_tc_first + adr_tc_child), (long*)(tc_child));
            dma_wait(&reply_get, 1);
            while (reply_get !=1 ) {}
            //* Compute open_bit
            if (open_bit == OPEN_BIT_IV) {
                open_bit = 0;
                for(i=0; i<N_CHILDREN; i++){
#ifdef PHI_R_TREE
                    // getPosForTree()
                    F32vec_ pos;
                    pos.x = tc_child[i].pos_phi;
                    pos.y = tc_child[i].pos_r;
                    pos.z = tc_child[i].pos.z;
                    //** Compute GetDistanceMinSqPeriodicX(pos_target_box, pos, len_peri_x);
                    F64vec_ cen0, len;
                    cen0.x = 0.5*(pos_target_box.high_.x + pos_target_box.low_.x);
                    cen0.y = 0.5*(pos_target_box.high_.y + pos_target_box.low_.y);
                    cen0.z = 0.5*(pos_target_box.high_.z + pos_target_box.low_.z);
                    len.x  = 0.5*(pos_target_box.high_.x - pos_target_box.low_.x);
                    len.y  = 0.5*(pos_target_box.high_.y - pos_target_box.low_.y);
                    len.z  = 0.5*(pos_target_box.high_.z - pos_target_box.low_.z);
                    dx = fabs(cen0.x - pos.x);
                    dx = (dx < (len_peri_x - dx)) ? dx : (len_peri_x - dx);
                    dx = (len.x < dx) ? (dx - len.x) : 0.0;
                    dy = (len.y < fabs(cen0.y - pos.y)) ? fabs(cen0.y - pos.y)-len.y : 0.0;
                    dz = (len.z < fabs(cen0.z - pos.z)) ? fabs(cen0.z - pos.z)-len.z : 0.0;
                    dis0 = dx*dx + dy*dy + dz*dz;
                    //** Compute GetDistanceMinSqPeriodicX(pos_target_box, tc_first[adr_tc+i].mom_.getVertexOut(), len_peri_x);
    #ifdef REMOVE_VERTEX
                    F64_ len_tmp = hlen1 + 2.0*r_search;
                    dx = fabs(cen0.x - cen1.x);
                    dx = (dx < (len_peri_x - dx)) ? dx : (len_peri_x - dx);
                    dx = (len_tmp < dx) ? (dx-len_tmp) : 0.0;
                    dy = (len_tmp < fabs(cen0.y - cen1.y)) ? fabs(cen0.y - cen1.y)-len_tmp : 0.0;
                    dz = (len_tmp < fabs(cen0.z - cen1.z)) ? fabs(cen0.z - cen1.z)-len_tmp : 0.0;
    #else
                    F64vec_ cen1;
                    cen1.x = 0.5*(tc_child[i].vertex_out_.high_.x + tc_child[i].vertex_out_.low_.x);
                    cen1.y = 0.5*(tc_child[i].vertex_out_.high_.y + tc_child[i].vertex_out_.low_.y);
                    cen1.z = 0.5*(tc_child[i].vertex_out_.high_.z + tc_child[i].vertex_out_.low_.z);
                    len.x  += 0.5*(tc_child[i].vertex_out_.high_.x - tc_child[i].vertex_out_.low_.x);
                    len.y  += 0.5*(tc_child[i].vertex_out_.high_.y - tc_child[i].vertex_out_.low_.y);
                    len.z  += 0.5*(tc_child[i].vertex_out_.high_.z - tc_child[i].vertex_out_.low_.z);
                    // Note that the calculation of `len` corresponds to the original implementation in ps_defs.hpp:
                    // const F64vec len = pos0.getHalfLength() + pos1.getHalfLength();
                    // Here, pos0.getHalfLength() is already stored in `len` before the start of the calculation.
                    // (see the calculation of `dis0`)
                    dx = fabs(cen0.x - cen1.x);
                    dx = (dx < (len_peri_x - dx)) ? dx : (len_peri_x - dx);
                    dx = (len.x < dx) ? (dx-len.x) : 0.0;
                    dy = (len.y < fabs(cen0.y - cen1.y)) ? fabs(cen0.y - cen1.y)-len.y : 0.0;
                    dz = (len.z < fabs(cen0.z - cen1.z)) ? fabs(cen0.z - cen1.z)-len.z : 0.0;
    #endif
                    dis1 = dx*dx + dy*dy + dz*dz;
                    //** Compute open_bit
                    open_bit |= ( (dis0 <= r_crit_sq_cur*0.25) || (dis1 <= 0.0) ) << i;
#else
                    F32vec_ pos = tc_child[i].pos;
                    //** Compute pos_target_box.getDistanceMinSQ(pos)
                    dx = (pos.x > pos_target_box.high_.x) ? (pos.x - pos_target_box.high_.x) : ( (pos.x < pos_target_box.low_.x) ? (pos_target_box.low_.x - pos.x) : 0.0 );
                    dy = (pos.y > pos_target_box.high_.y) ? (pos.y - pos_target_box.high_.y) : ( (pos.y < pos_target_box.low_.y) ? (pos_target_box.low_.y - pos.y) : 0.0 );
                    dz = (pos.z > pos_target_box.high_.z) ? (pos.z - pos_target_box.high_.z) : ( (pos.z < pos_target_box.low_.z) ? (pos_target_box.low_.z - pos.z) : 0.0 );
                    dr2 = dx*dx + dy*dy + dz*dz;
                    //** Compute pos_target_box.contained(tc_child[i].vertex_out_)
                    f_contained = (tc_child[i].vertex_out_.high_.x < pos_target_box.low_.x) || 
                                  (pos_target_box.high_.x < tc_child[i].vertex_out_.low_.x) ||
                                  (tc_child[i].vertex_out_.high_.y < pos_target_box.low_.y) ||
                                  (pos_target_box.high_.y < tc_child[i].vertex_out_.low_.y) ||
                                  (tc_child[i].vertex_out_.high_.z < pos_target_box.low_.z) ||
                                  (pos_target_box.high_.z < tc_child[i].vertex_out_.low_.z);
                    f_contained = f_contained ^ 0x1;
                    //** Compute open_bit
                    open_bit |= ( (dr2 <= r_crit_sq_cur*0.25) || f_contained ) << i;
#endif
                }
            }
            for(i=ista; i<N_CHILDREN; i++){
                if( tc_child[i].n_ptcl_ <= 0) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
                    //MakeListUsingTreeWithStack<TSM, TMLM, TCM, Ttc, Ttp, Tep, Tsp>
                    //    (tc_first, adr_tc_child+i, tp_first,
                    //     ep_first, ep_list, id_ep_list,
                    //     sp_first, sp_list, id_sp_list,
                    //     pos_target_box, r_crit_sq*0.25, n_leaf_limit, adr_tree_sp_first);
                    //* Save the current state to the stacks
                    stack_num++;
                    stack[stack_num-1].adr_tc = adr_tc_tgt;
                    stack[stack_num-1].r_crit_sq = r_crit_sq_cur;
                    stack[stack_num-1].open_bit = open_bit;
                    stack[stack_num-1].ista = i + 1;
#ifdef REMOVE_VERTEX
                    stack[stack_num-1].hlen_tree = hlen1;
                    stack[stack_num-1].cen_tree  =  cen1;
#endif // REMOVE_VERTEX
                    
                    //* Set the next target tc
                    stack_num++;
                    stack[stack_num-1].adr_tc = adr_tc_child + i;
                    stack[stack_num-1].r_crit_sq = r_crit_sq_cur * 0.25;
                    stack[stack_num-1].open_bit = OPEN_BIT_IV;
                    stack[stack_num-1].ista = 0;
#ifdef REMOVE_VERTEX
                    stack[stack_num-1].hlen_tree = hlen1 * 0.5;
                    stack[stack_num-1].cen_tree.x  = cen1.x + TREE_SHIFT_CENTER[i].x*hlen1;
                    stack[stack_num-1].cen_tree.y  = cen1.y + TREE_SHIFT_CENTER[i].y*hlen1;
                    stack[stack_num-1].cen_tree.z  = cen1.z + TREE_SHIFT_CENTER[i].z*hlen1;
#endif // REMOVE_VERTEX
                    
#if defined(DEBUG_MAKELIST)
                    //* Check
                    assert(stack_num < STACK_SIZE);
#endif
                    break;
                }
                else{ // far
                    //CopyInfoDistant(typename TSM::force_type(),
                    //                typename TCM::list_content_type(),
                    //                tc_first,
                    //                adr_tc_child+i,
                    //                adr_tc_child+i+adr_tree_sp_first,
                    //                sp_list,
                    //                id_sp_list);
                    S32_ adr_sp = adr_tc_child + i + adr_tree_sp_first;
                    id_sp_list[*size_of_sp_list] = adr_sp;
                    (*size_of_sp_list)++;
                }
            }
        }
        else{ //* leaf
            //CopyInfoClose(typename TSM::force_type(),
            //              typename TMLM::list_mode_type(),
            //              typename TCM::list_content_type(),
            //              tp_first,
            //              adr_ptcl_child,
            //              n_ptcl,
            //              ep_first,
            //              ep_list,
            //              id_ep_list,
            //              sp_first,
            //              sp_list,
            //              id_sp_list);
            S32_ cnt_adr_ptcl = adr_ptcl_child;
            //** Get tp
#if defined(DEBUG_MAKELIST)
            assert(n_ptcl < TP_ARRAY_SIZE); 
#endif
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tp * n_ptcl);
            dma(dma_get, (long*)((tpLM *)adr_tp_first + cnt_adr_ptcl), (long*)(tp));
            dma_wait(&reply_get, 1);
            while (reply_get !=1 ) {}
            //** Update interaction list
            for(ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
#ifdef REDUCE_MEMORY
                U32_ adr = tp[ip].adr_ptcl_;
                if( (adr>>31 & 0x1) == 0){
                    id_ep_list[*size_of_ep_list] = adr;
                    (*size_of_ep_list)++;
                }
                else{
                    id_sp_list[*size_of_sp_list] = (adr&0x7fffffff);
                    (*size_of_sp_list)++;
                }
#else
                if( ((tp[ip].adr_ptcl_>>31) & 0x1) == 0){
                    id_ep_list[*size_of_ep_list] = cnt_adr_ptcl;
                    (*size_of_ep_list)++;
                }
                else{
                    id_sp_list[*size_of_sp_list] = cnt_adr_ptcl;
                    (*size_of_sp_list)++;
                }
#endif
            }
        }
    }
    //* Release memory
    ldm_free(tc_cur,   bsize_tc );
    ldm_free(tc_child, bsize_tc_child ); 
    ldm_free(tp,       bsize_tp_array );
}

void MakeListUsingTree(void *args) {
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    //* Initialize local vars. for DMA comm.
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    //* Process the arguments
    int myrank                = (int)(((unsigned long*)args)[0]);
    int n_ipg                 = (int)(((unsigned long*)args)[1]);
    void *adr_ipg_            = (void *)((unsigned long*)args)[2];
    void *adr_tc_glb_         = (void *)((unsigned long*)args)[3];
    void *adr_tp_glb_         = (void *)((unsigned long*)args)[4];
    void *adr_id_epj          = (void *)((unsigned long*)args)[5];
    void *adr_id_spj          = (void *)((unsigned long*)args)[6];
    void *adr_n_disp_epj      = (void *)((unsigned long*)args)[7];
    void *adr_n_disp_spj      = (void *)((unsigned long*)args)[8];
    void *adr_n_epj           = (void *)((unsigned long*)args)[9];
    void *adr_n_spj           = (void *)((unsigned long*)args)[10];
    void *adr_r_crit_sq       = (void *)((unsigned long*)args)[11];
    int n_leaf_limit_         = (int)(((unsigned long*)args)[12]);
    int sp_tree_offset        = (int)(((unsigned long*)args)[13]);
#if PHI_R_TREE
    void *adr_len_peri_x      = (void *)((unsigned long*)args)[14];
#endif
    int ca_id_epj             = (int)(((unsigned long*)args)[15]);
    int ca_id_spj             = (int)(((unsigned long*)args)[16]);
    void *adr_n_epj_shortfall = (void *)((unsigned long*)args)[17];
    void *adr_n_spj_shortfall = (void *)((unsigned long*)args)[18];
#ifdef REMOVE_VERTEX
    void *adr_r_search = (void *)((unsigned long*)args)[19];
    void *adr_half_len_tree = (void *)((unsigned long*)args)[20];
    void *adr_cen_tree = (void *)((unsigned long*)args)[21];
#endif
    
#if 0
    //* Data Check
    if ((myrank == 0) && (my_id == 63)) {
        printf("n_ipg            = %d\n",n_ipg);
        printf("adr_ipg          = %lu\n",((unsigned long *)adr_ipg_));
        printf("adr_tc_glb_      = %lu\n",((unsigned long *)adr_tc_glb_));
        printf("adr_tp_glb_      = %lu\n",((unsigned long *)adr_tp_glb_));
        printf("adr_id_epj       = %lu\n",((unsigned long *)adr_id_epj));
        printf("adr_id_spj       = %lu\n",((unsigned long *)adr_id_spj));
        printf("adr_n_disp_epj   = %lu\n",((unsigned long *)adr_n_disp_epj));
        printf("adr_n_disp_spj   = %lu\n",((unsigned long *)adr_n_disp_spj));
        printf("adr_n_epj        = %lu\n",((unsigned long *)adr_n_epj));
        printf("adr_n_spj        = %lu\n",((unsigned long *)adr_n_spj));
        printf("adr_r_crit_sq    = %lu\n",((unsigned long *)adr_r_crit_sq));
        printf("n_leaf_limit_    = %d\n",n_leaf_limit_);
        printf("sp_tree_offset   = %d\n",sp_tree_offset);
        printf("adr_n_ipg_locals = %lu\n",((unsigned long *)adr_n_ipg_locals));
    }
#endif
    //* Memory allocation of local buffers
    enum {
        //N_EPJ_PER_GROUP_MAX = 10240,
        //N_SPJ_PER_GROUP_MAX = 3072,
        N_EPJ_PER_GROUP_MAX = 8000,
        N_SPJ_PER_GROUP_MAX = 5000,
    };
    size_t bsize_id_epj_array = sizeof(int) * N_EPJ_PER_GROUP_MAX;
    size_t bsize_id_spj_array = sizeof(int) * N_SPJ_PER_GROUP_MAX;
    size_t bsize_n_disp_epj = sizeof(int);
    size_t bsize_n_disp_spj = sizeof(int);
    size_t bsize_n_epj = sizeof(int);
    size_t bsize_n_spj = sizeof(int);
    size_t bsize_ipg   = sizeof(ipgLM);
    volatile int *id_epj_cpe_local = (int *) ldm_malloc( bsize_id_epj_array );
    volatile int *id_spj_cpe_local = (int *) ldm_malloc( bsize_id_spj_array );
    volatile int *n_disp_epj_cpe_local = (int *) ldm_malloc( bsize_n_disp_epj );
    volatile int *n_disp_spj_cpe_local = (int *) ldm_malloc( bsize_n_disp_spj );
    volatile int *n_epj_cpe_local = (int *) ldm_malloc( bsize_n_epj );
    volatile int *n_spj_cpe_local = (int *) ldm_malloc( bsize_n_spj );
    volatile int *n_epj_tbl = (int *) ldm_malloc( bsize_n_epj * NUMBER_OF_CPE);
    volatile int *n_spj_tbl = (int *) ldm_malloc( bsize_n_spj * NUMBER_OF_CPE);
    volatile ipgLM *ipg_ = (ipgLM *) ldm_malloc( bsize_ipg );
    volatile double *r_crit_sq = (double *) ldm_malloc( sizeof(double) );
    volatile double *len_peri_x = (double *) ldm_malloc( sizeof(double) );
    #ifdef REMOVE_VERTEX
    volatile double  *r_search = (double *) ldm_malloc( sizeof(double) );
    volatile double  *half_len_tree = (double *) ldm_malloc( sizeof(double) );
    volatile F64vec_ *cen_tree      = (double *) ldm_malloc( sizeof(double)*3 );
    #endif
    
    //* Main part of making interaction lists
    // (loop counters)
    int i,j,k,ig,id; 
    // (local vars. for making interaction list)
    int n_epj_local_1task,n_spj_local_1task;
    int n_epj_total_1task,n_spj_total_1task;
    int n_epj_total,n_spj_total;
    int n_disp_epj_local,n_disp_spj_local;
    int n_epj_shortfall,n_spj_shortfall;
    F64ort_ pos_target_box;
    //* Get r_crit_sq
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(double));
    dma(dma_get, (long*)((double *)adr_r_crit_sq), (long*)(r_crit_sq));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
    //* Get len_peri_x
#if PHI_R_TREE
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(double));
    dma(dma_get, (long*)((double *)adr_len_peri_x), (long*)(len_peri_x));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
#else
    // In this case, we set a dummy value to it.
    *len_peri_x = 0.0;
#endif

#ifdef REMOVE_VERTEX
    dma_descriptor_init(&dma_get, (unsigned *)&reply_get);
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(double));
    dma(dma_get, (long*)((double *)adr_r_search), (long*)(r_search));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
    
    dma_descriptor_init(&dma_get, (unsigned *)&reply_get);
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(double));
    dma(dma_get, (long*)((double *)adr_half_len_tree), (long*)(half_len_tree));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
    
    dma_descriptor_init(&dma_get, (unsigned *)&reply_get);
    reply_get = 0;
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, &reply_get);
    dma_set_size(&dma_get, sizeof(double)*3);
    dma(dma_get, (long*)((double *)adr_cen_tree), (long*)(cen_tree));
    dma_wait(&reply_get, 1);
    while (reply_get !=1 ) {}
#endif
    
#if 0
    if ((myrank == 0) && (my_id == 63)) {
        printf("r_crit_sq  = %e\n",*r_crit_sq);
#if PHI_R_TREE
        printf("len_peri_x = %e\n",*len_peri_x);
#endif
    }
#endif
    //* Calculate interaction list
    n_epj_total = 0;
    n_spj_total = 0;
    for(i=0; i<n_ipg; i+=NUMBER_OF_CPE){
        //** Set group id for this CPE
        ig = i + my_id;
        if (ig >= n_ipg) ig = -1;
        //** Reset list sizes
        n_epj_local_1task = 0;
        n_spj_local_1task = 0;
        if (ig != -1) {
            //** Get ipg_[ig]
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_ipg);
            dma(dma_get, (long*)((ipgLM *)adr_ipg_ + ig), (long*)(ipg_));
            dma_wait(&reply_get, 1);
            while (reply_get !=1 ) {}
            //** Set pos_target_box
            pos_target_box = (*ipg_).vertex_;
            //** Make interaction list for group `ig`

    #ifdef REMOVE_VERTEX
            MakeListUsingTreeWithStack(adr_tc_glb_, 0, adr_tp_glb_, 
                                       &n_epj_local_1task, id_epj_cpe_local, 
                                       &n_spj_local_1task, id_spj_cpe_local, 
                                       pos_target_box, *r_crit_sq, *len_peri_x,
                                       n_leaf_limit_, sp_tree_offset,
                                       *r_search, *half_len_tree, *cen_tree);
    #else
            MakeListUsingTreeWithStack(adr_tc_glb_, 0, adr_tp_glb_, 
                                       &n_epj_local_1task, id_epj_cpe_local, 
                                       &n_spj_local_1task, id_spj_cpe_local, 
                                       pos_target_box, *r_crit_sq, *len_peri_x,
                                       n_leaf_limit_, sp_tree_offset);
    #endif
            // Note that *len_peri_x is only used when the macro `PHI_R_TREE`
            // is defined.
        }
#if defined(DEBUG_MAKELIST)
        if ((n_epj_local_1task >= N_EPJ_PER_GROUP_MAX) ||
            (n_spj_local_1task >= N_SPJ_PER_GROUP_MAX)) {
            printf("n_epj_local_1task = %d\n",n_epj_local_1task);
            printf("n_spj_local_1task = %d\n",n_spj_local_1task);
        }
        assert(n_epj_local_1task < N_EPJ_PER_GROUP_MAX);
        assert(n_spj_local_1task < N_SPJ_PER_GROUP_MAX);
#endif
        //** Set n_epj_cpe_local & n_spj_cpe_local
        *n_epj_cpe_local = n_epj_local_1task;
        *n_spj_cpe_local = n_spj_local_1task;
        //** Broadcast n_epj, n_spj to compute offset
        n_epj_tbl[my_id] = n_epj_local_1task;
        n_spj_tbl[my_id] = n_spj_local_1task;
        n_epj_total_1task = 0;
        n_spj_total_1task = 0;
        for (id=0; id<NUMBER_OF_CPE; id++) {
            cpe_bcast_int32(id,&n_epj_tbl[id]);
            cpe_bcast_int32(id,&n_spj_tbl[id]);
            n_epj_total_1task += n_epj_tbl[id];
            n_spj_total_1task += n_spj_tbl[id];
        }
        //** Compute n_disp_epj_local, n_disp_spj_local
        n_disp_epj_local = 0;
        n_disp_spj_local = 0;
        for (id=0; id<my_id; id++) {
           n_disp_epj_local += n_epj_tbl[id];
           n_disp_spj_local += n_spj_tbl[id];
        }
        *n_disp_epj_cpe_local = n_epj_total + n_disp_epj_local;
        *n_disp_spj_cpe_local = n_spj_total + n_disp_spj_local;
        n_epj_total += n_epj_total_1task;
        n_spj_total += n_spj_total_1task;
        //** Check if there are sufficient memory 
        n_epj_shortfall = (n_epj_total > ca_id_epj) ? (n_epj_total - ca_id_epj) : 0;
        n_spj_shortfall = (n_spj_total > ca_id_spj) ? (n_spj_total - ca_id_spj) : 0;
        if ((n_epj_shortfall > 0) || (n_spj_shortfall > 0)) {
            if (my_id == 0) {
                // Put n_epj_shortfall
                reply_put = 0;
                dma_set_op(&dma_put, DMA_PUT);
                dma_set_mode(&dma_put, PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put, sizeof(int));
                dma(dma_put, (long*)(adr_n_epj_shortfall), (long*)(&n_epj_shortfall));
                dma_wait(&reply_put, 1);
                while (reply_put !=1 ) {}
                // Put n_spj_shortfall
                reply_put = 0;
                dma_set_op(&dma_put, DMA_PUT);
                dma_set_mode(&dma_put, PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put, sizeof(int));
                dma(dma_put, (long*)(adr_n_spj_shortfall), (long*)(&n_spj_shortfall));
                dma_wait(&reply_put, 1);
                while (reply_put !=1 ) {}
            }
            break;
        }
        //** Put the data of interaction lists
        if (ig != -1) {
            // (i) id_epj_cpe_local
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, sizeof(int) * (*n_epj_cpe_local));
            dma(dma_put, (long*)((int *)adr_id_epj + (*n_disp_epj_cpe_local)), (long*)(id_epj_cpe_local));
            dma_wait(&reply_put, 1);
            while (reply_put !=1 ) {}
            // (ii) id_spj_cpe_local
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, sizeof(int) * (*n_spj_cpe_local));
            dma(dma_put, (long*)((int *)adr_id_spj + (*n_disp_spj_cpe_local)), (long*)(id_spj_cpe_local));
            dma_wait(&reply_put, 1);
            while (reply_put !=1 ) {}
            // (iii) n_disp_epj_cpe_local
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, sizeof(int));
            dma(dma_put, (long*)((int *)adr_n_disp_epj + ig), (long*)(n_disp_epj_cpe_local));
            dma_wait(&reply_put, 1);
            while (reply_put !=1 ) {}
            // (iv) n_disp_spj_cpe_local
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, sizeof(int));
            dma(dma_put, (long*)((int *)adr_n_disp_spj + ig), (long*)(n_disp_spj_cpe_local));
            dma_wait(&reply_put, 1);
            while (reply_put !=1 ) {}
            // (v) n_epj_cpe_local
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, sizeof(int));
            dma(dma_put, (long*)((int *)adr_n_epj + ig), (long*)(n_epj_cpe_local));
            dma_wait(&reply_put, 1);
            while (reply_put !=1 ) {}
            // (vi) n_spj_cpe_local
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, sizeof(int));
            dma(dma_put, (long*)((int *)adr_n_spj + ig), (long*)(n_spj_cpe_local));
            dma_wait(&reply_put, 1);
            while (reply_put !=1 ) {}
        }
    }

    //* Release memory 
    ldm_free(id_epj_cpe_local, bsize_id_epj_array );
    ldm_free(id_spj_cpe_local, bsize_id_spj_array );
    ldm_free(n_disp_epj_cpe_local, bsize_n_disp_epj );
    ldm_free(n_disp_spj_cpe_local, bsize_n_disp_spj );
    ldm_free(n_epj_cpe_local, bsize_n_epj );
    ldm_free(n_spj_cpe_local, bsize_n_spj );
    ldm_free(n_epj_tbl, bsize_n_epj * NUMBER_OF_CPE );
    ldm_free(n_spj_tbl, bsize_n_spj * NUMBER_OF_CPE );
    ldm_free(ipg_, bsize_ipg );
    ldm_free(r_crit_sq, sizeof(double) );
    ldm_free(len_peri_x, sizeof(double) );

#ifdef REMOVE_VERTEX
    ldm_free(r_search, sizeof(double) );
    ldm_free(half_len_tree, sizeof(double) );
    ldm_free(cen_tree, sizeof(double)*3 );
#endif

}

