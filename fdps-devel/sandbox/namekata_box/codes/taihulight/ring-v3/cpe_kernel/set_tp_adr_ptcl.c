#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include"slave.h"
#include"dma.h"
#include"cpe_func.h"
#include"cpe_prof.h"

extern int MY_RANK_MPI;

static inline U32_ GetMsb(U32_ val){
    return (val>>31) & 0x1;
}

static inline U32_ SetMsb(U32_ val){
    return val | 0x80000000;
}

static inline U32_ ClearMsb(U32_ val){
    return val & 0x7fffffff;
}

static inline void CopyKeyOnly(tpLM * src, tpLM * dst){
    dst->key_hi_ = src->key_hi_;
#ifdef USE_96BIT_KEY
    dst->key_lo_ = src->key_lo_;
#endif
}

/* CPE communication functions */
static inline void cpe_bcast_int32(const int root_cpe_id, int* data) {
    int my_id = athread_get_id(-1);
    int my_col_id, my_row_id;
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
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

typedef int sumtype;
static void prefix_sum(sumtype val, sumtype *beg, sumtype *end, sumtype *psum){
	const int myid = athread_get_id(-1);
	const int jcol = myid % 8;
	const int irow = myid / 8;
	sumtype sbuf, rbuf, sum;
	sumtype v0 = val;

	sync_array_();
	// __shfl_up(1) 
	{
		sbuf = val;
		REG_PUTR(sbuf, (jcol+1)%8);
		REG_GETR(rbuf);
		if(jcol < 1){
			sbuf = rbuf;
			REG_PUTC(sbuf, (irow+1)%8);
			REG_GETC(rbuf);
		}
		if(myid >= 1){
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(2) 
	{
		sbuf = val;
		REG_PUTR(sbuf, (jcol+2)%8);
		REG_GETR(rbuf);
		if(jcol < 2){
			sbuf = rbuf;
			REG_PUTC(sbuf, (irow+1)%8);
			REG_GETC(rbuf);
		}
		if(myid >= 2){
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(4) 
	{
		sbuf = val;
		REG_PUTR(sbuf, (jcol+4)%8);
		REG_GETR(rbuf);
		if(jcol < 4){
			sbuf = rbuf;
			REG_PUTC(sbuf, (irow+1)%8);
			REG_GETC(rbuf);
		}
		if(myid >= 4){
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(8) 
	{
		if(myid < 56){
			sbuf = val;
			REG_PUTC(sbuf, (irow+1)%8);
		}
		if(myid >= 8){
			REG_GETC(rbuf);
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(16) 
	{
		if(myid < 48){
			sbuf = val;
			REG_PUTC(sbuf, (irow+2)%8);
		}
		if(myid >= 16){
			REG_GETC(rbuf);
			val += rbuf;
		}
	}
	sync_array_();

	// __shfl_up(32) 
	{
		if(myid < 32){
			sbuf = val;
			REG_PUTC(sbuf, (irow+4)%8);
		}
		if(myid >= 32){
			REG_GETC(rbuf);
			val += rbuf;
		}
	}
	sync_array_();

	// bcast the sum
	if(7 == irow){
		if(7 == jcol){
			sum = val;
			REG_PUTR(sum, 8);
		}else{
			REG_GETR(sum);
		}
		REG_PUTC(sum, 8);
	}else{
		REG_GETC(sum);
	}
	sync_array_();

	*beg  = val-v0;
	*end  = val;
	*psum = sum;
}

void SetTpAdrPtcl(void * args){
#if 1
    volatile int my_id, my_col_id, my_row_id;
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    //* Initialize local vars. for DMA comm.
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    //* Process the arguments
    int nptcl     = (int)(((unsigned long*)args)[0]);
    void *adr_tp  = (void *)((unsigned long*)args)[1];
    //* Local variables
    //-(loop counters)
    int i,j,k;
    //-(prefix sum calc.)
    volatile int beg,end,psum;
    //-(particle numbers)
    int n_epj_loc, n_spj_loc;
    int n_epj_tot, n_spj_tot;
    int n_disp_epj_loc, n_disp_spj_loc;
    //-(buffers)
    enum {
        CHUNK_SIZE = 64,
        CYCLE_SIZE = CHUNK_SIZE * NUMBER_OF_CPE,
    };
    size_t bsize_tp = sizeof(tpLM);
    size_t bsize_tp_chunk = bsize_tp * CHUNK_SIZE;
    volatile tpLM * tp = (tpLM *) ldm_malloc(bsize_tp_chunk);

    //* Update tp[]
    n_epj_tot = n_spj_tot = 0;
    for (i=0; i<nptcl; i+=CYCLE_SIZE) {
        int nrem = nptcl - i;
        int nn   = nrem < CYCLE_SIZE ? nrem : CYCLE_SIZE;
        //* Compute the task of each CPE
        int n_loc,my_offset;
        n_loc = nn/NUMBER_OF_CPE + ( (my_id < nn % NUMBER_OF_CPE) ? 1 : 0 );
        my_offset = i + (nn/NUMBER_OF_CPE)*my_id + ( (my_id < nn % NUMBER_OF_CPE) ? my_id : nn % NUMBER_OF_CPE );
        //* Get tp[] and counts the numbers of EPJ and SPJ
        n_epj_loc = n_spj_loc = 0;
        if (n_loc > 0) {
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tp * n_loc);
            dma(dma_get, (long*)((tpLM *)adr_tp + my_offset), (long*)(tp));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            for (j=0; j<n_loc; j++) {
                if (GetMsb(tp[j].adr_ptcl_) == 0) {
                    n_epj_loc++;
                } else {
                    n_spj_loc++;
                }
            }
        }
        //* Compute the prefix sum (Collective CPE communication)
        prefix_sum(n_epj_loc,&beg,&end,&psum);
        n_disp_epj_loc = n_epj_tot + beg;
        n_epj_tot += psum;
        prefix_sum(n_spj_loc,&beg,&end,&psum);
        n_disp_spj_loc = n_spj_tot + beg;
        n_spj_tot += psum; 
        //* Update tp[] and Put tp[]
        n_epj_loc = n_spj_loc = 0; // reset
        if (n_loc > 0) {
            for (j=0; j<n_loc; j++) {
                if (GetMsb(tp[j].adr_ptcl_) == 0) {
                    tp[j].adr_ptcl_ = n_disp_epj_loc + n_epj_loc;
                    n_epj_loc++;
                } else {
                    tp[j].adr_ptcl_ = SetMsb(n_disp_spj_loc + n_spj_loc);
                    n_spj_loc++;
                }
            }
            //* Put
            reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put, bsize_tp * n_loc);
            dma(dma_put, (long*)((tpLM *)adr_tp + my_offset), (long*)(tp));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }
    }

    //* Release memory
    ldm_free(tp, bsize_tp_chunk);
#else
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n_tot = (int)(((unsigned long*)args)[0]);
    void * adr_tp = (void *)((unsigned long*)args)[1];
    int my_n = n_tot/NUMBER_OF_CPE + ( (my_id < n_tot % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot % NUMBER_OF_CPE) ? my_id : n_tot % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tp = sizeof(tpLM);
    tpLM  * tp     = (tpLM *) ldm_malloc(bsize_tp*CHUNK_SIZE);
    tpLM  * tp_new = (tpLM *) ldm_malloc(bsize_tp*CHUNK_SIZE);
    int my_n_ep = 0;
    int my_n_sp = 0;
    // just count number to determine offset to write
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            if(GetMsb(tp[j].adr_ptcl_)==0){
                my_n_ep++;
            }
            else{
                my_n_sp++;
            }
        }
    }
    /*
    sync_array_();
    for(i=0;i<my_id*10000000; i++) NOP();
    if(MY_RANK_MPI==0){
        printf("A) my_id=%d, my_n_ep=%d , my_n_sp=%d \n",
               my_id, my_n_ep, my_n_sp);
    }
    sync_array_();
    */
    
    int n_ep_ar[NUMBER_OF_CPE];
    int n_sp_ar[NUMBER_OF_CPE];
    for (i=0; i<NUMBER_OF_CPE; i++) {
        n_ep_ar[i] = n_sp_ar[i] = 0;
    }
    n_ep_ar[my_id] = my_n_ep;
    n_sp_ar[my_id] = my_n_sp;
    for (i=0; i<NUMBER_OF_CPE; i++) {
        cpe_bcast_int32(i, &n_ep_ar[i]);
        cpe_bcast_int32(i, &n_sp_ar[i]);
    }
    int my_n_disp_ep = 0;
    int my_n_disp_sp = 0;
    for(i=0; i<my_id; i++){
        my_n_disp_ep += n_ep_ar[i];
        my_n_disp_sp += n_sp_ar[i];
    }
    /*
    sync_array_();
    for(i=0;i<my_id*10000000; i++) NOP();
    if(MY_RANK_MPI==0){
        printf("B) my_id=%d, my_n_disp_ep=%d , my_n_disp_sp=%d \n",
               my_id, my_n_disp_ep, my_n_disp_sp);
    }
    sync_array_();
    */
    my_n_ep = my_n_sp = 0;
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            CopyKeyOnly(&(tp[j]), &(tp_new[j]));
            if(GetMsb(tp[j].adr_ptcl_)==0){
                tp_new[j].adr_ptcl_ = my_n_disp_ep+my_n_ep;
                my_n_ep++;
            }
            else{
                tp_new[j].adr_ptcl_ = SetMsb(my_n_disp_sp+my_n_sp);
                my_n_sp++;
            }
        }
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tp * nn);
        dma(dma_put, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp_new));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    //* Release memory
    ldm_free(tp, bsize_tp*CHUNK_SIZE);
    ldm_free(tp_new, bsize_tp*CHUNK_SIZE);
#endif
}


void SetTpAdrPtclEpOnly(void * args){
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n_tot = (int)(((unsigned long*)args)[0]);
    void * adr_tp = (void *)((unsigned long*)args)[1];
    int my_n = n_tot/NUMBER_OF_CPE + ( (my_id < n_tot % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot % NUMBER_OF_CPE) ? my_id : n_tot % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tp = sizeof(tpLM);
    tpLM  * tp     = (tpLM *) ldm_malloc(bsize_tp*CHUNK_SIZE);
    tpLM  * tp_new = (tpLM *) ldm_malloc(bsize_tp*CHUNK_SIZE);
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            CopyKeyOnly(&(tp[j]), &(tp_new[j]));
            tp_new[j].adr_ptcl_ = my_offset+i+j;
        }
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tp * nn);
        dma(dma_put, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp_new));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    //* Release memory
    ldm_free(tp, bsize_tp*CHUNK_SIZE);
    ldm_free(tp_new, bsize_tp*CHUNK_SIZE);
}

void slave_SetTpLocAdrPtclFromTpCpe(void * args){
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n = (int)(((unsigned long*)args)[0]);
    void * adr_tp       = (void *)((unsigned long*)args)[1];
    void * adr_adr_ptcl = (void *)((unsigned long*)args)[2];
    int my_n = n/NUMBER_OF_CPE + ( (my_id < n % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n/NUMBER_OF_CPE)*my_id + ( (my_id < n % NUMBER_OF_CPE) ? my_id : n % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tp  = sizeof(tpLM);
    size_t bsize_adr = sizeof(U64_);
    tpLM  * tp      = (tpLM *) ldm_malloc(bsize_tp  * CHUNK_SIZE);
    U64_  * adr_ptcl = (U64_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            adr_ptcl[j] = tp[j].adr_ptcl_;
        }
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_adr * nn);
        dma(dma_put, (long*)((U64_ *)adr_adr_ptcl+my_offset+i), (long*)(adr_ptcl));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    //* Release memory
    ldm_free(tp, bsize_tp*CHUNK_SIZE);
    ldm_free(adr_ptcl, bsize_adr*CHUNK_SIZE);
}

void SetAdrGlb1stCpe(void * args){
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n = (int)(((unsigned long*)args)[0]);
    void * adr_tp       = (void *)((unsigned long*)args)[1];
    void * adr_adr_epj_org2glb = (void *)((unsigned long*)args)[2];
    void * adr_adr_epj_loc2glb = (void *)((unsigned long*)args)[3];
    int my_n = n / NUMBER_OF_CPE + ( (my_id < n % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n / NUMBER_OF_CPE)*my_id + ( (my_id < n % NUMBER_OF_CPE) ? my_id : n % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
#ifdef REMOVE_TP_LOC
    size_t bsize_tp = sizeof(U64_);
    U64_ * tp = (U64_ *)  ldm_malloc(bsize_tp * CHUNK_SIZE);
#else
    size_t bsize_tp = sizeof(tpLM);
    tpLM * tp = (tpLM *)  ldm_malloc(bsize_tp * CHUNK_SIZE);
#endif

    size_t bsize_adr = sizeof(S32_);
    S32_  * adr_epj_org2glb = (S32_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    S32_  * adr_epj_loc2glb = (S32_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    S32_ adr_org_tmp[CHUNK_SIZE];
    S32_ adr_glb_tmp[CHUNK_SIZE];
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        
        // get tp
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
#ifdef REMOVE_TP_LOC
        dma(dma_get, (long*)((U64_ *)adr_tp+my_offset+i), (long*)(tp));
#else
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
#endif
        while (reply_get != 1) {}
        dma_wait(&reply_get, 1);

        // prepare for getting adr_epj_org2glb
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_adr);
        for(j=0; j<nn; j++){
#ifdef REMOVE_TP_LOC
            adr_org_tmp[j] = (S32_)tp[j];
#else
            adr_org_tmp[j] = (S32_)tp[j].adr_ptcl_;
#endif
            dma(dma_get, (long*)((S32_ *)adr_adr_epj_org2glb+adr_org_tmp[j]),
                (long*)(&adr_glb_tmp[j]));
        }
        while (reply_get != nn) {}
        dma_wait(&reply_get, nn);

        // prepare for putting adr_epj_loc2glb
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_adr);
        for(j=0; j<nn; j++){
            dma(dma_put, (long*)((S32_ *)adr_adr_epj_loc2glb+my_offset+i+j),
                (long*)(&adr_glb_tmp[j]));
        }
        while (reply_put != nn) {}
        dma_wait(&reply_put, nn);
    }
    //* Release memory
    ldm_free(tp, bsize_tp*CHUNK_SIZE);
    ldm_free(adr_epj_org2glb, bsize_adr*CHUNK_SIZE);
    ldm_free(adr_epj_loc2glb, bsize_adr*CHUNK_SIZE);
}

void SetAdrGlb2ndCpe(void * args){
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n_epj_sorted = (int)(((unsigned long*)args)[0]);
    int n_adr = (int)(((unsigned long*)args)[1]);
    void * adr_adr_epj_glb2loc = (void *)((unsigned long*)args)[2];
    void * adr_adr_epj_loc2glb = (void *)((unsigned long*)args)[3];
    int my_n_epj_sorted = n_epj_sorted/NUMBER_OF_CPE
        + ( (my_id < n_epj_sorted % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset_epj_sorted = (n_epj_sorted / NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted % NUMBER_OF_CPE) ? my_id : n_epj_sorted % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_adr = sizeof(S32_);
    S32_  * adr_epj_glb2loc = (S32_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    S32_  * adr_epj_loc2glb = (S32_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    for (i=0; i<my_n_epj_sorted; i+=CHUNK_SIZE) {
        int nrem = my_n_epj_sorted - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_adr * nn);
        dma(dma_get, (long*)((S32_ *)adr_adr_epj_glb2loc+my_offset_epj_sorted+i),
            (long*)(adr_epj_glb2loc));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            adr_epj_glb2loc[j] = -1;
        }
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_adr * nn);
        dma(dma_put, (long*)((S32_ *)adr_adr_epj_glb2loc+my_offset_epj_sorted+i),
            (long*)(adr_epj_glb2loc));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    sync_array_();
    // buffer for DMA
    S32_ adr_tmp[CHUNK_SIZE]; 
    S32_ val_tmp[CHUNK_SIZE];
    int my_n_adr = n_adr / NUMBER_OF_CPE
        + ( (my_id < n_adr % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset_adr = (n_adr / NUMBER_OF_CPE)*my_id + ( (my_id < n_adr % NUMBER_OF_CPE) ? my_id : n_adr % NUMBER_OF_CPE );
    for (i=0; i<my_n_adr; i+=CHUNK_SIZE) {
        int nrem = my_n_adr - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        // get adr_epj_loc2glb
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_adr * nn);
        dma(dma_get, (long*)((S32_ *)adr_adr_epj_loc2glb+my_offset_adr+i),
            (long*)(adr_epj_loc2glb));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        // prepare for puting adr_epj_glb2loc
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_adr);
        for (j=0; j<nn; j++) {
            adr_tmp[j] = adr_epj_loc2glb[j];
            val_tmp[j] = my_offset_adr + i + j;
            dma(dma_put, (long*)((S32_ *)adr_adr_epj_glb2loc+adr_tmp[j]),
                (long*)(&val_tmp[j]));
        }

        dma_wait(&reply_put, nn);
        while (reply_put != nn) {}
    }

    //* Release memory
    ldm_free(adr_epj_glb2loc, bsize_adr*CHUNK_SIZE);
    ldm_free(adr_epj_loc2glb, bsize_adr*CHUNK_SIZE);
}


void SetPtclSortedGlbCpe(void * args){
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    dma_desc dma_put_adr_epj_org, dma_put_adr_epj_buf, dma_put_adr_spj_buf;
    volatile int reply_put_adr_epj_org = 0;
    volatile int reply_put_adr_epj_buf = 0;
    volatile int reply_put_adr_spj_buf = 0;
    dma_descriptor_init(&dma_put_adr_epj_org, &reply_put_adr_epj_org);
    dma_descriptor_init(&dma_put_adr_epj_buf, &reply_put_adr_epj_buf);
    dma_descriptor_init(&dma_put_adr_spj_buf, &reply_put_adr_spj_buf);
    
    dma_desc dma_get_epj_org, dma_get_spj_org;
    volatile int reply_get_epj_org = 0;
    volatile int reply_get_spj_org = 0;
    dma_descriptor_init(&dma_get_epj_org, &reply_get_epj_org);
    dma_descriptor_init(&dma_get_spj_org, &reply_get_spj_org);

    dma_desc dma_put_epj_sorted, dma_put_spj_sorted;
    volatile int reply_put_epj_sorted = 0;
    volatile int reply_put_spj_sorted = 0;
    dma_descriptor_init(&dma_put_epj_sorted, &reply_put_epj_sorted);
    dma_descriptor_init(&dma_put_spj_sorted, &reply_put_spj_sorted);
    
    int my_id = athread_get_id(-1);
    int n_glb_tot = (int)(((unsigned long*)args)[0]);
    int n_loc_tot = (int)(((unsigned long*)args)[1]);
    void * adr_tp = (void *)((unsigned long*)args)[2];
    void * adr_epj_sorted = (void *)((unsigned long*)args)[3];
    void * adr_epj_org = (void *)((unsigned long*)args)[4];
    void * adr_spj_sorted = (void *)((unsigned long*)args)[5];
    void * adr_spj_org = (void *)((unsigned long*)args)[6];
    void * adr_adr_epj_org2glb = (void *)((unsigned long*)args)[7];
    void * adr_adr_epj_buf2glb = (void *)((unsigned long*)args)[8];
    void * adr_adr_spj_buf2glb = (void *)((unsigned long*)args)[9];
    int my_n_glb = n_glb_tot/NUMBER_OF_CPE + ( (my_id < n_glb_tot % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset_glb = (n_glb_tot/NUMBER_OF_CPE)*my_id + ( (my_id < n_glb_tot % NUMBER_OF_CPE) ? my_id : n_glb_tot % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tp  = sizeof(tpLM);
    size_t bsize_adr = sizeof(S32_);
    size_t bsize_ep  = sizeof(epjLM);
    size_t bsize_sp  = sizeof(spjLM);
    tpLM tp[CHUNK_SIZE];
    tpLM tp_new[CHUNK_SIZE];
    volatile int my_n_epj = 0;
    volatile int my_n_spj = 0;
    // just count number to determine offset to write
    for (i=0; i<my_n_glb; i+=CHUNK_SIZE) {
        int nrem = my_n_glb - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset_glb+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            if(GetMsb(tp[j].adr_ptcl_)==0){
                my_n_epj++;
            }
            else{
                my_n_spj++;
            }
        }
    }
    sync_array_();

    // get offset to write by scan
    //int n_epj_ar[NUMBER_OF_CPE];
    //int n_spj_ar[NUMBER_OF_CPE];
    int * n_epj_ar = (int *) ldm_malloc(sizeof(int)*NUMBER_OF_CPE);
    int * n_spj_ar = (int *) ldm_malloc(sizeof(int)*NUMBER_OF_CPE);
    
    for (i=0; i<NUMBER_OF_CPE; i++) {
        n_epj_ar[i] = n_spj_ar[i] = 0;
    }
    n_epj_ar[my_id] = my_n_epj;
    n_spj_ar[my_id] = my_n_spj;
    sync_array_();
    for (i=0; i<NUMBER_OF_CPE; i++) {
        cpe_bcast_int32(i, &n_epj_ar[i]);
        cpe_bcast_int32(i, &n_spj_ar[i]);
    }
    sync_array_();
    volatile int my_n_disp_epj = 0;
    volatile int my_n_disp_spj = 0;
    for(i=0; i<my_id; i++){
        my_n_disp_epj += n_epj_ar[i];
        my_n_disp_spj += n_spj_ar[i];
    }
    ldm_free(n_epj_ar, sizeof(int)*NUMBER_OF_CPE);
    ldm_free(n_spj_ar, sizeof(int)*NUMBER_OF_CPE);
    sync_array_();

    /*
    EpjLM epj_org[CHUNK_SIZE];
    SpjLM spj_org[CHUNK_SIZE];
    S32_ my_n_epj_ar[CHUNK_SIZE];
    S32_ my_n_spj_ar[CHUNK_SIZE];
    S32_ adr_epj_ar[CHUNK_SIZE];
    S32_ adr_spj_ar[CHUNK_SIZE];
    */
    epjLM * epj_org = (epjLM *) ldm_malloc(bsize_ep * CHUNK_SIZE);
    spjLM * spj_org = (spjLM *) ldm_malloc(bsize_sp * CHUNK_SIZE);
    S32_ * my_n_epj_ar = (S32_*) ldm_malloc(sizeof(S32_)*CHUNK_SIZE);
    S32_ * my_n_spj_ar = (S32_*) ldm_malloc(sizeof(S32_)*CHUNK_SIZE);
    S32_ * adr_epj_ar  = (S32_*) ldm_malloc(sizeof(S32_)*CHUNK_SIZE);
    S32_ * adr_spj_ar  = (S32_*) ldm_malloc(sizeof(S32_)*CHUNK_SIZE);
    
    sync_array_();
    my_n_epj = my_n_spj = 0;
    for (i=0; i<my_n_glb; i+=CHUNK_SIZE) {
        int nrem = my_n_glb - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        // get tp
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset_glb+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        
        ////////////
        // to put adr
        // preparing for puting adr_epj_org
        reply_put_adr_epj_org = 0;
        dma_set_op(&dma_put_adr_epj_org, DMA_PUT);
        dma_set_mode(&dma_put_adr_epj_org, PE_MODE);
        dma_set_reply(&dma_put_adr_epj_org, &reply_put_adr_epj_org);
        dma_set_size(&dma_put_adr_epj_org, bsize_adr);

        // preparing for puting adr_epj_buf
        reply_put_adr_epj_buf = 0;
        dma_set_op(&dma_put_adr_epj_buf, DMA_PUT);
        dma_set_mode(&dma_put_adr_epj_buf, PE_MODE);
        dma_set_reply(&dma_put_adr_epj_buf, &reply_put_adr_epj_buf);
        dma_set_size(&dma_put_adr_epj_buf, bsize_adr);
        
        // preparing for puting adr_spj_buf
        reply_put_adr_spj_buf = 0;
        dma_set_op(&dma_put_adr_spj_buf, DMA_PUT);
        dma_set_mode(&dma_put_adr_spj_buf, PE_MODE);
        dma_set_reply(&dma_put_adr_spj_buf, &reply_put_adr_spj_buf);
        dma_set_size(&dma_put_adr_spj_buf, bsize_adr);
        // to put adr
        ////////////
        
        ////////////
        // to get epj_org and spj_org to put epj_sorted and spj_sorted
        // prepare for gettin epj_org
        reply_get_epj_org = 0;
        dma_set_op(&dma_get_epj_org, DMA_GET);
        dma_set_mode(&dma_get_epj_org, PE_MODE);
        dma_set_reply(&dma_get_epj_org, &reply_get_epj_org);
        dma_set_size(&dma_get_epj_org, bsize_ep);

        // prepare for gettin spj_org
        reply_get_spj_org = 0;
        dma_set_op(&dma_get_spj_org, DMA_GET);
        dma_set_mode(&dma_get_spj_org, PE_MODE);
        dma_set_reply(&dma_get_spj_org, &reply_get_spj_org);
        dma_set_size(&dma_get_spj_org, bsize_sp);
        // to get epj_org and spj_org to put epj_sorted and spj_sorted
        ////////////
        
        volatile int my_n_epj_tmp = 0;
        volatile int my_n_epj_org_tmp = 0;
        volatile int my_n_epj_buf_tmp = 0;
        volatile int my_n_spj_tmp = 0;
        for (j=0; j<nn; j++) {
            U32_ adr = tp[j].adr_ptcl_;
            CopyKeyOnly(&(tp[j]), &(tp_new[j]));
            if(GetMsb(adr)==0){
                tp_new[j].adr_ptcl_ = my_n_disp_epj + my_n_epj;
                adr_epj_ar[my_n_epj_tmp] = adr;
                my_n_epj_ar[my_n_epj_tmp] = my_n_disp_epj + my_n_epj;
                if(adr < n_loc_tot){ // particles are in own process
                    dma(dma_put_adr_epj_org,
                        (long*)((S32_*)adr_adr_epj_org2glb+adr_epj_ar[my_n_epj_tmp]),
                        (long*)((S32_*)my_n_epj_ar+my_n_epj_tmp));
                    my_n_epj_org_tmp++;
                }
                else{
                    dma(dma_put_adr_epj_buf,
                        (long*)((S32_*)adr_adr_epj_buf2glb+adr_epj_ar[my_n_epj_tmp]-n_loc_tot),
                        (long*)((S32_*)my_n_epj_ar+my_n_epj_tmp));
                    my_n_epj_buf_tmp++;
                }
                dma(dma_get_epj_org,
                    (long*)((epjLM*)adr_epj_org+adr_epj_ar[my_n_epj_tmp]),
                    (long*)((epjLM*)epj_org+my_n_epj_tmp));
                my_n_epj++;
                my_n_epj_tmp++;
            }
            else{
                tp_new[j].adr_ptcl_ = SetMsb(my_n_disp_spj + my_n_spj); // set tp
                adr_spj_ar[my_n_spj_tmp]  = ClearMsb(adr); // temporally store
                my_n_spj_ar[my_n_spj_tmp] = my_n_disp_spj + my_n_spj; // temporally store
                dma(dma_put_adr_spj_buf,
                    (long*)((S32_*)adr_adr_spj_buf2glb+adr_spj_ar[my_n_spj_tmp]),
                    (long*)((S32_*)my_n_spj_ar+my_n_spj_tmp)); //put adr_spj_buf2glb
                dma(dma_get_spj_org,
                    (long*)((spjLM*)adr_spj_org+adr_spj_ar[my_n_spj_tmp]),
                    (long*)((spjLM*)spj_org+my_n_spj_tmp));
                my_n_spj++;
                my_n_spj_tmp++;
            }
        }
        ///////////
        // put tp
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tp * nn);
        dma(dma_put, (long*)((tpLM *)adr_tp+my_offset_glb+i), (long*)(tp_new));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
        // put tp
        ///////////
        
        ///////////
        // put tp, adr_epj_org2glb, adr_epj_buf2glb, adr_spj_buf2glb
        // wait puting adr_epj_org2glb
        dma_wait(&reply_put_adr_epj_org, my_n_epj_org_tmp);
        while (reply_put_adr_epj_org != my_n_epj_org_tmp) {}
        
        // wait puting adr_epj_buf2glb
        dma_wait(&reply_put_adr_epj_buf, my_n_epj_buf_tmp);
        while (reply_put_adr_epj_buf !=  my_n_epj_buf_tmp){}
        
        // wait puting adr_spj_buf2glb
        dma_wait(&reply_put_adr_spj_buf, my_n_spj_tmp);
        while (reply_put_adr_spj_buf != my_n_spj_tmp) {}
        // put tp, adr_epj_org2glb, adr_epj_buf2glb, adr_spj_buf2glb
        ///////////

        /////////
        // get epj_org, spj_org to put epj_sorted and spj_sorted
        // wait getting epj_org
        dma_wait(&reply_get_epj_org, my_n_epj_tmp);
        while (reply_get_epj_org != my_n_epj_tmp) {}
        
        // wait getting spj_org
        dma_wait(&reply_get_spj_org, my_n_spj_tmp);
        while (reply_get_spj_org != my_n_spj_tmp) {}
        // get epj_org, spj_org
        /////////

        ///////////
        // put epj_sorted, spj_sorted
        // put epj_sorted
        if( my_n_epj_tmp > 0) {
            reply_put_epj_sorted = 0;
            dma_set_op(&dma_put_epj_sorted, DMA_PUT);
            dma_set_mode(&dma_put_epj_sorted, PE_MODE);
            dma_set_reply(&dma_put_epj_sorted, &reply_put_epj_sorted);
            dma_set_size(&dma_put_epj_sorted, bsize_ep*my_n_epj_tmp);
            dma(dma_put_epj_sorted,
                (long*)((epjLM*)adr_epj_sorted + (my_n_disp_epj + my_n_epj - my_n_epj_tmp)),
                (long*)(epj_org));
            dma_wait(&reply_put_epj_sorted, 1);
            while (reply_put_epj_sorted != 1) {}
        }

        // put spj_sorted
        if( my_n_spj_tmp > 0) {
            reply_put_spj_sorted = 0;
            dma_set_op(&dma_put_spj_sorted, DMA_PUT);
            dma_set_mode(&dma_put_spj_sorted, PE_MODE);
            dma_set_reply(&dma_put_spj_sorted, &reply_put_spj_sorted);
            dma_set_size(&dma_put_spj_sorted, bsize_sp*my_n_spj_tmp);
            dma(dma_put_spj_sorted,
                (long*)((spjLM*)adr_spj_sorted + (my_n_disp_spj + my_n_spj - my_n_spj_tmp)),
                (long*)((spjLM*)spj_org));
            dma_wait(&reply_put_spj_sorted, 1);
            while (reply_put_spj_sorted != 1) {}
        }
        // put epj_sorted, spj_sorted
        ///////////

    }
    
    ldm_free(epj_org, bsize_ep*CHUNK_SIZE);
    ldm_free(spj_org, bsize_sp*CHUNK_SIZE);
    ldm_free(my_n_epj_ar, sizeof(S32_)*CHUNK_SIZE);
    ldm_free(my_n_spj_ar, sizeof(S32_)*CHUNK_SIZE);
    ldm_free(adr_epj_ar, sizeof(S32_)*CHUNK_SIZE);
    ldm_free(adr_spj_ar, sizeof(S32_)*CHUNK_SIZE);

}

