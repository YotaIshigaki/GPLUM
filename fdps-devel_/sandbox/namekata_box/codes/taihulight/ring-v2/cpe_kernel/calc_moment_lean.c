/* Standard C headers */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
/* Headers specific to Sunway TaihuLight */
#include "slave.h"
#include "dma.h"
#include "cpe_func.h"
#include "cpe_prof.h"

void *ldm_malloc(size_t);
void ldm_free(void *, size_t);

#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))

#define tcLM error
// #define PROFILING

static inline void reduce_long(long *d) {
	int my_col, my_row;
	  get_col_id_(&my_col);
	  get_row_id_(&my_row);

	long t;
	unsigned int dest;
	if((my_row&1)==0) {
		dest = my_row+1;
		REG_PUTC(*d,dest);
	}
	else {
		REG_GETC(t);
		*d += t;

		if((my_row&2)==0) {
			dest = my_row+2;
			REG_PUTC(*d,dest);
		}
		else {
			REG_GETC(t);
			*d += t;

			if(my_row==3) {
				dest = 7;
				REG_PUTC(*d,dest);
			}
			else {
				REG_GETC(t);
				*d += t;

				if((my_col&1)==0) {
					dest = my_col+1;
					REG_PUTR(*d,dest);
				}
				else {
					REG_GETR(t);
					*d += t;

					if((my_col&2)==0) {
						dest = my_col+2;
						REG_PUTR(*d,dest);
					}
					else {
						REG_GETR(t);
						*d += t;

						if(my_col==3) {
							dest = 7;
							REG_PUTR(*d,dest);
						}
						else {
							REG_GETR(t);
							*d += t;
						}
					}
				}
			}
		}
	}
}


#ifdef PROFILING
typedef struct Profile{
	long get_tc;
	long get_tp;
	long get_epj;
	long get_spj;
	long get_chl;
	long put_tc;

	long num_get_tc;
	long num_get_tp;
	long num_get_epj;
	long num_get_spj;
	long num_get_chl;
	long num_put_tc;
} Profile, pProfile;
__thread_local Profile glb_prof;
#define BEG_PROF(name) do {long time ## name = rtc_();
#define END_PROF(name) glb_prof.name += rtc_() - time ## name; } while(0);
#define INC_CNT(name, cnt) glb_prof.num_ ##name += cnt
static void init_prof(){
	glb_prof.get_tc  = 0;
	glb_prof.get_tp  = 0;
	glb_prof.get_epj = 0;
	glb_prof.get_spj = 0;
	glb_prof.get_chl = 0;
	glb_prof.put_tc  = 0;
	glb_prof.num_get_tc  = 0;
	glb_prof.num_get_tp  = 0;
	glb_prof.num_get_epj = 0;
	glb_prof.num_get_spj = 0;
	glb_prof.num_get_chl = 0;
	glb_prof.num_put_tc  = 0;
}
static void print_prof(){
	printf("get_tc  : %e sec\n", glb_prof.get_tc / 1.45e9);
	printf("get_tp  : %e sec\n", glb_prof.get_tp / 1.45e9);
	printf("get_epj : %e sec\n", glb_prof.get_epj/ 1.45e9);
	printf("get_spj : %e sec\n", glb_prof.get_spj/ 1.45e9);
	printf("get_chl : %e sec\n", glb_prof.get_chl/ 1.45e9);
	printf("put_tc  : %e sec\n", glb_prof.put_tc / 1.45e9);
	printf("sum     : %e sec\n\n", 
			(glb_prof.get_tc +
			 glb_prof.get_tp             +
			 glb_prof.get_epj +
			 glb_prof.get_spj +
			 glb_prof.get_chl +
			 glb_prof.put_tc) / 1.45e9);
	double bw_get_tc  = ( glb_prof.num_get_tc  * sizeof(etcLM) ) / ( glb_prof.get_tc  / 1.45e9 );
	double bw_get_tp  = ( glb_prof.num_get_tp  * sizeof(tpLM)  ) / ( glb_prof.get_tp  / 1.45e9 );
	double bw_get_epj = ( glb_prof.num_get_epj * sizeof(epjLM) ) / ( glb_prof.get_epj / 1.45e9 );
	double bw_get_chl = ( glb_prof.num_get_chl * sizeof(etcLM) ) / ( glb_prof.get_chl / 1.45e9 );
	double bw_put_tc  = ( glb_prof.num_put_tc  * sizeof(etcLM) ) / ( glb_prof.put_tc  / 1.45e9 );

	printf("get_tc  : %e byte/sec\n", bw_get_tc);
	printf("get_tp  : %e byte/sec\n", bw_get_tp);
	printf("get_epj : %e byte/sec\n", bw_get_epj);
	// printf("get_spj : %e byte/sec\n", bw_get_spj);
	printf("get_chl : %e byte/sec\n", bw_get_chl);
	printf("put_tc  : %e byte/sec\n\n", bw_put_tc);

	printf("get_tc  : %ld byte\n", glb_prof.num_get_tc  * sizeof(etcLM));
	printf("get_tp  : %ld byte\n", glb_prof.num_get_tp  * sizeof(tpLM) );
	printf("get_epj : %ld byte\n", glb_prof.num_get_epj * sizeof(epjLM));
	printf("get_spj : %ld byte\n", glb_prof.num_get_spj * sizeof(spjLM));
	printf("get_chl : %ld byte\n", glb_prof.num_get_chl * sizeof(etcLM));
	printf("put_tc  : %ld byte\n", glb_prof.num_put_tc  * sizeof(etcLM));
	puts("\n");
}
static void print_bandwidth(double dt){
	long total = glb_prof.num_get_tc  * sizeof(etcLM) 
               + glb_prof.num_get_tp  * sizeof(tpLM ) 
               + glb_prof.num_get_epj * sizeof(epjLM)
               + glb_prof.num_get_spj * sizeof(spjLM)
               + glb_prof.num_get_chl * sizeof(etcLM) 
               + glb_prof.num_put_tc  * sizeof(etcLM);
	int my_id = athread_get_id(-1);
	// total = 1; // debug
	if(63 == my_id){
		printf("before : %ld (%ld M) byte\n", total, total>>20);
	}
	reduce_long(&total);
	if(63 == my_id){
		printf("after  : %ld (%ld M) byte\n", total, total>>20);
		printf("dt = %e\n", dt);
		double bw = total / dt;
		printf("Aggregate bandwidth: %e byte/sec\n", bw); 
	}
}
#else
#define BEG_PROF(name)
#define END_PROF(name)
#define init_prof()
#define print_prof()
#define INC_CNT(name, cnt)
#define print_bandwidth(dt)
#endif


static __thread_local dma_desc dma_get, dma_put;
static __thread_local volatile unsigned long reply_get;
static __thread_local volatile unsigned long reply_put;

__attribute__ ((noinline))
static void UpdateTc(etcLM *tc_tgt,
                     double rsrch,
                     int    n_leaf_limit,
                     etcLM *adr_tc_head,
                     epjLM *adr_epj_head,
                     int *ptr_epj_cache_beg){
    enum {
        CHUNK_SIZE = 64,
        EPJ_CACHE_SIZE = 32,
        //CHUNK_SIZE = 32,
        //EPJ_CACHE_SIZE = 16,
        //CHUNK_SIZE = 16,
        //EPJ_CACHE_SIZE = 8,
    };
    epjLM epj[CHUNK_SIZE];
    etcLM tc [N_CHILDREN];

    const int n_ptcl   = tc_tgt->n_ptcl_;
    const int adr_ptcl = tc_tgt->adr_ptcl_;

    int epj_cache_beg = *ptr_epj_cache_beg;
#if 0
    epjLM epj_cache[EPJ_CACHE_SIZE];
    int epj_cache_end = epj_cache_beg + EPJ_CACHE_SIZE;
    int epj_beg = adr_ptcl;
    int epj_end = epj_beg + n_ptcl;
#endif

    int i, j, k;
    
    double mass = 0.0;
    double pos_x = 0.0;
    double pos_y = 0.0;
    double pos_z = 0.0;
    /*
#ifdef PHI_R_TREE
    double pos_phi = 0.0;
    double pos_r = 0.0;
#endif
    */
    if (n_ptcl <= n_leaf_limit || tc_tgt->level_ == TREE_LEVEL_LIMIT){
#if 1
        for (j=0; j<n_ptcl; j+= CHUNK_SIZE) {
            int nrem = n_ptcl - j;
            int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
            //** Get epj
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, (unsigned *)&reply_get);
            dma_set_size(&dma_get, sizeof(epjLM) * nn);
            dma(dma_get, (long)((epjLM *)adr_epj_head + adr_ptcl + j), (long)(epj));
            dma_wait((unsigned *)&reply_get, 1);
            while (reply_get !=1 ) {}
            //** Moment calculation
            for (k=0; k<nn; k++){
                mass  += epj[k].mass;
                pos_x += epj[k].mass * epj[k].pos.x;
                pos_y += epj[k].mass * epj[k].pos.y;
                pos_z += epj[k].mass * epj[k].pos.z;
                /*
                  #ifdef PHI_R_TREE
                  pos_phi += epj[k].mass * epj[k].pos_phi;
                  pos_r   += epj[k].mass * epj[k].pos_r;
                  #endif // PHI_R_TREE
                */
            } // for(k)
        } // for(j)
#else
        /*
          while(epj_beg < epj_end){
			if(epj_beg >= epj_cache_end){
				reply_get = 0;
				dma_set_op(&dma_get, DMA_GET);
				dma_set_mode(&dma_get, PE_MODE);
				dma_set_reply(&dma_get, (unsigned *)&reply_get);
				dma_set_size(&dma_get, sizeof(epjLM) * EPJ_CACHE_SIZE);
				dma(dma_get, (long)((epjLM *)adr_epj_head + epj_beg), (long)(epj_cache));
				// dma_wait((unsigned *)&reply_get, 1);
				while (reply_get !=1 ) {}

				epj_cache_beg = epj_beg;
				epj_cache_end = epj_cache_beg + EPJ_CACHE_SIZE;
			}
			int nn = (epj_end > epj_cache_end) ? (epj_cache_end - epj_beg)
			                                   : (epj_end       - epj_beg);
			int addr = epj_beg - epj_cache_beg;
			assert(addr >= 0);
			assert(addr+nn <= EPJ_CACHE_SIZE);
			epjLM *epj = epj_cache + addr;
			for (k=0; k<nn; k++){
				mass  += epj[k].mass;
				pos_x += epj[k].mass * epj[k].pos.x;
				pos_y += epj[k].mass * epj[k].pos.y;
				pos_z += epj[k].mass * epj[k].pos.z;
			}
			if(-1 == athread_get_id(-1)){
				printf("epj_beg=%d, nn=%d, epj_end=%d\n", epj_beg, nn, epj_end);
			}
			epj_beg += nn;
		}
                */
#endif
        if(mass > 0.0){
            double minv = 1.0 / mass;
            pos_x *= minv;
            pos_y *= minv;
            pos_z *= minv;
                /*
            #ifdef PHI_R_TREE
                pos_phi *= minv;
                pos_r   *= minv;
            #endif
                */
        }
    }
    else { // if(!isLeaf)
        //** Get tc (childnodes)
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, (unsigned *)&reply_get);
        dma_set_size(&dma_get, sizeof(etcLM) * N_CHILDREN);
        dma(dma_get, (long)((etcLM *)adr_tc_head + tc_tgt->adr_tc_), (long)(tc));
        dma_wait((unsigned *)&reply_get, 1);
        while (reply_get !=1 ) {}
        //** Inline expansion of tc_tmp->mom_.accumulate( tc_tmp_tmp->mom_ );
        for (k=0; k<N_CHILDREN; k++){
            if (tc[k].n_ptcl_ == 0) continue;
            mass  += tc[k].mass;
            pos_x += tc[k].mass * tc[k].pos.x;
            pos_y += tc[k].mass * tc[k].pos.y;
            pos_z += tc[k].mass * tc[k].pos.z;
                /*
            #ifdef PHI_R_TREE
                pos_phi += tc[k].mass * tc[k].pos_phi;
                pos_r   += tc[k].mass * tc[k].pos_r;
            #endif //PHI_R_TREE
                */
        } // for(k)
        if(mass > 0.0){
            double minv = 1.0 / mass;
            pos_x *= minv;
            pos_y *= minv;
            pos_z *= minv;
                /*
            #ifdef PHI_R_TREE
                pos_phi *= minv;
                pos_r   *= minv;
            #endif //PHI_R_TREE
                */
        }
    } // end if(isLeaf)
    tc_tgt[0].pos.x = pos_x;
    tc_tgt[0].pos.y = pos_y;
    tc_tgt[0].pos.z = pos_z;
        /*
    #ifdef PHI_R_TREE
        tc_tgt[0].pos_phi = pos_phi;
        tc_tgt[0].pos_r   = pos_r;
    #endif
        */
    tc_tgt[0].mass  = mass;
    *ptr_epj_cache_beg = epj_cache_beg;
    return;
}

void CalcMomentLean(void *args) {
    int my_id = athread_get_id(-1);
    dma_descriptor_init(&dma_get, (unsigned *)&reply_get);
    dma_descriptor_init(&dma_put, (unsigned *)&reply_put);

    reply_get = 0;
    long largs[6];
    dma_set_op(&dma_get, DMA_GET);
    dma_set_mode(&dma_get, PE_MODE);
    dma_set_reply(&dma_get, (unsigned *)&reply_get);
    dma_set_size(&dma_get, sizeof(largs));
    dma(dma_get, (long)(args), (long)(largs));
    dma_wait((unsigned *)&reply_get, 1);
    while (reply_get !=1 ) {}

    void *adr_epj_head = (void *)largs[0];
    void *adr_tc_head  = (void *)largs[1];
    int n_leaf_limit   = (int   )largs[2];
    int head           = (int   )largs[3];
    int next           = (int   )largs[4];
    void *adr_rsrch    = (void *)largs[5];

    double rsrch = *(double *)adr_rsrch;

    int num_tc,num_tc_loc,my_offset;
    num_tc = next - head;
    num_tc_loc = num_tc/NUMBER_OF_CPE + ( (my_id < num_tc % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (num_tc/NUMBER_OF_CPE)*my_id + ( (my_id < num_tc % NUMBER_OF_CPE) ? my_id : num_tc % NUMBER_OF_CPE );

    enum {
        CHUNK_SIZE = 256,
        //CHUNK_SIZE = 128,
        //CHUNK_SIZE = 64,
    };
    etcLM tc_tgt[CHUNK_SIZE];

    int i,j,k; 
    int epj_cache_beg = -9999;
    for (i=0; i<num_tc_loc; i+=CHUNK_SIZE) {
        int nrem = num_tc_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;

        //* Get the target tc
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, (unsigned *)&reply_get);
        dma_set_size(&dma_get, nn * sizeof(etcLM));
        dma(dma_get, (long)((etcLM *)adr_tc_head + head + my_offset + i), (long)(tc_tgt));
        dma_wait((unsigned *)&reply_get, 1);
        while (reply_get !=1 ) {}

        reply_put = 0;
        int nreq_put = 0; 

        for(j=0; j<nn; j++){
            if (tc_tgt[j].n_ptcl_ == 0) {
                continue;
            }
            UpdateTc(tc_tgt + j,
                     rsrch, n_leaf_limit, 
                     adr_tc_head, adr_epj_head,
                     &epj_cache_beg);
            
            //* Put the target tc (LM -> MM) (overwrite the target tc on MPE)
            // reply_put = 0;
            dma_set_op(&dma_put, DMA_PUT);
            dma_set_mode(&dma_put, PE_MODE);
            dma_set_reply(&dma_put, (unsigned *)&reply_put);
            dma_set_size(&dma_put, sizeof(etcLM));
            dma(dma_put, (long)((etcLM *)adr_tc_head + head + my_offset + i + j), (long)(tc_tgt + j));
            //dma_wait((unsigned *)&reply_put, 1);
            //while (reply_put !=1 ) {}
            nreq_put++;
        }
        while (reply_put != nreq_put ) {}
    }
}	 

__attribute__ ((noinline))
static void UpdateTcGlobal(
		etcLM *tc_tgt,
		double rsrch,
		int    n_leaf_limit,
		etcLM *adr_tc_head,
		tpLM  *adr_tp_head,
		epjLM *adr_epj_head,
		spjLM *adr_spj_head)
{
	enum {
		CHUNK_SIZE = 32,
	};
	tpLM  tp [CHUNK_SIZE];
	epjLM epj[CHUNK_SIZE];
	spjLM spj[CHUNK_SIZE];
	etcLM tc [N_CHILDREN];
	int i, j, k;


	double mass = 0.0;
	double pos_x = 0.0;
	double pos_y = 0.0;
	double pos_z = 0.0;
        /*
#ifdef PHI_R_TREE
        double pos_phi = 0.0;
        double pos_r   = 0.0;
#endif
        */
	if (tc_tgt[0].n_ptcl_ == 0) {
		return;
	}
	if (tc_tgt[0].n_ptcl_ <= n_leaf_limit || tc_tgt[0].level_ == TREE_LEVEL_LIMIT){
		for (j=0; j<tc_tgt[0].n_ptcl_; j+= CHUNK_SIZE) {
			int nrem = tc_tgt[0].n_ptcl_ - j;
			int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;

			//* Get tp
#if 1
			BEG_PROF(get_tp)
				reply_get = 0;
			dma_set_op(&dma_get, DMA_GET);
			dma_set_mode(&dma_get, PE_MODE);
			dma_set_reply(&dma_get, (unsigned *)&reply_get);
			dma_set_size(&dma_get, nn * sizeof(tpLM));
			dma(dma_get, (long)((tpLM *)adr_tp_head + tc_tgt[0].adr_ptcl_ + j), (long)(tp));
			dma_wait((unsigned *)&reply_get, 1);
			while (reply_get !=1 ) {}
			INC_CNT(get_tp, nn);
			END_PROF(get_tp)

    #ifdef REDUCE_MEMORY
                        U32_ offset_epj = 0;
                        U32_ offset_spj = 0;
                        int has_epj=0, has_spj=0;
			for (k=0; k<nn; k++){
				if ( ((tp[k].adr_ptcl_>>31) & 0x1) == 0 && has_epj == 0) {// tp[k].adr_ptcl_ is U32_ type.
                                    offset_epj = tp[k].adr_ptcl_;
                                    has_epj = 1;
				}
                                else if ( ((tp[k].adr_ptcl_>>31) & 0x1) == 1 && has_spj == 0) {
                                    offset_spj = (tp[k].adr_ptcl_ & 0x7fffffff);
                                    has_spj = 1;
				}
			}
    #else
			int has_epj=0, has_spj=0;
			for (k=0; k<nn; k++){
				if ( ((tp[k].adr_ptcl_>>31) & 0x1) == 0) {// tp[k].adr_ptcl_ is U32_ type.
					has_epj = 1;
				}else{
					has_spj = 1;
				}
			}
    #endif
			reply_get = 0;
			if(has_epj){
				//* Get epj
				BEG_PROF(get_epj)
				reply_get = 0;
				dma_set_op(&dma_get, DMA_GET);
				dma_set_mode(&dma_get, PE_MODE);
				dma_set_reply(&dma_get, (unsigned *)&reply_get);
				dma_set_size(&dma_get, nn * sizeof(epjLM));
    #ifdef REDUCE_MEMORY
                                dma(dma_get, (long*)((epjLM *)adr_epj_head + offset_epj + j), (long*)(epj));
    #else                                
				dma(dma_get, (long)((epjLM *)adr_epj_head + tc_tgt[0].adr_ptcl_ + j), (long)(epj));
    #endif
				while (reply_get != 1 ) {}
				INC_CNT(get_epj, nn);
				END_PROF(get_epj)
			}
			if(has_spj){
				//* Get spj
				BEG_PROF(get_spj)
				reply_get = 0;
				dma_set_op(&dma_get, DMA_GET);
				dma_set_mode(&dma_get, PE_MODE);
				dma_set_reply(&dma_get, (unsigned *)&reply_get);
				dma_set_size(&dma_get, nn * sizeof(spjLM));
    #ifdef REDUCE_MEMORY
                                dma(dma_get, (long)((spjLM *)adr_spj_head + offset_spj + j), (long)(spj));
    #else
				dma(dma_get, (long)((spjLM *)adr_spj_head + tc_tgt[0].adr_ptcl_ + j), (long)(spj));
    #endif
				while (reply_get != 1 ) {}
				INC_CNT(get_spj, nn);
				END_PROF(get_spj)
			}
#else
			reply_get = 0;

			dma_set_size(&dma_get, nn * sizeof(tpLM));
			dma(dma_get, (long)((tpLM *)adr_tp_head + tc_tgt[0].adr_ptcl_ + j), (long)(tp));

			dma_set_size(&dma_get, nn * sizeof(epjLM));
			dma(dma_get, (long)((epjLM *)adr_epj_head + tc_tgt[0].adr_ptcl_ + j), (long)(epj));

			while (reply_get != 2 ) {}

			int has_epj=0, has_spj=0;
			for (k=0; k<nn; k++){
				if ( ((tp[k].adr_ptcl_>>31) & 0x1) == 0) {// tp[k].adr_ptcl_ is U32_ type.
					has_epj = 1;
				}else{
					has_spj = 1;
				}
			}
			if(has_spj){
				dma_set_size(&dma_get, nn * sizeof(spjLM));
				dma(dma_get, (long)((spjLM *)adr_spj_head + tc_tgt[0].adr_ptcl_ + j), (long)(spj));
				while (reply_get != 3 ) {}
			}
#endif

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
					mass  += epj[k].mass;
					pos_x += epj[k].mass * epj[k].pos.x;
					pos_y += epj[k].mass * epj[k].pos.y;
					pos_z += epj[k].mass * epj[k].pos.z;
                                        /*
                                   #ifdef PHI_R_TREE
					pos_phi += epj[k].mass * epj[k].pos_phi;
					pos_r   += epj[k].mass * epj[k].pos_r;
                                   #endif
                                        */
				}
                                else {
#ifdef REDUCE_MEMORY
                                    k = n_cnt_sp;
                                    n_cnt_sp++;
#endif                                    
					mass  += spj[k].mass;
					pos_x += spj[k].mass * spj[k].pos.x;
					pos_y += spj[k].mass * spj[k].pos.y;
					pos_z += spj[k].mass * spj[k].pos.z;
                                        /*
                                   #ifdef PHI_R_TREE
					pos_phi += spj[k].mass * spj[k].pos_phi;
					pos_r   += spj[k].mass * spj[k].pos_r;
                                   #endif
                                        */
				}
			} // for(k)
		} // for(j)
		if(mass > 0.0){
			double minv = 1.0 / mass;
			pos_x *= minv;
			pos_y *= minv;
			pos_z *= minv;
                        /*
                   #ifdef PHI_R_TREE
                        pos_phi *= minv;
                        pos_r   *= minv;
                   #endif
                        */
		}
	}else{
		//** Get tc (childnodes)
		BEG_PROF(get_chl)
			reply_get = 0;
		dma_set_op(&dma_get, DMA_GET);
		dma_set_mode(&dma_get, PE_MODE);
		dma_set_reply(&dma_get, (unsigned *)&reply_get);
		dma_set_size(&dma_get, sizeof(etcLM) * N_CHILDREN);
		dma(dma_get, (long)((etcLM *)adr_tc_head + tc_tgt[0].adr_tc_), (long)(tc));
		dma_wait((unsigned *)&reply_get, 1);
		while (reply_get !=1 ) {}
		INC_CNT(get_chl, N_CHILDREN);
		END_PROF(get_chl)
			for (k=0; k<N_CHILDREN; k++){
				if (tc[k].n_ptcl_ == 0) continue;
				mass  += tc[k].mass;
				pos_x += tc[k].mass * tc[k].pos.x;
				pos_y += tc[k].mass * tc[k].pos.y;
				pos_z += tc[k].mass * tc[k].pos.z;
                                /*
                            #ifdef PHI_R_TREE
                                pos_phi += tc[k].mass * tc[k].pos_phi;
                                pos_r   += tc[k].mass * tc[k].pos_r;
                            #endif
                                */
			}
		if(mass > 0.0){
			double minv = 1.0 / mass;
			pos_x *= minv;
			pos_y *= minv;
			pos_z *= minv;
                        /*
                   #ifdef PHI_R_TREE
                        pos_phi *= minv;
                        pos_r   *= minv;
                   #endif 
                        */                       
		}
	} // end if(isLeaf)
	tc_tgt[0].pos.x = pos_x;
	tc_tgt[0].pos.y = pos_y;
	tc_tgt[0].pos.z = pos_z;
	tc_tgt[0].mass  = mass;
        /*
    #ifdef PHI_R_TREE
        tc_tgt[0].pos_phi = pos_phi;
        tc_tgt[0].pos_r   = pos_r;
    #endif
        */
	return;
}

void CalcMomentLongGlobalTreeLean(void *args) {
	dma_descriptor_init(&dma_get, (unsigned *)&reply_get);
	dma_descriptor_init(&dma_put, (unsigned *)&reply_put);
	init_prof();
	int my_id = athread_get_id(-1);
	{
		int i;
		sync_array_();
		for(i=0; i<100*my_id; i++) {
			asm volatile ("nop");
		}
	}

	reply_get = 0;
	long largs[8];
	dma_set_op(&dma_get, DMA_GET);
	dma_set_mode(&dma_get, PE_MODE);
	dma_set_reply(&dma_get, (unsigned *)&reply_get);
	dma_set_size(&dma_get, sizeof(largs));
	dma(dma_get, (long)(args), (long)(largs));
	dma_wait((unsigned *)&reply_get, 1);
	while (reply_get !=1 ) {}

	sync_array_();
	long global_begin = rtc_();

	void *adr_tc_head  = (void *)largs[0];
	void *adr_tp_head  = (void *)largs[1];
	void *adr_epj_head = (void *)largs[2];
	void *adr_spj_head = (void *)largs[3];
	int n_leaf_limit   = (int   )largs[4];
	int head           = (int   )largs[5];
	int next           = (int   )largs[6];
	//void * adr_rsrch   = (void *)largs[7];
	double rsrch   = *(double *)&largs[7];

	//* Compute the task of each CPE
	int num_tc,num_tc_loc,my_offset;
	num_tc = next - head;
	num_tc_loc = num_tc/NUMBER_OF_CPE + ( (my_id < num_tc % NUMBER_OF_CPE) ? 1 : 0 );
	my_offset = (num_tc/NUMBER_OF_CPE)*my_id + ( (my_id < num_tc % NUMBER_OF_CPE) ? my_id : num_tc % NUMBER_OF_CPE );

	enum {
		CHUNK_SIZE = 256,
	};
	etcLM tc_tgt[CHUNK_SIZE];

	int i, j; 
	for (i=0; i<num_tc_loc; i+=CHUNK_SIZE) {
		int nrem = num_tc_loc - i;
		int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;

		//* Get the target tc
		BEG_PROF(get_tc)
			reply_get = 0;
		dma_set_op(&dma_get, DMA_GET);
		dma_set_mode(&dma_get, PE_MODE);
		dma_set_reply(&dma_get, (unsigned *)&reply_get);
		dma_set_size(&dma_get, nn * sizeof(etcLM));
		dma(dma_get, (long)((etcLM *)adr_tc_head + head + my_offset + i), (long)(tc_tgt));
		dma_wait((unsigned *)&reply_get, 1);
		while (reply_get !=1 ) {}
		INC_CNT(get_tc, nn);
		END_PROF(get_tc)

		reply_put = 0;
		int nreq_put = 0; 
		for(j=0; j<nn; j++){
			if (tc_tgt[j].n_ptcl_ == 0){
				continue;
			}
			UpdateTcGlobal(
					tc_tgt+j,
					rsrch, n_leaf_limit, 
					adr_tc_head, adr_tp_head, adr_epj_head, adr_spj_head);

			BEG_PROF(put_tc)
			dma_set_op(&dma_put, DMA_PUT);
			dma_set_mode(&dma_put, PE_MODE);
			dma_set_reply(&dma_put, (unsigned *)&reply_put);
			dma_set_size(&dma_put, sizeof(etcLM));
			dma(dma_put, (long)((etcLM *)adr_tc_head + head + my_offset + i + j), (long)(tc_tgt + j));
			// dma_wait((unsigned *)&reply_put, 1);
			// while (reply_put !=1 ) {}
			INC_CNT(put_tc, 1);
			END_PROF(put_tc)
			nreq_put++;
		}
		BEG_PROF(put_tc)
		while (reply_put != nreq_put ) {}
		END_PROF(put_tc)
	}

	sync_array_();
	double dt  = (rtc_() - global_begin) / 1.45e9;
	print_bandwidth(dt);
	fflush(stdout);
	sync_array_();

	if(0 == my_id){
		print_prof();
	}
}
