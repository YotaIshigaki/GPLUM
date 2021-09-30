#include"slave.h"
#include"simd.h"
#include"dma.h"
#include"cpe_prof.h"
#include"cpe_func.h"

#define ID_CHECK 63
#define REP8(A) A,A,A,A,A,A,A,A

#ifdef FULL_ASM
#define EPJ_COST 39  //54*1.2
#define SPJ_COST 30  //29*1.2
#define SAT_COST 50 //100*1.2
#else
#define EPJ_COST 54  //54*1.2
#define SPJ_COST 30  //29*1.2
#define SAT_COST 100 //100*1.2
#endif
#define PLANET_COST 32 //128/4
#define CONST_COST 237404 //((128*(8192*89+128*100)+128)/4*0.01
//#define CONST_COST 0

#define NW_LOCAL_MAX 512
#define NW_CPE_MAX 150

//void *ldm_malloc(size_t);
//void ldm_free(void *, size_t);
void abort();

static void final_insertion_sort_ul_decrease(unsigned long a[], int n) {
	int i, j;
	for (i=1; i<n; i++) {
		unsigned long  tmp = a[i];
		if (a[i - 1] < tmp) {
			j = i;
			do {
				a[j] = a[j-1];
				j--;
			} while (j>0 && a[j-1]<tmp);
			a[j] = tmp;
		}
	}
}

static unsigned long qsort_med3_ul(unsigned long x, unsigned long y, unsigned long z){
	return x<y ? y<z ? y
	                 : x<z ? z : x
	           : x<z ? x
	                 : y<z ? z : y;
}

__attribute__((noinline))
static void my_quicksort_ul_decrease(
		unsigned long a[], 
		int beg, // inclusive
		int end) // exclusive
{
	const int NSTACK = 64;
	int sp = 0;
	int stack[NSTACK][2];
//sort_body:
  while(1){
	int n = end - beg;
	if(n <= 8){
		final_insertion_sort_ul_decrease(a+beg, n);
	}else{
		int i = beg;
		int j = end - 1;
        int count=0;
        unsigned long piv = a[i+(j-i)/2]; 
		// T piv = a[(i+j)/2];
        //		unsigned long piv = qsort_med3_ul(a[i], a[i+(j-i)/2], a[j]);

		if(sp > 30) printf("sort [%d:%d), piv at %d, sp=%d, N=%d\n", beg, end, (i+j)/2, sp, end-beg);
bsplit:
		for(;;){
			while(a[i] > piv) ++i;
			while(a[j] < piv) --j;
			if(i >= j) break;
			unsigned long ai = a[i];
			unsigned long aj = a[j];
			a[i] = aj;
			a[j] = ai;
			i++;
			j--;
		}
        int right=end-(j+1);
        int left=i-beg;
        if(left==0||right==0) {
          if(count) {
            printf("failed \n");
            abort();
          }
          i=beg;
          j=end-1;
          piv = qsort_med3_ul(a[i], a[i+(j-i)/2], a[j]);
          count++;
          goto bsplit;
        }
		if( j+1 < end ){
			stack[sp][0] = j+1;
			stack[sp][1] = end;
			sp++;
		}
		if( beg < i ){
			end = i;
			// goto sort_body;
			continue;
		}
	}
	if(sp){
		sp--;
        //		assert(sp < NSTACK);
		beg = stack[sp][0];
		end = stack[sp][1];
		//goto sort_body;
		continue;
	}
	break;
  }
}


// root depth = 0, maximum num at level of depth = 2^depth
static inline int find_min_tree(unsigned long* trb, const int depth){
  int i,k=0,sf=1;
  int csf=sf;
  for(i=0; i<depth; i++) {
    k *= 2;
    int ksf = k+csf;
    k += trb[ksf]>trb[ksf+1];
    sf  *= 2;
    csf += sf;
  }
  return k;
}

// k range: (0,length-1)
static inline void insert_sum_tree(unsigned long* tree, const int k, const unsigned long val, const int length) {
  int sf=1;
  int csf=0;
  int i;
  for(i=length; i>0; i/=2) {
    tree[csf+k/i] += val;
    csf += sf;
    sf *= 2;
  }
};

// k range: (0,length-1)
static inline void insert_min_tree(unsigned long* tree, const int k, const unsigned long val, const int length) {
  int i,j=k;
  tree[length-1+k] += val;
  for(i=length; i>1; i/=2, j/=2) {
    //tree[csf+k/i] = val;
    int sj2 = i-1+((j%2==0)?(j+1):(j-1));
    int sj1 = i-1+j;
    tree[(i+j)/2-1] = (tree[sj2]>tree[sj1])?tree[sj1]:tree[sj2];
  }
};

void print_tree(unsigned long* tree, const int depth, const int mask) {
  int sf=1;
  int csf=0;
  int i;
  for(i=0; i<=depth; i++) {
    if(!mask||i==depth) {
      printf("CID           min");
      int k;
      for(k=0;k<sf;k++) printf("%12d",k);
      printf("\n");
      printf("L\%2d",i);
      unsigned long sum=0;
      for(k=0;k<sf;k++) sum+=tree[csf+k];
      printf("%14lu",sum);
      //unsigned long cmin = 9223372036854775807;
      //for(k=0;k<sf;k++) cmin = (cmin>tree[csf+k])?tree[csf+k]:cmin;
      //printf("%14lu",cmin);
      for(k=0;k<sf;k++) printf("%12lu",tree[csf+k]);
      printf("\n");
    }
    csf += sf;
    sf  *= 2;
  }
}

void balance_nwalk(void* arg) {
#ifdef PROFILE
  unsigned long t_count=rtc_();
#endif

//  unsigned int my_row,my_col,my_id;
//  get_col_id_(&my_col);
//  get_row_id_(&my_row);
//  my_id = my_row*8 + (my_row%2==0)?my_col:(7-my_col);
  

  Force_Kernel_Pars pars;
  volatile unsigned int par_reply = 0;
  athread_get(PE_MODE,arg,&pars,sizeof(Force_Kernel_Pars),(void*)&par_reply,0,0,0);
  while(par_reply!=1);

  if(pars.n_walk<=64) {
    //    if(my_id==0) {
      // int nw = (pars.n_walk+7)/8*8;
      int* nindx=(int*)(long)ldm_malloc(64*sizeof(int));
      int* nwcpe=(int*)(long)ldm_malloc(64*sizeof(int));
      int* ndspw=(int*)(long)ldm_malloc(64*sizeof(int));
      //      intv8 nifxt  = simd_set_intv8(0,1,2,3,4,5,6,7);
      //      intv8 nnext  = {REP8(8)};
      //      intv8 i1v8   = {REP8(1)};
      int i;
      for(i=0; i<pars.n_walk; i++) {
        nindx[i] = i;
        ndspw[i] = i;
        nwcpe[i] = 1;
        //        simd_store(nifxt, &nindx[i]);
        //        simd_store(nifxt, &ndspw[i]);
        //        simd_print_intv8(nnext);
        //        simd_print_intv8(nifxt);
        //        nifxt = nifxt + nnext;
        //        simd_print_intv8(nnext);
        //        simd_print_intv8(nifxt);
        //        printf("\n");
        //        simd_store(i1v8,  &nwcpe[i]);
      }
      for(i=pars.n_walk; i<64; i++) nwcpe[i] = 0;
      
      par_reply=0;
      athread_put(PE_MODE, nindx, pars.adr_n_walk,  64*sizeof(int), (void*)&par_reply, 0, 0);
      athread_put(PE_MODE, nwcpe, pars.n_walk_cpe,  64*sizeof(int), (void*)&par_reply, 0, 0);
      athread_put(PE_MODE, ndspw, pars.n_disp_walk, 64*sizeof(int), (void*)&par_reply, 0, 0);
      while(par_reply!=3);
      
      ldm_free(nindx,64*sizeof(int));
      ldm_free(nwcpe,64*sizeof(int));
      ldm_free(ndspw,64*sizeof(int));

      //    }
    return;
  }

  int nw_div=(pars.n_walk+NW_LOCAL_MAX-1)/NW_LOCAL_MAX;
  //  if (my_id<nw_div) return;
  int nlst_size = NW_LOCAL_MAX*sizeof(int);

  int* ni=(int*)(long)ldm_malloc(nlst_size);
  int* nj=(int*)(long)ldm_malloc(nlst_size);
  int* ns=(int*)(long)ldm_malloc(nlst_size);
  unsigned long* ncost=(unsigned long*)(long)ldm_malloc(NW_LOCAL_MAX*sizeof(unsigned long));
//  int nj[NW_LOCAL_MAX];
//  int ns[NW_LOCAL_MAX];
//  unsigned long ncost[NW_LOCAL_MAX];

  unsigned long* csum_tree=(unsigned long*)(long)ldm_malloc(128*sizeof(unsigned long));
  // unsigned long csum_tree[128]; // one more for alignment
  int* nw_cpe        =(int*)(long)ldm_malloc(64*sizeof(int));
  int* nw_cpe_cum    =(int*)(long)ldm_malloc(64*sizeof(int));
  int* nidx_cpe_data =(int*)(long)ldm_malloc(64*NW_CPE_MAX*sizeof(int));
  // int nw_cpe[64];
  // int nw_cpe_cum[64];
  // int nw_cpe_tot[64];
  int* nidx_cpe[64];
  int i;
  for(i=0; i<64; i++) nidx_cpe[i] = &nidx_cpe_data[NW_CPE_MAX*i];

#ifdef M_CHECK
  printf("Address position: ni=%ld nj=%ld ns=%ld ncost=%ld csum_tree=%ld nw_cpe=%ld nw_cpe_cum=%ld nidx_cpe[0]=%ld nidx_cpe[64]=%ld\n",(long)ni,(long)nj,(long)ns,(long)ncost,(long)csum_tree,(long)nw_cpe,(long)nw_cpe_cum,(long)nidx_cpe[0],(long)&nidx_cpe[63][NW_CPE_MAX]);
#endif

  // volatile intv8 i7v8 = {REP8(7)};
  // volatile intv8 i8v8 = {REP8(8)};
  
  int nsat = (pars.n_sat+3)/4;
  //intv8 nsatv8 = {REP8(nsat)};
  //intv8 epcost = {REP8(EPJ_COST)};
  //intv8 spcost = {REP8(SPJ_COST)};
  //intv8 satcost= {REP8(SAT_COST)};
  //intv8 pcost  = {REP8(PLANET_COST)};
  //intv8 ccost  = {REP8(CONST_COST)};

  dma_desc n_get;
  volatile unsigned int get_reply;
  dma_descriptor_init((dma_desc*)&n_get,(unsigned int*)&get_reply);
  dma_set_op(   (dma_desc*)&n_get, DMA_GET);
  dma_set_mode( (dma_desc*)&n_get, PE_MODE);
  dma_set_reply((dma_desc*)&n_get, &get_reply);
  
  int lk;
  for(lk=0; lk<nw_div; lk++) {
    int nw_lm = NW_LOCAL_MAX;
    int nw_shift = lk*NW_LOCAL_MAX;
    if (lk==nw_div-1) nw_lm = pars.n_walk%NW_LOCAL_MAX;
    if (nw_lm==0) nw_lm = NW_LOCAL_MAX;
  
    // int nw4 = ((nw_lm+7)/8*8);
    // int nlst_size = nw4*sizeof(int);
  
    //  int* ni = (int*)(long)ldm_malloc(nlst_size);
    //  int* nj = (int*)(long)ldm_malloc(nlst_size);
    //  int* ns = (int*)(long)ldm_malloc(nlst_size);

    int nw_start_index=lk*NW_LOCAL_MAX;
  
    dma_set_size( (dma_desc*)&n_get, nw_lm*sizeof(int));

    get_reply = 0;
    dma(n_get, (long)&pars.n_epi[nw_start_index], (long)ni); 
    dma(n_get, (long)&pars.n_epj[nw_start_index], (long)nj); 
    dma(n_get, (long)&pars.n_spj[nw_start_index], (long)ns); 
    while(get_reply<3);

#ifdef BW_DEBUG
    printf("Get ni[0]=%d,nj[0]=%d,ns[0]=%d\n",ni[0],nj[0],ns[0]);
#endif    
    // int* ncost=(int*)(long)ldm_malloc(nlst_size);
    // int* nindx=(int*)(long)ldm_malloc(nlst_size);

    //  intv8 nifxt  = simd_set_intv8(0,1,2,3,4,5,6,7);
    //  intv8 nnext  = REP8(8);

    int i;
//    for(i=0; i<nw_lm; i+=8) {
//      intv8 ki = *(intv8*)&ni[i];
//      intv8 kj = *(intv8*)&nj[i];
//      intv8 ks = *(intv8*)&ns[i];
//    
//      intv8 kcost = ((kj+i7v8)/i8v8 * epcost + 
//                     (ks+i7v8)/i8v8 * spcost + 
//                     nsatv8 * satcost + 
//                     pcost) * (ki+i7v8)/i8v8 + ccost;
//      simd_store(kcost, &ncost[i]);
//      //    simd_store(nifxt, &nindx[i]);
//      //    nifxt += nnext;
//    }
    for(i=0; i<nw_lm; i++) {
      ncost[i] = ((nj[i]+3)/4 * (unsigned long)EPJ_COST + 
                  (ns[i]+3)/4 * (unsigned long)SPJ_COST +
                  nsat        * (unsigned long)SAT_COST +
                  PLANET_COST)* (ni[i]+3)/4*4 + CONST_COST;
#ifdef BW_CHECK
      if((ncost[i]<<16>>16)!=ncost[i]) {
        printf("ncost[%d] = %lu overflow!\n",i,ncost[i]);
        abort();
      }
#endif
      ncost[i] = (ncost[i]<<16)|i;

    }

    //sorting
    my_quicksort_ul_decrease(ncost,0,nw_lm);

#ifdef BW_DEBUG
    for(i=0; i<nw_lm-1; i++) {
      if(ncost[i]<ncost[i+1]) {
        printf("Sorting error, loop=%d, i=%d,ncost[i]=%lu,ncost[i+1]=%lu\n",lk,i,ncost[i],ncost[i+1]);
        abort();
      }
    }
#endif
    
    int istart = 0;
    if(lk==0) {
      //intv8 i0v8 = {REP8(0)};
      for(i=0; i<64;  i++) nw_cpe[i] = 0;
        //simd_store(i0v8, &nw_cpe[i]);
      for(i=0; i<128; i++) csum_tree[i] = 0;
        //simd_store(i0v8, (int*)&csum_tree[i]);
#ifdef BW_DEBUG
      for(i=0; i<64; i++) {
        if(nw_cpe[i]!=0) printf("Error: nw_cpe[%d]=%d\n",i,nw_cpe[i]);
      }
      for(i=0; i<128; i++) {
        if(csum_tree[i]!=0) printf("Error: csum_tree[%d]=%d\n",i,csum_tree[i]);
      }
#endif
      // for(i=0; i<64; i+=8) simd_store(i0v8, &nw_cpe_tot[i]);
      // first 64
      for(i=0; i<64; i++) {
        nw_cpe[i]++;
        nidx_cpe[i][0]= (int)(ncost[i]&0xffff);
        insert_min_tree(csum_tree,i,(ncost[i]>>16),64);
      }
      istart = 64;
    }
//    else if(my_col!=7&&my_col!=0){
//    for(i=0; i<64; i+=8) {
//      REG_GETR(csum_tree);
//    }
//  }
  
    for(i=istart; i<nw_lm; i++) {
#ifdef BW_DEBUG_DEEP
      print_tree(csum_tree,6,1);
#endif
      int imin= find_min_tree(csum_tree,6);
      insert_min_tree(csum_tree,imin,(ncost[i]>>16),64);
      nidx_cpe[imin][nw_cpe[imin]]= (int)(ncost[i]&0xffff) + nw_shift;
      nw_cpe[imin]++;
#ifdef BW_DEBUG_DEEP
      printf("Insert i=%d, iw=%d, imin=%d, cost=%lu, nw[imin]=%d\n",i,(int)(ncost[i]&0xffff)+nw_shift,imin,(ncost[i]>>16),nw_cpe[imin]);
#endif
      if(nw_cpe[imin]>NW_CPE_MAX) {
        printf("Error! n_walk per cpe overflow! CPE ID=%d, maximum n_walk=%d, current nwalk index=%d, total n_walk=%d\n", imin, NW_CPE_MAX, i, nw_lm);
        abort();
      }
    }
  }

  //  nw_cpe_cum[0] = nw_cpe[0];
  nw_cpe_cum[0] = 0;
  for(lk=1; lk<64; lk++) nw_cpe_cum[lk] = nw_cpe_cum[lk-1] + nw_cpe[lk-1];
#ifdef BW_CHECK
  if(nw_cpe_cum[63]+nw_cpe[63]!=pars.n_walk) {
    printf("Error! nw_disp_walk not consistent with n_walk, nw_disp_walk[63]=%d nw_cpe[63]=%d n_walk=%d\n",nw_cpe_cum[63],nw_cpe[63],pars.n_walk);
    abort();
  }
#endif
  
  dma_desc n_put;
  volatile unsigned int put_reply;
  dma_descriptor_init((dma_desc*)&n_put,(unsigned int*)&put_reply);
  dma_set_op(   (dma_desc*)&n_put, DMA_PUT);
  dma_set_mode( (dma_desc*)&n_put, PE_MODE);
  dma_set_reply((dma_desc*)&n_put, &put_reply);
  
  put_reply=0;
  dma_set_size( (dma_desc*)&n_put, 64*sizeof(int));
  dma(n_put, (long)pars.n_walk_cpe,  (long)nw_cpe);
  dma(n_put, (long)pars.n_disp_walk, (long)nw_cpe_cum);
  while(put_reply<2);

  put_reply=0;
  for(lk=0; lk<64; lk++) {
    dma_set_size((dma_desc*)&n_put, nw_cpe[lk]*sizeof(int));
    dma(n_put, (long)&pars.adr_n_walk[nw_cpe_cum[lk]], (long)nidx_cpe[lk]);
    while(put_reply<=lk);
  }

#ifdef BW_CHECK
  print_tree(csum_tree,6,0);
//  printf("CPEID:");
//  for(lk=0; lk<64; lk++) printf("\t%d",lk);
//  printf("\nWORK: ");
//  for(lk=0; lk<64; lk++) printf("\t%d",csum_tree[63+lk]);
//  printf("\n");
  printf("last adr walk: n=%d\n",nw_cpe[63]);
  for(lk=0; lk<nw_cpe[63]; lk++) printf("\t%d",nidx_cpe[63][lk]);
  printf("\n");
#endif

  ldm_free(ni,nlst_size);
  ldm_free(nj,nlst_size);
  ldm_free(ns,nlst_size);
  ldm_free(ncost,NW_LOCAL_MAX*sizeof(unsigned long));
  ldm_free(csum_tree,   128*sizeof(unsigned long));
  ldm_free(nw_cpe,       64*sizeof(int));
  ldm_free(nw_cpe_cum,   64*sizeof(int));
  ldm_free(nidx_cpe_data,64*NW_CPE_MAX*sizeof(int));

#ifdef PROFILE
  unsigned long t_total = rtc_() - t_count;
  //  if(my_id==ID_CHECK) {
  printf("Balance_nwalk time=%e, cycles=%lu\n", (double)t_total/1.45e9, t_total);
    //  }
#endif

  return;
};
