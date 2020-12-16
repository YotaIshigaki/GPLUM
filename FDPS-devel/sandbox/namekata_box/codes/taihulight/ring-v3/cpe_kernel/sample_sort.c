#include <slave.h>
#include <simd.h>
#include <dma.h>
#include <assert.h>
#include "cpe_func.h"
//#include "common.h"



void sync_array_(){
	unsigned long sync_tmp=-1;
	asm volatile(
			"ldi %0, 0xff\n"
			"sync %0\n"
			"synr %0\n"
			:"=r"(sync_tmp)
			:
			:"memory");
}

void *ldm_malloc(size_t);
void ldm_free(void *, size_t);
void abort();

typedef struct KeyIndex {
    unsigned long key_hi_;
#ifdef USE_96BIT_KEY
    unsigned int key_lo_;
    unsigned int index;
#else
    unsigned int index;
    unsigned int pad;
#endif
} KeyIndex;

// 
extern unsigned long glb_parts[NSPLIT-1];

extern int glb_task_off[NSPLIT+1];
extern int glb_status[NSPLIT];

extern int glb_count[NTHREADS][NSPLIT];

extern int glb_enn;
extern void *glb_ptr_inp;
extern void *glb_ptr_work;
extern void *glb_ptr_out;
// 

static inline int comp_lt_ki(const KeyIndex * a, const KeyIndex * b){
#ifdef USE_96BIT_KEY
    return (a->key_hi_ != b->key_hi_) ? a->key_hi_ < b->key_hi_ : a->key_lo_ < b->key_lo_;
#else
    return a->key_hi_ < b->key_hi_;
#endif
}

static int where_it_goes_ki(const KeyIndex * val,
                            const KeyIndex parts[]){
    int loc = NSPLIT/2;
    int inc;
    for(inc = NSPLIT/4; inc; inc/=2){
        const KeyIndex piv = parts[loc-1];
        if(comp_lt_ki(val, &piv)){
            loc -= inc;
        }else{
            loc += inc;
        }
    }
    const KeyIndex piv = parts[loc-1];
    if(comp_lt_ki(val, &piv)){
        loc--;
    }
    return loc;
}

static void final_insertion_sort_ki(KeyIndex a[], int n) {
	int i, j;
	// Make it partially sorted
	for(i=0, j=n/2+n%2; i<n/2; i++, j++){
#if 1
		KeyIndex ai = a[i];
		KeyIndex aj = a[j];

		KeyIndex bi = comp_lt_ki(&ai, &aj) ? ai : aj;
		KeyIndex bj = comp_lt_ki(&ai, &aj) ? aj : ai;
		//KeyIndex bi = (ai.key_hi_ < aj.key_hi_) ? ai : aj;
		//KeyIndex bj = (ai.key_hi_ < aj.key_hi_) ? aj : ai;
                
		a[i] = bi;
		a[j] = bj;
#else
                // under construction
                // TODO: swap key_lo_
                unsigned long keyi = a[i].key_hi_;
                unsigned long keyj = a[j].key_hi_;
                unsigned int  idxi = a[i].index;
                unsigned int  idxj = a[j].index;
                
		unsigned long keyI = (keyi < keyj) ? keyi : keyj;
		unsigned long keyJ = (keyi < keyj) ? keyj : keyi;
		unsigned int idxI = (keyi < keyj) ? idxi : idxj;
		unsigned int idxJ = (keyi < keyj) ? idxj : idxi;

                a[i].key_hi_   = keyI;
		a[i].index = idxI;
                a[j].key_hi_   = keyJ;
		a[j].index = idxJ;

#endif
	}

	for (i=1; i<n; i++) {
#if 1
            KeyIndex tmp = a[i];
            if (comp_lt_ki(&tmp, &a[i-1])) {
                j = i;
                do {
                    a[j] = a[j-1];
                    j--;
                } while (j>0 && comp_lt_ki(&tmp, &a[j-1]));
                a[j] = tmp;
            }
                /*
#ifdef USE_96BIT_KEY
		if (a[i - 1].key_hi_ > tmp.key_hi_ || ( (a[i-1].key_hi_==tmp.key_hi_) && (a[i-1].key_lo_>tmp.key_lo_) ) ) {
			j = i;
			do {
				a[j] = a[j-1];
				j--;
			} while (j>0 && ((a[j-1].key_hi_ > tmp.key_hi_) || ( (a[i-1].key_hi_==tmp.key_hi_) && (a[i-1].key_lo_>tmp.key_lo_) ) ));
			a[j] = tmp;
		}
#else
		if (a[i - 1].key_hi_ > tmp.key_hi_) {
			j = i;
			do {
				a[j] = a[j-1];
				j--;
			} while (j>0 && a[j-1].key_hi_ > tmp.key_hi_);
			a[j] = tmp;
		}
#endif
                */
#else

                unsigned long tmpkey = a[i].key_hi_;
                unsigned int  tmpidx = a[i].index;
		if (a[i - 1].key_hi_ > tmpkey) {
			j = i;
			do {
				a[j] = a[j-1];
				j--;
			} while (j>0 && a[j-1].key_hi_ > tmpkey);
			a[j].key_hi_ = tmpkey;
			a[j].index = tmpidx;
		}
#endif
	}
}

static unsigned long qsort_med3_ul(unsigned long x, unsigned long y, unsigned long z){
	return x<y ? y<z ? y
	                 : x<z ? z : x
	           : x<z ? x
	                 : y<z ? z : y;
}

static unsigned int qsort_med3_ui(unsigned int x, unsigned int y, unsigned int z){
	return x<y ? y<z ? y
	                 : x<z ? z : x
	           : x<z ? x
	                 : y<z ? z : y;
}

static KeyIndex qsort_med3_ki(KeyIndex * x, KeyIndex * y, KeyIndex * z){
    return comp_lt_ki(x,y) ? comp_lt_ki(y,z) ? *y
        : comp_lt_ki(x,z) ? *z : *x
        : comp_lt_ki(x,z) ? *x
        : comp_lt_ki(y,z) ? *z : *y;
}


__attribute__((noinline))
static void my_quicksort_ki(
		KeyIndex a[], 
		int beg, // inclusive
		int end) // exclusive
{
	const int NSTACK = 64;
	int sp = 0;
	int stack[NSTACK][2];
//sort_body:
  while(1){
	int n = end - beg;
	if(n <= 32){
		final_insertion_sort_ki(a+beg, n);
	}else{
		int i = beg;
		int j = end - 1;
                int count=0;
                //unsigned long piv = qsort_med3_ul(a[i].key_hi_, a[(i+j)/2].key_hi_, a[j].key_hi_);
                KeyIndex piv = qsort_med3_ki(&a[i], &a[(i+j)/2], &a[j]);
		if(sp > 30) printf("sort [%d:%d), piv at %d, sp=%d, N=%d\n", beg, end, (i+j)/2, sp, end-beg);
split:                
		for(;;){
                    //while(a[i].key_hi_ < piv) ++i;
                    //while(a[j].key_hi_ > piv) --j;
                    while(comp_lt_ki(&a[i], &piv)) ++i;
                    while(comp_lt_ki(&piv, &a[j])) --j;
			if(i >= j) break;
			KeyIndex ai = a[i];
			KeyIndex aj = a[j];
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
                    piv = qsort_med3_ki(&a[i], &a[i+(j-i)/2], &a[j]);
                    /*
#ifdef USE_96BIT_KEY
                if( a[i].key_hi_==a[i+(j-i)/2].key_hi_ && a[i].key_hi_==a[j].key_hi_){
                    piv = qsort_med3_ui(a[i].key_lo_, a[i+(j-i)/2].key_lo_, a[j].key_lo_);
                }
                else{
                     piv = qsort_med3_ul(a[i].key_hi_, a[i+(j-i)/2].key_hi_, a[j].key_hi_);
                }
#else
                piv = qsort_med3_ul(a[i].key_hi_, a[i+(j-i)/2].key_hi_, a[j].key_hi_);
#endif
                    */
                    count++;
                    goto split;
                }
#if 0
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
#else
        int lbeg = beg, lend = i, rbeg = j+1, rend = end; 
        if (right < left) {
            // Swap (lbeg,lend) and (rbeg,rend) to store
            // the larger one to the stack.
            int tmp = rbeg;
            rbeg = lbeg;
            lbeg = tmp;
            tmp = rend;
            rend = lend;
            lend = tmp;
        }
        if( rbeg < rend ){
            stack[sp][0] = rbeg;
            stack[sp][1] = rend;
            sp++;
            if (sp > NSTACK) {
                printf("id= %d sp= %d piv.key_hi_ = %lu left = %d right = %d n = %d\n",
                       athread_get_id(-1), sp, piv.key_hi_, left, right, n);
            }
            assert(sp<=NSTACK);
        }
        if( lbeg < lend ){
            beg = lbeg;
            end = lend;
            // goto sort_body;
            continue;
        }                
#endif
	}
	if(sp){
		sp--;
		assert(sp < NSTACK);
		beg = stack[sp][0];
		end = stack[sp][1];
		//goto sort_body;
		continue;
	}
	break;
  }
}
 
static void final_insertion_sort_ul(unsigned long a[], int n) {
    int i, j;
    for (i=1; i<n; i++) {
        unsigned long  tmp = a[i];
        if (a[i - 1] > tmp) {
            j = i;
            do {
                a[j] = a[j-1];
                j--;
            } while (j>0 && a[j-1]>tmp);
            a[j] = tmp;
        }
    }
}

__attribute__((noinline))
static void my_quicksort_ul(
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
		final_insertion_sort_ul(a+beg, n);
	}else{
		int i = beg;
		int j = end - 1;
                int count=0;
		// T piv = a[(i+j)/2];
		unsigned long piv = qsort_med3_ul(a[i], a[(i+j)/2], a[j]);

		if(sp > 30) printf("sort [%d:%d), piv at %d, sp=%d, N=%d\n", beg, end, (i+j)/2, sp, end-beg);
                
split2:
		for(;;){
			while(a[i] < piv) ++i;
			while(a[j] > piv) --j;
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
            if (left==0||right==0) {
                if (count) {
                    printf("failed \n");
                    abort();
                }
                i=beg;
                j=end-1;
                piv = qsort_med3_ul(a[i], a[i+(j-i)/2], a[j]);
                count++;
                goto split2;
            }
#if 0
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
#else
            int lbeg = beg, lend = i, rbeg = j+1, rend = end; 
            if (right < left) {
                // Swap (lbeg,lend) and (rbeg,rend) to store
                // the larger one to the stack.
                int tmp = rbeg;
                rbeg = lbeg;
                lbeg = tmp;
                tmp = rend;
                rend = lend;
                lend = tmp;
            }
            if( rbeg < rend ){
                stack[sp][0] = rbeg;
                stack[sp][1] = rend;
                sp++;
                if (sp>NSTACK) {
                    printf("id= %d sp= %d piv = %lu left = %d right = %d n = %d\n",
                           athread_get_id(-1),sp,piv, left, right, n);
                }
                assert(sp<=NSTACK);
            }
            if ( lbeg < lend ){
                beg = lbeg;
                end = lend;
                // goto sort_body;
                continue;
            }
#endif
	}
	if(sp){
		sp--;
		assert(sp < NSTACK);
		beg = stack[sp][0];
		end = stack[sp][1];
		//goto sort_body;
		continue;
	}
	break;
  }
}

static inline int where_it_goes(unsigned long val,
                                const unsigned long parts[])
{
	int loc = NSPLIT/2;
	int inc;
	for(inc = NSPLIT/4; inc; inc/=2){
		const unsigned long piv = parts[loc-1];
		if(val < piv){
			loc -= inc;
		}else{
			loc += inc;
		}
	}
	const unsigned long piv = parts[loc-1];
	if(val < piv){
		loc--;
	}
	return loc;
}

__thread_local volatile unsigned long get_reply, put_reply;
// __thread_local unsigned long parts[NSPLIT-1];

static void inner_count_kernel(
		const KeyIndex      *inp,
		const int            n,
		const KeyIndex *parts,
		int                 *count)
{
	int i, ii;
#if 0
	for(i=0; i<n; i++){
		const int loc = where_it_goes(inp[i].key, parts);
		count[loc]++;
	}
#else
	KeyIndex buf[2048]; // 32 KiB in stack
	for(ii=0; ii<n; ii+=2048){
		int rem = n-ii;
		int size = rem<2048 ? rem : 2048;
		get_reply = 0;
		athread_get(PE_MODE, (void *)(inp+ii), buf, size*16, (unsigned long *)&get_reply, 0, 0, 0);
		while(get_reply < 1);
		for(i=0; i<size; i++){
                    //const int loc = where_it_goes(buf[i].key_hi_, parts);
                    const int loc = where_it_goes_ki(&buf[i], parts);
                    count[loc]++;
		}
	}
#endif
}



__attribute__((noinline))
static void inner_split_kernel(
		const KeyIndex      *inp,
		KeyIndex            *out,
		const int            n,
		const KeyIndex *parts,
		int                 *offset)
{
	int i;
	// naive version
#if 0
#  if 1
	for(i=0; i<n; i++){
            const int loc = where_it_goes(inp[i].key_hi_, parts);
            out[ offset[loc]++ ] = inp[i];
	}
#  else
	int ii;
	for(ii=0; ii<n; ii+=1024){
		int rem = n-ii;
		int size = rem<1024 ? rem : 1024;
		for(i=0; i<size; i++){
			// KeyIndex ki = inp[ii + i];
			//const int loc = where_it_goes(inp[ii+i].key, parts);
                    const int loc = where_it_goes(inp[ii+i].key_hi_, parts);
			out[offset[loc]++] = inp[ii+i];
		}
	}
#  endif
#else
	int ii;
	int count[NSPLIT] = {0, };

	KeyIndex work[NSPLIT][NCHUNK_SORT] __attribute__((aligned(256))); // 32 KiB in stack
	KeyIndex buf[1024] __attribute__((aligned(256))); // 16 KiB in stack

	int ilast = -1;

	dma_desc desc_put;
	dma_descriptor_init(&desc_put, (void *)&put_reply);
	dma_set_op         (&desc_put, DMA_PUT);
	// dma_set_mode       (&desc_put, PE_MODE);
	// dma_set_reply      (&desc_put, &reply_put);
	dma_set_size       (&desc_put, sizeof(work[0]));

	for(ii=0; ii<n; ii+=1024){
		int rem = n-ii;
		int size = rem<1024 ? rem : 1024;
		get_reply = 0;
                athread_get(PE_MODE, (void *)(inp+ii), buf, size*16, (unsigned long *)&get_reply, 0, 0, 0);
		while(get_reply < 1);
		for(i=0; i<size; i++){
			KeyIndex ki = buf[i];
                        const int loc = where_it_goes_ki(&ki, parts);
			if(loc == ilast){
				while(put_reply < 1);
				ilast = -1;
			}
			work[loc][count[loc] % NCHUNK_SORT ] = buf[i];
			count[loc]++;
			// DMA flush the work buffer
			if(0 == count[loc]%NCHUNK_SORT){ // buffer rounded
				if(-1 != ilast){
					while(put_reply < 1);
				}
				put_reply = 0;
#if 0
				athread_put(PE_MODE, work[loc], &out[offset[loc]+count[loc]-NCHUNK_SORT], sizeof(work[0]), (unsigned long *)&put_reply, 0, 0);
#else
				dma(desc_put, (long)&out[offset[loc]+count[loc]-NCHUNK_SORT], (long)work[loc]);
#endif
				// while(put_reply < 1);
				ilast = loc;
			}
		}
	}
	// return;
	if(-1 != ilast){
		while(put_reply < 1);
	}
	// Final flushing the DMA buffer
	int loc;
	int nreq = 0;
	put_reply = 0;
	for(loc=0; loc<NSPLIT; loc++){
		int nsend = count[loc] % NCHUNK_SORT;
		int sendoff = offset[loc] + count[loc]/NCHUNK_SORT*NCHUNK_SORT;
		if(nsend){
			nreq++;
#if 0
			athread_put(PE_MODE, work[loc], &out[sendoff], nsend*sizeof(KeyIndex), (unsigned long *)&put_reply, 0, 0);
#else
			dma_set_size(&desc_put, nsend * sizeof(KeyIndex));
			dma(desc_put, (long)&out[sendoff], (long)work[loc]);
#endif
		}
	}
	while(put_reply < nreq);
	for(loc=0; loc<NSPLIT; loc++){
		offset[loc] += count[loc];
	}
#endif
}


void final_kernel(){
    int i, j;
    int id0 = athread_get_id(-1);
    int myid;

    for(myid=id0; myid<NSPLIT; myid+=NTHREADS){
	const int beg = glb_task_off[myid  ];
	const int end = glb_task_off[myid+1];
	const int n = end - beg;
	if(n < 3000){
            //puts("small!!");
            KeyIndex *input  = (KeyIndex *)glb_ptr_out  + beg;
            KeyIndex *output = (KeyIndex *)glb_ptr_work + beg;
            KeyIndex *result = (KeyIndex *)glb_ptr_inp + beg;
            size_t buf_size = n * sizeof(KeyIndex);
            KeyIndex *buf = ldm_malloc(buf_size);

            get_reply = 0;
            athread_get(PE_MODE, &input[0], buf, buf_size, (unsigned long *)&get_reply, 0, 0, 0);
            while(get_reply < 1);

            my_quicksort_ki(buf, 0, n);
            
            put_reply = 0;
            athread_put(PE_MODE, buf, &output[0], buf_size, (unsigned long *)&put_reply, 0, 0);
            while(put_reply < 1);

            sync_array_(); // to avoid DMA bug

            get_reply = 0;
            athread_get(PE_MODE, &output[0], buf, buf_size, (unsigned long *)&get_reply, 0, 0, 0);
            while(get_reply < 1);

            my_quicksort_ki(buf, 0, n);

            put_reply = 0;
            //athread_put(PE_MODE, buf, &input[0], buf_size, (unsigned long *)&put_reply, 0, 0);
            athread_put(PE_MODE, buf, &result[0], buf_size, (unsigned long *)&put_reply, 0, 0);
            while(put_reply < 1);

            ldm_free(buf, buf_size);
            //return;
            continue;
	}
        //puts("big!!");
        //printf("id0=%d, n=%d \n", id0, n);
        
	const int oversample = 15;
        //const int oversample = 30;
	const int nsamples   = (oversample * NSPLIT) + NSPLIT - 1;

	//unsigned long *samples = ldm_malloc(8* nsamples);
        KeyIndex *samples = ldm_malloc(16* nsamples);

	const int stride = n / nsamples;
	j = (n%nsamples)/2;
	{
		KeyIndex *input  = (KeyIndex *)glb_ptr_out  + beg;
		for(i=0; i<nsamples; i++, j+=stride){
                    //samples[i] = input[j].key;
                    //samples[i] = input[j].key_hi_;
                    samples[i] = input[j];
		}
	}

        my_quicksort_ki(samples, 0, nsamples);
	KeyIndex parts[NSPLIT-1];
	for(i=oversample, j=0; i<nsamples; i+=oversample+1, j++){
		parts[j] = samples[i];
	}
	ldm_free(samples, 16 * nsamples);
        
        /*
          my_quicksort_ul(samples, 0, nsamples);
	unsigned long parts[NSPLIT-1];
	for(i=oversample, j=0; i<nsamples; i+=oversample+1, j++){
		parts[j] = samples[i];
	}
	ldm_free(samples, 8 * nsamples);
        */
#if 0
	for(ii=0; ii<64; ii++){
		if(myid == ii){
			for(i=0; i<NSPLIT-1; i++){
				printf("%ld, ", parts[i]);
			}
			puts("");
		}
		sync_array_();
	}
#endif
	int count[NSPLIT] = {0, };
	int offset[NSPLIT+1] = {0, };
	{
		KeyIndex *input  = (KeyIndex *)glb_ptr_out  + beg;
		inner_count_kernel(input, n, parts, count);
	}

	// Prefix sum
	for(i=0; i<NSPLIT; i++){
		offset[i+1] = offset[i] + count[i];
	}
	assert(offset[NSPLIT] == n);
#if 0
	for(ii=0; ii<64; ii++){
		if(myid == ii){
			for(i=0; i<NSPLIT+1; i++){
				printf("%d, ", offset[i]);
			}
			puts("");
		}
		sync_array_();
	}
#endif
	// out -> work
	{
		KeyIndex *input  = (KeyIndex *)glb_ptr_out  + beg;
		KeyIndex *output = (KeyIndex *)glb_ptr_work + beg;
		inner_split_kernel(input, output, n, parts, offset);
	}		 
	// DEBUG!!!!
	sync_array_();

	// final qsort
	// work -> out
	for(i=0; i<NSPLIT; i++){
		KeyIndex *input  = (KeyIndex *)glb_ptr_work + beg;
		//KeyIndex *output = (KeyIndex *)glb_ptr_out  + beg;
                KeyIndex *output = (KeyIndex *)glb_ptr_inp  + beg;

		const int ibeg = offset[i] - count[i];
		const int iend = offset[i];
		int num = iend - ibeg;
		if(num > 3000){
			//printf("Ignoring too big segment, id=%d, length=%d, in [%d:%d)\n",
         //       myid, num, beg+ibeg, beg+iend);	
			glb_status[myid] = -1;
		}else if(num > 0){ 
			size_t buf_size = num * sizeof(KeyIndex);
			KeyIndex *buf = ldm_malloc(buf_size);

			get_reply = 0;
			athread_get(PE_MODE, &input[ibeg], buf, buf_size, (unsigned long *)&get_reply, 0, 0, 0);
			while(get_reply < 1);

			my_quicksort_ki(buf, 0, num);

			put_reply = 0;
			athread_put(PE_MODE, buf, &output[ibeg], buf_size, (unsigned long *)&put_reply, 0, 0);
			while(put_reply < 1);

			ldm_free(buf, buf_size);
		}
	}
	// sync_array_();
  } // for(myid)

}

// inp -> out
void split_kernel(){
	int i, ii;
	int myid = athread_get_id(-1);
	unsigned long parts[NSPLIT-1];

	int off[NSPLIT] = {0, };
	int cnt[NSPLIT] = {0, };

	KeyIndex work[NSPLIT][NCHUNK_SORT]; // 32 KiB in stack
	// for(i=0; i<NSPLIT; i++) count[i] = 0;

	int beg = ((long)glb_enn * (myid+0)) / NTHREADS;
	int end = ((long)glb_enn * (myid+1)) / NTHREADS;

	// first fetch count
	get_reply = 0;
	athread_get(PE_MODE, glb_parts,       parts, sizeof(parts), (unsigned long *)&get_reply, 0, 0, 0);
	athread_get(PE_MODE, glb_count[myid], off,   sizeof(off),   (unsigned long *)&get_reply, 0, 0, 0);
	while(get_reply < 2);

	KeyIndex buf[1024]; // 16 KiB in stack
	int n = end - beg;

	int ilast = -1;

	dma_desc desc_put;
	dma_descriptor_init(&desc_put, (void *)&put_reply);
	dma_set_op         (&desc_put, DMA_PUT);
	// dma_set_mode       (&desc_put, PE_MODE);
	// dma_set_reply      (&desc_put, &reply_put);
	dma_set_size       (&desc_put, sizeof(work[0]));

	KeyIndex *ptr_inp  = (KeyIndex *)glb_ptr_inp;
	KeyIndex *ptr_out = (KeyIndex *)glb_ptr_out;

	for(ii=0; ii<n; ii+=1024){
		int rem = n-ii;
		int size = rem<1024 ? rem : 1024;
		get_reply = 0;
                athread_get(PE_MODE, ptr_inp+beg+ii, buf, size*16, (unsigned long *)&get_reply, 0, 0, 0);
		while(get_reply < 1);
		for(i=0; i<size; i++){
                    //const int loc = where_it_goes(buf[i].key, parts);
                    const int loc = where_it_goes(buf[i].key_hi_, parts);
			// glb_result[off[loc]++] = buf[i];
			if(loc == ilast){
				while(put_reply < 1);
				ilast = -1;
			}
			work[loc][cnt[loc] % NCHUNK_SORT ] = buf[i];
			cnt[loc]++;
			// DMA flush the work buffer
			if(0 == cnt[loc]%NCHUNK_SORT){ // buffer rounded
				if(-1 != ilast){
					while(put_reply < 1);
				}
				put_reply = 0;
#if 0
				athread_put(PE_MODE, work[loc], &glb_result[off[loc]+cnt[loc]-NCHUNK_SORT], sizeof(work[0]), (unsigned long *)&put_reply, 0, 0);
#else
				dma(desc_put, (long)&ptr_out[off[loc]+cnt[loc]-NCHUNK_SORT], (long)work[loc]);
#endif
				// while(put_reply < 1);
				ilast = loc;
			}
		}
	}
	if(-1 != ilast){
		while(put_reply < 1);
	}
	// Final flushing the DMA buffer
	int loc;
	int nreq = 0;
	put_reply = 0;
	for(loc=0; loc<NSPLIT; loc++){
		int nsend = cnt[loc] % NCHUNK_SORT;
		int sendoff = off[loc] + cnt[loc]/NCHUNK_SORT*NCHUNK_SORT;
		if(nsend){
			nreq++;
#if 0
			athread_put(PE_MODE, work[loc], &glb_result[sendoff], nsend*sizeof(KeyIndex), (unsigned long *)&put_reply, 0, 0);
#else
			dma_set_size(&desc_put, nsend*sizeof(KeyIndex) );
			dma(desc_put, (long)&ptr_out[sendoff], (long)work[loc]);
#endif
		}
	}
	while(put_reply < nreq);
}

void count_kernel(){
	int i, ii;
	int myid = athread_get_id(-1);
	unsigned long parts[NSPLIT-1];
	int count[NSPLIT] = {0, };
	// for(i=0; i<NSPLIT; i++) count[i] = 0;

	get_reply = 0;
	athread_get(PE_MODE, glb_parts, parts, sizeof(parts), (unsigned long *)&get_reply, 0, 0, 0);
	while(get_reply < 1);

	int beg = ((long)glb_enn * (myid+0)) / NTHREADS;
	int end = ((long)glb_enn * (myid+1)) / NTHREADS;
#if 0
	for(i=beg; i<end; i++){
            //const int loc = where_it_goes(glb_input[i].key, parts);
            const int loc = where_it_goes(glb_input[i].key_hi_, parts);
            count[loc]++;
	}
#else
	KeyIndex buf[2048]; // 32 KiB in stack
	int n = end - beg;
	for(ii=0; ii<n; ii+=2048){
		int rem = n-ii;
		int size = rem<2048 ? rem : 2048;
		get_reply = 0;
		KeyIndex *ptr_inp = (KeyIndex *)glb_ptr_inp;
                athread_get(PE_MODE, ptr_inp+beg+ii, buf, size*16, (unsigned long *)&get_reply, 0, 0, 0);
		while(get_reply < 1);
		for(i=0; i<size; i++){
                    //const int loc = where_it_goes(buf[i].key, parts);
                    const int loc = where_it_goes(buf[i].key_hi_, parts);
                    count[loc]++;
		}
	}
#endif

	put_reply = 0;
	athread_put(PE_MODE, count, glb_count[myid], sizeof(count), (unsigned long *)&put_reply, 0, 0);
	while(put_reply < 1);
}



