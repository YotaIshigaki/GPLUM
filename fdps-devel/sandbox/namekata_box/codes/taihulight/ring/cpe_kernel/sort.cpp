#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>
extern "C"{
#include <athread.h>
#include "cpe_func.h"
}

/*
#include <sys/time.h>
static double wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return 1.e-6*tv.tv_usec +tv.tv_sec;
}*/

// MPE cycle counts
static inline unsigned long rpcc_() {
	unsigned long rpcc;
	asm volatile ("rtc %0": "=r" (rpcc) : );
	return rpcc;
}

static inline int operator<(const KeyIndex &lhs, const KeyIndex &rhs){
	return lhs.key < rhs.key;
}

extern "C"{
//    unsigned long glb_parts[NSPLIT-1];
// 
//    int glb_task_off[NTHREADS+1];
//    int glb_status  [NTHREADS];
// 
//    int glb_count[NTHREADS][NSPLIT];
// 
//    int glb_enn;
//    void *glb_ptr_inp;
//    void *glb_ptr_work;
//    void *glb_ptr_out;

	extern void SLAVE_FUN(count_kernel)(void*);
	extern void SLAVE_FUN(split_kernel)(void*);
	extern void SLAVE_FUN(final_kernel)(void*);

}

/*
 * Usage: 
 *   ptr_inp, ptr_work, ptr_out should be (16*n)-byte buffer
 *   ptr_inp and ptr_work can point the same address when
 *   the data in ptr_inp can be discarded
 */

void samplesort(
		const int n,
		void *ptr_inp,
		void *ptr_work,
		void *ptr_out)
{
    sort_pars ipar;
    unsigned long glb_parts[NSPLIT-1];
    int glb_task_off[NTHREADS+1];
    int glb_status[NTHREADS];
    int glb_count[NTHREADS][NSPLIT];
    
    ipar.N = n;
    ipar.glb_input    = ptr_inp;
    ipar.glb_work     = ptr_work;
    ipar.glb_result   = ptr_out;

    ipar.glb_parts    = glb_parts;
    ipar.glb_task_off = glb_task_off;
    ipar.glb_status   = glb_status;
    ipar.glb_count    = glb_count[0];

//	::glb_enn = n;
//	::glb_ptr_inp  = ptr_inp;
//	::glb_ptr_work = ptr_work;
//	::glb_ptr_out  = ptr_out;

	const int nsamples = (OVERSAMPKLE * NSPLIT) + NSPLIT - 1;
	assert(nsamples < n);

	const KeyIndex *inp = (const KeyIndex *)ptr_inp;
	unsigned long *samples = (unsigned long *)ptr_out;

	const int stride = n / nsamples;
#ifdef SORT_DEBUG
	printf("nsamples = %d, stride = %d\n", nsamples, stride);
#endif
	for(int i=0, j = (n%nsamples)/2 ; i<nsamples; i++, j+=stride){
		assert(j < n);
		samples[i] = inp[j].key;
	}

	std::sort( samples, samples + nsamples );

	{
		int j=0;
		for(int i=OVERSAMPKLE; i<nsamples; i+=OVERSAMPKLE+1, j++){
			glb_parts[j] = samples[i];
		}
		assert(NSPLIT-1 == j);
	} 

	for(int j=0; j<NTHREADS; j++){
		for(int i=0; i<NSPLIT; i++){
			glb_count[j][i] = 0;
		}
	}

	__real_athread_spawn((void *)slave_count_kernel, (void *)&ipar);
	athread_join();

	// Prefix sum
	int sum = 0;
	for(int i=0; i<NSPLIT; i++){
		for(int j=0; j<NTHREADS; j++){
			int loc = sum;
			sum += glb_count[j][i];
			glb_count[j][i] = loc;
		}
	}
	assert(n == sum);

	int max=0, min=n;
	glb_task_off[0] = 0;
	for(int i=0; i<NSPLIT; i++){
		const int beg = glb_count[0][i];
		const int end = (i == NSPLIT-1) ? n : glb_count[0][i+1];
		// printf("%3d : %8d\n", i, end-beg);
		max = std::max(max, end-beg);
		min = std::min(min, end-beg);
		glb_task_off[i+1] = end;
		// printf("CPE %d : [%d:%d)\n", i, beg, end);
		glb_status[i] = 0;
	}
#ifdef SORT_DEBUG
	printf("max : %d, min : %d, avr : %d\n", max, min, n/NSPLIT);
	fflush(stdout);
#endif

	__real_athread_spawn((void *)slave_split_kernel, (void *)&ipar); // inp -> out
	athread_join();
	// puts("split done");

	__real_athread_spawn((void *)slave_final_kernel, (void *)&ipar); // out -> work -> out
	athread_join();
	// puts("final sort done");

	// error handling
	KeyIndex *out = (KeyIndex *)ptr_out;
	for(int i=0; i<NTHREADS; i++){
		if(glb_status[i]){
			int beg = glb_task_off[i+0];
			int end = glb_task_off[i+1];
			printf("error handling, [%d, %d)\n", beg, end);
			for(int i=beg; i<end; i++){
				out[i] = inp[i];
			}
			std::sort(out+beg, out+end);
		}
	}
	
#if 0
	// check result
	for(int i=1; i<n; i++){
		// is it really sorted?
		if(!(out[i].key >= out[i-1].key)){
			printf("key[%d]=%ld (index: %ld), key[%d]=%ld (index: %ld)\n",
					i-1, out[i-1].key, out[i-1].index,
					i-0, out[i-0].key, out[i-1].index);
			fflush(stdout);
		}
	}
	for(int i=1; i<n; i++){
		assert(out[i].key >= out[i-1].key);
	}
#endif
}
