#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>
#include<iostream>
extern "C"{
#include <athread.h>
#include "cpe_func.h"
}
//#include "common.h"
#include"mpi.h"

extern int     MY_RANK_MPI;

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

struct CompKey0{
    bool operator() (const KeyIndex & lhs,
                     const KeyIndex & rhs) {
       return lhs.key_hi_ < rhs.key_hi_;
   }
};


struct CompKey1{
    bool operator() (const KeyIndex & lhs,
                     const KeyIndex & rhs) {
#ifdef USE_96BIT_KEY
        return (lhs.key_hi_ != rhs.key_hi_) ? (lhs.key_hi_ < rhs.key_hi_) : (lhs.key_lo_ < rhs.key_lo_);
#else
        return lhs.key_hi_ < rhs.key_hi_;
#endif
   }
};


extern "C"{
    unsigned long glb_parts[NSPLIT-1];

    int glb_task_off[NSPLIT+1];
    int glb_status  [NSPLIT];

    int glb_count[NTHREADS][NSPLIT];

    int glb_enn;
    void *glb_ptr_inp; // inp and out
    void *glb_ptr_work; // work
    void *glb_ptr_out; // work2

    extern void SLAVE_FUN(count_kernel)();
    extern void SLAVE_FUN(split_kernel)();
    extern void SLAVE_FUN(final_kernel)();

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
		void *ptr_out){
    
	::glb_enn = n;
	::glb_ptr_inp  = ptr_inp;
	::glb_ptr_work = ptr_work;
	::glb_ptr_out  = ptr_out;

	const int nsamples = (OVERSAMPKLE * NSPLIT) + NSPLIT - 1;
	assert(nsamples < n);

	const KeyIndex *inp = (const KeyIndex *)ptr_inp;
	unsigned long *samples = (unsigned long *)ptr_out;

	const int stride = n / nsamples;


#if defined(DEBUG_SAMPLE_SORT)
        MPI_Barrier(MPI_COMM_WORLD);
        if(MY_RANK_MPI==0){
            std::cerr<<"sizeof(KeyIndex)= "<<sizeof(KeyIndex)<<std::endl;
            printf("nsamples = %d, stride = %d\n", nsamples, stride);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        
	for(int i=0, j = (n%nsamples)/2 ; i<nsamples; i++, j+=stride){
		assert(j < n);
                samples[i] = inp[j].key_hi_;
	}

	std::sort( samples, samples + nsamples);

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
        /*
        MPI_Barrier(MPI_COMM_WORLD);
        if(MY_RANK_MPI==0){
            for(int i=0; i<NSPLIT-1; i++){
                std::cerr<<"i= "<<i<<" glb_parts[i]= "<<glb_parts[i]<<std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        */
        
	__real_athread_spawn((void *)slave_count_kernel, (void *)0);
	athread_join();

	// Prefix sum
	int sum = 0;
	for(int i=0; i<NSPLIT; i++){
            for(int j=0; j<NTHREADS; j++){
#if defined(DEBUG_SAMPLE_SORT)
                if(MY_RANK_MPI==0){
                    std::cerr<<"j= "<<j<<" i= "<<i
                             <<" glb_count[j][i]= "<<glb_count[j][i]
                             <<std::endl;
                }
#endif
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

#if defined(DEBUG_SAMPLE_SORT)
        MPI_Barrier(MPI_COMM_WORLD);
        if(MY_RANK_MPI==0){
            std::cerr<<"sizeof(KeyIndex)= "<<sizeof(KeyIndex)<<std::endl;
            printf("max : %d, min : %d, avr : %d\n", max, min, n/NSPLIT);
            fflush(stdout);
            for(int i=0; i<NSPLIT-1; i++){
                std::cerr<<"i= "<<i<<" glb_parts[i]= "<<glb_parts[i]<<std::endl;
            }
            for(int i=0; i<NSPLIT; i++){
                std::cerr<<"i= "<<i<<" glb_task_off[i]= "<<glb_task_off[i]<<std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        
	__real_athread_spawn((void *)slave_split_kernel, (void *)0); // inp -> out(work2)
	athread_join();

#if defined(DEBUG_SAMPLE_SORT)
        for(int i=0; i<NSPLIT; i++){
            for(int j=glb_task_off[i]; j<glb_task_off[i+1]; j++){
                if(i==0)  assert(((KeyIndex *)glb_ptr_out+j)->key_hi_ <= glb_parts[i]);
                else if(i==NSPLIT-1) assert(((KeyIndex *)glb_ptr_out+j)->key_hi_ >= glb_parts[i-1]);
                else assert( ((KeyIndex *)glb_ptr_out+j)->key_hi_ <= glb_parts[i] && ((KeyIndex *)glb_ptr_out+j)->key_hi_ >= glb_parts[i-1]);
            }
        }
#endif
        
#if 0
        for(int i=0; i<NSPLIT; i++){
            //std::sort((KeyIndex *)glb_ptr_out+glb_task_off[i], (KeyIndex *)glb_ptr_out+glb_task_off[i+1]);
            std::sort((KeyIndex *)glb_ptr_out+glb_task_off[i], (KeyIndex *)glb_ptr_out+glb_task_off[i+1], CompKey1());
        }
#else
        __real_athread_spawn((void *)slave_final_kernel, (void *)0); // out(work2) -> work(work1) -> out(inp)
	athread_join();
#endif
        /*
        MPI_Barrier(MPI_COMM_WORLD);
        if(MY_RANK_MPI==0){
            printf("final sort done");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        */
        
	//KeyIndex *out = (KeyIndex *)ptr_out;
        KeyIndex *out = (KeyIndex *)ptr_inp;
        /*
        MPI_Barrier(MPI_COMM_WORLD);
        if(MY_RANK_MPI==0){
            printf("check A \n");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        */
#if 0
	fflush(stdout);
	for(int i=0; i<NSPLIT; i++){
		printf("%d ", glb_status[i]);
	}
	puts("");
#endif
	for(int i=0; i<NSPLIT; i++){
		if(glb_status[i]){
			int beg = glb_task_off[i+0];
			int end = glb_task_off[i+1];
			printf("error handling, [%d, %d)\n", beg, end);
			KeyIndex *work = (KeyIndex *)ptr_work;
			for(int i=beg; i<end; i++){
                            // out[i] = inp[i];
                            out[i] = work[i];
			}
			std::sort(out+beg, out+end, CompKey1());
		}
	}
        /*
        MPI_Barrier(MPI_COMM_WORLD);
        if(MY_RANK_MPI==0){
            printf("check B \n");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        */
#if defined(DEBUG_SAMPLE_SORT)
	// check result
	for(int i=1; i<n; i++){
		// is it really sorted?
            if(MY_RANK_MPI==0){
                if(!(inp[i].key_hi_ >= inp[i-1].key_hi_)){
                    printf("key[%d]=%ld (index: %ld), key[%d]=%ld (index: %ld)\n",
                           i-1, inp[i-1].key_hi_, inp[i-1].index,
                           i-0, inp[i-0].key_hi_, inp[i-0].index);
                    fflush(stdout);
		}
            }
	}
	for(int i=1; i<n; i++){
    #ifdef USE_96BIT_KEY
            assert( (inp[i].key_hi_ > inp[i-1].key_hi_)
                    || ( (inp[i].key_hi_ == inp[i-1].key_hi_) && (inp[i].key_lo_ >= inp[i-1].key_lo_)));
    #else
            assert( (inp[i].key_hi_ >= inp[i-1].key_hi_));
    #endif
	}
#endif
        /*
        MPI_Barrier(MPI_COMM_WORLD);
        if(MY_RANK_MPI==0){
            printf("check C \n");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        */

}
