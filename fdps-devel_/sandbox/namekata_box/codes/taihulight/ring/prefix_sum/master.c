#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <athread.h>

static inline unsigned long rpcc_() {
	unsigned long rpcc = 0;
	asm volatile("rtc %0": "=r" (rpcc) : );
	return rpcc;
}

extern void SLAVE_FUN(func)();

int main() {
	athread_init();
	{
		int j;
		long t0 = rpcc_();
		for(j=0; j<1; j++){
			athread_spawn(func, 0);
			athread_join();
		}
		long t1 = rpcc_();

		double dt = (t1-t0)/1.45e9/1000;
		printf("%e sec per fork-join\n", dt);
	}
	athread_halt();

	// reference
	unsigned long val[64], sum[64], i, j;
	for(i=0; i<64; i++){
		unsigned long v = i+11111;
		v *= v;
		v *= v;
		v *= v;
		v *= v;
		v %= 10;
		val[i] = v;
	}
	for(i=0; i<64; i++){
		unsigned long s = 0;
		for(j=0; j<=i; j++){
			s += val[j];
		}
		sum[i] = s;
	}
	puts("MPE reference (exclusive):");
	for(i=0; i<64; i++){
		printf("%8ld", sum[i] - val[i]);
		if(7 == i%8) printf("\n");
	}

	return 0;
}
