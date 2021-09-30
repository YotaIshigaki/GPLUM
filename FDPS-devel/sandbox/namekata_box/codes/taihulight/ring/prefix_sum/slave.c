#include <slave.h>
#include "cpe_prof.h"

void print_int_array(FILE *fp, int val, const char *msg){
	const int myid = athread_get_id(-1);
	int i;

	if(0 == myid){
		printf("%s\n", msg);
	}
	for(i=0; i<64; i++){
		if(myid == i){
			fprintf(fp,"%8d", val);
			if(7 == i%8) fprintf(fp, "\n");
			fflush(fp);
		}
		sync_array_();
	}
}

typedef unsigned long sumtype;
void prefix_sum(sumtype val, sumtype *beg, sumtype *end, sumtype *psum){
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


void func() {
	const int myid = athread_get_id(-1);
	const int jcol = myid % 8;
	const int irow = myid / 8;
	// int i;
	
	unsigned long val = myid+11111;
	val *= val;
	val *= val;
	val *= val;
	val *= val;
	val %= 10;

	unsigned long beg, end, sum;

	print_int_array(stdout, 10*(irow+1) + (jcol+1), "position:");
	print_int_array(stdout, val, "original:");

	prefix_sum(val, &beg, &end, &sum);
	print_int_array(stdout, beg, "exclusive:");
	print_int_array(stdout, end, "inclusive:");
	print_int_array(stdout, sum, "summation:");

#if 0
	unsigned long sbuf, rbuf;

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
	print_int_array(stdout, val, "after __shfl_up(1):");

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
	print_int_array(stdout, val, "after __shfl_up(2):");

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
	print_int_array(stdout, val, "after __shfl_up(4):");

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
	print_int_array(stdout, val, "after __shfl_up(8):");

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
	print_int_array(stdout, val, "after __shfl_up(16):");

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
	print_int_array(stdout, val, "after __shfl_up(32):");
	print_int_array(stdout, val-v0, "exclusive:");

	// bcast the sum
	unsigned long sum = val;
	if(7 == irow){
		if(7 == jcol){
			REG_PUTR(sum, 8);
		}else{
			REG_GETR(sum);
		}
		REG_PUTC(sum, 8);
	}else{
		REG_GETC(sum);
	}
	print_int_array(stdout, sum, "summation:");
#endif

	return;
}
