#include <slave.h>
#include <simd.h>

// extern __thread_local uint256 rsqrtd_magic;
// extern __thread_local uint256 rsqrtd_mask;
#define REP4(x) x,x,x,x
__thread_local uint256 rsqrtd_magic = { REP4(0x5fe6eb50c7b537a9) };
__thread_local uint256 rsqrtd_mask  = { REP4(0x7fffffffffffffff) };


void rsqrt_loop(const int N, void *work){
	asm volatile("stl $26,  8($31)");
	asm volatile("stl $29, 16($31)");
	asm volatile("stl $30, 24($31)");
	asm volatile("stl $15, 40($31)");
	asm volatile("stl $27, 48($31)");
	/*
	 * $16 : N
	 * $17 : work
	 * $18 : (i<N)
	 * 
	 * $9, $19,$29,$24 : x0, x1, x2, x3
	 * $10,$20,$30,$25 : y0, y1, y2, y3
	 * $11,$21,$1, $26 : h0, h1, h2, h3
	 * $12,$22,$14,$27 : s0, s1, s2, s3 (aka h^2 prevousely)
	 * $13,$23,$15,$28 : b0, b1, b2, b3
	 *
	 * $2 : c1
	 * $3 : c2
	 * $4 : c3
	 * $5 : c4
	 * $6 : c5
	 *
	 * $7 : rsqrtd_mask
	 * $8 : rsqrtd_magic
	 *
	 */
	/*
	 * &work[0].r2 = 128 + 0*192 = 128
	 * &work[1].r2 = 128 + 1*192 = 320
	 * &work[2].r2 = 128 + 2*192 = 512
	 * &work[3].r2 = 128 + 3*192 = 704
	 *
	 * &work[-4].rinv = 160 - 4*192 = -608
	 * &work[-3].rinv = 160 - 4*192 = -416
	 * &work[-2].rinv = 160 - 4*192 = -224
	 * &work[-1].rinv = 160 - 4*192 =  -32
	 */
	// int i;
#define WINC  "768"
#define LOFF0 "128"
#define LOFF1 "320"
#define LOFF2 "512"
#define LOFF3 "704"
#define SOFF0 "-608"
#define SOFF1 "-416"
#define SOFF2 "-224"
#define SOFF3 "-32"
	
#if 1
	// $2=1.0;
   	// $3=1./2.;
   	// $4=3./8.;
   	// $5=35./128.;
   	// $6=75./128.;
	asm volatile(
			"ldi   $6, -16384($31) \n\t"
			"ldi   $2, 1023($31)   \n\t"
			"ldi   $3, 511($31)    \n\t"
			"ldi   $4, 2043($31)   \n\t"
			"ldi   $5, 32675($31)  \n\t"
			"ldih  $6, 16355($6)   \n\t"
			"sll   $2, 52, $2      \n\t"
			"sll   $3, 53, $3      \n\t"
			"sll   $4, 51, $4      \n\t"
			"sll   $5, 47, $5      \n\t"
			"sll   $6, 32, $6      \n\t"
			"vshff $2, $2, $31, $2 \n\t"
		"vldd  $9,  " LOFF0 "($17)        \n\t"
			"vshff $3, $3, $31, $3 \n\t" 
		"vldd  $19, " LOFF1 "($17)        \n\t"
			"vshff $4, $4, $31, $4 \n\t" 
		"vldd  $29, " LOFF2 "($17)        \n\t"
			"vshff $5, $5, $31, $5 \n\t" 
		"vldd  $24, " LOFF3 "($17)        \n\t"
			"vshff $6, $6, $31, $6" );
#else
	doublev4 c1=1.0, c2=1./2., c3=3./8., c4=35./128., c5=75./128.;

	asm volatile(
		"vfmov %4, $6 \n\t"
		"vfmov %3, $5 \n\t"
		"vfmov %2, $4 \n\t"
		"vfmov %1, $3 \n\t"
		"vfmov %0, $2 \n\t"
		:
		: "r"(c1), "r"(c2), "r"(c3), "r"(c4), "r"(c5)
		);
#endif


	asm volatile(
		"ldih $7,rsqrtd_mask ($31) !tprelhi \n\t"
		"ldih $8,rsqrtd_magic($31) !tprelhi \n\t"
		"vldd $7,rsqrtd_mask ($7)  !tprello \n\t"
		"vldd $8,rsqrtd_magic($8)  !tprello "
		);

	asm volatile(
		// "vldd  $9,  128($17)        \n\t"
		// "vldd  $19, 320($17)        \n\t"
		// "vldd  $29, 512($17)        \n\t"
		// "vldd  $24, 704($17)        \n\t"

		"mov $31, $0 \n\t"
		".align 8 \n"
		".Lt_fdps_1:"
		);
	// for(i=0; i<N; i+=4){
		asm volatile(
			"srlow $9,  1,   $10      \n\t"
			"srlow $19, 1,   $20      \n\t"
			"srlow $29, 1,   $30      \n\t"
			"srlow $24, 1,   $25      \n\t"

			"vandw $10, $7,  $10      \n\t"
			"vandw $20, $7,  $20      \n\t"
			"vandw $30, $7,  $30      \n\t"
			"vandw $25, $7,  $25      \n\t"

			"vsubl $8,  $10, $10      \n\t"
			  "vstd  $12, " SOFF0 "($17)      \n\t"
			"vsubl $8,  $20, $20      \n\t"
			  "vstd  $22, " SOFF1 "($17)      \n\t"
			"vsubl $8,  $30, $30      \n\t"
			  "vstd  $14, " SOFF2 "($17)      \n\t"
			"vsubl $8,  $25, $25      \n\t"
			  "vstd  $27, " SOFF3 "($17)      \n\t"

			"vmuld $10, $10, $12      \n\t" // y, y, s
			"vmuld $20, $20, $22      \n\t" // y, y, s
			"vmuld $30, $30, $14      \n\t" // y, y, s
			"vmuld $25, $25, $27      \n\t" // y, y, s

			"vmuld $9,  $5,  $13      \n\t" // x, c4, b
			"vmuld $19, $5,  $23      \n\t" // x, c4, b
			"vmuld $29, $5,  $15      \n\t" // x, c4, b
			"vmuld $24, $5,  $28      \n\t" // x, c4, b

			"vnmad $9,  $12, $2, $11  \n\t" // x, s, c1, h
			"vnmad $19, $22, $2, $21  \n\t" // x, s, c1, h
			"vnmad $29, $14, $2, $1   \n\t" // x, s, c1, h
			"vnmad $24, $27, $2, $26  \n\t" // x, s, c1, h

			"vnmad $13, $12, $6, $13  \n\t" // b, s, c5, b
			"vnmad $23, $22, $6, $23  \n\t" // b, s, c5, b
			"vnmad $15, $14, $6, $15  \n\t" // b, s, c5, b
			"vnmad $28, $27, $6, $28  \n\t" // b, s, c5, b

			"vmuld $11, $11, $12      \n\t" // h, h, s
			"vmuld $21, $21, $22      \n\t" // h, h, s
			"vmuld $1,  $1,  $14      \n\t" // h, h, s
			"vmuld $26, $26, $27      \n\t" // h, h, s

			"vmad  $4,  $11, $3, $9   \n\t" // c3, h, c2, x
			"vmad  $4,  $21, $3, $19  \n\t" // c3, h, c2, x
			"vmad  $4,  $1,  $3, $29  \n\t" // c3, h, c2, x
			"vmad  $4,  $26, $3, $24  \n\t" // c3, h, c2, x
			
			"vmuld $11, $10, $11      \n\t" // h, y, h
			"vmuld $21, $20, $21      \n\t" // h, y, h
			"vmuld $1,  $30, $1       \n\t" // h, y, h
			"vmuld $26, $25, $26      \n\t" // h, y, h

			"vmad  $12, $13, $9,  $9  \n\t" // s, b, x, x
			// "ldi   $18, 128($18)       \n\t" // yout++;
			"vmad  $22, $23, $19, $19 \n\t" // s, b, x, x
			"ldi   $17, " WINC "($17)       \n\t" // work++;
			"vmad  $14, $15, $29, $29 \n\t" // s, b, x, x
			  "addw  $0, 4, $0   \n\t"
			"vmad  $27, $28, $24, $24 \n\t" // s, b, x, x
			  "cmplt $0, $16, $18 \n\t"

			"vmad  $11, $9,  $10, $12 \n\t" // h, x, y, s
			  "vldd  $9,  " LOFF0 "($17)     \n\t"
			"vmad  $21, $19, $20, $22 \n\t" // h, x, y, s
			  "vldd  $19, " LOFF1 "($17)     \n\t"
			"vmad  $1,  $29, $30, $14 \n\t" // h, x, y, s
			  "vldd  $29, " LOFF2 "($17)     \n\t"
			"vmad  $26, $24, $25, $27 \n\t" // h, x, y, s
			  "vldd  $24, " LOFF3 "($17)     \n\t"


			"bne   $18,.Lt_fdps_1"
			);
	// }
	asm volatile(
		"vstd  $12, " SOFF0 "($17)      \n\t"
		"vstd  $22, " SOFF1 "($17)      \n\t"
		"vstd  $14, " SOFF2 "($17)      \n\t"
		"vstd  $27, " SOFF3 "($17)      \n\t"
			);
	asm volatile("ldl $26,  8($31)");
	asm volatile("ldl $29, 16($31)");
	asm volatile("ldl $30, 24($31)");
	asm volatile("ldl $15, 40($31)");
	asm volatile("ldl $27, 48($31)");
	asm volatile("stl $31,  8($31)");
	asm volatile("stl $31, 16($31)");
	asm volatile("stl $31, 24($31)");
	asm volatile("stl $31, 40($31)");
	asm volatile("stl $31, 48($31)");
}
