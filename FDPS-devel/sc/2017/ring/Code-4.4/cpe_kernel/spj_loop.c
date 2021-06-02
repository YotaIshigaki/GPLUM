// Memory map
#define LDM_BASE "3200+"
#define S0 LDM_BASE "0($31)"
#define S1 LDM_BASE "8($31)"
#define S2 LDM_BASE "16($31)"
#define S3 LDM_BASE "24($31)"
#define V1 LDM_BASE "32($31)"
#define V2 LDM_BASE "64($31)"
#define V3 LDM_BASE "96($31)"
#define V4 LDM_BASE "128($31)"
#define V5 LDM_BASE "160($31)"
#define V6 LDM_BASE "192($31)"
#define V7 LDM_BASE "224($31)"
#define V8 LDM_BASE "256($31)"
#define V9 LDM_BASE "288($31)"
#define V10 LDM_BASE "320($31)"
#define V11 LDM_BASE "352($31)"
#define V12 LDM_BASE "384($31)"
#define V13 LDM_BASE "416($31)"
#define V14 LDM_BASE "448($31)"
#define V15 LDM_BASE "480($31)"
#define V16 LDM_BASE "512($31)"
#define V17 LDM_BASE "544($31)"
#define V18 LDM_BASE "576($31)"
#define V19 LDM_BASE "608($31)"
#define V20 LDM_BASE "640($31)"
#define V21 LDM_BASE "672($31)"
#define V22 LDM_BASE "704($31)"
#define V23 LDM_BASE "736($31)"
#define V24 LDM_BASE "768($31)"
#define V25 LDM_BASE "800($31)"
#define V26 LDM_BASE "832($31)"
#define V27 LDM_BASE "864($31)"
#define V28 LDM_BASE "896($31)"
#define V29 LDM_BASE "928($31)"
#define V30 LDM_BASE "960($31)"
#define V31 LDM_BASE "992($31)"

typedef double doublev4 __attribute__ ((__mode__(__V4DF__)));

__attribute__((noinline))
void sp_init_ldm_constants(){
	long rsqrtd_mask  = 0x7fffffffffffffff;
	asm volatile(
		"vshff %0, %0, $31, %0 \n\t" 
		"vstd %0, " V1 
		:: "r"(rsqrtd_mask)  );
	long rsqrtd_magic = 0x5fe6eb50c7b537a9;
	asm volatile(
		"vshff %0, %0, $31, %0 \n\t" 
		"vstd %0, " V2
		:: "r"(rsqrtd_magic)  );
	doublev4 c0 = 1.0;
	asm volatile ("vstd %0, " V3 :: "r"(c0) );
	doublev4 c1 = 1.0/2.0;
	asm volatile ("vstd %0, " V4 :: "r"(c1) );
	doublev4 c2 = 3.0/8.0;
	asm volatile ("vstd %0, " V5 :: "r"(c2) );
	doublev4 c3 = 5.0/16.0;
	asm volatile ("vstd %0, " V6 :: "r"(c3) );
	doublev4 c4 = 35.0/128.0;
	asm volatile ("vstd %0, " V7 :: "r"(c4) );
}
/*
 * $0-$7   : work1
 * $8-$15  : work2
 * $16-$23 : work3
 * $24-$30 : misc
 *
 *
 *
 *
 *
 *
 */

void CalcGrav_sp_loop(
		/*
		const int    n_epi, 
		const EpiMM *epi,     // double r[4], v[4]
		const EpjMM *epj,     // double r[4], v[4]
		ForceMM4    *force,   // double ax[4], ay[4], az[4], pot[4];
		*/
		)
{
	// save registers
	asm volatile(
		"stl $26,  8($31)\n\t"
		"stl $29, 16($31)\n\t"
		"stl $30, 24($31)\n\t"
		"stl $15, 40($31)\n\t"
		"stl $27, 48($31)\n\t"
	);
	asm volatile(
#if 0
		"addw $16, 3,   $16  \n\t"
		"srl  $16, 2,   $16  \n\t"
		"sll  $16, 7,   $16  \n\t"
#else
		"addw $16, 7,   $16  \n\t"
		"srl  $16, 3,   $16  \n\t"
		"sll  $16, 8,   $16  \n\t"
#endif
		"addl $16, $19, $16  \n\t" // fend
		"stl  $17, " S0 " \n\t" // epi
		"stl  $18, " S1 " \n\t" // epj
		"stl  $19, " S2 " \n\t" // force
		"stl  $16, " S3 " \n\t" // fend
		"mov $17, $28 \n\t" // epi
		"mov $19, $30 \n\t" // force
		// "mov $16, $30 \n\t" // fend
		// load epj
		"vldd $24,  0($18) \n\t"
		"vldd $25, 32($18) \n\t"
		"vldd $26, 64($18) \n\t"
		"vldd $27, 96($18) \n\t"
	);
	// loop body
	asm volatile(
		".align 6 \n"
		".Lt_fdps_4: \n\t"

		// load xi
		// shuffle xj
		"vshff $24, $24, $31, $8  \n\t"
		  "vldd $12,   0($28) \n\t"
		"vshff $25, $25, $31, $9  \n\t"
		  "vldd $13, 128($28) \n\t"
		"vshff $26, $26, $31, $10 \n\t"
		  "stl  $28, " S0 " \n\t" // epi
		"vshff $27, $27, $31, $11 \n\t"
		  "stl  $30, " S2 " \n\t" // force
		// dx
		"vsubd $8,  $12, $0 \n\t"
		  "vstd $16, -256($30) \n\t" // fx0
		"vsubd $9,  $12, $1 \n\t"
		  "vstd $17, -224($30) \n\t" // fy0
		"vsubd $10, $12, $2 \n\t"
		  "vstd $18, -192($30) \n\t" // fz0
		"vsubd $11, $12, $3 \n\t"
		  "vstd $19, -128($30) \n\t" // fx1
		"vsubd $8,  $13, $4 \n\t"
		  "vstd $20,  -96($30) \n\t" // fy1
		"vsubd $9,  $13, $5 \n\t"
		  "vstd $21,  -64($30) \n\t" // fz1
		"vsubd $10, $13, $6 \n\t"
		"vsubd $11, $13, $7 \n\t"
		// dx*dx
		"vmuld $0, $0, $8  \n\t" "vstd $0, " V8  " \n\t"
		"vmuld $1, $1, $9  \n\t" "vstd $1, " V9  " \n\t"
		"vmuld $2, $2, $10 \n\t" "vstd $2, " V10 " \n\t"
		"vmuld $3, $3, $11 \n\t" "vstd $3, " V11 " \n\t"
		"vmuld $4, $4, $12 \n\t" "vstd $4, " V12 " \n\t"
		"vmuld $5, $5, $13 \n\t" "vstd $5, " V13 " \n\t"
		  "vldd $16,  32($28) \n\t" // load yi
		"vmuld $6, $6, $14 \n\t" "vstd $6, " V14 " \n\t"
		  "vldd $17, 160($28) \n\t"
		"vmuld $7, $7, $15 \n\t" "vstd $7, " V15 " \n\t"
		// shuffle yj
		"ldi $30, 85($31) \n\t"
		"vshff $24, $24, $30, $18 \n\t"
		"vshff $25, $25, $30, $19 \n\t"
		"vshff $26, $26, $30, $20 \n\t"
		"vshff $27, $27, $30, $21 \n\t"
		// dy
		"vsubd $18, $16, $0 \n\t"
		"vsubd $19, $16, $1 \n\t"
		"vsubd $20, $16, $2 \n\t"
		"vsubd $21, $16, $3 \n\t"
		"vsubd $18, $17, $4 \n\t"
		"vsubd $19, $17, $5 \n\t"
		"vsubd $20, $17, $6 \n\t"
		"vsubd $21, $17, $7 \n\t"
		// dy*dy
		"vmad $0, $0, $8,  $8  \n\t" "vstd $0, " V16 " \n\t"
		"vmad $1, $1, $9,  $9  \n\t" "vstd $1, " V17 " \n\t"
		"vmad $2, $2, $10, $10 \n\t" "vstd $2, " V18 " \n\t"
		"vmad $3, $3, $11, $11 \n\t" "vstd $3, " V19 " \n\t"
		"vmad $4, $4, $12, $12 \n\t" "vstd $4, " V20 " \n\t"
		"vmad $5, $5, $13, $13 \n\t" "vstd $5, " V21 " \n\t"
		  "vldd $16,  64($28) \n\t"
		"vmad $6, $6, $14, $14 \n\t" "vstd $6, " V22 " \n\t"
		  "vldd $17, 192($28) \n\t"
		"vmad $7, $7, $15, $15 \n\t" "vstd $7, " V23 " \n\t"
		// load zi
		// shuffle zj
		"ldi $30, 170($31) \n\t"
		"vshff $24, $24, $30, $18 \n\t"
		"vshff $25, $25, $30, $19 \n\t"
		"vshff $26, $26, $30, $20 \n\t"
		"vshff $27, $27, $30, $21 \n\t"
		// dz
		"vsubd $18, $16, $0 \n\t"
		"vsubd $19, $16, $1 \n\t"
		"vsubd $20, $16, $2 \n\t"
		"vsubd $21, $16, $3 \n\t"
		"vsubd $18, $17, $4 \n\t"
		"vsubd $19, $17, $5 \n\t"
		"vsubd $20, $17, $6 \n\t"
		"vsubd $21, $17, $7 \n\t"
		// dz*dz
		"vmad $0, $0, $8,  $8  \n\t" "vstd $0, " V24 " \n\t"
		"vmad $1, $1, $9,  $9  \n\t" "vstd $1, " V25 " \n\t"
		"vmad $2, $2, $10, $10 \n\t" "vstd $2, " V26 " \n\t"
		"vmad $3, $3, $11, $11 \n\t" "vstd $3, " V27 " \n\t"
		"vmad $4, $4, $12, $12 \n\t" "vstd $4, " V28 " \n\t"
		"vmad $5, $5, $13, $13 \n\t" "vstd $5, " V29 " \n\t"
		"vmad $6, $6, $14, $14 \n\t" "vstd $6, " V30 " \n\t"
		  "vldd $16, " V1 " \n\t" // mask
		"vmad $7, $7, $15, $15 \n\t" "vstd $7, " V31 " \n\t"
		  "vldd $17, " V2 " \n\t" // magic

		// mask and magic
		// approx rsqrt
		"srlow $8,  1, $0 \n\t" // P0, 2 cycle
		"srlow $9,  1, $1 \n\t"
		"srlow $10, 1, $2 \n\t"
		"srlow $11, 1, $3 \n\t"
		"srlow $12, 1, $4 \n\t"
		"srlow $13, 1, $5 \n\t"
		"srlow $14, 1, $6 \n\t"
		"srlow $15, 1, $7 \n\t"
		"vandw $0, $16, $0 \n\t"
		"vandw $1, $16, $1 \n\t"
		"vandw $2, $16, $2 \n\t"
		"vandw $3, $16, $3 \n\t"
		"vandw $4, $16, $4 \n\t"
		"vandw $5, $16, $5 \n\t"
		"vandw $6, $16, $6 \n\t"
		"vandw $7, $16, $7 \n\t"
		"vsubl $17, $0, $0 \n\t"
		"vsubl $17, $1, $1 \n\t"
		"vsubl $17, $2, $2 \n\t"
		"vsubl $17, $3, $3 \n\t"
		"vsubl $17, $4, $4 \n\t"
		"vsubl $17, $5, $5 \n\t"
		"vsubl $17, $6, $6 \n\t"
		"vsubl $17, $7, $7 \n\t"
		// y^2
		"vmuld $0, $0, $16 \n\t"
		  "vldd $30, " V3 " \n\t" // c0 = 1.0
		"vmuld $1, $1, $17 \n\t"
		  "vldd $28, " V6 " \n\t" // c3 = 5/16
		"vmuld $2, $2, $18 \n\t"
		  "vldd $29, " V7 " \n\t" // c4 = 35/128
		"vmuld $3, $3, $19 \n\t"
		"vmuld $4, $4, $20 \n\t"
		"vmuld $5, $5, $21 \n\t"
		"vmuld $6, $6, $22 \n\t"
		"vmuld $7, $7, $23 \n\t"

		// h = 1.0 - x * y^2
		"vnmad $16, $8,  $30, $16 \n\t"
		"vnmad $17, $9,  $30, $17 \n\t"
		"vnmad $18, $10, $30, $18 \n\t"
		"vnmad $19, $11, $30, $19 \n\t"
		"vnmad $20, $12, $30, $20 \n\t"
		"vnmad $21, $13, $30, $21 \n\t"
		"vnmad $22, $14, $30, $22 \n\t"
		"vnmad $23, $15, $30, $23 \n\t"

		// h*c4 + c3
		"vmad $16, $29, $28, $8  \n\t"
		  "vldd $30, " V5 " \n\t" // c2 = 3/8
		"vmad $17, $29, $28, $9  \n\t"
		"vmad $18, $29, $28, $10 \n\t"
		"vmad $19, $29, $28, $11 \n\t"
		"vmad $20, $29, $28, $12 \n\t"
		"vmad $21, $29, $28, $13 \n\t"
		"vmad $22, $29, $28, $14 \n\t"
		"vmad $23, $29, $28, $15 \n\t"
		  "vldd $28, " V4 " \n\t" // c1 = 1/2

		// h*() + c2
		"vmad $16, $8,  $30, $8  \n\t"
		"vmad $17, $9,  $30, $9  \n\t"
		"vmad $18, $10, $30, $10 \n\t"
		"vmad $19, $11, $30, $11 \n\t"
		"vmad $20, $12, $30, $12 \n\t"
		"vmad $21, $13, $30, $13 \n\t"
		"vmad $22, $14, $30, $14 \n\t"
		"vmad $23, $15, $30, $15 \n\t"
		  "vldd $30, " V3 " \n\t" // c0 = 1.0
		// h*() + c1
		"vmad $16, $8,  $28, $8  \n\t"
		"vmad $17, $9,  $28, $9  \n\t"
		"vmad $18, $10, $28, $10 \n\t"
		"vmad $19, $11, $28, $11 \n\t"
		"vmad $20, $12, $28, $12 \n\t"
		"vmad $21, $13, $28, $13 \n\t"
		"vmad $22, $14, $28, $14 \n\t"
		"vmad $23, $15, $28, $15 \n\t"

		// h*() + c0
		"vmad $16, $8,  $30, $8  \n\t"
		"vmad $17, $9,  $30, $9  \n\t"
		"vmad $18, $10, $30, $10 \n\t"
		"vmad $19, $11, $30, $11 \n\t"
		"vmad $20, $12, $30, $12 \n\t"
		"vmad $21, $13, $30, $13 \n\t"
		"vmad $22, $14, $30, $14 \n\t"
		"vmad $23, $15, $30, $15 \n\t"

		// y *= poly(h)
		"vmuld $0, $8,  $0 \n\t"
		"vmuld $1, $9,  $1 \n\t"
		"vmuld $2, $10, $2 \n\t"
		"vmuld $3, $11, $3 \n\t"
		"vmuld $4, $12, $4 \n\t"
		"vmuld $5, $13, $5 \n\t"
		"vmuld $6, $14, $6 \n\t"
		"vmuld $7, $15, $7 \n\t"

		// shuffle mj
		"ldi $30, 255($31) \n\t"
		"vshff $24, $24, $30, $12 \n\t"
		"vshff $25, $25, $30, $13 \n\t"
		"vshff $26, $26, $30, $14 \n\t"
		"vshff $27, $27, $30, $15 \n\t"
		  "ldl  $30, " S2 " \n\t" // force

		// mj/r
		"vmuld $0, $12, $8  \n\t"
		"vmuld $1, $13, $9  \n\t"
		"vmuld $2, $14, $10 \n\t"
		"vmuld $3, $15, $11 \n\t"
		"vmuld $4, $12, $12  \n\t"
		"vmuld $5, $13, $13  \n\t"
		"vmuld $6, $14, $14 \n\t"
		"vmuld $7, $15, $15 \n\t"
		// ri^2
		"vmuld $0, $0, $0  \n\t"
		"vmuld $1, $1, $1  \n\t"
		"vmuld $2, $2, $2  \n\t"
		  "vldd $16,   0($30) \n\t" // fx0
		"vmuld $3, $3, $3  \n\t"
		  "vldd $17,  32($30) \n\t" // fy0
		"vmuld $4, $4, $4  \n\t"
		  "vldd $18,  64($30) \n\t" // fz0
		"vmuld $5, $5, $5  \n\t"
		  "vldd $19, 128($30) \n\t" // fx1
		"vmuld $6, $6, $6  \n\t"
		  "vldd $20, 160($30) \n\t" // fy1
		"vmuld $7, $7, $7  \n\t"
		  "vldd $21, 192($30) \n\t" // fz1
		// m/r^3
		"vmuld $8,  $0, $0  \n\t"
		  "vldd $8,  " V8 " \n\t" // load dx0, dy0, dz0
		"vmuld $9,  $1, $1  \n\t"
		  "vldd $9,  " V16 " \n\t"
		"vmuld $10, $2, $2  \n\t"
		  "vldd $10, " V24 " \n\t"
		"vmuld $11, $3, $3  \n\t"
		  "vldd $11, " V12 " \n\t" // load dx4, dy4, dz4
		"vmuld $12, $4, $4  \n\t"
		  "vldd $12, " V20 " \n\t"
		"vmuld $13, $5, $5  \n\t"
		  "vldd $13, " V28 " \n\t"
		"vmuld $14, $6, $6  \n\t"
		"vmuld $15, $7, $7  \n\t"
		
		// accumulate
		"vmad $8,  $0, $16, $16 \n\t"
		  "vldd $8,  " V9 " \n\t" // load dx1, dy1, dz1
		"vmad $9,  $0, $17, $17 \n\t"
		  "vldd $9,  " V17 " \n\t"
		"vmad $10, $0, $18, $18 \n\t"
		  "vldd $10, " V25 " \n\t"
		"vmad $11, $4, $19, $19 \n\t"
		  "vldd $11, " V13 " \n\t" // load dx5, dy5, dz5
		"vmad $12, $4, $20, $20 \n\t"
		  "vldd $12, " V21 " \n\t"
		"vmad $13, $4, $21, $21 \n\t"
		  "vldd $13, " V29 " \n\t"

		// accumulate
		"vmad $8,  $1, $16, $16 \n\t"
		  "vldd $8,  " V10 " \n\t" // load dx2, dy2, dz2
		"vmad $9,  $1, $17, $17 \n\t"
		  "vldd $9,  " V18 " \n\t"
		"vmad $10, $1, $18, $18 \n\t"
		  "vldd $10, " V26 " \n\t"
		"vmad $11, $5, $19, $19 \n\t"
		  "vldd $11, " V14 " \n\t" // load dx6, dy6, dz6
		"vmad $12, $5, $20, $20 \n\t"
		  "vldd $12, " V22 " \n\t"
		"vmad $13, $5, $21, $21 \n\t"
		  "vldd $13, " V30 " \n\t"

		// accumulate
		"vmad $8,  $2, $16, $16 \n\t"
		  "vldd $8,  " V11 " \n\t" // load dx3, dy3, dz3
		"vmad $9,  $2, $17, $17 \n\t"
		  "vldd $9,  " V19 " \n\t"
		"vmad $10, $2, $18, $18 \n\t"
		  "vldd $10, " V27 " \n\t"
		"vmad $11, $6, $19, $19 \n\t"
		  "vldd $11, " V15 " \n\t" // load dx7, dy7, dz7
		"vmad $12, $6, $20, $20 \n\t"
		  "vldd $12, " V23 " \n\t"
		"vmad $13, $6, $21, $21 \n\t"
		  "vldd $13, " V31 " \n\t"

		// accumulate
		"vmad $8,  $3, $16, $16 \n\t"
		  "ldl  $30, " S2 " \n\t" // force
		"vmad $9,  $3, $17, $17 \n\t"
		  "ldl  $28, " S0 " \n\t" // epi
		"vmad $10, $3, $18, $18 \n\t"
		  "ldl  $29, " S3 " \n\t" // fend
		"vmad $11, $7, $19, $19 \n\t"
		  "ldi $30, 256($30) \n\t" // force+=2
		"vmad $12, $7, $20, $20 \n\t"
		"vmad $13, $7, $21, $21 \n\t"
		  "ldi $28, 256($28) \n\t" // epi+=2

		"cmpule $29, $30, $0 \n\t"
		"beq $0,.Lt_fdps_4 \n\n\t"

		// store forces
		"vstd $16, -256($30) \n\t" // fx0
		"vstd $17, -224($30) \n\t" // fy0
		"vstd $18, -192($30) \n\t" // fz0
		"vstd $19, -128($30) \n\t" // fx1
		"vstd $20,  -96($30) \n\t" // fy1
		"vstd $21,  -64($30) \n\t" // fz1
	);
	// restore registers
	asm volatile(
		"ldl $26,  8($31)\n\t"
		"ldl $29, 16($31)\n\t"
		"ldl $30, 24($31)\n\t"
		"ldl $15, 40($31)\n\t"
		"ldl $27, 48($31)\n\t"
		"stl $31,  8($31)\n\t"
		"stl $31, 16($31)\n\t"
		"stl $31, 24($31)\n\t"
		"stl $31, 40($31)\n\t"
		"stl $31, 48($31)\n\t"
	);
}
