#ifdef FULL_ASM
// #define RSQRT_FULL_CONVERGENCE

// Memory map
// #define LDM_BASE "1536+"
#define LDM_BASE "0x600+"
#define S0 LDM_BASE "0($31)"
#define S1 LDM_BASE "8($31)"
#define S2 LDM_BASE "16($31)"
#define S3 LDM_BASE "24($31)"

#define V1 LDM_BASE "0x20($31)"
#define V2 LDM_BASE "0x40($31)"
#define V3 LDM_BASE "0x60($31)"
#define V4 LDM_BASE "0x80($31)"
#define V5 LDM_BASE "0xa0($31)"
#define V6 LDM_BASE "0xc0($31)"
#define V7 LDM_BASE "0xe0($31)"

#define DX0 LDM_BASE "0x100($31)"
#define DX1 LDM_BASE "0x120($31)"
#define DX2 LDM_BASE "0x140($31)"
#define DX3 LDM_BASE "0x160($31)"
#define DX4 LDM_BASE "0x180($31)"
#define DX5 LDM_BASE "0x1a0($31)"
#define DX6 LDM_BASE "0x1c0($31)"
#define DX7 LDM_BASE "0x1e0($31)"

#define DY0 LDM_BASE "0x200($31)"
#define DY1 LDM_BASE "0x220($31)"
#define DY2 LDM_BASE "0x240($31)"
#define DY3 LDM_BASE "0x260($31)"
#define DY4 LDM_BASE "0x280($31)"
#define DY5 LDM_BASE "0x2a0($31)"
#define DY6 LDM_BASE "0x2c0($31)"
#define DY7 LDM_BASE "0x2e0($31)"

#define DZ0 LDM_BASE "0x300($31)"
#define DZ1 LDM_BASE "0x320($31)"
#define DZ2 LDM_BASE "0x340($31)"
#define DZ3 LDM_BASE "0x360($31)"
#define DZ4 LDM_BASE "0x380($31)"
#define DZ5 LDM_BASE "0x3a0($31)"
#define DZ6 LDM_BASE "0x3c0($31)"
#define DZ7 LDM_BASE "0x3e0($31)"

#define XJ0  LDM_BASE "0x400($31)"
#define XJ1  LDM_BASE "0x420($31)"
#define XJ2  LDM_BASE "0x440($31)"
#define XJ3  LDM_BASE "0x460($31)"
#define YJ0  LDM_BASE "0x480($31)"
#define YJ1  LDM_BASE "0x4a0($31)"
#define YJ2  LDM_BASE "0x4c0($31)"
#define YJ3  LDM_BASE "0x4e0($31)"

#define ZJ0  LDM_BASE "0x500($31)"
#define ZJ1  LDM_BASE "0x520($31)"
#define ZJ2  LDM_BASE "0x540($31)"
#define ZJ3  LDM_BASE "0x560($31)"
#define MJ0  LDM_BASE "0x580($31)"
#define MJ1  LDM_BASE "0x5a0($31)"
#define MJ2  LDM_BASE "0x5c0($31)"
#define MJ3  LDM_BASE "0x5e0($31)"

#define VXJ0  LDM_BASE "0x600($31)"
#define VXJ1  LDM_BASE "0x620($31)"
#define VXJ2  LDM_BASE "0x640($31)"
#define VXJ3  LDM_BASE "0x660($31)"
#define VYJ0  LDM_BASE "0x680($31)"
#define VYJ1  LDM_BASE "0x6a0($31)"
#define VYJ2  LDM_BASE "0x6c0($31)"
#define VYJ3  LDM_BASE "0x6e0($31)"

#define VZJ0  LDM_BASE "0x700($31)"
#define VZJ1  LDM_BASE "0x720($31)"
#define VZJ2  LDM_BASE "0x740($31)"
#define VZJ3  LDM_BASE "0x760($31)"

#define V8    LDM_BASE "0x780($31)"
#define V9    LDM_BASE "0x7a0($31)"
#define V10   LDM_BASE "0x7c0($31)"
#define V11   LDM_BASE "0x7e0($31)"

//                      2048
#define RV0   LDM_BASE "0x800($31)"
#define RV1   LDM_BASE "0x820($31)"
#define RV2   LDM_BASE "0x840($31)"
#define RV3   LDM_BASE "0x860($31)"
#define RV4   LDM_BASE "0x880($31)"
#define RV5   LDM_BASE "0x8a0($31)"
#define RV6   LDM_BASE "0x8c0($31)"
#define RV7   LDM_BASE "0x8e0($31)"

#define V12   LDM_BASE "0x900($31)"

// satelite force
#define FSX0  LDM_BASE "0x920($31)"
#define FSX1  LDM_BASE "0x940($31)"
#define FSX2  LDM_BASE "0x960($31)"
#define FSX3  LDM_BASE "0x980($31)"
#define FSY0  LDM_BASE "0x9a0($31)"
#define FSY1  LDM_BASE "0x9c0($31)"
#define FSY2  LDM_BASE "0x9e0($31)"
#define FSY3  LDM_BASE "0xa00($31)"
#define FSZ0  LDM_BASE "0xa20($31)"
#define FSZ1  LDM_BASE "0xa40($31)"
#define FSZ2  LDM_BASE "0xa60($31)"
#define FSZ3  LDM_BASE "0xa80($31)"
//                      2688

#ifdef RSQRT_FULL_CONVERGENCE
#define RSQ0  LDM_BASE "0xaa0($31)"
#define RSQ1  LDM_BASE "0xac0($31)"
#define RSQ2  LDM_BASE "0xae0($31)"
#define RSQ3  LDM_BASE "0xb00($31)"
#define RSQ4  LDM_BASE "0xb20($31)"
#define RSQ5  LDM_BASE "0xb40($31)"
#define RSQ6  LDM_BASE "0xb60($31)"
#define RSQ7  LDM_BASE "0xb80($31)"
//                      2944
#endif

                           
                           
typedef double doublev4 __attribute__ ((__mode__(__V4DF__)));
                           
__attribute__((noinline))
void p_init_ldm_constants(void *pcoff){
	asm volatile(
		"vldd $1,  0($16) \n\t"
		"vldd $2, 32($16) \n\t"
		"vldd $3, 64($16) \n\t"
		"vldd $4, 96($16) \n\t"
		"vstd $1, " V9  " \n\t"
		"vstd $2, " V10 " \n\t"
		"vstd $3, " V11 " \n\t"
		"vstd $4, " V12 " \n\t"
			);
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
	long ltiny = 1L << 55;
	asm volatile(
		"vshff %0, %0, $31, %0 \n\t" 
		"vstd %0, " V8
		:: "r"(ltiny)  );

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

void CalcGrav_p(
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
		"addw $16, 7,   $16  \n\t"
		"srl  $16, 3,   $16  \n\t"
		"sll  $16, 8,   $16  \n\t"

		"addl $16, $19, $16  \n\t" // fend
		  "ldi $30,  85($31) \n\t"
		"stl  $17, " S0 " \n\t" // epi
		  "ldi $29, 170($31) \n\t"
		// "stl  $18, " S1 " \n\t" // epj
		"stl  $19, " S2 " \n\t" // force
		  "ldi $28, 255($31) \n\t"
		"stl  $16, " S3 " \n\t" // fend
		// load epj

		"vldd $20,   0($18) \n\t"
		"vldd $21,  32($18) \n\t"
		"vldd $22,  64($18) \n\t"
		"vldd $23,  96($18) \n\t"
		"vldd $24, 128($18) \n\t"
		"vldd $25, 160($18) \n\t"
		"vldd $26, 192($18) \n\t"
		"vldd $27, 224($18) \n\t"

		// shuffle epj
		"vshff $20, $20, $31, $0  \n\t"
		"vshff $20, $20, $30, $1  \n\t"
		  "vstd $0, "  XJ0 "\n\t"
		"vshff $20, $20, $29, $2  \n\t"
		  "vstd $1, "  YJ0 "\n\t"
		"vshff $20, $20, $28, $3  \n\t"
		  "vstd $2, "  ZJ0 "\n\t"
		"vshff $21, $21, $31, $4  \n\t"
		  "vstd $3, "  MJ0 "\n\t"
		"vshff $21, $21, $30, $5  \n\t"
		  "vstd $4, "  VXJ0 "\n\t"
		"vshff $21, $21, $29, $6  \n\t"
		  "vstd $5, "  VYJ0 "\n\t"
		// "vshff $21, $21, $28, $7  \n\t"

		"vshff $22, $22, $31, $8  \n\t"
		  "vstd $6, "  VZJ0 "\n\t"
		"vshff $22, $22, $30, $9  \n\t"
		  "vstd $8, "  XJ1 "\n\t"
		"vshff $22, $22, $29, $10 \n\t"
		  "vstd $9, "  YJ1 "\n\t"
		"vshff $22, $22, $28, $11 \n\t"
		  "vstd $10, " ZJ1 "\n\t"
		"vshff $23, $23, $31, $12 \n\t"
		  "vstd $11, " MJ1 "\n\t"
		"vshff $23, $23, $30, $13 \n\t"
		  "vstd $12, " VXJ1 "\n\t"
		"vshff $23, $23, $29, $14 \n\t"
		  "vstd $13, " VYJ1 "\n\t"
		// "vshff $23, $23, $28, $15 \n\t"


		"vshff $24, $24, $31, $0  \n\t"
		  "vstd $14, " VZJ1 "\n\t"
		"vshff $24, $24, $30, $1  \n\t"
		  "vstd $0, "  XJ2 "\n\t"
		"vshff $24, $24, $29, $2  \n\t"
		  "vstd $1, "  YJ2 "\n\t"
		"vshff $24, $24, $28, $3  \n\t"
		  "vstd $2, "  ZJ2 "\n\t"
		"vshff $25, $25, $31, $4  \n\t"
		  "vstd $3, "  MJ2 "\n\t"
		"vshff $25, $25, $30, $5  \n\t"
		  "vstd $4, "  VXJ2 "\n\t"
		"vshff $25, $25, $29, $6  \n\t"
		  "vstd $5, "  VYJ2 "\n\t"
		// "vshff $25, $25, $28, $7  \n\t"

		"vshff $26, $26, $31, $8  \n\t"
		  "vstd $6, "  VZJ2 "\n\t"
		"vshff $26, $26, $30, $9  \n\t"
		  "vstd $8, "  XJ3 "\n\t"
		"vshff $26, $26, $29, $10 \n\t"
		  "vstd $9, "  YJ3 "\n\t"
		"vshff $26, $26, $28, $11 \n\t"
		  "vstd $10, " ZJ3 "\n\t"
		"vshff $27, $27, $31, $12 \n\t"
		  "vstd $11, " MJ3 "\n\t"
		"vshff $27, $27, $30, $13 \n\t"
		  "vstd $12, " VXJ3 "\n\t"
		"vshff $27, $27, $29, $14 \n\t"
		  "vstd $13, " VYJ3 "\n\t"
		  "vstd $14, " VZJ3 "\n\t"
		// "vshff $27, $27, $28, $15 \n\t"

		"ldi $30, 0($17)    \n\t" // epi
		"vldd $25,   0($30) \n\t" // xi0
		"vldd $26,  32($30) \n\t" // vxi0
	);
	// loop body
	asm volatile(
		".align 6 \n"
		".Lt_fdps_5: \n\t"

		"ldl  $29, " S1 " \n\t" // epj
		"vldd $28, " V8 "\n\t"   // tiny

		/*
		 *  $20:$27 are free
		 */

		// DX,DVX

		"vldd $0, "  XJ0 "\n\t"
		"vldd $1, "  XJ1 "\n\t"
		"vldd $2, "  XJ2 "\n\t"
		"vldd $3, "  XJ3 "\n\t"

		// dx
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VXJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VXJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VXJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VXJ3 "\n\t"
		// dvx
		"vsubd $4,  $26, $4 \n\t"
		"vsubd $5,  $26, $5 \n\t"
		"vsubd $6,  $26, $6 \n\t"
		  "vldd $25,  64($30) \n\t" // yi0
		"vsubd $7,  $26, $7 \n\t"
		  "vldd $26,  96($30) \n\t" // vyi0
		// store dx
		// dx * dx
		"vmad $0, $0, $28, $8  \n\t"
		  "vstd $0, "  DX0 "\n\t"
		"vmad $1, $1, $28, $9  \n\t"
		  "vstd $1, "  DX1 "\n\t"
		"vmad $2, $2, $28, $10  \n\t"
		  "vstd $2, "  DX2 "\n\t"
		"vmad $3, $3, $28, $11  \n\t"
		  "vstd $3, "  DX3 "\n\t"
		// dx * dvx
		"vmuld $0, $4, $16  \n\t"
		  "vldd $0, "  YJ0 "\n\t"
		"vmuld $1, $5, $17  \n\t"
		  "vldd $1, "  YJ1 "\n\t"
		"vmuld $2, $6, $18  \n\t"
		  "vldd $2, "  YJ2 "\n\t"
		"vmuld $3, $7, $19  \n\t"
		  "vldd $3, "  YJ3 "\n\t"

		// DY,DVY

		// dy
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VYJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VYJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VYJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VYJ3 "\n\t"
		// dvy
		"vsubd $4,  $26, $4 \n\t"
		"vsubd $5,  $26, $5 \n\t"
		"vsubd $6,  $26, $6 \n\t"
		  "vldd $25, 128($30) \n\t" // zi0
		"vsubd $7,  $26, $7 \n\t"
		  "vldd $26, 160($30) \n\t" // vzi0
		// store dy
		// dy * dy
		"vmad $0, $0, $8,  $8  \n\t"
		  "vstd $0, "  DY0 "\n\t"
		"vmad $1, $1, $9,  $9  \n\t"
		  "vstd $1, "  DY1 "\n\t"
		"vmad $2, $2, $10, $10  \n\t"
		  "vstd $2, "  DY2 "\n\t"
		"vmad $3, $3, $11, $11  \n\t"
		  "vstd $3, "  DY3 "\n\t"
		// dy * dvy
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $0, "  ZJ0 "\n\t"
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $1, "  ZJ1 "\n\t"
		"vmad $2, $6, $18, $18  \n\t"
		  "vldd $2, "  ZJ2 "\n\t"
		"vmad $3, $7, $19, $19  \n\t"
		  "vldd $3, "  ZJ3 "\n\t"

		// DZ,DVZ

		// dz
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VZJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VZJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VZJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VZJ3 "\n\t"
		// dvz
		"vsubd $4,  $26, $4 \n\t"
		  "vldd $27, 320($30) \n\t" // yi0
		"vsubd $5,  $26, $5 \n\t"
		  "vldd $29, 352($30) \n\t" // yxi0
		"vsubd $6,  $26, $6 \n\t"
		  "vldd $25, 256($30) \n\t" // xi1
		"vsubd $7,  $26, $7 \n\t"
		  "vldd $26, 288($30) \n\t" // vxi1
		// store dz
		// dz * dz
		"vmad $0, $0, $8,  $8  \n\t"
		  "vstd $0, "  DZ0 "\n\t"
		"vmad $1, $1, $9,  $9  \n\t"
		  "vstd $1, "  DZ1 "\n\t"
		"vmad $2, $2, $10, $10  \n\t"
		  "vstd $2, "  DZ2 "\n\t"
		"vmad $3, $3, $11, $11  \n\t"
		  "vstd $3, "  DZ3 "\n\t"
		// dz * dvz
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $0, "  XJ0 "\n\t"
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $1, "  XJ1 "\n\t"
		"vmad $2, $6, $18, $18  \n\t"
		  "vldd $2, "  XJ2 "\n\t"
		"vmad $3, $7, $19, $19  \n\t"
		  "vldd $3, "  XJ3 "\n\t"
		// store rdotv

		// DX,DVX

		// dx
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VXJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VXJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VXJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VXJ3 "\n\t"
		// dvx
		"vsubd $4,  $26, $4 \n\t"
		  "vstd $16, "  RV0 "\n\t"
		"vsubd $5,  $26, $5 \n\t"
		  "vstd $17, "  RV1 "\n\t"
		"vsubd $6,  $26, $6 \n\t"
		  "vstd $18, "  RV2 "\n\t"
		"vsubd $7,  $26, $7 \n\t"
		  "vstd $19, "  RV3 "\n\t"
		// store dx
		// dx * dx
		"vmad $0, $0, $28, $12  \n\t"
		  "vstd $0, "  DX4 "\n\t"
		"vmad $1, $1, $28, $13  \n\t"
		  "vstd $1, "  DX5 "\n\t"
		"vmad $2, $2, $28, $14  \n\t"
		  "vstd $2, "  DX6 "\n\t"
		"vmad $3, $3, $28, $15  \n\t"
		  "vstd $3, "  DX7 "\n\t"
		// dx * dvx
		"vmuld $0, $4, $16  \n\t"
		  "vldd $0, "  YJ0 "\n\t"
		"vmuld $1, $5, $17  \n\t"
		  "vldd $1, "  YJ1 "\n\t"
		"vmuld $2, $6, $18  \n\t"
		  "vldd $2, "  YJ2 "\n\t"
		"vmuld $3, $7, $19  \n\t"
		  "vldd $3, "  YJ3 "\n\t"

		// DY,DVY

		// dy
		"vsubd $0,  $27, $0 \n\t"
		  "vldd $4, "  VYJ0 "\n\t"
		"vsubd $1,  $27, $1 \n\t"
		  "vldd $5, "  VYJ1 "\n\t"
		"vsubd $2,  $27, $2 \n\t"
		  "vldd $6, "  VYJ2 "\n\t"
		"vsubd $3,  $27, $3 \n\t"
		  "vldd $7, "  VYJ3 "\n\t"
		// dvy
		"vsubd $4,  $29, $4 \n\t"
		"vsubd $5,  $29, $5 \n\t"
		"vsubd $6,  $29, $6 \n\t"
		  "vldd $25, 384($30) \n\t" // zi0
		"vsubd $7,  $29, $7 \n\t"
		  "vldd $26, 416($30) \n\t" // zxi0
		// store dy
		// dy * dy
		"vmad $0, $0, $12, $12 \n\t"
		  "vstd $0, "  DY4 "\n\t"
		"vmad $1, $1, $13, $13 \n\t"
		  "vstd $1, "  DY5 "\n\t"
		"vmad $2, $2, $14, $14 \n\t"
		  "vstd $2, "  DY6 "\n\t"
		"vmad $3, $3, $15, $15 \n\t"
		  "vstd $3, "  DY7 "\n\t"
		// dy * dvy
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $0, "  ZJ0 "\n\t"
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $1, "  ZJ1 "\n\t"
		"vmad $2, $6, $18, $18  \n\t"
		  "vldd $2, "  ZJ2 "\n\t"
		"vmad $3, $7, $19, $19  \n\t"
		  "vldd $3, "  ZJ3 "\n\t"

		// DZ,DVZ

		// dz
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VZJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VZJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VZJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VZJ3 "\n\t"
		// dvz
		"vsubd $4,  $26, $4 \n\t"
		"vsubd $5,  $26, $5 \n\t"
		"vsubd $6,  $26, $6 \n\t"
		"vsubd $7,  $26, $7 \n\t"
		// store dz
		// dz * dz
		"vmad $0, $0, $12, $12 \n\t"
		  "vstd $0, "  DZ4 "\n\t"
		"vmad $1, $1, $13, $13 \n\t"
		  "vstd $1, "  DZ5 "\n\t"
		"vmad $2, $2, $14, $14 \n\t"
		  "vstd $2, "  DZ6 "\n\t"
		"vmad $3, $3, $15, $15 \n\t"
		  "vstd $3, "  DZ7 "\n\t"
		// dz * dvz
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $20, " V1 " \n\t" // mask
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $21, " V2 " \n\t" // magic
		"vmad $2, $6, $18, $18  \n\t"
		"vmad $3, $7, $19, $19  \n\t"

		// store rdotv

		// RSQRT-X5
		"srlow $8,  1, $0 \n\t" // P0, 2 cycle
		"srlow $9,  1, $1 \n\t"
		"srlow $10, 1, $2 \n\t"
		"srlow $11, 1, $3 \n\t"
		"srlow $12, 1, $4 \n\t"
		  "vstd $16, "  RV4 "\n\t"
		"srlow $13, 1, $5 \n\t"
		  "vstd $17, "  RV5 "\n\t"
		"srlow $14, 1, $6 \n\t"
		  "vstd $18, "  RV6 "\n\t"
		"srlow $15, 1, $7 \n\t"
		  "vstd $19, "  RV7 "\n\t"
		"vandw $0, $20, $0 \n\t"
		"vandw $1, $20, $1 \n\t"
		"vandw $2, $20, $2 \n\t"
		"vandw $3, $20, $3 \n\t"
		"vandw $4, $20, $4 \n\t"
		"vandw $5, $20, $5 \n\t"
		"vandw $6, $20, $6 \n\t"
		"vandw $7, $20, $7 \n\t"
		"vsubl $21, $0, $0 \n\t"
		"vsubl $21, $1, $1 \n\t"
		"vsubl $21, $2, $2 \n\t"
		"vsubl $21, $3, $3 \n\t"
		"vsubl $21, $4, $4 \n\t"
		"vsubl $21, $5, $5 \n\t"
		"vsubl $21, $6, $6 \n\t"
		"vsubl $21, $7, $7 \n\t"
#ifdef RSQRT_FULL_CONVERGENCE
		// store rsq
		  "vstd $8, "   RSQ0 "\n\t"
		  "vstd $9, "   RSQ1 "\n\t"
		  "vstd $10, "  RSQ2 "\n\t"
		  "vstd $11, "  RSQ3 "\n\t"
		  "vstd $12, "  RSQ4 "\n\t"
		  "vstd $13, "  RSQ5 "\n\t"
		  "vstd $14, "  RSQ6 "\n\t"
		  "vstd $15, "  RSQ7 "\n\t"
#endif
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
		  "vldd $12, "  MJ0 "\n\t"
		"vmuld $5, $13, $5 \n\t"
		  "vldd $13, "  MJ1 "\n\t"
		"vmuld $6, $14, $6 \n\t"
		  "vldd $14, "  MJ2 "\n\t"
		"vmuld $7, $15, $7 \n\t"
		  "vldd $15, "  MJ3 "\n\t"

#ifdef RSQRT_FULL_CONVERGENCE // for strict accuracy check
		// constants
		  "vldd $30, " V3 " \n\t" // c0 = 1.0
		  "vldd $29, " V4 " \n\t" // c1 = 1/2
		  "vldd $28, " V5 " \n\t" // c2 = 3/8
		// load x=rsq
		  "vldd $8, "   RSQ0 "\n\t"
		  "vldd $9, "   RSQ1 "\n\t"
		  "vldd $10, "  RSQ2 "\n\t"
		  "vldd $11, "  RSQ3 "\n\t"
		  "vldd $12, "  RSQ4 "\n\t"
		  "vldd $13, "  RSQ5 "\n\t"
		  "vldd $14, "  RSQ6 "\n\t"
		  "vldd $15, "  RSQ7 "\n\t"
		// y^2
		"vmuld $0, $0, $16 \n\t"
		"vmuld $1, $1, $17 \n\t"
		"vmuld $2, $2, $18 \n\t"
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
		// p = (1/2)h + (3/8)
		"vmad $16, $28, $29, $8   \n\t"
		"vmad $17, $28, $29, $9   \n\t"
		"vmad $18, $28, $29, $10  \n\t"
		"vmad $19, $28, $29, $11  \n\t"
		"vmad $20, $28, $29, $12  \n\t"
		"vmad $21, $28, $29, $13  \n\t"
		"vmad $22, $28, $29, $14  \n\t"
		"vmad $23, $28, $29, $15  \n\t"
		// y * h
		"vmuld $16, $0, $16 \n\t"
		"vmuld $17, $1, $17 \n\t"
		"vmuld $18, $2, $18 \n\t"
		"vmuld $19, $3, $19 \n\t"
		"vmuld $20, $4, $20 \n\t"
		"vmuld $21, $5, $21 \n\t"
		"vmuld $22, $6, $22 \n\t"
		"vmuld $23, $7, $23 \n\t"
		// y += (y*h)*p
		"vmad $16, $8,  $0, $0   \n\t"
		"vmad $17, $9,  $1, $1   \n\t"
		"vmad $18, $10, $2, $2   \n\t"
		"vmad $19, $11, $3, $3   \n\t"
		"vmad $20, $12, $4, $4   \n\t"
		  "vldd $12, "  MJ0 "\n\t"
		"vmad $21, $13, $5, $5   \n\t"
		  "vldd $13, "  MJ1 "\n\t"
		"vmad $22, $14, $6, $6   \n\t"
		  "vldd $14, "  MJ2 "\n\t"
		"vmad $23, $15, $7, $7   \n\t"
		  "vldd $15, "  MJ3 "\n\t"
#endif

		// a0, a1, b0, b1, c, d for $0:$23 (4 for each)
		// load mj
		// load constants

		// b = mj *a;
		"vmuld $12, $0,  $8 \n\t"
		  "vldd $24, "  V3  "\n\t" // 1.0
		"vmuld $13, $1,  $9 \n\t"
		  "vldd $28, "  V12 "\n\t" // r_collv4
		"vmuld $14, $2, $10 \n\t"
		  "vldd $25, "  V9  "\n\t" // m_r_coll_inv_qv4
		"vmuld $15, $3, $11 \n\t"
		  "vldd $26, "  V10 "\n\t" // k_v_miv4
		// c = a * a;
		"vmuld $0, $0, $16 \n\t"
		  "vldd $27, "  V11 "\n\t" // eta_v_miv4
		"vmuld $1, $1, $17 \n\t"
		"vmuld $2, $2, $18 \n\t"
		"vmuld $3, $3, $19 \n\t"
		// d = 1.0 - r_collv4 * a;
		"vnmad $0, $28, $24, $20 \n\t"
		  "vldd $0, " RV0 "\n\t"
		"vnmad $1, $28, $24, $21 \n\t"
		  "vldd $1, " RV1 "\n\t"
		"vnmad $2, $28, $24, $22 \n\t"
		  "vldd $2, " RV2 "\n\t"
		"vnmad $3, $28, $24, $23 \n\t"
		  "vldd $3, " RV3 "\n\t"
		// a = rdotv;
		// b *= c;
		"vmuld $8,  $16, $8  \n\t"
		"vmuld $9,  $17, $9  \n\t"
		"vmuld $10, $18, $10 \n\t"
		"vmuld $11, $19, $11 \n\t"
		// a *= c;
		"vmuld $0,  $16, $0  \n\t"
		"vmuld $1,  $17, $1  \n\t"
		"vmuld $2,  $18, $2  \n\t"
		"vmuld $3,  $19, $3  \n\t"
		// c = m_r_coll_inv_qv4 + k_v_miv4 * d;
		"vmad $20, $26, $25, $16 \n\t"
		"vmad $21, $26, $25, $17 \n\t"
		"vmad $22, $26, $25, $18 \n\t"
		"vmad $23, $26, $25, $19 \n\t"

			// b' = mj *a';
			"vmuld $12, $4, $12 \n\t"
			"vmuld $13, $5, $13 \n\t"
			"vmuld $14, $6, $14 \n\t"
			"vmuld $15, $7, $15 \n\t"

		// c += eta_v_miv4 * a;
		"vmad $27, $0, $16, $16 \n\t"
		"vmad $27, $1, $17, $17 \n\t"
		"vmad $27, $2, $18, $18 \n\t"
		"vmad $27, $3, $19, $19 \n\t"

			// c' = a' * a';
			"vmuld $4, $4, $0 \n\t"
			"vmuld $5, $5, $1 \n\t"
			"vmuld $6, $6, $2 \n\t"
			"vmuld $7, $7, $3 \n\t"

		// b = simd_vsellt(d, c, b);
		"vsellt $20, $16, $8,  $8  \n\t"
		"vsellt $21, $17, $9,  $9  \n\t"
		"vsellt $22, $18, $10, $10 \n\t"
		"vsellt $23, $19, $11, $11 \n\t"

		// DONE: c' for $0:$3 improves scheduling
		// d = 1.0 - r_collv4 * a';
		"vnmad $4, $28, $24, $20 \n\t"
		  "vldd $4, " RV4 "\n\t"
		"vnmad $5, $28, $24, $21 \n\t"
		  "vldd $5, " RV5 "\n\t"
		"vnmad $6, $28, $24, $22 \n\t"
		  "vldd $6, " RV6 "\n\t"
		"vnmad $7, $28, $24, $23 \n\t"
		  "vldd $7, " RV7 "\n\t"
		// a' = rdotv;
		// b' *= c';
		"vmuld $12, $0, $12  \n\t"
		    "ldl  $28, " S2 " \n\t" // force
		"vmuld $13, $1, $13  \n\t"
		"vmuld $14, $2, $14 \n\t"
		"vmuld $15, $3, $15 \n\t"
		// a' *= c';
		"vmuld $4,  $0, $4  \n\t"
		    "vldd $0,   " DX0 " \n\t" 
		"vmuld $5,  $1, $5  \n\t"
		    "vldd $1,    0($28) \n\t" // fx0
		"vmuld $6,  $2, $6  \n\t"
		    "vldd $2,   " DY0 " \n\t" 
		"vmuld $7,  $3, $7  \n\t"
		    "vldd $3,   32($28) \n\t" // fy0
		// c = m_r_coll_inv_qv4 + k_v_miv4 * d;
		"vmad $20, $26, $25, $16 \n\t"
		"vmad $21, $26, $25, $17 \n\t"
		"vmad $22, $26, $25, $18 \n\t"
		"vmad $23, $26, $25, $19 \n\t"

		// BUBBLE
		"vmad $8,  $0,  $1,  $1 \n\t"
		  "vldd $0,   " DX1 " \n\t" 
		"vmad $8,  $2,  $3,  $3 \n\t"
		  "vldd $2,   " DY1 " \n\t" 

		// c += eta_v_miv4 * a';
		"vmad $27, $4, $16, $16 \n\t"
		    "vldd $4,   " DZ0 " \n\t" 
		"vmad $27, $5, $17, $17 \n\t"
		    "vldd $5,   64($28) \n\t" // fz0
		"vmad $27, $6, $18, $18 \n\t"
		    "vldd $6,   " DX4 " \n\t" 
		"vmad $27, $7, $19, $19 \n\t"
		    "vldd $7,  128($28) \n\t" // fx1

		// BUBBLE
		"vmad $8,  $4,  $5,  $5 \n\t"
		  "vldd $4,   " DZ1 " \n\t" 

		// b' = simd_vsellt(d, c, b');
		"vsellt $20, $16, $12, $12  \n\t"
		    "vldd $20,  " DY4 " \n\t" 
		"vsellt $21, $17, $13, $13  \n\t"
		    "vldd $21, 160($28) \n\t" // fy1
		"vsellt $22, $18, $14, $14 \n\t"
		    "vldd $22,  " DZ4 " \n\t" 
		"vsellt $23, $19, $15, $15 \n\t"
		    "vldd $23, 192($28) \n\t" // fz1

		/*
		 * $0:  dx[0-3]
		 * $1:  fx0
		 * $2:  dy[0-3]
		 * $3:  fy0
		 * $4:  dz[0-3]
		 * $5:  fz0
		 * $6:  dx[4-7]
		 * $7:  fx1
		 * $20: dy[4-7]
		 * $21: fy1
		 * $22: dz[4-7]
		 * $23: fz1
		 */


		"vmad $12, $6,  $7,  $7 \n\t"
		  "vldd $6,   " DX5 " \n\t" 
		"vmad $12, $20, $21, $21 \n\t"
		  "vldd $20,  " DY5 " \n\t" 
		"vmad $12, $22, $23, $23 \n\t"
		  "vldd $22,  " DZ5 " \n\t" 


		"vmad $9,  $0,  $1,  $1 \n\t"
		  "vldd $0,   " DX2 " \n\t" 
		"vmad $9,  $2,  $3,  $3 \n\t"
		  "vldd $2,   " DY2 " \n\t" 
		"vmad $9,  $4,  $5,  $5 \n\t"
		  "vldd $4,   " DZ2 " \n\t" 
		"vmad $13, $6,  $7,  $7 \n\t"
		  "vldd $6,   " DX6 " \n\t" 
		"vmad $13, $20, $21, $21 \n\t"
		  "vldd $20,  " DY6 " \n\t" 
		"vmad $13, $22, $23, $23 \n\t"
		  "vldd $22,  " DZ6 " \n\t" 

		"ldl  $29, " S3 " \n\t" // fend
		"ldl  $30, " S0 " \n\t" // epi

		"vmad $10, $0,  $1,  $1 \n\t"
		  "vldd $0,   " DX3 " \n\t" 
		"vmad $10, $2,  $3,  $3 \n\t"
		  "vldd $2,   " DY3 " \n\t" 
		"vmad $10, $4,  $5,  $5 \n\t"
		  "vldd $4,   " DZ3 " \n\t" 
		"vmad $14, $6,  $7,  $7 \n\t"
		  "vldd $6,   " DX7 " \n\t" 
		"vmad $14, $20, $21, $21 \n\t"
		  "vldd $20,  " DY7 " \n\t" 
		"vmad $14, $22, $23, $23 \n\t"
		  "vldd $22,  " DZ7 " \n\t" 

		"vmad $11, $0,  $1,  $1 \n\t"
		  "ldi $28, 256($28) \n\t" // force+=2
		"vmad $11, $2,  $3,  $3 \n\t"
		  "ldi $30, 512($30) \n\t" // epi+=2
		"vmad $11, $4,  $5,  $5 \n\t"
		  "stl  $28, " S2 " \n\t"  // force
		"vmad $15, $6,  $7,  $7 \n\t"
		  "stl  $30, " S0 " \n\t"  // epi
		"vmad $15, $20, $21, $21 \n\t"
		  "vldd $25,   0($30) \n\t" // xi0
		"vmad $15, $22, $23, $23 \n\t"
		  "vldd $26,  32($30) \n\t" // vxi0

		// store forces

		"vstd $1,    0-256($28) \n\t"
		"vstd $3,   32-256($28) \n\t"
		"vstd $5,   64-256($28) \n\t"
		"vstd $7,  128-256($28) \n\t"
		"vstd $21, 160-256($28) \n\t"
		"vstd $23, 192-256($28) \n\t"

		"cmpule $29, $28, $0 \n\t"
		"beq $0,.Lt_fdps_5 \n\n\t"
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

void clear_satellite_force()
{
	asm volatile(
		"vstd $31, " FSX0 "\n\t"
		"vstd $31, " FSX1 "\n\t"
		"vstd $31, " FSX2 "\n\t"
		"vstd $31, " FSX3 "\n\t"
		"vstd $31, " FSY0 "\n\t"
		"vstd $31, " FSY1 "\n\t"
		"vstd $31, " FSY2 "\n\t"
		"vstd $31, " FSY3 "\n\t"
		"vstd $31, " FSZ0 "\n\t"
		"vstd $31, " FSZ1 "\n\t"
		"vstd $31, " FSZ2 "\n\t"
		"vstd $31, " FSZ3 "\n\t"
	);
}

void read_satellite_force(
		/*
		ForceMM4    *force,   // double ax[4], ay[4], az[4], pot[4];
		*/
		)
{
	asm volatile(
		"vldd $1,  " FSX0 "\n\t"
		"vldd $2,  " FSY0 "\n\t"
		"vldd $3,  " FSZ0 "\n\t"
		"vldd $4,  " FSX1 "\n\t"
		"vldd $5,  " FSY1 "\n\t"
		"vldd $6,  " FSZ1 "\n\t"
		"vldd $7,  " FSX2 "\n\t"
		"vldd $8,  " FSY2 "\n\t"
		"vldd $9,  " FSZ2 "\n\t"
		"vldd $10, " FSX3 "\n\t"
		"vldd $11, " FSY3 "\n\t"
		"vldd $12, " FSZ3 "\n\t"

		"vstd $1,    0($16) \n\t"
		"vstd $2,   32($16) \n\t"
		"vstd $3,   64($16) \n\t"
		"vstd $31,  96($16) \n\t"
		"vstd $4,  128($16) \n\t"
		"vstd $5,  160($16) \n\t"
		"vstd $6,  192($16) \n\t"
		"vstd $31, 224($16) \n\t"
		"vstd $7,  256($16) \n\t"
		"vstd $8,  288($16) \n\t"
		"vstd $9,  320($16) \n\t"
		"vstd $31, 352($16) \n\t"
		"vstd $10, 384($16) \n\t"
		"vstd $11, 416($16) \n\t"
		"vstd $12, 448($16) \n\t"
		"vstd $31, 480($16) \n\t"
	);
}

#if 1

// r 2890
void CalcGrav_ps(
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
		"addw $16, 7,   $16  \n\t"
		"srl  $16, 3,   $16  \n\t"
		"sll  $16, 8,   $16  \n\t"

		"addl $16, $19, $16  \n\t" // fend
		  "ldi $30,  85($31) \n\t"
		"stl  $17, " S0 " \n\t" // epi
		  "ldi $29, 170($31) \n\t"
		// "stl  $18, " S1 " \n\t" // epj
		"stl  $19, " S2 " \n\t" // force
		  "ldi $28, 255($31) \n\t"
		"stl  $16, " S3 " \n\t" // fend
		// load epj

		"vldd $20,   0($18) \n\t"
		"vldd $21,  32($18) \n\t"
		"vldd $22,  64($18) \n\t"
		"vldd $23,  96($18) \n\t"
		"vldd $24, 128($18) \n\t"
		"vldd $25, 160($18) \n\t"
		"vldd $26, 192($18) \n\t"
		"vldd $27, 224($18) \n\t"

		// shuffle epj
		"vshff $20, $20, $31, $0  \n\t"
		"vshff $20, $20, $30, $1  \n\t"
		  "vstd $0, "  XJ0 "\n\t"
		"vshff $20, $20, $29, $2  \n\t"
		  "vstd $1, "  YJ0 "\n\t"
		"vshff $20, $20, $28, $3  \n\t"
		  "vstd $2, "  ZJ0 "\n\t"
		"vshff $21, $21, $31, $4  \n\t"
		  "vstd $3, "  MJ0 "\n\t"
		"vshff $21, $21, $30, $5  \n\t"
		  "vstd $4, "  VXJ0 "\n\t"
		"vshff $21, $21, $29, $6  \n\t"
		  "vstd $5, "  VYJ0 "\n\t"
		// "vshff $21, $21, $28, $7  \n\t"

		"vshff $22, $22, $31, $8  \n\t"
		  "vstd $6, "  VZJ0 "\n\t"
		"vshff $22, $22, $30, $9  \n\t"
		  "vstd $8, "  XJ1 "\n\t"
		"vshff $22, $22, $29, $10 \n\t"
		  "vstd $9, "  YJ1 "\n\t"
		"vshff $22, $22, $28, $11 \n\t"
		  "vstd $10, " ZJ1 "\n\t"
		"vshff $23, $23, $31, $12 \n\t"
		  "vstd $11, " MJ1 "\n\t"
		"vshff $23, $23, $30, $13 \n\t"
		  "vstd $12, " VXJ1 "\n\t"
		"vshff $23, $23, $29, $14 \n\t"
		  "vstd $13, " VYJ1 "\n\t"
		// "vshff $23, $23, $28, $15 \n\t"


		"vshff $24, $24, $31, $0  \n\t"
		  "vstd $14, " VZJ1 "\n\t"
		"vshff $24, $24, $30, $1  \n\t"
		  "vstd $0, "  XJ2 "\n\t"
		"vshff $24, $24, $29, $2  \n\t"
		  "vstd $1, "  YJ2 "\n\t"
		"vshff $24, $24, $28, $3  \n\t"
		  "vstd $2, "  ZJ2 "\n\t"
		"vshff $25, $25, $31, $4  \n\t"
		  "vstd $3, "  MJ2 "\n\t"
		"vshff $25, $25, $30, $5  \n\t"
		  "vstd $4, "  VXJ2 "\n\t"
		"vshff $25, $25, $29, $6  \n\t"
		  "vstd $5, "  VYJ2 "\n\t"
		// "vshff $25, $25, $28, $7  \n\t"

		"vshff $26, $26, $31, $8  \n\t"
		  "vstd $6, "  VZJ2 "\n\t"
		"vshff $26, $26, $30, $9  \n\t"
		  "vstd $8, "  XJ3 "\n\t"
		"vshff $26, $26, $29, $10 \n\t"
		  "vstd $9, "  YJ3 "\n\t"
		"vshff $26, $26, $28, $11 \n\t"
		  "vstd $10, " ZJ3 "\n\t"
		"vshff $27, $27, $31, $12 \n\t"
		  "vstd $11, " MJ3 "\n\t"
		"vshff $27, $27, $30, $13 \n\t"
		  "vstd $12, " VXJ3 "\n\t"
		"vshff $27, $27, $29, $14 \n\t"
		  "vstd $13, " VYJ3 "\n\t"
		  "vstd $14, " VZJ3 "\n\t"
		// "vshff $27, $27, $28, $15 \n\t"

		"ldi $30, 0($17)    \n\t" // epi
		"vldd $25,   0($30) \n\t" // xi0
		"vldd $26,  32($30) \n\t" // vxi0
	);
	// loop body
	asm volatile(
		".align 6 \n"
		".Lt_fdps_6: \n\t"

		"ldl  $29, " S1 " \n\t" // epj
		"vldd $28, " V8 "\n\t"   // tiny

		/*
		 *  $20:$27 are free
		 */

		// DX,DVX

		"vldd $0, "  XJ0 "\n\t"
		"vldd $1, "  XJ1 "\n\t"
		"vldd $2, "  XJ2 "\n\t"
		"vldd $3, "  XJ3 "\n\t"

		// dx
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VXJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VXJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VXJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VXJ3 "\n\t"
		// dvx
		"vsubd $4,  $26, $4 \n\t"
		"vsubd $5,  $26, $5 \n\t"
		"vsubd $6,  $26, $6 \n\t"
		  "vldd $25,  64($30) \n\t" // yi0
		"vsubd $7,  $26, $7 \n\t"
		  "vldd $26,  96($30) \n\t" // vyi0
		// store dx
		// dx * dx
		"vmad $0, $0, $28, $8  \n\t"
		  "vstd $0, "  DX0 "\n\t"
		"vmad $1, $1, $28, $9  \n\t"
		  "vstd $1, "  DX1 "\n\t"
		"vmad $2, $2, $28, $10  \n\t"
		  "vstd $2, "  DX2 "\n\t"
		"vmad $3, $3, $28, $11  \n\t"
		  "vstd $3, "  DX3 "\n\t"
		// dx * dvx
		"vmuld $0, $4, $16  \n\t"
		  "vldd $0, "  YJ0 "\n\t"
		"vmuld $1, $5, $17  \n\t"
		  "vldd $1, "  YJ1 "\n\t"
		"vmuld $2, $6, $18  \n\t"
		  "vldd $2, "  YJ2 "\n\t"
		"vmuld $3, $7, $19  \n\t"
		  "vldd $3, "  YJ3 "\n\t"

		// DY,DVY

		// dy
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VYJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VYJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VYJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VYJ3 "\n\t"
		// dvy
		"vsubd $4,  $26, $4 \n\t"
		"vsubd $5,  $26, $5 \n\t"
		"vsubd $6,  $26, $6 \n\t"
		  "vldd $25, 128($30) \n\t" // zi0
		"vsubd $7,  $26, $7 \n\t"
		  "vldd $26, 160($30) \n\t" // vzi0
		// store dy
		// dy * dy
		"vmad $0, $0, $8,  $8  \n\t"
		  "vstd $0, "  DY0 "\n\t"
		"vmad $1, $1, $9,  $9  \n\t"
		  "vstd $1, "  DY1 "\n\t"
		"vmad $2, $2, $10, $10  \n\t"
		  "vstd $2, "  DY2 "\n\t"
		"vmad $3, $3, $11, $11  \n\t"
		  "vstd $3, "  DY3 "\n\t"
		// dy * dvy
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $0, "  ZJ0 "\n\t"
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $1, "  ZJ1 "\n\t"
		"vmad $2, $6, $18, $18  \n\t"
		  "vldd $2, "  ZJ2 "\n\t"
		"vmad $3, $7, $19, $19  \n\t"
		  "vldd $3, "  ZJ3 "\n\t"

		// DZ,DVZ

		// dz
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VZJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VZJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VZJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VZJ3 "\n\t"
		// dvz
		"vsubd $4,  $26, $4 \n\t"
		  "vldd $27, 320($30) \n\t" // yi0
		"vsubd $5,  $26, $5 \n\t"
		  "vldd $29, 352($30) \n\t" // yxi0
		"vsubd $6,  $26, $6 \n\t"
		  "vldd $25, 256($30) \n\t" // xi1
		"vsubd $7,  $26, $7 \n\t"
		  "vldd $26, 288($30) \n\t" // vxi1
		// store dz
		// dz * dz
		"vmad $0, $0, $8,  $8  \n\t"
		  "vstd $0, "  DZ0 "\n\t"
		"vmad $1, $1, $9,  $9  \n\t"
		  "vstd $1, "  DZ1 "\n\t"
		"vmad $2, $2, $10, $10  \n\t"
		  "vstd $2, "  DZ2 "\n\t"
		"vmad $3, $3, $11, $11  \n\t"
		  "vstd $3, "  DZ3 "\n\t"
		// dz * dvz
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $0, "  XJ0 "\n\t"
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $1, "  XJ1 "\n\t"
		"vmad $2, $6, $18, $18  \n\t"
		  "vldd $2, "  XJ2 "\n\t"
		"vmad $3, $7, $19, $19  \n\t"
		  "vldd $3, "  XJ3 "\n\t"
		// store rdotv

		// DX,DVX

		// dx
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VXJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VXJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VXJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VXJ3 "\n\t"
		// dvx
		"vsubd $4,  $26, $4 \n\t"
		  "vstd $16, "  RV0 "\n\t"
		"vsubd $5,  $26, $5 \n\t"
		  "vstd $17, "  RV1 "\n\t"
		"vsubd $6,  $26, $6 \n\t"
		  "vstd $18, "  RV2 "\n\t"
		"vsubd $7,  $26, $7 \n\t"
		  "vstd $19, "  RV3 "\n\t"
		// store dx
		// dx * dx
		"vmad $0, $0, $28, $12  \n\t"
		  "vstd $0, "  DX4 "\n\t"
		"vmad $1, $1, $28, $13  \n\t"
		  "vstd $1, "  DX5 "\n\t"
		"vmad $2, $2, $28, $14  \n\t"
		  "vstd $2, "  DX6 "\n\t"
		"vmad $3, $3, $28, $15  \n\t"
		  "vstd $3, "  DX7 "\n\t"
		// dx * dvx
		"vmuld $0, $4, $16  \n\t"
		  "vldd $0, "  YJ0 "\n\t"
		"vmuld $1, $5, $17  \n\t"
		  "vldd $1, "  YJ1 "\n\t"
		"vmuld $2, $6, $18  \n\t"
		  "vldd $2, "  YJ2 "\n\t"
		"vmuld $3, $7, $19  \n\t"
		  "vldd $3, "  YJ3 "\n\t"

		// DY,DVY

		// dy
		"vsubd $0,  $27, $0 \n\t"
		  "vldd $4, "  VYJ0 "\n\t"
		"vsubd $1,  $27, $1 \n\t"
		  "vldd $5, "  VYJ1 "\n\t"
		"vsubd $2,  $27, $2 \n\t"
		  "vldd $6, "  VYJ2 "\n\t"
		"vsubd $3,  $27, $3 \n\t"
		  "vldd $7, "  VYJ3 "\n\t"
		// dvy
		"vsubd $4,  $29, $4 \n\t"
		"vsubd $5,  $29, $5 \n\t"
		"vsubd $6,  $29, $6 \n\t"
		  "vldd $25, 384($30) \n\t" // zi0
		"vsubd $7,  $29, $7 \n\t"
		  "vldd $26, 416($30) \n\t" // zxi0
		// store dy
		// dy * dy
		"vmad $0, $0, $12, $12 \n\t"
		  "vstd $0, "  DY4 "\n\t"
		"vmad $1, $1, $13, $13 \n\t"
		  "vstd $1, "  DY5 "\n\t"
		"vmad $2, $2, $14, $14 \n\t"
		  "vstd $2, "  DY6 "\n\t"
		"vmad $3, $3, $15, $15 \n\t"
		  "vstd $3, "  DY7 "\n\t"
		// dy * dvy
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $0, "  ZJ0 "\n\t"
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $1, "  ZJ1 "\n\t"
		"vmad $2, $6, $18, $18  \n\t"
		  "vldd $2, "  ZJ2 "\n\t"
		"vmad $3, $7, $19, $19  \n\t"
		  "vldd $3, "  ZJ3 "\n\t"

		// DZ,DVZ

		// dz
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VZJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VZJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VZJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VZJ3 "\n\t"
		// dvz
		"vsubd $4,  $26, $4 \n\t"
		"vsubd $5,  $26, $5 \n\t"
		"vsubd $6,  $26, $6 \n\t"
		"vsubd $7,  $26, $7 \n\t"
		// store dz
		// dz * dz
		"vmad $0, $0, $12, $12 \n\t"
		  "vstd $0, "  DZ4 "\n\t"
		"vmad $1, $1, $13, $13 \n\t"
		  "vstd $1, "  DZ5 "\n\t"
		"vmad $2, $2, $14, $14 \n\t"
		  "vstd $2, "  DZ6 "\n\t"
		"vmad $3, $3, $15, $15 \n\t"
		  "vstd $3, "  DZ7 "\n\t"
		// dz * dvz
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $20, " V1 " \n\t" // mask
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $21, " V2 " \n\t" // magic
		"vmad $2, $6, $18, $18  \n\t"
		"vmad $3, $7, $19, $19  \n\t"

		// store rdotv

		// RSQRT-X5
		"srlow $8,  1, $0 \n\t" // P0, 2 cycle
		"srlow $9,  1, $1 \n\t"
		"srlow $10, 1, $2 \n\t"
		"srlow $11, 1, $3 \n\t"
		"srlow $12, 1, $4 \n\t"
		  "vstd $16, "  RV4 "\n\t"
		"srlow $13, 1, $5 \n\t"
		  "vstd $17, "  RV5 "\n\t"
		"srlow $14, 1, $6 \n\t"
		  "vstd $18, "  RV6 "\n\t"
		"srlow $15, 1, $7 \n\t"
		  "vstd $19, "  RV7 "\n\t"
		"vandw $0, $20, $0 \n\t"
		"vandw $1, $20, $1 \n\t"
		"vandw $2, $20, $2 \n\t"
		"vandw $3, $20, $3 \n\t"
		"vandw $4, $20, $4 \n\t"
		"vandw $5, $20, $5 \n\t"
		"vandw $6, $20, $6 \n\t"
		"vandw $7, $20, $7 \n\t"
		"vsubl $21, $0, $0 \n\t"
		"vsubl $21, $1, $1 \n\t"
		"vsubl $21, $2, $2 \n\t"
		"vsubl $21, $3, $3 \n\t"
		"vsubl $21, $4, $4 \n\t"
		"vsubl $21, $5, $5 \n\t"
		"vsubl $21, $6, $6 \n\t"
		"vsubl $21, $7, $7 \n\t"
#ifdef RSQRT_FULL_CONVERGENCE
		// store rsq
		  "vstd $8, "   RSQ0 "\n\t"
		  "vstd $9, "   RSQ1 "\n\t"
		  "vstd $10, "  RSQ2 "\n\t"
		  "vstd $11, "  RSQ3 "\n\t"
		  "vstd $12, "  RSQ4 "\n\t"
		  "vstd $13, "  RSQ5 "\n\t"
		  "vstd $14, "  RSQ6 "\n\t"
		  "vstd $15, "  RSQ7 "\n\t"
#endif
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
		  "vldd $12, "  MJ0 "\n\t"
		"vmuld $5, $13, $5 \n\t"
		  "vldd $13, "  MJ1 "\n\t"
		"vmuld $6, $14, $6 \n\t"
		  "vldd $14, "  MJ2 "\n\t"
		"vmuld $7, $15, $7 \n\t"
		  "vldd $15, "  MJ3 "\n\t"

#ifdef RSQRT_FULL_CONVERGENCE // for strict accuracy check
		// constants
		  "vldd $30, " V3 " \n\t" // c0 = 1.0
		  "vldd $29, " V4 " \n\t" // c1 = 1/2
		  "vldd $28, " V5 " \n\t" // c2 = 3/8
		// load x=rsq
		  "vldd $8, "   RSQ0 "\n\t"
		  "vldd $9, "   RSQ1 "\n\t"
		  "vldd $10, "  RSQ2 "\n\t"
		  "vldd $11, "  RSQ3 "\n\t"
		  "vldd $12, "  RSQ4 "\n\t"
		  "vldd $13, "  RSQ5 "\n\t"
		  "vldd $14, "  RSQ6 "\n\t"
		  "vldd $15, "  RSQ7 "\n\t"
		// y^2
		"vmuld $0, $0, $16 \n\t"
		"vmuld $1, $1, $17 \n\t"
		"vmuld $2, $2, $18 \n\t"
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
		// p = (1/2)h + (3/8)
		"vmad $16, $28, $29, $8   \n\t"
		"vmad $17, $28, $29, $9   \n\t"
		"vmad $18, $28, $29, $10  \n\t"
		"vmad $19, $28, $29, $11  \n\t"
		"vmad $20, $28, $29, $12  \n\t"
		"vmad $21, $28, $29, $13  \n\t"
		"vmad $22, $28, $29, $14  \n\t"
		"vmad $23, $28, $29, $15  \n\t"
		// y * h
		"vmuld $16, $0, $16 \n\t"
		"vmuld $17, $1, $17 \n\t"
		"vmuld $18, $2, $18 \n\t"
		"vmuld $19, $3, $19 \n\t"
		"vmuld $20, $4, $20 \n\t"
		"vmuld $21, $5, $21 \n\t"
		"vmuld $22, $6, $22 \n\t"
		"vmuld $23, $7, $23 \n\t"
		// y += (y*h)*p
		"vmad $16, $8,  $0, $0   \n\t"
		"vmad $17, $9,  $1, $1   \n\t"
		"vmad $18, $10, $2, $2   \n\t"
		"vmad $19, $11, $3, $3   \n\t"
		"vmad $20, $12, $4, $4   \n\t"
		  "vldd $12, "  MJ0 "\n\t"
		"vmad $21, $13, $5, $5   \n\t"
		  "vldd $13, "  MJ1 "\n\t"
		"vmad $22, $14, $6, $6   \n\t"
		  "vldd $14, "  MJ2 "\n\t"
		"vmad $23, $15, $7, $7   \n\t"
		  "vldd $15, "  MJ3 "\n\t"
#endif

		// a0, a1, b0, b1, c, d for $0:$23 (4 for each)
		// load mj
		// load constants

		// b = mj *a;
		"vmuld $12, $0,  $8 \n\t"
		  "vldd $24, "  V3  "\n\t" // 1.0
		"vmuld $13, $1,  $9 \n\t"
		  "vldd $28, "  V12 "\n\t" // r_collv4
		"vmuld $14, $2, $10 \n\t"
		  "vldd $25, "  V9  "\n\t" // m_r_coll_inv_qv4
		"vmuld $15, $3, $11 \n\t"
		  "vldd $26, "  V10 "\n\t" // k_v_miv4
		// c = a * a;
		"vmuld $0, $0, $16 \n\t"
		  "vldd $27, "  V11 "\n\t" // eta_v_miv4
		"vmuld $1, $1, $17 \n\t"
		"vmuld $2, $2, $18 \n\t"
		"vmuld $3, $3, $19 \n\t"
		// d = 1.0 - r_collv4 * a;
		"vnmad $0, $28, $24, $20 \n\t"
		  "vldd $0, " RV0 "\n\t"
		"vnmad $1, $28, $24, $21 \n\t"
		  "vldd $1, " RV1 "\n\t"
		"vnmad $2, $28, $24, $22 \n\t"
		  "vldd $2, " RV2 "\n\t"
		"vnmad $3, $28, $24, $23 \n\t"
		  "vldd $3, " RV3 "\n\t"
		// a = rdotv;
		// b *= c;
		"vmuld $8,  $16, $8  \n\t"
		"vmuld $9,  $17, $9  \n\t"
		"vmuld $10, $18, $10 \n\t"
		"vmuld $11, $19, $11 \n\t"
		// a *= c;
		"vmuld $0,  $16, $0  \n\t"
		"vmuld $1,  $17, $1  \n\t"
		"vmuld $2,  $18, $2  \n\t"
		"vmuld $3,  $19, $3  \n\t"
		// c = m_r_coll_inv_qv4 + k_v_miv4 * d;
		"vmad $20, $26, $25, $16 \n\t"
		"vmad $21, $26, $25, $17 \n\t"
		"vmad $22, $26, $25, $18 \n\t"
		"vmad $23, $26, $25, $19 \n\t"

			// b' = mj *a';
			"vmuld $12, $4, $12 \n\t"
			"vmuld $13, $5, $13 \n\t"
			"vmuld $14, $6, $14 \n\t"
			"vmuld $15, $7, $15 \n\t"

		// c += eta_v_miv4 * a;
		"vmad $27, $0, $16, $16 \n\t"
		"vmad $27, $1, $17, $17 \n\t"
		"vmad $27, $2, $18, $18 \n\t"
		"vmad $27, $3, $19, $19 \n\t"

			// c' = a' * a';
			"vmuld $4, $4, $0 \n\t"
			"vmuld $5, $5, $1 \n\t"
			"vmuld $6, $6, $2 \n\t"
			"vmuld $7, $7, $3 \n\t"

		// b = simd_vsellt(d, c, b);
		"vsellt $20, $16, $8,  $8  \n\t"
		"vsellt $21, $17, $9,  $9  \n\t"
		"vsellt $22, $18, $10, $10 \n\t"
		"vsellt $23, $19, $11, $11 \n\t"

		// DONE: c' for $0:$3 improves scheduling
		// d = 1.0 - r_collv4 * a';
		"vnmad $4, $28, $24, $20 \n\t"
		  "vldd $4, " RV4 "\n\t"
		"vnmad $5, $28, $24, $21 \n\t"
		  "vldd $5, " RV5 "\n\t"
		"vnmad $6, $28, $24, $22 \n\t"
		  "vldd $6, " RV6 "\n\t"
		"vnmad $7, $28, $24, $23 \n\t"
		  "vldd $7, " RV7 "\n\t"
		// a' = rdotv;
		// b' *= c';
		"vmuld $12, $0, $12  \n\t"
		    "ldl  $28, " S2 " \n\t" // force
		"vmuld $13, $1, $13  \n\t"
		"vmuld $14, $2, $14 \n\t"
		"vmuld $15, $3, $15 \n\t"
		// a' *= c';
		"vmuld $4,  $0, $4  \n\t"
		    "vldd $0,   " DX0 " \n\t" 
		"vmuld $5,  $1, $5  \n\t"
		    "vldd $1,    0($28) \n\t" // fx0
		"vmuld $6,  $2, $6  \n\t"
		    "vldd $2,   " DY0 " \n\t" 
		"vmuld $7,  $3, $7  \n\t"
		    "vldd $3,   32($28) \n\t" // fy0
		// c = m_r_coll_inv_qv4 + k_v_miv4 * d;
		"vmad $20, $26, $25, $16 \n\t"
		"vmad $21, $26, $25, $17 \n\t"
		"vmad $22, $26, $25, $18 \n\t"
		"vmad $23, $26, $25, $19 \n\t"

		// BUBBLE
		    "vldd $30, " FSX0 "\n\t"
		    "vldd $29, " FSY0 "\n\t"

		"vmad $8,  $0, $1,  $1 \n\t"
		"vmad $8,  $0, $30, $30 \n\t"
		  "vldd $0,   " DX1 " \n\t" 
		"vmad $8,  $2, $3,  $3 \n\t"
		"vmad $8,  $2, $29, $29 \n\t"
		  "vldd $2,   " DY1 " \n\t" 

		// c += eta_v_miv4 * a';
		"vmad $27, $4, $16, $16 \n\t"
		    "vldd $4,   " DZ0 " \n\t" 
		"vmad $27, $5, $17, $17 \n\t"
		    "vldd $5,   64($28) \n\t" // fz0
		"vmad $27, $6, $18, $18 \n\t"
		    "vldd $6,   " DX4 " \n\t" 
		"vmad $27, $7, $19, $19 \n\t"
		    "vldd $7,  128($28) \n\t" // fx1
		    "vldd $27, " FSZ0 "\n\t"

		// BUBBLE
		"vmad $8,  $4, $5,  $5 \n\t"
		"vmad $8,  $4, $27, $27 \n\t"
		  "vldd $4,   " DZ1 " \n\t" 

		// b' = simd_vsellt(d, c, b');
		"vsellt $20, $16, $12, $12  \n\t"
		    "vldd $20,  " DY4 " \n\t" 
		"vsellt $21, $17, $13, $13  \n\t"
		    "vldd $21, 160($28) \n\t" // fy1
		"vsellt $22, $18, $14, $14 \n\t"
		    "vldd $22,  " DZ4 " \n\t" 
		"vsellt $23, $19, $15, $15 \n\t"
		    "vldd $23, 192($28) \n\t" // fz1


		/*
		 * $0:  dx[0-3]
		 * $1:  fx0
		 * $2:  dy[0-3]
		 * $3:  fy0
		 * $4:  dz[0-3]
		 * $5:  fz0
		 * $6:  dx[4-7]
		 * $7:  fx1
		 * $20: dy[4-7]
		 * $21: fy1
		 * $22: dz[4-7]
		 * $23: fz1
		 */


		    "vldd $26, " FSX1 "\n\t"
		    "vldd $25, " FSY1 "\n\t"
		    "vldd $24, " FSZ1 "\n\t"
		"vmad $12, $6,  $7,  $7 \n\t"
		"vmad $12, $6,  $30, $30 \n\t"
		  "vldd $6,   " DX5 " \n\t" 
		"vmad $12, $20, $21, $21 \n\t"
		"vmad $12, $20, $29, $29 \n\t"
		  "vldd $20,  " DY5 " \n\t" 
		"vmad $12, $22, $23, $23 \n\t"
		"vmad $12, $22, $27, $27 \n\t"
		  "vldd $22,  " DZ5 " \n\t" 

		   "vstd $30, " FSX0 "\n\t"
		   "vstd $29, " FSY0 "\n\t"
		   "vstd $27, " FSZ0 "\n\t"
		   "vldd $30, " FSX2 "\n\t"
		   "vldd $29, " FSY2 "\n\t"
		   "vldd $27, " FSZ2 "\n\t"


		"vmad $9,  $0,  $1,  $1 \n\t"
		"vmad $9,  $0,  $26, $26 \n\t"
		  "vldd $0,   " DX2 " \n\t" 
		"vmad $9,  $2,  $3,  $3 \n\t"
		"vmad $9,  $2,  $25, $25 \n\t"
		  "vldd $2,   " DY2 " \n\t" 
		"vmad $9,  $4,  $5,  $5 \n\t"
		"vmad $9,  $4,  $24, $24 \n\t"
		  "vldd $4,   " DZ2 " \n\t" 
		"vmad $13, $6,  $7,  $7 \n\t"
		"vmad $13, $6,  $26, $26 \n\t"
		  "vldd $6,   " DX6 " \n\t" 
		"vmad $13, $20, $21, $21 \n\t"
		"vmad $13, $20, $25, $25 \n\t"
		  "vldd $20,  " DY6 " \n\t" 
		"vmad $13, $22, $23, $23 \n\t"
		"vmad $13, $22, $24, $24 \n\t"
		  "vldd $22,  " DZ6 " \n\t" 

		    "vstd $26, " FSX1 "\n\t"
		    "vstd $25, " FSY1 "\n\t"
		    "vstd $24, " FSZ1 "\n\t"


		"vmad $10, $0,  $1,  $1 \n\t"
		"vmad $10, $0,  $30, $30 \n\t"
		  "vldd $0,   " DX3 " \n\t" 
		"vmad $10, $2,  $3,  $3 \n\t"
		"vmad $10, $2,  $29, $29 \n\t"
		  "vldd $2,   " DY3 " \n\t" 
		"vmad $10, $4,  $5,  $5 \n\t"
		"vmad $10, $4,  $27, $27 \n\t"
		  "vldd $4,   " DZ3 " \n\t" 
		"vmad $14, $6,  $7,  $7 \n\t"
		"vmad $14, $6,  $30, $30 \n\t"
		  "vldd $6,   " DX7 " \n\t" 
		"vmad $14, $20, $21, $21 \n\t"
		"vmad $14, $20, $29, $29 \n\t"
		  "vldd $20,  " DY7 " \n\t" 
		"vmad $14, $22, $23, $23 \n\t"
		"vmad $14, $22, $27, $27 \n\t"
		  "vldd $22,  " DZ7 " \n\t" 

		    "vldd $26, " FSX3 "\n\t"
		    "vldd $25, " FSY3 "\n\t"
		    "vldd $24, " FSZ3 "\n\t"
		    "vstd $30, " FSX2 "\n\t"
		    "vstd $29, " FSY2 "\n\t"
		    "vstd $27, " FSZ2 "\n\t"

		"vmad $11, $0,  $1,  $1 \n\t"
		  "ldl  $30, " S0 " \n\t" // epi
		"vmad $11, $0,  $26, $26 \n\t"
		  "ldl  $29, " S3 " \n\t" // fend
		"vmad $11, $2,  $3,  $3 \n\t"
		  "ldi $28, 256($28) \n\t" // force+=2
		"vmad $11, $2,  $25, $25 \n\t"
		"vmad $11, $4,  $5,  $5 \n\t"
		"vmad $11, $4,  $24, $24 \n\t"
		  "stl  $28, " S2 " \n\t"  // force
		"vmad $15, $6,  $7,  $7 \n\t"
		"vmad $15, $6,  $26, $26 \n\t"
		"vmad $15, $20, $21, $21 \n\t"
		"vmad $15, $20, $25, $25 \n\t"
		  "ldi $30, 512($30) \n\t" // epi+=2
		"vmad $15, $22, $23, $23 \n\t"
		  "stl  $30, " S0 " \n\t"  // epi
		"vmad $15, $22, $24, $24 \n\t"

		    "vstd $26, " FSX3 "\n\t"
		    "vstd $25, " FSY3 "\n\t"
		    "vstd $24, " FSZ3 "\n\t"

		  "vldd $25,   0($30) \n\t" // xi0
		  "vldd $26,  32($30) \n\t" // vxi0

		// store forces

		"vstd $1,    0-256($28) \n\t"
		"vstd $3,   32-256($28) \n\t"
		"vstd $5,   64-256($28) \n\t"
		"vstd $7,  128-256($28) \n\t"
		"vstd $21, 160-256($28) \n\t"
		"vstd $23, 192-256($28) \n\t"

		"cmpule $29, $28, $0 \n\t"
		"beq $0,.Lt_fdps_6 \n\n\t"
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

#else
// original
void CalcGrav_ps(
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
		"addw $16, 7,   $16  \n\t"
		"srl  $16, 3,   $16  \n\t"
		"sll  $16, 8,   $16  \n\t"

		"addl $16, $19, $16  \n\t" // fend
		  "ldi $30,  85($31) \n\t"
		"stl  $17, " S0 " \n\t" // epi
		  "ldi $29, 170($31) \n\t"
		// "stl  $18, " S1 " \n\t" // epj
		"stl  $19, " S2 " \n\t" // force
		  "ldi $28, 255($31) \n\t"
		"stl  $16, " S3 " \n\t" // fend
		// load epj

		"vldd $20,   0($18) \n\t"
		"vldd $21,  32($18) \n\t"
		"vldd $22,  64($18) \n\t"
		"vldd $23,  96($18) \n\t"
		"vldd $24, 128($18) \n\t"
		"vldd $25, 160($18) \n\t"
		"vldd $26, 192($18) \n\t"
		"vldd $27, 224($18) \n\t"

		// shuffle epj
		"vshff $20, $20, $31, $0  \n\t"
		"vshff $20, $20, $30, $1  \n\t"
		  "vstd $0, "  XJ0 "\n\t"
		"vshff $20, $20, $29, $2  \n\t"
		  "vstd $1, "  YJ0 "\n\t"
		"vshff $20, $20, $28, $3  \n\t"
		  "vstd $2, "  ZJ0 "\n\t"
		"vshff $21, $21, $31, $4  \n\t"
		  "vstd $3, "  MJ0 "\n\t"
		"vshff $21, $21, $30, $5  \n\t"
		  "vstd $4, "  VXJ0 "\n\t"
		"vshff $21, $21, $29, $6  \n\t"
		  "vstd $5, "  VYJ0 "\n\t"
		// "vshff $21, $21, $28, $7  \n\t"

		"vshff $22, $22, $31, $8  \n\t"
		  "vstd $6, "  VZJ0 "\n\t"
		"vshff $22, $22, $30, $9  \n\t"
		  "vstd $8, "  XJ1 "\n\t"
		"vshff $22, $22, $29, $10 \n\t"
		  "vstd $9, "  YJ1 "\n\t"
		"vshff $22, $22, $28, $11 \n\t"
		  "vstd $10, " ZJ1 "\n\t"
		"vshff $23, $23, $31, $12 \n\t"
		  "vstd $11, " MJ1 "\n\t"
		"vshff $23, $23, $30, $13 \n\t"
		  "vstd $12, " VXJ1 "\n\t"
		"vshff $23, $23, $29, $14 \n\t"
		  "vstd $13, " VYJ1 "\n\t"
		// "vshff $23, $23, $28, $15 \n\t"


		"vshff $24, $24, $31, $0  \n\t"
		  "vstd $14, " VZJ1 "\n\t"
		"vshff $24, $24, $30, $1  \n\t"
		  "vstd $0, "  XJ2 "\n\t"
		"vshff $24, $24, $29, $2  \n\t"
		  "vstd $1, "  YJ2 "\n\t"
		"vshff $24, $24, $28, $3  \n\t"
		  "vstd $2, "  ZJ2 "\n\t"
		"vshff $25, $25, $31, $4  \n\t"
		  "vstd $3, "  MJ2 "\n\t"
		"vshff $25, $25, $30, $5  \n\t"
		  "vstd $4, "  VXJ2 "\n\t"
		"vshff $25, $25, $29, $6  \n\t"
		  "vstd $5, "  VYJ2 "\n\t"
		// "vshff $25, $25, $28, $7  \n\t"

		"vshff $26, $26, $31, $8  \n\t"
		  "vstd $6, "  VZJ2 "\n\t"
		"vshff $26, $26, $30, $9  \n\t"
		  "vstd $8, "  XJ3 "\n\t"
		"vshff $26, $26, $29, $10 \n\t"
		  "vstd $9, "  YJ3 "\n\t"
		"vshff $26, $26, $28, $11 \n\t"
		  "vstd $10, " ZJ3 "\n\t"
		"vshff $27, $27, $31, $12 \n\t"
		  "vstd $11, " MJ3 "\n\t"
		"vshff $27, $27, $30, $13 \n\t"
		  "vstd $12, " VXJ3 "\n\t"
		"vshff $27, $27, $29, $14 \n\t"
		  "vstd $13, " VYJ3 "\n\t"
		  "vstd $14, " VZJ3 "\n\t"
		// "vshff $27, $27, $28, $15 \n\t"

		"ldi $30, 0($17)    \n\t" // epi
		"vldd $25,   0($30) \n\t" // xi0
		"vldd $26,  32($30) \n\t" // vxi0
	);
	// loop body
	asm volatile(
		".align 6 \n"
		".Lt_fdps_6: \n\t"

		"ldl  $29, " S1 " \n\t" // epj
		"vldd $28, " V8 "\n\t"   // tiny

		/*
		 *  $20:$27 are free
		 */

		// DX,DVX

		"vldd $0, "  XJ0 "\n\t"
		"vldd $1, "  XJ1 "\n\t"
		"vldd $2, "  XJ2 "\n\t"
		"vldd $3, "  XJ3 "\n\t"

		// dx
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VXJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VXJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VXJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VXJ3 "\n\t"
		// dvx
		"vsubd $4,  $26, $4 \n\t"
		"vsubd $5,  $26, $5 \n\t"
		"vsubd $6,  $26, $6 \n\t"
		  "vldd $25,  64($30) \n\t" // yi0
		"vsubd $7,  $26, $7 \n\t"
		  "vldd $26,  96($30) \n\t" // vyi0
		// store dx
		// dx * dx
		"vmad $0, $0, $28, $8  \n\t"
		  "vstd $0, "  DX0 "\n\t"
		"vmad $1, $1, $28, $9  \n\t"
		  "vstd $1, "  DX1 "\n\t"
		"vmad $2, $2, $28, $10  \n\t"
		  "vstd $2, "  DX2 "\n\t"
		"vmad $3, $3, $28, $11  \n\t"
		  "vstd $3, "  DX3 "\n\t"
		// dx * dvx
		"vmuld $0, $4, $16  \n\t"
		  "vldd $0, "  YJ0 "\n\t"
		"vmuld $1, $5, $17  \n\t"
		  "vldd $1, "  YJ1 "\n\t"
		"vmuld $2, $6, $18  \n\t"
		  "vldd $2, "  YJ2 "\n\t"
		"vmuld $3, $7, $19  \n\t"
		  "vldd $3, "  YJ3 "\n\t"

		// DY,DVY

		// dy
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VYJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VYJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VYJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VYJ3 "\n\t"
		// dvy
		"vsubd $4,  $26, $4 \n\t"
		"vsubd $5,  $26, $5 \n\t"
		"vsubd $6,  $26, $6 \n\t"
		  "vldd $25, 128($30) \n\t" // zi0
		"vsubd $7,  $26, $7 \n\t"
		  "vldd $26, 160($30) \n\t" // vzi0
		// store dy
		// dy * dy
		"vmad $0, $0, $8,  $8  \n\t"
		  "vstd $0, "  DY0 "\n\t"
		"vmad $1, $1, $9,  $9  \n\t"
		  "vstd $1, "  DY1 "\n\t"
		"vmad $2, $2, $10, $10  \n\t"
		  "vstd $2, "  DY2 "\n\t"
		"vmad $3, $3, $11, $11  \n\t"
		  "vstd $3, "  DY3 "\n\t"
		// dy * dvy
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $0, "  ZJ0 "\n\t"
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $1, "  ZJ1 "\n\t"
		"vmad $2, $6, $18, $18  \n\t"
		  "vldd $2, "  ZJ2 "\n\t"
		"vmad $3, $7, $19, $19  \n\t"
		  "vldd $3, "  ZJ3 "\n\t"

		// DZ,DVZ

		// dz
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VZJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VZJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VZJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VZJ3 "\n\t"
		// dvz
		"vsubd $4,  $26, $4 \n\t"
		  "vldd $27, 320($30) \n\t" // yi0
		"vsubd $5,  $26, $5 \n\t"
		  "vldd $29, 352($30) \n\t" // yxi0
		"vsubd $6,  $26, $6 \n\t"
		  "vldd $25, 256($30) \n\t" // xi1
		"vsubd $7,  $26, $7 \n\t"
		  "vldd $26, 288($30) \n\t" // vxi1
		// store dz
		// dz * dz
		"vmad $0, $0, $8,  $8  \n\t"
		  "vstd $0, "  DZ0 "\n\t"
		"vmad $1, $1, $9,  $9  \n\t"
		  "vstd $1, "  DZ1 "\n\t"
		"vmad $2, $2, $10, $10  \n\t"
		  "vstd $2, "  DZ2 "\n\t"
		"vmad $3, $3, $11, $11  \n\t"
		  "vstd $3, "  DZ3 "\n\t"
		// dz * dvz
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $0, "  XJ0 "\n\t"
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $1, "  XJ1 "\n\t"
		"vmad $2, $6, $18, $18  \n\t"
		  "vldd $2, "  XJ2 "\n\t"
		"vmad $3, $7, $19, $19  \n\t"
		  "vldd $3, "  XJ3 "\n\t"
		// store rdotv

		// DX,DVX

		// dx
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VXJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VXJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VXJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VXJ3 "\n\t"
		// dvx
		"vsubd $4,  $26, $4 \n\t"
		  "vstd $16, "  RV0 "\n\t"
		"vsubd $5,  $26, $5 \n\t"
		  "vstd $17, "  RV1 "\n\t"
		"vsubd $6,  $26, $6 \n\t"
		  "vstd $18, "  RV2 "\n\t"
		"vsubd $7,  $26, $7 \n\t"
		  "vstd $19, "  RV3 "\n\t"
		// store dx
		// dx * dx
		"vmad $0, $0, $28, $12  \n\t"
		  "vstd $0, "  DX4 "\n\t"
		"vmad $1, $1, $28, $13  \n\t"
		  "vstd $1, "  DX5 "\n\t"
		"vmad $2, $2, $28, $14  \n\t"
		  "vstd $2, "  DX6 "\n\t"
		"vmad $3, $3, $28, $15  \n\t"
		  "vstd $3, "  DX7 "\n\t"
		// dx * dvx
		"vmuld $0, $4, $16  \n\t"
		  "vldd $0, "  YJ0 "\n\t"
		"vmuld $1, $5, $17  \n\t"
		  "vldd $1, "  YJ1 "\n\t"
		"vmuld $2, $6, $18  \n\t"
		  "vldd $2, "  YJ2 "\n\t"
		"vmuld $3, $7, $19  \n\t"
		  "vldd $3, "  YJ3 "\n\t"

		// DY,DVY

		// dy
		"vsubd $0,  $27, $0 \n\t"
		  "vldd $4, "  VYJ0 "\n\t"
		"vsubd $1,  $27, $1 \n\t"
		  "vldd $5, "  VYJ1 "\n\t"
		"vsubd $2,  $27, $2 \n\t"
		  "vldd $6, "  VYJ2 "\n\t"
		"vsubd $3,  $27, $3 \n\t"
		  "vldd $7, "  VYJ3 "\n\t"
		// dvy
		"vsubd $4,  $29, $4 \n\t"
		"vsubd $5,  $29, $5 \n\t"
		"vsubd $6,  $29, $6 \n\t"
		  "vldd $25, 384($30) \n\t" // zi0
		"vsubd $7,  $29, $7 \n\t"
		  "vldd $26, 416($30) \n\t" // zxi0
		// store dy
		// dy * dy
		"vmad $0, $0, $12, $12 \n\t"
		  "vstd $0, "  DY4 "\n\t"
		"vmad $1, $1, $13, $13 \n\t"
		  "vstd $1, "  DY5 "\n\t"
		"vmad $2, $2, $14, $14 \n\t"
		  "vstd $2, "  DY6 "\n\t"
		"vmad $3, $3, $15, $15 \n\t"
		  "vstd $3, "  DY7 "\n\t"
		// dy * dvy
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $0, "  ZJ0 "\n\t"
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $1, "  ZJ1 "\n\t"
		"vmad $2, $6, $18, $18  \n\t"
		  "vldd $2, "  ZJ2 "\n\t"
		"vmad $3, $7, $19, $19  \n\t"
		  "vldd $3, "  ZJ3 "\n\t"

		// DZ,DVZ

		// dz
		"vsubd $0,  $25, $0 \n\t"
		  "vldd $4, "  VZJ0 "\n\t"
		"vsubd $1,  $25, $1 \n\t"
		  "vldd $5, "  VZJ1 "\n\t"
		"vsubd $2,  $25, $2 \n\t"
		  "vldd $6, "  VZJ2 "\n\t"
		"vsubd $3,  $25, $3 \n\t"
		  "vldd $7, "  VZJ3 "\n\t"
		// dvz
		"vsubd $4,  $26, $4 \n\t"
		"vsubd $5,  $26, $5 \n\t"
		"vsubd $6,  $26, $6 \n\t"
		"vsubd $7,  $26, $7 \n\t"
		// store dz
		// dz * dz
		"vmad $0, $0, $12, $12 \n\t"
		  "vstd $0, "  DZ4 "\n\t"
		"vmad $1, $1, $13, $13 \n\t"
		  "vstd $1, "  DZ5 "\n\t"
		"vmad $2, $2, $14, $14 \n\t"
		  "vstd $2, "  DZ6 "\n\t"
		"vmad $3, $3, $15, $15 \n\t"
		  "vstd $3, "  DZ7 "\n\t"
		// dz * dvz
		"vmad $0, $4, $16, $16  \n\t"
		  "vldd $20, " V1 " \n\t" // mask
		"vmad $1, $5, $17, $17  \n\t"
		  "vldd $21, " V2 " \n\t" // magic
		"vmad $2, $6, $18, $18  \n\t"
		"vmad $3, $7, $19, $19  \n\t"

		// store rdotv

		// RSQRT-X5
		"srlow $8,  1, $0 \n\t" // P0, 2 cycle
		"srlow $9,  1, $1 \n\t"
		"srlow $10, 1, $2 \n\t"
		"srlow $11, 1, $3 \n\t"
		"srlow $12, 1, $4 \n\t"
		  "vstd $16, "  RV4 "\n\t"
		"srlow $13, 1, $5 \n\t"
		  "vstd $17, "  RV5 "\n\t"
		"srlow $14, 1, $6 \n\t"
		  "vstd $18, "  RV6 "\n\t"
		"srlow $15, 1, $7 \n\t"
		  "vstd $19, "  RV7 "\n\t"
		"vandw $0, $20, $0 \n\t"
		"vandw $1, $20, $1 \n\t"
		"vandw $2, $20, $2 \n\t"
		"vandw $3, $20, $3 \n\t"
		"vandw $4, $20, $4 \n\t"
		"vandw $5, $20, $5 \n\t"
		"vandw $6, $20, $6 \n\t"
		"vandw $7, $20, $7 \n\t"
		"vsubl $21, $0, $0 \n\t"
		"vsubl $21, $1, $1 \n\t"
		"vsubl $21, $2, $2 \n\t"
		"vsubl $21, $3, $3 \n\t"
		"vsubl $21, $4, $4 \n\t"
		"vsubl $21, $5, $5 \n\t"
		"vsubl $21, $6, $6 \n\t"
		"vsubl $21, $7, $7 \n\t"
#ifdef RSQRT_FULL_CONVERGENCE
		// store rsq
		  "vstd $8, "   RSQ0 "\n\t"
		  "vstd $9, "   RSQ1 "\n\t"
		  "vstd $10, "  RSQ2 "\n\t"
		  "vstd $11, "  RSQ3 "\n\t"
		  "vstd $12, "  RSQ4 "\n\t"
		  "vstd $13, "  RSQ5 "\n\t"
		  "vstd $14, "  RSQ6 "\n\t"
		  "vstd $15, "  RSQ7 "\n\t"
#endif
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
		  "vldd $12, "  MJ0 "\n\t"
		"vmuld $5, $13, $5 \n\t"
		  "vldd $13, "  MJ1 "\n\t"
		"vmuld $6, $14, $6 \n\t"
		  "vldd $14, "  MJ2 "\n\t"
		"vmuld $7, $15, $7 \n\t"
		  "vldd $15, "  MJ3 "\n\t"

#ifdef RSQRT_FULL_CONVERGENCE // for strict accuracy check
		// constants
		  "vldd $30, " V3 " \n\t" // c0 = 1.0
		  "vldd $29, " V4 " \n\t" // c1 = 1/2
		  "vldd $28, " V5 " \n\t" // c2 = 3/8
		// load x=rsq
		  "vldd $8, "   RSQ0 "\n\t"
		  "vldd $9, "   RSQ1 "\n\t"
		  "vldd $10, "  RSQ2 "\n\t"
		  "vldd $11, "  RSQ3 "\n\t"
		  "vldd $12, "  RSQ4 "\n\t"
		  "vldd $13, "  RSQ5 "\n\t"
		  "vldd $14, "  RSQ6 "\n\t"
		  "vldd $15, "  RSQ7 "\n\t"
		// y^2
		"vmuld $0, $0, $16 \n\t"
		"vmuld $1, $1, $17 \n\t"
		"vmuld $2, $2, $18 \n\t"
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
		// p = (1/2)h + (3/8)
		"vmad $16, $28, $29, $8   \n\t"
		"vmad $17, $28, $29, $9   \n\t"
		"vmad $18, $28, $29, $10  \n\t"
		"vmad $19, $28, $29, $11  \n\t"
		"vmad $20, $28, $29, $12  \n\t"
		"vmad $21, $28, $29, $13  \n\t"
		"vmad $22, $28, $29, $14  \n\t"
		"vmad $23, $28, $29, $15  \n\t"
		// y * h
		"vmuld $16, $0, $16 \n\t"
		"vmuld $17, $1, $17 \n\t"
		"vmuld $18, $2, $18 \n\t"
		"vmuld $19, $3, $19 \n\t"
		"vmuld $20, $4, $20 \n\t"
		"vmuld $21, $5, $21 \n\t"
		"vmuld $22, $6, $22 \n\t"
		"vmuld $23, $7, $23 \n\t"
		// y += (y*h)*p
		"vmad $16, $8,  $0, $0   \n\t"
		"vmad $17, $9,  $1, $1   \n\t"
		"vmad $18, $10, $2, $2   \n\t"
		"vmad $19, $11, $3, $3   \n\t"
		"vmad $20, $12, $4, $4   \n\t"
		  "vldd $12, "  MJ0 "\n\t"
		"vmad $21, $13, $5, $5   \n\t"
		  "vldd $13, "  MJ1 "\n\t"
		"vmad $22, $14, $6, $6   \n\t"
		  "vldd $14, "  MJ2 "\n\t"
		"vmad $23, $15, $7, $7   \n\t"
		  "vldd $15, "  MJ3 "\n\t"
#endif

		// a0, a1, b0, b1, c, d for $0:$23 (4 for each)
		// load mj
		// load constants

		// b = mj *a;
		"vmuld $12, $0,  $8 \n\t"
		  "vldd $24, "  V3  "\n\t" // 1.0
		"vmuld $13, $1,  $9 \n\t"
		  "vldd $28, "  V12 "\n\t" // r_collv4
		"vmuld $14, $2, $10 \n\t"
		  "vldd $25, "  V9  "\n\t" // m_r_coll_inv_qv4
		"vmuld $15, $3, $11 \n\t"
		  "vldd $26, "  V10 "\n\t" // k_v_miv4
		// c = a * a;
		"vmuld $0, $0, $16 \n\t"
		  "vldd $27, "  V11 "\n\t" // eta_v_miv4
		"vmuld $1, $1, $17 \n\t"
		"vmuld $2, $2, $18 \n\t"
		"vmuld $3, $3, $19 \n\t"
		// d = 1.0 - r_collv4 * a;
		"vnmad $0, $28, $24, $20 \n\t"
		  "vldd $0, " RV0 "\n\t"
		"vnmad $1, $28, $24, $21 \n\t"
		  "vldd $1, " RV1 "\n\t"
		"vnmad $2, $28, $24, $22 \n\t"
		  "vldd $2, " RV2 "\n\t"
		"vnmad $3, $28, $24, $23 \n\t"
		  "vldd $3, " RV3 "\n\t"
		// a = rdotv;
		// b *= c;
		"vmuld $8,  $16, $8  \n\t"
		"vmuld $9,  $17, $9  \n\t"
		"vmuld $10, $18, $10 \n\t"
		"vmuld $11, $19, $11 \n\t"
		// a *= c;
		"vmuld $0,  $16, $0  \n\t"
		"vmuld $1,  $17, $1  \n\t"
		"vmuld $2,  $18, $2  \n\t"
		"vmuld $3,  $19, $3  \n\t"
		// c = m_r_coll_inv_qv4 + k_v_miv4 * d;
		"vmad $20, $26, $25, $16 \n\t"
		"vmad $21, $26, $25, $17 \n\t"
		"vmad $22, $26, $25, $18 \n\t"
		"vmad $23, $26, $25, $19 \n\t"

			// b' = mj *a';
			"vmuld $12, $4, $12 \n\t"
			"vmuld $13, $5, $13 \n\t"
			"vmuld $14, $6, $14 \n\t"
			"vmuld $15, $7, $15 \n\t"

		// c += eta_v_miv4 * a;
		"vmad $27, $0, $16, $16 \n\t"
		"vmad $27, $1, $17, $17 \n\t"
		"vmad $27, $2, $18, $18 \n\t"
		"vmad $27, $3, $19, $19 \n\t"

			// c' = a' * a';
			"vmuld $4, $4, $0 \n\t"
			"vmuld $5, $5, $1 \n\t"
			"vmuld $6, $6, $2 \n\t"
			"vmuld $7, $7, $3 \n\t"

		// b = simd_vsellt(d, c, b);
		"vsellt $20, $16, $8,  $8  \n\t"
		"vsellt $21, $17, $9,  $9  \n\t"
		"vsellt $22, $18, $10, $10 \n\t"
		"vsellt $23, $19, $11, $11 \n\t"

		// DONE: c' for $0:$3 improves scheduling
		// d = 1.0 - r_collv4 * a';
		"vnmad $4, $28, $24, $20 \n\t"
		  "vldd $4, " RV4 "\n\t"
		"vnmad $5, $28, $24, $21 \n\t"
		  "vldd $5, " RV5 "\n\t"
		"vnmad $6, $28, $24, $22 \n\t"
		  "vldd $6, " RV6 "\n\t"
		"vnmad $7, $28, $24, $23 \n\t"
		  "vldd $7, " RV7 "\n\t"
		// a' = rdotv;
		// b' *= c';
		"vmuld $12, $0, $12  \n\t"
		    "ldl  $28, " S2 " \n\t" // force
		"vmuld $13, $1, $13  \n\t"
		"vmuld $14, $2, $14 \n\t"
		"vmuld $15, $3, $15 \n\t"
		// a' *= c';
		"vmuld $4,  $0, $4  \n\t"
		    "vldd $0,   " DX0 " \n\t" 
		"vmuld $5,  $1, $5  \n\t"
		    "vldd $1,    0($28) \n\t" // fx0
		"vmuld $6,  $2, $6  \n\t"
		    "vldd $2,   " DY0 " \n\t" 
		"vmuld $7,  $3, $7  \n\t"
		    "vldd $3,   32($28) \n\t" // fy0
		// c = m_r_coll_inv_qv4 + k_v_miv4 * d;
		"vmad $20, $26, $25, $16 \n\t"
		    "vldd $30, " FSX0 "\n\t"
		"vmad $21, $26, $25, $17 \n\t"
		    "vldd $29, " FSY0 "\n\t"
		"vmad $22, $26, $25, $18 \n\t"
		"vmad $23, $26, $25, $19 \n\t"

		// BUBBLE

		"vmad $8,  $0, $1,  $1 \n\t"
		"vmad $8,  $0, $30, $30 \n\t"
		  "vldd $0,   " DX1 " \n\t" 
		"vmad $8,  $2, $3,  $3 \n\t"
		"vmad $8,  $2, $29, $29 \n\t"
		  "vldd $2,   " DY1 " \n\t" 

		// c += eta_v_miv4 * a';
		"vmad $27, $4, $16, $16 \n\t"
		    "vldd $4,   " DZ0 " \n\t" 
		"vmad $27, $5, $17, $17 \n\t"
		    "vldd $5,   64($28) \n\t" // fz0
		"vmad $27, $6, $18, $18 \n\t"
		    "vldd $6,   " DX4 " \n\t" 
		"vmad $27, $7, $19, $19 \n\t"
		    "vldd $27, " FSZ0 "\n\t"

		// BUBBLE
		"vmad $8,  $4, $5,  $5 \n\t"
		    "vldd $7,  128($28) \n\t" // fx1
		"vmad $8,  $4, $27, $27 \n\t"
		  "vldd $4,   " DZ1 " \n\t" 

		// b' = simd_vsellt(d, c, b');
		"vsellt $20, $16, $12, $12  \n\t"
		    "vldd $20,  " DY4 " \n\t" 
		"vsellt $21, $17, $13, $13  \n\t"
		    "vldd $21, 160($28) \n\t" // fy1
		"vsellt $22, $18, $14, $14 \n\t"
		    "vldd $22,  " DZ4 " \n\t" 
		"vsellt $23, $19, $15, $15 \n\t"
		    "vldd $23, 192($28) \n\t" // fz1


		/*
		 * $0:  dx[0-3]
		 * $1:  fx0
		 * $2:  dy[0-3]
		 * $3:  fy0
		 * $4:  dz[0-3]
		 * $5:  fz0
		 * $6:  dx[4-7]
		 * $7:  fx1
		 * $20: dy[4-7]
		 * $21: fy1
		 * $22: dz[4-7]
		 * $23: fz1
		 */


		"vmad $12, $6,  $7,  $7 \n\t"
		    "vldd $26, " FSX1 "\n\t"
		"vmad $12, $6,  $30, $30 \n\t"
		  "vldd $6,   " DX5 " \n\t" 
		"vmad $12, $20, $21, $21 \n\t"
		    "vldd $25, " FSY1 "\n\t"
		"vmad $12, $20, $29, $29 \n\t"
		  "vldd $20,  " DY5 " \n\t" 
		"vmad $12, $22, $23, $23 \n\t"
		    "vldd $24, " FSZ1 "\n\t"
		"vmad $12, $22, $27, $27 \n\t"
		  "vldd $22,  " DZ5 " \n\t" 



		"vmad $9,  $0,  $1,  $1 \n\t"
		   "vstd $30, " FSX0 "\n\t"
		"vmad $9,  $0,  $26, $26 \n\t"
		  "vldd $0,   " DX2 " \n\t" 
		"vmad $9,  $2,  $3,  $3 \n\t"
		   "vstd $29, " FSY0 "\n\t"
		"vmad $9,  $2,  $25, $25 \n\t"
		  "vldd $2,   " DY2 " \n\t" 
		"vmad $9,  $4,  $5,  $5 \n\t"
		   "vstd $27, " FSZ0 "\n\t"
		"vmad $9,  $4,  $24, $24 \n\t"
		  "vldd $4,   " DZ2 " \n\t" 
		"vmad $13, $6,  $7,  $7 \n\t"
		   "vldd $30, " FSX2 "\n\t"
		"vmad $13, $6,  $26, $26 \n\t"
		  "vldd $6,   " DX6 " \n\t" 
		"vmad $13, $20, $21, $21 \n\t"
		   "vldd $29, " FSY2 "\n\t"
		"vmad $13, $20, $25, $25 \n\t"
		  "vldd $20,  " DY6 " \n\t" 
		"vmad $13, $22, $23, $23 \n\t"
		   "vldd $27, " FSZ2 "\n\t"
		"vmad $13, $22, $24, $24 \n\t"
		  "vldd $22,  " DZ6 " \n\t" 



		"vmad $10, $0,  $1,  $1 \n\t"
		    "vstd $26, " FSX1 "\n\t"
		"vmad $10, $0,  $30, $30 \n\t"
		  "vldd $0,   " DX3 " \n\t" 
		"vmad $10, $2,  $3,  $3 \n\t"
		    "vstd $25, " FSY1 "\n\t"
		"vmad $10, $2,  $29, $29 \n\t"
		  "vldd $2,   " DY3 " \n\t" 
		"vmad $10, $4,  $5,  $5 \n\t"
		    "vstd $24, " FSZ1 "\n\t"
		"vmad $10, $4,  $27, $27 \n\t"
		  "vldd $4,   " DZ3 " \n\t" 
		"vmad $14, $6,  $7,  $7 \n\t"
		    "vldd $26, " FSX3 "\n\t"
		"vmad $14, $6,  $30, $30 \n\t"
		  "vldd $6,   " DX7 " \n\t" 
		"vmad $14, $20, $21, $21 \n\t"
		    "vldd $25, " FSY3 "\n\t"
		"vmad $14, $20, $29, $29 \n\t"
		  "vldd $20,  " DY7 " \n\t" 
		"vmad $14, $22, $23, $23 \n\t"
		    "vldd $24, " FSZ3 "\n\t"
		"vmad $14, $22, $27, $27 \n\t"
		  "vldd $22,  " DZ7 " \n\t" 


		"vmad $11, $0,  $1,  $1 \n\t"
		  "ldl  $30, " S0 " \n\t" // epi
		"vmad $11, $0,  $26, $26 \n\t"
		  "ldl  $29, " S3 " \n\t" // fend
		"vmad $11, $2,  $3,  $3 \n\t"
		  "ldi $28, 256($28) \n\t" // force+=2
		"vmad $11, $2,  $25, $25 \n\t"
		  "stl  $28, " S2 " \n\t"  // force
		"vmad $11, $4,  $5,  $5 \n\t"
		    "vstd $30, " FSX2 "\n\t"
		"vmad $11, $4,  $24, $24 \n\t"
		    "vstd $29, " FSY2 "\n\t"
		"vmad $15, $6,  $7,  $7 \n\t"
		    "vstd $27, " FSZ2 "\n\t"
		"vmad $15, $6,  $26, $26 \n\t"
		"vmad $15, $20, $21, $21 \n\t"
		"vmad $15, $20, $25, $25 \n\t"
		  "ldi $30, 512($30) \n\t" // epi+=2
		"vmad $15, $22, $23, $23 \n\t"
		  "stl  $30, " S0 " \n\t"  // epi
		"vmad $15, $22, $24, $24 \n\t"

		    "vstd $26, " FSX3 "\n\t"
		    "vstd $25, " FSY3 "\n\t"
		    "vstd $24, " FSZ3 "\n\t"

		  "vldd $25,   0($30) \n\t" // xi0
		  "vldd $26,  32($30) \n\t" // vxi0

		// store forces

		"vstd $1,    0-256($28) \n\t"
		"vstd $3,   32-256($28) \n\t"
		"vstd $5,   64-256($28) \n\t"
		"vstd $7,  128-256($28) \n\t"
		"vstd $21, 160-256($28) \n\t"
		"vstd $23, 192-256($28) \n\t"

		"cmpule $29, $28, $0 \n\t"
		"beq $0,.Lt_fdps_6 \n\n\t"
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
#endif
#endif
