#if 0
#include <slave.h>
#include <simd.h>

typedef struct Work{
	doublev4 dx;
	doublev4 dy;
	doublev4 dz;
	doublev4 rdotv;
	doublev4 r2;
	doublev4 rinv;
} Work, *pWork;

typedef double REAL;
typedef struct{
  REAL r[4]; //x,y,z,m
  REAL v[4]; //vx,vy,vz,id
} EpiMM;

typedef EpiMM EpjMM;

extern __thread_local doublev4 glb_tiny;

void posvel_loop_C(
		const int    N,
		const EpiMM *epi,
		const EpjMM *epj,
		pWork        work)
{
	int i;
	doublev4 tiny = glb_tiny;

	doublev4 rj0 = *(doublev4 *)epj[0].r;
	doublev4 vj0 = *(doublev4 *)epj[0].v;
	doublev4 rj1 = *(doublev4 *)epj[1].r;
	doublev4 vj1 = *(doublev4 *)epj[1].v;
	doublev4 rj2 = *(doublev4 *)epj[2].r;
	doublev4 vj2 = *(doublev4 *)epj[2].v;
	doublev4 rj3 = *(doublev4 *)epj[3].r;
	doublev4 vj3 = *(doublev4 *)epj[3].v;

	for(i=0; i<N; i+=4){
		doublev4 xi  = *(doublev4 *)epi[i+0].r;
		doublev4 vxi = *(doublev4 *)epi[i+0].v;
		doublev4 yi  = *(doublev4 *)epi[i+1].r;
		doublev4 vyi = *(doublev4 *)epi[i+1].v;
		doublev4 zi  = *(doublev4 *)epi[i+2].r;
		doublev4 vzi = *(doublev4 *)epi[i+2].v;

		doublev4 dx, dy, dz, dvx, dvy, dvz, r2, rdotv;

		dx = xi - simd_vshff(rj0, rj0, 0x00); 
		dy = yi - simd_vshff(rj0, rj0, 0x55); 
		dz = zi - simd_vshff(rj0, rj0, 0xaa); 

		dvx = vxi - simd_vshff(vj0, vj0, 0x00); 
		dvy = vyi - simd_vshff(vj0, vj0, 0x55); 
		dvz = vzi - simd_vshff(vj0, vj0, 0xaa); 

		r2    = tiny + dx*dx + dy*dy + dz*dz;
		rdotv = dx*dvx + dy*dvy + dz*dvz; 

		work[i+0].dx    = dx;    //   0
		work[i+0].dy    = dy;    //  32
		work[i+0].dz    = dz;    //  64
		work[i+0].rdotv = rdotv; //  96
		work[i+0].r2    = r2;    // 128

		dx = xi - simd_vshff(rj1, rj1, 0x00); 
		dy = yi - simd_vshff(rj1, rj1, 0x55); 
		dz = zi - simd_vshff(rj1, rj1, 0xaa); 

		dvx = vxi - simd_vshff(vj1, vj1, 0x00); 
		dvy = vyi - simd_vshff(vj1, vj1, 0x55); 
		dvz = vzi - simd_vshff(vj1, vj1, 0xaa); 

		r2    = tiny + dx*dx + dy*dy + dz*dz;
		rdotv = dx*dvx + dy*dvy + dz*dvz; 

		work[i+1].dx    = dx;     // 192
		work[i+1].dy    = dy;     // 224
		work[i+1].dz    = dz;     // 256
		work[i+1].rdotv = rdotv;  // 288
		work[i+1].r2    = r2;     // 320

		dx = xi - simd_vshff(rj2, rj2, 0x00); 
		dy = yi - simd_vshff(rj2, rj2, 0x55); 
		dz = zi - simd_vshff(rj2, rj2, 0xaa); 

		dvx = vxi - simd_vshff(vj2, vj2, 0x00); 
		dvy = vyi - simd_vshff(vj2, vj2, 0x55); 
		dvz = vzi - simd_vshff(vj2, vj2, 0xaa); 

		r2    = tiny + dx*dx + dy*dy + dz*dz;
		rdotv = dx*dvx + dy*dvy + dz*dvz; 

		work[i+2].dx    = dx;    // 384
		work[i+2].dy    = dy;    // 416
		work[i+2].dz    = dz;    // 448
		work[i+2].rdotv = rdotv; // 480
		work[i+2].r2    = r2;    // 512

		dx = xi - simd_vshff(rj3, rj3, 0x00); 
		dy = yi - simd_vshff(rj3, rj3, 0x55); 
		dz = zi - simd_vshff(rj3, rj3, 0xaa); 

		dvx = vxi - simd_vshff(vj3, vj3, 0x00); 
		dvy = vyi - simd_vshff(vj3, vj3, 0x55); 
		dvz = vzi - simd_vshff(vj3, vj3, 0xaa); 

		r2    = tiny + dx*dx + dy*dy + dz*dz;
		rdotv = dx*dvx + dy*dvy + dz*dvz; 

		work[i+3].dx    = dx;    // 576
		work[i+3].dy    = dy;    // 608
		work[i+3].dz    = dz;    // 640
		work[i+3].rdotv = rdotv; // 672
		work[i+3].r2    = r2;    // 704
	}
}
#endif
#if 0
 $0  : (epi < N)
 $l  : dvy1
 $2  : rj0
 $3  : vj0
 $4  : rj1
 $5  : vj1
 $6  : rj2
 $7  : vj2
 $8  : rj2
 $9  : vj2
 $10 : xi
 $11 : vxi
 $12 : yi
 $13 : vyi
 $14 : zi
 $15 : vzi
 $16 : end = epi + N
 $17 : epi
 $18 : epj -> work
 $19 : dx0, dz0
 $20 : dx1, dz1
 $21 : dvx0, dvy0, dvz0
 $22 : dvx1, dvy1, dvz1
 $23 : dy0
 $24 : dy1
 $25 : dvy0
 $26 : rsq0
 $27 : rsq1
 $28 : rdotv0
 $29 : rdotv1
 $30 : tiny
#endif

#define RJ0 "$2"
#define RJ1 "$4"
#define RJ2 "$6"
#define RJ3 "$8"
#define VJ0 "$3"
#define VJ1 "$5"
#define VJ2 "$7"
#define VJ3 "$9"
#define XI  "$10"
#define YI  "$12"
#define ZI  "$14"
#define VXI "$11"
#define VYI "$13"
#define VZI "$15"

#define DX0  "$19"
#define DY0  "$23"
#define DZ0  "$19"
#define DX1  "$20"
#define DY1  "$24"
#define DZ1  "$20"
#define DVX0 "$21"
#define DVY0 "$25"
#define DVZ0 "$21"
#define DVX1 "$22"
#define DVY1 "$1"
#define DVZ1 "$22"
#define IMM55 DVY1
#define IMMaa DVZ1

#define RSQ0   "$26"
#define RSQ1   "$27"
#define RDOTV0 "$28"
#define RDOTV1 "$29"
#define DTINY  "$30"

#define VSHFF2(src, imm, dst) "vshff " src ", " src  ", " imm ", " dst " \n\t"
#define VSUBD(lhs, rhs, dst)  "vsubd " lhs ", " rhs ", " dst " \n\t"
#define VMULD(lhs, rhs, dst)  "vmuld " lhs ", " rhs ", " dst " \n\t"
#define VMAD(A, B, C, dst)   "vmad " A ", " B ", " C ", " dst " \n\t"
#define VLDD(reg, addr) "vldd " reg ", " addr "\n\t"
#define VSTD(reg, addr) "vstd " reg ", " addr "\n\t"
// #define VLDD(reg, addr) ""
// #define VSTD(reg, addr) ""


void posvel_loop(/*
		const int    N,
		const EpiMM *epi,
		const EpjMM *epj,
		pWork        work*/)
{
	// EpiMM *end = epi+N;
	asm volatile(
		"sll $16, 6, $16 \n\t"
		"addl $17, $16, $16 \n\t"
		);
	// save registers
	asm volatile(
		"stl $26,  8($31)\n\t"
		"stl $29, 16($31)\n\t"
		"stl $30, 24($31)\n\t"
		"stl $15, 40($31)\n\t"
		"stl $27, 48($31)\n\t"
	);

	asm volatile(
		  "ldih	$1,glb_tiny($31)         	!tprelhi \n\t"
		"vldd $2,   0($18) \n\t"
		  "vldd	$1,glb_tiny($1)          	!tprello \n\t"
		"vldd $3,  32($18) \n\t"
		"vldd $4,  64($18) \n\t"
		"vldd $5,  96($18) \n\t"
		"vldd $6, 128($18) \n\t"
		"vldd $7, 160($18) \n\t"
		"vldd $8, 192($18) \n\t"
		  "vfmov  $1, $30 \n\t"
		"vldd $9, 224($18) \n\t"
		"mov $19, $18"
		);

	asm volatile(
		VLDD(XI,  "  0($17)")
		VLDD(VXI, " 32($17)")
		VLDD(YI,  " 64($17)")
		VLDD(VYI, " 96($17)")
		VLDD(ZI,  "128($17)")
		VLDD(VZI, "160($17)")
		".align 8 \n"
		".Lt_fdps_2:"
		);
	//for(i=0; i<N; i+=4){
		asm volatile(
			  VMAD(DZ1, DZ1, RSQ1, RSQ1)  VSTD(DZ1, "-128($18)")
			  VMAD(DZ1, DVZ1, RDOTV1, RDOTV1)

			VSHFF2(RJ0, "$31", DX0)
			VSHFF2(VJ0, "$31", DVX0)
			VSUBD (XI,  DX0,  DX0)
			VSUBD (VXI, DVX0, DVX0)

			VSHFF2(RJ1, "$31", DX1)
			VSHFF2(VJ1, "$31", DVX1)
			"ldi " IMM55 ", 85($31) \n\t"
			VSUBD (XI,  DX1,  DX1)
			  VSTD(RSQ0,   "-256($18)")
			VSUBD (VXI, DVX1, DVX1)
			  VSTD(RDOTV0, "-288($18)")

			VSHFF2(RJ0, IMM55, DY0)
			  VSTD(RSQ1,   " -64($18)")
			VSHFF2(VJ0, IMM55, DVY0)
			  VSTD(RDOTV1, " -96($18)")
			
			VMAD(DX0, DX0, DTINY, RSQ0)  VSTD(DX0, "0($18)")
			VMULD(DX0, DVX0, RDOTV0)
			
			VSUBD (YI,  DY0,  DY0)
			VSUBD (VYI, DVY0, DVY0)

			VSHFF2(RJ1, IMM55, DY1)
			VSHFF2(VJ1, IMM55, DVY1)
			
			VMAD(DX1, DX1, DTINY, RSQ1)  VSTD(DX1, "192($18)")
			VMULD(DX1, DVX1, RDOTV1)

			VSUBD (YI,  DY1,  DY1)
			VSUBD (VYI, DVY1, DVY1)

			VMAD(DY0, DY0, RSQ0, RSQ0)  VSTD(DY0, "32($18)")
			VMAD(DY0, DVY0, RDOTV0, RDOTV0)

			"ldi " IMMaa ", 170($31) \n\t"
			VSHFF2(RJ0, IMMaa, DZ0)
			VSHFF2(VJ0, IMMaa, DVZ0)

			VSUBD (ZI,  DZ0,  DZ0)
			VSUBD (VZI, DVZ0, DVZ0)

			VMAD(DY1, DY1, RSQ1, RSQ1)  VSTD(DY1, "224($18)")
			VMAD(DY1, DVY1, RDOTV1, RDOTV1)

			VSHFF2(RJ1, IMMaa, DZ1)
			VSHFF2(VJ1, IMMaa, DVZ1)

			VSUBD (ZI,  DZ1,  DZ1)
			VSUBD (VZI, DVZ1, DVZ1)

			VMAD(DZ0, DZ0, RSQ0, RSQ0)  VSTD(DZ0, "64($18)")
			VMAD(DZ0, DVZ0, RDOTV0, RDOTV0)

			  VSHFF2(RJ2, "$31", DX0)
			  VSHFF2(VJ2, "$31", DVX0)
			  VSUBD (XI,  DX0,  DX0)
			  VSUBD (VXI, DVX0, DVX0)


			VMAD(DZ1, DZ1, RSQ1, RSQ1)  VSTD(DZ1, "256($18)")
			  "ldi " IMM55 ", 85($31) \n\t"
			VMAD(DZ1, DVZ1, RDOTV1, RDOTV1)
			"ldi  $17, 256($17)  \n\t"


			  VSHFF2(RJ3, "$31", DX1)
			VSTD(RSQ0,   "128($18)")
			  VSHFF2(VJ3, "$31", DVX1)
			VSTD(RDOTV0, " 96($18)")
			  VSUBD (XI,  DX1,  DX1)  VLDD(XI,  "  0($17)")
			  VSUBD (VXI, DVX1, DVX1) VLDD(VXI, " 32($17)")

			  VSHFF2(RJ2, IMM55, DY0)
			  VSHFF2(VJ2, IMM55, DVY0)
			
			  VMAD(DX0, DX0, DTINY, RSQ0)  VSTD(DX0, "384($18)")
			  VMULD(DX0, DVX0, RDOTV0)
			
			  VSUBD (YI,  DY0,  DY0)
			VSTD(RSQ1,   "320($18)")
			  VSUBD (VYI, DVY0, DVY0)
			VSTD(RDOTV1, "288($18)")

			  VSHFF2(RJ3, IMM55, DY1)
			  VSHFF2(VJ3, IMM55, DVY1)
			
			  VMAD(DX1, DX1, DTINY, RSQ1)  VSTD(DX1, "576($18)")
			  VMULD(DX1, DVX1, RDOTV1)

			  VSUBD (YI,  DY1,  DY1)  VLDD(YI,  " 64($17)")
			  VSUBD (VYI, DVY1, DVY1) VLDD(VYI, " 96($17)")

			  VMAD(DY0, DY0, RSQ0, RSQ0)  VSTD(DY0, "416($18)")
			  VMAD(DY0, DVY0, RDOTV0, RDOTV0)

			  "ldi " IMMaa ", 170($31) \n\t"
			  VSHFF2(RJ2, IMMaa, DZ0)
			  VSHFF2(VJ2, IMMaa, DVZ0)

			  VSUBD (ZI,  DZ0,  DZ0)
			  VSUBD (VZI, DVZ0, DVZ0)

			  VMAD(DY1, DY1, RSQ1, RSQ1)  VSTD(DY1, "608($18)")
			  VMAD(DY1, DVY1, RDOTV1, RDOTV1)

			  VSHFF2(RJ3, IMMaa, DZ1)
			"ldi  $18, 768($18)  \n\t"
			  VSHFF2(VJ3, IMMaa, DVZ1)
			"cmpule $16, $17, $0 \n\t"

			  VSUBD (ZI,  DZ1,  DZ1)  VLDD(ZI,  "128($17)")
			  VSUBD (VZI, DVZ1, DVZ1) VLDD(VZI, "160($17)")

			  VMAD(DZ0, DZ0, RSQ0, RSQ0)  VSTD(DZ0, "-320($18)")
			  VMAD(DZ0, DVZ0, RDOTV0, RDOTV0)
			"beq $0,.Lt_fdps_2 \n\n\t"

			  VMAD(DZ1, DZ1, RSQ1, RSQ1)  VSTD(DZ1, "-128($18)")
			  VMAD(DZ1, DVZ1, RDOTV1, RDOTV1)

			  VSTD(RSQ0,   "-256($18)")
			  VSTD(RDOTV0, "-288($18)")
			  VSTD(RSQ1,   " -64($18)")
			  VSTD(RDOTV1, " -96($18)")

				);
	// }

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
