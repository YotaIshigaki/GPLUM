#if 0
#include <slave.h>
#include <simd.h>

typedef struct Work{
	doublev4 dx;    //   0
	doublev4 dy;    //  32
	doublev4 dz;    //  64
	doublev4 rdotv; //  96
	doublev4 r2;
	doublev4 rinv;  // 160
} Work, *pWork;

typedef double REAL;
typedef struct{
  REAL r[4]; //x,y,z,m
  REAL v[4]; //vx,vy,vz,id
} EpiMM;

typedef EpiMM EpjMM;

typedef struct{
  REAL ax[4],ay[4],az[4],pot[4];
} ForceMM4;

typedef struct{
  REAL ax,ay,az,pot;
} ForceMM;

typedef struct{
  doublev4 m_r_coll_inv_qv4;
  doublev4 k_v_miv4;
  doublev4 eta_v_miv4;   
  doublev4 r_collv4;   
} Force_coff;

void force_loop_C(
		const int         n_epi, 
		// const EpiMM      *epi,     // double r[4], v[4]
		const EpjMM      *epj,     // double r[4], v[4]
		ForceMM4         *force,   // double ax[4], ay[4], az[4], pot[4];
		const Force_coff *pcoff,   // doublev4 m_r_coll_inv_qv4, k_v_miv4, eta_v_miv4;   
		Work             *work     // doublev4 dx, dy, dz, rdotv, r2, rinv
		)
{
	int i, j;
	const doublev4 m_r_coll_inv_qv4 = pcoff->m_r_coll_inv_qv4;
	const doublev4 k_v_miv4         = pcoff->k_v_miv4; 
	const doublev4 eta_v_miv4       = pcoff->eta_v_miv4;   
	const doublev4 r_collv4         = pcoff->r_collv4;   

	const double mjs[4] = {epj[0].r[3], epj[1].r[3], epj[2].r[3], epj[3].r[3],};

	for(i=0; i<n_epi; i+=4){
		doublev4 fx = *(doublev4 *)force[i/4].ax;
		doublev4 fy = *(doublev4 *)force[i/4].ay;
		doublev4 fz = *(doublev4 *)force[i/4].az;
		doublev4 fp = *(doublev4 *)force[i/4].pot;

		for(j=0; j<4; j++){
			// const double mj = epj[j].r[3];
			const double mj = mjs[j];

			const doublev4 dx    = work[i+j].dx;
			const doublev4 dy    = work[i+j].dy;
			const doublev4 dz    = work[i+j].dz;
			const doublev4 rdotv = work[i+j].rdotv;
			// const doublev4 r_sq  = work[i+j].r2;
			const doublev4 r_inv = work[i+j].rinv;

			const doublev4 r_inv_sq = r_inv * r_inv;
			const doublev4 xij      = 1.0 - r_collv4 * r_inv;
			const doublev4 rvji     = rdotv * r_inv_sq; 

			const doublev4 fsd = m_r_coll_inv_qv4 + k_v_miv4 * xij + eta_v_miv4 * rvji;
			const doublev4 pij = mj * r_inv;
			const doublev4 mri3 = simd_vsellt(xij, fsd, r_inv_sq * pij);

			fx -= mri3 * dx;
			fy -= mri3 * dy;
			fz -= mri3 * dz;
			fp -= pij;
		}

		*(doublev4 *)force[i/4].ax  = fx;
		*(doublev4 *)force[i/4].ay  = fy;
		*(doublev4 *)force[i/4].az  = fz;
		*(doublev4 *)force[i/4].pot = fp;
	}
}

void force_loop_C2(
		const int         n_epi,   // $16
		const EpjMM      *epj,     // $18
		ForceMM4         *force,   // $19
		const Force_coff *pcoff,   // $20
		Work             *work     // $21
		)
{
	const doublev4 m_r_coll_inv_qv4 = pcoff->m_r_coll_inv_qv4;
	const doublev4 k_v_miv4         = pcoff->k_v_miv4; 
	const doublev4 eta_v_miv4       = pcoff->eta_v_miv4;   
	const doublev4 r_collv4         = pcoff->r_collv4;   

	const doublev4 mj = {epj[0].r[3], epj[1].r[3], epj[2].r[3], epj[3].r[3]};
	// const doublev4 mj = epj[0].r[3];
	// const double mjs[4] = {epj[0].r[3], epj[1].r[3], epj[2].r[3], epj[3].r[3],};

	ForceMM4 *fend = force + ((n_epi+3)/4);

	while((long)force < (long)fend){
	// while(i < n_epi){
		int j;
		doublev4 fx = *(doublev4 *)force[0].ax;
		doublev4 fy = *(doublev4 *)force[0].ay;
		doublev4 fz = *(doublev4 *)force[0].az;
		// doublev4 fp = *(doublev4 *)force[0].pot;

		for(j=0; j<4; j++){
			// const double mj = epj[j].r[3];
			const doublev4 dx    = work[j].dx;
			const doublev4 dy    = work[j].dy;
			const doublev4 dz    = work[j].dz;
#if 0
			const doublev4 rdotv = work[j].rdotv;
			const doublev4 r_inv = work[j].rinv;

			const doublev4 r_inv_sq = r_inv * r_inv;
			const doublev4 xij      = 1.0 - r_collv4 * r_inv;
			const doublev4 rvji     = rdotv * r_inv_sq; 

			const doublev4 fsd = m_r_coll_inv_qv4 + k_v_miv4 * xij + eta_v_miv4 * rvji;
			const doublev4 pij = mj * r_inv;
			const doublev4 mri3 = simd_vsellt(xij, fsd, r_inv_sq * pij);

			fx -= mri3 * dx;
			fy -= mri3 * dy;
			fz -= mri3 * dz;
			// fp -= pij;
#else
			doublev4 a, b, c, d;
			a = work[j].rinv;
			b = mj *a;
			c = a * a;
			d = 1.0 - r_collv4 * a;

			a = work[j].rdotv;
			b *= c;
			a *= c;
			c = m_r_coll_inv_qv4 + k_v_miv4 * d;
			c += eta_v_miv4 * a;
			b = simd_vsellt(d, c, b);

			fx -= b * dx;
			fy -= b * dy;
			fz -= b * dz; // 11 operations
#endif
		}

		*(doublev4 *)force[0].ax  = fx;
		*(doublev4 *)force[0].ay  = fy;
		*(doublev4 *)force[0].az  = fz;
		// *(doublev4 *)force[0].pot = fp;

		force++;
		work+=4;
	}
}
#endif

#if 0
 $0  : a0
 $1  : a1
 $2  : a2
 $3  : a3
 $4  : b0
 $5  : b1
 $6  : b2
 $7  : b3
 $8  : c0
 $9  : c1
 $10 : c2
 $11 : c3
 $12 : d0
 $13 : d1
 $14 : d2
 $15 : d3
 $16 : n_epi -> fend
 $17 : epj -> mj
 $18 : force
 $19 : pcoff -> one
 $20 : work
 $21 : mrcol
 $22 : kv
 $23 : eta
 $24 : rcoll
 $25 : dx
 $26 : dy
 $27 : dz
 $28 : fx 
 $29 : fy
 $30 : fz
#endif

#define A0 "$0"
#define A1 "$1"
#define A2 "$2"
#define A3 "$3"
#define B0 "$4"
#define B1 "$5"
#define B2 "$6"
#define B3 "$7"
#define C0 "$8"
#define C1 "$9"
#define C2 "$10"
#define C3 "$11"
#define D0 "$12"
#define D1 "$13"
#define D2 "$14"
#define D3 "$15"

#define ONE    "$19"
#define MRCOL  "$21"
#define KV     "$22"
#define ETA    "$23"
#define RCOLL  "$24"

#define DX  "$25"
#define DY  "$26"
#define DZ  "$27"
#define FX  "$28"
#define FY  "$29"
#define FZ  "$30"

#define VMULD(lhs, rhs, dst)  "vmuld " lhs ", " rhs ", " dst " \n\t"
#define VMAD(A, B, C, dst)    "vmad  "  A ", " B ", " C ", " dst " \n\t"
#define VNMAD(A, B, C, dst)   "vnmad"   A ", " B ", " C ", " dst " \n\t"
#define VSELLT(A, B, C, dst)  "vsellt " A ", " B ", " C ", " dst " \n\t"
#define VLDD(reg, addr) "vldd " reg ", " addr "\n\t"
#define VSTD(reg, addr) "vstd " reg ", " addr "\n\t"

void force_loop(){
	// save registers
	asm volatile(
		"stl $26,  8($31)\n\t"
		"stl $29, 16($31)\n\t"
		"stl $30, 24($31)\n\t"
		"stl $15, 40($31)\n\t"
		"stl $27, 48($31)\n\t"
	);
	// initialize
	asm volatile(
		// fend = force + ((n_epi+3)/4);
		"addw $16, 3,   $16 \n\t"
		"srl  $16, 2,   $16 \n\t"
		"sll  $16, 7,   $16 \n\t"
		"addl $16, $18, $16 \n\t" 

		// epj (sorry for the pollution)
		"ldl $0,   24($17) \n\t" 
		"ldl $1,   88($17) \n\t" 
		"ldl $2,  152($17) \n\t" 
		"ldl $3,  216($17) \n\t" 
		"vshff $0, $0, $31, $0 \n\t"
		"vshff $1, $1, $31, $1 \n\t"
		"vshff $2, $2, $31, $2 \n\t"
		"vshff $3, $3, $31, $3 \n\t"
		VSTD("$0", "128($31)")
		VSTD("$1", "160($31)")
		VSTD("$2", "192($31)")
		VSTD("$3", "224($31)")

		// pcoff
		VLDD(MRCOL,  " 0($19)" )
		VLDD(KV,     "32($19)" )
		VLDD(ETA,    "64($19)" )
		VLDD(RCOLL,  "96($19)" )
		"ldi   $19, 1023($31)     \n\t"
		"sll   $19, 52, $19       \n\t"
		"vshff $19, $19, $31, $19 \n\t"

		VLDD(A0, "160($20)" )
		VLDD(A1, "352($20)" )
		VLDD(A2, "544($20)" )
		VLDD(A3, "736($20)" )
	);
	// loop body
	asm volatile(
		".align 8 \n"
		".Lt_fdps_3:"

		VMULD(A0, A0, C0) VLDD(D0, "128($31)" )
		VMULD(A1, A1, C1) VLDD(D1, "160($31)" )
		VMULD(A2, A2, C2) VLDD(D2, "192($31)" )
		    VNMAD(B2, DX, FX, FX) VLDD(DX, "-192($20)" )
		    VNMAD(B2, DY, FY, FY) VLDD(DY, "-160($20)" )
		    VNMAD(B2, DZ, FZ, FZ) VLDD(DZ, "-128($20)" )
		VMULD(A3, A3, C3) VLDD(D3, "224($31)" )


		VMULD(A0, D0, B0)
		VMULD(A1, D1, B1)
		VMULD(A2, D2, B2)
		  VNMAD(B3, DX, FX, FX)
		  VNMAD(B3, DY, FY, FY)
		  VNMAD(B3, DZ, FZ, FZ)
		VMULD(A3, D3, B3)

		VNMAD(A0, RCOLL, ONE, D0)
		  VLDD(A0, " 96($20)" )
		VNMAD(A1, RCOLL, ONE, D1)
		  VLDD(A1, "288($20)" )
		VNMAD(A2, RCOLL, ONE, D2)
		  VLDD(A2, "480($20)" )
		VNMAD(A3, RCOLL, ONE, D3)
		  VLDD(A3, "672($20)" )

		VMULD(B0, C0, B0)
		VMULD(B1, C1, B1)
		VMULD(B2, C2, B2)
		VMULD(B3, C3, B3)

		VMULD(A0, C0, A0)
		VMULD(A1, C1, A1)
		VMULD(A2, C2, A2)
		VMULD(A3, C3, A3)

		VMAD(D0, KV, MRCOL, C0)
		  VSTD(FX, "-128($18)" )
		VMAD(D1, KV, MRCOL, C1)
		  VSTD(FY,  "-96($18)" )
		VMAD(D2, KV, MRCOL, C2)
		  VSTD(FZ,  "-64($18)" )
		VMAD(D3, KV, MRCOL, C3)

		// BUBBLE

		VMAD(A0, ETA, C0, C0)
		  VLDD(FX,  "0($18)" )
		VMAD(A1, ETA, C1, C1)
		  VLDD(FY, "32($18)" )
		VMAD(A2, ETA, C2, C2)
		  VLDD(FZ, "64($18)" )
		VMAD(A3, ETA, C3, C3)

		// BUBBLE
		VLDD(A0, "928($20)" )
		  "ldi $18, 128($18) \n\t"
		VLDD(A1, "1120($20)" )
		VLDD(A2, "1312($20)" )

		VSELLT(D0, C0, B0, B0)
		  VLDD(DX, "  0($20)" )
		VSELLT(D1, C1, B1, B1)
		  VLDD(DY, " 32($20)" )
		VSELLT(D2, C2, B2, B2)
		  VLDD(DZ, " 64($20)" )
		VSELLT(D3, C3, B3, B3)
		  VLDD(A3, "1504($20)" )
		  "cmpule $16, $18, $8 \n\t"

		// BUBBLE

		VNMAD(B0, DX, FX, FX)
		  VLDD(DX, "192($20)" )
		VNMAD(B0, DY, FY, FY)
		  VLDD(DY, "224($20)" )
		VNMAD(B0, DZ, FZ, FZ)
		  VLDD(DZ, "256($20)" )

		// BUBBLE

		VNMAD(B1, DX, FX, FX)
		  VLDD(DX, "384($20)" )
		VNMAD(B1, DY, FY, FY)
		  VLDD(DY, "416($20)" )
		VNMAD(B1, DZ, FZ, FZ)
		  VLDD(DZ, "448($20)" )

		"ldi $20, 768($20) \n\t"
		"beq $8,.Lt_fdps_3 \n\n\t"
	);
	// post loop
	asm volatile(
		VNMAD(B2, DX, FX, FX)
		  VLDD(DX, "-192($20)" )
		VNMAD(B2, DY, FY, FY)
		  VLDD(DY, "-160($20)" )
		VNMAD(B2, DZ, FZ, FZ)
		  VLDD(DZ, "-128($20)" )

		VNMAD(B3, DX, FX, FX)
		VNMAD(B3, DY, FY, FY)
		VNMAD(B3, DZ, FZ, FZ)

		// BUBBLE

		VSTD(FX, "-128($18)" )
		VSTD(FY,  "-96($18)" )
		VSTD(FZ,  "-64($18)" )
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
