#include <cassert>
#include "v2r8.h"

class PhantomGrapeQuad{
public:
	enum{
        NIMAX = 2048,
        NJMAX = 262144,
	};
private:
	double xibuf  [NIMAX/2]  [3][2]; // x, y, z
	double accpbuf[NIMAX/2]  [4][2]; // ax, ay, az, pot
	double epjbuf [NJMAX]    [4];    // x, y, z, m
	double spjbuf [NJMAX]    [10];   // x, y, z, m, xx, yy, xy, yz, zx, tr
	double eps2;
	static double get_a_NaN(){
		union{ long   l; double d; } m;
		m.l = -1;
		return m.d;
	}
public:
	PhantomGrapeQuad() : eps2(get_a_NaN()) {} // default NaN
	void set_eps2(const double _eps2){ this->eps2 = _eps2; }
	void set_epj_one(const int addr, const double x, const double y, const double z, const double m) {
		epjbuf[addr][0] = x;
		epjbuf[addr][1] = y;
		epjbuf[addr][2] = z;
		epjbuf[addr][3] = m;
	}
	template <typename EPJ_t>
	void set_epj(const int nj, const EPJ_t epj[]);
	void set_spj_one(
			const int addr, 
			const double x,   const double y,   const double z,   const double m,
			const double qxx, const double qyy, const double qzz,
			const double qxy, const double qyz, const double qzx)
	{
		const double tr = qxx + qyy + qzz;
		spjbuf[addr][0] = x;
		spjbuf[addr][1] = y;
		spjbuf[addr][2] = z;
		spjbuf[addr][3] = m;
		spjbuf[addr][4] = 3.0 * qxx - tr;
		spjbuf[addr][5] = 3.0 * qyy - tr;
		spjbuf[addr][6] = 3.0 * qxy;
		spjbuf[addr][7] = 3.0 * qyz;
		spjbuf[addr][8] = 3.0 * qzx;
		spjbuf[addr][9] = -(eps2 * tr);
	}
	template <typename SPJ_t>
	void set_spj(const int nj, const SPJ_t spj[]);
	void set_xi_one(const int addr, const double x, const double y, const double z){
		const int ah = addr / 2;
		const int al = addr % 2;
		xibuf[ah][0][al] = x;
		xibuf[ah][1][al] = y;
		xibuf[ah][2][al] = z;
	}
	template <typename EPI_t>
	void set_xi(const int ni, const EPI_t spj[]);
	template <typename real_t>
	void get_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot){
		const int ah = addr / 2;
		const int al = addr % 2;
		ax  = accpbuf[ah][0][al];
		ay  = accpbuf[ah][1][al];
		az  = accpbuf[ah][2][al];
		pot = accpbuf[ah][3][al];
	}
	template <typename real_t>
	void accum_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot){
		const int ah = addr / 2;
		const int al = addr % 2;
		ax  += accpbuf[ah][0][al];
		ay  += accpbuf[ah][1][al];
		az  += accpbuf[ah][2][al];
		pot += accpbuf[ah][3][al];
	}

	void run_epj(const int ni, const int nj){
	    if(ni > NIMAX || nj > NJMAX){
		std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
	    }
	    assert(ni <= NIMAX);
	    assert(nj <= NJMAX);
	    kernel_epj_unroll4(ni, nj);
	}
	void run_spj(const int ni, const int nj){
	    if(ni > NIMAX || nj > NJMAX){
		std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
	    }
	    assert(ni <= NIMAX);
	    assert(nj <= NJMAX);
	    kernel_spj_unroll4(ni, nj);
	}
private:
	__attribute__ ((always_inline))
	void inline_kernel_epj(
			const v2r8 eps2,
			const v2r8 jp_xy,
			const v2r8 jp_zm,
			const v2r8 xi,
			const v2r8 yi,
			const v2r8 zi,
			v2r8 &ax,
			v2r8 &ay,
			v2r8 &az,
			v2r8 &pot){
	    const v2r8_bcl xj(jp_xy);
	    const v2r8_bch yj(jp_xy);
	    const v2r8_bcl zj(jp_zm);
	    const v2r8_bch mj(jp_zm);
	    const v2r8 dx = xj - xi;
	    const v2r8 dy = yj - yi;
	    const v2r8 dz = zj - zi;
	    const v2r8 r2   = ((eps2 + dx*dx) +  dy*dy) + dz*dz;
	    const v2r8 ri1  = r2.rsqrta_x3();
	    const v2r8 ri2  = ri1  * ri1;
	    const v2r8 mri1 = mj   * ri1;
	    const v2r8 mri3 = mri1 * ri2;
	    pot = mj.nmsub(ri1, pot);
	    ax += mri3 * dx;
	    ay += mri3 * dy;
	    az += mri3 * dz;
	}
	__attribute__ ((noinline))
	void kernel_epj_unroll4(const int ni, const int nj){
		for(int i=0; i<ni; i+=8){
			const v2r8 veps2 = v2r8(eps2);
			const v2r8 xi0 = v2r8::load(xibuf[0 + i/2][0]);
			const v2r8 yi0 = v2r8::load(xibuf[0 + i/2][1]);
			const v2r8 zi0 = v2r8::load(xibuf[0 + i/2][2]);
			const v2r8 xi1 = v2r8::load(xibuf[1 + i/2][0]);
			const v2r8 yi1 = v2r8::load(xibuf[1 + i/2][1]);
			const v2r8 zi1 = v2r8::load(xibuf[1 + i/2][2]);
			const v2r8 xi2 = v2r8::load(xibuf[2 + i/2][0]);
			const v2r8 yi2 = v2r8::load(xibuf[2 + i/2][1]);
			const v2r8 zi2 = v2r8::load(xibuf[2 + i/2][2]);
			const v2r8 xi3 = v2r8::load(xibuf[3 + i/2][0]);
			const v2r8 yi3 = v2r8::load(xibuf[3 + i/2][1]);
			const v2r8 zi3 = v2r8::load(xibuf[3 + i/2][2]);
			v2r8 ax0(0.0);
			v2r8 ay0(0.0);
			v2r8 az0(0.0);
			v2r8 po0(0.0);
			v2r8 ax1(0.0);
			v2r8 ay1(0.0);
			v2r8 az1(0.0);
			v2r8 po1(0.0);
			v2r8 ax2(0.0);
			v2r8 ay2(0.0);
			v2r8 az2(0.0);
			v2r8 po2(0.0);
			v2r8 ax3(0.0);
			v2r8 ay3(0.0);
			v2r8 az3(0.0);
			v2r8 po3(0.0);
			for(int j=0; j<nj; j++){
				const v2r8 jp_xy = v2r8::load(epjbuf[j] + 0);
				const v2r8 jp_zm = v2r8::load(epjbuf[j] + 2);
				inline_kernel_epj(veps2, jp_xy, jp_zm, xi0, yi0, zi0, ax0, ay0, az0, po0);
				inline_kernel_epj(veps2, jp_xy, jp_zm, xi1, yi1, zi1, ax1, ay1, az1, po1);
				inline_kernel_epj(veps2, jp_xy, jp_zm, xi2, yi2, zi2, ax2, ay2, az2, po2);
				inline_kernel_epj(veps2, jp_xy, jp_zm, xi3, yi3, zi3, ax3, ay3, az3, po3);
			}
			ax0.store(accpbuf[0 + i/2][0]);
			ay0.store(accpbuf[0 + i/2][1]);
			az0.store(accpbuf[0 + i/2][2]);
			po0.store(accpbuf[0 + i/2][3]);
			ax1.store(accpbuf[1 + i/2][0]);
			ay1.store(accpbuf[1 + i/2][1]);
			az1.store(accpbuf[1 + i/2][2]);
			po1.store(accpbuf[1 + i/2][3]);
			ax2.store(accpbuf[2 + i/2][0]);
			ay2.store(accpbuf[2 + i/2][1]);
			az2.store(accpbuf[2 + i/2][2]);
			po2.store(accpbuf[2 + i/2][3]);
			ax3.store(accpbuf[3 + i/2][0]);
			ay3.store(accpbuf[3 + i/2][1]);
			az3.store(accpbuf[3 + i/2][2]);
			po3.store(accpbuf[3 + i/2][3]);
		}
	}
	__attribute__ ((always_inline))
	void inline_kernel_spj(
			const v2r8 eps2,
			const v2r8 jp_xy,
			const v2r8 jp_zm,
			const v2r8 jp_xx_yy,
			const v2r8 jp_xy_yz,
			const v2r8 jp_zx_tr,
			const v2r8 xi,
			const v2r8 yi,
			const v2r8 zi,
			v2r8 &ax,
			v2r8 &ay,
			v2r8 &az,
			v2r8 &pot){
		const v2r8_bcl xj(jp_xy);
		const v2r8_bch yj(jp_xy);
		const v2r8_bcl zj(jp_zm);
		const v2r8 dx = xj - xi;
		const v2r8 dy = yj - yi;
		const v2r8 dz = zj - zi;
		const v2r8 r2   = ((eps2 + dx*dx) +  dy*dy) + dz*dz;
		const v2r8 ri1  = r2.rsqrta_x3();
		const v2r8 ri2  = ri1  * ri1;
		const v2r8 ri3  = ri1  * ri2;
		const v2r8 ri4  = ri2  * ri2;
		const v2r8 ri5  = ri2  * ri3;
		const v2r8_bcl q_xx(jp_xx_yy);
		const v2r8_bch q_yy(jp_xx_yy);
		const v2r8     q_zz = jp_xx_yy.hadd();  // -q_zz = q_xx + q_yy
		const v2r8_bcl q_xy(jp_xy_yz);
		const v2r8_bch q_yz(jp_xy_yz);
		const v2r8_bcl q_zx(jp_zx_tr);
		const v2r8 q_tr = v2r8::unpckh(jp_zx_tr, jp_zx_tr); // -eps2 * tr(q)
		const v2r8 mj  = v2r8::unpckh(jp_zm, jp_zm);
		v2r8 qr_x = q_xx * dx;
		qr_x = q_xy.madd(dy, qr_x);
		qr_x = q_zx.madd(dz, qr_x);
		v2r8 qr_y = q_xy * dx;
		qr_y = q_yy.madd(dy, qr_y);
		qr_y = q_yz.madd(dz, qr_y);
		v2r8 qr_z = q_zx * dx;
		qr_z = q_yz.madd(dy, qr_z);
		qr_z = qr_z - q_zz * dz;
	   	ax = v2r8::nmsub(ri5, qr_x, ax);
	   	ay = v2r8::nmsub(ri5, qr_y, ay);
	   	az = v2r8::nmsub(ri5, qr_z, az);
		const v2r8 rqr = ((q_tr + qr_x*dx) + qr_y*dy) + qr_z*dz;
		const v2r8 rqr_ri4 = rqr * ri4;
		const v2r8 cc(0.5, 2.5);
		const v2r8 madd05 = v2r8_bcl(cc).madd(rqr_ri4, mj);
		const v2r8 madd25 = v2r8_bch(cc).madd(rqr_ri4, mj) * ri3;
		pot -= madd05 * ri1;
		ax = v2r8::madd(madd25, dx, ax);
		ay = v2r8::madd(madd25, dy, ay);
		az = v2r8::madd(madd25, dz, az);
	}
	__attribute__ ((noinline))
	void kernel_spj_unroll4(const int ni, const int nj){
		for(int i=0; i<ni; i+=8){
			const v2r8 veps2 = v2r8(eps2);
			const v2r8 xi0 = v2r8::load(xibuf[0 + i/2][0]);
			const v2r8 yi0 = v2r8::load(xibuf[0 + i/2][1]);
			const v2r8 zi0 = v2r8::load(xibuf[0 + i/2][2]);
			const v2r8 xi1 = v2r8::load(xibuf[1 + i/2][0]);
			const v2r8 yi1 = v2r8::load(xibuf[1 + i/2][1]);
			const v2r8 zi1 = v2r8::load(xibuf[1 + i/2][2]);
			const v2r8 xi2 = v2r8::load(xibuf[2 + i/2][0]);
			const v2r8 yi2 = v2r8::load(xibuf[2 + i/2][1]);
			const v2r8 zi2 = v2r8::load(xibuf[2 + i/2][2]);
			const v2r8 xi3 = v2r8::load(xibuf[3 + i/2][0]);
			const v2r8 yi3 = v2r8::load(xibuf[3 + i/2][1]);
			const v2r8 zi3 = v2r8::load(xibuf[3 + i/2][2]);
			v2r8 ax0(0.0);
			v2r8 ay0(0.0);
			v2r8 az0(0.0);
			v2r8 po0(0.0);
			v2r8 ax1(0.0);
			v2r8 ay1(0.0);
			v2r8 az1(0.0);
			v2r8 po1(0.0);
			v2r8 ax2(0.0);
			v2r8 ay2(0.0);
			v2r8 az2(0.0);
			v2r8 po2(0.0);
			v2r8 ax3(0.0);
			v2r8 ay3(0.0);
			v2r8 az3(0.0);
			v2r8 po3(0.0);
			for(int j=0; j<nj; j++){
				const v2r8 jp_xy    = v2r8::load(spjbuf[j] + 0);
				const v2r8 jp_zm    = v2r8::load(spjbuf[j] + 2);
				const v2r8 jp_xx_yy = v2r8::load(spjbuf[j] + 4);
				const v2r8 jp_xy_yz = v2r8::load(spjbuf[j] + 6);
				const v2r8 jp_zx_tr = v2r8::load(spjbuf[j] + 8);
				inline_kernel_spj(veps2, jp_xy, jp_zm, jp_xx_yy, jp_xy_yz, jp_zx_tr, xi0, yi0, zi0, ax0, ay0, az0, po0);
				inline_kernel_spj(veps2, jp_xy, jp_zm, jp_xx_yy, jp_xy_yz, jp_zx_tr, xi1, yi1, zi1, ax1, ay1, az1, po1);
				inline_kernel_spj(veps2, jp_xy, jp_zm, jp_xx_yy, jp_xy_yz, jp_zx_tr, xi2, yi2, zi2, ax2, ay2, az2, po2);
				inline_kernel_spj(veps2, jp_xy, jp_zm, jp_xx_yy, jp_xy_yz, jp_zx_tr, xi3, yi3, zi3, ax3, ay3, az3, po3);
			}
			ax0.store(accpbuf[0 + i/2][0]);
			ay0.store(accpbuf[0 + i/2][1]);
			az0.store(accpbuf[0 + i/2][2]);
			po0.store(accpbuf[0 + i/2][3]);
			ax1.store(accpbuf[1 + i/2][0]);
			ay1.store(accpbuf[1 + i/2][1]);
			az1.store(accpbuf[1 + i/2][2]);
			po1.store(accpbuf[1 + i/2][3]);
			ax2.store(accpbuf[2 + i/2][0]);
			ay2.store(accpbuf[2 + i/2][1]);
			az2.store(accpbuf[2 + i/2][2]);
			po2.store(accpbuf[2 + i/2][3]);
			ax3.store(accpbuf[3 + i/2][0]);
			ay3.store(accpbuf[3 + i/2][1]);
			az3.store(accpbuf[3 + i/2][2]);
			po3.store(accpbuf[3 + i/2][3]);
		}
	}
} __attribute__ ((aligned(128)));
