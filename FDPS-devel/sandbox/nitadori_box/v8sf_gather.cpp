#include <cstdio>
#include <cstdlib>
#include <cstring>

typedef float v4sf  __attribute__((vector_size(16)));
typedef float v8sf  __attribute__((vector_size(32)));
typedef int   v8si  __attribute__((vector_size(32)));
typedef int   v4si  __attribute__((vector_size(16)));

void v8si_mask_scatter(v8si val, v8si msk, int *base, v8si idx){
	v4si vlo = __builtin_ia32_vextractf128_si256(val, 0);
	v4si vhi = __builtin_ia32_vextractf128_si256(val, 1);
	v4si mlo = __builtin_ia32_vextractf128_si256(msk, 0);
	v4si mhi = __builtin_ia32_vextractf128_si256(msk, 1);

	v4si take0 = {-1, 0, 0, 0};
	v4si take1 = {0, -1, 0, 0};
	v4si take2 = {0, 0, -1, 0};
	v4si take3 = {0, 0, 0, -1};

	v4si m0 = mlo & take0;
	v4si m1 = mlo & take1;
	v4si m2 = mlo & take2;
	v4si m3 = mlo & take3;
	v4si m4 = mhi & take0;
	v4si m5 = mhi & take1;
	v4si m6 = mhi & take2;
	v4si m7 = mhi & take3;

	int off[8] __attribute__((aligned(32)));
	*(v8si *)off = idx;

	int *adr0 = base + off[0] - 0;
	int *adr1 = base + off[1] - 1;
	int *adr2 = base + off[2] - 2;
	int *adr3 = base + off[3] - 3;
	int *adr4 = base + off[4] - 0;
	int *adr5 = base + off[5] - 1;
	int *adr6 = base + off[6] - 2;
	int *adr7 = base + off[7] - 3;
	// avxintrin.h:  __builtin_ia32_maskstoreps ((__v4sf *)__P, (__v4si)__M, (__v4sf)__A);
	// avx2intrin.h:  __builtin_ia32_maskstored ((__v4si *)__X, (__v4si)__M, (__v4si)__Y);
	__builtin_ia32_maskstored((v4si *)adr0, m0, vlo);
	__builtin_ia32_maskstored((v4si *)adr1, m1, vlo);
	__builtin_ia32_maskstored((v4si *)adr2, m2, vlo);
	__builtin_ia32_maskstored((v4si *)adr3, m3, vlo);

	__builtin_ia32_maskstored((v4si *)adr4, m4, vhi);
	__builtin_ia32_maskstored((v4si *)adr5, m5, vhi);
	__builtin_ia32_maskstored((v4si *)adr6, m6, vhi);
	__builtin_ia32_maskstored((v4si *)adr7, m7, vhi);
}

template<typename PackType, typename MembType>
v8sf v8sf_struct_load(const PackType &, const MembType &base){
	const char *p = (char *)(&base);
	const size_t stride = sizeof(PackType);

	const float f0 = (float)( *(MembType *)(p + 0 * stride) );
	const float f1 = (float)( *(MembType *)(p + 1 * stride) );
	const float f2 = (float)( *(MembType *)(p + 2 * stride) );
	const float f3 = (float)( *(MembType *)(p + 3 * stride) );
	const float f4 = (float)( *(MembType *)(p + 4 * stride) );
	const float f5 = (float)( *(MembType *)(p + 5 * stride) );
	const float f6 = (float)( *(MembType *)(p + 6 * stride) );
	const float f7 = (float)( *(MembType *)(p + 7 * stride) );

	return (v8sf){f0, f1, f2, f3, f4, f5, f6, f7};
}

template<typename PackType, typename MembType>
v8sf v8sf_struct_load_indirect(const PackType &, const MembType &base, const v8si vaddr){
	const char *p = (char *)(&base);
	const size_t stride = sizeof(PackType);

	int addr[8] __attribute__((aligned(32)));
	*(v8si *)addr = vaddr;

	const float f0 = (float)( *(MembType *)(p + addr[0] * stride) );
	const float f1 = (float)( *(MembType *)(p + addr[1] * stride) );
	const float f2 = (float)( *(MembType *)(p + addr[2] * stride) );
	const float f3 = (float)( *(MembType *)(p + addr[3] * stride) );
	const float f4 = (float)( *(MembType *)(p + addr[4] * stride) );
	const float f5 = (float)( *(MembType *)(p + addr[5] * stride) );
	const float f6 = (float)( *(MembType *)(p + addr[6] * stride) );
	const float f7 = (float)( *(MembType *)(p + addr[7] * stride) );

	return (v8sf){f0, f1, f2, f3, f4, f5, f6, f7};
}

template<typename PackType, typename MembType>
v8si v8si_struct_load_indirect(const PackType &, const MembType &base, const v8si vaddr){
	const char *p = (char *)(&base);
	const size_t stride = sizeof(PackType);

	int addr[8] __attribute__((aligned(32)));
	*(v8si *)addr = vaddr;

	const int i0 = (float)( *(MembType *)(p + addr[0] * stride) );
	const int i1 = (float)( *(MembType *)(p + addr[1] * stride) );
	const int i2 = (float)( *(MembType *)(p + addr[2] * stride) );
	const int i3 = (float)( *(MembType *)(p + addr[3] * stride) );
	const int i4 = (float)( *(MembType *)(p + addr[4] * stride) );
	const int i5 = (float)( *(MembType *)(p + addr[5] * stride) );
	const int i6 = (float)( *(MembType *)(p + addr[6] * stride) );
	const int i7 = (float)( *(MembType *)(p + addr[7] * stride) );

	return (v8si){i0, i1, i2, i3, i4, i5, i6, i7};
}

static void print_v8sf(const v8sf v){
	float f[8] __attribute__((aligned(32)));
	*(v8sf *)f = v;
	printf("{%f, %f, %f, %f}, {%f, %f, %f, %f}\n",
			f[0], f[1], f[2], f[3],
			f[4], f[5], f[6], f[7]);
}

static void print_v8si(const v8si v){
	int i[8] __attribute__((aligned(32)));
	*(v8si *)i = v;
	printf("{%d, %d, %d, %d}, {%d, %d, %d, %d}\n",
			i[0], i[1], i[2], i[3],
			i[4], i[5], i[6], i[7]);
}

struct Pack{
	int    i;
	float  f;
	double d;
};

int main(){
	static Pack pack[64];

	for(int i=0; i<64; i++){
		pack[i].i = i;
		pack[i].f = 0.1 * i;
		pack[i].d = 10. * i;
	}

	v8sf v0 = v8sf_struct_load(pack[0], pack[0].f);
	print_v8sf(v0);

	v8si pi = {3,1,4,1,5,9,2,6};
	v8sf v1 = v8sf_struct_load_indirect(pack[0], pack[0].f, pi);
	print_v8sf(v1);

	v8sf v2 = v8sf_struct_load_indirect(pack[0], pack[0].i, pi);
	print_v8sf(v2);

	v8sf v3 = v8sf_struct_load_indirect(pack[0], pack[0].d, pi);
	print_v8sf(v3);

	v8si v4 = v8si_struct_load_indirect(pack[0], pack[0].i, pi);
	print_v8si(v4);

	// test for mask scatter

	for(int n=0; n<4; n++){
		static int buf[256];
		memset(buf, -1, 1024);
		int val [8] __attribute__((aligned(32)));
		int loc [8] __attribute__((aligned(32)));
		int mask[8] __attribute__((aligned(32)));

		for(int i=0; i<8; i++){
			val [i] = i;
			loc [i] = lrand48() % 256;
			mask[i] = (lrand48()>>16)%4 ? -1 : 0;
		}

		v8si_mask_scatter(*(v8si *)val, *(v8si *)mask, buf, *(v8si *)loc);

		puts("");
		for(int i=0; i<8; i++){
			printf("i=%d, val=%d, loc=%3d, mask=%d, buf=%d\n",
					i, val[i], loc[i], mask[i], buf[loc[i]]);
		}
	}

	return 0;
}

#if 0
{
	const S32 ith = Comm::getThreadNum();
	const __m256  inv_theta_sq_ps   = _mm256_set1_ps( (float)(1.0/(theta_*theta_)) );
	const __m256  half_ps           = _mm256_set1_ps(0.5f);
	const __m256  zero_ps           = _mm256_set1_ps(0.0f);
	const __m256i zero_pi32         = _mm256_set1_epi32(0);
	const __m256i one_pi32          = _mm256_set1_epi32(1);
	const __m256i minus_one_pi32    = _mm256_set1_epi32(-1);
	const __m256  mask_abs_ps       = (__m256)_mm256_set1_epi32(0x7fffffff);
	const __m256i n_leaf_limit_pi32 = _mm256_set1_epi32(n_leaf_limit_+1); // plus one is needed
	const __m256i level_limit_pi32  = _mm256_set1_epi32(TREE_LEVEL_LIMIT);// plus one is not needed
	const __m256  all_bit_ps        = (__m256)_mm256_set1_epi32(0xffffffff);
	const __m256i all_bit_pi32      = _mm256_set1_epi32(0xffffffff);
	const __m256i msb_pi32          = _mm256_set1_epi32(0x80000000);
	S32 iw_head = (n_walk/n_thread)*(ith+0) + std::min((ith+0), n_walk%n_thread);
	S32 iw_tail = (n_walk/n_thread)*(ith+1) + std::min((ith+1), n_walk%n_thread);

	for(int iw=iw_head; iw<iw_tail; iw++){
		const S32 first_adr_ip = ipg_[walk_grp_head+iw].adr_ptcl_; 
		n_epi_array[iw] = ipg_[walk_grp_head+iw].n_ptcl_;
		epi_array  [iw] = epi_sorted_.getPointer(first_adr_ip);
		force_array[iw] = force_sorted_.getPointer(first_adr_ip);
	}
	//if(ith==1) std::cerr<<"ith= "<<ith<<" iw_head= "<<iw_head<<" iw_tail= "<<iw_tail<<" n_walk= "<<n_walk<<std::endl;
	//std::cerr<<"wg= "<<wg<<" iw_head= "<<iw_head<<" iw_tail= "<<iw_tail<<" n_walk= "<<n_walk<<std::endl;
	__m256i n_sp_pi32 = _mm256_set1_epi32(0);
	__m256i n_ep_pi32 = _mm256_set1_epi32(0);
	for(int iw=iw_head; iw<iw_tail; iw += size_vec){
		S32 id_ipg_base = iw + walk_grp_head;
		//if(ith==1) std::cerr<<"CHECK-1: iw= "<<iw<<std::endl;
		for(S32 k=0; k<size_vec; k++){
			epj_for_force_tmp[ith][k].clearSize();
			spj_for_force_tmp[ith][k].clearSize();
			adr_cell_recorder[ith][k].clearSize();
			record_type      [ith][k].clearSize();
		}

		const F64ort &vtx0 = ipg_[id_ipg_base+0].vertex_;
		const F64ort &vtx1 = ipg_[id_ipg_base+1].vertex_;
		const F64ort &vtx2 = ipg_[id_ipg_base+2].vertex_;
		const F64ort &vtx3 = ipg_[id_ipg_base+3].vertex_;
		const F64ort &vtx4 = ipg_[id_ipg_base+4].vertex_;
		const F64ort &vtx5 = ipg_[id_ipg_base+5].vertex_;
		const F64ort &vtx6 = ipg_[id_ipg_base+6].vertex_;
		const F64ort &vtx7 = ipg_[id_ipg_base+7].vertex_;
		
		const __m256 ipg_low_x = {
			(float)vtx0.low_.x, (float)vtx1.low_.x, (float)vtx2.low_.x, (float)vtx3.low_.x,
			(float)vtx4.low_.x, (float)vtx5.low_.x, (float)vtx6.low_.x, (float)vtx7.low_.x };
		const __m256 ipg_low_y = {
			(float)vtx0.low_.y, (float)vtx1.low_.y, (float)vtx2.low_.y, (float)vtx3.low_.y,
			(float)vtx4.low_.y, (float)vtx5.low_.y, (float)vtx6.low_.y, (float)vtx7.low_.y };
		const __m256 ipg_low_z = {
			(float)vtx0.low_.z, (float)vtx1.low_.z, (float)vtx2.low_.z, (float)vtx3.low_.z,
			(float)vtx4.low_.z, (float)vtx5.low_.z, (float)vtx6.low_.z, (float)vtx7.low_.z };

		const __m256 ipg_high_x = {
			(float)vtx0.high_.x, (float)vtx1.high_.x, (float)vtx2.high_.x, (float)vtx3.high_.x,
			(float)vtx4.high_.x, (float)vtx5.high_.x, (float)vtx6.high_.x, (float)vtx7.high_.x };
		const __m256 ipg_high_y = {
			(float)vtx0.high_.y, (float)vtx1.high_.y, (float)vtx2.high_.y, (float)vtx3.high_.y,
			(float)vtx4.high_.y, (float)vtx5.high_.y, (float)vtx6.high_.y, (float)vtx7.high_.y };
		const __m256 ipg_high_z = {
			(float)vtx0.high_.z, (float)vtx1.high_.z, (float)vtx2.high_.z, (float)vtx3.high_.z,
			(float)vtx4.high_.z, (float)vtx5.high_.z, (float)vtx6.high_.z, (float)vtx7.high_.z };

		__m256 cx = _mm256_mul_ps(_mm256_add_ps(ipg_high_x, ipg_low_x), half_ps);
		__m256 cy = _mm256_mul_ps(_mm256_add_ps(ipg_high_y, ipg_low_y), half_ps);
		__m256 cz = _mm256_mul_ps(_mm256_add_ps(ipg_high_z, ipg_low_z), half_ps);
		__m256 lx = _mm256_mul_ps(_mm256_sub_ps(ipg_high_x, ipg_low_x), half_ps);
		__m256 ly = _mm256_mul_ps(_mm256_sub_ps(ipg_high_y, ipg_low_y), half_ps);
		__m256 lz = _mm256_mul_ps(_mm256_sub_ps(ipg_high_z, ipg_low_z), half_ps);

		__m256i adr_cell_pi32 = zero_pi32;
		U32 finished[size_vec];
		for(S32 i=0; i<size_vec; i++) finished[i] = 0xffffffff;
		S32 size_vec_new = std::min((iw_tail-iw), size_vec);
		for(S32 i=0; i<size_vec_new; i++) finished[i] = 0;
		//__m256i is_finished = _mm256_set_epi32(finished[7], finished[6], finished[5], finished[4],
		//		finished[3], finished[2], finished[1], finished[0]);
		__m256i is_finished = _mm256_loadu_si256((__m256i const *)finished);

		int flag = 0;
		int n_rem;
		int adr_cell   [8] __attribute__(aligned(32));
		int ep_loc     [8] __attribute__(aligned(32));
		int sp_loc     [8] __attribute__(aligned(32));
		int adr_ptcl_tc[8] __attribute__(aligned(32));
		int inc_sp     [8] __attribute__(aligned(32));
		int n_epj_tmp  [8] __attribute__(aligned(32));
		int n_spj_tmp  [8] __attribute__(aligned(32));
		int rem_id     [8] __attribute__(aligned(32));
		__m256i n_epj_tmp_pi32 = zero_pi32;
		__m256i n_spj_tmp_pi32 = zero_pi32;
		for(S32 k=0; k<8; k++) n_epj_tmp[k] = n_spj_tmp[k] = 0;

		do{
			*(v8si*) adr_cell = (v8si) adr_cell_pi32;
			__m256 fulle_len_sq_ps = {
				tc_glb_[adr_cell[0]].full_len_sq_, tc_glb_[adr_cell[1]].full_len_sq_, 
				tc_glb_[adr_cell[2]].full_len_sq_, tc_glb_[adr_cell[3]].full_len_sq_, 
				tc_glb_[adr_cell[4]].full_len_sq_, tc_glb_[adr_cell[5]].full_len_sq_, 
				tc_glb_[adr_cell[6]].full_len_sq_, tc_glb_[adr_cell[7]].full_len_sq_}; 
			__m256 r_crit_sq_ps = _mm256_mul_ps(fulle_len_sq_ps, inv_theta_sq_ps);

			//__asm__("# comment: jx, jy, jz");
#if 1

			__m256 jx = {(float)tc_glb_[adr_cell[0]].mom_.pos.x, (float)tc_glb_[adr_cell[1]].mom_.pos.x,
				(float)tc_glb_[adr_cell[2]].mom_.pos.x, (float)tc_glb_[adr_cell[3]].mom_.pos.x,
				(float)tc_glb_[adr_cell[4]].mom_.pos.x, (float)tc_glb_[adr_cell[5]].mom_.pos.x,
				(float)tc_glb_[adr_cell[6]].mom_.pos.x, (float)tc_glb_[adr_cell[7]].mom_.pos.x};
			__m256 jy = {(float)tc_glb_[adr_cell[0]].mom_.pos.y, (float)tc_glb_[adr_cell[1]].mom_.pos.y,
				(float)tc_glb_[adr_cell[2]].mom_.pos.y, (float)tc_glb_[adr_cell[3]].mom_.pos.y,
				(float)tc_glb_[adr_cell[4]].mom_.pos.y, (float)tc_glb_[adr_cell[5]].mom_.pos.y,
				(float)tc_glb_[adr_cell[6]].mom_.pos.y, (float)tc_glb_[adr_cell[7]].mom_.pos.y};
			__m256 jz = {(float)tc_glb_[adr_cell[0]].mom_.pos.z, (float)tc_glb_[adr_cell[1]].mom_.pos.z,
				(float)tc_glb_[adr_cell[2]].mom_.pos.z, (float)tc_glb_[adr_cell[3]].mom_.pos.z,
				(float)tc_glb_[adr_cell[4]].mom_.pos.z, (float)tc_glb_[adr_cell[5]].mom_.pos.z,
				(float)tc_glb_[adr_cell[6]].mom_.pos.z, (float)tc_glb_[adr_cell[7]].mom_.pos.z};
#else
			const F32vec & mom0 = tc_glb_[adr_cell[0]].mom_.getPos();
			const F32vec & mom1 = tc_glb_[adr_cell[1]].mom_.getPos();
			const F32vec & mom2 = tc_glb_[adr_cell[2]].mom_.getPos();
			const F32vec & mom3 = tc_glb_[adr_cell[3]].mom_.getPos();
			const F32vec & mom4 = tc_glb_[adr_cell[4]].mom_.getPos();
			const F32vec & mom5 = tc_glb_[adr_cell[5]].mom_.getPos();
			const F32vec & mom6 = tc_glb_[adr_cell[6]].mom_.getPos();
			const F32vec & mom7 = tc_glb_[adr_cell[7]].mom_.getPos();

			__m256 jx = _mm256_set_m128_forgcc(_mm_loadu_ps(&mom4.x), 
					_mm_loadu_ps(&mom0.x));
			__m256 jy = _mm256_set_m128_forgcc(_mm_loadu_ps(&mom5.x), 
					_mm_loadu_ps(&mom1.x));
			__m256 jz = _mm256_set_m128_forgcc(_mm_loadu_ps(&mom6.x), 
					_mm_loadu_ps(&mom2.x));
			__m256 jw = _mm256_set_m128_forgcc(_mm_loadu_ps(&mom7.x), 
					_mm_loadu_ps(&mom3.x));

			transpose_4ymm(jx, jy, jz, jw);
#endif

			//__asm__("# comment: calc dr");

			__m256 dx = _mm256_and_ps(mask_abs_ps, (cx - jx));
			dx = dx - lx;
			dx = _mm256_max_ps(dx, zero_ps);
			__m256 dy = _mm256_and_ps(mask_abs_ps, (cy - jy));
			dy = dy - ly;
			dy = _mm256_max_ps(dy, zero_ps);
			__m256 dz = _mm256_and_ps(mask_abs_ps, (cz - jz));
			dz = dz - lz;
			dz = _mm256_max_ps(dz, zero_ps);
			__m256 dr_sq_ps = _mm256_fmadd_ps(dz, dz, _mm256_fmadd_ps(dy, dy, (dx*dx)));

			__m256 is_close  = _mm256_cmp_ps(dr_sq_ps, r_crit_sq_ps, 2); // 2: <=
			__m256 is_far    = _mm256_xor_ps(is_close, all_bit_ps);

			//__asm__("# comment: adr_next");
#if 1
			__m256i adr_next_pi32 = _mm256_set_epi32(tc_glb_[adr_cell[7]].next_adr_, tc_glb_[adr_cell[6]].next_adr_,
					tc_glb_[adr_cell[5]].next_adr_, tc_glb_[adr_cell[4]].next_adr_,
					tc_glb_[adr_cell[3]].next_adr_, tc_glb_[adr_cell[2]].next_adr_,
					tc_glb_[adr_cell[1]].next_adr_, tc_glb_[adr_cell[0]].next_adr_);

			__m256i adr_more_pi32 = _mm256_set_epi32(tc_glb_[adr_cell[7]].more_adr_, tc_glb_[adr_cell[6]].more_adr_,
					tc_glb_[adr_cell[5]].more_adr_, tc_glb_[adr_cell[4]].more_adr_,
					tc_glb_[adr_cell[3]].more_adr_, tc_glb_[adr_cell[2]].more_adr_,
					tc_glb_[adr_cell[1]].more_adr_, tc_glb_[adr_cell[0]].more_adr_);

			__m256i n_ptcl_pi32 = _mm256_set_epi32(tc_glb_[adr_cell[7]].n_ptcl_, tc_glb_[adr_cell[6]].n_ptcl_,
					tc_glb_[adr_cell[5]].n_ptcl_, tc_glb_[adr_cell[4]].n_ptcl_,
					tc_glb_[adr_cell[3]].n_ptcl_, tc_glb_[adr_cell[2]].n_ptcl_,
					tc_glb_[adr_cell[1]].n_ptcl_, tc_glb_[adr_cell[0]].n_ptcl_);

			__m256i level_pi32 = _mm256_set_epi32(tc_glb_[adr_cell[7]].level_, tc_glb_[adr_cell[6]].level_,
					tc_glb_[adr_cell[5]].level_, tc_glb_[adr_cell[4]].level_,
					tc_glb_[adr_cell[3]].level_, tc_glb_[adr_cell[2]].level_,
					tc_glb_[adr_cell[1]].level_, tc_glb_[adr_cell[0]].level_);
#else

			__m128i * next0 = (__m128i*)(&(tc_glb_[adr_cell[0]].next_adr_));
			__m128i * next1 = (__m128i*)(&(tc_glb_[adr_cell[1]].next_adr_));
			__m128i * next2 = (__m128i*)(&(tc_glb_[adr_cell[2]].next_adr_));
			__m128i * next3 = (__m128i*)(&(tc_glb_[adr_cell[3]].next_adr_));
			__m128i * next4 = (__m128i*)(&(tc_glb_[adr_cell[4]].next_adr_));
			__m128i * next5 = (__m128i*)(&(tc_glb_[adr_cell[5]].next_adr_));
			__m128i * next6 = (__m128i*)(&(tc_glb_[adr_cell[6]].next_adr_));
			__m128i * next7 = (__m128i*)(&(tc_glb_[adr_cell[7]].next_adr_));

			__m256i adr_next_pi32 = _mm256_set_m128i_forgcc(_mm_loadu_si128(next4), 
					_mm_loadu_si128(next0));
			__m256i adr_more_pi32 = _mm256_set_m128i_forgcc(_mm_loadu_si128(next5), 
					_mm_loadu_si128(next1));
			__m256i n_ptcl_pi32   = _mm256_set_m128i_forgcc(_mm_loadu_si128(next6), 
					_mm_loadu_si128(next2));
			__m256i level_pi32    = _mm256_set_m128i_forgcc(_mm_loadu_si128(next7), 
					_mm_loadu_si128(next3));
			transpose_4ymm(adr_next_pi32, adr_more_pi32, n_ptcl_pi32, level_pi32);
#endif

			//__asm__("# comment: calc adr");
			__m256 is_leaf   = (__m256)_mm256_or_si256(_mm256_cmpgt_epi32(n_leaf_limit_pi32, n_ptcl_pi32), _mm256_cmpgt_epi32(level_pi32, level_limit_pi32));
			__m256i go_next  = (__m256i)_mm256_or_ps(_mm256_and_ps(is_close, is_leaf), is_far);
			__m256i go_more  = _mm256_xor_si256(go_next, all_bit_pi32);

#if 0
			int is_far_mask = _mm256_movemask_ps(is_far);
			int is_close_mask = _mm256_movemask_ps( _mm256_and_ps(is_close, is_leaf) );

			for(int k=0; k<size_vec_new;  k++){
				if( (flag>>k & 0x1) != 1){
					if( (is_far_mask>>k & 0x1) ){
						spj_for_force_tmp[ith][k].increaseSize();
						spj_for_force_tmp[ith][k].back().copyFromMoment(tc_glb_[adr_cell[k]].mom_);
					}
					else if(is_close_mask>>k & 0x1){
						U32 adr_ptcl = tc_glb_[adr_cell[k]].adr_ptcl_;
						S32 n_ptcl = tc_glb_[adr_cell[k]].n_ptcl_;
						//std::cerr<<"n_ptcl= "<<n_ptcl<<std::endl;
						epj_for_force_tmp[ith][k].reserveEmptyAreaAtLeast( n_ptcl );
						spj_for_force_tmp[ith][k].reserveEmptyAreaAtLeast( n_ptcl );              
						for(S32 ip=0; ip<n_ptcl; ip++, adr_ptcl++){
							if( !GetMSB(tp_glb_[adr_ptcl].adr_ptcl_) ){
								epj_for_force_tmp[ith][k].push_back( epj_sorted_[adr_ptcl] );
							}
							else{
								spj_for_force_tmp[ith][k].push_back( spj_sorted_[adr_ptcl] );
							}
						}
					}
				}
			}
#else
			__m256i set_ptcl_in_cell = _mm256_andnot_si256(is_finished, (__m256i)is_far);
			__m256i sp_loc_offset_pi32 = _mm256_andnot_si256(set_ptcl_in_cell, one_pi32); // set -> loc_offset = 0, not set -> loc_offset = 1
			__m256i inc_sp_pi32 = _mm256_and_si256(set_ptcl_in_cell, one_pi32); // set -> loc_offset = 1, not set -> loc_offset = 0
			n_spj_tmp_pi32 = _mm256_add_epi32(n_spj_tmp_pi32, inc_sp_pi32); 

			/*
			   float jx_tmp[8];
			 *(v8sf*)jx_tmp = jx;

			 float jy_tmp[8];
			 *(v8sf*)jy_tmp = jy;

			 float jz_tmp[8];
			 *(v8sf*)jz_tmp = jz;

			 float jm_tmp[8];
			 *(v8sf*)jm_tmp = jm;
			 */
			__m256i sp_loc_pi32 = _mm256_add_epi32(n_sp_pi32, sp_loc_offset_pi32); // offset_pi32=1 -> sp is not substitue, 
			n_sp_pi32 = _mm256_add_epi32(n_sp_pi32, inc_sp_pi32);
			*(v8si*) sp_loc = (v8si)sp_loc_pi32;

			int mask_set_ptcl_in_cell = _mm256_movemask_ps((__m256)set_ptcl_in_cell);
			//__asm__("# comment: copy cell");
			for(int k=0; k<8;  k++){
				if( (mask_set_ptcl_in_cell>>k & 0x1) ){
					//__asm__("# comment: copy cell (copy sp)");
					spj_for_force_tmp[ith][k][sp_loc[k]].copyFromMoment(tc_glb_[adr_cell[k]].mom_);
				}
			}

			__m256i n_ptcl_tmp_pi32 = n_ptcl_pi32;
			__m256i loop_flag_in_pi32 = _mm256_andnot_si256(is_finished, (__m256i)_mm256_and_ps(is_close, is_leaf));
			__m256i adr_ptcl_tc_pi32 = _mm256_set_epi32( tc_glb_[adr_cell[7]].adr_ptcl_,  tc_glb_[adr_cell[6]].adr_ptcl_,
					tc_glb_[adr_cell[5]].adr_ptcl_,  tc_glb_[adr_cell[4]].adr_ptcl_,
					tc_glb_[adr_cell[3]].adr_ptcl_,  tc_glb_[adr_cell[2]].adr_ptcl_,
					tc_glb_[adr_cell[1]].adr_ptcl_,  tc_glb_[adr_cell[0]].adr_ptcl_);
			int loop_flag_in = _mm256_movemask_ps((__m256)loop_flag_in_pi32);

			//__asm__("# comment: copy leaf begin");

			//hozon
			while(loop_flag_in != 0){
				__m256i adr_ptcl_tp_pi32 = _mm256_set_epi32( tp_glb_[tc_glb_[adr_cell[7]].adr_ptcl_].adr_ptcl_,  tp_glb_[tc_glb_[adr_cell[6]].adr_ptcl_].adr_ptcl_,
						tp_glb_[tc_glb_[adr_cell[5]].adr_ptcl_].adr_ptcl_,  tp_glb_[tc_glb_[adr_cell[4]].adr_ptcl_].adr_ptcl_,
						tp_glb_[tc_glb_[adr_cell[3]].adr_ptcl_].adr_ptcl_,  tp_glb_[tc_glb_[adr_cell[2]].adr_ptcl_].adr_ptcl_,
						tp_glb_[tc_glb_[adr_cell[1]].adr_ptcl_].adr_ptcl_,  tp_glb_[tc_glb_[adr_cell[0]].adr_ptcl_].adr_ptcl_); 

				__m256i is_ep = _mm256_and_si256(_mm256_cmpeq_epi32(_mm256_and_si256(adr_ptcl_tp_pi32, msb_pi32), zero_pi32), loop_flag_in_pi32);
				__m256i is_sp = _mm256_and_si256(_mm256_xor_si256(is_ep, all_bit_pi32), loop_flag_in_pi32);
				__m256i ep_loc_offset_pi32 = _mm256_andnot_si256( _mm256_and_si256(is_ep, loop_flag_in_pi32),    one_pi32);
				__m256i ep_loc_pi32 = _mm256_add_epi32(n_ep_pi32, ep_loc_offset_pi32);
				__m256i inc_ep_pi32 = _mm256_and_si256( _mm256_and_si256(is_ep, loop_flag_in_pi32),    one_pi32);
				n_ep_pi32 = _mm256_add_epi32(n_ep_pi32, inc_ep_pi32);
				sp_loc_offset_pi32 = _mm256_andnot_si256( _mm256_andnot_si256(is_sp, loop_flag_in_pi32), one_pi32);
				sp_loc_pi32 = _mm256_add_epi32(n_sp_pi32, sp_loc_offset_pi32);
				inc_sp_pi32 = _mm256_and_si256( _mm256_and_si256(is_sp, loop_flag_in_pi32),    one_pi32);
				n_sp_pi32 = _mm256_add_epi32(n_sp_pi32, inc_sp_pi32);

				n_epj_tmp_pi32 = _mm256_add_epi32(n_epj_tmp_pi32, inc_ep_pi32); 
				n_spj_tmp_pi32 = _mm256_add_epi32(n_spj_tmp_pi32, inc_sp_pi32); 

				*(v8si*)ep_loc      = (v8si)ep_loc_pi32;
				*(v8si*)sp_loc      = (v8si)sp_loc_pi32;
				*(v8si*)adr_ptcl_tc = (v8si)adr_ptcl_tc_pi32;
				*(v8si*)inc_sp      = (v8si)inc_sp_pi32;
				//__asm__("# comment: copy leaf");
				for(int k=0; k<size_vec_new;  k++){
					if( (loop_flag_in>>k & 0x1) ){
						//__asm__("# comment: copy leaf (ep, sp)");
						if(inc_sp[k] == 1){
							spj_for_force_tmp[ith][k][sp_loc[k]] = spj_sorted_[adr_ptcl_tc[k]];
						}
						else{
							epj_for_force_tmp[ith][k][ep_loc[k]] = epj_sorted_[adr_ptcl_tc[k]];
						}
					}
				}
				adr_ptcl_tc_pi32 = _mm256_add_epi32(adr_ptcl_tc_pi32, one_pi32);
				n_ptcl_tmp_pi32  = _mm256_sub_epi32(n_ptcl_tmp_pi32,  one_pi32 );
				loop_flag_in_pi32 = _mm256_and_si256(loop_flag_in_pi32, _mm256_cmpgt_epi32(n_ptcl_tmp_pi32, zero_pi32));
				loop_flag_in = _mm256_movemask_ps((__m256)loop_flag_in_pi32);
			} // end of while loop

			/*
			   rem_id[0] = 0;
			   rem_id[1] = 1;
			   rem_id[2] = 2;
			   rem_id[3] = 3;
			   rem_id[4] = 4;
			   rem_id[5] = 5;
			   rem_id[6] = 6;
			   rem_id[7] = 7;
			   n_rem = 8;
			   while(loop_flag_in != 0){
			   __m256i adr_ptcl_tp_pi32 = _mm256_set_epi32( tp_glb_[tc_glb_[adr_cell[7]].adr_ptcl_].adr_ptcl_,  tp_glb_[tc_glb_[adr_cell[6]].adr_ptcl_].adr_ptcl_,
			   tp_glb_[tc_glb_[adr_cell[5]].adr_ptcl_].adr_ptcl_,  tp_glb_[tc_glb_[adr_cell[4]].adr_ptcl_].adr_ptcl_,
			   tp_glb_[tc_glb_[adr_cell[3]].adr_ptcl_].adr_ptcl_,  tp_glb_[tc_glb_[adr_cell[2]].adr_ptcl_].adr_ptcl_,
			   tp_glb_[tc_glb_[adr_cell[1]].adr_ptcl_].adr_ptcl_,  tp_glb_[tc_glb_[adr_cell[0]].adr_ptcl_].adr_ptcl_); 
			   __m256i is_ep = _mm256_and_si256(_mm256_cmpeq_epi32(_mm256_and_si256(adr_ptcl_tp_pi32, msb_pi32), zero_pi32), loop_flag_in_pi32);
			   __m256i is_sp = _mm256_and_si256(_mm256_xor_si256(is_ep, all_bit_pi32), loop_flag_in_pi32);
			   __m256i ep_loc_offset_pi32 = _mm256_andnot_si256( _mm256_and_si256(is_ep, loop_flag_in_pi32),    one_pi32);
			   __m256i ep_loc_pi32 = _mm256_add_epi32(n_ep_pi32, ep_loc_offset_pi32);
			   __m256i inc_ep_pi32 = _mm256_and_si256( _mm256_and_si256(is_ep, loop_flag_in_pi32),    one_pi32);
			   n_ep_pi32 = _mm256_add_epi32(n_ep_pi32, inc_ep_pi32);
			   sp_loc_offset_pi32 = _mm256_andnot_si256( _mm256_andnot_si256(is_sp, loop_flag_in_pi32), one_pi32);
			   sp_loc_pi32 = _mm256_add_epi32(n_sp_pi32, sp_loc_offset_pi32);
			   inc_sp_pi32 = _mm256_and_si256( _mm256_and_si256(is_sp, loop_flag_in_pi32),    one_pi32);
			   n_sp_pi32 = _mm256_add_epi32(n_sp_pi32, inc_sp_pi32);

			   n_epj_tmp_pi32 = _mm256_add_epi32(n_epj_tmp_pi32, inc_ep_pi32); 
			   n_spj_tmp_pi32 = _mm256_add_epi32(n_spj_tmp_pi32, inc_sp_pi32); 

			 *(v8si*)ep_loc      = (v8si)ep_loc_pi32;
			 *(v8si*)sp_loc      = (v8si)sp_loc_pi32;
			 *(v8si*)adr_ptcl_tc = (v8si)adr_ptcl_tc_pi32;
			 *(v8si*)inc_sp      = (v8si)inc_sp_pi32;
			//__asm__("# comment: copy leaf");
			int n_rem_new = 0;
			for(int k=0; k<n_rem;  k++){
			if( (loop_flag_in>>rem_id[k] & 0x1) ){
			//__asm__("# comment: copy leaf (ep, sp)");
			if(inc_sp[k] == 1){
			spj_for_force_tmp[ith][k][sp_loc[k]] = spj_sorted_[adr_ptcl_tc[k]];
			}
			else{
			epj_for_force_tmp[ith][k][ep_loc[k]] = epj_sorted_[adr_ptcl_tc[k]];
			}
			rem_id[n_rem_new] = rem_id[k];
			n_rem_new++;
			}
			}
			n_rem = n_rem_new;
			adr_ptcl_tc_pi32 = _mm256_add_epi32(adr_ptcl_tc_pi32, one_pi32);
			n_ptcl_tmp_pi32  = _mm256_sub_epi32(n_ptcl_tmp_pi32,  one_pi32 );
			loop_flag_in_pi32 = _mm256_and_si256(loop_flag_in_pi32, _mm256_cmpgt_epi32(n_ptcl_tmp_pi32, zero_pi32));
			loop_flag_in = _mm256_movemask_ps((__m256)loop_flag_in_pi32);
			} // end of while loop
			 */


#endif
			go_next = _mm256_andnot_si256(is_finished, go_next);
			go_more = _mm256_andnot_si256(is_finished, go_more);
			adr_cell_pi32 = _mm256_add_epi32(_mm256_and_si256(adr_next_pi32, go_next), _mm256_and_si256(adr_more_pi32, go_more));
			is_finished = _mm256_or_si256(is_finished, _mm256_cmpeq_epi32(adr_cell_pi32, minus_one_pi32));
			adr_cell_pi32 = _mm256_andnot_si256(is_finished, adr_cell_pi32);
			flag = _mm256_movemask_ps((__m256)is_finished);

		}while(flag != 0xff);

		*(v8si*)n_epj_tmp = (v8si)n_epj_tmp_pi32;
		*(v8si*)n_spj_tmp = (v8si)n_spj_tmp_pi32;
		for(S32 k=0; k<size_vec_new; k++){
			/*
			   S32 ep_loc = n_epj_disp_thread[ith][cnt_thread[ith]];
			   S32 sp_loc = n_spj_disp_thread[ith][cnt_thread[ith]];
			   for(S32 k2=0; k2<n_epj_tmp[k]; k2++, ep_loc++){
			   epj_for_force_[ith][ep_loc] = epj_for_force_tmp[ith][k][k2];
			   }
			   for(S32 k2=0; k2<n_spj_tmp[k]; k2++, sp_loc++){
			   spj_for_force_[ith][sp_loc] = spj_for_force_tmp[ith][k][k2];
			   }
			 */
			epj_array[iw+k] = epj_for_force_tmp[ith][k].getPointer();
			spj_array[iw+k] = spj_for_force_tmp[ith][k].getPointer();
			n_epj_array[iw+k] = n_epj_tmp[k];
			n_spj_array[iw+k] = n_spj_tmp[k];
			n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_epj_array[iw+k] + n_epj_disp_thread[ith][cnt_thread[ith]];
			n_spj_disp_thread[ith][cnt_thread[ith]+1] = n_spj_array[iw+k] + n_spj_disp_thread[ith][cnt_thread[ith]];
			n_interaction_ep_ep_array[ith] += ((S64)n_epj_array[iw+k]*(S64)n_epi_array[iw+k]);
			n_interaction_ep_sp_array[ith] += ((S64)n_spj_array[iw+k]*(S64)n_epi_array[iw+k]);
			iw2ith[iw+k] = ith;
			iw2cnt[iw+k] = cnt_thread[ith];
			cnt_thread[ith]++;
		}
	}
}
#endif
