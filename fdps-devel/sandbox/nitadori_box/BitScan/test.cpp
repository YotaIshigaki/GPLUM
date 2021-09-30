#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <sys/time.h>

static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return double(tv.tv_sec) + 1.e-6*double(tv.tv_usec);
}

struct MyBitArray0{
	size_t N;
	bool *array;

	void initialize(const size_t _N){
		N = _N;
		array = new bool[N];
	}
	~MyBitArray0(){
		if(array) delete [] array;
	}
	void set_true_at(size_t loc){
		array[loc] = true;
	}

	__attribute__((noinline))
	void scan_all(std::vector<int> &list) const {
		size_t N = this->N;
		int ntrue = 0;
		for(size_t i=0; i<N; i++){
			if(array[i]) list[ntrue++] = i;
		}
	}
};

struct MyBitArray1{
	std::vector<bool> array;

	void initialize(const size_t N){
		array.clear();
		array.resize(N, false);
	}

	void set_true_at(size_t loc){
		array[loc] = true;
	}

	__attribute__((noinline))
	void scan_all(std::vector<int> &list) const {
		size_t N = array.size();
#if 0
		list.clear();
		for(size_t i=0; i<N; i++){
			if(array[i]) list.push_back(i);
		}
#else
		int ntrue = 0;
		for(size_t i=0; i<N; i++){
			if(array[i]) list[ntrue++] = i;
		}
#endif
	}
};

template<bool USE_POPCNT>
struct MyBitArray2{
	typedef unsigned long u64;
	std::vector<u64> array;

	void initialize(const size_t N){
		assert(8 == sizeof(u64));
		size_t len = N/8 + (N%8 ? 1: 0);
		array.clear();
		array.resize(len, 0);
	}

	void set_true_at(size_t loc){
		size_t ah = loc / 64;
		size_t al = loc % 64;
		array[ah] |= 1L << al;
	}

	__attribute__((noinline))
	void scan_all(std::vector<int> &list) const {
		int len = array.size();
		int ntrue = 0;
		for(int i=0, base=0; i<len; i++, base+=64){
			if(array[i]){
				u64 tmp = array[i];
				if(USE_POPCNT){
					const int nbit = __builtin_popcountl(tmp);
					for(int j=0; j<nbit; j++){
						int ctz = __builtin_ctzl(tmp);

						list[ntrue++] = base + ctz;

						// tmp ^= 1L<<ctz;
						tmp &= (tmp-1); // BLSR (reset lower set bit) instruction
					}
				}else{
					for(int j=0; j<8; j++){
						u64 bit8 = 0xff & tmp;
						if(bit8){
							const int loc = base + 8*j;
							if(bit8 & (1<<0)) list[ntrue++] = loc+0;
							if(bit8 & (1<<1)) list[ntrue++] = loc+1;
							if(bit8 & (1<<2)) list[ntrue++] = loc+2;
							if(bit8 & (1<<3)) list[ntrue++] = loc+3;
							if(bit8 & (1<<4)) list[ntrue++] = loc+4;
							if(bit8 & (1<<5)) list[ntrue++] = loc+5;
							if(bit8 & (1<<6)) list[ntrue++] = loc+6;
							if(bit8 & (1<<7)) list[ntrue++] = loc+7;
						}
						tmp >>= 8;
					}
				}
			}
		}
	}
};

template <typename BitArray>
void benchmark(
		const unsigned long     NBITS,
		const int               NLOOP,
		const std::vector<int> &inp,
		std::vector<int>       &out)
{
	BitArray ba;
	ba.initialize(NBITS);

	const int Ntrue = inp.size();
	for(int i=0; i<Ntrue; i++){
	   	ba.set_true_at(inp[i]);
	}

	const double t0 = get_wtime();
	for(int k=0; k<NLOOP; k++){
		ba.scan_all(out);
	}
	const double t1 = get_wtime();

	assert(inp.size() == out.size());
	for(int i=0; i<Ntrue; i++){
		// printf("(%d, %d)\n", inp[i], out[i]);
		assert(inp[i] == out[i]);
	}

	// printf("%e sec\n", t1-t0);
	double scan_per_sec = double(NLOOP) * double(NBITS) / (t1-t0);
	double elem_per_sec = double(NLOOP) * double(Ntrue) / (t1-t0);
	printf("%s:\n", typeid(BitArray).name());
	printf("%e bit scan per second\n", scan_per_sec);
	printf("%e active bits per second\n", elem_per_sec);
}

int main(){
	enum{
		NBITS = 1<<24,
		NTRUE = 1<<20,
		NLOOP = 10,
	};

	srand48(20150426);
	std::vector<int> inp, out;
	inp.resize(NTRUE);
	for(int i=0; i<NTRUE; i++){
		inp[i] = lrand48() % NBITS;
	}
	std::sort(inp.begin(), inp.end());
	inp.erase(
			std::unique(inp.begin(), inp.end()),
			inp.end());
	printf("New size = %d\n", (int)inp.size());
	const int Ntrue = inp.size();

	// out.reserve(NTRUE);
	out.resize(Ntrue);

	benchmark<MyBitArray0> (NBITS, NLOOP, inp, out);
	benchmark<MyBitArray1> (NBITS, NLOOP, inp, out);
	benchmark<MyBitArray2<true>  > (NBITS, NLOOP, inp, out);
	benchmark<MyBitArray2<false> > (NBITS, NLOOP, inp, out);

	return 0;
}
