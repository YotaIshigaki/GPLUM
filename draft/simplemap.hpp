#pragma once
#include <unistd.h>
//#include "samplesortlib.hpp"
//
// simplemap.hpp
//
// a simple (but OMP-friendly) replacement of std::map
//
// provide:
//
//#define SIMPLEMAP_RANGE_CHECK	    
namespace SimpleMapLib{

    const int slack=1024;
    template <typename KeyT, typename ValT>
    class Map{
    public:

	class KeyValuePair{
	public:
	    KeyT key;
	    ValT value;
	};
	int32_t size;
	int32_t memsize;
	KeyValuePair * p;
	bool sorted;
	void initialize(ValT s){
	    if (s>0){
		p = new KeyValuePair[s+slack];
		size=s;
		memsize = s+slack;
	    }else{
		p= NULL;
		size=0;
		memsize=0;
	    }
	    sorted = false;
	}
    
	Map(ValT s){ initialize(s);	}
	Map(){ initialize(0);	}
	~Map(){
	    if (p!=NULL)  delete[] p;
	}
	void resize(ValT s){
	    sorted = false;
	    size=s;
	    if (s >memsize){
		if (p!=NULL) delete[] p;
		p = new KeyValuePair[s+slack];
		memsize = s+slack;
	    }
	}
	void clear(){
	    sorted=false;
	    size=0;
	}
		
		
	void makemap()
	{
#if 0	    
	    printf("data before sort\n");
	    for(auto i=0;i<size; i++){
		printf("it=%d i= %d  key=%d val=%d\n",
		       omp_get_thread_num(),
		       i, (int) p[i].key, (int) p[i].value);
	    }
#endif	    
		
	    SampleSortLib::samplesort(p, size,[](KeyValuePair & l)
				->auto{return l.key;});
	    sorted=true;
	}
	void set(KeyT  k, ValT v)
	{
#ifdef SIMPLEMAP_RANGE_CHECKE	    
	    if (v >= size){
		printf("simplemap::set failed, too large index:%d size:%d %d\n",
		       v, size, memsize);
		exit(-1);
	    }
#endif	    
	    p[v].key=k;
	    p[v].value=v;
	}

	ValT find(KeyT k){
	    if (!sorted)makemap();
		
	    if (k < p[0].key || k > p[size-1].key ){
		return (ValT)(size);
	    }
	    ValT low=0;
	    ValT high=size-1;
	    auto mid = (low + high)/2;
	    while (low < high-1){
		if (p[mid].key < k){
		    low=mid;
		}else{
		    high=mid;
		}
		mid = (low + high)/2;
		if (p[mid].key == k){
		    return  (ValT) mid;
		}
	    }
	    if (p[high].key == k){
		return  (ValT)high;
	    }else if (p[low].key == k){
		return  (ValT)low;
	    }
	    return (ValT)(size);

	}
	ValT value(KeyT k){return p[find(k)].value;}
	ValT at(KeyT k){return p[find(k)].value;}
	ValT second(ValT v){return p[v].value;}
	ValT end(){return (ValT)size;}
	ValT operator [](KeyT  k){return at(k);}
	void dump()
	{
	    printf("Size: %d\n", size);
	    for(auto i=0;i<size; i++){
		printf("%lld: %d\n", (int64_t)p[i].key, p[i].value);
	    }
	}
    };
}
	    

    
