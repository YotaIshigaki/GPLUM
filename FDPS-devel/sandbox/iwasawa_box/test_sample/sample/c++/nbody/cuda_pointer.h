#pragma once
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#if 0
#  include <cutil.h>
#else
#  include <helper_cuda.h>
#  define CUDA_SAFE_CALL checkCudaErrors
#endif

template <typename T>
struct cudaPointer{
	T *dev_pointer;
	T *host_pointer;
	size_t size;
	cudaPointer(){
		dev_pointer = NULL;
		host_pointer = NULL;
		size = 0;
	}
//        ~cudaPointer(){
//         free();
//      }
	void allocate(int _size){
		size = _size;
		void *p;
		CUDA_SAFE_CALL(cudaMalloc(&p, size * sizeof(T)));
		assert(p);
		dev_pointer = (T*)p;
		CUDA_SAFE_CALL(cudaMallocHost(&p, size * sizeof(T)));
		assert(p);
		host_pointer = (T*)p;
	}
	void free(){
		CUDA_SAFE_CALL(cudaFree(dev_pointer));
		CUDA_SAFE_CALL(cudaFreeHost(host_pointer));
		dev_pointer = NULL;
		host_pointer = NULL;
		size = 0;
	}

	void htod(int count){
#ifdef SANITY_CHECK_CUDA_POINTER
            assert(size_t(count) < size);
            //std::cerr<<"count= "<<count<<" size= "<<size<<std::endl;
#endif
		CUDA_SAFE_CALL(cudaMemcpy(dev_pointer, host_pointer, count * sizeof(T), cudaMemcpyHostToDevice));
	}
	void htod(){
		this->htod(size);
	}
	void dtoh(int count){
#ifdef SANITY_CHECK_CUDA_POINTER
            assert(size_t(count) < size);
            //std::cerr<<"count= "<<count<<" size= "<<size<<std::endl;
#endif
		CUDA_SAFE_CALL(cudaMemcpy(host_pointer, dev_pointer, count * sizeof(T), cudaMemcpyDeviceToHost));
	}
	void dtoh(){
		this->dtoh(size);
	}
    
    // async version
    void htod(int count, cudaStream_t & stream){
#ifdef SANITY_CHECK_CUDA_POINTER
        assert(size_t(count) < size);
        //std::cerr<<"count= "<<count<<" size= "<<size<<std::endl;
#endif
        CUDA_SAFE_CALL(cudaMemcpyAsync(dev_pointer, host_pointer, count * sizeof(T), cudaMemcpyHostToDevice, stream));
    }
    void htod(cudaStream_t & stream){
        this->htod(size, stream);
    }
    void dtoh(int count, cudaStream_t & stream){
#ifdef SANITY_CHECK_CUDA_POINTER
        assert(size_t(count) < size);
        //std::cerr<<"count= "<<count<<" size= "<<size<<std::endl;
#endif
        CUDA_SAFE_CALL(cudaMemcpyAsync(host_pointer, dev_pointer, count * sizeof(T), cudaMemcpyDeviceToHost, stream));
    }
    void dtoh(cudaStream_t & stream){
        this->dtoh(size, stream);
    }

    // async version
    void htod(int first_loc, int count, cudaStream_t & stream){
#ifdef SANITY_CHECK_CUDA_POINTER
        assert(size_t(count) < size);
#endif
        CUDA_SAFE_CALL(cudaMemcpyAsync(dev_pointer+first_loc, host_pointer+first_loc, 
                                       count * sizeof(T), cudaMemcpyHostToDevice, stream));
    }
    void dtoh(int first_loc, int count, cudaStream_t & stream){
#ifdef SANITY_CHECK_CUDA_POINTER
        assert(size_t(count) < size);
#endif
        CUDA_SAFE_CALL(cudaMemcpyAsync(host_pointer+first_loc, dev_pointer+first_loc, 
                                       count * sizeof(T), cudaMemcpyDeviceToHost, stream));
    }

	T &operator[] (int i){
		return host_pointer[i];
	}
	operator T* (){
		return dev_pointer;
	}
};
