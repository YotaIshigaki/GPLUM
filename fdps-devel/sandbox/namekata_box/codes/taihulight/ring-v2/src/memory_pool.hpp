#include<iostream>
#include<cstdlib>
#include<cassert>

namespace ParticleSimulator{
    class MemoryPool{
    private:
        enum{
            ALLIGE_SIZE       = 8,
            N_SEGMENT_LIMIT   = 10000,
        };
        MemoryPool(){}
        ~MemoryPool(){}
        MemoryPool(const MemoryPool & mem);
        MemoryPool & operator = (const MemoryPool & mem);
        void * head_;
        void * top_;
        size_t cap_;
        size_t size_;
        size_t n_segment_;
        size_t cap_per_seg_[N_SEGMENT_LIMIT];
        bool   used_per_seg_[N_SEGMENT_LIMIT];
        static MemoryPool & getInstance(){
            static MemoryPool inst;
            return inst;
        }
        static size_t getAllignSize(const size_t _size){
            return (((_size-1)/ALLIGE_SIZE)+1)*ALLIGE_SIZE;
        }
    public:
        static void initialize(const size_t _cap){
            size_t cap_tmp = getInstance().getAllignSize(_cap);
            getInstance().cap_ = _cap;
            getInstance().size_ = 0;
            getInstance().head_ = malloc(getInstance().cap_);
            getInstance().top_ = getInstance().head_;
            getInstance().n_segment_ = 0;
            for(size_t i=0; i<N_SEGMENT_LIMIT; i++){
                getInstance().cap_per_seg_[i] = 0;
                getInstance().used_per_seg_[i] = false;
            }
        }
        static void * alloc(const size_t _size, int & _id_mpool){

            const size_t size_allign = getAllignSize(_size);
            size_t cap_cum = 0;
            for(int i=0; i<getInstance().n_segment_; i++){
                if( !getInstance().used_per_seg_[i] && getInstance().cap_per_seg_[i] >= size_allign){
                    // insert to middle 
                    getInstance().used_per_seg_[i] = true;
                    _id_mpool = i;
                    return (void*)((char*)getInstance().head_ + cap_cum);
                }
                cap_cum += getInstance().cap_per_seg_[i];
            }
            // add to tail
            assert(N_SEGMENT_LIMIT > getInstance().n_segment_);
            void * top_prev = getInstance().top_;
            assert(getInstance().cap_ > getInstance().size_ + size_allign);
            getInstance().top_   = ((char*)getInstance().top_) + size_allign;
            getInstance().size_ += size_allign;
            getInstance().cap_per_seg_[getInstance().n_segment_] = size_allign;
            getInstance().used_per_seg_[getInstance().n_segment_] = true;
            _id_mpool = getInstance().n_segment_;
            getInstance().n_segment_++;
            return top_prev;
        }
        static void free(const int id_seg){
            getInstance().used_per_seg_[id_seg] = false;
            if(id_seg == getInstance().n_segment_-1){
                for(int i=id_seg; i>=0; i--){
                    if(getInstance().used_per_seg_[i] == true) break;
                    getInstance().size_ -= getInstance().cap_per_seg_[i];
                    getInstance().cap_per_seg_[i] = 0;
                    getInstance().n_segment_--;
                }
            }
            getInstance().top_ = ((char*)getInstance().head_) + getInstance().size_;
        }
        static void dump(){
            std::cerr<<"head_= "<<getInstance().head_<<std::endl;
            std::cerr<<"top_= "<<getInstance().top_<<std::endl;
            std::cerr<<"cap_= "<<getInstance().cap_<<std::endl;
            std::cerr<<"size_= "<<getInstance().size_<<std::endl;
            std::cerr<<"n_segment= "<<getInstance().n_segment_<<std::endl;
            for(int i=0; i<getInstance().n_segment_; i++){
                std::cerr<<"i= "<<i
                         <<" cap= "<<getInstance().cap_per_seg_[i]
                         <<" used= "<<getInstance().used_per_seg_[i]
                         <<std::endl;
            }
        }
    };

    template<class T>
    class TempArray{
        TempArray(const TempArray &);
        TempArray & operator = (const TempArray &);
        T * data_;
        size_t size_;
        size_t capacity_;
        int id_mpool_;
    public:
        TempArray() : data_(NULL), size_(0), capacity_(0) {}
        void reserve(const size_t _cap){
            assert(data_ == NULL);
            capacity_ = _cap;
            size_ = 0;
            data_ = (T*)MemoryPool::alloc(sizeof(T)*capacity_, id_mpool_);
            //MemoryPool::dump();
        }
        void free(){
            //std::cerr<<"check a"<<std::endl;
            if(data_ == NULL) return;
            capacity_ = 0;
            size_     = 0;
            data_     = NULL;
            //std::cerr<<"check b"<<std::endl;
            MemoryPool::free(id_mpool_);
            //std::cerr<<"check c"<<std::endl;
        }

        const T & operator [] (const int i) const {
#ifdef SANITY_CHECK_TEMP_ARRAY
            if(i < 0 || capacity_ <= i || size_<= i ){
                dumpImpl();
                std::cerr<<"i= "<<i<<" capacity_= "<<capacity_<<" size_= "<<size_<<std::endl;
                std::cerr<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
                std::cerr<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
            }
            assert(i >= 0);
            assert(capacity_ > i);
            assert(size_ > i);
#endif
            return data_[i];
        }

        T & operator [] (const int i){
#ifdef SANITY_CHECK_TEMP_ARRAY
            if(i < 0 || capacity_ <= i || size_<= i ){
                dumpImpl();
                std::cerr<<"i= "<<i<<" capacity_= "<<capacity_<<" size_= "<<size_<<std::endl;
                std::cerr<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
                std::cerr<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
            }
            assert(i >= 0);
            assert(capacity_ > i);
            assert(size_ > i);
#endif
            return data_[i];
        }

        T * getPointer(const int i=0) const { return data_+i; }

        void resizeNoInitialize (const size_t n){
            assert(data_ == NULL);
            capacity_ = (n + (n+3)/3);
            reserve(capacity_);
            size_ = n;
        }
        size_t capacity() const {
            return capacity_;
        }
        size_t size() const {
            return size_;
        }

        size_t getMemSize() const { return capacity_ * sizeof(T); }
        
        void dump(){
            std::cerr<<"capacity_= "<<capacity_<<std::endl;
            std::cerr<<"size_= "<<size_<<std::endl;
            std::cerr<<"id_mpool_= "<<id_mpool_<<std::endl;
            std::cerr<<"data_= "<<data_<<std::endl;
            MemoryPool::dump();
        }
    };
}
