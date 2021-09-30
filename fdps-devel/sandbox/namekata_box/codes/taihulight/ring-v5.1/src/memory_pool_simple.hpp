#include<iostream>
#include<cstdlib>
#include<cassert>

namespace ParticleSimulator{
    class MemoryPool{
    private:
        enum{
            ALLIGE_SIZE = 8,
        };
        MemoryPool(){}
        ~MemoryPool(){}
        MemoryPool(const MemoryPool & mem);
        MemoryPool & operator = (const MemoryPool & mem);
        void * bottom_;
        void * top_;
        size_t cap_;
        size_t size_;
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
            getInstance().cap_ = cap_tmp;
            getInstance().size_ = 0;
            getInstance().bottom_ = malloc(getInstance().cap_);
            getInstance().top_ = getInstance().bottom_;
        }
        static void * alloc(const size_t _size){
            const size_t size_allign = getAllignSize(_size);
            void * top_prev = getInstance().top_;
            // add to tail
            if(getInstance().cap_ <= (getInstance().size_ + size_allign)){
                std::cerr<<"cap_= "<<getInstance().cap_
                         <<" size_= "<<getInstance().size_
                         <<" size_allign= "<<size_allign
                         <<std::endl;
            }
            assert(getInstance().cap_ > getInstance().size_ + size_allign);
            getInstance().top_   = ((char*)getInstance().top_) + size_allign;
            getInstance().size_ += size_allign;
            return top_prev;
        }
        static void free(){
            getInstance().size_ = 0;
            getInstance().top_ = getInstance().bottom_;
        }
        static void dump(){
            std::cerr<<"bottom_= "<<getInstance().bottom_<<std::endl;
            std::cerr<<"top_= "<<getInstance().top_<<std::endl;
            std::cerr<<"cap_= "<<getInstance().cap_<<std::endl;
            std::cerr<<"size_= "<<getInstance().size_<<std::endl;
        }

        static size_t getLeftSize(){
            return getInstance().cap_ - getInstance().size_;
        }
        
    };

    template<class T>
    class TempArray{
        TempArray(const TempArray &);
        TempArray & operator = (const TempArray &);
        //T * data_;
        size_t size_;
        size_t capacity_;
    public:
        T * data_;
        TempArray() : data_(NULL), size_(0), capacity_(0) {}
        void reserve(const size_t _cap){
            assert(data_ == NULL);
            capacity_ = _cap;
            size_ = 0;
            data_ = (T*)MemoryPool::alloc(sizeof(T)*capacity_);
        }
        void free(){

            if(data_ == NULL) return;
            capacity_ = 0;
            size_     = 0;
            data_     = NULL;
            MemoryPool::free();
        }

        const T & operator [] (const int i) const {
#ifdef SANITY_CHECK_TEMP_ARRAY
            if(i < 0 || capacity_ <= i || size_<= i ){
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

        T * getPointer(const int i=0) const {
#ifdef SANITY_CHECK_TEMP_ARRAY
            if(i < 0 || capacity_ <= i || size_<= i ){
                std::cerr<<"i= "<<i<<" capacity_= "<<capacity_<<" size_= "<<size_<<std::endl;
                std::cerr<<"function: "<<__FUNCTION__<<", line: "<<__LINE__<<", file: "<<__FILE__<<std::endl;
                std::cerr<<"typeid(T).name()"<<typeid(T).name()<<std::endl;
            }
            assert(i >= 0);
            assert(capacity_ > i);
            assert(size_ > i);
#endif
            return data_+i;
        }

        void resizeNoInitialize (const size_t n){
            assert(data_ == NULL);
            capacity_ = n;
            reserve(capacity_);
            size_ = n;
        }

        void pushBackNoCheck(const T & val){
            data_[size_++] = val;
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
            std::cerr<<"data_= "<<data_<<std::endl;
            MemoryPool::dump();
        }
        
        // special function
        void setSize (const size_t n){
            assert(n < capacity_);
            size_ = n;
        }
        size_t getLeftSizeOfMemoryPool(){
            return (MemoryPool::getLeftSize())/sizeof(T);
        }
    };
}
