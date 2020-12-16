#pragma once

#include<cassert>
#include<fstream>

#include "MT.hpp"
#include"ps_defs.hpp"

namespace ParticleSimulator{
    template<class Tptcl>
    class ParticleSystem{
    private:
        CountT n_ptcl_send_;
        CountT n_ptcl_recv_;
        TimeProfile time_profile_;
        const S32 n_smp_ave_ = 30;
        ReallocatableArray<Tptcl> ptcl_;
        ReallocatableArray<S32> idx_remove_ptcl_; // add 2016/09/14
        S32 n_smp_ptcl_tot_;
        bool first_call_by_initialize;
        bool first_call_by_setAverageTargetNumberOfSampleParticlePerProcess;
        bool first_call_by_DomainInfo_collect_sample_particle;
        CountT n_proc_exch_;
        inline bool determineWhetherParticleIsInDomain(const F64vec & pos,
                                                       const F64ort & domain) {
            bool ret = true;
            ret = ret && (domain.low_.x <= pos.x) && (pos.x < domain.high_.x);
            ret = ret && (domain.low_.y <= pos.y) && (pos.y < domain.high_.y);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            ret = ret && (domain.low_.z <= pos.z) && (pos.z < domain.high_.z);
#endif
            return ret;
        }

        S32 searchWhichDomainParticleGoTo(const F64vec & pos,
                                          const S32 n_domain [],
                                          const F64ort domain []) {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
            S32 idomain = 0;
            const S32 ny = n_domain[1];
            while(domain[idomain].high_.x <= pos.x)
                idomain += ny;
            while(domain[idomain].high_.y <= pos.y)
                idomain++;
            return idomain;
#else
            S32 idomain = 0;
            const S32 nynz = n_domain[1] * n_domain[2];
            while(domain[idomain].high_.x <= pos.x)
                idomain += nynz;
            const S32 nz   = n_domain[2];
            while(domain[idomain].high_.y <= pos.y)
                idomain += nz;
            while(domain[idomain].high_.z <= pos.z)
                idomain++;            
            return idomain;
#endif
        }

    public:
#ifdef TEST_VARIADIC_TEMPLATE
        void ParticleSystemDummyFunc(){}
#endif
        TimeProfile getTimeProfile() const {
            return time_profile_;
        }
        void clearTimeProfile(){
            time_profile_.clear();
        }
        ParticleSystem() {
            first_call_by_setAverageTargetNumberOfSampleParticlePerProcess = true;
            first_call_by_initialize = true;
            first_call_by_DomainInfo_collect_sample_particle = true;
            //n_smp_ptcl_tot_ = n_smp_ave_ * Comm::getNumberOfProc();
        }

        void initialize(const S32 cap=10000) {
            assert(first_call_by_initialize);
            first_call_by_initialize = false;
            n_smp_ptcl_tot_ = n_smp_ave_ * Comm::getNumberOfProc();
            n_ptcl_send_ = n_ptcl_recv_ = 0;
            MT::init_genrand(Comm::getRank()); // 2017.11.01
            createParticle(cap);
            //first_call_by_DomainInfo_collect_sample_particle = true;
            //n_smp_ptcl_tot_ = n_smp_ave_ * Comm::getNumberOfProc();
        }

        void setAverageTargetNumberOfSampleParticlePerProcess(const S32 &nsampleperprocess) {
            assert(first_call_by_setAverageTargetNumberOfSampleParticlePerProcess);
            first_call_by_setAverageTargetNumberOfSampleParticlePerProcess = false;
            n_smp_ptcl_tot_ = nsampleperprocess * Comm::getNumberOfProc();
        }

        S32 getTargetNumberOfSampleParticle() {
            return n_smp_ptcl_tot_;
        }

        bool getFirstCallByDomainInfoCollectSampleParticle() {
            if(first_call_by_DomainInfo_collect_sample_particle) {
                first_call_by_DomainInfo_collect_sample_particle = false;
                return true;
            } else {
                return false;
            }
        }

        void createParticle(const S32 n_limit, bool clear=true){
            //n_ptcl_limit_ = n_limit;
            //ptcl_ = new Tptcl[n_ptcl_limit_];
            ptcl_.reserve(n_limit);
            ptcl_.resizeNoInitialize(0);
        }


        //void setNumberOfParticleLocal(const S32 n){ n_ptcl_ = n; }
        void setNumberOfParticleLocal(const S32 n){
            //15/02/20 Hosono bug(?) fix.
            //ptcl_.reserve(n*3+1000);
            ptcl_.resizeNoInitialize(n);
        }
        ////////////////
        // 05/01/30 Hosono From
        ////////////////
        //dummy class for the case if User does NOT define the file header.
        //TO BE PRIVATE
        struct DummyHeader{
            void writeAscii(FILE* fp) const{
            }
            S32 readAscii (FILE* fp){
                return -1;
            }
            void writeBinary(FILE* fp) const{
            }
            S32 readBinary (FILE* fp){
                return -1;
            }
        };

        // fp must points to the head address of the first partticle
        S32 countPtclAscii(FILE * fp){
            S32 n = 0;
            for(int c ; (c = getc(fp)) != EOF ; n += '\n' == c ? 1 : 0){}
            return n;
        }
        S32 countPtclBinary(FILE * fp,
                            void (Tptcl::*pFuncPtcl)(FILE*)){
            S32 n = 0;
            long end_of_header = ftell(fp);
            Tptcl ptcl_tmp;
            (ptcl_tmp.*pFuncPtcl)(fp);
            long size_ptcl = ftell(fp) - end_of_header;
            fseek(fp, 0, SEEK_END);
            n = (ftell(fp)-end_of_header) / size_ptcl;
            return n;
        }
        //2016_11_03 modified IO functions to handle multiple file formats
        //read
        template <class Theader>
        void readParticleImpl(const char * const filename,
                              const char * const format,
                              Theader * const header,
                              void (Tptcl::*pFuncPtcl)(FILE*),
                              S32 (Theader::*pFuncHead)(FILE*),
                              const char * open_format){
            if(format == NULL){//Read from single file
                if(Comm::getRank() == 0){
                    FILE* fp = fopen(filename, open_format);
                    if(fp == NULL){
                        PARTICLE_SIMULATOR_PRINT_ERROR("can not open input file ");
                        std::cerr<<"filename: "<<filename<<std::endl;
                        Abort(-1);
                    }
                    S32 n_ptcl_ = (header->*pFuncHead)(fp);
                    while('\n' == getc(fp));
                    fseek(fp, -1, SEEK_CUR);
                    if(n_ptcl_ < 0){//User does NOT return # of ptcl
                        //n_ptcl_ = 0;
                        if(strcmp(open_format, "rb")==0){
                            n_ptcl_ = countPtclBinary(fp, pFuncPtcl);
                        }
                        else{
                            //count # of lines
                            //KN
                            n_ptcl_ = countPtclAscii(fp);
                        }
                        fseek(fp, 0, SEEK_SET);
                        (header->*pFuncHead)(fp);
                        while('\n' == getc(fp));
                        fseek(fp, -1, SEEK_CUR);
                    }
                    //Inform the # of ptcl for each process.
                    const S32 n_proc = Comm::getNumberOfProc();
                    S32 *n_ptcl = new S32[n_proc];
                    for(S32 i = 0 ; i < n_proc ; ++ i){
                        n_ptcl[i] = n_ptcl_ / n_proc;
                    }
                    n_ptcl[0] += n_ptcl_ % n_proc;
                    #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    MPI_Scatter(n_ptcl, 1, GetDataType<S32>(), &n_ptcl_, 1, GetDataType<S32>(), 0, MPI_COMM_WORLD);
                    #endif
                    //allocate ptcl.
                    //First of all, Rank 0 reads its own particle.
                    this->createParticle(n_ptcl_ << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_);
                    for(int i = 0 ; i < n_ptcl_ ; ++ i){
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    //Read remaining data to buffer and send them to appropriate process.
                    #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    for(S32 rank = 1 ; rank < n_proc ; ++ rank){
                        Tptcl * buffer = new Tptcl[n_ptcl[rank]];
                        for(int i = 0 ; i < n_ptcl[rank] ; ++ i){
                            (buffer[i].*pFuncPtcl)(fp);
                        }
                        MPI_Send(buffer, n_ptcl[rank], GetDataType<Tptcl>(), rank, 0, MPI_COMM_WORLD);
                        delete [] buffer;
                    }
                    #endif
                    //End.
                    delete [] n_ptcl;
                    fclose(fp);
                }else{
                    #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    //Receive the # of ptcl from Rank 0
                    S32 n_ptcl_loc;
                    S32 *n_ptcl = new S32[Comm::getNumberOfProc()];
                    MPI_Scatter(n_ptcl, 1, GetDataType<S32>(), &n_ptcl_loc, 1, GetDataType<S32>(), 0, MPI_COMM_WORLD);
                    delete [] n_ptcl;
                    //allocate ptcl.
                    this->createParticle(n_ptcl_loc << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_loc);
                    MPI_Status stat;
                    MPI_Recv(ptcl_.getPointer(), ptcl_.size(), GetDataType<Tptcl>(), 0, 0, MPI_COMM_WORLD, &stat);
                    #endif
                }
            }else{//Read from multiple file
               char input[256];
                sprintf(input, format, filename, Comm::getNumberOfProc(), Comm::getRank());
                FILE* fp = fopen(input, open_format);
                if(fp == NULL){
                    PARTICLE_SIMULATOR_PRINT_ERROR("can not open input file");
                    std::cerr<<"filename: "<<input<<std::endl;
                    Abort(-1);
                }
                S32 n_ptcl_ = (header->*pFuncHead)(fp);
                while('\n' == getc(fp));
                fseek(fp, -1, SEEK_CUR);
                if(n_ptcl_ < 0){//User does NOT return # of ptcl
                    if(strcmp(open_format, "rb")==0){
                        n_ptcl_ = countPtclBinary(fp, pFuncPtcl);
                    }
                    else{
                        n_ptcl_ = countPtclAscii(fp);
                    }
                    this->createParticle(n_ptcl_ << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_);
                    fseek(fp, 0, SEEK_SET);
                    //S32 tmp = (header->*pFuncHead)(fp);
                    while('\n' == getc(fp));
                    fseek(fp, -1, SEEK_CUR);                    
                    for(S32 i = 0 ; i < ptcl_.size() ; ++ i){
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    fclose(fp);
                }
                else{
                    //User returns # of ptcl.
                    this->createParticle(n_ptcl_ << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_);
                    for(S32 i = 0 ; i < n_ptcl_ ; ++ i){
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    fclose(fp);
                }
            }
        }        

        template <class Theader>
        void readParticleAscii(const char * const filename, const char * const format, Theader& header){
            readParticleImpl(filename, format, &header, &Tptcl::readAscii, &Theader::readAscii, "r");
        }
        template <class Theader>
        void readParticleAscii(const char * const filename, Theader& header){
            readParticleImpl(filename, NULL, &header, &Tptcl::readAscii, &Theader::readAscii, "r");
        }
        void readParticleAscii(const char * const filename, const char * const format){
            readParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::readAscii, &DummyHeader::readAscii, "r");
        }
        void readParticleAscii(const char * const filename){
            readParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::readAscii, &DummyHeader::readAscii, "r");
        }

        template <class Theader>
        void readParticleAscii(const char * const filename, const char * const format, Theader& header, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl(filename, format, &header, pFunc, &Theader::readAscii, "r");
        }
        template <class Theader>
        void readParticleAscii(const char * const filename, Theader& header, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl(filename, NULL, &header, pFunc, &Theader::readAscii, "r");
        }
        void readParticleAscii(const char * const filename, const char * const format, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::readAscii, "r");
        }
        void readParticleAscii(const char * const filename, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::readAscii, "r");
        }

        template <class Theader>
        void readParticleBinary(const char * const filename, const char * const format, Theader& header){
            readParticleImpl(filename, format, &header, &Tptcl::readBinary, &Theader::readBinary, "rb");
        }
        template <class Theader>
        void readParticleBinary(const char * const filename, Theader& header){
            readParticleImpl(filename, NULL, &header, &Tptcl::readBinary, &Theader::readBinary, "rb");
        }
        void readParticleBinary(const char * const filename, const char * const format){
            readParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::readBinary, &DummyHeader::readBinary, "rb");
        }
        void readParticleBinary(const char * const filename){
            readParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::readBinary, &DummyHeader::readBinary, "rb");
        }
        
        template <class Theader>
        void readParticleBinary(const char * const filename, const char * const format, Theader& header, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl(filename, format, &header, pFunc, &Theader::readBinary, "rb");
        }
        template <class Theader>
        void readParticleBinary(const char * const filename, Theader& header, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl(filename, NULL, &header, pFunc, &Theader::readBinary, "rb");
        }
        void readParticleBinary(const char * const filename, const char * const format, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::readBinary, "rb");
        }
        void readParticleBinary(const char * const filename, void (Tptcl::*pFunc)(FILE*)){
            readParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::readBinary, "rb");
        }

        template <class Theader>
        void writeParticleImpl(const char * const filename,
                               const char * const format,
                               const Theader * const header,
                               void (Tptcl::*pFuncPtcl)(FILE*)const,
                               void (Theader::*pFuncHead)(FILE*)const,
                               const char * open_format){
            if(format == NULL){
                #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                //declare local # of ptcl.
                const S32 n_ptcl_ = ptcl_.size();
                //get # of process.
                const S32 n_proc = Comm::getNumberOfProc();
                //get # of ptcls in each process.
                S32 *n_ptcl = new S32[n_proc];
                //Gather # of particles.
                //MPI_Allgather(&n_ptcl_, 1, GetDataType<S32>(), n_ptcl, 1, GetDataType<S32>(), MPI_COMM_WORLD);
                MPI_Allgather(const_cast<S32*>(&n_ptcl_), 1, GetDataType<S32>(), n_ptcl, 1, GetDataType<S32>(), MPI_COMM_WORLD);
                //set displacement
                S32 *n_ptcl_displs = new S32[n_proc+1];
                n_ptcl_displs[0] = 0;
                for(S32 i = 0 ; i < n_proc ; ++ i){
                    n_ptcl_displs[i+1] = n_ptcl_displs[i] + n_ptcl[i];
                }
                const S32 n_tot = n_ptcl_displs[n_proc];
                Tptcl *ptcl = new Tptcl[n_tot];
                //gather data
                MPI_Gatherv(ptcl_.getPointer(), n_ptcl_, GetDataType<Tptcl>(), ptcl, n_ptcl, n_ptcl_displs, GetDataType<Tptcl>(), 0, MPI_COMM_WORLD);
                if(Comm::getRank() == 0){
                    FILE* fp = fopen(filename, open_format);
                    if(fp == NULL){
                        PARTICLE_SIMULATOR_PRINT_ERROR("can not open output file");
                        std::cerr<<"output file: "<<filename<<std::endl;
                        Abort(-1);
                    }
                    //header->writeAscii(fp);
                    (header->*pFuncHead)(fp);
                    for(S32 i = 0 ; i < n_tot ; ++ i){
                        //ptcl[i].writeAscii(fp);
                        (ptcl[i].*pFuncPtcl)(fp);
                    }
                    fclose(fp);
                }
                delete [] n_ptcl;
                delete [] n_ptcl_displs;
                delete [] ptcl;
                #else
                const S32 n_tot = ptcl_.size();
                if(Comm::getRank() == 0){
                    FILE* fp = fopen(filename, open_format);
                    //header->writeAscii(fp);
                    (header->*pFuncHead)(fp);
                    for(S32 i = 0 ; i < n_tot ; ++ i){
                        //ptcl_[i].writeAscii(fp);
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    fclose(fp);
                }
                #endif
            }else{
                char output[256];
                sprintf(output, format, filename, Comm::getNumberOfProc(), Comm::getRank());
                FILE* fp = fopen(output, open_format);
                if(fp == NULL){
                    PARTICLE_SIMULATOR_PRINT_ERROR("can not open output file");
                    std::cerr<<"output file: "<<output<<std::endl;
                    Abort(-1);
                }
                //header->writeAscii(fp);
                (header->*pFuncHead)(fp);
                for(S32 i = 0 ; i < ptcl_.size() ; ++ i){
                    //ptcl_[i].writeAscii(fp);
                    (ptcl_[i].*pFuncPtcl)(fp);
                }
                fclose(fp);
            }
        }
        //write
        template <class Theader>
        void writeParticleAscii(const char * const filename, const char * const format, const Theader& header){
            writeParticleImpl(filename, format, &header, &Tptcl::writeAscii, &Theader::writeAscii, "w");
        }
        template <class Theader>
        void writeParticleAscii(const char * const filename, const Theader& header){
            writeParticleImpl(filename, NULL, &header, &Tptcl::writeAscii, &Theader::writeAscii, "w");
        }
        void writeParticleAscii(const char * const filename, const char * format){
            writeParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::writeAscii, &DummyHeader::writeAscii, "w");
        }
        void writeParticleAscii(const char * const filename){
            writeParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::writeAscii, &DummyHeader::writeAscii, "w");
        }

        template <class Theader>
        void writeParticleAscii(const char * const filename, const char * const format, const Theader& header, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl(filename, format, &header, pFunc, &Theader::writeAscii, "w");
        }
        template <class Theader>
        void writeParticleAscii(const char * const filename, const Theader& header, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl(filename, NULL, &header, pFunc, &Theader::writeAscii, "w");
        }
        void writeParticleAscii(const char * const filename, const char * format, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::writeAscii, "w");
        }
        void writeParticleAscii(const char * const filename, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::writeAscii, "w");
        }


        template <class Theader>
        void writeParticleBinary(const char * const filename, const char * const format, const Theader& header){
            writeParticleImpl(filename, format, &header, &Tptcl::writeBinary, &Theader::writeBinary, "wb");
        }
        template <class Theader>
        void writeParticleBinary(const char * const filename, const Theader& header){
            writeParticleImpl(filename, NULL, &header, &Tptcl::writeBinary, &Theader::writeBinary, "wb");
        }
        void writeParticleBinary(const char * const filename, const char * format){
            writeParticleImpl<DummyHeader>(filename, format, NULL, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
        }
        void writeParticleBinary(const char * const filename){
            writeParticleImpl<DummyHeader>(filename, NULL, NULL, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
        }

        template <class Theader>
        void writeParticleBinary(const char * const filename, const char * const format, const Theader& header, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl(filename, format, &header, pFunc, &Theader::writeBinary, "wb");
        }
        template <class Theader>
        void writeParticleBinary(const char * const filename, const Theader& header, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl(filename, NULL, &header, pFunc, &Theader::writeBinary, "wb");
        }
        void writeParticleBinary(const char * const filename, const char * format, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &DummyHeader::writeBinary, "wb");
        }
        void writeParticleBinary(const char * const filename, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &DummyHeader::writeBinary, "wb");
        }
        

        

        Tptcl & operator [] (const S32 id) {return ptcl_[id];}
        const Tptcl & operator [] (const S32 id) const {return ptcl_[id];}
        Tptcl & getParticle(const S32 id=0) {return ptcl_[id];}
        const Tptcl & getParticle(const S32 id=0) const {return ptcl_[id];}
        Tptcl * getParticlePointer(const S32 id=0) const {return ptcl_+id;}
        //S32 getNumberOfParticleLocal() const {return n_ptcl_;}
        S32 getNumberOfParticleLocal() const {return ptcl_.size();}
        ////////////////
        // 05/02/04 Hosono From
        // 11th jul MI modified
        ////////////////
        S64 getNumberOfParticleGlobal() const {
            #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const S64 n_ptcl_loc = ptcl_.size();
            return Comm::getSum(n_ptcl_loc);
            #else
            return ptcl_.size();
            #endif
        }
        ////////////////
        // 05/02/04 Hosono To
        ////////////////

        F64 getHalfLength(const F64vec & center=0.0){
            F64 hl_max_loc = (ptcl_[0].getPos() - center).applyEach(Abs<F64>()).getMax();
            const S32 n_ptcl = ptcl_.size();
            for(size_t i=1; i<n_ptcl; i++){
                F64 hl_tmp = (ptcl_[i].getPos() - center).applyEach(Abs<F64>()).getMax();
                hl_max_loc = (hl_max_loc > hl_tmp) ? hl_max_loc : hl_tmp;
            }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const F64 hl_max_glb = Comm::getMaxValue(hl_max_loc);
            return hl_max_glb;
#else
            return hl_max_loc;
#endif
        }

// *******************************************************************        
// ************** This can be replaced with MT method. ***************
// *******************************************************************        
        inline F64 frand() {
//                return (double) rand() / ((double) RAND_MAX + 1.);
            //return genrand_res53();
            return MT::genrand_res53();
        }
// *******************************************************************        
// *******************************************************************        

        inline S32 getUniformDistributionFromArg1ToArg2(S32 arg1,
                                                        S32 arg2) {
            S32 random_number;
            
            random_number = (S32)((arg2 - arg1 + 1) * frand()) + arg1;
            
            return random_number;
        }

        /* AT_DEBU
        void getSampleParticle(S32 & number_of_sample_particle,
                               F32vec pos_sample[],
                               const F32 weight=1.0) {
        */

        S32 getNumberOfSampleParticleLocal(const F32 weight){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            S64 nloc = (S64)ptcl_.size();
            const auto weight_all = Comm::getSum(weight);
            S32 number_of_sample_particle = (weight * n_smp_ptcl_tot_) / weight_all;
#if 0
            // modified to limit # of sample particles by M.I. 
            const F32 coef_limitter = 0.2;
            S64 nglb = Comm::getSum( (S64)nloc );
            S32 number_of_sample_particle_limit = ((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter)); // lower limit
            number_of_sample_particle = (number_of_sample_particle > number_of_sample_particle_limit) ? number_of_sample_particle : number_of_sample_particle_limit;
#endif
            number_of_sample_particle = (number_of_sample_particle < nloc) ? number_of_sample_particle : nloc;
            
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            std::cout<<"weight="<<weight<<" weight_all="<<weight_all<<std::endl;
            std::cout<<"n_smp_ptcl_tot_="<<n_smp_ptcl_tot_<<std::endl;
            std::cout<<"(weight * n_smp_ptcl_tot_) / weight_all="<<(weight * n_smp_ptcl_tot_) / weight_all<<std::endl;
            std::cout<<"((S64)nloc * n_smp_ptcl_tot_)="<< ((S64)nloc * n_smp_ptcl_tot_)<<std::endl;
            std::cout<<"number_of_sample_particle(final)="<<number_of_sample_particle<<std::endl;
#endif
            return number_of_sample_particle;
#else
            return 0;
#endif

        }

        void getSampleParticle2(const S32 & number_of_sample_particle,
                                F64vec pos_sample[]){
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            RandomSampling(&ptcl_[0], pos_sample, ptcl_.size(), number_of_sample_particle, [](const Tptcl & src, F64vec & dst){dst = src.getPos();});
            return;
#endif
        }

        
        void getSampleParticle(S32 & number_of_sample_particle,
                               F64vec pos_sample[],
                               const F32 weight) {

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            S64 nloc = (S64)ptcl_.size();
            F32 weight_all = Comm::getSum(weight);
            number_of_sample_particle = (weight * n_smp_ptcl_tot_) / weight_all;
#if 0
// modified to limit # of sample particles by M.I. 
            const F32 coef_limitter = 0.2;
            S64 nglb = Comm::getSum( (S64)nloc );
            S32 number_of_sample_particle_limit = ((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter)); // lower limit
            number_of_sample_particle = (number_of_sample_particle > number_of_sample_particle_limit) ? number_of_sample_particle : number_of_sample_particle_limit;
#endif
            number_of_sample_particle = (number_of_sample_particle < nloc) ? number_of_sample_particle : nloc;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            std::cout<<"weight="<<weight<<" weight_all="<<weight_all<<std::endl;
            std::cout<<"n_smp_ptcl_tot_="<<n_smp_ptcl_tot_<<std::endl;
            std::cout<<"(weight * n_smp_ptcl_tot_) / weight_all="<<(weight * n_smp_ptcl_tot_) / weight_all<<std::endl;
            std::cout<<"((S64)nloc * n_smp_ptcl_tot_)="<< ((S64)nloc * n_smp_ptcl_tot_)<<std::endl;
            std::cout<<"number_of_sample_particle(final)="<<number_of_sample_particle<<std::endl;
#endif
            S32 *record = new S32[number_of_sample_particle];
            for(S32 i = 0; i < number_of_sample_particle; i++) {
                S32 j = getUniformDistributionFromArg1ToArg2(i, nloc-1);
                Tptcl hold = ptcl_[j];
                ptcl_[j]   = ptcl_[i];
                ptcl_[i]   = hold;
                record[i]  = j;
            }
            for(S32 i = 0; i < number_of_sample_particle; i++) {
                pos_sample[i] = ptcl_[i].getPos();
            }
            for(S32 i = number_of_sample_particle - 1; i >= 0; i--) {
                S32 j = record[i];
                Tptcl hold = ptcl_[j];
                ptcl_[j]   = ptcl_[i];
                ptcl_[i]   = hold;
            }
            delete [] record;
            return;
#endif
        }

        // new version (with switch)
        // for search, use tree with level 3(x, y, z) 
        // must be consistend with geometry of domains.
        template<class Tdinfo>
        void exchangeParticle(Tdinfo & dinfo,
                              const bool flag_serialize=false) {
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            F64 time_offset = GetWtime();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const S32 nloc  = ptcl_.size();
            const S32 rank  = Comm::getRank();
            const S32 nproc = Comm::getNumberOfProc();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort thisdomain = dinfo.getPosDomain(rank);

            ReallocatableArray<S32> nsend(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<S32> nsend_disp(nproc+1, nproc+1, MemoryAllocMode::Pool);
            ReallocatableArray<S32> nrecv(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<S32> nrecv_disp(nproc+1, nproc+1, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Request> req_send(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Request> req_recv(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Status> status(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<Tptcl> ptcl_send(0, 0, MemoryAllocMode::Pool);
            
            for(S32 i = 0; i < nproc; i++) {
                nsend[i] = nsend_disp[i] = nrecv[i] = nrecv_disp[i] = 0;
            }
            nsend_disp[nproc] = nrecv_disp[nproc] = 0;
            F64 time_offset_inner = GetWtime();
            for(S32 ip = 0; ip < nloc; ip++) {
                if( dinfo.getPosRootDomain().notOverlapped(ptcl_[ip].getPos()) ){
                    PARTICLE_SIMULATOR_PRINT_ERROR("A particle is out of root domain");
                    std::cerr<<"position of the particle="<<ptcl_[ip].getPos()<<std::endl;
                    std::cerr<<"position of the root domain="<<dinfo.getPosRootDomain()<<std::endl;
                    Abort(-1);
                }
                if(!determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    nsend[srank]++;
                }
            }

            nsend_disp[0] = 0;
            for(S32 i = 0; i < nproc; i++) {
                nsend_disp[i+1] += nsend_disp[i] + nsend[i];
            }
            ptcl_send.resizeNoInitialize( nsend_disp[nproc] );
            // ****************************************************
            // *** align send particles on ptcl_send_ *************
            for(S32 i = 0; i < nproc; i++) nsend[i] = 0;
            S32 iloc = 0;
            for(S32 ip = 0; ip < nloc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    iloc++;
                } else {
                    S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    S32 jloc = nsend[srank] + nsend_disp[srank];
                    ptcl_send[jloc] = ptcl_[ip];
                    nsend[srank]++;
                }
            }
            ptcl_.resizeNoInitialize(iloc);
            //time_profile_.exchange_particle__find_particle = GetWtime() - time_offset_inner;
            time_profile_.exchange_particle__find_particle += GetWtime() - time_offset_inner;

            // ****************************************************
            // *** receive the number of receive particles ********
            time_offset_inner = GetWtime();
            Comm::allToAll(nsend.getPointer(), 1, nrecv.getPointer());

            nrecv_disp[0] = 0;
            for(S32 i=0; i<nproc; i++){
                nrecv_disp[i+1] = nrecv_disp[i] + nrecv[i];
            }
            //ptcl_recv_.resizeNoInitialize( nrecv_disp[nproc] );
            ptcl_.resizeNoInitialize( nrecv_disp[nproc] + iloc);
            const S32 n_proc_comm_limit = 500;
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0) n_proc_send++;
                if(nrecv[i] > 0) n_proc_recv++;
            }
            bool flag_one_to_one_comm = ( n_proc_send < n_proc_comm_limit && n_proc_recv < n_proc_comm_limit) ? true : false;

            if( Comm::synchronizeConditionalBranchAND(flag_one_to_one_comm) ){
                n_proc_send = n_proc_recv = 0;
                for(S32 ib = 1; ib < nproc; ib++) {
                    S32 idsend = (ib + rank) % nproc;
                    if(nsend[idsend] > 0) {
                        S32 adrsend = nsend_disp[idsend];
                        S32 tagsend = (rank < idsend) ? rank : idsend;
                        MPI_Isend(ptcl_send.getPointer(adrsend), nsend[idsend], GetDataType<Tptcl>(), idsend, tagsend, MPI_COMM_WORLD, &req_send[n_proc_send++]);
                    }
                    S32 idrecv = (nproc + rank - ib) % nproc;
                    if(nrecv[idrecv] > 0) {
                        //S32 adrrecv = nrecv_disp[idrecv];
                        S32 adrrecv = nrecv_disp[idrecv] + iloc;
                        S32 tagrecv = (rank < idrecv) ? rank : idrecv;
                        //MPI_Irecv(ptcl_recv_.getPointer(adrrecv), nrecv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv, MPI_COMM_WORLD, &req_recv[n_proc_recv++]);
                        MPI_Irecv(ptcl_.getPointer(adrrecv), nrecv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv, MPI_COMM_WORLD, &req_recv[n_proc_recv++]);                        
                    }
                }
                MPI_Waitall(n_proc_send, req_send.getPointer(), status.getPointer());
                MPI_Waitall(n_proc_recv, req_recv.getPointer(), status.getPointer());
            }
            else{
                //Comm::allToAllV(ptcl_send_.getPointer(), nsend, nsend_disp, ptcl_recv_.getPointer(), nrecv, nrecv_disp);
                Comm::allToAllV(ptcl_send.getPointer(), nsend.getPointer(), nsend_disp.getPointer(), ptcl_.getPointer(iloc), nrecv.getPointer(), nrecv_disp.getPointer());
            }

            //const S32 nrecv_tot = nrecv_disp[nproc];
            //ptcl_.reserveEmptyAreaAtLeast( nrecv_tot );
            //for(S32 ip = 0; ip < nrecv_tot; ip++) {
            //    ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            //}
            // **************************************************** 
            n_ptcl_send_ += nsend_disp[nproc];
            n_ptcl_recv_ += nrecv_disp[nproc];
            time_profile_.exchange_particle__exchange_particle += GetWtime() - time_offset_inner;

            //delete [] nsend;
            //delete [] nsend_disp;
            //delete [] nrecv;
            //delete [] nrecv_disp;
            //delete [] req_send;
            //delete [] req_recv;
            //delete [] status;

            nsend.freeMem();
            nsend_disp.freeMem();
            nrecv.freeMem();
            nrecv_disp.freeMem();
            req_send.freeMem();
            req_recv.freeMem();
            status.freeMem();
            ptcl_send.freeMem();
#else
            n_ptcl_send_ = 0;
            n_ptcl_recv_ = 0;
#endif
            //const auto wtime_end = GetWtime();
            time_profile_.exchange_particle += GetWtime() - time_offset;
            //std::cerr<<"wtime_end-time_offset= "<<wtime_end-time_offset
            //         <<std::endl;
        }
        
        template<class Tcomp>
        void sortParticle(Tcomp comp){
            const S32 n = ptcl_.size();
            std::sort(ptcl_.getPointer(), ptcl_.getPointer(n), comp);
        }

        template<class Tcomp>
        void sortParticleStable(Tcomp comp){
            const S32 n = ptcl_.size();
            std::stable_sort(ptcl_.getPointer(), ptcl_.getPointer(n), comp);
        }

        // for DEBUG functions
        template<class Treal, class Tvec>
        void calcCMDirect(Treal & mass_cm, Tvec & pos_cm){
            mass_cm = 0.0;
            pos_cm = 0.0;
            const S32 n_ptcl = ptcl_.size();
            for(S32 i=0; i<n_ptcl; i++){
                mass_cm += ptcl_[i].mass;
                pos_cm += ptcl_[i].mass * ptcl_[i].pos;
            }
            pos_cm /= mass_cm;
        }

        bool checkExchangeParticleAllParticleInside(DomainInfo & dinfo);
        bool checkExchangeParticleSumOfNumberOfParticle(DomainInfo & dinfo,
                                                        S32 ntot_init);

        void adjustPositionIntoRootDomain(const DomainInfo & dinfo){
            const F64ort pos_root = dinfo.getPosRootDomain();
            const F64vec len_root = pos_root.getFullLength();
            const S32 n = ptcl_.size();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n; i++){
#if __cplusplus < 201103L
                F64vec pos_new = ptcl_[i].getPos() ;
#else
                auto pos_new = ptcl_[i].getPos() ;
#endif
                //if( pos_root.notOverlapped(pos_new) ){
                    while(pos_new.x < pos_root.low_.x){
                        pos_new.x += len_root.x;
                    }
                    while(pos_new.x >= pos_root.high_.x){
                        pos_new.x -= len_root.x;
                    }
                    if(pos_new.x == pos_root.high_.x){
                        pos_new.x = pos_root.low_.x;
                    }
                    while(pos_new.y < pos_root.low_.y){
                        pos_new.y += len_root.y;
                    }
                    while(pos_new.y >= pos_root.high_.y){
                        pos_new.y -= len_root.y;
                    }
                    if(pos_new.y == pos_root.high_.y){
                        pos_new.y = pos_root.low_.y;
                    }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                    while(pos_new.z < pos_root.low_.z){
                        pos_new.z += len_root.z;
                    }
                    while(pos_new.z >= pos_root.high_.z){
                        pos_new.z -= len_root.z;
                    }
                    if(pos_new.z == pos_root.high_.z){
                        pos_new.z = pos_root.low_.z;
                    }
#endif
                    //}
                ptcl_[i].setPos(pos_new);
            }
        }

        size_t getMemSizeUsed() const {
            return ptcl_.getMemSize();
        }
        size_t getUsedMemorySize() const {
            return getMemSizeUsed();
        }
        CountT getNumberOfParticleSendLocal() const { return (CountT)n_ptcl_send_; }
        CountT getNumberOfParticleRecvLocal() const { return (CountT)n_ptcl_recv_; }
        CountT getNumberOfParticleSendGlobal() const { return Comm::getSum((CountT)n_ptcl_send_); }
        CountT getNumberOfParticleRecvGlobal() const { return Comm::getSum((CountT)n_ptcl_recv_); }
        CountT getNumberOfProcExchangeLocal()  const { return (CountT)n_proc_exch_; }
        CountT getNumberOfProcExchangeGlobal() const { return Comm::getSum((CountT)n_proc_exch_); }
        void clearCounterAll(){
            n_ptcl_send_ = n_ptcl_recv_ = 0;
            n_proc_exch_ = 0;
            time_profile_.clear();
        }

        //////////
        // add and remove particles
        void addOneParticle(const Tptcl & fp){
            ptcl_.push_back(fp);
        }
        void removeParticle(const S32 * idx, const S32 n_remove){
            idx_remove_ptcl_.resizeNoInitialize(n_remove);
            for(S32 i=0; i<n_remove; i++){
                idx_remove_ptcl_[i] = idx[i];
            }
            std::sort(idx_remove_ptcl_.getPointer(), idx_remove_ptcl_.getPointer(n_remove));
            S32 * ptr_end = std::unique(idx_remove_ptcl_.getPointer(), idx_remove_ptcl_.getPointer(n_remove));
            const S32 n_remove_tmp = ptr_end - idx_remove_ptcl_.getPointer();
            const S32 n_prev = ptcl_.size();
            S32 i_loc = n_prev-1;
            for(S32 i=n_remove_tmp-1; i>=0; i--){
                std::swap(ptcl_[idx_remove_ptcl_[i]], ptcl_[i_loc]);
                i_loc--;
            }
            ptcl_.resizeNoInitialize(i_loc+1);
            /*
            // original
            const S32 n_prev = ptcl_.size();
            flag_remove.resizeNoInitialize(n_prev);
            for(S32 i=0; i<n_prev; i++){
                flag_remove[i] = false;
            }
            for(S32 i=0; i<n_remove; i++){
                S32 idx_tmp = idx[i];
                flag_remove[idx_tmp] = true;
            }
            S32 i_loc = n_prev-1;
            for(S32 i=n_prev-1; i>=0; i--){
                if(flag_remove[i] == true){
                    std::swap(ptcl_[i], ptcl_[i_loc]);
                    i_loc--;
                }
            }
            ptcl_.resizeNoInitialize(i_loc+1);
            */
        }

        void dumpMemSizeUsed(std::ostream & fout){
            F64 sum_loc,sum_max;
            sum_loc = (double)(ptcl_.getMemSize()
                              +idx_remove_ptcl_.getMemSize()) * 1e-9;
                               //+ptcl_send_.getMemSize()
                               //+ptcl_recv_.getMemSize()) * 1e-9;
            if (Comm::getRank() == 0) {
                fout<<"ptcl_.getMemSize()= "<<ptcl_.getMemSize()<<std::endl;
                fout<<"    ptcl_.size()= "<<ptcl_.size()<<std::endl;
                fout<<"    sizeof(Tptcl)= "<<sizeof(Tptcl)<<std::endl;
                fout<<"idx_remove_ptcl_.getMemSize()= "<<idx_remove_ptcl_.getMemSize()<<std::endl;
                fout<<"sum[GB]= " << sum_loc << std::endl;
            }
            sum_max = Comm::getMaxValue(sum_loc);
            if (Comm::getRank() == 0) {
               fout << "sum[GB](psys,max.) = " << sum_max << std::endl;
            }
        }

        S32 searchWhichDomainParticleGoToOpen(const F64vec & pos,
                                              const S32 n_domain [],
                                              const F64ort domain [],
                                              const S32 my_rank){
            auto rank = my_rank;
            auto shift = n_domain[1]*n_domain[2];
            if(pos.x < domain[rank].low_.x){
                do{
                    rank -= shift;
                }while(pos.x < domain[rank].low_.x);
            }
            else if(pos.x >= domain[rank].high_.x) {
                do{
                    rank += shift;
                }while(pos.x >= domain[rank].high_.x);
            }
            shift = n_domain[2];
            if(pos.y < domain[rank].low_.y){
                do{
                    rank -= shift;
                }while(pos.y < domain[rank].low_.y);
            }
            else if(pos.y >= domain[rank].high_.y) {
                do{
                    rank += shift;
                }while(pos.y >= domain[rank].high_.y);
            }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            shift = 1;
            if(pos.z < domain[rank].low_.z){
                do{
                    rank -= shift;
                }while(pos.z < domain[rank].low_.z);
            }
            else if(pos.z >= domain[rank].high_.z) {
                do{
                    rank += shift;
                }while(pos.z >= domain[rank].high_.z);
            }
#endif
            return rank;
        }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        S32 searchWhichDomainParticleGoTo(const F64vec & pos,
                                          const S32 n_domain [],
                                          const F64ort domain [],
                                          const S32 my_rank,
                                          const F64vec len_peri,
                                          const bool pa[],
                                          S32vec & dnp){
            const auto rank_x_org = my_rank / (n_domain[1]*n_domain[2]);
            const auto rank_y_org = (my_rank / n_domain[2]) % n_domain[1];
            const auto rank_z_org = my_rank % n_domain[2];
            auto rank_x = rank_x_org;
            auto rank_y = rank_y_org;
            auto rank_z = rank_z_org;
            auto rank = my_rank;
            // x direction
            dnp[0] = 0;
            auto shift = n_domain[1]*n_domain[2];
            if(pa[0]){
                if(domain[rank].low_.x > pos.x || pos.x >= domain[rank].high_.x){
                    // out of self domain
                    auto x_tmp = (pos.x < domain[rank].low_.x) ?
                        (domain[rank].low_.x-pos.x < pos.x+len_peri.x-domain[rank].high_.x ? pos.x : pos.x+len_peri.x) :
                        (pos.x-domain[rank].high_.x < domain[rank].low_.x-(pos.x-len_peri.x) ? pos.x : pos.x-len_peri.x);
                    if(x_tmp < domain[rank].low_.x){
                        if(x_tmp != pos.x)
                            rank = (n_domain[0]-1)*n_domain[1]*n_domain[2] + rank_y*n_domain[2] + rank_z;
                        else
                            rank -= shift;
                        while(pos.x < domain[rank].low_.x){
                            rank -= shift;
                        };
                    }
                    else if(x_tmp >= domain[rank].high_.x){
                        if(x_tmp != pos.x) rank = rank_y*n_domain[2] + rank_z;
                        else rank += shift;
                        while(pos.x >= domain[rank].high_.x){
                            rank += shift;
                        };
                    }
                }
                rank_x = rank / (n_domain[1]*n_domain[2]);
                const auto dx = std::abs(rank_x-rank_x_org);
                dnp[0] = std::min(dx, n_domain[0]-dx);
            }
            else{
                if(pos.x < domain[rank].low_.x){
                    do{
                        rank -= shift;
                    }while(pos.x < domain[rank].low_.x);
                }
                else if(pos.x >= domain[rank].high_.x) {
                    do{
                        rank += shift;
                    }while(pos.x >= domain[rank].high_.x);
                }
                rank_x = rank / (n_domain[1]*n_domain[2]);
                dnp[0] = std::abs(rank_x-rank_x_org);
            }
            //rank_x = rank / (n_domain[1]*n_domain[2]);
            //dnp[0] = std::abs(rank_x-rank_x_org);

            // y direction
            dnp[1] = 0;
            shift = n_domain[2];
            if(pa[1]){
                if(domain[rank].low_.y > pos.y || pos.y >= domain[rank].high_.y){
                    auto y_tmp = (pos.y < domain[rank].low_.y) ?
                        (domain[rank].low_.y-pos.y < pos.y+len_peri.y-domain[rank].high_.y ? pos.y : pos.y+len_peri.y) :
                        (pos.y-domain[rank].high_.y < domain[rank].low_.y-(pos.y-len_peri.y) ? pos.y : pos.y-len_peri.y);
                    if(y_tmp < domain[rank].low_.y){
                        if(y_tmp != pos.y) rank = rank_x*n_domain[1]*n_domain[2] + (n_domain[1]-1)*n_domain[2] + rank_z;
                        else rank -= shift;
                        while(pos.y < domain[rank].low_.y){
                            rank -= shift;
                        };
                    }
                    else if(y_tmp >= domain[rank].high_.y){
                        if(y_tmp != pos.y) rank = rank_x*n_domain[1]*n_domain[2] + rank_z;
                        else rank += shift;
                        while(pos.y >= domain[rank].high_.y){
                            rank += shift;
                        };
                    }
                }
                rank_y = (rank / n_domain[2]) % n_domain[1];
                const auto dy = std::abs(rank_y-rank_y_org);
                dnp[1] = std::min(dy, n_domain[1]-dy);
            }
            else{
                if(pos.y < domain[rank].low_.y){
                    do{
                        rank -= shift;
                    }while(pos.y < domain[rank].low_.y);
                }
                else if(pos.y >= domain[rank].high_.y) {
                    do{
                        rank += shift;
                    }while(pos.y >= domain[rank].high_.y);
                }
                rank_y = (rank / n_domain[2]) % n_domain[1];
                dnp[1] = std::abs(rank_y-rank_y_org);                
            }
            //rank_y = (rank / n_domain[2]) % n_domain[1];
            //dnp[1] = std::abs(rank_y-rank_y_org);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            dnp[2] = 0;
            shift = 1;
            if(pa[2]){
                if(domain[rank].low_.z > pos.z || pos.z >= domain[rank].high_.z){
                    auto z_tmp = (pos.z < domain[rank].low_.z) ?
                        (domain[rank].low_.z-pos.z < pos.z+len_peri.z-domain[rank].high_.z ? pos.z : pos.z+len_peri.z) :
                        (pos.z-domain[rank].high_.z < domain[rank].low_.z-(pos.z-len_peri.z) ? pos.z : pos.z-len_peri.z);
                    if(z_tmp < domain[rank].low_.z){
                        if(z_tmp != pos.z) rank = rank_x*n_domain[1]*n_domain[2] + rank_y*n_domain[2] + n_domain[2]-1;
                        else rank -= shift;
                        while(pos.z < domain[rank].low_.z){
                            rank -= shift;
                        };
                    }
                    else if(z_tmp >= domain[rank].high_.z){
                        if(z_tmp != pos.z) rank = rank_x*n_domain[1]*n_domain[2] + rank_y*n_domain[2];
                        else rank += shift;
                        while(pos.z >= domain[rank].high_.z){
                            rank += shift;
                        };
                    }
                }
                rank_z = rank % n_domain[2];
                const auto dz = std::abs(rank_z-rank_z_org);
                dnp[2] = std::min(dz, n_domain[2]-dz);
            }
            else{
                if(pos.z < domain[rank].low_.z){
                    do{
                        rank -= shift;
                    }while(pos.z < domain[rank].low_.z);
                }
                else if(pos.z >= domain[rank].high_.z) {
                    do{
                        rank += shift;
                    }while(pos.z >= domain[rank].high_.z);
                }
                rank_z = rank % n_domain[2];
                dnp[2] = std::abs(rank_z-rank_z_org);
            }
            //rank_z = rank % n_domain[2];
            //dnp[2] = std::abs(rank_z-rank_z_org);
#endif
            return rank;
        }

        void exchangeNumberOfPtcl(ReallocatableArray<S32> & n_send, // const
                                  ReallocatableArray<S32> & n_recv,
                                  const S32 n_domain[],
                                  const S32vec & dnp_max){
            const auto my_rank = Comm::getRank();
            const auto n_proc  = Comm::getNumberOfProc();
            ReallocatableArray<MPI_Request> req_recv(n_proc, n_proc, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Status> status(n_proc, n_proc, MemoryAllocMode::Pool);
            const auto rank_x_org = my_rank / (n_domain[1]*n_domain[2]);
            const auto rank_y_org = (my_rank / n_domain[2]) % n_domain[1];
            const auto rank_z_org = my_rank % n_domain[2];
            S32 n_cnt = 0;
            S32 dnp_x_head = -dnp_max[0];
            S32 dnp_x_tail = dnp_max[0];
            if(dnp_max[0]*2+1 > n_domain[0]){
                dnp_x_head = -n_domain[0] / 2;
                dnp_x_tail = (n_domain[0]%2==0) ? (n_domain[0]/2 - 1) : n_domain[0]/2;
            }
            S32 dnp_y_head = -dnp_max[1];
            S32 dnp_y_tail = dnp_max[1];
            if(dnp_max[1]*2+1 > n_domain[1]){
                dnp_y_head = -n_domain[1] / 2;
                dnp_y_tail = (n_domain[1]%2==0) ? (n_domain[1]/2 - 1) : n_domain[1]/2;
            }
            S32 dnp_z_head = -dnp_max[2];
            S32 dnp_z_tail = dnp_max[2];
            if(dnp_max[2]*2+1 > n_domain[2]){
                dnp_z_head = -n_domain[2] / 2;
                dnp_z_tail = (n_domain[2]%2==0) ? (n_domain[2]/2 - 1) : n_domain[2]/2;
            }            
            for(S32 i=dnp_x_head; i<=dnp_x_tail; i++){
                auto rank_x = (rank_x_org+n_domain[0]+i) % n_domain[0];
                for(S32 j=dnp_y_head; j<=dnp_y_tail; j++){
                    auto rank_y = (rank_y_org+n_domain[1]+j) % n_domain[1];
                    for(S32 k=dnp_z_head; k<=dnp_z_tail; k++){
                        auto rank_z = (rank_z_org+n_domain[2]+k) % n_domain[2];
                        const auto rank = rank_x*(n_domain[1]*n_domain[2]) + rank_y*n_domain[2] + rank_z;
                        MPI_Irecv(&n_recv[rank], 1, GetDataType<S32>(), rank, 0, MPI_COMM_WORLD, &req_recv[n_cnt]);
                        n_cnt++;
                    }
                }
            }
            n_cnt = 0;
            for(S32 i=dnp_x_head; i<=dnp_x_tail; i++){
                auto rank_x = (rank_x_org+n_domain[0]+i) % n_domain[0];
                for(S32 j=dnp_y_head; j<=dnp_y_tail; j++){
                    auto rank_y = (rank_y_org+n_domain[1]+j) % n_domain[1];
                    for(S32 k=dnp_z_head; k<=dnp_z_tail; k++){
                        auto rank_z = (rank_z_org+n_domain[2]+k) % n_domain[2];
                        const auto rank = rank_x*(n_domain[1]*n_domain[2]) + rank_y*n_domain[2] + rank_z;
                        MPI_Send(&n_send[rank], 1, GetDataType<S32>(), rank, 0, MPI_COMM_WORLD);
                        n_cnt++;
                    }
                }
            }
            MPI_Waitall(n_cnt, &req_recv[0], &status[0]);
        }
#endif
        
        template<class Tdinfo>
        void exchangeParticle2(Tdinfo & dinfo,
                               const bool flag_serialize=false) {
            if (flag_serialize == true) {
                PARTICLE_SIMULATOR_PRINT_ERROR("serialization is not yet supported.");
                Abort(-1);
            }
            F64 time_offset = GetWtime();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const S32 nloc  = ptcl_.size();
            const S32 rank  = Comm::getRank();
            const S32 nproc = Comm::getNumberOfProc();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            ReallocatableArray<S32> nsend(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<S32> nsend_disp(nproc+1, nproc+1, MemoryAllocMode::Pool);
            ReallocatableArray<S32> nrecv(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<S32> nrecv_disp(nproc+1, nproc+1, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Request> req_send(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Request> req_recv(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<MPI_Status> status(nproc, nproc, MemoryAllocMode::Pool);
            ReallocatableArray<Tptcl> ptcl_send(0, 0, MemoryAllocMode::Pool);
            ReallocatableArray<S32> rank_send(nloc, nloc, MemoryAllocMode::Pool);
            for(S32 i = 0; i < nproc; i++) {
                nsend[i] = nsend_disp[i] = nrecv[i] = nrecv_disp[i] = 0;
            }
            nsend_disp[nproc] = nrecv_disp[nproc] = 0;
            F64 time_offset_inner = GetWtime();
            F64ort pos_root_domain = dinfo.getPosRootDomain();
            F64vec len_peri = pos_root_domain.high_ - pos_root_domain.low_;
            bool pa[DIMENSION_LIMIT];
            dinfo.getPeriodicAxis(pa);
            S32vec dnp_max_loc = 0;
            S32 iloc = nloc-1;
            for(S32 ip = nloc-1; ip >= 0; ip--) {
                if( dinfo.getPosRootDomain().notOverlapped(ptcl_[ip].getPos()) ){
                    PARTICLE_SIMULATOR_PRINT_ERROR("A particle is out of root domain");
                    std::cerr<<"position of the particle="<<ptcl_[ip].getPos()<<std::endl;
                    std::cerr<<"position of the root domain="<<dinfo.getPosRootDomain()<<std::endl;
                    Abort(-1);
                }
                S32vec dnp;
                S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain, Comm::getRank(), len_peri, pa, dnp);
                rank_send[ip] = rank;
                if(srank == rank) continue;
                std::swap(ptcl_[ip], ptcl_[iloc]);
                rank_send[iloc] = srank;
                iloc--;
                dnp_max_loc.x = std::max(dnp_max_loc.x, dnp.x);
                dnp_max_loc.y = std::max(dnp_max_loc.y, dnp.y);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                dnp_max_loc.z = std::max(dnp_max_loc.z, dnp.z);
#endif
                nsend[srank]++;
            }
            const auto n_remain = iloc+1;
            S32vec dnp_max_glb = Comm::getMaxValue(dnp_max_loc);
            nsend_disp[0] = 0;
            for(S32 i = 0; i < nproc; i++) {
                nsend_disp[i+1] += nsend_disp[i] + nsend[i];
            }
            ptcl_send.resizeNoInitialize( nsend_disp[nproc] );
            // ****************************************************
            // *** align send particles on ptcl_send_ *************
            for(S32 i = 0; i < nproc; i++) nsend[i] = 0;
            for(S32 ip = n_remain; ip < nloc; ip++) {
                const auto srank = rank_send[ip];
                S32 jloc = nsend[srank] + nsend_disp[srank];
                ptcl_send[jloc] = ptcl_[ip];
                nsend[srank]++;
            }
            ptcl_.resizeNoInitialize(n_remain);
            time_profile_.exchange_particle__find_particle += GetWtime() - time_offset_inner;
            
            // ****************************************************
            // *** receive the number of receive particles ********
            time_offset_inner = GetWtime();
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
            if(2*dnp_max_glb[0]*2*dnp_max_glb[1] < nproc/4){
#else
            if(2*dnp_max_glb[0]*2*dnp_max_glb[1]*2*dnp_max_glb[2] < nproc/8){
#endif
                exchangeNumberOfPtcl(nsend, nrecv, n_domain, dnp_max_glb);
            }
            else{
                //exchangeNumberOfPtcl(nsend, nrecv, n_domain, dnp_max_glb);
                Comm::allToAll(nsend.getPointer(), 1, nrecv.getPointer());
            }

            nrecv_disp[0] = 0;
            for(S32 i=0; i<nproc; i++){
                nrecv_disp[i+1] = nrecv_disp[i] + nrecv[i];
            }
            ptcl_.resizeNoInitialize( nrecv_disp[nproc] + n_remain);

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
            if(2*dnp_max_glb[0]*2*dnp_max_glb[1] < nproc/4){
#else
            if(2*dnp_max_glb[0]*2*dnp_max_glb[1]*2*dnp_max_glb[2] < nproc/8){
#endif
                S32 n_cnt = 0;
                for(S32 ib = 1; ib < nproc; ib++) {
                    S32 idrecv = (nproc + rank - ib) % nproc;
                    if(nrecv[idrecv] > 0) {
                        S32 adrrecv = nrecv_disp[idrecv] + n_remain;
                        S32 tagrecv = (rank < idrecv) ? rank : idrecv;
                        MPI_Irecv(ptcl_.getPointer(adrrecv), nrecv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv, MPI_COMM_WORLD, &req_recv[n_cnt]);
                        n_cnt++;
                    }
                }
                for(S32 ib = 1; ib < nproc; ib++) {
                    S32 idsend = (ib + rank) % nproc;
                    if(nsend[idsend] > 0) {
                        S32 adrsend = nsend_disp[idsend];
                        S32 tagsend = (rank < idsend) ? rank : idsend;
                        MPI_Send(ptcl_send.getPointer(adrsend), nsend[idsend], GetDataType<Tptcl>(), idsend, tagsend, MPI_COMM_WORLD);
                    }
                }
                MPI_Waitall(n_cnt, &req_recv[0], &status[0]);
            }
            else{
                Comm::allToAllV(ptcl_send.getPointer(), nsend.getPointer(), nsend_disp.getPointer(), ptcl_.getPointer(n_remain), nrecv.getPointer(), nrecv_disp.getPointer());
            }
            n_ptcl_send_ += nsend_disp[nproc];
            n_ptcl_recv_ += nrecv_disp[nproc];
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
            n_proc_exch_ += 2*dnp_max_glb[0]*2*dnp_max_glb[1];
#else
            n_proc_exch_ += 2*dnp_max_glb[0]*2*dnp_max_glb[1]*2*dnp_max_glb[2];
#endif            
            time_profile_.exchange_particle__exchange_particle += GetWtime() - time_offset_inner;
            nsend.freeMem();
            nsend_disp.freeMem();
            nrecv.freeMem();
            nrecv_disp.freeMem();
            req_send.freeMem();
            req_recv.freeMem();
            status.freeMem();
            ptcl_send.freeMem();
#else
            n_ptcl_send_ = 0;
            n_ptcl_recv_ = 0;
#endif
            time_profile_.exchange_particle += GetWtime() - time_offset;
        }
    };
}
