#pragma once

#include<cassert>
#include<fstream>
#include<time.h>

#include "MT.hpp"
#include"ps_defs.hpp"

#ifdef SUNWAY
extern "C"{
    #include <athread.h>
    #include <cpe_func.h>
    void SLAVE_FUN(CheckIsInDomain)(void *);
    void SLAVE_FUN(FindDomainParticleGoTo)(void *);
    void SLAVE_FUN(SetParticleToSendBuffer)(void *);
    void SLAVE_FUN(GetCylCoord)(void *);
    void SLAVE_FUN(FindDestinationRank)(void *);
    void SLAVE_FUN(MakeSendBuffers)(void *);
}
#endif

namespace ParticleSimulator{
    long int DISP_WALK_EX_PTCL = 0;
    long int N_NOT_MOVE = 0;
    long int N_LOC_ORG = 0;
    double   DISP_WALK_EX_PTCL_AVE = 0;
    double HIT_RATIO = 0.0;

    
    template<class Tptcl>
    class ParticleSystem{
    private:
        CountT n_ptcl_send_;
        CountT n_ptcl_recv_;
        TimeProfile time_profile_;
        static const S32 n_smp_ave_ = 30;
        ReallocatableArray<Tptcl> ptcl_;
        ReallocatableArray<S32> idx_remove_ptcl_; // add 2016/09/14
        ReallocatableArray<Tptcl> ptcl_send_;
        ReallocatableArray<Tptcl> * ptcl_send_buf_;
        ReallocatableArray<Tptcl> ptcl_recv_;
        S32 n_smp_ptcl_tot_;
        bool first_call_by_initialize;
        bool first_call_by_setAverageTargetNumberOfSampleParticlePerProcess;
        bool first_call_by_DomainInfo_collect_sample_particle;
        inline bool determineWhetherParticleIsInDomain(const F64vec & pos,
                                                       const F64ort & domain) {
            bool ret = true;
	    ret *= (domain.low_.x <= pos.x) * (pos.x < domain.high_.x);
	    ret *= (domain.low_.y <= pos.y) * (pos.y < domain.high_.y);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	    ret *= (domain.low_.z <= pos.z) * (pos.z < domain.high_.z);
#endif
            return ret;
        }


        S32 searchWhichDomainParticleGoToPeriodicX(const F64vec & pos,
                                                   const S32 n_domain [],
                                                   const F64ort domain [],
                                                   const F64 peri_len_x,
                                                   const S32 previous_rank=0){
            S32 idomain = previous_rank;
            if(determineWhetherParticleIsInDomain(pos, domain[idomain])){
                return idomain;
            }
            S32 shift = n_domain[1] * n_domain[2];
            F64 box_cen_x = (domain[idomain].high_.x + domain[idomain].low_.x)*0.5;
            F64 box_len_x = (domain[idomain].high_.x - domain[idomain].low_.x)*0.5;
            F64 dx = pos.x-box_cen_x;
            if(dx > box_len_x){
                if(dx > 0.0 && dx > peri_len_x*0.5){
                    F64 rank_y = (idomain/shift) % n_domain[1];
                    idomain = rank_y*n_domain[2];
                }
                else if(dx < 0.0 && -dx > peri_len_x*0.5){
                    F64 rank_y = (idomain/shift) % n_domain[1];
                    idomain = rank_y*n_domain[2] + (n_domain[0]-1)*shift;
                }
            }
            // x direction
            if(domain[idomain].high_.x <= pos.x){
                do{
                    idomain += shift;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].high_.x <= pos.x);
            }
            else if(domain[idomain].low_.x > pos.x){
                do{
                    idomain -= shift;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].low_.x > pos.x);
            }
            
            // y direction
            shift = n_domain[2];
            if(domain[idomain].high_.y <= pos.y){
                do{
                    idomain += shift;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].high_.y <= pos.y);
            }
            else if(domain[idomain].low_.y > pos.y){
                do{
                    idomain -= shift;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].low_.y > pos.y);
            }
            
            // z direction
            if(domain[idomain].high_.z <= pos.z){
                do{
                    idomain++;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].high_.z <= pos.z);
            }
            else if(domain[idomain].low_.z > pos.z){
                do{
                    idomain--;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].low_.z > pos.z);
            }
            return idomain;
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
            while(domain[idomain].high_.x <= pos.x){
                idomain += nynz;
                DISP_WALK_EX_PTCL++;
            }
            const S32 nz   = n_domain[2];
            while(domain[idomain].high_.y <= pos.y){
                idomain += nz;
                DISP_WALK_EX_PTCL++;
            }
            while(domain[idomain].high_.z <= pos.z){
                idomain++;
                DISP_WALK_EX_PTCL++;
            }
            return idomain;
#endif
        }
        
        S32 searchWhichDomainParticleGoTo2(const F64vec & pos,
                                           const S32 n_domain [],
                                           const F64ort domain [],
                                           const S32 previous_rank=0){
            S32 idomain = previous_rank;
            S32 shift = n_domain[1] * n_domain[2];
            //S32 disp_x = 0;
            //S32 disp_y = 0;
            //S32 disp_z = 0;
            // x direction
            if(domain[idomain].high_.x <= pos.x){
                do{
                    idomain += shift;
                    DISP_WALK_EX_PTCL++;
                    //disp_x++;
                }while(domain[idomain].high_.x <= pos.x);
            }
            else if(domain[idomain].low_.x > pos.x){
                do{
                    idomain -= shift;
                    DISP_WALK_EX_PTCL++;
                    //disp_x++;
                }while(domain[idomain].low_.x > pos.x);
            }
            
            // y direction
            shift = n_domain[2];
            if(domain[idomain].high_.y <= pos.y){
                do{
                    idomain += shift;
                    DISP_WALK_EX_PTCL++;
                    //disp_y++;
                }while(domain[idomain].high_.y <= pos.y);
            }
            else if(domain[idomain].low_.y > pos.y){
                do{
                    idomain -= shift;
                    DISP_WALK_EX_PTCL++;
                    //disp_y++;
                }while(domain[idomain].low_.y > pos.y);
            }
            
            // z direction
            if(domain[idomain].high_.z <= pos.z){
                do{
                    idomain++;
                    DISP_WALK_EX_PTCL++;
                    //disp_z++;
                }while(domain[idomain].high_.z <= pos.z);
            }
            else if(domain[idomain].low_.z > pos.z){
                do{
                    idomain--;
                    DISP_WALK_EX_PTCL++;
                    //disp_z++;
                }while(domain[idomain].low_.z > pos.z);
            }
            /*
            S32 disp_sun = disp_x + disp_y + disp_z;
            //if(disp_sun > 5 && previous_rank != 0){
            if(disp_sun > 5){
                std::cerr<<"previous_rank= "<<previous_rank
                         <<" my_rank= "<<Comm::getRank()
                         <<" idomain= "<<idomain
                         <<std::endl;
                std::cerr<<"disp_x= "<<disp_x
                         <<" disp_y= "<<disp_y
                         <<" disp_z= "<<disp_z
                         <<" sum= "<<disp_x+disp_y+disp_z
                         <<std::endl;
                std::cerr<<"pos.x= "<<pos.x
                         <<" pos.y= "<<pos.y
                         <<" pos.z= "<<pos.z
                         <<std::endl;
            }
            */
            return idomain;
        }

        S32 searchWhichDomainParticleGoTo3(const F64vec & pos,
                                           const S32 n_domain [],
                                           const F64ort domain [],
                                           const S32 my_rank,
                                           const S32 previous_rank=0){
            S32 idomain = previous_rank;
            S32 shift = n_domain[1] * n_domain[2];
            U32 flag = 0;
            if(domain[idomain].high_.x <= pos.x){
                do{
                    idomain += shift;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].high_.x <= pos.x);
            }
            else if(domain[idomain].low_.x > pos.x){
                do{
                    idomain -= shift;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].low_.x > pos.x);
            }
            else{
                flag++;
            }
            
            // y direction
            shift = n_domain[2];
            if(domain[idomain].high_.y <= pos.y){
                do{
                    idomain += shift;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].high_.y <= pos.y);
            }
            else if(domain[idomain].low_.y > pos.y){
                do{
                    idomain -= shift;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].low_.y > pos.y);
            }
            else{
                flag++;
            }
            
            // z direction
            if(domain[idomain].high_.z <= pos.z){
                do{
                    idomain++;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].high_.z <= pos.z);
            }
            else if(domain[idomain].low_.z > pos.z){
                do{
                    idomain--;
                    DISP_WALK_EX_PTCL++;
                }while(domain[idomain].low_.z > pos.z);
            }
            else{
                flag++;
            }
            if(flag == 3) return my_rank;
            else return idomain;
        }
        

    public:
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
	
        void initialize() {
            assert(first_call_by_initialize);
            first_call_by_initialize = false;
	    n_smp_ptcl_tot_ = n_smp_ave_ * Comm::getNumberOfProc();
            n_ptcl_send_ = n_ptcl_recv_ = 0;
            ptcl_send_buf_ = new ReallocatableArray<Tptcl>[Comm::getNumberOfProc()];
            
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
                    //S32 n_ptcl_ = header->readAscii(fp);
                    S32 n_ptcl_ = (header->*pFuncHead)(fp);
                    while('\n' == getc(fp));
                    fseek(fp, -1, SEEK_CUR);
                    if(n_ptcl_ < 0){//User does NOT return # of ptcl
                        //count # of lines
                        n_ptcl_ = 0;
                        //KN
                        for(int c ; (c = getc(fp)) != EOF ; n_ptcl_ += '\n' == c ? 1 : 0){}
                        fclose(fp);
                        fp = fopen(filename, "r");
                        //header->readAscii(fp);
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
                    MPI::COMM_WORLD.Scatter(n_ptcl, 1, GetDataType<S32>(), &n_ptcl_, 1, GetDataType<S32>(), 0);
                    #endif
                    //allocate ptcl.
                    //First of all, Rank 0 reads its own particle.
                    this->createParticle(n_ptcl_ << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_);
                    for(int i = 0 ; i < n_ptcl_ ; ++ i){
                        //ptcl_[i].readAscii(fp);
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    //Read remaining data to buffer and send them to appropriate process.
                    #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    for(S32 rank = 1 ; rank < n_proc ; ++ rank){
                        Tptcl * buffer = new Tptcl[n_ptcl[rank]];
                        for(int i = 0 ; i < n_ptcl[rank] ; ++ i){
                            //buffer[i].readAscii(fp);
                            (buffer[i].*pFuncPtcl)(fp);
                        }
                        MPI::COMM_WORLD.Send(buffer, n_ptcl[rank], GetDataType<Tptcl>(), rank, 0);
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
                    MPI::COMM_WORLD.Scatter(n_ptcl, 1, GetDataType<S32>(), &n_ptcl_loc, 1, GetDataType<S32>(), 0);
                    delete [] n_ptcl;
                    //allocate ptcl.
                    this->createParticle(n_ptcl_loc << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_loc);
                    MPI::COMM_WORLD.Recv(ptcl_.getPointer(), ptcl_.size(), GetDataType<Tptcl>(), 0, 0);
                    #endif
                }
            }else{//Read from multiple file
                char input[256];
                sprintf(input, format, filename, Comm::getNumberOfProc(), Comm::getRank());
                FILE* fp = fopen(input, "r");
                if(fp == NULL){
                    PARTICLE_SIMULATOR_PRINT_ERROR("can not open input file");
                    std::cerr<<"filename: "<<input<<std::endl;
                    Abort(-1);
                }
                //S32 n_ptcl_ = header->readAscii(fp);
                S32 n_ptcl_ = (header->*pFuncHead)(fp);
                while('\n' == getc(fp));
                fseek(fp, -1, SEEK_CUR);
                if(n_ptcl_ >= 0){
                    //User returns # of ptcl.
                    this->createParticle(n_ptcl_ << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_);
                    for(S32 i = 0 ; i < n_ptcl_ ; ++ i){
                        //ptcl_[i].readAscii(fp);
                        (ptcl_[i].*pFuncPtcl)(fp);
                    }
                    fclose(fp);
                }else{//User does NOT return # of ptcl
                    //count # of lines
                    n_ptcl_ = 0;
                    for(int c ; (c = getc(fp)) != EOF ; n_ptcl_ += c == '\n' ? 1 : 0){}
                    fclose(fp);
                    //
                    FILE* fp = fopen(input, "r");
                    //header->readAscii(fp);
                    (header->*pFuncHead)(fp);
                    while('\n' == getc(fp));
		    fseek(fp, -1, SEEK_CUR);
                    this->createParticle(n_ptcl_ << 2);//Magic shift
                    ptcl_.resizeNoInitialize(n_ptcl_);
                    for(S32 i = 0 ; i < ptcl_.size() ; ++ i){
                        //ptcl_[i].readAscii(fp);
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
                MPI::COMM_WORLD.Allgather(&n_ptcl_, 1, GetDataType<S32>(), n_ptcl, 1, GetDataType<S32>());
                //set displacement
                S32 *n_ptcl_displs = new S32[n_proc+1];
                n_ptcl_displs[0] = 0;
                for(S32 i = 0 ; i < n_proc ; ++ i){
                    n_ptcl_displs[i+1] = n_ptcl_displs[i] + n_ptcl[i];
                }
                const S32 n_tot = n_ptcl_displs[n_proc];
                Tptcl *ptcl = new Tptcl[n_tot];
                //gather data
                MPI::COMM_WORLD.Gatherv(ptcl_.getPointer(), n_ptcl_, GetDataType<Tptcl>(), ptcl, n_ptcl, n_ptcl_displs, GetDataType<Tptcl>(), 0);
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
            writeParticleImpl<DummyHeader>(filename, format, NULL, pFunc, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
        }
        void writeParticleBinary(const char * const filename, void (Tptcl::*pFunc)(FILE*)const){
            writeParticleImpl<DummyHeader>(filename, NULL, NULL, pFunc, &Tptcl::writeBinary, &DummyHeader::writeBinary, "wb");
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
/*
            #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            //get # of process.
            const S32 n_proc = Comm::getNumberOfProc();
            //get # of ptcls in each process.
            S64 *n_ptcl = new S64[n_proc];
            //Gather # of particles.
            S64 n_ptcl_ = ptcl_.size();
            MPI::COMM_WORLD.Allgather(&n_ptcl_, 1, GetDataType<S64>(), n_ptcl, 1, GetDataType<S64>());
            //set displacement
            S32 *n_ptcl_displs = new S32[n_proc+1];
            n_ptcl_displs[0] = 0;
            for(S32 i = 0 ; i < n_proc ; ++ i){
                n_ptcl_displs[i+1] = n_ptcl_displs[i] + n_ptcl[i];
            }
            return n_ptcl_displs[n_proc];
            #else
            return ptcl_.size();
            #endif
*/
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
            F64 hl_max_glb;
            MPI::COMM_WORLD.Allreduce(&hl_max_loc, &hl_max_glb, 1, GetDataType<F64>(), MPI::MAX);
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
        void getSampleParticle(S32 & number_of_sample_particle,
                               F64vec pos_sample[],
                               const F32 weight) {

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL

            S64 nloc = (S64)ptcl_.size();

            F32 weight_all = Comm::getSum(weight);
            number_of_sample_particle = (S32)((weight * n_smp_ptcl_tot_) / weight_all);
#if 0
// modified to limit # of sample particles by M.I. 
            const F32 coef_limitter = 0.2;
            S64 nglb = Comm::getSum( (S64)nloc );
            S32 number_of_sample_particle_limit = ((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter)); // lower limit
            number_of_sample_particle = (number_of_sample_particle > number_of_sample_particle_limit) ? number_of_sample_particle : number_of_sample_particle_limit;
#endif
            number_of_sample_particle = (number_of_sample_particle < nloc) ? number_of_sample_particle : nloc;
            //std::cout<<"number_of_sample_particle= "<<number_of_sample_particle<<std::endl;
            //std::cout<<"weight= "<<weight<<" weight_all= "<<weight_all<<std::endl;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            std::cout<<"weight="<<weight<<" weight_all="<<weight_all<<std::endl;
            std::cout<<"n_smp_ptcl_tot_="<<n_smp_ptcl_tot_<<std::endl;
            //std::cout<<"nglb="<<nglb<<" nloc="<<nloc<<std::endl;
            std::cout<<"(weight * n_smp_ptcl_tot_) / weight_all="<<(weight * n_smp_ptcl_tot_) / weight_all<<std::endl;
            std::cout<<"((S64)nloc * n_smp_ptcl_tot_)="<< ((S64)nloc * n_smp_ptcl_tot_)<<std::endl;
            //std::cout<<"((F32)nglb * (1.0 + coef_limitter))="<<((F32)nglb * (1.0 + coef_limitter))<<std::endl;
            //std::cout<<"((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter))="<<((S64)nloc * n_smp_ptcl_tot_) / ((F32)nglb * (1.0 + coef_limitter))<<std::endl;
            std::cout<<"number_of_sample_particle(final)="<<number_of_sample_particle<<std::endl;
#endif

            const F64 MY_PI = 4.0*atan(1.0);
            S32 *record = new S32[number_of_sample_particle];
            for(S32 i = 0; i < number_of_sample_particle; i++) {
                    S32 j = getUniformDistributionFromArg1ToArg2(i, nloc-1);
                    Tptcl hold = ptcl_[j];
                    ptcl_[j]   = ptcl_[i];
                    ptcl_[i]   = hold;
                    record[i]  = j;
                }
                for(S32 i = 0; i < number_of_sample_particle; i++) {
    #ifdef PHI_R_TREE
                    const F64 x_tmp = ptcl_[i].getPos().x + 1.0e-12*(frand()-0.5);
                    const F64 y_tmp = ptcl_[i].getPos().y + 1.0e-12*(frand()-0.5);
                    F64 phi = atan2(y_tmp, x_tmp);
                    if(phi < 0.0) phi += 2.0*MY_PI;
                    else if(phi >= 2.0*MY_PI) phi -= 2.0*MY_PI;
                    pos_sample[i].x = phi;
                    pos_sample[i].y = sqrt(x_tmp*x_tmp + y_tmp*y_tmp);
                    pos_sample[i].z = ptcl_[i].getPos().z;
    #else
                    pos_sample[i] = ptcl_[i].getPos();
    #endif
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


        template<class Tdinfo>
        void exchangeParticle2(Tdinfo & dinfo) {
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK0 @exchangeParticle2"<<std::endl;
            //for(S32 ip = 0; ip < nloc; ip++) {
            //    F64vec pos = ptcl_[ip].getPos();
            //    if ((std::isfinite(pos.x) != true) ||
            //        (std::isfinite(pos.y) != true) ||
            //        (std::isfinite(pos.z) != true)) {
            //        std::cout << "[before] nan detected: " << rank << " "
            //                  << nloc << " " << ip << std::endl;
            //    }
            //}
#endif
            F64 time_offset = GetWtime();
            const S32 nloc  = ptcl_.size();
            const S32 rank  = MPI::COMM_WORLD.Get_rank();
            const S32 nproc = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();

            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort thisdomain = dinfo.getPosDomain(rank);

            //* Original
            S32 * nsend  = new S32[nproc];
            S32 * nsend_disp  = new S32[nproc+1];
            S32 * nrecv  = new S32[nproc];
            S32 * nrecv_disp  = new S32[nproc+1];
            MPI::Request * req_send = new MPI::Request[nproc];
            MPI::Request * req_recv = new MPI::Request[nproc];
            for(S32 i = 0; i < nproc; i++) {
                nsend[i] = nsend_disp[i] = nrecv[i] = nrecv_disp[i] = 0;
            }
            nsend_disp[nproc] = nrecv_disp[nproc] = 0;
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK1 @exchangeParticle2"<<std::endl;
#endif


            F64 time_offset_inner = GetWtime();
            for(S32 ip = 0; ip < nloc; ip++) {
                //F64vec pos = ptcl_[ip].getPos();
                //if ((std::isfinite(pos.x) != true) ||
                //    (std::isfinite(pos.y) != true) ||
                //    (std::isfinite(pos.z) != true)) {
                //    std::cout << "[after] nan detected: " << rank << " "
                //              << nloc << " " << ip << std::endl;
                //}
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
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            //Finalize();
            //std::exit(0);
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK2 @exchangeParticle2"<<std::endl;
#endif

            nsend_disp[0] = 0;
            for(S32 i = 0; i < nproc; i++) {
                nsend_disp[i+1] += nsend_disp[i] + nsend[i];
            }
            ptcl_send_.resizeNoInitialize( nsend_disp[nproc] );
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK3 @exchangeParticle2"<<std::endl;
#endif
            // ****************************************************
            // *** align send particles on ptcl_send_ *************
            for(S32 i = 0; i < nproc; i++) nsend[i] = 0;
            S32 iloc = 0;
            DISP_WALK_EX_PTCL = 0;
            N_NOT_MOVE = 0;
            N_LOC_ORG = nloc;
            for(S32 ip = 0; ip < nloc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    N_NOT_MOVE++;
                    iloc++;
                } else {
                    S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    S32 jloc = nsend[srank] + nsend_disp[srank];
                    ptcl_send_[jloc] = ptcl_[ip];
                    nsend[srank]++;
                    /*
                    if(rank == srank){
                        std::cerr<<"rank= "<<rank
                                 <<" thisdomain= "<<thisdomain
                                 <<" ptcl_[ip].getPos()= "<<ptcl_[ip].getPos()
                                 <<std::endl;
                    }
                    */
                }
            }
            DISP_WALK_EX_PTCL_AVE = (double)DISP_WALK_EX_PTCL / (N_LOC_ORG-N_NOT_MOVE);
            std::cerr<<"rank= "<<rank<<" DISP_WALK_EX_PTCL= "<<DISP_WALK_EX_PTCL
                     <<" DISP_WALK_EX_PTCL_AVE= "<<DISP_WALK_EX_PTCL_AVE
                     <<std::endl;
            ptcl_.resizeNoInitialize(iloc);
            time_profile_.exchange_particle__find_particle += GetWtime() - time_offset_inner;

            // ****************************************************
            // *** receive the number of receive particles ********
            time_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK4 @exchangeParticle2"<<std::endl;
#endif
#if 1
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 rank_x = dinfo.getRank1d(0);
            const S32 rank_y = dinfo.getRank1d(1);
            S32 dx_max_loc = 0;
            S32 dy_max_loc = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0){
                    dx_max_loc = abs( (i / n_proc_y) - rank_x) > dx_max_loc ? abs( (i / n_proc_y) - rank_x) : dx_max_loc;
                    dy_max_loc = abs( (i % n_proc_y) - rank_y) > dy_max_loc ? abs( (i % n_proc_y) - rank_y) : dy_max_loc;
                }
            }
            S32 dx_max_glb = Comm::getMaxValue(dx_max_loc);
            S32 dy_max_glb = Comm::getMaxValue(dy_max_loc);
            for(S32 ix=-dx_max_glb; ix<=dx_max_glb; ix++){
                if(rank_x+ix<0 || rank_x+ix>=n_proc_x) continue;
                for(S32 iy=-dy_max_glb; iy<=dy_max_glb; iy++){
                    if(rank_y+iy<0 || rank_y+iy>=n_proc_y) continue;
                    S32 rank = (rank_x+ix)*n_proc_y + (rank_y+iy);
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(nsend+rank, 1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(nrecv+rank, 1, GetDataType<S32>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK5 @exchangeParticle2"<<std::endl;
#endif
            /*
            // for debug
            S32 * nrecv_tmp = new S32[nproc];
            Comm::allToAll(nsend, 1, nrecv_tmp);
            for(S32 i=0; i<nproc; i++){
                if(Comm::getRank()==0){
                    std::cerr<<"rank= "<<rank
                             <<" nrecv_tmp[i]= "<<nrecv_tmp[i]
                             <<" nrecv[i]= "<<nrecv[i]
                             <<std::endl;
                }
                //nrecv[i] = nrecv_tmp[i];
            }
            delete [] nrecv_tmp;
            */
#else
            Comm::allToAll(nsend, 1, nrecv);
#endif
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK6 @exchangeParticle2"<<std::endl;
#endif
            nrecv_disp[0] = 0;
            for(S32 i=0; i<nproc; i++){
                nrecv_disp[i+1] = nrecv_disp[i] + nrecv[i];
            }
            ptcl_recv_.resizeNoInitialize( nrecv_disp[nproc] );

            //const S32 n_proc_comm_limit = 500;
            //const S32 n_proc_comm_limit = nproc / 2;
            const S32 n_proc_comm_limit = nproc + 1;
            n_proc_send = n_proc_recv = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0) n_proc_send++;
                if(nrecv[i] > 0) n_proc_recv++;
            }
            bool flag_one_to_one_comm = ( n_proc_send < n_proc_comm_limit && n_proc_recv < n_proc_comm_limit) ? true : false;
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK7 @exchangeParticle2"<<std::endl;
#endif
            if( Comm::synchronizeConditionalBranchAND(flag_one_to_one_comm) ){
                n_proc_send = n_proc_recv = 0;
                for(S32 ib = 1; ib < nproc; ib++) {
                    S32 idsend = (ib + rank) % nproc;
                    if(nsend[idsend] > 0) {
                        S32 adrsend = nsend_disp[idsend];
                        S32 tagsend = (rank < idsend) ? rank : idsend;
                        req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), nsend[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    }
                    S32 idrecv = (nproc + rank - ib) % nproc;
                    if(nrecv[idrecv] > 0) {
                        S32 adrrecv = nrecv_disp[idrecv];
                        S32 tagrecv = (rank < idrecv) ? rank : idrecv;
                        req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), nrecv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    }
                }
                MPI::Request::Waitall(n_proc_send, req_send);
                MPI::Request::Waitall(n_proc_recv, req_recv);
            }
            else{
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (0)"<<std::endl;
                CommForAllToAll<Tptcl, 3> comm_a2a_3d;
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (1)"<<std::endl;
                comm_a2a_3d.executeV(ptcl_send_, ptcl_recv_, nsend, nrecv);
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (2)"<<std::endl;
                /*
                CommForAllToAll<Tptcl, 2> comm_a2a_2d;
                comm_a2a_2d.executeV(ptcl_send_, ptcl_recv_, nsend, nrecv);
                */
                /*
                Comm::allToAllV(ptcl_send_.getPointer(), nsend, nsend_disp,
                                ptcl_recv_.getPointer(), nrecv, nrecv_disp);
                */
            }
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK8 @exchangeParticle2"<<std::endl;
#endif
            const S32 nrecv_tot = nrecv_disp[nproc];
            ptcl_.reserveEmptyAreaAtLeast( nrecv_tot );
            for(S32 ip = 0; ip < nrecv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK9 @exchangeParticle2"<<std::endl;
#endif
            // **************************************************** 
            n_ptcl_send_ += nsend_disp[nproc];
            n_ptcl_recv_ += nrecv_disp[nproc];
            time_profile_.exchange_particle__exchange_particle += GetWtime() - time_offset_inner;

            delete [] nsend;
            delete [] nsend_disp;
            delete [] nrecv;
            delete [] nrecv_disp;
            delete [] req_send;
            delete [] req_recv;
            time_profile_.exchange_particle += GetWtime() - time_offset;
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK10 @exchangeParticle2"<<std::endl;
#endif
        }

#if 1
      // new version (with switch)
      // for search, use tree with level 3(x, y, z) 
      // must be consistend with geometry of domains.
        template<class Tdinfo>
        void exchangeParticle(Tdinfo & dinfo) {
            F64 time_offset = GetWtime();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const S32 nloc  = ptcl_.size();
            const S32 rank  = MPI::COMM_WORLD.Get_rank();
            const S32 nproc = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();

            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort thisdomain = dinfo.getPosDomain(rank);
    
            Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK1 @exchangeParticle()" << std::endl;

            S32 * nsend  = new S32[nproc];
            S32 * nsend_disp  = new S32[nproc+1];
            S32 * nrecv  = new S32[nproc];
            S32 * nrecv_disp  = new S32[nproc+1];
            MPI::Request * req_send = new MPI::Request[nproc];
            MPI::Request * req_recv = new MPI::Request[nproc];
            for(S32 i = 0; i < nproc; i++) {
                nsend[i] = nsend_disp[i] = nrecv[i] = nrecv_disp[i] = 0;
            }
            nsend_disp[nproc] = nrecv_disp[nproc] = 0;
            Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK2 @exchangeParticle()" << std::endl;
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
            Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK3 @exchangeParticle()" << std::endl;
            //Finalize();
            //std::exit(0);

            nsend_disp[0] = 0;
            for(S32 i = 0; i < nproc; i++) {
                nsend_disp[i+1] += nsend_disp[i] + nsend[i];
            }
            ptcl_send_.resizeNoInitialize( nsend_disp[nproc] );
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
                    ptcl_send_[jloc] = ptcl_[ip];
                    nsend[srank]++;
                }
            }
            ptcl_.resizeNoInitialize(iloc);
            //time_profile_.exchange_particle__find_particle = GetWtime() - time_offset_inner;
            time_profile_.exchange_particle__find_particle += GetWtime() - time_offset_inner;

            Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK4 @exchangeParticle()" << std::endl;

            // ****************************************************
            // *** receive the number of receive particles ********
            time_offset_inner = GetWtime();
            //Comm::allToAll(nsend, 1, nrecv);
            //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (0)"<<std::endl;
            CommForAllToAll<S32, 3> comm_a2a_3d;
            //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (1)"<<std::endl;
            comm_a2a_3d.execute(nsend, 1, nrecv);
            //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (2)"<<std::endl;
            
            nrecv_disp[0] = 0;
            for(S32 i=0; i<nproc; i++){
                nrecv_disp[i+1] = nrecv_disp[i] + nrecv[i];
            }
            ptcl_recv_.resizeNoInitialize( nrecv_disp[nproc] );
            Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK5 @exchangeParticle()" << std::endl;

            const S32 n_proc_comm_limit = 500;
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0) n_proc_send++;
                if(nrecv[i] > 0) n_proc_recv++;
            }
            bool flag_one_to_one_comm = ( n_proc_send < n_proc_comm_limit && n_proc_recv < n_proc_comm_limit) ? true : false;
            Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK6 @exchangeParticle()" << std::endl;

            if( Comm::synchronizeConditionalBranchAND(flag_one_to_one_comm) ){
                n_proc_send = n_proc_recv = 0;
                for(S32 ib = 1; ib < nproc; ib++) {
                    S32 idsend = (ib + rank) % nproc;
                    if(nsend[idsend] > 0) {
                        S32 adrsend = nsend_disp[idsend];
                        S32 tagsend = (rank < idsend) ? rank : idsend;
                        req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), nsend[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    }
                    S32 idrecv = (nproc + rank - ib) % nproc;
                    if(nrecv[idrecv] > 0) {
                        S32 adrrecv = nrecv_disp[idrecv];
                        S32 tagrecv = (rank < idrecv) ? rank : idrecv;
                        req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), nrecv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    }
                }
                MPI::Request::Waitall(n_proc_send, req_send);
                MPI::Request::Waitall(n_proc_recv, req_recv);
            }
            else{
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoallv (0)"<<std::endl;
                CommForAllToAll<Tptcl, 3> comm_a2av_3d;
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoallv (1)"<<std::endl;
                comm_a2av_3d.executeV(ptcl_send_, ptcl_recv_, nsend, nrecv);
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoallv (2)"<<std::endl;
                /*
                Comm::allToAllV(ptcl_send_.getPointer(), nsend, nsend_disp,
                                ptcl_recv_.getPointer(), nrecv, nrecv_disp);
                */
            }
            Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK7 @exchangeParticle()" << std::endl;

            const S32 nrecv_tot = nrecv_disp[nproc];
            ptcl_.reserveEmptyAreaAtLeast( nrecv_tot );
            for(S32 ip = 0; ip < nrecv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
            Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK8 @exchangeParticle()" << std::endl;
            // **************************************************** 
            n_ptcl_send_ += nsend_disp[nproc];
            n_ptcl_recv_ += nrecv_disp[nproc];
            time_profile_.exchange_particle__exchange_particle += GetWtime() - time_offset_inner;

            delete [] nsend;
            delete [] nsend_disp;
            delete [] nrecv;
            delete [] nrecv_disp;
            delete [] req_send;
            delete [] req_recv;
#else
            n_ptcl_send_ = 0;
            n_ptcl_recv_ = 0;
#endif
            time_profile_.exchange_particle += GetWtime() - time_offset;
        }
#else

        template<class Tdinfo>
        void exchangeParticle(Tdinfo & dinfo) {
            F64 time_offset = GetWtime();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            //const S32 nloc  = n_ptcl_;
            const S32 nloc  = ptcl_.size();
            const S32 rank  = MPI::COMM_WORLD.Get_rank();
            const S32 nproc = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();

            /* AT_DEBUG
            F32ort * pos_domain = dinfo.getPointerOfPosDomain();
            F32ort thisdomain = dinfo.getPosDomain(rank);
            */
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort thisdomain = dinfo.getPosDomain(rank);

            S32 nsendtot = 0;
            S32 *nsend0  = new S32[nproc];
            S32 *nsend1  = new S32[nproc];
            S32 nrecvtot = 0;
            S32 *nrecv0  = new S32[nproc];
            S32 *nrecv1  = new S32[nproc];
            for(S32 i = 0; i < nproc; i++) {
                nsend0[i] = 0;
                nsend1[i] = 0;
                nrecv0[i] = 0;
                nrecv1[i] = 0;
            }
            MPI::Request * req_send = new MPI::Request[nproc];
            MPI::Request * req_recv = new MPI::Request[nproc];
            // *** count the number of send particles preliminary *
            for(S32 ip = 0; ip < nloc; ip++) {
                if( dinfo.getPosRootDomain().notOverlapped(ptcl_[ip].getPos()) ){
                    PARTICLE_SIMULATOR_PRINT_ERROR("A particle is out of root domain");
                    std::cerr<<"position of the particle="<<ptcl_[ip].getPos()<<std::endl;
                    std::cerr<<"position of the root domain="<<dinfo.getPosRootDomain()<<std::endl;
                    Abort(-1);
                }
                if(!determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    nsend1[srank]++;
                    nsendtot++;
                }
            }
            nsend0[0] = 0;
            for(S32 i = 1; i < nproc; i++) {
                nsend0[i] += nsend0[i-1] + nsend1[i-1];
            }
            //ptcl_send_ = new Tptcl[nsendtot];
            ptcl_send_.resizeNoInitialize(nsendtot);
            // ****************************************************
            // *** align send particles on ptcl_send_ *************
            //std::cerr<<"check 0"<<std::endl;
            for(S32 i = 0; i < nproc; i++)
                nsend1[i] = 0;
            S32 iloc = 0;
            for(S32 ip = 0; ip < nloc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    iloc++;
                } else {
                    S32 srank = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    S32 jloc = nsend0[srank] + nsend1[srank];
                    ptcl_send_[jloc] = ptcl_[ip];
                    nsend1[srank]++;
                }
            }           
            //n_ptcl_ = iloc;
	    //std::cerr<<"check 1"<<std::endl;
	    ptcl_.resizeNoInitialize(iloc);
	    //std::cerr<<"ptcl_.size()="<<ptcl_.size()<<std::endl;
/*
            if(rank == 0) {
                char filename[1024];
                FILE *fp;
                sprintf(filename, "out/send_%04d_%04d.txt", rank, rank);
                fp = fopen(filename, "w");
                for(S32 j = 0; j < n_ptcl_; j++)
                    fprintf(fp, "%+e %+e %+e\n",
                            ptcl_[j].getPos()[0],
                            ptcl_[j].getPos()[1],
                            ptcl_[j].getPos()[2]);
                fclose(fp);
                for(S32 i = 1; i < nproc; i++) {
                    S32 srank = (rank + i) % nproc;
                    sprintf(filename, "out/send_%04d_%04d.txt", rank, srank);
                    fp = fopen(filename, "w");
                    S32 next = (srank + 1 < nproc) ? nsend0[srank+1] : nsendtot;
                    for(S32 j = nsend0[srank]; j < next; j++)
                        fprintf(fp, "%+e %+e %+e\n",
                                ptcl_send_[j].getPos()[0],
                                ptcl_send_[j].getPos()[1],
                                ptcl_send_[j].getPos()[2]);
                    fclose(fp);
                }
            }
*/
            // ****************************************************
            // *** receive the number of receive particles ********
            MPI::COMM_WORLD.Alltoall(nsend1, 1, GetDataType<S32>(), nrecv1, 1, GetDataType<S32>());
	    //std::cerr<<"check 2"<<std::endl;
            for(S32 i = 0; i < nproc; i++)
                nrecvtot += nrecv1[i];
            //assert(n_ptcl_ + nrecvtot <= n_ptcl_limit_);
	    //ptcl_.reserve(n_ptcl_ + nrecvtot);
	    //ptcl_.resizeNoInitialize(n_ptcl_ + nrecvtot);
            //ptcl_recv_ = new Tptcl[nrecvtot];
	    ptcl_recv_.resizeNoInitialize(nrecvtot);
	    //std::cerr<<"check 3"<<std::endl;
/*
            {
                char filename[1024];
                FILE *fp;
                sprintf(filename, "out/send_%04d.txt", rank);
                fp = fopen(filename, "w");
                fprintf(fp, "%d", rank);
                for(S32 i = 0; i < nproc; i++) {
                    fprintf(fp, "%6d", nsend1[i]);
                }
                fprintf(fp, "\n");
                fclose(fp);
                sprintf(filename, "out/recv_%04d.txt", rank);
                fp = fopen(filename, "w");
                fprintf(fp, "%d", rank);
                for(S32 i = 0; i < nproc; i++) {
                    fprintf(fp, "%6d", nrecv[i]);
                }
                fprintf(fp, "\n");
                fclose(fp);
            }
*/
            // ****************************************************
            // *** send and receive particles *********************
            nrecv0[0] = 0;
            for(S32 i = 1; i < nproc; i++) {
                nrecv0[i] += nrecv0[i-1] + nrecv1[i-1];
            }
            S32 nsendnode = 0;
            S32 nrecvnode = 0;
	    //std::cerr<<"check 4"<<std::endl;
            for(S32 ib = 1; ib < nproc; ib++) {
                S32 idsend = (ib + rank) % nproc;
                if(nsend1[idsend] > 0) {
                    S32 adrsend = nsend0[idsend];                    
                    S32 tagsend = (rank < idsend) ? rank : idsend;                    
                    //req_send[nsendnode] = MPI::COMM_WORLD.Isend(ptcl_send_+adrsend, nsend1[idsend]*sizeof(Tptcl), MPI::BYTE, idsend, tagsend);
		    req_send[nsendnode] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), nsend1[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    nsendnode++;
                }                
                S32 idrecv = (nproc + rank - ib) % nproc;
                if(nrecv1[idrecv] > 0) {
                    S32 adrrecv = nrecv0[idrecv];
                    S32 tagrecv = (rank < idrecv) ? rank : idrecv;
                    //req_recv[nrecvnode] = MPI::COMM_WORLD.Irecv(ptcl_recv_+adrrecv, nrecv1[idrecv]*sizeof(Tptcl), MPI::BYTE, idrecv, tagrecv);
		    req_recv[nrecvnode] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), nrecv1[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    nrecvnode++;
                }
            }
            MPI::Request::Waitall(nsendnode, req_send);
            MPI::Request::Waitall(nrecvnode, req_recv);
            //std::cerr<<"check 5"<<std::endl;
            // ****************************************************            
            // *** align particles ********************************
            /*
              for(S32 ip = 0; ip < nrecvtot; ip++) {
              ptcl_[n_ptcl_] = ptcl_recv_[ip];
              n_ptcl_++;
            }
	    */
            ptcl_.reserve( ptcl_.size()+nrecvtot );
            //ptcl_.dump("dump ptcl");
            for(S32 ip = 0; ip < nrecvtot; ip++) {
                ptcl_.pushBackNoCheck(ptcl_recv_[ip]);
            }
	    //std::cerr<<"ptcl_.size()="<<ptcl_.size()<<std::endl;
            // ****************************************************            

            delete [] nsend0;
            delete [] nsend1;
            delete [] nrecv0;
            delete [] nrecv1;
            delete [] req_send;
            delete [] req_recv;
            //delete [] ptcl_send_;
            //delete [] ptcl_recv_;
#endif
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            std::cout<<"ptcl_.size()="<<ptcl_.size()<<std::endl;
            if(ptcl_.size() > 0){
                std::cout<<"ptcl_[0].getPos()="<<ptcl_[0].getPos()<<std::endl;
                std::cout<<"ptcl_[0].getRSearch()="<<ptcl_[0].getRSearch()<<std::endl;
            }
#endif
            //time_profile_.exchange_particle = GetWtime() - time_offset;
            time_profile_.exchange_particle += GetWtime() - time_offset;
        }
#endif

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
                F64vec pos_new = ptcl_[i].getPos() ;
                //if( pos_root.notOverlapped(pos_new) ){
                    while(pos_new.x < pos_root.low_.x){
                        pos_new.x += len_root.x;
                    }
                    while(pos_new.x > pos_root.high_.x){
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
            return ptcl_.getMemSize() + ptcl_send_.getMemSize() + ptcl_recv_.getMemSize();
        }
	/*
        CountT getNumberOfParticleSendLocal() const { return (CountT)ptcl_send_.size(); }
        CountT getNumberOfParticleRecvLocal() const { return (CountT)ptcl_recv_.size(); }
        CountT getNumberOfParticleSendGlobal() const { return Comm::getSum((CountT)ptcl_send_.size()); }
        CountT getNumberOfParticleRecvGlobal() const { return Comm::getSum((CountT)ptcl_recv_.size()); }
	*/
        CountT getNumberOfParticleSendLocal() const { return (CountT)n_ptcl_send_; }
        CountT getNumberOfParticleRecvLocal() const { return (CountT)n_ptcl_recv_; }
        CountT getNumberOfParticleSendGlobal() const { return Comm::getSum((CountT)n_ptcl_send_); }
        CountT getNumberOfParticleRecvGlobal() const { return Comm::getSum((CountT)n_ptcl_recv_); }
        void clearCounterAll(){
            n_ptcl_send_ = n_ptcl_recv_ = 0;
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

        template<class Tdinfo>
        void exchangeParticle3(Tdinfo & dinfo) {
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK0 @exchangeParticle3"<<std::endl;
#endif
            F64 time_offset = GetWtime();
            const S32 nloc  = ptcl_.size();
            const S32 my_rank  = MPI::COMM_WORLD.Get_rank();
            const S32 nproc = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort thisdomain = dinfo.getPosDomain(my_rank);
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 my_rank_x = dinfo.getRank1d(0);
            const S32 my_rank_y = dinfo.getRank1d(1);

#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK1 @exchangeParticle3"<<std::endl;
#endif
            
            S32 * nsend  = new S32[nproc];
            S32 * nsend_disp  = new S32[nproc+1];
            S32 * nrecv  = new S32[nproc];
            S32 * nrecv_disp  = new S32[nproc+1];
            MPI::Request * req_send = new MPI::Request[nproc];
            MPI::Request * req_recv = new MPI::Request[nproc];
            for(S32 i = 0; i < nproc; i++) {
                nsend[i] = nsend_disp[i] = nrecv[i] = nrecv_disp[i] = 0;
            }
            nsend_disp[nproc] = nrecv_disp[nproc] = 0;
            F64 time_offset_inner = GetWtime();

#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK2 @exchangeParticle3"<<std::endl;
#endif

#if 1
            S32 iloc = 0;
            S32 n_send_tot = 0;
            S32 * id_to_rank = new S32[nloc];
            bool break_flag = false;
            for(S32 ip = 0; ip < nloc; ip++) {
                break_flag = false;
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    id_to_rank[ip] = my_rank;
                    iloc++;
                }
                else{
                    for(S32 lev=1;; lev++){
                        S32 lev_x = lev;
                        S32 lev_y = lev;
                        for(S32 ix=-lev_x; ix<=lev_x; ix++){
                            for(S32 iy=-lev_y; iy<=lev_y; iy++){
                                if( std::abs(ix) !=lev_x &&  std::abs(iy) != lev_y) continue;
                                if( (my_rank_x+ix) < 0 || (my_rank_x+ix) >= n_proc_x || (my_rank_y+iy) < 0 || (my_rank_y+iy) >= n_proc_y ) continue;
                                S32 target_rank = (my_rank_x+ix)*n_proc_y + (my_rank_y+iy);
                                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), pos_domain[target_rank])){
                                    id_to_rank[ip] = target_rank;
                                    nsend[target_rank]++;
                                    break_flag = true;
                                }
                                if(break_flag) break;
                            }
                            if(break_flag) break;
                        }
                        if(break_flag) break;
                    }
                }
            }
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK3 @exchangeParticle3"<<std::endl;
#endif
            ptcl_.resizeNoInitialize(iloc);
            iloc = 0;
            nsend_disp[0] = 0;
            for(S32 i = 0; i < nproc; i++) {
                nsend_disp[i+1] += nsend_disp[i] + nsend[i];
            }
            ptcl_send_.resizeNoInitialize( nsend_disp[nproc] );
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK4 @exchangeParticle3"<<std::endl;
#endif
            for(S32 ip=0; ip<nloc; ip++){
                int rank_tmp = id_to_rank[ip];
                if( rank_tmp == my_rank){
                    ptcl_[iloc] = ptcl_[ip];
                    iloc++;                    
                }
                else{
                    ptcl_send_[nsend_disp[rank_tmp]] = ptcl_[ip];
                    nsend_disp[rank_tmp]++; // temporally nsend_disp increase here.
                }
            }
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK5 @exchangeParticle3"<<std::endl;
#endif
            nsend_disp[0] = 0;
            for(S32 i = 0; i < nproc; i++) {
                nsend_disp[i+1] = nsend_disp[i] + nsend[i];
            }
            delete [] id_to_rank;
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK6 @exchangeParticle3"<<std::endl;
#endif
#else
            //original
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
            ptcl_send_.resizeNoInitialize( nsend_disp[nproc] );
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 2"<<std::endl;
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
                    ptcl_send_[jloc] = ptcl_[ip];
                    nsend[srank]++;
                }
            }
            ptcl_.resizeNoInitialize(iloc);            
#endif
            time_profile_.exchange_particle__find_particle += GetWtime() - time_offset_inner;
            
            // ****************************************************
            // *** receive the number of receive particles ********
            time_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 3"<<std::endl;
#if 1
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK7 @exchangeParticle3"<<std::endl;
#endif
            S32 dx_max_loc = 0;
            S32 dy_max_loc = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0){
                    dx_max_loc = abs( (i / n_proc_y) - my_rank_x) > dx_max_loc ? abs( (i / n_proc_y) - my_rank_x) : dx_max_loc;
                    dy_max_loc = abs( (i % n_proc_y) - my_rank_y) > dy_max_loc ? abs( (i % n_proc_y) - my_rank_y) : dy_max_loc;
                }
            }
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK8 @exchangeParticle3"<<std::endl;
#endif
            S32 dx_max_glb = Comm::getMaxValue(dx_max_loc);
            S32 dy_max_glb = Comm::getMaxValue(dy_max_loc);
            for(S32 ix=-dx_max_glb; ix<=dx_max_glb; ix++){
                if(my_rank_x+ix<0 || my_rank_x+ix>=n_proc_x) continue;
                for(S32 iy=-dy_max_glb; iy<=dy_max_glb; iy++){
                    if(my_rank_y+iy<0 || my_rank_y+iy>=n_proc_y) continue;
                    S32 target_rank = (my_rank_x+ix)*n_proc_y + (my_rank_y+iy);
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(nsend+target_rank, 1, GetDataType<S32>(), target_rank, 0);
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(nrecv+target_rank, 1, GetDataType<S32>(), target_rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK9 @exchangeParticle3"<<std::endl;
#endif
            /*
            // for debug
            S32 * nrecv_tmp = new S32[nproc];
            Comm::allToAll(nsend, 1, nrecv_tmp);
            for(S32 i=0; i<nproc; i++){
                if(Comm::getRank()==0){
                    std::cerr<<"rank= "<<rank
                             <<" nrecv_tmp[i]= "<<nrecv_tmp[i]
                             <<" nrecv[i]= "<<nrecv[i]
                             <<std::endl;
                }
                //nrecv[i] = nrecv_tmp[i];
            }
            delete [] nrecv_tmp;
            */
#else
            Comm::allToAll(nsend, 1, nrecv);
#endif
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK10 @exchangeParticle3"<<std::endl;
#endif
            nrecv_disp[0] = 0;
            for(S32 i=0; i<nproc; i++){
                nrecv_disp[i+1] = nrecv_disp[i] + nrecv[i];
            }
            ptcl_recv_.resizeNoInitialize( nrecv_disp[nproc] );

            //const S32 n_proc_comm_limit = 500;
            //const S32 n_proc_comm_limit = nproc / 2;
            //const S32 n_proc_comm_limit = nproc + 1;
            const S32 n_proc_comm_limit = 1;
            n_proc_send = n_proc_recv = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0) n_proc_send++;
                if(nrecv[i] > 0) n_proc_recv++;
            }
            bool flag_one_to_one_comm = ( n_proc_send < n_proc_comm_limit && n_proc_recv < n_proc_comm_limit) ? true : false;
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK11 @exchangeParticle3"<<std::endl;
#endif
            if( Comm::synchronizeConditionalBranchAND(flag_one_to_one_comm) ){
                n_proc_send = n_proc_recv = 0;
                for(S32 ib = 1; ib < nproc; ib++) {
                    S32 idsend = (ib + my_rank) % nproc;
                    if(nsend[idsend] > 0) {
                        S32 adrsend = nsend_disp[idsend];
                        S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                        req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), nsend[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    }
                    S32 idrecv = (nproc + my_rank - ib) % nproc;
                    if(nrecv[idrecv] > 0) {
                        S32 adrrecv = nrecv_disp[idrecv];
                        S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                        req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), nrecv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    }
                }
                MPI::Request::Waitall(n_proc_send, req_send);
                MPI::Request::Waitall(n_proc_recv, req_recv);
            }
            else{
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (0)"<<std::endl;
                CommForAllToAll<Tptcl, 3> comm_a2a_3d;
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (1)"<<std::endl;
                comm_a2a_3d.executeV(ptcl_send_, ptcl_recv_, nsend, nrecv);
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (2)"<<std::endl;
                /*
                CommForAllToAll<Tptcl, 2> comm_a2a_2d;
                comm_a2a_2d.executeV(ptcl_send_, ptcl_recv_, nsend, nrecv);
                */
                /*
                Comm::allToAllV(ptcl_send_.getPointer(), nsend, nsend_disp,
                                ptcl_recv_.getPointer(), nrecv, nrecv_disp);
                */
            }
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK12 @exchangeParticle3"<<std::endl;
#endif
            const S32 nrecv_tot = nrecv_disp[nproc];
            ptcl_.reserveEmptyAreaAtLeast( nrecv_tot );
            for(S32 ip = 0; ip < nrecv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK13 @exchangeParticle3"<<std::endl;
#endif
            // **************************************************** 
            n_ptcl_send_ += nsend_disp[nproc];
            n_ptcl_recv_ += nrecv_disp[nproc];
            time_profile_.exchange_particle__exchange_particle += GetWtime() - time_offset_inner;

            /*
            /////////////
            // FOR DEBUG
            for(S32 ip=0; ip<ptcl_.size(); ip++) {
                if(!determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)){
                    std::cerr<<"ptcl_[ip].getPos()= "<<ptcl_[ip].getPos()
                             <<" thisdomain= "<<thisdomain
                             <<std::endl;
                    assert(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain));
                }
            }
            // FOR DEBUG
            /////////////
            */
            
            delete [] nsend;
            delete [] nsend_disp;
            delete [] nrecv;
            delete [] nrecv_disp;
            delete [] req_send;
            delete [] req_recv;

            time_profile_.exchange_particle += GetWtime() - time_offset;
#ifdef PARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK14 @exchangeParticle3"<<std::endl;
#endif
        }

        template<class Tdinfo>
        void exchangeParticle4(Tdinfo & dinfo) {
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 1"<<std::endl;
            F64 time_offset = GetWtime();
            const S32 nloc  = ptcl_.size();
            const S32 rank  = MPI::COMM_WORLD.Get_rank();
            const S32 nproc = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort thisdomain = dinfo.getPosDomain(rank);
            S32 * nsend  = new S32[nproc];
            S32 * nsend_disp  = new S32[nproc+1];
            S32 * nrecv  = new S32[nproc];
            S32 * nrecv_disp  = new S32[nproc+1];
            MPI::Request * req_send = new MPI::Request[nproc];
            MPI::Request * req_recv = new MPI::Request[nproc];
            for(S32 i = 0; i < nproc; i++) {
                nsend[i] = nsend_disp[i] = nrecv[i] = nrecv_disp[i] = 0;
            }
            nsend_disp[nproc] = nrecv_disp[nproc] = 0;
            F64 time_offset_inner = GetWtime();
            S32 previous_rank = rank;
            for(S32 ip = 0; ip < nloc; ip++) {
                if( dinfo.getPosRootDomain().notOverlapped(ptcl_[ip].getPos()) ){
                    PARTICLE_SIMULATOR_PRINT_ERROR("A particle is out of root domain");
                    std::cerr<<"position of the particle="<<ptcl_[ip].getPos()<<std::endl;
                    std::cerr<<"position of the root domain="<<dinfo.getPosRootDomain()<<std::endl;
                    Abort(-1);
                }
                if(!determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    S32 srank = searchWhichDomainParticleGoTo2(ptcl_[ip].getPos(), n_domain, pos_domain, previous_rank);
                    //S32 srank_tmp = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    //assert(srank == srank_tmp);
                    nsend[srank]++;
                    previous_rank = srank;
                }
            }
            nsend_disp[0] = 0;
            for(S32 i = 0; i < nproc; i++) {
                nsend_disp[i+1] += nsend_disp[i] + nsend[i];
            }
            ptcl_send_.resizeNoInitialize( nsend_disp[nproc] );
            // ****************************************************
            // *** align send particles on ptcl_send_ *************
            for(S32 i = 0; i < nproc; i++) nsend[i] = 0;
            S32 iloc = 0;
            DISP_WALK_EX_PTCL = 0;
            N_NOT_MOVE = 0;
            N_LOC_ORG = nloc;
            previous_rank = rank;
            HIT_RATIO = 0.0;
            for(S32 ip = 0; ip < nloc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    previous_rank = rank;
                    N_NOT_MOVE++;
                    iloc++;
                } else {
                    S32 srank = searchWhichDomainParticleGoTo2(ptcl_[ip].getPos(), n_domain, pos_domain, previous_rank);
                    //S32 srank_tmp = searchWhichDomainParticleGoTo(ptcl_[ip].getPos(), n_domain, pos_domain);
                    //assert(srank == srank_tmp);
                    S32 jloc = nsend[srank] + nsend_disp[srank];
                    ptcl_send_[jloc] = ptcl_[ip];
                    nsend[srank]++;
                    if(srank == previous_rank) HIT_RATIO += 1.0;
                    previous_rank = srank;
                }
            }
#if 0
            DISP_WALK_EX_PTCL_AVE = (double)DISP_WALK_EX_PTCL / (N_LOC_ORG-N_NOT_MOVE);
            HIT_RATIO /= (N_LOC_ORG-N_NOT_MOVE);
#endif
            /*
            std::cerr<<"rank= "<<rank<<" DISP_WALK_EX_PTCL= "<<DISP_WALK_EX_PTCL
                     <<" DISP_WALK_EX_PTCL_AVE= "<<DISP_WALK_EX_PTCL_AVE
                     <<" HIT_RATIO= "<<HIT_RATIO
                     <<std::endl;            
            */
            ptcl_.resizeNoInitialize(iloc);
            time_profile_.exchange_particle__find_particle += GetWtime() - time_offset_inner;

            // ****************************************************
            // *** receive the number of receive particles ********
            time_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 3"<<std::endl;
#if 1
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 rank_x = dinfo.getRank1d(0);
            const S32 rank_y = dinfo.getRank1d(1);
            S32 dx_max_loc = 0;
            S32 dy_max_loc = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0){
                    dx_max_loc = abs( (i / n_proc_y) - rank_x) > dx_max_loc ? abs( (i / n_proc_y) - rank_x) : dx_max_loc;
                    dy_max_loc = abs( (i % n_proc_y) - rank_y) > dy_max_loc ? abs( (i % n_proc_y) - rank_y) : dy_max_loc;
                }
            }
            S32 dx_max_glb = Comm::getMaxValue(dx_max_loc);
            S32 dy_max_glb = Comm::getMaxValue(dy_max_loc);
            for(S32 ix=-dx_max_glb; ix<=dx_max_glb; ix++){
                if(rank_x+ix<0 || rank_x+ix>=n_proc_x) continue;
                for(S32 iy=-dy_max_glb; iy<=dy_max_glb; iy++){
                    if(rank_y+iy<0 || rank_y+iy>=n_proc_y) continue;
                    S32 rank = (rank_x+ix)*n_proc_y + (rank_y+iy);
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(nsend+rank, 1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(nrecv+rank, 1, GetDataType<S32>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
            /*
            // for debug
            S32 * nrecv_tmp = new S32[nproc];
            Comm::allToAll(nsend, 1, nrecv_tmp);
            for(S32 i=0; i<nproc; i++){
                if(Comm::getRank()==0){
                    std::cerr<<"rank= "<<rank
                             <<" nrecv_tmp[i]= "<<nrecv_tmp[i]
                             <<" nrecv[i]= "<<nrecv[i]
                             <<std::endl;
                }
                //nrecv[i] = nrecv_tmp[i];
            }
            delete [] nrecv_tmp;
            */
#else
            Comm::allToAll(nsend, 1, nrecv);
#endif
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 4"<<std::endl;
            nrecv_disp[0] = 0;
            for(S32 i=0; i<nproc; i++){
                nrecv_disp[i+1] = nrecv_disp[i] + nrecv[i];
            }
            ptcl_recv_.resizeNoInitialize( nrecv_disp[nproc] );

            //const S32 n_proc_comm_limit = 500;
            //const S32 n_proc_comm_limit = nproc / 2;
            const S32 n_proc_comm_limit = nproc+1;
            n_proc_send = n_proc_recv = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0) n_proc_send++;
                if(nrecv[i] > 0) n_proc_recv++;
            }
            bool flag_one_to_one_comm = ( n_proc_send < n_proc_comm_limit && n_proc_recv < n_proc_comm_limit) ? true : false;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 5"<<std::endl;
            if( Comm::synchronizeConditionalBranchAND(flag_one_to_one_comm) ){
                n_proc_send = n_proc_recv = 0;
                for(S32 ib = 1; ib < nproc; ib++) {
                    S32 idsend = (ib + rank) % nproc;
                    if(nsend[idsend] > 0) {
                        S32 adrsend = nsend_disp[idsend];
                        S32 tagsend = (rank < idsend) ? rank : idsend;
                        req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), nsend[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    }
                    S32 idrecv = (nproc + rank - ib) % nproc;
                    if(nrecv[idrecv] > 0) {
                        S32 adrrecv = nrecv_disp[idrecv];
                        S32 tagrecv = (rank < idrecv) ? rank : idrecv;
                        req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), nrecv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    }
                }
                MPI::Request::Waitall(n_proc_send, req_send);
                MPI::Request::Waitall(n_proc_recv, req_recv);
            }
            else{
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (0)"<<std::endl;
                CommForAllToAll<Tptcl, 3> comm_a2a_3d;
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (1)"<<std::endl;
                comm_a2a_3d.executeV(ptcl_send_, ptcl_recv_, nsend, nrecv);
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (2)"<<std::endl;
            }
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 6"<<std::endl;
            const S32 nrecv_tot = nrecv_disp[nproc];
            ptcl_.reserveEmptyAreaAtLeast( nrecv_tot );
            for(S32 ip = 0; ip < nrecv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
            // **************************************************** 
            n_ptcl_send_ += nsend_disp[nproc];
            n_ptcl_recv_ += nrecv_disp[nproc];
            time_profile_.exchange_particle__exchange_particle += GetWtime() - time_offset_inner;

            delete [] nsend;
            delete [] nsend_disp;
            delete [] nrecv;
            delete [] nrecv_disp;
            delete [] req_send;
            delete [] req_recv;
            time_profile_.exchange_particle += GetWtime() - time_offset;
            //std::cerr<<"time_profile_.exchange_particle= "<<time_profile_.exchange_particle<<std::endl;
        }

        template<class Tdinfo>
        void exchangeParticle5(Tdinfo & dinfo) {
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK -04 @exchangeParticle5"<<std::endl;
            F64 time_offset = GetWtime();
            const S32 nloc  = ptcl_.size();
            const S32 rank  = MPI::COMM_WORLD.Get_rank();
            const S32 nproc = MPI::COMM_WORLD.Get_size();
            std::cerr<<"nproc= "<<nproc<<std::endl;
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort thisdomain = dinfo.getPosDomain(rank);
            Comm::barrier(); if(Comm::getRank()==0) std::cerr<<"OK -03 @exchangeParticle5"<<std::endl;
            S32 * nsend  = new S32[nproc];
            Comm::barrier(); if(Comm::getRank()==0) std::cerr<<"OK -02.9 @exchangeParticle5"<<std::endl;
            S32 * nsend_disp  = new S32[nproc+1];
            Comm::barrier(); if(Comm::getRank()==0) std::cerr<<"OK -02.8 @exchangeParticle5"<<std::endl;
            S32 * nrecv  = new S32[nproc];
            Comm::barrier(); if(Comm::getRank()==0) std::cerr<<"OK -02.7 @exchangeParticle5"<<std::endl;
            S32 * nrecv_disp  = new S32[nproc+1];
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK -02.5 @exchangeParticle5"<<std::endl;
            MPI::Request * req_send = new MPI::Request[nproc];
            MPI::Request * req_recv = new MPI::Request[nproc];
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK -02 @exchangeParticle5"<<std::endl;
            for(S32 i = 0; i < nproc; i++) {
                nsend[i] = nsend_disp[i] = nrecv[i] = nrecv_disp[i] = 0;
            }
            nsend_disp[nproc] = nrecv_disp[nproc] = 0;
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK -01 @exchangeParticle5"<<std::endl;
            F64 time_offset_inner = GetWtime();
            S32 previous_rank = rank;
            for(S32 i = 0; i < nproc; i++) nsend[i] = nrecv[i] = 0;
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 00 @exchangeParticle5"<<std::endl;
            S32 iloc = 0;
            DISP_WALK_EX_PTCL = 0;
            N_NOT_MOVE = 0;
            N_LOC_ORG = nloc;
            previous_rank = rank;
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 01 @exchangeParticle5"<<std::endl;
            for(S32 ip = 0; ip < nloc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), thisdomain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    previous_rank = rank;
                    N_NOT_MOVE++;
                    iloc++;
                } else {
                    S32 srank = searchWhichDomainParticleGoTo2(ptcl_[ip].getPos(), n_domain, pos_domain, previous_rank);
                    ptcl_send_buf_[srank].push_back(ptcl_[ip]);
                    if(srank == previous_rank) HIT_RATIO += 1.0;
                    previous_rank = srank;
                    nsend[srank]++;
                }
            }
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 02 @exchangeParticle5"<<std::endl;
#if 0
            DISP_WALK_EX_PTCL_AVE = (double)DISP_WALK_EX_PTCL / (N_LOC_ORG-N_NOT_MOVE);
            HIT_RATIO /= (N_LOC_ORG-N_NOT_MOVE);
#endif
            /*
            std::cerr<<"rank= "<<rank<<" DISP_WALK_EX_PTCL= "<<DISP_WALK_EX_PTCL
                     <<" DISP_WALK_EX_PTCL_AVE= "<<DISP_WALK_EX_PTCL_AVE
                     <<" HIT_RATIO= "<<HIT_RATIO
                     <<std::endl;
            */
            nsend_disp[0] = 0;
            for(S32 i = 0; i < nproc; i++) {
                nsend_disp[i+1] = nsend_disp[i] + nsend[i];
            }
            ptcl_send_.resizeNoInitialize( nsend_disp[nproc] );
            ptcl_.resizeNoInitialize(iloc);
            S32 n_cnt = 0;
            for(S32 i=0; i<nproc; i++) {
                const S32 n_tmp = nsend[i];
                for(S32 j=0; j<n_tmp; j++){
                    ptcl_send_[n_cnt++] = ptcl_send_buf_[i][j];
                }
            }
            assert(n_cnt == nloc-iloc);
            time_profile_.exchange_particle__find_particle += GetWtime() - time_offset_inner;

            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 03 @exchangeParticle5"<<std::endl;
            
            // ****************************************************
            // *** receive the number of receive particles ********
            time_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 3"<<std::endl;
#if 1
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 rank_x = dinfo.getRank1d(0);
            const S32 rank_y = dinfo.getRank1d(1);
            S32 dx_max_loc = 0;
            S32 dy_max_loc = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0){
                    dx_max_loc = abs( (i / n_proc_y) - rank_x) > dx_max_loc ? abs( (i / n_proc_y) - rank_x) : dx_max_loc;
                    dy_max_loc = abs( (i % n_proc_y) - rank_y) > dy_max_loc ? abs( (i % n_proc_y) - rank_y) : dy_max_loc;
                }
            }
            S32 dx_max_glb = Comm::getMaxValue(dx_max_loc);
            S32 dy_max_glb = Comm::getMaxValue(dy_max_loc);
            for(S32 ix=-dx_max_glb; ix<=dx_max_glb; ix++){
                if(rank_x+ix<0 || rank_x+ix>=n_proc_x) continue;
                for(S32 iy=-dy_max_glb; iy<=dy_max_glb; iy++){
                    if(rank_y+iy<0 || rank_y+iy>=n_proc_y) continue;
                    S32 rank = (rank_x+ix)*n_proc_y + (rank_y+iy);
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(nsend+rank, 1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(nrecv+rank, 1, GetDataType<S32>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 04 @exchangeParticle5"<<std::endl;
#else
            Comm::allToAll(nsend, 1, nrecv);
#endif
            nrecv_disp[0] = 0;
            for(S32 i=0; i<nproc; i++){
                nrecv_disp[i+1] = nrecv_disp[i] + nrecv[i];
            }
            ptcl_recv_.resizeNoInitialize( nrecv_disp[nproc] );

            //const S32 n_proc_comm_limit = 1;
            const S32 n_proc_comm_limit = nproc+1;
            n_proc_send = n_proc_recv = 0;
            for(S32 i=0; i<nproc; i++){
                if(nsend[i] > 0) n_proc_send++;
                if(nrecv[i] > 0) n_proc_recv++;
            }
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 05 @exchangeParticle5"<<std::endl;            
            bool flag_one_to_one_comm = ( n_proc_send < n_proc_comm_limit && n_proc_recv < n_proc_comm_limit) ? true : false;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 5"<<std::endl;
            if( Comm::synchronizeConditionalBranchAND(flag_one_to_one_comm) ){
                n_proc_send = n_proc_recv = 0;
                for(S32 ib = 1; ib < nproc; ib++) {
                    S32 idsend = (ib + rank) % nproc;
                    if(nsend[idsend] > 0) {
                        S32 adrsend = nsend_disp[idsend];
                        S32 tagsend = (rank < idsend) ? rank : idsend;
                        req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), nsend[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    }
                    S32 idrecv = (nproc + rank - ib) % nproc;
                    if(nrecv[idrecv] > 0) {
                        S32 adrrecv = nrecv_disp[idrecv];
                        S32 tagrecv = (rank < idrecv) ? rank : idrecv;
                        req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), nrecv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    }
                }
                MPI::Request::Waitall(n_proc_send, req_send);
                MPI::Request::Waitall(n_proc_recv, req_recv);
            }
            else{
                if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (0)"<<std::endl;
                CommForAllToAll<Tptcl, 3> comm_a2a_3d;
                if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (1)"<<std::endl;
                comm_a2a_3d.executeV(ptcl_send_, ptcl_recv_, nsend, nrecv);
                if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (2)"<<std::endl;
            }
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 06 @exchangeParticle5"<<std::endl;
            const S32 nrecv_tot = nrecv_disp[nproc];
            ptcl_.reserveEmptyAreaAtLeast( nrecv_tot );
            for(S32 ip = 0; ip < nrecv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
            // **************************************************** 
            n_ptcl_send_ += nsend_disp[nproc];
            n_ptcl_recv_ += nrecv_disp[nproc];
            time_profile_.exchange_particle__exchange_particle += GetWtime() - time_offset_inner;

            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 07 @exchangeParticle5"<<std::endl;
            
            delete [] nsend;
            delete [] nsend_disp;
            delete [] nrecv;
            delete [] nrecv_disp;
            delete [] req_send;
            delete [] req_recv;
            time_profile_.exchange_particle += GetWtime() - time_offset;
            //std::cerr<<"time_profile_.exchange_particle= "<<time_profile_.exchange_particle<<std::endl;
        }        

        template<class Tdinfo>
        void exchangeParticle6(Tdinfo & dinfo) {
            F64 time_offset = GetWtime();
            const F64 peri_len_x = 2.0 * 4.0 * atan(1.0);
            const S32 n_loc   = ptcl_.size();
            const S32 my_rank = MPI::COMM_WORLD.Get_rank();
            const S32 n_proc  = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
    #ifdef DEBUG_PRINT_EX_PTCL6
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 00 @exchangeParticle6"<<std::endl;
    #endif
            S32 * n_send  = new S32[n_proc];
            S32 * n_disp_send = new S32[n_proc+1];
            S32 * n_recv  = new S32[n_proc];
            S32 * n_disp_recv = new S32[n_proc+1];
            MPI::Request * req_send = new MPI::Request[n_proc];
            MPI::Request * req_recv = new MPI::Request[n_proc];
    #ifdef DEBUG_PRINT_EX_PTCL6
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 01 @exchangeParticle6"<<std::endl;
    #endif
            for(S32 i=0; i<n_proc; i++) {
                n_send[i] = n_disp_send[i] = n_recv[i] = n_disp_recv[i] = 0;
                ptcl_send_buf_[i].clearSize();
            }
            n_disp_send[n_proc] = n_disp_recv[n_proc] = 0;
            F64 wtime_offset_inner = GetWtime();
            S32 previous_rank = my_rank;
            for(S32 i=0; i<n_proc; i++) n_send[i] = n_recv[i] = 0;
            S32 iloc = 0;
            DISP_WALK_EX_PTCL = 0;
            N_NOT_MOVE = 0;
            N_LOC_ORG = n_loc;
    //#ifdef SUNWAY
    #if 0
    #else //SUNWAY
            previous_rank = my_rank;
            for(S32 ip=0; ip<n_loc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), pos_my_domain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    previous_rank = my_rank;
                    N_NOT_MOVE++;
                    iloc++;
                } else {
        #if 1
                    S32 srank = searchWhichDomainParticleGoToPeriodicX(ptcl_[ip].getPos(), n_domain, pos_domain, peri_len_x, previous_rank);
        #else
                    S32 srank = searchWhichDomainParticleGoTo2(ptcl_[ip].getPos(), n_domain, pos_domain, previous_rank);
        #endif
                    ptcl_send_buf_[srank].push_back(ptcl_[ip]);
                    if(srank == previous_rank) HIT_RATIO += 1.0;
                    previous_rank = srank;
                    n_send[srank]++;
                }
            }
    #endif //SUNWAY
            
    #ifdef DEBUG_PRINT_EX_PTCL6
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 02 @exchangeParticle6"<<std::endl;
    #endif
    #if 0
            DISP_WALK_EX_PTCL_AVE = (double)DISP_WALK_EX_PTCL / (N_LOC_ORG-N_NOT_MOVE);
            HIT_RATIO /= (N_LOC_ORG-N_NOT_MOVE);
    #endif
            n_disp_send[0] = 0;
            for(S32 i=0; i<n_proc; i++) {
                n_disp_send[i+1] = n_disp_send[i] + n_send[i];
            }
            ptcl_send_.resizeNoInitialize( n_disp_send[n_proc] );
            ptcl_.resizeNoInitialize(iloc);
            S32 n_cnt = 0;
            for(S32 i=0; i<n_proc; i++) {
                const S32 n_tmp = n_send[i];
                for(S32 j=0; j<n_tmp; j++){
                    ptcl_send_[n_cnt++] = ptcl_send_buf_[i][j];
                }
            }
            assert(n_cnt == n_loc-iloc);
            assert(n_cnt == n_disp_send[n_proc]);
            time_profile_.exchange_particle__find_particle += GetWtime() - wtime_offset_inner;
    #ifdef DEBUG_PRINT_EX_PTCL6
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 03 @exchangeParticle6"<<std::endl;
    #endif
            // ****************************************************
            // *** receive the number of receive particles ********
            wtime_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 3"<<std::endl;
        #if 1
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 rank_x = dinfo.getRank1d(0);
            const S32 rank_y = dinfo.getRank1d(1);
            S32 dx_max_loc = 0;
            S32 dy_max_loc = 0;
            for(S32 i=0; i<n_proc; i++){
                if(n_send[i] > 0){
                    S32 dx = abs( (i / n_proc_y) - rank_x);
                    if(dx > n_proc_x / 2) dx = n_domain[0] - dx;
                    dx_max_loc =  dx > dx_max_loc ? dx : dx_max_loc;
                    dy_max_loc = abs( (i % n_proc_y) - rank_y) > dy_max_loc ? abs( (i % n_proc_y) - rank_y) : dy_max_loc;
                }
            }
            S32 dx_max_glb = Comm::getMaxValue(dx_max_loc);
            S32 dy_max_glb = Comm::getMaxValue(dy_max_loc);
            #ifdef DEBUG_PRINT_EX_PTCL6
            Comm::barrier(); 
            if(Comm::getRank()==0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<std::endl;
                std::cerr<<"dy_max_glb= "<<dy_max_glb<<std::endl;
            }
            #endif
            assert(dx_max_glb*2+1 < n_domain[0]); // TO DO; have to consider this case.
            for(S32 ix=-dx_max_glb; ix<=dx_max_glb; ix++){
                S32 rank_x_target = rank_x+ix;
                if(rank_x_target < 0) rank_x_target = rank_x_target + n_proc_x;
                else if(rank_x_target >= n_proc_x) rank_x_target = rank_x_target - n_proc_x;
            #ifdef DEBUG_PRINT_EX_PTCL6                
                if(rank_x_target < 0 || rank_x_target >= n_proc_x){
                    std::cerr<<"rank= "<<Comm::getRank()
                             <<" rank_x_target= "<<rank_x_target
                             <<" ix= "<<ix
                             <<" dx_max_glb= "<<dx_max_glb
                             <<" n_proc_x= "<<n_proc_x
                             <<" rank_x= "<<rank_x
                             <<std::endl;
                }
             #endif
                assert(rank_x_target >= 0 && rank_x_target < n_proc_x);
                for(S32 iy=-dy_max_glb; iy<=dy_max_glb; iy++){
                    S32 rank_y_target = rank_y+iy;
                    if(rank_y_target < 0 || rank_y_target >= n_proc_y) continue;
                    S32 rank = (rank_x_target)*n_proc_y + (rank_y_target);
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(n_send+rank, 1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(n_recv+rank, 1, GetDataType<S32>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
            #ifdef DEBUG_PRINT_EX_PTCL6
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 04 @exchangeParticle6"<<std::endl;
            #endif
        #else
            Comm::allToAll(nsend, 1, nrecv);
        #endif
            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            ptcl_recv_.resizeNoInitialize( n_disp_recv[n_proc] );
        #ifdef DEBUG_PRINT_EX_PTCL6
            if(Comm::getRank() == 0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<" dy_max_glb= "<<dy_max_glb<<std::endl;
                for(S32 i=0; i<n_proc; i++){
                    std::cerr<<"n_send[i]="<<n_send[i]<<" n_recv[i]="<<n_recv[i]<<std::endl;
                }
            }
        #endif
            const S32 n_proc_comm_limit = n_proc+1;
            n_proc_send = n_proc_recv = 0;
            for(S32 i=0; i<n_proc; i++){
                if(n_send[i] > 0) n_proc_send++;
                if(n_recv[i] > 0) n_proc_recv++;
            }
        #ifdef DEBUG_PRINT_EX_PTCL6
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 05 @exchangeParticle6"<<std::endl;
        #endif
        #if 1
            n_proc_send = n_proc_recv = 0;
            for(S32 ib = 1; ib < n_proc; ib++) {
                S32 idsend = (ib + my_rank) % n_proc;
                if(n_send[idsend] > 0) {
                    S32 adrsend = n_disp_send[idsend];
                    S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), n_send[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                }
                S32 idrecv = (n_proc + my_rank - ib) % n_proc;
                if(n_recv[idrecv] > 0) {
                    S32 adrrecv = n_disp_recv[idrecv];
                    S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), n_recv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);            
        #else
            bool flag_one_to_one_comm = ( n_proc_send < n_proc_comm_limit && n_proc_recv < n_proc_comm_limit) ? true : false;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 5"<<std::endl;
            if( Comm::synchronizeConditionalBranchAND(flag_one_to_one_comm) ){
                n_proc_send = n_proc_recv = 0;
                for(S32 ib = 1; ib < n_proc; ib++) {
                    S32 idsend = (ib + my_rank) % n_proc;
                    if(n_send[idsend] > 0) {
                        S32 adrsend = n_disp_send[idsend];
                        S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                        req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), n_send[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    }
                    S32 idrecv = (n_proc + my_rank - ib) % n_proc;
                    if(n_recv[idrecv] > 0) {
                        S32 adrrecv = n_disp_recv[idrecv];
                        S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                        req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), n_recv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    }
                }
                MPI::Request::Waitall(n_proc_send, req_send);
                MPI::Request::Waitall(n_proc_recv, req_recv);
            }
            else{
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (0)"<<std::endl;
                CommForAllToAll<Tptcl, 3> comm_a2a_3d;
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (1)"<<std::endl;
                comm_a2a_3d.executeV(ptcl_send_, ptcl_recv_, n_send, n_recv);
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (2)"<<std::endl;
            }
        #endif
    #ifdef DEBUG_PRINT_EX_PTCL6            
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 06 @exchangeParticle6"<<std::endl;
    #endif
            const S32 n_recv_tot = n_disp_recv[n_proc];
            ptcl_.reserveEmptyAreaAtLeast( n_recv_tot );
            for(S32 ip = 0; ip < n_recv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
            // **************************************************** 
            n_ptcl_send_ += n_disp_send[n_proc];
            n_ptcl_recv_ += n_disp_recv[n_proc];
            time_profile_.exchange_particle__exchange_particle += GetWtime() - wtime_offset_inner;
        #ifdef DEBUG_PRINT_EX_PTCL6
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 07 @exchangeParticle6"<<std::endl;
        #endif            
            delete [] n_send;
            delete [] n_disp_send;
            delete [] n_recv;
            delete [] n_disp_recv;
            delete [] req_send;
            delete [] req_recv;
            time_profile_.exchange_particle += GetWtime() - time_offset;
            //std::cerr<<"time_profile_.exchange_particle= "<<time_profile_.exchange_particle<<std::endl;
        }


        template<class Tdinfo>
        void exchangeParticle7(Tdinfo & dinfo) {
            F64 time_offset = GetWtime();
            const F64 peri_len_x = 2.0 * 4.0 * atan(1.0);
            const S32 n_loc   = ptcl_.size();
            const S32 my_rank = MPI::COMM_WORLD.Get_rank();
            const S32 n_proc  = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
    #ifdef DEBUG_PRINT_EX_PTCL7
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 00 @exchangeParticle7"<<std::endl;
    #endif
            S32 * n_send  = new S32[n_proc];
            S32 * n_disp_send = new S32[n_proc+1];
            S32 * n_recv  = new S32[n_proc];
            S32 * n_disp_recv = new S32[n_proc+1];
            MPI::Request * req_send = new MPI::Request[n_proc];
            MPI::Request * req_recv = new MPI::Request[n_proc];
    #ifdef DEBUG_PRINT_EX_PTCL7
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 01 @exchangeParticle7"<<std::endl;
    #endif
            for(S32 i=0; i<n_proc; i++) {
                n_send[i] = n_disp_send[i] = n_recv[i] = n_disp_recv[i] = 0;
                ptcl_send_buf_[i].clearSize();
            }
            n_disp_send[n_proc] = n_disp_recv[n_proc] = 0;
            F64 wtime_offset_inner = GetWtime();
            S32 previous_rank = my_rank;
            for(S32 i=0; i<n_proc; i++) n_send[i] = n_recv[i] = 0;
            S32 iloc = 0;
            DISP_WALK_EX_PTCL = 0;
            N_NOT_MOVE = 0;
            N_LOC_ORG = n_loc;

    #ifdef DEBUG_PRINT_EX_PTCL7
            S32 n_out = 0;
            for(S32 i=0; i<n_loc; i++){
                if(pos_my_domain.notContained(ptcl_[i].pos)){
                    n_out++;
                }
            }
            int * adr_ex = new int[n_out];
            n_out = 0;
            for(S32 i=0; i<n_loc; i++){
                if(pos_my_domain.notContained(ptcl_[i].pos)){
                    adr_ex[n_out++] = i;
                }
            }
    #endif
            
    #ifdef SUNWAY
            Comm::barrier();
            F64 MY_PI = 4.0*atan(1.0);
            F64ort pos_my_domain_ret;
            enum{
                N_ADR_EXCHANGE_LIMIT = 100000,
            };
            S32 adr_exchange[N_ADR_EXCHANGE_LIMIT];
            S32 n_exchange = 0;
            unsigned long args[5];
            args[0] = (unsigned long)(n_loc);
            args[1] = (unsigned long)(&pos_my_domain);
            args[2] = (unsigned long) ptcl_.getPointer();
            args[3] = (unsigned long)(&adr_exchange);
            args[4] = (unsigned long)(&n_exchange);
            __real_athread_spawn((void*)slave_CheckIsInDomain, args);
            athread_join();
        #ifdef DEBUG_PRINT_EX_PTCL7
            std::cerr<<"my_rank= "<<my_rank
                     <<" n_exchange= "<<n_exchange
                     <<" n_out= "<<n_out
                     <<std::endl;            
            assert(n_out == n_exchange);
            for(S32 i=0; i<n_exchange; i++){
                assert(adr_ex[i] == adr_exchange[i]);
            }
            delete [] adr_ex;
            //exit(1);
        #endif
            iloc = n_loc; // to be consistent with original difinition (not -1)
            for(S32 ip=n_exchange-1; ip>=0; ip--){
                S32 adr = adr_exchange[ip];
                const F64 x_tmp = ptcl_[adr].getPos().x;
                const F64 y_tmp = ptcl_[adr].getPos().y;
                const F64 z_tmp = ptcl_[adr].getPos().z;
                F64 phi = atan2(y_tmp, x_tmp);
                if(phi < 0.0) phi += 2.0*MY_PI;
                F64 r = sqrt(x_tmp*x_tmp + y_tmp*y_tmp);
                S32 srank = searchWhichDomainParticleGoToPeriodicX(F64vec(phi, r, z_tmp), n_domain, pos_domain, peri_len_x, previous_rank);
                //S32 srank = searchWhichDomainParticleGoToPeriodicX(ptcl_[adr].getPos(), n_domain, pos_domain, peri_len_x, previous_rank);
                ptcl_send_buf_[srank].push_back(ptcl_[adr]);
                if(srank == previous_rank) HIT_RATIO += 1.0;
                previous_rank = srank;
                n_send[srank]++;
		std::swap(ptcl_[adr], ptcl_[iloc-1]);
                iloc--;
            }
        #ifdef DEBUG_PRINT_EX_PTCL7
            for(S32 ip=0; ip<iloc; ip++){
                assert(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), pos_my_domain));
            }
            for(S32 ip=iloc; ip<n_loc; ip++){
                assert(!determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), pos_my_domain));
            }
            for(S32 ib=0; ib<n_proc; ib++){
                for(S32 ip=0; ip<n_send[ib]; ip++){
                    assert( determineWhetherParticleIsInDomain(ptcl_send_buf_[ib][ip].getPos(), pos_domain[ib]) );
                }
            }
        #endif
    #else //SUNWAY
            previous_rank = my_rank;
            for(S32 ip=0; ip<n_loc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), pos_my_domain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    previous_rank = my_rank;
                    N_NOT_MOVE++;
                    iloc++;
                }
                else {
        #if 1
                    S32 srank = searchWhichDomainParticleGoToPeriodicX(ptcl_[ip].getPos(), n_domain, pos_domain, peri_len_x, previous_rank);
        #else
                    S32 srank = searchWhichDomainParticleGoTo2(ptcl_[ip].getPos(), n_domain, pos_domain, previous_rank);
        #endif
                    ptcl_send_buf_[srank].push_back(ptcl_[ip]);
                    if(srank == previous_rank) HIT_RATIO += 1.0;
                    previous_rank = srank;
                    n_send[srank]++;
                }
            }
    #endif //SUNWAY
            
    #ifdef DEBUG_PRINT_EX_PTCL7
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 02 @exchangeParticle7"<<std::endl;
    #endif
    #if 0
            DISP_WALK_EX_PTCL_AVE = (double)DISP_WALK_EX_PTCL / (N_LOC_ORG-N_NOT_MOVE);
            HIT_RATIO /= (N_LOC_ORG-N_NOT_MOVE);
    #endif
            n_disp_send[0] = 0;
            for(S32 i=0; i<n_proc; i++) {
                n_disp_send[i+1] = n_disp_send[i] + n_send[i];
            }
            ptcl_send_.resizeNoInitialize( n_disp_send[n_proc] );
            ptcl_.resizeNoInitialize(iloc);
            S32 n_cnt = 0;
            for(S32 i=0; i<n_proc; i++) {
                const S32 n_tmp = n_send[i];
                for(S32 j=0; j<n_tmp; j++){
                    ptcl_send_[n_cnt++] = ptcl_send_buf_[i][j];
                }
            }
            assert(n_cnt == n_loc-iloc);
            assert(n_cnt == n_disp_send[n_proc]);
            time_profile_.exchange_particle__find_particle += GetWtime() - wtime_offset_inner;
    #ifdef DEBUG_PRINT_EX_PTCL7
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 03 @exchangeParticle7"<<std::endl;
    #endif
            // ****************************************************
            // *** receive the number of receive particles ********
            wtime_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 3"<<std::endl;
        #if 1
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 rank_x = dinfo.getRank1d(0);
            const S32 rank_y = dinfo.getRank1d(1);
            S32 dx_max_loc = 0;
            S32 dy_max_loc = 0;
            for(S32 i=0; i<n_proc; i++){
                if(n_send[i] > 0){
                    S32 dx = abs( (i / n_proc_y) - rank_x);
                    if(dx > n_proc_x / 2) dx = n_domain[0] - dx;
                    dx_max_loc =  dx > dx_max_loc ? dx : dx_max_loc;
                    dy_max_loc = abs( (i % n_proc_y) - rank_y) > dy_max_loc ? abs( (i % n_proc_y) - rank_y) : dy_max_loc;
                }
            }
            S32 dx_max_glb = Comm::getMaxValue(dx_max_loc);
            S32 dy_max_glb = Comm::getMaxValue(dy_max_loc);
            #ifdef DEBUG_PRINT_EX_PTCL7
            Comm::barrier(); 
            if(Comm::getRank()==0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<std::endl;
                std::cerr<<"dy_max_glb= "<<dy_max_glb<<std::endl;
            }
            #endif
            assert(dx_max_glb*2+1 < n_domain[0]); // TO DO; have to consider this case.
            for(S32 ix=-dx_max_glb; ix<=dx_max_glb; ix++){
                S32 rank_x_target = rank_x+ix;
                if(rank_x_target < 0) rank_x_target = rank_x_target + n_proc_x;
                else if(rank_x_target >= n_proc_x) rank_x_target = rank_x_target - n_proc_x;
            #ifdef DEBUG_PRINT_EX_PTCL7
                if(rank_x_target < 0 || rank_x_target >= n_proc_x){
                    std::cerr<<"rank= "<<Comm::getRank()
                             <<" rank_x_target= "<<rank_x_target
                             <<" ix= "<<ix
                             <<" dx_max_glb= "<<dx_max_glb
                             <<" n_proc_x= "<<n_proc_x
                             <<" rank_x= "<<rank_x
                             <<std::endl;
                }
             #endif
                assert(rank_x_target >= 0 && rank_x_target < n_proc_x);
                for(S32 iy=-dy_max_glb; iy<=dy_max_glb; iy++){
                    S32 rank_y_target = rank_y+iy;
                    if(rank_y_target < 0 || rank_y_target >= n_proc_y) continue;
                    S32 rank = (rank_x_target)*n_proc_y + (rank_y_target);
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(n_send+rank, 1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(n_recv+rank, 1, GetDataType<S32>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
            #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 04 @exchangeParticle7"<<std::endl;
            #endif
        #else
            Comm::allToAll(nsend, 1, nrecv);
        #endif
            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            ptcl_recv_.resizeNoInitialize( n_disp_recv[n_proc] );
        #ifdef DEBUG_PRINT_EX_PTCL8
            if(Comm::getRank() == 0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<" dy_max_glb= "<<dy_max_glb<<std::endl;
                for(S32 i=0; i<n_proc; i++){
                    std::cerr<<"n_send[i]="<<n_send[i]<<" n_recv[i]="<<n_recv[i]<<std::endl;
                }
            }
        #endif
            const S32 n_proc_comm_limit = n_proc+1;
            n_proc_send = n_proc_recv = 0;
            for(S32 i=0; i<n_proc; i++){
                if(n_send[i] > 0) n_proc_send++;
                if(n_recv[i] > 0) n_proc_recv++;
            }
        #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 05 @exchangeParticle7"<<std::endl;
        #endif
        #if 1
            n_proc_send = n_proc_recv = 0;
            for(S32 ib = 1; ib < n_proc; ib++) {
                S32 idsend = (ib + my_rank) % n_proc;
                if(n_send[idsend] > 0) {
                    S32 adrsend = n_disp_send[idsend];
                    S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), n_send[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                }
                S32 idrecv = (n_proc + my_rank - ib) % n_proc;
                if(n_recv[idrecv] > 0) {
                    S32 adrrecv = n_disp_recv[idrecv];
                    S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), n_recv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);            
        #else
            bool flag_one_to_one_comm = ( n_proc_send < n_proc_comm_limit && n_proc_recv < n_proc_comm_limit) ? true : false;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 5"<<std::endl;
            if( Comm::synchronizeConditionalBranchAND(flag_one_to_one_comm) ){
                n_proc_send = n_proc_recv = 0;
                for(S32 ib = 1; ib < n_proc; ib++) {
                    S32 idsend = (ib + my_rank) % n_proc;
                    if(n_send[idsend] > 0) {
                        S32 adrsend = n_disp_send[idsend];
                        S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                        req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), n_send[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                    }
                    S32 idrecv = (n_proc + my_rank - ib) % n_proc;
                    if(n_recv[idrecv] > 0) {
                        S32 adrrecv = n_disp_recv[idrecv];
                        S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                        req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), n_recv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                    }
                }
                MPI::Request::Waitall(n_proc_send, req_send);
                MPI::Request::Waitall(n_proc_recv, req_recv);
            }
            else{
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (0)"<<std::endl;
                CommForAllToAll<Tptcl, 3> comm_a2a_3d;
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (1)"<<std::endl;
                comm_a2a_3d.executeV(ptcl_send_, ptcl_recv_, n_send, n_recv);
                //if(Comm::getRank() == 0) std::cerr<<"ex ptcl using 3d alltoall (2)"<<std::endl;
            }
        #endif
    #ifdef DEBUG_PRINT_EX_PTCL7
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 06 @exchangeParticle7"<<std::endl;
    #endif
            const S32 n_recv_tot = n_disp_recv[n_proc];
            ptcl_.reserveEmptyAreaAtLeast( n_recv_tot );
            for(S32 ip = 0; ip < n_recv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
            // **************************************************** 
            n_ptcl_send_ += n_disp_send[n_proc];
            n_ptcl_recv_ += n_disp_recv[n_proc];
            time_profile_.exchange_particle__exchange_particle += GetWtime() - wtime_offset_inner;
        #ifdef DEBUG_PRINT_EX_PTCL7
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 07 @exchangeParticle7"<<std::endl;
        #endif            
            delete [] n_send;
            delete [] n_disp_send;
            delete [] n_recv;
            delete [] n_disp_recv;
            delete [] req_send;
            delete [] req_recv;
            time_profile_.exchange_particle += GetWtime() - time_offset;
            //std::cerr<<"time_profile_.exchange_particle= "<<time_profile_.exchange_particle<<std::endl;
        }



        template<class Tdinfo>
        void exchangeParticle8(Tdinfo & dinfo, bool flag_reuse=false) {
            F64 time_offset = GetWtime();
            const F64 peri_len_x = 2.0 * 4.0 * atan(1.0);
            const S32 n_loc   = ptcl_.size();
            const S32 my_rank = MPI::COMM_WORLD.Get_rank();
            const S32 n_proc  = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
            #ifdef DEBUG_PRINT_EX_PTCL8
                Comm::barrier();
                if(Comm::getRank()==0)std::cerr<<"OK 00 @exchangeParticle8"<<std::endl;
            #endif
            static S32 __n_rem = 0; // not move
            static bool first = true;
            static S32 * __n_send = NULL;
            if(first){
                __n_send  = new S32[n_proc];
                first = false;
            }
            //S32 * rank_send = new S32[n_loc];
            S32 * n_disp_send = new S32[n_proc+1];
            S32 * n_recv  = new S32[n_proc];
            S32 * n_disp_recv = new S32[n_proc+1];
            MPI::Request * req_send = new MPI::Request[n_proc];
            MPI::Request * req_recv = new MPI::Request[n_proc];
            #ifdef DEBUG_PRINT_EX_PTCL8
                Comm::barrier();
                if(Comm::getRank()==0)std::cerr<<"OK 01 @exchangeParticle8"<<std::endl;
            #endif
            F64 wtime_offset_inner = GetWtime();
            S32 previous_rank = my_rank;
            S32 iloc = 0;
    #ifdef SUNWAY
            Comm::barrier();
            if(flag_reuse){
#if 0                
    #ifdef DEBUG_PRINT_EX_PTCL8
                if(my_rank == 0){
                    std::cerr<<"new loop: __n_rem= "<<__n_rem<<std::endl;
                    for(S32 i=0; i<n_proc; i++){
                        std::cerr<<"i= "<<i<<" __n_send[i]= "<<__n_send[i]<<std::endl;
                    }
                }
    #endif
                const S32 n_div = 3;
                const S32 n_crit = n_loc/n_div + 1;
                S32 n_proc_target = 0;
                S32 rank_target[n_div];
                S32 n_ptcl_target[n_div]; // used only for size estimate.
                for(S32 i=0; i<n_div; i++){
                    rank_target[i] = 0;
                    n_ptcl_target[i] = 0;
                }
                for(S32 i=0; i<n_proc; i++){
                    if( (i==my_rank && __n_rem > n_crit)
                        || (i!=my_rank && __n_send[i] > n_crit) ){
                        assert(n_proc_target < n_div);
                        rank_target[n_proc_target] = i;
                        n_ptcl_target[n_proc_target] = __n_send[i];
                        n_proc_target++;                        
                    }
                }

    #ifdef DEBUG_PRINT_EX_PTCL8
                if(my_rank == 0){
                    std::cerr<<"rank= "<<my_rank
                             <<" n_crit= "<<n_crit
                             <<" __n_rem= "<<__n_rem
                             <<" n_proc_target= "<<n_proc_target
                             <<" rank_target[0]= "<<rank_target[0]
                             <<" rank_target[1]= "<<rank_target[1]
                             <<" rank_target[2]= "<<rank_target[2]
                             <<" n_ptcl_target[0]= "<<n_ptcl_target[0]
                             <<" n_ptcl_target[1]= "<<n_ptcl_target[1]
                             <<" n_ptcl_target[2]= "<<n_ptcl_target[2]
                             <<std::endl;
                }
    #endif
                
                for(S32 i=0; i<n_proc_target; i++){
                    if(rank_target[i] == my_rank){
                        ptcl_send_buf_[rank_target[i]].reserve( (__n_rem + __n_rem/3) + 1000);
                    }
                    else{
                        ptcl_send_buf_[rank_target[i]].reserve( (n_ptcl_target[i] + n_ptcl_target[i]/3) + 1000);
                    }
                }
                for(S32 i=0; i<n_proc; i++) {
                    __n_send[i] = n_disp_send[i] = n_recv[i] = n_disp_recv[i] = 0;
                    ptcl_send_buf_[i].clearSize();
                }
                n_disp_send[n_proc] = n_disp_recv[n_proc] = 0;
                

                S32 * adr_no_set = new S32[n_loc];
                //for(S32 i=0; i<n_loc; i++) adr_no_set[i] = -2;
                S32 n_no_set = 0;
                unsigned long args[20];

                args[0] = (unsigned long)(n_loc);
                args[1] = (unsigned long)(my_rank);
                args[2] = (unsigned long)(n_domain[0]);
                args[3] = (unsigned long)(n_domain[1]);
                args[4] = (unsigned long) &peri_len_x;
                args[5] = (unsigned long) ptcl_.getPointer();
                args[6] = (unsigned long) pos_domain;
                args[7] = (unsigned long) rank_send;
                args[8] = (unsigned long)(n_proc_target);
                args[9] = (unsigned long)(adr_no_set);
                args[10] = (unsigned long)(&n_no_set);
                
                args[11] = (unsigned long)(ptcl_send_buf_[rank_target[0]].getPointer()); // OUT
                args[12] = (unsigned long)(ptcl_send_buf_[rank_target[1]].getPointer());
                args[13] = (unsigned long)(ptcl_send_buf_[rank_target[2]].getPointer());
                
                args[14] = (unsigned long)(__n_send+rank_target[0]); //OUT
                args[15] = (unsigned long)(__n_send+rank_target[1]);
                args[16] = (unsigned long)(__n_send+rank_target[2]);

                args[17] = (unsigned long)(rank_target[0]); //IN
                args[18] = (unsigned long)(rank_target[1]);
                args[19] = (unsigned long)(rank_target[2]);

                __real_athread_spawn((void*)slave_SetParticleToSendBuffer, args);
                athread_join();
                time_profile_.exchange_particle__find_particle = GetWtime() - wtime_offset_inner;

                S32 n_loc_tmp = n_no_set;
                for(S32 i=0; i<n_proc_target; i++) n_loc_tmp += __n_send[rank_target[i]];
    #ifdef DEBUG_PRINT_EX_PTCL8
                if(my_rank == 0){
                    std::cerr<<"rank= "<<my_rank
                             <<" n_crit= "<<n_crit
                             <<" __n_rem= "<<__n_rem
                             <<" n_proc_target= "<<n_proc_target
                             <<" rank_target[0]= "<<rank_target[0]
                             <<" rank_target[1]= "<<rank_target[1]
                             <<" rank_target[2]= "<<rank_target[2]
                             <<" __n_send[rank_target[0]]= "<<__n_send[rank_target[0]]
                             <<" __n_send[rank_target[1]]= "<<__n_send[rank_target[1]]
                             <<" __n_send[rank_target[2]]= "<<__n_send[rank_target[2]]
                             <<" n_ptcl_target[0]= "<<n_ptcl_target[0]
                             <<" n_ptcl_target[1]= "<<n_ptcl_target[1]
                             <<" n_ptcl_target[2]= "<<n_ptcl_target[2]
                             <<" n_no_set= "<<n_no_set
                             <<" n_loc= "<<n_loc
                             <<" n_loc_tmp= "<<n_loc_tmp
                             <<std::endl;
                }
    #endif
                assert(n_loc == n_loc_tmp);

    #ifdef DEBUG_PRINT_EX_PTCL8
                for(S32 i=0; i<n_proc_target; i++){
                    for(S32 j=0; j<n_ptcl_target[i]; j++){
                        F64 pos_x = ptcl_send_buf_[rank_target[i]][j].pos.x;
                        F64 pos_y = ptcl_send_buf_[rank_target[i]][j].pos.y;
                        F64 pos_z = ptcl_send_buf_[rank_target[i]][j].pos.z;
                        F64 pos_phi = atan2(pos_y, pos_x);
                        if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                        F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                        assert( determineWhetherParticleIsInDomain(F64vec(pos_phi, pos_r, pos_z), pos_domain[rank_target[i]]) );
                    }
                    for(S32 j=0; j<n_no_set; j++){
                        F64 pos_x = ptcl_[adr_no_set[j]].pos.x;
                        F64 pos_y = ptcl_[adr_no_set[j]].pos.y;
                        F64 pos_z = ptcl_[adr_no_set[j]].pos.z;
                        F64 pos_phi = atan2(pos_y, pos_x);
                        if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                        F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                        assert( !determineWhetherParticleIsInDomain( F64vec(pos_phi, pos_r, pos_z), pos_domain[rank_target[i]]) );
                    }
                }
    #endif
                //exit(1);
                
                for(S32 i=0; i<n_proc_target; i++){
                    ptcl_send_buf_[rank_target[i]].resizeNoInitialize(__n_send[rank_target[i]]);
                }
                
                int flag_my_rank_copy = 0;
                for(S32 i=0; i<n_proc_target; i++){
                    if(my_rank == rank_target[i]){
                        flag_my_rank_copy = 1;
                        break;
                    }
                }

                iloc = 0;
                if(flag_my_rank_copy){
                    for(S32 i=0; i<n_no_set; i++){
                        S32 adr_tmp = adr_no_set[i];
                        int rank_tmp = rank_send[adr_tmp];
                        ptcl_send_buf_[rank_tmp].push_back(ptcl_[adr_tmp]);
                        __n_send[rank_tmp]++;
    #ifdef DEBUG_PRINT_EX_PTCL8
                        F64 pos_x = ptcl_[adr_tmp].pos.x;
                        F64 pos_y = ptcl_[adr_tmp].pos.y;
                        F64 pos_z = ptcl_[adr_tmp].pos.z;
                        F64 pos_phi = atan2(pos_y, pos_x);
                        if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                        F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                        assert( determineWhetherParticleIsInDomain( F64vec(pos_phi, pos_r, pos_z), pos_domain[rank_tmp]) );
    #endif
                    }
                    for(S32 i=0; i<__n_send[my_rank]; i++, iloc++){
                        ptcl_[iloc] = ptcl_send_buf_[my_rank][i];
    #ifdef DEBUG_PRINT_EX_PTCL8
                        F64 pos_x = ptcl_[iloc].pos.x;
                        F64 pos_y = ptcl_[iloc].pos.y;
                        F64 pos_z = ptcl_[iloc].pos.z;
                        F64 pos_phi = atan2(pos_y, pos_x);
                        if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                        F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                        assert( determineWhetherParticleIsInDomain( F64vec(pos_phi, pos_r, pos_z), pos_domain[my_rank]) );
    #endif
                    }
                }
                else{
                    for(S32 i=0; i<n_no_set; i++){
                        S32 adr_tmp = adr_no_set[i];
                        int rank_tmp = rank_send[adr_tmp];
                        if(rank_tmp == my_rank){
                            ptcl_[iloc++] = ptcl_[adr_tmp];
                        }
                        else{
                            ptcl_send_buf_[rank_tmp].push_back(ptcl_[adr_tmp]);
                            __n_send[rank_tmp]++;
                        }
                    }
                }
    #ifdef DEBUG_PRINT_EX_PTCL8
                for(S32 i=0; i<iloc; i++){
                    F64 pos_x = ptcl_[i].pos.x;
                    F64 pos_y = ptcl_[i].pos.y;
                    F64 pos_z = ptcl_[i].pos.z;
                    F64 pos_phi = atan2(pos_y, pos_x);
                    if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                    F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                    assert( determineWhetherParticleIsInDomain( F64vec(pos_phi, pos_r, pos_z), pos_domain[my_rank]) );
                }

                for(S32 i=0; i<n_proc; i++){
                    if(i == my_rank) continue;
                    for(S32 j=0; j<__n_send[i]; j++){
                        F64 pos_x = ptcl_send_buf_[i][j].pos.x;
                        F64 pos_y = ptcl_send_buf_[i][j].pos.y;
                        F64 pos_z = ptcl_send_buf_[i][j].pos.z;
                        F64 pos_phi = atan2(pos_y, pos_x);
                        if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                        F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);                        
                        assert( determineWhetherParticleIsInDomain( F64vec(pos_phi, pos_r, pos_z), pos_domain[i]) );
                    }
                }
    #endif
                __n_rem = iloc;
                delete [] adr_no_set;
                //exit(1);
#endif // #if 0 if sentence
            }
            else{
                #ifdef DEBUG_PRINT_EX_PTCL8
                    Comm::barrier();
                    if(Comm::getRank()==0)std::cerr<<"OK 02 @exchangeParticle8"<<std::endl;
                #endif
                for(S32 i=0; i<n_proc; i++) {
                    __n_send[i] = n_disp_send[i] = n_recv[i] = n_disp_recv[i] = 0;
                    ptcl_send_buf_[i].clearSize();
                }
                n_disp_send[n_proc] = n_disp_recv[n_proc] = 0;
#if 1
                #ifdef DEBUG_PRINT_EX_PTCL8
                    Comm::barrier();
                    if(Comm::getRank()==0)std::cerr<<"OK 03 @exchangeParticle8"<<std::endl;
                #endif
                double * pos_phi = new double[n_loc];
                double * pos_r = new double[n_loc];
                unsigned long args[5];
                args[0] = (unsigned long)(n_loc);
                args[1] = (unsigned long) &peri_len_x;
                args[2] = (unsigned long) ptcl_.getPointer();
                args[3] = (unsigned long) pos_phi;
                args[4] = (unsigned long) pos_r;
                __real_athread_spawn((void*)slave_GetCylCoord, args);
                athread_join();
                #ifdef DEBUG_PRINT_EX_PTCL8
                    Comm::barrier();
                    if(Comm::getRank()==0)std::cerr<<"OK 04 @exchangeParticle8"<<std::endl;
                #endif
                iloc = 0;
                previous_rank = my_rank;
                for(S32 ip=0; ip<n_loc; ip++) {
                    F64 pos_z = ptcl_[ip].getPos().z;
                    F64vec pos_tmp = F64vec(pos_phi[ip], pos_r[ip], pos_z);
                    S32 srank = searchWhichDomainParticleGoToPeriodicX(pos_tmp, n_domain, pos_domain, peri_len_x, previous_rank);
                    if(srank == my_rank){
                        ptcl_[iloc++] = ptcl_[ip];
                    }
                    else{
                        ptcl_send_buf_[srank].push_back(ptcl_[ip]);
                        __n_send[srank]++;
                    }
                    previous_rank = srank;
                }
                __n_rem = iloc;
                delete [] pos_phi;
                delete [] pos_r;
                #ifdef DEBUG_PRINT_EX_PTCL8
                    Comm::barrier();
                    if(Comm::getRank()==0)std::cerr<<"OK 05 @exchangeParticle8"<<std::endl;
                #endif
#elif 0
                unsigned long args[8];
                args[0] = (unsigned long)(n_loc);
                args[1] = (unsigned long)(my_rank);
                args[2] = (unsigned long)(n_domain[0]);
                args[3] = (unsigned long)(n_domain[1]);
                args[4] = (unsigned long) &peri_len_x;
                args[5] = (unsigned long) ptcl_.getPointer();
                args[6] = (unsigned long) pos_domain;
                args[7] = (unsigned long) rank_send;
                __real_athread_spawn((void*)slave_FindDomainParticleGoTo, args);
                athread_join();
#else
                for(S32 ip=0; ip<n_loc; ip++) {
                    F64 pos_x = ptcl_[ip].getPos().x;
                    F64 pos_y = ptcl_[ip].getPos().y;
                    F64 pos_z = ptcl_[ip].getPos().z;
                    F64 pos_phi = atan2(pos_y, pos_x);
                    if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                    F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                    F64vec pos_tmp = F64vec(pos_phi, pos_r, pos_z);
                    if(determineWhetherParticleIsInDomain(pos_tmp, pos_my_domain)) {
                        rank_send[ip] = my_rank;
                    }
                    else{
                        S32 srank = searchWhichDomainParticleGoToPeriodicX(pos_tmp, n_domain, pos_domain, peri_len_x, previous_rank);
                        previous_rank = srank;
                        rank_send[ip] = srank;
                    }
                }
#endif
                time_profile_.exchange_particle__find_particle = GetWtime() - wtime_offset_inner;

        #ifdef DEBUG_PRINT_EX_PTCL8
                //int * n_send_tmp = new int[n_proc];
                //int * rank_send_tmp = new int[n_loc];
                //for(int i=0; i<n_proc; i++) n_send_tmp[i] = 0;
                //for(S32 ip=0; ip<n_loc; ip++) {
                //    F64 pos_x = ptcl_[ip].getPos().x;
                //    F64 pos_y = ptcl_[ip].getPos().y;
                //    F64 pos_z = ptcl_[ip].getPos().z;
                //    F64 pos_phi = atan2(pos_y, pos_x);
                //    if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                //    F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                //    F64vec pos_tmp = F64vec(pos_phi, pos_r, pos_z);
                //    if(determineWhetherParticleIsInDomain(pos_tmp, pos_my_domain)) {
                //        n_send_tmp[my_rank]++;
                //        rank_send_tmp[ip] = my_rank;
                //    }
                //    else{
                //        S32 srank = searchWhichDomainParticleGoToPeriodicX(pos_tmp, n_domain, pos_domain, peri_len_x, previous_rank);
                //        previous_rank = srank;
                //        n_send_tmp[srank]++;
                //        rank_send_tmp[ip] = srank;
                //    }
                //}
                //for(S32 ip=0; ip<n_loc; ip++){
                //    assert(rank_send_tmp[ip] == rank_send[ip]);
                //}
                //delete [] n_send_tmp;
                //delete [] rank_send_tmp;
        #endif

                /*
                iloc = 0;
                for(S32 i=0; i<n_loc; i++){
                    int rank_tmp = rank_send[i];
                    if(rank_tmp == my_rank){
                        ptcl_[iloc++] = ptcl_[i];
                    }
                    else{
                        ptcl_send_buf_[rank_tmp].push_back(ptcl_[i]);
                        __n_send[rank_tmp]++;
                    }
                }
                __n_rem = iloc;
                */
            }
            n_disp_send[0] = 0;
            for(S32 i=0; i<n_proc; i++) {
                n_disp_send[i+1] = n_disp_send[i] + __n_send[i];
            }
            ptcl_.resizeNoInitialize(iloc);



        #ifdef DEBUG_PRINT_EX_PTCL8
            for(S32 ip=0; ip<iloc; ip++){
                F64 pos_x = ptcl_[ip].getPos().x;
                F64 pos_y = ptcl_[ip].getPos().y;
                F64 pos_z = ptcl_[ip].getPos().z;
                F64 pos_phi = atan2(pos_y, pos_x);
                if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                F64vec pos_tmp = F64vec(pos_phi, pos_r, pos_z);                
                assert(determineWhetherParticleIsInDomain(pos_tmp, pos_my_domain));
            }
            for(S32 ib=0; ib<n_proc; ib++){
                for(S32 ip=0; ip<__n_send[ib]; ip++){
                    F64 pos_x = ptcl_send_buf_[ib][ip].getPos().x;
                    F64 pos_y = ptcl_send_buf_[ib][ip].getPos().y;
                    F64 pos_z = ptcl_send_buf_[ib][ip].getPos().z;
                    F64 pos_phi = atan2(pos_y, pos_x);
                    if(pos_phi < 0.0) pos_phi += 2.0*4.0*atan(1.0);
                    F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                    F64vec pos_tmp = F64vec(pos_phi, pos_r, pos_z);                
                    assert( determineWhetherParticleIsInDomain(pos_tmp, pos_domain[ib]) );
                }
            }
        #endif
    #else //SUNWAY
            assert(false);
            /*
            previous_rank = my_rank;
            for(S32 ip=0; ip<n_loc; ip++) {
                if(determineWhetherParticleIsInDomain(ptcl_[ip].getPos(), pos_my_domain)) {
                    ptcl_[iloc] = ptcl_[ip];
                    previous_rank = my_rank;
                    N_NOT_MOVE++;
                    iloc++;
                }
                else {
        #if 1
                    S32 srank = searchWhichDomainParticleGoToPeriodicX(ptcl_[ip].getPos(), n_domain, pos_domain, peri_len_x, previous_rank);
        #else
                    S32 srank = searchWhichDomainParticleGoTo2(ptcl_[ip].getPos(), n_domain, pos_domain, previous_rank);
        #endif
                    ptcl_send_buf_[srank].push_back(ptcl_[ip]);
                    previous_rank = srank;
                    __n_send[srank]++;
                }
            }
            */
    #endif //SUNWAY
            
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 06 @exchangeParticle8"<<std::endl;
    #endif
            //time_profile_.exchange_particle__find_particle += GetWtime() - wtime_offset_inner;
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 07 @exchangeParticle8"<<std::endl;
    #endif
            // ****************************************************
            // *** receive the number of receive particles ********
            wtime_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 3"<<std::endl;
    #if 1
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 rank_x = dinfo.getRank1d(0);
            const S32 rank_y = dinfo.getRank1d(1);
            S32 dx_max_loc = 0;
            S32 dy_max_loc = 0;
            for(S32 i=0; i<n_proc; i++){
                if(__n_send[i] > 0){
                    S32 dx = abs( (i / n_proc_y) - rank_x);
                    if(dx > n_proc_x / 2) dx = n_domain[0] - dx;
                    dx_max_loc =  dx > dx_max_loc ? dx : dx_max_loc;
                    dy_max_loc = abs( (i % n_proc_y) - rank_y) > dy_max_loc ? abs( (i % n_proc_y) - rank_y) : dy_max_loc;
                }
            }
            S32 dx_max_glb = Comm::getMaxValue(dx_max_loc);
            S32 dy_max_glb = Comm::getMaxValue(dy_max_loc);
        #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); 
            if(Comm::getRank()==0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<std::endl;
                std::cerr<<"dy_max_glb= "<<dy_max_glb<<std::endl;
            }
        #endif
            assert(dx_max_glb*2+1 < n_domain[0]); // TO DO; have to consider this case.
            for(S32 ix=-dx_max_glb; ix<=dx_max_glb; ix++){
                S32 rank_x_target = rank_x+ix;
                if(rank_x_target < 0) rank_x_target = rank_x_target + n_proc_x;
                else if(rank_x_target >= n_proc_x) rank_x_target = rank_x_target - n_proc_x;
        #ifdef DEBUG_PRINT_EX_PTCL8
                if(rank_x_target < 0 || rank_x_target >= n_proc_x){
                    std::cerr<<"rank= "<<Comm::getRank()
                             <<" rank_x_target= "<<rank_x_target
                             <<" ix= "<<ix
                             <<" dx_max_glb= "<<dx_max_glb
                             <<" n_proc_x= "<<n_proc_x
                             <<" rank_x= "<<rank_x
                             <<std::endl;
                }
         #endif
                assert(rank_x_target >= 0 && rank_x_target < n_proc_x);
                for(S32 iy=-dy_max_glb; iy<=dy_max_glb; iy++){
                    S32 rank_y_target = rank_y+iy;
                    if(rank_y_target < 0 || rank_y_target >= n_proc_y) continue;
                    S32 rank = (rank_x_target)*n_proc_y + (rank_y_target);
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(__n_send+rank, 1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(n_recv+rank, 1, GetDataType<S32>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
        #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 08 @exchangeParticle8"<<std::endl;
        #endif
    #else
            Comm::allToAll(nsend, 1, nrecv);
    #endif
            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            ptcl_recv_.resizeNoInitialize( n_disp_recv[n_proc] );
    #ifdef DEBUG_PRINT_EX_PTCL8
            if(Comm::getRank() == 0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<" dy_max_glb= "<<dy_max_glb<<std::endl;
                for(S32 i=0; i<n_proc; i++){
                    std::cerr<<"__n_send[i]="<<__n_send[i]<<" n_recv[i]="<<n_recv[i]<<std::endl;
                }
            }
    #endif
            const S32 n_proc_comm_limit = n_proc+1;
            n_proc_send = n_proc_recv = 0;
            for(S32 i=0; i<n_proc; i++){
                if(__n_send[i] > 0) n_proc_send++;
                if(n_recv[i] > 0) n_proc_recv++;
            }
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 09 @exchangeParticle7"<<std::endl;
    #endif
            n_proc_send = n_proc_recv = 0;
            for(S32 ib = 1; ib < n_proc; ib++) {
                S32 idsend = (ib + my_rank) % n_proc;
                if(__n_send[idsend] > 0) {
                    S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_buf_[idsend].getPointer(), __n_send[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                }
                S32 idrecv = (n_proc + my_rank - ib) % n_proc;
                if(n_recv[idrecv] > 0) {
                    S32 adrrecv = n_disp_recv[idrecv];
                    S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), n_recv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 10 @exchangeParticle8"<<std::endl;
    #endif
            const S32 n_recv_tot = n_disp_recv[n_proc];
            ptcl_.reserveEmptyAreaAtLeast( n_recv_tot );
            for(S32 ip = 0; ip < n_recv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
            // **************************************************** 
            n_ptcl_send_ += n_disp_send[n_proc];
            n_ptcl_recv_ += n_disp_recv[n_proc];
            time_profile_.exchange_particle__exchange_particle = GetWtime() - wtime_offset_inner;
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 11 @exchangeParticle8"<<std::endl;
    #endif
            //delete [] rank_send;
            delete [] n_disp_send;
            delete [] n_recv;
            delete [] n_disp_recv;
            delete [] req_send;
            delete [] req_recv;
            time_profile_.exchange_particle = GetWtime() - time_offset;
        }
        
        template<class Tdinfo>
        void exchangeParticle8v2(Tdinfo & dinfo, const S32 op_mode, const S32 debug_flag) {
            F64 time_offset = GetWtime();
            const F64 peri_len_x = 2.0 * 4.0 * atan(1.0);
            const S32 n_loc   = ptcl_.size();
            const S32 my_rank = MPI::COMM_WORLD.Get_rank();
            const S32 n_proc  = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK1 @exchangeParticle8v2"<<std::endl;
#endif
            static S32 __n_rem = 0; // not move
            static bool first = true;
            static S32 * __n_send = NULL;
            if(first){
                __n_send  = new S32[n_proc];
                first = false;
            }
            //S32 * rank_send = new S32[n_loc];
            S32 * n_disp_send = new S32[n_proc+1];
            S32 * n_recv  = new S32[n_proc];
            S32 * n_disp_recv = new S32[n_proc+1];
            MPI::Request * req_send = new MPI::Request[n_proc];
            MPI::Request * req_recv = new MPI::Request[n_proc];
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK2 @exchangeParticle8v2"<<std::endl;
#endif
            F64 wtime_offset_inner = GetWtime();
            S32 previous_rank = my_rank;
            S32 iloc = 0;
#ifdef SUNWAY
            for(S32 i=0; i<n_proc; i++) {
                __n_send[i] = n_disp_send[i] = n_recv[i] = n_disp_recv[i] = 0;
                ptcl_send_buf_[i].clearSize();
            }
            n_disp_send[n_proc] = n_disp_recv[n_proc] = 0;
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK3 @exchangeParticle8v2"<<std::endl;
            Comm::barrier();
            const double start_time = MPI_Wtime();
#endif
            if (op_mode == 1) {
                //* [Original implementation]
                double * pos_phi = new double[n_loc];
                double * pos_r = new double[n_loc];
                unsigned long args[16];
                args[0] = (unsigned long)(n_loc);
                args[1] = (unsigned long) &peri_len_x;
                args[2] = (unsigned long) ptcl_.getPointer();
                args[3] = (unsigned long) pos_phi;
                args[4] = (unsigned long) pos_r;
                __real_athread_spawn((void*)slave_GetCylCoord, args);
                athread_join();

                iloc = 0;
                previous_rank = my_rank;
                for(S32 ip=0; ip<n_loc; ip++) {
                    F64 pos_z = ptcl_[ip].getPos().z;
                    F64vec pos_tmp = F64vec(pos_phi[ip], pos_r[ip], pos_z);
                    S32 srank = searchWhichDomainParticleGoToPeriodicX(pos_tmp, n_domain, pos_domain, peri_len_x, previous_rank);
                    if(srank == my_rank){
                        ptcl_[iloc++] = ptcl_[ip];
                    }
                    else{
                        ptcl_send_buf_[srank].push_back(ptcl_[ip]);
                        __n_send[srank]++;
                    }
                    previous_rank = srank;
                    //[DEBUG]
                    if ((my_rank == 22) && ((ip == 6091427) || (ip == 6091435))) {
                        printf("[org] ip = %d, id = %lu\n",ip,ptcl_[ip].id);
                    }
                }
                __n_rem = iloc;
                delete [] pos_phi;
                delete [] pos_r;
            } else {
                //* [New implementation]
                // Examine the distination rank for each particle and
                // count the number of particles to be sent for each rank.
                volatile S32 * dest_rank = new S32[n_loc];
                volatile S32 * n_dest_rank_cpe = new S32[NUMBER_OF_CPE];
                volatile S32 * dest_rank_list_cpe = new S32[NUMBER_OF_CPE * n_proc];
                volatile S32 * n_send_cpe = new S32[NUMBER_OF_CPE * n_proc];
                volatile S32 * n_disp_send_cpe = new S32[NUMBER_OF_CPE * n_proc];
                volatile U64 * adr_ptcl_send_buf_ = new U64[n_proc];
                unsigned long args[16];
                args[0]  = (unsigned long) my_rank;
                args[1]  = (unsigned long) n_proc;
                args[2]  = (unsigned long) n_domain;
                args[3]  = (unsigned long) n_loc;
                args[4]  = (unsigned long) ptcl_.getPointer();
                args[5]  = (unsigned long) pos_domain;
                args[6]  = (unsigned long) dest_rank;
                args[7]  = (unsigned long) n_dest_rank_cpe;
                args[8]  = (unsigned long) dest_rank_list_cpe;
                args[9]  = (unsigned long) n_send_cpe;
                args[10] = (unsigned long) debug_flag;
                __real_athread_spawn((void*)slave_FindDestinationRank, args);
                athread_join();
                // Compute __n_send[]
                for (S32 i=0; i<n_loc; i++) {
                   const S32 srank = dest_rank[i];
                   assert((0 <= srank) && (srank < n_proc));
                   __n_send[srank]++;
                }
                __n_send[my_rank]=0;
                // Compute the prefix sum (n_disp_send_cpe)
                //-(0th CPE)
                for (S32 i=0; i<n_dest_rank_cpe[0]; i++) n_disp_send_cpe[i] = 0;
                //-(the other CPEs) 
                for (S32 id_tgt=1; id_tgt<NUMBER_OF_CPE; id_tgt++) {
                    const S32 offset_tgt = id_tgt * n_proc;
                    for (S32 i=0; i<n_dest_rank_cpe[id_tgt]; i++) {
                        const S32 srank_tgt = dest_rank_list_cpe[offset_tgt + i];
                        n_disp_send_cpe[offset_tgt + i] = 0;
                        for (S32 id=0; id<id_tgt; id++) {
                            const S32 offset = id * n_proc;
                            for (S32 j=0; j<n_dest_rank_cpe[id]; j++) {
                                const S32 srank = dest_rank_list_cpe[offset + j];
                                if (srank == srank_tgt) {
                                    n_disp_send_cpe[offset_tgt + i] += n_send_cpe[offset + j];
                                }
                            }
                        }
                    }
                }
#if 0
                // [CHECK]
                Comm::barrier();
                if (Comm::getRank() == 0) 
                    std::cout << "outputing ex_ptcl_tbl..." << std::endl;
                std::ofstream ofs;
                std::stringstream filenum;
                std::string filename;
                filenum << std::setfill('0') << std::setw(5) << my_rank;
                filename = "./ex_ptcl_tbl" + filenum.str() + ".txt";
                ofs.open(filename.c_str(), std::ios::trunc);
                S32 n_send_tot = 0;
                for (S32 cpe_id=0; cpe_id<NUMBER_OF_CPE; cpe_id++) {
                    ofs << "##################################" << std::endl;
                    ofs << "cpe_id = " << cpe_id << " "
                        << "n_dest_rank = " << n_dest_rank_cpe[cpe_id] << " "
                        << std::endl;
                    ofs << "---- the content of the lists ----" << std::endl;
                    const S32 offset = cpe_id * n_proc;
                    for (S32 j=0; j<n_dest_rank_cpe[cpe_id]; j++) {
                        ofs << "dest_rank = " << dest_rank_list_cpe[offset + j] << " "
                            << "n_send = " << n_send_cpe[offset + j] << " "
                            << "n_disp_send = " << n_disp_send_cpe[offset + j] << " "
                            << std::endl;
                        n_send_tot += n_send_cpe[offset + j];
                    }
                }
                ofs << "(n_loc = " << n_loc << " n_send_tot = " << n_send_tot << ")" << std::endl;
                ofs.close();
                Comm::barrier();
                if (Comm::getRank() == 0) 
                    std::cout << "output of ex_ptcl_tbl is completed!" << std::endl;
#endif
                // Memory allocation for send buffers
                S32 nrank_send_loc = 0, nrank_send_max;
                for (S32 rank=0; rank<n_proc; rank++) {
                    if (rank != my_rank) {
                        if (__n_send[rank] > 0) {
                            ptcl_send_buf_[rank].reserve(__n_send[rank]);
                            nrank_send_loc++;
                        }
                        ptcl_send_buf_[rank].resizeNoInitialize(__n_send[rank]);
                        adr_ptcl_send_buf_[rank] = (unsigned long) ptcl_send_buf_[rank].getPointer();
                    } else {
                        adr_ptcl_send_buf_[rank] = (unsigned long) ptcl_.getPointer();
                    }
                }
                //MPI_Allreduce(&nrank_send_loc,&nrank_send_max, 
                //              1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
                //if (my_rank == 0) std::cout << "nrank_send_max = " << nrank_send_max << std::endl;
                // Pack particles into the send buffers
                args[0] = (unsigned long) my_rank;
                args[1] = (unsigned long) n_proc;
                args[2] = (unsigned long) n_loc;
                args[3] = (unsigned long) ptcl_.getPointer();
                args[4] = (unsigned long) dest_rank;
                args[5] = (unsigned long) n_dest_rank_cpe;
                args[6] = (unsigned long) dest_rank_list_cpe;
                args[7] = (unsigned long) n_disp_send_cpe;
                args[8] = (unsigned long) adr_ptcl_send_buf_; 
                args[9] = (unsigned long) debug_flag;
                __real_athread_spawn((void*)slave_MakeSendBuffers, args);
                athread_join();
#if 0
                // Shrink ptcl_
                iloc = 0; // Note that iloc is used in the later part of the code.
                for (S32 i=0; i<n_loc; i++) {
                    const S32 srank = dest_rank[i];
                    if (srank == my_rank) {
                        ptcl_[iloc++] = ptcl_[i];
                    }
                }
#endif
                // Release memory 
                delete [] dest_rank;
                delete [] n_dest_rank_cpe;
                delete [] dest_rank_list_cpe;
                delete [] n_send_cpe;
                delete [] n_disp_send_cpe;
                delete [] adr_ptcl_send_buf_;
            }
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            const double end_time = MPI_Wtime();
            if(Comm::getRank()==0)std::cerr<<"OK6 @exchangeParticle8v2"<<std::endl;
            if (Comm::getRank() == 0)
                std::cout << "wtime(srank) = " << end_time - start_time << " [s]" << std::endl;

#if 0
            // [CHECK]
            if (debug_flag) {
            Comm::barrier();
            if (Comm::getRank() == 0) 
                std::cout << "outputing ptcl_send_buf_..." << std::endl;
            //if ((1 <= Comm::getRank()) && (Comm::getRank() <= 3)) {// tmp
            //if (Comm::getRank() == 22) {// tmp
            std::ofstream ofs;
            std::stringstream filenum;
            std::string filename;
            filenum << std::setfill('0') << std::setw(5) << my_rank;
            //filename = "./ptcl_send_buf_-old-" + filenum.str() + ".dat";
            filename = "./ptcl_send_buf_-new-" + filenum.str() + ".dat";
            //ofs.open(filename.c_str(), std::ios::binary|std::ios::trunc);
            ofs.open(filename.c_str(), std::ios::trunc);
            for (S32 rank=0; rank<n_proc; rank++) 
                if (ptcl_send_buf_[rank].size() > 0) 
                    for (S32 i=0; i<ptcl_send_buf_[rank].size(); i++)
                        ofs << ptcl_send_buf_[rank][i].pos << " "
                            << ptcl_send_buf_[rank][i].mass << " "
                            << ptcl_send_buf_[rank][i].vel << " "
                            << ptcl_send_buf_[rank][i].id << " "
                            << std::endl;
            ofs.close();
            //}// tmp
            Comm::barrier();
            if (Comm::getRank() == 0) 
                std::cout << "output of ptcl_send_buf_ is completed!" << std::endl;
            }
#endif
#endif
#if 1
            //* For a test run
            Comm::barrier();
            if (Comm::getRank() == 0) std::cout << "OK!" << std::endl;
            athread_halt();
            Finalize();
            std::exit(0);
#endif
            time_profile_.exchange_particle__find_particle = GetWtime() - wtime_offset_inner;

            n_disp_send[0] = 0;
            for(S32 i=0; i<n_proc; i++) {
                n_disp_send[i+1] = n_disp_send[i] + __n_send[i];
            }
            ptcl_.resizeNoInitialize(iloc);


#else //SUNWAY
            assert(false);
#endif //SUNWAY
            
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK7 @exchangeParticle8v2"<<std::endl;
#endif
            //time_profile_.exchange_particle__find_particle += GetWtime() - wtime_offset_inner;
#ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK8 @exchangeParticle8v2"<<std::endl;
#endif
            // ****************************************************
            // *** receive the number of receive particles ********
            wtime_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"in exptcl 3"<<std::endl;
#if 1
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 rank_x = dinfo.getRank1d(0);
            const S32 rank_y = dinfo.getRank1d(1);
            S32 dx_max_loc = 0;
            S32 dy_max_loc = 0;
            for(S32 i=0; i<n_proc; i++){
                if(__n_send[i] > 0){
                    S32 dx = abs( (i / n_proc_y) - rank_x);
                    if(dx > n_proc_x / 2) dx = n_domain[0] - dx;
                    dx_max_loc =  dx > dx_max_loc ? dx : dx_max_loc;
                    dy_max_loc = abs( (i % n_proc_y) - rank_y) > dy_max_loc ? abs( (i % n_proc_y) - rank_y) : dy_max_loc;
                }
            }
            S32 dx_max_glb = Comm::getMaxValue(dx_max_loc);
            S32 dy_max_glb = Comm::getMaxValue(dy_max_loc);
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier(); 
            if(Comm::getRank()==0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<std::endl;
                std::cerr<<"dy_max_glb= "<<dy_max_glb<<std::endl;
            }
#endif
            assert(dx_max_glb*2+1 < n_domain[0]); // TO DO; have to consider this case.
            for(S32 ix=-dx_max_glb; ix<=dx_max_glb; ix++){
                S32 rank_x_target = rank_x+ix;
                if(rank_x_target < 0) rank_x_target = rank_x_target + n_proc_x;
                else if(rank_x_target >= n_proc_x) rank_x_target = rank_x_target - n_proc_x;
#ifdef DEBUG_PRINT_EX_PTCL8_V2
                if(rank_x_target < 0 || rank_x_target >= n_proc_x){
                    std::cerr<<"rank= "<<Comm::getRank()
                             <<" rank_x_target= "<<rank_x_target
                             <<" ix= "<<ix
                             <<" dx_max_glb= "<<dx_max_glb
                             <<" n_proc_x= "<<n_proc_x
                             <<" rank_x= "<<rank_x
                             <<std::endl;
                }
#endif
                assert(rank_x_target >= 0 && rank_x_target < n_proc_x);
                for(S32 iy=-dy_max_glb; iy<=dy_max_glb; iy++){
                    S32 rank_y_target = rank_y+iy;
                    if(rank_y_target < 0 || rank_y_target >= n_proc_y) continue;
                    S32 rank = (rank_x_target)*n_proc_y + (rank_y_target);
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(__n_send+rank, 1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(n_recv+rank, 1, GetDataType<S32>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
#ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK9 @exchangeParticle8v2"<<std::endl;
#endif
#else
            Comm::allToAll(nsend, 1, nrecv);
#endif
            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            ptcl_recv_.resizeNoInitialize( n_disp_recv[n_proc] );
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            if(Comm::getRank() == 0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<" dy_max_glb= "<<dy_max_glb<<std::endl;
                for(S32 i=0; i<n_proc; i++){
                    std::cerr<<"__n_send[i]="<<__n_send[i]<<" n_recv[i]="<<n_recv[i]<<std::endl;
                }
            }
#endif
            const S32 n_proc_comm_limit = n_proc+1;
            n_proc_send = n_proc_recv = 0;
            for(S32 i=0; i<n_proc; i++){
                if(__n_send[i] > 0) n_proc_send++;
                if(n_recv[i] > 0) n_proc_recv++;
            }
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK10 @exchangeParticle8v2"<<std::endl;
#endif
            n_proc_send = n_proc_recv = 0;
            for(S32 ib = 1; ib < n_proc; ib++) {
                S32 idsend = (ib + my_rank) % n_proc;
                if(__n_send[idsend] > 0) {
                    S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_buf_[idsend].getPointer(), __n_send[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                }
                S32 idrecv = (n_proc + my_rank - ib) % n_proc;
                if(n_recv[idrecv] > 0) {
                    S32 adrrecv = n_disp_recv[idrecv];
                    S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), n_recv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK11 @exchangeParticle8v2"<<std::endl;
#endif
            const S32 n_recv_tot = n_disp_recv[n_proc];
            ptcl_.reserveEmptyAreaAtLeast( n_recv_tot );
            for(S32 ip = 0; ip < n_recv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
            // **************************************************** 
            n_ptcl_send_ += n_disp_send[n_proc];
            n_ptcl_recv_ += n_disp_recv[n_proc];
            time_profile_.exchange_particle__exchange_particle = GetWtime() - wtime_offset_inner;
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK12 @exchangeParticle8v2"<<std::endl;
#endif
            //delete [] rank_send;
            delete [] n_disp_send;
            delete [] n_recv;
            delete [] n_disp_recv;
            delete [] req_send;
            delete [] req_recv;
            time_profile_.exchange_particle = GetWtime() - time_offset;
#if 0
            //* For a test run
            Comm::barrier();
            if (Comm::getRank() == 0) {
                std::cout << "exchangeParticle8v2() completed!" << std::endl;
                std::cout << "etime(exchPtcl8v2) = " << time_profile_.exchange_particle << " [s]" << std::endl;
            }
            athread_halt();
            Finalize();
            std::exit(0);
#endif
        }

        

        void dumpMemSizeUsed(std::ostream & fout){
            F64 norm_fac=1e-9;
            std::string str_unit=" [GB]";
            U64 size_of_ptcl_send_buf_ = 0;
            for (S32 rank=0; rank<Comm::getNumberOfProc(); rank++)
                size_of_ptcl_send_buf_ += ptcl_send_buf_[rank].getMemSize();
            F64 sum = (double)(ptcl_.getMemSize() 
                             + idx_remove_ptcl_.getMemSize()
                             + ptcl_send_.getMemSize()
                             + size_of_ptcl_send_buf_
                             + ptcl_recv_.getMemSize())
                             * norm_fac;
            F64 sum_max;
            MPI::COMM_WORLD.Allreduce(&sum,&sum_max,1,MPI::DOUBLE,MPI::MAX);

            if (Comm::getRank() == 0) {
                fout<<"*** ParticleSystem class ***" << std::endl;
                fout<<"ptcl_             = " << ptcl_.getMemSize()*norm_fac             << str_unit << std::endl;
                fout<<"    ptcl_.size()  = " << ptcl_.size()                                        << std::endl;
                fout<<"    sizeof(Tptcl) = " << sizeof(Tptcl)                           << " [B]"   << std::endl;
                fout<<"idx_remove_ptcl_  = " << idx_remove_ptcl_.getMemSize()*norm_fac  << str_unit << std::endl;
                fout<<"ptcl_send_        = " << ptcl_send_.getMemSize()*norm_fac        << str_unit << std::endl;
                fout<<"ptcl_send_buf_    = " << size_of_ptcl_send_buf_*norm_fac         << str_unit << std::endl;
                fout<<"ptcl_recv_        = " << ptcl_recv_.getMemSize()*norm_fac        << str_unit << std::endl;
                fout<<"sum               = " << sum                                     << str_unit << std::endl;
                fout<<"sum (psys,max.)   = " << sum_max                                 << str_unit << std::endl;
            }
        }

        void freeCommunicationBuffer(){
            ptcl_send_.freeMem();
            ptcl_recv_.freeMem();
        }

        void allGatherParticle(Tptcl *& ptcl_glb,
                               S64 & n_ptcl_glb) const {
            S32 n_proc = Comm::getNumberOfProc();
            S32 * n_loc_ar = new S32[n_proc];
            S32 * n_loc_disp = new S32[n_proc+1];
            S32 n_loc = this->getNumberOfParticleLocal();
            Comm::allGather(&n_loc, 1, n_loc_ar);
            n_loc_disp[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_loc_disp[i+1] = n_loc_disp[i] + n_loc_ar[i];
            }
            n_ptcl_glb = n_loc_disp[n_proc];
            ptcl_glb = new Tptcl[n_ptcl_glb];
            Comm::allGatherV(this->ptcl_.getPointer(), n_loc,
                             ptcl_glb, n_loc_ar, n_loc_disp);
            delete [] n_loc_ar;
            delete [] n_loc_disp;
        }
        
    };
}
