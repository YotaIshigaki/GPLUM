#pragma once

#include<cassert>
#include<fstream>

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
    
    bool IsInBox1d(const F64vec & pos, const F64ort & domain, const S32 cid){
            return (pos[cid] < domain.high_[cid]) && (pos[cid] >= domain.low_[cid]);
        }
    
    
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
        //ReallocatableArray<Tptcl> * ptcl_send_buf_;
        TempArray<Tptcl> * ptcl_send_buf_;
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


        
        S32 searchWhichDomainParticleGoToPeriodicX2(const F64vec & pos,
                                                    const S32 n_domain [],
                                                    const F64ort domain [],
                                                    const S32 previous_rank=0){
            S32 n_proc = n_domain[0] * n_domain[1] * n_domain[2]; 
            S32 idomain = previous_rank;
            if(determineWhetherParticleIsInDomain(pos, domain[idomain])){
                return idomain;
            }
            S32 shift = n_domain[1] * n_domain[2];
            // x direction
            for(S32 i=0; i<n_domain[0]; i++){
                S32 idomain_tmp = ( (idomain+i*shift) < n_proc) ? (idomain+i*shift) : (idomain+i*shift) - n_proc;
                if(IsInBox1d(pos, domain[idomain_tmp], 0)){
                    idomain = idomain_tmp;
                    break;
                }
                idomain_tmp = ( (idomain-i*shift) >= 0) ? (idomain-i*shift) : n_proc + (idomain-i*shift);
                if(IsInBox1d(pos, domain[idomain_tmp], 0)){
                    idomain = idomain_tmp;
                    break;
                }
            }
            // y direction
            shift = n_domain[2];
            if(domain[idomain].high_.y <= pos.y){
                do{
                    idomain += shift;
                }while(domain[idomain].high_.y <= pos.y);
            }
            else if(domain[idomain].low_.y > pos.y){
                do{
                    idomain -= shift;
                }while(domain[idomain].low_.y > pos.y);
            }
            // z direction
            if(domain[idomain].high_.z <= pos.z){
                do{
                    idomain++;
                }while(domain[idomain].high_.z <= pos.z);
            }
            else if(domain[idomain].low_.z > pos.z){
                do{
                    idomain--;
                }while(domain[idomain].low_.z > pos.z);
            }
            return idomain;
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
            //ptcl_send_buf_ = new ReallocatableArray<Tptcl>[Comm::getNumberOfProc()];
            ptcl_send_buf_ = new TempArray<Tptcl>[Comm::getNumberOfProc()];
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
        void calcDxDyMax(const Tdinfo & dinfo, const S32 n_send[], const S32 n_proc, S32 & dx_max, S32 & dy_max){
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 my_rank_x = dinfo.getRank1d(0);
            const S32 my_rank_y = dinfo.getRank1d(1);
            S32 dx_max_loc = 0;
            S32 dy_max_loc = 0;
            for(S32 i=0; i<n_proc; i++) {
                if(n_send[i] > 0){
                    S32 dx = abs( (i / n_proc_y) - my_rank_x);
                    if(dx > n_proc_x / 2) dx = n_proc_x - dx;
                    S32 dy = abs( (i % n_proc_y) - my_rank_y);
                    dx_max_loc =  dx > dx_max_loc ? dx : dx_max_loc;
                    dy_max_loc =  dy > dy_max_loc ? dy : dy_max_loc;
                }
            }
            dx_max = Comm::getMaxValue(dx_max_loc);
            dy_max = Comm::getMaxValue(dy_max_loc);
        }

        template<class Tdinfo>
        void exchangeNumberOfPtcl(const Tdinfo & dinfo,
                                  const S32 dx_max_glb, const S32 dy_max_glb,
                                  S32 * n_send, S32 * n_recv,
                                  MPI::Request * req_send, MPI::Request * req_recv){
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 my_rank_x = dinfo.getRank1d(0);
            const S32 my_rank_y = dinfo.getRank1d(1);
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            for(S32 ix=-dx_max_glb; ix<=dx_max_glb; ix++){
                S32 rank_x_target = my_rank_x+ix;
                if(rank_x_target < 0) rank_x_target = rank_x_target + n_proc_x;
                else if(rank_x_target >= n_proc_x) rank_x_target = rank_x_target - n_proc_x;
#if 0
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
                    S32 rank_y_target = my_rank_y+iy;
                    if(rank_y_target < 0 || rank_y_target >= n_proc_y) continue;
                    S32 rank = (rank_x_target)*n_proc_y + (rank_y_target);
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(n_send+rank, 1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(n_recv+rank, 1, GetDataType<S32>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
        }
        
        template<class Tdinfo>
        void exchangeParticle8(Tdinfo & dinfo,
                               bool flag_reuse=false) {
            F64 time_offset = GetWtime();
            const F64 peri_len_x = 2.0 * 4.0 * atan(1.0);
            const S32 n_loc   = ptcl_.size();
            const S32 my_rank = MPI::COMM_WORLD.Get_rank();
            const S32 n_proc  = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 00 @exchangeParticle8"<<std::endl;
    #endif
            TempArray<S32> n_send;
            TempArray<S32> n_disp_send;
            TempArray<S32> n_recv;
            TempArray<S32> n_disp_recv;
            TempArray<MPI::Request> req_send;
            TempArray<MPI::Request> req_recv;
            TempArray<double> pos_phi;
            TempArray<double> pos_r;
            n_send.resizeNoInitialize(n_proc);
            n_recv.resizeNoInitialize(n_proc);
            n_disp_send.resizeNoInitialize(n_proc+1);
            n_disp_recv.resizeNoInitialize(n_proc+1);
            req_send.resizeNoInitialize(n_proc);
            req_recv.resizeNoInitialize(n_proc);
            pos_phi.resizeNoInitialize(n_loc);
            pos_r.resizeNoInitialize(n_loc);
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 01 @exchangeParticle8"<<std::endl;
    #endif
            F64 wtime_offset_inner = GetWtime();
            S32 previous_rank = my_rank;
            S32 i_loc = 0;

    #ifdef SUNWAY
            if(flag_reuse){
                assert(false);
            }
            else{
                for(S32 i=0; i<n_proc; i++) {
                    n_send[i] = n_disp_send[i] = n_recv[i] = n_disp_recv[i] = 0;
                }
                n_disp_send[n_proc] = n_disp_recv[n_proc] = 0;
                #if 0
                unsigned long args[5];
                args[0] = (unsigned long)(n_loc);
                args[1] = (unsigned long) &peri_len_x;
                args[2] = (unsigned long) ptcl_.getPointer();
                args[3] = (unsigned long) pos_phi.getPointer();
                args[4] = (unsigned long) pos_r.getPointer();
                __real_athread_spawn((void*)slave_GetCylCoord, args);
                athread_join();
                #else
                for(S32 i=0; i<n_loc; i++) {
                    F64 pos_x = ptcl_[i].pos.x;
                    F64 pos_y = ptcl_[i].pos.y;
                    F64 pos_z = ptcl_[i].pos.z;
                    pos_r[i]   = sqrt(pos_x*pos_x + pos_y*pos_y);
                    pos_phi[i] = atan2(pos_y, pos_x);
                    if(pos_phi[i] < 0.0) pos_phi[i] += 8.0*atan(1.0);
                }
                #endif
                i_loc = 0;
                previous_rank = my_rank;
                for(S32 ip=0; ip<n_loc; ip++) {
                    F64 pos_z = ptcl_[ip].getPos().z;
                    F64vec pos_tmp = F64vec(pos_phi[ip], pos_r[ip], pos_z);
                    S32 srank = searchWhichDomainParticleGoToPeriodicX2(pos_tmp, n_domain, pos_domain, previous_rank);
                    if(srank != my_rank){
                        // just count
                        n_send[srank]++;
                    }
                    previous_rank = srank;
                }

                for(S32 ib=0; ib<n_proc; ib++){
                    if(ib==my_rank) continue;
                    ptcl_send_buf_[ib].reserve(n_send[ib]);
                }
                previous_rank = my_rank;
                for(S32 ip=0; ip<n_loc; ip++) {
                    F64 pos_z = ptcl_[ip].getPos().z;
                    F64vec pos_tmp = F64vec(pos_phi[ip], pos_r[ip], pos_z);
                    S32 srank = searchWhichDomainParticleGoToPeriodicX2(pos_tmp, n_domain, pos_domain, previous_rank);
                    if(srank == my_rank){
                        ptcl_[i_loc++] = ptcl_[ip];
                    }
                    else{
                        ptcl_send_buf_[srank].pushBackNoCheck(ptcl_[ip]);
                    }
                    previous_rank = srank;
                }
                time_profile_.exchange_particle__find_particle = GetWtime() - wtime_offset_inner;
            }
            n_disp_send[0] = 0;
            for(S32 i=0; i<n_proc; i++) {
                n_disp_send[i+1] = n_disp_send[i] + n_send[i];
            }
    #else //SUNWAY
            assert(false);
    #endif //SUNWAY
            
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 02 @exchangeParticle8"<<std::endl;
    #endif
            //time_profile_.exchange_particle__find_particle += GetWtime() - wtime_offset_inner;
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 03 @exchangeParticle8"<<std::endl;
    #endif
            // ****************************************************
            // *** receive the number of receive particles ********
            wtime_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            S32 dx_max_glb, dy_max_glb;
            calcDxDyMax(dinfo, n_send.getPointer(), n_proc, dx_max_glb, dy_max_glb);
        #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<std::endl;
                std::cerr<<"dy_max_glb= "<<dy_max_glb<<std::endl;
            }
        #endif
            assert(dx_max_glb*2+1 < n_domain[0]); // TO DO; have to consider this case.
            exchangeNumberOfPtcl(dinfo, dx_max_glb, dy_max_glb,
                                 n_send.getPointer(), n_recv.getPointer(),
                                 req_send.getPointer(), req_recv.getPointer());
        #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 04 @exchangeParticle8"<<std::endl;
        #endif

            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            ptcl_.resizeNoInitialize( i_loc + n_disp_recv[n_proc] );
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
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 05 @exchangeParticle8"<<std::endl;
    #endif
            n_proc_send = n_proc_recv = 0;
            for(S32 ib = 1; ib < n_proc; ib++) {
                S32 idsend = (ib + my_rank) % n_proc;
                if(n_send[idsend] > 0) {
                    S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_buf_[idsend].getPointer(), n_send[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                }
                S32 idrecv = (n_proc + my_rank - ib) % n_proc;
                if(n_recv[idrecv] > 0) {
                    S32 adrrecv = n_disp_recv[idrecv];
                    S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_.getPointer(i_loc+adrrecv), n_recv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send.getPointer());
            MPI::Request::Waitall(n_proc_recv, req_recv.getPointer());
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 06 @exchangeParticle8"<<std::endl;
    #endif
            // **************************************************** 
            n_ptcl_send_ += n_disp_send[n_proc];
            n_ptcl_recv_ += n_disp_recv[n_proc];
            time_profile_.exchange_particle__exchange_particle = GetWtime() - wtime_offset_inner;
    #ifdef DEBUG_PRINT_EX_PTCL8
            Comm::barrier(); if(Comm::getRank()==0)std::cerr<<"OK 07 @exchangeParticle8"<<std::endl;
    #endif
            for(S32 i=0; i<n_proc; i++) {
                ptcl_send_buf_[i].free();
            }
            pos_r.free();
            pos_phi.free();
            n_send.free();
            n_disp_send.free();
            n_recv.free();
            n_disp_recv.free();
            req_send.free();
            req_recv.free();
            time_profile_.exchange_particle = GetWtime() - time_offset;
        }

        
        template<class Tdinfo>
        void exchangeParticle8v2(Tdinfo & dinfo, const S32 op_mode=2, const S32 debug_flag=0) {
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
                //ptcl_send_buf_[i].clearSize();
            }
            n_disp_send[n_proc] = n_disp_recv[n_proc] = 0;
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK3 @exchangeParticle8v2"<<std::endl;
            Comm::barrier();
            const double start_time = MPI_Wtime();
#endif
            if (op_mode == 1) {
                assert(0);
                /*
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
                }
                __n_rem = iloc;
                delete [] pos_phi;
                delete [] pos_r;
                */
            } else {
                //* [New implementation]
                // Examine the distination rank for each particle and
                // count the number of particles to be sent for each rank.
                volatile S32 * dest_rank = new S32[n_loc];
                volatile S32 * n_dest_rank_cpe = new S32[NUMBER_OF_CPE];
                volatile S32 * dest_rank_list_cpe = new S32[NUMBER_OF_CPE * n_proc];
                volatile U64 * adr_ptcl_send_buf_ = new U64[n_proc];
                unsigned long args[16];
                args[0] = (unsigned long) my_rank;
                args[1] = (unsigned long) n_proc;
                args[2] = (unsigned long) n_domain;
                args[3] = (unsigned long) n_loc;
                args[4] = (unsigned long) ptcl_.getPointer();
                args[5] = (unsigned long) pos_domain;
                args[6] = (unsigned long) dest_rank;
                args[7] = (unsigned long) n_dest_rank_cpe;
                args[8] = (unsigned long) dest_rank_list_cpe;
                args[9] = (unsigned long) debug_flag;
                __real_athread_spawn((void*)slave_FindDestinationRank, args);
                athread_join();
                // Calculate __n_send[]
                for (S32 i=0; i<n_loc; i++) {
                   const S32 srank = dest_rank[i];
                   assert((0 <= srank) && (srank < n_proc));
                   __n_send[srank]++;
                }
                iloc = __n_send[my_rank];
                __n_send[my_rank]=0;
                // Make dest_rank_list
                int n_dest_rank = 0;
                std::vector<int> dest_rank_list;
                dest_rank_list.reserve(2048);
                for (S32 id=0; id<NUMBER_OF_CPE; id++) {
                    const S32 offset = id * n_proc;
                    for (S32 i=0; i<n_dest_rank_cpe[id]; i++) {
                        const S32 srank = dest_rank_list_cpe[offset + i];
                        if (n_dest_rank > 0) {
                           S32 is_exist = 0;
                           for (S32 j=0; j<n_dest_rank; j++) {
                               if (srank == dest_rank_list[j]) {
                                  is_exist = 1;
                                  break;
                               }
                           }
                           if (!is_exist) {
                               dest_rank_list.push_back(srank);
                               n_dest_rank++;
                           }
                        } else {
                            dest_rank_list.push_back(srank);
                            n_dest_rank++;
                        }
                    }
                }
                // Memory allocation for send buffers
                for (S32 rank=0; rank<n_proc; rank++) {
                    if (rank != my_rank) {
                        /*
                        if (__n_send[rank] > 0) {
                            ptcl_send_buf_[rank].reserve(__n_send[rank]);
                        }
                        ptcl_send_buf_[rank].resizeNoInitialize(__n_send[rank]);
                        */
                        ptcl_send_buf_[rank].resizeNoInitialize(__n_send[rank]);
                        adr_ptcl_send_buf_[rank] = (unsigned long) ptcl_send_buf_[rank].getPointer();
                    } else {
                        adr_ptcl_send_buf_[rank] = (unsigned long) ptcl_.getPointer();
                    }
                }
                // Pack particles into the send buffers
                args[0] = (unsigned long) my_rank;
                args[1] = (unsigned long) n_proc;
                args[2] = (unsigned long) n_loc;
                args[3] = (unsigned long) ptcl_.getPointer();
                args[4] = (unsigned long) dest_rank;
                args[5] = (unsigned long) n_dest_rank;
                args[6] = (unsigned long) &dest_rank_list[0];
                args[7] = (unsigned long) adr_ptcl_send_buf_; 
                args[8] = (unsigned long) debug_flag;
                __real_athread_spawn((void*)slave_MakeSendBuffers, args);
                athread_join();
                // Release memory 
                delete [] dest_rank;
                delete [] n_dest_rank_cpe;
                delete [] dest_rank_list_cpe;
                delete [] adr_ptcl_send_buf_;
            }
#ifdef DEBUG_PRINT_EX_PTCL8_V2
            Comm::barrier();
            const double end_time = MPI_Wtime();
            if(Comm::getRank()==0)std::cerr<<"OK6 @exchangeParticle8v2"<<std::endl;
            if (Comm::getRank() == 0)
                std::cout << "wtime(srank) = " << end_time - start_time << " [s]" << std::endl;

#endif
            time_profile_.exchange_particle__find_particle = GetWtime() - wtime_offset_inner;

            n_disp_send[0] = 0;
            for(S32 i=0; i<n_proc; i++) {
                n_disp_send[i+1] = n_disp_send[i] + __n_send[i];
            }
            ptcl_.resizeNoInitialize(iloc);

#if 0
            // [CHECK]
            //if (debug_flag) {
            Comm::barrier();
            if (Comm::getRank() == 0) 
                std::cout << "outputing ptcl_send_buf_..." << std::endl;
            if (my_rank == 4) {
            std::ofstream ofs;
            std::stringstream filenum;
            std::string filename;
            filenum << std::setfill('0') << std::setw(5) << my_rank;
            //filename = "./ptcl_send_buf_-old-" + filenum.str() + ".dat";
            filename = "./ptcl_send_buf_-new-" + filenum.str() + ".dat";
            //ofs.open(filename.c_str(), std::ios::binary|std::ios::trunc);
            ofs.open(filename.c_str(), std::ios::trunc);
#if 0
            for (S32 rank=0; rank<n_proc; rank++) {
                if (ptcl_send_buf_[rank].size() > 0) {
                    std::cout << "rank = " << rank << " size = " << ptcl_send_buf_[rank].size() << std::endl;
                    for (S32 i=0; i<ptcl_send_buf_[rank].size(); i++) {
                        ofs << ptcl_send_buf_[rank][i].pos.x << " "
                            << ptcl_send_buf_[rank][i].pos.y << " "
                            << ptcl_send_buf_[rank][i].pos.z << " "
                            << ptcl_send_buf_[rank][i].mass << " "
                            << ptcl_send_buf_[rank][i].vel.x << " "
                            << ptcl_send_buf_[rank][i].vel.y << " "
                            << ptcl_send_buf_[rank][i].vel.z << " "
                            << ptcl_send_buf_[rank][i].id << " "
                            << std::endl;
                    }
                }
            }
#else
            std::cout << "ptcl_.size() = " << ptcl_.size() << std::endl;
            for (S32 i=0; i<ptcl_.size(); i++)
                ofs << ptcl_[i].pos << " "
                    << ptcl_[i].mass << " "
                    << ptcl_[i].vel << " "
                    << ptcl_[i].id << " "
                    << std::endl;
#endif
            ofs.close();
            }
            Comm::barrier();
            if (Comm::getRank() == 0) 
                std::cout << "output of ptcl_send_buf_ is completed!" << std::endl;
#endif
#if 0
            //* For a test run
            if (debug_flag) {
            Comm::barrier();
            if (Comm::getRank() == 0) std::cout << "OK!" << std::endl;
            athread_halt();
            Finalize();
            std::exit(0);
            }
#endif

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

            for(S32 i=0; i<n_proc; i++) {
                ptcl_send_buf_[i].free();
                if(i==0){
                    std::cerr<<"CHECK A"<<std::endl;
                }
                assert(ptcl_send_buf_[i].data_ == NULL);
            }
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
        

        template<class Tdinfo>
        void exchangeParticle8v3(Tdinfo & dinfo, const S32 op_mode=2, const S32 debug_flag=0) {
            F64 time_offset = GetWtime();
            const F64 peri_len_x = 2.0 * 4.0 * atan(1.0);
            const S32 n_loc   = ptcl_.size();
            const S32 my_rank = MPI::COMM_WORLD.Get_rank();
            const S32 n_proc  = MPI::COMM_WORLD.Get_size();
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK1 @exchangeParticle8v3"<<std::endl;
#endif
            TempArray<S32> n_send;
            TempArray<S32> n_disp_send;
            TempArray<S32> n_recv;
            TempArray<S32> n_disp_recv;
            TempArray<MPI::Request> req_send;
            TempArray<MPI::Request> req_recv;
            n_send.resizeNoInitialize(n_proc);
            n_recv.resizeNoInitialize(n_proc);
            n_disp_send.resizeNoInitialize(n_proc+1);
            n_disp_recv.resizeNoInitialize(n_proc+1);
            req_send.resizeNoInitialize(n_proc);
            req_recv.resizeNoInitialize(n_proc);
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK2 @exchangeParticle8v3"<<std::endl;
#endif
            F64 wtime_offset_inner = GetWtime();
            S32 iloc = 0;
#ifdef SUNWAY
            for(S32 i=0; i<n_proc; i++) {
                n_send[i] = n_disp_send[i] = n_recv[i] = n_disp_recv[i] = 0;
                //ptcl_send_buf_[i].clearSize();
            }
            n_disp_send[n_proc] = n_disp_recv[n_proc] = 0;
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK3 @exchangeParticle8v3"<<std::endl;
            Comm::barrier();
            const double start_time = MPI_Wtime();
#endif
            if (op_mode == 1) {
                assert(0);
                /*
                //* [Original implementation]
                S32 previous_rank = my_rank;
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
                        n_send[srank]++;
                    }
                    previous_rank = srank;
                }
                delete [] pos_phi;
                delete [] pos_r;
                */
            } else {
                //* [New implementation]
                // Examine the distination rank for each particle and
                // count the number of particles to be sent for each rank.
                TempArray<S32> dest_rank;
                TempArray<S32> n_dest_rank_cpe;
                TempArray<S32> dest_rank_list_cpe;
                TempArray<U64> adr_ptcl_send_buf_;
                dest_rank.resizeNoInitialize(n_loc);
                n_dest_rank_cpe.resizeNoInitialize(NUMBER_OF_CPE);
                dest_rank_list_cpe.resizeNoInitialize(NUMBER_OF_CPE * n_proc);
                adr_ptcl_send_buf_.resizeNoInitialize(n_proc);
                unsigned long args[16];
                args[0] = (unsigned long) my_rank;
                args[1] = (unsigned long) n_proc;
                args[2] = (unsigned long) n_domain;
                args[3] = (unsigned long) n_loc;
                args[4] = (unsigned long) ptcl_.getPointer();
                args[5] = (unsigned long) pos_domain;
                args[6] = (unsigned long) dest_rank.getPointer();
                args[7] = (unsigned long) n_dest_rank_cpe.getPointer();
                args[8] = (unsigned long) dest_rank_list_cpe.getPointer();
                args[9] = (unsigned long) debug_flag;
                __real_athread_spawn((void*)slave_FindDestinationRank, args);
                athread_join();
                // Calculate __n_send[]
                for (S32 i=0; i<n_loc; i++) {
                   const S32 srank = dest_rank[i];
#if 0
                   if( (0 > srank) || (srank >= n_proc) ){
                       std::cerr<<"my_rank= "<<Comm::getRank()
                                <<" i= "<<i
                                <<" ptcl_[i].id= "<<ptcl_[i].id
                                <<" ptcl_[i].pos= "<<ptcl_[i].pos
                                <<std::endl;
                   }
#endif
                   assert((0 <= srank) && (srank < n_proc));
                   n_send[srank]++;
                }
                iloc = n_send[my_rank];
                n_send[my_rank]=0;
                // Make dest_rank_list
                int n_dest_rank = 0;
                std::vector<int> dest_rank_list;
                dest_rank_list.reserve(2048);
                for (S32 id=0; id<NUMBER_OF_CPE; id++) {
                    const S32 offset = id * n_proc;
                    for (S32 i=0; i<n_dest_rank_cpe[id]; i++) {
                        const S32 srank = dest_rank_list_cpe[offset + i];
                        if (n_dest_rank > 0) {
                           S32 is_exist = 0;
                           for (S32 j=0; j<n_dest_rank; j++) {
                               if (srank == dest_rank_list[j]) {
                                  is_exist = 1;
                                  break;
                               }
                           }
                           if (!is_exist) {
                               dest_rank_list.push_back(srank);
                               n_dest_rank++;
                           }
                        } else {
                            dest_rank_list.push_back(srank);
                            n_dest_rank++;
                        }
                    }
                }
                // Memory allocation for send buffers
                for (S32 rank=0; rank<n_proc; rank++) {
                    if (rank != my_rank) {
                        ptcl_send_buf_[rank].resizeNoInitialize(n_send[rank]);
                        adr_ptcl_send_buf_[rank] = (unsigned long) ptcl_send_buf_[rank].getPointer();
                    } else {
                        adr_ptcl_send_buf_[rank] = (unsigned long) ptcl_.getPointer();
                    }
                }
                // Pack particles into the send buffers
                args[0] = (unsigned long) my_rank;
                args[1] = (unsigned long) n_proc;
                args[2] = (unsigned long) n_loc;
                args[3] = (unsigned long) ptcl_.getPointer();
                args[4] = (unsigned long) dest_rank.getPointer();
                args[5] = (unsigned long) n_dest_rank;
                args[6] = (unsigned long) &dest_rank_list[0];
                args[7] = (unsigned long) adr_ptcl_send_buf_.getPointer(); 
                args[8] = (unsigned long) debug_flag;
                __real_athread_spawn((void*)slave_MakeSendBuffers, args);
                athread_join();
                // Release memory 
                dest_rank.free();
                n_dest_rank_cpe.free();
                dest_rank_list_cpe.free();
                adr_ptcl_send_buf_.free();
            }
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            const double end_time = MPI_Wtime();
            if(Comm::getRank()==0)std::cerr<<"OK6 @exchangeParticle8v3"<<std::endl;
            if (Comm::getRank() == 0)
                std::cout << "wtime(srank) = " << end_time - start_time << " [s]" << std::endl;

#endif
            time_profile_.exchange_particle__find_particle = GetWtime() - wtime_offset_inner;

            n_disp_send[0] = 0;
            for(S32 i=0; i<n_proc; i++) {
                n_disp_send[i+1] = n_disp_send[i] + n_send[i];
            }

#else //SUNWAY
            assert(false);
#endif //SUNWAY
            
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK7 @exchangeParticle8v3"<<std::endl;
#endif
            //time_profile_.exchange_particle__find_particle += GetWtime() - wtime_offset_inner;

            // ****************************************************
            // *** receive the number of receive particles ********
            wtime_offset_inner = GetWtime();
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            S32 dx_max_glb, dy_max_glb;
            calcDxDyMax(dinfo, n_send.getPointer(), n_proc, dx_max_glb, dy_max_glb);
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"dx_max_glb= "<<dx_max_glb<<std::endl;
                std::cerr<<"dy_max_glb= "<<dy_max_glb<<std::endl;
            }
#endif
            assert(dx_max_glb*2+1 < n_domain[0]); // TO DO; have to consider this case.
            exchangeNumberOfPtcl(dinfo, dx_max_glb, dy_max_glb,
                                 n_send.getPointer(), n_recv.getPointer(),
                                 req_send.getPointer(), req_recv.getPointer());
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK 04 @exchangeParticle8v3"<<std::endl;
#endif

            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            ptcl_.resizeNoInitialize( iloc + n_disp_recv[n_proc] );
#ifdef DEBUG_PRINT_EX_PTCL8_V3
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
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK 05 @exchangeParticle8v3"<<std::endl;
#endif
            n_proc_send = n_proc_recv = 0;
            for(S32 ib = 1; ib < n_proc; ib++) {
                S32 idsend = (ib + my_rank) % n_proc;
                if(n_send[idsend] > 0) {
                    S32 tagsend = (my_rank < idsend) ? my_rank : idsend;
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_buf_[idsend].getPointer(), n_send[idsend], GetDataType<Tptcl>(), idsend, tagsend);
                }
                S32 idrecv = (n_proc + my_rank - ib) % n_proc;
                if(n_recv[idrecv] > 0) {
                    S32 adrrecv = n_disp_recv[idrecv];
                    S32 tagrecv = (my_rank < idrecv) ? my_rank : idrecv;
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_.getPointer(iloc+adrrecv), n_recv[idrecv], GetDataType<Tptcl>(), idrecv, tagrecv);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send.getPointer());
            MPI::Request::Waitall(n_proc_recv, req_recv.getPointer());
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK 06 @exchangeParticle8v3"<<std::endl;
#endif
            // **************************************************** 
            n_ptcl_send_ += n_disp_send[n_proc];
            n_ptcl_recv_ += n_disp_recv[n_proc];
            time_profile_.exchange_particle__exchange_particle = GetWtime() - wtime_offset_inner;
#ifdef DEBUG_PRINT_EX_PTCL8_V3
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK12 @exchangeParticle8v3"<<std::endl;
#endif

            for(S32 i=0; i<n_proc; i++) {
                ptcl_send_buf_[i].free();
                assert(ptcl_send_buf_[i].data_ == NULL);
            }
            n_send.free();
            n_disp_send.free();
            n_recv.free();
            n_disp_recv.free();
            req_send.free();
            req_recv.free();
            time_profile_.exchange_particle = GetWtime() - time_offset;
#if 0
            // [CHECK]
            Comm::barrier();
            if (Comm::getRank() == 0) 
                std::cout << "outputing ptcl_ after exchangeParticle..." << std::endl;
            if (my_rank == 4) {
            std::ofstream ofs;
            std::stringstream filenum;
            std::string filename;
            filenum << std::setfill('0') << std::setw(5) << my_rank;
            filename = "./ptcl_-v3-" + filenum.str() + ".dat";
            //ofs.open(filename.c_str(), std::ios::binary|std::ios::trunc);
            ofs.open(filename.c_str(), std::ios::trunc);
            std::cout << "ptcl_.size() = " << ptcl_.size() << std::endl;
            for (S32 i=0; i<ptcl_.size(); i++)
                ofs << ptcl_[i].pos << " "
                    << ptcl_[i].mass << " "
                    << ptcl_[i].vel << " "
                    << ptcl_[i].id << " "
                    << std::endl;
            ofs.close();
            }
            Comm::barrier();
            if (Comm::getRank() == 0) 
                std::cout << "output of ptcl_ is completed!" << std::endl;
#endif
#if 0
            //* For a test run
            Comm::barrier();
            if (Comm::getRank() == 0) std::cout << "OK!" << std::endl;
            athread_halt();
            Finalize();
            std::exit(0);
#endif
        }
        
        template<class Tdinfo>
        void exchangeParticleSw(Tdinfo & dinfo,
                                   bool flag_reuse_prev_result=false) {
            F64 time_offset = GetWtime();
            const F64 peri_len_x = 2.0 * 4.0 * atan(1.0);
            const S32 n_loc   = ptcl_.size();
            const S32 my_rank = MPI::COMM_WORLD.Get_rank();
            const S32 my_rank_x = dinfo.getRank1d(0);
            const S32 my_rank_y = dinfo.getRank1d(1);
            const S32 n_proc  = MPI::COMM_WORLD.Get_size();
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
            S32 * n_send = new S32[n_proc];
            S32 * rank_send_of_ptcl = new S32[n_loc];
            S32 n_send_tot = 0;
            S32 iloc = 0;
            F64 wtime_offset_inner = GetWtime();
            if(flag_reuse_prev_result){
                // under construction
            }
            else{
                //////////////////////
                // search send rank //
                // calc send rank on MPE (for debug)
                S32 previous_rank = my_rank;
                for(S32 ip=0; ip<n_loc; ip++) {
                    F64 pos_x = ptcl_[ip].getPos().x;
                    F64 pos_y = ptcl_[ip].getPos().y;
                    F64 pos_z = ptcl_[ip].getPos().z;
                    F64 pos_phi = atan2(pos_y, pos_x);
                    if(pos_phi < 0.0) pos_phi += peri_len_x;
                    F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                    F64vec pos_tmp = F64vec(pos_phi, pos_r, pos_z);
                    if(determineWhetherParticleIsInDomain(pos_tmp, pos_my_domain)) {
                        rank_send_of_ptcl[ip] = my_rank;
                    }
                    else{
                        S32 srank = searchWhichDomainParticleGoToPeriodicX(pos_tmp, n_domain, pos_domain, peri_len_x, previous_rank);
                        previous_rank = srank;
                        rank_send_of_ptcl[ip] = srank;
                    }
                    assert( determineWhetherParticleIsInDomain(pos_tmp, pos_domain[rank_send_of_ptcl[ip]]) );
                }
                for(S32 i=0; i<n_proc; i++){
                    n_send[i] = 0;
                    ptcl_send_buf_[i].clearSize();
                }
                for(S32 ip=0; ip<n_loc; ip++) {
                    const S32 rank_tmp = rank_send_of_ptcl[ip];
                    if(rank_tmp == my_rank){
                        ptcl_[iloc++] = ptcl_[ip];
                    }
                    else{
                        ptcl_send_buf_[rank_tmp].push_back(ptcl_[ip]);
                        n_send_tot++;
                        n_send[rank_tmp]++;
                    }
                }
                ptcl_.resizeNoInitialize(iloc);
            }
            time_profile_.exchange_particle__find_particle = GetWtime() - wtime_offset_inner;
            S32 dx_max = 0;
            S32 dy_max = 0;
            calcDxDyMax(dinfo, n_send, n_proc, dx_max, dy_max);
            assert(2*dx_max+1 < n_proc_x); // TO DO; have to consider this case.
            // ****************************************************
            // *** receive the number of exchange particles ********
            const S32 n_proc_comm_limit = (2*dx_max+1)*(2*dy_max+1);
            S32 * n_recv  = new S32[n_proc_comm_limit];
            S32 * n_disp_recv = new S32[n_proc_comm_limit+1];
            S32 * rank_comm = new S32[n_proc_comm_limit];
            MPI::Request * req_send = new MPI::Request[n_proc_comm_limit];
            MPI::Request * req_recv = new MPI::Request[n_proc_comm_limit];
            wtime_offset_inner = GetWtime();
            S32 n_proc_comm = 0;
            for(S32 ix=-dx_max; ix<=dx_max; ix++){
                S32 rank_x_target = my_rank_x+ix;
                if(rank_x_target < 0) rank_x_target = rank_x_target + n_proc_x;
                else if(rank_x_target >= n_proc_x) rank_x_target = rank_x_target - n_proc_x;
                assert(rank_x_target >= 0 && rank_x_target < n_proc_x);
                for(S32 iy=-dy_max; iy<=dy_max; iy++){
                    S32 rank_y_target = my_rank_y+iy;
                    if(rank_y_target < 0 || rank_y_target >= n_proc_y) continue;
                    S32 rank = (rank_x_target)*n_proc_y + (rank_y_target);
                    rank_comm[n_proc_comm] = rank;
                    req_send[n_proc_comm] = MPI::COMM_WORLD.Isend(n_send+rank,        1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_comm] = MPI::COMM_WORLD.Irecv(n_recv+n_proc_comm, 1, GetDataType<S32>(), rank, 0);
                    n_proc_comm++;
                }
            }
            assert(n_proc_comm <= n_proc_comm_limit);
            MPI::Request::Waitall(n_proc_comm, req_send);
            MPI::Request::Waitall(n_proc_comm, req_recv);
            
            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc_comm; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            ptcl_recv_.resizeNoInitialize( n_disp_recv[n_proc_comm] );

            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            for(S32 i=0; i<n_proc_comm; i++){
                S32 rank = rank_comm[i];
                if(n_send[rank] > 0) {
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_buf_[rank].getPointer(), n_send[rank], GetDataType<Tptcl>(), rank, 0);
                }
                if(n_recv[i] > 0) {
                    S32 adrrecv = n_disp_recv[i];
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), n_recv[i], GetDataType<Tptcl>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
            
            const S32 n_recv_tot = n_disp_recv[n_proc_comm];
            ptcl_.reserveEmptyAreaAtLeast( n_recv_tot );
            for(S32 ip=0; ip<n_recv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }

#if defined(DEBUG_PRINT_EX_PTCL_SW)
            Comm::barrier();
            if(my_rank == 0){
                std::cerr<<"ptcl_.size()= "<<ptcl_.size()<<std::endl;
                std::cerr<<"n_send_tot= "<<n_send_tot
                         <<" n_disp_recv[n_proc_comm]= "<<n_disp_recv[n_proc_comm]
                         <<" n_proc_comm= "<<n_proc_comm
                         <<std::endl;
            }
            for(S32 ip=0; ip<ptcl_.size(); ip++){
                F64 pos_x = ptcl_[ip].getPos().x;
                F64 pos_y = ptcl_[ip].getPos().y;
                F64 pos_z = ptcl_[ip].getPos().z;
                F64 pos_phi = atan2(pos_y, pos_x);
                if(pos_phi < 0.0) pos_phi += peri_len_x;
                F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                F64vec pos_tmp = F64vec(pos_phi, pos_r, pos_z);
                if(my_rank == 0){
                    if( ip % ((ptcl_.size()/100)+1) == 0){
                        std::cerr<<"pos_tmp= "<<pos_tmp
                                 <<" pos_my_domain)= "<<pos_my_domain
                                 <<std::endl;
                    }
                }
                assert(determineWhetherParticleIsInDomain(pos_tmp, pos_my_domain));
            }
            Comm::barrier();
#endif 
            // **************************************************** 
            n_ptcl_send_ += n_send_tot;
            n_ptcl_recv_ += n_disp_recv[n_proc_comm];
            time_profile_.exchange_particle__exchange_particle = GetWtime() - wtime_offset_inner;
            delete [] n_send;
            delete [] rank_send_of_ptcl;
            delete [] n_recv;
            delete [] n_disp_recv;
            delete [] rank_comm;
            delete [] req_send;
            delete [] req_recv;
            time_profile_.exchange_particle = GetWtime() - time_offset;
        }

        template<class Tdinfo>
        void exchangeParticleSw2(Tdinfo & dinfo,
                                 const F64 delta_ax){
            F64 time_offset = GetWtime();
            const F64 peri_len_x = 2.0 * 4.0 * atan(1.0);
            const S32 n_loc   = ptcl_.size();
            const S32 my_rank = MPI::COMM_WORLD.Get_rank();
            const S32 my_rank_x = dinfo.getRank1d(0);
            const S32 my_rank_y = dinfo.getRank1d(1);
            const S32 n_proc  = MPI::COMM_WORLD.Get_size();
            const S32 n_proc_x = dinfo.getNDomain(0);
            const S32 n_proc_y = dinfo.getNDomain(1);
            const S32 * n_domain = dinfo.getPointerOfNDomain();
            const F64ort * pos_domain = dinfo.getPointerOfPosDomain();
            const F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
            S32 * n_send = new S32[n_proc];
            S32 * rank_send_of_ptcl = new S32[n_loc];
            S32 n_send_tot = 0;
            S32 iloc = 0;
            F64 dx = peri_len_x / n_proc_x;
            F64 dy = delta_ax   / n_proc_y;
            F64 y_in  = 1.0 - (0.5*delta_ax);
            F64 y_out = 1.0 + (0.5*delta_ax);
            F64 wtime_offset_inner = GetWtime();
            ////////////////
            // initialize //
            for(S32 i=0; i<n_proc; i++){
                n_send[i] = 0;
                ptcl_send_buf_[i].clearSize();
            }            
            //////////////////////
            // search send rank //
            for(S32 ip=0; ip<n_loc; ip++) {
                F64 pos_x = ptcl_[ip].getPos().x;
                F64 pos_y = ptcl_[ip].getPos().y;
                F64 pos_z = ptcl_[ip].getPos().z;
                F64 pos_phi = atan2(pos_y, pos_x);
                if(pos_phi < 0.0) pos_phi += peri_len_x;
                F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                S32 rank_x_send = pos_phi / dx;
                S32 rank_y_send = (pos_r-y_in)   / dy;
                if(pos_r < y_in) rank_y_send = 0;
                else if(pos_r >= y_out) rank_y_send = n_proc_y-1;
                S32 rank_send = rank_x_send*n_proc_y + rank_y_send;
                rank_send_of_ptcl[ip] = rank_send;
                if(rank_send == my_rank){
                    ptcl_[iloc++] = ptcl_[ip];
                }
                else{
                    ptcl_send_buf_[rank_send].push_back(ptcl_[ip]);
                    n_send_tot++;
                    n_send[rank_send]++;
                }
                assert(rank_send < n_proc);
                if(rank_x_send != 0)          assert(pos_phi >= pos_domain[rank_send].low_.x);
                if(rank_x_send != n_proc_x-1) assert(pos_phi <  pos_domain[rank_send].high_.x);
                if(rank_y_send != 0)          assert(pos_r >= pos_domain[rank_send].low_.y);
                if(rank_y_send != n_proc_y-1) assert(pos_r <  pos_domain[rank_send].high_.y);
            }
            ptcl_.resizeNoInitialize(iloc);
            time_profile_.exchange_particle__find_particle = GetWtime() - wtime_offset_inner;
            S32 dx_max = 0;
            S32 dy_max = 0;
            calcDxDyMax(dinfo, n_send, n_proc, dx_max, dy_max);
            assert(2*dx_max+1 < n_proc_x); // TO DO; have to consider this case.
            // ****************************************************
            // *** receive the number of exchange particles ********
            const S32 n_proc_comm_limit = (2*dx_max+1)*(2*dy_max+1);
            S32 * n_recv  = new S32[n_proc_comm_limit];
            S32 * n_disp_recv = new S32[n_proc_comm_limit+1];
            S32 * rank_comm = new S32[n_proc_comm_limit];
            MPI::Request * req_send = new MPI::Request[n_proc_comm_limit];
            MPI::Request * req_recv = new MPI::Request[n_proc_comm_limit];
            wtime_offset_inner = GetWtime();
            S32 n_proc_comm = 0;
            for(S32 ix=-dx_max; ix<=dx_max; ix++){
                S32 rank_x_target = my_rank_x+ix;
                if(rank_x_target < 0) rank_x_target = rank_x_target + n_proc_x;
                else if(rank_x_target >= n_proc_x) rank_x_target = rank_x_target - n_proc_x;
                assert(rank_x_target >= 0 && rank_x_target < n_proc_x);
                for(S32 iy=-dy_max; iy<=dy_max; iy++){
                    S32 rank_y_target = my_rank_y+iy;
                    if(rank_y_target < 0 || rank_y_target >= n_proc_y) continue;
                    S32 rank = (rank_x_target)*n_proc_y + (rank_y_target);
                    rank_comm[n_proc_comm] = rank;
                    req_send[n_proc_comm] = MPI::COMM_WORLD.Isend(n_send+rank,        1, GetDataType<S32>(), rank, 0);
                    req_recv[n_proc_comm] = MPI::COMM_WORLD.Irecv(n_recv+n_proc_comm, 1, GetDataType<S32>(), rank, 0);
                    n_proc_comm++;
                }
            }
            assert(n_proc_comm <= n_proc_comm_limit);
            MPI::Request::Waitall(n_proc_comm, req_send);
            MPI::Request::Waitall(n_proc_comm, req_recv);
            n_disp_recv[0] = 0;
            for(S32 i=0; i<n_proc_comm; i++){
                n_disp_recv[i+1] = n_disp_recv[i] + n_recv[i];
            }
            ptcl_recv_.resizeNoInitialize( n_disp_recv[n_proc_comm] );
            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            for(S32 i=0; i<n_proc_comm; i++){
                S32 rank = rank_comm[i];
                if(n_send[rank] > 0) {
                    req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_buf_[rank].getPointer(), n_send[rank], GetDataType<Tptcl>(), rank, 0);
                }
                if(n_recv[i] > 0) {
                    S32 adrrecv = n_disp_recv[i];
                    req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), n_recv[i], GetDataType<Tptcl>(), rank, 0);
                }
            }
            MPI::Request::Waitall(n_proc_send, req_send);
            MPI::Request::Waitall(n_proc_recv, req_recv);
            
            const S32 n_recv_tot = n_disp_recv[n_proc_comm];
            ptcl_.reserveEmptyAreaAtLeast( n_recv_tot );
            for(S32 ip=0; ip<n_recv_tot; ip++) {
                ptcl_.pushBackNoCheck( ptcl_recv_[ip] );
            }
#if defined(DEBUG_PRINT_EX_PTCL_SW2)
            Comm::barrier();
            if(my_rank == 0){
                std::cerr<<"ptcl_.size()= "<<ptcl_.size()<<std::endl;
                std::cerr<<"n_send_tot= "<<n_send_tot
                         <<" n_disp_recv[n_proc_comm]= "<<n_disp_recv[n_proc_comm]
                         <<" n_proc_comm= "<<n_proc_comm
                         <<std::endl;
            }
            for(S32 ip=0; ip<ptcl_.size(); ip++){
                F64 pos_x = ptcl_[ip].getPos().x;
                F64 pos_y = ptcl_[ip].getPos().y;
                F64 pos_z = ptcl_[ip].getPos().z;
                F64 pos_phi = atan2(pos_y, pos_x);
                if(pos_phi < 0.0) pos_phi += peri_len_x;
                F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                F64vec pos_tmp = F64vec(pos_phi, pos_r, pos_z);
                if(my_rank == 0){
                    if( ip % ((ptcl_.size()/100)+1) == 0){
                        std::cerr<<"pos_tmp= "<<pos_tmp
                                 <<" pos_my_domain)= "<<pos_my_domain
                                 <<std::endl;
                    }
                }
                assert(determineWhetherParticleIsInDomain(pos_tmp, pos_my_domain));
            }
            Comm::barrier();
#endif 
            // **************************************************** 
            n_ptcl_send_ += n_send_tot;
            n_ptcl_recv_ += n_disp_recv[n_proc_comm];
            time_profile_.exchange_particle__exchange_particle = GetWtime() - wtime_offset_inner;
            delete [] n_send;
            delete [] rank_send_of_ptcl;
            delete [] n_recv;
            delete [] n_disp_recv;
            delete [] rank_comm;
            delete [] req_send;
            delete [] req_recv;
            time_profile_.exchange_particle = GetWtime() - time_offset;
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
            /*
            for(S32 i=0; i<Comm::getNumberOfProc(); i++){
                ptcl_send_buf_[i].freeMem();
            }
            */
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
