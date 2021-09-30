#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>

bool FLAG_PARA;

class FileHeader{
public:
    PS::F64 time;
    PS::S32 n_loc;
    PS::S32 n_glb;
    FileHeader(const PS::F64 t, const PS::S64 nl, const PS::S64 ng): time(t), n_loc(nl), n_glb(ng){}
    PS::S32 readAscii(FILE * fp){
	fscanf(fp, "%lf%d%d", &time, &n_loc, &n_glb);
	if(FLAG_PARA)
	    return n_loc;
	else
	    return n_glb;
    }
    void writeAscii(FILE* fp){
        fprintf(fp, "%lf   %d   %d\n", time, n_loc, n_glb);
    }
    PS::S32 readBinary(FILE * fp){
        fread(&time, sizeof(time), 1, fp);
        fread(&n_loc, sizeof(n_loc), 1, fp);
	fread(&n_glb, sizeof(n_glb), 1, fp);
	if(FLAG_PARA)
	    return n_loc;
	else
	    return n_glb;
    }
    void writeBinary(FILE * fp){
        fwrite(&time, sizeof(time), 1, fp);
        fwrite(&n_loc, sizeof(n_loc), 1, fp);
	fwrite(&n_glb, sizeof(n_glb), 1, fp);
    }
};

class FP{
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    void readAscii(FILE * fp){
	fscanf(fp, "%lld %lf %lf %lf %lf", &id, &mass, &pos.x, &pos.y, &pos.z);
    }
    void writeAscii(FILE * fp) const {
	fprintf(fp, "%lld  %lf  %lf  %lf  %lf\n", id, mass, pos.x, pos.y, pos.z);
    }
    void readBinary(FILE * fp){
	fread(&id,   sizeof(id), 1, fp);
	fread(&mass, sizeof(mass), 1, fp);
	fread(&pos.x, sizeof(pos.x), 1, fp);
	fread(&pos.y, sizeof(pos.y), 1, fp);
	fread(&pos.z, sizeof(pos.z), 1, fp);
    }    

    void writeBinary(FILE * fp){
        fwrite(&id, sizeof(id), 1, fp);
        fwrite(&mass, sizeof(mass), 1, fp);
        fwrite(&pos.x, sizeof(pos.x), 1, fp);
	fwrite(&pos.y, sizeof(pos.y), 1, fp);
	fwrite(&pos.z, sizeof(pos.z), 1, fp);
    }


    void readAscii2(FILE * fp){
	fscanf(fp, "%lld %lf %lf %lf %lf", &id, &mass, &pos.x, &pos.y, &pos.z);
    }
    void writeAscii2(FILE * fp) const {
	fprintf(fp, "%lld  %lf  %lf  %lf  %lf\n", id, mass, pos.x, pos.y, pos.z);
    }
    void readBinary2(FILE * fp){
	fread(&id,   sizeof(id), 1, fp);
	fread(&mass, sizeof(mass), 1, fp);
	fread(&pos.x, sizeof(pos.x), 1, fp);
	fread(&pos.y, sizeof(pos.y), 1, fp);
	fread(&pos.z, sizeof(pos.z), 1, fp);
    }    
    void writeBinary2(FILE * fp){
        fwrite(&id, sizeof(id), 1, fp);
        fwrite(&mass, sizeof(mass), 1, fp);
        fwrite(&pos.x, sizeof(pos.x), 1, fp);
	fwrite(&pos.y, sizeof(pos.y), 1, fp);
	fwrite(&pos.z, sizeof(pos.z), 1, fp);
    }
    
    PS::F64vec getPos() const {
	return pos;
    }
};

template<typename Tsys0, typename Tsys1>
void check(Tsys0 & sys0, Tsys1 & sys1){
    const auto n0 = sys0.getNumberOfParticleLocal();
    const auto n1 = sys1.getNumberOfParticleLocal();
    assert(n0 == n1);
    for(auto i=0; i<n0; i++){
	auto dr = sys1[i].pos - sys0[i].pos;
	assert(fabs(dr.x) < 1e-13);
	assert(fabs(dr.y) < 1e-13);
	assert(fabs(dr.z) < 1e-13);
    }
    for(auto i=0; i<n0; i++){
	sys1[i].id = -10;
	sys1[i].pos = -100;
    }
}


int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);
    PS::S32 my_rank = PS::Comm::getRank();

    PS::ParticleSystem<FP> system;
    PS::S32 n_loc = 5;
    system.setNumberOfParticleLocal(n_loc);
    PS::S32 n_glb = system.getNumberOfParticleGlobal();
    FileHeader header_out(0.0, n_loc, n_glb);
    FileHeader header_in(-1.0, -1, -10);
    for(auto i=0; i<n_loc; i++){
	system[i].id = my_rank*10 + i;
	system[i].mass = 1.0 / n_loc;
	system[i].pos = PS::F64vec((PS::F64)(i+my_rank), (PS::F64)(-2*i+my_rank), (PS::F64)(3*i+my_rank));
    }
    
    const char fmt[] = "%s_%05d_%05d.dat";

    const char base_ascii_0[] = "./result/base_ascii_0";
    FLAG_PARA = true;
    system.writeParticleAscii(base_ascii_0, fmt, header_out);
    PS::ParticleSystem<FP> system_in;
    system_in.readParticleAscii(base_ascii_0, fmt, header_in);
    check(system, system_in);
			
    system.writeParticleAscii(base_ascii_0, fmt, header_out, &FP::writeAscii2);
    system_in.readParticleAscii(base_ascii_0, fmt, header_in, &FP::readAscii2);
    check(system, system_in);			
       


    const char base_ascii_1[] = "./result/base_ascii_1";
    FLAG_PARA = true;
    system.writeParticleAscii(base_ascii_1, fmt);
    system_in.readParticleAscii(base_ascii_1, fmt);
    check(system, system_in);

    system.writeParticleAscii(base_ascii_1, fmt, &FP::writeAscii2);
    system_in.readParticleAscii(base_ascii_1, fmt, &FP::readAscii2);
    check(system, system_in);


    char siri_ascii_0[] = "./result/siri_ascii_0.dat";
    FLAG_PARA = false;
    system.writeParticleAscii(siri_ascii_0, header_out);
    system_in.readParticleAscii(siri_ascii_0, header_in);
    check(system, system_in);

    system.writeParticleAscii(siri_ascii_0, header_out, &FP::writeAscii2);
    system_in.readParticleAscii(siri_ascii_0, header_in, &FP::readAscii2);
    check(system, system_in);

    
    char siri_ascii_1[] = "./result/siri_ascii_1.dat";
    FLAG_PARA = false;
    system.writeParticleAscii(siri_ascii_1);
    system_in.readParticleAscii(siri_ascii_1);
    check(system, system_in);

    system.writeParticleAscii(siri_ascii_1, &FP::writeAscii2);
    system_in.readParticleAscii(siri_ascii_1, &FP::readAscii2);
    check(system, system_in);    


    const char base_binary_0[] = "./result/base_binary_0";
    FLAG_PARA = true;
    system.writeParticleBinary(base_binary_0, fmt, header_out);
    system_in.readParticleBinary(base_binary_0, fmt, header_in);
    check(system, system_in);

    system.writeParticleBinary(base_binary_0, fmt, header_out, &FP::writeBinary2);
    system_in.readParticleBinary(base_binary_0, fmt, header_in, &FP::readBinary2);
    check(system, system_in);
    

    
    const char base_binary_1[] = "./result/base_binary_1";
    FLAG_PARA = true;
    system.writeParticleBinary(base_binary_1, fmt);
    system_in.readParticleBinary(base_binary_1, fmt);
    check(system, system_in);

    system.writeParticleBinary(base_binary_1, fmt, &FP::writeBinary2);
    system_in.readParticleBinary(base_binary_1, fmt, &FP::readBinary2);
    check(system, system_in);
    
    
    char siri_binary_0[] = "./result/siri_binary_0.dat";
    FLAG_PARA = false;
    system.writeParticleBinary(siri_binary_0, header_out);
    system_in.readParticleBinary(siri_binary_0, header_in);
    check(system, system_in);

    system.writeParticleBinary(siri_binary_0, header_out, &FP::writeBinary2);
    system_in.readParticleBinary(siri_binary_0, header_in, &FP::readBinary2);
    check(system, system_in);
    

    
    char siri_binary_1[] = "./result/siri_binary_1.dat";
    FLAG_PARA = false;
    system.writeParticleBinary(siri_binary_1);
    system_in.readParticleBinary(siri_binary_1);
    check(system, system_in);

    system.writeParticleBinary(siri_binary_1, &FP::writeBinary2);
    system_in.readParticleBinary(siri_binary_1, &FP::readBinary2);
    check(system, system_in);


    PS::Finalize();
}
