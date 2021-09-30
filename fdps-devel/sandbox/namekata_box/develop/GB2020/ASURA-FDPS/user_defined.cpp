/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "mathematical_constants.h"
#include "debug_utilities.h"
#include "run_parameters.h"
#include "SPH_kernel.h"
#include "user_defined.h"
/* CELib header */
#include "CELib.h"

/* Implementations of user-defined classes */
//** File header class (used for IO)
PS::S32 FileHeader::readAscii(FILE* fp) {
    fscanf(fp, "%lf\n", &this->time);
    fscanf(fp, "%ld\n", &this->numPtcl);
    return this->numPtcl;
}
void FileHeader::writeAscii(FILE* fp) const {
    fprintf(fp, "%e\n",  this->time);
    fprintf(fp, "%ld\n", this->numPtcl);
}
PS::S32 FileHeader::readBinary(FILE* fp) {
    fread(&this->time, sizeof(this->time), 1, fp);
    fread(&this->numPtcl, sizeof(this->numPtcl), 1, fp);
    return this->numPtcl;
}
void FileHeader::writeBinary(FILE* fp) const {
    fwrite(&this->time, sizeof(this->time), 1, fp);
    fwrite(&this->numPtcl, sizeof(this->numPtcl), 1, fp);
}

PS::S32 TestHeader::readAscii(FILE* fp) {
    fscanf(fp, "%lf\n", &this->time);
    return -1;
}
void TestHeader::writeAscii(FILE* fp) const {
    fprintf(fp, "%e\n",  this->time);
}


//** Force class for gravity calculation
void Force_grav::clear() {
    this->acc = 0.0;
    this->pot = 0.0;
}

//** Force classes for SPH calculation
void Force_knl_sz::clear() {
    this->flag = 0;
    this->rad = 0.0;
}

void Force_dens::clear(){
    this->dens  = 0.0;
    this->pres  = 0.0;
    this->gradh = 0.0;
    this->divv  = 0.0;
    this->rotv  = 0.0;
}

void Force_hydro::clear(){
    this->acc = 0.0;
    this->eng_dot = 0.0;
}

//** Full Particle Classes
//**** FP_dm
PS::S64 FP_dm::getId() const {
    return this->id;
}
PS::F64vec FP_dm::getPos() const {
    return this->pos;
}
PS::F64 FP_dm::getCharge() const {
    return this->mass;
}
void FP_dm::setPos(const PS::F64vec& pos){
   this->pos = pos;
}
void FP_dm::copyFromForce(const Force_grav & f) {
    this->acc = f.acc;
    this->pot = f.pot;
}
void FP_dm::writeAscii(FILE* fp) const {
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    fprintf(fp,
            "%lld\t%e\t"
            "%e\t%e\t"
            "%e\t%e\n", 
            this->id, this->mass,
            this->pos.x, this->pos.y,
            this->vel.x, this->vel.y);
#else
    fprintf(fp,
            "%lld\t%e\t"
            "%e\t%e\t%e\t"
            "%e\t%e\t%e\n", 
            this->id, this->mass,
            this->pos.x, this->pos.y, this->pos.z,
            this->vel.x, this->vel.y, this->vel.z);
#endif
}
void FP_dm::readAscii(FILE* fp) {
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    fscanf(fp,
           "%ld\t%lf\t"
           "%lf\t%lf\t"
           "%lf\t%lf\n", 
           &this->id, &this->mass,
           &this->pos.x, &this->pos.y,
           &this->vel.x, &this->vel.y);
#else
    fscanf(fp,
           "%ld\t%lf\t"
           "%lf\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\n", 
           &this->id, &this->mass,
           &this->pos.x, &this->pos.y, &this->pos.z,
           &this->vel.x, &this->vel.y, &this->vel.z);
#endif
}
void FP_dm::writeRestartData(FILE* fp) const {
    fwrite(&this->id, sizeof(this->id), 1, fp);
    fwrite(&this->mass, sizeof(this->mass), 1, fp);
    fwrite(&this->pos, sizeof(this->pos), 1, fp);
    fwrite(&this->vel, sizeof(this->vel), 1, fp);
}
void FP_dm::readRestartData(FILE* fp) {
    fread(&this->id, sizeof(this->id), 1, fp);
    fread(&this->mass, sizeof(this->mass), 1, fp);
    fread(&this->pos, sizeof(this->pos), 1, fp);
    fread(&this->vel, sizeof(this->vel), 1, fp);
}
void FP_dm::clearGravitationalForce() {
    this->acc = 0.0;
    this->pot = 0.0;
}
void FP_dm::dump(const std::string msg_id,
                    const PS::S32 arr_idx,
                    const std::string caller_name,
                    std::ostream & fout) const {
    fout << msg_id
         << " (FP_dm)" 
         << " rank = " << PS::Comm::getRank()
         << " arr_idx = " << arr_idx
         << " id = " << this->id
         << " mass = " << this->mass
         << " pos = " << this->pos
         << " vel = " << this->vel
         << " acc = " << this->acc
         << " pot = " << this->pot
         << " @" << caller_name
         << std::endl;
}

//**** FP_star
FP_star::FP_star(): FBcnt(0) {} 
PS::S64 FP_star::getId() const {
    return this->id;
}
PS::F64vec FP_star::getPos() const {
    return this->pos;
}
PS::F64 FP_star::getCharge() const {
    return this->mass;
}
PS::F64 FP_star::getRSearch() const {
    return this->FBrad;
}
void FP_star::setPos(const PS::F64vec& pos){
   this->pos = pos;
}
void FP_star::copyFromForce(const Force_grav & f) {
    this->acc = f.acc;
    this->pot = f.pot;
}
void FP_star::copyFromForce(const Force_knl_sz & f) {
    this->flag  = f.flag;
    this->FBrad = f.rad;
}
void FP_star::copyFromForce(const Force_dens & f) {
    this->dens  = f.dens;
}
void FP_star::writeAscii(FILE* fp) const {
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    fprintf(fp,
            "%lld\t"
            "%e\t%e\t"
            "%e\t%e\t"
            "%e\t%e\t", 
            this->id,
            this->mass0, this->mass,
            this->pos.x, this->pos.y,
            this->vel.x, this->vel.y);
#else
    fprintf(fp,
            "%lld\t"
            "%e\t%e\t"
            "%e\t%e\t%e\t"
            "%e\t%e\t%e\t", 
            this->id,
            this->mass0, this->mass,
            this->pos.x, this->pos.y, this->pos.z,
            this->vel.x, this->vel.y, this->vel.z);
#endif
    for (PS::S32 k = 0; k < CELibYield_Number; k++)
        fprintf(fp, "%e\t", this->mabn[k]);
    fprintf(fp,
            "%e\t%llu\t"
            "%e\t%e\t%e\t%e\t"
            "%e\n",
            this->t_form, this->FBcnt,
            this->t_SNII, this->t_SNIa, this->t_AGB, this->t_NSM,
            this->FBrad);
}
void FP_star::readAscii(FILE* fp) {
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    fscanf(fp,
           "%ld\t"
           "%lf\t%lf\t"
           "%lf\t%lf\t"
           "%lf\t%lf\t", 
           &this->id,
           &this->mass0, &this->mass,
           &this->pos.x, &this->pos.y,
           &this->vel.x, &this->vel.y);
#else
    fscanf(fp,
           "%ld\t"
           "%lf\t%lf\t"
           "%lf\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\t", 
           &this->id,
           &this->mass0, &this->mass,
           &this->pos.x, &this->pos.y, &this->pos.z,
           &this->vel.x, &this->vel.y, &this->vel.z);
#endif
    for (PS::S32 k = 0; k < CELibYield_Number; k++)
        fscanf(fp, "%lf\t", &this->mabn[k]);
    fscanf(fp,
           "%lf\t%lu\t"
           "%lf\t%lf\t%lf\t%lf\t"
           "%lf\n",
           &this->t_form, &this->FBcnt,
           &this->t_SNII, &this->t_SNIa, &this->t_AGB, &this->t_NSM,
           &this->FBrad);
}
void FP_star::writeRestartData(FILE* fp) const {
    fwrite(&this->id, sizeof(this->id), 1, fp);
    fwrite(&this->mass0, sizeof(this->mass0), 1, fp);
    fwrite(&this->mass, sizeof(this->mass), 1, fp);
    fwrite(&this->pos, sizeof(this->pos), 1, fp);
    fwrite(&this->vel, sizeof(this->vel), 1, fp);
    fwrite(&this->mabn[0], sizeof(this->mabn[0]), CELibYield_Number, fp);
    fwrite(&this->t_form, sizeof(this->t_form), 1, fp);
    fwrite(&this->FBcnt, sizeof(this->FBcnt), 1, fp);
    fwrite(&this->t_SNII, sizeof(this->t_SNII), 1, fp);
    fwrite(&this->t_SNIa, sizeof(this->t_SNIa), 1, fp);
    fwrite(&this->t_AGB, sizeof(this->t_AGB), 1, fp);
    fwrite(&this->t_NSM, sizeof(this->t_NSM), 1, fp);
    fwrite(&this->FBrad, sizeof(this->FBrad), 1, fp);
}
void FP_star::readRestartData(FILE* fp) {
    fread(&this->id, sizeof(this->id), 1, fp);
    fread(&this->mass0, sizeof(this->mass0), 1, fp);
    fread(&this->mass, sizeof(this->mass), 1, fp);
    fread(&this->pos, sizeof(this->pos), 1, fp);
    fread(&this->vel, sizeof(this->vel), 1, fp);
    fread(&this->mabn[0], sizeof(this->mabn[0]), CELibYield_Number, fp);
    fread(&this->t_form, sizeof(this->t_form), 1, fp);
    fread(&this->FBcnt, sizeof(this->FBcnt), 1, fp);
    fread(&this->t_SNII, sizeof(this->t_SNII), 1, fp);
    fread(&this->t_SNIa, sizeof(this->t_SNIa), 1, fp);
    fread(&this->t_AGB, sizeof(this->t_AGB), 1, fp);
    fread(&this->t_NSM, sizeof(this->t_NSM), 1, fp);
    fread(&this->FBrad, sizeof(this->FBrad), 1, fp);
}
void FP_star::setId(const PS::S64 cid, const PS::S64 pid) {
    assert(0 <= cid && cid < (2<<nbit));
    assert(0 <= pid && pid < (0xffffffffffffffff >> (nbit+1)));
    this->id = pid | (cid << (64-(nbit+1)));
}
PS::S64 FP_star::getPid() const {
    PS::S64 id = this->id;
    return ((id << (nbit+1))>>(nbit+1));
}
PS::S64 FP_star::getCid() const {
    PS::S64 id = this->id;
    return static_cast<PS::S32>(id >> (64-(nbit+1)));
}
void FP_star::copyAbundanceFrom(const PS::F64 mabn[]) {
    for (PS::S32 i = 0; i < CELibYield_Number; i++) this->mabn[i] = mabn[i];
}
PS::F64 FP_star::getMetallicity() const {
    PS::F64 Z {0.0};
    for (PS::S32 i = 2; i < CELibYield_Number; i++) Z += this->mabn[i];
    return Z;
}
PS::U32 FP_star::getSNIICount() const {
    return (this->FBcnt >> 48);
}
void FP_star::setSNIICount(const PS::U64 cnt) {
    this->FBcnt &= 0x0000ffffffffffff; // clear
    this->FBcnt |= (cnt << 48);
}
PS::U32 FP_star::getSNIaCount() const {
    return ((this->FBcnt >> 32) & 0xffff);
}
void FP_star::setSNIaCount(const PS::U64 cnt) {
    this->FBcnt &= 0xffff0000ffffffff; // clear
    this->FBcnt |= (cnt << 32);
}
PS::U32 FP_star::getAGBCount() const {
    return ((this->FBcnt >> 16) & 0xffff);
}
void FP_star::setAGBCount(const PS::U64 cnt) {
    this->FBcnt &= 0xffffffff0000ffff; // clear
    this->FBcnt |= (cnt << 16);
}
PS::U32 FP_star::getNSMCount() const {
    return (this->FBcnt & 0xffff);
}
void FP_star::setNSMCount(const PS::U64 cnt) {
    this->FBcnt &= 0xffffffffffff0000; // clear
    this->FBcnt |= cnt;
}
bool FP_star::feedbackAsSNII(const PS::F64 t, const PS::F64 dt) const {
    return (t <= this->t_SNII && this->t_SNII < t+dt);
}
bool FP_star::feedbackAsSNIa(const PS::F64 t, const PS::F64 dt) const {
    return (t <= this->t_SNIa && this->t_SNIa < t+dt);
}
bool FP_star::feedbackAsAGB(const PS::F64 t, const PS::F64 dt) const {
    return (t <= this->t_AGB  && this->t_AGB  < t+dt);
}
bool FP_star::feedbackAsNSM(const PS::F64 t, const PS::F64 dt) const {
    return (t <= this->t_NSM  && this->t_NSM  < t+dt);
}
void FP_star::clearGravitationalForce() {
    this->acc = 0.0;
    this->pot = 0.0;
}
void FP_star::dump(const std::string msg_id,
                   const PS::S32 arr_idx,
                   const std::string caller_name,
                   std::ostream & fout) const {
    fout << msg_id
         << " (FP_star)" 
         << " rank = " << PS::Comm::getRank()
         << " arr_idx = " << arr_idx
         << " id = " << this->id
         << " mass0 = " << this->mass0
         << " mass = " << this->mass
         << " pos = " << this->pos
         << " vel = " << this->vel
         << " acc = " << this->acc
         << " pot = " << this->pot
         << " @" << caller_name
         << std::endl;
}
//*** FP_gas
FP_gas::FP_gas() {
    this->h_dot_prev = 0.0;
    this->n_stars = 0;
}
PS::S64 FP_gas::getId() const {
    return this->id;
}
PS::F64 FP_gas::getCharge() const {
    return this->mass;
}
PS::F64vec FP_gas::getPos() const {
    return this->pos;
}
PS::F64 FP_gas::getRSearch() const {
    return this->h;
}
void FP_gas::setPos(const PS::F64vec& pos){
   this->pos = pos;
}
void FP_gas::copyFromForce(const Force_grav & f) {
    this->acc_grav = f.acc;
    this->pot_grav = f.pot;
}
void FP_gas::copyFromForce(const Force_knl_sz & f) {
    this->flag = f.flag;
    this->h = f.rad;
}
void FP_gas::copyFromForce(const Force_dens & f){
    this->dens  = f.dens;
    this->pres  = f.pres;
    this->gradh = f.gradh;
    this->divv  = f.divv;
    this->rotv  = f.rotv;
}
void FP_gas::copyFromForce(const Force_hydro& f){
    this->acc_hydro = f.acc;
    this->eng_dot   = f.eng_dot;
    this->dt        = f.dt;
}
void FP_gas::writeAscii(FILE* fp) const {
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    fprintf(fp,
            "%lld\t%e\t%e\t"
            "%e\t%e\t"
            "%e\t%e\t"
            "%e\t%e\t%e\t",
            this->id, this->mass0, this->mass,
            this->pos.x, this->pos.y,
            this->vel.x, this->vel.y,
            this->eng, this->dens, this->h);
#else
    fprintf(fp,
            "%lld\t%e\t%e\t"
            "%e\t%e\t%e\t"
            "%e\t%e\t%e\t"
            "%e\t%e\t%e\t",
            this->id, this->mass0, this->mass,
            this->pos.x, this->pos.y, this->pos.z,
            this->vel.x, this->vel.y, this->vel.z,
            this->eng, this->dens, this->h);
#endif
    for (PS::S32 k = 0; k < CELibYield_Number; k++)
        fprintf(fp, "%e\t", this->mabn[k]);
    fprintf(fp, "%d\n", this->n_stars);
}
void FP_gas::readAscii(FILE* fp){
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    fscanf(fp,
           "%ld\t%lf\t%lf\t"
           "%lf\t%lf\t"
           "%lf\t%lf\t"
           "%lf\t%lf\t%lf\t",
           &this->id, &this->mass0, &this->mass,
           &this->pos.x, &this->pos.y,
           &this->vel.x, &this->vel.y,
           &this->eng, &this->dens, &this->h);
#else
    fscanf(fp,
           "%ld\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\t",
           &this->id, &this->mass0, &this->mass,
           &this->pos.x, &this->pos.y, &this->pos.z,
           &this->vel.x, &this->vel.y, &this->vel.z,
           &this->eng, &this->dens, &this->h);
#endif
    for (PS::S32 k = 0; k < CELibYield_Number; k++)
        fscanf(fp, "%lf\t", &this->mabn[k]);
    fscanf(fp, "%d\n", &this->n_stars);
}
void FP_gas::writeRestartData(FILE* fp) const {
    fwrite(&this->id, sizeof(this->id), 1, fp);
    fwrite(&this->mass0, sizeof(this->mass0), 1, fp);
    fwrite(&this->mass, sizeof(this->mass), 1, fp);
    fwrite(&this->pos, sizeof(this->pos), 1, fp);
    fwrite(&this->vel, sizeof(this->vel), 1, fp);
    fwrite(&this->eng, sizeof(this->eng), 1, fp);
    fwrite(&this->h, sizeof(this->h), 1, fp);
    fwrite(&this->alpha, sizeof(this->alpha), 1, fp);
    fwrite(&this->mabn[0], sizeof(this->mabn[0]), CELibYield_Number, fp);
    fwrite(&this->n_stars, sizeof(this->n_stars), 1, fp);
}
void FP_gas::readRestartData(FILE* fp) {
    fread(&this->id, sizeof(this->id), 1, fp);
    fread(&this->mass0, sizeof(this->mass0), 1, fp);
    fread(&this->mass, sizeof(this->mass), 1, fp);
    fread(&this->pos, sizeof(this->pos), 1, fp);
    fread(&this->vel, sizeof(this->vel), 1, fp);
    fread(&this->eng, sizeof(this->eng), 1, fp);
    fread(&this->h, sizeof(this->h), 1, fp);
    fread(&this->alpha, sizeof(this->alpha), 1, fp);
    fread(&this->mabn[0], sizeof(this->mabn[0]), CELibYield_Number, fp);
    fread(&this->n_stars, sizeof(this->n_stars), 1, fp);
}
void FP_gas::writeGlassData(FILE* fp) const {
    fwrite(&this->pos, sizeof(this->pos), 1, fp);
}
void FP_gas::readGlassData(FILE* fp) {
    fread(&this->pos, sizeof(this->pos), 1, fp);
}
PS::F64 FP_gas::getMass() const {
    return this->mass;
}
PS::F64 FP_gas::getKernelSupportRadius() const {
    return this->h;
}
PS::F64 FP_gas::getMassDensity() const {
    return this->dens;
}
PS::F64 FP_gas::getTemperature() const {
    return run_param::ism::mu * this->pres / this->dens;
}
PS::F64 FP_gas::getInternalEnergyDensity() const {
    return this->pres / (run_param::ism::gamma - 1.0);
}
PS::F64 FP_gas::getMetallicity() const {
    PS::F64 Z {0.0};
    for (PS::S32 i = 2; i < CELibYield_Number; i++) Z += this->mabn[i];
    return Z;
}
void FP_gas::calcSoundSpeed(){
    this->snds = std::sqrt(run_param::ism::gamma * this->pres / this->dens);
    this->BalSW = std::fabs(this->divv)
                / ( std::fabs(this->divv)
                  + std::sqrt(this->rotv * this->rotv)
                  + 1.0e-4 * this->snds / this->h); 
}
void FP_gas::applyTemperatureLimits() {
#ifdef ASURA_FDPS_USE_LOWER_TEMPERATURE_LIMIT
    {
        const PS::F64 T = (run_param::sph::T_lower_limit/run_param::unit::temp);
        const PS::F64 eng = T/(run_param::ism::mu*(run_param::ism::gamma-1.0));
        this->eng = std::max(this->eng, eng);
        this->eng_half = std::max(this->eng_half, eng);
    }
#endif
#ifdef ASURA_FDPS_USE_UPPER_TEMPERATURE_LIMIT
    {
        const PS::F64 T = (run_param::sph::T_upper_limit/run_param::unit::temp);
        const PS::F64 eng = T/(run_param::ism::mu*(run_param::ism::gamma-1.0));
        this->eng = std::min(this->eng, eng);
        this->eng_half = std::min(this->eng_half, eng);
    }
#endif
}
void FP_gas::clearGravitationalForce() {
    this->acc_grav = 0.0;
    this->pot_grav = 0.0;
}
void FP_gas::dump(const std::string msg_id,
                  const PS::S32 arr_idx,
                  const std::string caller_name,
                  std::ostream & fout) const {
    fout << msg_id
         << " (FP_gas)" 
         << " rank = " << PS::Comm::getRank()
         << " arr_idx = " << arr_idx
         << " id = " << this->id
         << " mass0 = " << this->mass0
         << " mass = " << this->mass
         << " pos = " << this->pos
         << " vel = " << this->vel
         << " vel_half = " << this->vel_half
         << " acc_grav = " << this->acc_grav
         << " acc_hydro = " << this->acc_hydro
         << " eng = " << this->eng
         << " eng_half = " << this->eng_half
         << " eng_dot = " << this->eng_dot
         << " dens = " << this->dens
         << " pres = " << this->pres
         << " h = " << this->h
         << " h_prev = " << this->h_prev
         << " h_dot_prev = " << this->h_dot_prev
         << " alpha = " << this->alpha
         << " @" << caller_name
         << std::endl;
}

//** Essential Particle Class
//*** EP_grav
PS::S64 EP_grav::getId() const {
    return this->id;
}
PS::F64 EP_grav::getCharge() const {
    return this->mass;
}
PS::F64 EP_grav::getEps2() const {
    return SQ(this->eps);
}
PS::F64vec EP_grav::getPos() const {
    return this->pos;
}
void EP_grav::copyFromFP(const FP_dm& fp) {
    this->id   = fp.id;
    this->mass = fp.mass;
    this->pos  = fp.pos;
    this->eps  = run_param::grav::soft::eps_dm;
}
void EP_grav::copyFromFP(const FP_star& fp) {
    this->id   = fp.id;
    this->mass = fp.mass;
    this->pos  = fp.pos;
    this->eps  = run_param::grav::soft::eps_star;
}
void EP_grav::copyFromFP(const FP_gas& fp) {
    this->id   = fp.id;
    this->mass = fp.mass;
    this->pos  = fp.pos;
    this->eps  = run_param::grav::soft::eps_gas;
}
void EP_grav::writeAscii(std::ostream & fout) const { /* do nothing */ }
void EP_grav::dump(const std::string msg_id,
                   const PS::S32 arr_idx,
                   const std::string caller_name,
                   std::ostream & fout) const {
    fout << msg_id
         << " (EP_grav)" 
         << " rank = " << PS::Comm::getRank()
         << " arr_idx = " << arr_idx
         << " id = " << this->id
         << " mass = " << this->mass
         << " pos = " << this->pos
         << " eps = " << this->eps
         << " @" << caller_name
         << std::endl;
}

//*** EPI_knl_sz
PS::S64 EPI_knl_sz::getId() const {
    return this->id;
}
PS::F64vec EPI_knl_sz::getPos() const {
    return this->pos;
}
PS::F64 EPI_knl_sz::getRSearch() const {
    return this->rad;
}
PS::F64 EPI_knl_sz::getRPhysical() const {
    return this->rad;
}
void EPI_knl_sz::setPos(const PS::F64vec& pos) {
    this->pos = pos;
}
void EPI_knl_sz::copyFromFP(const FP_gas& fp) {
    this->type = ParticleType::Gas;
    this->flag = fp.flag;
    this->id   = fp.id;
    this->pos  = fp.pos;
    this->rad  = fp.h;
}
void EPI_knl_sz::copyFromFP(const FP_star& fp) {
    this->type = ParticleType::Star;
    this->flag = fp.flag;
    this->id   = fp.id;
    this->pos  = fp.pos;
    this->rad  = fp.FBrad;
}
void EPI_knl_sz::writeAscii(std::ostream & fout) const {/* do nothing */}
void EPI_knl_sz::dump(const std::string msg_id,
                      const PS::S32 arr_idx,
                      const std::string caller_name,
                      std::ostream & fout) const {
    fout << msg_id
         << " (EPI_knl_sz)" 
         << " rank = " << PS::Comm::getRank()
         << " arr_idx = " << arr_idx
         << " id = " << this->id
         << " pos = " << this->pos
         << " rad = " << this->rad
         << " @" << caller_name
         << std::endl;
}

//*** EPJ_knl_sz
PS::F64vec EPJ_knl_sz::getPos() const {
    return this->pos;
}
void EPJ_knl_sz::setPos(const PS::F64vec& pos) {
    this->pos = pos;
}
void EPJ_knl_sz::copyFromFP(const FP_gas& fp) {
    this->type = ParticleType::Gas;
    this->pos  = fp.pos;
}
void EPJ_knl_sz::copyFromFP(const FP_star& fp) {
    this->type = ParticleType::Star;
    this->pos  = fp.pos;
}
void EPJ_knl_sz::writeAscii(std::ostream & fout) const { /* do nothing */ }
void EPJ_knl_sz::dump(const std::string msg_id,
                     const PS::S32 arr_idx,
                     const std::string caller_name,
                     std::ostream & fout) const {
    fout << msg_id
         << " (EPJ_knl_sz)" 
         << " rank = " << PS::Comm::getRank()
         << " arr_idx = " << arr_idx
         << " pos = " << this->pos
         << " @" << caller_name
         << std::endl;
}

//*** EP_hydro
PS::S64 EP_hydro::getId() const {
    return this->id;
}
PS::F64vec EP_hydro::getPos() const {
    return this->pos;
}
PS::F64 EP_hydro::getRSearch() const {
    return this->h;
}
PS::F64 EP_hydro::getRPhysical() const {
    return this->h;
}
void EP_hydro::setPos(const PS::F64vec& pos) {
    this->pos = pos;
}
void EP_hydro::copyFromFP(const FP_gas& fp){
    this->type  = ParticleType::Gas; 
    this->rank  = PS::Comm::getRank();
    this->idx   = fp.idx;
    this->id    = fp.id;
    this->pos   = fp.pos;
    this->vel   = fp.vel;
    this->mass  = fp.mass;
    this->eng   = fp.eng;
    this->h     = fp.h;
    this->dens  = fp.dens;
    this->pres  = fp.pres;
    this->gradh = fp.gradh;
    this->snds  = fp.snds;
    this->BalSW = fp.BalSW;
    this->alpha = fp.alpha;
}
void EP_hydro::copyFromFP(const FP_star& fp) {
    this->type  = ParticleType::Star; 
    this->rank  = -1; // not used
    this->idx   = -1; // not used
    this->id    = fp.id;
    this->pos   = fp.pos;
    this->vel   = 0.0; // not used
    this->mass  = 0.0; // treated as mass-less particle
    this->eng   = 0.0; // treated as energy-less particle
    this->h     = fp.FBrad;
    this->dens  = 0.0; // not used
    this->pres  = 0.0; // not used
    this->gradh = 0.0; // not used
    this->snds  = 0.0; // not used
    this->BalSW = 0.0; // not used
    this->alpha = 0.0; // not used
}
void EP_hydro::writeAscii(std::ostream & fout) const {
    const std::string delim = "   ";
    fout << this->id << delim
         << this->pos.x << delim
         << this->pos.y << delim
         << this->pos.z << delim
         << this->vel.x << delim
         << this->vel.y << delim
         << this->vel.z << delim
         << this->mass << delim
         << this->eng << delim
         << this->h << delim
         << this->dens << delim
         << this->pres << delim
         << this->gradh << delim
         << this->snds << delim
         << this->BalSW << delim
         << this->alpha << std::endl;
}
void EP_hydro::dump(const std::string msg_id,
                    const PS::S32 arr_idx,
                    const std::string caller_name,
                    std::ostream & fout) const {
    fout << msg_id
         << " (EP_hydro)" 
         << " rank = " << PS::Comm::getRank()
         << " arr_idx = " << arr_idx
         << " idx = " << this->idx
         << " id = " << this->id
         << " pos = " << this->pos
         << " vel = " << this->vel
         << " mass = " << this->mass
         << " eng = " << this->eng
         << " h = " << this->h
         << " dens = " << this->dens
         << " pres = " << this->pres
         << " gradh = " << this->gradh
         << " snds = " << this->snds
         << " BalSW = " << this->BalSW
         << " alpha = " << this->alpha
         << " @" << caller_name
         << std::endl;
}

bool output_flag {false};
PS::F64 et_rij_calc {0.0};
PS::F64 et_bisec {0.0};

PS::F64 getPredictedKernelSize(const PS::S32 N_ngb_cur,
                               const PS::F64 h_cur) {
    // Using the current number of neighbors and the current
    // kernel size, this function predicts a kernel size 
    // in which the specified number of neighbors is contained.
    // See Eq.(10),(11) in Thacker et al.(2000)[MNRAS,319,619].
    // Some modifications are made by Daisuke Namekata.
    constexpr PS::F64 limiter = 10.0;
    const PS::F64 N_ngb = run_param::sph::N_ngb; // casted to F64
#if ASURA_FDPS_KERNEL_SIZE_PRED_METHOD == ASURA_FDPS_THACKER_2000
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    const PS::F64 s = std::sqrt(N_ngb/std::max(N_ngb_cur,1));
#else
    const PS::F64 s = std::cbrt(N_ngb/std::max(N_ngb_cur,1));
#endif
    PS::F64 a;
    if (s < 1.0) {
        a = 0.2 * (1.0 + SQ(s));
    } else {
        const PS::F64 sinv = 1.0/s;
        a = 0.2 * (1.0 + CUBE(sinv));
    }
    return h_cur * std::min(limiter, (1.0 - a + a*s));
#elif ASURA_FDPS_KERNEL_SIZE_PRED_METHOD == ASURA_FDPS_THACKER_2000_MODIFIED
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    assert(N_ngb > 4.0);
    const PS::F64 sm = std::sqrt(N_ngb/(N_ngb + 4.0));
    const PS::F64 sp = std::sqrt(N_ngb/(N_ngb - 4.0));
    const PS::F64  s = std::sqrt(N_ngb/std::max(N_ngb_cur, 1));
#else
    assert(N_ngb > 6.0);
    const PS::F64 sm = std::cbrt(N_ngb/(N_ngb + 6.0));
    const PS::F64 sp = std::cbrt(N_ngb/(N_ngb - 6.0));
    const PS::F64  s = std::cbrt(N_ngb/std::max(N_ngb_cur, 1));
#endif
    PS::F64 a,s0;
    if (s < 1.0) {
        a = 0.2 * (1.0 + SQ(s));
        s0 = (1.0 - sm)/std::sqrt(std::log(1.5));
    } else {
        const PS::F64 sinv = 1.0/s;
        a = 0.2 * (1.0 + CUBE(sinv));
        s0 = (sp - 1.0)/std::sqrt(std::log(1.5));
    }
    const PS::F64 p = 0.5 * std::exp(-SQ((s-1.0)/s0));
    return h_cur * std::min(limiter, (1.0 - p) * (1.0 - a + a*s) + p);
#endif // ASURA_FDPS_KERNEL_SIZE_PRED_METHOD
}


/* Interaction functions */
void CalcKernelSize::operator () (const EPI_knl_sz * ep_i,
                                  const PS::S32 n_ip,
                                  const EPJ_knl_sz * ep_j,
                                  const PS::S32 n_jp,
                                  Force_knl_sz * force) {
    // for debug
    //if (output_flag) dbg_utils::fout << n_jp << std::endl;
    // Allocate caches
    std::vector<PS::F64> rij;
    rij.resize(n_jp);
    // Density calculation
    for (PS::S32 i = 0; i < n_ip; i++) {
#if defined(ASURA_FDPS_USE_FIXED_KERNEL_SIZE)
        force[i].flag = 1;
        force[i].rad = ep_i[i].rad;
        continue;
#endif
        // Skip this particle if its radius has been already determined.
        if (ep_i[i].flag == 1) {
            force[i].flag = 1;
            force[i].rad = ep_i[i].rad;
            continue;
        }
        // Cache rij[]
        PS::F64 h = ep_i[i].rad; // maximum allowance because of rsrch=h
        assert(h > 0.0);
        PS::F64 start_time = PS::GetWtime();
        PS::S32 N_ngb = 0;
        for (PS::S32 j = 0; j < n_jp; j++) {
            const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
            if (ep_j[j].type == ParticleType::Gas) rij[j] = std::sqrt(dr * dr); // gas particle
            else rij[j] = 2.0 * h; // star particle
            if (rij[j] < h) N_ngb++;
        }
        et_rij_calc += PS::GetWtime() - start_time;
        // Set the allowed range of the number of neighbor particles
        PS::S32 N_ngb_tgt, N_ngb_min, N_ngb_max;
        if (ep_i[i].type == ParticleType::Gas) {
            N_ngb_tgt = run_param::sph::N_ngb;
            N_ngb_min = run_param::sph::N_ngb - run_param::sph::N_ngb_mgn;
            N_ngb_max = run_param::sph::N_ngb + run_param::sph::N_ngb_mgn;
        } else if (ep_i[i].type == ParticleType::Star) {
            N_ngb_tgt = run_param::fb::N_ngb;
            N_ngb_min = run_param::fb::N_ngb - run_param::fb::N_ngb_mgn;
            N_ngb_max = run_param::fb::N_ngb + run_param::fb::N_ngb_mgn;
        }
        // Adjust h so that the number of neighbors is in the specified range
        start_time = PS::GetWtime();
#if defined(ASURA_FDPS_USE_NTH_ELEMENT_TO_FIND_KERNEL_SIZE)
        if (N_ngb < N_ngb_tgt) {
            // In this case, we need to increase the search radius and
            // perform calcForceAllAndWriteBack again. At this time,
            // we skip this particle forcibly.
            const PS::F64 h_next = run_param::sph::h2h_next * h;
            force[i].flag = 0;
            force[i].rad = std::max(h_next, getPredictedKernelSize(N_ngb, h)); 
        } else {
            // Find an appropriate h by using std::nth_element
            constexpr PS::F64 eps = 1.0e-6;
            std::nth_element(rij.begin(), rij.begin() + N_ngb_tgt - 1, rij.end());
            force[i].flag = 1;
            force[i].rad = (1.0 + eps) * rij[N_ngb_tgt - 1];
            // eps is introduced to include the (N_ngb_tgt-1)th particle as
            // a neighbor.
        }
#else
        if (N_ngb < N_ngb_min) {
            // In this case, we need to increase the search radius and
            // perform calcForceAllAndWriteBack again. At this time,
            // we skip this particle forcibly.
            const PS::F64 h_next = run_param::sph::h2h_next * h;
            force[i].flag = 0;
            force[i].rad = std::max(h_next, getPredictedKernelSize(N_ngb, h)); 
        } else {
            // Find an appropriate h by the bisection method.
            PS::F64 h_L = 0.0;
            PS::F64 h_U = h;
            PS::F64 h_bisec = 0.75 * h;
            PS::S32 iter = 1;
            for (;;) {
                // Calculate number of neighbor particles
                N_ngb = 0;
                for (PS::S32 j = 0; j < n_jp; j++) if (rij[j] < h_bisec) N_ngb++;
                // Check if the termination condition is satisfied or not
                if ((N_ngb_min <= N_ngb) && (N_ngb <= N_ngb_max)) {
                    h = h_bisec;
                    break;
                }
                // Calculate the next h_bisec and update h_L & h_U
                if ((N_ngb < N_ngb_min) && (h_L < h_bisec)) h_L = h_bisec;
                if ((N_ngb > N_ngb_max) && (h_bisec < h_U)) h_U = h_bisec;
                h_bisec = 0.5*(h_L + h_U);
                iter++;
            }
            force[i].flag = 1;
            force[i].rad = h;
        }
#endif // ASURA_FDPS_USE_NTH_ELEMENT_TO_FIND_KERNEL_SIZE
        et_bisec += PS::GetWtime() - start_time;
    }
}

void CalcDensityAndPressure::operator () (const EP_hydro * ep_i,
                               const PS::S32 n_ip,
                               const EP_hydro * ep_j,
                               const PS::S32 n_jp,
                               Force_dens * force){
    // for debug
    if (output_flag) dbg_utils::fout << "n_jp = " << n_jp << std::endl;
    // Local constants
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    constexpr PS::F64 DIM = 2.0;
#else
    constexpr PS::F64 DIM = 3.0;
#endif
    constexpr PS::F64 eps = 1.0e-6; // introduced to prevent the grad-h term from diverging
    // Allocate caches
    PS::F64 * rij = (PS::F64 *)malloc(sizeof(PS::F64) * n_jp);
    // Density calculation
    for (PS::S32 i = 0; i < n_ip; i++) {
        //==================
        //   Gas particle
        //==================
        if (ep_i[i].type == ParticleType::Gas) {
            // Cache rij[]
            const PS::F64 h = ep_i[i].h;
            assert(h > 0.0);
            for (PS::S32 j = 0; j < n_jp; j++) {
                const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
                if (ep_j[j].type == ParticleType::Gas) rij[j] = std::sqrt(dr * dr); // gas particle
                else rij[j] = 2.0 * h; // To exclude star particles from the summations below
            }
            // Calculate the smoothed density and pressure (U)
            PS::F64 dens = 0.0;
            PS::F64 pres = 0.0;
            for (PS::S32 j = 0; j < n_jp; j++) {
                const PS::F64 mj = ep_j[j].mass;
                const PS::F64 uj = ep_j[j].eng;
                const PS::F64 Wij = W(rij[j], h);
                dens += mj * Wij;
                pres += (run_param::ism::gamma -1.0) * mj * uj * Wij;
            }
            force[i].dens = dens;
            force[i].pres = pres;
#if ASURA_FDPS_GRAD_H_TERM == ASURA_FDPS_EQ15_IN_HOPKINS_2013
            // Calculate the grad-h term in Eq.(59) in Saitoh & Makino (2013)
            // or Eq.(15) in Hopkins (2013).
            PS::F64 dPdh = 0.0;
            for (PS::S32 j = 0; j < n_jp; j++) {
                const PS::F64 mj = ep_j[j].mass;
                const PS::F64 uj = ep_j[j].eng;
                const PS::F64 dWdh_ij = dWdh(rij[j], h);
                dPdh += (run_param::ism::gamma - 1.0) * mj * uj* dWdh_ij;
            }
            force[i].gradh = 1.0 / (eps + 1.0 + (h * dPdh)/(DIM * pres));
#elif ASURA_FDPS_GRAD_H_TERM == ASURA_FDPS_EQ18_IN_HOPKINS_2013
            // Calculate the grad-h term in Eq.(18) in Hopkins (2013)
            // Note that here we only calculate A PART OF the grad-h term,
            // which depends only on i-particles in Eq.(18).
            PS::F64 n = 0.0;
            PS::F64 dndh = 0.0;
            PS::F64 dPdh = 0.0;
            for (PS::S32 j = 0; j < n_jp; j++) {
                const PS::F64 mj = ep_j[j].mass;
                const PS::F64 uj = ep_j[j].eng;
                const PS::F64 W_cstr_ij = W_cstr(rij[j], h);
                const PS::F64 dWdh_cstr_ij = dWdh_cstr(rij[j], h);
                const PS::F64 dWdh_ij = dWdh(rij[j], h);
                n += W_cstr_ij;
                dndh += dWdh_cstr_ij;
                dPdh += (run_param::ism::gamma - 1.0) * mj * uj* dWdh_ij;
            }
            force[i].gradh = ((h * dPdh)/(DIM * n)) / (eps + 1.0 + (h * dndh)/(DIM * n));
#endif
            // Check the grad-h term
#if defined(ENABLE_NAN_CHECK)
            if (std::isnan(force[i].gradh) ||
                std::isinf(force[i].gradh)) {
                ep_i[i].dump("[nan(i)]", i, __func__, dbg_utils::fout);
                assert(false);
            }
#endif
            // Calculate \div v & \rot v for Balsara switch
            for (PS::S32 j = 0; j < n_jp; j++) {
                const PS::F64 mj = ep_j[j].mass;
                const PS::F64 uj = ep_j[j].eng;
                const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
                const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
                force[i].divv -= (run_param::ism::gamma - 1.0) * mj * uj
                                * dv * gradW(dr, h);
                force[i].rotv -= (run_param::ism::gamma - 1.0) * mj * uj
                                * dv ^ gradW(dr, h);
            }
            force[i].divv /= force[i].pres;
            force[i].rotv /= force[i].pres;
        //==================
        //   FB particle
        //==================
        } else if (ep_i[i].type == ParticleType::Star) {
            const PS::F64vec pos_i = ep_i[i].pos;
            const PS::F64 h_i = ep_i[i].h;
            PS::F64 dens = 0.0;
            for (PS::S32 j = 0; j < n_jp; j++) {
                const PS::F64 mj = ep_j[j].mass;
                const PS::F64vec dr = ep_j[j].pos - pos_i;
                const PS::F64 rij = std::sqrt(dr * dr);
                dens += mj * W(rij, h_i);
            }
            force[i].dens = dens;
        }
    }
    // Release memory
    free(rij);
}

void CalcHydroForce::operator () (const EP_hydro * ep_i,
                                  const PS::S32 n_ip,
                                  const EP_hydro * ep_j,
                                  const PS::S32 n_jp,
                                  Force_hydro * force) {
    const PS::F64 gamma = run_param::ism::gamma;
    for (PS::S32 i = 0; i < n_ip; i++){
        const PS::F64vec pos_i = ep_i[i].pos;
        const PS::F64vec vel_i = ep_i[i].vel;
        const PS::F64 mi = ep_i[i].mass;
        const PS::F64 hi = ep_i[i].h;
        const PS::F64 rhoi = ep_i[i].dens;
        const PS::F64 ui = ep_i[i].eng;
        const PS::F64 Pi = ep_i[i].pres;
        const PS::F64 fi = ep_i[i].gradh;
        const PS::F64 ci = ep_i[i].snds;
        PS::F64 v_sig_max = 0.0;
        for (PS::S32 j = 0; j < n_jp; j++){
            const PS::F64vec dr = pos_i - ep_j[j].pos;
            const PS::F64vec dv = vel_i - ep_j[j].vel;
            const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / std::sqrt(dr * dr) : 0;
            const PS::F64 v_sig = ci + ep_j[j].snds - 3.0 * w_ij;
            v_sig_max = std::max(v_sig_max, v_sig);
            const PS::F64 alpha_ij = 0.5 * (ep_i[i].alpha + ep_j[j].alpha);
            const PS::F64 AV = - 0.5 * alpha_ij * v_sig * w_ij
                               / (0.5 * (rhoi + ep_j[j].dens))
                               * 0.5 * (ep_i[i].BalSW + ep_j[j].BalSW);
            const PS::F64vec gradW_i  = gradW(dr, hi);
            const PS::F64vec gradW_j  = gradW(dr, ep_j[j].h);
            const PS::F64vec gradW_ij = 0.5 * (gradW_i + gradW_j);
            const PS::F64 mj = ep_j[j].mass;
            const PS::F64 uj = ep_j[j].eng;
            const PS::F64 Pj = ep_j[j].pres;
            const PS::F64 fj = ep_j[j].gradh;
#if ASURA_FDPS_GRAD_H_TERM == ASURA_FDPS_EQ15_IN_HOPKINS_2013
            force[i].acc -= SQ(gamma - 1.0) * mj * ui * uj
                            * ( fi * gradW_i / Pi
                              + fj * gradW_j / Pj)
                          + mj * AV * gradW_ij;
            force[i].eng_dot += SQ(gamma - 1.0) * mj * ui * uj * fi * gradW_i * dv / Pi
                              + mj * 0.5 * AV * gradW_ij * dv;
#elif ASURA_FDPS_GRAD_H_TERM == ASURA_FDPS_EQ18_IN_HOPKINS_2013
            const PS::F64 fij = 1.0 - fi/((gamma - 1.0) * mj * uj);
            const PS::F64 fji = 1.0 - fj/((gamma - 1.0) * mi * ui);
            force[i].acc -= SQ(gamma - 1.0) * mj * ui * uj
                            * ( fij * gradW_i / Pi
                              + fji * gradW_j / Pj)
                          + mj * AV * gradW_ij;
            force[i].eng_dot += SQ(gamma - 1.0) * mj * ui * uj * fij * gradW_i * dv / Pi
                              + mj * 0.5 * AV * gradW_ij * dv;
#endif

#if defined(ENABLE_NAN_CHECK)
            if (force[i].acc.isnan() ||
                force[i].acc.isinf() ||
                std::isnan(force[i].eng_dot) ||
                std::isinf(force[i].eng_dot)) {
                ep_i[i].dump("[nan(i)]", i, __func__, dbg_utils::fout);
                ep_j[j].dump("[nan(j)]", j, __func__, dbg_utils::fout);
                assert(false);
            }
#endif
        }
        force[i].dt = run_param::sph::CFL * 2.0 * ep_i[i].h / v_sig_max;
    }
}
