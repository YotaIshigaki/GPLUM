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
void Force_dens::clear(){
    this->flag  = 0;
    this->dens  = 0.0;
    this->gradh = 0.0;
    this->divv  = 0.0;
    this->rotv  = 0.0;
}
void Force_hydro::clear(){
    this->acc = 0.0;
    this->eng_dot = 0.0;
    this->ent_dot = 0.0;
}

//** Full Particle Classes
//**** FP_nbody
PS::S64 FP_nbody::getId() const {
    return this->id;
}
PS::F64vec FP_nbody::getPos() const {
    return this->pos;
}
PS::F64 FP_nbody::getCharge() const {
    return this->mass;
}
void FP_nbody::setPos(const PS::F64vec& pos){
   this->pos = pos;
}
void FP_nbody::copyFromForce(const Force_grav & f) {
    this->acc = f.acc;
    this->pot = f.pot;
}
void FP_nbody::writeAscii(FILE* fp) const {
    fprintf(fp,
            "%lld\t%e\t"
            "%e\t%e\t%e\t"
            "%e\t%e\t%e\n", 
            this->id, this->mass,
            this->pos.x, this->pos.y, this->pos.z,
            this->vel.x, this->vel.y, this->vel.z);
}
void FP_nbody::readAscii(FILE* fp) {
    fscanf(fp,
           "%ld\t%lf\t"
           "%lf\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\n", 
           &this->id, &this->mass,
           &this->pos.x, &this->pos.y, &this->pos.z,
           &this->vel.x, &this->vel.y, &this->vel.z);
}
void FP_nbody::writeRestartData(FILE* fp) const {
    fwrite(&this->id, sizeof(this->id), 1, fp);
    fwrite(&this->mass, sizeof(this->mass), 1, fp);
    fwrite(&this->pos, sizeof(this->pos), 1, fp);
    fwrite(&this->vel, sizeof(this->vel), 1, fp);
}
void FP_nbody::readRestartData(FILE* fp) {
    fread(&this->id, sizeof(this->id), 1, fp);
    fread(&this->mass, sizeof(this->mass), 1, fp);
    fread(&this->pos, sizeof(this->pos), 1, fp);
    fread(&this->vel, sizeof(this->vel), 1, fp);
}
void FP_nbody::dump(const std::string msg_id,
                    const PS::S32 arr_idx,
                    const std::string caller_name,
                    std::ostream & fout) const {
    fout << msg_id
         << " (FP_nbody)" 
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
void FP_star::copyFromForce(const Force_dens & f) {
    this->flag  = f.flag;
    this->dens  = f.dens;
    this->FBrad = f.smth;
}
void FP_star::writeAscii(FILE* fp) const {
    fprintf(fp,
            "%lld\t%lld\t"
            "%e\t%e\t"
            "%e\t%e\t%e\t"
            "%e\t%e\t%e\t", 
            this->pid, this->id,
            this->mass0, this->mass,
            this->pos.x, this->pos.y, this->pos.z,
            this->vel.x, this->vel.y, this->vel.z);
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
    fscanf(fp,
           "%ld\t%ld\t"
           "%lf\t%lf\t"
           "%lf\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\t", 
           &this->pid, &this->id,
           &this->mass0, &this->mass,
           &this->pos.x, &this->pos.y, &this->pos.z,
           &this->vel.x, &this->vel.y, &this->vel.z);
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
    fwrite(&this->pid, sizeof(this->pid), 1, fp);
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
    fread(&this->pid, sizeof(this->pid), 1, fp);
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
void FP_star::dump(const std::string msg_id,
                   const PS::S32 arr_idx,
                   const std::string caller_name,
                   std::ostream & fout) const {
    fout << msg_id
         << " (FP_star)" 
         << " rank = " << PS::Comm::getRank()
         << " arr_idx = " << arr_idx
         << " pid = " << this->pid
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
FP_gas::FP_gas(): n_stars(0) {}
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
    return this->smth;
}
void FP_gas::setPos(const PS::F64vec& pos){
   this->pos = pos;
}
void FP_gas::copyFromForce(const Force_grav& f) {
    this->acc_grav = f.acc;
    this->pot_grav = f.pot;
}
void FP_gas::copyFromForce(const Force_dens& f){
    this->flag  = f.flag;
    this->dens  = f.dens;
    this->smth  = f.smth;
    this->gradh = f.gradh;
    this->divv  = f.divv;
    this->rotv  = f.rotv;
    
}
void FP_gas::copyFromForce(const Force_hydro& f){
    this->acc_hydro = f.acc;
    this->eng_dot   = f.eng_dot;
    this->ent_dot   = f.ent_dot;
    this->dt        = f.dt;
}
void FP_gas::writeAscii(FILE* fp) const {
    fprintf(fp,
            "%lld\t%e\t%e\t"
            "%e\t%e\t%e\t"
            "%e\t%e\t%e\t"
            "%e\t%e\t%e\t",
            this->id, this->mass0, this->mass,
            this->pos.x, this->pos.y, this->pos.z,
            this->vel.x, this->vel.y, this->vel.z,
            this->eng, this->dens, this->smth);
    for (PS::S32 k = 0; k < CELibYield_Number; k++)
        fprintf(fp, "%e\t", this->mabn[k]);
    fprintf(fp, "%d\n", this->n_stars);
}
void FP_gas::readAscii(FILE* fp){
    fscanf(fp,
           "%ld\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\t"
           "%lf\t%lf\t%lf\t",
           &this->id, &this->mass0, &this->mass,
           &this->pos.x, &this->pos.y, &this->pos.z,
           &this->vel.x, &this->vel.y, &this->vel.z,
           &this->eng, &this->dens, &this->smth);
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
    fwrite(&this->dens, sizeof(this->dens), 1, fp);
    fwrite(&this->smth, sizeof(this->smth), 1, fp);
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
    fread(&this->dens, sizeof(this->dens), 1, fp);
    fread(&this->smth, sizeof(this->smth), 1, fp);
    fread(&this->mabn[0], sizeof(this->mabn[0]), CELibYield_Number, fp);
    fread(&this->n_stars, sizeof(this->n_stars), 1, fp);
}
PS::F64 FP_gas::getMass() const {
    return this->mass;
}
PS::F64 FP_gas::getKernelSupportRadius() const {
    return this->smth;
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
void FP_gas::setEntropy(){
    this->ent = (run_param::ism::gamma - 1.0) * this->eng 
              / std::pow(this->dens, run_param::ism::gamma - 1.0);
}
void FP_gas::setPressure(){
#if defined(ISOTHERMAL_EOS)
    // In this case, eng = const.
    this->pres = (run_param::ism::gamma - 1.0) * this->dens * this->eng;
    this->ent  = this->pres / std::pow(this->dens, run_param::ism::gamma);
#else
#if defined(USE_ENTROPY)
    this->pres = this->ent * std::pow(this->dens, run_param::ism::gamma);
    this->eng  = this->pres / ((run_param::ism::gamma - 1.0) * this->dens);
#else
    this->pres = (run_param::ism::gamma - 1.0) * this->dens * this->eng;
    this->ent  = this->pres / std::pow(this->dens, run_param::ism::gamma);
#endif
#endif
    this->snds = std::sqrt(run_param::ism::gamma * this->pres / this->dens);
#if defined(USE_BALSARA_SWITCH)
    this->BalSW = std::fabs(this->divv)
                / ( std::fabs(this->divv)
                  + std::sqrt(this->rotv * this->rotv)
                  + 1.0e-4 * this->snds / this->smth); 
#else
    this->BalSW = 1.0;
#endif
}
void FP_gas::setPressureFromSpecificInternalEnergy(const PS::F64 eng_new) {
    this->eng  = eng_new;
    this->pres = (run_param::ism::gamma - 1.0) * this->dens * this->eng;
    this->ent  = this->pres / std::pow(this->dens, run_param::ism::gamma);
    this->snds = std::sqrt(run_param::ism::gamma * this->pres / this->dens);
#if defined(USE_BALSARA_SWITCH)
    this->BalSW = std::fabs(this->divv)
                / ( std::fabs(this->divv) 
                  + std::sqrt(this->rotv * this->rotv)
                  + 1.0e-4 * this->snds / this->smth); 
#else
    this->BalSW = 1.0;
#endif
}
void FP_gas::setPressureFromInternalEnergyDensity(const PS::F64 U) {
    this->pres = (run_param::ism::gamma - 1.0) * U;
    this->eng  = this->pres / ((run_param::ism::gamma - 1.0) * this->dens);
    this->ent  = this->pres / std::pow(this->dens, run_param::ism::gamma);
    this->snds = std::sqrt(run_param::ism::gamma * this->pres / this->dens);
#if defined(USE_BALSARA_SWITCH)
    this->BalSW = std::fabs(this->divv) 
                / ( std::fabs(this->divv)
                  + std::sqrt(this->rotv * this->rotv)
                  + 1.0e-4 * this->snds / this->smth); 
#else
    this->BalSW = 1.0;
#endif
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
         << " ent = " << this->ent
         << " ent_half = " << this->ent_half
         << " ent_dot = " << this->ent_dot
         << " dens = " << this->dens
         << " pres = " << this->pres
         << " smth = " << this->smth
         << " @" << caller_name
         << std::endl;
}

//** Essential Particle Class
//*** EP_grav
PS::F64 EP_grav::eps;
PS::S64 EP_grav::getId() const {
    return this->id;
}
PS::F64 EP_grav::getCharge() const {
    return this->mass;
}
PS::F64vec EP_grav::getPos() const {
    return this->pos;
}
void EP_grav::copyFromFP(const FP_nbody& fp) {
    this->id   = fp.id;
    this->mass = fp.mass;
    this->pos  = fp.pos;
}
void EP_grav::copyFromFP(const FP_star& fp) {
    this->id   = fp.id;
    this->mass = fp.mass;
    this->pos  = fp.pos;
}
void EP_grav::copyFromFP(const FP_gas& fp) {
    this->id   = fp.id;
    this->mass = fp.mass;
    this->pos  = fp.pos;
}
//*** EP_hydro
PS::S64 EP_hydro::getId() const {
    return this->id;
}
PS::F64vec EP_hydro::getPos() const {
    return this->pos;
}
PS::F64 EP_hydro::getRSearch() const {
    return run_param::sim::SCF_smth * this->smth;
}
void EP_hydro::setPos(const PS::F64vec& pos) {
    this->pos = pos;
}
void EP_hydro::copyFromFP(const FP_gas& fp){
    this->type  = 0; // to indicate `gas` particle
    this->rank  = PS::Comm::getRank();
    this->idx   = fp.idx;
    this->id    = fp.id;
    this->pos   = fp.pos;
    this->vel   = fp.vel;
    this->mass  = fp.mass;
    this->smth  = fp.smth;
    this->dens  = fp.dens;
    this->pres  = fp.pres;
    this->gradh = fp.gradh;
    this->snds  = fp.snds;
    this->BalSW = fp.BalSW;
}
void EP_hydro::copyFromFP(const FP_star& fp) {
    this->type  = 1; // to indicate `star` particle
    this->rank  = -1; // not used
    this->idx   = -1; // not used
    this->id    = fp.id;
    this->pos   = fp.pos;
    this->vel   = 0.0; // not used
    this->mass  = 0.0; // treated as mass-less particle
    this->smth  = fp.FBrad;
    this->dens  = 0.0; // not used
    this->pres  = 0.0; // not used
    this->gradh = 0.0; // not used
    this->snds  = 0.0; // not used
    this->BalSW = 0.0; // not used
}

/* Interaction functions */
void CalcDensity::operator () (const EP_hydro * ep_i,
                               const PS::S32 n_ip,
                               const EP_hydro * ep_j,
                               const PS::S32 n_jp,
                               Force_dens * force){
    const PS::F64 eps = 1.0e-6;
    const PS::F64 M_trgt = run_param::sim::mass_avg
                         * run_param::sim::N_neighbor;
    const PS::S32 mgn = 1;
    const PS::F64 coeff = 4.0 * math_const::pi / 3.0;
    for (PS::S32 i = 0; i < n_ip; i++) {
        //==================
        //   Gas particle
        //==================
        if (ep_i[i].type == 0) {
            bool iter_flag {true};
            PS::F64 h = ep_i[i].smth;
            const PS::F64 h_max_alw = run_param::sim::SCF_smth * h; // maximum allowance
            PS::F64 h_L = 0.0;
            PS::F64 h_U = h_max_alw;
            // Software caches
            PS::F64 * mj  = (PS::F64 *)malloc(sizeof(PS::F64) * n_jp);
            PS::F64 * rij = (PS::F64 *)malloc(sizeof(PS::F64) * n_jp);
            for (PS::S32 j = 0; j < n_jp; j++) {
                mj[j] = ep_j[j].mass;
                const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
                rij[j] = std::sqrt(dr * dr);
            }
            // Check if the range [h_L, h_U] brackets the root or not
            PS::F64 dens_U = 0.0;
            for (PS::S32 j = 0; j < n_jp; j++) {
               dens_U += mj[j] * W(rij[j], h_U);
            }
            const PS::F64 M_U = coeff * CUBE(h_U) * dens_U;
            const PS::F64 f_U = M_U/M_trgt - 1.0;
            if (f_U <= 0.0) {
                // In this case, [h_L, h_U] does not brackets the root.
                // Hence, we skip this particle forcibly.
                // In order to determine consistently the density
                // and the smoothing length for this particle,
                // we must re-perform calcForceAllAndWriteBack().
                force[i].flag = 0;
                force[i].dens = dens_U;
                force[i].smth = h_max_alw;
                iter_flag = false;
            }
            // Find the root using both the bisection search and the Newton method
            if (iter_flag) {
                const PS::S32 iter_max=100;
                const PS::S32 iter_crit=80;
                const PS::S32 iter_switch=10;
                PS::S32 iter = 0;

                PS::F64 h_bisec = h;
                PS::F64 h_newton = h;
                for (;;) {
#if !defined(USE_NEWTON_METHOD_FROM_THE_BEGINNING)
                    if (iter < iter_switch) {
                        //------------------------------------------------------------------
                        // In this case, we use the bisection search only to find the root.
                        //------------------------------------------------------------------
                        // Calculate density
                        PS::F64 dens_bisec {0.0};
                        for (PS::S32 j = 0; j < n_jp; j++) {
                            dens_bisec  += mj[j] * W(rij[j], h_bisec);
                        }
                        // Check if the current value of the smoohting length satisfies 
                        // Eq.(5) in Springel (2005).
                        const PS::F64 M_bisec = coeff * CUBE(h_bisec) * dens_bisec;
                        const PS::F64 f_bisec = M_bisec/M_trgt - 1.0;
                        if (std::abs(f_bisec) < eps) {
                            // In this case, Eq.(5) holds within a specified accuracy.
                            force[i].flag = 1;
                            force[i].dens = dens_bisec;
                            force[i].smth = h_bisec;
                            break;
                        }
                        // Calculate the next h_bisec and update h_L & h_U
                        if ((f_bisec < 0.0) && (h_L < h_bisec)) h_L = h_bisec;
                        if ((f_bisec > 0.0) && (h_bisec < h_U)) h_U = h_bisec;
                        h_bisec = 0.5*(h_L + h_U);
                    } else {
                        //------------------------------------------------------------------
                        // In this case, we use the bisection search and the Newton method
                        // to find the root.
                        //------------------------------------------------------------------
                        if (iter == iter_switch) h_newton = h_bisec; // initialize h_newton
#endif // !USE_NEWTON_METHOD_FROM_THE_BEGINNING
                        // Calculate density
                        PS::F64 dens_bisec {0.0};
                        PS::F64 dens_newton {0.0};
                        PS::F64 d1dens_newton {0.0};
                        for (PS::S32 j = 0; j < n_jp; j++) {
                            dens_bisec  += mj[j] * W(rij[j], h_bisec);
                            dens_newton += mj[j] * W(rij[j], h_newton);
                            d1dens_newton += mj[j] * dWdh(rij[j], h_newton);
                        }
                        // Check if the current value of the smoohting length satisfies 
                        // Eq.(5) in Springel (2005).
                        const PS::F64 M_bisec  = coeff * CUBE(h_bisec) * dens_bisec;
                        const PS::F64 f_bisec  = M_bisec/M_trgt - 1.0;
                        const PS::F64 M_newton = coeff * CUBE(h_newton) * dens_newton;
                        const PS::F64 f_newton = M_newton/M_trgt - 1.0;
                        if ((std::abs(f_bisec) < eps) || (std::abs(f_newton) < eps)) {
                            // In this case, Eq.(5) holds within a specified accuracy.
                            if (std::abs(f_newton) < std::abs(f_bisec)) {
                                force[i].flag = 1;
                                force[i].dens = dens_newton;
                                force[i].smth = h_newton;
                            } else {
                                force[i].flag = 1;
                                force[i].dens = dens_bisec;
                                force[i].smth = h_bisec;
                            }
                            break;
                        }
                        // Calculate the next h_bisec and update h_L & h_U
                        if ((f_bisec < 0.0) && (h_L < h_bisec)) h_L = h_bisec;
                        if ((f_bisec > 0.0) && (h_bisec < h_U)) h_U = h_bisec;
                        h_bisec = 0.5*(h_L + h_U);
                        // Calculate the next h_newton
                        const PS::F64 h_newton_prev = h_newton;
                        const PS::F64 df_newton = coeff * SQ(h_newton) * (3.0 * dens_newton + h_newton * d1dens_newton);
                        if (df_newton != 0.0) {
                            h_newton -= f_newton/df_newton;
                            if (h_newton <= h_L || h_U <= h_newton) h_newton = h_bisec;
                            // In this case, h_newton is too far from the true solution,
                            // and it seems that a first-order Talyor series approximation
                            // used in Newton's method is not good. Hence, we set h_bisec
                            // to h_newton as we expect h_bisec is a better guess.
                        } else {
                            // In this case, this particle interacts itself only.
                            // Because W(0,hi) = hi * dWdh(0,hi) holds, df_newton = 0.
                            // We do not use Newton's method in this case.
                            h_newton = h_bisec;
                        }
#if 0
                        // for debugging [start]
                        if (h_bisec <= 0.0 || h_newton <= 0.0) {
                            if (PS::Comm::getRank() == 31 && PS::Comm::getThreadNum() == 0) {
                            //if (PS::Comm::getRank() == 31) {
                            dbg_utils::fout << "[h <= 0] " << std::endl;
                            dbg_utils::fout << " iter = " << iter << std::endl;
                            dbg_utils::fout << " rank = " << PS::Comm::getRank() << std::endl;
                            dbg_utils::fout << " th = " << PS::Comm::getThreadNum() << std::endl;
                            dbg_utils::fout << " id = " << ep_i[i].id << std::endl;
                            dbg_utils::fout << " dens_bisec = " << dens_bisec << std::endl;
                            dbg_utils::fout << " dens_newton = " << dens_newton << std::endl;
                            dbg_utils::fout << " d1dens_newton = " << d1dens_newton << std::endl;
                            dbg_utils::fout << " h_bisec = " << h_bisec << std::endl;
                            dbg_utils::fout << " f_bisec = " << f_bisec << std::endl;
                            dbg_utils::fout << " h_newton = " << h_newton << std::endl;
                            dbg_utils::fout << " h_newton_prev = " << h_newton_prev << std::endl;
                            dbg_utils::fout << " f_newton = " << f_newton << std::endl;
                            dbg_utils::fout << " df_newton = " << df_newton << std::endl;
                            dbg_utils::fout << " f/df = " << f_newton/df_newton << std::endl;
                            dbg_utils::fout << " h_L = " << h_L << std::endl;
                            dbg_utils::fout << " h_U = " << h_U << std::endl;
                            assert(false);
                            }
                        }
                        // for debugging [end]
#endif

#if defined(SET_LIMIT_TO_ITERATION_IN_CALC_DENSITY)
                        // Check
                        if (iter == iter_max) {
                            std::cout << "Too many iteration in " << __func__ << std::endl;
                            assert(false);
                        }
                        if (iter >= iter_crit) {
                            dbg_utils::fout << "iter = " << iter 
                                            << " rank = " << PS::Comm::getRank()
                                            << " id = " << ep_i[i].id 
                                            << " dens_bisec = " << dens_bisec
                                            << " dens_newton = " << dens_newton
                                            << " h_bisec = " << h_bisec
                                            << " h_newton = " << h_newton
                                            << " h_L = " << h_L 
                                            << " h_U = " << h_U 
                                            << " dh = " << (h_U - h_L)
                                            << " f_bisec = " << f_bisec
                                            << " f_newton = " << f_newton
                                            << std::endl;
                        }
#endif // SET_LIMIT_TO_ITERATION_IN_CALC_DENSITY
#if !defined(USE_NEWTON_METHOD_FROM_THE_BEGINNING)
                    }
#endif // !USE_NEWTON_METHOD_FROM_THE_BEGINNING

                    // Update iter
                    iter++;

                }
            }
            // Calculate grad-h term
            if (force[i].flag == 1) {
                const PS::F64 h = force[i].smth;
                const PS::F64 dens = force[i].dens;
                PS::F64 drho_dh = 0.0;
                for (PS::S32 j = 0; j < n_jp; j++) {
                   drho_dh += mj[j] * dWdh(rij[j], h);
                }
                force[i].gradh = 1.0 / (1.0 + (h * drho_dh) / (3.0 * dens));
            } 
            else {
                force[i].gradh = 1.0; // dummy value
            }
#if defined(USE_BALSARA_SWITCH)
            // Compute \div v & \rot v for Balsara switch
            for (PS::S32 j = 0; j < n_jp; j++) {
               const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
               const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;
               force[i].divv -= mj[j] * dv * gradW(dr, force[i].smth);
               force[i].rotv -= mj[j] * dv ^ gradW(dr, force[i].smth);
            }
            force[i].divv /= force[i].dens;
            force[i].rotv /= force[i].dens;
#endif
            // Release memory
            free(mj);
            free(rij);
        //==================
        //   FB particle
        //==================
        } else if (ep_i[i].type == 1) {
            bool iter_flag {true};
            PS::F64 h = ep_i[i].smth;
            const PS::F64 h_max_alw = run_param::sim::SCF_smth * h; // maximum allowance
            PS::F64 h_L = 0.0;
            PS::F64 h_U = h_max_alw;
            // Software caches
            PS::F64 * mj  = (PS::F64 *)malloc(sizeof(PS::F64) * n_jp);
            PS::F64 * rij = (PS::F64 *)malloc(sizeof(PS::F64) * n_jp);
            for (PS::S32 j = 0; j < n_jp; j++) {
                mj[j] = ep_j[j].mass;
                const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
                if (ep_j[j].type == 0) rij[j] = std::sqrt(dr * dr); // gas particle
                else rij[j] = 2.0 * h_max_alw; // star particle
            }
            // Check if the range [h2_L, h2_U] brackets the root or not
            PS::S32 num_U {0};
            for (PS::S32 j = 0; j < n_jp; j++) if (rij[j] < h_U) num_U++;
            if (num_U < run_param::sim::N_neighbor_FB) {
                // In this case, [h2_L, h2_U] does not brackets the root.
                // Hence, we skip this particle forcibly.
                // In order to determine the feedback radius for this particle,
                // we must re-perform calcForceAllAndWriteBack().
                force[i].flag = 0;
                force[i].dens = 1.0; // dummy value
                force[i].smth = h_max_alw;
                iter_flag = false;
            }
            // Find the root using both the bisection search and the Newton method
            if (iter_flag) {
                PS::F64 h_bisec = h;
                for (;;) {
                    //------------------------------------------------------------------
                    // In this case, we use the bisection search only to find the root.
                    //------------------------------------------------------------------
                    // Calculate number of particles
                    PS::F64 num {0};
                    for (PS::S32 j = 0; j < n_jp; j++) if (rij[j] < h_bisec) num++;
                    //std::cout << "num = " << num << std::endl;
                    // Check if the termination condition is satisfied or not
                    if ((run_param::sim::N_neighbor_FB - mgn <= num) &&
                        (num <= run_param::sim::N_neighbor_FB + mgn)) {
                        force[i].flag = 1;
                        force[i].dens = 0.0;
                        for (PS::S32 j = 0; j < n_jp; j++) {
                            force[i].dens += mj[j] * W(rij[j], h_bisec);
                        }
                        force[i].smth = h_bisec;
                        break;
                    }
                    // Calculate the next h_bisec and update h_L & h_U
                    if ((num < run_param::sim::N_neighbor_FB - mgn) && (h_L < h_bisec)) h_L = h_bisec;
                    if ((num > run_param::sim::N_neighbor_FB + mgn) && (h_bisec < h_U)) h_U = h_bisec;
                    h_bisec = 0.5*(h_L + h_U);

                }
            }
            // Release memory
            free(mj);
            free(rij);
        }
    }
}

void CalcHydroForce::operator () (const EP_hydro * ep_i,
                                  const PS::S32 n_ip,
                                  const EP_hydro * ep_j,
                                  const PS::S32 n_jp,
                                  Force_hydro * force) {
    for (PS::S32 i = 0; i < n_ip; i++){
        const PS::F64vec pos_i = ep_i[i].pos;
        const PS::F64vec vel_i = ep_i[i].vel;
        const PS::F64 smth_i   = ep_i[i].smth;
        const PS::F64 dens_i   = ep_i[i].dens;
        const PS::F64 pres_i   = ep_i[i].pres;
        const PS::F64 f_i      = ep_i[i].gradh;
        const PS::F64 snds_i   = ep_i[i].snds;
        const PS::F64 povrho2_i = pres_i / (dens_i * dens_i);
        PS::F64 v_sig_max = 0.0;
        for (PS::S32 j = 0; j < n_jp; j++){
            const PS::F64vec dr = pos_i - ep_j[j].pos;
            const PS::F64vec dv = vel_i - ep_j[j].vel;
            const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / std::sqrt(dr * dr) : 0;
            const PS::F64 v_sig = snds_i + ep_j[j].snds - 3.0 * w_ij;
            v_sig_max = std::max(v_sig_max, v_sig);
            const PS::F64 AV = - 0.5 * run_param::sim::alpha_AV * v_sig * w_ij / (0.5 * (dens_i + ep_j[j].dens))
                               * 0.5 * (ep_i[i].BalSW + ep_j[j].BalSW);
            const PS::F64vec gradW_i  = gradW(dr, smth_i);
            const PS::F64vec gradW_j  = gradW(dr, ep_j[j].smth);
            const PS::F64vec gradW_ij = 0.5 * (gradW_i + gradW_j);
            const PS::F64 povrho2_j = ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens);
            const PS::F64 f_j = ep_j[j].gradh;
            force[i].acc     -= ep_j[j].mass * (f_i * povrho2_i * gradW_i
                                               +f_j * povrho2_j * gradW_j
                                               +AV * gradW_ij);
            force[i].eng_dot += ep_j[j].mass * (f_i * povrho2_i * gradW_i
                                               +0.5 * AV * gradW_ij) * dv;
            force[i].ent_dot += 0.5 * ep_j[j].mass * AV * gradW_ij * dv;
            // for debugging [start]
            if (force[i].acc.isnan() ||
                force[i].acc.isinf() ||
                std::isnan(force[i].eng_dot) ||
                std::isinf(force[i].eng_dot) ||
                std::isnan(force[i].ent_dot) ||
                std::isinf(force[i].ent_dot)) {
                dbg_utils::fout << "[hyd] rank = " << PS::Comm::getRank()
                                << " id = " << ep_i[i].id 
                                << " pos_i = " << pos_i 
                                << " vel_i = " << vel_i
                                << " smth_i = " << smth_i
                                << " dens_i = " << dens_i
                                << " pres_i = " << pres_i
                                << " f_i    = " << f_i
                                << " snds_i = " << snds_i
                                << " BalSW_i = " << ep_i[i].BalSW
                                << std::endl;
                dbg_utils::fout << "pos_j = " << ep_j[j].pos
                                << " vel_j = " << ep_j[j].vel
                                << " smth_j = " << ep_j[j].smth
                                << " dens_j = " << ep_j[j].dens
                                << " pres_j = " << ep_j[j].pres
                                << " f_j = " << f_j
                                << " snds_j = " << ep_j[j].snds
                                << " BalSW_j = " << ep_j[j].BalSW
                                << std::endl;
                dbg_utils::fout << "w_ij  = " << w_ij
                                << " v_sig = " << v_sig
                                << " AV = " << AV
                                << " gradW_i = " << gradW_i
                                << " gradW_j = " << gradW_j
                                << std::endl;
                assert(false);
            }
            // for debugging [end]
        }
        const PS::F64 p = run_param::ism::gamma - 1.0;
        force[i].ent_dot *= p/std::pow(dens_i, p);
        force[i].dt = run_param::sim::CFL_hydro * 2.0 * ep_i[i].smth / v_sig_max;
    }
}
