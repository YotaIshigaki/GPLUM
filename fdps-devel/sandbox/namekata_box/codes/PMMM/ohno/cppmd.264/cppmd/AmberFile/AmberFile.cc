#include "AmberFile.h"

#include <cstdlib>

#include "AmberPrmtop.h"
#include "AmberRestrt.h"
#include "UnitParameter.h"
#include "LJAmber94.h"
#include "Timer.h"

namespace {

typedef std::vector<amberfile::AmberPrmtop::int_type> intvec_t;

void PushBackBond(const amberfile::AmberPrmtop& prmtop,
                  const intvec_t& v, 
                  CovalentBondList* bond_list) {
  for (intvec_t::size_type i = 0; i < v.size(); i += 3) {
    CovalentBondInfo::Bond b;
    b.id_of_atom[0] = v[i]/3;
    b.id_of_atom[1] = v[i+1]/3;
    b.typeofbond = v[i+2] - 1;
    bond_list->bond.push_back(b);
  }
}
void PushBackBondParameter(const amberfile::AmberPrmtop& prmtop,
                           CovalentBondParameterList* parameter_list) {
  parameter_list->bond.reserve(prmtop.numbnd);
  for (amberfile::AmberPrmtop::int_type i = 0; i < prmtop.numbnd; ++i) {
    CovalentBondInfo::BondParameter p;
    p.force_constant = UnitParameter::normalizeEnergy_kcal_mol(prmtop.rk[i]);
    //    p.force_constant = prmtop.rk[i];
    p.equilibrium_length = prmtop.req[i];
    parameter_list->bond.push_back(p);
  }
}

void PushBackAngle(const amberfile::AmberPrmtop& prmtop,
                   const intvec_t& v,
                   CovalentBondList* bond_list) {
  for (intvec_t::size_type i = 0; i < v.size(); i += 4) {
    CovalentBondInfo::Angle a;
    a.id_of_atom[0] = v[i]/3;
    a.id_of_atom[1] = v[i+1]/3;
    a.id_of_atom[2] = v[i+2]/3;
    a.typeofangle = v[i+3] - 1;
    bond_list->angle.push_back(a);
  }
}
void PushBackAngleParameter(const amberfile::AmberPrmtop& prmtop,
                            CovalentBondParameterList* parameter_list) {
  parameter_list->angle.reserve(prmtop.numang);
  for (amberfile::AmberPrmtop::int_type i = 0; i < prmtop.numang; ++i) {
    CovalentBondInfo::AngleParameter p;
    p.force_constant = UnitParameter::normalizeEnergy_kcal_mol(prmtop.tk[i]);
    p.equilibrium_angle = prmtop.teq[i];
    parameter_list->angle.push_back(p);
  }
}

void PushBackTorsionAndImproper(const amberfile::AmberPrmtop& prmtop,
                                const intvec_t& v,
                                CovalentBondList* bond_list,
                                CovalentBondParameterList* parameter_list) {
  for (intvec_t::size_type i = 0; i < v.size(); i += 5) {
    if (v[i+3] < 0) {  // v[n+3]<0: improper, otherwise: torsion
      CovalentBondInfo::Improper imp;
      imp.id_of_atom[0] = v[i]/3;
      imp.id_of_atom[1] = v[i+1]/3;   // v[i+2]<0 means end group to be ignored
      imp.id_of_atom[2] = std::abs(static_cast<int>(v[i+2]))/3;
      imp.id_of_atom[3] = -v[i+3]/3;  // v[i+3] is always negative here
      imp.typeofimproper = parameter_list->improper.size();   // v[i+4]-1;
      bond_list->improper.push_back(imp);
      intvec_t::size_type ic = v[i+4]-1;
      do {
        CovalentBondInfo::ImproperParameter p;
        p.force_constant = UnitParameter::normalizeEnergy_kcal_mol(prmtop.pk[ic]);
        p.periodicity = prmtop.pn[ic];
        p.phase = prmtop.phase[ic];
        parameter_list->improper.push_back(p);
      } while (prmtop.pn[ic++] < 0);
    } else {
      CovalentBondInfo::Torsion tor;
      tor.id_of_atom[0] = v[i]/3;
      tor.id_of_atom[1] = v[i+1]/3;
      tor.id_of_atom[2] = abs(v[i+2])/3;  // v[i+2]<0 means !!!
      tor.id_of_atom[3] = v[i+3]/3;  // v[i+3] is always >=0 here
      tor.typeoftorsion = parameter_list->torsion.size();     // v[i+4]-1
      bond_list->torsion.push_back(tor);
      intvec_t::size_type ic = v[i+4]-1;
      do {
        CovalentBondInfo::TorsionParameter p;
        p.force_constant = UnitParameter::normalizeEnergy_kcal_mol(prmtop.pk[ic]);
        p.periodicity = prmtop.pn[ic];
        p.phase = prmtop.phase[ic];
        parameter_list->torsion.push_back(p);
      } while (prmtop.pn[ic++] < 0);
    }
  }
}

}  // namespace

template<class PA>
int ReadAmberFile(const std::string& prmtop_file,
                  const std::string& restrt_file,
                  PA* particle_array,
                  CovalentBondList* bond_list,
                  CovalentBondParameterList* parameter_list,
                  SpaceVector<double>* box,
                  SpaceVector<double>* skew,
                  int restart,
                  double *start_time) {
  using namespace amberfile;

  //  double read_top_start = getrealtime();
  AmberPrmtop prmtop;
  ReadFromFile(prmtop_file, &prmtop);
  //  double read_top_end = getrealtime();
  //  printf("read prmtop %g sec\n",read_top_end-read_top_start);

  //  double read_rst_start = getrealtime();
  AmberRestrt restrt;
  ReadFromFile(restrt_file, &restrt, restart);
  //  double read_rst_end = getrealtime();
  //  printf("read restrt %g sec\n",read_rst_end-read_rst_start);

  *start_time = restrt.current_time;

  // construct particle array
  particle_array->clear();
  particle_array->reserve(prmtop.natom);
  for (amberfile::AmberPrmtop::int_type i = 0; i < prmtop.natom; ++i) {
    Particle p;
    p.position.x = restrt.coordinates[3*i    ];
    p.position.y = restrt.coordinates[3*i + 1];
    p.position.z = restrt.coordinates[3*i + 2];
    p.velocity.x = UnitParameter::normalizeVelocity(restrt.velocities[3*i    ]*2045.5);
    p.velocity.y = UnitParameter::normalizeVelocity(restrt.velocities[3*i + 1]*2045.5);
    p.velocity.z = UnitParameter::normalizeVelocity(restrt.velocities[3*i + 2]*2045.5);
    p.force = Force();
    p.charge = prmtop.chrg[i] / 18.2223;
    p.mass = prmtop.amass[i] / UnitParameter::CarbonWeight;
    p.inv_mass = 1.0 / p.mass;
    intvec_t::size_type at = 0;
    {
      do{
        if(prmtop.isymbl[i]==ljamber94parameters[at].name)break;
        at++;
      }while(ljamber94parameters[at].name!=LJAmberParameterEP.name);
    }
    //! rescan unmached atom name by 1st character
    /*!
      ex. "CY" is not in ljamber94parameters, it match to "C".
      If "CY" must match other "C?", append "CY" with same as "C?" 
      to LJAmber94.h
      */
    if(ljamber94parameters[at].name==LJAmberParameterEP.name){
      at = 0;
      std::string name1(prmtop.isymbl[i].substr(0,1));
      do{
        if(name1==ljamber94parameters[at].name)break;
        at++;
      }while(ljamber94parameters[at].name!=LJAmberParameterEP.name);
      if(ljamber94parameters[at].name==LJAmberParameterEP.name){
        std::cout << " Warning : atom name " << prmtop.isymbl[i] << " not found in ljamber94parameters" << std::endl;
      }else{
        if(DebugLog::verbose > 1){
          std::cout << " atom name " << prmtop.isymbl[i] << " match as " << name1 << std::endl;
        }
      }
    }
    p.atomtype = at;
    p.atomid = i;
    particle_array->push_back(p);
  }

  // DEBUG code  
#if 0
  {
    int ljsize = sizeof(ljamber94parameters)/sizeof(LJAmberParameter);
    std::vector<int> atomtype_dist(ljsize,0);
    for(int p=0;p<particle_array->size();p++){
      atomtype_dist[getatomtype((*particle_array),p)]++;
    }
    for(int a=0;a<ljsize;a++){
      if(atomtype_dist[a]>0){
        std::cout << " found atomtype " << a << " " << ljamber94parameters[a].name << " " << atomtype_dist[a] << std::endl;
      }
    }
  }
#endif

  // construct covalent bond list
  bond_list->clear();
  parameter_list->clear();
  //// bond
  bond_list->bond.reserve(prmtop.bh.size() + prmtop.b.size());
  PushBackBond(prmtop, prmtop.bh, bond_list);
  PushBackBond(prmtop, prmtop.b,  bond_list);
  PushBackBondParameter(prmtop, parameter_list);
  //// angle
  bond_list->angle.reserve(prmtop.th.size() + prmtop.t.size());
  PushBackAngle(prmtop, prmtop.th, bond_list);
  PushBackAngle(prmtop, prmtop.t,  bond_list);
  PushBackAngleParameter(prmtop, parameter_list);
  //// torsion & improper
  PushBackTorsionAndImproper(prmtop, prmtop.ph, bond_list, parameter_list);
  PushBackTorsionAndImproper(prmtop, prmtop.p,  bond_list, parameter_list);
  //// box
  box->x = restrt.box[0]; box->y = restrt.box[1]; box->z = restrt.box[2];
  //// skew
  skew->x = restrt.box[3]; skew->y = restrt.box[4]; skew->z = restrt.box[5];


  //! set bond_list.torsion.calc14interaction
  intvec_t exclude_offset;
  intvec_t::value_type extotal=0;
  for(typename PA::size_type i=0;i<particle_array->size();i++){
    exclude_offset.push_back(extotal);
    extotal += prmtop.numex[i];
  }
  if(DebugLog::verbose > 1){
    std::cout << " number of exclude pair " << extotal << std::endl;
  }
  std::vector<int> excluded_flag(extotal,0);
  intvec_t::value_type assigned_ex=0;
  //! check excluded by bond
  for(std::vector<CovalentBondInfo::Bond>::size_type b=0;
      b<bond_list->bond.size();b++){
    amberfile::AmberPrmtop::int_type b0 = bond_list->bond[b].id_of_atom[0];
    amberfile::AmberPrmtop::int_type b1 = bond_list->bond[b].id_of_atom[1];
    if(b0>b1)std::swap(b0,b1);
    intvec_t::value_type ep;
    for(ep=exclude_offset[b0];ep<exclude_offset[b0]+prmtop.numex[b0];ep++){
      if(b1==(prmtop.natex[ep]-1)){
        excluded_flag[ep] = 1;
        assigned_ex++;
        break;
      }
    }
    if(ep==exclude_offset[b0]+prmtop.numex[b0]){
      std::cout << " Warnig : (" << b0 << "," << b1 <<
          ") is not found in Amber Exclude list." << std::endl;
    }
  }
  //! check excluded by angle
  for(std::vector<CovalentBondInfo::Angle>::size_type a=0;
      a<bond_list->angle.size();a++){
    amberfile::AmberPrmtop::int_type a0 = bond_list->angle[a].id_of_atom[0];
    amberfile::AmberPrmtop::int_type a2 = bond_list->angle[a].id_of_atom[2];
    if(a0>a2)std::swap(a0,a2);
    intvec_t::value_type ep;
    for(ep=exclude_offset[a0];ep<exclude_offset[a0]+prmtop.numex[a0];ep++){
      if(a2==(prmtop.natex[ep]-1)){
        if(excluded_flag[ep]!=0){
          std::cout << " Warnig : exclude (" << a0 << "," << a2 <<
              ") duplicate by angle. 3 or 4-membered ring?" << std::endl;
        }else{
          assigned_ex++;
        }
        excluded_flag[ep] = 1;
        break;
      }
    }
    if(ep==exclude_offset[a0]+prmtop.numex[a0]){
      std::cout << " Warnig : (" << a0 << "," << a2 << ") is not found in Amber Exclude list." << std::endl;
    }
  }
  //! check excluded by torsion and set calc14interaction
  for(std::vector<CovalentBondInfo::Torsion>::size_type t=0;
      t<bond_list->torsion.size();t++){
    amberfile::AmberPrmtop::int_type t0 = bond_list->torsion[t].id_of_atom[0];
    amberfile::AmberPrmtop::int_type t3 = bond_list->torsion[t].id_of_atom[3];
    if(t0>t3)std::swap(t0,t3);
    intvec_t::value_type ep;
    for(ep=exclude_offset[t0];ep<exclude_offset[t0]+prmtop.numex[t0];ep++){
      if(t3==(prmtop.natex[ep]-1)){
        if(excluded_flag[ep]!=0){
          bond_list->torsion[t].calc14interaction = false;
        }else{
          bond_list->torsion[t].calc14interaction = true;
          assigned_ex++;
        }
        excluded_flag[ep] = 1;
        break;
      }
    }
    if(ep==exclude_offset[t0]+prmtop.numex[t0]){
      std::cout << " Warnig : (" << t0 << "," << t3 <<
          ") is not found in Amber Exclude list(" << exclude_offset[t0] <<
          "+" << prmtop.numex[t0] << ")." << std::endl;
    }
  }
  if(DebugLog::verbose > 2){
    std::cout << " number of bond assigned to exclude pair ";
    for(std::vector<int>::size_type i=0;i<excluded_flag.size();i++){
      std::cout << " " << excluded_flag[i];
    }
    std::cout << std::endl;
  }

  return 0;   // 0: success
}

#ifdef OLDPARTICLE
template
int ReadAmberFile(const std::string& prmtop_file,
                  const std::string& restrt_file,
                  ParticleArray* particle_array,
                  CovalentBondList* bond_list,
                  CovalentBondParameterList* parameter_list,
                  SpaceVector<double>* box,
                  SpaceVector<double>* skew,
                  int restart,
                  double *start_time);
#else
template
int ReadAmberFile(const std::string& prmtop_file,
                  const std::string& restrt_file,
                  CombinedParticleArray* particle_array,
                  CovalentBondList* bond_list,
                  CovalentBondParameterList* parameter_list,
                  SpaceVector<double>* box,
                  SpaceVector<double>* skew,
                  int restart,
                  double *start_time);
#endif
