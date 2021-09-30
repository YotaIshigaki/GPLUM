#include "HDFSnapshotData.h"

using namespace HDFSnapshot;

HDFDump::HDFDump(const int _mode)
  : snapshot(NULL), mode(_mode), num_atom(0)
{
}

HDFDump::HDFDump(const std::string& filename, const int _mode)
  : snapshot(NULL), mode(_mode), num_atom(0)
{
  init(filename);
}

HDFDump::~HDFDump()
{
  if(snapshot!=NULL)delete snapshot;
}

void HDFDump::dumpmode()
{
  mode = HDFDUMP;
}

void HDFDump::restoremode()
{
  mode = HDFRESTORE;
}

int HDFDump::init(const std::string& filename)
{
  snapshot = new HDFSnapshotFile(filename);
  std::cout << "HDFSnapshotFile " << filename << std::endl;
  return 1-snapshot->exists();
}

int HDFDump::close()
{
  if(snapshot!=NULL)delete snapshot;
  return 0;
}

void HDFDump::getParameter(double& deltat, int& step,
			   int& reduction_counter, int& print_counter,
			   int& crd_counter, int& rst_counter,
			   AtomID& size_of_particle)
{
  SnapshotDatatype *sd = snapshot->getSnapshotDatatype();
  deltat = sd->TimeStepInterval;
  step = sd->CurrentTimeStep;
  reduction_counter = sd->CurrentReductionCounter;
  print_counter = sd->CurrentPrintCounter;
  crd_counter = sd->CurrentCRDCounter;
  rst_counter = sd->CurrentRSTCounter;
  size_of_particle = sd->AtomArrayLength;
}

void HDFDump::setParameter(double deltat, int step,
			   int reduction_counter, int print_counter,
			   int crd_counter, int rst_counter,
			   AtomID size_of_particle)
{
  SnapshotGroup *sg = snapshot->getSnapshotGroup();
  sg->setParam(deltat,(long long)step,
	       (long long)reduction_counter,(long long)print_counter,
	       (long long)crd_counter,(long long)rst_counter,
	       (long long)size_of_particle);
}

void HDFDump::getInfo(int& rank, int& num_rank, string& version, AtomID& size_of_particle, HDFSnapshot::HDFSnapshotFile::SnapshotType& snapshot_type)
{
  rank = snapshot->getRank();
  num_rank = snapshot->getNumberOfRank();
  version = snapshot->getProgramVersion();
  size_of_particle = snapshot->getNumberOfAtom();
  snapshot_type = snapshot->getSnapshotType();
}

void HDFDump::setInfo(int rank, int num_rank, string* version, AtomID size_of_particle, HDFSnapshot::HDFSnapshotFile::SnapshotType snapshot_type)
{
  snapshot->setRank(rank,num_rank);
  snapshot->setProgramVersion(version);
  snapshot->setNumberOfAtom(size_of_particle);
  snapshot->setSnapshotType(snapshot_type);
}

void HDFDump::getIntegration(double& LJEnergyCorrection,
			     double& KineticEnergy,
			     double& PotentialEnergy,
			     double& TotalEnergy,
			     double& Virial,
			     double& EtaPosition,
			     double& EtaVelocity,
			     double& EtaForce,
			     double& LogVPosition,
			     double& LogVVelocity,
			     double& LogVForce,
			     double& Volume )
{
  IntegrationDatatype *integ = snapshot->getIntegrationDatatype();
  LJEnergyCorrection = integ->LJEnergyCorrection;
  KineticEnergy      = integ->KineticEnergy     ;
  PotentialEnergy    = integ->PotentialEnergy   ;
  TotalEnergy        = integ->TotalEnergy       ;
  Virial             = integ->Virial            ;
  EtaPosition        = integ->EtaPosition       ;
  EtaVelocity        = integ->EtaVelocity       ;
  EtaForce           = integ->EtaForce          ;
  LogVPosition       = integ->LogVPosition      ;
  LogVVelocity       = integ->LogVVelocity      ;
  LogVForce          = integ->LogVForce         ;
  Volume             = integ->Volume            ;
}

void HDFDump::setIntegration(double LJEnergyCorrection,
			     double KineticEnergy,
			     double PotentialEnergy,
			     double TotalEnergy,
			     double Virial,
			     double EtaPosition,
			     double EtaVelocity,
			     double EtaForce,
			     double LogVPosition,
			     double LogVVelocity,
			     double LogVForce,
			     double Volume )
{
  IntegrationDatatype *integ = snapshot->getIntegrationDatatype();
  integ->LJEnergyCorrection = LJEnergyCorrection;
  integ->KineticEnergy = KineticEnergy;
  integ->PotentialEnergy = PotentialEnergy;
  integ->TotalEnergy = TotalEnergy;
  integ->Virial = Virial;
  integ->EtaPosition = EtaPosition;
  integ->EtaVelocity = EtaVelocity;
  integ->EtaForce = EtaForce;
  integ->LogVPosition = LogVPosition;
  integ->LogVVelocity = LogVVelocity;
  integ->LogVForce = LogVForce;
  integ->Volume = Volume;
}

void HDFDump::setBoxDatatype(SpaceVector<double> boxsize, int type, SpaceVector<double> angle)
{
  BoxDatatype *boxdatatype = snapshot->getBoxDatatype();
  boxdatatype->SizeX = boxsize.x;
  boxdatatype->SizeY = boxsize.y;
  boxdatatype->SizeZ = boxsize.z;
  boxdatatype->Type = type;
  boxdatatype->BetaAngle1 = angle.x;
  boxdatatype->BetaAngle2 = angle.y;
  boxdatatype->BetaAngle3 = angle.z;
}


void HDFDump::getBoxDatatype(SpaceVector<double>& boxsize, int& type, SpaceVector<double>& angle)
{
  BoxDatatype *boxdatatype = snapshot->getBoxDatatype();
  boxsize.x = boxdatatype->SizeX;
  boxsize.y = boxdatatype->SizeY;
  boxsize.z = boxdatatype->SizeZ;
  type = boxdatatype->Type;
  angle.x = boxdatatype->BetaAngle1;
  angle.y = boxdatatype->BetaAngle2;
  angle.z = boxdatatype->BetaAngle3;
}

template<class PA>
void HDFDump::getParticle(PA& particlearray,
			  std::vector<TypeRange>& typerangearray)
{
  TyperangeDataset *typerangeset = snapshot->getTyperangeDataset();
  long long num_tr = typerangeset->getNumberOfTyperange();
  typerangearray.resize(num_tr);
  TyperangeDatatype *typeranges = typerangeset->getTyperangeDatatype();
  AtomID num_atom = 0;
  AtomID max_index = 0;
  AtomID min_index = typeranges[0].Begin;
  for(size_t tr=0;tr<num_tr;tr++){
    typerangearray[tr].begin = typeranges[tr].Begin;
    typerangearray[tr].end = typeranges[tr].End;
    typerangearray[tr].lj.begin = typeranges[tr].LJBegin;
    typerangearray[tr].lj.end = typeranges[tr].LJEnd;
    typerangearray[tr].ljcoulomb.begin = typeranges[tr].LJCoulombBegin;
    typerangearray[tr].ljcoulomb.end = typeranges[tr].LJCoulombEnd;
    typerangearray[tr].coulomb.begin = typeranges[tr].CoulombBegin;
    typerangearray[tr].coulomb.end = typeranges[tr].CoulombEnd;
    num_atom += (typeranges[tr].End - typeranges[tr].Begin);
    if(typeranges[tr].End>max_index){
      max_index = typeranges[tr].End;
    }
    if(typeranges[tr].Begin<min_index){
      min_index = typeranges[tr].Begin;
    }
    //    std::cout << "TypeRange " << tr << " : " <<  typerangearray[tr].begin << ":" << typerangearray[tr].end << std::endl;
  }
  std::cout << "Number of atom (from typerange) " << num_atom << std::endl;
  std::cout << "Max index of particle array (from typerange) " << max_index << std::endl;

  if(max_index>particlearray.size()){
    particlearray.resize(max_index);
  }
  AtomDataset *atomset = snapshot->getAtomDataset();
  //  atomset->setNumberOfAtom(num_atom);
  AtomDatatype *atoms = atomset->getAtomDatatype();
  long long na=0;
  for(size_t tr=0;tr<num_tr;tr++){
    for(int i=typerangearray[tr].begin;i<typerangearray[tr].end;i++){
      getatomid(particlearray,i)     = atoms[na].AtomId     ;
      getpos(particlearray,i).x      = atoms[na].PositionX  ;
      getpos(particlearray,i).y      = atoms[na].PositionY  ;
      getpos(particlearray,i).z      = atoms[na].PositionZ  ;
      getcharge(particlearray,i)     = atoms[na].Charge     ;
      getvelocity(particlearray,i).x = atoms[na].VelocityX  ;
      getvelocity(particlearray,i).y = atoms[na].VelocityY  ;
      getvelocity(particlearray,i).z = atoms[na].VelocityZ  ;
      getforce(particlearray,i).x    = atoms[na].ForceX     ;
      getforce(particlearray,i).y    = atoms[na].ForceY     ;
      getforce(particlearray,i).z    = atoms[na].ForceZ     ;
      getmass(particlearray,i)       = atoms[na].Mass       ;
      getinvmass(particlearray,i)    = atoms[na].InverseMass;
      getatomtype(particlearray,i)   = atoms[na].AtomType   ;
      na++;
    }
  }
  //  std::cout << "First(" << min_index << ") atom " << getatomid(particlearray,min_index) << " " << getpos(particlearray,min_index) << std::endl;
  //  std::cout << "Second(" << min_index+1 << ") atom " << getatomid(particlearray,min_index+1) << " " << getpos(particlearray,min_index+1) << std::endl;
  //  std::cout << "Last(" << max_index-1 << ") atom " << getatomid(particlearray,max_index-1) << " " << getpos(particlearray,max_index-1) << std::endl;
}

template
void HDFDump::getParticle(CombinedParticleArray& particlearray,
			  std::vector<TypeRange>& typerangearray);

template<class PA>
void HDFDump::setParticle(const PA& particlearray,
			  const std::vector<TypeRange>& typerangearray)
{
  TyperangeDataset *typerangeset = snapshot->getTyperangeDataset();
  size_t num_tr = typerangearray.size();
  typerangeset->setNumberOfTyperange((long long)num_tr);

  TyperangeDatatype *typeranges = typerangeset->getTyperangeDatatype();
  long long num_atom=0;
  for(size_t tr=0;tr<num_tr;tr++){
    typeranges[tr].Begin = typerangearray[tr].begin;
    typeranges[tr].End = typerangearray[tr].end;
    typeranges[tr].LJBegin = typerangearray[tr].lj.begin;
    typeranges[tr].LJEnd = typerangearray[tr].lj.end;
    typeranges[tr].LJCoulombBegin = typerangearray[tr].ljcoulomb.begin;
    typeranges[tr].LJCoulombEnd = typerangearray[tr].ljcoulomb.end;
    typeranges[tr].CoulombBegin = typerangearray[tr].coulomb.begin;
    typeranges[tr].CoulombEnd = typerangearray[tr].coulomb.end;
    num_atom += (typeranges[tr].End - typeranges[tr].Begin);
  }

  AtomDataset *atomset = snapshot->getAtomDataset();
  atomset->setNumberOfAtom(num_atom);
  AtomDatatype *atoms = atomset->getAtomDatatype();
  long long na=0;
  for(size_t tr=0;tr<num_tr;tr++){
    for(int i=typerangearray[tr].begin;i<typerangearray[tr].end;i++){
      atoms[na].AtomId = getatomid(particlearray,i);
      atoms[na].PositionX = getpos(particlearray,i).x;
      atoms[na].PositionY = getpos(particlearray,i).y;
      atoms[na].PositionZ = getpos(particlearray,i).z;
      atoms[na].Charge = getcharge(particlearray,i);
      atoms[na].VelocityX = getvelocity(particlearray,i).x;
      atoms[na].VelocityY = getvelocity(particlearray,i).y;
      atoms[na].VelocityZ = getvelocity(particlearray,i).z;
      atoms[na].ForceX = getforce(particlearray,i).x;
      atoms[na].ForceY = getforce(particlearray,i).y;
      atoms[na].ForceZ = getforce(particlearray,i).z;
      atoms[na].Mass = getmass(particlearray,i);
      atoms[na].InverseMass = getinvmass(particlearray,i);
      atoms[na].AtomType = getatomtype(particlearray,i);
      na++;
    }
  }
}

template
void HDFDump::setParticle(const CombinedParticleArray& particlearray,
			  const std::vector<TypeRange>& typerangearray);

void HDFDump::getBondlists(std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  BondlistarrayGroup *bondlists = snapshot->getBondlistarrayGroup();
  size_t num_bondlist = bondlists->getNumberOfBondlist();
  bondlistarray.resize(num_bondlist);
  //  std::cout << "bondlist " << num_bondlist << std::endl;
  for(int b=0;b<num_bondlist;b++){
    BondlistGroup *bondlist = bondlists->getBondlistGroup(b);
    BondNumberDatatype *bondnumberDatatype = bondlist->getBondNumberDataset()->getBondNumberDatatype();
    size_t num_bond     = bondnumberDatatype->NumberOfBond    ;
    size_t num_angle    = bondnumberDatatype->NumberOfAngle   ;
    size_t num_torsion  = bondnumberDatatype->NumberOfTrosion ;
    size_t num_improper = bondnumberDatatype->NumberOfImproper;
    bondlistarray[b].BondArray.resize(num_bond);	
    bondlistarray[b].AngleArray.resize(num_angle);	
    bondlistarray[b].TorsionArray.resize(num_torsion);	
    bondlistarray[b].ImproperArray.resize(num_improper);
    
    BondDataset *bondset = bondlist->getBondDataset();
    if(num_bond>0){
      BondDatatype *bond = bondset->getBondDatatype();
      for(size_t i=0;i<num_bond;i++){
	(bondlistarray[b].BondArray[i].id_of_atom[0]) = bond[i].AtomId0   ;
	(bondlistarray[b].BondArray[i].id_of_atom[1]) = bond[i].AtomId1   ;
	(bondlistarray[b].BondArray[i].typeofbond)    = bond[i].TypeOfBond;
	(bondlistarray[b].BondArray[i].shake)         = bond[i].Shake     ;
      }
    }
    AngleDataset *angleset = bondlist->getAngleDataset();
    if(num_angle>0){
      AngleDatatype *angle = angleset->getAngleDatatype();
      for(size_t i=0;i<num_angle;i++){
	(bondlistarray[b].AngleArray[i].id_of_atom[0]) = angle[i].AtomId0    ;
	(bondlistarray[b].AngleArray[i].id_of_atom[1]) = angle[i].AtomId1    ;
	(bondlistarray[b].AngleArray[i].id_of_atom[2]) = angle[i].AtomId2    ;
	(bondlistarray[b].AngleArray[i].typeofangle)   = angle[i].TypeOfAngle;
      }
    }
    TorsionDataset *torsionset = bondlist->getTorsionDataset();
    if(num_torsion>0){
      TorsionDatatype *torsion = torsionset->getTorsionDatatype();
      for(size_t i=0;i<num_torsion;i++){
	(bondlistarray[b].TorsionArray[i].id_of_atom[0]) = torsion[i].AtomId0      ;
	(bondlistarray[b].TorsionArray[i].id_of_atom[1]) = torsion[i].AtomId1      ;
	(bondlistarray[b].TorsionArray[i].id_of_atom[2]) = torsion[i].AtomId2      ;
	(bondlistarray[b].TorsionArray[i].id_of_atom[3]) = torsion[i].AtomId3      ;
	(bondlistarray[b].TorsionArray[i].typeoftorsion) = torsion[i].TypeOfTorsion;
	(bondlistarray[b].TorsionArray[i].calc14interaction) = (torsion[i].Calc14==1);
      }
    }
    ImproperDataset *improperset = bondlist->getImproperDataset();
    if(num_improper>0){
      ImproperDatatype *improper = improperset->getImproperDatatype();
      for(size_t i=0;i<num_improper;i++){
	(bondlistarray[b].ImproperArray[i].id_of_atom[0]) = improper[i].AtomId0;
	(bondlistarray[b].ImproperArray[i].id_of_atom[1]) = improper[i].AtomId1;
	(bondlistarray[b].ImproperArray[i].id_of_atom[2]) = improper[i].AtomId2;
	(bondlistarray[b].ImproperArray[i].id_of_atom[3]) = improper[i].AtomId3;
	(bondlistarray[b].ImproperArray[i].typeofimproper) = improper[i].TypeOfImproper;
      }
    }
  }
}

void HDFDump::setBondlists(const std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  BondlistarrayGroup *bondlists = snapshot->getBondlistarrayGroup();
  size_t num_bondlist = bondlistarray.size();
  bondlists->setNumberOfBondlist((long long)num_bondlist);
  //  std::cout << "bondlist " << num_bondlist << std::endl;
  for(int b=0;b<num_bondlist;b++){
    BondlistGroup *bondlist = bondlists->getBondlistGroup(b);
    size_t num_bond =bondlistarray[b].BondArray.size();
    size_t num_angle = bondlistarray[b].AngleArray.size();
    size_t num_torsion = bondlistarray[b].TorsionArray.size();
    size_t num_improper = bondlistarray[b].ImproperArray.size();
    BondNumberDatatype *bondnumberDatatype = bondlist->getBondNumberDataset()->getBondNumberDatatype();
    bondnumberDatatype->NumberOfBond = num_bond;
    bondnumberDatatype->NumberOfAngle = num_angle;
    bondnumberDatatype->NumberOfTrosion = num_torsion;
    bondnumberDatatype->NumberOfImproper = num_improper;
    BondDataset *bondset = bondlist->getBondDataset();
    bondset->setNumberOfBond((int)num_bond);
    if(num_bond>0){
      BondDatatype *bond = bondset->getBondDatatype();
      for(size_t i=0;i<num_bond;i++){
	bond[i].AtomId0 = (bondlistarray[b].BondArray[i].id_of_atom[0]);
	bond[i].AtomId1 = (bondlistarray[b].BondArray[i].id_of_atom[1]);
	bond[i].TypeOfBond = (bondlistarray[b].BondArray[i].typeofbond);
	bond[i].Shake = (bondlistarray[b].BondArray[i].shake);
      }
    }
    AngleDataset *angleset = bondlist->getAngleDataset();
    angleset->setNumberOfAngle((int)num_angle);
    AngleDatatype *angle = angleset->getAngleDatatype();
    for(size_t i=0;i<num_angle;i++){
      angle[i].AtomId0 = (bondlistarray[b].AngleArray[i].id_of_atom[0]);
      angle[i].AtomId1 = (bondlistarray[b].AngleArray[i].id_of_atom[1]);
      angle[i].AtomId2 = (bondlistarray[b].AngleArray[i].id_of_atom[2]);
      angle[i].TypeOfAngle = (bondlistarray[b].AngleArray[i].typeofangle);
    }
    TorsionDataset *torsionset = bondlist->getTorsionDataset();
    torsionset->setNumberOfTorsion((int)num_torsion);
    TorsionDatatype *torsion = torsionset->getTorsionDatatype();
    for(size_t i=0;i<num_torsion;i++){
      torsion[i].AtomId0 = (bondlistarray[b].TorsionArray[i].id_of_atom[0]);
      torsion[i].AtomId1 = (bondlistarray[b].TorsionArray[i].id_of_atom[1]);
      torsion[i].AtomId2 = (bondlistarray[b].TorsionArray[i].id_of_atom[2]);
      torsion[i].AtomId3 = (bondlistarray[b].TorsionArray[i].id_of_atom[3]);
      torsion[i].TypeOfTorsion = (bondlistarray[b].TorsionArray[i].typeoftorsion);
      int c14 = (bondlistarray[b].TorsionArray[i].calc14interaction ? 1 : 0);
      torsion[i].Calc14 = (c14);
    }
    ImproperDataset *improperset = bondlist->getImproperDataset();
    improperset->setNumberOfImproper((int)num_improper);
    ImproperDatatype *improper = improperset->getImproperDatatype();
    for(size_t i=0;i<num_improper;i++){
      improper[i].AtomId0 = (bondlistarray[b].ImproperArray[i].id_of_atom[0]);
      improper[i].AtomId1 = (bondlistarray[b].ImproperArray[i].id_of_atom[1]);
      improper[i].AtomId2 = (bondlistarray[b].ImproperArray[i].id_of_atom[2]);
      improper[i].AtomId3 = (bondlistarray[b].ImproperArray[i].id_of_atom[3]);
      improper[i].TypeOfImproper = (bondlistarray[b].ImproperArray[i].typeofimproper);
    }
  }
}

void HDFDump::getWaterlist(WaterList& waterlist)
{
  WaterDataset *waterset = snapshot->getWaterDataset();
  size_t num_water = waterset->getNumberOfWater();
  WaterDatatype *waters = waterset->getWaterDatatype();
  for(size_t i=0;i<num_water;i++){
    waterlist.add_water(waters[i].OWId0, waters[i].HWId1, waters[i].HWId2);
  }
}

void HDFDump::setWaterlist(const WaterList& waterlist)
{
  size_t num_water = waterlist.size();
  WaterDataset *waterset = snapshot->getWaterDataset();
  waterset->setNumberOfWater(num_water);
  WaterDatatype *waters = waterset->getWaterDatatype();
  size_t i=0;
  for(WaterList::const_iterator w=waterlist.begin();w!=waterlist.end();++w){
    waters[i].OWId0 = w->first;
    waters[i].HWId1 = w->second.h1;
    waters[i].HWId2 = w->second.h2;
    i++;
  }
}

void HDFDump::getShakelist(ShakeList& shakelist)
{
  ShakeDataset *shakeset = snapshot->getShakeDataset();
  size_t num_shake = shakeset->getNumberOfShake();
  ShakeDatatype *shakes = shakeset->getShakeDatatype();
  for(size_t i=0;i<num_shake;i++){
    int haid = shakes[i].HAId0;
    H1List h1list;
    h1list.nh1 = shakes[i].Nh1;
    h1list.h1[0] = shakes[i].H1Id1;
    h1list.h1[1] = shakes[i].H1Id2;
    h1list.h1[2] = shakes[i].H1Id3;
    h1list.h1[3] = shakes[i].H1Id4;
    h1list.bondtype[0] = shakes[i].BondType1;
    h1list.bondtype[1] = shakes[i].BondType2;
    h1list.bondtype[2] = shakes[i].BondType3;
    h1list.bondtype[3] = shakes[i].BondType4;
    shakelist.add_shake(haid,h1list);
  }
}

void HDFDump::setShakelist(const ShakeList& shakelist)
{
  size_t num_shake = shakelist.size();
  ShakeDataset *shakeset = snapshot->getShakeDataset();
  shakeset->setNumberOfShake(num_shake);
  ShakeDatatype *shakes = shakeset->getShakeDatatype();
  size_t i=0;
  for(ShakeList::const_iterator s=shakelist.begin();s!=shakelist.end();++s){
    shakes[i].HAId0 = s->first;
    shakes[i].Nh1 = s->second.nh1;
    shakes[i].H1Id1 = s->second.h1[0];
    shakes[i].H1Id2 = s->second.h1[1];
    shakes[i].H1Id3 = s->second.h1[2];
    shakes[i].H1Id4 = s->second.h1[3];
    shakes[i].BondType1 = s->second.bondtype[0];
    shakes[i].BondType2 = s->second.bondtype[1];
    shakes[i].BondType3 = s->second.bondtype[2];
    shakes[i].BondType4 = s->second.bondtype[3];
    i++;
  }
}

void HDFDump::dump()
{
  snapshot->write();
}

void HDFDump::restore(size_t step)
{
  //  snapshot->outputInfo();
  snapshot->read(step);
}
