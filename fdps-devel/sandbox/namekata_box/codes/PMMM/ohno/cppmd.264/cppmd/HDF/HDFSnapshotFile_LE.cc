/* @file HDFSnapshotFile.cpp
 * HDF Snapshot
 * @version 1.0.0
 */
#include "HDFSnapshotFile.h"
//#include "ErrorPos.h"
#include <sstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>

// from ErrorPos.h
#include <string>
#include <sstream>
#ifndef NDEBUG
#define errorPos(x) makeErrorPos(x,__FILE__,__LINE__,__FUNCTION__)
#else
#define errorPos(x) (x)
#endif
namespace {
const std::string makeErrorPos(const std::string s, char const *file, int line, char const* func)
{
  std::ostringstream os;
  os << s << ": " << file << ":" << line << ":" << func;
  return os.str();
}
}


using namespace std;

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

//#define SNAPSHOT_DEBUG

namespace HDFSnapshot {

void writeAttributeString( Group *group , string attrName , string attrValue ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### writeAttributeString() start #####" << endl;
#endif
    //
    DataSpace space( H5S_SCALAR );
    StrType datatype( PredType::C_S1 , H5T_VARIABLE );
    Attribute *attr = new Attribute( group->createAttribute( attrName , datatype , space ) );
    attr->write( datatype , attrValue );
    attr->close();
    delete attr;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### writeAttributeString() end   #####" << endl;
#endif
}

void updateAttributeString( Group *group , string attrName , string attrValue ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### updateAttributeString() start #####" << endl;
#endif
    //
    StrType datatype( PredType::C_S1 , H5T_VARIABLE );
    Attribute *attr = new Attribute( group->openAttribute( attrName ) );
    attr->write( datatype , attrValue );
    attr->close();
    delete attr;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### updateAttributeString() end   #####" << endl;
#endif
}

void readAttributeString( Group *group , string attrName , string *attrValue ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### readAttributeString() start #####" << endl;
#endif
    //
    StrType datatype( PredType::C_S1 , H5T_VARIABLE );
    Attribute *attr = new Attribute( group->openAttribute( attrName ) );
    attr->read( datatype , *attrValue );
    //cout << attrName << " : " << *attrValue << endl;
    attr->close();
    delete attr;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### readAttributeString() end   #####" << endl;
#endif
}

void writeAttribute( Group *group , string attrName , void* attrValue , PredType type ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### writeAttribute() start #####" << endl;
#endif
    //
    DataSpace space( H5S_SCALAR );
    StrType datatype( type );
    Attribute *attr = new Attribute( group->createAttribute( attrName , datatype , space ) );
    attr->write( datatype , attrValue );
    attr->close();
    delete attr;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### writeAttribute() end   #####" << endl;
#endif
}

void updateAttribute( Group *group , string attrName , void* attrValue , PredType type ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### updateAttribute() start #####" << endl;
#endif
    //
    StrType datatype( type );
    Attribute *attr = new Attribute( group->openAttribute( attrName ) );
    attr->write( datatype , attrValue );
    attr->close();
    delete attr;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### updateAttribute() end   #####" << endl;
#endif
}

void readAttribute( Group *group , string attrName , void* attrValue , PredType type ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### readAttributeDouble() start #####" << endl;
#endif
    //
    StrType datatype( type );
    Attribute *attr = new Attribute( group->openAttribute( attrName ) );
    attr->read( datatype , attrValue );
    attr->close();
    delete attr;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### readAttributeDouble() end   #####" << endl;
#endif
}

BoxDataset::BoxDataset() {
}

CompType* BoxDataset::createCompType() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### BoxDataset::createCompType() start ####" << endl;
#endif
    //
    CompType* datatype = new CompType( sizeof( BoxDatatype ) );
    //
    datatype->insertMember( "SizeX" , HOFFSET(BoxDatatype,SizeX) , PredType::IEEE_F64LE );
    datatype->insertMember( "SizeY" , HOFFSET(BoxDatatype,SizeY) , PredType::IEEE_F64LE );
    datatype->insertMember( "SizeZ" , HOFFSET(BoxDatatype,SizeZ) , PredType::IEEE_F64LE );
    datatype->insertMember( "Type" , HOFFSET(BoxDatatype,Type) , PredType::STD_I64LE );
    datatype->insertMember( "BetaAngle1" , HOFFSET(BoxDatatype,BetaAngle1) , PredType::IEEE_F64LE );
    datatype->insertMember( "BetaAngle2" , HOFFSET(BoxDatatype,BetaAngle2) , PredType::IEEE_F64LE );
    datatype->insertMember( "BetaAngle3" , HOFFSET(BoxDatatype,BetaAngle3) , PredType::IEEE_F64LE );
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### BoxDataset::createCompType() end   ####" << endl;
#endif
    //
    return datatype;
}

void BoxDataset::write( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### BoxDataset::write() start ####" << endl;
#endif
    //
    DataSpace space;
    CompType* datatype = createCompType();
    //
    DataSet* dataset = new DataSet( group->createDataSet( "BoxDataset" , *datatype , space ) );
    dataset->write( &boxDatatype , *datatype );
    delete dataset;
    space.close();
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### BoxDataset::write() end   ####" << endl;
#endif
}

void BoxDataset::read( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### BoxDataset::read() start ####" << endl;
#endif
    //
    CompType* datatype = createCompType();
    //
    DataSet* dataset = new DataSet( group->openDataSet( "BoxDataset" ) );
    //
    dataset->read( &boxDatatype , *datatype );
    delete dataset;
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### BoxDataset::read() end   ####" << endl;
#endif
}

AtomDataset::AtomDataset() : atomDatatype(0),numberOfAtom(0) {
}

CompType* AtomDataset::createCompType() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### AtomDataset::createCompType() start ####" << endl;
#endif
    //
    CompType* datatype = new CompType( sizeof( AtomDatatype ) );
    //
    datatype->insertMember( "AtomId"    , HOFFSET(AtomDatatype,AtomId)    , PredType::STD_I64LE  );
    datatype->insertMember( "PositionX" , HOFFSET(AtomDatatype,PositionX) , PredType::IEEE_F64LE );
    datatype->insertMember( "PositionY" , HOFFSET(AtomDatatype,PositionY) , PredType::IEEE_F64LE );
    datatype->insertMember( "PositionZ" , HOFFSET(AtomDatatype,PositionZ) , PredType::IEEE_F64LE );
    datatype->insertMember( "Charge"    , HOFFSET(AtomDatatype,Charge)    , PredType::IEEE_F64LE );
    datatype->insertMember( "VelocityX" , HOFFSET(AtomDatatype,VelocityX) , PredType::IEEE_F64LE );
    datatype->insertMember( "VelocityY" , HOFFSET(AtomDatatype,VelocityY) , PredType::IEEE_F64LE );
    datatype->insertMember( "VelocityZ" , HOFFSET(AtomDatatype,VelocityZ) , PredType::IEEE_F64LE );
    datatype->insertMember( "ForceX"    , HOFFSET(AtomDatatype,ForceX)    , PredType::IEEE_F64LE );
    datatype->insertMember( "ForceY"    , HOFFSET(AtomDatatype,ForceY)    , PredType::IEEE_F64LE );
    datatype->insertMember( "ForceZ"    , HOFFSET(AtomDatatype,ForceZ)    , PredType::IEEE_F64LE );
    datatype->insertMember( "Mass"      , HOFFSET(AtomDatatype,Mass)      , PredType::IEEE_F64LE );
    datatype->insertMember( "InverseMass" , HOFFSET(AtomDatatype,InverseMass) , PredType::IEEE_F64LE );
    datatype->insertMember( "AtomType"  , HOFFSET(AtomDatatype,AtomType)  , PredType::STD_I64LE  );
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### AtomDataset::createCompType() end   ####" << endl;
#endif
    //
    return datatype;
}

void AtomDataset::write( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### AtomDataset::write() start ####" << endl;
#endif
    //
    if( numberOfAtom > 0 ) {
        hsize_t dim[] = {numberOfAtom};
        DataSpace* space = new DataSpace(1,dim);
        CompType* datatype = createCompType();
        //
        DataSet* dataset = new DataSet( group->createDataSet( "AtomDataset" , *datatype , *space ) );
        //
        dataset->write( atomDatatype , *datatype );
        delete dataset;
        delete space;
        delete datatype;
    }
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### AtomDataset::write() end   ####" << endl;
#endif
}

void AtomDataset::read( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### AtomDataset::read() start ####" << endl;
#endif
    //
    CompType* datatype = createCompType();
    //
    DataSet* dataset = new DataSet( group->openDataSet( "AtomDataset" ) );
    DataSpace space = dataset->getSpace();
    //cout << "numberOfAtom : " << space.getSelectNpoints() << endl;
    setNumberOfAtom( space.getSelectNpoints() );
    space.close();
    //
    dataset->read( atomDatatype , *datatype );
    delete dataset;
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### AtomDataset::read() end   ####" << endl;
#endif
}

TyperangeDataset::TyperangeDataset() : typerangeDatatype(0),numberOfTyperange(0) {
}

CompType* TyperangeDataset::createCompType() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### TyperangeDataset::createCompType() start ####" << endl;
#endif
    //
    CompType* datatype = new CompType( sizeof( TyperangeDatatype ) );
    //
    datatype->insertMember( "Begin"   , HOFFSET(TyperangeDatatype,Begin)   , PredType::STD_I64LE  );
    datatype->insertMember( "End"     , HOFFSET(TyperangeDatatype,End)     , PredType::STD_I64LE  );
    datatype->insertMember( "LJBegin" , HOFFSET(TyperangeDatatype,LJBegin) , PredType::STD_I64LE  );
    datatype->insertMember( "LJEnd"   , HOFFSET(TyperangeDatatype,LJEnd)   , PredType::STD_I64LE  );
    datatype->insertMember( "LJCoulombBegin", HOFFSET(TyperangeDatatype,LJCoulombBegin), PredType::STD_I64LE  );
    datatype->insertMember( "LJCoulombEnd", HOFFSET(TyperangeDatatype,LJCoulombEnd), PredType::STD_I64LE  );
    datatype->insertMember( "CoulombBegin", HOFFSET(TyperangeDatatype,CoulombBegin), PredType::STD_I64LE  );
    datatype->insertMember( "CoulombEnd", HOFFSET(TyperangeDatatype,CoulombEnd), PredType::STD_I64LE  );
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### TyperangeDataset::createCompType() end   ####" << endl;
#endif
    //
    return datatype;
}

void TyperangeDataset::write( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### TyperangeDataset::write() start ####" << endl;
#endif
    //
    if( numberOfTyperange > 0 ) {
        hsize_t dim[] = {numberOfTyperange};
        DataSpace* space = new DataSpace(1,dim);
        CompType* datatype = createCompType();
        //
        DataSet* dataset = new DataSet( group->createDataSet( "TyperangeDataset" , *datatype , *space ) );
        //
        dataset->write( typerangeDatatype , *datatype );
        delete dataset;
        delete space;
        delete datatype;
    }
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### TyperangeDataset::write() end   ####" << endl;
#endif
}

void TyperangeDataset::read( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### TyperangeDataset::read() start ####" << endl;
#endif
    //
    CompType* datatype = createCompType();
    //
    DataSet* dataset = new DataSet( group->openDataSet( "TyperangeDataset" ) );
    DataSpace space = dataset->getSpace();
    //cout << "numberOfTyperange : " << space.getSelectNpoints() << endl;
    setNumberOfTyperange( space.getSelectNpoints() );
    space.close();
    //
    dataset->read( typerangeDatatype , *datatype );
    delete dataset;
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### TyperangeDataset::read() end   ####" << endl;
#endif
}

BondDataset::BondDataset() : bondDatatype(0), numberOfBond(0) {
}

CompType* BondDataset::createMemType() {
    CompType* datatype = new CompType( sizeof( BondDatatype ) );
    datatype->insertMember( "AtomId0" , HOFFSET(BondDatatype,AtomId0) , PredType::STD_I64LE );
    datatype->insertMember( "AtomId1" , HOFFSET(BondDatatype,AtomId1) , PredType::STD_I64LE );
    datatype->insertMember( "TypeOfBond" , HOFFSET(BondDatatype,TypeOfBond) , PredType::STD_I64LE );
    datatype->insertMember( "Shake" , HOFFSET(BondDatatype,Shake) , PredType::STD_I64LE );
    return datatype;
}

void BondDataset::read( Group *structureGroup ) {
    CompType* memtype = createMemType();
    DataSet* dataset = new DataSet( structureGroup->openDataSet( "BondDataset" ) );
    DataSpace space = dataset->getSpace();
    setNumberOfBond( space.getSelectNpoints() );
    space.close();
    //
    dataset->read( bondDatatype , *memtype );
    //
    delete dataset;
    delete memtype;
}

void BondDataset::write( Group *structureGroup ) {
    if (numberOfBond > 0) {
        hsize_t dim[] = {numberOfBond};
        DataSpace* space = new DataSpace(1,dim);
        CompType* memtype = createMemType();
        //
        DataSet* dataset = new DataSet( structureGroup->createDataSet( "BondDataset" , *memtype , *space ) );
        dataset->write( bondDatatype , *memtype );
        delete dataset;
        delete space;
        delete memtype;
    }
    //
}

AngleDataset::AngleDataset() {
    angleDatatype = 0;
    numberOfAngle = 0;
    //
}

CompType* AngleDataset::createMemType() {
    CompType* datatype = new CompType( sizeof( AngleDatatype ) );
    datatype->insertMember( "AtomId0" , HOFFSET(AngleDatatype,AtomId0) , PredType::STD_I64LE );
    datatype->insertMember( "AtomId1" , HOFFSET(AngleDatatype,AtomId1) , PredType::STD_I64LE );
    datatype->insertMember( "AtomId2" , HOFFSET(AngleDatatype,AtomId2) , PredType::STD_I64LE );
    datatype->insertMember( "TypeOfAngle" , HOFFSET(AngleDatatype,TypeOfAngle) , PredType::STD_I64LE );
    //
    return datatype;
}

void AngleDataset::read( Group *structureGroup ) {
    CompType* memtype = createMemType();
    DataSet* dataset = new DataSet( structureGroup->openDataSet( "AngleDataset" ) );
    DataSpace space = dataset->getSpace();
    setNumberOfAngle( space.getSelectNpoints() );
    space.close();
    //
    dataset->read( angleDatatype , *memtype );
    //
    delete dataset;
    delete memtype;
}

void AngleDataset::write( Group *structureGroup ) {
    //
    if (numberOfAngle > 0) {
        hsize_t dim[] = {numberOfAngle};
        DataSpace* space = new DataSpace(1,dim);
        CompType* memtype = createMemType();
        //
        DataSet* dataset = new DataSet( structureGroup->createDataSet( "AngleDataset" , *memtype , *space ) );
        dataset->write( angleDatatype , *memtype );
        delete dataset;
        delete space;
        delete memtype;
    }
    //
}


TorsionDataset::TorsionDataset() {
    torsionDatatype = 0;
    numberOfTorsion = 0;
}

CompType* TorsionDataset::createMemType() {
    //
    CompType* datatype = new CompType( sizeof( TorsionDatatype ) );
    datatype->insertMember( "AtomId0" , HOFFSET(TorsionDatatype,AtomId0) , PredType::STD_I64LE );
    datatype->insertMember( "AtomId1" , HOFFSET(TorsionDatatype,AtomId1) , PredType::STD_I64LE );
    datatype->insertMember( "AtomId2" , HOFFSET(TorsionDatatype,AtomId2) , PredType::STD_I64LE );
    datatype->insertMember( "AtomId3" , HOFFSET(TorsionDatatype,AtomId3) , PredType::STD_I64LE );
    datatype->insertMember( "TypeOfTorsion" , HOFFSET(TorsionDatatype,TypeOfTorsion) , PredType::STD_I64LE );
    datatype->insertMember( "Calc14" , HOFFSET(TorsionDatatype,Calc14) , PredType::STD_I64LE );
    //
    return datatype;
}

void TorsionDataset::read( Group *dihedralGroup ) {
    //
    CompType* memtype = createMemType();
    DataSet* dataset = new DataSet( dihedralGroup->openDataSet( "TorsionDataset" ) );
    DataSpace space = dataset->getSpace();
    setNumberOfTorsion( space.getSelectNpoints() );
    space.close();
    //
    dataset->read( torsionDatatype , *memtype );
    //
    delete dataset;
    delete memtype;
    //
}

void TorsionDataset::write( Group *dihedralGroup ) {
    //
    if (numberOfTorsion > 0) {
        hsize_t dim[] = {numberOfTorsion};
        DataSpace* space = new DataSpace(1,dim);
        CompType* memtype = createMemType();
        //
        DataSet* dataset = new DataSet( dihedralGroup->createDataSet( "TorsionDataset" , *memtype , *space ) );
        dataset->write( torsionDatatype , *memtype );
        delete dataset;
        delete space;
        delete memtype;
    }
    //
}

//

ImproperDataset::ImproperDataset() {
    //
    improperDatatype = 0;
    numberOfImproper = 0;
    //
}

CompType* ImproperDataset::createMemType() {
    //
    CompType* datatype = new CompType( sizeof( ImproperDatatype ) );
    datatype->insertMember( "AtomId0" , HOFFSET(ImproperDatatype,AtomId0) , PredType::STD_I64LE );
    datatype->insertMember( "AtomId1" , HOFFSET(ImproperDatatype,AtomId1) , PredType::STD_I64LE );
    datatype->insertMember( "AtomId2" , HOFFSET(ImproperDatatype,AtomId2) , PredType::STD_I64LE );
    datatype->insertMember( "AtomId3" , HOFFSET(ImproperDatatype,AtomId3) , PredType::STD_I64LE );
    datatype->insertMember( "TypeOfImproper" , HOFFSET(ImproperDatatype,TypeOfImproper) , PredType::STD_I64LE );
    //
    return datatype;
}

void ImproperDataset::read( Group *dihedralGroup ) {
    //
    CompType* memtype = createMemType();
    DataSet* dataset = new DataSet( dihedralGroup->openDataSet( "ImproperDataset" ) );
    DataSpace space = dataset->getSpace();
    setNumberOfImproper( space.getSelectNpoints() );
    space.close();
    //
    dataset->read( improperDatatype , *memtype );
    //
    delete dataset;
    delete memtype;
    //
}

void ImproperDataset::write( Group *dihedralGroup ) {
    //
    if (numberOfImproper > 0) {
        hsize_t dim[] = {numberOfImproper};
        DataSpace* space = new DataSpace(1,dim);
        CompType* memtype = createMemType();
        //
        DataSet* dataset = new DataSet( dihedralGroup->createDataSet( "ImproperDataset" , *memtype , *space ) );
        dataset->write( improperDatatype , *memtype );
        delete dataset;
        delete space;
        delete memtype;
    }
    //
}

BondlistGroup::BondlistGroup() {
  bondDataset = new BondDataset();
  angleDataset = new AngleDataset();
  torsionDataset = new TorsionDataset();
  improperDataset = new ImproperDataset();
}

void BondlistGroup::read(Group *bondlistarrayGroup){

  stringstream s;
  s << "BondlistGroup_" << index;
  Group* bondlistGroup = NULL;
  try{
    bondlistGroup = new Group( bondlistarrayGroup->openGroup(s.str().c_str()));
  } catch(GroupIException e ){
    e.printError();
    string msg = "Can't find group [/SnapshotGroup/BondlistarrayGroup/"+s.str()+"] on HDF snapshot file";
    throw std::runtime_error(errorPos(msg));
  }
  try {
    bondDataset->read(bondlistGroup);
  } catch( Exception e ) {
    e.printError();
    string msg = "Can't read dataset [BondDataset] on HDF restart file";
    throw std::runtime_error(errorPos(msg));
  }
  try {
    angleDataset->read(bondlistGroup);
  } catch( Exception e ) {
    e.printError();
    string msg = "Can't read dataset [AngleDataset] on HDF restart file";
    throw std::runtime_error(errorPos(msg));
  }
  try {
    torsionDataset->read(bondlistGroup);
  } catch( Exception e ) {
    e.printError();
    string msg = "Can't read dataset [TorsionDataset] on HDF restart file";
    throw std::runtime_error(errorPos(msg));
  }
  try {
    improperDataset->read(bondlistGroup);
  } catch( Exception e ) {
    e.printError();
    string msg = "Can't read dataset [ImproperDataset] on HDF restart file";
    throw std::runtime_error(errorPos(msg));
  }
  delete bondlistGroup;
}

void BondlistGroup::write(Group *bondlistarrayGroup){
  stringstream s;
  s << "BondlistGroup_" << index;
  Group* bondlistGroup = new Group( bondlistarrayGroup->createGroup(s.str().c_str()));
  bondDataset->write(bondlistGroup);
  angleDataset->write(bondlistGroup);
  torsionDataset->write(bondlistGroup);
  improperDataset->write(bondlistGroup);
  delete bondlistGroup;
}

BondlistarrayGroup::BondlistarrayGroup(){
  bondlistGroup = NULL;
  numberOfBondlist = 0;
}

void BondlistarrayGroup::read( Group *group){

  Group* bondlistarrayGroup = NULL;
  try {
    bondlistarrayGroup = new Group(group->openGroup("BondlistarrayGroup"));
  } catch (Exception e ) {
    e.printError();
    string msg = "Can't find group [/SnapshotGroup/BondlistarrayGroup] on HDF restart file";
    throw std::runtime_error(errorPos(msg));
  }
  int num = bondlistarrayGroup->getNumObjs();
  setNumberOfBondlist(num);
  try {
    for(int n=0;n<numberOfBondlist;n++){
      bondlistGroup[n].setIndex(n);
      bondlistGroup[n].read(bondlistarrayGroup);
    }
  } catch (Exception e ) {
    e.printError();
    string msg = "Can't read group [BondlistGroup] on HDF restart file";
    throw std::runtime_error(errorPos(msg));
  }
  delete bondlistarrayGroup;
}

void BondlistarrayGroup::write( Group *group){
  Group* bondlistarrayGroup = new Group(group->createGroup("BondlistarrayGroup"));
  try {
    for(int n=0;n<numberOfBondlist;n++){
      bondlistGroup[n].write(bondlistarrayGroup);
    }
  } catch( Exception e ) {
    e.printError();
    string msg = "Can't write group [BondlistGroup] on HDF restart file";
    throw std::runtime_error(errorPos(msg));
  }
  delete bondlistarrayGroup;
}


WaterDataset::WaterDataset() : waterDatatype(0),numberOfWater(0) {
}

CompType* WaterDataset::createCompType() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### WaterDataset::createCompType() start ####" << endl;
#endif
    //
    CompType* datatype = new CompType( sizeof( WaterDatatype ) );
    //
    datatype->insertMember( "OWId0" , HOFFSET(WaterDatatype,OWId0) , PredType::STD_I64LE  );
    datatype->insertMember( "HWId1" , HOFFSET(WaterDatatype,HWId1) , PredType::STD_I64LE  );
    datatype->insertMember( "HWId2" , HOFFSET(WaterDatatype,HWId2) , PredType::STD_I64LE  );
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### WaterDataset::createCompType() end   ####" << endl;
#endif
    //
    return datatype;
}

void WaterDataset::write( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### WaterDataset::write() start ####" << endl;
#endif
    //
    if( numberOfWater > 0 ) {
        hsize_t dim[] = {numberOfWater};
        DataSpace* space = new DataSpace(1,dim);
        CompType* datatype = createCompType();
        //
        DataSet* dataset = new DataSet( group->createDataSet( "WaterDataset" , *datatype , *space ) );
        //
        dataset->write( waterDatatype , *datatype );
        delete dataset;
        delete space;
        delete datatype;
    }
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### WaterDataset::write() end   ####" << endl;
#endif
}

void WaterDataset::read( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### WaterDataset::read() start ####" << endl;
#endif
    //
    CompType* datatype = createCompType();
    //
    DataSet* dataset = new DataSet( group->openDataSet( "WaterDataset" ) );
    DataSpace space = dataset->getSpace();
    //cout << "numberOfWater : " << space.getSelectNpoints() << endl;
    setNumberOfWater( space.getSelectNpoints() );
    space.close();
    //
    dataset->read( waterDatatype , *datatype );
    delete dataset;
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### WaterDataset::read() end   ####" << endl;
#endif
}



ShakeDataset::ShakeDataset() : shakeDatatype(0),numberOfShake(0) {
}

CompType* ShakeDataset::createCompType() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### ShakeDataset::createCompType() start ####" << endl;
#endif
    //
    CompType* datatype = new CompType( sizeof( ShakeDatatype ) );
    //
    datatype->insertMember( "HAId0"     , HOFFSET(ShakeDatatype,HAId0)     , PredType::STD_I64LE  );
    datatype->insertMember( "Nh1"       , HOFFSET(ShakeDatatype,Nh1)       , PredType::STD_I64LE  );
    datatype->insertMember( "H1Id1"     , HOFFSET(ShakeDatatype,H1Id1)     , PredType::STD_I64LE  );
    datatype->insertMember( "H1Id2"     , HOFFSET(ShakeDatatype,H1Id2)     , PredType::STD_I64LE  );
    datatype->insertMember( "H1Id3"     , HOFFSET(ShakeDatatype,H1Id3)     , PredType::STD_I64LE  );
    datatype->insertMember( "H1Id4"     , HOFFSET(ShakeDatatype,H1Id4)     , PredType::STD_I64LE  );
    datatype->insertMember( "BondType1" , HOFFSET(ShakeDatatype,BondType1) , PredType::STD_I64LE  );
    datatype->insertMember( "BondType2" , HOFFSET(ShakeDatatype,BondType2) , PredType::STD_I64LE  );
    datatype->insertMember( "BondType3" , HOFFSET(ShakeDatatype,BondType3) , PredType::STD_I64LE  );
    datatype->insertMember( "BondType4" , HOFFSET(ShakeDatatype,BondType4) , PredType::STD_I64LE  );
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### ShakeDataset::createCompType() end   ####" << endl;
#endif
    //
    return datatype;
}

void ShakeDataset::write( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### ShakeDataset::write() start ####" << endl;
#endif
    //
    if( numberOfShake > 0 ) {
        hsize_t dim[] = {numberOfShake};
        DataSpace* space = new DataSpace(1,dim);
        CompType* datatype = createCompType();
        //
        DataSet* dataset = new DataSet( group->createDataSet( "ShakeDataset" , *datatype , *space ) );
        //
        dataset->write( shakeDatatype , *datatype );
        delete dataset;
        delete space;
        delete datatype;
    }
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### ShakeDataset::write() end   ####" << endl;
#endif
}

void ShakeDataset::read( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### ShakeDataset::read() start ####" << endl;
#endif
    //
    CompType* datatype = createCompType();
    //
    DataSet* dataset = new DataSet( group->openDataSet( "ShakeDataset" ) );
    DataSpace space = dataset->getSpace();
    //cout << "numberOfShake : " << space.getSelectNpoints() << endl;
    setNumberOfShake( space.getSelectNpoints() );
    space.close();
    //
    dataset->read( shakeDatatype , *datatype );
    delete dataset;
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### ShakeDataset::read() end   ####" << endl;
#endif
}


IntegrationDataset::IntegrationDataset() {
}

CompType* IntegrationDataset::createCompType() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### IntegrationDataset::createCompType() start ####" << endl;
#endif
    //
    CompType* datatype = new CompType( sizeof( IntegrationDatatype ) );
    //
    datatype->insertMember( "LJEnergyCorrection" , HOFFSET(IntegrationDatatype,LJEnergyCorrection) , PredType::IEEE_F64LE );
    datatype->insertMember( "KineticEnergy" , HOFFSET(IntegrationDatatype,KineticEnergy) , PredType::IEEE_F64LE );
    datatype->insertMember( "PotentialEnergy" , HOFFSET(IntegrationDatatype,PotentialEnergy) , PredType::IEEE_F64LE );
    datatype->insertMember( "TotalEnergy" , HOFFSET(IntegrationDatatype,TotalEnergy) , PredType::IEEE_F64LE );
    datatype->insertMember( "Virial" , HOFFSET(IntegrationDatatype,Virial) , PredType::IEEE_F64LE );
    datatype->insertMember( "EtaPosition" , HOFFSET(IntegrationDatatype,EtaPosition) , PredType::IEEE_F64LE );
    datatype->insertMember( "EtaVelocity" , HOFFSET(IntegrationDatatype,EtaVelocity) , PredType::IEEE_F64LE );
    datatype->insertMember( "EtaForce" , HOFFSET(IntegrationDatatype,EtaForce) , PredType::IEEE_F64LE );
    datatype->insertMember( "LogVPosition" , HOFFSET(IntegrationDatatype,LogVPosition) , PredType::IEEE_F64LE );
    datatype->insertMember( "LogVVelocity" , HOFFSET(IntegrationDatatype,LogVVelocity) , PredType::IEEE_F64LE );
    datatype->insertMember( "LogVForce" , HOFFSET(IntegrationDatatype,LogVForce) , PredType::IEEE_F64LE );
    datatype->insertMember( "Volume" , HOFFSET(IntegrationDatatype,Volume) , PredType::IEEE_F64LE );
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### IntegrationDataset::createCompType() end   ####" << endl;
#endif
    //
    return datatype;
}

void IntegrationDataset::write( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### IntegrationDataset::write() start ####" << endl;
#endif
    //
    CompType* datatype = createCompType();
    //
    DataSpace space;
    DataSet*  dataset = new DataSet( group->createDataSet( "IntegrationDataset" , *datatype , space ) );
    dataset->write( &integrationDatatype , *datatype );
    delete dataset;
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### IntegrationDataset::write() end   ####" << endl;
#endif
}

void IntegrationDataset::read( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### IntegrationDataset::read() start ####" << endl;
#endif
    //
    CompType* datatype = createCompType();
    //
    DataSet* dataset = new DataSet( group->openDataSet( "IntegrationDataset" ) );
    dataset->read( &integrationDatatype , *datatype );
    delete dataset;
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### IntegrationDataset::read() end   ####" << endl;
#endif
}



SnapshotDataset::SnapshotDataset() {
}

CompType* SnapshotDataset::createCompType() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotDataset::createCompType() start ####" << endl;
#endif
    //
    CompType* datatype = new CompType( sizeof( SnapshotDatatype ) );
    //
    datatype->insertMember( "TimeStepInterval" , HOFFSET(SnapshotDatatype,TimeStepInterval) , PredType::IEEE_F64LE );
    datatype->insertMember( "CurrentTimeStep"  , HOFFSET(SnapshotDatatype,CurrentTimeStep)  , PredType::STD_I64LE );
    datatype->insertMember( "CurrentReductionCounter"  , HOFFSET(SnapshotDatatype,CurrentReductionCounter)  , PredType::STD_I64LE );
    datatype->insertMember( "CurrentPrintCounter"  , HOFFSET(SnapshotDatatype,CurrentPrintCounter)  , PredType::STD_I64LE );
    datatype->insertMember( "CurrentCRDCounter"  , HOFFSET(SnapshotDatatype,CurrentCRDCounter)  , PredType::STD_I64LE );
    datatype->insertMember( "CurrentRSTCounter"  , HOFFSET(SnapshotDatatype,CurrentRSTCounter)  , PredType::STD_I64LE );
    datatype->insertMember( "AtomArrayLength"  , HOFFSET(SnapshotDatatype,AtomArrayLength)  , PredType::STD_I64LE );
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotDataset::createCompType() end   ####" << endl;
#endif
    //
    return datatype;
}

void SnapshotDataset::write( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotDataset::write() start ####" << endl;
#endif
    //
    CompType* datatype = createCompType();
    //
    DataSpace space;
    DataSet*  dataset = new DataSet( group->createDataSet( "SnapshotDataset" , *datatype , space ) );
    dataset->write( &snapshotDatatype , *datatype );
    delete dataset;
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotDataset::write() end   ####" << endl;
#endif
}

void SnapshotDataset::read( Group *group ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotDataset::read() start ####" << endl;
#endif
    //
    CompType* datatype = createCompType();
    //
    DataSet* dataset = new DataSet( group->openDataSet( "SnapshotDataset" ) );
    dataset->read( &snapshotDatatype , *datatype );
    delete dataset;
    delete datatype;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotDataset::read() end   ####" << endl;
#endif
}

//-------------------------- SnapshotGroup -----------------------------

SnapshotGroup::SnapshotGroup() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotGroup() start #####" << endl;
#endif
    //
    snapshotDataset = new SnapshotDataset();
    integrationDataset = new IntegrationDataset();
    boxDataset = new BoxDataset();
    atomDataset = new AtomDataset();
    typerangeDataset = new TyperangeDataset();
    bondlistarrayGroup = new BondlistarrayGroup();
    waterDataset = new WaterDataset();
    shakeDataset = new ShakeDataset();
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotGroup() end #####" << endl;
#endif
}

void SnapshotGroup::write( Group* rootGroup ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotGroup::write() start #####" << endl;
#endif
    //
    stringstream s;
    s << snapshotDataset->getSnapshotDatatype()->CurrentTimeStep;
    //cout << "group name:" << s.str() << ":" << endl;
    Group* group = new Group( rootGroup->createGroup( s.str().c_str() ) );
    //
    try {
        snapshotDataset->write( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't write dataset [SnapshotDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
        integrationDataset->write( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't write dataset [IntegrationDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
        boxDataset->write( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't write dataset [BoxDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
        atomDataset->write( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't write dataset [AtomDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
        typerangeDataset->write( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't write dataset [TyperangeDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
      bondlistarrayGroup->write( group);
    } catch( Exception e ) {
      e.printError();
      string msg = "Can't write dataset [BondlistarrayGroup] on HDF snapshot file";
      throw std::runtime_error(errorPos(msg));
    }
    try {
        waterDataset->write( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't write dataset [WaterDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
        shakeDataset->write( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't write dataset [ShakeDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    //
    if( group != NULL ) delete group;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotGroup::write() end #####" << endl;
#endif
}

void SnapshotGroup::read( Group* rootGroup , size_t step ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotGroup::read() start #####" << endl;
#endif
    //
    stringstream s;
    s << step;
    //cout << "group name:" << s.str() << ":" << endl;
    Group* group = NULL;
    try {
        group = new Group( rootGroup->openGroup( s.str().c_str() ) );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't find group [SnapshotGroup:"+s.str()+"] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
        snapshotDataset->read( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't read dataset [SnapshotDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
        integrationDataset->read( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't read dataset [IntegrationDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
        boxDataset->read( group );
    } catch( Exception e ) {
        e.printError();
        string msg = "Can't read dataset [BoxDataset] on HDF snapshot file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
        atomDataset->read( group );
    } catch( Exception e ) {
    }
    try {
        typerangeDataset->read( group );
    } catch( Exception e ) {
    }
    try {
      bondlistarrayGroup->read( group);
    } catch( Exception e) {
      e.printError();
        string msg = "Can't read group [BondlistarrayGroup] on HDF restart file";
        throw std::runtime_error(errorPos(msg));
    }
    try {
      waterDataset->read( group );
    } catch( Exception e ) {
    }
    try {
      shakeDataset->read( group );
    } catch( Exception e ) {
    }
    //
    delete group;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotGroup::read() end #####" << endl;
#endif
}

void SnapshotGroup::setParam( double timeStepInterval, long long currentTimeStep, 
			      long long currentReductionCounter, long long currentPrintCounter,
			      long long currentCRDCounter, long long currentRSTCounter,
			      long long atomArrayLength ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotGroup::setParam() start #####" << endl;
#endif
    //cout << "currentTimeStep : " << currentTimeStep << endl;
    //
    snapshotDataset->getSnapshotDatatype()->TimeStepInterval = timeStepInterval;
    snapshotDataset->getSnapshotDatatype()->CurrentTimeStep = currentTimeStep;
    snapshotDataset->getSnapshotDatatype()->CurrentReductionCounter = currentReductionCounter;
    snapshotDataset->getSnapshotDatatype()->CurrentPrintCounter = currentPrintCounter;
    snapshotDataset->getSnapshotDatatype()->CurrentCRDCounter = currentCRDCounter;
    snapshotDataset->getSnapshotDatatype()->CurrentRSTCounter = currentRSTCounter;
    snapshotDataset->getSnapshotDatatype()->AtomArrayLength = atomArrayLength;
    //
    atomDataset->setNumberOfAtom( atomArrayLength );
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### SnapshotGroup::setParam() end #####" << endl;
#endif
}

void SnapshotGroup::outputInfo( Group* rootGroup ) {
    cout << "##### SnapshotGroup::outputInfo() start #####" << endl;
    //
    hsize_t numberOfStep = rootGroup->getNumObjs();
    cout << "Number Of Step : " << numberOfStep << endl;
    vector<long> steps;
    for( int n = 0 ; n < numberOfStep ; ++n ) {
        H5std_string name = rootGroup->getObjnameByIdx( n );
        long step = atol(name.c_str());
        steps.push_back(step);
        //cout << "n : " << n << " , " << step << endl;
    }
    std::sort(steps.begin(),steps.end());
    for( int n = 0 ; n < numberOfStep ; ++n ) {
        cout << "n : " << n << " , " << steps[n] << endl;
    }
    //
    cout << "##### SnapshotGroup::outputInfo() end #####" << endl;
}


//-------------------------- HDFSnapshotFile ---------------------------

HDFSnapshotFile::HDFSnapshotFile( string snapshotFilename ) {
#ifdef SNAPSHOT_DEBUG
    cout << "##### HDFSnapshotFile() start #####" << endl;
#endif
    Exception::dontPrint();
    this->snapshotFilename = snapshotFilename;
    //cout << "snapshotFilename : " << snapshotFilename << endl;
    rank = -1;
    numberOfRank = -1;
    numberOfAtom = 0;
    snapshotGroup = new SnapshotGroup();
#ifdef SNAPSHOT_DEBUG
    cout << "##### HDFSnapshotFile() end   #####" << endl;
#endif
}

int HDFSnapshotFile::exists() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### HDFSnapshotFile::exists() start #####" << endl;
#endif
    //
    int isHDF5 = H5Fis_hdf5( snapshotFilename.c_str() );
    if( isHDF5 <= 0 ) {
        return 0;
    } else {
        return 1;
    }
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### HDFSnapshotFile::exists() end #####" << endl;
#endif
}

void HDFSnapshotFile::write() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### HDFSnapshotFile::write() start #####" << endl;
#endif
    //
    int isHDF5 = exists();
    //cout << "isHDF5 : " << isHDF5 << endl;
    H5File* snapshotFile;
    if( isHDF5 <= 0) {
        try {
            snapshotFile = new H5File( snapshotFilename , H5F_ACC_TRUNC );
        } catch( Exception e ) {
            cout << "H5F_ACC_TRUNC" << endl;
            e.printError();
            string msg = "Can't create H5File on HDF snapshot file(H5F_ACC_TRUNC)";
            throw std::runtime_error(errorPos(msg));
        }
    } else {
        try {
            snapshotFile = new H5File( snapshotFilename , H5F_ACC_RDWR );
        } catch( Exception e ) {
            cout << "H5F_ACC_RDWR" << endl;
            e.printError();
            string msg = "Can't open H5File on HDF snapshot file(H5F_ACC_RDWR)";
            throw std::runtime_error(errorPos(msg));
        }
    }
    Group* rootGroup = new Group( snapshotFile->openGroup( "/" ) );
    //
    time_t fileTimeStamp = time(NULL);
    string timeStamp = getFileTimeStampString( fileTimeStamp );
    //cout << "fileTimeStamp : " << timeStamp << endl;
    if( isHDF5 <= 0) {
        //cout << "programVersion : " << *programVersion << endl;
        try {
            writeAttributeString( rootGroup , "ProgramVersion" , *programVersion );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't write attribute [ProgramVersion] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
        //cout << "restartFileVersion : " << SNAPSHOT_FILE_VERSION << endl;
        try {
            writeAttributeString( rootGroup , "SnapshotFileVersion" , SNAPSHOT_FILE_VERSION );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't write attribute [SnapshotFileVersion] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
        //
        try {
            writeAttribute( rootGroup , "NumberOfAtom" , &numberOfAtom , PredType::IEEE_F64LE );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't write attribute [NumberOfAtom] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
        try {
            writeAttribute( rootGroup , "Rank" , &rank , PredType::STD_I32LE );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't write attribute [Rank] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
        try {
            writeAttribute( rootGroup , "NumberOfRank" , &numberOfRank , PredType::STD_I32LE );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't write attribute [NumberOfRank] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
        //
        try {
            writeAttributeString( rootGroup , "FileTimeStamp" , timeStamp );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't write attribute [FileTimeStamp] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
        //
        long long numberOfTimeStep = 1;
        try {
            writeAttribute( rootGroup , "NumberOfTimeStep" , &numberOfTimeStep , PredType::STD_I64LE );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't write attribute [NumberOfTimeStep] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
    } else {
        try {
            updateAttributeString( rootGroup , "FileTimeStamp" , timeStamp );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't update attribute [FileTimeStamp] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
        long long numberOfTimeStep;
        try {
            readAttribute( rootGroup , "NumberOfTimeStep" , &numberOfTimeStep , PredType::STD_I64LE );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't write attribute [NumberOfTimeStep] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
        ++numberOfTimeStep;
        try {
            updateAttribute( rootGroup , "NumberOfTimeStep" , &numberOfTimeStep , PredType::STD_I64LE );
        } catch( Exception e ) {
            e.printError();
            string msg = "Can't update attribute [NumberOfTimeStep] on HDF snapshot file";
            throw std::runtime_error(errorPos(msg));
        }
    }
    //
    snapshotGroup->write( rootGroup );
    snapshotFile->close();
    //
    if( rootGroup != NULL ) delete rootGroup;
    if( snapshotFile != NULL ) delete snapshotFile;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### HDFSnapshotFile::write() end   #####" << endl;
#endif
}

void HDFSnapshotFile::read() {
#ifdef SNAPSHOT_DEBUG
    cout << "##### HDFSnapshotFile::read() start #####" << endl;
#endif
    //
    H5File* snapshotFile = new H5File( snapshotFilename , H5F_ACC_RDONLY );
    Group* rootGroup = new Group( snapshotFile->openGroup( "/" ) );
    //
    readAttribute( rootGroup , "NumberOfRank" , &numberOfRank , PredType::STD_I32LE );
    //
    snapshotFile->close();
    //
    if( rootGroup != NULL ) delete rootGroup;
    if( snapshotFile != NULL ) delete snapshotFile;
    //
#ifdef SNAPSHOT_DEBUG
    cout << "##### HDFSnapshotFile::read() end   #####" << endl;
#endif
}

void HDFSnapshotFile::read(size_t step) {
    cout << "##### HDFSnapshotFile::read() start #####" << endl;
    //
    H5File* snapshotFile = new H5File( snapshotFilename , H5F_ACC_RDONLY );
    Group* rootGroup = new Group( snapshotFile->openGroup( "/" ) );
    //
    snapshotGroup->read( rootGroup , step);
    //
    SnapshotDatatype* snapshotDatatype = getSnapshotDatatype();
    cout << "SnapshotDatatype ----------" << endl;
    cout << "   TimeStepInterval : " << snapshotDatatype->TimeStepInterval << endl;
    cout << "   CurrentTimeStep  : " << snapshotDatatype->CurrentTimeStep << endl;
    cout << "   CurrentReductionCounter  : " << snapshotDatatype->CurrentReductionCounter << endl;
    cout << "   CurrentPrintCounter  : " << snapshotDatatype->CurrentPrintCounter << endl;
    cout << "   CurrentCRDCounter  : " << snapshotDatatype->CurrentCRDCounter << endl;
    cout << "   CurrentRSTCounter  : " << snapshotDatatype->CurrentRSTCounter << endl;
    cout << "   AtomArrayLength  : " << snapshotDatatype->AtomArrayLength << endl;
    //
    IntegrationDatatype* integrationDatatype = getIntegrationDatatype();
    cout << "IntegrationDatatype ----------" << endl;
    cout << "   LJEnergyCorrection : " << integrationDatatype->LJEnergyCorrection << endl;
    cout << "   KineticEnergy : " << integrationDatatype->KineticEnergy << endl;
    cout << "   PotentialEnergy : " << integrationDatatype->PotentialEnergy << endl;
    cout << "   TotalEnergy : " << integrationDatatype->TotalEnergy << endl;
    cout << "   Virial : " << integrationDatatype->Virial << endl;
    cout << "   EtaPosition : " << integrationDatatype->EtaPosition << endl;
    cout << "   EtaVelocity : " << integrationDatatype->EtaVelocity << endl;
    cout << "   EtaForce : " << integrationDatatype->EtaForce << endl;
    cout << "   LogVPosition : " << integrationDatatype->LogVPosition << endl;
    cout << "   LogVVelocity : " << integrationDatatype->LogVVelocity << endl;
    cout << "   LogVForce : " << integrationDatatype->LogVForce << endl;
    cout << "   Volume : " << integrationDatatype->Volume << endl;
    //
    BoxDatatype* boxDatatype = getBoxDatatype();
    cout << "BoxDatatype ----------" << endl;
    cout << "   SizeX      : " << boxDatatype->SizeX << endl;
    cout << "   SizeY      : " << boxDatatype->SizeY << endl;
    cout << "   SizeZ      : " << boxDatatype->SizeZ << endl;
    cout << "   Type       : " << boxDatatype->Type << endl;
    cout << "   BetaAngle1 : " << boxDatatype->BetaAngle1 << endl;
    cout << "   BetaAngle2 : " << boxDatatype->BetaAngle2 << endl;
    cout << "   BetaAngle3 : " << boxDatatype->BetaAngle3 << endl;
    //
    AtomDatatype* atomDatatype = getAtomDatatype();
    cout << "AtomDatatype ----------" << endl;
    long long numberOfAtom = getAtomDataset()->getNumberOfAtom();
    cout << "   NumberOfAtom : " << numberOfAtom << endl;
    for( int n = 0 ; n < 10 ; ++n ) {
        cout << "   " << n << "-th pos : " << atomDatatype[n].PositionX << " , " << atomDatatype[n].PositionY << " , " << atomDatatype[n].PositionZ << endl;
    }
    //
    if( rootGroup != NULL ) delete rootGroup;
    if( snapshotFile != NULL ) delete snapshotFile;
    //
    cout << "##### HDFSnapshotFile::read() end   #####" << endl;
}

void HDFSnapshotFile::outputInfo() {
    cout << "##### HDFSnapshotFile::outputInfo() start #####" << endl;
    //
    H5File* snapshotFile = new H5File( snapshotFilename , H5F_ACC_RDONLY );
    Group* rootGroup = new Group( snapshotFile->openGroup( "/" ) );
    //
    string programVersion;
    readAttributeString( rootGroup , "ProgramVersion" , &programVersion );
    cout << "ProgramVersion      : " << programVersion << endl;
    string snapshotFileVersion;
    readAttributeString( rootGroup , "SnapshotFileVersion" , &snapshotFileVersion );
    cout << "SnapshotFileVersion : " << snapshotFileVersion << endl;
    string fileTimeStamp;
    readAttributeString( rootGroup , "FileTimeStamp" , &fileTimeStamp );
    cout << "FileTimeStamp       : " << fileTimeStamp << endl;
    long long numberOfTimeStep;
    readAttribute( rootGroup , "NumberOfTimeStep" , &numberOfTimeStep , PredType::STD_I64LE );
    cout << "NumberOfTimeStep    : " << numberOfTimeStep << endl;
    hsize_t size = snapshotFile->getFileSize();
    cout << "Snapshot file size  : " << size << endl;
    //
    //snapshotGroup->outputInfo( rootGroup );
    //
    if( rootGroup != NULL ) delete rootGroup;
    if( snapshotFile != NULL ) delete snapshotFile;
    //
    cout << "##### HDFSnapshotFile::outputInfo() end   #####" << endl;
}

string HDFSnapshotFile::getFileTimeStampString( time_t fileTimeStamp ) {
    char timeStamp[20];
    struct tm* ts = localtime( &fileTimeStamp );
    sprintf( timeStamp , "%04d/%02d/%02d %02d:%02d:%02d",
          ts->tm_year+1900, ts->tm_mon+1, ts->tm_mday, ts->tm_hour, ts->tm_min, ts->tm_sec);
    return timeStamp;
}

}

