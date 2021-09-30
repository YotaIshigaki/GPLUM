#ifndef HDFSNAPSHOTFILE_H
#define HDFSNAPSHOTFILE_H
#include "H5Cpp.h"
#include <vector>
#include <iostream>

using namespace std;

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

namespace HDFSnapshot {

#define SNAPSHOT_FILE_VERSION "1.00.00"

typedef struct BoxDatatype {
    double SizeX;
    double SizeY;
    double SizeZ;
    long long Type;
    double BetaAngle1;
    double BetaAngle2;
    double BetaAngle3;
  /*
    BoxDatatype() : SizeX(0),SizeY(0),SizeZ(0),Type(1),
                    BetaAngle1(0),BetaAngle2(0),BetaAngle3(0)
    {
    }
  */
} BoxDatatype;

class BoxDataset {
public:
    BoxDataset();
    BoxDatatype* getBoxDatatype() { return &boxDatatype; }
    void read( Group *group );
    void write( Group *group );
private:
    CompType* createCompType();
    BoxDatatype boxDatatype;
};

typedef struct AtomDatatype {
    long long AtomId;
    double PositionX;
    double PositionY;
    double PositionZ;
    double Charge;
    double VelocityX;
    double VelocityY;
    double VelocityZ;
    double ForceX;
    double ForceY;
    double ForceZ;
    double Mass;
    double InverseMass;
    long long AtomType;
  /*
    AtomDatatype() : AtomId(0),PositionX(0),PositionY(0),PositionZ(0),Charge(0),
                           VelocityX(0),VelocityY(0),VelocityZ(0),
      ForceX(0),ForceY(0),ForceZ(0),
      Mass(0),InverseMass(0),AtomType(0)
    {
    }
  */
} AtomDatatype;

class AtomDataset {
public:
    AtomDataset();
    ~AtomDataset() {
        if( atomDatatype != NULL ) delete[] atomDatatype;
    }
    AtomDatatype* getAtomDatatype() { return atomDatatype; }
    void read( Group *group );
    void write( Group *group );
    void setNumberOfAtom( long long no ) {
        numberOfAtom = no;
        if( atomDatatype != NULL ) delete[] atomDatatype;
        atomDatatype = new AtomDatatype[numberOfAtom];
    }
    long long getNumberOfAtom() { return numberOfAtom; }
private:
    CompType* createCompType();
    AtomDatatype *atomDatatype;
    long long numberOfAtom;
};

typedef struct TyperangeDatatype {
    long long Begin;
    long long End;
    long long LJBegin;
    long long LJEnd;
    long long LJCoulombBegin;
    long long LJCoulombEnd;
    long long CoulombBegin;
    long long CoulombEnd;
  /*
    TyperangeDatatype() : Begin(0),End(0),
                          LJBegin(0),LJEnd(0),
                          LJCoulombBegin(0),LJCoulombEnd(0),
                          CoulombBegin(0),CoulombEnd(0)
    {
    }
  */
} TyperangeDatatype;

class TyperangeDataset {
public:
    TyperangeDataset();
    ~TyperangeDataset() {
        if( typerangeDatatype != NULL ) delete[] typerangeDatatype;
    }
    TyperangeDatatype* getTyperangeDatatype() { return typerangeDatatype; }
    void read( Group *group );
    void write( Group *group );
    void setNumberOfTyperange( long long no ) {
        numberOfTyperange = no;
        if( typerangeDatatype != NULL ) delete[] typerangeDatatype;
        typerangeDatatype = new TyperangeDatatype[numberOfTyperange];
    }
    long long getNumberOfTyperange() { return numberOfTyperange; }
private:
    CompType* createCompType();
    TyperangeDatatype *typerangeDatatype;
    long long numberOfTyperange;
};

typedef struct BondDatatype {
  long long AtomId0;
  long long AtomId1;
  long long TypeOfBond;
  long long Shake;
} BondDatatype;

class BondDataset {
 public:
  BondDataset();
  ~BondDataset(){
    if( bondDatatype != NULL ) delete[] bondDatatype;
  }
  BondDatatype *getBondDatatype() { return bondDatatype; }
  void read( Group *structureGroup );
  void write( Group *structureGroup );
  int getNumberOfBond() { return numberOfBond; }
  void setNumberOfBond( int no ) {
    numberOfBond = no;
    if( bondDatatype != NULL ) delete[] bondDatatype;
    if( numberOfBond > 0 ) bondDatatype = new BondDatatype[numberOfBond];
  }
private:
  CompType* createMemType();
  BondDatatype *bondDatatype;
  int numberOfBond;
};

typedef struct AngleDatatype {
  long long AtomId0;
  long long AtomId1;
  long long AtomId2;
  long long TypeOfAngle;
} AngleDatatype;

class AngleDataset {
public:
    AngleDataset();
    ~AngleDataset() {
        if( angleDatatype != NULL ) delete[] angleDatatype;
    }
    AngleDatatype* getAngleDatatype() { return angleDatatype; }
    void read( Group *structureGroup );
    void write( Group *structureGroup );
    int getNumberOfAngle() { return numberOfAngle; }
    void setNumberOfAngle( int no ) {
        numberOfAngle = no;
        if( angleDatatype != NULL ) delete[] angleDatatype;
        if( numberOfAngle > 0 ) angleDatatype = new AngleDatatype[numberOfAngle];
   }
private:
    CompType* createMemType();
    AngleDatatype *angleDatatype;
    int numberOfAngle;
};

typedef struct TorsionDatatype {
  long long AtomId0;
  long long AtomId1;
  long long AtomId2;
  long long AtomId3;
  long long TypeOfTorsion;
  long long Calc14;
} TorsionDatatype;

class TorsionDataset {
public:
    TorsionDataset();
    ~TorsionDataset() {
        if( torsionDatatype != NULL ) delete[] torsionDatatype;
    }
    TorsionDatatype* getTorsionDatatype() { return torsionDatatype; }
    void read( Group *dihedralGroup );
    void write( Group *dihedralGroup );
    int getNumberOfTorsion() { return numberOfTorsion; }
    void setNumberOfTorsion( int no ) {
        numberOfTorsion = no;
        if( torsionDatatype != NULL ) delete[] torsionDatatype;
        if( numberOfTorsion > 0 ) torsionDatatype = new TorsionDatatype[numberOfTorsion];
   }
private:
    CompType* createMemType();
    TorsionDatatype *torsionDatatype;
    int numberOfTorsion;
};

typedef struct ImproperDatatype {
  long long AtomId0;
  long long AtomId1;
  long long AtomId2;
  long long AtomId3;
  long long TypeOfImproper;
} ImproperDatatype;

class ImproperDataset {
public:
    ImproperDataset();
    ~ImproperDataset() {
        if( improperDatatype != NULL ) delete[] improperDatatype;
    }
    ImproperDatatype* getImproperDatatype() { return improperDatatype; }
    void read( Group *dihedralGroup );
    void write( Group *dihedralGroup );
    int getNumberOfImproper() { return numberOfImproper; }
    void setNumberOfImproper( int no ) {
        numberOfImproper = no;
        if( improperDatatype != NULL ) delete[] improperDatatype;
        if( numberOfImproper > 0 ) improperDatatype = new ImproperDatatype[numberOfImproper];
   }
private:
    CompType* createMemType();
    ImproperDatatype *improperDatatype;
    int numberOfImproper;
};

typedef struct BondNumberDatatype{
  long long NumberOfBond;
  long long NumberOfAngle;
  long long NumberOfTrosion;
  long long NumberOfImproper;
} BondNumberDatatype;

class BondNumberDataset {
public:
    BondNumberDataset();
    BondNumberDatatype* getBondNumberDatatype() { return &bondnumberDatatype; }
    void read( Group *group );
    void write( Group *group );
private:
    CompType* createCompType();
    BondNumberDatatype bondnumberDatatype;
};

class BondlistGroup {
 public:
  BondlistGroup();
  ~BondlistGroup() {
    if(bondDataset!=NULL)delete bondDataset;
    if(angleDataset!=NULL)delete angleDataset;
    if(torsionDataset!=NULL)delete torsionDataset;
    if(improperDataset!=NULL)delete improperDataset;
  }
  void setIndex(int index){
    this->index = index;
  }
  BondNumberDataset *getBondNumberDataset() {return bondnumberDataset;}
  BondDataset* getBondDataset() {return bondDataset;}
  AngleDataset* getAngleDataset() {return angleDataset;}
  TorsionDataset* getTorsionDataset() {return torsionDataset;}
  ImproperDataset* getImproperDataset() {return improperDataset;}
  void read (Group *bondlistarrayGroup);
  void write (Group *bondlistarrayGroup);

 private:
  BondNumberDataset *bondnumberDataset;
  BondDataset* bondDataset;
  AngleDataset* angleDataset;
  TorsionDataset* torsionDataset;
  ImproperDataset* improperDataset;
  int index;
};

class BondlistarrayGroup {
 public:
  BondlistarrayGroup();
  ~BondlistarrayGroup(){
    if(bondlistGroup!=NULL)delete[] bondlistGroup;
  }
  void setNumberOfBondlist(int numberOfBondlist) {
    this->numberOfBondlist = numberOfBondlist;
    if(bondlistGroup!=NULL)delete[] bondlistGroup;
    bondlistGroup = new BondlistGroup[numberOfBondlist];
    for(int n=0;n<numberOfBondlist; n++){
      bondlistGroup[n].setIndex(n);
    }
  }
  BondlistGroup* getBondlistGroup(int bondlistIndex){
    if((bondlistIndex<0) || (bondlistIndex>=numberOfBondlist)){
      return NULL;
    }
    return &bondlistGroup[bondlistIndex];
  }
  int getNumberOfBondlist(){return numberOfBondlist;}
  void read( Group *group);
  void write( Group *group);

 private:
  BondlistGroup* bondlistGroup;
  int numberOfBondlist;
};

typedef struct WaterDatatype {
  long long OWId0;
  long long HWId1;
  long long HWId2;
} WaterDatatype;

class WaterDataset {
public:
    WaterDataset();
    ~WaterDataset() {
        if( waterDatatype != NULL ) delete[] waterDatatype;
    }
    WaterDatatype* getWaterDatatype() { return waterDatatype; }
    void read( Group *structureGroup );
    void write( Group *structureGroup );
    int getNumberOfWater() { return numberOfWater; }
    void setNumberOfWater( int no ) {
        numberOfWater = no;
        if( waterDatatype != NULL ) delete[] waterDatatype;
        if( numberOfWater > 0 ) waterDatatype = new WaterDatatype[numberOfWater];
    }
 private:
    CompType* createCompType();
    WaterDatatype *waterDatatype;
    int numberOfWater;
};

#ifndef MAX_HBOND
#define MAX_HBOND 4    /// ParticleInfo.h
#endif
typedef struct ShakeDatatype {
  long long HAId0;
  long long Nh1;
  long long H1Id1;
  long long H1Id2;
  long long H1Id3;
  long long H1Id4;
  long long BondType1;
  long long BondType2;
  long long BondType3;
  long long BondType4;
} ShakeDatatype;

class ShakeDataset {
public:
    ShakeDataset();
    ~ShakeDataset() {
        if( shakeDatatype != NULL ) delete[] shakeDatatype;
    }
    ShakeDatatype* getShakeDatatype() { return shakeDatatype; }
    void read( Group *structureGroup );
    void write( Group *structureGroup );
    int getNumberOfShake() { return numberOfShake; }
    void setNumberOfShake( int no ) {
        numberOfShake = no;
        if( shakeDatatype != NULL ) delete[] shakeDatatype;
        if( numberOfShake > 0 ) shakeDatatype = new ShakeDatatype[numberOfShake];
    }
 private:
    CompType* createCompType();
    ShakeDatatype *shakeDatatype;
    int numberOfShake;
};

typedef struct IntegrationDatatype {
  double LJEnergyCorrection;
  double KineticEnergy;
  double PotentialEnergy;
  double TotalEnergy;
  double Virial;
  double EtaPosition;
  double EtaVelocity;
  double EtaForce;
  double LogVPosition;
  double LogVVelocity;
  double LogVForce;
  double Volume;
  /*
  IntegrationDatatype() : LJEnergyCorrection(0),
                      KineticEnergy(0),PotentialEnergy(0),TotalEnergy(0),Virial(0),
                      EtaPosition(0),EtaVelocity(0),EtaForce(0),
                      LogVPosition(0),LogVVelocity(0),LogVForce(0),Volume(0)
  {}
  */
} IntegrationDatatype;

class IntegrationDataset {
public:
    IntegrationDataset();
    ~IntegrationDataset(){};
    IntegrationDatatype* getIntegrationDatatype() { return &integrationDatatype; }
    void read( Group *group );
    void write( Group *group );
private:
    CompType* createCompType();
    IntegrationDatatype integrationDatatype;
};

typedef struct SnapshotDatatype {
    double TimeStepInterval;
    long long CurrentTimeStep;
    long long CurrentReductionCounter;
    long long CurrentPrintCounter;
    long long CurrentCRDCounter;
    long long CurrentRSTCounter;
    long long AtomArrayLength;
  /*
    SnapshotDatatype() : TimeStepInterval(0),CurrentTimeStep(0),
                         CurrentReductionCounter(0),CurrentPrintCounter(0),
                         CurrentCRDCounter(0),CurrentRSTCounter(0),
                         AtomArrayLength(0)
    {
    }
  */
} SnapshotDatatype;

class SnapshotDataset {
public:
    SnapshotDataset();
    ~SnapshotDataset() {
    }
    void write( Group* group );
    void read( Group* group );
    SnapshotDatatype* getSnapshotDatatype() { return &snapshotDatatype; }
private:
    SnapshotDatatype snapshotDatatype;
    CompType* createCompType();
};

//------------------------- SnapshotGroup -------------------------

class SnapshotGroup {
public:
    SnapshotGroup();
    ~SnapshotGroup() {
        cout << "##### ~SnapshotGroup start #####" << endl;
        if( snapshotDataset != NULL ) delete snapshotDataset;
        if( integrationDataset != NULL ) delete integrationDataset;
        if( boxDataset != NULL ) delete boxDataset;
        if( atomDataset != NULL ) delete atomDataset;
        if( typerangeDataset != NULL ) delete typerangeDataset;
	if( bondlistarrayGroup != NULL) delete bondlistarrayGroup;
	if( waterDataset !=NULL ) delete waterDataset;
	if( shakeDataset !=NULL ) delete shakeDataset;
        cout << "##### ~SnapshotGroup end #####" << endl;
    }
    void write( Group* rootGroup );
    void read( Group* rootGroup , size_t step, int verbose=0 );

    void setParam( double timeStepInterval , long long currentTimeStep, 
		   long long currentReductionCounter, long long currentPrintCounter,
		   long long currentCRDCounter, long long currentRSTCounter,
		   long long atomArrayLength );
    SnapshotDataset* getSnapshotDataset() { return snapshotDataset; }
    IntegrationDataset* getIntegrationDataset() { return integrationDataset; }
    BoxDataset* getBoxDataset() { return boxDataset; }
    AtomDataset* getAtomDataset() { return atomDataset; }
    TyperangeDataset* getTyperangeDataset() { return typerangeDataset; }
    BondlistarrayGroup* getBondlistarrayGroup() {return bondlistarrayGroup; }
    WaterDataset* getWaterDataset() {return waterDataset; }
    ShakeDataset* getShakeDataset() {return shakeDataset; }
    void outputInfo( Group* rootGroup );
    int getsteps( Group* rootGroup, std::vector<long>& steps);

private:
    SnapshotDataset* snapshotDataset;
    IntegrationDataset* integrationDataset;
    BoxDataset* boxDataset;
    AtomDataset* atomDataset;
    TyperangeDataset* typerangeDataset;
    BondlistarrayGroup* bondlistarrayGroup;
    WaterDataset* waterDataset;
    ShakeDataset* shakeDataset;
};

//------------------------- HDFSnapshotFile -----------------------

class HDFSnapshotFile {
public:
    enum SnapshotType {
      Single,
      Parallel,
      SingleDump,
      ParallelDump
    };

    HDFSnapshotFile( string snapshotFilename );
    int exists();
    void write();
    void read();
    void read( size_t step, int verbose=0);

    SnapshotGroup* getSnapshotGroup() { return snapshotGroup; }
    SnapshotDatatype* getSnapshotDatatype() { return snapshotGroup->getSnapshotDataset()->getSnapshotDatatype(); }
    IntegrationDatatype* getIntegrationDatatype() { return snapshotGroup->getIntegrationDataset()->getIntegrationDatatype(); }
    BoxDatatype* getBoxDatatype() { return snapshotGroup->getBoxDataset()->getBoxDatatype(); }
    AtomDataset* getAtomDataset() { return snapshotGroup->getAtomDataset(); }
    AtomDatatype* getAtomDatatype() { return getAtomDataset()->getAtomDatatype(); }
    TyperangeDataset* getTyperangeDataset() { return snapshotGroup->getTyperangeDataset(); }
    TyperangeDatatype* getTyperangeDatatype() { return getTyperangeDataset()->getTyperangeDatatype(); }
    BondlistarrayGroup* getBondlistarrayGroup() {return snapshotGroup->getBondlistarrayGroup();}
    WaterDataset *getWaterDataset() {return snapshotGroup->getWaterDataset();}
    WaterDatatype *getWaterDatatype() {return getWaterDataset()->getWaterDatatype();}
    ShakeDataset *getShakeDataset() {return snapshotGroup->getShakeDataset();}
    ShakeDatatype *getShakeDatatype() {return getShakeDataset()->getShakeDatatype();}


    void setProgramVersion( string* programVersion ) {
        this->programVersion = programVersion;
    }
    void setNumberOfAtom( long long numberOfAtom ) { this->numberOfAtom = numberOfAtom;}
    void setRank( int rank , int numberOfRank ) {
        //cout << "##### setRank start #####" << endl;
        this->rank = rank;
        this->numberOfRank = numberOfRank;
        //cout << "##### setRank end #####" << endl;
    }
    void setSnapshotType(SnapshotType snapshotytype){ this->snapshotType = snapshotytype; }

    string getProgramVersion() {return *programVersion; }
    long long getNumberOfAtom() {return numberOfAtom; }
    int getRank() { return rank; }
    int getNumberOfRank() { return numberOfRank; }
    SnapshotType getSnapshotType() { return snapshotType; }

    void outputInfo();


private:
    string snapshotFilename;
    string *programVersion;

    int count;
    int rank;
    int numberOfRank;
    long long numberOfAtom;
    SnapshotType snapshotType;

    SnapshotGroup* snapshotGroup;

    string getFileTimeStampString( time_t fileTimeStamp );
    time_t getFileTimeStampTime( string fileTimeStampString );
};

}
#endif // HDFSNAPSHOTFILE_H
