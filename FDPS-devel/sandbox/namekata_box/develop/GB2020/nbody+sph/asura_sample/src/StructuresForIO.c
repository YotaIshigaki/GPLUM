#include "config.h"
#include "DataStructures.h"
#include "StructuresForIO.h"

// This file was generated at 2019-10-09 16:39:07 +0900

struct StructPbodyIOCompact CopyPbodyToTemporalStructureCompact(const int index){

	struct StructPbodyIOCompact TempPbodyIOCompact;

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    TempPbodyIOCompact.GlobalID = Pbody[index]->GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
    TempPbodyIOCompact.Active = Pbody[index]->Active;

// omit int InteractionList; // <TMP> <LEAN>

    TempPbodyIOCompact.Pos[0] = (float)Pbody[index]->Pos[0];
    TempPbodyIOCompact.Pos[1] = (float)Pbody[index]->Pos[1];
    TempPbodyIOCompact.Pos[2] = (float)Pbody[index]->Pos[2];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    TempPbodyIOCompact.Vel[0] = (float)Pbody[index]->Vel[0];
    TempPbodyIOCompact.Vel[1] = (float)Pbody[index]->Vel[1];
    TempPbodyIOCompact.Vel[2] = (float)Pbody[index]->Vel[2];
    TempPbodyIOCompact.Velh[0] = (float)Pbody[index]->Velh[0];
    TempPbodyIOCompact.Velh[1] = (float)Pbody[index]->Velh[1];
    TempPbodyIOCompact.Velh[2] = (float)Pbody[index]->Velh[2];
    TempPbodyIOCompact.Acc[0] = (float)Pbody[index]->Acc[0];
    TempPbodyIOCompact.Acc[1] = (float)Pbody[index]->Acc[1];
    TempPbodyIOCompact.Acc[2] = (float)Pbody[index]->Acc[2];
    TempPbodyIOCompact.AccOld[0] = (float)Pbody[index]->AccOld[0];
    TempPbodyIOCompact.AccOld[1] = (float)Pbody[index]->AccOld[1];
    TempPbodyIOCompact.AccOld[2] = (float)Pbody[index]->AccOld[2];
    TempPbodyIOCompact.Pot = (float)Pbody[index]->Pot;
    TempPbodyIOCompact.Mass = (float)Pbody[index]->Mass;
    TempPbodyIOCompact.Eps = (float)Pbody[index]->Eps;

    TempPbodyIOCompact.Type = Pbody[index]->Type;
// omit short     k;            // <TMP> <LEAN>
    TempPbodyIOCompact.dt = (float)Pbody[index]->dt;
// omit double    EraLocal;  // <TMP> <LEAN>

	return TempPbodyIOCompact;

}

StructPbody CopyTemporalStructureCompactToPbodyCompact(struct StructPbodyIOCompact PbodyIOCompact){

	StructPbody TemporalPbody;
	memset(&TemporalPbody,0,sizeof(StructPbody));

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    TemporalPbody.GlobalID = PbodyIOCompact.GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
    TemporalPbody.Active = PbodyIOCompact.Active;

// omit int InteractionList; // <TMP> <LEAN>

    TemporalPbody.Pos[0] = PbodyIOCompact.Pos[0];
    TemporalPbody.Pos[1] = PbodyIOCompact.Pos[1];
    TemporalPbody.Pos[2] = PbodyIOCompact.Pos[2];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    TemporalPbody.Vel[0] = PbodyIOCompact.Vel[0];
    TemporalPbody.Vel[1] = PbodyIOCompact.Vel[1];
    TemporalPbody.Vel[2] = PbodyIOCompact.Vel[2];
    TemporalPbody.Velh[0] = PbodyIOCompact.Velh[0];
    TemporalPbody.Velh[1] = PbodyIOCompact.Velh[1];
    TemporalPbody.Velh[2] = PbodyIOCompact.Velh[2];
    TemporalPbody.Acc[0] = PbodyIOCompact.Acc[0];
    TemporalPbody.Acc[1] = PbodyIOCompact.Acc[1];
    TemporalPbody.Acc[2] = PbodyIOCompact.Acc[2];
    TemporalPbody.AccOld[0] = PbodyIOCompact.AccOld[0];
    TemporalPbody.AccOld[1] = PbodyIOCompact.AccOld[1];
    TemporalPbody.AccOld[2] = PbodyIOCompact.AccOld[2];
    TemporalPbody.Pot = PbodyIOCompact.Pot;
    TemporalPbody.Mass = PbodyIOCompact.Mass;
    TemporalPbody.Eps = PbodyIOCompact.Eps;

    TemporalPbody.Type = PbodyIOCompact.Type;
// omit short     k;            // <TMP> <LEAN>
    TemporalPbody.dt = PbodyIOCompact.dt;
// omit double    EraLocal;  // <TMP> <LEAN>

	return TemporalPbody;

}

struct StructPbodyIOCompactDouble CopyPbodyToTemporalStructureCompactDouble(const int index){

	struct StructPbodyIOCompactDouble TempPbodyIOCompactDouble;

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    TempPbodyIOCompactDouble.GlobalID = Pbody[index]->GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
    TempPbodyIOCompactDouble.Active = Pbody[index]->Active;

// omit int InteractionList; // <TMP> <LEAN>

    TempPbodyIOCompactDouble.Pos[0] = (double)Pbody[index]->Pos[0];
    TempPbodyIOCompactDouble.Pos[1] = (double)Pbody[index]->Pos[1];
    TempPbodyIOCompactDouble.Pos[2] = (double)Pbody[index]->Pos[2];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    TempPbodyIOCompactDouble.Vel[0] = (double)Pbody[index]->Vel[0];
    TempPbodyIOCompactDouble.Vel[1] = (double)Pbody[index]->Vel[1];
    TempPbodyIOCompactDouble.Vel[2] = (double)Pbody[index]->Vel[2];
    TempPbodyIOCompactDouble.Velh[0] = (double)Pbody[index]->Velh[0];
    TempPbodyIOCompactDouble.Velh[1] = (double)Pbody[index]->Velh[1];
    TempPbodyIOCompactDouble.Velh[2] = (double)Pbody[index]->Velh[2];
    TempPbodyIOCompactDouble.Acc[0] = (double)Pbody[index]->Acc[0];
    TempPbodyIOCompactDouble.Acc[1] = (double)Pbody[index]->Acc[1];
    TempPbodyIOCompactDouble.Acc[2] = (double)Pbody[index]->Acc[2];
    TempPbodyIOCompactDouble.AccOld[0] = (double)Pbody[index]->AccOld[0];
    TempPbodyIOCompactDouble.AccOld[1] = (double)Pbody[index]->AccOld[1];
    TempPbodyIOCompactDouble.AccOld[2] = (double)Pbody[index]->AccOld[2];
    TempPbodyIOCompactDouble.Pot = (double)Pbody[index]->Pot;
    TempPbodyIOCompactDouble.Mass = (double)Pbody[index]->Mass;
    TempPbodyIOCompactDouble.Eps = (double)Pbody[index]->Eps;

    TempPbodyIOCompactDouble.Type = Pbody[index]->Type;
// omit short     k;            // <TMP> <LEAN>
    TempPbodyIOCompactDouble.dt = (double)Pbody[index]->dt;
// omit double    EraLocal;  // <TMP> <LEAN>

	return TempPbodyIOCompactDouble;

}

StructPbody CopyTemporalStructureCompactDoubleToPbodyCompactDouble(struct StructPbodyIOCompactDouble PbodyIOCompactDouble){

	StructPbody TemporalPbody;
	memset(&TemporalPbody,0,sizeof(StructPbody));

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    TemporalPbody.GlobalID = PbodyIOCompactDouble.GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
    TemporalPbody.Active = PbodyIOCompactDouble.Active;

// omit int InteractionList; // <TMP> <LEAN>

    TemporalPbody.Pos[0] = PbodyIOCompactDouble.Pos[0];
    TemporalPbody.Pos[1] = PbodyIOCompactDouble.Pos[1];
    TemporalPbody.Pos[2] = PbodyIOCompactDouble.Pos[2];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    TemporalPbody.Vel[0] = PbodyIOCompactDouble.Vel[0];
    TemporalPbody.Vel[1] = PbodyIOCompactDouble.Vel[1];
    TemporalPbody.Vel[2] = PbodyIOCompactDouble.Vel[2];
    TemporalPbody.Velh[0] = PbodyIOCompactDouble.Velh[0];
    TemporalPbody.Velh[1] = PbodyIOCompactDouble.Velh[1];
    TemporalPbody.Velh[2] = PbodyIOCompactDouble.Velh[2];
    TemporalPbody.Acc[0] = PbodyIOCompactDouble.Acc[0];
    TemporalPbody.Acc[1] = PbodyIOCompactDouble.Acc[1];
    TemporalPbody.Acc[2] = PbodyIOCompactDouble.Acc[2];
    TemporalPbody.AccOld[0] = PbodyIOCompactDouble.AccOld[0];
    TemporalPbody.AccOld[1] = PbodyIOCompactDouble.AccOld[1];
    TemporalPbody.AccOld[2] = PbodyIOCompactDouble.AccOld[2];
    TemporalPbody.Pot = PbodyIOCompactDouble.Pot;
    TemporalPbody.Mass = PbodyIOCompactDouble.Mass;
    TemporalPbody.Eps = PbodyIOCompactDouble.Eps;

    TemporalPbody.Type = PbodyIOCompactDouble.Type;
// omit short     k;            // <TMP> <LEAN>
    TemporalPbody.dt = PbodyIOCompactDouble.dt;
// omit double    EraLocal;  // <TMP> <LEAN>

	return TemporalPbody;

}

struct StructPbodyIOCompactMix CopyPbodyToTemporalStructureCompactMix(const int index){

	struct StructPbodyIOCompactMix TempPbodyIOCompactMix;

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    TempPbodyIOCompactMix.GlobalID = Pbody[index]->GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
    TempPbodyIOCompactMix.Active = Pbody[index]->Active;

// omit int InteractionList; // <TMP> <LEAN>

    TempPbodyIOCompactMix.Pos[0] = (float)Pbody[index]->Pos[0];
    TempPbodyIOCompactMix.Pos[1] = (float)Pbody[index]->Pos[1];
    TempPbodyIOCompactMix.Pos[2] = (float)Pbody[index]->Pos[2];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    TempPbodyIOCompactMix.Vel[0] = (float)Pbody[index]->Vel[0];
    TempPbodyIOCompactMix.Vel[1] = (float)Pbody[index]->Vel[1];
    TempPbodyIOCompactMix.Vel[2] = (float)Pbody[index]->Vel[2];
    TempPbodyIOCompactMix.Velh[0] = (float)Pbody[index]->Velh[0];
    TempPbodyIOCompactMix.Velh[1] = (float)Pbody[index]->Velh[1];
    TempPbodyIOCompactMix.Velh[2] = (float)Pbody[index]->Velh[2];
    TempPbodyIOCompactMix.Acc[0] = (double)Pbody[index]->Acc[0];
    TempPbodyIOCompactMix.Acc[1] = (double)Pbody[index]->Acc[1];
    TempPbodyIOCompactMix.Acc[2] = (double)Pbody[index]->Acc[2];
    TempPbodyIOCompactMix.AccOld[0] = (double)Pbody[index]->AccOld[0];
    TempPbodyIOCompactMix.AccOld[1] = (double)Pbody[index]->AccOld[1];
    TempPbodyIOCompactMix.AccOld[2] = (double)Pbody[index]->AccOld[2];
    TempPbodyIOCompactMix.Pot = (double)Pbody[index]->Pot;
    TempPbodyIOCompactMix.Mass = (float)Pbody[index]->Mass;
    TempPbodyIOCompactMix.Eps = (float)Pbody[index]->Eps;

    TempPbodyIOCompactMix.Type = Pbody[index]->Type;
// omit short     k;            // <TMP> <LEAN>
    TempPbodyIOCompactMix.dt = (double)Pbody[index]->dt;
// omit double    EraLocal;  // <TMP> <LEAN>

	return TempPbodyIOCompactMix;

}

StructPbody CopyTemporalStructureCompactMixToPbodyCompactMix(struct StructPbodyIOCompactMix PbodyIOCompactMix){

	StructPbody TemporalPbody;
	memset(&TemporalPbody,0,sizeof(StructPbody));

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    TemporalPbody.GlobalID = PbodyIOCompactMix.GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
    TemporalPbody.Active = PbodyIOCompactMix.Active;

// omit int InteractionList; // <TMP> <LEAN>

    TemporalPbody.Pos[0] = PbodyIOCompactMix.Pos[0];
    TemporalPbody.Pos[1] = PbodyIOCompactMix.Pos[1];
    TemporalPbody.Pos[2] = PbodyIOCompactMix.Pos[2];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    TemporalPbody.Vel[0] = PbodyIOCompactMix.Vel[0];
    TemporalPbody.Vel[1] = PbodyIOCompactMix.Vel[1];
    TemporalPbody.Vel[2] = PbodyIOCompactMix.Vel[2];
    TemporalPbody.Velh[0] = PbodyIOCompactMix.Velh[0];
    TemporalPbody.Velh[1] = PbodyIOCompactMix.Velh[1];
    TemporalPbody.Velh[2] = PbodyIOCompactMix.Velh[2];
    TemporalPbody.Acc[0] = PbodyIOCompactMix.Acc[0];
    TemporalPbody.Acc[1] = PbodyIOCompactMix.Acc[1];
    TemporalPbody.Acc[2] = PbodyIOCompactMix.Acc[2];
    TemporalPbody.AccOld[0] = PbodyIOCompactMix.AccOld[0];
    TemporalPbody.AccOld[1] = PbodyIOCompactMix.AccOld[1];
    TemporalPbody.AccOld[2] = PbodyIOCompactMix.AccOld[2];
    TemporalPbody.Pot = PbodyIOCompactMix.Pot;
    TemporalPbody.Mass = PbodyIOCompactMix.Mass;
    TemporalPbody.Eps = PbodyIOCompactMix.Eps;

    TemporalPbody.Type = PbodyIOCompactMix.Type;
// omit short     k;            // <TMP> <LEAN>
    TemporalPbody.dt = PbodyIOCompactMix.dt;
// omit double    EraLocal;  // <TMP> <LEAN>

	return TemporalPbody;

}

struct StructPbodyIOLean CopyPbodyToTemporalStructureLean(const int index){

	struct StructPbodyIOLean TempPbodyIOLean;

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    TempPbodyIOLean.GlobalID = Pbody[index]->GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
// omit bool    Active; // <LEAN> This is used for the individual time step.

// omit int InteractionList; // <TMP> <LEAN>

    TempPbodyIOLean.Pos[0] = (float)Pbody[index]->Pos[0];
    TempPbodyIOLean.Pos[1] = (float)Pbody[index]->Pos[1];
    TempPbodyIOLean.Pos[2] = (float)Pbody[index]->Pos[2];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    TempPbodyIOLean.Vel[0] = (float)Pbody[index]->Vel[0];
    TempPbodyIOLean.Vel[1] = (float)Pbody[index]->Vel[1];
    TempPbodyIOLean.Vel[2] = (float)Pbody[index]->Vel[2];
// omit double    Velh[3];      // <MIX> <LEAN> The \"half-step slided\" velocity of the particle. For leapfrog integration.
    TempPbodyIOLean.Acc[0] = (float)Pbody[index]->Acc[0];
    TempPbodyIOLean.Acc[1] = (float)Pbody[index]->Acc[1];
    TempPbodyIOLean.Acc[2] = (float)Pbody[index]->Acc[2];
// omit double    AccOld[3];    // <LEAN> The acceleration of the particle.
    TempPbodyIOLean.Pot = (float)Pbody[index]->Pot;
    TempPbodyIOLean.Mass = (float)Pbody[index]->Mass;
    TempPbodyIOLean.Eps = (float)Pbody[index]->Eps;

    TempPbodyIOLean.Type = Pbody[index]->Type;
// omit short     k;            // <TMP> <LEAN>
// omit double    dt;           // <LEAN>
// omit double    EraLocal;  // <TMP> <LEAN>

	return TempPbodyIOLean;

}

StructPbody CopyTemporalStructureLeanToPbodyLean(struct StructPbodyIOLean PbodyIOLean){

	StructPbody TemporalPbody;
	memset(&TemporalPbody,0,sizeof(StructPbody));

// omit StructPbodyptr Next; // <TMP> <LEAN>
// omit void    (*Baryon);   // <TMP> <LEAN>

    TemporalPbody.GlobalID = PbodyIOLean.GlobalID;
// omit unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
// omit unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

// omit bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false.
// omit bool    Active; // <LEAN> This is used for the individual time step.

// omit int InteractionList; // <TMP> <LEAN>

    TemporalPbody.Pos[0] = PbodyIOLean.Pos[0];
    TemporalPbody.Pos[1] = PbodyIOLean.Pos[1];
    TemporalPbody.Pos[2] = PbodyIOLean.Pos[2];
// omit double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    TemporalPbody.Vel[0] = PbodyIOLean.Vel[0];
    TemporalPbody.Vel[1] = PbodyIOLean.Vel[1];
    TemporalPbody.Vel[2] = PbodyIOLean.Vel[2];
// omit double    Velh[3];      // <MIX> <LEAN> The \"half-step slided\" velocity of the particle. For leapfrog integration.
    TemporalPbody.Acc[0] = PbodyIOLean.Acc[0];
    TemporalPbody.Acc[1] = PbodyIOLean.Acc[1];
    TemporalPbody.Acc[2] = PbodyIOLean.Acc[2];
// omit double    AccOld[3];    // <LEAN> The acceleration of the particle.
    TemporalPbody.Pot = PbodyIOLean.Pot;
    TemporalPbody.Mass = PbodyIOLean.Mass;
    TemporalPbody.Eps = PbodyIOLean.Eps;

    TemporalPbody.Type = PbodyIOLean.Type;
// omit short     k;            // <TMP> <LEAN>
// omit double    dt;           // <LEAN>
// omit double    EraLocal;  // <TMP> <LEAN>

	return TemporalPbody;

}

struct StructPhydroIOCompact CopyPhydroToTemporalStructureCompact(const int index){

	struct StructPhydroIOCompact TempPhydroIOCompact;

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    TempPhydroIOCompact.Nlist = Phydro[index]->Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


    TempPhydroIOCompact.Leaf = Phydro[index]->Leaf;

    TempPhydroIOCompact.GravityKickFlag = Phydro[index]->GravityKickFlag;
    TempPhydroIOCompact.k_hydro = Phydro[index]->k_hydro;
    TempPhydroIOCompact.dt_hydro = (float)Phydro[index]->dt_hydro;
    TempPhydroIOCompact.EraLocal_hydro = (float)Phydro[index]->EraLocal_hydro;

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    TempPhydroIOCompact.k_hydro_localmin = Phydro[index]->k_hydro_localmin;
    TempPhydroIOCompact.k_hydro_localmin_old = Phydro[index]->k_hydro_localmin_old;
    TempPhydroIOCompact.NextUpdateEra = (float)Phydro[index]->NextUpdateEra;
    TempPhydroIOCompact.dt_hydro_localmin = (float)Phydro[index]->dt_hydro_localmin;


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TempPhydroIOCompact.Rho = (float)Phydro[index]->Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TempPhydroIOCompact.EnergyDensity = (float)Phydro[index]->EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TempPhydroIOCompact.Kernel = (float)Phydro[index]->Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

    TempPhydroIOCompact.Gradh = (float)Phydro[index]->Gradh;

    TempPhydroIOCompact.NumberDensity = (float)Phydro[index]->NumberDensity;
    TempPhydroIOCompact.GradN = (float)Phydro[index]->GradN;
    TempPhydroIOCompact.fij = (float)Phydro[index]->fij;


    TempPhydroIOCompact.Div = (float)Phydro[index]->Div;
    TempPhydroIOCompact.Rot[0] = (float)Phydro[index]->Rot[0];
    TempPhydroIOCompact.Rot[1] = (float)Phydro[index]->Rot[1];
    TempPhydroIOCompact.Rot[2] = (float)Phydro[index]->Rot[2];
    TempPhydroIOCompact.F = (float)Phydro[index]->F;
    TempPhydroIOCompact.HydroAcc[0] = (float)Phydro[index]->HydroAcc[0];
    TempPhydroIOCompact.HydroAcc[1] = (float)Phydro[index]->HydroAcc[1];
    TempPhydroIOCompact.HydroAcc[2] = (float)Phydro[index]->HydroAcc[2];
    TempPhydroIOCompact.U = (float)Phydro[index]->U;
// omit double    UPred;      // <TMP> <LEAN>
    TempPhydroIOCompact.Du = (float)Phydro[index]->Du;
    TempPhydroIOCompact.DuPrev = (float)Phydro[index]->DuPrev;
    TempPhydroIOCompact.DuCooling = (float)Phydro[index]->DuCooling;


    TempPhydroIOCompact.PseudoDensity = (float)Phydro[index]->PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    TempPhydroIOCompact.Zw = (float)Phydro[index]->Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    TempPhydroIOCompact.DZw = (float)Phydro[index]->DZw;
    TempPhydroIOCompact.Ddif = (float)Phydro[index]->Ddif;


    TempPhydroIOCompact.DQheat = (float)Phydro[index]->DQheat;
    TempPhydroIOCompact.dMass = (float)Phydro[index]->dMass;
    TempPhydroIOCompact.Z = (float)Phydro[index]->Z;
    TempPhydroIOCompact.ZII = (float)Phydro[index]->ZII;
    TempPhydroIOCompact.ZIa = (float)Phydro[index]->ZIa;
    TempPhydroIOCompact.dZII = (float)Phydro[index]->dZII;
    TempPhydroIOCompact.dZIa = (float)Phydro[index]->dZIa;
    TempPhydroIOCompact.Vsig = (float)Phydro[index]->Vsig;
    TempPhydroIOCompact.Alpha = (float)Phydro[index]->Alpha;
    TempPhydroIOCompact.DAlpha = (float)Phydro[index]->DAlpha;


    TempPhydroIOCompact.DZdiff = (float)Phydro[index]->DZdiff;
    TempPhydroIOCompact.Sxy = (float)Phydro[index]->Sxy;
    TempPhydroIOCompact.Sxz = (float)Phydro[index]->Sxz;
    TempPhydroIOCompact.Syx = (float)Phydro[index]->Syx;
    TempPhydroIOCompact.Syz = (float)Phydro[index]->Syz;
    TempPhydroIOCompact.Szx = (float)Phydro[index]->Szx;
    TempPhydroIOCompact.Szy = (float)Phydro[index]->Szy;

#if VISCOSITY_TYPE == 1 //{
    TempPhydroIOCompact.Bxx = (float)Phydro[index]->Bxx;
    TempPhydroIOCompact.Bxy = (float)Phydro[index]->Bxy;
    TempPhydroIOCompact.Bxz = (float)Phydro[index]->Bxz;
    TempPhydroIOCompact.Byx = (float)Phydro[index]->Byx;
    TempPhydroIOCompact.Byy = (float)Phydro[index]->Byy;
    TempPhydroIOCompact.Byz = (float)Phydro[index]->Byz;
    TempPhydroIOCompact.Bzx = (float)Phydro[index]->Bzx;
    TempPhydroIOCompact.Bzy = (float)Phydro[index]->Bzy;
    TempPhydroIOCompact.Bzz = (float)Phydro[index]->Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TempPhydroIOCompact.G0 = (float)Phydro[index]->G0;
    TempPhydroIOCompact.fH2 = (float)Phydro[index]->fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    TempPhydroIOCompact.SpawnTimes = Phydro[index]->SpawnTimes;
    TempPhydroIOCompact.SpawnMass = (float)Phydro[index]->SpawnMass;
#endif
#ifdef USE_CELIB
    TempPhydroIOCompact.Elements[0] = (float)Phydro[index]->Elements[0];
    TempPhydroIOCompact.Elements[1] = (float)Phydro[index]->Elements[1];
    TempPhydroIOCompact.Elements[2] = (float)Phydro[index]->Elements[2];
    TempPhydroIOCompact.Elements[3] = (float)Phydro[index]->Elements[3];
    TempPhydroIOCompact.Elements[4] = (float)Phydro[index]->Elements[4];
    TempPhydroIOCompact.Elements[5] = (float)Phydro[index]->Elements[5];
    TempPhydroIOCompact.Elements[6] = (float)Phydro[index]->Elements[6];
    TempPhydroIOCompact.Elements[7] = (float)Phydro[index]->Elements[7];
    TempPhydroIOCompact.Elements[8] = (float)Phydro[index]->Elements[8];
    TempPhydroIOCompact.Elements[9] = (float)Phydro[index]->Elements[9];
    TempPhydroIOCompact.Elements[10] = (float)Phydro[index]->Elements[10];
    TempPhydroIOCompact.Elements[11] = (float)Phydro[index]->Elements[11];
    TempPhydroIOCompact.Elements[12] = (float)Phydro[index]->Elements[12];
    TempPhydroIOCompact.Elements[13] = (float)Phydro[index]->Elements[13];
    TempPhydroIOCompact.Elements[14] = (float)Phydro[index]->Elements[14];
    TempPhydroIOCompact.Elements[15] = (float)Phydro[index]->Elements[15];
    TempPhydroIOCompact.Elements[16] = (float)Phydro[index]->Elements[16];
    TempPhydroIOCompact.Elements[17] = (float)Phydro[index]->Elements[17];
#endif // USE_CELIB

    TempPhydroIOCompact.GradRho = (float)Phydro[index]->GradRho;
    TempPhydroIOCompact.G0thin = (float)Phydro[index]->G0thin;
    TempPhydroIOCompact.G0thinLocal = (float)Phydro[index]->G0thinLocal;
    TempPhydroIOCompact.G0extLocal = (float)Phydro[index]->G0extLocal;
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TempPhydroIOCompact.MomentumRP[0] = (float)Phydro[index]->MomentumRP[0];
    TempPhydroIOCompact.MomentumRP[1] = (float)Phydro[index]->MomentumRP[1];
    TempPhydroIOCompact.MomentumRP[2] = (float)Phydro[index]->MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPhydroIOCompact.MomentumFB[0] = (float)Phydro[index]->MomentumFB[0];
    TempPhydroIOCompact.MomentumFB[1] = (float)Phydro[index]->MomentumFB[1];
    TempPhydroIOCompact.MomentumFB[2] = (float)Phydro[index]->MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TempPhydroIOCompact.Uhot = (float)Phydro[index]->Uhot;
    TempPhydroIOCompact.Rhohot = (float)Phydro[index]->Rhohot;
    TempPhydroIOCompact.Mhot = (float)Phydro[index]->Mhot;
    TempPhydroIOCompact.MultiphaseEndTime = (float)Phydro[index]->MultiphaseEndTime;
    TempPhydroIOCompact.MultiphaseFlag = Phydro[index]->MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TempPhydroIOCompact.SplitFlag = Phydro[index]->SplitFlag;
    TempPhydroIOCompact.SplitGeneration = Phydro[index]->SplitGeneration;
    TempPhydroIOCompact.SplitNthChildren = Phydro[index]->SplitNthChildren;
    TempPhydroIOCompact.AncestorGlobalID = Phydro[index]->AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TempPhydroIOCompact.Pot = (float)Phydro[index]->Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TempPhydroIOCompact.Tag = Phydro[index]->Tag;
#endif

	return TempPhydroIOCompact;

}

struct StructPhydroIOCompact CopyPhydroToTemporalStructureCompactElement(StructPhydroptr const Ph){

	struct StructPhydroIOCompact TempPhydroIOCompact;

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    TempPhydroIOCompact.Nlist = Ph->Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


    TempPhydroIOCompact.Leaf = Ph->Leaf;

    TempPhydroIOCompact.GravityKickFlag = Ph->GravityKickFlag;
    TempPhydroIOCompact.k_hydro = Ph->k_hydro;
    TempPhydroIOCompact.dt_hydro = Ph->dt_hydro;
    TempPhydroIOCompact.EraLocal_hydro = Ph->EraLocal_hydro;

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    TempPhydroIOCompact.k_hydro_localmin = Ph->k_hydro_localmin;
    TempPhydroIOCompact.k_hydro_localmin_old = Ph->k_hydro_localmin_old;
    TempPhydroIOCompact.NextUpdateEra = Ph->NextUpdateEra;
    TempPhydroIOCompact.dt_hydro_localmin = Ph->dt_hydro_localmin;


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TempPhydroIOCompact.Rho = Ph->Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TempPhydroIOCompact.EnergyDensity = Ph->EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TempPhydroIOCompact.Kernel = Ph->Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

    TempPhydroIOCompact.Gradh = Ph->Gradh;

    TempPhydroIOCompact.NumberDensity = Ph->NumberDensity;
    TempPhydroIOCompact.GradN = Ph->GradN;
    TempPhydroIOCompact.fij = Ph->fij;


    TempPhydroIOCompact.Div = Ph->Div;
    TempPhydroIOCompact.Rot[0] = Ph->Rot[0];
    TempPhydroIOCompact.Rot[1] = Ph->Rot[1];
    TempPhydroIOCompact.Rot[2] = Ph->Rot[2];
    TempPhydroIOCompact.F = Ph->F;
    TempPhydroIOCompact.HydroAcc[0] = Ph->HydroAcc[0];
    TempPhydroIOCompact.HydroAcc[1] = Ph->HydroAcc[1];
    TempPhydroIOCompact.HydroAcc[2] = Ph->HydroAcc[2];
    TempPhydroIOCompact.U = Ph->U;
// omit double    UPred;      // <TMP> <LEAN>
    TempPhydroIOCompact.Du = Ph->Du;
    TempPhydroIOCompact.DuPrev = Ph->DuPrev;
    TempPhydroIOCompact.DuCooling = Ph->DuCooling;


    TempPhydroIOCompact.PseudoDensity = Ph->PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    TempPhydroIOCompact.Zw = Ph->Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    TempPhydroIOCompact.DZw = Ph->DZw;
    TempPhydroIOCompact.Ddif = Ph->Ddif;


    TempPhydroIOCompact.DQheat = Ph->DQheat;
    TempPhydroIOCompact.dMass = Ph->dMass;
    TempPhydroIOCompact.Z = Ph->Z;
    TempPhydroIOCompact.ZII = Ph->ZII;
    TempPhydroIOCompact.ZIa = Ph->ZIa;
    TempPhydroIOCompact.dZII = Ph->dZII;
    TempPhydroIOCompact.dZIa = Ph->dZIa;
    TempPhydroIOCompact.Vsig = Ph->Vsig;
    TempPhydroIOCompact.Alpha = Ph->Alpha;
    TempPhydroIOCompact.DAlpha = Ph->DAlpha;


    TempPhydroIOCompact.DZdiff = Ph->DZdiff;
    TempPhydroIOCompact.Sxy = Ph->Sxy;
    TempPhydroIOCompact.Sxz = Ph->Sxz;
    TempPhydroIOCompact.Syx = Ph->Syx;
    TempPhydroIOCompact.Syz = Ph->Syz;
    TempPhydroIOCompact.Szx = Ph->Szx;
    TempPhydroIOCompact.Szy = Ph->Szy;

#if VISCOSITY_TYPE == 1 //{
    TempPhydroIOCompact.Bxx = Ph->Bxx;
    TempPhydroIOCompact.Bxy = Ph->Bxy;
    TempPhydroIOCompact.Bxz = Ph->Bxz;
    TempPhydroIOCompact.Byx = Ph->Byx;
    TempPhydroIOCompact.Byy = Ph->Byy;
    TempPhydroIOCompact.Byz = Ph->Byz;
    TempPhydroIOCompact.Bzx = Ph->Bzx;
    TempPhydroIOCompact.Bzy = Ph->Bzy;
    TempPhydroIOCompact.Bzz = Ph->Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TempPhydroIOCompact.G0 = Ph->G0;
    TempPhydroIOCompact.fH2 = Ph->fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    TempPhydroIOCompact.SpawnTimes = Ph->SpawnTimes;
    TempPhydroIOCompact.SpawnMass = Ph->SpawnMass;
#endif
#ifdef USE_CELIB
    TempPhydroIOCompact.Elements[0] = Ph->Elements[0];
    TempPhydroIOCompact.Elements[1] = Ph->Elements[1];
    TempPhydroIOCompact.Elements[2] = Ph->Elements[2];
    TempPhydroIOCompact.Elements[3] = Ph->Elements[3];
    TempPhydroIOCompact.Elements[4] = Ph->Elements[4];
    TempPhydroIOCompact.Elements[5] = Ph->Elements[5];
    TempPhydroIOCompact.Elements[6] = Ph->Elements[6];
    TempPhydroIOCompact.Elements[7] = Ph->Elements[7];
    TempPhydroIOCompact.Elements[8] = Ph->Elements[8];
    TempPhydroIOCompact.Elements[9] = Ph->Elements[9];
    TempPhydroIOCompact.Elements[10] = Ph->Elements[10];
    TempPhydroIOCompact.Elements[11] = Ph->Elements[11];
    TempPhydroIOCompact.Elements[12] = Ph->Elements[12];
    TempPhydroIOCompact.Elements[13] = Ph->Elements[13];
    TempPhydroIOCompact.Elements[14] = Ph->Elements[14];
    TempPhydroIOCompact.Elements[15] = Ph->Elements[15];
    TempPhydroIOCompact.Elements[16] = Ph->Elements[16];
    TempPhydroIOCompact.Elements[17] = Ph->Elements[17];
#endif // USE_CELIB

    TempPhydroIOCompact.GradRho = Ph->GradRho;
    TempPhydroIOCompact.G0thin = Ph->G0thin;
    TempPhydroIOCompact.G0thinLocal = Ph->G0thinLocal;
    TempPhydroIOCompact.G0extLocal = Ph->G0extLocal;
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TempPhydroIOCompact.MomentumRP[0] = Ph->MomentumRP[0];
    TempPhydroIOCompact.MomentumRP[1] = Ph->MomentumRP[1];
    TempPhydroIOCompact.MomentumRP[2] = Ph->MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPhydroIOCompact.MomentumFB[0] = Ph->MomentumFB[0];
    TempPhydroIOCompact.MomentumFB[1] = Ph->MomentumFB[1];
    TempPhydroIOCompact.MomentumFB[2] = Ph->MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TempPhydroIOCompact.Uhot = Ph->Uhot;
    TempPhydroIOCompact.Rhohot = Ph->Rhohot;
    TempPhydroIOCompact.Mhot = Ph->Mhot;
    TempPhydroIOCompact.MultiphaseEndTime = Ph->MultiphaseEndTime;
    TempPhydroIOCompact.MultiphaseFlag = Ph->MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TempPhydroIOCompact.SplitFlag = Ph->SplitFlag;
    TempPhydroIOCompact.SplitGeneration = Ph->SplitGeneration;
    TempPhydroIOCompact.SplitNthChildren = Ph->SplitNthChildren;
    TempPhydroIOCompact.AncestorGlobalID = Ph->AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TempPhydroIOCompact.Pot = Ph->Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TempPhydroIOCompact.Tag = Ph->Tag;
#endif

	return TempPhydroIOCompact;

}

StructPhydro CopyTemporalStructureCompactToPhydroCompact(struct StructPhydroIOCompact PhydroIOCompact){

	StructPhydro TemporalPhydro;
	memset(&TemporalPhydro,0,sizeof(StructPhydro));

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    TemporalPhydro.Nlist = PhydroIOCompact.Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


    TemporalPhydro.Leaf = PhydroIOCompact.Leaf;

    TemporalPhydro.GravityKickFlag = PhydroIOCompact.GravityKickFlag;
    TemporalPhydro.k_hydro = PhydroIOCompact.k_hydro;
    TemporalPhydro.dt_hydro = PhydroIOCompact.dt_hydro;
    TemporalPhydro.EraLocal_hydro = PhydroIOCompact.EraLocal_hydro;

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    TemporalPhydro.k_hydro_localmin = PhydroIOCompact.k_hydro_localmin;
    TemporalPhydro.k_hydro_localmin_old = PhydroIOCompact.k_hydro_localmin_old;
    TemporalPhydro.NextUpdateEra = PhydroIOCompact.NextUpdateEra;
    TemporalPhydro.dt_hydro_localmin = PhydroIOCompact.dt_hydro_localmin;


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TemporalPhydro.Rho = PhydroIOCompact.Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TemporalPhydro.EnergyDensity = PhydroIOCompact.EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TemporalPhydro.Kernel = PhydroIOCompact.Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

    TemporalPhydro.Gradh = PhydroIOCompact.Gradh;

    TemporalPhydro.NumberDensity = PhydroIOCompact.NumberDensity;
    TemporalPhydro.GradN = PhydroIOCompact.GradN;
    TemporalPhydro.fij = PhydroIOCompact.fij;


    TemporalPhydro.Div = PhydroIOCompact.Div;
    TemporalPhydro.Rot[0] = PhydroIOCompact.Rot[0];
    TemporalPhydro.Rot[1] = PhydroIOCompact.Rot[1];
    TemporalPhydro.Rot[2] = PhydroIOCompact.Rot[2];
    TemporalPhydro.F = PhydroIOCompact.F;
    TemporalPhydro.HydroAcc[0] = PhydroIOCompact.HydroAcc[0];
    TemporalPhydro.HydroAcc[1] = PhydroIOCompact.HydroAcc[1];
    TemporalPhydro.HydroAcc[2] = PhydroIOCompact.HydroAcc[2];
    TemporalPhydro.U = PhydroIOCompact.U;
// omit double    UPred;      // <TMP> <LEAN>
    TemporalPhydro.Du = PhydroIOCompact.Du;
    TemporalPhydro.DuPrev = PhydroIOCompact.DuPrev;
    TemporalPhydro.DuCooling = PhydroIOCompact.DuCooling;


    TemporalPhydro.PseudoDensity = PhydroIOCompact.PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    TemporalPhydro.Zw = PhydroIOCompact.Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    TemporalPhydro.DZw = PhydroIOCompact.DZw;
    TemporalPhydro.Ddif = PhydroIOCompact.Ddif;


    TemporalPhydro.DQheat = PhydroIOCompact.DQheat;
    TemporalPhydro.dMass = PhydroIOCompact.dMass;
    TemporalPhydro.Z = PhydroIOCompact.Z;
    TemporalPhydro.ZII = PhydroIOCompact.ZII;
    TemporalPhydro.ZIa = PhydroIOCompact.ZIa;
    TemporalPhydro.dZII = PhydroIOCompact.dZII;
    TemporalPhydro.dZIa = PhydroIOCompact.dZIa;
    TemporalPhydro.Vsig = PhydroIOCompact.Vsig;
    TemporalPhydro.Alpha = PhydroIOCompact.Alpha;
    TemporalPhydro.DAlpha = PhydroIOCompact.DAlpha;


    TemporalPhydro.DZdiff = PhydroIOCompact.DZdiff;
    TemporalPhydro.Sxy = PhydroIOCompact.Sxy;
    TemporalPhydro.Sxz = PhydroIOCompact.Sxz;
    TemporalPhydro.Syx = PhydroIOCompact.Syx;
    TemporalPhydro.Syz = PhydroIOCompact.Syz;
    TemporalPhydro.Szx = PhydroIOCompact.Szx;
    TemporalPhydro.Szy = PhydroIOCompact.Szy;

#if VISCOSITY_TYPE == 1 //{
    TemporalPhydro.Bxx = PhydroIOCompact.Bxx;
    TemporalPhydro.Bxy = PhydroIOCompact.Bxy;
    TemporalPhydro.Bxz = PhydroIOCompact.Bxz;
    TemporalPhydro.Byx = PhydroIOCompact.Byx;
    TemporalPhydro.Byy = PhydroIOCompact.Byy;
    TemporalPhydro.Byz = PhydroIOCompact.Byz;
    TemporalPhydro.Bzx = PhydroIOCompact.Bzx;
    TemporalPhydro.Bzy = PhydroIOCompact.Bzy;
    TemporalPhydro.Bzz = PhydroIOCompact.Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TemporalPhydro.G0 = PhydroIOCompact.G0;
    TemporalPhydro.fH2 = PhydroIOCompact.fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    TemporalPhydro.SpawnTimes = PhydroIOCompact.SpawnTimes;
    TemporalPhydro.SpawnMass = PhydroIOCompact.SpawnMass;
#endif
#ifdef USE_CELIB
    TemporalPhydro.Elements[0] = PhydroIOCompact.Elements[0];
    TemporalPhydro.Elements[1] = PhydroIOCompact.Elements[1];
    TemporalPhydro.Elements[2] = PhydroIOCompact.Elements[2];
    TemporalPhydro.Elements[3] = PhydroIOCompact.Elements[3];
    TemporalPhydro.Elements[4] = PhydroIOCompact.Elements[4];
    TemporalPhydro.Elements[5] = PhydroIOCompact.Elements[5];
    TemporalPhydro.Elements[6] = PhydroIOCompact.Elements[6];
    TemporalPhydro.Elements[7] = PhydroIOCompact.Elements[7];
    TemporalPhydro.Elements[8] = PhydroIOCompact.Elements[8];
    TemporalPhydro.Elements[9] = PhydroIOCompact.Elements[9];
    TemporalPhydro.Elements[10] = PhydroIOCompact.Elements[10];
    TemporalPhydro.Elements[11] = PhydroIOCompact.Elements[11];
    TemporalPhydro.Elements[12] = PhydroIOCompact.Elements[12];
    TemporalPhydro.Elements[13] = PhydroIOCompact.Elements[13];
    TemporalPhydro.Elements[14] = PhydroIOCompact.Elements[14];
    TemporalPhydro.Elements[15] = PhydroIOCompact.Elements[15];
    TemporalPhydro.Elements[16] = PhydroIOCompact.Elements[16];
    TemporalPhydro.Elements[17] = PhydroIOCompact.Elements[17];
#endif // USE_CELIB

    TemporalPhydro.GradRho = PhydroIOCompact.GradRho;
    TemporalPhydro.G0thin = PhydroIOCompact.G0thin;
    TemporalPhydro.G0thinLocal = PhydroIOCompact.G0thinLocal;
    TemporalPhydro.G0extLocal = PhydroIOCompact.G0extLocal;
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TemporalPhydro.MomentumRP[0] = PhydroIOCompact.MomentumRP[0];
    TemporalPhydro.MomentumRP[1] = PhydroIOCompact.MomentumRP[1];
    TemporalPhydro.MomentumRP[2] = PhydroIOCompact.MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TemporalPhydro.MomentumFB[0] = PhydroIOCompact.MomentumFB[0];
    TemporalPhydro.MomentumFB[1] = PhydroIOCompact.MomentumFB[1];
    TemporalPhydro.MomentumFB[2] = PhydroIOCompact.MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TemporalPhydro.Uhot = PhydroIOCompact.Uhot;
    TemporalPhydro.Rhohot = PhydroIOCompact.Rhohot;
    TemporalPhydro.Mhot = PhydroIOCompact.Mhot;
    TemporalPhydro.MultiphaseEndTime = PhydroIOCompact.MultiphaseEndTime;
    TemporalPhydro.MultiphaseFlag = PhydroIOCompact.MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPhydro.SplitFlag = PhydroIOCompact.SplitFlag;
    TemporalPhydro.SplitGeneration = PhydroIOCompact.SplitGeneration;
    TemporalPhydro.SplitNthChildren = PhydroIOCompact.SplitNthChildren;
    TemporalPhydro.AncestorGlobalID = PhydroIOCompact.AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TemporalPhydro.Pot = PhydroIOCompact.Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TemporalPhydro.Tag = PhydroIOCompact.Tag;
#endif

	return TemporalPhydro;

}

struct StructPhydroIOCompactDouble CopyPhydroToTemporalStructureCompactDouble(const int index){

	struct StructPhydroIOCompactDouble TempPhydroIOCompactDouble;

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    TempPhydroIOCompactDouble.Nlist = Phydro[index]->Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


    TempPhydroIOCompactDouble.Leaf = Phydro[index]->Leaf;

    TempPhydroIOCompactDouble.GravityKickFlag = Phydro[index]->GravityKickFlag;
    TempPhydroIOCompactDouble.k_hydro = Phydro[index]->k_hydro;
    TempPhydroIOCompactDouble.dt_hydro = (double)Phydro[index]->dt_hydro;
    TempPhydroIOCompactDouble.EraLocal_hydro = (double)Phydro[index]->EraLocal_hydro;

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    TempPhydroIOCompactDouble.k_hydro_localmin = Phydro[index]->k_hydro_localmin;
    TempPhydroIOCompactDouble.k_hydro_localmin_old = Phydro[index]->k_hydro_localmin_old;
    TempPhydroIOCompactDouble.NextUpdateEra = (double)Phydro[index]->NextUpdateEra;
    TempPhydroIOCompactDouble.dt_hydro_localmin = (double)Phydro[index]->dt_hydro_localmin;


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TempPhydroIOCompactDouble.Rho = (double)Phydro[index]->Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TempPhydroIOCompactDouble.EnergyDensity = (double)Phydro[index]->EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TempPhydroIOCompactDouble.Kernel = (double)Phydro[index]->Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

    TempPhydroIOCompactDouble.Gradh = (double)Phydro[index]->Gradh;

    TempPhydroIOCompactDouble.NumberDensity = (double)Phydro[index]->NumberDensity;
    TempPhydroIOCompactDouble.GradN = (double)Phydro[index]->GradN;
    TempPhydroIOCompactDouble.fij = (double)Phydro[index]->fij;


    TempPhydroIOCompactDouble.Div = (double)Phydro[index]->Div;
    TempPhydroIOCompactDouble.Rot[0] = (double)Phydro[index]->Rot[0];
    TempPhydroIOCompactDouble.Rot[1] = (double)Phydro[index]->Rot[1];
    TempPhydroIOCompactDouble.Rot[2] = (double)Phydro[index]->Rot[2];
    TempPhydroIOCompactDouble.F = (double)Phydro[index]->F;
    TempPhydroIOCompactDouble.HydroAcc[0] = (double)Phydro[index]->HydroAcc[0];
    TempPhydroIOCompactDouble.HydroAcc[1] = (double)Phydro[index]->HydroAcc[1];
    TempPhydroIOCompactDouble.HydroAcc[2] = (double)Phydro[index]->HydroAcc[2];
    TempPhydroIOCompactDouble.U = (double)Phydro[index]->U;
// omit double    UPred;      // <TMP> <LEAN>
    TempPhydroIOCompactDouble.Du = (double)Phydro[index]->Du;
    TempPhydroIOCompactDouble.DuPrev = (double)Phydro[index]->DuPrev;
    TempPhydroIOCompactDouble.DuCooling = (double)Phydro[index]->DuCooling;


    TempPhydroIOCompactDouble.PseudoDensity = (double)Phydro[index]->PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    TempPhydroIOCompactDouble.Zw = (double)Phydro[index]->Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    TempPhydroIOCompactDouble.DZw = (double)Phydro[index]->DZw;
    TempPhydroIOCompactDouble.Ddif = (double)Phydro[index]->Ddif;


    TempPhydroIOCompactDouble.DQheat = (double)Phydro[index]->DQheat;
    TempPhydroIOCompactDouble.dMass = (double)Phydro[index]->dMass;
    TempPhydroIOCompactDouble.Z = (double)Phydro[index]->Z;
    TempPhydroIOCompactDouble.ZII = (double)Phydro[index]->ZII;
    TempPhydroIOCompactDouble.ZIa = (double)Phydro[index]->ZIa;
    TempPhydroIOCompactDouble.dZII = (double)Phydro[index]->dZII;
    TempPhydroIOCompactDouble.dZIa = (double)Phydro[index]->dZIa;
    TempPhydroIOCompactDouble.Vsig = (double)Phydro[index]->Vsig;
    TempPhydroIOCompactDouble.Alpha = (double)Phydro[index]->Alpha;
    TempPhydroIOCompactDouble.DAlpha = (double)Phydro[index]->DAlpha;


    TempPhydroIOCompactDouble.DZdiff = (double)Phydro[index]->DZdiff;
    TempPhydroIOCompactDouble.Sxy = (double)Phydro[index]->Sxy;
    TempPhydroIOCompactDouble.Sxz = (double)Phydro[index]->Sxz;
    TempPhydroIOCompactDouble.Syx = (double)Phydro[index]->Syx;
    TempPhydroIOCompactDouble.Syz = (double)Phydro[index]->Syz;
    TempPhydroIOCompactDouble.Szx = (double)Phydro[index]->Szx;
    TempPhydroIOCompactDouble.Szy = (double)Phydro[index]->Szy;

#if VISCOSITY_TYPE == 1 //{
    TempPhydroIOCompactDouble.Bxx = (double)Phydro[index]->Bxx;
    TempPhydroIOCompactDouble.Bxy = (double)Phydro[index]->Bxy;
    TempPhydroIOCompactDouble.Bxz = (double)Phydro[index]->Bxz;
    TempPhydroIOCompactDouble.Byx = (double)Phydro[index]->Byx;
    TempPhydroIOCompactDouble.Byy = (double)Phydro[index]->Byy;
    TempPhydroIOCompactDouble.Byz = (double)Phydro[index]->Byz;
    TempPhydroIOCompactDouble.Bzx = (double)Phydro[index]->Bzx;
    TempPhydroIOCompactDouble.Bzy = (double)Phydro[index]->Bzy;
    TempPhydroIOCompactDouble.Bzz = (double)Phydro[index]->Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TempPhydroIOCompactDouble.G0 = (double)Phydro[index]->G0;
    TempPhydroIOCompactDouble.fH2 = (double)Phydro[index]->fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    TempPhydroIOCompactDouble.SpawnTimes = Phydro[index]->SpawnTimes;
    TempPhydroIOCompactDouble.SpawnMass = (double)Phydro[index]->SpawnMass;
#endif
#ifdef USE_CELIB
    TempPhydroIOCompactDouble.Elements[0] = (double)Phydro[index]->Elements[0];
    TempPhydroIOCompactDouble.Elements[1] = (double)Phydro[index]->Elements[1];
    TempPhydroIOCompactDouble.Elements[2] = (double)Phydro[index]->Elements[2];
    TempPhydroIOCompactDouble.Elements[3] = (double)Phydro[index]->Elements[3];
    TempPhydroIOCompactDouble.Elements[4] = (double)Phydro[index]->Elements[4];
    TempPhydroIOCompactDouble.Elements[5] = (double)Phydro[index]->Elements[5];
    TempPhydroIOCompactDouble.Elements[6] = (double)Phydro[index]->Elements[6];
    TempPhydroIOCompactDouble.Elements[7] = (double)Phydro[index]->Elements[7];
    TempPhydroIOCompactDouble.Elements[8] = (double)Phydro[index]->Elements[8];
    TempPhydroIOCompactDouble.Elements[9] = (double)Phydro[index]->Elements[9];
    TempPhydroIOCompactDouble.Elements[10] = (double)Phydro[index]->Elements[10];
    TempPhydroIOCompactDouble.Elements[11] = (double)Phydro[index]->Elements[11];
    TempPhydroIOCompactDouble.Elements[12] = (double)Phydro[index]->Elements[12];
    TempPhydroIOCompactDouble.Elements[13] = (double)Phydro[index]->Elements[13];
    TempPhydroIOCompactDouble.Elements[14] = (double)Phydro[index]->Elements[14];
    TempPhydroIOCompactDouble.Elements[15] = (double)Phydro[index]->Elements[15];
    TempPhydroIOCompactDouble.Elements[16] = (double)Phydro[index]->Elements[16];
    TempPhydroIOCompactDouble.Elements[17] = (double)Phydro[index]->Elements[17];
#endif // USE_CELIB

    TempPhydroIOCompactDouble.GradRho = (double)Phydro[index]->GradRho;
    TempPhydroIOCompactDouble.G0thin = (double)Phydro[index]->G0thin;
    TempPhydroIOCompactDouble.G0thinLocal = (double)Phydro[index]->G0thinLocal;
    TempPhydroIOCompactDouble.G0extLocal = (double)Phydro[index]->G0extLocal;
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TempPhydroIOCompactDouble.MomentumRP[0] = (double)Phydro[index]->MomentumRP[0];
    TempPhydroIOCompactDouble.MomentumRP[1] = (double)Phydro[index]->MomentumRP[1];
    TempPhydroIOCompactDouble.MomentumRP[2] = (double)Phydro[index]->MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPhydroIOCompactDouble.MomentumFB[0] = (double)Phydro[index]->MomentumFB[0];
    TempPhydroIOCompactDouble.MomentumFB[1] = (double)Phydro[index]->MomentumFB[1];
    TempPhydroIOCompactDouble.MomentumFB[2] = (double)Phydro[index]->MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TempPhydroIOCompactDouble.Uhot = (double)Phydro[index]->Uhot;
    TempPhydroIOCompactDouble.Rhohot = (double)Phydro[index]->Rhohot;
    TempPhydroIOCompactDouble.Mhot = (double)Phydro[index]->Mhot;
    TempPhydroIOCompactDouble.MultiphaseEndTime = (double)Phydro[index]->MultiphaseEndTime;
    TempPhydroIOCompactDouble.MultiphaseFlag = Phydro[index]->MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TempPhydroIOCompactDouble.SplitFlag = Phydro[index]->SplitFlag;
    TempPhydroIOCompactDouble.SplitGeneration = Phydro[index]->SplitGeneration;
    TempPhydroIOCompactDouble.SplitNthChildren = Phydro[index]->SplitNthChildren;
    TempPhydroIOCompactDouble.AncestorGlobalID = Phydro[index]->AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TempPhydroIOCompactDouble.Pot = (double)Phydro[index]->Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TempPhydroIOCompactDouble.Tag = Phydro[index]->Tag;
#endif

	return TempPhydroIOCompactDouble;

}

struct StructPhydroIOCompactDouble CopyPhydroToTemporalStructureCompactDoubleElement(StructPhydroptr const Ph){

	struct StructPhydroIOCompactDouble TempPhydroIOCompactDouble;

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    TempPhydroIOCompactDouble.Nlist = Ph->Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


    TempPhydroIOCompactDouble.Leaf = Ph->Leaf;

    TempPhydroIOCompactDouble.GravityKickFlag = Ph->GravityKickFlag;
    TempPhydroIOCompactDouble.k_hydro = Ph->k_hydro;
    TempPhydroIOCompactDouble.dt_hydro = Ph->dt_hydro;
    TempPhydroIOCompactDouble.EraLocal_hydro = Ph->EraLocal_hydro;

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    TempPhydroIOCompactDouble.k_hydro_localmin = Ph->k_hydro_localmin;
    TempPhydroIOCompactDouble.k_hydro_localmin_old = Ph->k_hydro_localmin_old;
    TempPhydroIOCompactDouble.NextUpdateEra = Ph->NextUpdateEra;
    TempPhydroIOCompactDouble.dt_hydro_localmin = Ph->dt_hydro_localmin;


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TempPhydroIOCompactDouble.Rho = Ph->Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TempPhydroIOCompactDouble.EnergyDensity = Ph->EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TempPhydroIOCompactDouble.Kernel = Ph->Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

    TempPhydroIOCompactDouble.Gradh = Ph->Gradh;

    TempPhydroIOCompactDouble.NumberDensity = Ph->NumberDensity;
    TempPhydroIOCompactDouble.GradN = Ph->GradN;
    TempPhydroIOCompactDouble.fij = Ph->fij;


    TempPhydroIOCompactDouble.Div = Ph->Div;
    TempPhydroIOCompactDouble.Rot[0] = Ph->Rot[0];
    TempPhydroIOCompactDouble.Rot[1] = Ph->Rot[1];
    TempPhydroIOCompactDouble.Rot[2] = Ph->Rot[2];
    TempPhydroIOCompactDouble.F = Ph->F;
    TempPhydroIOCompactDouble.HydroAcc[0] = Ph->HydroAcc[0];
    TempPhydroIOCompactDouble.HydroAcc[1] = Ph->HydroAcc[1];
    TempPhydroIOCompactDouble.HydroAcc[2] = Ph->HydroAcc[2];
    TempPhydroIOCompactDouble.U = Ph->U;
// omit double    UPred;      // <TMP> <LEAN>
    TempPhydroIOCompactDouble.Du = Ph->Du;
    TempPhydroIOCompactDouble.DuPrev = Ph->DuPrev;
    TempPhydroIOCompactDouble.DuCooling = Ph->DuCooling;


    TempPhydroIOCompactDouble.PseudoDensity = Ph->PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    TempPhydroIOCompactDouble.Zw = Ph->Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    TempPhydroIOCompactDouble.DZw = Ph->DZw;
    TempPhydroIOCompactDouble.Ddif = Ph->Ddif;


    TempPhydroIOCompactDouble.DQheat = Ph->DQheat;
    TempPhydroIOCompactDouble.dMass = Ph->dMass;
    TempPhydroIOCompactDouble.Z = Ph->Z;
    TempPhydroIOCompactDouble.ZII = Ph->ZII;
    TempPhydroIOCompactDouble.ZIa = Ph->ZIa;
    TempPhydroIOCompactDouble.dZII = Ph->dZII;
    TempPhydroIOCompactDouble.dZIa = Ph->dZIa;
    TempPhydroIOCompactDouble.Vsig = Ph->Vsig;
    TempPhydroIOCompactDouble.Alpha = Ph->Alpha;
    TempPhydroIOCompactDouble.DAlpha = Ph->DAlpha;


    TempPhydroIOCompactDouble.DZdiff = Ph->DZdiff;
    TempPhydroIOCompactDouble.Sxy = Ph->Sxy;
    TempPhydroIOCompactDouble.Sxz = Ph->Sxz;
    TempPhydroIOCompactDouble.Syx = Ph->Syx;
    TempPhydroIOCompactDouble.Syz = Ph->Syz;
    TempPhydroIOCompactDouble.Szx = Ph->Szx;
    TempPhydroIOCompactDouble.Szy = Ph->Szy;

#if VISCOSITY_TYPE == 1 //{
    TempPhydroIOCompactDouble.Bxx = Ph->Bxx;
    TempPhydroIOCompactDouble.Bxy = Ph->Bxy;
    TempPhydroIOCompactDouble.Bxz = Ph->Bxz;
    TempPhydroIOCompactDouble.Byx = Ph->Byx;
    TempPhydroIOCompactDouble.Byy = Ph->Byy;
    TempPhydroIOCompactDouble.Byz = Ph->Byz;
    TempPhydroIOCompactDouble.Bzx = Ph->Bzx;
    TempPhydroIOCompactDouble.Bzy = Ph->Bzy;
    TempPhydroIOCompactDouble.Bzz = Ph->Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TempPhydroIOCompactDouble.G0 = Ph->G0;
    TempPhydroIOCompactDouble.fH2 = Ph->fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    TempPhydroIOCompactDouble.SpawnTimes = Ph->SpawnTimes;
    TempPhydroIOCompactDouble.SpawnMass = Ph->SpawnMass;
#endif
#ifdef USE_CELIB
    TempPhydroIOCompactDouble.Elements[0] = Ph->Elements[0];
    TempPhydroIOCompactDouble.Elements[1] = Ph->Elements[1];
    TempPhydroIOCompactDouble.Elements[2] = Ph->Elements[2];
    TempPhydroIOCompactDouble.Elements[3] = Ph->Elements[3];
    TempPhydroIOCompactDouble.Elements[4] = Ph->Elements[4];
    TempPhydroIOCompactDouble.Elements[5] = Ph->Elements[5];
    TempPhydroIOCompactDouble.Elements[6] = Ph->Elements[6];
    TempPhydroIOCompactDouble.Elements[7] = Ph->Elements[7];
    TempPhydroIOCompactDouble.Elements[8] = Ph->Elements[8];
    TempPhydroIOCompactDouble.Elements[9] = Ph->Elements[9];
    TempPhydroIOCompactDouble.Elements[10] = Ph->Elements[10];
    TempPhydroIOCompactDouble.Elements[11] = Ph->Elements[11];
    TempPhydroIOCompactDouble.Elements[12] = Ph->Elements[12];
    TempPhydroIOCompactDouble.Elements[13] = Ph->Elements[13];
    TempPhydroIOCompactDouble.Elements[14] = Ph->Elements[14];
    TempPhydroIOCompactDouble.Elements[15] = Ph->Elements[15];
    TempPhydroIOCompactDouble.Elements[16] = Ph->Elements[16];
    TempPhydroIOCompactDouble.Elements[17] = Ph->Elements[17];
#endif // USE_CELIB

    TempPhydroIOCompactDouble.GradRho = Ph->GradRho;
    TempPhydroIOCompactDouble.G0thin = Ph->G0thin;
    TempPhydroIOCompactDouble.G0thinLocal = Ph->G0thinLocal;
    TempPhydroIOCompactDouble.G0extLocal = Ph->G0extLocal;
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TempPhydroIOCompactDouble.MomentumRP[0] = Ph->MomentumRP[0];
    TempPhydroIOCompactDouble.MomentumRP[1] = Ph->MomentumRP[1];
    TempPhydroIOCompactDouble.MomentumRP[2] = Ph->MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPhydroIOCompactDouble.MomentumFB[0] = Ph->MomentumFB[0];
    TempPhydroIOCompactDouble.MomentumFB[1] = Ph->MomentumFB[1];
    TempPhydroIOCompactDouble.MomentumFB[2] = Ph->MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TempPhydroIOCompactDouble.Uhot = Ph->Uhot;
    TempPhydroIOCompactDouble.Rhohot = Ph->Rhohot;
    TempPhydroIOCompactDouble.Mhot = Ph->Mhot;
    TempPhydroIOCompactDouble.MultiphaseEndTime = Ph->MultiphaseEndTime;
    TempPhydroIOCompactDouble.MultiphaseFlag = Ph->MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TempPhydroIOCompactDouble.SplitFlag = Ph->SplitFlag;
    TempPhydroIOCompactDouble.SplitGeneration = Ph->SplitGeneration;
    TempPhydroIOCompactDouble.SplitNthChildren = Ph->SplitNthChildren;
    TempPhydroIOCompactDouble.AncestorGlobalID = Ph->AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TempPhydroIOCompactDouble.Pot = Ph->Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TempPhydroIOCompactDouble.Tag = Ph->Tag;
#endif

	return TempPhydroIOCompactDouble;

}

StructPhydro CopyTemporalStructureCompactDoubleToPhydroCompactDouble(struct StructPhydroIOCompactDouble PhydroIOCompactDouble){

	StructPhydro TemporalPhydro;
	memset(&TemporalPhydro,0,sizeof(StructPhydro));

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    TemporalPhydro.Nlist = PhydroIOCompactDouble.Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


    TemporalPhydro.Leaf = PhydroIOCompactDouble.Leaf;

    TemporalPhydro.GravityKickFlag = PhydroIOCompactDouble.GravityKickFlag;
    TemporalPhydro.k_hydro = PhydroIOCompactDouble.k_hydro;
    TemporalPhydro.dt_hydro = PhydroIOCompactDouble.dt_hydro;
    TemporalPhydro.EraLocal_hydro = PhydroIOCompactDouble.EraLocal_hydro;

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    TemporalPhydro.k_hydro_localmin = PhydroIOCompactDouble.k_hydro_localmin;
    TemporalPhydro.k_hydro_localmin_old = PhydroIOCompactDouble.k_hydro_localmin_old;
    TemporalPhydro.NextUpdateEra = PhydroIOCompactDouble.NextUpdateEra;
    TemporalPhydro.dt_hydro_localmin = PhydroIOCompactDouble.dt_hydro_localmin;


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TemporalPhydro.Rho = PhydroIOCompactDouble.Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TemporalPhydro.EnergyDensity = PhydroIOCompactDouble.EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TemporalPhydro.Kernel = PhydroIOCompactDouble.Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

    TemporalPhydro.Gradh = PhydroIOCompactDouble.Gradh;

    TemporalPhydro.NumberDensity = PhydroIOCompactDouble.NumberDensity;
    TemporalPhydro.GradN = PhydroIOCompactDouble.GradN;
    TemporalPhydro.fij = PhydroIOCompactDouble.fij;


    TemporalPhydro.Div = PhydroIOCompactDouble.Div;
    TemporalPhydro.Rot[0] = PhydroIOCompactDouble.Rot[0];
    TemporalPhydro.Rot[1] = PhydroIOCompactDouble.Rot[1];
    TemporalPhydro.Rot[2] = PhydroIOCompactDouble.Rot[2];
    TemporalPhydro.F = PhydroIOCompactDouble.F;
    TemporalPhydro.HydroAcc[0] = PhydroIOCompactDouble.HydroAcc[0];
    TemporalPhydro.HydroAcc[1] = PhydroIOCompactDouble.HydroAcc[1];
    TemporalPhydro.HydroAcc[2] = PhydroIOCompactDouble.HydroAcc[2];
    TemporalPhydro.U = PhydroIOCompactDouble.U;
// omit double    UPred;      // <TMP> <LEAN>
    TemporalPhydro.Du = PhydroIOCompactDouble.Du;
    TemporalPhydro.DuPrev = PhydroIOCompactDouble.DuPrev;
    TemporalPhydro.DuCooling = PhydroIOCompactDouble.DuCooling;


    TemporalPhydro.PseudoDensity = PhydroIOCompactDouble.PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    TemporalPhydro.Zw = PhydroIOCompactDouble.Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    TemporalPhydro.DZw = PhydroIOCompactDouble.DZw;
    TemporalPhydro.Ddif = PhydroIOCompactDouble.Ddif;


    TemporalPhydro.DQheat = PhydroIOCompactDouble.DQheat;
    TemporalPhydro.dMass = PhydroIOCompactDouble.dMass;
    TemporalPhydro.Z = PhydroIOCompactDouble.Z;
    TemporalPhydro.ZII = PhydroIOCompactDouble.ZII;
    TemporalPhydro.ZIa = PhydroIOCompactDouble.ZIa;
    TemporalPhydro.dZII = PhydroIOCompactDouble.dZII;
    TemporalPhydro.dZIa = PhydroIOCompactDouble.dZIa;
    TemporalPhydro.Vsig = PhydroIOCompactDouble.Vsig;
    TemporalPhydro.Alpha = PhydroIOCompactDouble.Alpha;
    TemporalPhydro.DAlpha = PhydroIOCompactDouble.DAlpha;


    TemporalPhydro.DZdiff = PhydroIOCompactDouble.DZdiff;
    TemporalPhydro.Sxy = PhydroIOCompactDouble.Sxy;
    TemporalPhydro.Sxz = PhydroIOCompactDouble.Sxz;
    TemporalPhydro.Syx = PhydroIOCompactDouble.Syx;
    TemporalPhydro.Syz = PhydroIOCompactDouble.Syz;
    TemporalPhydro.Szx = PhydroIOCompactDouble.Szx;
    TemporalPhydro.Szy = PhydroIOCompactDouble.Szy;

#if VISCOSITY_TYPE == 1 //{
    TemporalPhydro.Bxx = PhydroIOCompactDouble.Bxx;
    TemporalPhydro.Bxy = PhydroIOCompactDouble.Bxy;
    TemporalPhydro.Bxz = PhydroIOCompactDouble.Bxz;
    TemporalPhydro.Byx = PhydroIOCompactDouble.Byx;
    TemporalPhydro.Byy = PhydroIOCompactDouble.Byy;
    TemporalPhydro.Byz = PhydroIOCompactDouble.Byz;
    TemporalPhydro.Bzx = PhydroIOCompactDouble.Bzx;
    TemporalPhydro.Bzy = PhydroIOCompactDouble.Bzy;
    TemporalPhydro.Bzz = PhydroIOCompactDouble.Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TemporalPhydro.G0 = PhydroIOCompactDouble.G0;
    TemporalPhydro.fH2 = PhydroIOCompactDouble.fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    TemporalPhydro.SpawnTimes = PhydroIOCompactDouble.SpawnTimes;
    TemporalPhydro.SpawnMass = PhydroIOCompactDouble.SpawnMass;
#endif
#ifdef USE_CELIB
    TemporalPhydro.Elements[0] = PhydroIOCompactDouble.Elements[0];
    TemporalPhydro.Elements[1] = PhydroIOCompactDouble.Elements[1];
    TemporalPhydro.Elements[2] = PhydroIOCompactDouble.Elements[2];
    TemporalPhydro.Elements[3] = PhydroIOCompactDouble.Elements[3];
    TemporalPhydro.Elements[4] = PhydroIOCompactDouble.Elements[4];
    TemporalPhydro.Elements[5] = PhydroIOCompactDouble.Elements[5];
    TemporalPhydro.Elements[6] = PhydroIOCompactDouble.Elements[6];
    TemporalPhydro.Elements[7] = PhydroIOCompactDouble.Elements[7];
    TemporalPhydro.Elements[8] = PhydroIOCompactDouble.Elements[8];
    TemporalPhydro.Elements[9] = PhydroIOCompactDouble.Elements[9];
    TemporalPhydro.Elements[10] = PhydroIOCompactDouble.Elements[10];
    TemporalPhydro.Elements[11] = PhydroIOCompactDouble.Elements[11];
    TemporalPhydro.Elements[12] = PhydroIOCompactDouble.Elements[12];
    TemporalPhydro.Elements[13] = PhydroIOCompactDouble.Elements[13];
    TemporalPhydro.Elements[14] = PhydroIOCompactDouble.Elements[14];
    TemporalPhydro.Elements[15] = PhydroIOCompactDouble.Elements[15];
    TemporalPhydro.Elements[16] = PhydroIOCompactDouble.Elements[16];
    TemporalPhydro.Elements[17] = PhydroIOCompactDouble.Elements[17];
#endif // USE_CELIB

    TemporalPhydro.GradRho = PhydroIOCompactDouble.GradRho;
    TemporalPhydro.G0thin = PhydroIOCompactDouble.G0thin;
    TemporalPhydro.G0thinLocal = PhydroIOCompactDouble.G0thinLocal;
    TemporalPhydro.G0extLocal = PhydroIOCompactDouble.G0extLocal;
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TemporalPhydro.MomentumRP[0] = PhydroIOCompactDouble.MomentumRP[0];
    TemporalPhydro.MomentumRP[1] = PhydroIOCompactDouble.MomentumRP[1];
    TemporalPhydro.MomentumRP[2] = PhydroIOCompactDouble.MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TemporalPhydro.MomentumFB[0] = PhydroIOCompactDouble.MomentumFB[0];
    TemporalPhydro.MomentumFB[1] = PhydroIOCompactDouble.MomentumFB[1];
    TemporalPhydro.MomentumFB[2] = PhydroIOCompactDouble.MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TemporalPhydro.Uhot = PhydroIOCompactDouble.Uhot;
    TemporalPhydro.Rhohot = PhydroIOCompactDouble.Rhohot;
    TemporalPhydro.Mhot = PhydroIOCompactDouble.Mhot;
    TemporalPhydro.MultiphaseEndTime = PhydroIOCompactDouble.MultiphaseEndTime;
    TemporalPhydro.MultiphaseFlag = PhydroIOCompactDouble.MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPhydro.SplitFlag = PhydroIOCompactDouble.SplitFlag;
    TemporalPhydro.SplitGeneration = PhydroIOCompactDouble.SplitGeneration;
    TemporalPhydro.SplitNthChildren = PhydroIOCompactDouble.SplitNthChildren;
    TemporalPhydro.AncestorGlobalID = PhydroIOCompactDouble.AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TemporalPhydro.Pot = PhydroIOCompactDouble.Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TemporalPhydro.Tag = PhydroIOCompactDouble.Tag;
#endif

	return TemporalPhydro;

}

struct StructPhydroIOCompactMix CopyPhydroToTemporalStructureCompactMix(const int index){

	struct StructPhydroIOCompactMix TempPhydroIOCompactMix;

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    TempPhydroIOCompactMix.Nlist = Phydro[index]->Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


    TempPhydroIOCompactMix.Leaf = Phydro[index]->Leaf;

    TempPhydroIOCompactMix.GravityKickFlag = Phydro[index]->GravityKickFlag;
    TempPhydroIOCompactMix.k_hydro = Phydro[index]->k_hydro;
    TempPhydroIOCompactMix.dt_hydro = (double)Phydro[index]->dt_hydro;
    TempPhydroIOCompactMix.EraLocal_hydro = (double)Phydro[index]->EraLocal_hydro;

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    TempPhydroIOCompactMix.k_hydro_localmin = Phydro[index]->k_hydro_localmin;
    TempPhydroIOCompactMix.k_hydro_localmin_old = Phydro[index]->k_hydro_localmin_old;
    TempPhydroIOCompactMix.NextUpdateEra = (double)Phydro[index]->NextUpdateEra;
    TempPhydroIOCompactMix.dt_hydro_localmin = (double)Phydro[index]->dt_hydro_localmin;


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TempPhydroIOCompactMix.Rho = (double)Phydro[index]->Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TempPhydroIOCompactMix.EnergyDensity = (double)Phydro[index]->EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TempPhydroIOCompactMix.Kernel = (double)Phydro[index]->Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

    TempPhydroIOCompactMix.Gradh = (double)Phydro[index]->Gradh;

    TempPhydroIOCompactMix.NumberDensity = (double)Phydro[index]->NumberDensity;
    TempPhydroIOCompactMix.GradN = (double)Phydro[index]->GradN;
    TempPhydroIOCompactMix.fij = (double)Phydro[index]->fij;


    TempPhydroIOCompactMix.Div = (double)Phydro[index]->Div;
    TempPhydroIOCompactMix.Rot[0] = (double)Phydro[index]->Rot[0];
    TempPhydroIOCompactMix.Rot[1] = (double)Phydro[index]->Rot[1];
    TempPhydroIOCompactMix.Rot[2] = (double)Phydro[index]->Rot[2];
    TempPhydroIOCompactMix.F = (double)Phydro[index]->F;
    TempPhydroIOCompactMix.HydroAcc[0] = (double)Phydro[index]->HydroAcc[0];
    TempPhydroIOCompactMix.HydroAcc[1] = (double)Phydro[index]->HydroAcc[1];
    TempPhydroIOCompactMix.HydroAcc[2] = (double)Phydro[index]->HydroAcc[2];
    TempPhydroIOCompactMix.U = (double)Phydro[index]->U;
// omit double    UPred;      // <TMP> <LEAN>
    TempPhydroIOCompactMix.Du = (double)Phydro[index]->Du;
    TempPhydroIOCompactMix.DuPrev = (double)Phydro[index]->DuPrev;
    TempPhydroIOCompactMix.DuCooling = (double)Phydro[index]->DuCooling;


    TempPhydroIOCompactMix.PseudoDensity = (double)Phydro[index]->PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    TempPhydroIOCompactMix.Zw = (double)Phydro[index]->Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    TempPhydroIOCompactMix.DZw = (double)Phydro[index]->DZw;
    TempPhydroIOCompactMix.Ddif = (double)Phydro[index]->Ddif;


    TempPhydroIOCompactMix.DQheat = (double)Phydro[index]->DQheat;
    TempPhydroIOCompactMix.dMass = (double)Phydro[index]->dMass;
    TempPhydroIOCompactMix.Z = (double)Phydro[index]->Z;
    TempPhydroIOCompactMix.ZII = (double)Phydro[index]->ZII;
    TempPhydroIOCompactMix.ZIa = (double)Phydro[index]->ZIa;
    TempPhydroIOCompactMix.dZII = (double)Phydro[index]->dZII;
    TempPhydroIOCompactMix.dZIa = (double)Phydro[index]->dZIa;
    TempPhydroIOCompactMix.Vsig = (double)Phydro[index]->Vsig;
    TempPhydroIOCompactMix.Alpha = (double)Phydro[index]->Alpha;
    TempPhydroIOCompactMix.DAlpha = (double)Phydro[index]->DAlpha;


    TempPhydroIOCompactMix.DZdiff = (double)Phydro[index]->DZdiff;
    TempPhydroIOCompactMix.Sxy = (double)Phydro[index]->Sxy;
    TempPhydroIOCompactMix.Sxz = (double)Phydro[index]->Sxz;
    TempPhydroIOCompactMix.Syx = (double)Phydro[index]->Syx;
    TempPhydroIOCompactMix.Syz = (double)Phydro[index]->Syz;
    TempPhydroIOCompactMix.Szx = (double)Phydro[index]->Szx;
    TempPhydroIOCompactMix.Szy = (double)Phydro[index]->Szy;

#if VISCOSITY_TYPE == 1 //{
    TempPhydroIOCompactMix.Bxx = (double)Phydro[index]->Bxx;
    TempPhydroIOCompactMix.Bxy = (double)Phydro[index]->Bxy;
    TempPhydroIOCompactMix.Bxz = (double)Phydro[index]->Bxz;
    TempPhydroIOCompactMix.Byx = (double)Phydro[index]->Byx;
    TempPhydroIOCompactMix.Byy = (double)Phydro[index]->Byy;
    TempPhydroIOCompactMix.Byz = (double)Phydro[index]->Byz;
    TempPhydroIOCompactMix.Bzx = (double)Phydro[index]->Bzx;
    TempPhydroIOCompactMix.Bzy = (double)Phydro[index]->Bzy;
    TempPhydroIOCompactMix.Bzz = (double)Phydro[index]->Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TempPhydroIOCompactMix.G0 = (double)Phydro[index]->G0;
    TempPhydroIOCompactMix.fH2 = (double)Phydro[index]->fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    TempPhydroIOCompactMix.SpawnTimes = Phydro[index]->SpawnTimes;
    TempPhydroIOCompactMix.SpawnMass = (double)Phydro[index]->SpawnMass;
#endif
#ifdef USE_CELIB
    TempPhydroIOCompactMix.Elements[0] = (float)Phydro[index]->Elements[0];
    TempPhydroIOCompactMix.Elements[1] = (float)Phydro[index]->Elements[1];
    TempPhydroIOCompactMix.Elements[2] = (float)Phydro[index]->Elements[2];
    TempPhydroIOCompactMix.Elements[3] = (float)Phydro[index]->Elements[3];
    TempPhydroIOCompactMix.Elements[4] = (float)Phydro[index]->Elements[4];
    TempPhydroIOCompactMix.Elements[5] = (float)Phydro[index]->Elements[5];
    TempPhydroIOCompactMix.Elements[6] = (float)Phydro[index]->Elements[6];
    TempPhydroIOCompactMix.Elements[7] = (float)Phydro[index]->Elements[7];
    TempPhydroIOCompactMix.Elements[8] = (float)Phydro[index]->Elements[8];
    TempPhydroIOCompactMix.Elements[9] = (float)Phydro[index]->Elements[9];
    TempPhydroIOCompactMix.Elements[10] = (float)Phydro[index]->Elements[10];
    TempPhydroIOCompactMix.Elements[11] = (float)Phydro[index]->Elements[11];
    TempPhydroIOCompactMix.Elements[12] = (float)Phydro[index]->Elements[12];
    TempPhydroIOCompactMix.Elements[13] = (float)Phydro[index]->Elements[13];
    TempPhydroIOCompactMix.Elements[14] = (float)Phydro[index]->Elements[14];
    TempPhydroIOCompactMix.Elements[15] = (float)Phydro[index]->Elements[15];
    TempPhydroIOCompactMix.Elements[16] = (float)Phydro[index]->Elements[16];
    TempPhydroIOCompactMix.Elements[17] = (float)Phydro[index]->Elements[17];
#endif // USE_CELIB

    TempPhydroIOCompactMix.GradRho = (double)Phydro[index]->GradRho;
    TempPhydroIOCompactMix.G0thin = (double)Phydro[index]->G0thin;
    TempPhydroIOCompactMix.G0thinLocal = (double)Phydro[index]->G0thinLocal;
    TempPhydroIOCompactMix.G0extLocal = (double)Phydro[index]->G0extLocal;
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TempPhydroIOCompactMix.MomentumRP[0] = (double)Phydro[index]->MomentumRP[0];
    TempPhydroIOCompactMix.MomentumRP[1] = (double)Phydro[index]->MomentumRP[1];
    TempPhydroIOCompactMix.MomentumRP[2] = (double)Phydro[index]->MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPhydroIOCompactMix.MomentumFB[0] = (double)Phydro[index]->MomentumFB[0];
    TempPhydroIOCompactMix.MomentumFB[1] = (double)Phydro[index]->MomentumFB[1];
    TempPhydroIOCompactMix.MomentumFB[2] = (double)Phydro[index]->MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TempPhydroIOCompactMix.Uhot = (double)Phydro[index]->Uhot;
    TempPhydroIOCompactMix.Rhohot = (double)Phydro[index]->Rhohot;
    TempPhydroIOCompactMix.Mhot = (double)Phydro[index]->Mhot;
    TempPhydroIOCompactMix.MultiphaseEndTime = (double)Phydro[index]->MultiphaseEndTime;
    TempPhydroIOCompactMix.MultiphaseFlag = Phydro[index]->MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TempPhydroIOCompactMix.SplitFlag = Phydro[index]->SplitFlag;
    TempPhydroIOCompactMix.SplitGeneration = Phydro[index]->SplitGeneration;
    TempPhydroIOCompactMix.SplitNthChildren = Phydro[index]->SplitNthChildren;
    TempPhydroIOCompactMix.AncestorGlobalID = Phydro[index]->AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TempPhydroIOCompactMix.Pot = (double)Phydro[index]->Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TempPhydroIOCompactMix.Tag = Phydro[index]->Tag;
#endif

	return TempPhydroIOCompactMix;

}

struct StructPhydroIOCompactMix CopyPhydroToTemporalStructureCompactMixElement(StructPhydroptr const Ph){

	struct StructPhydroIOCompactMix TempPhydroIOCompactMix;

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    TempPhydroIOCompactMix.Nlist = Ph->Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


    TempPhydroIOCompactMix.Leaf = Ph->Leaf;

    TempPhydroIOCompactMix.GravityKickFlag = Ph->GravityKickFlag;
    TempPhydroIOCompactMix.k_hydro = Ph->k_hydro;
    TempPhydroIOCompactMix.dt_hydro = Ph->dt_hydro;
    TempPhydroIOCompactMix.EraLocal_hydro = Ph->EraLocal_hydro;

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    TempPhydroIOCompactMix.k_hydro_localmin = Ph->k_hydro_localmin;
    TempPhydroIOCompactMix.k_hydro_localmin_old = Ph->k_hydro_localmin_old;
    TempPhydroIOCompactMix.NextUpdateEra = Ph->NextUpdateEra;
    TempPhydroIOCompactMix.dt_hydro_localmin = Ph->dt_hydro_localmin;


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TempPhydroIOCompactMix.Rho = Ph->Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TempPhydroIOCompactMix.EnergyDensity = Ph->EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TempPhydroIOCompactMix.Kernel = Ph->Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

    TempPhydroIOCompactMix.Gradh = Ph->Gradh;

    TempPhydroIOCompactMix.NumberDensity = Ph->NumberDensity;
    TempPhydroIOCompactMix.GradN = Ph->GradN;
    TempPhydroIOCompactMix.fij = Ph->fij;


    TempPhydroIOCompactMix.Div = Ph->Div;
    TempPhydroIOCompactMix.Rot[0] = Ph->Rot[0];
    TempPhydroIOCompactMix.Rot[1] = Ph->Rot[1];
    TempPhydroIOCompactMix.Rot[2] = Ph->Rot[2];
    TempPhydroIOCompactMix.F = Ph->F;
    TempPhydroIOCompactMix.HydroAcc[0] = Ph->HydroAcc[0];
    TempPhydroIOCompactMix.HydroAcc[1] = Ph->HydroAcc[1];
    TempPhydroIOCompactMix.HydroAcc[2] = Ph->HydroAcc[2];
    TempPhydroIOCompactMix.U = Ph->U;
// omit double    UPred;      // <TMP> <LEAN>
    TempPhydroIOCompactMix.Du = Ph->Du;
    TempPhydroIOCompactMix.DuPrev = Ph->DuPrev;
    TempPhydroIOCompactMix.DuCooling = Ph->DuCooling;


    TempPhydroIOCompactMix.PseudoDensity = Ph->PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    TempPhydroIOCompactMix.Zw = Ph->Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    TempPhydroIOCompactMix.DZw = Ph->DZw;
    TempPhydroIOCompactMix.Ddif = Ph->Ddif;


    TempPhydroIOCompactMix.DQheat = Ph->DQheat;
    TempPhydroIOCompactMix.dMass = Ph->dMass;
    TempPhydroIOCompactMix.Z = Ph->Z;
    TempPhydroIOCompactMix.ZII = Ph->ZII;
    TempPhydroIOCompactMix.ZIa = Ph->ZIa;
    TempPhydroIOCompactMix.dZII = Ph->dZII;
    TempPhydroIOCompactMix.dZIa = Ph->dZIa;
    TempPhydroIOCompactMix.Vsig = Ph->Vsig;
    TempPhydroIOCompactMix.Alpha = Ph->Alpha;
    TempPhydroIOCompactMix.DAlpha = Ph->DAlpha;


    TempPhydroIOCompactMix.DZdiff = Ph->DZdiff;
    TempPhydroIOCompactMix.Sxy = Ph->Sxy;
    TempPhydroIOCompactMix.Sxz = Ph->Sxz;
    TempPhydroIOCompactMix.Syx = Ph->Syx;
    TempPhydroIOCompactMix.Syz = Ph->Syz;
    TempPhydroIOCompactMix.Szx = Ph->Szx;
    TempPhydroIOCompactMix.Szy = Ph->Szy;

#if VISCOSITY_TYPE == 1 //{
    TempPhydroIOCompactMix.Bxx = Ph->Bxx;
    TempPhydroIOCompactMix.Bxy = Ph->Bxy;
    TempPhydroIOCompactMix.Bxz = Ph->Bxz;
    TempPhydroIOCompactMix.Byx = Ph->Byx;
    TempPhydroIOCompactMix.Byy = Ph->Byy;
    TempPhydroIOCompactMix.Byz = Ph->Byz;
    TempPhydroIOCompactMix.Bzx = Ph->Bzx;
    TempPhydroIOCompactMix.Bzy = Ph->Bzy;
    TempPhydroIOCompactMix.Bzz = Ph->Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TempPhydroIOCompactMix.G0 = Ph->G0;
    TempPhydroIOCompactMix.fH2 = Ph->fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    TempPhydroIOCompactMix.SpawnTimes = Ph->SpawnTimes;
    TempPhydroIOCompactMix.SpawnMass = Ph->SpawnMass;
#endif
#ifdef USE_CELIB
    TempPhydroIOCompactMix.Elements[0] = Ph->Elements[0];
    TempPhydroIOCompactMix.Elements[1] = Ph->Elements[1];
    TempPhydroIOCompactMix.Elements[2] = Ph->Elements[2];
    TempPhydroIOCompactMix.Elements[3] = Ph->Elements[3];
    TempPhydroIOCompactMix.Elements[4] = Ph->Elements[4];
    TempPhydroIOCompactMix.Elements[5] = Ph->Elements[5];
    TempPhydroIOCompactMix.Elements[6] = Ph->Elements[6];
    TempPhydroIOCompactMix.Elements[7] = Ph->Elements[7];
    TempPhydroIOCompactMix.Elements[8] = Ph->Elements[8];
    TempPhydroIOCompactMix.Elements[9] = Ph->Elements[9];
    TempPhydroIOCompactMix.Elements[10] = Ph->Elements[10];
    TempPhydroIOCompactMix.Elements[11] = Ph->Elements[11];
    TempPhydroIOCompactMix.Elements[12] = Ph->Elements[12];
    TempPhydroIOCompactMix.Elements[13] = Ph->Elements[13];
    TempPhydroIOCompactMix.Elements[14] = Ph->Elements[14];
    TempPhydroIOCompactMix.Elements[15] = Ph->Elements[15];
    TempPhydroIOCompactMix.Elements[16] = Ph->Elements[16];
    TempPhydroIOCompactMix.Elements[17] = Ph->Elements[17];
#endif // USE_CELIB

    TempPhydroIOCompactMix.GradRho = Ph->GradRho;
    TempPhydroIOCompactMix.G0thin = Ph->G0thin;
    TempPhydroIOCompactMix.G0thinLocal = Ph->G0thinLocal;
    TempPhydroIOCompactMix.G0extLocal = Ph->G0extLocal;
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TempPhydroIOCompactMix.MomentumRP[0] = Ph->MomentumRP[0];
    TempPhydroIOCompactMix.MomentumRP[1] = Ph->MomentumRP[1];
    TempPhydroIOCompactMix.MomentumRP[2] = Ph->MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPhydroIOCompactMix.MomentumFB[0] = Ph->MomentumFB[0];
    TempPhydroIOCompactMix.MomentumFB[1] = Ph->MomentumFB[1];
    TempPhydroIOCompactMix.MomentumFB[2] = Ph->MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TempPhydroIOCompactMix.Uhot = Ph->Uhot;
    TempPhydroIOCompactMix.Rhohot = Ph->Rhohot;
    TempPhydroIOCompactMix.Mhot = Ph->Mhot;
    TempPhydroIOCompactMix.MultiphaseEndTime = Ph->MultiphaseEndTime;
    TempPhydroIOCompactMix.MultiphaseFlag = Ph->MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TempPhydroIOCompactMix.SplitFlag = Ph->SplitFlag;
    TempPhydroIOCompactMix.SplitGeneration = Ph->SplitGeneration;
    TempPhydroIOCompactMix.SplitNthChildren = Ph->SplitNthChildren;
    TempPhydroIOCompactMix.AncestorGlobalID = Ph->AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TempPhydroIOCompactMix.Pot = Ph->Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TempPhydroIOCompactMix.Tag = Ph->Tag;
#endif

	return TempPhydroIOCompactMix;

}

StructPhydro CopyTemporalStructureCompactMixToPhydroCompactMix(struct StructPhydroIOCompactMix PhydroIOCompactMix){

	StructPhydro TemporalPhydro;
	memset(&TemporalPhydro,0,sizeof(StructPhydro));

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    TemporalPhydro.Nlist = PhydroIOCompactMix.Nlist;
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


    TemporalPhydro.Leaf = PhydroIOCompactMix.Leaf;

    TemporalPhydro.GravityKickFlag = PhydroIOCompactMix.GravityKickFlag;
    TemporalPhydro.k_hydro = PhydroIOCompactMix.k_hydro;
    TemporalPhydro.dt_hydro = PhydroIOCompactMix.dt_hydro;
    TemporalPhydro.EraLocal_hydro = PhydroIOCompactMix.EraLocal_hydro;

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    TemporalPhydro.k_hydro_localmin = PhydroIOCompactMix.k_hydro_localmin;
    TemporalPhydro.k_hydro_localmin_old = PhydroIOCompactMix.k_hydro_localmin_old;
    TemporalPhydro.NextUpdateEra = PhydroIOCompactMix.NextUpdateEra;
    TemporalPhydro.dt_hydro_localmin = PhydroIOCompactMix.dt_hydro_localmin;


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TemporalPhydro.Rho = PhydroIOCompactMix.Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TemporalPhydro.EnergyDensity = PhydroIOCompactMix.EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TemporalPhydro.Kernel = PhydroIOCompactMix.Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

    TemporalPhydro.Gradh = PhydroIOCompactMix.Gradh;

    TemporalPhydro.NumberDensity = PhydroIOCompactMix.NumberDensity;
    TemporalPhydro.GradN = PhydroIOCompactMix.GradN;
    TemporalPhydro.fij = PhydroIOCompactMix.fij;


    TemporalPhydro.Div = PhydroIOCompactMix.Div;
    TemporalPhydro.Rot[0] = PhydroIOCompactMix.Rot[0];
    TemporalPhydro.Rot[1] = PhydroIOCompactMix.Rot[1];
    TemporalPhydro.Rot[2] = PhydroIOCompactMix.Rot[2];
    TemporalPhydro.F = PhydroIOCompactMix.F;
    TemporalPhydro.HydroAcc[0] = PhydroIOCompactMix.HydroAcc[0];
    TemporalPhydro.HydroAcc[1] = PhydroIOCompactMix.HydroAcc[1];
    TemporalPhydro.HydroAcc[2] = PhydroIOCompactMix.HydroAcc[2];
    TemporalPhydro.U = PhydroIOCompactMix.U;
// omit double    UPred;      // <TMP> <LEAN>
    TemporalPhydro.Du = PhydroIOCompactMix.Du;
    TemporalPhydro.DuPrev = PhydroIOCompactMix.DuPrev;
    TemporalPhydro.DuCooling = PhydroIOCompactMix.DuCooling;


    TemporalPhydro.PseudoDensity = PhydroIOCompactMix.PseudoDensity;
// omit double    PseudoDensityPred;      // <TMP> The predictor of y
    TemporalPhydro.Zw = PhydroIOCompactMix.Zw;
// omit double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    TemporalPhydro.DZw = PhydroIOCompactMix.DZw;
    TemporalPhydro.Ddif = PhydroIOCompactMix.Ddif;


    TemporalPhydro.DQheat = PhydroIOCompactMix.DQheat;
    TemporalPhydro.dMass = PhydroIOCompactMix.dMass;
    TemporalPhydro.Z = PhydroIOCompactMix.Z;
    TemporalPhydro.ZII = PhydroIOCompactMix.ZII;
    TemporalPhydro.ZIa = PhydroIOCompactMix.ZIa;
    TemporalPhydro.dZII = PhydroIOCompactMix.dZII;
    TemporalPhydro.dZIa = PhydroIOCompactMix.dZIa;
    TemporalPhydro.Vsig = PhydroIOCompactMix.Vsig;
    TemporalPhydro.Alpha = PhydroIOCompactMix.Alpha;
    TemporalPhydro.DAlpha = PhydroIOCompactMix.DAlpha;


    TemporalPhydro.DZdiff = PhydroIOCompactMix.DZdiff;
    TemporalPhydro.Sxy = PhydroIOCompactMix.Sxy;
    TemporalPhydro.Sxz = PhydroIOCompactMix.Sxz;
    TemporalPhydro.Syx = PhydroIOCompactMix.Syx;
    TemporalPhydro.Syz = PhydroIOCompactMix.Syz;
    TemporalPhydro.Szx = PhydroIOCompactMix.Szx;
    TemporalPhydro.Szy = PhydroIOCompactMix.Szy;

#if VISCOSITY_TYPE == 1 //{
    TemporalPhydro.Bxx = PhydroIOCompactMix.Bxx;
    TemporalPhydro.Bxy = PhydroIOCompactMix.Bxy;
    TemporalPhydro.Bxz = PhydroIOCompactMix.Bxz;
    TemporalPhydro.Byx = PhydroIOCompactMix.Byx;
    TemporalPhydro.Byy = PhydroIOCompactMix.Byy;
    TemporalPhydro.Byz = PhydroIOCompactMix.Byz;
    TemporalPhydro.Bzx = PhydroIOCompactMix.Bzx;
    TemporalPhydro.Bzy = PhydroIOCompactMix.Bzy;
    TemporalPhydro.Bzz = PhydroIOCompactMix.Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TemporalPhydro.G0 = PhydroIOCompactMix.G0;
    TemporalPhydro.fH2 = PhydroIOCompactMix.fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    TemporalPhydro.SpawnTimes = PhydroIOCompactMix.SpawnTimes;
    TemporalPhydro.SpawnMass = PhydroIOCompactMix.SpawnMass;
#endif
#ifdef USE_CELIB
    TemporalPhydro.Elements[0] = PhydroIOCompactMix.Elements[0];
    TemporalPhydro.Elements[1] = PhydroIOCompactMix.Elements[1];
    TemporalPhydro.Elements[2] = PhydroIOCompactMix.Elements[2];
    TemporalPhydro.Elements[3] = PhydroIOCompactMix.Elements[3];
    TemporalPhydro.Elements[4] = PhydroIOCompactMix.Elements[4];
    TemporalPhydro.Elements[5] = PhydroIOCompactMix.Elements[5];
    TemporalPhydro.Elements[6] = PhydroIOCompactMix.Elements[6];
    TemporalPhydro.Elements[7] = PhydroIOCompactMix.Elements[7];
    TemporalPhydro.Elements[8] = PhydroIOCompactMix.Elements[8];
    TemporalPhydro.Elements[9] = PhydroIOCompactMix.Elements[9];
    TemporalPhydro.Elements[10] = PhydroIOCompactMix.Elements[10];
    TemporalPhydro.Elements[11] = PhydroIOCompactMix.Elements[11];
    TemporalPhydro.Elements[12] = PhydroIOCompactMix.Elements[12];
    TemporalPhydro.Elements[13] = PhydroIOCompactMix.Elements[13];
    TemporalPhydro.Elements[14] = PhydroIOCompactMix.Elements[14];
    TemporalPhydro.Elements[15] = PhydroIOCompactMix.Elements[15];
    TemporalPhydro.Elements[16] = PhydroIOCompactMix.Elements[16];
    TemporalPhydro.Elements[17] = PhydroIOCompactMix.Elements[17];
#endif // USE_CELIB

    TemporalPhydro.GradRho = PhydroIOCompactMix.GradRho;
    TemporalPhydro.G0thin = PhydroIOCompactMix.G0thin;
    TemporalPhydro.G0thinLocal = PhydroIOCompactMix.G0thinLocal;
    TemporalPhydro.G0extLocal = PhydroIOCompactMix.G0extLocal;
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TemporalPhydro.MomentumRP[0] = PhydroIOCompactMix.MomentumRP[0];
    TemporalPhydro.MomentumRP[1] = PhydroIOCompactMix.MomentumRP[1];
    TemporalPhydro.MomentumRP[2] = PhydroIOCompactMix.MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TemporalPhydro.MomentumFB[0] = PhydroIOCompactMix.MomentumFB[0];
    TemporalPhydro.MomentumFB[1] = PhydroIOCompactMix.MomentumFB[1];
    TemporalPhydro.MomentumFB[2] = PhydroIOCompactMix.MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TemporalPhydro.Uhot = PhydroIOCompactMix.Uhot;
    TemporalPhydro.Rhohot = PhydroIOCompactMix.Rhohot;
    TemporalPhydro.Mhot = PhydroIOCompactMix.Mhot;
    TemporalPhydro.MultiphaseEndTime = PhydroIOCompactMix.MultiphaseEndTime;
    TemporalPhydro.MultiphaseFlag = PhydroIOCompactMix.MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPhydro.SplitFlag = PhydroIOCompactMix.SplitFlag;
    TemporalPhydro.SplitGeneration = PhydroIOCompactMix.SplitGeneration;
    TemporalPhydro.SplitNthChildren = PhydroIOCompactMix.SplitNthChildren;
    TemporalPhydro.AncestorGlobalID = PhydroIOCompactMix.AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TemporalPhydro.Pot = PhydroIOCompactMix.Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TemporalPhydro.Tag = PhydroIOCompactMix.Tag;
#endif

	return TemporalPhydro;

}

struct StructPhydroIOLean CopyPhydroToTemporalStructureLean(const int index){

	struct StructPhydroIOLean TempPhydroIOLean;

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
// omit unsigned int Nlist;      // <LEAN>
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


// omit int Leaf; // <LEAN> The ID to the HydroTree leaves.

// omit bool      GravityKickFlag; // <LEAN> The kick flag of gravitataional acceleration.
// omit int       k_hydro;         // <LEAN>
// omit double    dt_hydro;        // <LEAN>
// omit double    EraLocal_hydro;  // <LEAN> The local time in a current block.

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
// omit int       k_hydro_localmin;  // <LEAN>
// omit int       k_hydro_localmin_old;  // <LEAN>
// omit double    NextUpdateEra;     // <LEAN> The next update timing in Era. This variable is used when the HydroTimeStepLimiterFlag is true.
// omit double    dt_hydro_localmin; // <LEAN>


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TempPhydroIOLean.Rho = (float)Phydro[index]->Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TempPhydroIOLean.EnergyDensity = (float)Phydro[index]->EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TempPhydroIOLean.Kernel = (float)Phydro[index]->Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

// omit double    Gradh; // <LEAN> Grad-h term

// omit double    NumberDensity; // <LEAN> Number density
// omit double    GradN;         // <LEAN> Grad-N term
// omit double    fij;           // <LEAN> Part of fij. fij = 1-fij/(Uj)


// omit double    Div;        // <LEAN>
// omit double    Rot[3];     // <LEAN>
// omit double    F;          // <LEAN>
    TempPhydroIOLean.HydroAcc[0] = (float)Phydro[index]->HydroAcc[0];
    TempPhydroIOLean.HydroAcc[1] = (float)Phydro[index]->HydroAcc[1];
    TempPhydroIOLean.HydroAcc[2] = (float)Phydro[index]->HydroAcc[2];
    TempPhydroIOLean.U = (float)Phydro[index]->U;
// omit double    UPred;      // <TMP> <LEAN>
    TempPhydroIOLean.Du = (float)Phydro[index]->Du;
// omit double    DuPrev;     // <LEAN> Du/Dt.
// omit double    DuCooling;  // <LEAN> The amount of the energy loss by the radiative cooling.


    TempPhydroIOLean.PseudoDensity = (float)Phydro[index]->PseudoDensity;
    TempPhydroIOLean.PseudoDensityPred = (float)Phydro[index]->PseudoDensityPred;
    TempPhydroIOLean.Zw = (float)Phydro[index]->Zw;
    TempPhydroIOLean.ZwPred = (float)Phydro[index]->ZwPred;
    TempPhydroIOLean.DZw = (float)Phydro[index]->DZw;
    TempPhydroIOLean.Ddif = (float)Phydro[index]->Ddif;


// omit double    DQheat;     // <LEAN> The amount of the heating energy by SNe, which injects during dt.
// omit double    dMass;      // <LEAN> The additional mass by SNe.
    TempPhydroIOLean.Z = (float)Phydro[index]->Z;
    TempPhydroIOLean.ZII = (float)Phydro[index]->ZII;
    TempPhydroIOLean.ZIa = (float)Phydro[index]->ZIa;
    TempPhydroIOLean.dZII = (float)Phydro[index]->dZII;
    TempPhydroIOLean.dZIa = (float)Phydro[index]->dZIa;
// omit double    Vsig;       // <LEAN> Signal Velocity  or Maximum Velocity difference.
    TempPhydroIOLean.Alpha = (float)Phydro[index]->Alpha;
// omit double    DAlpha;     // <LEAN>


    TempPhydroIOLean.DZdiff = (float)Phydro[index]->DZdiff;
    TempPhydroIOLean.Sxy = (float)Phydro[index]->Sxy;
    TempPhydroIOLean.Sxz = (float)Phydro[index]->Sxz;
    TempPhydroIOLean.Syx = (float)Phydro[index]->Syx;
    TempPhydroIOLean.Syz = (float)Phydro[index]->Syz;
    TempPhydroIOLean.Szx = (float)Phydro[index]->Szx;
    TempPhydroIOLean.Szy = (float)Phydro[index]->Szy;

#if VISCOSITY_TYPE == 1 //{
    TempPhydroIOLean.Bxx = (float)Phydro[index]->Bxx;
    TempPhydroIOLean.Bxy = (float)Phydro[index]->Bxy;
    TempPhydroIOLean.Bxz = (float)Phydro[index]->Bxz;
    TempPhydroIOLean.Byx = (float)Phydro[index]->Byx;
    TempPhydroIOLean.Byy = (float)Phydro[index]->Byy;
    TempPhydroIOLean.Byz = (float)Phydro[index]->Byz;
    TempPhydroIOLean.Bzx = (float)Phydro[index]->Bzx;
    TempPhydroIOLean.Bzy = (float)Phydro[index]->Bzy;
    TempPhydroIOLean.Bzz = (float)Phydro[index]->Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TempPhydroIOLean.G0 = (float)Phydro[index]->G0;
    TempPhydroIOLean.fH2 = (float)Phydro[index]->fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
// omit short   SpawnTimes; // <LEAN>
// omit double  SpawnMass;  // <LEAN>
#endif
#ifdef USE_CELIB
    TempPhydroIOLean.Elements[0] = (float)Phydro[index]->Elements[0];
    TempPhydroIOLean.Elements[1] = (float)Phydro[index]->Elements[1];
    TempPhydroIOLean.Elements[2] = (float)Phydro[index]->Elements[2];
    TempPhydroIOLean.Elements[3] = (float)Phydro[index]->Elements[3];
    TempPhydroIOLean.Elements[4] = (float)Phydro[index]->Elements[4];
    TempPhydroIOLean.Elements[5] = (float)Phydro[index]->Elements[5];
    TempPhydroIOLean.Elements[6] = (float)Phydro[index]->Elements[6];
    TempPhydroIOLean.Elements[7] = (float)Phydro[index]->Elements[7];
    TempPhydroIOLean.Elements[8] = (float)Phydro[index]->Elements[8];
    TempPhydroIOLean.Elements[9] = (float)Phydro[index]->Elements[9];
    TempPhydroIOLean.Elements[10] = (float)Phydro[index]->Elements[10];
    TempPhydroIOLean.Elements[11] = (float)Phydro[index]->Elements[11];
    TempPhydroIOLean.Elements[12] = (float)Phydro[index]->Elements[12];
    TempPhydroIOLean.Elements[13] = (float)Phydro[index]->Elements[13];
    TempPhydroIOLean.Elements[14] = (float)Phydro[index]->Elements[14];
    TempPhydroIOLean.Elements[15] = (float)Phydro[index]->Elements[15];
    TempPhydroIOLean.Elements[16] = (float)Phydro[index]->Elements[16];
    TempPhydroIOLean.Elements[17] = (float)Phydro[index]->Elements[17];
#endif // USE_CELIB

// omit double GradRho;               // <LEAN>
// omit double G0thin;                // <LEAN>
// omit double G0thinLocal;           // <LEAN>
// omit double G0extLocal;            // <LEAN>
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TempPhydroIOLean.MomentumRP[0] = (float)Phydro[index]->MomentumRP[0];
    TempPhydroIOLean.MomentumRP[1] = (float)Phydro[index]->MomentumRP[1];
    TempPhydroIOLean.MomentumRP[2] = (float)Phydro[index]->MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPhydroIOLean.MomentumFB[0] = (float)Phydro[index]->MomentumFB[0];
    TempPhydroIOLean.MomentumFB[1] = (float)Phydro[index]->MomentumFB[1];
    TempPhydroIOLean.MomentumFB[2] = (float)Phydro[index]->MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TempPhydroIOLean.Uhot = (float)Phydro[index]->Uhot;
    TempPhydroIOLean.Rhohot = (float)Phydro[index]->Rhohot;
    TempPhydroIOLean.Mhot = (float)Phydro[index]->Mhot;
    TempPhydroIOLean.MultiphaseEndTime = (float)Phydro[index]->MultiphaseEndTime;
    TempPhydroIOLean.MultiphaseFlag = Phydro[index]->MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TempPhydroIOLean.SplitFlag = Phydro[index]->SplitFlag;
    TempPhydroIOLean.SplitGeneration = Phydro[index]->SplitGeneration;
    TempPhydroIOLean.SplitNthChildren = Phydro[index]->SplitNthChildren;
    TempPhydroIOLean.AncestorGlobalID = Phydro[index]->AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TempPhydroIOLean.Pot = (float)Phydro[index]->Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TempPhydroIOLean.Tag = Phydro[index]->Tag;
#endif

	return TempPhydroIOLean;

}

struct StructPhydroIOLean CopyPhydroToTemporalStructureLeanElement(StructPhydroptr const Ph){

	struct StructPhydroIOLean TempPhydroIOLean;

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
// omit unsigned int Nlist;      // <LEAN>
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


// omit int Leaf; // <LEAN> The ID to the HydroTree leaves.

// omit bool      GravityKickFlag; // <LEAN> The kick flag of gravitataional acceleration.
// omit int       k_hydro;         // <LEAN>
// omit double    dt_hydro;        // <LEAN>
// omit double    EraLocal_hydro;  // <LEAN> The local time in a current block.

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
// omit int       k_hydro_localmin;  // <LEAN>
// omit int       k_hydro_localmin_old;  // <LEAN>
// omit double    NextUpdateEra;     // <LEAN> The next update timing in Era. This variable is used when the HydroTimeStepLimiterFlag is true.
// omit double    dt_hydro_localmin; // <LEAN>


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TempPhydroIOLean.Rho = Ph->Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TempPhydroIOLean.EnergyDensity = Ph->EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TempPhydroIOLean.Kernel = Ph->Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

// omit double    Gradh; // <LEAN> Grad-h term

// omit double    NumberDensity; // <LEAN> Number density
// omit double    GradN;         // <LEAN> Grad-N term
// omit double    fij;           // <LEAN> Part of fij. fij = 1-fij/(Uj)


// omit double    Div;        // <LEAN>
// omit double    Rot[3];     // <LEAN>
// omit double    F;          // <LEAN>
    TempPhydroIOLean.HydroAcc[0] = Ph->HydroAcc[0];
    TempPhydroIOLean.HydroAcc[1] = Ph->HydroAcc[1];
    TempPhydroIOLean.HydroAcc[2] = Ph->HydroAcc[2];
    TempPhydroIOLean.U = Ph->U;
// omit double    UPred;      // <TMP> <LEAN>
    TempPhydroIOLean.Du = Ph->Du;
// omit double    DuPrev;     // <LEAN> Du/Dt.
// omit double    DuCooling;  // <LEAN> The amount of the energy loss by the radiative cooling.


    TempPhydroIOLean.PseudoDensity = Ph->PseudoDensity;
    TempPhydroIOLean.PseudoDensityPred = Ph->PseudoDensityPred;
    TempPhydroIOLean.Zw = Ph->Zw;
    TempPhydroIOLean.ZwPred = Ph->ZwPred;
    TempPhydroIOLean.DZw = Ph->DZw;
    TempPhydroIOLean.Ddif = Ph->Ddif;


// omit double    DQheat;     // <LEAN> The amount of the heating energy by SNe, which injects during dt.
// omit double    dMass;      // <LEAN> The additional mass by SNe.
    TempPhydroIOLean.Z = Ph->Z;
    TempPhydroIOLean.ZII = Ph->ZII;
    TempPhydroIOLean.ZIa = Ph->ZIa;
    TempPhydroIOLean.dZII = Ph->dZII;
    TempPhydroIOLean.dZIa = Ph->dZIa;
// omit double    Vsig;       // <LEAN> Signal Velocity  or Maximum Velocity difference.
    TempPhydroIOLean.Alpha = Ph->Alpha;
// omit double    DAlpha;     // <LEAN>


    TempPhydroIOLean.DZdiff = Ph->DZdiff;
    TempPhydroIOLean.Sxy = Ph->Sxy;
    TempPhydroIOLean.Sxz = Ph->Sxz;
    TempPhydroIOLean.Syx = Ph->Syx;
    TempPhydroIOLean.Syz = Ph->Syz;
    TempPhydroIOLean.Szx = Ph->Szx;
    TempPhydroIOLean.Szy = Ph->Szy;

#if VISCOSITY_TYPE == 1 //{
    TempPhydroIOLean.Bxx = Ph->Bxx;
    TempPhydroIOLean.Bxy = Ph->Bxy;
    TempPhydroIOLean.Bxz = Ph->Bxz;
    TempPhydroIOLean.Byx = Ph->Byx;
    TempPhydroIOLean.Byy = Ph->Byy;
    TempPhydroIOLean.Byz = Ph->Byz;
    TempPhydroIOLean.Bzx = Ph->Bzx;
    TempPhydroIOLean.Bzy = Ph->Bzy;
    TempPhydroIOLean.Bzz = Ph->Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TempPhydroIOLean.G0 = Ph->G0;
    TempPhydroIOLean.fH2 = Ph->fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
// omit short   SpawnTimes; // <LEAN>
// omit double  SpawnMass;  // <LEAN>
#endif
#ifdef USE_CELIB
    TempPhydroIOLean.Elements[0] = Ph->Elements[0];
    TempPhydroIOLean.Elements[1] = Ph->Elements[1];
    TempPhydroIOLean.Elements[2] = Ph->Elements[2];
    TempPhydroIOLean.Elements[3] = Ph->Elements[3];
    TempPhydroIOLean.Elements[4] = Ph->Elements[4];
    TempPhydroIOLean.Elements[5] = Ph->Elements[5];
    TempPhydroIOLean.Elements[6] = Ph->Elements[6];
    TempPhydroIOLean.Elements[7] = Ph->Elements[7];
    TempPhydroIOLean.Elements[8] = Ph->Elements[8];
    TempPhydroIOLean.Elements[9] = Ph->Elements[9];
    TempPhydroIOLean.Elements[10] = Ph->Elements[10];
    TempPhydroIOLean.Elements[11] = Ph->Elements[11];
    TempPhydroIOLean.Elements[12] = Ph->Elements[12];
    TempPhydroIOLean.Elements[13] = Ph->Elements[13];
    TempPhydroIOLean.Elements[14] = Ph->Elements[14];
    TempPhydroIOLean.Elements[15] = Ph->Elements[15];
    TempPhydroIOLean.Elements[16] = Ph->Elements[16];
    TempPhydroIOLean.Elements[17] = Ph->Elements[17];
#endif // USE_CELIB

// omit double GradRho;               // <LEAN>
// omit double G0thin;                // <LEAN>
// omit double G0thinLocal;           // <LEAN>
// omit double G0extLocal;            // <LEAN>
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TempPhydroIOLean.MomentumRP[0] = Ph->MomentumRP[0];
    TempPhydroIOLean.MomentumRP[1] = Ph->MomentumRP[1];
    TempPhydroIOLean.MomentumRP[2] = Ph->MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPhydroIOLean.MomentumFB[0] = Ph->MomentumFB[0];
    TempPhydroIOLean.MomentumFB[1] = Ph->MomentumFB[1];
    TempPhydroIOLean.MomentumFB[2] = Ph->MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TempPhydroIOLean.Uhot = Ph->Uhot;
    TempPhydroIOLean.Rhohot = Ph->Rhohot;
    TempPhydroIOLean.Mhot = Ph->Mhot;
    TempPhydroIOLean.MultiphaseEndTime = Ph->MultiphaseEndTime;
    TempPhydroIOLean.MultiphaseFlag = Ph->MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TempPhydroIOLean.SplitFlag = Ph->SplitFlag;
    TempPhydroIOLean.SplitGeneration = Ph->SplitGeneration;
    TempPhydroIOLean.SplitNthChildren = Ph->SplitNthChildren;
    TempPhydroIOLean.AncestorGlobalID = Ph->AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TempPhydroIOLean.Pot = Ph->Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TempPhydroIOLean.Tag = Ph->Tag;
#endif

	return TempPhydroIOLean;

}

StructPhydro CopyTemporalStructureLeanToPhydroLean(struct StructPhydroIOLean PhydroIOLean){

	StructPhydro TemporalPhydro;
	memset(&TemporalPhydro,0,sizeof(StructPhydro));

// omit StructPhydroptr Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>


// omit bool Use;                // <TMP> <LEAN>
// omit unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
// omit unsigned int Nlist;      // <LEAN>
// omit unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).

// omit bool    Active;          // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.
// omit bool    CoolingFlag;     // <TMP> <LEAN> if \"OFF\", the particle skips the cooling routine.


// omit int Leaf; // <LEAN> The ID to the HydroTree leaves.

// omit bool      GravityKickFlag; // <LEAN> The kick flag of gravitataional acceleration.
// omit int       k_hydro;         // <LEAN>
// omit double    dt_hydro;        // <LEAN>
// omit double    EraLocal_hydro;  // <LEAN> The local time in a current block.

// omit bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
// omit int       k_hydro_localmin;  // <LEAN>
// omit int       k_hydro_localmin_old;  // <LEAN>
// omit double    NextUpdateEra;     // <LEAN> The next update timing in Era. This variable is used when the HydroTimeStepLimiterFlag is true.
// omit double    dt_hydro_localmin; // <LEAN>


// omit double    Mass;         // <TMP> <LEAN>
// omit double    PosP[3];      // <TMP> <LEAN>
// omit double    VelP[3];      // <TMP> <LEAN>

    TemporalPhydro.Rho = PhydroIOLean.Rho;
// omit double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    TemporalPhydro.EnergyDensity = PhydroIOLean.EnergyDensity;
// omit double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
// omit double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    TemporalPhydro.Kernel = PhydroIOLean.Kernel;
// omit double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2

// omit double    Gradh; // <LEAN> Grad-h term

// omit double    NumberDensity; // <LEAN> Number density
// omit double    GradN;         // <LEAN> Grad-N term
// omit double    fij;           // <LEAN> Part of fij. fij = 1-fij/(Uj)


// omit double    Div;        // <LEAN>
// omit double    Rot[3];     // <LEAN>
// omit double    F;          // <LEAN>
    TemporalPhydro.HydroAcc[0] = PhydroIOLean.HydroAcc[0];
    TemporalPhydro.HydroAcc[1] = PhydroIOLean.HydroAcc[1];
    TemporalPhydro.HydroAcc[2] = PhydroIOLean.HydroAcc[2];
    TemporalPhydro.U = PhydroIOLean.U;
// omit double    UPred;      // <TMP> <LEAN>
    TemporalPhydro.Du = PhydroIOLean.Du;
// omit double    DuPrev;     // <LEAN> Du/Dt.
// omit double    DuCooling;  // <LEAN> The amount of the energy loss by the radiative cooling.


    TemporalPhydro.PseudoDensity = PhydroIOLean.PseudoDensity;
    TemporalPhydro.PseudoDensityPred = PhydroIOLean.PseudoDensityPred;
    TemporalPhydro.Zw = PhydroIOLean.Zw;
    TemporalPhydro.ZwPred = PhydroIOLean.ZwPred;
    TemporalPhydro.DZw = PhydroIOLean.DZw;
    TemporalPhydro.Ddif = PhydroIOLean.Ddif;


// omit double    DQheat;     // <LEAN> The amount of the heating energy by SNe, which injects during dt.
// omit double    dMass;      // <LEAN> The additional mass by SNe.
    TemporalPhydro.Z = PhydroIOLean.Z;
    TemporalPhydro.ZII = PhydroIOLean.ZII;
    TemporalPhydro.ZIa = PhydroIOLean.ZIa;
    TemporalPhydro.dZII = PhydroIOLean.dZII;
    TemporalPhydro.dZIa = PhydroIOLean.dZIa;
// omit double    Vsig;       // <LEAN> Signal Velocity  or Maximum Velocity difference.
    TemporalPhydro.Alpha = PhydroIOLean.Alpha;
// omit double    DAlpha;     // <LEAN>


    TemporalPhydro.DZdiff = PhydroIOLean.DZdiff;
    TemporalPhydro.Sxy = PhydroIOLean.Sxy;
    TemporalPhydro.Sxz = PhydroIOLean.Sxz;
    TemporalPhydro.Syx = PhydroIOLean.Syx;
    TemporalPhydro.Syz = PhydroIOLean.Syz;
    TemporalPhydro.Szx = PhydroIOLean.Szx;
    TemporalPhydro.Szy = PhydroIOLean.Szy;

#if VISCOSITY_TYPE == 1 //{
    TemporalPhydro.Bxx = PhydroIOLean.Bxx;
    TemporalPhydro.Bxy = PhydroIOLean.Bxy;
    TemporalPhydro.Bxz = PhydroIOLean.Bxz;
    TemporalPhydro.Byx = PhydroIOLean.Byx;
    TemporalPhydro.Byy = PhydroIOLean.Byy;
    TemporalPhydro.Byz = PhydroIOLean.Byz;
    TemporalPhydro.Bzx = PhydroIOLean.Bzx;
    TemporalPhydro.Bzy = PhydroIOLean.Bzy;
    TemporalPhydro.Bzz = PhydroIOLean.Bzz;
#endif // VISCOSITY_TYPE == 1 //}




#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    TemporalPhydro.G0 = PhydroIOLean.G0;
    TemporalPhydro.fH2 = PhydroIOLean.fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
// omit short   SpawnTimes; // <LEAN>
// omit double  SpawnMass;  // <LEAN>
#endif
#ifdef USE_CELIB
    TemporalPhydro.Elements[0] = PhydroIOLean.Elements[0];
    TemporalPhydro.Elements[1] = PhydroIOLean.Elements[1];
    TemporalPhydro.Elements[2] = PhydroIOLean.Elements[2];
    TemporalPhydro.Elements[3] = PhydroIOLean.Elements[3];
    TemporalPhydro.Elements[4] = PhydroIOLean.Elements[4];
    TemporalPhydro.Elements[5] = PhydroIOLean.Elements[5];
    TemporalPhydro.Elements[6] = PhydroIOLean.Elements[6];
    TemporalPhydro.Elements[7] = PhydroIOLean.Elements[7];
    TemporalPhydro.Elements[8] = PhydroIOLean.Elements[8];
    TemporalPhydro.Elements[9] = PhydroIOLean.Elements[9];
    TemporalPhydro.Elements[10] = PhydroIOLean.Elements[10];
    TemporalPhydro.Elements[11] = PhydroIOLean.Elements[11];
    TemporalPhydro.Elements[12] = PhydroIOLean.Elements[12];
    TemporalPhydro.Elements[13] = PhydroIOLean.Elements[13];
    TemporalPhydro.Elements[14] = PhydroIOLean.Elements[14];
    TemporalPhydro.Elements[15] = PhydroIOLean.Elements[15];
    TemporalPhydro.Elements[16] = PhydroIOLean.Elements[16];
    TemporalPhydro.Elements[17] = PhydroIOLean.Elements[17];
#endif // USE_CELIB

// omit double GradRho;               // <LEAN>
// omit double G0thin;                // <LEAN>
// omit double G0thinLocal;           // <LEAN>
// omit double G0extLocal;            // <LEAN>
// omit double G0thinNextUpdateTime;  // <TMP> <LEAN>


#ifdef USE_RADIATION_PRESSURE //{
    TemporalPhydro.MomentumRP[0] = PhydroIOLean.MomentumRP[0];
    TemporalPhydro.MomentumRP[1] = PhydroIOLean.MomentumRP[1];
    TemporalPhydro.MomentumRP[2] = PhydroIOLean.MomentumRP[2];
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    TemporalPhydro.MomentumFB[0] = PhydroIOLean.MomentumFB[0];
    TemporalPhydro.MomentumFB[1] = PhydroIOLean.MomentumFB[1];
    TemporalPhydro.MomentumFB[2] = PhydroIOLean.MomentumFB[2];
#endif //USE_MOMENTUM_FEEDBACK //{


    TemporalPhydro.Uhot = PhydroIOLean.Uhot;
    TemporalPhydro.Rhohot = PhydroIOLean.Rhohot;
    TemporalPhydro.Mhot = PhydroIOLean.Mhot;
    TemporalPhydro.MultiphaseEndTime = PhydroIOLean.MultiphaseEndTime;
    TemporalPhydro.MultiphaseFlag = PhydroIOLean.MultiphaseFlag;


#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPhydro.SplitFlag = PhydroIOLean.SplitFlag;
    TemporalPhydro.SplitGeneration = PhydroIOLean.SplitGeneration;
    TemporalPhydro.SplitNthChildren = PhydroIOLean.SplitNthChildren;
    TemporalPhydro.AncestorGlobalID = PhydroIOLean.AncestorGlobalID;
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    TemporalPhydro.Pot = PhydroIOLean.Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    TemporalPhydro.Tag = PhydroIOLean.Tag;
#endif

	return TemporalPhydro;

}

struct StructPstarIOCompact CopyPstarToTemporalStructureCompact(const int index){

	struct StructPstarIOCompact TempPstarIOCompact;

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPstarIOCompact.TypeII = Pstar[index]->TypeII;
    TempPstarIOCompact.TypeIa = Pstar[index]->TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TempPstarIOCompact.TypeIIProb = Pstar[index]->TypeIIProb;
#endif

    TempPstarIOCompact.IMFTYPE = Pstar[index]->IMFTYPE;
    TempPstarIOCompact.NthChildren = Pstar[index]->NthChildren;
    TempPstarIOCompact.ParentGlobalID = Pstar[index]->ParentGlobalID;

    TempPstarIOCompact.InitialMass = (float)Pstar[index]->InitialMass;
    TempPstarIOCompact.Mass = (float)Pstar[index]->Mass;
    TempPstarIOCompact.FormationTime = (float)Pstar[index]->FormationTime;
    TempPstarIOCompact.Z = (float)Pstar[index]->Z;
    TempPstarIOCompact.ZII = (float)Pstar[index]->ZII;
    TempPstarIOCompact.ZIa = (float)Pstar[index]->ZIa;
    TempPstarIOCompact.TempForm = (float)Pstar[index]->TempForm;
    TempPstarIOCompact.RhoForm = (float)Pstar[index]->RhoForm;
    TempPstarIOCompact.StromgrenRadius = (float)Pstar[index]->StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TempPstarIOCompact.SNIICount = Pstar[index]->SNIICount;
    TempPstarIOCompact.SNIaCount = Pstar[index]->SNIaCount;
    TempPstarIOCompact.EventTimeSNII = (float)Pstar[index]->EventTimeSNII;
    TempPstarIOCompact.EventTimeSNIa = (float)Pstar[index]->EventTimeSNIa;

    TempPstarIOCompact.AGBCount = Pstar[index]->AGBCount;
    TempPstarIOCompact.EventTimeAGB = (float)Pstar[index]->EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TempPstarIOCompact.NSMCount = Pstar[index]->NSMCount;
    TempPstarIOCompact.EventTimeNSM = (float)Pstar[index]->EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TempPstarIOCompact.Elements[0] = (float)Pstar[index]->Elements[0];
    TempPstarIOCompact.Elements[1] = (float)Pstar[index]->Elements[1];
    TempPstarIOCompact.Elements[2] = (float)Pstar[index]->Elements[2];
    TempPstarIOCompact.Elements[3] = (float)Pstar[index]->Elements[3];
    TempPstarIOCompact.Elements[4] = (float)Pstar[index]->Elements[4];
    TempPstarIOCompact.Elements[5] = (float)Pstar[index]->Elements[5];
    TempPstarIOCompact.Elements[6] = (float)Pstar[index]->Elements[6];
    TempPstarIOCompact.Elements[7] = (float)Pstar[index]->Elements[7];
    TempPstarIOCompact.Elements[8] = (float)Pstar[index]->Elements[8];
    TempPstarIOCompact.Elements[9] = (float)Pstar[index]->Elements[9];
    TempPstarIOCompact.Elements[10] = (float)Pstar[index]->Elements[10];
    TempPstarIOCompact.Elements[11] = (float)Pstar[index]->Elements[11];
    TempPstarIOCompact.Elements[12] = (float)Pstar[index]->Elements[12];
    TempPstarIOCompact.Elements[13] = (float)Pstar[index]->Elements[13];
    TempPstarIOCompact.Elements[14] = (float)Pstar[index]->Elements[14];
    TempPstarIOCompact.Elements[15] = (float)Pstar[index]->Elements[15];
    TempPstarIOCompact.Elements[16] = (float)Pstar[index]->Elements[16];
    TempPstarIOCompact.Elements[17] = (float)Pstar[index]->Elements[17];
#endif // USE_CELIB //}


    TempPstarIOCompact.LFUV = (float)Pstar[index]->LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TempPstarIOCompact.BolometricLuminosity = (float)Pstar[index]->BolometricLuminosity;
    TempPstarIOCompact.RadiusRP = (float)Pstar[index]->RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TempPstarIOCompact.StellarWindEnergy = (float)Pstar[index]->StellarWindEnergy;
    TempPstarIOCompact.RadiusSW = (float)Pstar[index]->RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPstarIOCompact.nH_ave = (float)Pstar[index]->nH_ave;
    TempPstarIOCompact.Z_ave = (float)Pstar[index]->Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TempPstarIOCompact.SplitGeneration = Pstar[index]->SplitGeneration;
    TempPstarIOCompact.SplitNthChildren = Pstar[index]->SplitNthChildren;
    TempPstarIOCompact.GrandParentGlobalID = Pstar[index]->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPstarIOCompact;

}

struct StructPstarIOCompact CopyPstarToTemporalStructureCompactElement(StructPstarptr const Ps){

	struct StructPstarIOCompact TempPstarIOCompact;

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPstarIOCompact.TypeII = Ps->TypeII;
    TempPstarIOCompact.TypeIa = Ps->TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TempPstarIOCompact.TypeIIProb = Ps->TypeIIProb;
#endif

    TempPstarIOCompact.IMFTYPE = Ps->IMFTYPE;
    TempPstarIOCompact.NthChildren = Ps->NthChildren;
    TempPstarIOCompact.ParentGlobalID = Ps->ParentGlobalID;

    TempPstarIOCompact.InitialMass = Ps->InitialMass;
    TempPstarIOCompact.Mass = Ps->Mass;
    TempPstarIOCompact.FormationTime = Ps->FormationTime;
    TempPstarIOCompact.Z = Ps->Z;
    TempPstarIOCompact.ZII = Ps->ZII;
    TempPstarIOCompact.ZIa = Ps->ZIa;
    TempPstarIOCompact.TempForm = Ps->TempForm;
    TempPstarIOCompact.RhoForm = Ps->RhoForm;
    TempPstarIOCompact.StromgrenRadius = Ps->StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TempPstarIOCompact.SNIICount = Ps->SNIICount;
    TempPstarIOCompact.SNIaCount = Ps->SNIaCount;
    TempPstarIOCompact.EventTimeSNII = Ps->EventTimeSNII;
    TempPstarIOCompact.EventTimeSNIa = Ps->EventTimeSNIa;

    TempPstarIOCompact.AGBCount = Ps->AGBCount;
    TempPstarIOCompact.EventTimeAGB = Ps->EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TempPstarIOCompact.NSMCount = Ps->NSMCount;
    TempPstarIOCompact.EventTimeNSM = Ps->EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TempPstarIOCompact.Elements[0] = Ps->Elements[0];
    TempPstarIOCompact.Elements[1] = Ps->Elements[1];
    TempPstarIOCompact.Elements[2] = Ps->Elements[2];
    TempPstarIOCompact.Elements[3] = Ps->Elements[3];
    TempPstarIOCompact.Elements[4] = Ps->Elements[4];
    TempPstarIOCompact.Elements[5] = Ps->Elements[5];
    TempPstarIOCompact.Elements[6] = Ps->Elements[6];
    TempPstarIOCompact.Elements[7] = Ps->Elements[7];
    TempPstarIOCompact.Elements[8] = Ps->Elements[8];
    TempPstarIOCompact.Elements[9] = Ps->Elements[9];
    TempPstarIOCompact.Elements[10] = Ps->Elements[10];
    TempPstarIOCompact.Elements[11] = Ps->Elements[11];
    TempPstarIOCompact.Elements[12] = Ps->Elements[12];
    TempPstarIOCompact.Elements[13] = Ps->Elements[13];
    TempPstarIOCompact.Elements[14] = Ps->Elements[14];
    TempPstarIOCompact.Elements[15] = Ps->Elements[15];
    TempPstarIOCompact.Elements[16] = Ps->Elements[16];
    TempPstarIOCompact.Elements[17] = Ps->Elements[17];
#endif // USE_CELIB //}


    TempPstarIOCompact.LFUV = Ps->LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TempPstarIOCompact.BolometricLuminosity = Ps->BolometricLuminosity;
    TempPstarIOCompact.RadiusRP = Ps->RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TempPstarIOCompact.StellarWindEnergy = Ps->StellarWindEnergy;
    TempPstarIOCompact.RadiusSW = Ps->RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPstarIOCompact.nH_ave = Ps->nH_ave;
    TempPstarIOCompact.Z_ave = Ps->Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TempPstarIOCompact.SplitGeneration = Ps->SplitGeneration;
    TempPstarIOCompact.SplitNthChildren = Ps->SplitNthChildren;
    TempPstarIOCompact.GrandParentGlobalID = Ps->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPstarIOCompact;

}

StructPstar CopyTemporalStructureCompactToPstarCompact(struct StructPstarIOCompact PstarIOCompact){

	StructPstar TemporalPstar;
	memset(&TemporalPstar,0,sizeof(StructPstar));

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TemporalPstar.TypeII = PstarIOCompact.TypeII;
    TemporalPstar.TypeIa = PstarIOCompact.TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TemporalPstar.TypeIIProb = PstarIOCompact.TypeIIProb;
#endif

    TemporalPstar.IMFTYPE = PstarIOCompact.IMFTYPE;
    TemporalPstar.NthChildren = PstarIOCompact.NthChildren;
    TemporalPstar.ParentGlobalID = PstarIOCompact.ParentGlobalID;

    TemporalPstar.InitialMass = PstarIOCompact.InitialMass;
    TemporalPstar.Mass = PstarIOCompact.Mass;
    TemporalPstar.FormationTime = PstarIOCompact.FormationTime;
    TemporalPstar.Z = PstarIOCompact.Z;
    TemporalPstar.ZII = PstarIOCompact.ZII;
    TemporalPstar.ZIa = PstarIOCompact.ZIa;
    TemporalPstar.TempForm = PstarIOCompact.TempForm;
    TemporalPstar.RhoForm = PstarIOCompact.RhoForm;
    TemporalPstar.StromgrenRadius = PstarIOCompact.StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TemporalPstar.SNIICount = PstarIOCompact.SNIICount;
    TemporalPstar.SNIaCount = PstarIOCompact.SNIaCount;
    TemporalPstar.EventTimeSNII = PstarIOCompact.EventTimeSNII;
    TemporalPstar.EventTimeSNIa = PstarIOCompact.EventTimeSNIa;

    TemporalPstar.AGBCount = PstarIOCompact.AGBCount;
    TemporalPstar.EventTimeAGB = PstarIOCompact.EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TemporalPstar.NSMCount = PstarIOCompact.NSMCount;
    TemporalPstar.EventTimeNSM = PstarIOCompact.EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TemporalPstar.Elements[0] = PstarIOCompact.Elements[0];
    TemporalPstar.Elements[1] = PstarIOCompact.Elements[1];
    TemporalPstar.Elements[2] = PstarIOCompact.Elements[2];
    TemporalPstar.Elements[3] = PstarIOCompact.Elements[3];
    TemporalPstar.Elements[4] = PstarIOCompact.Elements[4];
    TemporalPstar.Elements[5] = PstarIOCompact.Elements[5];
    TemporalPstar.Elements[6] = PstarIOCompact.Elements[6];
    TemporalPstar.Elements[7] = PstarIOCompact.Elements[7];
    TemporalPstar.Elements[8] = PstarIOCompact.Elements[8];
    TemporalPstar.Elements[9] = PstarIOCompact.Elements[9];
    TemporalPstar.Elements[10] = PstarIOCompact.Elements[10];
    TemporalPstar.Elements[11] = PstarIOCompact.Elements[11];
    TemporalPstar.Elements[12] = PstarIOCompact.Elements[12];
    TemporalPstar.Elements[13] = PstarIOCompact.Elements[13];
    TemporalPstar.Elements[14] = PstarIOCompact.Elements[14];
    TemporalPstar.Elements[15] = PstarIOCompact.Elements[15];
    TemporalPstar.Elements[16] = PstarIOCompact.Elements[16];
    TemporalPstar.Elements[17] = PstarIOCompact.Elements[17];
#endif // USE_CELIB //}


    TemporalPstar.LFUV = PstarIOCompact.LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TemporalPstar.BolometricLuminosity = PstarIOCompact.BolometricLuminosity;
    TemporalPstar.RadiusRP = PstarIOCompact.RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TemporalPstar.StellarWindEnergy = PstarIOCompact.StellarWindEnergy;
    TemporalPstar.RadiusSW = PstarIOCompact.RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TemporalPstar.nH_ave = PstarIOCompact.nH_ave;
    TemporalPstar.Z_ave = PstarIOCompact.Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPstar.SplitGeneration = PstarIOCompact.SplitGeneration;
    TemporalPstar.SplitNthChildren = PstarIOCompact.SplitNthChildren;
    TemporalPstar.GrandParentGlobalID = PstarIOCompact.GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TemporalPstar;

}

struct StructPstarIOCompactDouble CopyPstarToTemporalStructureCompactDouble(const int index){

	struct StructPstarIOCompactDouble TempPstarIOCompactDouble;

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPstarIOCompactDouble.TypeII = Pstar[index]->TypeII;
    TempPstarIOCompactDouble.TypeIa = Pstar[index]->TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TempPstarIOCompactDouble.TypeIIProb = Pstar[index]->TypeIIProb;
#endif

    TempPstarIOCompactDouble.IMFTYPE = Pstar[index]->IMFTYPE;
    TempPstarIOCompactDouble.NthChildren = Pstar[index]->NthChildren;
    TempPstarIOCompactDouble.ParentGlobalID = Pstar[index]->ParentGlobalID;

    TempPstarIOCompactDouble.InitialMass = (double)Pstar[index]->InitialMass;
    TempPstarIOCompactDouble.Mass = (double)Pstar[index]->Mass;
    TempPstarIOCompactDouble.FormationTime = (double)Pstar[index]->FormationTime;
    TempPstarIOCompactDouble.Z = (double)Pstar[index]->Z;
    TempPstarIOCompactDouble.ZII = (double)Pstar[index]->ZII;
    TempPstarIOCompactDouble.ZIa = (double)Pstar[index]->ZIa;
    TempPstarIOCompactDouble.TempForm = (double)Pstar[index]->TempForm;
    TempPstarIOCompactDouble.RhoForm = (double)Pstar[index]->RhoForm;
    TempPstarIOCompactDouble.StromgrenRadius = (double)Pstar[index]->StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TempPstarIOCompactDouble.SNIICount = Pstar[index]->SNIICount;
    TempPstarIOCompactDouble.SNIaCount = Pstar[index]->SNIaCount;
    TempPstarIOCompactDouble.EventTimeSNII = (double)Pstar[index]->EventTimeSNII;
    TempPstarIOCompactDouble.EventTimeSNIa = (double)Pstar[index]->EventTimeSNIa;

    TempPstarIOCompactDouble.AGBCount = Pstar[index]->AGBCount;
    TempPstarIOCompactDouble.EventTimeAGB = (double)Pstar[index]->EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TempPstarIOCompactDouble.NSMCount = Pstar[index]->NSMCount;
    TempPstarIOCompactDouble.EventTimeNSM = (double)Pstar[index]->EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TempPstarIOCompactDouble.Elements[0] = (double)Pstar[index]->Elements[0];
    TempPstarIOCompactDouble.Elements[1] = (double)Pstar[index]->Elements[1];
    TempPstarIOCompactDouble.Elements[2] = (double)Pstar[index]->Elements[2];
    TempPstarIOCompactDouble.Elements[3] = (double)Pstar[index]->Elements[3];
    TempPstarIOCompactDouble.Elements[4] = (double)Pstar[index]->Elements[4];
    TempPstarIOCompactDouble.Elements[5] = (double)Pstar[index]->Elements[5];
    TempPstarIOCompactDouble.Elements[6] = (double)Pstar[index]->Elements[6];
    TempPstarIOCompactDouble.Elements[7] = (double)Pstar[index]->Elements[7];
    TempPstarIOCompactDouble.Elements[8] = (double)Pstar[index]->Elements[8];
    TempPstarIOCompactDouble.Elements[9] = (double)Pstar[index]->Elements[9];
    TempPstarIOCompactDouble.Elements[10] = (double)Pstar[index]->Elements[10];
    TempPstarIOCompactDouble.Elements[11] = (double)Pstar[index]->Elements[11];
    TempPstarIOCompactDouble.Elements[12] = (double)Pstar[index]->Elements[12];
    TempPstarIOCompactDouble.Elements[13] = (double)Pstar[index]->Elements[13];
    TempPstarIOCompactDouble.Elements[14] = (double)Pstar[index]->Elements[14];
    TempPstarIOCompactDouble.Elements[15] = (double)Pstar[index]->Elements[15];
    TempPstarIOCompactDouble.Elements[16] = (double)Pstar[index]->Elements[16];
    TempPstarIOCompactDouble.Elements[17] = (double)Pstar[index]->Elements[17];
#endif // USE_CELIB //}


    TempPstarIOCompactDouble.LFUV = (double)Pstar[index]->LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TempPstarIOCompactDouble.BolometricLuminosity = (double)Pstar[index]->BolometricLuminosity;
    TempPstarIOCompactDouble.RadiusRP = (double)Pstar[index]->RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TempPstarIOCompactDouble.StellarWindEnergy = (double)Pstar[index]->StellarWindEnergy;
    TempPstarIOCompactDouble.RadiusSW = (double)Pstar[index]->RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPstarIOCompactDouble.nH_ave = (double)Pstar[index]->nH_ave;
    TempPstarIOCompactDouble.Z_ave = (double)Pstar[index]->Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TempPstarIOCompactDouble.SplitGeneration = Pstar[index]->SplitGeneration;
    TempPstarIOCompactDouble.SplitNthChildren = Pstar[index]->SplitNthChildren;
    TempPstarIOCompactDouble.GrandParentGlobalID = Pstar[index]->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPstarIOCompactDouble;

}

struct StructPstarIOCompactDouble CopyPstarToTemporalStructureCompactDoubleElement(StructPstarptr const Ps){

	struct StructPstarIOCompactDouble TempPstarIOCompactDouble;

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPstarIOCompactDouble.TypeII = Ps->TypeII;
    TempPstarIOCompactDouble.TypeIa = Ps->TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TempPstarIOCompactDouble.TypeIIProb = Ps->TypeIIProb;
#endif

    TempPstarIOCompactDouble.IMFTYPE = Ps->IMFTYPE;
    TempPstarIOCompactDouble.NthChildren = Ps->NthChildren;
    TempPstarIOCompactDouble.ParentGlobalID = Ps->ParentGlobalID;

    TempPstarIOCompactDouble.InitialMass = Ps->InitialMass;
    TempPstarIOCompactDouble.Mass = Ps->Mass;
    TempPstarIOCompactDouble.FormationTime = Ps->FormationTime;
    TempPstarIOCompactDouble.Z = Ps->Z;
    TempPstarIOCompactDouble.ZII = Ps->ZII;
    TempPstarIOCompactDouble.ZIa = Ps->ZIa;
    TempPstarIOCompactDouble.TempForm = Ps->TempForm;
    TempPstarIOCompactDouble.RhoForm = Ps->RhoForm;
    TempPstarIOCompactDouble.StromgrenRadius = Ps->StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TempPstarIOCompactDouble.SNIICount = Ps->SNIICount;
    TempPstarIOCompactDouble.SNIaCount = Ps->SNIaCount;
    TempPstarIOCompactDouble.EventTimeSNII = Ps->EventTimeSNII;
    TempPstarIOCompactDouble.EventTimeSNIa = Ps->EventTimeSNIa;

    TempPstarIOCompactDouble.AGBCount = Ps->AGBCount;
    TempPstarIOCompactDouble.EventTimeAGB = Ps->EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TempPstarIOCompactDouble.NSMCount = Ps->NSMCount;
    TempPstarIOCompactDouble.EventTimeNSM = Ps->EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TempPstarIOCompactDouble.Elements[0] = Ps->Elements[0];
    TempPstarIOCompactDouble.Elements[1] = Ps->Elements[1];
    TempPstarIOCompactDouble.Elements[2] = Ps->Elements[2];
    TempPstarIOCompactDouble.Elements[3] = Ps->Elements[3];
    TempPstarIOCompactDouble.Elements[4] = Ps->Elements[4];
    TempPstarIOCompactDouble.Elements[5] = Ps->Elements[5];
    TempPstarIOCompactDouble.Elements[6] = Ps->Elements[6];
    TempPstarIOCompactDouble.Elements[7] = Ps->Elements[7];
    TempPstarIOCompactDouble.Elements[8] = Ps->Elements[8];
    TempPstarIOCompactDouble.Elements[9] = Ps->Elements[9];
    TempPstarIOCompactDouble.Elements[10] = Ps->Elements[10];
    TempPstarIOCompactDouble.Elements[11] = Ps->Elements[11];
    TempPstarIOCompactDouble.Elements[12] = Ps->Elements[12];
    TempPstarIOCompactDouble.Elements[13] = Ps->Elements[13];
    TempPstarIOCompactDouble.Elements[14] = Ps->Elements[14];
    TempPstarIOCompactDouble.Elements[15] = Ps->Elements[15];
    TempPstarIOCompactDouble.Elements[16] = Ps->Elements[16];
    TempPstarIOCompactDouble.Elements[17] = Ps->Elements[17];
#endif // USE_CELIB //}


    TempPstarIOCompactDouble.LFUV = Ps->LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TempPstarIOCompactDouble.BolometricLuminosity = Ps->BolometricLuminosity;
    TempPstarIOCompactDouble.RadiusRP = Ps->RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TempPstarIOCompactDouble.StellarWindEnergy = Ps->StellarWindEnergy;
    TempPstarIOCompactDouble.RadiusSW = Ps->RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPstarIOCompactDouble.nH_ave = Ps->nH_ave;
    TempPstarIOCompactDouble.Z_ave = Ps->Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TempPstarIOCompactDouble.SplitGeneration = Ps->SplitGeneration;
    TempPstarIOCompactDouble.SplitNthChildren = Ps->SplitNthChildren;
    TempPstarIOCompactDouble.GrandParentGlobalID = Ps->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPstarIOCompactDouble;

}

StructPstar CopyTemporalStructureCompactDoubleToPstarCompactDouble(struct StructPstarIOCompactDouble PstarIOCompactDouble){

	StructPstar TemporalPstar;
	memset(&TemporalPstar,0,sizeof(StructPstar));

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TemporalPstar.TypeII = PstarIOCompactDouble.TypeII;
    TemporalPstar.TypeIa = PstarIOCompactDouble.TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TemporalPstar.TypeIIProb = PstarIOCompactDouble.TypeIIProb;
#endif

    TemporalPstar.IMFTYPE = PstarIOCompactDouble.IMFTYPE;
    TemporalPstar.NthChildren = PstarIOCompactDouble.NthChildren;
    TemporalPstar.ParentGlobalID = PstarIOCompactDouble.ParentGlobalID;

    TemporalPstar.InitialMass = PstarIOCompactDouble.InitialMass;
    TemporalPstar.Mass = PstarIOCompactDouble.Mass;
    TemporalPstar.FormationTime = PstarIOCompactDouble.FormationTime;
    TemporalPstar.Z = PstarIOCompactDouble.Z;
    TemporalPstar.ZII = PstarIOCompactDouble.ZII;
    TemporalPstar.ZIa = PstarIOCompactDouble.ZIa;
    TemporalPstar.TempForm = PstarIOCompactDouble.TempForm;
    TemporalPstar.RhoForm = PstarIOCompactDouble.RhoForm;
    TemporalPstar.StromgrenRadius = PstarIOCompactDouble.StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TemporalPstar.SNIICount = PstarIOCompactDouble.SNIICount;
    TemporalPstar.SNIaCount = PstarIOCompactDouble.SNIaCount;
    TemporalPstar.EventTimeSNII = PstarIOCompactDouble.EventTimeSNII;
    TemporalPstar.EventTimeSNIa = PstarIOCompactDouble.EventTimeSNIa;

    TemporalPstar.AGBCount = PstarIOCompactDouble.AGBCount;
    TemporalPstar.EventTimeAGB = PstarIOCompactDouble.EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TemporalPstar.NSMCount = PstarIOCompactDouble.NSMCount;
    TemporalPstar.EventTimeNSM = PstarIOCompactDouble.EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TemporalPstar.Elements[0] = PstarIOCompactDouble.Elements[0];
    TemporalPstar.Elements[1] = PstarIOCompactDouble.Elements[1];
    TemporalPstar.Elements[2] = PstarIOCompactDouble.Elements[2];
    TemporalPstar.Elements[3] = PstarIOCompactDouble.Elements[3];
    TemporalPstar.Elements[4] = PstarIOCompactDouble.Elements[4];
    TemporalPstar.Elements[5] = PstarIOCompactDouble.Elements[5];
    TemporalPstar.Elements[6] = PstarIOCompactDouble.Elements[6];
    TemporalPstar.Elements[7] = PstarIOCompactDouble.Elements[7];
    TemporalPstar.Elements[8] = PstarIOCompactDouble.Elements[8];
    TemporalPstar.Elements[9] = PstarIOCompactDouble.Elements[9];
    TemporalPstar.Elements[10] = PstarIOCompactDouble.Elements[10];
    TemporalPstar.Elements[11] = PstarIOCompactDouble.Elements[11];
    TemporalPstar.Elements[12] = PstarIOCompactDouble.Elements[12];
    TemporalPstar.Elements[13] = PstarIOCompactDouble.Elements[13];
    TemporalPstar.Elements[14] = PstarIOCompactDouble.Elements[14];
    TemporalPstar.Elements[15] = PstarIOCompactDouble.Elements[15];
    TemporalPstar.Elements[16] = PstarIOCompactDouble.Elements[16];
    TemporalPstar.Elements[17] = PstarIOCompactDouble.Elements[17];
#endif // USE_CELIB //}


    TemporalPstar.LFUV = PstarIOCompactDouble.LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TemporalPstar.BolometricLuminosity = PstarIOCompactDouble.BolometricLuminosity;
    TemporalPstar.RadiusRP = PstarIOCompactDouble.RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TemporalPstar.StellarWindEnergy = PstarIOCompactDouble.StellarWindEnergy;
    TemporalPstar.RadiusSW = PstarIOCompactDouble.RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TemporalPstar.nH_ave = PstarIOCompactDouble.nH_ave;
    TemporalPstar.Z_ave = PstarIOCompactDouble.Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPstar.SplitGeneration = PstarIOCompactDouble.SplitGeneration;
    TemporalPstar.SplitNthChildren = PstarIOCompactDouble.SplitNthChildren;
    TemporalPstar.GrandParentGlobalID = PstarIOCompactDouble.GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TemporalPstar;

}

struct StructPstarIOCompactMix CopyPstarToTemporalStructureCompactMix(const int index){

	struct StructPstarIOCompactMix TempPstarIOCompactMix;

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPstarIOCompactMix.TypeII = Pstar[index]->TypeII;
    TempPstarIOCompactMix.TypeIa = Pstar[index]->TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TempPstarIOCompactMix.TypeIIProb = Pstar[index]->TypeIIProb;
#endif

    TempPstarIOCompactMix.IMFTYPE = Pstar[index]->IMFTYPE;
    TempPstarIOCompactMix.NthChildren = Pstar[index]->NthChildren;
    TempPstarIOCompactMix.ParentGlobalID = Pstar[index]->ParentGlobalID;

    TempPstarIOCompactMix.InitialMass = (float)Pstar[index]->InitialMass;
    TempPstarIOCompactMix.Mass = (float)Pstar[index]->Mass;
    TempPstarIOCompactMix.FormationTime = (double)Pstar[index]->FormationTime;
    TempPstarIOCompactMix.Z = (double)Pstar[index]->Z;
    TempPstarIOCompactMix.ZII = (double)Pstar[index]->ZII;
    TempPstarIOCompactMix.ZIa = (double)Pstar[index]->ZIa;
    TempPstarIOCompactMix.TempForm = (double)Pstar[index]->TempForm;
    TempPstarIOCompactMix.RhoForm = (double)Pstar[index]->RhoForm;
    TempPstarIOCompactMix.StromgrenRadius = (double)Pstar[index]->StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TempPstarIOCompactMix.SNIICount = Pstar[index]->SNIICount;
    TempPstarIOCompactMix.SNIaCount = Pstar[index]->SNIaCount;
    TempPstarIOCompactMix.EventTimeSNII = (double)Pstar[index]->EventTimeSNII;
    TempPstarIOCompactMix.EventTimeSNIa = (double)Pstar[index]->EventTimeSNIa;

    TempPstarIOCompactMix.AGBCount = Pstar[index]->AGBCount;
    TempPstarIOCompactMix.EventTimeAGB = (double)Pstar[index]->EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TempPstarIOCompactMix.NSMCount = Pstar[index]->NSMCount;
    TempPstarIOCompactMix.EventTimeNSM = (double)Pstar[index]->EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TempPstarIOCompactMix.Elements[0] = (float)Pstar[index]->Elements[0];
    TempPstarIOCompactMix.Elements[1] = (float)Pstar[index]->Elements[1];
    TempPstarIOCompactMix.Elements[2] = (float)Pstar[index]->Elements[2];
    TempPstarIOCompactMix.Elements[3] = (float)Pstar[index]->Elements[3];
    TempPstarIOCompactMix.Elements[4] = (float)Pstar[index]->Elements[4];
    TempPstarIOCompactMix.Elements[5] = (float)Pstar[index]->Elements[5];
    TempPstarIOCompactMix.Elements[6] = (float)Pstar[index]->Elements[6];
    TempPstarIOCompactMix.Elements[7] = (float)Pstar[index]->Elements[7];
    TempPstarIOCompactMix.Elements[8] = (float)Pstar[index]->Elements[8];
    TempPstarIOCompactMix.Elements[9] = (float)Pstar[index]->Elements[9];
    TempPstarIOCompactMix.Elements[10] = (float)Pstar[index]->Elements[10];
    TempPstarIOCompactMix.Elements[11] = (float)Pstar[index]->Elements[11];
    TempPstarIOCompactMix.Elements[12] = (float)Pstar[index]->Elements[12];
    TempPstarIOCompactMix.Elements[13] = (float)Pstar[index]->Elements[13];
    TempPstarIOCompactMix.Elements[14] = (float)Pstar[index]->Elements[14];
    TempPstarIOCompactMix.Elements[15] = (float)Pstar[index]->Elements[15];
    TempPstarIOCompactMix.Elements[16] = (float)Pstar[index]->Elements[16];
    TempPstarIOCompactMix.Elements[17] = (float)Pstar[index]->Elements[17];
#endif // USE_CELIB //}


    TempPstarIOCompactMix.LFUV = (double)Pstar[index]->LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TempPstarIOCompactMix.BolometricLuminosity = (double)Pstar[index]->BolometricLuminosity;
    TempPstarIOCompactMix.RadiusRP = (double)Pstar[index]->RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TempPstarIOCompactMix.StellarWindEnergy = (double)Pstar[index]->StellarWindEnergy;
    TempPstarIOCompactMix.RadiusSW = (double)Pstar[index]->RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPstarIOCompactMix.nH_ave = (double)Pstar[index]->nH_ave;
    TempPstarIOCompactMix.Z_ave = (double)Pstar[index]->Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TempPstarIOCompactMix.SplitGeneration = Pstar[index]->SplitGeneration;
    TempPstarIOCompactMix.SplitNthChildren = Pstar[index]->SplitNthChildren;
    TempPstarIOCompactMix.GrandParentGlobalID = Pstar[index]->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPstarIOCompactMix;

}

struct StructPstarIOCompactMix CopyPstarToTemporalStructureCompactMixElement(StructPstarptr const Ps){

	struct StructPstarIOCompactMix TempPstarIOCompactMix;

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPstarIOCompactMix.TypeII = Ps->TypeII;
    TempPstarIOCompactMix.TypeIa = Ps->TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TempPstarIOCompactMix.TypeIIProb = Ps->TypeIIProb;
#endif

    TempPstarIOCompactMix.IMFTYPE = Ps->IMFTYPE;
    TempPstarIOCompactMix.NthChildren = Ps->NthChildren;
    TempPstarIOCompactMix.ParentGlobalID = Ps->ParentGlobalID;

    TempPstarIOCompactMix.InitialMass = Ps->InitialMass;
    TempPstarIOCompactMix.Mass = Ps->Mass;
    TempPstarIOCompactMix.FormationTime = Ps->FormationTime;
    TempPstarIOCompactMix.Z = Ps->Z;
    TempPstarIOCompactMix.ZII = Ps->ZII;
    TempPstarIOCompactMix.ZIa = Ps->ZIa;
    TempPstarIOCompactMix.TempForm = Ps->TempForm;
    TempPstarIOCompactMix.RhoForm = Ps->RhoForm;
    TempPstarIOCompactMix.StromgrenRadius = Ps->StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TempPstarIOCompactMix.SNIICount = Ps->SNIICount;
    TempPstarIOCompactMix.SNIaCount = Ps->SNIaCount;
    TempPstarIOCompactMix.EventTimeSNII = Ps->EventTimeSNII;
    TempPstarIOCompactMix.EventTimeSNIa = Ps->EventTimeSNIa;

    TempPstarIOCompactMix.AGBCount = Ps->AGBCount;
    TempPstarIOCompactMix.EventTimeAGB = Ps->EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TempPstarIOCompactMix.NSMCount = Ps->NSMCount;
    TempPstarIOCompactMix.EventTimeNSM = Ps->EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TempPstarIOCompactMix.Elements[0] = Ps->Elements[0];
    TempPstarIOCompactMix.Elements[1] = Ps->Elements[1];
    TempPstarIOCompactMix.Elements[2] = Ps->Elements[2];
    TempPstarIOCompactMix.Elements[3] = Ps->Elements[3];
    TempPstarIOCompactMix.Elements[4] = Ps->Elements[4];
    TempPstarIOCompactMix.Elements[5] = Ps->Elements[5];
    TempPstarIOCompactMix.Elements[6] = Ps->Elements[6];
    TempPstarIOCompactMix.Elements[7] = Ps->Elements[7];
    TempPstarIOCompactMix.Elements[8] = Ps->Elements[8];
    TempPstarIOCompactMix.Elements[9] = Ps->Elements[9];
    TempPstarIOCompactMix.Elements[10] = Ps->Elements[10];
    TempPstarIOCompactMix.Elements[11] = Ps->Elements[11];
    TempPstarIOCompactMix.Elements[12] = Ps->Elements[12];
    TempPstarIOCompactMix.Elements[13] = Ps->Elements[13];
    TempPstarIOCompactMix.Elements[14] = Ps->Elements[14];
    TempPstarIOCompactMix.Elements[15] = Ps->Elements[15];
    TempPstarIOCompactMix.Elements[16] = Ps->Elements[16];
    TempPstarIOCompactMix.Elements[17] = Ps->Elements[17];
#endif // USE_CELIB //}


    TempPstarIOCompactMix.LFUV = Ps->LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TempPstarIOCompactMix.BolometricLuminosity = Ps->BolometricLuminosity;
    TempPstarIOCompactMix.RadiusRP = Ps->RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TempPstarIOCompactMix.StellarWindEnergy = Ps->StellarWindEnergy;
    TempPstarIOCompactMix.RadiusSW = Ps->RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPstarIOCompactMix.nH_ave = Ps->nH_ave;
    TempPstarIOCompactMix.Z_ave = Ps->Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TempPstarIOCompactMix.SplitGeneration = Ps->SplitGeneration;
    TempPstarIOCompactMix.SplitNthChildren = Ps->SplitNthChildren;
    TempPstarIOCompactMix.GrandParentGlobalID = Ps->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPstarIOCompactMix;

}

StructPstar CopyTemporalStructureCompactMixToPstarCompactMix(struct StructPstarIOCompactMix PstarIOCompactMix){

	StructPstar TemporalPstar;
	memset(&TemporalPstar,0,sizeof(StructPstar));

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TemporalPstar.TypeII = PstarIOCompactMix.TypeII;
    TemporalPstar.TypeIa = PstarIOCompactMix.TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TemporalPstar.TypeIIProb = PstarIOCompactMix.TypeIIProb;
#endif

    TemporalPstar.IMFTYPE = PstarIOCompactMix.IMFTYPE;
    TemporalPstar.NthChildren = PstarIOCompactMix.NthChildren;
    TemporalPstar.ParentGlobalID = PstarIOCompactMix.ParentGlobalID;

    TemporalPstar.InitialMass = PstarIOCompactMix.InitialMass;
    TemporalPstar.Mass = PstarIOCompactMix.Mass;
    TemporalPstar.FormationTime = PstarIOCompactMix.FormationTime;
    TemporalPstar.Z = PstarIOCompactMix.Z;
    TemporalPstar.ZII = PstarIOCompactMix.ZII;
    TemporalPstar.ZIa = PstarIOCompactMix.ZIa;
    TemporalPstar.TempForm = PstarIOCompactMix.TempForm;
    TemporalPstar.RhoForm = PstarIOCompactMix.RhoForm;
    TemporalPstar.StromgrenRadius = PstarIOCompactMix.StromgrenRadius;
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TemporalPstar.SNIICount = PstarIOCompactMix.SNIICount;
    TemporalPstar.SNIaCount = PstarIOCompactMix.SNIaCount;
    TemporalPstar.EventTimeSNII = PstarIOCompactMix.EventTimeSNII;
    TemporalPstar.EventTimeSNIa = PstarIOCompactMix.EventTimeSNIa;

    TemporalPstar.AGBCount = PstarIOCompactMix.AGBCount;
    TemporalPstar.EventTimeAGB = PstarIOCompactMix.EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TemporalPstar.NSMCount = PstarIOCompactMix.NSMCount;
    TemporalPstar.EventTimeNSM = PstarIOCompactMix.EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TemporalPstar.Elements[0] = PstarIOCompactMix.Elements[0];
    TemporalPstar.Elements[1] = PstarIOCompactMix.Elements[1];
    TemporalPstar.Elements[2] = PstarIOCompactMix.Elements[2];
    TemporalPstar.Elements[3] = PstarIOCompactMix.Elements[3];
    TemporalPstar.Elements[4] = PstarIOCompactMix.Elements[4];
    TemporalPstar.Elements[5] = PstarIOCompactMix.Elements[5];
    TemporalPstar.Elements[6] = PstarIOCompactMix.Elements[6];
    TemporalPstar.Elements[7] = PstarIOCompactMix.Elements[7];
    TemporalPstar.Elements[8] = PstarIOCompactMix.Elements[8];
    TemporalPstar.Elements[9] = PstarIOCompactMix.Elements[9];
    TemporalPstar.Elements[10] = PstarIOCompactMix.Elements[10];
    TemporalPstar.Elements[11] = PstarIOCompactMix.Elements[11];
    TemporalPstar.Elements[12] = PstarIOCompactMix.Elements[12];
    TemporalPstar.Elements[13] = PstarIOCompactMix.Elements[13];
    TemporalPstar.Elements[14] = PstarIOCompactMix.Elements[14];
    TemporalPstar.Elements[15] = PstarIOCompactMix.Elements[15];
    TemporalPstar.Elements[16] = PstarIOCompactMix.Elements[16];
    TemporalPstar.Elements[17] = PstarIOCompactMix.Elements[17];
#endif // USE_CELIB //}


    TemporalPstar.LFUV = PstarIOCompactMix.LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TemporalPstar.BolometricLuminosity = PstarIOCompactMix.BolometricLuminosity;
    TemporalPstar.RadiusRP = PstarIOCompactMix.RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TemporalPstar.StellarWindEnergy = PstarIOCompactMix.StellarWindEnergy;
    TemporalPstar.RadiusSW = PstarIOCompactMix.RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TemporalPstar.nH_ave = PstarIOCompactMix.nH_ave;
    TemporalPstar.Z_ave = PstarIOCompactMix.Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPstar.SplitGeneration = PstarIOCompactMix.SplitGeneration;
    TemporalPstar.SplitNthChildren = PstarIOCompactMix.SplitNthChildren;
    TemporalPstar.GrandParentGlobalID = PstarIOCompactMix.GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TemporalPstar;

}

struct StructPstarIOLean CopyPstarToTemporalStructureLean(const int index){

	struct StructPstarIOLean TempPstarIOLean;

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPstarIOLean.TypeII = Pstar[index]->TypeII;
    TempPstarIOLean.TypeIa = Pstar[index]->TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TempPstarIOLean.TypeIIProb = Pstar[index]->TypeIIProb;
#endif

    TempPstarIOLean.IMFTYPE = Pstar[index]->IMFTYPE;
    TempPstarIOLean.NthChildren = Pstar[index]->NthChildren;
    TempPstarIOLean.ParentGlobalID = Pstar[index]->ParentGlobalID;

// omit double  InitialMass;    // <MIX> <LEAN> Initial Mass
    TempPstarIOLean.Mass = (float)Pstar[index]->Mass;
    TempPstarIOLean.FormationTime = (float)Pstar[index]->FormationTime;
    TempPstarIOLean.Z = (float)Pstar[index]->Z;
    TempPstarIOLean.ZII = (float)Pstar[index]->ZII;
    TempPstarIOLean.ZIa = (float)Pstar[index]->ZIa;
// omit double  TempForm;       // <LEAN> Temperature(t=FormationTime)
// omit double  RhoForm;        // <LEAN> Density(t=FormationTime) // UnitMass/UnitLength^3
// omit double  StromgrenRadius;// <LEAN>
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TempPstarIOLean.SNIICount = Pstar[index]->SNIICount;
    TempPstarIOLean.SNIaCount = Pstar[index]->SNIaCount;
    TempPstarIOLean.EventTimeSNII = (float)Pstar[index]->EventTimeSNII;
    TempPstarIOLean.EventTimeSNIa = (float)Pstar[index]->EventTimeSNIa;

    TempPstarIOLean.AGBCount = Pstar[index]->AGBCount;
    TempPstarIOLean.EventTimeAGB = (float)Pstar[index]->EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TempPstarIOLean.NSMCount = Pstar[index]->NSMCount;
    TempPstarIOLean.EventTimeNSM = (float)Pstar[index]->EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TempPstarIOLean.Elements[0] = (float)Pstar[index]->Elements[0];
    TempPstarIOLean.Elements[1] = (float)Pstar[index]->Elements[1];
    TempPstarIOLean.Elements[2] = (float)Pstar[index]->Elements[2];
    TempPstarIOLean.Elements[3] = (float)Pstar[index]->Elements[3];
    TempPstarIOLean.Elements[4] = (float)Pstar[index]->Elements[4];
    TempPstarIOLean.Elements[5] = (float)Pstar[index]->Elements[5];
    TempPstarIOLean.Elements[6] = (float)Pstar[index]->Elements[6];
    TempPstarIOLean.Elements[7] = (float)Pstar[index]->Elements[7];
    TempPstarIOLean.Elements[8] = (float)Pstar[index]->Elements[8];
    TempPstarIOLean.Elements[9] = (float)Pstar[index]->Elements[9];
    TempPstarIOLean.Elements[10] = (float)Pstar[index]->Elements[10];
    TempPstarIOLean.Elements[11] = (float)Pstar[index]->Elements[11];
    TempPstarIOLean.Elements[12] = (float)Pstar[index]->Elements[12];
    TempPstarIOLean.Elements[13] = (float)Pstar[index]->Elements[13];
    TempPstarIOLean.Elements[14] = (float)Pstar[index]->Elements[14];
    TempPstarIOLean.Elements[15] = (float)Pstar[index]->Elements[15];
    TempPstarIOLean.Elements[16] = (float)Pstar[index]->Elements[16];
    TempPstarIOLean.Elements[17] = (float)Pstar[index]->Elements[17];
#endif // USE_CELIB //}


    TempPstarIOLean.LFUV = (float)Pstar[index]->LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TempPstarIOLean.BolometricLuminosity = (float)Pstar[index]->BolometricLuminosity;
    TempPstarIOLean.RadiusRP = (float)Pstar[index]->RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TempPstarIOLean.StellarWindEnergy = (float)Pstar[index]->StellarWindEnergy;
    TempPstarIOLean.RadiusSW = (float)Pstar[index]->RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPstarIOLean.nH_ave = (float)Pstar[index]->nH_ave;
    TempPstarIOLean.Z_ave = (float)Pstar[index]->Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TempPstarIOLean.SplitGeneration = Pstar[index]->SplitGeneration;
    TempPstarIOLean.SplitNthChildren = Pstar[index]->SplitNthChildren;
    TempPstarIOLean.GrandParentGlobalID = Pstar[index]->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPstarIOLean;

}

struct StructPstarIOLean CopyPstarToTemporalStructureLeanElement(StructPstarptr const Ps){

	struct StructPstarIOLean TempPstarIOLean;

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPstarIOLean.TypeII = Ps->TypeII;
    TempPstarIOLean.TypeIa = Ps->TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TempPstarIOLean.TypeIIProb = Ps->TypeIIProb;
#endif

    TempPstarIOLean.IMFTYPE = Ps->IMFTYPE;
    TempPstarIOLean.NthChildren = Ps->NthChildren;
    TempPstarIOLean.ParentGlobalID = Ps->ParentGlobalID;

// omit double  InitialMass;    // <MIX> <LEAN> Initial Mass
    TempPstarIOLean.Mass = Ps->Mass;
    TempPstarIOLean.FormationTime = Ps->FormationTime;
    TempPstarIOLean.Z = Ps->Z;
    TempPstarIOLean.ZII = Ps->ZII;
    TempPstarIOLean.ZIa = Ps->ZIa;
// omit double  TempForm;       // <LEAN> Temperature(t=FormationTime)
// omit double  RhoForm;        // <LEAN> Density(t=FormationTime) // UnitMass/UnitLength^3
// omit double  StromgrenRadius;// <LEAN>
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TempPstarIOLean.SNIICount = Ps->SNIICount;
    TempPstarIOLean.SNIaCount = Ps->SNIaCount;
    TempPstarIOLean.EventTimeSNII = Ps->EventTimeSNII;
    TempPstarIOLean.EventTimeSNIa = Ps->EventTimeSNIa;

    TempPstarIOLean.AGBCount = Ps->AGBCount;
    TempPstarIOLean.EventTimeAGB = Ps->EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TempPstarIOLean.NSMCount = Ps->NSMCount;
    TempPstarIOLean.EventTimeNSM = Ps->EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TempPstarIOLean.Elements[0] = Ps->Elements[0];
    TempPstarIOLean.Elements[1] = Ps->Elements[1];
    TempPstarIOLean.Elements[2] = Ps->Elements[2];
    TempPstarIOLean.Elements[3] = Ps->Elements[3];
    TempPstarIOLean.Elements[4] = Ps->Elements[4];
    TempPstarIOLean.Elements[5] = Ps->Elements[5];
    TempPstarIOLean.Elements[6] = Ps->Elements[6];
    TempPstarIOLean.Elements[7] = Ps->Elements[7];
    TempPstarIOLean.Elements[8] = Ps->Elements[8];
    TempPstarIOLean.Elements[9] = Ps->Elements[9];
    TempPstarIOLean.Elements[10] = Ps->Elements[10];
    TempPstarIOLean.Elements[11] = Ps->Elements[11];
    TempPstarIOLean.Elements[12] = Ps->Elements[12];
    TempPstarIOLean.Elements[13] = Ps->Elements[13];
    TempPstarIOLean.Elements[14] = Ps->Elements[14];
    TempPstarIOLean.Elements[15] = Ps->Elements[15];
    TempPstarIOLean.Elements[16] = Ps->Elements[16];
    TempPstarIOLean.Elements[17] = Ps->Elements[17];
#endif // USE_CELIB //}


    TempPstarIOLean.LFUV = Ps->LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TempPstarIOLean.BolometricLuminosity = Ps->BolometricLuminosity;
    TempPstarIOLean.RadiusRP = Ps->RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TempPstarIOLean.StellarWindEnergy = Ps->StellarWindEnergy;
    TempPstarIOLean.RadiusSW = Ps->RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TempPstarIOLean.nH_ave = Ps->nH_ave;
    TempPstarIOLean.Z_ave = Ps->Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TempPstarIOLean.SplitGeneration = Ps->SplitGeneration;
    TempPstarIOLean.SplitNthChildren = Ps->SplitNthChildren;
    TempPstarIOLean.GrandParentGlobalID = Ps->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPstarIOLean;

}

StructPstar CopyTemporalStructureLeanToPstarLean(struct StructPstarIOLean PstarIOLean){

	StructPstar TemporalPstar;
	memset(&TemporalPstar,0,sizeof(StructPstar));

// omit StructPstarptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TemporalPstar.TypeII = PstarIOLean.TypeII;
    TemporalPstar.TypeIa = PstarIOLean.TypeIa;
// omit bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    TemporalPstar.TypeIIProb = PstarIOLean.TypeIIProb;
#endif

    TemporalPstar.IMFTYPE = PstarIOLean.IMFTYPE;
    TemporalPstar.NthChildren = PstarIOLean.NthChildren;
    TemporalPstar.ParentGlobalID = PstarIOLean.ParentGlobalID;

// omit double  InitialMass;    // <MIX> <LEAN> Initial Mass
    TemporalPstar.Mass = PstarIOLean.Mass;
    TemporalPstar.FormationTime = PstarIOLean.FormationTime;
    TemporalPstar.Z = PstarIOLean.Z;
    TemporalPstar.ZII = PstarIOLean.ZII;
    TemporalPstar.ZIa = PstarIOLean.ZIa;
// omit double  TempForm;       // <LEAN> Temperature(t=FormationTime)
// omit double  RhoForm;        // <LEAN> Density(t=FormationTime) // UnitMass/UnitLength^3
// omit double  StromgrenRadius;// <LEAN>
// omit double  Density;        // <TMP> <LEAN>
#ifdef USE_CELIB //{
    TemporalPstar.SNIICount = PstarIOLean.SNIICount;
    TemporalPstar.SNIaCount = PstarIOLean.SNIaCount;
    TemporalPstar.EventTimeSNII = PstarIOLean.EventTimeSNII;
    TemporalPstar.EventTimeSNIa = PstarIOLean.EventTimeSNIa;

    TemporalPstar.AGBCount = PstarIOLean.AGBCount;
    TemporalPstar.EventTimeAGB = PstarIOLean.EventTimeAGB;

#ifdef USE_CELIB_NSM //{
    TemporalPstar.NSMCount = PstarIOLean.NSMCount;
    TemporalPstar.EventTimeNSM = PstarIOLean.EventTimeNSM;
#endif // USE_CELIB_NSM //}
    TemporalPstar.Elements[0] = PstarIOLean.Elements[0];
    TemporalPstar.Elements[1] = PstarIOLean.Elements[1];
    TemporalPstar.Elements[2] = PstarIOLean.Elements[2];
    TemporalPstar.Elements[3] = PstarIOLean.Elements[3];
    TemporalPstar.Elements[4] = PstarIOLean.Elements[4];
    TemporalPstar.Elements[5] = PstarIOLean.Elements[5];
    TemporalPstar.Elements[6] = PstarIOLean.Elements[6];
    TemporalPstar.Elements[7] = PstarIOLean.Elements[7];
    TemporalPstar.Elements[8] = PstarIOLean.Elements[8];
    TemporalPstar.Elements[9] = PstarIOLean.Elements[9];
    TemporalPstar.Elements[10] = PstarIOLean.Elements[10];
    TemporalPstar.Elements[11] = PstarIOLean.Elements[11];
    TemporalPstar.Elements[12] = PstarIOLean.Elements[12];
    TemporalPstar.Elements[13] = PstarIOLean.Elements[13];
    TemporalPstar.Elements[14] = PstarIOLean.Elements[14];
    TemporalPstar.Elements[15] = PstarIOLean.Elements[15];
    TemporalPstar.Elements[16] = PstarIOLean.Elements[16];
    TemporalPstar.Elements[17] = PstarIOLean.Elements[17];
#endif // USE_CELIB //}


    TemporalPstar.LFUV = PstarIOLean.LFUV;


#ifdef USE_RADIATION_PRESSURE //{
    TemporalPstar.BolometricLuminosity = PstarIOLean.BolometricLuminosity;
    TemporalPstar.RadiusRP = PstarIOLean.RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    TemporalPstar.StellarWindEnergy = PstarIOLean.StellarWindEnergy;
    TemporalPstar.RadiusSW = PstarIOLean.RadiusSW;
#endif //USE_STELLAR_WIND //{

#ifdef USE_FEEDBACK_TIMESTEP //{
// omit double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    TemporalPstar.nH_ave = PstarIOLean.nH_ave;
    TemporalPstar.Z_ave = PstarIOLean.Z_ave;
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPstar.SplitGeneration = PstarIOLean.SplitGeneration;
    TemporalPstar.SplitNthChildren = PstarIOLean.SplitNthChildren;
    TemporalPstar.GrandParentGlobalID = PstarIOLean.GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TemporalPstar;

}

struct StructPsinkIOCompact CopyPsinkToTemporalStructureCompact(const int index){

	struct StructPsinkIOCompact TempPsinkIOCompact;

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPsinkIOCompact.ParentGlobalID = Psink[index]->ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    TempPsinkIOCompact.dt_localmin = (float)Psink[index]->dt_localmin;

    TempPsinkIOCompact.NumberofAbsorbedParticles = Psink[index]->NumberofAbsorbedParticles;
    TempPsinkIOCompact.AccretionRadius = (float)Psink[index]->AccretionRadius;
    TempPsinkIOCompact.MergingDistance = (float)Psink[index]->MergingDistance;
    TempPsinkIOCompact.AM[0] = (float)Psink[index]->AM[0];
    TempPsinkIOCompact.AM[1] = (float)Psink[index]->AM[1];
    TempPsinkIOCompact.AM[2] = (float)Psink[index]->AM[2];
    TempPsinkIOCompact.FormationTime = (float)Psink[index]->FormationTime;
    TempPsinkIOCompact.Z = (float)Psink[index]->Z;
    TempPsinkIOCompact.ZII = (float)Psink[index]->ZII;
    TempPsinkIOCompact.ZIa = (float)Psink[index]->ZIa;
    TempPsinkIOCompact.AccretionMass = (float)Psink[index]->AccretionMass;
    TempPsinkIOCompact.AccretionMassGas = (float)Psink[index]->AccretionMassGas;
    TempPsinkIOCompact.AccretionMassStar = (float)Psink[index]->AccretionMassStar;
    TempPsinkIOCompact.AccretionMassToBH = (float)Psink[index]->AccretionMassToBH;

#ifdef USE_CELIB
    TempPsinkIOCompact.Elements[0] = (float)Psink[index]->Elements[0];
    TempPsinkIOCompact.Elements[1] = (float)Psink[index]->Elements[1];
    TempPsinkIOCompact.Elements[2] = (float)Psink[index]->Elements[2];
    TempPsinkIOCompact.Elements[3] = (float)Psink[index]->Elements[3];
    TempPsinkIOCompact.Elements[4] = (float)Psink[index]->Elements[4];
    TempPsinkIOCompact.Elements[5] = (float)Psink[index]->Elements[5];
    TempPsinkIOCompact.Elements[6] = (float)Psink[index]->Elements[6];
    TempPsinkIOCompact.Elements[7] = (float)Psink[index]->Elements[7];
    TempPsinkIOCompact.Elements[8] = (float)Psink[index]->Elements[8];
    TempPsinkIOCompact.Elements[9] = (float)Psink[index]->Elements[9];
    TempPsinkIOCompact.Elements[10] = (float)Psink[index]->Elements[10];
    TempPsinkIOCompact.Elements[11] = (float)Psink[index]->Elements[11];
    TempPsinkIOCompact.Elements[12] = (float)Psink[index]->Elements[12];
    TempPsinkIOCompact.Elements[13] = (float)Psink[index]->Elements[13];
    TempPsinkIOCompact.Elements[14] = (float)Psink[index]->Elements[14];
    TempPsinkIOCompact.Elements[15] = (float)Psink[index]->Elements[15];
    TempPsinkIOCompact.Elements[16] = (float)Psink[index]->Elements[16];
    TempPsinkIOCompact.Elements[17] = (float)Psink[index]->Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TempPsinkIOCompact.SplitGeneration = Psink[index]->SplitGeneration;
    TempPsinkIOCompact.SplitNthChildren = Psink[index]->SplitNthChildren;
    TempPsinkIOCompact.GrandParentGlobalID = Psink[index]->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPsinkIOCompact;

}

struct StructPsinkIOCompact CopyPsinkToTemporalStructureCompactElement(StructPsinkptr const Psk){

	struct StructPsinkIOCompact TempPsinkIOCompact;

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPsinkIOCompact.ParentGlobalID = Psk->ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    TempPsinkIOCompact.dt_localmin = Psk->dt_localmin;

    TempPsinkIOCompact.NumberofAbsorbedParticles = Psk->NumberofAbsorbedParticles;
    TempPsinkIOCompact.AccretionRadius = Psk->AccretionRadius;
    TempPsinkIOCompact.MergingDistance = Psk->MergingDistance;
    TempPsinkIOCompact.AM[0] = Psk->AM[0];
    TempPsinkIOCompact.AM[1] = Psk->AM[1];
    TempPsinkIOCompact.AM[2] = Psk->AM[2];
    TempPsinkIOCompact.FormationTime = Psk->FormationTime;
    TempPsinkIOCompact.Z = Psk->Z;
    TempPsinkIOCompact.ZII = Psk->ZII;
    TempPsinkIOCompact.ZIa = Psk->ZIa;
    TempPsinkIOCompact.AccretionMass = Psk->AccretionMass;
    TempPsinkIOCompact.AccretionMassGas = Psk->AccretionMassGas;
    TempPsinkIOCompact.AccretionMassStar = Psk->AccretionMassStar;
    TempPsinkIOCompact.AccretionMassToBH = Psk->AccretionMassToBH;

#ifdef USE_CELIB
    TempPsinkIOCompact.Elements[0] = Psk->Elements[0];
    TempPsinkIOCompact.Elements[1] = Psk->Elements[1];
    TempPsinkIOCompact.Elements[2] = Psk->Elements[2];
    TempPsinkIOCompact.Elements[3] = Psk->Elements[3];
    TempPsinkIOCompact.Elements[4] = Psk->Elements[4];
    TempPsinkIOCompact.Elements[5] = Psk->Elements[5];
    TempPsinkIOCompact.Elements[6] = Psk->Elements[6];
    TempPsinkIOCompact.Elements[7] = Psk->Elements[7];
    TempPsinkIOCompact.Elements[8] = Psk->Elements[8];
    TempPsinkIOCompact.Elements[9] = Psk->Elements[9];
    TempPsinkIOCompact.Elements[10] = Psk->Elements[10];
    TempPsinkIOCompact.Elements[11] = Psk->Elements[11];
    TempPsinkIOCompact.Elements[12] = Psk->Elements[12];
    TempPsinkIOCompact.Elements[13] = Psk->Elements[13];
    TempPsinkIOCompact.Elements[14] = Psk->Elements[14];
    TempPsinkIOCompact.Elements[15] = Psk->Elements[15];
    TempPsinkIOCompact.Elements[16] = Psk->Elements[16];
    TempPsinkIOCompact.Elements[17] = Psk->Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TempPsinkIOCompact.SplitGeneration = Psk->SplitGeneration;
    TempPsinkIOCompact.SplitNthChildren = Psk->SplitNthChildren;
    TempPsinkIOCompact.GrandParentGlobalID = Psk->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPsinkIOCompact;

}

StructPsink CopyTemporalStructureCompactToPsinkCompact(struct StructPsinkIOCompact PsinkIOCompact){

	StructPsink TemporalPsink;
	memset(&TemporalPsink,0,sizeof(StructPsink));

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TemporalPsink.ParentGlobalID = PsinkIOCompact.ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    TemporalPsink.dt_localmin = PsinkIOCompact.dt_localmin;

    TemporalPsink.NumberofAbsorbedParticles = PsinkIOCompact.NumberofAbsorbedParticles;
    TemporalPsink.AccretionRadius = PsinkIOCompact.AccretionRadius;
    TemporalPsink.MergingDistance = PsinkIOCompact.MergingDistance;
    TemporalPsink.AM[0] = PsinkIOCompact.AM[0];
    TemporalPsink.AM[1] = PsinkIOCompact.AM[1];
    TemporalPsink.AM[2] = PsinkIOCompact.AM[2];
    TemporalPsink.FormationTime = PsinkIOCompact.FormationTime;
    TemporalPsink.Z = PsinkIOCompact.Z;
    TemporalPsink.ZII = PsinkIOCompact.ZII;
    TemporalPsink.ZIa = PsinkIOCompact.ZIa;
    TemporalPsink.AccretionMass = PsinkIOCompact.AccretionMass;
    TemporalPsink.AccretionMassGas = PsinkIOCompact.AccretionMassGas;
    TemporalPsink.AccretionMassStar = PsinkIOCompact.AccretionMassStar;
    TemporalPsink.AccretionMassToBH = PsinkIOCompact.AccretionMassToBH;

#ifdef USE_CELIB
    TemporalPsink.Elements[0] = PsinkIOCompact.Elements[0];
    TemporalPsink.Elements[1] = PsinkIOCompact.Elements[1];
    TemporalPsink.Elements[2] = PsinkIOCompact.Elements[2];
    TemporalPsink.Elements[3] = PsinkIOCompact.Elements[3];
    TemporalPsink.Elements[4] = PsinkIOCompact.Elements[4];
    TemporalPsink.Elements[5] = PsinkIOCompact.Elements[5];
    TemporalPsink.Elements[6] = PsinkIOCompact.Elements[6];
    TemporalPsink.Elements[7] = PsinkIOCompact.Elements[7];
    TemporalPsink.Elements[8] = PsinkIOCompact.Elements[8];
    TemporalPsink.Elements[9] = PsinkIOCompact.Elements[9];
    TemporalPsink.Elements[10] = PsinkIOCompact.Elements[10];
    TemporalPsink.Elements[11] = PsinkIOCompact.Elements[11];
    TemporalPsink.Elements[12] = PsinkIOCompact.Elements[12];
    TemporalPsink.Elements[13] = PsinkIOCompact.Elements[13];
    TemporalPsink.Elements[14] = PsinkIOCompact.Elements[14];
    TemporalPsink.Elements[15] = PsinkIOCompact.Elements[15];
    TemporalPsink.Elements[16] = PsinkIOCompact.Elements[16];
    TemporalPsink.Elements[17] = PsinkIOCompact.Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPsink.SplitGeneration = PsinkIOCompact.SplitGeneration;
    TemporalPsink.SplitNthChildren = PsinkIOCompact.SplitNthChildren;
    TemporalPsink.GrandParentGlobalID = PsinkIOCompact.GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TemporalPsink;

}

struct StructPsinkIOCompactDouble CopyPsinkToTemporalStructureCompactDouble(const int index){

	struct StructPsinkIOCompactDouble TempPsinkIOCompactDouble;

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPsinkIOCompactDouble.ParentGlobalID = Psink[index]->ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    TempPsinkIOCompactDouble.dt_localmin = (double)Psink[index]->dt_localmin;

    TempPsinkIOCompactDouble.NumberofAbsorbedParticles = Psink[index]->NumberofAbsorbedParticles;
    TempPsinkIOCompactDouble.AccretionRadius = (double)Psink[index]->AccretionRadius;
    TempPsinkIOCompactDouble.MergingDistance = (double)Psink[index]->MergingDistance;
    TempPsinkIOCompactDouble.AM[0] = (double)Psink[index]->AM[0];
    TempPsinkIOCompactDouble.AM[1] = (double)Psink[index]->AM[1];
    TempPsinkIOCompactDouble.AM[2] = (double)Psink[index]->AM[2];
    TempPsinkIOCompactDouble.FormationTime = (double)Psink[index]->FormationTime;
    TempPsinkIOCompactDouble.Z = (double)Psink[index]->Z;
    TempPsinkIOCompactDouble.ZII = (double)Psink[index]->ZII;
    TempPsinkIOCompactDouble.ZIa = (double)Psink[index]->ZIa;
    TempPsinkIOCompactDouble.AccretionMass = (double)Psink[index]->AccretionMass;
    TempPsinkIOCompactDouble.AccretionMassGas = (double)Psink[index]->AccretionMassGas;
    TempPsinkIOCompactDouble.AccretionMassStar = (double)Psink[index]->AccretionMassStar;
    TempPsinkIOCompactDouble.AccretionMassToBH = (double)Psink[index]->AccretionMassToBH;

#ifdef USE_CELIB
    TempPsinkIOCompactDouble.Elements[0] = (double)Psink[index]->Elements[0];
    TempPsinkIOCompactDouble.Elements[1] = (double)Psink[index]->Elements[1];
    TempPsinkIOCompactDouble.Elements[2] = (double)Psink[index]->Elements[2];
    TempPsinkIOCompactDouble.Elements[3] = (double)Psink[index]->Elements[3];
    TempPsinkIOCompactDouble.Elements[4] = (double)Psink[index]->Elements[4];
    TempPsinkIOCompactDouble.Elements[5] = (double)Psink[index]->Elements[5];
    TempPsinkIOCompactDouble.Elements[6] = (double)Psink[index]->Elements[6];
    TempPsinkIOCompactDouble.Elements[7] = (double)Psink[index]->Elements[7];
    TempPsinkIOCompactDouble.Elements[8] = (double)Psink[index]->Elements[8];
    TempPsinkIOCompactDouble.Elements[9] = (double)Psink[index]->Elements[9];
    TempPsinkIOCompactDouble.Elements[10] = (double)Psink[index]->Elements[10];
    TempPsinkIOCompactDouble.Elements[11] = (double)Psink[index]->Elements[11];
    TempPsinkIOCompactDouble.Elements[12] = (double)Psink[index]->Elements[12];
    TempPsinkIOCompactDouble.Elements[13] = (double)Psink[index]->Elements[13];
    TempPsinkIOCompactDouble.Elements[14] = (double)Psink[index]->Elements[14];
    TempPsinkIOCompactDouble.Elements[15] = (double)Psink[index]->Elements[15];
    TempPsinkIOCompactDouble.Elements[16] = (double)Psink[index]->Elements[16];
    TempPsinkIOCompactDouble.Elements[17] = (double)Psink[index]->Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TempPsinkIOCompactDouble.SplitGeneration = Psink[index]->SplitGeneration;
    TempPsinkIOCompactDouble.SplitNthChildren = Psink[index]->SplitNthChildren;
    TempPsinkIOCompactDouble.GrandParentGlobalID = Psink[index]->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPsinkIOCompactDouble;

}

struct StructPsinkIOCompactDouble CopyPsinkToTemporalStructureCompactDoubleElement(StructPsinkptr const Psk){

	struct StructPsinkIOCompactDouble TempPsinkIOCompactDouble;

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPsinkIOCompactDouble.ParentGlobalID = Psk->ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    TempPsinkIOCompactDouble.dt_localmin = Psk->dt_localmin;

    TempPsinkIOCompactDouble.NumberofAbsorbedParticles = Psk->NumberofAbsorbedParticles;
    TempPsinkIOCompactDouble.AccretionRadius = Psk->AccretionRadius;
    TempPsinkIOCompactDouble.MergingDistance = Psk->MergingDistance;
    TempPsinkIOCompactDouble.AM[0] = Psk->AM[0];
    TempPsinkIOCompactDouble.AM[1] = Psk->AM[1];
    TempPsinkIOCompactDouble.AM[2] = Psk->AM[2];
    TempPsinkIOCompactDouble.FormationTime = Psk->FormationTime;
    TempPsinkIOCompactDouble.Z = Psk->Z;
    TempPsinkIOCompactDouble.ZII = Psk->ZII;
    TempPsinkIOCompactDouble.ZIa = Psk->ZIa;
    TempPsinkIOCompactDouble.AccretionMass = Psk->AccretionMass;
    TempPsinkIOCompactDouble.AccretionMassGas = Psk->AccretionMassGas;
    TempPsinkIOCompactDouble.AccretionMassStar = Psk->AccretionMassStar;
    TempPsinkIOCompactDouble.AccretionMassToBH = Psk->AccretionMassToBH;

#ifdef USE_CELIB
    TempPsinkIOCompactDouble.Elements[0] = Psk->Elements[0];
    TempPsinkIOCompactDouble.Elements[1] = Psk->Elements[1];
    TempPsinkIOCompactDouble.Elements[2] = Psk->Elements[2];
    TempPsinkIOCompactDouble.Elements[3] = Psk->Elements[3];
    TempPsinkIOCompactDouble.Elements[4] = Psk->Elements[4];
    TempPsinkIOCompactDouble.Elements[5] = Psk->Elements[5];
    TempPsinkIOCompactDouble.Elements[6] = Psk->Elements[6];
    TempPsinkIOCompactDouble.Elements[7] = Psk->Elements[7];
    TempPsinkIOCompactDouble.Elements[8] = Psk->Elements[8];
    TempPsinkIOCompactDouble.Elements[9] = Psk->Elements[9];
    TempPsinkIOCompactDouble.Elements[10] = Psk->Elements[10];
    TempPsinkIOCompactDouble.Elements[11] = Psk->Elements[11];
    TempPsinkIOCompactDouble.Elements[12] = Psk->Elements[12];
    TempPsinkIOCompactDouble.Elements[13] = Psk->Elements[13];
    TempPsinkIOCompactDouble.Elements[14] = Psk->Elements[14];
    TempPsinkIOCompactDouble.Elements[15] = Psk->Elements[15];
    TempPsinkIOCompactDouble.Elements[16] = Psk->Elements[16];
    TempPsinkIOCompactDouble.Elements[17] = Psk->Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TempPsinkIOCompactDouble.SplitGeneration = Psk->SplitGeneration;
    TempPsinkIOCompactDouble.SplitNthChildren = Psk->SplitNthChildren;
    TempPsinkIOCompactDouble.GrandParentGlobalID = Psk->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPsinkIOCompactDouble;

}

StructPsink CopyTemporalStructureCompactDoubleToPsinkCompactDouble(struct StructPsinkIOCompactDouble PsinkIOCompactDouble){

	StructPsink TemporalPsink;
	memset(&TemporalPsink,0,sizeof(StructPsink));

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TemporalPsink.ParentGlobalID = PsinkIOCompactDouble.ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    TemporalPsink.dt_localmin = PsinkIOCompactDouble.dt_localmin;

    TemporalPsink.NumberofAbsorbedParticles = PsinkIOCompactDouble.NumberofAbsorbedParticles;
    TemporalPsink.AccretionRadius = PsinkIOCompactDouble.AccretionRadius;
    TemporalPsink.MergingDistance = PsinkIOCompactDouble.MergingDistance;
    TemporalPsink.AM[0] = PsinkIOCompactDouble.AM[0];
    TemporalPsink.AM[1] = PsinkIOCompactDouble.AM[1];
    TemporalPsink.AM[2] = PsinkIOCompactDouble.AM[2];
    TemporalPsink.FormationTime = PsinkIOCompactDouble.FormationTime;
    TemporalPsink.Z = PsinkIOCompactDouble.Z;
    TemporalPsink.ZII = PsinkIOCompactDouble.ZII;
    TemporalPsink.ZIa = PsinkIOCompactDouble.ZIa;
    TemporalPsink.AccretionMass = PsinkIOCompactDouble.AccretionMass;
    TemporalPsink.AccretionMassGas = PsinkIOCompactDouble.AccretionMassGas;
    TemporalPsink.AccretionMassStar = PsinkIOCompactDouble.AccretionMassStar;
    TemporalPsink.AccretionMassToBH = PsinkIOCompactDouble.AccretionMassToBH;

#ifdef USE_CELIB
    TemporalPsink.Elements[0] = PsinkIOCompactDouble.Elements[0];
    TemporalPsink.Elements[1] = PsinkIOCompactDouble.Elements[1];
    TemporalPsink.Elements[2] = PsinkIOCompactDouble.Elements[2];
    TemporalPsink.Elements[3] = PsinkIOCompactDouble.Elements[3];
    TemporalPsink.Elements[4] = PsinkIOCompactDouble.Elements[4];
    TemporalPsink.Elements[5] = PsinkIOCompactDouble.Elements[5];
    TemporalPsink.Elements[6] = PsinkIOCompactDouble.Elements[6];
    TemporalPsink.Elements[7] = PsinkIOCompactDouble.Elements[7];
    TemporalPsink.Elements[8] = PsinkIOCompactDouble.Elements[8];
    TemporalPsink.Elements[9] = PsinkIOCompactDouble.Elements[9];
    TemporalPsink.Elements[10] = PsinkIOCompactDouble.Elements[10];
    TemporalPsink.Elements[11] = PsinkIOCompactDouble.Elements[11];
    TemporalPsink.Elements[12] = PsinkIOCompactDouble.Elements[12];
    TemporalPsink.Elements[13] = PsinkIOCompactDouble.Elements[13];
    TemporalPsink.Elements[14] = PsinkIOCompactDouble.Elements[14];
    TemporalPsink.Elements[15] = PsinkIOCompactDouble.Elements[15];
    TemporalPsink.Elements[16] = PsinkIOCompactDouble.Elements[16];
    TemporalPsink.Elements[17] = PsinkIOCompactDouble.Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPsink.SplitGeneration = PsinkIOCompactDouble.SplitGeneration;
    TemporalPsink.SplitNthChildren = PsinkIOCompactDouble.SplitNthChildren;
    TemporalPsink.GrandParentGlobalID = PsinkIOCompactDouble.GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TemporalPsink;

}

struct StructPsinkIOCompactMix CopyPsinkToTemporalStructureCompactMix(const int index){

	struct StructPsinkIOCompactMix TempPsinkIOCompactMix;

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPsinkIOCompactMix.ParentGlobalID = Psink[index]->ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    TempPsinkIOCompactMix.dt_localmin = (double)Psink[index]->dt_localmin;

    TempPsinkIOCompactMix.NumberofAbsorbedParticles = Psink[index]->NumberofAbsorbedParticles;
    TempPsinkIOCompactMix.AccretionRadius = (double)Psink[index]->AccretionRadius;
    TempPsinkIOCompactMix.MergingDistance = (double)Psink[index]->MergingDistance;
    TempPsinkIOCompactMix.AM[0] = (double)Psink[index]->AM[0];
    TempPsinkIOCompactMix.AM[1] = (double)Psink[index]->AM[1];
    TempPsinkIOCompactMix.AM[2] = (double)Psink[index]->AM[2];
    TempPsinkIOCompactMix.FormationTime = (double)Psink[index]->FormationTime;
    TempPsinkIOCompactMix.Z = (double)Psink[index]->Z;
    TempPsinkIOCompactMix.ZII = (double)Psink[index]->ZII;
    TempPsinkIOCompactMix.ZIa = (double)Psink[index]->ZIa;
    TempPsinkIOCompactMix.AccretionMass = (double)Psink[index]->AccretionMass;
    TempPsinkIOCompactMix.AccretionMassGas = (double)Psink[index]->AccretionMassGas;
    TempPsinkIOCompactMix.AccretionMassStar = (double)Psink[index]->AccretionMassStar;
    TempPsinkIOCompactMix.AccretionMassToBH = (double)Psink[index]->AccretionMassToBH;

#ifdef USE_CELIB
    TempPsinkIOCompactMix.Elements[0] = (float)Psink[index]->Elements[0];
    TempPsinkIOCompactMix.Elements[1] = (float)Psink[index]->Elements[1];
    TempPsinkIOCompactMix.Elements[2] = (float)Psink[index]->Elements[2];
    TempPsinkIOCompactMix.Elements[3] = (float)Psink[index]->Elements[3];
    TempPsinkIOCompactMix.Elements[4] = (float)Psink[index]->Elements[4];
    TempPsinkIOCompactMix.Elements[5] = (float)Psink[index]->Elements[5];
    TempPsinkIOCompactMix.Elements[6] = (float)Psink[index]->Elements[6];
    TempPsinkIOCompactMix.Elements[7] = (float)Psink[index]->Elements[7];
    TempPsinkIOCompactMix.Elements[8] = (float)Psink[index]->Elements[8];
    TempPsinkIOCompactMix.Elements[9] = (float)Psink[index]->Elements[9];
    TempPsinkIOCompactMix.Elements[10] = (float)Psink[index]->Elements[10];
    TempPsinkIOCompactMix.Elements[11] = (float)Psink[index]->Elements[11];
    TempPsinkIOCompactMix.Elements[12] = (float)Psink[index]->Elements[12];
    TempPsinkIOCompactMix.Elements[13] = (float)Psink[index]->Elements[13];
    TempPsinkIOCompactMix.Elements[14] = (float)Psink[index]->Elements[14];
    TempPsinkIOCompactMix.Elements[15] = (float)Psink[index]->Elements[15];
    TempPsinkIOCompactMix.Elements[16] = (float)Psink[index]->Elements[16];
    TempPsinkIOCompactMix.Elements[17] = (float)Psink[index]->Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TempPsinkIOCompactMix.SplitGeneration = Psink[index]->SplitGeneration;
    TempPsinkIOCompactMix.SplitNthChildren = Psink[index]->SplitNthChildren;
    TempPsinkIOCompactMix.GrandParentGlobalID = Psink[index]->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPsinkIOCompactMix;

}

struct StructPsinkIOCompactMix CopyPsinkToTemporalStructureCompactMixElement(StructPsinkptr const Psk){

	struct StructPsinkIOCompactMix TempPsinkIOCompactMix;

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPsinkIOCompactMix.ParentGlobalID = Psk->ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    TempPsinkIOCompactMix.dt_localmin = Psk->dt_localmin;

    TempPsinkIOCompactMix.NumberofAbsorbedParticles = Psk->NumberofAbsorbedParticles;
    TempPsinkIOCompactMix.AccretionRadius = Psk->AccretionRadius;
    TempPsinkIOCompactMix.MergingDistance = Psk->MergingDistance;
    TempPsinkIOCompactMix.AM[0] = Psk->AM[0];
    TempPsinkIOCompactMix.AM[1] = Psk->AM[1];
    TempPsinkIOCompactMix.AM[2] = Psk->AM[2];
    TempPsinkIOCompactMix.FormationTime = Psk->FormationTime;
    TempPsinkIOCompactMix.Z = Psk->Z;
    TempPsinkIOCompactMix.ZII = Psk->ZII;
    TempPsinkIOCompactMix.ZIa = Psk->ZIa;
    TempPsinkIOCompactMix.AccretionMass = Psk->AccretionMass;
    TempPsinkIOCompactMix.AccretionMassGas = Psk->AccretionMassGas;
    TempPsinkIOCompactMix.AccretionMassStar = Psk->AccretionMassStar;
    TempPsinkIOCompactMix.AccretionMassToBH = Psk->AccretionMassToBH;

#ifdef USE_CELIB
    TempPsinkIOCompactMix.Elements[0] = Psk->Elements[0];
    TempPsinkIOCompactMix.Elements[1] = Psk->Elements[1];
    TempPsinkIOCompactMix.Elements[2] = Psk->Elements[2];
    TempPsinkIOCompactMix.Elements[3] = Psk->Elements[3];
    TempPsinkIOCompactMix.Elements[4] = Psk->Elements[4];
    TempPsinkIOCompactMix.Elements[5] = Psk->Elements[5];
    TempPsinkIOCompactMix.Elements[6] = Psk->Elements[6];
    TempPsinkIOCompactMix.Elements[7] = Psk->Elements[7];
    TempPsinkIOCompactMix.Elements[8] = Psk->Elements[8];
    TempPsinkIOCompactMix.Elements[9] = Psk->Elements[9];
    TempPsinkIOCompactMix.Elements[10] = Psk->Elements[10];
    TempPsinkIOCompactMix.Elements[11] = Psk->Elements[11];
    TempPsinkIOCompactMix.Elements[12] = Psk->Elements[12];
    TempPsinkIOCompactMix.Elements[13] = Psk->Elements[13];
    TempPsinkIOCompactMix.Elements[14] = Psk->Elements[14];
    TempPsinkIOCompactMix.Elements[15] = Psk->Elements[15];
    TempPsinkIOCompactMix.Elements[16] = Psk->Elements[16];
    TempPsinkIOCompactMix.Elements[17] = Psk->Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TempPsinkIOCompactMix.SplitGeneration = Psk->SplitGeneration;
    TempPsinkIOCompactMix.SplitNthChildren = Psk->SplitNthChildren;
    TempPsinkIOCompactMix.GrandParentGlobalID = Psk->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPsinkIOCompactMix;

}

StructPsink CopyTemporalStructureCompactMixToPsinkCompactMix(struct StructPsinkIOCompactMix PsinkIOCompactMix){

	StructPsink TemporalPsink;
	memset(&TemporalPsink,0,sizeof(StructPsink));

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TemporalPsink.ParentGlobalID = PsinkIOCompactMix.ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
    TemporalPsink.dt_localmin = PsinkIOCompactMix.dt_localmin;

    TemporalPsink.NumberofAbsorbedParticles = PsinkIOCompactMix.NumberofAbsorbedParticles;
    TemporalPsink.AccretionRadius = PsinkIOCompactMix.AccretionRadius;
    TemporalPsink.MergingDistance = PsinkIOCompactMix.MergingDistance;
    TemporalPsink.AM[0] = PsinkIOCompactMix.AM[0];
    TemporalPsink.AM[1] = PsinkIOCompactMix.AM[1];
    TemporalPsink.AM[2] = PsinkIOCompactMix.AM[2];
    TemporalPsink.FormationTime = PsinkIOCompactMix.FormationTime;
    TemporalPsink.Z = PsinkIOCompactMix.Z;
    TemporalPsink.ZII = PsinkIOCompactMix.ZII;
    TemporalPsink.ZIa = PsinkIOCompactMix.ZIa;
    TemporalPsink.AccretionMass = PsinkIOCompactMix.AccretionMass;
    TemporalPsink.AccretionMassGas = PsinkIOCompactMix.AccretionMassGas;
    TemporalPsink.AccretionMassStar = PsinkIOCompactMix.AccretionMassStar;
    TemporalPsink.AccretionMassToBH = PsinkIOCompactMix.AccretionMassToBH;

#ifdef USE_CELIB
    TemporalPsink.Elements[0] = PsinkIOCompactMix.Elements[0];
    TemporalPsink.Elements[1] = PsinkIOCompactMix.Elements[1];
    TemporalPsink.Elements[2] = PsinkIOCompactMix.Elements[2];
    TemporalPsink.Elements[3] = PsinkIOCompactMix.Elements[3];
    TemporalPsink.Elements[4] = PsinkIOCompactMix.Elements[4];
    TemporalPsink.Elements[5] = PsinkIOCompactMix.Elements[5];
    TemporalPsink.Elements[6] = PsinkIOCompactMix.Elements[6];
    TemporalPsink.Elements[7] = PsinkIOCompactMix.Elements[7];
    TemporalPsink.Elements[8] = PsinkIOCompactMix.Elements[8];
    TemporalPsink.Elements[9] = PsinkIOCompactMix.Elements[9];
    TemporalPsink.Elements[10] = PsinkIOCompactMix.Elements[10];
    TemporalPsink.Elements[11] = PsinkIOCompactMix.Elements[11];
    TemporalPsink.Elements[12] = PsinkIOCompactMix.Elements[12];
    TemporalPsink.Elements[13] = PsinkIOCompactMix.Elements[13];
    TemporalPsink.Elements[14] = PsinkIOCompactMix.Elements[14];
    TemporalPsink.Elements[15] = PsinkIOCompactMix.Elements[15];
    TemporalPsink.Elements[16] = PsinkIOCompactMix.Elements[16];
    TemporalPsink.Elements[17] = PsinkIOCompactMix.Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPsink.SplitGeneration = PsinkIOCompactMix.SplitGeneration;
    TemporalPsink.SplitNthChildren = PsinkIOCompactMix.SplitNthChildren;
    TemporalPsink.GrandParentGlobalID = PsinkIOCompactMix.GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TemporalPsink;

}

struct StructPsinkIOLean CopyPsinkToTemporalStructureLean(const int index){

	struct StructPsinkIOLean TempPsinkIOLean;

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPsinkIOLean.ParentGlobalID = Psink[index]->ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
// omit double  dt_localmin;    // <LEAN> The local minimum time-steps.

    TempPsinkIOLean.NumberofAbsorbedParticles = Psink[index]->NumberofAbsorbedParticles;
    TempPsinkIOLean.AccretionRadius = (float)Psink[index]->AccretionRadius;
    TempPsinkIOLean.MergingDistance = (float)Psink[index]->MergingDistance;
    TempPsinkIOLean.AM[0] = (float)Psink[index]->AM[0];
    TempPsinkIOLean.AM[1] = (float)Psink[index]->AM[1];
    TempPsinkIOLean.AM[2] = (float)Psink[index]->AM[2];
    TempPsinkIOLean.FormationTime = (float)Psink[index]->FormationTime;
    TempPsinkIOLean.Z = (float)Psink[index]->Z;
    TempPsinkIOLean.ZII = (float)Psink[index]->ZII;
    TempPsinkIOLean.ZIa = (float)Psink[index]->ZIa;
    TempPsinkIOLean.AccretionMass = (float)Psink[index]->AccretionMass;
    TempPsinkIOLean.AccretionMassGas = (float)Psink[index]->AccretionMassGas;
    TempPsinkIOLean.AccretionMassStar = (float)Psink[index]->AccretionMassStar;
    TempPsinkIOLean.AccretionMassToBH = (float)Psink[index]->AccretionMassToBH;

#ifdef USE_CELIB
    TempPsinkIOLean.Elements[0] = (float)Psink[index]->Elements[0];
    TempPsinkIOLean.Elements[1] = (float)Psink[index]->Elements[1];
    TempPsinkIOLean.Elements[2] = (float)Psink[index]->Elements[2];
    TempPsinkIOLean.Elements[3] = (float)Psink[index]->Elements[3];
    TempPsinkIOLean.Elements[4] = (float)Psink[index]->Elements[4];
    TempPsinkIOLean.Elements[5] = (float)Psink[index]->Elements[5];
    TempPsinkIOLean.Elements[6] = (float)Psink[index]->Elements[6];
    TempPsinkIOLean.Elements[7] = (float)Psink[index]->Elements[7];
    TempPsinkIOLean.Elements[8] = (float)Psink[index]->Elements[8];
    TempPsinkIOLean.Elements[9] = (float)Psink[index]->Elements[9];
    TempPsinkIOLean.Elements[10] = (float)Psink[index]->Elements[10];
    TempPsinkIOLean.Elements[11] = (float)Psink[index]->Elements[11];
    TempPsinkIOLean.Elements[12] = (float)Psink[index]->Elements[12];
    TempPsinkIOLean.Elements[13] = (float)Psink[index]->Elements[13];
    TempPsinkIOLean.Elements[14] = (float)Psink[index]->Elements[14];
    TempPsinkIOLean.Elements[15] = (float)Psink[index]->Elements[15];
    TempPsinkIOLean.Elements[16] = (float)Psink[index]->Elements[16];
    TempPsinkIOLean.Elements[17] = (float)Psink[index]->Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TempPsinkIOLean.SplitGeneration = Psink[index]->SplitGeneration;
    TempPsinkIOLean.SplitNthChildren = Psink[index]->SplitNthChildren;
    TempPsinkIOLean.GrandParentGlobalID = Psink[index]->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPsinkIOLean;

}

struct StructPsinkIOLean CopyPsinkToTemporalStructureLeanElement(StructPsinkptr const Psk){

	struct StructPsinkIOLean TempPsinkIOLean;

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TempPsinkIOLean.ParentGlobalID = Psk->ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
// omit double  dt_localmin;    // <LEAN> The local minimum time-steps.

    TempPsinkIOLean.NumberofAbsorbedParticles = Psk->NumberofAbsorbedParticles;
    TempPsinkIOLean.AccretionRadius = Psk->AccretionRadius;
    TempPsinkIOLean.MergingDistance = Psk->MergingDistance;
    TempPsinkIOLean.AM[0] = Psk->AM[0];
    TempPsinkIOLean.AM[1] = Psk->AM[1];
    TempPsinkIOLean.AM[2] = Psk->AM[2];
    TempPsinkIOLean.FormationTime = Psk->FormationTime;
    TempPsinkIOLean.Z = Psk->Z;
    TempPsinkIOLean.ZII = Psk->ZII;
    TempPsinkIOLean.ZIa = Psk->ZIa;
    TempPsinkIOLean.AccretionMass = Psk->AccretionMass;
    TempPsinkIOLean.AccretionMassGas = Psk->AccretionMassGas;
    TempPsinkIOLean.AccretionMassStar = Psk->AccretionMassStar;
    TempPsinkIOLean.AccretionMassToBH = Psk->AccretionMassToBH;

#ifdef USE_CELIB
    TempPsinkIOLean.Elements[0] = Psk->Elements[0];
    TempPsinkIOLean.Elements[1] = Psk->Elements[1];
    TempPsinkIOLean.Elements[2] = Psk->Elements[2];
    TempPsinkIOLean.Elements[3] = Psk->Elements[3];
    TempPsinkIOLean.Elements[4] = Psk->Elements[4];
    TempPsinkIOLean.Elements[5] = Psk->Elements[5];
    TempPsinkIOLean.Elements[6] = Psk->Elements[6];
    TempPsinkIOLean.Elements[7] = Psk->Elements[7];
    TempPsinkIOLean.Elements[8] = Psk->Elements[8];
    TempPsinkIOLean.Elements[9] = Psk->Elements[9];
    TempPsinkIOLean.Elements[10] = Psk->Elements[10];
    TempPsinkIOLean.Elements[11] = Psk->Elements[11];
    TempPsinkIOLean.Elements[12] = Psk->Elements[12];
    TempPsinkIOLean.Elements[13] = Psk->Elements[13];
    TempPsinkIOLean.Elements[14] = Psk->Elements[14];
    TempPsinkIOLean.Elements[15] = Psk->Elements[15];
    TempPsinkIOLean.Elements[16] = Psk->Elements[16];
    TempPsinkIOLean.Elements[17] = Psk->Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TempPsinkIOLean.SplitGeneration = Psk->SplitGeneration;
    TempPsinkIOLean.SplitNthChildren = Psk->SplitNthChildren;
    TempPsinkIOLean.GrandParentGlobalID = Psk->GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TempPsinkIOLean;

}

StructPsink CopyTemporalStructureLeanToPsinkLean(struct StructPsinkIOLean PsinkIOLean){

	StructPsink TemporalPsink;
	memset(&TemporalPsink,0,sizeof(StructPsink));

// omit StructPsinkptr  Next;   // <TMP> <LEAN>
// omit StructPbodyptr  Body;   // <TMP> <LEAN>

// omit bool    Use;            // <TMP> <LEAN>
    TemporalPsink.ParentGlobalID = PsinkIOLean.ParentGlobalID;

// omit double  PosP[3];        // <TMP> <LEAN>
// omit double  VelP[3];        // <TMP> <LEAN>
// omit double  dt_localmin;    // <LEAN> The local minimum time-steps.

    TemporalPsink.NumberofAbsorbedParticles = PsinkIOLean.NumberofAbsorbedParticles;
    TemporalPsink.AccretionRadius = PsinkIOLean.AccretionRadius;
    TemporalPsink.MergingDistance = PsinkIOLean.MergingDistance;
    TemporalPsink.AM[0] = PsinkIOLean.AM[0];
    TemporalPsink.AM[1] = PsinkIOLean.AM[1];
    TemporalPsink.AM[2] = PsinkIOLean.AM[2];
    TemporalPsink.FormationTime = PsinkIOLean.FormationTime;
    TemporalPsink.Z = PsinkIOLean.Z;
    TemporalPsink.ZII = PsinkIOLean.ZII;
    TemporalPsink.ZIa = PsinkIOLean.ZIa;
    TemporalPsink.AccretionMass = PsinkIOLean.AccretionMass;
    TemporalPsink.AccretionMassGas = PsinkIOLean.AccretionMassGas;
    TemporalPsink.AccretionMassStar = PsinkIOLean.AccretionMassStar;
    TemporalPsink.AccretionMassToBH = PsinkIOLean.AccretionMassToBH;

#ifdef USE_CELIB
    TemporalPsink.Elements[0] = PsinkIOLean.Elements[0];
    TemporalPsink.Elements[1] = PsinkIOLean.Elements[1];
    TemporalPsink.Elements[2] = PsinkIOLean.Elements[2];
    TemporalPsink.Elements[3] = PsinkIOLean.Elements[3];
    TemporalPsink.Elements[4] = PsinkIOLean.Elements[4];
    TemporalPsink.Elements[5] = PsinkIOLean.Elements[5];
    TemporalPsink.Elements[6] = PsinkIOLean.Elements[6];
    TemporalPsink.Elements[7] = PsinkIOLean.Elements[7];
    TemporalPsink.Elements[8] = PsinkIOLean.Elements[8];
    TemporalPsink.Elements[9] = PsinkIOLean.Elements[9];
    TemporalPsink.Elements[10] = PsinkIOLean.Elements[10];
    TemporalPsink.Elements[11] = PsinkIOLean.Elements[11];
    TemporalPsink.Elements[12] = PsinkIOLean.Elements[12];
    TemporalPsink.Elements[13] = PsinkIOLean.Elements[13];
    TemporalPsink.Elements[14] = PsinkIOLean.Elements[14];
    TemporalPsink.Elements[15] = PsinkIOLean.Elements[15];
    TemporalPsink.Elements[16] = PsinkIOLean.Elements[16];
    TemporalPsink.Elements[17] = PsinkIOLean.Elements[17];
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    TemporalPsink.SplitGeneration = PsinkIOLean.SplitGeneration;
    TemporalPsink.SplitNthChildren = PsinkIOLean.SplitNthChildren;
    TemporalPsink.GrandParentGlobalID = PsinkIOLean.GrandParentGlobalID;
#endif // USE_PARTICLE_SPLITTING //}

	return TemporalPsink;

}


