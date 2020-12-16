#pragma once
#define EffSAVecSize (7)

#include "StellarFeedback.h"

void CalcEffectiveSurfaceAreaPrev(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticleOriginal);
void CalcEffectiveSurfaceAreaSum(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticleOriginal);
void CalcEffectiveSurfaceAreaVec(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticleOriginal);
void CalcEffectiveSurfaceAreaVecCheck(const int NExplosion, struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticleOriginal);
bool CheckExporFlagGatherScatter(const int Index, const int SendRank);
