#include "SetPairList.h"

void generateParticleRangeSets(const SetsRange& i_sets, const SetsRange& j_sets,
                               const std::vector<PairList>& jsetlist_for_isets,
                               std::vector<PairRangeList>& jranges_for_ljisets,
                               std::vector<PairRangeList>& jranges_for_ljcoulombisets,
                               std::vector<PairRangeList>& jranges_for_coulombisets)
{
  PairRangeList prlist;
  jranges_for_ljisets.clear();
  for(int iset=0;iset<i_sets.number_of_sets;iset++){
    prlist.clear();
    addrange(i_sets.lj[iset],jsetlist_for_isets[iset],j_sets.lj,
             prlist);
    addrange(i_sets.lj[iset],jsetlist_for_isets[iset],j_sets.ljcoulomb,
             prlist);
    jranges_for_ljisets.push_back(prlist);
  }
  for(int iset=0;iset<i_sets.number_of_sets;iset++){
    prlist.clear();
    addrange(i_sets.ljcoulomb[iset],jsetlist_for_isets[iset],j_sets.lj,
             prlist);
    jranges_for_ljisets.push_back(prlist);
  }
  jranges_for_ljcoulombisets.clear();
  for(int iset=0;iset<i_sets.number_of_sets;iset++){
    prlist.clear();
    addrange(i_sets.ljcoulomb[iset],jsetlist_for_isets[iset],j_sets.ljcoulomb,
             prlist);
    jranges_for_ljcoulombisets.push_back(prlist);
  }
  jranges_for_coulombisets.clear();
  for(int iset=0;iset<i_sets.number_of_sets;iset++){
    prlist.clear();
    addrange(i_sets.ljcoulomb[iset],jsetlist_for_isets[iset],j_sets.coulomb,
             prlist);
    jranges_for_coulombisets.push_back(prlist);
  }
  for(int iset=0;iset<i_sets.number_of_sets;iset++){
    prlist.clear();
    addrange(i_sets.coulomb[iset],jsetlist_for_isets[iset],j_sets.ljcoulomb,
             prlist);
    addrange(i_sets.coulomb[iset],jsetlist_for_isets[iset],j_sets.coulomb,
             prlist);
    jranges_for_coulombisets.push_back(prlist);
  }
}
