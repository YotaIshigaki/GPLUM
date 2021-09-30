#ifndef SETPAIRLIST_H
#define SETPAIRLIST_H
#include <vector>
#include "ParticleInfo.h"

//! index ranges of subsets for each forcetype @n
/*! |--number_of_lj--|--number_of_ljcoulomb--|--number_of_coulomb--| @n
    |  lj particles  |  ljcoulomb particles  |  coulomb particles  | @n
    |0|1|2|.|i |.|N-1|0|1|2|....|i|......|N-1|0|1|2|....|i|....|N-1| @n
            |  | @n
    lj[i].begin lj[i].end         @n
    N : number_of_sets @n
*/
typedef struct {
  int number_of_sets;
  int number_of_lj;
  int number_of_ljcoulomb;
  int number_of_coulomb;
  std::vector<ParticleRange> lj;
  std::vector<ParticleRange> ljcoulomb;
  std::vector<ParticleRange> coulomb;
} SetsRange;

typedef std::vector<ParticleRange> PairRangeList;

void clearPairRangeList(PairRangeList& prlist);

void clearPairRangeListArray(std::vector<PairRangeList> prlistarray);


inline bool checkrange(ParticleRange range){
  return (range.end>range.begin);
}

inline void addrange(const ParticleRange& ir, const PairList& jsetlist, 
                     const std::vector<ParticleRange>& jrlist, PairRangeList& jrange)
{
  if(checkrange(ir)){   // check iset  exists
    for(PairList::const_iterator it = jsetlist.begin();
        it!=jsetlist.end(); ++it){
      int jset = *it;
      if(checkrange(jrlist[jset])){ // check jset exists
        if(jrange.empty()){
          jrange.push_back(jrlist[jset]);
        }else{
          if(jrlist[jset].begin==jrange[jrange.size()-1].end){
            jrange[jrange.size()-1].end = jrlist[jset].end;
          }else{
            jrange.push_back(jrlist[jset]);
          }
        }
      }
    }
  }
}

void generateParticleRangeSets(const SetsRange& i_sets, const SetsRange& j_sets,
                               const std::vector<PairList>& jsetlist_for_isets,
                               std::vector<PairRangeList>& jranges_for_ljisets,
                               std::vector<PairRangeList>& jranges_for_ljcoulombisets,
                               std::vector<PairRangeList>& jranges_for_coulombisets);

#endif
