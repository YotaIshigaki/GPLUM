#include <iostream>
#include "SetPairList.h"

using namespace std;

void nums_to_sets(int nums[][3], int size, SetsRange& sets)
{
  int pn;
  sets.number_of_sets = size;
  sets.number_of_lj = 0;
  sets.number_of_ljcoulomb = 0;
  sets.number_of_coulomb = 0;
  pn = 0;
  for(int n=0;n<sets.number_of_sets;n++){
    ParticleRange pr = {pn,pn};
    sets.number_of_lj += nums[n][0];
    pn += nums[n][0];
    pr.end = pn;
    sets.lj.push_back(pr);
  }
  for(int n=0;n<sets.number_of_sets;n++){
    ParticleRange pr = {pn,pn};
    sets.number_of_ljcoulomb += nums[n][1];
    pn += nums[n][1];
    pr.end = pn;
    sets.ljcoulomb.push_back(pr);
  }
  for(int n=0;n<sets.number_of_sets;n++){
    ParticleRange pr = {pn,pn};
    sets.number_of_coulomb += nums[n][2];
    pn += nums[n][2];
    pr.end = pn;
    sets.coulomb.push_back(pr);
  }
}

int
main()
{
  int i_nums[][3] = {{1,12,3},{1,8,7}};
  int j_nums[][3] = {{1,10,5},{1,12,3},{1,8,7},{1,5,10}};

  int nicell=2, njcell=4;
  std::vector<PairList> jsetlists;
  SetsRange isets;
  SetsRange jsets;
  std::vector<PairRangeList> jranges_for_ljisets;
  std::vector<PairRangeList> jranges_for_ljcoulombisets;
  std::vector<PairRangeList> jranges_for_coulombisets;
  PairList jl;

  int tmp[20];
  
  jl.push_back(0);
  jl.push_back(1);
  jl.push_back(3);
  jsetlists.push_back(jl);
  jl.clear();
  jl.push_back(1);
  jl.push_back(2);
  jsetlists.push_back(jl);
  
  nums_to_sets(i_nums,sizeof(i_nums)/sizeof(i_nums[0]),isets);
  nums_to_sets(j_nums,sizeof(j_nums)/sizeof(j_nums[0]),jsets);

  for(int n=0;n<jsets.lj.size();n++){
    cout << "LJ " << n << " " << jsets.lj[n].begin << " - " << jsets.lj[n].end << endl;
  }
  for(int n=0;n<jsets.ljcoulomb.size();n++){
    cout << "LJCoulomb " << n << " " << jsets.ljcoulomb[n].begin << " - " << jsets.ljcoulomb[n].end << endl;
  }
  for(int n=0;n<jsets.coulomb.size();n++){
    cout << "Coulomb " << n << " " << jsets.coulomb[n].begin << " - " << jsets.coulomb[n].end << endl;
  }

  generateParticleRangeSets(isets,jsets,jsetlists,
                            jranges_for_ljisets,
                            jranges_for_ljcoulombisets,
                            jranges_for_coulombisets);
  
  for(int n=0;n<jranges_for_ljisets.size();n++){
    cout << "i-set " << n;
    for(int j=0;j<jranges_for_ljisets[n].size();j++){
      cout << " " << jranges_for_ljisets[n][j].begin << "--" << jranges_for_ljisets[n][j].end ;
    }
    cout << endl;
  }
  for(int n=0;n<jranges_for_ljcoulombisets.size();n++){
    cout << "i-set " << n + isets.number_of_sets;
    for(int j=0;j<jranges_for_ljcoulombisets[n].size();j++){
      cout << " " << jranges_for_ljcoulombisets[n][j].begin << "--" << jranges_for_ljcoulombisets[n][j].end ;
    }
    cout << endl;
  }
  for(int n=0;n<jranges_for_coulombisets.size();n++){
    cout << "i-set " << n + isets.number_of_sets;
    for(int j=0;j<jranges_for_coulombisets[n].size();j++){
      cout << " " << jranges_for_coulombisets[n][j].begin << "--" << jranges_for_coulombisets[n][j].end ;
    }
    cout << endl;
  }
}
