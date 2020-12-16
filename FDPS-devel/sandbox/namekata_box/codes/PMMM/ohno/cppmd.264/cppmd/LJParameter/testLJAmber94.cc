#include <iostream>
#include <vector>
#include <boost/multi_array.hpp>
#include <boost/timer.hpp>
#include "LJAmber94.h"

using namespace std;

int
main(int argc, char **argv)
{
  boost::timer tm;
  double time;

  LJAmber94 ljamber94;

  int num_lj = ljamber94.num_lj;


  LJMixparameterArray ljma94;
  ljamber94.convertLJMixparameterArray(ljma94);

  vector<LJParameter> ljps;
  for(int i=0;i<ljamber94.num_lj;i++){
    /*
    ljps.push_back(LJParameter(ljamber94.ljamberparameters[i].name,
                               ljamber94.ljamberparameters[i].r*pow(2.0,5.0/6.0),
                               ljamber94.ljamberparameters[i].epsilon));
    */
    ljps.push_back(ljamber94.getLJParameter(i));
  }

  for(int i=0;i<ljamber94.num_lj;i++){
    cout << i;
    cout << " " << ljamber94.ljamberparameters[i].name;
    cout << " " << ljamber94.ljamberparameters[i].r;
    cout << " " << ljamber94.ljamberparameters[i].epsilon;
    cout << " : " << ljps[i].getName();
    cout << " " << ljps[i].getSigma();
    cout << " " << ljps[i].getEpsilon();
    cout << endl;
  }

  boost::multi_array<LJParameter,2> ljaa(boost::extents[num_lj][num_lj]);
  boost::multi_array<LJParameter,2> ljaad(boost::extents[num_lj][num_lj]);
  for(int i=0;i<num_lj;i++){
    for(int j=0;j<num_lj;j++){
      ljaa[i][j] = LJParameterStuff::mixParameter(ljps[i],ljps[j]);
      ljaad[i][j] = ljamber94.getmixParameter(i,j);
    }
  }

  for(int i=0;i<ljamber94.num_lj;i++){
    cout << i ;
    cout << " " << ljaa[0][i].getName();
    cout << " " << ljma94[0][i].potFirst-ljaa[0][i].getPotentialFirst();
    cout << " " << ljma94[0][i].potSecond-ljaa[0][i].getPotentialSecond();
    cout << " " << ljma94[0][i].forceFirst-ljaa[0][i].getForceFirst();
    cout << " " << ljma94[0][i].forceSecond-ljaa[0][i].getForceSecond();
    cout << " " << ljaad[0][i].getPotentialFirst()-ljaa[0][i].getPotentialFirst();
    cout << " " << ljaad[0][i].getPotentialSecond()-ljaa[0][i].getPotentialSecond();
    cout << " " << ljaad[0][i].getForceFirst()-ljaa[0][i].getForceFirst();
    cout << " " << ljaad[0][i].getForceSecond()-ljaa[0][i].getForceSecond();
    cout << endl;
  }

}
