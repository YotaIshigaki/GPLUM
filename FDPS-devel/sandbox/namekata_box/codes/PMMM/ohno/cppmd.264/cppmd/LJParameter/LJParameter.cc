
#include "LJParameter.h"

using namespace std;

LJParameter::LJParameter(const string _name, const double _sigma, const double _epsilon)
  : name(_name), sigma(_sigma), epsilon(_epsilon)
{
  setLocalParm();
}
void LJParameter::setLocalParm()
{
  potFirst =  4.0*epsilon*pow(sigma,12);
  potSecond = 4.0*epsilon*pow(sigma,6);
  forceFirst =  48.0*epsilon*pow(sigma,12);
  forceSecond = 24.0*epsilon*pow(sigma,6);
}
namespace LJParameterStuff {
  const LJParameter mixParameter(const LJParameter& p1, const LJParameter& p2){
    string name(p1.getName());
    name += "-";
    name += p2.getName();
    return LJParameter(name,
                       (p1.getSigma()+p2.getSigma())*.5,
                       sqrt(p1.getEpsilon()*p2.getEpsilon()));
  }
}

