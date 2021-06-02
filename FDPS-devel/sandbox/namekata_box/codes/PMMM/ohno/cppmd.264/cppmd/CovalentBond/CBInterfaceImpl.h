#ifndef CBINTERFACEIMPL_H
#define CBINTERFACEIMPL_H

#include <stdexcept>
#include <cstdio>
#include "CBInterface.h"
#include "ErrorPos.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace CBModule {
#ifdef ORG_BOND
template<>
CBInterface::CBBondIterator CBInterface::begin<CBInterface::Bond>
                                              (Context pContext)
{
  return pContext->pBond->begin();
}

template<>
CBInterface::CBBondIterator CBInterface::end<CBInterface::Bond>
                                              (Context pContext)
{
  return pContext->pBond->end();
}
#endif
#ifdef  ORG_ANGLE 
template<>
CBInterface::CBAngleIterator CBInterface::begin<CBInterface::Angle>
                                              (Context pContext)
{
  return pContext->pAngle->begin();
}

template<>
CBInterface::CBAngleIterator CBInterface::end<CBInterface::Angle>
                                              (Context pContext)
{
  return pContext->pAngle->end();
}
#endif
#ifdef ORG_TORSION
template<>
CBInterface::CBTorsionIterator CBInterface::begin<CBInterface::Torsion>
                                              (Context pContext)
{
  return pContext->pTorsion->begin();
}

template<>
CBInterface::CBTorsionIterator CBInterface::end<CBInterface::Torsion>
                                              (Context pContext)
{
  return pContext->pTorsion->end();
}

template<>
CBInterface::CBImproperIterator CBInterface::begin<CBInterface::Improper>
                                              (Context pContext)
{
  return pContext->pImproper->begin();
}

template<>
CBInterface::CBImproperIterator CBInterface::end<CBInterface::Improper>
                                              (Context pContext)
{
  return pContext->pImproper->end();
}
 
#endif
}



#endif
