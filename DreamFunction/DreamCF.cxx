/*
 * DreamCF.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#include "DreamCF.h"

DreamCF::DreamCF()
:fPair(nullptr)
,fPairShifted(nullptr)
,fPairLast(nullptr)
{
}

DreamCF::~DreamCF()
{
}

void DreamCF::ShiftForEmpty() {
  int SEkMin=fPair->GetSEDist()->FindFirstBinAbove(0);
  if (!fPairLast)
  {
    fPairLast=new DreamPair(fPair,"tmp");
  }

}
