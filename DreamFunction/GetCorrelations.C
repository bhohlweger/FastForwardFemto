#include "TROOT.h"

void GetCorrelations(const char* filename, const char* prefix) {
  ReadDreamFile* DreamFile=new ReadDreamFile(6,6);
  DreamFile->SetAnalysisFile(filename,prefix);

  DreamCF* pp=new DreamCF();
//  pp->SetPairOne(DreamFile->GetPairDistributions(0,0,"0"));
//  pp->SetPairTwo(DreamFile->GetPairDistributions(1,1,"0"));
}
