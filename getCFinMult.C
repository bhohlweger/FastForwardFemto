void getCFinMult(const char* filename) {
  TFile *input=TFile::Open(filename);
  TString FolderName = "MBResults";
  TH2F* CentVsRefMult=nullptr;
  TDirectoryFile *dir=(TDirectoryFile*)(input->FindObjectAny(FolderName.Data()));
  TString centOrMult = "Mult";
  std::vector<int> binLimits;
  if (centOrMult == "Mult") {
    binLimits.push_back(-1);
    binLimits.push_back(27);
    binLimits.push_back(54);
    binLimits.push_back(200);
  }
  const int nMultBins = 3;
  TH2F* SEIntPart;
  TH2F* SEIntAPart;
  TH2F* MEIntPart;
  TH2F* MEIntAPart;
  TH1F* SEPartBin[nMultBins];
  TH1F* SEAPartBin[nMultBins];
  TH1F* MEPartBin[nMultBins];
  TH1F* MEAPartBin[nMultBins];

  if (dir) {
    TList *tmp=dir->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *OutList;
    dir->GetObject(name,OutList);
    if (!OutList) {
      std::cout << "No QAResult\n";
    } else {
      TList* CFPartList = (TList*)OutList->FindObject("Particle0_Particle0");
      SEIntPart =
          (TH2F*)CFPartList->FindObject(Form("SE%sDist_Particle0_Particle0",
                                             centOrMult.Data()));
      if (!SEIntPart) {
        std::cout << "No SEIntPart found \n";
        return;
      }

      MEIntPart =
          (TH2F*)CFPartList->FindObject(Form("ME%sDist_Particle0_Particle0",
                                             centOrMult.Data()));
      if (!MEIntPart) {
        std::cout << "No MEIntPart found \n";
        return;
      }
    }
  }

}
