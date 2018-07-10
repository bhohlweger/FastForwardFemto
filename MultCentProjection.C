Color_t fFillColors[8] = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9,
    kOrange-9, kCyan-8, kYellow-7};
Color_t fColors[8]  = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1,
    kOrange-1, kCyan+2, kYellow+2};
Style_t fMarkers[10] = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare,
    kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette, bool title)
{
  gStyle->Reset("Plain");
  gStyle->SetCanvasPreferGL(1);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  auto ci = 1179;
  auto color = new TColor(ci, 1, 0, 0, " ", 0.);
  gStyle->SetStatColor(ci);
  gStyle->SetTitleColor(ci);
  gStyle->SetCanvasColor(ci);
  gStyle->SetPadColor(ci);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //  gStyle->SetPadBottomMargin(0.15);
  //  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHisto(TH1 *histo, int marker, int color)
{
  histo->GetXaxis()->SetLabelSize(0.06);
  histo->GetXaxis()->SetTitleSize(0.07);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.0);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.06);
  histo->GetYaxis()->SetTitleSize(0.07);
  histo->GetYaxis()->SetLabelOffset(0.02);
  histo->GetYaxis()->SetTitleOffset(1.0);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}


void MultCentProjection(const char *file) {
  SetStyle(false,true);
  TFile* ExpFile = TFile::Open(file);
  TString FolderName = "MBEvtCuts";
  TH2F* CentVsRefMult=nullptr;
  TDirectoryFile *dir=(TDirectoryFile*)(ExpFile->FindObjectAny(FolderName.Data()));
  std::vector<int> binLimits ={-1,27,54,200};
  const int nMultBins = 3;
  TH1F* CentDistPerMultBin[nMultBins];
  if (dir) {
    TList *tmp=dir->GetListOfKeys();
    TString name=tmp->At(0)->GetName();
    TList *OutList;
    dir->GetObject(name,OutList);
    if (!OutList) {
      std::cout << "No QAResult\n";
    } else {
      CentVsRefMult=(TH2F*)OutList->FindObject("CentvsRefMult");
//      CentVsRefMult->DrawCopy();
      if (!CentVsRefMult) {
        std::cout << "No Histogram found \n";
        return;
      }
    }
  }
  TCanvas *MultDist= new TCanvas("MultDist","MultDist",0,0,1500,1100);
  auto* leg= new TLegend(0.18, 0.7, 0.55, 0.85);
  for (int iMult = 0;iMult < nMultBins;++iMult) {
    TString HistName = Form("MultBin_%d",iMult+1);
    int firstBin=CentVsRefMult->GetYaxis()->FindBin(binLimits.at(iMult)+1);
    int lastBin=CentVsRefMult->GetYaxis()->FindBin(binLimits.at(iMult+1));
    CentDistPerMultBin[iMult]=(TH1F*)CentVsRefMult->ProjectionX(HistName.Data(),firstBin,lastBin);
    SetStyleHisto(CentDistPerMultBin[iMult],1,iMult+1);
    CentDistPerMultBin[iMult]->SetTitle(Form(";Centrality (%%);Counts"));
    MultDist->cd();
    if (iMult ==0 ) CentDistPerMultBin[iMult]->DrawCopy();
    else CentDistPerMultBin[iMult]->DrawCopy("SAME");
    leg->AddEntry(CentDistPerMultBin[iMult],
                  Form("Multiplicity from %d to %d",
                       binLimits.at(iMult)+1,binLimits.at(iMult+1)), "pe");
  }
  MultDist->cd();
  leg->Draw("SAME");
}
