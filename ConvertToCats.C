#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"

std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-3, kYellow-7};
std::vector<int> fColors     = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2};
std::vector<int> fMarkers    = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette=false, bool title=false)
{
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
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

  const int NRGBs = 6;
  Double_t stops[NRGBs];
  for(int i=0; i<NRGBs; ++i) stops[i] = float(i)/(NRGBs-1);

  Double_t red[NRGBs]   = { 1.,  29./255., 25./255., 27./255., 32./255., 24./255.};
  Double_t green[NRGBs] = { 1., 221./255., 160./255., 113./255., 74./255., 37./255.};
  Double_t blue[NRGBs] = {  1., 221./255., 184./255., 154./255., 129./255., 98./255.};
  TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHisto(TH1 *histo, int marker, int color)
{
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerSize(1.5);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleGraph(TGraph *histo, int marker, int color)
{
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerSize(1.5);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* Calculate_CF(TH1F* histRE_relK,TH1F* histME_relK, TString CFname,Double_t normleft,Double_t normright, const char* folder, float spinningDepth)
{
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  if(strcmp(folder, "") == 0) {
    Double_t norm_relK = histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) / histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright));
    std::cout << histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) << '\t' <<
        histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright)) << '\t' <<
        norm_relK << std::endl;
    Hist_CF->Divide(histRE_relK,histME_relK,1,norm_relK);
  }
  else {
    histME_relK->Scale(1.f/spinningDepth);
    Hist_CF->Divide(histRE_relK,histME_relK);
  }

  return Hist_CF;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* add_CF(TH1F* hist_CF1, TH1F* hist_CF2,TString HistName)
{
  //Calculate CFs with error weighting
  TH1F* hist_CF_sum = (TH1F*)hist_CF1->Clone(HistName.Data());

  Int_t NBins = hist_CF_sum->GetNbinsX();

  for(Int_t i=0;i<NBins;i++)
  {
    double CF1_val = hist_CF1->GetBinContent(i+1);
    double CF1_err = hist_CF1->GetBinError(i+1);
    double CF2_val = hist_CF2->GetBinContent(i+1);
    double CF2_err = hist_CF2->GetBinError(i+1);

    //average for bin i:
    if(CF1_val != 0. && CF2_val != 0.)
    {
      double CF1_err_weight = 1./TMath::Power(CF1_err,2.);
      double CF2_err_weight = 1./TMath::Power(CF2_err,2.);

      double CF_sum_average = (CF1_err_weight*CF1_val + CF2_err_weight*CF2_val)/(CF1_err_weight+CF2_err_weight);
      double CF_sum_err = 1./TMath::Sqrt(CF1_err_weight+CF2_err_weight);

      hist_CF_sum->SetBinContent(i+1,CF_sum_average);
      hist_CF_sum->SetBinError(i+1,CF_sum_err);
    }
    else if(CF1_val == 0. && CF2_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF2_val);
      hist_CF_sum->SetBinError(i+1,CF2_err);
    }
    else if(CF2_val == 0 && CF1_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF1_val);
      hist_CF_sum->SetBinError(i+1,CF1_err);
    }

  }
  return hist_CF_sum;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *DrawSystematicError(TH1F* histexp,TH1F *histerr,double errorwidth)
{
  //Input are the experimental histogram with statistical errors only and a histogram containing the errors

  const int histbins = histexp->GetNbinsX();

  TGraphErrors *ge_SysError_C2 = new TGraphErrors();

  for(int i=0;i<histbins;i++)
  {
    if(histexp->GetBinCenter(i+1) > 0.2) continue;
    ge_SysError_C2->SetPoint(i, histexp->GetBinCenter(i+1), histexp->GetBinContent(i+1));
    ge_SysError_C2->SetPointError(i, errorwidth, histerr->GetBinContent(i+1));
  }

  return ge_SysError_C2;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *FemtoModelFitBands(TGraph *grMedian1, TGraph *grLower, TGraph *grUpper)
{
  TGraphErrors *grFemtoModel = new TGraphErrors();
  grFemtoModel->SetName(grMedian1->GetName());
  double x, yM1, yLo, yUp;
  int count = 0;
  for(int i=0; i< grMedian1->GetN(); ++i) {
    grMedian1->GetPoint(i, x, yM1);
    grLower->GetPoint(i, x, yLo);
    grUpper->GetPoint(i, x, yUp);
    std::vector<float> yAll;
    yAll.push_back(yM1);
    yAll.push_back(yLo);
    yAll.push_back(yUp);
    std::sort(yAll.begin(), yAll.end());
    grFemtoModel->SetPoint(count, x/1000.f, (yAll[2]+yAll[0])/2.f);
    grFemtoModel->SetPointError(count++, 0, (yAll[2]+yAll[0])/2.f - yAll[0]);
  }
  return grFemtoModel;
}


TH1F* ConvertToOtherUnit(TH1F* HistCF,int Scale,TString name) {
  int nBins = HistCF->GetNbinsX();
  double kMin=HistCF->GetXaxis()->GetXmin();
  double kMax=HistCF->GetXaxis()->GetXmax();
  TH1F* HistScaled = new TH1F(name.Data(),name.Data(),nBins,kMin*Scale,kMax*Scale);
  for (int iBin = 1; iBin <= nBins; ++iBin ) {
    HistScaled->SetBinContent(iBin,HistCF->GetBinContent(iBin));
    HistScaled->SetBinError(iBin,HistCF->GetBinError(iBin));
  }
  return HistScaled;
}



void ConvertToCats(const char *fileExp, TString partPair) {
  TFile* _file0=TFile::Open(fileExp,"READ");
  const float normleft = 200;
  const float normright = 400;

  TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny("MBResults"));
  TList *Results;
  dirResults->GetObject("MBResults",Results);
  TList* tmpFolder;
  if (partPair.Contains("pp")) {
    TH1F *ppApApCF[3];
    tmpFolder=(TList*)Results->FindObject("Particle0_Particle0");
    TH1F* ppSE = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle0")),
        1000,"hPartSe");
    if (!ppSE ) {
      std::cout << "No pp SE Hist \n";
    }
    TH1F* ppME = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle0")),
        1000,"hPartMe");
    ppApApCF[0] = Calculate_CF(ppSE,ppME,"hCkPartNorm",normleft,normright, "", 0);
    if (!ppME ) {
      std::cout << "No pp ME Hist \n";
    }

    tmpFolder=(TList*)Results->FindObject("Particle1_Particle1");
    TH1F* ApApSE = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle1")),
        1000,"hAntiPartSe");
    if (!ApApSE ) {
      std::cout << "No ApAp SE Hist \n";
    }
    TH1F* ApApME = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle1")),
        1000,"hAntiPartMe");
    if (!ApApME ) {
      std::cout << "No ApAp ME Hist \n";
    }
    ppApApCF[1] = Calculate_CF(ApApSE,ApApME,"hCkAntiPartNorm",normleft,normright, "", 0);
    ppApApCF[2] = add_CF(ppApApCF[0],ppApApCF[1],"hCkTotNormWeight");

    TFile *output=
        TFile::Open("AnalysisResults_AndiBernieSystME_CkProtonProton_0.root","RECREATE");
    ppSE->Write();
    ppME->Write();
    ApApSE->Write();
    ApApME->Write();
    ppApApCF[0]->Write();
    ppApApCF[1]->Write();
    ppApApCF[2]->Write();

    output->Close();
  }

  if (partPair.Contains("pL")) {
    TH1F *pLApALCF[3];
    tmpFolder=(TList*)Results->FindObject("Particle0_Particle2");
    TH1F* pLSE = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle2")),
        1000,"hPartSe");
    TH1F* pLME = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle2")),
        1000,"hPartMe");
    pLApALCF[0] = Calculate_CF(pLSE,pLME,"hCkPartNorm",normleft,normright, "", 0);

    tmpFolder=(TList*)Results->FindObject("Particle1_Particle3");
    TH1F* ApALSE = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle3")),
        1000,"hAntiPartSe");
    TH1F* ApALME = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle3")),
        1000,"hAntiPartMe");
    pLApALCF[1] = Calculate_CF(ApALSE,ApALME,"hCkAntiPartNorm",normleft,normright, "", 0);
    pLApALCF[2] = add_CF(pLApALCF[0],pLApALCF[1],"hCkTotNormWeight");

    TFile *output=
        TFile::Open("AnalysisResults_AndiBernieSystME_CkProtonLambda_0.root","RECREATE");
    pLSE->Write();
    pLME->Write();
    ApALSE->Write();
    ApALME->Write();
    pLApALCF[0]->Write();
    pLApALCF[1]->Write();
    pLApALCF[2]->Write();

    output->Close();
  }

  if (partPair.Contains("LL")) {
    TH1F *LLALALCF[3];
    tmpFolder=(TList*)Results->FindObject("Particle2_Particle2");
    TH1F* LLSE = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("SEDist_Particle2_Particle2")),
        1000,"hPartSe");
    TH1F* LLME = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("MEDist_Particle2_Particle2")) ,
        1000,"hPartMe");
    LLALALCF[0] = Calculate_CF(LLSE,LLME,"hCkPartNorm",normleft,normright, "", 0);

    tmpFolder=(TList*)Results->FindObject("Particle3_Particle3");
    TH1F* ALALSE = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("SEDist_Particle3_Particle3")),
        1000,"hAntiPartSe");
    TH1F* ALALME = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("MEDist_Particle3_Particle3")),
        1000,"hAntiPartMe");
    LLALALCF[1] = Calculate_CF(ALALSE,ALALME,"hCkAntiPartNorm",normleft,normright, "", 0);
    LLALALCF[2] = add_CF(LLALALCF[0],LLALALCF[1],"hCkTotNormWeight");

    TFile *output=
        TFile::Open("AnalysisResults_AndiBernieSystME_CkLambdaLambda_0.root","RECREATE");
    LLSE->Write();
    LLME->Write();
    ALALSE->Write();
    ALALME->Write();
    LLALALCF[0]->Write();
    LLALALCF[1]->Write();
    LLALALCF[2]->Write();

    output->Close();
  }

  if (partPair.Contains("pXi")) {
    TH1F *pXiApAXiCF[3];
    tmpFolder=(TList*)Results->FindObject("Particle0_Particle4");
    TH1F* pXiSE = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle4")),
        1000,"hPartSe");
    TH1F* pXiME = ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle4")),
        1000,"hPartMe");
    pXiApAXiCF[0] = Calculate_CF(pXiSE,pXiME,"hCkPartNorm",normleft,normright, "", 0);

    tmpFolder=(TList*)Results->FindObject("Particle1_Particle5");
    TH1F *ApAXiSE= ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle5")),
        1000,"hAntiPartSe");
    TH1F *ApAXiME=ConvertToOtherUnit(
        ((TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle5")),
        1000,"hAntiPartMe");
    pXiApAXiCF[1] = Calculate_CF(ApAXiSE,ApAXiME,"hCkAntiPartNorm",normleft,normright, "", 0);
    pXiApAXiCF[2] = add_CF(pXiApAXiCF[0],pXiApAXiCF[1],"hCkTotNormWeight");

    TFile *output=
        TFile::Open("AnalysisResults_AndiBernieSystME_CkProtonXim_0.root","RECREATE");
    pXiSE->Write();
    pXiME->Write();
    ApAXiSE->Write();
    ApAXiME->Write();
    pXiApAXiCF[0]->Write();
    pXiApAXiCF[1]->Write();
    pXiApAXiCF[2]->Write();

    output->Close();
  }


  return;
}
