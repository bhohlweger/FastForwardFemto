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

std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7};
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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* Calculate_CF(TH1F* histRE_relK,TH1F* histME_relK, TString CFname,Double_t normleft,Double_t normright, const char* folder, float spinningDepth)
{
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  if(strcmp(folder, "") == 0) {
    Double_t norm_relK = histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) / histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright));
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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraph *convertInGev(TGraph *gr) {
  TGraph *grOut = new TGraph();
  grOut->SetName(gr->GetName());
  double x,y;
  for(int i=0; i< gr->GetN(); ++i) {
    gr->GetPoint(i, x, y);
    grOut->SetPoint(i, x/1000.f, y);
  }
  return grOut;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *convertHistoInGev(TH1F *gr) {
  TGraphErrors *grOut = new TGraphErrors();
  grOut->SetName(gr->GetName());
  for(int i=0; i< gr->GetNbinsX(); ++i) {
    grOut->SetPoint(i, gr->GetBinCenter(i)/1000.f, gr->GetBinContent(i));
    grOut->SetPointError(i, 0, gr->GetBinError(i));
  }
  return grOut;
}

void Closure(TString expfile) {
  SetStyle(false,true);
  const float normleft = 0.2;
  const float normright = 0.4;
  const float spinningDepth = 10.f;

  const float right = 0.025;
  const float top = 0.025;
  int minBin = 1;
  int divisionBin = 25;
  int maxBin = 26;
    int RebinFac = 2;
     int RebinFacPP = 1;
  const char* prefix="MB";
  const char* addon="";
  TFile* _file0=TFile::Open(expfile);
  TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon),Results);


  TList* tmpFolder=(TList*)Results->FindObject("Particle0_Particle0");
  TH1F* histRE_relK_pp = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle0");
    histRE_relK_pp->Rebin(RebinFacPP);
  TH1F* histME_relK_pp = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle0");
    histME_relK_pp->Rebin(RebinFacPP);
//  TH2F* histRE_relK_ppMult = (TH1F*)tmpFolder->FindObject("SEMultDist_Particle0_Particle0");
//  TH1F *histRE_relK_LpBinMultLow = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_NonCent",minBin,divisionBin);
//  TH1F *histRE_relK_LpBinMultHigh = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_Cent",divisionBin+1,maxBin);
//
//  TH2F* histME_relK_ppMult = (TH1F*)tmpFolder->FindObject("MEMultDist_Particle0_Particle0");


  tmpFolder=(TList*)Results->FindObject("Particle1_Particle1");
  TH1F* histRE_relK_ApAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle1");
    histRE_relK_ApAp->Rebin(RebinFacPP);
  TH1F* histME_relK_ApAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle1");
    histME_relK_ApAp->Rebin(RebinFacPP);
//  TH2F* histRE_relK_ApApMult = (TH1F*)tmpFolder->FindObject("SEMultDist_Particle1_Particle1");
//  TH1F *histRE_relK_LpBinMultLow = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_NonCent",minBin,divisionBin);
//  TH1F *histRE_relK_LpBinMultHigh = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_Cent",divisionBin+1,maxBin);
//
//  TH2F* histME_relK_ApApMult = (TH1F*)tmpFolder->FindObject("MEMultDist_Particle1_Particle1");


  tmpFolder=(TList*)Results->FindObject("Particle0_Particle2");
  TH1F* histRE_relK_Lp = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle2");
  TH1F* histME_relK_Lp = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle2");
  histRE_relK_Lp->Rebin(RebinFac);
      histME_relK_Lp->Rebin(RebinFac);
//  TH2F* histRE_relK_LpMult = (TH1F*)tmpFolder->FindObject("SEMultDist_Particle0_Particle2");
//  TH1F *histRE_relK_LpBinMultLow = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_NonCent",minBin,divisionBin);
//  TH1F *histRE_relK_LpBinMultHigh = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_Cent",divisionBin+1,maxBin);
//
//  TH2F* histME_relK_LpMult = (TH1F*)tmpFolder->FindObject("MEMultDist_Particle0_Particle2");
//  TH1F *histME_relK_LpBinMultLow = (TH1F*)histME_relK_LpMult->ProjectionX("MELP_NonCent",minBin,divisionBin);
//  TH1F *histME_relK_LpBinMultHigh = (TH1F*)histME_relK_LpMult->ProjectionX("MELP_Cent",divisionBin+1,maxBin);


  tmpFolder=(TList*)Results->FindObject("Particle1_Particle3");
  TH1F* histRE_relK_ALAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle3");
  TH1F* histME_relK_ALAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle3");
    histRE_relK_ALAp->Rebin(RebinFac);
    histME_relK_ALAp->Rebin(RebinFac);
//  TH2F* histRE_relK_ALApMult = (TH1F*)tmpFolder->FindObject("SEMultDist_Particle1_Particle3");
//  TH1F *histRE_relK_LpBinMultLow = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_NonCent",minBin,divisionBin);
//  TH1F *histRE_relK_LpBinMultHigh = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_Cent",divisionBin+1,maxBin);
//
//  TH2F* histME_relK_ALApMult = (TH1F*)tmpFolder->FindObject("MEMultDist_Particle1_Particle3");


  tmpFolder=(TList*)Results->FindObject("Particle2_Particle2");
  TH1F* histRE_relK_LL = (TH1F*)tmpFolder->FindObject("SEDist_Particle2_Particle2");
  TH1F* histME_relK_LL = (TH1F*)tmpFolder->FindObject("MEDist_Particle2_Particle2");
    histRE_relK_LL->Rebin(RebinFac);
    histME_relK_LL->Rebin(RebinFac);
//  TH2F* histRE_relK_LLMult = (TH1F*)tmpFolder->FindObject("SEMultDist_Particle2_Particle2");
//  TH1F *histRE_relK_LpBinMultLow = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_NonCent",minBin,divisionBin);
//  TH1F *histRE_relK_LpBinMultHigh = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_Cent",divisionBin+1,maxBin);
//
//  TH2F* histME_relK_LLMult = (TH1F*)tmpFolder->FindObject("MEMultDist_Particle2_Particle2");


  tmpFolder=(TList*)Results->FindObject("Particle3_Particle3");
  TH1F* histRE_relK_ALAL = (TH1F*)tmpFolder->FindObject("SEDist_Particle3_Particle3");
  TH1F* histME_relK_ALAL = (TH1F*)tmpFolder->FindObject("MEDist_Particle3_Particle3");
    histRE_relK_ALAL->Rebin(RebinFac);
    histME_relK_ALAL->Rebin(RebinFac);
//  TH2F* histRE_relK_ALALMult = (TH1F*)tmpFolder->FindObject("SEMultDist_Particle3_Particle3");
//  TH1F *histRE_relK_LpBinMultLow = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_NonCent",minBin,divisionBin);
//  TH1F *histRE_relK_LpBinMultHigh = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_Cent",divisionBin+1,maxBin);
//
//  TH2F* histME_relK_ALALMult = (TH1F*)tmpFolder->FindObject("MEMultDist_Particle3_Particle3");


  tmpFolder=(TList*)Results->FindObject("Particle0_Particle4");
  TH1F* histRE_relK_Xip = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle4");
  TH1F* histME_relK_Xip = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle4");
    histRE_relK_Xip->Rebin(RebinFac);
    histME_relK_Xip->Rebin(RebinFac);
//  TH2F* histRE_relK_XipMult = (TH1F*)tmpFolder->FindObject("SEMultDist_Particle0_Particle4");
//  TH1F *histRE_relK_LpBinMultLow = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_NonCent",minBin,divisionBin);
//  TH1F *histRE_relK_LpBinMultHigh = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_Cent",divisionBin+1,maxBin);
//
//  TH2F* histME_relK_XipMult = (TH1F*)tmpFolder->FindObject("MEMultDist_Particle0_Particle4");


  tmpFolder=(TList*)Results->FindObject("Particle1_Particle5");
  TH1F* histRE_relK_AXiAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle5");
  TH1F* histME_relK_AXiAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle5");
    histRE_relK_AXiAp->Rebin(RebinFac);
    histME_relK_AXiAp->Rebin(RebinFac);
//  TH2F* histRE_relK_AXiApMult = (TH1F*)tmpFolder->FindObject("SEMultDist_Particle1_Particle5");
//  TH1F *histRE_relK_LpBinMultLow = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_NonCent",minBin,divisionBin);
//  TH1F *histRE_relK_LpBinMultHigh = (TH1F*)histRE_relK_LpMult->ProjectionX("SELP_Cent",divisionBin+1,maxBin);
//
//  TH2F* histME_relK_AXiApMult = (TH1F*)tmpFolder->FindObject("MEMultDist_Particle1_Particle5");

  TH1F *hist_CF_Lp_ALAp_exp[3];
  TH1F *hist_CF_LL_ALAL_exp[3];
  TH1F *hist_CF_LL_ALAL_expMult[2][3];
  TH1F *hist_CF_pp_ApAp_exp[3];
  TH1F *hist_CF_pp_ApAp_expMult[2][3];
  TH1F *hist_CF_pXi_ApAXi_exp[3];
  TH1F *hist_CF_pXi_ApAXi_expMult[2][3];

  hist_CF_Lp_ALAp_exp[0] = Calculate_CF(histRE_relK_Lp,histME_relK_Lp,"hist_CF_Lp_exp",normleft,normright, addon, spinningDepth);
  hist_CF_Lp_ALAp_exp[1] = Calculate_CF(histRE_relK_ALAp,histME_relK_ALAp,"hist_CF_ALAp_exp",normleft,normright, addon, spinningDepth);
  hist_CF_Lp_ALAp_exp[2] = add_CF(hist_CF_Lp_ALAp_exp[0],hist_CF_Lp_ALAp_exp[1],"hist_CF_Lp_ALAp_exp_sum");
  hist_CF_LL_ALAL_exp[0] = Calculate_CF(histRE_relK_LL,histME_relK_LL,"hist_CF_LL_exp",normleft,normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[1] = Calculate_CF(histRE_relK_ALAL,histME_relK_ALAL,"hist_CF_LL_exp",normleft,normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[2] = add_CF(hist_CF_LL_ALAL_exp[0],hist_CF_LL_ALAL_exp[1],"hist_CF_LL_ALAL_exp_sum");
  hist_CF_pp_ApAp_exp[0] = Calculate_CF(histRE_relK_pp,histME_relK_pp,"hist_CF_pp",normleft,normright, addon, spinningDepth);
  hist_CF_pp_ApAp_exp[1] = Calculate_CF(histRE_relK_ApAp,histME_relK_ApAp,"hist_CF_ApAp",normleft,normright, addon, spinningDepth);
  hist_CF_pp_ApAp_exp[2] = add_CF(hist_CF_pp_ApAp_exp[0],hist_CF_pp_ApAp_exp[1],"hist_CF_pp_ApAp_exp_sum");
  hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(histRE_relK_Xip,histME_relK_Xip,"hist_CF_pXi",normleft,normright, addon, spinningDepth);
  hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(histRE_relK_AXiAp,histME_relK_AXiAp,"hist_CF_ApAXi",normleft,normright, addon, spinningDepth);
  hist_CF_pXi_ApAXi_exp[2] = add_CF(hist_CF_pXi_ApAXi_exp[0],hist_CF_pXi_ApAXi_exp[1],"hist_CF_pXi_ApAXi_exp_sum");

  SetStyleHisto(hist_CF_pp_ApAp_exp[0], 1,1);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[0], 1,1);
  SetStyleHisto(hist_CF_LL_ALAL_exp[0], 1,1);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[0], 1,1);
  SetStyleHisto(hist_CF_pp_ApAp_exp[1], 0,2);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[1], 0,2);
  SetStyleHisto(hist_CF_LL_ALAL_exp[1], 0,2);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[1], 0,2);
  SetStyleHisto(hist_CF_pp_ApAp_exp[2],   0,3);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[2],   0,3);
  SetStyleHisto(hist_CF_LL_ALAL_exp[2],   0,3);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[2], 0,3);


  //  SetStyleHisto(histRE_relK_pp, 0,1);
  //  SetStyleHisto(histME_relK_pp, 0,2);
//======================================================================================================================================
  TCanvas *c1 = new TCanvas("pp-barpbarp","pp-barpbarp",0,0,2000,1100);
  c1->Divide(2,3);
  c1->cd(1);
  c1->cd(1)->SetRightMargin(right);
  c1->cd(1)->SetTopMargin(top);

  histRE_relK_pp->Scale(1./(double)
                        (histRE_relK_pp->Integral(histRE_relK_pp->FindBin(normleft),histRE_relK_pp->FindBin(normright))));
  histME_relK_pp->Scale(1./(double)
                        (histME_relK_pp->Integral(histME_relK_pp->FindBin(normleft),histME_relK_pp->FindBin(normright))));

  histRE_relK_pp->SetMarkerColor(fColors[0]);
  histRE_relK_pp->SetLineColor(fColors[0]);
  histRE_relK_pp->SetFillColor(fColors[0]);
  histRE_relK_pp->SetTitle("p-p; k* (GeV/#it{c}); #it{C}(k*)");
  histRE_relK_pp->GetXaxis()->SetRangeUser(0, 2);
  //histRE_relK_pp->GetYaxis()->SetRangeUser(0, 0.05);
  histRE_relK_pp->Draw("pe");
  histME_relK_pp->SetMarkerColor(fColors[1]);
  histME_relK_pp->SetLineColor(fColors[1]);
  histME_relK_pp->SetFillColor(fColors[1]);
  histME_relK_pp->Draw("pe same");

  auto* leg= new TLegend(0.65, 0.7, 0.95, 0.95);
  leg->AddEntry(histRE_relK_pp, "SE", "pe");
  leg->AddEntry(histME_relK_pp, "ME", "pe");
  leg->Draw("same");

  histRE_relK_ApAp->Scale(1./(double)
                          (histRE_relK_ApAp->Integral(histRE_relK_ApAp->FindBin(normleft),histRE_relK_ApAp->FindBin(normright))));
  histME_relK_ApAp->Scale(1./(double)
                          (histME_relK_ApAp->Integral(histME_relK_ApAp->FindBin(normleft),histME_relK_ApAp->FindBin(normright))));

  c1->cd(2);
  histRE_relK_ApAp->SetMarkerColor(fColors[0]);
  histRE_relK_ApAp->SetLineColor(fColors[0]);
  histRE_relK_ApAp->SetFillColor(fColors[0]);
  histRE_relK_ApAp->SetTitle("#bar{p}-#bar{p}; k* (GeV/#it{c}); #it{C}(k*)");
  histRE_relK_ApAp->GetXaxis()->SetRangeUser(0, 2);
  //histRE_relK_ApAp->GetYaxis()->SetRangeUser(0, 0.05);
  histRE_relK_ApAp->Draw("pe");
  histME_relK_ApAp->SetMarkerColor(fColors[1]);
  histME_relK_ApAp->SetLineColor(fColors[1]);
  histME_relK_ApAp->SetFillColor(fColors[1]);
  histME_relK_ApAp->Draw("pe same");

  c1->cd(3);
  hist_CF_pp_ApAp_exp[0]->Draw("pe");
  hist_CF_pp_ApAp_exp[0]->GetXaxis()->SetRangeUser(0, 2);
  hist_CF_pp_ApAp_exp[0]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_pp_ApAp_exp[0]->SetTitle("p-p; k* (GeV/#it{c}); #it{C}(k*)");


  c1->cd(4);
  hist_CF_pp_ApAp_exp[1]->Draw("pe");
  hist_CF_pp_ApAp_exp[1]->GetXaxis()->SetRangeUser(0, 2);
    hist_CF_pp_ApAp_exp[1]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_pp_ApAp_exp[1]->SetTitle("#bar{p}-#bar{p}; k* (GeV/#it{c}); #it{C}(k*)");

  c1->cd(5);
  hist_CF_pp_ApAp_exp[2]->Draw("pe");
  hist_CF_pp_ApAp_exp[2]->GetXaxis()->SetRangeUser(0, 2);
    hist_CF_pp_ApAp_exp[2]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_pp_ApAp_exp[2]->SetTitle("Sum p-p; k* (GeV/#it{c}); #it{C}(k*)");

  //======================================================================================================================================
  TCanvas *c2 = new TCanvas("pL-barpbarL","pL-barpbarL",0,0,2000,1100);
  c2->Divide(2,3);
  c2->cd(1);
  c2->cd(1)->SetRightMargin(right);
  c2->cd(1)->SetTopMargin(top);

  histRE_relK_Lp->Scale(1./(double)
                        (histRE_relK_Lp->Integral(histRE_relK_Lp->FindBin(normleft),histRE_relK_Lp->FindBin(normright))));
  histME_relK_Lp->Scale(1./(double)
                        (histME_relK_Lp->Integral(histME_relK_Lp->FindBin(normleft),histME_relK_Lp->FindBin(normright))));


  histRE_relK_Lp->SetMarkerColor(fColors[0]);
  histRE_relK_Lp->SetLineColor(fColors[0]);
  histRE_relK_Lp->SetFillColor(fColors[0]);
  histRE_relK_Lp->SetTitle("p-L; k* (GeV/#it{c}); #it{C}(k*)");
  histRE_relK_Lp->GetXaxis()->SetRangeUser(0, 2);
  //histRE_relK_Lp->GetYaxis()->SetRangeUser(0, 0.3);
  histRE_relK_Lp->Draw("pe");
  histME_relK_Lp->SetMarkerColor(fColors[1]);
  histME_relK_Lp->SetLineColor(fColors[1]);
  histME_relK_Lp->SetFillColor(fColors[1]);
  histME_relK_Lp->Draw("pe same");

  leg->Draw("same");

  histRE_relK_ALAp->Scale(1./(double)
                          (histRE_relK_ALAp->Integral(histRE_relK_ALAp->FindBin(normleft),histRE_relK_ALAp->FindBin(normright))));
  histME_relK_ALAp->Scale(1./(double)
                          (histME_relK_ALAp->Integral(histME_relK_ALAp->FindBin(normleft),histME_relK_ALAp->FindBin(normright))));

  c2->cd(2);
  histRE_relK_ALAp->SetMarkerColor(fColors[0]);
  histRE_relK_ALAp->SetLineColor(fColors[0]);
  histRE_relK_ALAp->SetFillColor(fColors[0]);
  histRE_relK_ALAp->SetTitle("#bar{p}-#bar{L}; k* (GeV/#it{c}); #it{C}(k*)");
  histRE_relK_ALAp->GetXaxis()->SetRangeUser(0, 2);
  //histRE_relK_ALAp->GetYaxis()->SetRangeUser(0, 0.3);
  histRE_relK_ALAp->Draw("pe");
  histME_relK_ALAp->SetMarkerColor(fColors[1]);
  histME_relK_ALAp->SetLineColor(fColors[1]);
  histME_relK_ALAp->SetFillColor(fColors[1]);
  histME_relK_ALAp->Draw("pe same");

  c2->cd(3);
  hist_CF_Lp_ALAp_exp[0]->Draw("pe");
  hist_CF_Lp_ALAp_exp[0]->GetXaxis()->SetRangeUser(0, 2);
    hist_CF_Lp_ALAp_exp[0]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_Lp_ALAp_exp[0]->SetTitle("p-L; k* (GeV/#it{c}); #it{C}(k*)");


  c2->cd(4);
  hist_CF_Lp_ALAp_exp[1]->Draw("pe");
  hist_CF_Lp_ALAp_exp[1]->GetXaxis()->SetRangeUser(0, 2);
    hist_CF_Lp_ALAp_exp[1]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_Lp_ALAp_exp[1]->SetTitle("#bar{p}-#bar{L}; k* (GeV/#it{c}); #it{C}(k*)");

  c2->cd(5);
  hist_CF_Lp_ALAp_exp[2]->Draw("pe");
  hist_CF_Lp_ALAp_exp[2]->GetXaxis()->SetRangeUser(0, 2);
    hist_CF_Lp_ALAp_exp[2]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_Lp_ALAp_exp[2]->SetTitle("Sum p-L; k* (GeV/#it{c}); #it{C}(k*)");

  //======================================================================================================================================
  TCanvas *c3 = new TCanvas("LL-barLbarL","LL-barLbarL",0,0,2000,1100);
  c3->Divide(2,3);
  c3->cd(1);
  c3->cd(1)->SetRightMargin(right);
  c3->cd(1)->SetTopMargin(top);

  histRE_relK_LL->Scale(1./(double)
                          (histRE_relK_LL->Integral(histRE_relK_LL->FindBin(normleft),histRE_relK_LL->FindBin(normright))));
  histME_relK_LL->Scale(1./(double)
                          (histME_relK_LL->Integral(histME_relK_LL->FindBin(normleft),histME_relK_LL->FindBin(normright))));

  histRE_relK_LL->SetMarkerColor(fColors[0]);
  histRE_relK_LL->SetLineColor(fColors[0]);
  histRE_relK_LL->SetFillColor(fColors[0]);
  histRE_relK_LL->SetTitle("L-L; k* (GeV/#it{c}); #it{C}(k*)");
  histRE_relK_LL->GetXaxis()->SetRangeUser(0, 2);
//    histRE_relK_LL->GetYaxis()->SetRangeUser(0, 0.3);
  histRE_relK_LL->Draw("pe");
  histME_relK_LL->SetMarkerColor(fColors[1]);
  histME_relK_LL->SetLineColor(fColors[1]);
  histME_relK_LL->SetFillColor(fColors[1]);
  histME_relK_LL->Draw("pe same");

  leg->Draw("same");

  histRE_relK_ALAL->Scale(1./(double)
                          (histRE_relK_ALAL->Integral(histRE_relK_ALAL->FindBin(normleft),histRE_relK_ALAL->FindBin(normright))));
  histME_relK_ALAL->Scale(1./(double)
                          (histME_relK_ALAL->Integral(histME_relK_ALAL->FindBin(normleft),histME_relK_ALAL->FindBin(normright))));

  c3->cd(2);
  histRE_relK_ALAL->SetMarkerColor(fColors[0]);
  histRE_relK_ALAL->SetLineColor(fColors[0]);
  histRE_relK_ALAL->SetFillColor(fColors[0]);
  histRE_relK_ALAL->SetTitle("#bar{L}-#bar{L}; k* (GeV/#it{c}); #it{C}(k*)");
  histRE_relK_ALAL->GetXaxis()->SetRangeUser(0, 2);
 // histRE_relK_ALAL->GetYaxis()->SetRangeUser(0, .3);
  histRE_relK_ALAL->Draw("pe");
  histME_relK_ALAL->SetMarkerColor(fColors[1]);
  histME_relK_ALAL->SetLineColor(fColors[1]);
  histME_relK_ALAL->SetFillColor(fColors[1]);
  histME_relK_ALAL->Draw("pe same");

  c3->cd(3);
  hist_CF_LL_ALAL_exp[0]->Draw("pe");
  hist_CF_LL_ALAL_exp[0]->GetXaxis()->SetRangeUser(0, 2);
    hist_CF_LL_ALAL_exp[0]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_LL_ALAL_exp[0]->SetTitle("L-L; k* (GeV/#it{c}); #it{C}(k*)");


  c3->cd(4);
  hist_CF_LL_ALAL_exp[1]->Draw("pe");
  hist_CF_LL_ALAL_exp[1]->GetXaxis()->SetRangeUser(0, 2);
    hist_CF_LL_ALAL_exp[1]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_LL_ALAL_exp[1]->SetTitle("#bar{L}-#bar{L}; k* (GeV/#it{c}); #it{C}(k*)");

  c3->cd(5);
  hist_CF_LL_ALAL_exp[2]->Draw("pe");
  hist_CF_LL_ALAL_exp[2]->GetXaxis()->SetRangeUser(0, 2);
    hist_CF_LL_ALAL_exp[2]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_LL_ALAL_exp[2]->SetTitle("Sum L-L; k* (GeV/#it{c}); #it{C}(k*)");

  TCanvas *c4 = new TCanvas("pXi-barpbarXi","pXi-barpbarXi",0,0,2000,1100);
  c4->Divide(2,3);
  c4->cd(1);
  c4->cd(1)->SetRightMargin(right);
  c4->cd(1)->SetTopMargin(top);


  histRE_relK_Xip->Scale(1./(double)
                          (histRE_relK_Xip->Integral(histRE_relK_Xip->FindBin(normleft),histRE_relK_Xip->FindBin(normright))));
  histME_relK_Xip->Scale(1./(double)
                          (histME_relK_Xip->Integral(histME_relK_Xip->FindBin(normleft),histME_relK_Xip->FindBin(normright))));


  histRE_relK_Xip->SetMarkerColor(fColors[0]);
  histRE_relK_Xip->SetLineColor(fColors[0]);
  histRE_relK_Xip->SetFillColor(fColors[0]);
  histRE_relK_Xip->SetTitle("p-Xi; k* (GeV/#it{c}); #it{C}(k*)");
  histRE_relK_Xip->GetXaxis()->SetRangeUser(0, 2);
//  histRE_relK_Xip->GetYaxis()->SetRangeUser(0, 0.3);
  histRE_relK_Xip->Draw("pe");
  histME_relK_Xip->SetMarkerColor(fColors[1]);
  histME_relK_Xip->SetLineColor(fColors[1]);
  histME_relK_Xip->SetFillColor(fColors[1]);
  histME_relK_Xip->Draw("pe same");

  leg->Draw("same");

  histRE_relK_AXiAp->Scale(1./(double)
                          (histRE_relK_AXiAp->Integral(histRE_relK_AXiAp->FindBin(normleft),histRE_relK_AXiAp->FindBin(normright))));
  histME_relK_AXiAp->Scale(1./(double)
                          (histME_relK_AXiAp->Integral(histME_relK_AXiAp->FindBin(normleft),histME_relK_AXiAp->FindBin(normright))));

  c4->cd(2);
  histRE_relK_AXiAp->SetMarkerColor(fColors[0]);
  histRE_relK_AXiAp->SetLineColor(fColors[0]);
  histRE_relK_AXiAp->SetFillColor(fColors[0]);
  histRE_relK_AXiAp->SetTitle("#bar{p}-#bar{Xi}; k* (GeV/#it{c}); #it{C}(k*)");
  histRE_relK_AXiAp->GetXaxis()->SetRangeUser(0, 2);
//  histRE_relK_AXiAp->GetYaxis()->SetRangeUser(0, 0.3);
  histRE_relK_AXiAp->Draw("pe");
  histME_relK_AXiAp->SetMarkerColor(fColors[1]);
  histME_relK_AXiAp->SetLineColor(fColors[1]);
  histME_relK_AXiAp->SetFillColor(fColors[1]);
  histME_relK_AXiAp->Draw("pe same");

  c4->cd(3);
  hist_CF_pXi_ApAXi_exp[0]->Draw("pe");
  hist_CF_pXi_ApAXi_exp[0]->GetXaxis()->SetRangeUser(0, 2);
    hist_CF_pXi_ApAXi_exp[0]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_pXi_ApAXi_exp[0]->SetTitle("p-Xi; k* (GeV/#it{c}); #it{C}(k*)");


  c4->cd(4);
  hist_CF_pXi_ApAXi_exp[1]->Draw("pe");
  hist_CF_pXi_ApAXi_exp[1]->GetXaxis()->SetRangeUser(0, 2);
        hist_CF_pXi_ApAXi_exp[1]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_pXi_ApAXi_exp[1]->SetTitle("#bar{p}-#bar{Xi}; k* (GeV/#it{c}); #it{C}(k*)");

  c4->cd(5);
  hist_CF_pXi_ApAXi_exp[2]->Draw("pe");
  hist_CF_pXi_ApAXi_exp[2]->GetXaxis()->SetRangeUser(0, 2);
        hist_CF_pXi_ApAXi_exp[2]->GetYaxis()->SetRangeUser(0.8, 1.3);
  hist_CF_pXi_ApAXi_exp[2]->SetTitle("Sum p-Xi; k* (GeV/#it{c}); #it{C}(k*)");


}
