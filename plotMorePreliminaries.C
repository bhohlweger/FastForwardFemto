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
  histo->SetMarkerSize(1.25);
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
TGraphErrors *convertInGev(TGraphErrors *gr) {
 TGraphErrors *grOut = new TGraphErrors();
 grOut->SetName(gr->GetName());
 double x,y;
 for(int i=0; i< gr->GetN(); ++i) {
   gr->GetPoint(i, x, y);
   grOut->SetPoint(i, x/1000.f, y);
   grOut->SetPointError(i, gr->GetErrorX(i)/1000.f, gr->GetErrorY(i));
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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void plotMorePreliminaries(const char *expfile = "~/Results/LHC17p_fast/AnalysisResults.root", const char *CATSfileGaus = "", const char *CATSfileEPOS = "", const char *MoritaFile = "")
{
  gStyle->SetCanvasPreferGL(1);
  const float right = 0.025;
  const float top = 0.025;

  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;
  const float spinningDepth = 10.f;

  // for pythia comparison
  const float rebin = 2;

  const int energy = 13; // TeV
  const char *system = "pp";

  const char *prefix = "MB";
  const char *addon = "";

  float r = 1.188;
  float rErr = 0.009;
  float rSystErrUp = 0.016;
  float rSystErrDown = 0.009;

  //  EPOS
  float rEPOS = 1.476;
  float rEPOSErr = 0.015;
  float rEPOSSystErrUp = 0.007;
  float rEPOSSystErrDown = 0.014;

  SetStyle();

  // EXP DATA
  TFile* _file0=TFile::Open(expfile);
  TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon),Results);
  TList* tmpFolder=(TList*)Results->FindObject("Particle0_Particle0");
  TH1F* histRE_relK_pp = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle0");
  TH1F* histME_relK_pp = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle0");
  tmpFolder=(TList*)Results->FindObject("Particle1_Particle1");
  TH1F* histRE_relK_ApAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle1");
  TH1F* histME_relK_ApAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle1");
  tmpFolder=(TList*)Results->FindObject("Particle0_Particle2");
  TH1F* histRE_relK_Lp = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle2");
  TH1F* histME_relK_Lp = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle2");
  tmpFolder=(TList*)Results->FindObject("Particle1_Particle3");
  TH1F* histRE_relK_ALAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle3");
  TH1F* histME_relK_ALAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle3");
  tmpFolder=(TList*)Results->FindObject("Particle2_Particle2");
  TH1F* histRE_relK_LL = (TH1F*)tmpFolder->FindObject("SEDist_Particle2_Particle2");
  TH1F* histME_relK_LL = (TH1F*)tmpFolder->FindObject("MEDist_Particle2_Particle2");
  tmpFolder=(TList*)Results->FindObject("Particle3_Particle3");
  TH1F* histRE_relK_ALAL = (TH1F*)tmpFolder->FindObject("SEDist_Particle3_Particle3");
  TH1F* histME_relK_ALAL = (TH1F*)tmpFolder->FindObject("MEDist_Particle3_Particle3");
  tmpFolder=(TList*)Results->FindObject("Particle0_Particle4");
  TH1F* histRE_relK_Xip = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle4");
  TH1F* histME_relK_Xip = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle4");
  tmpFolder=(TList*)Results->FindObject("Particle1_Particle5");
  TH1F* histRE_relK_AXiAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle5");
  TH1F* histME_relK_AXiAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle5");
  TH1F *hist_CF_Lp_ALAp_exp[3];
  TH1F *hist_CF_LL_ALAL_exp[3];
  TH1F *hist_CF_pp_ApAp_exp[3];
  TH1F *hist_CF_pXi_ApAXi_exp[3];

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
  SetStyleHisto(hist_CF_pp_ApAp_exp[2], 0,0);
  SetStyleHisto(hist_CF_Lp_ALAp_exp[2], 0,0);
  SetStyleHisto(hist_CF_LL_ALAL_exp[2], 0,0);
  SetStyleHisto(hist_CF_pXi_ApAXi_exp[2], 0,0);
  hist_CF_LL_ALAL_exp[2]->SetMarkerSize(1.5);

  // SYSTEMATIC UNCERTAINTIES
  TString input_sys_pp = "C2totalsysPP.root";
  TString input_sys_pL = "C2totalsysPL.root";
  TString input_sys_LL = "C2totalsysLL.root";
  TString input_sys_pXi = "C2totalsysPXi.root";
  TFile* file_sys_pp = new TFile(input_sys_pp.Data());
  TFile* file_sys_pL = new TFile(input_sys_pL.Data());
  TFile* file_sys_pXi = new TFile(input_sys_pXi.Data());
  TFile* file_sys_LL = new TFile(input_sys_LL.Data());
  TH1F* hist_sys_pp = (TH1F*)file_sys_pp->Get("C2totalsysPP");
  hist_sys_pp->SetLineWidth(2.0);
  TH1F* hist_sys_pL = (TH1F*)file_sys_pL->Get("C2totalsysPL");
  hist_sys_pL->SetLineWidth(2.0);
  TH1F* hist_sys_LL = (TH1F*)file_sys_LL->Get("C2totalsysLL");
  hist_sys_LL->SetLineWidth(2.0);
  TH1F* hist_sys_pXi = (TH1F*)file_sys_pXi->Get("C2totalsysPXi");
  hist_sys_pXi->SetLineWidth(2.0);

  // UNCERTAINTIES ON THE FIT
  TFile *systFitGaus = TFile::Open(CATSfileGaus);
  TGraph *grpLDefaultLOGaus = nullptr;
  TGraph *grpLLowLOGaus = nullptr;
  TGraph *grpLUpLOGaus = nullptr;
  TGraphErrors *grFemtopLLOGaus = nullptr;
  grpLDefaultLOGaus = (TGraph*)systFitGaus->Get("pLamGraphDefault_LO");
  grpLLowLOGaus = (TGraph*)systFitGaus->Get("Copy_pLamGraphLowerLim_LO");
  grpLUpLOGaus = (TGraph*)systFitGaus->Get("pLamGraphUpperLim_LO");
  grFemtopLLOGaus = FemtoModelFitBands(grpLDefaultLOGaus, grpLLowLOGaus, grpLUpLOGaus);
  grFemtopLLOGaus->SetFillColor(fColors[3]);
  grFemtopLLOGaus->SetLineColor(fColors[3]);

  TFile *systFitEPOS = TFile::Open(CATSfileEPOS);
  TGraph *grpLDefaultNLOEPOS = nullptr;
  TGraph *grpLLowNLOEPOS = nullptr;
  TGraph *grpLUpNLOEPOS = nullptr;
  TGraphErrors *grFemtopLNLOEPOS = nullptr;
  grpLDefaultNLOEPOS = (TGraph*)systFitEPOS->Get("pLamGraphDefault_NLO");
  grpLLowNLOEPOS = (TGraph*)systFitEPOS->Get("Copy_pLamGraphLowerLim_NLO");
  grpLUpNLOEPOS = (TGraph*)systFitEPOS->Get("pLamGraphUpperLim_NLO");
  grFemtopLNLOEPOS = FemtoModelFitBands(grpLDefaultNLOEPOS, grpLLowNLOEPOS, grpLUpNLOEPOS);
  grFemtopLNLOEPOS->SetFillColor(fColors[1]);
  grFemtopLNLOEPOS->SetLineColor(fColors[1]);

  TGraph *grFakeSyst = new TGraph();
  grFakeSyst->SetFillColor(fFillColors[0]);
  grFakeSyst->SetLineColor(kBlack);
  grFakeSyst->SetMarkerStyle(20);
  grFakeSyst->SetMarkerSize(1.5);

  TGraph *grFakePP = new TGraph();
  grFakePP->SetLineColor(fColors[7]);
  grFakePP->SetLineWidth(4);
  TGraph *grFakePLnlo = new TGraph();
  grFakePLnlo->SetLineColor(fColors[1]);
  grFakePLnlo->SetLineWidth(4);
  TGraph *grFakePLlo = new TGraph();
  grFakePLlo->SetLineColor(fColors[3]);
  grFakePLlo->SetLineWidth(4);
  TGraph *grFakeLL = new TGraph();
  grFakeLL->SetLineColor(fColors[5]);
  grFakeLL->SetLineWidth(4);
  TGraph *grFakeLLStar = new TGraph();
  grFakeLLStar->SetLineColor(fColors[6]);
  grFakeLLStar->SetLineWidth(4);
  TGraph *grFakepXi= new TGraph();
  grFakepXi->SetLineColor(fColors[6]);
  grFakepXi->SetLineWidth(4);
  TGraph *grFakeXiCoulomb= new TGraph();
  grFakeXiCoulomb->SetLineColor(fColors[7]);
  grFakeXiCoulomb->SetLineWidth(4);

  TGraphErrors *Tgraph_syserror_pp_ApAp = DrawSystematicError(hist_CF_pp_ApAp_exp[2], hist_sys_pp, 0.002);
  TGraphErrors *Tgraph_syserror_pL_ApAL = DrawSystematicError(hist_CF_Lp_ALAp_exp[2], hist_sys_pL, 0.005);
  TGraphErrors *Tgraph_syserror_LL_ALAL = DrawSystematicError(hist_CF_LL_ALAL_exp[2], hist_sys_LL, 0.005);
  TGraphErrors *Tgraph_syserror_pXi_ApAXi = DrawSystematicError(hist_CF_pXi_ApAXi_exp[2], hist_sys_pXi, 0.005);

  auto* morita = TFile::Open(MoritaFile);
  auto* grnd46 = convertInGev((TGraphErrors*)morita->Get("gCkTheory_ND46"));
  grnd46->SetFillColor(fColors[2]);
  grnd46->SetLineColor(fColors[2]);
  grnd46->SetLineWidth(3);
  grnd46->SetLineStyle(3);
  auto* grnsc97f = convertInGev((TGraphErrors*)morita->Get("gCkTheory_NSC97f"));
  grnsc97f->SetFillColor(fColors[5]);
  grnsc97f->SetLineColor(fColors[5]);
  grnsc97f->SetLineWidth(3);
  grnsc97f->SetLineStyle(2);
  auto* grESC08 = convertInGev((TGraphErrors*)morita->Get("gCkTheory_ESC08"));
  grESC08->SetFillColor(fColors[4]);
  grESC08->SetLineColor(fColors[4]);
  grESC08->SetLineWidth(3);
  auto* grEhime = convertInGev((TGraphErrors*)morita->Get("gCkTheory_Ehime"));
  grEhime->SetFillColor(fColors[2]);
  grEhime->SetLineColor(fColors[2]);
  grEhime->SetLineWidth(3);
  auto* grHKMYY = convertInGev((TGraphErrors*)morita->Get("gCkTheory_HKMYY"));
  grHKMYY->SetFillColor(fColors[5]);
  grHKMYY->SetLineColor(fColors[5]);
  grHKMYY->SetLineWidth(3);
  auto* grNoSI = convertInGev((TGraphErrors*)morita->Get("gCkTheory_NoSI"));
  grNoSI->SetFillColor(kGray+1);
  grNoSI->SetLineColor(kGray+1);
  grNoSI->SetLineWidth(3);
  grNoSI->SetLineStyle(3);

  gStyle->SetErrorX(0.001);


  // PRELIMINARY PLOTS
  TCanvas *Can_CF_pL = new TCanvas("pL","pL", 0,0,650,550);
  Can_CF_pL->SetRightMargin(right);
  Can_CF_pL->SetTopMargin(top);
  Tgraph_syserror_pL_ApAL->SetLineColor(kWhite);
  Tgraph_syserror_pL_ApAL->Draw("ap");
  Tgraph_syserror_pL_ApAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_pL_ApAL->GetYaxis()->SetRangeUser(0.8, 2.);
  grFemtopLLOGaus->Draw("l3 same");
  grFemtopLNLOEPOS->Draw("l3 same");
//  if(!EPOS && grFemtopLLO) grFemtopLLO->Draw("l3 same");
  Tgraph_syserror_pL_ApAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_pL_ApAL->Draw("2 same");
  hist_CF_Lp_ALAp_exp[2]->Draw("pe same");
  TLegend *legLp2 = new TLegend(0.475,0.495,0.8,0.81);
  legLp2->SetBorderSize(0);
  legLp2->SetTextFont(42);
  legLp2->SetTextSize(gStyle->GetTextSize()*0.75);
  legLp2->AddEntry(hist_CF_Lp_ALAp_exp[2], "p#Lambda #oplus #bar{p}#bar{#Lambda} pairs", "pe");
  legLp2->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legLp2->AddEntry(grFemtopLLOGaus,"Femtoscopic fit (#chiEFT LO)","l");
  legLp2->AddEntry((TObject*)0x0, Form("Gaussian source #it{r_{0}} = %.3f fm", r),"");
  legLp2->AddEntry(grFemtopLNLOEPOS,"Femtoscopic fit (#chiEFT NLO)","l");
  legLp2->AddEntry((TObject*)0x0,Form("EPOS source #it{N}_{R} = %.3f", rEPOS),"");
  legLp2->Draw("same");
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.495, 0.875, "ALICE Preliminary");
  BeamText.DrawLatex(0.495, 0.825, Form("%s #sqrt{#it{s}} = %i TeV", system, energy));
  TLatex ref;
  ref.SetTextSize(gStyle->GetTextSize()*0.6);
  ref.SetNDC(kTRUE);
  ref.DrawLatex(0.5575, 0.465, "Nucl. Phys. A915 (2013) 24.");
  Can_CF_pL->Print("ANplot/CF_pL_Gaus-EPOS_prelim.pdf");

  TCanvas *Can_CF_LL = new TCanvas("LL","LL", 0,0,650,550);
  Can_CF_LL->SetRightMargin(right);
  Can_CF_LL->SetTopMargin(top);
  Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);
  Tgraph_syserror_LL_ALAL->Draw("ap");
  Tgraph_syserror_LL_ALAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_LL_ALAL->GetYaxis()->SetRangeUser(0.4, 2.3);
  grnd46->Draw("l same");
  grnsc97f->Draw("l same");
  grESC08->Draw("le3 same");
  grEhime->Draw("le3 same");
  grHKMYY->Draw("le3 same");
  grNoSI->Draw("l same");
  Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_LL_ALAL->Draw("2 same");
  hist_CF_LL_ALAL_exp[2]->Draw("pe same");
  TLegend *legLL2 = new TLegend(0.21,0.795,0.53,0.86);
  legLL2->SetBorderSize(0);
  legLL2->SetTextFont(42);
  legLL2->SetTextSize(gStyle->GetTextSize()*0.75);
  legLL2->AddEntry(grFakeSyst, "#Lambda#Lambda #oplus #bar{#Lambda}#bar{#Lambda} pairs", "pef");
  legLL2->Draw("same");
  TLegend *legLL = new TLegend(0.5,0.66,0.95,0.86);
  legLL->SetBorderSize(0);
  legLL->SetTextFont(42);
  legLL->SetTextSize(gStyle->GetTextSize()*0.75);
  legLL->SetNColumns(2);
  legLL->AddEntry(grnd46, "ND46", "l");
  legLL->AddEntry(grnsc97f, "NSC97f", "l");
  legLL->AddEntry(grESC08, "ESC08", "l");
  legLL->AddEntry(grHKMYY, "HKMYY", "l");
  legLL->AddEntry(grEhime, "Ehime", "l");
  legLL->AddEntry(grNoSI, "Quantum statistics", "l");
  legLL->Draw("same");
  ref.DrawLatex(0.555, 0.63, "PRC 91 (2015) 024916.");
  BeamText.DrawLatex(0.23, 0.875, Form("ALICE %s #sqrt{#it{s}} = %i TeV", system, energy));
  Can_CF_LL->Print("ANplot/CF_LL_Models.pdf");
}
