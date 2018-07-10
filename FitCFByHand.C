#include "TPad.h"
Color_t fFillColors[8] = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9,
    kOrange-9, kCyan-8, kYellow-7};
Color_t fColors[8]  = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1,
    kOrange-1, kCyan+2, kYellow+2};
Style_t fMarkers[10] = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare,
    kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a and c are the weighting entities, i.e. wMean = (weightA*A + weightB*d) / (a+weightB)
float weightedMean(float weightA, float A, float weightB, float B)
{
  return (weightA*A + weightB*B)/(weightA+weightB);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a and weightB are the weighting entities, i.e. wMean = (weightA*A + weightB*B) / (a+weightB)
float weightedMeanError(float weightA, float A, float weightB, float B,
                        float weightAErr, float AErr, float weightBErr, float BErr)
{
  return std::sqrt(weightAErr*weightAErr*std::pow(((weightB*(A -
      B))/((weightA + weightB)*(weightA + weightB))), 2) + AErr*AErr*
                   std::pow(weightA/(weightA + weightB), 2)  +
                   weightBErr*weightBErr*std::pow((weightA* (-A + B))/((weightA +
                       weightB)*(weightA + weightB)) , 2) + BErr*BErr*
                       std::pow(weightB/(weightA + weightB), 2));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *getSignalHisto(TF1 *function, TH1F *histo, float
                     rangeLow, float rangeHigh, const char *name)
{
  const int firstBin = histo->FindBin(rangeLow);
  const int lastBin  = histo->FindBin(rangeHigh);
  TH1F *result = new TH1F(Form("result_%.2f_%.2f_%s", rangeLow, rangeHigh,
                               name), "", histo->GetNbinsX(), histo->GetXaxis()->GetXmin(),
                          histo->GetXaxis()->GetXmax());
  for(int i = firstBin; i<lastBin; ++i) {
    float weight = histo->GetBinContent(i) - function->Eval(histo->GetBinCenter(i));
    result->Fill(histo->GetBinCenter(i), weight);
    result->SetBinError(i, histo->GetBinError(i));
  }
  result->SetFillColor(kGray+1);
  result->SetLineColor(kGray+1);
  return result;
}
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
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.02);
  histo->GetYaxis()->SetTitleOffset(1.2);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// second order polynomial + double gaus to Lambda peak
void FitLambda(TH1F* histo, float lowerBound,float upperBound,TCanvas *cout, int iCD)
{
  //    if (x[0] > 50 && x[0] < 100) {
  std::vector<int> fBKGExlusion = {
      40,40,30,
      30,40,40,
      50,50,50
  };
  TF1 *fBackground = new TF1("fBackground", [&](double *x, double *p) {
    if ((x[0] < fBKGExlusion.at(iCD-1)  || x[0] > 55) && x[0] < 95) {
      TF1::RejectPoint();
      return (double)0; }
    return p[0]+TMath::Exp(p[1]+p[2]*x[0]);
  }, histo->GetXaxis()->GetXmin(), 350 , 3);
  TFitResultPtr backgroundR = histo->Fit("fBackground", "SR0", "",histo->GetXaxis()->GetXmin(), 220);

  // parse then to proper TF1
  TF1 *fBackground2 = new TF1("fBackground2","pol0(0)+expo(1)", 0, 350 );

  fBackground2->SetLineStyle(4);

  fBackground2->SetParameter(0, fBackground->GetParameter(0));
  fBackground2->SetParameter(1, fBackground->GetParameter(1));
  fBackground2->SetParameter(2, fBackground->GetParameter(2));
  //   remove background from signal
  TH1F *signalOnly = getSignalHisto(fBackground2, histo, 40,120,Form("%s_signal_only", histo->GetName()));
  signalOnly->GetXaxis()->SetRangeUser(0,500);
  //  signalOnly->DrawCopy();

  // fit signal only
  TF1 *fSignalSingleGauss = new TF1("fSignalSingleGauss","gaus(0)",40,100);
  fSignalSingleGauss->SetParLimits(1,55,80);
//  fSignalSingleGauss->SetParLimits(2,0,20);
  TFitResultPtr r = signalOnly->Fit("fSignalSingleGauss","SR0", "",50,95);
  //  fSignalSingleGauss->Draw("SAME");

  // Extract signal as integral
  float signal = fSignalSingleGauss->Integral(lowerBound, upperBound);
  //  float signalErr =
  //      fSignalSingleGauss->IntegralError(
  //          lowerBound, upperBound,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())
  //          /double(histo->GetBinWidth(1));
  std::cout << "Signal: " <<signal << std::endl;
  TF1 *fLambda = new TF1("fLambda","pol0(0)+expo(1)+gaus(3)", 0,300);
  fLambda->SetParameter(0, fBackground->GetParameter(0));
  fLambda->SetParameter(1, fBackground->GetParameter(1));
  fLambda->SetParameter(2, fBackground->GetParameter(2));
  fLambda->SetParameter(3, fSignalSingleGauss->GetParameter(0));
  fLambda->SetParameter(4, fSignalSingleGauss->GetParameter(1));
  fLambda->SetParameter(5, fSignalSingleGauss->GetParameter(2));
  fLambda->SetLineColor(fColors[2]);
  cout->cd(iCD);
  SetStyleHisto(histo,fMarkers[1],fColors[1]);
  histo->SetMinimum(0.6);
  histo->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  histo->DrawCopy();
  fBackground2->SetLineColor(fColors[2]);
  fBackground2->SetLineStyle(2);
  fBackground2->Draw("SAME");
  fLambda->Draw("SAME");
  float mean1 = fLambda->GetParameter(4);
  float width1 = fLambda->GetParameter(5);
  float background = fBackground2->Integral(lowerBound, upperBound);
  background-=(upperBound-lowerBound);
  std::cout << "Background: " << background << std::endl;
  //  float backgroundErr = fBackground2->IntegralError(lowerBound,
  //                                              upperBound, backgroundR->GetParams(),
  //                                              backgroundR->GetCovarianceMatrix().GetMatrixArray())
  //                          /double(histo->GetBinWidth(1));
  TLine* line = new TLine(lowerBound, 0.6, lowerBound, 1.6);
  line->Draw("Same");
  line->SetLineStyle(3);
  line->SetLineColorAlpha(fColors[0],0.7);
  TLine* line2 = new TLine(upperBound, 0.6, upperBound, 1.6);
  line2->Draw("Same");
  line2->SetLineStyle(3);
  line2->SetLineColorAlpha(fColors[0],0.7);
  TLine* line3 = new TLine(8, 1, 270, 1);
  line3->SetLineStyle(3);
  line3->SetLineColorAlpha(fColors[3],0.7);
  line3->Draw("Same");
  TLatex LambdaLabel;
  LambdaLabel.SetNDC(kTRUE);
  LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
  LambdaLabel.DrawLatex(gPad->GetUxmax()-0.6, gPad->GetUymax()-0.35,
                        Form("#splitline{"
                            "#splitline{< k* > = %.1fMeV/#it{c}^{2}}"
                            "{#sigma= %.1fMeV/#it{c}^{2}}}"
                            "{signal/(signal+bkg) = %.1f %%}",
                            mean1,
                            width1,
                            signal/(signal+background)*100.f));

}


void FitCFByHand() {
  SetStyle(false,false);
  TCanvas *cout = new TCanvas("HandMade","HandMade",0,0,1500,1100);
  cout->Divide(3,3);
  float lowerBound=55;
  float upperBound=95;
  int HistCounter = 1;
  float binWidth = 4;
  for (int i=1; i<6;++i) {
    if ( i == 4 ) continue;
    TFile *file = TFile::Open(Form("CFOutput_pXim_Rebin_%d.root",i),"READ");
    TH1F* CFInput =  (TH1F*)file->Get("hCkTotNormWeight_Shifted");
    if (!CFInput) {
      std::cout << "no Input hist";
      return;
    }
    CFInput->GetXaxis()->SetRangeUser(0,270);
    TLatex LambdaLabel;
    FitLambda(CFInput, lowerBound,upperBound,cout,HistCounter);
    cout->cd(HistCounter++);
    LambdaLabel.SetNDC(kTRUE);
    LambdaLabel.SetTextSize(gStyle->GetTextSize()*1.1);
    LambdaLabel.DrawLatex(gPad->GetUxmax()-0.6, gPad->GetUymax()-0.17,Form("Bin Width: %.0f MeV/#it{c}^{2}",binWidth*i));
  }
  float confLevel = 0.1;
  for (int iCon = 0; iCon < 5; ++iCon ) {
    TFile *file = TFile::Open(Form("CFOutput_pXim_RelFac_%.2f.root",confLevel),"READ");
    TH1F* CFInput =  (TH1F*)file->Get("hCkTotNormWeightIndepBinning");
    if (!CFInput) {
      std::cout << "no Input hist";
      return;
    }
    CFInput->GetXaxis()->SetRangeUser(0,270);
    FitLambda(CFInput, lowerBound,upperBound,cout,HistCounter);
    cout->cd(HistCounter++);
    TLatex LambdaLabel;
    LambdaLabel.SetNDC(kTRUE);
    LambdaLabel.SetTextSize(gStyle->GetTextSize()*1.1);
    LambdaLabel.DrawLatex(gPad->GetUxmax()-0.6, gPad->GetUymax()-0.17,Form("Confidence Level per Bin: %.2f", confLevel));
    confLevel+=0.05;
  }
}
