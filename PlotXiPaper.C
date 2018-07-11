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
    float weight = histo->GetBinContent(i) -
        function->Eval(histo->GetBinCenter(i));
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
void FitLambda(TH1F* histo, float &signal, float
               &signalErr, float &background, float &backgroundErr, float lowerBound,
               float upperBound)
{
  // Fit Background with second order polynomial, excluding Mlambda +/- 10 MeV
  TF1 *fBackground = new TF1("fBackground", [&](double *x, double *p) {
    if (x[0] > 1.31 && x[0] < 1.335) {TF1::RejectPoint(); return
        (double)0; } return p[0] + p[1]*x[0] + p[2]*x[0]*x[0]; }, 1.29, 1.35, 3);
  TFitResultPtr backgroundR = histo->Fit("fBackground", "SRQ0", "",1.295,1.349);

  // parse then to proper TF1
  TF1 *fBackground2 = new TF1("fBackground2","pol2", 0, 1.4);
  fBackground2->SetParameter(0, fBackground->GetParameter(0));
  fBackground2->SetParameter(1, fBackground->GetParameter(1));
  fBackground2->SetParameter(2, fBackground->GetParameter(2));

  // remove background from signal
  TH1F *signalOnly = getSignalHisto(fBackground2, histo, 1.3, 1.34,
                                    Form("%s_signal_only", histo->GetName()));

  // fit signal only
  TF1 *fSignalSingleGauss = new TF1("fSignalSingleGauss", "gaus(0)",1.31,1.34);
  //  signalOnly->DrawCopy();
  signalOnly->Fit("fSignalSingleGauss");
  TF1 *fSignalGauss = new TF1("fSignalGauss", "gaus(0) + gaus(3)", 1.3,
                              1.4);
  fSignalGauss->SetParameter(0, 0.75 * histo->GetMaximum());
  fSignalGauss->SetParameter(1, fSignalSingleGauss->GetParameter(1));
  fSignalGauss->SetParameter(2, 2.f*fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParLimits(2, 0.5*fSignalSingleGauss->GetParameter(2),1e2*2.f*fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParameter(3, 0.2 * histo->GetMaximum());
  fSignalGauss->SetParameter(4, fSignalSingleGauss->GetParameter(1));
  fSignalGauss->SetParLimits(4,
                             fSignalSingleGauss->GetParameter(1)-fSignalSingleGauss->GetParameter(2),
                             fSignalSingleGauss->GetParameter(1)+fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParameter(5, 0.5*fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParLimits(5, 0.5*fSignalSingleGauss->GetParameter(2),1e2*2.f*fSignalSingleGauss->GetParameter(2));
  TFitResultPtr r = signalOnly->Fit("fSignalGauss", "SRQ0", "", 1.29,
                                    1.38);

  // Extract signal as integral
  signal = fSignalGauss->Integral(lowerBound, upperBound)
                    /double(histo->GetBinWidth(1));
  signalErr = fSignalGauss->IntegralError(lowerBound, upperBound,
                                          r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())
                    /double(histo->GetBinWidth(1));

  TF1 *fLambda = new TF1("fLambda", "fBackground2 + fSignalGauss", 1.25,
                         1.4);
  fLambda->SetNpx(1000);
  fLambda->SetParameter(3, 0.75 * histo->GetMaximum());
  fLambda->SetParameter(4, fSignalGauss->GetParameter((1)));
  fLambda->SetParameter(5, fSignalGauss->GetParameter((2)));
  fLambda->SetParameter(6, 0.2 * histo->GetMaximum());
  fLambda->SetParameter(7, fSignalGauss->GetParameter((4)));
  fLambda->SetParameter(8, fSignalGauss->GetParameter((5)));
  fLambda->SetLineColor(fColors[2]);
  histo->Fit("fLambda", "SRQ", "", 1.28, 1.4);

  TF1 *fLambda_background = new TF1("fLambda_background", "pol2(0)",
                                    1.28, 1.45);
  fLambda_background->SetParameter(0, fLambda->GetParameter(0));
  fLambda_background->SetParameter(1, fLambda->GetParameter(1));
  fLambda_background->SetParameter(2, fLambda->GetParameter(2));
  fLambda_background->SetLineStyle(3);
  fLambda_background->SetLineColor(fColors[2]);

  background = fLambda_background->Integral(lowerBound, upperBound)
                    /double(histo->GetBinWidth(1));
  backgroundErr = fLambda_background->IntegralError(lowerBound,
                                                    upperBound, backgroundR->GetParams(),
                                                    backgroundR->GetCovarianceMatrix().GetMatrixArray())
                    /double(histo->GetBinWidth(1));

  histo->GetListOfFunctions()->Add(fLambda_background);
//  TH1F* FitValues = new TH1F("results","results",10,0,10);
//  FitValues->SetBinContent(1,signal);
//  FitValues->SetBinContent(2,signalErr);
//  FitValues->SetBinContent(1,background);
//  FitValues->SetBinContent(2,backgroundErr);
//
//  histo->GetListOfFunctions()->Add(FitValues);

  delete signalOnly;
  delete fSignalGauss;
  delete fSignalSingleGauss;
}

void executeFit(int iBin, TCanvas *PtXi, TH1F *xiPt, TH1F *Purity,
                double ptmin, double ptmax,TString period,float &signalOut,
                float &backgroundOut, float &massMeanOut, float &widthMeanOut) {

  TPad *padPt;
  if (iBin >= 0 ) {
    padPt=(TPad*)PtXi->cd(1+iBin);
    padPt->SetRightMargin(0.0);
    padPt->SetTopMargin(0.008);
    padPt->SetBottomMargin(0.15);
    padPt->SetLeftMargin(0.07);
    padPt->Draw();
  } else {
    PtXi->cd();
  }
  float signal=0;
  float signalErr=0;
  float background=0;
  float backgroundErr=0;
  float lowerBound=1.322-0.005;
  float upperBound=1.322+0.005;
  FitLambda(xiPt, signal, signalErr, background, backgroundErr, lowerBound,upperBound);

  TList *funListLambda = xiPt->GetListOfFunctions();
  TF1 *fLambdaTotal = (TF1*)funListLambda->FindObject("fLambda");

  float amp1 = fLambdaTotal->GetParameter(3);
  float amp2 = fLambdaTotal->GetParameter(6);
  float mean1 = fLambdaTotal->GetParameter(4);
  float mean2 = fLambdaTotal->GetParameter(7);
  float width1 = fLambdaTotal->GetParameter(5);
  float width2 = fLambdaTotal->GetParameter(8);

  float meanMass = weightedMean(amp1, mean1, amp2, mean2);
  float meanWidth = weightedMean(amp1, width1, amp2, width2);
  signalOut = signal;
  backgroundOut = background;
  massMeanOut = meanMass;
  widthMeanOut = meanWidth;

  if (iBin >= 0) {
    TLatex LambdaLabel;
    LambdaLabel.SetNDC(kTRUE);
    LambdaLabel.SetTextSize(gStyle->GetTextSize()*1.2);
    LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,
                          Form("#splitline{#splitline{#splitline{#splitline{#splitline{%.2f < p_{T} < %.2f}{#Xi^{-}: %.0f}}{m_{#Xi} = %.1fMeV/#it{c}^{2}}}{#sigma_{#Xi}= %.1fMeV/#it{c}^{2}}}{Purity = %.1f %%}}{Sig/Bkg = %.1f}",
                               ptmin,ptmax,signal, meanMass*1000.f, meanWidth*1000.f,
                               signal/(signal+background)*100.f,
                               signal/background));//#splitline{#splitline{
  } else {
    TLatex LambdaLabel;
    LambdaLabel.SetNDC(kTRUE);
    LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
    if (iBin == -1) {
      LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.42,
                          Form("#splitline{#splitline{#splitline{#splitline{#Xi^{-}: %.0f}{m_{#Xi} = %.1fMeV/#it{c}^{2}}}{#sigma_{#Xi}= %.1fMeV/#it{c}^{2}}}{Purity = %.1f %%}}{Sig/Bkg = %.1f}",
                               signal, meanMass*1000.f, meanWidth*1000.f,
                               signal/(signal+background)*100.f,
                               signal/background));
    } else if (iBin == -2 ) {
      LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.42,
                          Form("#splitline{#splitline{#splitline{#splitline{#bar{#Xi}^{-}: %.0f}{m_{#Xi} = %.1fMeV/#it{c}^{2}}}{#sigma_{#Xi}= %.1fMeV/#it{c}^{2}}}{Purity = %.1f %%}}{Sig/Bkg = %.1f}",
                               signal, meanMass*1000.f, meanWidth*1000.f,
                               signal/(signal+background)*100.f,
                               signal/background));
    }
    TLatex BeamTextLambda;
    BeamTextLambda.SetNDC(kTRUE);
    BeamTextLambda.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.2,
                             period.Data());
  }

  TLine line;
  line.SetLineColor(kGray);
  line.SetLineStyle(9);
  line.SetLineWidth(3);
  line.DrawLine(1.317,0,1.317,5e4);
  line.DrawLine(1.327,0,1.327,5e4);
  if(Purity) {
    Purity->SetBinContent(iBin+1,signal/(signal+background)*100.f);
  }
  return;
}

void plotBothPeaks(TCanvas *c, TH1F *Xi, TH1F *AXi,
                   double Xisignal1,double Xibackground1,
                   double XimassMean1, double XiwidthMean1,
                   double Xisignal2,double Xibackground2,
                   double XimassMean2, double XiwidthMean2) {
  c->cd();
  double iMmin=1.305;
  double iMmax=1.335;
  TH1F *dummy2 = new TH1F("dummy2","dummy2",10,iMmin,iMmax);
  dummy2->GetXaxis()->SetTitle("M_{#Lambda#pi} (GeV/#it{c}^{2})");
  dummy2->GetXaxis()->SetTitleOffset(0.8);
  dummy2->GetYaxis()->SetTitle("d#it{N}/d#it{M} (GeV/#it{c}^2)^{-1}");
  dummy2->SetMaximum(Xi->GetMaximum()*1.1);
  dummy2->GetXaxis()->SetLabelSize(0.045);
  dummy2->GetXaxis()->SetTitleSize(0.05);
  dummy2->GetXaxis()->SetLabelOffset(0.01);
  dummy2->GetXaxis()->SetTitleOffset(1.2);
  dummy2->GetXaxis()->SetLabelFont(42);
  dummy2->GetYaxis()->SetLabelSize(0.045);
  dummy2->GetYaxis()->SetTitleSize(0.05);
  dummy2->GetYaxis()->SetLabelOffset(0.02);
  dummy2->GetYaxis()->SetTitleOffset(0.8);
  dummy2->DrawCopy();

  TList *funListLambda =AXi->GetListOfFunctions();

  TF1 *fLambdaTotal = (TF1*)funListLambda->FindObject("fLambda");
  fLambdaTotal->SetFillColor(fColors[3]);
  fLambdaTotal->SetMarkerColor(fColors[3]);
  fLambdaTotal->SetLineColor(fColors[3]);

  TF1 *fLambdaBackground = (TF1*)funListLambda->FindObject("fLambda_background");
  fLambdaBackground->SetFillColor(fColors[3]);
  fLambdaBackground->SetMarkerColor(fColors[3]);
  fLambdaBackground->SetLineColor(fColors[3]);
//  Xi->GetXaxis()->SetRangeUser(1.305,1.335);
  Xi->SetFillColor(fColors[2]);
  Xi->SetMarkerColor(fColors[2]);
  Xi->SetLineColor(fColors[2]);

  TLegend *leg = new TLegend(0.65, 0.6, 0.75, 0.75);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(4000);
  leg->SetTextSize(gStyle->GetTextSize()*0.85);
  leg->AddEntry(Xi, "#Xi^{-} Candidates: M_{#Lambda#pi^{-}}", "PL");
  leg->AddEntry(AXi, "#bar{#Xi}^{-} Candidates: M_{#bar{#Lambda}#pi^{+}}", "PL");
  Xi->DrawCopy("SAME");

  AXi->SetFillColor(fColors[3]);
  AXi->SetMarkerColor(fColors[3]);
  AXi->SetLineColor(fColors[3]);

  AXi->DrawCopy("same");
  TLatex LambdaLabel;
  LambdaLabel.SetNDC(kTRUE);
  LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
  LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.42,
                         Form("#splitline{#splitline{#splitline{#splitline{#Xi^{-}(#bar{#Xi}^{+}): %.0f (%0.f)}{m_{#Xi}(m_{#bar{#Xi}}) = %.1f (%.1f) MeV/#it{c}^{2}}}{#sigma_{#Xi} (#sigma_{#bar{#Xi}})= %.1f (%.1f) MeV/#it{c}^{2}}}{Purity = %.1f(%.1f) %%}}{Sig/Bkg = %.1f(%.1f)}",
                              Xisignal1,Xisignal2, XimassMean1*1000.f,XimassMean2*1000.f, XiwidthMean1*1000.f,XiwidthMean2*1000.f,
                              Xisignal1/(Xisignal1+Xibackground1)*100.f,Xisignal2/(Xisignal2+Xibackground2)*100.f,
                              Xisignal1/Xibackground1,Xisignal2/Xibackground2));
  LambdaLabel.DrawLatex(0.65,0.76,"p-Pb #sqrt{#it{s_{NN}}} = 5.02 TeV");
  LambdaLabel.DrawLatex(0.780946,0.065,"M_{#Lambda#pi} (GeV/#it{c}^{2})");

  TLine line;
  line.SetLineColor(kGray);
  line.SetLineStyle(9);
  line.SetLineWidth(3);
  line.DrawLine(1.317,0,1.317,5e4);
  line.DrawLine(1.327,0,1.327,5e4);
  leg->Draw("Same");
}

void plotBothPeaks2(TCanvas *c, TH1F *Xi, TH1F *AXi) {
  c->cd();
  double iMmin=1.3219 -0.03;//1.2919
  double iMmax=1.3219 +0.03;//1.3519
  TH1F *dummy = new TH1F("dummy","dummy",10,iMmin,iMmax);
  dummy->GetXaxis()->SetTitle("M_{#Lambda #pi} (GeV/#it{c}^{2})");
  dummy->GetYaxis()->SetTitle("d#it{N}/d#it{M} (GeV/#it{c}^{2})^{-1}");
  dummy->SetMaximum(Xi->GetMaximum()*1.1);
  dummy->GetXaxis()->SetLabelSize(0.045);
  dummy->GetXaxis()->SetTitleSize(0.05);
  dummy->GetXaxis()->SetLabelOffset(0.01);
  dummy->GetXaxis()->SetTitleOffset(1.2);
  dummy->GetXaxis()->SetLabelFont(42);
  dummy->GetYaxis()->SetLabelSize(0.045);
  dummy->GetYaxis()->SetTitleSize(0.05);
  dummy->GetYaxis()->SetLabelOffset(0.02);
  dummy->GetYaxis()->SetTitleOffset(0.8);
  dummy->DrawCopy("");


//  Xi->SetFillColor(fColors[2]);
  Xi->SetMarkerColor(fColors[2]);
  Xi->SetLineColor(fColors[2]);
  Xi->SetLineWidth(3);

//  AXi->SetFillColor(fColors[3]);
  AXi->SetMarkerColor(fColors[3]);
  AXi->SetLineColor(fColors[3]);
  AXi->SetLineWidth(3);

  TLegend *leg = new TLegend(0.65, 0.6, 0.75, 0.75);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(4000);
  leg->SetTextSize(gStyle->GetTextSize()*0.85);
  leg->AddEntry(Xi, "#Xi^{-} Candidates: M_{#Lambda#pi^{-}}", "PL");
  leg->AddEntry(AXi, "#bar{#Xi}^{+} Candidates: M_{#bar{#Lambda}#pi^{+}}", "PL");
  Xi->DrawCopy("HIST same");

  AXi->DrawCopy("HIST same");

  TLine line;
  line.SetLineColor(kGray);
  line.SetLineStyle(9);
  line.SetLineWidth(3);
  line.DrawLine(1.317,0,1.317,5e4);
  line.DrawLine(1.327,0,1.327,5e4);

  TLatex LambdaLabel;
  LambdaLabel.SetNDC(kTRUE);
  LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
//  LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.42,
//                         Form("#splitline{#splitline{#splitline{#splitline{#Xi^{-}(#bar{#Xi}^{+}): %.0f (%0.f)}{m_{#Xi}(m_{#bar{#Xi}}) = %.1f (%.1f) MeV/#it{c}^{2}}}{#sigma_{#Xi} (#sigma_{#bar{#Xi}})= %.1f (%.1f) MeV/#it{c}^{2}}}{Purity = %.1f(%.1f) %%}}{Sig/Bkg = %.1f(%.1f)}",
//                              Xisignal1,Xisignal2, XimassMean1*1000.f,XimassMean2*1000.f, XiwidthMean1*1000.f,XiwidthMean2*1000.f,
//                              Xisignal1/(Xisignal1+Xibackground1)*100.f,Xisignal2/(Xisignal2+Xibackground2)*100.f,
//                              Xisignal1/Xibackground1,Xisignal2/Xibackground2));
  LambdaLabel.DrawLatex(0.65,0.76,"p-Pb #sqrt{#it{s_{NN}}} = 5.02 TeV");
  LambdaLabel.DrawLatex(0.780946,0.065,"M_{#Lambda#pi} (GeV/#it{c}^{2})");
  leg->Draw("Same");
}

void PlotXiPaper(TString fileName, const char *prefix) {
  TGaxis::SetMaxDigits(3);
  double iMmin=1.285;
  double iMmax=1.345;
  SetStyle(false,false);
  TH1F* xiPt[13];
  TH1F* axiPt[13];

  TH2F* xiMass2D;
  TH2F* axiMass2D;
  TH1F* xiMass;
  TH1F* axiMass;
  TFile *file=TFile::Open(fileName.Data());
  //  file->ls();
  TString dirNameXi="";
  TString dirNameAXi="";
  dirNameXi+=prefix;
  dirNameAXi+=prefix;
  dirNameXi+="CascadeCuts";
  dirNameAXi+="AntiCascadeCuts";
  TDirectoryFile *dirExpXi=(TDirectoryFile*)(file->FindObjectAny(dirNameXi.Data()));
  TDirectoryFile *dirExpAXi=(TDirectoryFile*)(file->FindObjectAny(dirNameAXi.Data()));

  TCanvas *Xi= new TCanvas("Xi","Xi",0,0,1500,1100);
  TCanvas *AXi= new TCanvas("AXi","AXi",0,0,1500,1100);
  TCanvas *XiAXi= new TCanvas("XiAXi","XiAXi",0,0,1500,1100);
  TCanvas *XiAXiClone= new TCanvas("XiAXiClone","XiAXiClone",0,0,1500,1100);
  TCanvas *PtXi= new TCanvas("PtXi","PtXi",0,0,1500,1100);
  TCanvas *PtAXi= new TCanvas("PtAXi","PtAXi",0,0,1500,1100);
  TH1F *Purity=0;
  TH1F *aPurity=0;
  PtXi->Divide(4,4);
  PtAXi->Divide(4,4);
  if (dirExpXi&&dirExpAXi) {
    TList *tmpXi=dirExpXi->GetListOfKeys();
    TList *tmpAXi=dirExpAXi->GetListOfKeys();

    TString nameXi=tmpXi->At(0)->GetName();
    TString nameAXi=tmpAXi->At(0)->GetName();

    TList *outputXi;
    TList *outputAXi;
    dirExpXi->GetObject(nameXi,outputXi);
    dirExpAXi->GetObject(nameAXi,outputAXi);
    if (!outputXi && !outputAXi) {
      std::cout << "No Output \n";
    }
    //output->ls();
    TList *Cascade=(TList*)outputXi->FindObject("Cascade");
    TList *AntiCascade=(TList*)outputAXi->FindObject("Cascade");
    if(!Cascade) {
      std::cout << "No Cascade \n";
    }
    if(!AntiCascade) {
      std::cout << "No Anti Cascade \n";
    }
    xiMass2D=(TH2F*)Cascade->FindObject("InvMassXiPt");
    xiMass2D->SetNameTitle("XiMassPt","XiMassPt");
    xiMass=(TH1F*)xiMass2D->ProjectionY("InvMassXi");
    TH1F* XiMassClone = (TH1F*)xiMass->Clone("CloneXiMass");
    Purity=new TH1F("Purity","Purity",13,xiMass2D->GetXaxis()->GetBinLowEdge(1),xiMass2D->GetXaxis()->GetBinUpEdge(13));

    axiMass2D=(TH2F*)AntiCascade->FindObject("InvMassXiPt");
    axiMass=(TH1F*)axiMass2D->ProjectionY("InvMassaXi");
    TH1F* aXiMassClone =(TH1F*)axiMass->Clone("CloneAXiMass");
    axiMass->SetNameTitle("aXiMassPt","aXiMassPt");
    aPurity=new TH1F("aPurity","aPurity",13,axiMass2D->GetXaxis()->GetBinLowEdge(1),axiMass2D->GetXaxis()->GetBinUpEdge(13));
    float signal, background,mass,width;

    for (int iBin = 0; iBin<13; ++iBin) {
      TString XiNamePt=Form("XiPt_%i",iBin);
      xiPt[iBin]=(TH1F*)xiMass2D->ProjectionY(XiNamePt.Data(),iBin+1,iBin+2);
      SetStyleHisto(xiPt[iBin],2,1);
      xiPt[iBin]->GetXaxis()->SetRangeUser(iMmin,iMmax);
      xiPt[iBin]->GetXaxis()->SetTitle("M_{#pi#Lambda} (GeV/#it{c}^{2})");

      executeFit(iBin, PtXi, xiPt[iBin],Purity,xiMass2D->GetXaxis()->GetBinCenter(iBin+1),
                 xiMass2D->GetXaxis()->GetBinCenter(iBin+2),"",signal,background,mass,width);

      TString aXiNamePt=Form("aXiPt_%i",iBin);
      axiPt[iBin]=(TH1F*)axiMass2D->ProjectionY(aXiNamePt.Data(),iBin+1,iBin+2);
      SetStyleHisto(axiPt[iBin],2,1);
      axiPt[iBin]->GetXaxis()->SetRangeUser(iMmin,iMmax);
      axiPt[iBin]->GetXaxis()->SetTitle("M_{#bar{#pi}#bar{#Lambda}} (GeV/#it{c}^{2})");

      executeFit(iBin, PtAXi, axiPt[iBin],aPurity,axiMass2D->GetXaxis()->GetBinCenter(iBin+1),
                 axiMass2D->GetXaxis()->GetBinCenter(iBin+2),"",signal,background,mass,width);
    }

    SetStyleHisto(xiMass,20,1);
    xiMass->GetXaxis()->SetRangeUser(iMmin,iMmax);
    xiMass->GetXaxis()->SetTitle("M_{#pi#Lambda} (GeV/#it{c}^{2})");
    xiMass->GetYaxis()->SetTitle("d#it{N}/d#it{M} (GeV/#it{c}^{2})^{-1}");
    xiMass->GetXaxis()->SetTitleOffset(0.9);
    xiMass->GetYaxis()->SetTitleOffset(0.8);
    if (!xiMass) {
      std::cout << "Missing Hist" << std::endl;
    }
    float signal1, background1,mass1,width1;
    float signal2, background2,mass2,width2;

    executeFit(-1,Xi,xiMass,nullptr,0.3,6.5,"p-Pb (2016) #sqrt{#it{s_{NN}}} = 5.02 TeV",
               signal1,background1,mass1,width1);

    SetStyleHisto(axiMass,20,1);
    axiMass->GetXaxis()->SetRangeUser(iMmin,iMmax);
    axiMass->GetXaxis()->SetTitle("M_{#bar{#pi}#bar{#Lambda}} (GeV/#it{c}^{2})");
    axiMass->GetYaxis()->SetTitle("d#it{N}/d#it{M} (GeV/#it{c}^{2})^{-1}");
    axiMass->GetXaxis()->SetTitleOffset(0.9);
    axiMass->GetYaxis()->SetTitleOffset(0.8);

    if (!axiMass) {
      std::cout << "Missing Hist" << std::endl;
    }
    executeFit(-2,AXi,axiMass,nullptr,0.3,6.5,"p-Pb (2016) #sqrt{#it{s_{NN}}} = 5.02 TeV",
               signal2,background2,mass2,width2);

    plotBothPeaks(
        XiAXi,xiMass,axiMass,signal1,background1,mass1,width1,signal2,background2,mass2,width2);
    plotBothPeaks2(
        XiAXiClone,XiMassClone,aXiMassClone);

    TCanvas *cPurity=new TCanvas("cPurity","cPurity");
    cPurity->cd();
    SetStyleHisto(Purity,21,2);
    Purity->GetYaxis()->SetRangeUser(80,100);
    Purity->SetTitle(";P_{T} (GeV/#it{c});Purity");
    Purity->SetStats(0);
    Purity->GetYaxis()->SetTitleOffset(1.);
    Purity->GetXaxis()->SetTitleOffset(.9);
    Purity->DrawCopy("p");
  } else {
    std::cout << "dirExp missing \n";
  }
}
