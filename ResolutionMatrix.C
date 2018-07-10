
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
  gStyle->SetPalette(kRainBow);

//  const int NRGBs = 6;
//  Double_t stops[NRGBs];
//  for(int i=0; i<NRGBs; ++i) stops[i] = float(i)/(NRGBs-1);
//
//  Double_t red[NRGBs]   = { 1.,  29./255., 25./255., 27./255., 32./255., 24./255.};
//  Double_t green[NRGBs] = { 1., 221./255., 160./255., 113./255., 74./255., 37./255.};
//  Double_t blue[NRGBs] = {  1., 221./255., 184./255., 154./255., 129./255., 98./255.};
//  TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHisto(TH2 *histo)
{
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
}

void ResolutionMatrix(const char* filename) {
  SetStyle(false,true);
  TFile *reso=TFile::Open(filename);
  TH2F* ResoHist[4];
  ResoHist[0] = (TH2F*)reso->Get("hSigmaMeV_Proton_Proton");
  ResoHist[1] = (TH2F*)reso->Get("hSigmaMeV_Proton_Lambda");
  ResoHist[2] = (TH2F*)reso->Get("hSigmaMeV_Lambda_Lambda");
  ResoHist[3] = (TH2F*)reso->Get("hSigmaMeV_Proton_Xim");
  TCanvas* c1 = new TCanvas("c1","c1",0,0,2000,1100);
  TString PairName[4] = {"p-p", "p-#Lambda", "#Lambda-#Lambda","p-#Xi"};
  c1->Divide(2,2);
  c1->cd(1);
  for (int iBin = 1; iBin < 5; ++iBin) {
    TPad *padPt=(TPad*)c1->cd(iBin);
    padPt->SetRightMargin(0.11);
      padPt->SetTopMargin(0.08);
      padPt->SetBottomMargin(0.13);
      padPt->SetLeftMargin(0.12);
    padPt->Draw();
    SetStyleHisto(ResoHist[iBin-1]);
    ResoHist[iBin-1]->SetTitle(Form("Resolution %s;k*_{Generated} (MeV/c);k*_{Reconstructed} (MeV/c)",PairName[iBin-1].Data()));
    ResoHist[iBin-1]->DrawCopy("COLZ");
  }
}
