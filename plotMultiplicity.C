
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
histo->SetMarkerStyle(fMarkers[marker]);
histo->SetMarkerColor(fColors[color]);
histo->SetLineColor(fColors[color]);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* Calculate_CF(TH1F* histRE_relK,TH1F* histME_relK, TString CFname,Double_t normleft,Double_t normright)
{
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  Double_t norm_relK = histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) / histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright));

  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  Hist_CF->Divide(histRE_relK,histME_relK,1,norm_relK);

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
void plotMultiplicity(const char* file="~/Results/LHC17p_fast/AnalysisResults.root")
{


  const float right = 0.025;
  const float top = 0.025;

  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;

  TFile* _file0=TFile::Open(file);
  TList *listTP=0;
  _file0->GetObject("/PWGCF_PLFemto_0/TP_dir_0",listTP);
  if(!listTP) _file0->GetObject("/PWGCF_PLFemto_0/TPdir_0",listTP);
  TH1F* histRE_relK_Lp[12];
  TH1F* histME_relK_Lp[12];
  TH1F* histRE_relK_ALAp[12];
  TH1F* histME_relK_ALAp[12];
  TH1F* histRE_relK_LL[12];
  TH1F* histME_relK_LL[12];
  TH1F* histRE_relK_ALAL[12];
  TH1F* histME_relK_ALAL[12];
  TH1F* histRE_relK_pp[12];
  TH1F* histME_relK_pp[12];
  TH1F* histRE_relK_ApAp[12];
  TH1F* histME_relK_ApAp[12];
  TH1F *hist_CF_Lp_ALAp_exp[3][12];
  TH1F *hist_CF_LL_ALAL_exp[3][12];
  TH1F *hist_CF_pp_ApAp_exp[3][12];

  TH1F* histRE_relK_Lp_merge[3];
  TH1F* histME_relK_Lp_merge[3];
  TH1F* histRE_relK_ALAp_merge[3];
  TH1F* histME_relK_ALAp_merge[3];
  TH1F* histRE_relK_LL_merge[3];
  TH1F* histME_relK_LL_merge[3];
  TH1F* histRE_relK_ALAL_merge[3];
  TH1F* histME_relK_ALAL_merge[3];
  TH1F* histRE_relK_pp_merge[3];
  TH1F* histME_relK_pp_merge[3];
  TH1F* histRE_relK_ApAp_merge[3];
  TH1F* histME_relK_ApAp_merge[3];
  TH1F *hist_CF_Lp_ALAp_exp_merge[3][3];
  TH1F *hist_CF_LL_ALAL_exp_merge[3][3];
  TH1F *hist_CF_pp_ApAp_exp_merge[3][3];

  auto *c = new TCanvas("Can_CF_fitting","Can_CF_fitting",0,0,1500,550);
  c->Divide(3,1);
  for(int i =0; i<3; ++i) {
    c->cd(i+1)->SetRightMargin(right);
    c->cd(i+1)->SetTopMargin(top);
  }

  auto* dummy1 = new TGraph();
  dummy1->SetPoint(0, 0, 0);
  auto* dummy2 = new TGraph();
  dummy2->SetPoint(0, 0, 0);
  auto* dummy3 = new TGraph();
  dummy3->SetPoint(0, 0, 0);

  c->cd(1);
  dummy1->Draw("ape");
  dummy1->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  dummy1->GetXaxis()->SetRangeUser(0, 0.125);
  dummy1->GetXaxis()->SetNdivisions(505);
  dummy1->GetYaxis()->SetRangeUser(0, 8);

  c->cd(2);
  dummy2->Draw("ape");
  dummy2->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  dummy2->GetXaxis()->SetRangeUser(0, 0.2);
  dummy2->GetXaxis()->SetNdivisions(505);
  dummy2->GetYaxis()->SetRangeUser(0., 8);

  c->cd(3);
  dummy3->Draw("ape");
  dummy3->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  dummy3->GetXaxis()->SetRangeUser(0, 0.2);
  dummy3->GetXaxis()->SetNdivisions(505);
  dummy3->GetYaxis()->SetRangeUser(0., 8);

  auto*d = new TCanvas("ppMulti", "ppMulti", 1500,1500);
  d->Divide(4,3);
  auto*e = new TCanvas("pLMulti", "pLMulti", 1500,1500);
  e->Divide(4,3);
  auto*f = new TCanvas("LLMulti", "LLMulti", 1500,1500);
  f->Divide(4,3);

  int counter = 0;


  for(int i=0; i<12; ++i) {
    histRE_relK_Lp[i] = (TH1F*)listTP->FindObject(Form("fProtonLambdaRelKMulti_%i", i));
    histME_relK_Lp[i] = (TH1F*)listTP->FindObject(Form("fProtonLambdaRelKMEMulti_%i", i));
    histRE_relK_ALAp[i] = (TH1F*)listTP->FindObject(Form("fAntiProtonAntiLambdaRelKMulti_%i", i));
    histME_relK_ALAp[i] = (TH1F*)listTP->FindObject(Form("fAntiProtonAntiLambdaRelKMEMulti_%i", i));
    histRE_relK_LL[i] = (TH1F*)listTP->FindObject(Form("fLambdaLambdaRelKMulti_%i", i));
    histME_relK_LL[i] = (TH1F*)listTP->FindObject(Form("fLambdaLambdaRelKMEMulti_%i", i));
    histRE_relK_ALAL[i] = (TH1F*)listTP->FindObject(Form("fAntiLambdaAntiLambdaRelKMulti_%i", i));
    histME_relK_ALAL[i] = (TH1F*)listTP->FindObject(Form("fAntiLambdaAntiLambdaRelKMEMulti_%i", i));
    histRE_relK_pp[i] = (TH1F*)listTP->FindObject(Form("fProtonProtonRelKMulti_%i", i));
    histME_relK_pp[i] = (TH1F*)listTP->FindObject(Form("fProtonProtonRelKMEMulti_%i", i));
    histRE_relK_ApAp[i] = (TH1F*)listTP->FindObject(Form("fAntiProtonAntiProtonRelKMulti_%i", i));
    histME_relK_ApAp[i] = (TH1F*)listTP->FindObject(Form("fAntiProtonAntiProtonRelKMEMulti_%i", i));
    hist_CF_Lp_ALAp_exp[0][i] = Calculate_CF(histRE_relK_Lp[i],histME_relK_Lp[i],Form("hist_CF_Lp_exp_%i", i),normleft,normright);
    hist_CF_Lp_ALAp_exp[1][i] = Calculate_CF(histRE_relK_ALAp[i],histME_relK_ALAp[i],Form("hist_CF_ALAp_exp_%i", i),normleft,normright);
    hist_CF_Lp_ALAp_exp[2][i] = add_CF(hist_CF_Lp_ALAp_exp[0][i],hist_CF_Lp_ALAp_exp[1][i],Form("hist_CF_Lp_ALAp_exp_sum_%i", i));
    hist_CF_LL_ALAL_exp[0][i] = Calculate_CF(histRE_relK_LL[i],histME_relK_LL[i],Form("hist_CF_LL_exp_%i", i),normleft,normright);
    hist_CF_LL_ALAL_exp[1][i] = Calculate_CF(histRE_relK_ALAL[i],histME_relK_ALAL[i],Form("hist_CF_LL_exp_%i", i),normleft,normright);
    hist_CF_LL_ALAL_exp[2][i] = add_CF(hist_CF_LL_ALAL_exp[0][i],hist_CF_LL_ALAL_exp[1][i],Form("hist_CF_LL_ALAL_exp_sum_%i", i));
    hist_CF_pp_ApAp_exp[0][i] = Calculate_CF(histRE_relK_pp[i],histME_relK_pp[i],Form("hist_CF_pp_%i", i),normleft,normright);
    hist_CF_pp_ApAp_exp[1][i] = Calculate_CF(histRE_relK_ApAp[i],histME_relK_ApAp[i],Form("hist_CF_ApAp_%i", i),normleft,normright);
    hist_CF_pp_ApAp_exp[2][i] = add_CF(hist_CF_pp_ApAp_exp[0][i],hist_CF_pp_ApAp_exp[1][i],Form("hist_CF_pp_ApAp_exp_sum_%i", i));
    SetStyleHisto(hist_CF_pp_ApAp_exp[2][i], i%8,i%8);
    SetStyleHisto(hist_CF_Lp_ALAp_exp[2][i], i%8,i%8);
    SetStyleHisto(hist_CF_LL_ALAL_exp[2][i], i%8,i%8);

    if(i==0 || i ==4 || i ==6) {
      histRE_relK_Lp_merge[counter] = histRE_relK_Lp[i];
      histME_relK_Lp_merge[counter] = histME_relK_Lp[i];
      histRE_relK_ALAp_merge[counter] = histRE_relK_ALAp[i];
      histME_relK_ALAp_merge[counter] = histME_relK_ALAp[i];
      histRE_relK_LL_merge[counter] = histRE_relK_LL[i];
      histME_relK_LL_merge[counter] = histME_relK_LL[i];
      histRE_relK_ALAL_merge[counter] = histRE_relK_ALAL[i];
      histME_relK_ALAL_merge[counter] = histME_relK_ALAL[i];
      histRE_relK_pp_merge[counter] = histRE_relK_pp[i];
      histME_relK_pp_merge[counter] = histME_relK_pp[i];
      histRE_relK_ApAp_merge[counter] = histRE_relK_ApAp[i];
      histME_relK_ApAp_merge[counter] = histME_relK_ApAp[i];
    }
    else {
      histRE_relK_Lp_merge[counter]->Add(histRE_relK_Lp[i]);
      histME_relK_Lp_merge[counter]->Add(histME_relK_Lp[i]);
      histRE_relK_ALAp_merge[counter]->Add(histRE_relK_ALAp[i]);
      histME_relK_ALAp_merge[counter]->Add(histME_relK_ALAp[i]);
      histRE_relK_LL_merge[counter]->Add(histRE_relK_LL[i]);
      histME_relK_LL_merge[counter]->Add(histME_relK_LL[i]);
      histRE_relK_ALAL_merge[counter]->Add(histRE_relK_ALAL[i]);
      histME_relK_ALAL_merge[counter]->Add(histME_relK_ALAL[i]);
      histRE_relK_pp_merge[counter]->Add(histRE_relK_pp[i]);
      histME_relK_pp_merge[counter]->Add(histME_relK_pp[i]);
      histRE_relK_ApAp_merge[counter]->Add(histRE_relK_ApAp[i]);
      histME_relK_ApAp_merge[counter]->Add(histME_relK_ApAp[i]);
    }
    c->cd(1);
    hist_CF_pp_ApAp_exp[2][i]->Draw("pe same");
    c->cd(2);
    hist_CF_Lp_ALAp_exp[2][i]->Draw("pe same");
    c->cd(3);
    hist_CF_LL_ALAL_exp[2][i]->Draw("pe same");

    d->cd(i+1);
    dummy1->Draw("ape");
    hist_CF_pp_ApAp_exp[2][i]->Draw("pe same");
    TLatex protonCount;
    protonCount.SetNDC(kTRUE);
    protonCount.DrawLatex(0.2, 0.8, Form("p: %.0f", histRE_relK_pp[i]->Integral(histRE_relK_pp[i]->FindBin(0), histRE_relK_pp[i]->FindBin(0.2))));
    protonCount.DrawLatex(0.2, 0.725, Form("#bar{p}: %.0f", histRE_relK_ApAp[i]->Integral(histRE_relK_ApAp[i]->FindBin(0), histRE_relK_ApAp[i]->FindBin(0.2))));

    e->cd(i+1);
    dummy2->Draw("ape");
    hist_CF_Lp_ALAp_exp[2][i]->Draw("pe same");
    protonCount.DrawLatex(0.2, 0.8, Form("p: %.0f", histRE_relK_Lp[i]->Integral(histRE_relK_Lp[i]->FindBin(0), histRE_relK_Lp[i]->FindBin(0.2))));
    protonCount.DrawLatex(0.2, 0.725, Form("#bar{p}: %.0f", histRE_relK_ALAp[i]->Integral(histRE_relK_ALAp[i]->FindBin(0), histRE_relK_ALAp[i]->FindBin(0.2))));

    f->cd(i+1);
    dummy3->Draw("ape");
    hist_CF_LL_ALAL_exp[2][i]->Draw("pe same");
    protonCount.DrawLatex(0.2, 0.8, Form("p: %.0f", histRE_relK_LL[i]->Integral(histRE_relK_LL[i]->FindBin(0), histRE_relK_LL[i]->FindBin(0.2))));
    protonCount.DrawLatex(0.2, 0.725, Form("#bar{p}: %.0f", histRE_relK_ALAL[i]->Integral(histRE_relK_ALAL[i]->FindBin(0), histRE_relK_ALAL[i]->FindBin(0.2))));

    if(i == 3) ++counter;
    if(i == 5) ++counter;
    std::cout << i << " " << counter << "\n";
  }


  auto *g = new TCanvas("can","can",0,0,1500,550);
  g->Divide(3,1);
  for(int i =0; i<3; ++i) {
    g->cd(i+1)->SetRightMargin(right);
    g->cd(i+1)->SetTopMargin(top);
    g->cd(i+1);

    hist_CF_Lp_ALAp_exp_merge[0][i] = Calculate_CF(histRE_relK_Lp_merge[i],histME_relK_Lp_merge[i],Form("hist_CF_Lp_exp_%i", i),normleft,normright);
    hist_CF_Lp_ALAp_exp_merge[1][i] = Calculate_CF(histRE_relK_ALAp_merge[i],histME_relK_ALAp_merge[i],Form("hist_CF_ALAp_exp_%i", i),normleft,normright);
    hist_CF_Lp_ALAp_exp_merge[2][i] = add_CF(hist_CF_Lp_ALAp_exp_merge[0][i],hist_CF_Lp_ALAp_exp_merge[1][i],Form("hist_CF_Lp_ALAp_exp_sum_%i", i));
    hist_CF_LL_ALAL_exp_merge[0][i] = Calculate_CF(histRE_relK_LL_merge[i],histME_relK_LL_merge[i],Form("hist_CF_LL_exp_%i", i),normleft,normright);
    hist_CF_LL_ALAL_exp_merge[1][i] = Calculate_CF(histRE_relK_ALAL_merge[i],histME_relK_ALAL_merge[i],Form("hist_CF_LL_exp_%i", i),normleft,normright);
    hist_CF_LL_ALAL_exp_merge[2][i] = add_CF(hist_CF_LL_ALAL_exp_merge[0][i],hist_CF_LL_ALAL_exp_merge[1][i],Form("hist_CF_LL_ALAL_exp_sum_%i", i));
    hist_CF_pp_ApAp_exp_merge[0][i] = Calculate_CF(histRE_relK_pp_merge[i],histME_relK_pp_merge[i],Form("hist_CF_pp_%i", i),normleft,normright);
    hist_CF_pp_ApAp_exp_merge[1][i] = Calculate_CF(histRE_relK_ApAp_merge[i],histME_relK_ApAp_merge[i],Form("hist_CF_ApAp_%i", i),normleft,normright);
    hist_CF_pp_ApAp_exp_merge[2][i] = add_CF(hist_CF_pp_ApAp_exp_merge[0][i],hist_CF_pp_ApAp_exp_merge[1][i],Form("hist_CF_pp_ApAp_exp_sum_%i", i));
    SetStyleHisto(hist_CF_pp_ApAp_exp_merge[2][i], i%8,i%8);
    SetStyleHisto(hist_CF_Lp_ALAp_exp_merge[2][i], i%8,i%8);
    SetStyleHisto(hist_CF_LL_ALAL_exp_merge[2][i], i%8,i%8);

    g->cd(1);
    hist_CF_pp_ApAp_exp_merge[2][i]->Draw("pe same");
    g->cd(2);
    hist_CF_Lp_ALAp_exp_merge[2][i]->Draw("pe same");
    g->cd(3);
    hist_CF_LL_ALAL_exp_merge[2][i]->Draw("pe same");
  }


}
