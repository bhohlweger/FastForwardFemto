// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette=false, bool title=true)
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
gStyle->SetPalette(57);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHistoMany(TH2 *histo)
{
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.055);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleSize(0.055);
  histo->GetYaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetTitleOffset(1.25);
}

void plotMC(const char *simfile = "~/Results/LHC17l3b/AnalysisResults.root")
{
  SetStyle();
  TFile* _file0sim=TFile::Open(simfile);
  TList *listTPsim=0;
  _file0sim->GetObject("/PWGCF_PLFemto_0/TP_dir_0",listTPsim);
  if(!listTPsim) _file0sim->GetObject("/PWGCF_PLFemto_0/TPdir_0",listTPsim);
  TH2F* fDeltaEtaDeltaPhiProton[9];
  TH2F* fDeltaEtaDeltaPhiProtonME[9];
  TH2F* fDeltaEtaDeltaPhiProtonV0proton[9];
  TH2F* fDeltaEtaDeltaPhiProtonV0protonME[9];
  TH2F* fDeltaEtaDeltaPhiProtonV0pion[9];
  TH2F* fDeltaEtaDeltaPhiProtonV0pionME[9];

  std::vector<int> r = {85,105,125,145,165,185,205,225,245};

  auto* c = new TCanvas("pp", "pp", 1500, 1000);
  c->Divide(3,3);
  for(int i=0; i<9; ++i) {
    c->cd(i+1);
    c->cd(i+1)->SetLogz();
    fDeltaEtaDeltaPhiProton[i] = (TH2F*)listTPsim->FindObject(Form("fDeltaEtaDeltaPhiTPCradpp%i", i));
    fDeltaEtaDeltaPhiProtonME[i] = (TH2F*)listTPsim->FindObject(Form("fDeltaEtaDeltaPhiTPCradppME%i", i));
    fDeltaEtaDeltaPhiProton[i]->Divide(fDeltaEtaDeltaPhiProtonME[i]);
    SetStyleHistoMany(fDeltaEtaDeltaPhiProton[i]);
    fDeltaEtaDeltaPhiProton[i]->Draw("colz");
    fDeltaEtaDeltaPhiProton[i]->SetTitle(Form("p-p -- TPC radius #it{r} = %i cm; |#Delta#eta|; |#Delta#phi*| (rad)", r[i]));
    fDeltaEtaDeltaPhiProton[i]->GetXaxis()->SetRangeUser(0, 0.5);
    fDeltaEtaDeltaPhiProton[i]->GetYaxis()->SetRangeUser(0, 0.5);
  }

  auto* d = new TCanvas("pp(L)", "pp(L)", 1500, 1000);
  d->Divide(3,3);
  for(int i=0; i<9; ++i) {
    d->cd(i+1);
    d->cd(i+1)->SetLogz();
    fDeltaEtaDeltaPhiProtonV0proton[i] = (TH2F*)listTPsim->FindObject(Form("fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProton%i", i));
    fDeltaEtaDeltaPhiProtonV0protonME[i] = (TH2F*)listTPsim->FindObject(Form("fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughProtonME%i", i));
    fDeltaEtaDeltaPhiProtonV0proton[i]->Divide(fDeltaEtaDeltaPhiProtonV0protonME[i]);
    SetStyleHistoMany(fDeltaEtaDeltaPhiProtonV0proton[i]);
    fDeltaEtaDeltaPhiProtonV0proton[i]->Draw("colz");
    fDeltaEtaDeltaPhiProtonV0proton[i]->SetTitle(Form("p-p_{#Lambda} -- TPC radius #it{r} = %i cm; |#Delta#eta|; |#Delta#phi*| (rad)", r[i]));
    fDeltaEtaDeltaPhiProtonV0proton[i]->GetXaxis()->SetRangeUser(0, 0.5);
    fDeltaEtaDeltaPhiProtonV0proton[i]->GetYaxis()->SetRangeUser(0, 0.5);
  }

  auto* e = new TCanvas("ppi(L)", "ppi(L)", 1500, 1000);
  e->Divide(3,3);
  for(int i=0; i<9; ++i) {
    e->cd(i+1);
    e->cd(i+1)->SetLogz();
    fDeltaEtaDeltaPhiProtonV0pion[i] = (TH2F*)listTPsim->FindObject(Form("fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPion%i", i));
    fDeltaEtaDeltaPhiProtonV0pionME[i] = (TH2F*)listTPsim->FindObject(Form("fDeltaEtaDeltaPhiTPCradLpPrimProtonDaughPionME%i", i));
    fDeltaEtaDeltaPhiProtonV0pion[i]->Divide(fDeltaEtaDeltaPhiProtonV0pionME[i]);
    SetStyleHistoMany(fDeltaEtaDeltaPhiProtonV0pion[i]);
    fDeltaEtaDeltaPhiProtonV0pion[i]->Draw("colz");
    fDeltaEtaDeltaPhiProtonV0pion[i]->SetTitle(Form("p-#pi_{#Lambda} -- TPC radius #it{r} = %i cm; |#Delta#eta|; |#Delta#phi*| (rad)", r[i]));
    fDeltaEtaDeltaPhiProtonV0pion[i]->GetXaxis()->SetRangeUser(0, 0.5);
    fDeltaEtaDeltaPhiProtonV0pion[i]->GetYaxis()->SetRangeUser(0, 0.5);
  }
}
