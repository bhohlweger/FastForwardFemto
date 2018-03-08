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
  gStyle->SetPalette(57);
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
  histo->SetLineWidth(2);
  histo->Sumw2();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHisto(TH2 *histo)
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
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void plotAN(const char* file="~/Results/LHC17p_fast/AnalysisResults.root")
{
  SetStyle();

  TFile* _file0=TFile::Open(file);
  TList *listSP=0;
  TList *listPID=0;
  TList *listTP=0;
  TList *listAliEvent=0;
  _file0->GetObject("/PWGCF_PLFemto_0/SPdir_0",listSP);
  _file0->GetObject("/PWGCF_PLFemto_0/PIDdir_0",listPID);
  _file0->GetObject("/PWGCF_PLFemto_0/TP_dir_0",listTP);
  _file0->GetObject("/PWGCF_PLFemto_0/AliEventCuts_0",listAliEvent);

  // Event variables
  auto* histCutStats = (TH1F*)listAliEvent->FindObject("fCutStats");
  auto* histZvertex = (TH1F*)listAliEvent->FindObject("Vtz_selected");
  auto* histNtracklets = (TH1F*)listSP->FindObject("fNTracklets");

  // PID variables
  auto* protonTPC = (TH2F*)listPID->FindObject("fProtonNSigmaTPC");
  auto* protonCombined = (TH2F*)listPID->FindObject("fProtonNSigmaCombined");

  // Single particle plots
  auto* protonDCAxy = (TH1F*)listSP->FindObject("fProtonDCAxy");
  auto* protonDCAxyCutZ = (TH1F*)listSP->FindObject("fProtonDCAxyCutz");
  auto* protonPt = (TH1F*)listSP->FindObject("fProtonPt");
  auto* protonPhi = (TH1F*)listSP->FindObject("fProtonPhi");
  auto* protonEta = (TH1F*)listSP->FindObject("fProtonEta");
  auto* lambdaPt = (TH1F*)listSP->FindObject("fLambdaPt");
  auto* lambdaPhi = (TH1F*)listSP->FindObject("fLambdaPhi");
  auto* lambdaEta = (TH1F*)listSP->FindObject("fLambdaEta");
  auto* lambdaDaugherDCA = (TH1F*)listSP->FindObject("fLambdaDCADaughterTracks");
  auto* lambdaPosDaughterDCA = (TH1F*)listSP->FindObject("fLambdaDCAPosdaughPrimVertex");
  auto* lambdaNegDaughterDCA = (TH1F*)listSP->FindObject("fLambdaDCANegdaughPrimVertex");
  auto* lambdaTransRadius = (TH1F*)listSP->FindObject("fLambdaTransverseRadius");
  auto* lambdaInvMass = (TH1F*)listSP->FindObject("fInvMassLambdawCuts");
  auto* antilambdaInvMass = (TH1F*)listSP->FindObject("fInvMassAntiLambdawCuts");
  auto* lambdaK0InvMass = (TH1F*)listSP->FindObject("fInvMassMissIDK0s");
  auto* lambdaK0InvMassCut = (TH1F*)listSP->FindObject("fInvMassMissIDK0swCuts");

  // Particle pair plots
  auto* V0V0shared = (TH1F*)listSP->FindObject("fNV0TrackSharing");
  auto* aV0aV0shared = (TH1F*)listSP->FindObject("fNAntiV0TrackSharing");
  auto* V0pshared = (TH1F*)listSP->FindObject("fNV0protonSharedTracks");
  auto* aV0apshared = (TH1F*)listSP->FindObject("fNAntiV0AntiprotonSharedTracks");

  auto *cCutStats = new TCanvas();
  SetStyleHisto(histCutStats, 0, 1);
  histCutStats->Draw();
  cCutStats->Print("ANplot/CutStats.pdf");

  auto *cEventProps= new TCanvas();
  cEventProps->Divide(2,1);
  cEventProps->cd(1);
  cEventProps->cd(1)->SetLogy();
  SetStyleHisto(histNtracklets, 0, 1);
  histNtracklets->Draw("hist");
  cEventProps->cd(2);
  SetStyleHisto(histZvertex, 0, 1);
  histZvertex->Draw("hist");
  cEventProps->Print("ANplot/eventProperties.pdf");

  auto *cPIDproton= new TCanvas();
  SetStyleHisto(protonTPC);
  protonTPC->SetTitle("; #it{p} (GeV/#it{c}); |n#sigma_{TPC}|");
  SetStyleHisto(protonCombined);
  protonCombined->SetTitle("; #it{p} (GeV/#it{c}); n#sigma_{comb} = #sqrt{n#sigma_{TPC}^{2} + n#sigma_{TOF}^{2}}");
  cPIDproton->Divide(2,1);
  cPIDproton->cd(1);
  protonTPC->Draw("colz");
  auto *cutLineA = new TLine(0, 3, 0.75,3);
  cutLineA->SetLineColor(kGray);
  cutLineA->SetLineWidth(2);
  cutLineA->SetLineStyle(2);
  cutLineA->Draw("same");
  cPIDproton->cd(2);
  protonCombined->Draw("colz");
  auto *cutLineB = new TLine(0.75, 3, 4.05,3);
  cutLineB->SetLineColor(kGray);
  cutLineB->SetLineWidth(2);
  cutLineB->SetLineStyle(2);
  cutLineB->Draw("same");
  cPIDproton->Print("ANplot/ProtonPID.pdf");

  auto *cQAproton= new TCanvas("cQAproton", "cQAproton", 1250,1000);
  SetStyleHisto(protonPt, 0, 1);
  protonPt->SetTitle("; #it{p}_{T} (GeV/#it{c}); Entries");
  SetStyleHisto(protonPhi, 0, 1);
  protonPhi->SetTitle("; #phi (rad); Entries");
  SetStyleHisto(protonEta, 0, 1);
  protonEta->SetTitle("; #eta; Entries");
  SetStyleHisto(protonDCAxy, 0, 1);
  protonDCAxy->SetTitle("; DCA_{xy} (cm); Entries");
  SetStyleHisto(protonDCAxyCutZ, 0, 2);
  cQAproton->Divide(2,2);
  cQAproton->cd(1);
  cQAproton->cd(1)->SetLogy();
  protonPt->Draw("hist");
  cQAproton->cd(2);
  protonPhi->Draw("hist");
  cQAproton->cd(3);
  protonEta->Draw("hist");
  cQAproton->cd(4);
  cQAproton->cd(4)->SetLogy();
  protonDCAxy->Draw("hist");
  protonDCAxyCutZ->Draw("hist same");
  auto *cutLineC = new TLine(-0.2, 0, -0.2, protonDCAxy->GetMaximum());
  cutLineC->SetLineColor(kGray+3);
  cutLineC->SetLineWidth(2);
  cutLineC->SetLineStyle(2);
  cutLineC->Draw("same");
  auto *cutLineD = new TLine(0.2, 0, 0.2, protonDCAxy->GetMaximum());
  cutLineD->SetLineColor(kGray+3);
  cutLineD->SetLineWidth(2);
  cutLineD->SetLineStyle(2);
  cutLineD->Draw("same");
  cQAproton->Print("ANplot/ProtonQA.pdf");

  auto *cQAlambda = new TCanvas("cQAlambda", "cQAlambda", 1250,1000);
  SetStyleHisto(lambdaPt, 0, 1);
  lambdaPt->SetTitle("; #it{p}_{T} (GeV/#it{c}); Entries");
  SetStyleHisto(lambdaPhi, 0, 1);
  lambdaPhi->SetTitle("; #phi (rad); Entries");
  SetStyleHisto(lambdaEta, 0, 1);
  lambdaEta->SetTitle("; #eta; Entries");
  cQAlambda->Divide(2,2);
  cQAlambda->cd(1);
  cQAlambda->cd(1)->SetLogy();
  lambdaPt->Draw("hist");
  cQAlambda->cd(2);
  lambdaPhi->Draw("hist");
  cQAlambda->cd(3);
  lambdaEta->Draw("hist");
  cQAlambda->Print("ANplot/LambdaQA.pdf");

  auto *cQAlambdaDaughter = new TCanvas("cQAlambdaDaughter", "cQAlambdaDaughter", 1250,1000);
  SetStyleHisto(lambdaDaugherDCA, 0, 1);
  lambdaDaugherDCA->SetTitle("; DCA(p, #pi) (cm); Entries");
  SetStyleHisto(lambdaPosDaughterDCA, 0, 1);
  lambdaPosDaughterDCA->SetTitle("; DCA(p, PV) (cm); Entries");
  SetStyleHisto(lambdaNegDaughterDCA, 0, 1);
  lambdaNegDaughterDCA->SetTitle("; DCA(#pi, PV) (cm); Entries");
  SetStyleHisto(lambdaTransRadius, 0,1);
  lambdaTransRadius->SetTitle("; #it{r}_{xy} (cm); Entries");
  cQAlambdaDaughter->Divide(2,2);
  cQAlambdaDaughter->cd(1);
  cQAlambdaDaughter->cd(1)->SetLogy();
  lambdaDaugherDCA->Draw("hist");
  cQAlambdaDaughter->cd(2);
  cQAlambdaDaughter->cd(2)->SetLogy();
  lambdaPosDaughterDCA->Draw("hist");
  cQAlambdaDaughter->cd(3);
  cQAlambdaDaughter->cd(3)->SetLogy();
  lambdaNegDaughterDCA->Draw("hist");
  cQAlambdaDaughter->cd(4);
  lambdaTransRadius->Draw("hist");
  cQAlambdaDaughter->Print("ANplot/LambdaDaughterQA.pdf");

  auto* cLambdaMass = new TCanvas();
  SetStyleHisto(lambdaInvMass, 0, 1);
  SetStyleHisto(antilambdaInvMass, 1, 2);
  lambdaInvMass->SetTitle("; M_{p#pi^{-}} (GeV/#it{c}^{2}}); Entries");
  SetStyleHisto(lambdaK0InvMass, 0, 1);
  lambdaK0InvMass->SetTitle("; M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}}); Entries");
  SetStyleHisto(lambdaK0InvMassCut, 0, 2);
  cLambdaMass->Divide(2,1);
  cLambdaMass->cd(1);
  lambdaInvMass->Draw("hist");
  lambdaInvMass->GetXaxis()->SetRangeUser(1.1, 1.13);
  // antilambdaInvMass->Draw("hist same");
  cLambdaMass->cd(2);
  lambdaK0InvMass->Draw("hist");
  lambdaK0InvMass->GetXaxis()->SetRangeUser(0.4, 0.6);
  lambdaK0InvMassCut->Draw("hist same");
  cLambdaMass->Print("ANplot/LambdaInvMass.pdf");

  auto* cShared = new TCanvas("cShared", "cShared", 1250,1000);
  SetStyleHisto(V0V0shared, 0, 1);
  V0V0shared->SetTitle("; # V0 pairs with shared tracks / event; Entries");
  SetStyleHisto(aV0aV0shared, 0, 1);
  aV0aV0shared->SetTitle("; # #bar{V0} pairs with shared tracks / event; Entries");
  SetStyleHisto(V0pshared, 0, 1);
  V0pshared->SetTitle("; # p-V0 pairs with shared tracks / event; Entries");
  SetStyleHisto(aV0apshared, 0,1);
  aV0apshared->SetTitle("; # #bar{p}-#bar{V0} pairs with shared tracks / event; Entries");
  cShared->Divide(2,2);
  cShared->cd(1);
  cShared->cd(1)->SetLogy();
  V0V0shared->Draw("hist");
  cShared->cd(2);
  cShared->cd(2)->SetLogy();
  aV0aV0shared->Draw("hist");
  cShared->cd(3);
  cShared->cd(3)->SetLogy();
  V0pshared->Draw("hist");
  cShared->cd(4);
  cShared->cd(4)->SetLogy();
  aV0apshared->Draw("hist");
  cShared->Print("ANplot/SharedTracks.pdf");
}
