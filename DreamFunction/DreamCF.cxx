/*
 * DreamCF.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#include "DreamCF.h"
#include "TCanvas.h"
#include <iostream>
DreamCF::DreamCF()
:fPair(nullptr)
,fPairShifted(nullptr)
,fPairLast(nullptr)
{
}

DreamCF::~DreamCF()
{
}

void DreamCF::ShiftForEmpty() {
  fPairShifted=new DreamPair("fted");
  if (!fPairLast)
  {
    fPairLast=fPairShifted;
  }
  TH1F* SE=fPair->GetSEDist();
  TH1F* ME=fPair->GetMEDist();

  TH2F* SEMult=fPair->GetSEMultDist();
  TH2F* MEMult=fPair->GetMEMultDist();

  int SEkMinBin=SE->FindFirstBinAbove(0);
  //+1 since we start counting bins at 1
  int nBinsEffSE=SE->GetNbinsX()-SEkMinBin+1;
  float SEkMin=SE->GetXaxis()->GetBinLowEdge(SEkMinBin);
  float SEkMax=SE->GetXaxis()->GetBinUpEdge(SE->GetNbinsX());
  int multBins=SEMult->GetNbinsY();
  int multMax=SEMult->GetYaxis()->GetBinUpEdge(multBins);

  const char* SEHistName=Form("%sShi",SE->GetName());
  TH1F* SEShifted=new TH1F(SEHistName,SEHistName,nBinsEffSE,SEkMin,SEkMax);
  SEShifted->Sumw2();
  const char* MEHistName=Form("%sShi",ME->GetName());
  TH1F* MEShifted=new TH1F(MEHistName,MEHistName,nBinsEffSE,SEkMin,SEkMax);
  MEShifted->Sumw2();
  const char* SEMultHistName=Form("%sShi",SEMult->GetName());
  TH2F* SEMultShifted=new TH2F(SEMultHistName,SEMultHistName,
                               nBinsEffSE,SEkMin,SEkMax,
                               multBins,0,multMax);
  SEMultShifted->Sumw2();
  const char* MEMultHistName=Form("%sShi",MEMult->GetName());
  TH2F* MEMultShifted=new TH2F(MEMultHistName,MEMultHistName,
                               nBinsEffSE,SEkMin,SEkMax,
                               multBins,0,multMax);
  MEMultShifted->Sumw2();
  int ckBin=1;
  for (int ikBin=SEkMinBin;ikBin<=SE->GetNbinsX();++ikBin) {
    SEShifted->SetBinContent(ckBin,SE->GetBinContent(ikBin));
    SEShifted->SetBinError(ckBin,SE->GetBinError(ikBin));

    MEShifted->SetBinContent(ckBin,ME->GetBinContent(ikBin));
    MEShifted->SetBinError(ckBin,ME->GetBinError(ikBin));
//    for (int iMult=1;iMul)

    ckBin++;
  }
  fPairShifted->SetSEDist(SEShifted);
  fPairShifted->SetSEMultDist(SEMultShifted);
  fPairShifted->SetMEDist(MEShifted);
  fPairShifted->SetMEMultDist(MEMultShifted);
//  const char* canName=Form("c%s",SE->GetName());
//  TCanvas *c1 = new TCanvas(canName,canName,2000,1000);
//  c1->Divide(2,1);
//  c1->cd(1);
//  SE->Draw();
//  SEShifted->SetLineColor(2);
//  SEShifted->Draw("same");
//  c1->cd(2);
//  ME->Draw();
//  MEShifted->SetLineColor(2);
//  MEShifted->Draw("same");
}
