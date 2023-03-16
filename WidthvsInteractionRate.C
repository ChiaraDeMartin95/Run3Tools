// This macro was originally written by:
// chiara.de.martin@cern.ch

#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TPad.h"
#include "TSpline.h"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);
}

void StyleCanvas(TCanvas *canvas, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(LMargin);
  canvas->SetRightMargin(RMargin);
  canvas->SetTopMargin(TMargin);
  canvas->SetBottomMargin(BMargin);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  // gStyle->SetPalette(55, 0);
}

void StylePad(TPad *pad, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  pad->SetFillColor(0);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(LMargin);
  pad->SetRightMargin(RMargin);
  pad->SetTopMargin(TMargin);
  pad->SetBottomMargin(BMargin);
}

TSpline3 *sp3;
Double_t spline(Double_t *x, Double_t *p)
{
  Double_t xx = x[0];
  return sp3->Eval(xx);
}

const Int_t numPart = 7;
const Int_t numChoice = 4; // mean, sigma, purity, yield
Float_t ParticleMassPDG[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
TString Spart[numPart] = {"K0S", "Lam", "ALam", "XiMin", "XiPlu", "OmMin", "OmPlu"};

const Int_t numIntRate = 4;
const Int_t numRecoTypes = 3; // pass2, pass3, MC
Float_t InteractionRate[numIntRate] = {10, 100, 500, 1000};
TString SIntRate[numIntRate] = {"22q 6/15 kHz", "22r 100 kHz", "22m 500 kHz", "22o 1 MHz"};
Float_t YLowSigma[numPart] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002};
Float_t YUpSigma[numPart] = {0.018, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015};
// periods: 22q, 22r, 22m, 22o
// datasets pass2: LHC22q_pass2, LHC22r_pass2_low_rate_subset, LHC22m_pass2_subset, LHC22o_pass2_triggersel
// datasets pass3: LHC22q_pass3, LHC22r_pass3_low_rate_subset, LHC22m_pass3_relval_cpu2, TO BE DEFINED
// runs pass2: 12 runs (all), 529066 + 529067,  ??, 526712
// runs pass3: 12 runs (all), 529066 + 529067, 523308, ??
// Trains pass2: 60333, 70105 (ongoing), to be done, 58982
// Trains pass3: 64267, 70106 (ongoing), 63492, to be done
// status pass2: done, ongoing, NA, done
// status pass3: done, ongoing, done, NA

void WidthvsInteractionRate(TString OutputDir = "",
                            Float_t PtValue = 1,
                            Int_t ParticleType = 0)
{

  TString SFileInPass2[numIntRate] = {"LHC22q_pass2/PostProcess_LHC22q_pass2_529039_treno60333.root", "", "", "LHC22o_pass2/PostProcess_qa_LHC22o_pass2_triggersel_Train58982.root"};
  TString SFileInPass3[numIntRate] = {"", "", "LHC22m_pass3/PostProcess_qa_LHC22m_pass3_relval_cpu2_Train63492.root", ""};
  TString SFileInMC[numIntRate] = {"LHC21k6_MC_pp/PostProcess_qa_LHC21k6_Train58981.root", "", "", ""};
  Float_t Sigma[numIntRate] = {0};
  Float_t ErrSigma[numIntRate] = {0};
  Int_t Color[numIntRate] = {881, 628, 418, kBlue};
  Int_t Style[numIntRate] = {33, 20, 21, 27};
  TGraphErrors *gvsIR[numRecoTypes];
  TH1F *histoOut[numRecoTypes];
  TH1F *histoIn;
  TFile *fileIn;
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 500);
  StyleCanvas(canvas, 0.15, 0.05, 0.05, 0.15);
  TLegend *legend = new TLegend(0.75, 0.7, 0.9, 0.9);

  TString SPass[numRecoTypes] = {"pass2", "pass3", "MC"};

  Float_t MinPt = 0;
  Float_t MaxPt = 0;

  for (Int_t recoType = 0; recoType < numRecoTypes; recoType++)
  {
    histoOut[recoType] = new TH1F(Form("vsIRPass3_%i", recoType), Form("vsIRPass3_%i", recoType), 4, 0, 4);
    for (Int_t i = 0; i < numIntRate; i++)
    { // loop over interaction rates
      // if (recoType != 2 && i != 0)
      // continue;
      if (recoType == 0 && (i == 1 || i == 2))
        continue;
      if (recoType == 1 && (i == 0 || i == 1 || i == 3))
        continue;
      if (recoType == 2 && i != 0)
        continue;
      if (recoType == 0)
        fileIn = new TFile("../Run3QA/" + SFileInPass2[i], "");
      else if (recoType == 1)
        fileIn = new TFile("../Run3QA/" + SFileInPass3[i], "");
      else
        fileIn = new TFile("../Run3QA/" + SFileInMC[i], "");
      if (!fileIn)
        return;
      histoIn = (TH1F *)fileIn->Get("Sigma_" + Spart[ParticleType]);
      Sigma[i] = histoIn->GetBinContent(histoIn->GetXaxis()->FindBin(PtValue + 0.001));
      ErrSigma[i] = histoIn->GetBinError(histoIn->GetXaxis()->FindBin(PtValue + 0.001));
      MinPt = histoIn->GetXaxis()->GetBinLowEdge(histoIn->GetXaxis()->FindBin(PtValue + 0.001));
      MaxPt = histoIn->GetXaxis()->GetBinUpEdge(histoIn->GetXaxis()->FindBin(PtValue + 0.001));
      histoOut[recoType]->SetBinContent(i + 1, Sigma[i]);
      histoOut[recoType]->SetBinError(i + 1, ErrSigma[i]);
      histoOut[recoType]->GetXaxis()->SetBinLabel(i + 1, SIntRate[i]);
    }

    gvsIR[recoType] = new TGraphErrors(numIntRate, InteractionRate, Sigma, 0, ErrSigma);
    StyleHisto(histoOut[recoType], YLowSigma[ParticleType], YUpSigma[ParticleType], Color[recoType], Style[recoType], "IR", "Sigma (GeV/c^{2})", "", 0, 0, 0, 1.5, 1.5, 2);
    histoOut[recoType]->Draw("same");
    legend->AddEntry(histoOut[recoType], SPass[recoType], "pl");
  }

  TLegend *legendParticle = new TLegend(0.2, 0.7, 0.35, 0.9);
  legendParticle->SetMargin(0);
  legendParticle->SetTextSize(0.04);
  legendParticle->AddEntry("", "pp 13 TeV", "");
  legendParticle->AddEntry("", Spart[ParticleType] + Form(" %.1f < p_{T} < %.1f ", MinPt, MaxPt), "");
  legend->Draw("");
  legendParticle->Draw("");
  canvas->SaveAs(OutputDir + "Sigma" + Spart[ParticleType] + "_vsIR_ " + Form("PtInterval%.1f-%.1f", MinPt, MaxPt) + ".pdf");
}
