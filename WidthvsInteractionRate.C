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
const Int_t numChoice = 5; // mean, sigma, purity, S over B, yield
Float_t ParticleMassPDG[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
TString Spart[numPart] = {"K0S", "Lam", "ALam", "XiMin", "XiPlu", "OmMin", "OmPlu"};
TString SpartLeg[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};

const Int_t numIntRate = 4;
const Int_t numRecoTypes = 3; // pass2, pass3, MC
Float_t InteractionRate[numIntRate] = {10, 100, 500, 1000};
TString SIntRate[numIntRate] = {"22q 6/15 kHz", "22r 100 kHz", "22m 500 kHz", "22o 1 MHz"};

Float_t YLowMean[numPart] = {0.485, 1.110, 1.110, 1.316, 1.316, 1.664, 1.664};
Float_t YUpMean[numPart] = {0.51, 1.130, 1.130, 1.327, 1.327, 1.68, 1.68};
Float_t YLowSigma[numPart] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002};
Float_t YUpSigma[numPart] = {0.03, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015};
Float_t YLowPurity[numPart] = {0.4, 0.4, 0.4, 0, 0, 0, 0};
Float_t YUpPurity[numPart] = {1, 1, 1, 1, 1, 1, 1};
Float_t YLowSB[numPart] = {0, 0, 0, 0, 0, 0, 0};
Float_t YUpSB[numPart] = {10, 10, 10, 5, 5, 5, 5};
Float_t YLowYield[4] = {0.95, 0.5, 0.9, 0.8};
Float_t YUpYield[4] = {1.05, 1.5, 1.1, 1.2};

// periods: 22q, 22r, 22m, 22o
// datasets pass2: LHC22q_pass2, LHC22r_pass2_low_rate_subset, LHC22m_pass2_subset, LHC22o_pass2_triggersel
// datasets pass3: LHC22q_pass3, LHC22r_pass3_low_rate_subset, LHC22m_pass3_relval_cpu2, LHC22o_pass3_HIR_small
// runs pass2: 12 runs (all), 529066 + 529067, 523308, 526712
// runs pass3: 12 runs (all), 529066 + 529067, 523308, 526712
// Trains pass2: 60333, 70105, 70287, 58982
// Trains pass3: 64263, 70106, 63492, 70598
// status pass2: done, done, done, done
// status pass3: done, done, done, done
// date of post process pass2: 22/03, 17/03, 19/03, 22/03
// date of post process pass3: 20/03, 17/03, 22/03, 21/03
//NOTE: for the V0s, all topological selections are the same in all periods, in data and in MC

void WidthvsInteractionRate(Int_t PlotType = 0,
                            TString OutputDir = "",
                            Float_t PtValue = 1,
                            Int_t ParticleType = 0)
{

  TString SFileInPass2[numIntRate] = {"Periods/LHC22q_pass2/PostProcess_qa_LHC22q_pass2_Train60333.root",
                                      "Periods/LHC22r_pass2/PostProcess_qa_LHC22r_pass2_low_rate_subset_Train70105.root",
                                      "Periods/LHC22m_pass2/PostProcess_qa_LHC22m_pass2_test_Train70287.root",
                                      "Periods/LHC22o_pass2/PostProcess_qa_LHC22o_pass2_triggersel_Train58982.root"};
  TString SFileInPass3[numIntRate] = {"Periods/LHC22q_pass3/PostProcess_qa_LHC22q_pass3_Train64263.root",
                                      "Periods/LHC22r_pass3/PostProcess_qa_LHC22r_pass3_low_rate_subset_Train70106.root",
                                      "Periods/LHC22m_pass3/PostProcess_qa_LHC22m_pass3_relval_cpu2_Train63492.root",
                                      "Periods/LHC22o_pass3/PostProcess_qa_LHC22o_pass3_HIR_small_Train70598.root"};
  TString SFileInMC[numIntRate] = {"Periods/LHC21k6_MC_pp/PostProcess_qaNew_LHC21k6_Train58981.root",
                                   "",
                                   "",
                                   ""};
  Float_t Sigma[numIntRate] = {0};
  Float_t ErrSigma[numIntRate] = {0};
  Int_t Color[numIntRate] = {881, 628, 418, kBlue};
  Int_t Style[numIntRate] = {33, 20, 21, 27};
  TString SPlotType[numChoice] = {"Mean", "Sigma", "Purity", "SOverB", "Yield"};
  TString TitleY[numChoice] = {"Mean (GeV/c^{2})", "Sigma (GeV/c^{2})", "S/(S+B)", "S/B", "1/N_{evt} dN/dp_{T} [GeV/c]^{-1}"};
  TGraphErrors *gvsIR[numRecoTypes];
  TH1F *histoOut[numRecoTypes];
  TF1 *lineMC = new TF1("pol0", "pol0", 0, 4);
  TH1F *histoIn;
  TFile *fileIn;
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  StyleCanvas(canvas, 0.15, 0.05, 0.05, 0.12);
  TLegend *legend;
  if (PlotType == 0)
    legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  else if (PlotType == 1)
    legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  else if (PlotType == 2)
    legend = new TLegend(0.6, 0.1, 0.9, 0.3);
  else
    legend = new TLegend(0.6, 0.1, 0.9, 0.3);

  TString SPass[numRecoTypes] = {"pass2", "pass3", "MC"};

  Float_t MinPt = 0;
  Float_t MaxPt = 0;

  for (Int_t recoType = 0; recoType < numRecoTypes; recoType++) // pass2, pass3, MC
  {
    histoOut[recoType] = new TH1F(Form("vsIRPass3_%i", recoType), Form("vsIRPass3_%i", recoType), 4, 0, 4);
    for (Int_t i = 0; i < numIntRate; i++)
    { // loop over interaction rates
      // if (recoType != 2 && i != 0)
      // continue;
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
      histoIn = (TH1F *)fileIn->Get(SPlotType[PlotType] + "_" + Spart[ParticleType]);
      Sigma[i] = histoIn->GetBinContent(histoIn->GetXaxis()->FindBin(PtValue + 0.001));
      ErrSigma[i] = histoIn->GetBinError(histoIn->GetXaxis()->FindBin(PtValue + 0.001));
      MinPt = histoIn->GetXaxis()->GetBinLowEdge(histoIn->GetXaxis()->FindBin(PtValue + 0.001));
      MaxPt = histoIn->GetXaxis()->GetBinUpEdge(histoIn->GetXaxis()->FindBin(PtValue + 0.001));
      if (recoType == 2)
        lineMC->SetParameter(0, Sigma[i]);
      histoOut[recoType]->SetBinContent(i + 1, Sigma[i]);
      histoOut[recoType]->SetBinError(i + 1, ErrSigma[i]);
      histoOut[recoType]->GetXaxis()->SetBinLabel(i + 1, SIntRate[i]);
    }

    gvsIR[recoType] = new TGraphErrors(numIntRate, InteractionRate, Sigma, 0, ErrSigma);
    Float_t YLow = 0;
    Float_t YUp = 0;
    if (PlotType == 0)
    {
      YLow = YLowMean[ParticleType];
      YUp = YUpMean[ParticleType];
    }
    else if (PlotType == 1)
    {
      YLow = YLowSigma[ParticleType];
      YUp = YUpSigma[ParticleType];
    }
    else if (PlotType == 2)
    {
      YLow = YLowPurity[ParticleType];
      YUp = 1;
    }
    else 
    {
      YLow = 0;
      YUp = 1.2 * histoOut[recoType]->GetBinContent(histoOut[recoType]->GetMaximumBin());
    }
    StyleHisto(histoOut[recoType], YLow, YUp, Color[recoType], Style[recoType], "IR", TitleY[PlotType], "", 0, 0, 0, 1.2, 1.5, 2);
    if (recoType != 2)
    {
      legend->AddEntry(histoOut[recoType], SPass[recoType], "pl");
      if (PlotType != 4 && PlotType != 3)
        histoOut[recoType]->Draw("same");
    }
    else
    {
      lineMC->SetLineColor(Color[recoType]);
      legend->AddEntry(lineMC, SPass[recoType], "l");
      if (PlotType == 3 || PlotType==4)
      {
        histoOut[0]->GetYaxis()->SetRangeUser(YLow, YUp);
        histoOut[0]->Draw("same");
        histoOut[1]->Draw("same");
      }
      lineMC->Draw("same");
    }
  }

  TLegend *legendParticle = new TLegend(0.2, 0.7, 0.35, 0.9);
  legendParticle->SetMargin(0);
  legendParticle->SetTextSize(0.04);
  legendParticle->AddEntry("", "pp 13 TeV", "");
  legendParticle->AddEntry("", SpartLeg[ParticleType] + Form(", %.1f < p_{T} < %.1f GeV/#it{c}", MinPt, MaxPt), "");
  // legend->Draw("");
  // legendParticle->Draw("");
  canvas->SaveAs(OutputDir + SPlotType[PlotType] + "_" + Spart[ParticleType] + "_vsIR_" + Form("PtInterval%.1f-%.1f", MinPt, MaxPt) + ".pdf");
  canvas->SaveAs(OutputDir + SPlotType[PlotType] + "_" + Spart[ParticleType] + "_vsIR_" + Form("PtInterval%.1f-%.1f", MinPt, MaxPt) + ".png");
}
