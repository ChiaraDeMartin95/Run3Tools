// This macro provides the comparison of the following variables in two different periods / passes :
// 0) mean of gaussian fits
// 1) sigma of gaussian fits
// 2) purity of selected sample
// 3) raw yields
// It takes in input the output of the post processing macro PostProcessV0AndCascQA_AO2D.C
// ========================================================================
//
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
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TPad.h"
#include "TSpline.h"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRan\
ge,
                Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
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
const Int_t numChoice = 5; // mean, sigma, purity, yield, efficiency for MC
Float_t ParticleMassPDG[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};

// histo0 -> num
// histo1 -> denom
void CompareSigmaWidthPurity(TString year0 = "LHC24l_run553588_skimmed",
                             TString year1 = "LHC24l_run553588",
                             TString yearRatioToPub = "",
                             TString Sfilein0 = "../TriggerForRun3/EventFiltering2024_skimmedDATA/Yields_Omega_LHC24al_run553588_skimmed_OneGaussFit_TrackedCascades.root",
                             TString Sfilein1 = "../TriggerForRun3/EventFiltering2024_skimmedDATA/Yields_Omega_LHC24al_run553588_OneGaussFit_TrackedCascades.root",
                             TString OutputDir = "../TriggerForRun3/EventFiltering2024_skimmedDATA/",
                             Bool_t isPseudoEfficiency = 0,
                             Bool_t isOnlyPseudoEfficiency = 0,
                             TString SPublishedYieldForPseudoEff =
                                 "../PublishedYield13TeV/HEPData-ins1748157-v1-Table", // directory where published yields are stored
                             Bool_t ispp = 1,
                             Bool_t isYieldFromInvMassPostProcess = 1,
                             Bool_t isMC = 0)
{
  // isYieldFromInvMassPostProcess = 1 if files in input are the oputput of the macro Yields_from_invmass.C
  Int_t Choice = 0;
  Int_t ChosenType = -1;
  TString TypeHisto[numChoice] = {"Mean", "Sigma", "Purity", "Yield", "EfficiencyvsPt"};
  TString Spart[numPart] = {"K0S", "Lam", "ALam", "XiMin", "XiPlu", "OmMin", "OmPlu"};
  TString NamePart[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};
  TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "Efficiency"};
  TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";

  Float_t YLowMean[numPart] = {0.485, 1.110, 1.110, 1.316, 1.316, 1.664, 1.664};
  Float_t YUpMean[numPart] = {0.51, 1.130, 1.130, 1.327, 1.327, 1.68, 1.68};
  Float_t YLowSigma[numPart] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002};
  Float_t YUpSigma[numPart] = {0.025, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015};
  Float_t YLowPurity[numPart] = {0, 0, 0, 0, 0, 0, 0};

  Float_t YLow[numPart] = {0};
  Float_t YUp[numPart] = {0};
  Float_t YLowRatio[numChoice] = {0.95, 0.8, 0.9, 0.8, 0};
  Float_t YUpRatio[numChoice] = {1.05, 1.2, 1.1, 1.2, 2};

  Int_t color0 = kRed + 2;
  Int_t color1 = kBlue + 2;

  TH1F *histo0[numPart];
  TH1F *histo1[numPart];
  TH1F *histoRatio[numPart];
  TCanvas *canvas[numPart];
  TPad *pad1[numPart];
  TPad *pad2[numPart];

  TString Sfileout = "";

  if (isOnlyPseudoEfficiency)
    isPseudoEfficiency = 1;
  if (isYieldFromInvMassPostProcess)
  {
    cout << "Do you want to study Xi (=3) or Omega (=5)?" << endl;
    cin >> ChosenType;
  }
  Bool_t isGen = 0;
  if (!isOnlyPseudoEfficiency)
  {
    cout << "Do you want to compare Mean (=0), Sigma (=1), Purity (=2) or Yield per event (=3) ot efficiency (=4, only MC)?" << endl;
    cin >> Choice;
    if (!isMC && Choice == 4)
      return;
    if (isMC && Choice == 3)
    {
      cout << "Do you want the reco yield (0) or gen yield (1) ? " << endl;
      cin >> isGen;
      if (isGen)
        TypeHisto[Choice] = "GeneratedParticles_Rintegrated";
    }
    Bool_t YieldNorm = 0;
    if (Choice == 3)
    {
      cout << "Do you want the yield normalised per event or not?";
      cin >> YieldNorm;
      if (YieldNorm == 1)
        TypeHisto[Choice] = "Yield";
      else
        TypeHisto[Choice] = "YieldNotNorm";
    }
    cout << Choice << " " << TypeHisto[Choice] << endl;
    if (Choice > (numChoice - 1))
    {
      cout << "Option not implemented" << endl;
      return;
    }
    if (isPseudoEfficiency == 1 && Choice != 3)
    {
      cout << "To compute pseudoefficiency you must compare Yields" << endl;
      return;
    }

    Float_t ScalingFactorUnskimmed = 1;
    if (Choice == 3)
    {
      cout << "Do you want to scale the denominator of the ratio? Write down the scalign factor (float)" << endl;
      cout << "Usage case: you are comparing omega yields in skimmed / unskimmed files and you need to multiply the skimmed ratio by the downscaling factor (in 2024 is 0.2)" << endl;
      cin >> ScalingFactorUnskimmed;
    }

    TFile *filein0 = new TFile(Sfilein0, "");
    if (!filein0)
    {
      cout << "No input file n.0" << endl;
      return;
    }
    TFile *filein1 = new TFile(Sfilein1, "");
    if (!filein1)
    {
      cout << "No input file n.1" << endl;
      return;
    }

    // Sfileout = OutputDir + "Compare" + TypeHisto[Choice] + "_" + Spart[ChosenType] + "_" + year0 + "vs" + year1;
    Sfileout = OutputDir + "Compare" + TypeHisto[Choice] + "_";
    if (isYieldFromInvMassPostProcess)
      Sfileout += Spart[ChosenType] + "_";
    Sfileout += year0 + "vs" + year1;
    //Sfileout += "_IsTrackedCascades";
    cout << "Output file: " << Sfileout << endl;

    if (year0 == "LHC22f_pass2")
      color0 = kOrange + 2;

    for (Int_t part = 0; part < numPart; part++)
    {
      if (isYieldFromInvMassPostProcess && part != ChosenType)
        continue;
      cout << "\n\e[35mParticle:\e[39m " << Spart[part] << endl;
      if (Choice == 0)
      {
        YLow[part] = YLowMean[part];
        YUp[part] = YUpMean[part];
      }
      else if (Choice == 1)
      {
        YLow[part] = YLowSigma[part];
        YUp[part] = YUpSigma[part];
      }
      else if (Choice == 2)
      {
        YLow[part] = YLowPurity[part];
        YUp[part] = 1;
      }
      else if (Choice == 4)
      {
        YLow[part] = 0;
        YUp[part] = 1;
      }

      TString inputName = TypeHisto[Choice] + "_" + Spart[part];
      if (isYieldFromInvMassPostProcess)
        inputName = "histo" + TypeHisto[Choice];
      histo0[part] = (TH1F *)filein0->Get(inputName);
      if (!histo0[part])
      {
        cout << "No histo name: " << inputName << " in file0" << endl;
        return;
      }
      histo0[part]->SetName(inputName + "_file0");
      histo1[part] = (TH1F *)filein1->Get(inputName);
      if (!histo1[part])
      {
        cout << "No histo name: " << inputName << " in file1" << endl;
        return;
      }
      histo1[part]->SetName(inputName + "_file1");

      // Ratios
      histo1[part]->Scale(ScalingFactorUnskimmed);
      histoRatio[part] = (TH1F *)histo0[part]->Clone(inputName + "_Ratio");
      if (histo0[part]->GetNbinsX() != histo1[part]->GetNbinsX())
      {
        cout << "The number of bins of the two histograms are different " << endl;
        return;
      }
      histoRatio[part]->Divide(histo1[part]);

      for (Int_t b = 1; b <= histoRatio[part]->GetNbinsX(); b++)
      {
        // cout << "Num: " << histo0[part]->GetBinContent(b) << endl;
        // cout << "Denom " << histo1[part]->GetBinContent(b) << endl;
        // cout << "Ratio " << histoRatio[part]->GetBinContent(b) << endl;
      }

      if (Choice == 3)
      {
        YLow[part] = 0;
        if (histo1[part]->GetBinContent(histo1[part]->GetMaximumBin()) > histo0[part]->GetBinContent(histo1[part]->GetMaximumBin()))
        {
          YUp[part] = 1.2 * histo1[part]->GetBinContent(histo1[part]->GetMaximumBin());
        }
        else
        {
          YUp[part] = 1.2 * histo0[part]->GetBinContent(histo0[part]->GetMaximumBin());
        }
      }
      if (Choice != 3)
      {
        YLow[part] += 10e-5;
        if (Choice != 2)
          YUp[part] -= 10e-5;
      }

      canvas[part] = new TCanvas("canvas" + Spart[part], "canvas" + Spart[part], 1000, 800);
      StyleCanvas(canvas[part], 0.15, 0.05, 0.05, 0.15);
      pad1[part] = new TPad("pad1" + Spart[part], "pad1" + Spart[part], 0, 0.36, 1, 1);
      pad2[part] = new TPad("pad2" + Spart[part], "pad2" + Spart[part], 0, 0.01, 1, 0.35);
      StylePad(pad1[part], 0.15, 0.05, 0.05, 0.01);
      StylePad(pad2[part], 0.15, 0.05, 0.03, 0.2);

      TLegend *legend;
      if (Spart[part] == "XiPlu" && Choice == 0)
        legend = new TLegend(0.5, 0.25, 0.8, 0.45);
      // else if (Spart[part] == "K0S" && Choice == 1)
      // legend = new TLegend(0.5, 0.5, 0.8, 0.7);
      else if (part <= 2 && Choice == 2)
        legend = new TLegend(0.5, 0.25, 0.8, 0.45);
      else
        legend = new TLegend(0.5, 0.75, 0.8, 0.9);
      legend->AddEntry("", NamePart[part], "");

      StyleHisto(histo0[part], YLow[part], YUp[part], color0, 33, "", TitleY[Choice], "", 0, 0, 0, 1.5, 1.5, 2);
      StyleHisto(histo1[part], YLow[part], YUp[part], color1, 33, "", TitleY[Choice], "", 0, 0, 0, 1.5, 1.5, 2);
      StyleHisto(histoRatio[part], YLowRatio[Choice], YUpRatio[Choice], color0, 33, TitleXPt, "Ratio to " + year1, "", 0, 0, 0, 1.5, 1.5, 2);
      histoRatio[part]->GetXaxis()->SetLabelSize(0.08);
      histoRatio[part]->GetXaxis()->SetTitleSize(0.08);
      histoRatio[part]->GetXaxis()->SetTitleOffset(1.2);
      histoRatio[part]->GetYaxis()->SetLabelSize(0.08);
      histoRatio[part]->GetYaxis()->SetTitleSize(0.08);
      histoRatio[part]->GetYaxis()->SetTitleOffset(0.8);

      TF1 *lineMass = new TF1("lineMass", "pol0", 0, 8);
      lineMass->SetParameter(0, ParticleMassPDG[part]);
      lineMass->SetLineColor(kBlack);
      lineMass->SetLineStyle(7);

      TF1 *lineAt1 = new TF1("lineAt1", "pol0", 0, 8);
      lineAt1->SetParameter(0, 1);
      lineAt1->SetLineColor(kBlack);
      lineAt1->SetLineStyle(7);

      canvas[part]->cd();
      pad1[part]->Draw();
      pad1[part]->cd();
      histo0[part]->Draw("same");
      histo1[part]->Draw("same");
      if (Choice == 0)
        lineMass->DrawClone("same");
      legend->AddEntry(histo0[part], year0, "pl");
      legend->AddEntry(histo1[part], year1, "pl");
      legend->Draw("");

      canvas[part]->cd();
      pad2[part]->Draw();
      pad2[part]->cd();
      histoRatio[part]->Draw("same");
      lineAt1->Draw("same");

      canvas[part]->SaveAs(Sfileout + ".pdf");
      canvas[part]->SaveAs(Sfileout + ".png");
    }
  }

  // Pseudoefficiency
  TH1F *histoPub[numPart];
  TH1F *histoPubError[numPart];
  TH1F *histoNum[numPart];
  TH1F *histoRatioToPub[numPart];
  TDirectoryFile *dir;
  TDirectoryFile *dir1;
  TDirectoryFile *dir2;
  TString SfilePub = "";
  TString FileName[numPart] = {"1", "2", "2", "3", "3", "4", "4"};
  TString HistoNumber[numPart] = {"11", "11", "11", "11", "11", "6", "6"};
  TSpline3 *splinePub[numPart];
  TF1 *fsplinePub[numPart];
  TCanvas *canvasRatioToPub[numPart];
  TString SfileoutRatioToPub = OutputDir + "PseudoEfficiencyQA_" + yearRatioToPub + Spart[ChosenType];
  if (!ispp)
    SfileoutRatioToPub = OutputDir + "PseudoEfficiency_";
  // SfileoutRatioToPub += Spart[3];
  Float_t YLowRatioToPub = 0;
  // DATA
  Float_t YUpRatioToPub[numPart] = {0.2, 0.2, 0.2, 0.05, 0.05, 0.05, 0.05};
  // MC
  // Float_t YUpRatioToPub[numPart] = {0.4, 0.3, 0.3, 0.2, 0.2, 0.1, 0.1};
  Float_t LowPub[numPart] = {0, 0.4, 0.4, 0.6, 0.6, 0.9, 0.9}; // lower values of published yields
  Float_t UpPub[numPart] = {10, 6.5, 6.5, 5.5, 5.5, 5.5, 5.5}; // upper value of published yields
  Float_t Low[numPart] = {0, 0.4, 0.4, 0.6, 0.6, 0.9, 0.9};
  Float_t Up[numPart] = {0, 4, 4, 4, 4, 4, 4};

  TFile *fileTemp;
  TString SfileTemp = "";

  if (isPseudoEfficiency)
  {
    cout << "\n\n\e[35mComputation of pseudoefficiency \e[39m" << endl;
    for (Int_t part = 0; part < numPart; part++)
    {
      if (isOnlyPseudoEfficiency) // I take numerator from external file
      {
        SfileTemp = "../Run3QA/Periods/LHC22m_pass3/PostProcess_qa_LHC22m_pass3_relval_cpu2_Train78545.root";
        fileTemp = new TFile(SfileTemp, "");
        if (!fileTemp)
          return;
      }

      cout << "\n\e[35mParticle:\e[39m " << Spart[part] << endl;
      if (isOnlyPseudoEfficiency)
        cout << "File temp: " << SfileTemp << endl;
      if (ispp)
        SfilePub = SPublishedYieldForPseudoEff + "_" + FileName[part] + ".root";
      else
        SfilePub = SPublishedYieldForPseudoEff + ".root";

      TFile *filePub = new TFile(SfilePub, "");
      if (!filePub)
      {
        cout << "File " << SfilePub << " not available " << endl;
        return;
      }
      cout << "Input file with published yields: " << SfilePub << endl;

      if (ispp)
        dir = (TDirectoryFile *)filePub->Get("Table " + FileName[part]);
      else if (!ispp && part == 0)
        dir = (TDirectoryFile *)filePub->Get("Table 4");
      else if (!ispp && part == 3)
      {
        dir = (TDirectoryFile *)filePub->Get("chists_Xi_Tot_variations");
        if (!dir)
          return;
        dir1 = (TDirectoryFile *)dir->Get("Default");
        if (!dir1)
          return;
        dir2 = (TDirectoryFile *)dir1->Get("Average");
        if (!dir2)
          return;
      }
      else
      {
        cout << "Not implemented" << endl;
        return;
      }
      if (!dir)
      {
        cout << "Input dir not available " << endl;
        return;
      }
      if (ispp)
        histoPub[part] = (TH1F *)dir->Get("Hist1D_y" + HistoNumber[part]);
      else if (part == 0)
        histoPub[part] = (TH1F *)dir->Get("Hist1D_y1");
      else if (part == 3)
      {
        histoPub[part] = (TH1F *)dir2->Get("Yield_stat_9");
      }
      else
      {
        cout << "Not implemented" << endl;
        return;
      }
      if (!histoPub[part])
      {
        cout << "Published histo not found" << endl;
        return;
      }
      histoPub[part]->SetName("histoYieldPub" + Spart[part]);

      if (ispp)
        histoPubError[part] = (TH1F *)dir->Get("Hist1D_y" + HistoNumber[part] + "_e1");
      else if (part == 0)
        histoPubError[part] = (TH1F *)dir->Get("Hist1D_y1_e1");
      else if (part == 3)
      {
        histoPubError[part] = (TH1F *)histoPub[part]->Clone("histoYieldPubError" + Spart[part]);
        for (Int_t b = 1; b <= histoPubError[part]->GetNbinsX(); b++)
        {
          histoPubError[part]->SetBinContent(b, histoPubError[part]->GetBinError(b));
        }
      }
      else
      {
        cout << "Not implemented" << endl;
        return;
      }
      if (!histoPubError[part])
      {
        cout << "Published histo of STAT. ERRORS not found" << endl;
        return;
      }
      histoPubError[part]->SetName("histoYieldPubError" + Spart[part]);

      if (part != 0 && !isYieldFromInvMassPostProcess)
        histoPub[part]->Scale(1. / 2); // particle and antiparticle yields are summed
      splinePub[part] = new TSpline3(histoPub[part], "Spline" + Spart[part]);
      sp3 = (TSpline3 *)splinePub[part]->Clone("SplineClone" + Spart[part]);
      fsplinePub[part] = new TF1("fSpline" + Spart[part], spline, 0, 10);

      if (isOnlyPseudoEfficiency)
        if (isYieldFromInvMassPostProcess)
          histoNum[part] = (TH1F *)fileTemp->Get("histoYield");
        else
          histoNum[part] = (TH1F *)fileTemp->Get("Yield_" + Spart[part]);
      else
        histoNum[part] = (TH1F *)histo0[part]->Clone("histoNum" + Spart[part]);

      if (!histoNum[part])
      {
        cout << "histoNum not found " << endl;
        return;
      }
      for (Int_t b = 1; b <= histoPub[part]->GetNbinsX(); b++)
      {
        // cout << fsplinePub[part]->Eval(histoPub[part]->GetBinCenter(b)) << " vs "
        //<< histoPub[part]->GetBinContent(b) << endl;
      }
      histoRatioToPub[part] = (TH1F *)histoNum[part]->Clone("RatioToPub" + Spart[part]);
      for (Int_t b = 1; b <= histoNum[part]->GetNbinsX(); b++)
      {
        Float_t ALow = histoNum[part]->GetXaxis()->GetBinLowEdge(b);
        Float_t AUp = histoNum[part]->GetXaxis()->GetBinUpEdge(b);
        if (fsplinePub[part]->Integral(ALow, AUp) <= 0)
          continue;
        Float_t Numerator = histoNum[part]->GetBinContent(b) * histoNum[part]->GetBinWidth(b);

        Float_t SplineIntegralError = 0;
        if (part != 0)
          histoPubError[part]->Scale(1. / 2);
        for (Int_t b = 1; b <= histoPubError[part]->GetNbinsX(); b++)
        {
          if (histoPubError[part]->GetBinCenter(b) < ALow)
            continue;
          if (histoPubError[part]->GetBinCenter(b) > AUp)
            continue;
          SplineIntegralError += pow(histoPubError[part]->GetBinContent(b), 2);
        }
        SplineIntegralError = sqrt(SplineIntegralError);

        histoRatioToPub[part]->SetBinContent(b, Numerator / fsplinePub[part]->Integral(ALow, AUp));
        histoRatioToPub[part]->SetBinError(b, sqrt(
                                                  pow(histoNum[part]->GetBinError(b) / histoNum[part]->GetBinContent(b), 2) +
                                                  pow(SplineIntegralError / fsplinePub[part]->Integral(ALow, AUp), 2)) *
                                                  histoRatioToPub[part]->GetBinContent(b));

        cout << histoNum[part]->GetBinCenter(b) << " Alow: " << ALow << " AUp " << AUp << endl;
        cout << histoNum[part]->GetBinContent(b) << " +- " << histoNum[part]->GetBinError(b) << endl;
        cout << fsplinePub[part]->Integral(ALow, AUp) << " +- " << SplineIntegralError << endl;
        cout << histoRatioToPub[part]->GetBinContent(b) << " +- " << histoRatioToPub[part]->GetBinError(b) << endl;
      }
      canvasRatioToPub[part] = new TCanvas("canvasRatioToPub" + Spart[part], "canvasRatioToPub" + Spart[part], 1000, 800);
      StyleCanvas(canvasRatioToPub[part], 0.15, 0.05, 0.05, 0.15);

      TLegend *legendRatioToPub;
      if (part == 1 || part == 2)
        legendRatioToPub = new TLegend(0.4, 0.4, 0.8, 0.6);
      else
        legendRatioToPub = new TLegend(0.4, 0.6, 0.8, 0.8);
      legendRatioToPub->AddEntry("", NamePart[part], "");
      // legendRatioToPub->AddEntry("", year0, "");

      StyleHisto(histoRatioToPub[part], YLowRatioToPub, YUpRatioToPub[part], color0, 33, TitleXPt, "Ratio to published yield" + Spart[part], "", 0, 0, 0, 1.5, 1.5, 2);
      if (part == 2 || part == 4 || part == 6)
      {
        histoRatioToPub[part]->SetLineColor(kBlue + 3);
        histoRatioToPub[part]->SetMarkerColor(kBlue + 3);
      }
      histoRatioToPub[part]->GetXaxis()->SetRangeUser(Low[part], Up[part]);
      histoRatioToPub[part]->Draw("same");
      if (part == 2 || part == 4 || part == 6)
      {
        canvasRatioToPub[part - 1]->cd();
        histoRatioToPub[part]->Draw("same");
      }

      // legendRatioToPub->Draw("");
      if (!ispp)
        canvasRatioToPub[part]->SaveAs(SfileoutRatioToPub + ".pdf");
      else
      {
        if (part == 0)
          canvasRatioToPub[part]->SaveAs(SfileoutRatioToPub + ".pdf(");
        else if (part == numPart - 1)
        {
          canvasRatioToPub[part]->SaveAs(SfileoutRatioToPub + ".pdf");
          canvasRatioToPub[part - 1]->SaveAs(SfileoutRatioToPub + ".pdf)");
        }
        else
        {
          canvasRatioToPub[part]->SaveAs(SfileoutRatioToPub + ".pdf");
          if (part == 2 || part == 4)
          {
            canvasRatioToPub[part - 1]->SaveAs(SfileoutRatioToPub + ".pdf");
          }
        }
      }
    }
  }

  if (!isOnlyPseudoEfficiency)
  {
    cout << "\nI started from the files: " << endl;
    cout << Sfilein0 << "\n"
         << Sfilein1 << endl;

    cout << "\nI created the file: " << endl;
    cout << Sfileout << endl;
  }
  if (isPseudoEfficiency)
    cout << "\nPseudoefficiency is stored here: " << SfileoutRatioToPub << endl;
  cout << "Numerator and denominator spectra are obtained from these files:";
  if (isOnlyPseudoEfficiency)
    cout << "\nNumerator: " << SfileTemp << endl;
  else
    cout << "\nNumerator: " << Sfilein0 << endl;
  cout << "Denominator: " << SPublishedYieldForPseudoEff << endl;
}
