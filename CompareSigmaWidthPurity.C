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
const Int_t numChoice = 4; // mean, sigma, purity, yield
Float_t ParticleMassPDG[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};

// histo0 -> num
// histo1 -> denom
void CompareSigmaWidthPurity(TString year0 = "LHC22o_pass2_triggsel" /*"LHC22m_pass2" /*"LHC22f_pass2" /*"LHC22m_pass1"*/,
                             TString year1 = "LHC22m_pass2" /*"LHC22s_PbPb" /*"LHC22f_pass2"  "LHC22m_pass2"*/,
                             TString Sfilein0 = "../Run3QA/LHC22o_pass2/Yields_Xi_LHC22o_triggsel_Train57049.root" /*"../Run3QA/LHC22f_pass2/PostProcessing_Train47130_22f.root" /* "../Run3QA/LHC22m_pass1/PostProcessing_Train44413_22m.root"*/,
                             TString Sfilein1 = "../Run3QA/LHC22m_pass2/Yields_Xi_LHC22m_pass2_Train52781.root" /*"../Run3QA/LHC22f_pass2/PostProcessing_Train47130_22f.root"// "../LHC16k_pass3/PostProcessLHC16k_pass3.root"*/,
                             TString OutputDir = "../Run3QA/LHC22o_pass2/" /*"../Run3QA/LHC22m_pass2/"*/,
                             Bool_t isPseudoEfficiency = 0,
                             Bool_t isOnlyPseudoEfficiency = 0,
                             TString SPublishedYieldForPseudoEff =                 /*"../PbPbYields5TeV/SpectraAlessandro" /*"../PbPbYields5TeV/SpectraK0sPublished276TeV_2040Cent" */
                             "../PublishedYield13TeV/HEPData-ins1748157-v1-Table", // directory where published yields are stored
                             Bool_t ispp = 1,
                             Bool_t isYieldFromInvMassPostProcess = 1)
{
  // isYieldFromInvMassPostProcess = 1 if files in input are the oputput of the macro Yields_from_invmass.C
  Int_t Choice = 0;
  Int_t ChosenType = -1;
  TString TypeHisto[4] = {"Mean", "Sigma", "Purity", "Yield"};
  TString Spart[numPart] = {"K0S", "Lam", "ALam", "XiMin", "XiPlu", "OmMin", "OmPlu"};
  TString NamePart[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};
  TString TitleY[4] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}"};
  TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";

  Float_t YLowMean[numPart] = {0.485, 1.110, 1.110, 1.316, 1.316, 1.664, 1.664};
  Float_t YUpMean[numPart] = {0.51, 1.130, 1.130, 1.327, 1.327, 1.68, 1.68};
  Float_t YLowSigma[numPart] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002};
  Float_t YUpSigma[numPart] = {0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015};
  Float_t YLowPurity[numPart] = {0, 0, 0, 0, 0, 0, 0};

  Float_t YLow[numPart] = {0};
  Float_t YUp[numPart] = {0};
  /* default
  Float_t YLowRatio[4] = {0.95, 0.9, 0.01, 0};
  Float_t YUpRatio[4] = {1.05, 3.0, 1, 0.1};
  */
  Float_t YLowRatio[4] = {0.95, 0.5, 1, 2};
  Float_t YUpRatio[4] = {1.05, 1.5, 3, 20};

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
  if (!isOnlyPseudoEfficiency)
  {
    cout << "Do you want to compare Mean (=0), Sigma (=1), Purity (=2) or Yield per event (=3)?" << endl;
    cin >> Choice;
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

    Sfileout = OutputDir + "Compare" + TypeHisto[Choice] + "_" + year0 + "vs" + year1;

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
        } else {
          YUp[part] = 1.2 * histo0[part]->GetBinContent(histo0[part]->GetMaximumBin());
        }

        if (part > 2)
          // YUpRatio[Choice] = 0.01; 22m_pass1 vs LHC16k
          // YUpRatio[Choice] = 0.2;
          YUpRatio[Choice] = 20;
        else
          YUpRatio[Choice] = 0.1; // for 22m vs 22f comparison
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
      else if (Spart[part] == "K0S" && Choice == 1)
        legend = new TLegend(0.5, 0.5, 0.8, 0.7);
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

      TF1 *lineMass = new TF1("pol0", "pol0", 0, 8);
      lineMass->SetParameter(0, ParticleMassPDG[part]);
      lineMass->SetLineColor(kBlack);
      lineMass->SetLineStyle(7);

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

      if (part == 0)
        canvas[part]->SaveAs(Sfileout + ".pdf(");
      else if (part == numPart - 1)
        canvas[part]->SaveAs(Sfileout + ".pdf)");
      else
        canvas[part]->SaveAs(Sfileout + ".pdf");
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
  TString SfileoutRatioToPub = OutputDir + "PseudoEfficiency_" + year1 + "_";
  if (!ispp)
    SfileoutRatioToPub = OutputDir + "PseudoEfficiency_";
  SfileoutRatioToPub += Spart[3];
  Float_t YLowRatioToPub = 0;
  Float_t YUpRatioToPub = 0.01;
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
      // if (part != 0 && part != 3 && !ispp)
      if (part != 3 && !ispp && isOnlyPseudoEfficiency)
        continue;
      if (part != 3 && ispp && isOnlyPseudoEfficiency)
        continue;
      if (isOnlyPseudoEfficiency) // I take numerator from external file
      {
        if (part == 0)
          SfileTemp = "../Run3QA/LHC22s_PbPb/Yields_K0s_LHC22s_PbPb_FrancescaBin.root";
        else if (part == 3)
          SfileTemp = "../Run3QA/LHC22s_PbPb/Yields_XiNeg_LHC22s_PbPb.root";
        if (ispp && part == 3)
          // SfileTemp = "../Run3QA/LHC22m_pass2/Yields_XiNeg_LHC22m_pass2_Train52782_hasTOF.root";
          // SfileTemp = "../Run3QA/LHC21k6_MC_pp/Yields_Xi_LHC21k6_Train54362.root",
          SfileTemp = "../Run3QA/LHC22o_pass2/Yields_XiNeg_LHC22o_triggsel_Train57049.root",
          fileTemp = new TFile(SfileTemp, "");
        if (!fileTemp)
          return;
      }

      cout << "\n\e[35mParticle:\e[39m " << Spart[part] << endl;
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

      // if (part != 0)
      // histoPub[part]->Scale(1. / 2); // particle and antiparticle yields are summed
      splinePub[part] = new TSpline3(histoPub[part], "Spline" + Spart[part]);
      sp3 = (TSpline3 *)splinePub[part]->Clone("SplineClone" + Spart[part]);
      fsplinePub[part] = new TF1("fSpline" + Spart[part], spline, 0, 10);

      if (isOnlyPseudoEfficiency)
        histoNum[part] = (TH1F *)fileTemp->Get("histoYield");
      else
        histoNum[part] = (TH1F *)histo0[part]->Clone("histoNum" + Spart[part]);

      if (!histoNum[part])
        return;
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
          histoPubError[part]->Scale(1.2);
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

      if (part > 2)
      {
        // YUpRatioToPub = 0.001; //this is OK for 22m_pass1
        // YUpRatioToPub = 0.01; // this is OK for 22f_pass2
        YUpRatioToPub = 0.1; // this is OK for MC
      }
      // YUpRatioToPub = 0.5; //this is OK for Run 2
      StyleHisto(histoRatioToPub[part], YLowRatioToPub, YUpRatioToPub, color0, 33, TitleXPt, "Ratio to published yield", "", 0, 0, 0, 1.5, 1.5, 2);
      histoRatioToPub[part]->GetXaxis()->SetRangeUser(Low[part], Up[part]);
      histoRatioToPub[part]->Draw("");
      // legendRatioToPub->Draw("");
      if (!ispp)
        canvasRatioToPub[part]->SaveAs(SfileoutRatioToPub + ".pdf");
      else
      {
        if (part == 0)
          canvasRatioToPub[part]->SaveAs(SfileoutRatioToPub + ".pdf(");
        else if (part == numPart - 1)
          canvasRatioToPub[part]->SaveAs(SfileoutRatioToPub + ".pdf)");
        else
          canvasRatioToPub[part]->SaveAs(SfileoutRatioToPub + ".pdf");
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
