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

Float_t YLowMean[numPart] = {0.485, 1.110, 1.110, 1.316, 1.316, 1.664, 1.664};
Float_t YUpMean[numPart] = {0.51, 1.130, 1.130, 1.327, 1.327, 1.68, 1.68};
Float_t YLowSigma[numPart] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002};
Float_t YUpSigma[numPart] = {0.025, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015};
Float_t YLowPurity[numPart] = {0, 0, 0, 0, 0, 0, 0};

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0};

Float_t YLowRatio[numChoice] = {0.95, 0, 0.9, 0.8, 0};
Float_t YUpRatio[numChoice] = {1.05, 1.2, 1.1, 1.2, 2};

void CompareSigmaWidthPurity_MultipleFiles(
    Int_t Chosenfile = 0,
    TString filename = "list_QC.txt",
    TString OutputDir = "../Run3QA/Periods/PbPb2023/LHC23zx_cpass0/",
    Bool_t isMC = 0)
{

  std::vector<std::string> name;
  std::vector<std::string> nameLegend;
  std::ifstream file(Form("%s", filename.Data()));
  std::string remove = "../Run3QA/Periods/PbPb2023/LHC23zx_cpass0/PostProcessing_";
  std::string remove2 = ".root";

  cout << filename.Data() << endl;
  if (file.is_open())
  {
    std::string line;

    while (std::getline(file, line))
    {
      name.push_back(line);
      cout << line << endl;
      size_t pos = line.find(remove);
      if (pos != std::string::npos)
      {
        line.erase(pos, remove.length());
      }
      size_t pos2 = line.find(remove2);
      if (pos2 != std::string::npos)
      {
        line.erase(pos2, remove2.length());
      }
      nameLegend.push_back(line);
    }
    file.close();
  }
  else
  {
    std::cerr << "Unable to open file: " << std::endl;
  }

  const Int_t numFiles = name.size();
  Int_t color[14] = {634, 628, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
  TString Sfilein[numFiles];

  Int_t Choice = 0;
  Int_t ChosenType = -1;
  TString TypeHisto[numChoice] = {"Mean", "Sigma", "Purity", "Yield", "EfficiencyvsPt"};
  TString Spart[numPart] = {"K0S", "Lam", "ALam", "XiMin", "XiPlu", "OmMin", "OmPlu"};
  TString NamePart[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};
  TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "Efficiency"};
  TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";

  TH1F *histo[numPart];
  TH1F *histoDenom;
  TH1F *histoRatio[numPart];
  TCanvas *canvas[numPart];
  TPad *pad1[numPart];
  TPad *pad2[numPart];

  TFile *filein[numFiles];

  cout << "Do you want to compare Mean (=0), Sigma (=1), Purity (=2) or Yield per event (=3) ot efficiency (=4, only MC)?" << endl;
  cin >> Choice;
  if (!isMC && Choice == 4)
    return;
  Bool_t isGen = 0;
  if (isMC && Choice == 3)
  {
    cout << "Do you want the reco yield (0) or gen yield (1) ? " << endl;
    cin >> isGen;
    if (isGen)
      TypeHisto[Choice] = "GeneratedParticles_Rintegrated";
  }
  cout << Choice << " " << TypeHisto[Choice] << endl;
  if (Choice > (numChoice - 1))
  {
    cout << "Option not implemented" << endl;
    return;
  }

  for (Int_t part = 0; part < numPart; part++)
  {
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

    canvas[part] = new TCanvas("canvas" + Spart[part], "canvas" + Spart[part], 1000, 800);
    StyleCanvas(canvas[part], 0.15, 0.05, 0.05, 0.15);
    pad1[part] = new TPad("pad1" + Spart[part], "pad1" + Spart[part], 0, 0.36, 1, 1);
    pad2[part] = new TPad("pad2" + Spart[part], "pad2" + Spart[part], 0, 0.01, 1, 0.35);
    StylePad(pad1[part], 0.15, 0.05, 0.05, 0.01);
    StylePad(pad2[part], 0.15, 0.05, 0.03, 0.2);

    TLegend *legend = new TLegend(0.5, 0.75, 0.8, 0.9);
    legend->AddEntry("", NamePart[part], "");

    TF1 *lineMass = new TF1("lineMass", "pol0", 0, 8);
    lineMass->SetParameter(0, ParticleMassPDG[part]);
    lineMass->SetLineColor(kBlack);
    lineMass->SetLineStyle(7);

    TF1 *lineAt1 = new TF1("lineAt1", "pol0", 0, 8);
    lineAt1->SetParameter(0, 1);
    lineAt1->SetLineColor(kBlack);
    lineAt1->SetLineStyle(7);

    for (Int_t ifile = 0; ifile < numFiles; ifile++)
    {
      filein[ifile] = TFile::Open(Form("%s", name[ifile].c_str()));
      histo[ifile] = (TH1F *)filein[ifile]->Get(inputName);
      if (!histo[ifile])
      {
        cout << "No histo name: " << inputName << " in file" << ifile << endl;
        return;
      }
      histo[ifile]->SetName(inputName + Form("%i", ifile));

      // Denom
      if (ifile == 0)
      {
        histoDenom = (TH1F *)filein[Chosenfile]->Get(inputName);
        histoDenom->SetName("histoDenom");
      }

      // Ratios
      histoRatio[ifile] = (TH1F *)histo[ifile]->Clone(inputName + Form("_Ratio%i", ifile));
      histoRatio[ifile]->Divide(histoDenom);

      for (Int_t b = 1; b <= histoRatio[ifile]->GetNbinsX(); b++)
      {
        // cout << "Num: " << histo[ifile]->GetBinContent(b) << endl;
        // cout << "Denom " << histoDenom->GetBinContent(b) << endl;
        // cout << "Ratio " << histoRatio[ifile]->GetBinContent(b) << endl;
      }

      StyleHisto(histo[ifile], YLow[part], YUp[part], color[ifile], 33, "", TitleY[Choice], "", 0, 0, 0, 1.5, 1.5, 2);
      StyleHisto(histoRatio[ifile], YLowRatio[Choice], YUpRatio[Choice], color[ifile], 33, TitleXPt, Form("Ratio to %s", nameLegend[Chosenfile].c_str()), "", 0, 0, 0, 1.5, 1.5, 2);
      histoRatio[ifile]->GetXaxis()->SetLabelSize(0.08);
      histoRatio[ifile]->GetXaxis()->SetTitleSize(0.08);
      histoRatio[ifile]->GetXaxis()->SetTitleOffset(1.2);
      histoRatio[ifile]->GetYaxis()->SetLabelSize(0.08);
      histoRatio[ifile]->GetYaxis()->SetTitleSize(0.08);
      histoRatio[ifile]->GetYaxis()->SetTitleOffset(0.8);

      legend->AddEntry(histo[ifile], Form("%s", nameLegend[ifile].c_str()), "pl");

      canvas[part]->cd();
      pad1[part]->Draw();
      pad1[part]->cd();
      histo[ifile]->Draw("same");
      if (Choice == 0)
        lineMass->DrawClone("same");

      canvas[part]->cd();
      pad2[part]->Draw();
      pad2[part]->cd();
      histoRatio[ifile]->Draw("same");
      lineAt1->Draw("same");

    } // end loop on files

    canvas[part]->cd();
    pad1[part]->cd();
    legend->Draw("same");

    TString Sfileout = OutputDir + "Compare" + TypeHisto[Choice];
    cout << "Output file: " << Sfileout << endl;
    if (part == 0)
      canvas[part]->SaveAs(Sfileout + ".pdf(");
    else if (part == numPart - 1)
      canvas[part]->SaveAs(Sfileout + ".pdf)");
    else
      canvas[part]->SaveAs(Sfileout + ".pdf");
  }
}
