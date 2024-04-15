#include "Riostream.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
#include <TLine.h>
#include <TSpline.h>
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
// #include "Constants.h"
#include "ErrRatioCorr.C"
// #include "InputVar.h"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(1.5);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->SetTitle(title);
}

void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset); // 1.2
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

void SetFont(TH1F *histo)
{
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetLabelFont(43);
}
void SetTickLength(TH1F *histo, Float_t TickLengthX, Float_t TickLengthY)
{
  histo->GetXaxis()->SetTickLength(TickLengthX);
  histo->GetYaxis()->SetTickLength(TickLengthY);
}

void SetHistoTextSize(TH1F *histo, Float_t XSize, Float_t XLabelSize, Float_t XOffset, Float_t XLabelOffset, Float_t YSize, Float_t YLabelSize, Float_t YOffset, Float_t YLabelOffset)
{
  histo->GetXaxis()->SetTitleSize(XSize);
  histo->GetXaxis()->SetLabelSize(XLabelSize);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetXaxis()->SetLabelOffset(XLabelOffset);
  histo->GetYaxis()->SetTitleSize(YSize);
  histo->GetYaxis()->SetLabelSize(YLabelSize);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetYaxis()->SetLabelOffset(YLabelOffset);
}

void StyleCanvas(TCanvas *canvas, Float_t TopMargin, Float_t BottomMargin, Float_t LeftMargin, Float_t RightMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(TopMargin);
  gPad->SetLeftMargin(LeftMargin);
  gPad->SetBottomMargin(BottomMargin);
  gPad->SetRightMargin(RightMargin);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
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

// take spectra in input
// produces ratio of spectra wrt 0-100% multiplciity class

const Int_t numPart = 9;
const Int_t numChoice = 8;

Int_t ColorMult[] = {634, 628, 807, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
Float_t SizeMult[] = {2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8};
Float_t SizeMultRatio[] = {1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8};
// Int_t MarkerMult[] = {20, 21, 33, 34, 29, 24, 27, 28, 25, 25, 25, 20, 21, 33};
Int_t MarkerMult[] = {20, 21, 33, 34, 29, 20, 21, 33, 34, 29, 20, 21, 33, 34, 29};

const Float_t massParticle[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.32171, 1.67245, 1.67245, 1.67245};
TString Spart[numPart] = {"K0s", "Lambda", "AntiLambda", "XiNeg", "XiPos", "Xi", "OmegaNeg", "OmegaPos", "Omega"};
TString SpartType[numPart] = {"K0s", "Lambda", "Lambda", "Xi", "Xi", "Xi", "Omega", "Omega", "Omega"};
TString NamePart[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Xi^{#pm}", "#Omega^{-}", "#Omega^{+}", "#Omega^{#pm}"};

TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};

TString SIsBkgParab[2] = {"_BkgRetta", "_BkgParab"};

TString TypeHisto[numChoice] = {"Mean", "Sigma", "Purity", "Yield", "Significance", "Efficiency", "YieldCorr", "hNSigmaSummaryGaus_hTPCNsigma"};
TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "Yield/#sigma_{Yield}", "Efficiency x acceptance", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "NSigmaTPC"};

TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";
TString TitleYYield = "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}";
TString TitleYYieldPtInt = "1/#it{N}_{evt} d#it{N}/d#it{y}";
TString TitleYYieldPtIntToMB = "(1/#it{N}_{evt} d#it{N}/d#it{y}) / (1/#it{N}_{evt} d#it{N}/d#it{y})_{MB}";
TString TitleXMult = "#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
TString TitleXMultToMB = "#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5} / #LTd#it{N}_{ch}/d#it{#eta}#GT^{MB}_{|#it{#eta}|<0.5}";

Float_t YLowMean[numPart] = {0.485, 1.110, 1.110, 1.316, 1.316, 1.316, 1.66, 1.66, 1.66};
Float_t YUpMean[numPart] = {0.51, 1.130, 1.130, 1.327, 1.327, 1.327, 1.68, 1.68, 1.68};
Float_t YLowSigma[numPart] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0, 0.0, 0.0};
Float_t YUpSigma[numPart] = {0.025, 0.015, 0.015, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008};
Float_t YLowPurity[numPart] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0};

Float_t YLowRatio[numChoice] = {0.99, 0.4, 0.5, 0, 0.8, 0.8, 0.8, 0.8};
Float_t YUpRatio[numChoice] = {1.01, 1.6, 1.5, 2, 1.2, 1.2, 1.2, 1.2};

const Int_t numPeriods = 7;
// TString Srun[numPeriods] = {"LHC23f", "LHC23i", "LHC23g", "LHC23j", "LHC23k", "LHC23l", "LHC23m", "LHC23n", "LHC23o", "LHC23q", "LHC23r", "LHC23y", "LHC23z"};
// TString Srun[numPeriods] = {"LHC23p", "LHC23s", "LHC23t", "LHC23u", "LHC23za", "LHC23zc", "LHC23zj", "LHC23zl", "LHC23l"};
// TString Srun[numPeriods] = {"LHC23v", "LHC23ze", "LHC23zf", "LHC23zg", "LHC23zh", "LHC23zi", "LHC23zk", "LHC23zm", "LHC23zn","LHC23zq","LHC23zr","LHC23zt", "LHC23l"};

// 500 kHz periods (11) N.B. zs contains too few events
// TString Srun[numPeriods] = {"LHC23l", "LHC23p", "LHC23v", "LHC23za", "LHC23zc", "LHC23zj", "LHC23zl","LHC23zt"};
// low IR periods
// TString Srun[numPeriods] = {"LHC23i", "LHC23o", "LHC23zg", "LHC23zm", "LHC23n", "LHC23ze", "LHC23zf", "LHC23zk"};
// TString IRLegend[numPeriods] = {"50 kHz", "50 kHz", "50 kHz", "50 kHz", "10 kHz", "23 kHz", "10 kHz", "10 kHz"};
// Mixed sample
// TString Srun[numPeriods] = {"LHC23l", "LHC23o", "LHC23i", "LHC23zk", "LHC23zf", "LHC23m", "LHC23t" };
// TString IRLegend[numPeriods] = {"500 kHz", "50 kHz", "50 kHz", "10 kHz", "10 kHz", "250 kHz", "1000 kHz"};
TString Srun[numPeriods] = {"537505", "537511", "537531", "537546", "537547", "537551", "537553"};
TString IRLegend[numPeriods] = {""};

void MeanSigmaPurityQAFilter(Int_t part = 8,
                             Int_t ChosenRun = 0,
                             Int_t Choice = 0,
                             TString OutputDir = "../TriggerForRun3/EventFiltering2023/LHC23v_pass4/",
                             TString year = "")
{

  gStyle->SetOptStat(0);
  if (ChosenRun > (numPeriods - 1))
  {
    cout << "Chosen Mult outside of available range" << endl;
    return;
  }
  cout << Choice << " " << TypeHisto[Choice] << endl;
  if (Choice > (numChoice - 1))
  {
    cout << "Option not implemented" << endl;
    return;
  }
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
  else if (Choice == 3)
  {
    if (part == 6 || part == 7 || part == 8)
    {
      YUp[part] = 1e-3;
      YLow[part] = 1e-9;
    }
    else if (part == 3 || part == 4 || part == 5)
    {
      YUp[part] = 1e-2;
      YLow[part] = 1e-7;
    }
  }
  else if (Choice == 7)
  {
    YLow[part] = -5;
    YUp[part] = 5;
    TitleXPt = "p(GeV/c)";
  }

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TFile *fileIn[numPeriods];

  // fileout name
  Float_t MaxTPC = 1.2;
  Int_t dauparticle = 0;
  TString TypeHistoSuffix = "";
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "PlotRatios_" + year;
  stringout += "_" + TypeHisto[Choice];
  if (Choice == 7)
  {
    cout << "which particle?" << endl;
    if (part == 5)
      cout << "0: V0PiPlus, 1:V0PiMinus, 2:V0Proton, 3:V0AntiProton " << endl;
    else if (part == 8)
      cout << "0: BachKaMinus, 1:BachKaPlus" << endl;
    cin >> dauparticle;
    if (part == 5)
    {
      if (dauparticle == 0)
        TypeHistoSuffix = "XiV0PiPlus";
      else if (dauparticle == 1)
        TypeHistoSuffix = "XiV0PiMinus";
      else if (dauparticle == 2)
        TypeHistoSuffix = "XiV0Proton";
      else if (dauparticle == 3)
        TypeHistoSuffix = "XiV0AntiProton";
    }
    else if (part == 8)
    {
      MaxTPC = 0.7;
      if (dauparticle == 0)
        TypeHistoSuffix = "OmegaBachKaMinus";
      else if (dauparticle == 1)
        TypeHistoSuffix = "OmegaBachKaPlus";
    }
    TypeHisto[Choice] += TypeHistoSuffix;
    stringout += TypeHistoSuffix;
  }
  stringout += "_" + Spart[part];
  stringoutpdf = stringout;
  stringout += ".root";
  TFile *fileout = new TFile(stringout, "RECREATE");

  // canvases
  TCanvas *canvasPtSpectra = new TCanvas("canvasPtSpectra", "canvasPtSpectra", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  if (Choice == 7)
    LLUpperPad = 0.;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.); // L, R, T, B
  if (Choice == 7)
    StylePad(pad1, 0.18, 0.02, 0.03, 0.15);
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  TH1F *fHistSpectrum[numPeriods];
  TH1F *fHistSpectrumScaled[numPeriods];
  TString sScaleFactorFinal[numPeriods];
  TH1F *fHistSpectrumMultRatio[numPeriods];

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendAllMult = new TLegend(0.22, 0.03, 0.9, 0.28);
  if (Choice==7) legendAllMult = new TLegend(0.22, 0.2, 0.9, 0.4);
  legendAllMult->SetHeader("Runs");
  legendAllMult->SetNColumns(3);
  legendAllMult->SetFillStyle(0);
  TLegendEntry *lheaderAllMult = (TLegendEntry *)legendAllMult->GetListOfPrimitives()->First();
  lheaderAllMult->SetTextSize(0.04);

  TLegend *LegendTitle;
  if (part == 5 && Choice == 2)
    LegendTitle = new TLegend(0.54, 0.55, 0.95, 0.72);
  else
    LegendTitle = new TLegend(0.54, 0.75, 0.95, 0.92);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.04);
  // LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  LegendTitle->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  LegendTitle->AddEntry("", NamePart[part] + ", |y| < 0.5", "");

  TLine *lineat1Mult = new TLine(0, 1, 4, 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  TLine *lineat08Mult = new TLine(0, 0.8, 4, 0.8);
  lineat08Mult->SetLineColor(1);
  lineat08Mult->SetLineStyle(2);

  TLine *lineat12Mult = new TLine(0, 1.2, 4, 1.2);
  lineat12Mult->SetLineColor(1);
  lineat12Mult->SetLineStyle(2);

  // get spectra in multiplicity classes
  for (Int_t m = numPeriods - 1; m >= 0; m--)
  {
    PathIn = "../TriggerForRun3/EventFiltering2023/LHC23v_pass4/Yields_";
    PathIn += Spart[part];
    PathIn += "_" + year;
    // PathIn += "_Run" + Srun[m];
    PathIn += Srun[m];
    PathIn += "_OneGaussFit";
    PathIn += ".root";
    if (Choice == 7)
    {
      PathIn = "../TriggerForRun3/EventFiltering2023/LHC23v_pass4/FilterPostProcessing_";
      PathIn += Srun[m] + ".root";
    }

    cout << "Path in : " << PathIn << endl;

    fileIn[m] = TFile::Open(PathIn);
    if (Choice == 7)
      fHistSpectrum[m] = (TH1F *)fileIn[m]->Get(TypeHisto[Choice]);
    else
      fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histo" + TypeHisto[Choice]);
    fHistSpectrum[m]->SetName("histoSpectrum_" + Srun[m]);
    if (!fHistSpectrum[m])
    {
      cout << " no hist " << endl;
      return;
    }
  } // end loop on mult

  // draw spectra in multiplicity classes
  Float_t xTitle = 15;
  Float_t xOffset = 4;
  Float_t yTitle = 30;
  Float_t yOffset = 2;

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.05;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.042;

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 8);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasPtSpectra->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXPt, TitleY[Choice], "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, 4);
  if (Choice == 7)
    hDummy->GetXaxis()->SetRangeUser(0, MaxTPC);
  pad1->Draw();
  pad1->cd();
  if (Choice == 3)
    gPad->SetLogy();
  hDummy->Draw("same");

  for (Int_t m = numPeriods - 1; m >= 0; m--)
  {
    cout << "Period: " << m << " " << Srun[m] << endl;
    if (Srun[m] == "LHC23p" && part == 8)
      continue;
    fHistSpectrumScaled[m] = (TH1F *)fHistSpectrum[m]->Clone("fHistSpectrumScaled_" + Srun[m]);
    for (Int_t b = 1; b <= fHistSpectrum[m]->GetNbinsX(); b++)
    {
      cout << "bin " << b << " " << fHistSpectrum[m]->GetBinContent(b) << "+-" << fHistSpectrum[m]->GetBinError(b) << endl;
      cout << "bin " << b << " " << fHistSpectrumScaled[m]->GetBinContent(b) << "+-" << fHistSpectrumScaled[m]->GetBinError(b) << endl;
    }
    SetFont(fHistSpectrumScaled[m]);
    fHistSpectrumScaled[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumScaled[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumScaled[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumScaled[m]->SetMarkerSize(0.6 * SizeMult[m]);
    fHistSpectrumScaled[m]->GetYaxis()->SetRangeUser(YLow[part], YUp[part]);
    fHistSpectrumScaled[m]->Draw("same e0x0");
    sScaleFactorFinal[m] = "";
    legendAllMult->AddEntry(fHistSpectrumScaled[m], Srun[m] + sScaleFactorFinal[m] + " " + IRLegend[m] + " ", "pef");
  } // end loop on mult
  LegendTitle->Draw("");
  legendAllMult->Draw("");

  // Compute and draw spectra ratios
  Float_t LimSupMultRatio = 5.1;
  Float_t LimInfMultRatio = 1e-2;
  Float_t YoffsetSpectraRatio = 1.1;
  Float_t xTitleR = 35;
  Float_t xOffsetR = 1;
  Float_t yTitleR = 30;
  Float_t yOffsetR = 2;

  Float_t xLabelR = 25;
  Float_t yLabelR = 25;
  Float_t xLabelOffsetR = 0.02;
  Float_t yLabelOffsetR = 0.04;

  TString TitleYSpectraRatio = "Ratio to " + Srun[ChosenRun];
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, 0, 8);
  for (Int_t i = 1; i <= hDummyRatio->GetNbinsX(); i++)
    hDummyRatio->SetBinContent(i, 1e-12);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, YLowRatio[Choice], YUpRatio[Choice], 1, 1, TitleXPt, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickX, tickY);
  hDummyRatio->GetXaxis()->SetRangeUser(0, 4);
  if (Choice == 7)
    hDummyRatio->GetXaxis()->SetRangeUser(0, MaxTPC);
  if (Choice != 7)
  {
    canvasPtSpectra->cd();
    padL1->Draw();
    padL1->cd();
    hDummyRatio->Draw("same");
  }

  for (Int_t m = numPeriods - 1; m >= 0; m--)
  {
    if (Srun[m] == "LHC23p" && part == 8)
      continue;
    fHistSpectrumMultRatio[m] = (TH1F *)fHistSpectrum[m]->Clone("fHistSpectrumMultRatio_" + Srun[m]);
    fHistSpectrumMultRatio[m]->Divide(fHistSpectrum[ChosenRun]);
    ErrRatioCorr(fHistSpectrum[m], fHistSpectrum[ChosenRun], fHistSpectrumMultRatio[m], 0);
    for (Int_t b = 1; b <= fHistSpectrum[m]->GetNbinsX(); b++)
    {
      cout << "binR " << b << " " << fHistSpectrum[m]->GetBinContent(b) << endl;
      cout << "binR " << b << " " << fHistSpectrum[ChosenRun]->GetBinContent(b) << endl;
      cout << "binR " << b << " " << fHistSpectrumMultRatio[m]->GetBinContent(b) << endl;
    }
    fHistSpectrumMultRatio[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumMultRatio[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumMultRatio[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumMultRatio[m]->SetMarkerSize(SizeMultRatio[m]);

    if (m != ChosenRun && Choice != 7)
    {
      fHistSpectrumMultRatio[m]->Draw("same e0x0");
      lineat1Mult->Draw("same");
      lineat08Mult->Draw("same");
      lineat12Mult->Draw("same");
    }

  } // end loop on mult
  fileout->Close();
  canvasPtSpectra->SaveAs(stringoutpdf + ".pdf");
  canvasPtSpectra->SaveAs(stringoutpdf + ".png");
  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}
