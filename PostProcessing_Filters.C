#include <Riostream.h>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TCutG.h>
#include "TFitResult.h"
#include "TLegend.h"

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

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX,
                TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
{
  if (XRange)
    histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->GetYaxis()->SetRangeUser(Low, Up);
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

Bool_t reject;
Double_t fparab(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[3] == 0)
  {
    LimInf = 0.474;
    LimSup = 0.520;
  }
  else if (par[3] == 4 || par[3] == 5 || par[3] == 8)
  {
    LimInf = 1.310;
    LimSup = 1.335;
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

Double_t fretta(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[2] == 0)
  {
    LimInf = 0.47;
    LimSup = 0.530;
  }
  else if (par[2] == 4 || par[2] == 5 || par[2] == 8)
  {
    LimInf = 1.310;
    LimSup = 1.335;
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0];
}

const Int_t numPart = 2;
TString TitleInvMass[numPart] = {"(#Lambda, #pi^{-}) & (#bar{#Lambda}, #pi^{+}) invariant mass (GeV/#it{c}^{2})", "(#Lambda, K^{-}) & (#bar{#Lambda}, K^{+}) invariant mass (GeV/#it{c}^{2})"};
TString namehisto[numPart] = {"", ""};
Float_t LowLimitMass[numPart] = {1.29, 1.62};
Float_t UpLimitMass[numPart] = {1.35, 1.72};
Float_t LowMassRange[numPart] = {1.29, 1.62};
Float_t UpMassRange[numPart] = {1.35, 1.72};

Float_t min_range_signal[numPart] = {1.31, 1.66}; // estremi region fit segnale (gaussiane)
Float_t max_range_signal[numPart] = {1.334, 1.685};
Float_t min_histo[numPart] = {1.30, 1.62}; // estremi del range degli istogrammi
Float_t max_histo[numPart] = {1.342, 1.72};
Float_t liminf[numPart] = {1.30, 1.66}; // estremi regione fit del bkg e total
Float_t limsup[numPart] = {1.342, 1.685};

Float_t lim_inf_mean[numPart] = {1.31, 1.66};
Float_t lim_sup_mean[numPart] = {1.33, 1.685};
Float_t lim_inf_sigma[numPart] = {0};
Float_t lim_sup_sigma[numPart] = {0.008, 0.008};
Float_t lim_inf_errmean[numPart] = {0};
Float_t lim_sup_errmean[numPart] = {10, 10}; // loooooose
Float_t lim_inf_errsigma[numPart] = {0};
Float_t lim_sup_errsigma[numPart] = {10, 10}; // loose

const Float_t massParticle[numPart] = {1.32171, 1.67245};
// TString Spart[numPart+2] = {"XiNeg", "XiPos", "OmegaNeg", "OmegaPlus"};
TString Spart[numPart] = {"Xi", "Omega"};

void PostProcessing_Filters(TString year = "23h_PIDpass3",
                            TString SPathIn = "../TriggerForRun3/EventFiltering2023/LHC23h/AnalysisResults_LHC23h_apass4_PIDapass3_CEFP.root",
                            TString OutputDir = "../TriggerForRun3/EventFiltering2023/LHC23h/",
                            Float_t ptthr = 7,
                            Bool_t isOldVersionBf2701 = 0)
{

  cout << "Input file: " << SPathIn << endl;
  TFile *filein = new TFile(SPathIn, "");
  if (!filein)
  {
    cout << "FileIn not available" << endl;
    return;
  }

  TDirectoryFile *dirCascB;
  dirCascB = (TFile *)filein->Get("cascade-builder");
  if (!dirCascB)
  {
    cout << "cascade-builder not available" << endl;
    return;
  }

  TDirectoryFile *dir;
  dir = (TFile *)filein->Get("lf-strangeness-filter");
  if (!dir)
  {
    cout << "dir not available" << endl;
    return;
  }

  TDirectoryFile *dirQA;
  dirQA = (TDirectoryFile *)dir->Get("QAHistos");
  if (!dirQA)
  {
    cout << "Directory QAHistos not available" << endl;
    return;
  }

  TH1F *hEvents = (TH1F *)dir->Get("hProcessedEvents");
  if (!hEvents)
  {
    cout << "hProcessedEvents not available " << endl;
    return;
  }

  // Nevents vs trigger type
  Double_t NEvents = hEvents->GetBinContent(1);
  cout << "NEvents " << NEvents << endl;
  hEvents->Scale(1. / NEvents);
  hEvents->GetYaxis()->SetRangeUser(10e-9, 1);
  hEvents->SetTitle("Rejection factors");
  hEvents->GetXaxis()->SetTitle("");
  hEvents->GetYaxis()->SetTitle("Rejection factor");

  TCanvas *canvasNEvents = new TCanvas("canvasNEvents", "canvasNEvents", 800, 500);
  gPad->SetLogy();
  hEvents->Draw("text");

  // Z vertex distribution of selected events
  TH1F *hVtxZ = (TH1F *)dirQA->Get("hVtxZ");
  if (!hVtxZ)
  {
    cout << "hVtxZ not available" << endl;
    return;
  }
  TCanvas *canvasVtxZ = new TCanvas("canvasVtxZ", "canvasVtxZ", 800, 500);
  hVtxZ->SetTitle("Z-vertex distribution of selected events");
  hVtxZ->Draw("");
  TLegend *legendNumberEvents = new TLegend(0.6, 0.8, 0.9, 0.9);
  legendNumberEvents->SetBorderSize(0);
  legendNumberEvents->AddEntry("", Form("Number of selected events: %.0f", NEvents));
  legendNumberEvents->Draw();

  // cascade candidates

  TH1F *hCascCandidatesB = (TH1F *)dirCascB->Get("hCascadeCriteria");
  if (!hCascCandidatesB)
  {
    cout << "hCascadeCriteria in cascade-builder not available " << endl;
    return;
  }

  Int_t NInitialCascB = hCascCandidatesB->GetBinContent(1);
  hCascCandidatesB->Scale(1. / NInitialCascB);
  hCascCandidatesB->GetXaxis()->SetBinLabel(1, "All candidates");
  hCascCandidatesB->GetXaxis()->SetBinLabel(2, "Lambda mass");
  hCascCandidatesB->GetXaxis()->SetBinLabel(3, "Bach TPC refit");
  hCascCandidatesB->GetXaxis()->SetBinLabel(4, "Bach TPC crossed rows");
  hCascCandidatesB->GetXaxis()->SetBinLabel(5, "Bach DCAxy");
  hCascCandidatesB->GetXaxis()->SetBinLabel(6, "Casc DCA dau");
  hCascCandidatesB->GetXaxis()->SetBinLabel(7, "Casc Cos PA");
  hCascCandidatesB->GetXaxis()->SetBinLabel(8, "Casc Radius");
  TCanvas *canvasCascB = new TCanvas("canvasCascB", "canvasCascB", 800, 500);
  hCascCandidatesB->Draw("text");

  TH1F *hCascCandidates = (TH1F *)dir->Get("hCandidate");
  if (!hCascCandidates)
  {
    cout << "hCandidate not available " << endl;
    return;
  }

  Int_t NInitialCasc = hCascCandidates->GetBinContent(1);
  hCascCandidates->Scale(1. / NInitialCasc);
  TCanvas *canvasCasc = new TCanvas("canvasCasc", "canvasCasc", 800, 500);
  hCascCandidates->Draw("text");

  // Track QA
  TDirectoryFile *dirQATrack;
  dirQATrack = (TDirectoryFile *)dir->Get("QAHistosTriggerParticles");
  if (!dirQATrack)
  {
    cout << "Directory QAHistosTriggerParticles not available" << endl;
    return;
  }

  TCanvas *canvasTrackQA = new TCanvas("canvasTrackQA", "canvasTrackQA", 800, 500);
  canvasTrackQA->Divide(3, 2);

  TCanvas *canvasTrackQA2D = new TCanvas("canvasTrackQA2D", "canvasTrackQA2D", 800, 500);
  canvasTrackQA2D->Divide(2, 2);

  const Int_t numQATrackHistos = 5;
  TH1F *hTrackQAPt = (TH1F *)dirQATrack->Get("hPtTriggerAllEv");
  canvasTrackQA->cd(1);
  hTrackQAPt->Draw("");

  TH2F *hTrackQA[numQATrackHistos - 1];
  TH1F *hTrackQA1D[numQATrackHistos - 1];
  TString hSTrackQA[numQATrackHistos - 1] = {"hEtaTriggerAllEv", "hPhiTriggerAllEv", "hDCAxyTriggerAllEv", "hDCAzTriggerAllEv"};

  const Int_t numPtTrigg = 6;
  // Float_t binptTrigg[numPtTrigg + 1] = {0, 0.5, 1, 1.5, 2, 3, 10};
  // Float_t binptTrigg[numPtTrigg + 1] = {5, 5.5, 6.0, 6.5, 7, 8, 10};
  Float_t binptTrigg[numPtTrigg + 1] = {7, 7.5, 8, 8.5, 9, 9.5, 10};
  Int_t ColorPtTrigg[numPtTrigg + 1] = {634, 628, 797, 815, 418, 429, 867};
  Int_t MarkerPtTrigg[numPtTrigg + 1] = {20, 21, 22, 33, 20, 21, 22};

  TLegend *legendPtTrigg = new TLegend(0.1, 0.1, 0.9, 0.9);
  legendPtTrigg->SetBorderSize(0);

  for (Int_t i = 0; i < numQATrackHistos - 1; i++)
  {
    hTrackQA[i] = (TH2F *)dirQATrack->Get(hSTrackQA[i]);
    if (!hTrackQA[i])
    {
      cout << hSTrackQA[i] << " not available" << endl;
      return;
    }

    canvasTrackQA2D->cd(i + 1);
    hTrackQA[i]->Draw("colz");

    canvasTrackQA->cd(i + 2);
    for (Int_t pt = 0; pt < numPtTrigg; pt++)
    {
      // if (binptTrigg[pt] > 9)
      // continue;

      hTrackQA1D[i] = (TH1F *)hTrackQA[i]->ProjectionX(Form("hTrackQA_var%i_pt%i", i, pt), hTrackQA[i]->GetYaxis()->FindBin(binptTrigg[pt] + 0.001), hTrackQA[i]->GetYaxis()->FindBin(binptTrigg[pt + 1] - 0.001));
      hTrackQA1D[i]->SetLineColor(ColorPtTrigg[pt]);
      hTrackQA1D[i]->SetMarkerColor(ColorPtTrigg[pt]);
      hTrackQA1D[i]->SetMarkerSize(MarkerPtTrigg[pt]);
      hTrackQA1D[i]->SetMarkerSize(1.5);
      if (hSTrackQA[i] == "hDCAxyTriggerAllEv" || hSTrackQA[i] == "hDCAzTriggerAllEv")
        hTrackQA1D[i]->Scale(1. / hTrackQA1D[i]->Integral(hTrackQA1D[i]->FindBin(-0.1), hTrackQA1D[i]->FindBin(0.1)));
      else
        hTrackQA1D[i]->Scale(1. / hTrackQA1D[i]->GetEntries());

      // if (hSTrackQA[i] == "hPhiTriggerAllEv") hTrackQA1D[i]->Rebin(2);
      // hTrackQA1D[i]->Rebin(4);
      if (hSTrackQA[i] == "hDCAxyTriggerAllEv")
        hTrackQA1D[i]->GetXaxis()->SetRangeUser(-0.05, 0.05);
      hTrackQA1D[i]->GetYaxis()->SetRangeUser(0, 2 * hTrackQA1D[i]->GetMaximum());

      if (i == 0)
        legendPtTrigg->AddEntry(hTrackQA1D[i], Form("%.1f < p_{T} < %.1f GeV/c", binptTrigg[pt], binptTrigg[pt + 1]), "pl");
      hTrackQA1D[i]->Draw("same e l");
    }
  }
  canvasTrackQA->cd(6);
  legendPtTrigg->Draw("");

  TH1F *hTriggerParticles = (TH1F *)dirQATrack->Get("hTriggeredParticlesAllEv");
  if (!hTriggerParticles)
  {
    return;
  }
  hTriggerParticles->GetXaxis()->SetRangeUser(0, 8);
  hTriggerParticles->SetTitle("Number of trigger particles (p_{T} > p_{T}^{th}) per event");
  TCanvas *canvasTrigger = new TCanvas("canvasTrigger", "canvasTrigger", 800, 500);
  gPad->SetLogy();
  hTriggerParticles->Draw("");

  TH1F *hRejFactorshXi = (TH1F *)dir->Get("hEvtvshMinPt");
  if (!hRejFactorshXi)
  {
    return;
  }
  TCanvas *canvasRejFactorshXi = new TCanvas("canvasRejFactorshXi", "canvasRejFactorshXi", 800, 500);
  hRejFactorshXi->Scale(1. / NEvents);
  hRejFactorshXi->GetXaxis()->SetRangeUser(ptthr, 11);
  hRejFactorshXi->SetTitle("Rejection factors of hXi events for p_{T} > p_{T}^{th}");
  hRejFactorshXi->GetXaxis()->SetTitle("p_{T}^{th} (GeV/c)");
  hRejFactorshXi->GetYaxis()->SetTitle("Rejection factor");
  hRejFactorshXi->Draw("text");

  // TPC nsigma distributions for all cascade candidates after selections
  TH2F *hNSigmaTPCvsPt[12];
  TString SNSigmaTPCvsPt[12] = {"hTPCNsigmaXiBachPiPlus", "hTPCNsigmaXiV0PiPlus", "hTPCNsigmaXiV0AntiProton",
                                "hTPCNsigmaXiBachPiMinus", "hTPCNsigmaXiV0PiMinus", "hTPCNsigmaXiV0Proton",
                                "hTPCNsigmaOmegaBachKaPlus", "hTPCNsigmaOmegaV0PiPlus", "hTPCNsigmaOmegaV0AntiProton",
                                "hTPCNsigmaOmegaBachKaMinus", "hTPCNsigmaOmegaV0PiMinus", "hTPCNsigmaOmegaV0Proton"};
  if (isOldVersionBf2701)
  {
    SNSigmaTPCvsPt[0] = "hTPCNsigmaPi";
    SNSigmaTPCvsPt[1] = "hTPCNsigmaPr";
    SNSigmaTPCvsPt[2] = "hTPCNsigmaBachPi";
    SNSigmaTPCvsPt[3] = "hTPCNsigmaBachKa";
  }
  for (Int_t dau = 0; dau < 12; dau++)
  {
    if (isOldVersionBf2701 && dau > 3)
      continue;
    hNSigmaTPCvsPt[dau] = (TH2F *)dirQA->Get(SNSigmaTPCvsPt[dau]);
    if (!hNSigmaTPCvsPt[dau])
    {
      cout << SNSigmaTPCvsPt[dau] << " not available" << endl;
      return;
    }
  }
  const Int_t numPtDau = 6; // six pt intervals

  Float_t binptDauPi[numPtDau + 1] = {0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.2};
  Float_t binptDauPr[numPtDau + 1] = {0.1, 0.3, 0.5, 0.7, 0.8, 1.2, 1.6};
  Float_t binptDauBachPi[numPtDau + 1] = {0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.2};
  Float_t binptDauBachKa[numPtDau + 1] = {0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.2};
  Float_t binptDau[numPtDau + 1] = {0, 0, 0, 0, 0, 0, 0};

  TString SPtDau[numPtDau] = {""};
  TH1F *hNSigmaTPC[12][numPtDau];
  TCanvas *canvasTPC[12];
  TH1F *hNSigmaSummary[12];
  TH1F *hNSigmaSummaryGaus[12];
  TCanvas *cNSigmaSummary = new TCanvas("cNSigmaSummary", "cNSigmaSummary", 800, 500);
  cNSigmaSummary->Divide(6, 2);

  for (Int_t dau = 0; dau < 12; dau++)
  {
    if (isOldVersionBf2701 && dau > 3)
      continue;

    canvasTPC[dau] = new TCanvas(Form("canvasTPC_dau%i", dau), Form("canvasTPC_dau%i", dau), 1800, 1400);
    canvasTPC[dau]->Divide(numPtDau / 2, 2);
    StyleCanvas(canvasTPC[dau], 0.15, 0.05, 0.05, 0.15);

    for (Int_t pt = 0; pt < numPtDau + 1; pt++)
    {
      if (dau == 0 || dau == 3)
        binptDau[pt] = binptDauBachPi[pt];
      else if (dau == 1 || dau == 4 || dau == 7 || dau == 10)
        binptDau[pt] = binptDauPi[pt];
      else if (dau == 2 || dau == 5 || dau == 8 || dau == 11)
        binptDau[pt] = binptDauPr[pt];
      else if (dau == 6 || dau == 9)
        binptDau[pt] = binptDauBachKa[pt];

      if (isOldVersionBf2701)
      {
        if (dau == 0)
          binptDau[pt] = binptDauPi[pt];
        else if (dau == 1)
          binptDau[pt] = binptDauPr[pt];
        else if (dau == 2)
          binptDau[pt] = binptDauBachPi[pt];
        else if (dau == 6)
          binptDau[pt] = binptDauBachKa[pt];
      }
    }

    hNSigmaSummary[dau] = new TH1F("hNSigmaSummary_" + SNSigmaTPCvsPt[dau], "hNSigmaSummary_" + SNSigmaTPCvsPt[dau], numPtDau, binptDau);
    hNSigmaSummaryGaus[dau] = new TH1F("hNSigmaSummaryGaus_" + SNSigmaTPCvsPt[dau], "hNSigmaSummaryGaus_" + SNSigmaTPCvsPt[dau], numPtDau, binptDau);

    for (Int_t pt = 0; pt < numPtDau; pt++)
    {
      SPtDau[pt] = Form("%.1f < p_{T} < %.1f", binptDau[pt], binptDau[pt + 1]);

      hNSigmaTPC[dau][pt] = (TH1F *)hNSigmaTPCvsPt[dau]->ProjectionX(Form("hNSigmaTPC_dau%i_pt%i", dau, pt), hNSigmaTPCvsPt[dau]->GetYaxis()->FindBin(binptDau[pt] + 0.001), hNSigmaTPCvsPt[dau]->GetYaxis()->FindBin(binptDau[pt + 1] - 0.001));
      hNSigmaTPC[dau][pt]->Rebin(2);
      StyleHisto(hNSigmaTPC[dau][pt], 0, 1.2 * hNSigmaTPC[dau][pt]->GetBinContent(hNSigmaTPC[dau][pt]->GetMaximumBin()), 1, 20, "Nsigma", "Counts", SNSigmaTPCvsPt[dau] + " " + SPtDau[pt], 1, -6, 6, 1.4, 1.4, 1.2);
      hNSigmaTPC[dau][pt]->GetXaxis()->SetRangeUser(-6, 6);
      hNSigmaTPC[dau][pt]->GetYaxis()->SetRangeUser(0, 1.2 * hNSigmaTPC[dau][pt]->GetBinContent(hNSigmaTPC[dau][pt]->GetMaximumBin()));
      TF1* gaussTPC = new TF1("gaussTPC", "gaus", -2.5, 2.5);
      gaussTPC->SetLineColor(2);
      hNSigmaTPC[dau][pt]->Fit("gaussTPC", "R+");
      canvasTPC[dau]->cd(pt + 1);
      gPad->SetTopMargin(0.08);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.1);
      gPad->SetBottomMargin(0.2);
      hNSigmaTPC[dau][pt]->DrawClone("e same");
      hNSigmaSummary[dau]->SetBinContent(pt + 1, hNSigmaTPC[dau][pt]->GetMean());
      hNSigmaSummary[dau]->SetBinError(pt + 1, hNSigmaTPC[dau][pt]->GetRMS());
      hNSigmaSummaryGaus[dau]->SetBinContent(pt + 1, gaussTPC->GetParameter(1));
      hNSigmaSummaryGaus[dau]->SetBinError(pt + 1, gaussTPC->GetParameter(2));
    }
    // summary plot
    cNSigmaSummary->cd(dau + 1);
    hNSigmaSummary[dau]->SetTitle(SNSigmaTPCvsPt[dau]);
    hNSigmaSummary[dau]->GetYaxis()->SetRangeUser(-4, 4);
    hNSigmaSummary[dau]->Draw("same");
  }

  // Cascade topological variables (after all topo sleections, for the time being)
  TDirectoryFile *dirCascVar;
  dirCascVar = (TDirectoryFile *)dir->Get("QAHistosTopologicalVariables");
  if (!dirCascVar)
  {
    cout << "Directory QAHistosTopologicalVariables not available" << endl;
    return;
  }

  const int NTopCascVar = 11;

  TH1F *hTopVar[NTopCascVar];
  TCanvas *canvasTopology[2][2];
  canvasTopology[0][0] = new TCanvas("canvasTopologyXi1", "canvasTopologyXi1", 1800, 1400);
  canvasTopology[0][0]->Divide(3, 2);
  canvasTopology[0][1] = new TCanvas("canvasTopologyXi2", "canvasTopologyXi2", 1800, 1400);
  canvasTopology[0][1]->Divide(3, 2);
  canvasTopology[1][0] = new TCanvas("canvasTopologyOmega1", "canvasTopologyOmega1", 1800, 1400);
  canvasTopology[1][0]->Divide(3, 2);
  canvasTopology[1][1] = new TCanvas("canvasTopologyOmega2", "canvasTopologyOmega2", 1800, 1400);
  canvasTopology[1][1]->Divide(3, 2);

  TString TopVarCascInput[NTopCascVar] = {"CascCosPA", "V0CosPA", "CascRadius", "V0Radius",
                                          "InvMassLambda", "DCACascDaughters", "DCAV0Daughters", "DCABachToPV",
                                          "DCAV0ToPV", "DCAPosToPV", "DCANegToPV"};
  TString TopVarCasc[NTopCascVar] = {"Casc #it{cos}#theta_{PA}", "V0 #it{cos}#theta_{PA}", "Casc #it{R}", "V0 #it{R}",
                                     "#it{m}_{inv} #Lambda Daughter", "DCA Casc Daughters", "DCA V0 Daughters", "DCA Bach. To PV",
                                     "DCA V0 To PV", "DCA Pos. To PV", "DCA Neg. To PV"};
  TString TopVarCascUnit[NTopCascVar] = {"", "", "(cm)", "(cm)",
                                         "(GeV/#it{c}^2)", "(cm)", "(cm)", "(cm)",
                                         "(cm)", "(cm)", "(cm)"};
  Double_t TopVarCascCuts[NTopCascVar] = {0.95, 0.95, 0.3, 1,
                                          0, 2, 2, 0.05,
                                          0, 0.05, 0.05};
  TLine *lineCascCuts[NTopCascVar];
  TString SXiorOmega[2] = {"Xi", "Omega"};
  TString prefix = "h";

  for (Int_t XiorOmega = 0; XiorOmega < 2; XiorOmega++)
  {
    if (isOldVersionBf2701)
    {
      SXiorOmega[0] = "";
      prefix = "";
      if (XiorOmega == 1)
        continue;
    }
    for (Int_t var = 0; var < NTopCascVar; var++)
    {
      hTopVar[var] = (TH1F *)dirCascVar->Get(prefix + TopVarCascInput[var] + SXiorOmega[XiorOmega]);
      if (!hTopVar[var])
      {
        cout << prefix << TopVarCascInput[var] << SXiorOmega[XiorOmega] << " not available" << endl;
        return;
      }
      hTopVar[var]->Scale(1. / NEvents);
      hTopVar[var]->GetYaxis()->SetRangeUser(0.1 * hTopVar[var]->GetMinimum(1.e-10), 10 * hTopVar[var]->GetMaximum());
      hTopVar[var]->SetTitle(TopVarCasc[var]);
      hTopVar[var]->GetXaxis()->SetTitle(TopVarCasc[var] + " " + TopVarCascUnit[var]);
      hTopVar[var]->GetYaxis()->SetTitle("1/N_{ev} Counts");
      if (var < 6)
        canvasTopology[XiorOmega][0]->cd(var + 1);
      else
        canvasTopology[XiorOmega][1]->cd(var + 1 - 6);
      gPad->SetLogy();
      hTopVar[var]->DrawCopy("hist");
      lineCascCuts[var] = new TLine(TopVarCascCuts[var], hTopVar[var]->GetMinimum(), TopVarCascCuts[var], hTopVar[var]->GetMaximum());
      // lineCascCuts[var] = new TLine(TopVarCascCuts[var], 0, TopVarCascCuts[var], 100);
      lineCascCuts[var]->SetLineColor(867);
      lineCascCuts[var]->Draw("");
    }
  }
  // repidity, eta, ctau and pt
  const int NKineCascVar = 8;

  TH1F *hKineVar[NKineCascVar];
  TCanvas *canvasKine[2];
  canvasKine[0] = new TCanvas("canvasKine1", "canvasKine1", 1800, 1400);
  canvasKine[0]->Divide(2, 2);
  canvasKine[1] = new TCanvas("canvasKine2", "canvasKine2", 1800, 1400);
  canvasKine[1]->Divide(2, 2);

  TString KineVarCascInput[NKineCascVar] = {"hPtXi", "hEtaXi", "hRapXi", "hProperLifetimeXi", "hPtOmega", "hEtaOmega", "hRapOmega", "hProperLifetimeOmega"};
  TString KineVarCasc[NKineCascVar] = {"p_{T} Xi", "#eta Xi", "y Xi", "ctau Xi", "p_{T} Omega", "#eta Omega", "y Omega", "ctau Omega"};
  TString KineVarCascUnit[NKineCascVar] = {"GeV/c", "", "", "(cm)", "GeV/c", "", "", "(cm)"};

  for (Int_t var = 0; var < NKineCascVar; var++)
  {
    if (isOldVersionBf2701 && var > 5)
      continue;
    if (isOldVersionBf2701 || (KineVarCascInput[var] != "hProperLifetimeXi" && KineVarCascInput[var] != "hProperLifetimeOmega"))
      hKineVar[var] = (TH1F *)dirQA->Get(KineVarCascInput[var]);
    else
      hKineVar[var] = (TH1F *)dirCascVar->Get(KineVarCascInput[var]);
    if (!hKineVar[var])
    {
      cout << KineVarCascInput[var] << " not available" << endl;
      return;
    }
    hKineVar[var]->Scale(1. / NEvents);
    hKineVar[var]->GetYaxis()->SetRangeUser(0.1 * hKineVar[var]->GetMinimum(1.e-10), 1.2 * hKineVar[var]->GetMaximum());
    hKineVar[var]->SetTitle(KineVarCasc[var]);
    hKineVar[var]->GetXaxis()->SetTitle(KineVarCasc[var] + " " + KineVarCascUnit[var]);
    hKineVar[var]->GetYaxis()->SetTitle("1/N_{ev} Counts");
    if (var < 4)
      canvasKine[0]->cd(var + 1);
    else
      canvasKine[1]->cd(var + 1 - 4);
    // gPad->SetLogy();
    hKineVar[var]->DrawCopy("hist");
  }

  // inv mass distributions of Xi (all selected xis) and Omegas (all selected omegas) vs pT
  TH2F *hMassvsPt[2];
  hMassvsPt[0] = (TH2F *)dirQA->Get("hMassXiAfterSelvsPt");
  hMassvsPt[1] = (TH2F *)dirQA->Get("hMassOmegaAfterSelvsPt");
  if (!hMassvsPt[0])
  {
    cout << "hMassXiAfterSelvsPt not available" << endl;
    return;
  }
  if (!hMassvsPt[1])
  {
    cout << "hMassOmegaAfterSelvsPt not available" << endl;
    return;
  }
  const Int_t numPt = 6; // six pt intervals
  // Xi
  Float_t binpt[numPt + 1] = {1.0, 1.5, 1.8, 2.1, 2.4, 2.6, 6.0};

  // Omega
  // Float_t binpt[numPt + 1] = {1.0, 1.5, 2.0, 3.0, 6.0};

  TString SPt[numPt] = {""};
  TH1F *hInvMass[numPart][numPt];

  TCanvas *canvas[numPart];

  TH1F *histoCountsPerEvent[numPart];
  TH1F *histoYield[numPart];

  Float_t counts = 0;
  Float_t errcount = 0;

  for (Int_t part = 0; part < 2; part++)
  {
    canvas[part] = new TCanvas(Form("canvas_part%i", part), Form("canvas_part%i", part), 1800, 1400);
    canvas[part]->Divide(numPt / 2, 2);
    StyleCanvas(canvas[part], 0.15, 0.05, 0.05, 0.15);

    histoCountsPerEvent[part] = new TH1F(Form("histoCountsPerEvent_part%i", part), Form("histoCountsPerEvent_part%i", part), numPt, binpt);
    histoYield[part] = new TH1F(Form("histoYield_part%i", part), Form("histoYield_part%i", part), numPt, binpt);

    for (Int_t pt = 0; pt < numPt; pt++)
    {
      SPt[pt] = Form("%.1f < p_{T} < %.1f", binpt[pt], binpt[pt + 1]);

      hInvMass[part][pt] = (TH1F *)hMassvsPt[part]->ProjectionX(Form("hInvMass_part%i_pt%i", part, pt), hMassvsPt[part]->GetYaxis()->FindBin(binpt[pt] + 0.001), hMassvsPt[part]->GetYaxis()->FindBin(binpt[pt + 1] - 0.001));
      hInvMass[part][pt]->Rebin(2);
      StyleHisto(hInvMass[part][pt], 0, 1.2 * hInvMass[part][pt]->GetBinContent(hInvMass[part][pt]->GetMaximumBin()), 1, 20, TitleInvMass[part], "Counts", SPt[pt], 1, LowLimitMass[part], UpLimitMass[part], 1.4, 1.4, 1.2);
      canvas[part]->cd(pt + 1);
      gPad->SetTopMargin(0.08);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.1);
      gPad->SetBottomMargin(0.2);
      hInvMass[part][pt]->Draw("e same");

      counts = 0;
      errcount = 0;
      for (Int_t bmass = hInvMass[part][pt]->GetXaxis()->FindBin(LowMassRange[part]); bmass <= hInvMass[part][pt]->GetXaxis()->FindBin(UpMassRange[part]); bmass++)
      {
        counts += hInvMass[part][pt]->GetBinContent(bmass);
      }
      errcount = sqrt(counts);
      histoCountsPerEvent[part]->SetBinContent(pt + 1, counts / NEvents / histoCountsPerEvent[part]->GetBinWidth(pt + 1));
      histoCountsPerEvent[part]->SetBinError(pt + 1, errcount / NEvents / histoCountsPerEvent[part]->GetBinWidth(pt + 1));
    }
  }

  // h-Xi events vs pt,thr

  TString Soutputfile = OutputDir + "FilterPostProcessing_" + year;
  canvasNEvents->SaveAs(Soutputfile + ".pdf(");
  canvasRejFactorshXi->SaveAs(Soutputfile + ".pdf");
  canvasVtxZ->SaveAs(Soutputfile + ".pdf");
  canvasTrackQA2D->SaveAs(Soutputfile + ".pdf");
  canvasTrackQA->SaveAs(Soutputfile + ".pdf");
  canvasTrigger->SaveAs(Soutputfile + ".pdf");
  for (Int_t dau = 0; dau < 12; dau++)
  {
    if (isOldVersionBf2701 && dau > 3)
      continue;
    canvasTPC[dau]->SaveAs(Soutputfile + ".pdf");
  }
  canvasCasc->SaveAs(Soutputfile + ".pdf");
  canvasCascB->SaveAs(Soutputfile + ".pdf");
  canvasTopology[0][0]->SaveAs(Soutputfile + ".pdf");
  canvasTopology[0][1]->SaveAs(Soutputfile + ".pdf");
  if (!isOldVersionBf2701)
  {
    canvasTopology[1][0]->SaveAs(Soutputfile + ".pdf");
    canvasTopology[1][1]->SaveAs(Soutputfile + ".pdf");
  }
  canvasKine[0]->SaveAs(Soutputfile + ".pdf");
  canvasKine[1]->SaveAs(Soutputfile + ".pdf");
  canvas[0]->SaveAs(Soutputfile + ".pdf");
  canvas[1]->SaveAs(Soutputfile + ".pdf)");

  TFile *outputfile = new TFile(Soutputfile + ".root", "RECREATE");
  outputfile->WriteTObject(hEvents);
  outputfile->WriteTObject(hCascCandidates);
  outputfile->WriteTObject(hTriggerParticles);
  outputfile->WriteTObject(canvasTrackQA);

  for (Int_t dau = 0; dau < 12; dau++)
  {
    outputfile->WriteTObject(hNSigmaSummary[dau]);
    outputfile->WriteTObject(hNSigmaSummaryGaus[dau]);
    for (Int_t pt = 0; pt < numPtDau; pt++)
    {
      outputfile->WriteTObject(hNSigmaTPC[dau][pt]);
    }
  }

  outputfile->Close();
  cout << "Ho creato il file: " << Soutputfile << " (.pdf and .root)" << endl;

  cout << "\nTotal number of processed events " << NEvents << endl;
}
