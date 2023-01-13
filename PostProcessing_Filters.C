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

void PostProcessing_Filters(TString year = "LHC22m_pass2_NoGlobalTracks",
                            TString SPathIn = /*"../TriggerForRun3/AnalysisResults_FinalTOT_NoTOF.root"*/
                            "../TriggerForRun3/AnalysisResults_22mpass2_New_8_NoGlobalTrackHC_Eta0.9.root",
                            Int_t part = 0,
                            Bool_t UseTwoGauss = 1,
                            Bool_t isBkgParab = 0,
                            Bool_t isMeanFixedPDG = 0,
                            Float_t sigmacentral = 4)
{

  TFile *filein = new TFile(SPathIn, "");
  if (!filein)
  {
    cout << "FileIn not available" << endl;
    return;
  }
  TDirectoryFile *dir;
  dir = (TFile *)filein->Get("lf-strangeness-filter");
  if (!dir)
  {
    cout << "dir not available" << endl;
    return;
  }

  TH1F *hEvents = (TH1F *)dir->Get("hProcessedEvents");
  if (!hEvents)
  {
    cout << "hProcessedEvents not available " << endl;
    return;
  }

  Int_t NEvents = hEvents->GetBinContent(1);
  hEvents->Scale(1. / NEvents);
  hEvents->GetYaxis()->SetRangeUser(10e-8, 1);

  TCanvas *canvasNEvents = new TCanvas("canvasNEvents", "canvasNEvents", 800, 500);
  gPad->SetLogy();
  hEvents->Draw("");

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

  const Int_t numPtTrigg = 7;
  Float_t binptTrigg[numPtTrigg + 1] = {0, 0.5, 1, 1.5, 2, 3, 10};
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
      if (binptTrigg[pt] > 2.5)
        continue;

      hTrackQA1D[i] = (TH1F *)hTrackQA[i]->ProjectionX(Form("hTrackQA_var%i_pt%i", i, pt), hTrackQA[i]->GetYaxis()->FindBin(binptTrigg[pt] + 0.001), hTrackQA[i]->GetYaxis()->FindBin(binptTrigg[pt + 1] - 0.001));
      hTrackQA1D[i]->SetLineColor(ColorPtTrigg[pt]);
      hTrackQA1D[i]->SetMarkerColor(ColorPtTrigg[pt]);
      hTrackQA1D[i]->SetMarkerSize(MarkerPtTrigg[pt]);
      hTrackQA1D[i]->SetMarkerSize(1.5);
      if ( hSTrackQA[i] == "hDCAxyTriggerAllEv" || hSTrackQA[i] == "hDCAzTriggerAllEv") hTrackQA1D[i]->Scale(1. / hTrackQA1D[i]->Integral(hTrackQA1D[i]->FindBin(-0.1),hTrackQA1D[i]->FindBin(0.1)));
      else  hTrackQA1D[i]->Scale(1. / hTrackQA1D[i]->GetEntries());

      // if (hSTrackQA[i] == "hPhiTriggerAllEv") hTrackQA1D[i]->Rebin(2);
      hTrackQA1D[i]->Rebin(4);
      if (hSTrackQA[i] == "hDCAxyTriggerAllEv")
        hTrackQA1D[i]->GetXaxis()->SetRangeUser(-0.1, 0.1);
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
  TCanvas *canvasTrigger = new TCanvas("canvasTrigger", "canvasTrigger", 800, 500);
  hTriggerParticles->Draw("");

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
  TCanvas *canvasTopology[2];
  canvasTopology[0] = new TCanvas("canvasTopology1", "canvasTopology1", 1800, 1400);
  canvasTopology[0]->Divide(3, 2);
  canvasTopology[1] = new TCanvas("canvasTopology2", "canvasTopology2", 1800, 1400);
  canvasTopology[1]->Divide(3, 2);

  TString TopVarCascInput[NTopCascVar] = {"CascCosPA", "V0CosPA", "CascRadius", "V0Radius",
                                          "InvMassLambda", "DCACascDaughters", "DCAV0Daughters", "DCABachToPV",
                                          "DCAV0ToPV", "DCAPosToPV", "DCANegToPV"};
  // TString TopVarCascInput[NTopCascVar] = {"hXiCascCosPA", "hXiV0CosPA", "hXiCascRadius", "hXiV0Radius",
  //                                         "hXiInvMassLambda", "hXiDCACascDaughters", "hXiDCAV0Daughters", "hXiDCABachToPV",
  //                                         "hXiDCAV0ToPV", "hXiDCAPosToPV", "hXiDCANegToPV"};
  TString TopVarCasc[NTopCascVar] = {"Casc #it{cos}#theta_{PA}", "V0 #it{cos}#theta_{PA}", "Casc #it{R}", "V0 #it{R}",
                                     "#it{m}_{inv} #Lambda Daughter", "DCA Casc Daughters", "DCA V0 Daughters", "DCA Bach. To PV",
                                     "DCA V0 To PV", "DCA Pos. To PV", "DCA Neg. To PV"};
  TString TopVarCascUnit[NTopCascVar] = {"", "", "(cm)", "(cm)",
                                         "(GeV/#it{c}^2)", "(cm)", "(cm)", "(cm)",
                                         "(cm)", "(cm)", "(cm)"};

  for (Int_t var = 0; var < NTopCascVar; var++)
  {
    hTopVar[var] = (TH1F *)dirCascVar->Get(TopVarCascInput[var]);
    if (!hTopVar[var])
    {
      cout << TopVarCascInput[var] << " not available" << endl;
      return;
    }
    hTopVar[var]->Scale(1. / NEvents);
    hTopVar[var]->GetYaxis()->SetRangeUser(0.1 * hTopVar[var]->GetMinimum(1.e-10), 10 * hTopVar[var]->GetMaximum());
    hTopVar[var]->SetTitle(TopVarCasc[var]);
    hTopVar[var]->GetXaxis()->SetTitle(TopVarCasc[var] + " " + TopVarCascUnit[var]);
    hTopVar[var]->GetYaxis()->SetTitle("1/N_{ev} Counts");
    if (var < 6)
      canvasTopology[0]->cd(var + 1);
    else
      canvasTopology[1]->cd(var + 1 - 6);
    gPad->SetLogy();
    hTopVar[var]->DrawCopy("hist");
  }

  // inv mass distributions of Xi (all selected xis) and Omegas (all selected omegas) vs pT
  TDirectoryFile *dirCasc;
  dirCasc = (TDirectoryFile *)dir->Get("QAHistos");
  if (!dirCasc)
  {
    cout << "Directory QAHistos not available" << endl;
    return;
  }

  TH2F *hMassvsPt[2];
  hMassvsPt[0] = (TH2F *)dirCasc->Get("hMassXiAfterSelvsPt");
  hMassvsPt[1] = (TH2F *)dirCasc->Get("hMassOmegaAfterSelvsPt");
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

  TString Soutputfile = "../TriggerForRun3/" + year;
  canvasNEvents->SaveAs(Soutputfile + ".pdf(");
  canvasTrackQA2D->SaveAs(Soutputfile + ".pdf");
  canvasTrackQA->SaveAs(Soutputfile + ".pdf");
  canvasTrigger->SaveAs(Soutputfile + ".pdf");
  canvasCasc->SaveAs(Soutputfile + ".pdf");
  canvasTopology[0]->SaveAs(Soutputfile + ".pdf");
  canvasTopology[1]->SaveAs(Soutputfile + ".pdf");
  canvas[0]->SaveAs(Soutputfile + ".pdf");
  canvas[1]->SaveAs(Soutputfile + ".pdf)");

  TFile *outputfile = new TFile(Soutputfile + ".root", "RECREATE");
  outputfile->WriteTObject(hEvents);
  outputfile->WriteTObject(hCascCandidates);
  outputfile->WriteTObject(hTriggerParticles);
  outputfile->WriteTObject(canvasTrackQA);
  outputfile->Close();
  cout << "Ho creato il file: " << Soutputfile << " (.pdf and .root)" << endl;

  cout << "\nTotal number of processed events " << NEvents << endl;
}
