// Adaptation of Post Processing macro used for Strangeness QA of Run2 AODs
// ========================================================================
//
//    This macro was originally written for Run2 Strangeness QA by:
//    alessandro.balbino@cern.ch
//    chiara.de.martin@cern.ch

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

Bool_t reject = kTRUE;
Double_t fparab(Double_t *x, Double_t *par)
{
  const Int_t numPart = 7;
  Float_t liminf[numPart] = {0.48, 1.11, 1.11, 1.31, 1.31, 1.665, 1.665};
  Float_t limsup[numPart] = {0.54, 1.13, 1.13, 1.335, 1.335, 1.685, 1.685};
  Int_t part = par[3];
  if (reject && x[0] > liminf[part] && x[0] < limsup[part])
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

void checkExactLimit(TH1F *lHist, Float_t limit, bool isMinLimit, TLatex labels[]);
void checkIntervalWidth(TH1F *lHist, float width, TLatex labels[]);
void checkBoundaries(TH1F *lHist, Float_t lower_bound, Float_t upper_bound, TLatex labels[]);
void checkBoundaries(Float_t value1, Float_t lower_bound1, Float_t upper_bound1, Float_t value2, Float_t lower_bound2, Float_t upper_bound2, TLatex labels[]);
void setPadOptions(bool useLogY);
void setHistGraphics(TH1F *lHist, bool isOldPass);
float findMaximum(TH1F *lHist, float width);
float findMaxValue(TH1F *lHist1, TH1F *lHist2);
Double_t SetEfficiencyError(Int_t k, Int_t n);

void PostProcessV0AndCascQA_AO2D(TString CollType = "pp", Bool_t isMC = true, Int_t RebinTPC = 1,
                                 Int_t SkipCascFits = 0,                                                                                                                                                                                                                                                                                                      // 0 = don't skip, 1 = skip partc, 2 = skip all cascades
                                 Bool_t TopologyOnly = false,                                                                                                                                                                                                                                                                                                 // true = only topology analysis, false = complete analysis
                                 TString PathIn = "/Users/mbp-cdm-01/Desktop/AssegnoRicerca/Run3Analyses/OmegavsMult/EnrichedMonteCarlos/AnalysisResults_TestGapTriggeredMCBis_NoTPC.root"/*"../Run3QA/Periods/LHC23e1b/AnalysisResults_qatask_LHC23e1b_Train103380.root" /*"../Run3QA/Periods/LHC23e1_Test/AnalysisResults_qatask_LHC23e1_Test2.root" /*"../Run3QA/Periods/LHC22q_pass3/AnalysisResults_qatask_LHC22q_pass3_Train64263.root"*/, // input file name
                                 TString PathOut = "/Users/mbp-cdm-01/Desktop/AssegnoRicerca/Run3Analyses/OmegavsMult/EnrichedMonteCarlos/PostProcessing_TestGapTriggeredMCBis_NoTPC"/* "../Run3QA/Periods/LHC23e1b/PostProcess_qa_LHC23e1b_Train103380"/*"../Run3QA/Periods/LHC23e1_Test/PostProcess_qa_LHC23e1_Test2"*/,                                                                                                          // output file name
                                 Bool_t CheckOldPass = false,                                                                                                                                                                                                                                                                                                 // true to compare two passes
                                 TString OldPassPath = "..",                                                                                                                                                                                                                                                                                                  // input/output file name (old pass to be compared with)
                                 Bool_t isMassvsRadiusPlots = 0)
{
  // Define pass names
  TString pass_names[2] = {"pass1", "pass3"};
  // Define text
  TLatex cutCheckLabels[8] = {TLatex(0.35, 0.8, "#scale[0.8]{#color[3]{BOUNDARIES: OK!!}}"), TLatex(0.2, 0.7, "#scale[0.8]{#color[2]{TIGHTER CUT WRT EXPECTED!!}}"), TLatex(0.2, 0.7, "#scale[0.8]{#color[42]{LOOSER CUT WRT EXPECTED!!}}"), TLatex(0.2, 0.63, "#scale[0.8]{#color[2]{PROBLEM FOR ANALYSIS}}"), TLatex(0.2, 0.63, "#scale[0.8]{#color[42]{NOT AN ISSUE FOR ANALYSIS}}"), TLatex(0.3, 0.8, "#scale[0.8]{#color[42]{OUTSIDE BOUNDARIES!!}}"), TLatex(0.3, 0.75, "#scale[0.8]{#color[42]{PLEASE CHECK.}}"), TLatex(0.35, 0.8, "#scale[0.8]{#color[3]{OK!!}}")};
  for (int i = 0; i < 8; i++)
  {
    cutCheckLabels[i].SetNDC();
  }
  TFile *f = new TFile(PathIn, "");
  if (!f)
  {
    printf("\nFILE NOT FOUND!\n");
    return;
  }
  TDirectoryFile *dir = (TDirectoryFile *)f->Get("v0cascades-q-a");
  if (!dir)
  {
    printf("\nMAIN DIRECTORY NOT FOUND!\n");
    return;
  }
  TDirectoryFile *dirEvt = (TDirectoryFile *)dir->Get("histos-eve");
  if (!dirEvt)
  {
    printf("\nEVENTS DIRECTORY NOT FOUND!\n");
    return;
  }
  TDirectoryFile *dirCasc = (TDirectoryFile *)dir->Get("histos-Casc");
  if (!dirCasc)
  {
    printf("\nCASC DIRECTORY NOT FOUND!\n");
    return;
  }
  TDirectoryFile *dirV0 = (TDirectoryFile *)dir->Get("histos-V0");
  if (!dirV0)
  {
    printf("\nV0 DIRECTORY NOT FOUND!\n");
    return;
  }

  // gStyle->SetOptStat(kFALSE);

  const Int_t numPart = 7;
  TString NamehistoInvMass[numPart] = {"InvMassK0S", "InvMassLambda", "InvMassAntiLambda", "InvMassXiPlus", "InvMassXiMinus", "InvMassOmegaPlus", "InvMassOmegaMinus"};
  TString NamePart[numPart] = {"K0S", "Lam", "ALam", "XiPlu", "XiMin", "OmPlu", "OmMin"};
  Float_t MassPart[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
  Float_t MassPartLimits[numPart] = {0.005, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
  Float_t WidthPartLimits[2][numPart] = {{0.003, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001}, {0.01, 0.005, 0.005, 0.005, 0.005, 0.008, 0.008}};

  // histo of generated particles (MC)
  TCanvas *canvasMC[numPart];
  TCanvas *canvasMC_Rdependence[numPart];
  TH3F *fHistGen3D = {nullptr};
  TH2F *fHistGen2D_Rintegrated = {nullptr};
  TH2F *fHistGen2D_ptintegrated = {nullptr};
  TH1F *fHistGen_Rintegrated[numPart];
  TH1F *fHistGen_Rintegrated_Unscaled[numPart];
  TH1F *fHistGen_ptintegrated[numPart];
  TH1F *fHistGen_ptintegrated_Unscaled[numPart];
  TH1F *fHistEffvsPt[numPart];
  TH1F *fHistEffvsRadius[numPart];
  Float_t CountsGenPt[numPart] = {0};
  TH1F *histoPartPerEvent = new TH1F("histoPartPerEvent", "histoPartPerEvent", numPart, 0, numPart);

  cout << "Processing Events" << endl;
  TH1F *henum = (TH1F *)dirEvt->Get("hEventCounter");
  if (!henum)
  {
    cout << "henum histo not found! " << endl;
    return;
  }
  Double_t NEvents = henum->GetEntries();

  if (isMC)
  {
    fHistGen3D = (TH3F *)dirEvt->Get("GeneratedParticles");
    if (!fHistGen3D)
      return;
  }

  TFile *f_old = CheckOldPass ? new TFile(OldPassPath, "") : nullptr;
  if (isMassvsRadiusPlots)
    PathOut += "_MassvsRadius";
  TFile *f_out = new TFile(PathOut + ".root", "recreate");

  if (isMC)
  {
    // Radius integrated
    fHistGen2D_Rintegrated = (TH2F *)fHistGen3D->Project3D("yxo");
    fHistGen2D_Rintegrated->SetName("GeneratedParticles2D_RadiusIntegrated");
    fHistGen2D_Rintegrated->Write();

    // pt integrated
    fHistGen2D_ptintegrated = (TH2F *)fHistGen3D->Project3D("zxo");
    fHistGen2D_ptintegrated->SetName("GeneratedParticles2D_ptIntegrated");
    fHistGen2D_ptintegrated->Write();
  }

  // V0 topology
  TCanvas *canvasTopologyV0 = new TCanvas("canvasTopologyV0", "", 1000, 700);
  canvasTopologyV0->Divide(3, 2);
  const Int_t NTopV0variables = 5;
  TString TopVarV0Input[NTopV0variables] = {"V0DCAPosToPV", "V0DCANegToPV", "V0DCAV0Daughters", "CosPA", "V0Radius"};
  TString TopVarV0[NTopV0variables] = {"V0 DCA Pos. To PV", "V0 DCA Neg. To PV", "V0 DCA V0 Daughters", "#it{cos}#theta_{PA}", "#it{R}"};
  TString TopVarV0Unit[NTopV0variables] = {"(cm)", "(cm)", "(#sigma)", "", "(cm)"};
  float TopVarV0CutsRun2[NTopV0variables];
  float TopVarV0Cuts[NTopV0variables];
  const float TopVarV0CutsPbPbRun2[NTopV0variables] = {0.1, 0.1, 1.4, 0.95, 0.9}; // Run 2
  const float TopVarV0CutsPbPb[NTopV0variables] = {0.05, 0.05, 2, 0.9, 0.5};      // Run 3 --> Fix third element (DCAV0Daughters) I do not find value of parameter in SVertexerParams
  const float TopVarV0CutsppRun2[NTopV0variables] = {0.03, 0.03, 2, 0.95, 0.9};   // Run 2
  const float TopVarV0Cutspp[NTopV0variables] = {0.05, 0.05, 2, 0.9, 0.5};        // Run 3
  TLine *SelLineRun2;
  TLine *SelLine;
  for (Int_t i = 0; i < NTopV0variables; i++)
  {
    if (CollType == "PbPb")
      TopVarV0CutsRun2[i] = TopVarV0CutsPbPbRun2[i];
    else
      TopVarV0CutsRun2[i] = TopVarV0CutsppRun2[i];
    if (CollType == "PbPb")
      TopVarV0Cuts[i] = TopVarV0CutsPbPb[i];
    else
      TopVarV0Cuts[i] = TopVarV0Cutspp[i];
  }
  TH1F *fHistTopV0[NTopV0variables];
  for (int var = 0; var < NTopV0variables; var++)
  {
    fHistTopV0[var] = (TH1F *)dirV0->Get(TopVarV0Input[var]);
    if (!fHistTopV0[var])
    {
      cout << "Histo " << TopVarV0Input[var] << " with topo variables of V0 not found" << endl;
      return;
    }
    fHistTopV0[var]->Scale(1. / NEvents);
    fHistTopV0[var]->GetYaxis()->SetRangeUser(0.1 * fHistTopV0[var]->GetMinimum(1.e-10), 10 * fHistTopV0[var]->GetMaximum());
    fHistTopV0[var]->SetTitle(TopVarV0[var]);
    fHistTopV0[var]->GetXaxis()->SetTitle(TopVarV0[var] + " " + TopVarV0Unit[var]);
    fHistTopV0[var]->GetYaxis()->SetTitle("1/N_{ev} Counts");
    canvasTopologyV0->cd(var + 1);
    fHistTopV0[var]->DrawCopy("hist");
    SelLineRun2 = new TLine(TopVarV0CutsRun2[var], fHistTopV0[var]->GetMinimum(), TopVarV0CutsRun2[var], fHistTopV0[var]->GetMaximum());
    SelLineRun2->SetLineColor(kGray);
    SelLineRun2->DrawClone("same");
    SelLine = new TLine(TopVarV0Cuts[var], fHistTopV0[var]->GetMinimum(), TopVarV0Cuts[var], fHistTopV0[var]->GetMaximum());
    SelLine->DrawClone("same");
    checkExactLimit(fHistTopV0[var], TopVarV0Cuts[var], (var != 2) ? true : false, cutCheckLabels);
    setPadOptions(true);
    fHistTopV0[var]->Write();
  }
  canvasTopologyV0->SaveAs(PathOut + ".pdf(");

  // Cascade topology
  cout << "Processing Cascade topology" << endl;
  const int kNCanvases = 5;
  TCanvas *canvasTopologyCasc[kNCanvases];
  const int kNvariblesPerCanvas[kNCanvases] = {2, 6, 6, 4, 4};
  const int kCumvariblesPerCanvas[kNCanvases + 1] = {0, 2, 8, 14, 18, 22};
  const int NTopCascInput = 18;
  const int NTopCascVar = 22;
  TString TopVarCascInput[NTopCascInput] = {"CascPt", "V0Ctau", "CascCosPA", "V0CosPA", "V0CosPAToXi", "CascRadius", "V0Radius", "InvMassLambdaDaughter", "DcaCascDaughters", "DcaV0Daughters", "DcaBachToPV", "DcaV0ToPV", "DcaPosToPV", "DcaNegToPV", "CascyXi", "CascCtauXi", "CascyOmega", "CascCtauOmega"};
  TString TopVarCasc[NTopCascVar] = {"Casc #it{p}_{T}", "V0 #it{c}#tau", "Casc #it{cos}#theta_{PA}", "V0 #it{cos}#theta_{PA}", "V0 #it{cos}#theta_{PA} To Casc", "Casc #it{R}", "V0 #it{R}", "#it{m}_{inv} #Lambda Daughter", "DCA Casc Daughters", "DCA V0 Daughters", "DCA Bach. To PV", "DCA V0 To PV", "DCA Pos. To PV", "DCA Neg. To PV", "#it{y}_{#Xi^{-}}", "#it{c}#tau_{ #Xi^{-}}", "#it{y}_{#Omega^{-}}", "#it{c}#tau_{ #Omega^{-}}", "#it{y}_{#Xi^{+}}", "#it{c}#tau_{ #Xi^{+}}", "#it{y}_{#Omega^{+}}", "#it{c}#tau_{ #Omega^{+}}"};
  TString TopVarCascUnit[NTopCascVar] = {"(GeV/#it{c})", "(cm)", "", "", "", "(cm)", "(cm)", "(GeV/#it{c}^2)", "(cm)", "(cm)", "(cm)", "(cm)", "(cm)", "(cm)", "", "(cm)", "", "(cm)", "", "(cm)", "", "(cm)"};
  float TopVarCascCutsRun2[NTopCascVar];
  float TopVarCascCuts[NTopCascVar];
  const float TopVarCascCutsPbPbRun2[NTopCascVar] = {-100., -100., 0.95, 0.95, -100., 0.5, 0.9, -100., 1.4, 1.4, 0.02, 0.05, 0.1, 0.1, -100., -100., -100., -100., -100., -100., -100., -100.};
  const float TopVarCascCutsppRun2[NTopCascVar] = {-100., -100., 0.95, 0.95, -100., 0.5, 0.9, -100., 2.0, 2.0, 0.05, 0.05, 0.03, 0.03, -100., -100., -100., -100., -100., -100., -100., -100.};
  const float TopVarCascCutsPbPb[NTopCascVar] = {-100., -100., 0.7, 0.8, -100., 0, 0.5, -100., 2.0, 2.0, 0.05, 0, 0.05, 0.05, -100., -100., -100., -100., -100., -100., -100., 100.};
  const float TopVarCascCutspp[NTopCascVar] = {-100., -100., 0.7, 0.8, -100., 0, 0.5, -100., 2.0, 2.0, 0.05, 0, 0.05, 0.05, -100., -100., -100., -100., -100., -100., -100., 100.};
  for (Int_t i = 0; i < NTopCascVar; i++)
  {
    if (CollType == "PbPb")
      TopVarCascCutsRun2[i] = TopVarCascCutsPbPbRun2[i];
    else
      TopVarCascCutsRun2[i] = TopVarCascCutsppRun2[i];
    if (CollType == "PbPb")
      TopVarCascCuts[i] = TopVarCascCutsPbPb[i];
    else
      TopVarCascCuts[i] = TopVarCascCutspp[i];
  }
  Float_t InvMassWindow = 0.012;
  if (CollType == "PbPb")
    InvMassWindow = 0.012;
  else
    InvMassWindow = 0.008;

  const int TopVarCascCutsCheckLimit[NTopCascVar] = {0, 0, 2, 2, 0, 2, 2, 0, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0}; // 0 = limit not checked, 1 = max limit, 2 = min limit
  TH2F *fHistTopCasc2D[NTopCascInput];
  TH1F *fHistTopCasc1D[NTopCascVar];

  for (int i_Canv = 0; i_Canv < kNCanvases; i_Canv++)
  {
    canvasTopologyCasc[i_Canv] = new TCanvas(Form("canvasTopologyCasc_%i", i_Canv), "", 1000, 700);
    if (!i_Canv)
    {
      canvasTopologyCasc[i_Canv]->Divide(2, 1);
    }
    else
    {
      (i_Canv > 2) ? canvasTopologyCasc[i_Canv]->Divide(2, 2) : canvasTopologyCasc[i_Canv]->Divide(3, 2);
    }
    for (int var = kCumvariblesPerCanvas[i_Canv]; var < kCumvariblesPerCanvas[i_Canv + 1]; var++)
    {
      if (var < 18)
      {
        fHistTopCasc2D[var] = (TH2F *)dirCasc->Get(TopVarCascInput[var]);
        if (!fHistTopCasc2D[var])
        {
          cout << "Histo " << TopVarCascInput[var] << " with topo variables of Cascades not found" << endl;
          return;
        }
        fHistTopCasc1D[var] = (var < 14) ? (TH1F *)fHistTopCasc2D[var]->ProjectionX(TopVarCascInput[var]) : (TH1F *)fHistTopCasc2D[var]->ProjectionX(TopVarCascInput[var], 1, 1);
      }
      else
      {
        fHistTopCasc1D[var] = (TH1F *)fHistTopCasc2D[var - 4]->ProjectionX(TopVarCascInput[var - 4], 2, 2);
      }

      canvasTopologyCasc[i_Canv]->cd(var + 1 - kCumvariblesPerCanvas[i_Canv]);
      fHistTopCasc1D[var]->Scale(1. / NEvents);
      if (var > 1 && var < 14 && var != 7)
        fHistTopCasc1D[var]->GetYaxis()->SetRangeUser(0.1 * fHistTopCasc1D[var]->GetMinimum(1.e-10), 10 * fHistTopCasc1D[var]->GetMaximum());
      fHistTopCasc1D[var]->GetXaxis()->SetTitle(TopVarCasc[var] + " " + TopVarCascUnit[var]);
      fHistTopCasc1D[var]->GetYaxis()->SetTitle("1/N_{ev} Counts");
      fHistTopCasc1D[var]->SetTitle(TopVarCasc[var]);
      fHistTopCasc1D[var]->DrawCopy("hist");
      SelLineRun2 = new TLine(TopVarCascCutsRun2[var], fHistTopCasc1D[var]->GetMinimum(), TopVarCascCutsRun2[var], fHistTopCasc1D[var]->GetMaximum());
      SelLineRun2->SetLineColor(kGray);
      SelLineRun2->DrawClone("same");
      SelLine = new TLine(TopVarCascCuts[var], fHistTopCasc1D[var]->GetMinimum(), TopVarCascCuts[var], fHistTopCasc1D[var]->GetMaximum());
      SelLine->SetLineColor(kBlack);
      SelLine->DrawClone("same");
      if (TopVarCascCutsCheckLimit[var])
      {
        checkExactLimit(fHistTopCasc1D[var], TopVarCascCuts[var], TopVarCascCutsCheckLimit[var] - 1, cutCheckLabels);
      }
      if (var == 7)
      {
        checkIntervalWidth(fHistTopCasc1D[var], InvMassWindow, cutCheckLabels);
      }
      if (var > 1 && var != 7 && var <= 13)
      {
        setPadOptions(true);
      }
      else
      {
        setPadOptions(false);
      }
      fHistTopCasc1D[var]->Write();
    }
    (TopologyOnly && i_Canv == (kNCanvases - 1)) ? canvasTopologyCasc[i_Canv]->SaveAs(PathOut + ".pdf)") : canvasTopologyCasc[i_Canv]->SaveAs(PathOut + ".pdf");
  }
  if (TopologyOnly)
    return;

  // Define pt binning
  const Int_t num = 9;
  Int_t numPt = 8;
  Int_t numRadius = 9;
  Int_t numEff = 0;
  Int_t numEffPart[numPart] = {0};
  if (isMassvsRadiusPlots)
    numEff = numRadius;
  else
    numEff = numPt;

  // Float_t NPt[num][2] = {{0., 10.}, {0., 1.}, {1., 2.}, {2., 3.}, {3., 4.}, {4., 6.}, {6., 10.}};
  // TString SPt[num] = {"0.0-10.0", "0.0-1.0", "1.0-2.0", "2.0-3.0", "3.0-4.0", "4.0-6.0", "6.0-10.0"};
  // Float_t PtVector[num] = {0., 1., 2., 3., 4., 6., 10.};
  Float_t NPt[num][2] = {{0., 6.}, {0., 0.7}, {0.7, 1.}, {1., 1.5}, {1.5, 2}, {2., 3.}, {3., 4.}, {4., 6.}};
  TString SPt[num] = {"0.0-6.0", "0.0-0.7", "0.7-1.0", "1.0-1.5", "1.5-2.0", "2.0-3.0", "3.0-4.0", "4.0-6.0"};
  Float_t PtVector[num] = {0., 0.7, 1., 1.5, 2., 3., 4., 6.};

  Float_t NRadius[num][2] = {{0., 50.}, {0., 3.}, {3., 6.}, {6., 9.}, {9., 12.}, {12., 15.}, {15., 20.}, {20., 30.}, {30., 50.}};
  TString SRadius[num] = {"0.0-50.0", "0.0-3.0", "3.0-6.0", "6.0-9.0", "9.0-12.0", "12.0-15.0", "15.0-20.0", "20.0-30.0", "30.0-50.0"};
  Float_t RadiusVector[num] = {0., 3., 6., 9., 12., 15., 20., 30, 50};

  if (isMassvsRadiusPlots)
  {
    for (Int_t i = 0; i < numRadius; i++)
    {
      PtVector[i] = RadiusVector[i];
      SPt[i] = SRadius[i];
      for (Int_t j = 0; j < 2; j++)
      {
        NPt[i][j] = NRadius[i][j];
      }
    }
  }

  // dE/dx measurements of lambda daughters
  cout << "Processing lambda daughters" << endl;

  const int nDaughters = 2;
  TString LambdaDaugters[nDaughters] = {"#pi", "p"};
  TString LambdaDaugtersLong[nDaughters] = {"Pion", "Proton"};
  TH2F *fhistodEdxDaughters2D[nDaughters];
  TCanvas *canvasDaughtersdEdx[nDaughters];
  TH1F *fhistodEdxDaughters1D[num][nDaughters];
  TH1F *fHistPeak[nDaughters];
  TCanvas *canvasLambdaDaughtersvsPt = new TCanvas("canvasLambdaDaughtersvsPt", "canvasLambdaDaughtersvsPt", 1200, 800);
  canvasLambdaDaughtersvsPt->Divide(2, 1);
  for (int iDaughter = 0; iDaughter < nDaughters; iDaughter++)
  {
    fhistodEdxDaughters2D[iDaughter] = (TH2F *)dirV0->Get("Response" + LambdaDaugtersLong[iDaughter] + "FromLambda");
    canvasDaughtersdEdx[iDaughter] = new TCanvas("canvas" + LambdaDaugtersLong[iDaughter] + "dEdx", "canvas" + LambdaDaugtersLong[iDaughter] + "dEdx", 1200, 800);
    canvasDaughtersdEdx[iDaughter]->Divide(3, 2);
    fHistPeak[iDaughter] = new TH1F("TPCNsigmaPeak_" + LambdaDaugtersLong[iDaughter], "TPCNsigmaPeak_" + LambdaDaugtersLong[iDaughter], numPt - 1, PtVector);
    for (Int_t pt = 0; pt < numPt; pt++)
    {
      if (pt == 0)
        continue;
      fhistodEdxDaughters1D[pt - 1][iDaughter] = (TH1F *)fhistodEdxDaughters2D[iDaughter]->ProjectionY(Form(LambdaDaugtersLong[iDaughter] + "dEdx_1D_Pt%i", pt), fhistodEdxDaughters2D[iDaughter]->GetXaxis()->FindBin(NPt[pt][0] + 0.0001), fhistodEdxDaughters2D[iDaughter]->GetXaxis()->FindBin(NPt[pt][1] - 0.0001));
      if (RebinTPC)
        fhistodEdxDaughters1D[pt - 1][iDaughter]->Rebin(RebinTPC);
      fhistodEdxDaughters1D[pt - 1][iDaughter]->SetTitle("TPC n_{#sigma} " + LambdaDaugters[iDaughter] + " from #Lambda " + SPt[pt] + " GeV/c");
      fhistodEdxDaughters1D[pt - 1][iDaughter]->GetXaxis()->SetTitle("n_{#sigma}");
      fhistodEdxDaughters1D[pt - 1][iDaughter]->GetYaxis()->SetTitle("Counts");
      canvasDaughtersdEdx[iDaughter]->cd(pt);
      fhistodEdxDaughters1D[pt - 1][iDaughter]->Draw();
      setPadOptions(false);
      fHistPeak[iDaughter]->SetBinContent(pt, findMaximum(fhistodEdxDaughters1D[pt - 1][iDaughter], 3.));
      fHistPeak[iDaughter]->SetBinError(pt, 0.);
    }
    canvasDaughtersdEdx[iDaughter]->SaveAs(PathOut + ".pdf");
    canvasLambdaDaughtersvsPt->cd(iDaughter + 1);
    fHistPeak[iDaughter]->GetYaxis()->SetRangeUser(-3, 3);
    fHistPeak[iDaughter]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistPeak[iDaughter]->GetYaxis()->SetTitle("N_{#sigma} Peak");
    fHistPeak[iDaughter]->Draw("e");
    setPadOptions(false);
  }
  canvasLambdaDaughtersvsPt->SaveAs(PathOut + ".pdf");

  // V0 and cascade invariant mass
  cout << "Processing V0 and cascade invariant masses" << endl;
  TCanvas *canvasMassvsPt[numPart];
  TCanvas *canvasEta[numPart];
  TCanvas *canvasPhi[numPart];
  TCanvas *canvasPhi_2D[numPart];
  TCanvas *canvasITSNLayer[numPart];
  TCanvas *canvasMass[numPart];
  TCanvas *canvasResultsvsPt[numPart];
  TH3F *fhistoInvMass3DTrue[numPart];
  TH2F *fhistoInvMass2D[numPart];
  TH2F *fhistoInvMass2DRadius[numPart];
  TH2F *fhistoInvMass2DTrue_Rintegrated[numPart];
  TH2F *fhistoInvMass2DTrue_ptintegrated[numPart];
  TH1F *fhistoInvMass1D[numPart][num];
  TH1F *fhistoPt1DTrue[numPart];
  TH1F *fhistoR1DTrue[numPart];

  TH3F *fhistoInvMass3D[numPart];
  TH2F *fhistoMassvsEta1_Eta2Neg[numPart];
  TH2F *fhistoMassvsEta1_Eta2Pos[numPart];
  TH1F *fhistoMassEtaPos[numPart];
  TH1F *fhistoMassEtaNeg[numPart];
  TH1F *fhistoMassEta1Pos_Eta2Neg[numPart];
  TH1F *fhistoMassEta1Neg_Eta2Pos[numPart];

  TH3F *fhistoInvMass3D_Phi[numPart];
  TH2F *fhistoInvMass2D_Phi[numPart];
  TH2F *fhistoMassvsPhi1_Phi2Neg[numPart];
  TH2F *fhistoMassvsPhi1_Phi2Pos[numPart];
  TH1F *fhistoMassPhiPos[numPart];
  TH1F *fhistoMassPhiNeg[numPart];
  TH1F *fhistoMassPhi1Pos_Phi2Neg[numPart];
  TH1F *fhistoMassPhi1Neg_Phi2Pos[numPart];

  TH2F *fhistoInvMass2D_ITSNL[numPart];
  TH3F *fhistoInvMass3D_ITSNL[numPart];

  Float_t min_range_signal[numPart] = {0.485, 1.112, 1.112, 1.315, 1.315, 1.668, 1.668};
  Float_t max_range_signal[numPart] = {0.51, 1.12, 1.12, 1.328, 1.328, 1.677, 1.677};
  Float_t liminf[numPart] = {0.44, 1.09, 1.09, 1.30, 1.30, 1.65, 1.65};
  Float_t limsup[numPart] = {0.58, 1.15, 1.15, 1.35, 1.35, 1.70, 1.70};
  TF1 *gauss[numPart][num];
  TF1 *bkgparab[numPart][num];
  TF1 *total[numPart][num];
  TF1 *totalbis[numPart][num];

  Double_t parGaussParab[numPart][num][6];
  TFitResultPtr fFitResultPtr0[numPart][num];
  TFitResultPtr fFitResultPtr1[numPart][num];

  Double_t Mean[numPart][num];
  Double_t ErrMean[numPart][num];
  Double_t Sigma[numPart][num];
  Double_t ErrSigma[numPart][num];
  Double_t sigmas[numPart][num];
  Double_t sigmab[numPart][num];
  Double_t YieldS[numPart][num];
  Double_t YieldSBis[numPart][num];
  Double_t YieldB[numPart][num];
  Float_t YieldSNotScaled[numPart][num];
  Float_t YieldSBisNotScaled[numPart][num];
  Float_t YieldBNotScaled[numPart][num];
  Double_t IntegralS[numPart][num];
  Double_t IntegralB[numPart][num];
  Double_t ErrYieldS[numPart][num];
  Double_t ErrYieldSBis[numPart][num];
  Double_t ErrYieldB[numPart][num];
  Double_t SSB[numPart][num];
  Double_t ErrSSB[numPart][num];
  Double_t SB[numPart][num];
  Double_t ErrSB[numPart][num];
  Double_t SPlusB[numPart][num] = {0};
  Float_t BinWidth[numPart][num] = {0};
  TString SPtPart[numPart][num];

  TH1F *fHistMean[numPart];
  TH1F *fHistSigma[numPart];
  TH1F *fHistYield[numPart];
  TH1F *fHistSSB[numPart];
  TH1F *fHistSB[numPart];

  for (Int_t part = 0; part < numPart; part++)
  {

    if (part == 0)
    {
      NPt[0][0] = 0.;
      NPt[1][0] = 0.;
      PtVector[0] = 0.;
    }
    else if (part == 1 || part == 2)
    {
      NPt[0][0] = 0.4;
      NPt[1][0] = 0.4;
      PtVector[0] = 0.4;
      SPt[0] = "0.4-6.0";
      SPt[1] = "0.4-0.7";
    }
    else if (part == 3 || part == 4)
    {
      NPt[0][0] = 0.3;
      NPt[1][0] = 0.3;
      PtVector[0] = 0.3;
      SPt[0] = "0.3-6.0";
      SPt[1] = "0.3-0.7";
      numEff = numPt - 1;
    }
    else
    {
      NPt[0][0] = 0.7; // lower limit
      NPt[1][0] = 0.7;
      NPt[2][0] = 1.3;
      NPt[3][0] = 1.8;
      NPt[1][1] = 1.3;
      NPt[2][1] = 1.8;
      NPt[3][1] = 3;
      SPt[0] = "0.7-6.0";
      SPt[1] = "0.7-1.3";
      SPt[2] = "1.3-1.8";
      SPt[3] = "1.8-3.0";
      PtVector[0] = 0.7;
      PtVector[1] = 1.3;
      PtVector[2] = 1.8;
      PtVector[3] = 3;
      numEff = 4;
    }
    numEffPart[part] = numEff;

    // yields of generated particles and selected true particles + efficiencis (only in the case of MC)
    if (isMC)
    {

      Int_t partGen = part;

      canvasMC[part] = new TCanvas(Form("canvasMC%i", part), Form("canvasMC%i", part), 1200, 800);
      canvasMC[part]->Divide(3, 1);
      canvasMC[part]->SetTopMargin(0.1);

      canvasMC_Rdependence[part] = new TCanvas(Form("canvasMC_Rdependence%i", part), Form("canvasMC_Rdependence%i", part), 1200, 800);
      canvasMC_Rdependence[part]->Divide(3, 1);
      canvasMC_Rdependence[part]->SetTopMargin(0.1);

      if (part == 0)
        partGen = 0; // K0s
      else if (part == 1)
        partGen = 2; // Lambda
      else if (part == 2)
        partGen = 4; // Alambda
      else if (part == 3)
        partGen = 8; // Xi+
      else if (part == 4)
        partGen = 6; // Xi-
      else if (part == 5)
        partGen = 12; // Omega+
      else if (part == 6)
        partGen = 10; // Omega-

      fHistGen_Rintegrated[part] = (TH1F *)fHistGen2D_Rintegrated->ProjectionY("GeneratedParticles_Rintegrated_" + NamePart[part], fHistGen2D_Rintegrated->GetXaxis()->FindBin(partGen + 0.5), fHistGen2D_Rintegrated->GetXaxis()->FindBin(partGen + 0.5), "e");
      fHistGen_Rintegrated[part]->Rebin(4);
      if (part == 0)
        fHistGen_Rintegrated[part]->Write();

      fHistGen_ptintegrated[part] = (TH1F *)fHistGen2D_ptintegrated->ProjectionY("GeneratedParticles_ptintegrated_" + NamePart[part], fHistGen2D_ptintegrated->GetXaxis()->FindBin(partGen + 0.5), fHistGen2D_ptintegrated->GetXaxis()->FindBin(partGen + 0.5), "e");
      fHistGen_ptintegrated[part]->Rebin(8);
      if (part == 0)
        fHistGen_ptintegrated[part]->Write();

      if (part < 3)
      {
        fhistoInvMass3DTrue[part] = (TH3F *)dirV0->Get(NamehistoInvMass[part] + "True");
      }
      else
      {
        fhistoInvMass3DTrue[part] = (TH3F *)dirCasc->Get(NamehistoInvMass[part] + "True");
      }
      if (!fhistoInvMass3DTrue[part])
        continue;

      fhistoInvMass2DTrue_Rintegrated[part] = (TH2F *)fhistoInvMass3DTrue[part]->Project3D("zxo");
      fhistoInvMass2DTrue_Rintegrated[part]->SetName(NamehistoInvMass[part] + "True_Rintegrated");
      if (part == 0)
        fhistoInvMass2DTrue_Rintegrated[part]->Write();
      fhistoPt1DTrue[part] = (TH1F *)fhistoInvMass2DTrue_Rintegrated[part]->ProjectionX(NamehistoInvMass[part] + "True_PtProj", 0, -1, "E");
      fhistoPt1DTrue[part]->Rebin(4);
      if (part == 0)
        fhistoPt1DTrue[part]->Write();

      fhistoInvMass2DTrue_ptintegrated[part] = (TH2F *)fhistoInvMass3DTrue[part]->Project3D("zyo");
      fhistoInvMass2DTrue_ptintegrated[part]->SetName(NamehistoInvMass[part] + "True_ptintegrated");
      if (part == 0)
        fhistoInvMass2DTrue_ptintegrated[part]->Write();
      fhistoR1DTrue[part] = (TH1F *)fhistoInvMass2DTrue_ptintegrated[part]->ProjectionX(NamehistoInvMass[part] + "True_RProj", 0, -1, "E");
      fhistoR1DTrue[part]->Rebin(8);
      if (part == 0)
        fhistoR1DTrue[part]->Write();

      fHistEffvsPt[part] = (TH1F *)fhistoPt1DTrue[part]->Clone("EfficiencyvsPt_" + NamePart[part]);
      fHistEffvsPt[part]->Divide(fHistGen_Rintegrated[part]);
      for (Int_t b = 1; b <= fhistoPt1DTrue[part]->GetNbinsX(); b++)
      {
        fHistEffvsPt[part]->SetBinError(b, SetEfficiencyError(fhistoPt1DTrue[part]->GetBinContent(b), fHistGen_Rintegrated[part]->GetBinContent(b)));
      }

      fHistEffvsRadius[part] = (TH1F *)fhistoR1DTrue[part]->Clone("EfficiencyvsRadius_" + NamePart[part]);
      fHistEffvsRadius[part]->Divide(fHistGen_ptintegrated[part]);
      for (Int_t b = 1; b <= fhistoR1DTrue[part]->GetNbinsX(); b++)
      {
        fHistEffvsRadius[part]->SetBinError(b, SetEfficiencyError(fhistoR1DTrue[part]->GetBinContent(b), fHistGen_Rintegrated[part]->GetBinContent(b)));
      }

      //fHistEffvsPt[part]->GetXaxis()->SetRangeUser(PtVector[1], PtVector[numPt - 1]);
      fHistEffvsRadius[part]->GetXaxis()->SetRangeUser(0, 50);
      fhistoPt1DTrue[part]->Scale(1. / NEvents / fhistoPt1DTrue[part]->GetBinWidth(1));
      fhistoR1DTrue[part]->Scale(1. / NEvents / fhistoR1DTrue[part]->GetBinWidth(1));
      fHistGen_Rintegrated_Unscaled[part] = (TH1F*)fHistGen_Rintegrated[part]->Clone("fHistGen_Rintegrated_Unscaled_" + NamePart[part]);
      fHistGen_ptintegrated_Unscaled[part] = (TH1F*)fHistGen_ptintegrated[part]->Clone("fHistGen_ptintegrated_Unscaled_" + NamePart[part]);
      fHistGen_Rintegrated[part]->Scale(1. / NEvents / fHistGen_Rintegrated[part]->GetBinWidth(1));
      fHistGen_ptintegrated[part]->Scale(1. / NEvents / fHistGen_ptintegrated[part]->GetBinWidth(1));

      canvasMC[part]->cd(1);
      canvasMC[part]->cd(1)->SetLeftMargin(0.25);
      canvasMC[part]->cd(1)->SetRightMargin(0.1);
      canvasMC[part]->cd(1)->SetBottomMargin(0.15);
      canvasMC[part]->cd(1)->SetTopMargin(0.1);
      setPadOptions(false);
      setHistGraphics(fHistGen_Rintegrated[part], false);
      fHistGen_Rintegrated[part]->GetYaxis()->SetTitle("1/N_{ev} dN/dp_{T}");
      fHistGen_Rintegrated[part]->GetYaxis()->SetTitleOffset(1.);
      fHistGen_Rintegrated[part]->SetTitle("Generated " + NamePart[part] + " Yield");
      fHistGen_Rintegrated[part]->Draw("e0");
      fHistGen_Rintegrated[part]->Write();
      fHistGen_Rintegrated_Unscaled[part]->Write();

      canvasMC[part]->cd(2);
      canvasMC[part]->cd(2)->SetLeftMargin(0.25);
      canvasMC[part]->cd(2)->SetRightMargin(0.1);
      canvasMC[part]->cd(2)->SetBottomMargin(0.15);
      canvasMC[part]->cd(2)->SetTopMargin(0.1);
      setPadOptions(false);
      setHistGraphics(fhistoPt1DTrue[part], false);
      fhistoPt1DTrue[part]->GetYaxis()->SetTitle("1/N_{ev} dN/dp_{T}");
      fhistoPt1DTrue[part]->GetYaxis()->SetTitleOffset(1.2);
      fhistoPt1DTrue[part]->SetTitle("Selected true " + NamePart[part] + " Yield");
      fhistoPt1DTrue[part]->Draw("e0");
      fhistoPt1DTrue[part]->Write();

      canvasMC[part]->cd(3);
      canvasMC[part]->cd(3)->SetLeftMargin(0.25);
      canvasMC[part]->cd(3)->SetRightMargin(0.1);
      canvasMC[part]->cd(3)->SetBottomMargin(0.15);
      canvasMC[part]->cd(3)->SetTopMargin(0.1);
      canvasMC[part]->cd(3)->SetGridy();
      setPadOptions(false);
      setHistGraphics(fHistEffvsPt[part], false);
      fHistEffvsPt[part]->GetYaxis()->SetRangeUser(0, 1.);
      fHistEffvsPt[part]->SetTitle("Efficiency " + NamePart[part]);
      fHistEffvsPt[part]->Draw("e0");
      fHistEffvsPt[part]->Write();

      canvasMC[part]->SaveAs(PathOut + ".pdf");

      cout << "\n\e[35mParticle:\e[39m " << NamePart[part] << endl;
      cout << "Number of events " << NEvents << endl;
      cout << "#particl/event: " << fHistGen_Rintegrated[part]->Integral("width") << " " << fHistGen_ptintegrated[part]->Integral("width") << " " << endl;
      //histoPartPerEvent->SetBinContent(part + 1, fHistGen_Rintegrated[part]->Integral("width"));
      histoPartPerEvent->SetBinContent(part + 1, fHistGen_Rintegrated_Unscaled[part]->Integral("")/NEvents);
      histoPartPerEvent->SetBinError(part + 1, sqrt(fHistGen_Rintegrated_Unscaled[part]->Integral())/NEvents);
      histoPartPerEvent->GetXaxis()->SetBinLabel(part + 1, NamePart[part]);

      canvasMC_Rdependence[part]->cd(1);
      canvasMC_Rdependence[part]->cd(1)->SetLeftMargin(0.25);
      canvasMC_Rdependence[part]->cd(1)->SetRightMargin(0.1);
      canvasMC_Rdependence[part]->cd(1)->SetBottomMargin(0.15);
      canvasMC_Rdependence[part]->cd(1)->SetTopMargin(0.1);
      setPadOptions(false);
      setHistGraphics(fHistGen_ptintegrated[part], false);
      fHistGen_ptintegrated[part]->GetYaxis()->SetTitle("1/N_{ev} dN/dR");
      fHistGen_ptintegrated[part]->GetXaxis()->SetTitle("Radius (cm)");
      fHistGen_ptintegrated[part]->SetTitle("Generated " + NamePart[part] + " Yield");
      fHistGen_ptintegrated[part]->Draw("e0");
      fHistGen_ptintegrated[part]->Write();
      fHistGen_ptintegrated_Unscaled[part]->Write();

      canvasMC_Rdependence[part]->cd(2);
      canvasMC_Rdependence[part]->cd(2)->SetLeftMargin(0.25);
      canvasMC_Rdependence[part]->cd(2)->SetRightMargin(0.1);
      canvasMC_Rdependence[part]->cd(2)->SetBottomMargin(0.15);
      canvasMC_Rdependence[part]->cd(2)->SetTopMargin(0.1);
      setPadOptions(false);
      setHistGraphics(fhistoR1DTrue[part], false);
      fhistoR1DTrue[part]->GetYaxis()->SetTitle("1/N_{ev} dN/dR");
      fhistoR1DTrue[part]->GetXaxis()->SetTitle("Radius (cm)");
      fhistoR1DTrue[part]->SetTitle("Selected true " + NamePart[part] + " Yield");
      fhistoR1DTrue[part]->Draw("e0");
      fhistoR1DTrue[part]->Write();

      canvasMC_Rdependence[part]->cd(3);
      canvasMC_Rdependence[part]->cd(3)->SetLeftMargin(0.25);
      canvasMC_Rdependence[part]->cd(3)->SetRightMargin(0.1);
      canvasMC_Rdependence[part]->cd(3)->SetBottomMargin(0.15);
      canvasMC_Rdependence[part]->cd(3)->SetTopMargin(0.1);
      canvasMC_Rdependence[part]->cd(3)->SetGridy();
      setPadOptions(false);
      setHistGraphics(fHistEffvsRadius[part], false);
      fHistEffvsRadius[part]->GetYaxis()->SetRangeUser(0, 1.);
      fHistEffvsRadius[part]->GetXaxis()->SetTitle("Radius (cm)");
      fHistEffvsRadius[part]->SetTitle("Efficiency vs decay radius " + NamePart[part]);
      fHistEffvsRadius[part]->Draw("e0");
      fHistEffvsRadius[part]->Write();

      canvasMC_Rdependence[part]->SaveAs(PathOut + ".pdf");
    }

    fHistMean[part] = new TH1F("Mean_" + NamePart[part], "Mean_" + NamePart[part], numEff - 1, PtVector);
    fHistSigma[part] = new TH1F("Sigma_" + NamePart[part], "Sigma_" + NamePart[part], numEff - 1, PtVector);
    fHistYield[part] = new TH1F("Yield_" + NamePart[part], "Yield_" + NamePart[part], numEff - 1, PtVector);
    fHistSSB[part] = new TH1F("Purity_" + NamePart[part], "Purity_" + NamePart[part], numEff - 1, PtVector);
    fHistSB[part] = new TH1F("SOverB_" + NamePart[part], "SOverB_" + NamePart[part], numEff - 1, PtVector);
    fHistYield[part]->Sumw2();

    canvasMass[part] = new TCanvas(Form("canvasMass%i", part), NamehistoInvMass[part], 1200, 800);
    canvasMassvsPt[part] = new TCanvas(Form("canvasMassPt%i", part), NamehistoInvMass[part], 1200, 800);
    if (numEff > 7)
      canvasMassvsPt[part]->Divide(3, 3);
    else
      canvasMassvsPt[part]->Divide(3, 2);
    canvasResultsvsPt[part] = new TCanvas(Form("canvasResults%i", part), NamePart[part], 1200, 800);
    canvasResultsvsPt[part]->Divide(2, 2);

    TString NameHisto = "";
    if (isMassvsRadiusPlots)
      NameHisto = "_Radius";
    else
      NameHisto = "";

    if (part < 3)
    {
      fhistoInvMass2D[part] = (TH2F *)dirV0->Get(NamehistoInvMass[part] + NameHisto);
    }
    else
    {
      fhistoInvMass2D[part] = (TH2F *)dirCasc->Get(NamehistoInvMass[part] + NameHisto);
    }
    if (!fhistoInvMass2D[part])
      continue;

    // Fit loop
    cout << "Processing fits" << endl;
    for (Int_t pt = 0; pt < numEff; pt++)
    {
      fhistoInvMass1D[part][pt] = (TH1F *)fhistoInvMass2D[part]->ProjectionY(NamehistoInvMass[part] + Form("1D_Pt%i", pt), fhistoInvMass2D[part]->GetXaxis()->FindBin(NPt[pt][0] + 0.0001), fhistoInvMass2D[part]->GetXaxis()->FindBin(NPt[pt][1] - 0.0001));
      if (isMassvsRadiusPlots)
        fhistoInvMass1D[part][pt]->SetTitle(NamehistoInvMass[part] + " " + SPt[pt] + " cm");
      else
        fhistoInvMass1D[part][pt]->SetTitle(NamehistoInvMass[part] + " " + SPt[pt] + " GeV/c");

      (pt) ? canvasMassvsPt[part]->cd(pt) : canvasMass[part]->cd();
      setPadOptions(false);
      if (pt && ((SkipCascFits == 1 && part > 4) || (SkipCascFits == 2 && part > 2)))
      {
        // fhistoInvMass1D[part][pt]->Draw();
        continue;
      }

      fhistoInvMass1D[part][pt]->GetXaxis()->SetTitle("inv. mass (GeV/c^{2})");
      fhistoInvMass1D[part][pt]->GetYaxis()->SetTitle("Counts");
      fhistoInvMass1D[part][pt]->SetMarkerColor(1);
      fhistoInvMass1D[part][pt]->SetMarkerStyle(21);
      fhistoInvMass1D[part][pt]->SetLineColor(1);
      fhistoInvMass1D[part][pt]->Rebin(2);
      fhistoInvMass1D[part][pt]->Draw("same ep");

      // Fit the invariant mass with gaussian + pol
      gauss[part][pt] = new TF1("gauss_" + NamePart[part] + Form("_Pt%i", pt), "gaus", min_range_signal[part], max_range_signal[part]);
      gauss[part][pt]->SetLineColor(kRed);
      gauss[part][pt]->SetParName(0, "norm");
      gauss[part][pt]->SetParName(1, "mean");
      gauss[part][pt]->SetParName(2, "sigma");
      gauss[part][pt]->SetParameters(1., MassPart[part], 0.0002);
      gauss[part][pt]->SetParLimits(0, 0., 1.1 * fhistoInvMass1D[part][pt]->GetBinContent(fhistoInvMass1D[part][pt]->GetMaximumBin()));
      gauss[part][pt]->SetParLimits(1, min_range_signal[part], max_range_signal[part]);
      gauss[part][pt]->SetParLimits(2, 0.0001, 0.05);

      bkgparab[part][pt] = new TF1("parab_" + NamePart[part] + Form("_Pt%i", pt), fparab, liminf[part], limsup[part], 4);
      bkgparab[part][pt]->SetLineColor(kGreen);
      bkgparab[part][pt]->FixParameter(3, part);

      total[part][pt] = new TF1("total_" + NamePart[part] + Form("_Pt%i", pt), "gaus(0)+pol2(3)", liminf[part], limsup[part]);
      total[part][pt]->SetLineColor(881);
      total[part][pt]->SetLineWidth(1);
      total[part][pt]->SetParName(0, "norm");
      total[part][pt]->SetParName(1, "mean");
      total[part][pt]->SetParName(2, "sigma");
      total[part][pt]->SetParName(3, "p0");
      total[part][pt]->SetParName(4, "p1");
      total[part][pt]->SetParName(5, "p2");

      total[part][pt]->SetParameters(1., MassPart[part], 0.0002);
      total[part][pt]->SetParLimits(0, 0., 1.1 * fhistoInvMass1D[part][pt]->GetBinContent(fhistoInvMass1D[part][pt]->GetMaximumBin()));
      total[part][pt]->SetParLimits(1, min_range_signal[part], max_range_signal[part]);
      total[part][pt]->SetParLimits(2, 0.0001, 0.05);
      fhistoInvMass1D[part][pt]->Fit(gauss[part][pt], "R0Q");
      fhistoInvMass1D[part][pt]->Fit(bkgparab[part][pt], "R0Q");

      gauss[part][pt]->GetParameters(&parGaussParab[part][pt][0]);
      bkgparab[part][pt]->GetParameters(&parGaussParab[part][pt][3]);
      total[part][pt]->SetParameters(parGaussParab[part][pt]);

      fFitResultPtr0[part][pt] = fhistoInvMass1D[part][pt]->Fit(total[part][pt], "SRB+Q");
      if (fFitResultPtr0[part][pt] < 0.)
        continue; // AIMERICcondition in case fFitResultPtr0[part][pt] is empty

      if (!pt)
      {
        TLatex parValues;
        parValues.DrawLatexNDC(0.585, 0.8, Form("#scale[0.65]{#mu = (%.7f #pm %.7f) GeV/c^{2}}", total[part][pt]->GetParameter(1), total[part][pt]->GetParError(1)));
        parValues.DrawLatexNDC(0.585, 0.75, Form("#scale[0.65]{#sigma = (%.7f #pm %.7f) GeV/c^{2}}", total[part][pt]->GetParameter(2), total[part][pt]->GetParError(2)));
        checkBoundaries(total[part][pt]->GetParameter(1), (1. - MassPartLimits[part]) * MassPart[part], (1. + MassPartLimits[part]) * MassPart[part], total[part][pt]->GetParameter(2), WidthPartLimits[0][part], WidthPartLimits[1][part], cutCheckLabels);

        continue;
      }

      totalbis[part][pt] = (TF1 *)total[part][pt]->Clone("totalbis_" + NamePart[part] + Form("_Pt%i", pt));
      fFitResultPtr1[part][pt] = fFitResultPtr0[part][pt];
      total[part][pt]->FixParameter(3, 0);
      total[part][pt]->FixParameter(4, 0);
      total[part][pt]->FixParameter(5, 0);
      totalbis[part][pt]->FixParameter(0, 0);
      totalbis[part][pt]->FixParameter(1, 0);
      totalbis[part][pt]->FixParameter(2, 0);

      Mean[part][pt] = total[part][pt]->GetParameter(1);
      ErrMean[part][pt] = total[part][pt]->GetParError(1);
      Sigma[part][pt] = total[part][pt]->GetParameter(2);
      ErrSigma[part][pt] = total[part][pt]->GetParError(2);

      // Comput yields and Purity
      IntegralS[part][pt] = total[part][pt]->Integral(Mean[part][pt] - 3 * Sigma[part][pt], Mean[part][pt] + 3 * Sigma[part][pt]);
      IntegralB[part][pt] = totalbis[part][pt]->Integral(Mean[part][pt] - 3 * Sigma[part][pt], Mean[part][pt] + 3 * Sigma[part][pt]);
      sigmas[part][pt] = total[part][pt]->IntegralError(Mean[part][pt] - 3 * Sigma[part][pt], Mean[part][pt] + 3 * Sigma[part][pt], fFitResultPtr0[part][pt]->GetParams(), (fFitResultPtr0[part][pt]->GetCovarianceMatrix()).GetMatrixArray());
      sigmab[part][pt] = totalbis[part][pt]->IntegralError(Mean[part][pt] - 3 * Sigma[part][pt], Mean[part][pt] + 3 * Sigma[part][pt], fFitResultPtr1[part][pt]->GetParams(), (fFitResultPtr1[part][pt]->GetCovarianceMatrix()).GetMatrixArray());

      YieldS[part][pt] = IntegralS[part][pt] / fhistoInvMass1D[part][pt]->GetBinWidth(1);
      YieldB[part][pt] = IntegralB[part][pt] / fhistoInvMass1D[part][pt]->GetBinWidth(1);
      ErrYieldS[part][pt] = sigmas[part][pt] / fhistoInvMass1D[part][pt]->GetBinWidth(1);
      ErrYieldB[part][pt] = sigmab[part][pt] / fhistoInvMass1D[part][pt]->GetBinWidth(1);

      for (Int_t b = fhistoInvMass1D[part][pt]->GetXaxis()->FindBin(Mean[part][pt] - 3 * Sigma[part][pt]); b <= fhistoInvMass1D[part][pt]->GetXaxis()->FindBin(Mean[part][pt] + 3 * Sigma[part][pt]); b++)
      {
        SPlusB[part][pt] += fhistoInvMass1D[part][pt]->GetBinContent(b);
      }

      YieldSBis[part][pt] = SPlusB[part][pt] - YieldB[part][pt];
      ErrYieldSBis[part][pt] = sqrt(SPlusB[part][pt] + pow(ErrYieldB[part][pt], 2));

      SB[part][pt] = (SPlusB[part][pt] - YieldB[part][pt]) / YieldB[part][pt];
      ErrSB[part][pt] = SB[part][pt] * sqrt(pow(ErrYieldB[part][pt] / YieldB[part][pt], 2) + 1. / SPlusB[part][pt]);

      // SSB[part][pt] = YieldS[part][pt] / SPlusB[part][pt];
      // ErrSSB[part][pt] = SSB[part][pt] * sqrt(pow(ErrYieldS[part][pt] / YieldS[part][pt], 2) + 1. / SPlusB[part][pt]);
      SSB[part][pt] = (SPlusB[part][pt] - YieldB[part][pt]) / SPlusB[part][pt];
      ErrSSB[part][pt] = SSB[part][pt] * sqrt(1. / SPlusB[part][pt] + pow(ErrYieldB[part][pt] / YieldB[part][pt], 2));

      YieldSNotScaled[part][pt] = YieldS[part][pt];
      YieldSBisNotScaled[part][pt] = YieldSBis[part][pt];
      YieldBNotScaled[part][pt] = YieldB[part][pt];

      YieldS[part][pt] = YieldS[part][pt] / fHistYield[part]->GetBinWidth(pt) / NEvents;
      YieldSBis[part][pt] = YieldSBis[part][pt] / fHistYield[part]->GetBinWidth(pt) / NEvents;
      YieldB[part][pt] = YieldB[part][pt] / fHistYield[part]->GetBinWidth(pt) / NEvents;
      ErrYieldS[part][pt] = ErrYieldS[part][pt] / fHistYield[part]->GetBinWidth(pt) / NEvents;
      ErrYieldSBis[part][pt] = ErrYieldSBis[part][pt] / fHistYield[part]->GetBinWidth(pt) / NEvents;
      ErrYieldB[part][pt] = ErrYieldB[part][pt] / fHistYield[part]->GetBinWidth(pt) / NEvents;

      BinWidth[part][pt] = fHistYield[part]->GetBinWidth(pt);
      SPtPart[part][pt] = SPt[pt];

      fHistMean[part]->SetBinContent(pt, Mean[part][pt]);
      fHistMean[part]->SetBinError(pt, ErrMean[part][pt]);
      fHistSigma[part]->SetBinContent(pt, Sigma[part][pt]);
      fHistSigma[part]->SetBinError(pt, ErrSigma[part][pt]);
      fHistYield[part]->SetBinContent(pt, YieldSBis[part][pt]);
      fHistYield[part]->SetBinError(pt, ErrYieldSBis[part][pt]);
      fHistSSB[part]->SetBinContent(pt, SSB[part][pt]);
      fHistSSB[part]->SetBinError(pt, ErrSSB[part][pt]);
      fHistSB[part]->SetBinContent(pt, SB[part][pt]);
      fHistSB[part]->SetBinError(pt, ErrSB[part][pt]);
    }

    canvasMass[part]->SaveAs(PathOut + ".pdf");
    canvasMassvsPt[part]->SaveAs(PathOut + ".pdf");
    if ((SkipCascFits == 1 && part > 4) || (SkipCascFits == 2 && part > 2))
      continue;

    fHistMean[part]->Write();
    fHistSigma[part]->Write();
    fHistYield[part]->Write();
    fHistSSB[part]->Write();
    fHistSB[part]->Write();

    // Plot mass mean value, width and yields and check their values
    setHistGraphics(fHistMean[part], false);
    if (isMassvsRadiusPlots)
      fHistMean[part]->GetXaxis()->SetTitle("R (cm)");
    fHistMean[part]->GetYaxis()->SetTitle("#mu (GeV/c^{2})");
    setHistGraphics(fHistSigma[part], false);
    if (isMassvsRadiusPlots)
      fHistSigma[part]->GetXaxis()->SetTitle("R (cm)");
    fHistSigma[part]->GetYaxis()->SetTitle("#sigma (GeV/c^{2})");
    setHistGraphics(fHistYield[part], false);
    if (isMassvsRadiusPlots)
      fHistYield[part]->GetXaxis()->SetTitle("R (cm)");
    fHistYield[part]->GetYaxis()->SetTitle("1/N_{ev} dN/dp_{T}");
    setHistGraphics(fHistSSB[part], false);
    if (isMassvsRadiusPlots)
      fHistSSB[part]->GetXaxis()->SetTitle("R (cm)");
    fHistSSB[part]->GetYaxis()->SetTitle("S/(S+B)");
    if (isMassvsRadiusPlots)
      fHistSB[part]->GetXaxis()->SetTitle("R (cm)");
    fHistSB[part]->GetYaxis()->SetTitle("S/B");

    canvasResultsvsPt[part]->cd(1);
    setPadOptions(false);
    if (!CheckOldPass)
    {
      // fHistMean[part]->GetYaxis()->SetRangeUser((1.-2*MassPartLimits[part])*MassPart[part], (1.+2*MassPartLimits[part])*MassPart[part]);
      fHistMean[part]->GetYaxis()->SetRangeUser((1. - 5 * MassPartLimits[part]) * MassPart[part], (1. + 5 * MassPartLimits[part]) * MassPart[part]);
      fHistMean[part]->Draw("e");
    }
    else
    {
      TH1F *fHistMean_old = (TH1F *)f_old->Get("Mean_" + NamePart[part]);
      setHistGraphics(fHistMean_old, true);
      TRatioPlot *rp = new TRatioPlot(fHistMean[part], fHistMean_old);
      rp->SetH1DrawOpt("e");
      rp->Draw();
      rp->GetLowerRefGraph()->GetYaxis()->SetRangeUser(0.98, 1.02);
      ((TH1F *)rp->GetUpperRefObject())->GetYaxis()->SetRangeUser((1. - 2 * MassPartLimits[part]) * MassPart[part], (1. + 2 * MassPartLimits[part]) * MassPart[part]);
      rp->GetUpperPad()->cd();
    }
    TLine PDG_mass = TLine(NPt[0][0], MassPart[part], NPt[0][1], MassPart[part]);
    TLine PDG_mass_Limits[2] = {TLine(NPt[0][0], (1. - MassPartLimits[part]) * MassPart[part], NPt[0][1], (1. - MassPartLimits[part]) * MassPart[part]), TLine(NPt[0][0], (1. + MassPartLimits[part]) * MassPart[part], NPt[0][1], (1. + MassPartLimits[part]) * MassPart[part])};
    PDG_mass.SetLineColor(kBlack);
    PDG_mass.SetLineStyle(2);
    PDG_mass.Draw("same");
    for (int iLim = 0; iLim < 2; iLim++)
    {
      PDG_mass_Limits[iLim].SetLineColor(kBlack);
      PDG_mass_Limits[iLim].SetLineStyle(1);
      PDG_mass_Limits[iLim].Draw("same");
    }
    checkBoundaries(fHistMean[part], (1. - MassPartLimits[part]) * MassPart[part], (1. + MassPartLimits[part]) * MassPart[part], cutCheckLabels);

    canvasResultsvsPt[part]->cd(2);
    setPadOptions(false);
    if (!CheckOldPass)
    {
      fHistSigma[part]->GetYaxis()->SetRangeUser(0., WidthPartLimits[1][part] + 0.005);
      fHistSigma[part]->Draw("e");
    }
    else
    {
      TH1F *fHistSigma_old = (TH1F *)f_old->Get("Sigma_" + NamePart[part]);
      setHistGraphics(fHistSigma_old, true);
      TRatioPlot *rp = new TRatioPlot(fHistSigma[part], fHistSigma_old);
      rp->SetH1DrawOpt("e");
      rp->Draw();
      rp->GetLowerRefGraph()->GetYaxis()->SetRangeUser(0.8, 1.2);
      ((TH1F *)rp->GetUpperRefObject())->GetYaxis()->SetRangeUser(0., WidthPartLimits[1][part] + 0.005);
      rp->GetUpperPad()->cd();
      TLegend *leg = new TLegend(0.50, 0.63, 0.8, 0.73);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(fHistSigma[part], pass_names[1], "PE");
      leg->AddEntry(fHistSigma_old, pass_names[0], "PE");
      leg->Draw();
    }

    TLine Width_Limits[2] = {TLine(NPt[0][0], WidthPartLimits[0][part], NPt[0][1], WidthPartLimits[0][part]), TLine(NPt[0][0], WidthPartLimits[1][part], NPt[0][1], WidthPartLimits[1][part])};
    for (int iLim = 0; iLim < 2; iLim++)
    {
      Width_Limits[iLim].SetLineColor(kBlack);
      Width_Limits[iLim].SetLineStyle(1);
      Width_Limits[iLim].Draw("same");
    }
    checkBoundaries(fHistSigma[part], WidthPartLimits[0][part], WidthPartLimits[1][part], cutCheckLabels);

    canvasResultsvsPt[part]->cd(3);
    setPadOptions(false);
    if (!CheckOldPass)
    {
      fHistYield[part]->Draw("e");
    }
    else
    {
      TH1F *fHistYield_old = (TH1F *)f_old->Get("Yield_" + NamePart[part]);
      setHistGraphics(fHistYield_old, true);
      TRatioPlot *rp = new TRatioPlot(fHistYield[part], fHistYield_old);
      rp->SetH1DrawOpt("e");
      std::vector<double> lines = {0.5, 1., 1.5, 2.};
      rp->SetGridlines(lines);
      rp->Draw();
      rp->GetLowerRefGraph()->GetYaxis()->SetRangeUser(0.5, 2.);
      ((TH1F *)rp->GetUpperRefObject())->GetYaxis()->SetRangeUser(0., 1.1 * findMaxValue(fHistYield[part], fHistYield_old));
      rp->GetUpperPad()->cd();
    }

    canvasResultsvsPt[part]->cd(4);
    setPadOptions(false);
    if (!CheckOldPass)
    {
      fHistSSB[part]->GetYaxis()->SetRangeUser(0., 1.);
      fHistSSB[part]->Draw("e");
    }
    else
    {
      TH1F *fHistSSB_old = (TH1F *)f_old->Get("SSB_" + NamePart[part]);
      setHistGraphics(fHistSSB_old, true);
      TRatioPlot *rp = new TRatioPlot(fHistSSB[part], fHistSSB_old);
      rp->SetH1DrawOpt("e");
      std::vector<double> lines = {0.5, 1., 1.5, 2.};
      rp->SetGridlines(lines);
      rp->Draw();
      rp->GetLowerRefGraph()->GetYaxis()->SetRangeUser(0.5, 2.);
      ((TH1F *)rp->GetUpperRefObject())->GetYaxis()->SetRangeUser(0., 1.1 * findMaxValue(fHistSSB[part], fHistSSB_old));
      rp->GetUpperPad()->cd();
    }

    canvasResultsvsPt[part]->SaveAs(PathOut + ".pdf");

    // mass vs eta daughters
    if (isMassvsRadiusPlots)
    {
      if (part > 2)
        continue;
      canvasEta[part] = new TCanvas(Form("canvasEta%i", part), NamehistoInvMass[part], 1200, 800);
      canvasEta[part]->Divide(2, 2);
      fhistoInvMass3D[part] = (TH3F *)dirV0->Get(NamehistoInvMass[part] + "_EtaDaughters");
      if (!fhistoInvMass3D[part])
        continue;
      // eta1 = pos daughter, eta2 = neg daughter
      fhistoInvMass3D[part]->GetYaxis()->SetRange(fhistoInvMass3D[part]->GetYaxis()->FindBin(0.001), fhistoInvMass3D[part]->GetYaxis()->FindBin(0.999));
      fhistoMassvsEta1_Eta2Pos[part] = (TH2F *)fhistoInvMass3D[part]->Project3D("zxo"); // mass vs Eta1 (only pos Eta2)
      fhistoMassvsEta1_Eta2Pos[part]->SetName(NamehistoInvMass[part] + "_Eta2Pos");
      fhistoMassvsEta1_Eta2Pos[part]->SetTitle("Mass vs Eta of positive daughter (only for positive eta value of the other daughter)");
      fhistoMassEtaPos[part] = (TH1F *)fhistoMassvsEta1_Eta2Pos[part]->ProjectionY("hInvMass_EtaPos" + NamePart[part], fhistoMassvsEta1_Eta2Pos[part]->GetXaxis()->FindBin(0.001), fhistoMassvsEta1_Eta2Pos[part]->GetXaxis()->FindBin(0.999));
      fhistoMassEta1Neg_Eta2Pos[part] = (TH1F *)fhistoMassvsEta1_Eta2Pos[part]->ProjectionY("hInvMass_Eta1NegEta2Pos" + NamePart[part], fhistoMassvsEta1_Eta2Pos[part]->GetXaxis()->FindBin(-0.999), fhistoMassvsEta1_Eta2Pos[part]->GetXaxis()->FindBin(-0.001));

      fhistoInvMass3D[part]->GetYaxis()->SetRange(fhistoInvMass3D[part]->GetYaxis()->FindBin(-0.999), fhistoInvMass3D[part]->GetYaxis()->FindBin(-0.001));
      fhistoMassvsEta1_Eta2Neg[part] = (TH2F *)fhistoInvMass3D[part]->Project3D("zxo"); // mass vs Eta1 (only neg Eta2)
      fhistoMassvsEta1_Eta2Neg[part]->SetName(NamehistoInvMass[part] + "_Eta2Neg");
      fhistoMassvsEta1_Eta2Neg[part]->SetTitle("Mass vs Eta of positive daughter (only for negative eta value of the other daughter)");
      fhistoMassEtaNeg[part] = (TH1F *)fhistoMassvsEta1_Eta2Neg[part]->ProjectionY("hInvMass_EtaNeg" + NamePart[part], fhistoMassvsEta1_Eta2Pos[part]->GetXaxis()->FindBin(-0.999), fhistoMassvsEta1_Eta2Pos[part]->GetXaxis()->FindBin(-0.001));
      fhistoMassEta1Pos_Eta2Neg[part] = (TH1F *)fhistoMassvsEta1_Eta2Neg[part]->ProjectionY("hInvMass_Eta1PosEta2Neg" + NamePart[part], fhistoMassvsEta1_Eta2Pos[part]->GetXaxis()->FindBin(0.001), fhistoMassvsEta1_Eta2Pos[part]->GetXaxis()->FindBin(0.999));

      fhistoMassEta1Pos_Eta2Neg[part]->SetTitle("Invariant mass " + NamePart[part] + " opposite sign eta daughters (#eta_{+} > 0, #eta_{-} < 0)");
      fhistoMassEta1Neg_Eta2Pos[part]->SetTitle("Invariant mass " + NamePart[part] + " opposite sign eta daughters (#eta_{+} < 0, #eta_{-} > 0)");
      fhistoMassEtaPos[part]->SetTitle("Invariant mass " + NamePart[part] + " same sign eta daughters (#eta_{+} > 0, #eta_{-} > 0)");
      fhistoMassEtaNeg[part]->SetTitle("Invariant mass " + NamePart[part] + " same sign eta daughters (#eta_{+} < 0, #eta_{-} < 0)");

      canvasEta[part]->cd(1);
      fhistoMassEtaPos[part]->Draw();
      canvasEta[part]->cd(2);
      fhistoMassEtaNeg[part]->Draw();
      canvasEta[part]->cd(3);
      fhistoMassEta1Pos_Eta2Neg[part]->Draw();
      canvasEta[part]->cd(4);
      fhistoMassEta1Neg_Eta2Pos[part]->Draw();

      canvasEta[part]->SaveAs(PathOut + ".pdf");

      // mass vs phi daughters
      canvasPhi[part] = new TCanvas(Form("canvasPhi%i", part), NamehistoInvMass[part], 1200, 800);
      canvasPhi[part]->Divide(2, 2);

      fhistoInvMass3D_Phi[part] = (TH3F *)dirV0->Get(NamehistoInvMass[part] + "_PhiDaughters");
      if (!fhistoInvMass3D_Phi[part])
        continue;
      // phi1 = pos daughter, phi2 = neg daughter
      fhistoInvMass3D_Phi[part]->GetYaxis()->SetRange(fhistoInvMass3D_Phi[part]->GetYaxis()->FindBin(0.001), fhistoInvMass3D_Phi[part]->GetYaxis()->FindBin(TMath::Pi()));
      fhistoMassvsPhi1_Phi2Pos[part] = (TH2F *)fhistoInvMass3D_Phi[part]->Project3D("zxo"); // mass vs Phi1 (only pos Phi2)
      fhistoMassPhiPos[part] = (TH1F *)fhistoMassvsPhi1_Phi2Pos[part]->ProjectionY("hInvMass_PhiPos" + NamePart[part], fhistoMassvsPhi1_Phi2Pos[part]->GetXaxis()->FindBin(0.001), fhistoMassvsPhi1_Phi2Pos[part]->GetXaxis()->FindBin(TMath::Pi()));
      fhistoMassPhi1Neg_Phi2Pos[part] = (TH1F *)fhistoMassvsPhi1_Phi2Pos[part]->ProjectionY("hInvMass_Phi1NegPhi2Pos" + NamePart[part], fhistoMassvsPhi1_Phi2Pos[part]->GetXaxis()->FindBin(TMath::Pi()), fhistoMassvsPhi1_Phi2Pos[part]->GetXaxis()->FindBin(2 * TMath::Pi()));

      fhistoInvMass3D_Phi[part]->GetYaxis()->SetRange(fhistoInvMass3D_Phi[part]->GetYaxis()->FindBin(TMath::Pi()), fhistoInvMass3D_Phi[part]->GetYaxis()->FindBin(2 * TMath::Pi()));
      fhistoMassvsPhi1_Phi2Neg[part] = (TH2F *)fhistoInvMass3D_Phi[part]->Project3D("zxo"); // mass vs Phi1 (only neg Phi2)
      fhistoMassPhiNeg[part] = (TH1F *)fhistoMassvsPhi1_Phi2Neg[part]->ProjectionY("hInvMass_PhiNeg" + NamePart[part], fhistoMassvsPhi1_Phi2Pos[part]->GetXaxis()->FindBin(TMath::Pi()), fhistoMassvsPhi1_Phi2Pos[part]->GetXaxis()->FindBin(2 * TMath::Pi()));
      fhistoMassPhi1Pos_Phi2Neg[part] = (TH1F *)fhistoMassvsPhi1_Phi2Neg[part]->ProjectionY("hInvMass_Phi1PosPhi2Neg" + NamePart[part], fhistoMassvsPhi1_Phi2Pos[part]->GetXaxis()->FindBin(0.001), fhistoMassvsPhi1_Phi2Pos[part]->GetXaxis()->FindBin(TMath::Pi()));

      fhistoMassPhi1Pos_Phi2Neg[part]->SetTitle("Invariant mass " + NamePart[part] + " opposite sign Phi daughters ( 0 < #phi_{+} < #pi, #pi < #phi_{-} < 2#pi)");
      fhistoMassPhi1Neg_Phi2Pos[part]->SetTitle("Invariant mass " + NamePart[part] + " opposite sign Phi daughters (#pi< #phi_{+} < 2#pi, 0 < #phi_{-} < #pi)");
      fhistoMassPhiPos[part]->SetTitle("Invariant mass " + NamePart[part] + " same sign Phi daughters ( 0 < #phi_{+} < #pi, 0 < #phi_{-} < #pi)");
      fhistoMassPhiNeg[part]->SetTitle("Invariant mass " + NamePart[part] + " same sign Phi daughters (#pi< #phi_{+} < 2#pi, #pi < #phi_{-} < 2#pi)");

      canvasPhi[part]->cd(1);
      fhistoMassPhiPos[part]->Draw();
      canvasPhi[part]->cd(2);
      fhistoMassPhiNeg[part]->Draw();
      canvasPhi[part]->cd(3);
      fhistoMassPhi1Pos_Phi2Neg[part]->Draw();
      canvasPhi[part]->cd(4);
      fhistoMassPhi1Neg_Phi2Pos[part]->Draw();

      canvasPhi[part]->SaveAs(PathOut + ".pdf");

      canvasPhi_2D[part] = new TCanvas(Form("canvasPhi_2D%i", part), NamehistoInvMass[part], 1200, 800);
      fhistoInvMass3D_Phi[part]->GetYaxis()->SetRange(1, fhistoInvMass3D_Phi[part]->GetNbinsX());
      fhistoInvMass2D_Phi[part] = (TH2F *)fhistoInvMass3D_Phi[part]->Project3D("yx"); // mass vs Phi1 (only pos Phi2)

      fhistoInvMass2D_Phi[part]->GetYaxis()->SetTitle("#phi Neg Daughter");
      fhistoInvMass2D_Phi[part]->GetXaxis()->SetTitle("#phi Pos Daughter");

      canvasPhi_2D[part]->cd(1);
      fhistoInvMass2D_Phi[part]->Draw("colz");
      canvasPhi_2D[part]->SaveAs(PathOut + ".pdf");

      // Number of ITS layers of daughters
      canvasITSNLayer[part] = new TCanvas(Form("canvasITSNLayer%i", part), NamehistoInvMass[part], 1200, 800);
      fhistoInvMass3D_ITSNL[part] = (TH3F *)dirV0->Get(NamehistoInvMass[part] + "_ITSMapDaughters");

      fhistoInvMass2D_ITSNL[part] = (TH2F *)fhistoInvMass3D_ITSNL[part]->Project3D("yx");
      fhistoInvMass2D_ITSNL[part]->GetYaxis()->SetTitle("# ITS Layers Neg Daughter");
      fhistoInvMass2D_ITSNL[part]->GetXaxis()->SetTitle("# ITS Layers Pos Daughter");

      canvasITSNLayer[part]->cd(1);
      fhistoInvMass2D_ITSNL[part]->Draw("colz");
      canvasITSNLayer[part]->SaveAs(PathOut + ".pdf");
    }
  }

  // draw ratio between particle and antiparticle yields
  TCanvas *canvasRatioYieldsPAntiP;
  canvasRatioYieldsPAntiP = new TCanvas("canvasRatioYieldsPAntiP", "Ratio particle/antiparticle yields", 1200, 800);
  if (SkipCascFits != 2)
    canvasRatioYieldsPAntiP->Divide(3 - SkipCascFits, 1);
  TH1F *fHistRatioYields[3];
  TString NameRatioPart[3] = {"Lambda", "Xi", "Omega"};
  Float_t RatioLowerBound[3] = {0.75, 0.85, 0.85};
  Float_t RatioUpperBound[3] = {1.1, 1.1, 1.1};
  TLine *LineOne = new TLine(0.5, 1, 10, 1);

  for (Int_t i = 0; i < 3; i++)
  {
    if ((SkipCascFits == 1 && i > 1) || (SkipCascFits == 2 && i > 0))
      continue;
    if (i != 0)
    {
      fHistRatioYields[i] = (TH1F *)fHistYield[2 * i + 1]->Clone("fHistYieldRatio_" + NameRatioPart[i]);
      fHistRatioYields[i]->Divide(fHistYield[2 * i + 2]);
      fHistRatioYields[i]->SetTitle("yield ratio " + NamePart[2 * i + 1] + "/" + NamePart[2 * i + 2]);
    }
    else
    {
      fHistRatioYields[i] = (TH1F *)fHistYield[2 * i + 2]->Clone("fHistYieldRatio_" + NameRatioPart[i]);
      fHistRatioYields[i]->Divide(fHistYield[2 * i + 1]);
      fHistRatioYields[i]->SetTitle("yield ratio " + NamePart[2 * i + 2] + "/" + NamePart[2 * i + 1]);
    }
    canvasRatioYieldsPAntiP->cd(i + 1);
    fHistRatioYields[i]->SetLineColor(628);
    fHistRatioYields[i]->SetMarkerColor(628);
    fHistRatioYields[i]->GetYaxis()->SetRangeUser(0.7, 1.2);
    fHistRatioYields[i]->GetYaxis()->SetTitle("Ratio");
    fHistRatioYields[i]->Draw("e");
    //  checkBoundaries(fHistRatioYields[i], RatioLowerBound[i], RatioUpperBound[i], cutCheckLabels);
    LineOne->Draw("same");
    setPadOptions(false);
  }
  if (isMC)
    canvasRatioYieldsPAntiP->SaveAs(PathOut + ".pdf");
  else
    canvasRatioYieldsPAntiP->SaveAs(PathOut + ".pdf)");

  f_out->Close();
  cout << "Mission accomplished! " << endl;

  Float_t RelUnc = 0.05;
  cout << "Estimates! " << endl;
  cout << "How many events do we need for a small statistical uncertainty " << endl;
  cout << "on the spectrum in the largest pt interval(the one with less counts)? " << endl;
  cout << " Which uncertainty do you want? " << endl;
  // cin >> RelUnc;
  Float_t EstimateNevt = 0;
  Float_t EstimateNevtCorrect = 0;
  Float_t SvsRelUnc = 0;
  for (Int_t part = 0; part < numPart; part++)
  {
    // if (part != 3 && part != 4 && part!=0 )
    if (part == 5 || part == 6)
      continue;
    cout << "\nParticle: " << NamePart[part] << endl;
    for (Int_t pt = 0; pt < numEffPart[part]; pt++)
    {
      if (pt == 0)
        continue;
      SvsRelUnc = (1. + 2. / SB[part][pt]) / pow(RelUnc, 2);
      EstimateNevtCorrect = SvsRelUnc * 1. / BinWidth[part][pt] * 1. / YieldSBis[part][pt];
      EstimateNevt = 1. / pow(RelUnc, 2) * SB[part][pt] / (SB[part][pt] + 1) * 1. / BinWidth[part][pt] * 1. / YieldSBis[part][pt];
      cout << "\nPt: " << SPtPart[part][pt] << " (" << BinWidth[part][pt] << ")" << endl;
      cout << "SigmaN/N " << 1. / sqrt(SPlusB[part][pt]) << " SigmaS/S (correct) " << sqrt(YieldBNotScaled[part][pt] + SPlusB[part][pt]) / YieldSBisNotScaled[part][pt] << endl;
      cout << "Rel error on yield: " << ErrYieldSBis[part][pt] / YieldSBis[part][pt] << endl;
      cout << "Ratio (approx / correct) " << 1. / sqrt(SPlusB[part][pt]) / (sqrt(YieldBNotScaled[part][pt] + SPlusB[part][pt]) / YieldSBisNotScaled[part][pt]) << endl;
      cout << "Number of events needed for " << RelUnc << " uncertainty: " << EstimateNevt << " better: " << EstimateNevtCorrect << endl;
      cout << "Actual number of events in this period: " << NEvents << endl;
    }
  }

  TCanvas *canvasPartPerEvent = new TCanvas("canvasPartPerEvent", "canvasPartPerEvent", 800, 600);
  canvasPartPerEvent->cd();

  histoPartPerEvent->SetLineColor(628);
  histoPartPerEvent->SetMarkerColor(628);
  histoPartPerEvent->GetYaxis()->SetTitle("# particles / event");
  if (isMC)
  {
    histoPartPerEvent->GetYaxis()->SetRangeUser(0, 0.5);
    histoPartPerEvent->Draw("e");
    canvasPartPerEvent->SaveAs(PathOut + ".pdf)");
  }
}

void checkExactLimit(TH1F *lHist, Float_t limit, bool isMinLimit, TLatex labels[])
{
  Int_t lastbin = lHist->GetNbinsX();
  Float_t binEdge = 0;
  if (isMinLimit)
  {
    binEdge = 0; // should be zero if the minum is underflow
    for (Int_t i = 1; i < lastbin; i++)
    {
      if (lHist->GetBinContent(i) != 0)
        break;
      binEdge = lHist->GetBinLowEdge(i + 1); // minimum edge of histogram
    }
  }
  else
  {
    binEdge = lHist->GetBinLowEdge(lastbin); // should be lastbin if the maximum is overflow
    for (Int_t i = lastbin; i >= 0; i--)
    {
      if (lHist->GetBinContent(i) != 0)
        break;
      binEdge = lHist->GetBinLowEdge(i);
    }
  }
  if (binEdge < limit)
  {
    isMinLimit ? labels[2].Draw() : labels[1].Draw();
    isMinLimit ? labels[4].Draw() : labels[3].Draw();
  }
  else if (binEdge > limit)
  {
    isMinLimit ? labels[1].Draw() : labels[2].Draw();
    isMinLimit ? labels[3].Draw() : labels[4].Draw();
  }
  else
  {
    labels[0].Draw();
  }
}
void checkIntervalWidth(TH1F *lHist, float width, TLatex labels[])
{
  Int_t lastbin = lHist->GetNbinsX();
  Float_t binEdges[2] = {0., 0.};
  for (Int_t i = 1; i < lastbin; i++)
  {
    if (lHist->GetBinContent(i) != 0)
      break;
    binEdges[0] = lHist->GetBinLowEdge(i + 1);
  }
  for (Int_t i = lastbin; i >= 0; i--)
  {
    if (lHist->GetBinContent(i) != 0)
      break;
    binEdges[1] = lHist->GetBinLowEdge(i);
  }
  if (binEdges[1] - binEdges[0] > width + lHist->GetBinWidth(1) + 0.000001)
  {
    labels[2].Draw();
    labels[4].Draw();
  }
  else if (binEdges[1] - binEdges[0] < width)
  {
    labels[1].Draw();
    labels[3].Draw();
  }
  else
  {
    labels[0].Draw();
  }
}
void checkBoundaries(TH1F *lHist, Float_t lower_bound, Float_t upper_bound, TLatex labels[])
{
  bool isOk = true;
  Int_t lastbin = lHist->GetNbinsX();
  Float_t content = 0;
  for (Int_t i = 1; i <= lastbin; i++)
  {
    content = lHist->GetBinContent(i);
    if (content >= lower_bound && content <= upper_bound)
      continue;
    isOk = false;
    break;
  }

  if (isOk)
  {
    labels[0].Draw();
  }
  else
  {
    labels[5].Draw();
    labels[6].Draw();
  }
}
void checkBoundaries(Float_t value1, Float_t lower_bound1, Float_t upper_bound1, Float_t value2, Float_t lower_bound2, Float_t upper_bound2, TLatex labels[])
{
  if (value1 > lower_bound1 && value1 < upper_bound1 && value2 > lower_bound2 && value2 < upper_bound2)
  {
    labels[7].Draw();
  }
  else
  {
    labels[5].Draw();
  }
}
void setPadOptions(bool useLogY = false)
{
  if (useLogY)
    gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
}
void setHistGraphics(TH1F *lHist, bool isOldPass)
{
  isOldPass ? lHist->SetLineColor(kBlue) : lHist->SetLineColor(kRed);
  isOldPass ? lHist->SetMarkerColor(kBlue) : lHist->SetMarkerColor(kRed);
  lHist->SetMarkerStyle(33);
  lHist->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  lHist->GetXaxis()->SetTitleSize(0.05);
  lHist->GetYaxis()->SetTitleSize(0.05);
}
float findMaximum(TH1F *lHist, float width = 3.)
{
  lHist->GetXaxis()->SetRange(lHist->GetXaxis()->FindBin(-width), lHist->GetXaxis()->FindBin(width));
  float max = lHist->GetXaxis()->GetBinCenter(lHist->GetMaximumBin());
  TLine *ar1 = new TLine(0.5, 1000000, 0.5, 0.);
  TLine *lineMax = new TLine(max, 0., max, lHist->GetMaximum());
  lineMax->SetLineStyle(kDashed);
  lineMax->Draw();
  return max;
}
float findMaxValue(TH1F *lHist1, TH1F *lHist2)
{
  float max1 = lHist1->GetMaximum();
  float max2 = lHist2->GetMaximum();
  return max1 > max2 ? max1 : max2;
}

Double_t SetEfficiencyError(Int_t k, Int_t n)
{
  return sqrt(((Double_t)k + 1) * ((Double_t)k + 2) / (n + 2) / (n + 3) - pow((Double_t)(k + 1), 2) / pow(n + 2, 2));
}
