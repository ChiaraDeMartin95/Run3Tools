//----- Selectivity OSS errors Selectivity = n/N sigma = rad(n)/N; ratio 2025/2024 error = function, with FullCorr=1
#include <iostream>
#include <TAttLine.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TVirtualPad.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TSystem.h>
#include <string>
#include <TROOT.h>

const Float_t StrLimit = 4.3E-5;
const Float_t LFLimit = 5E-5;
const Int_t numTriggers = 12; // Interesting trigger
Int_t iscale = 0;             // 0 if input file has no mistakes. At some point axis labels where shifted by one bin and in those cases iscale = 1 should be set.

void ErrRatioCorr(TH1F *hNum, TH1F *hDenom, TH1F *hRatio, Int_t FullCorr);

void QAplots_runbyrun(string period = "LHC26ac_batch2")
{

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  Long_t *dummy1 = 0, *dummy2 = 0, *dummy3 = 0, *dummy4 = 0;

  TStopwatch time;
  time.Start();

  char pass_name[2][30] = {"apass1", "apass1_skimmed"}; // AnalysisResults_fullrun_LHC25ah_ctf_skim_full_564704

  // LHC26ac_batch1
  // const int nruns_2025 = 20;
  // int runnumber_2025[nruns_2025] = {569945, 569947, 569978, 569980, 569981, 569982, 570011, 570012, 570025, 570036, 570049, 570050, 570051, 570054, 570065, 570066, 570077, 570079, 570091, 570102};
  // LHC26ac_batch2
  const int nruns_2025 = 14;
  int runnumber_2025[nruns_2025] = {570128, 570143, 570158, 570159, 570160, 570162, 570163, 570165, 570166, 570191, 570205, 570243, 570279, 570292};

  //------------ read files
  TFile *file_in2025[nruns_2025] = {0x0};
  // TFile *file_in_skimmed[nruns] = {0x0};

  for (int irun = 0; irun < nruns_2025; irun++)
  {
    if (gSystem->GetPathInfo(Form("../TriggerForRun3/EventFiltering2026/LHC26ac_batch2/AnalysisResults_fullrun_LHC26ac_apass1_%d.root", runnumber_2025[irun]), dummy1, dummy2, dummy3, dummy4) != 0)
      cout << "File Not Found! Try again" << endl;
    file_in2025[irun] = new TFile(Form("../TriggerForRun3/EventFiltering2026/LHC26ac_batch2/AnalysisResults_fullrun_LHC26ac_apass1_%d.root", runnumber_2025[irun]), "read");
    printf("Open File: %s\n", file_in2025[irun]->GetName());
  }

  // color palette
  const int NRGBs = 5;
  double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red[NRGBs] = {0.00, 0.00, 0.78, 1.00, 0.51};
  double green[NRGBs] = {0.00, 0.81, 0.90, 0.20, 0.00};
  double blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};
  int FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, nruns_2025);
  int color;

  //------------ #events for each filter
  TH1F *hEvSel2025[nruns_2025] = {0x0};
  TH1F *hEvSelFinal2025[nruns_2025] = {0x0};
  TH1F *hMultFT0MNorm[nruns_2025] = {0x0};
  TH1F *hMultFT0MNormRatio[nruns_2025] = {0x0};
  Double_t TotEvtMax2025 = 0;

  for (int i = 0; i < nruns_2025; i++)
  {
    hEvSel2025[i] = (TH1F *)file_in2025[i]->Get("lf-strangeness-filter/hProcessedEvents");
    hEvSel2025[i]->SetName(Form("hEvSel2025_%d", i));
    if (hEvSel2025[i]->GetBinContent(1) > TotEvtMax2025)
      TotEvtMax2025 = hEvSel2025[i]->GetBinContent(1);
    hMultFT0MNorm[i] = (TH1F *)file_in2025[i]->Get("lf-strangeness-filter/EventsvsMultiplicity/AllEventsvsMultiplicityFT0MNorm");
    hMultFT0MNorm[i]->Scale(1. / hEvSel2025[i]->GetBinContent(1));
    hMultFT0MNorm[i]->Rebin(4);
  }

  TCanvas *cMultFT0MNorm = new TCanvas("cMultFT0MNorm", "cMultFT0MNorm", 800, 800);
  cMultFT0MNorm->cd();
  cMultFT0MNorm->SetMargin(0.15, 0.15, 0.13, 0.08); // left,right,bottom,top
  cMultFT0MNorm->SetLogy();
  for (int i = 0; i < nruns_2025; i++)
  {
    color = nruns_2025 - i - 1 + FI;
    hMultFT0MNorm[i]->SetLineColor(color);
    hMultFT0MNorm[i]->SetTitle("T0M distribution normalised to processed events");
    hMultFT0MNorm[i]->Draw("hist same");
  }
  TLegend *legendMultFT0MNorm = new TLegend(0.2, 0.55, 0.51, 0.9);
  legendMultFT0MNorm->SetFillStyle(0);
  legendMultFT0MNorm->SetBorderSize(0);
  legendMultFT0MNorm->SetTextAlign(12);
  legendMultFT0MNorm->SetTextSize(0.025);
  for (int i = 0; i < nruns_2025; i++)
  {
    legendMultFT0MNorm->AddEntry(hMultFT0MNorm[i], Form("Run %d", runnumber_2025[i]), "l");
  }
  legendMultFT0MNorm->Draw();
  cMultFT0MNorm->SaveAs(Form("cMultFT0MNorm_%s.pdf", period.c_str()));
  cMultFT0MNorm->SaveAs(Form("cMultFT0MNorm_%s.png", period.c_str()));
  cMultFT0MNorm->SaveAs(Form("cMultFT0MNorm_%s.eps", period.c_str()));

  TCanvas *cMultFT0MNormRatio = new TCanvas("cMultFT0MNormRatio", "cMultFT0MNormRatio", 800, 800);
  cMultFT0MNormRatio->cd();
  cMultFT0MNormRatio->SetMargin(0.15, 0.15, 0.13, 0.08); // left,right,bottom,top
  for (int i = 0; i < nruns_2025; i++)
  {
    hMultFT0MNormRatio[i] = (TH1F *)hMultFT0MNorm[i]->Clone(Form("hMultFT0MNormRatio_%d", i));
    hMultFT0MNormRatio[i]->Divide(hMultFT0MNorm[0]);
    hMultFT0MNormRatio[i]->Draw("hist same");
  }
  legendMultFT0MNorm->Draw();
  cMultFT0MNormRatio->SaveAs(Form("cMultFT0MNormRatio_%s.pdf", period.c_str()));
  cMultFT0MNormRatio->SaveAs(Form("cMultFT0MNormRatio_%s.png", period.c_str()));
  cMultFT0MNormRatio->SaveAs(Form("cMultFT0MNormRatio_%s.eps", period.c_str()));

  for (int i = 0; i < nruns_2025; i++)
  {
    color = nruns_2025 - i - 1 + FI;
    hEvSelFinal2025[i] = (TH1F *)hEvSel2025[0]->Clone(Form("hEvSelFinal2025_%d", i));
    for (int j = 1; j <= hEvSelFinal2025[i]->GetNbinsX(); j++)
    {
      hEvSelFinal2025[i]->SetBinContent(j, hEvSel2025[i]->GetBinContent(j));
      hEvSelFinal2025[i]->SetBinError(j, TMath::Sqrt(hEvSel2025[i]->GetBinContent(j)));
    }
    hEvSelFinal2025[i]->SetLineColor(color);
  }

  // name of each filter/column
  TString binlabel[21] = {"All events", "Processed events", "Events w/ high-pT hadron", "Omegas", "h-Omega", "2Xi", "3Xi", "4Xi", "Xi-N", "Omega large R", "Xi", "Tracked Xi", "Tracked Omega", "HighMultFT0M+Omega", "Double Omega", "Omega+Xi", "Lam+Lam", "HighMultTrack+Omega", "HighMultFT0M", "HighMultTrack", "sigma-p"};

  for (int i = 0; i < nruns_2025; i++)
  {
    for (Int_t j = 1; j <= hEvSelFinal2025[i]->GetNbinsX(); j++)
    {
      hEvSelFinal2025[i]->GetXaxis()->SetBinLabel(j, binlabel[j - 1]);
    }
  }

  TH1F *hEvSel_filters2025[nruns_2025] = {0x0};

  for (int i = 0; i < nruns_2025; i++)
  {
    hEvSel_filters2025[i] = new TH1F(Form("hEvSel_filters%d", i), ";;Events", numTriggers - 1, 0, numTriggers - 1);
    hEvSel_filters2025[i]->SetBinContent(1, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Omegas") + iscale));
    hEvSel_filters2025[i]->SetBinContent(2, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("h-Omega") + iscale));
    hEvSel_filters2025[i]->SetBinContent(3, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Xi-N") + iscale));
    hEvSel_filters2025[i]->SetBinContent(4, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("HighMultFT0M+Omega") + iscale));
    hEvSel_filters2025[i]->SetBinContent(5, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Tracked Omega") + iscale));
    hEvSel_filters2025[i]->SetBinContent(6, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Double Omega") + iscale));
    hEvSel_filters2025[i]->SetBinContent(7, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Omega+Xi") + iscale));
    hEvSel_filters2025[i]->SetBinContent(8, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Lam+Lam") + iscale));
    hEvSel_filters2025[i]->SetBinContent(9, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("HighMultTrack+Omega") + iscale));
    hEvSel_filters2025[i]->SetBinContent(10, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("HighMultFT0M") + iscale));
    hEvSel_filters2025[i]->SetBinContent(11, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("HighMultTrack") + iscale));
    hEvSel_filters2025[i]->SetBinContent(12, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("sigma-p") + iscale));
  }

  cout << "hleoo" << endl;
  //------------ #events for each filter
  TCanvas *cfull = new TCanvas("cfull", "cfull", 800, 800);
  cfull->cd();
  cfull->SetMargin(0.15, 0.15, 0.13, 0.08); // left,right,bottom,top
  cfull->SetLogy();

  TH1F *h = new TH1F(Form("h%d", 0), ";;Events", 21, -1, 20);
  h->SetTitle(Form("%s", period.c_str()));
  h->GetXaxis()->SetBinLabel(1, "All events");
  h->GetXaxis()->SetBinLabel(2, "Processed events");
  h->GetXaxis()->SetBinLabel(3, "Events w/ high-pT hadron");
  h->GetXaxis()->SetBinLabel(4, "Omega");
  h->GetXaxis()->SetBinLabel(5, "h-Omega");
  h->GetXaxis()->SetBinLabel(6, "2Xi");
  h->GetXaxis()->SetBinLabel(7, "3Xi");
  h->GetXaxis()->SetBinLabel(8, "4Xi");
  h->GetXaxis()->SetBinLabel(9, "Xi-N");
  h->GetXaxis()->SetBinLabel(10, "Omega large R");
  h->GetXaxis()->SetBinLabel(11, "Xi");
  h->GetXaxis()->SetBinLabel(12, "Tracked Xi");
  h->GetXaxis()->SetBinLabel(13, "Tracked Omega");
  h->GetXaxis()->SetBinLabel(14, "HighMultFT0M+Omega");
  h->GetXaxis()->SetBinLabel(15, "Double Omega");
  h->GetXaxis()->SetBinLabel(16, "Omega+Xi");
  h->GetXaxis()->SetBinLabel(17, "Lam+Lam");
  h->GetXaxis()->SetBinLabel(18, "HighMultTrack+Omega");
  h->GetXaxis()->SetBinLabel(19, "HighMultFT0M");
  h->GetXaxis()->SetBinLabel(20, "HighMultTrack");
  h->GetXaxis()->SetBinLabel(21, "sigma-p");
  h->SetStats(0);
  h->GetYaxis()->SetRangeUser(0.5, 1.5 * TotEvtMax2025);
  h->Draw();
  TLegend *leg = new TLegend(0.4, 0.8, 0.85, 0.91);
  leg->SetTextSize(0.02);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  cout << "hleoo" << endl;

  for (int i = 0; i < nruns_2025; i++)
  {
    color = nruns_2025 - i - 1 + FI;
    for (Int_t j = 1; j <= hEvSelFinal2025[i]->GetNbinsX(); j++)
      hEvSelFinal2025[i]->GetXaxis()->SetBinLabel(j, h->GetXaxis()->GetBinLabel(j));

    hEvSelFinal2025[i]->GetYaxis()->SetRangeUser(0.5, 1.5 * TotEvtMax2025);
    hEvSelFinal2025[i]->SetLineColor(color);
    hEvSelFinal2025[i]->SetMarkerColor(color);
    hEvSelFinal2025[i]->SetMarkerSize(1);
    hEvSelFinal2025[i]->SetMarkerStyle(20);
    hEvSelFinal2025[i]->DrawCopy("SAME e");
    leg->AddEntry(hEvSelFinal2025[i], "config2025", "pl");
    leg->Draw("same");
  }
  // cfull->SaveAs(Form("cfull_%s_total.png",period.c_str()));

  //------------ selectivity: #events for each filter/#total events (all processed events)
  TH1F *hSelectivity2025[numTriggers][nruns_2025];
  Double_t tot_filters2025[numTriggers];
  Double_t tot_processed2025 = 0.;
  TH1F *hSelectivity_period2025[numTriggers];
  //"Omegas", "h-Omega", "Xi-N", "HighMultFT0M+Omega", "Tracked Omega", "Double Omega", "Omega+Xi",
  //"Lam+Lam", "HighMultTrack+Omega", "HighMultFT0M", "HighMultTrack", "sigma-p"
  Float_t LowLimit[numTriggers] = {2e-6, 6e-6, 1.3e-6, 6e-6, 5e-6, 3e-7, 1.5e-6, 2e-7, 1e-6, 2e-5, 1e-6, 4.5e-6};
  Float_t UpLimit[numTriggers] = {3e-4, 1.5e-5, 2e-6, 10e-6, 9e-6, 6e-7, 3.5e-6, 6e-7, 4e-6, 10e-5, 7e-6, 6.5e-6};
  for (int ifilter = 0; ifilter < numTriggers; ifilter++)
  {
    for (int irun = 0; irun < nruns_2025; irun++)
    {
      hSelectivity2025[ifilter][irun] = new TH1F(Form("Selectivity2025_filter_%d_run_%d", ifilter, runnumber_2025[irun]), Form("Selectivity2026_filter_%d_run_%d", ifilter, runnumber_2025[irun]), nruns_2025, 0, nruns_2025);

      hSelectivity2025[ifilter][irun]->SetBinContent(irun + 1, hEvSel_filters2025[irun]->GetBinContent(ifilter + 1) / hEvSelFinal2025[irun]->GetBinContent(hEvSelFinal2025[irun]->GetXaxis()->FindBin("Processed events")));
      tot_filters2025[ifilter] += hEvSel_filters2025[irun]->GetBinContent(ifilter + 1);
      tot_processed2025 += hEvSelFinal2025[irun]->GetBinContent(hEvSelFinal2025[irun]->GetXaxis()->FindBin("Processed events"));
    }

    hSelectivity_period2025[ifilter] = new TH1F(Form("Selectivity2026_filter_%d_period_%s", ifilter, period.c_str()), Form("Selectivity2026_filter_%d_period_%s", ifilter, period.c_str()), numTriggers - 1, 0, numTriggers - 1);
    hSelectivity_period2025[ifilter]->SetBinContent(ifilter + 1, tot_filters2025[ifilter] / tot_processed2025);
  }

  TH1F *hSel_allfilters2025 = new TH1F(Form("hSel2026_allfilters_run_%d", runnumber_2025[0]), Form("hSel2026_allfilters_run_%d", runnumber_2025[0]), numTriggers, 0, numTriggers);
  for (int ifilter = 0; ifilter < numTriggers; ifilter++)
  {
    hSel_allfilters2025->SetBinContent(ifilter + 1, hSelectivity2025[ifilter][0]->GetBinContent(1));
    hSel_allfilters2025->SetBinError(ifilter + 1, TMath::Sqrt(hSelectivity2025[ifilter][0]->GetBinContent(1)) / tot_processed2025);
  }

  TString *xlabel[numTriggers];
  for (int ifilter = 0; ifilter < numTriggers; ifilter++)
  {
    if (ifilter == 0)
      xlabel[ifilter] = new TString("Omega");
    if (ifilter == 1)
      xlabel[ifilter] = new TString("hOmega");
    if (ifilter == 2)
      xlabel[ifilter] = new TString("Xi-N");
    if (ifilter == 3)
      xlabel[ifilter] = new TString("HighMultFT0M+Omega");
    if (ifilter == 4)
      xlabel[ifilter] = new TString("Tracked Omega");
    if (ifilter == 5)
      xlabel[ifilter] = new TString("Double Omega");
    if (ifilter == 6)
      xlabel[ifilter] = new TString("Omega+Xi");
    if (ifilter == 7)
      xlabel[ifilter] = new TString("Lam+Lam");
    if (ifilter == 8)
      xlabel[ifilter] = new TString("HighMultTrack+Omega");
    if (ifilter == 9)
      xlabel[ifilter] = new TString("HighMultFT0M");
    if (ifilter == 10)
      xlabel[ifilter] = new TString("HighMultTrack");
    if (ifilter == 11)
      xlabel[ifilter] = new TString("sigma-p");
  }

  TCanvas *cselallfilters = new TCanvas("cselallfilters", "cselallfilters", 1400, 1400);
  cselallfilters->cd();
  cselallfilters->SetMargin(0.15, 0.05, 0.13, 0.01); // left,right,bottom,top
  cselallfilters->SetLogy();
  TLegend *leg_sel = new TLegend(0.25, 0.9, 0.85, 0.95);
  leg_sel->SetTextSize(0.02);
  leg_sel->SetTextFont(42);
  leg_sel->SetBorderSize(0);
  leg_sel->SetNColumns(2);

  TH1D *frame = new TH1D("frame", "", 8, 0, 8);
  frame->GetYaxis()->SetTitle("Selectivity");
  frame->SetStats(0);
  for (int ifilter = 1; ifilter < 9; ifilter++)
  {
    frame->GetXaxis()->SetBinLabel(ifilter, xlabel[ifilter - 1]->Data());
    frame->GetYaxis()->SetRangeUser(hSel_allfilters2025->GetMinimum() * 0.2, hSel_allfilters2025->GetMaximum() * 1.5); // LowLimit[ifilter],UpLimit[ifilter]);
  }
  frame->Draw("");

  TLine *linelimit = new TLine(0, 5e-5, 8, 5e-5);
  linelimit->SetLineColor(kBlack);
  linelimit->SetLineWidth(2);
  linelimit->SetLineStyle(2);
  linelimit->Draw("same");
  hSel_allfilters2025->SetTitle("Selectivity");
  hSel_allfilters2025->GetYaxis()->SetRangeUser(hSel_allfilters2025->GetMinimum() * 0.2, hSel_allfilters2025->GetMaximum() * 1.5);
  hSel_allfilters2025->SetLineColor(kRed);
  hSel_allfilters2025->SetMarkerColor(kRed);
  hSel_allfilters2025->SetMarkerSize(3);
  hSel_allfilters2025->SetMarkerStyle(20);
  hSel_allfilters2025->DrawCopy("PSAME");
  leg_sel->AddEntry(hSel_allfilters2025, "2026 ac batch2", "pl");
  leg_sel->Draw("same");
  cselallfilters->SaveAs("cselallfilters_total_2026ac_batch2.png");

  //------------ selectivity for each filter for each run
  TCanvas *cselectivity[numTriggers];
  TH1F *hselectivity_2025[numTriggers];
  TFile *fileout = new TFile(Form("QAplots_runbyrun_%s.root", period.c_str()), "recreate");

  for (int ifilter = 0; ifilter < numTriggers; ifilter++)
  {
    cselectivity[ifilter] = new TCanvas(Form("cselectivity_filter%d", ifilter), Form("cselectivity_filter%d", ifilter), 800, 800);
    cselectivity[ifilter]->cd();
    cselectivity[ifilter]->SetMargin(0.13, 0.01, 0.13, 0.08); // left,right,bottom,top
    // cselectivity[ifilter]->SetLogy();
    cselectivity[ifilter]->SetGridy();
    cselectivity[ifilter]->SetGridx();

    hselectivity_2025[ifilter] = new TH1F(Form("hselectivity_filter%d", ifilter), ";;Trigger Selectivity", nruns_2025, 0, nruns_2025);
    if (ifilter == 0)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: Omega", period.c_str()));
    if (ifilter == 1)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: hOmega", period.c_str()));
    if (ifilter == 2)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: Xi-N", period.c_str()));
    if (ifilter == 3)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: HighMultFT0M+Omega", period.c_str()));
    if (ifilter == 4)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: Tracked Omega", period.c_str()));
    if (ifilter == 5)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: Double Omega", period.c_str()));
    if (ifilter == 6)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: Omega+Xi", period.c_str()));
    if (ifilter == 7)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: Lam+Lam", period.c_str()));
    if (ifilter == 8)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: HighMultTrack+Omega", period.c_str()));
    if (ifilter == 9)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: HighMultFT0M", period.c_str()));
    if (ifilter == 10)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: HighMultTrack", period.c_str()));
    if (ifilter == 11)
      hselectivity_2025[ifilter]->SetTitle(Form("%s: sigma-p", period.c_str()));

    for (int irun = 0; irun < nruns_2025; irun++)
      hselectivity_2025[ifilter]->GetXaxis()->SetBinLabel(irun + 1, Form("%d", runnumber_2025[irun]));
    hselectivity_2025[ifilter]->GetYaxis()->SetRangeUser(LowLimit[ifilter], UpLimit[ifilter]);
    hselectivity_2025[ifilter]->Draw();

    cout << "Selectivity for filter: " << hselectivity_2025[ifilter]->GetTitle() << endl;
    float runAvg = 0;
    int runCounter = 0;
    for (int irun = 0; irun < nruns_2025; irun++)
    {
      color = nruns_2025 - irun - 1 + FI;
      cout << "For run: " << runnumber_2025[irun] << " " << hSelectivity2025[ifilter][irun]->GetBinContent(irun + 1) << endl;
      runAvg += hSelectivity2025[ifilter][irun]->GetBinContent(irun + 1);
      runCounter++;
      hSelectivity2025[ifilter][irun]->GetYaxis()->SetRangeUser(LowLimit[ifilter], UpLimit[ifilter]);
      hSelectivity2025[ifilter][irun]->SetLineColor(color);
      hSelectivity2025[ifilter][irun]->SetMarkerColor(color);
      hSelectivity2025[ifilter][irun]->SetMarkerSize(1.3);
      hSelectivity2025[ifilter][irun]->SetMarkerStyle(20);
      hSelectivity2025[ifilter][irun]->DrawCopy("PSAME");
    }
    runAvg /= runCounter;
    cout << "Average selectivity for filter " << ifilter << ": " << runAvg << endl;
    TF1 *fAvg = new TF1(Form("fAvg_filter%d", ifilter), "[0]", 0, nruns_2025);
    fAvg->SetParameter(0, runAvg);
    fAvg->SetLineColor(kMagenta + 2);
    fAvg->SetLineWidth(3);
    fAvg->SetLineStyle(2);
    fAvg->Draw("same");
    TF1 *fAvgPlus = new TF1(Form("fAvgPlus_filter%d", ifilter), "[0]", 0, nruns_2025);
    fAvgPlus->SetParameter(0, runAvg + 0.1 * runAvg);
    fAvgPlus->SetLineColor(kBlack);
    fAvgPlus->SetLineWidth(2);
    fAvgPlus->SetLineStyle(2);
    fAvgPlus->Draw("same");
    TF1 *fAvgMinus = new TF1(Form("fAvgMinus_filter%d", ifilter), "[0]", 0, nruns_2025);
    fAvgMinus->SetParameter(0, runAvg - 0.1 * runAvg);
    fAvgMinus->SetLineColor(kBlack);
    fAvgMinus->SetLineWidth(2);
    fAvgMinus->SetLineStyle(2);
    fAvgMinus->Draw("same");
    if (ifilter == 0)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Omega.png", period.c_str()));
    if (ifilter == 1)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_hOmega.png", period.c_str()));
    if (ifilter == 2)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Xi-N.png", period.c_str()));
    if (ifilter == 3)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_HighMultFT0M+Omega.png", period.c_str()));
    if (ifilter == 4)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_OmegaTracked.png", period.c_str()));
    if (ifilter == 5)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_OmegaDouble.png", period.c_str()));
    if (ifilter == 6)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Omega+Xi.png", period.c_str()));
    if (ifilter == 7)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Lam+Lam.png", period.c_str()));
    if (ifilter == 8)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_HighMultTrack+Omega.png", period.c_str()));
    if (ifilter == 9)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_HighMultFT0M.png", period.c_str()));
    if (ifilter == 10)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_HighMultTrack.png", period.c_str()));
    if (ifilter == 11)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_sigma-p.png", period.c_str()));

    fileout->cd();
    fAvg->Write();
  }
  fileout->Close();
  cout << "I stored the average values for each filter in the file: " << endl;
  cout << fileout->GetName() << endl;

  // cout << hSelectivity_period[0]->GetBinContent(1) << endl;
}

void ErrRatioCorr(TH1F *hNum, TH1F *hDenom, TH1F *hRatio, Int_t FullCorr)
{
  // FullCorr == 1 means ro = 1; //full correlation
  // FullCorr == 0 means ro = sigmaDenom/sigmaNum ; //data sample at numerator is a subset of denominator
  // FullCorr == 2 :another possibility (followed by Fiorella's according to Roger Barlow): I think this is an approximation of Barlow's prescription, and my procedure is better. I found out that errors calculated in this way vary within +-10% from those calculated by me.

  Float_t Err1 = 0;
  Float_t Err2 = 0;
  Float_t ErrC = 0;
  Float_t Err = 0;

  for (Int_t b = 1; b <= hNum->GetNbinsX(); b++)
  {
    if (hNum->GetBinContent(b) == 0 || hDenom->GetBinContent(b) == 0)
    {
      hRatio->SetBinError(b, 0);
      continue;
    }
    Err1 = pow(hNum->GetBinError(b) / hNum->GetBinContent(b), 2);
    Err2 = pow(hDenom->GetBinError(b) / hDenom->GetBinContent(b), 2);
    if (FullCorr == 0)
    {
      if (hDenom->GetBinError(b) < hNum->GetBinError(b))
        ErrC = pow(hDenom->GetBinError(b), 2) / (hNum->GetBinContent(b) * hDenom->GetBinContent(b)); // Num is a subsample of denom
      else
        ErrC = pow(hNum->GetBinError(b), 2) / (hNum->GetBinContent(b) * hDenom->GetBinContent(b)); // Denom is a subsample of num
    }
    else if (FullCorr == 1)
    {
      ErrC = hDenom->GetBinError(b) * hNum->GetBinError(b) / (hNum->GetBinContent(b) * hDenom->GetBinContent(b));
    }
    Err = sqrt(Err1 + Err2 - 2 * ErrC); // are we sure there's the 2?
    if (Err1 + Err2 - ErrC < 0)
    {
      cout << "Error not defined! (NaN)" << endl;
      Err = sqrt(Err1 + Err2);
    }
    hRatio->SetBinError(b, Err * hRatio->GetBinContent(b));
    if (FullCorr == 2)
    {
      hRatio->SetBinError(b, sqrt(TMath::Abs(pow(hNum->GetBinError(b), 2) - pow(hDenom->GetBinError(b), 2))) / hDenom->GetBinContent(b));
    }
    //  cout << "bin: " << b << " err 1: " << Err*hRatio->GetBinContent(b) << " err 2: " <<  hRatio->GetBinError(b) << " 2/1: " << hRatio->GetBinError(b) / Err*hRatio->GetBinContent(b) << endl;
  }
}
