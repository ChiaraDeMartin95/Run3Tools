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
Int_t numBins = 8; // Interesting trigger
Int_t iscale = 0;  // 0 if input file has no mistakes. At some point axis labels where shifted by one bin and in those cases iscale = 1 should be set.

void QAplots(string period = "LHC25am")
{

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  Long_t *dummy1 = 0, *dummy2 = 0, *dummy3 = 0, *dummy4 = 0;

  TStopwatch time;
  time.Start();

  char pass_name[2][30] = {"apass1", "apass1_skimmed"};
  const int nruns = 9;
  int runnumber[nruns] = {566848, 566859, 566871, 566893, 566905, 566906, 566907, 566921, 566935};
  // int runnumber[nruns] = {566947, 566980, 566995, 567005, 567006, 567014, 567017, 567023, 567039, 567040, 567050, 567063, 567078, 567109, 567121, 567122, 567136, 567147, 567148, 567149, 567160, 567174, 567186};

  //------------ read files
  TFile *file_in2025[nruns] = {0x0};
  TFile *file_in_skimmed[nruns] = {0x0};

  // file_in2025[0] = new TFile("AnalysisResults_HY_newtask_v2.root","read");
  // printf("Open File: %s\n",file_in2025[0]->GetName());
  //
  TFile *file_in2024[nruns] = {0x0};
  // file_in2024[0] = new TFile("AnalysisResults_HY_oldtask.root","read");
  //  printf("Open File: %s\n",file_in2024[0]->GetName());

  for (int irun = 0; irun < nruns; irun++)
  {
    // if(gSystem->GetPathInfo(Form("AnalysisResults_fullrun_%s_%s_%d.root",period.c_str(),pass_name[0],runnumber[irun]),dummy1,dummy2,dummy3,dummy4) != 0) cout << "File Not Found! Try again" << endl;
    file_in2024[irun] = new TFile(Form("../TriggerForRun3/EventFiltering2025/AnalysisResults_fullrun_skimmed_LHC25am_batch3_%d.root", runnumber[irun]), "read");
    file_in2025[irun] = new TFile(Form("../TriggerForRun3/EventFiltering2025/AnalysisResults_fullrun_LHC25am_batch3_%d.root", runnumber[irun]), "read");
    // printf("Open File: %s\n",file_in[irun]->GetName());
    //
    file_in_skimmed[irun] = new TFile(Form("../TriggerForRun3/EventFiltering2025/AnalysisResults_fullrun_skimmed_LHC25am_batch3_%d.root", runnumber[irun]), "read");
  }

  // color palette
  const int NRGBs = 5;
  double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red[NRGBs] = {0.00, 0.00, 0.78, 1.00, 0.51};
  double green[NRGBs] = {0.00, 0.81, 0.90, 0.20, 0.00};
  double blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};
  int FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, nruns);
  int color;

  //------------ #events for each filter
  TH1F *hEvSel2025[nruns] = {0x0};
  TH1F *hEvSelFinal2025[nruns] = {0x0};
  Double_t TotEvtMax2025 = 0;

  for (int i = 0; i < nruns; i++)
  {
    hEvSel2025[i] = (TH1F *)file_in2025[i]->Get("lf-strangeness-filter/hProcessedEvents");
    hEvSel2025[i]->SetName(Form("hEvSel2025_%d", i));
    if (hEvSel2025[i]->GetBinContent(1) > TotEvtMax2025)
      TotEvtMax2025 = hEvSel2025[i]->GetBinContent(1);
  }

  for (int i = 0; i < nruns; i++)
  {
    color = nruns - i - 1 + FI;
    hEvSelFinal2025[i] = (TH1F *)hEvSel2025[0]->Clone(Form("hEvSelFinal2025_%d", i));
    for (int j = 1; j <= hEvSelFinal2025[i]->GetNbinsX(); j++)
    {
      hEvSelFinal2025[i]->SetBinContent(j, hEvSel2025[i]->GetBinContent(j));
      hEvSelFinal2025[i]->SetBinError(j, TMath::Sqrt(hEvSel2025[i]->GetBinContent(j)));
    }
    hEvSelFinal2025[i]->SetLineColor(color);
  }

  TH1F *hEvSel2024[nruns] = {0x0};
  TH1F *hEvSelFinal2024[nruns] = {0x0};
  Double_t TotEvtMax2024 = 0;

  for (int i = 0; i < nruns; i++)
  {
    hEvSel2024[i] = (TH1F *)file_in2024[i]->Get("lf-strangeness-filter/hProcessedEvents");
    hEvSel2024[i]->SetName(Form("hEvSel2024_%d", i));
    if (hEvSel2024[i]->GetBinContent(1) > TotEvtMax2024)
      TotEvtMax2024 = hEvSel2024[i]->GetBinContent(1);
  }

  for (int i = 0; i < nruns; i++)
  {
    color = nruns - i - 1 + FI;
    hEvSelFinal2024[i] = (TH1F *)hEvSel2024[0]->Clone(Form("hEvSelFinal2024_%d", i));
    for (int j = 1; j <= hEvSelFinal2024[i]->GetNbinsX(); j++)
    {
      hEvSelFinal2024[i]->SetBinContent(j, hEvSel2024[i]->GetBinContent(j));
      hEvSelFinal2024[i]->SetBinError(j, TMath::Sqrt(hEvSel2024[i]->GetBinContent(j)));
    }
    hEvSelFinal2024[i]->SetLineColor(color);
  }

  //------------ #events for each filter from skimmed data
  TH1F *hEvSel_skimmed[nruns] = {0x0};
  TH1F *hEvSelFinal_skimmed[nruns] = {0x0};

  for (int i = 0; i < nruns; i++)
  {
    hEvSel_skimmed[i] = (TH1F *)file_in_skimmed[i]->Get("lf-strangeness-filter/hProcessedEvents");
    hEvSel_skimmed[i]->SetName(Form("hEvSel_skimmed_%d", i));
  }

  for (int i = 0; i < nruns; i++)
  {
    color = nruns - i - 1 + FI;
    hEvSelFinal_skimmed[i] = (TH1F *)hEvSel_skimmed[0]->Clone(Form("hEvSelFinal_skimmed_%d", i));
    for (int j = 1; j <= hEvSelFinal_skimmed[i]->GetNbinsX(); j++)
    {
      hEvSelFinal_skimmed[i]->SetBinContent(j, hEvSel_skimmed[i]->GetBinContent(j));
      hEvSelFinal_skimmed[i]->SetBinError(j, TMath::Sqrt(hEvSel_skimmed[i]->GetBinContent(j)));
      cout << "skimmed run: " << runnumber[i] << " |bin: " << j << " val: " << hEvSelFinal_skimmed[i]->GetBinContent(j) << endl;
      // if(runnumber[i]==551759) cout << "skimmed run: " << runnumber[i] << " |bin: " << j << " val: " << hEvSelFinal_skimmed[i]->GetBinContent(j) << endl;
    }
    hEvSelFinal_skimmed[i]->SetLineColor(color);
  }

  // name of each filter/column
  char binlabel[18][30] = {"All events", "Processed events", "Events w/ high-pT hadron", "Omegas", "h-Omega", "2Xi", "3Xi", "4Xi", "Xi-N", "Omega large R", "Xi", "Tracked Xi", "Tracked Omega", "HighMult+Omega", "Double Omega", "Omega+Xi", "Lam+Lam"};

  for (int i = 0; i < nruns; i++)
  {
    for (Int_t j = 1; j <= hEvSelFinal2025[i]->GetNbinsX(); j++)
    {
      hEvSelFinal2025[i]->GetXaxis()->SetBinLabel(j, binlabel[j - 1]);
    }

    for (Int_t j = 1; j <= hEvSelFinal_skimmed[i]->GetNbinsX(); j++){
      hEvSelFinal_skimmed[i]->GetXaxis()->SetBinLabel(j,binlabel[j-1]);
    }
  }

  for (int i = 0; i < nruns; i++)
  {
    for (Int_t j = 1; j <= hEvSelFinal2024[i]->GetNbinsX(); j++)
    {
      hEvSelFinal2024[i]->GetXaxis()->SetBinLabel(j, binlabel[j - 1]);
    }

    // for (Int_t j = 1; j <= hEvSelFinal_skimmed[i]->GetNbinsX(); j++){
    //   hEvSelFinal_skimmed[i]->GetXaxis()->SetBinLabel(j,binlabel[j-1]);
    // }
  }

  TH1F *hEvSel_filters2025[nruns] = {0x0};
  TH1F *hEvSel_skimmed_filters[nruns] = {0x0};

  for (int i = 0; i < nruns; i++)
  {
    hEvSel_filters2025[i] = new TH1F(Form("hEvSel_filters%d", i), ";;Events", numBins - 1, 0, numBins - 1);
    hEvSel_filters2025[i]->SetBinContent(1, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Omegas") + iscale));
    hEvSel_filters2025[i]->SetBinContent(2, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("h-Omega") + iscale));
    hEvSel_filters2025[i]->SetBinContent(3, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Xi-N") + iscale));
    hEvSel_filters2025[i]->SetBinContent(4, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("HighMult+Omega") + iscale));
    hEvSel_filters2025[i]->SetBinContent(5, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Tracked Omega") + iscale));
    hEvSel_filters2025[i]->SetBinContent(6, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Double Omega") + iscale));
    hEvSel_filters2025[i]->SetBinContent(7, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Omega+Xi") + iscale));
    hEvSel_filters2025[i]->SetBinContent(8, hEvSelFinal2025[i]->GetBinContent(hEvSelFinal2025[i]->GetXaxis()->FindBin("Lam+Lam") + iscale));
    hEvSel_skimmed_filters[i] = new TH1F(Form("hEvSel_skimmed_filters%d", i), ";;Events", numBins - 1, 0, numBins - 1);
    hEvSel_skimmed_filters[i]->SetBinContent(1, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Omegas") + iscale));
    hEvSel_skimmed_filters[i]->SetBinContent(2, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("h-Omega") + iscale));
    hEvSel_skimmed_filters[i]->SetBinContent(3, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Xi-N") + iscale));
    hEvSel_skimmed_filters[i]->SetBinContent(4, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("HighMult+Omega") + iscale));
    hEvSel_skimmed_filters[i]->SetBinContent(5, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Tracked Omega") + iscale));
    hEvSel_skimmed_filters[i]->SetBinContent(6, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Double Omega") + iscale));
    hEvSel_skimmed_filters[i]->SetBinContent(7, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Omega+Xi") + iscale));
    hEvSel_skimmed_filters[i]->SetBinContent(8, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Lam+Lam") + iscale));
  }

  TH1F *hEvSel_filters2024[nruns] = {0x0};
  // TH1F *hEvSel_skimmed_filters[nruns] = {0x0};

  for (int i = 0; i < nruns; i++)
  {
    hEvSel_filters2024[i] = new TH1F(Form("hEvSel_filters%d", i), ";;Events", numBins - 1, 0, numBins - 1);
    hEvSel_filters2024[i]->SetBinContent(1, hEvSelFinal2024[i]->GetBinContent(hEvSelFinal2024[i]->GetXaxis()->FindBin("Omegas") + iscale));
    hEvSel_filters2024[i]->SetBinContent(2, hEvSelFinal2024[i]->GetBinContent(hEvSelFinal2024[i]->GetXaxis()->FindBin("h-Omega") + iscale));
    hEvSel_filters2024[i]->SetBinContent(3, hEvSelFinal2024[i]->GetBinContent(hEvSelFinal2024[i]->GetXaxis()->FindBin("Xi-N") + iscale));
    hEvSel_filters2024[i]->SetBinContent(4, hEvSelFinal2024[i]->GetBinContent(hEvSelFinal2024[i]->GetXaxis()->FindBin("HighMult+Omega") + iscale));
    hEvSel_filters2024[i]->SetBinContent(5, hEvSelFinal2024[i]->GetBinContent(hEvSelFinal2024[i]->GetXaxis()->FindBin("Tracked Omega") + iscale));
    hEvSel_filters2024[i]->SetBinContent(6, hEvSelFinal2024[i]->GetBinContent(hEvSelFinal2024[i]->GetXaxis()->FindBin("Double Omega") + iscale));
    hEvSel_filters2024[i]->SetBinContent(7, hEvSelFinal2024[i]->GetBinContent(hEvSelFinal2024[i]->GetXaxis()->FindBin("Omega+Xi") + iscale));
    // hEvSel_skimmed_filters[i] = new TH1F(Form("hEvSel_skimmed_filters%d", i), ";;Events", numBins - 1, 0, numBins - 1);
    // hEvSel_skimmed_filters[i]->SetBinContent(1, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Omegas") + iscale));
    // hEvSel_skimmed_filters[i]->SetBinContent(2, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("h-Omega") + iscale));
    // hEvSel_skimmed_filters[i]->SetBinContent(3, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Xi-N") + iscale));
    // hEvSel_skimmed_filters[i]->SetBinContent(4, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("HighMult+Omega") + iscale));
    // hEvSel_skimmed_filters[i]->SetBinContent(5, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Tracked Omega") + iscale));
    // hEvSel_skimmed_filters[i]->SetBinContent(6, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Double Omega") + iscale));
    // hEvSel_skimmed_filters[i]->SetBinContent(7, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Omega+Xi") + iscale));
  }

  //------------ ratio skimmed/unskimmed

  TH1F *hRatio[nruns];
  for (int irun = 0; irun < nruns; irun++)
  {
    hRatio[irun] = new TH1F(Form("RatioSkim_UnSkim_run_%d", runnumber[irun]), Form("RatioSkim_UnSkim_run_%d", runnumber[irun]), hEvSel_skimmed_filters[irun]->GetNbinsX(), 0, hEvSel_skimmed_filters[irun]->GetNbinsX());
    for (int ibin = 0; ibin < hEvSel_skimmed_filters[irun]->GetNbinsX(); ibin++)
    {
      hRatio[irun]->SetBinContent(ibin + 1, hEvSel_skimmed_filters[irun]->GetBinContent(ibin + 1) / hEvSel_filters2025[irun]->GetBinContent(ibin + 1));
      cout << "run: " << runnumber[irun] << " |bin: " << ibin + 1 << " skimmed: " << hEvSel_skimmed_filters[irun]->GetBinContent(ibin + 1) << " unskimmed: " << hEvSel_filters2025[irun]->GetBinContent(ibin + 1) << " ratio: " << hRatio[irun]->GetBinContent(ibin + 1) << endl;
    }
    // hRatio[irun]->Sumw2();
  }

  TCanvas *cratio = new TCanvas("cratio", "cratio", 800, 800);
  cratio->cd();
  cratio->SetMargin(0.15, 0.1, 0.15, 0.1);
  TH1F *hratio = new TH1F(Form("hratio%d", 0), ";;Ratio Skimmed/Unskimmed", numBins, 0, numBins);
  hratio->SetTitle("LHC24am batch3");
  hratio->GetXaxis()->SetBinLabel(1, "OmegaDueToOtherFilters");
  hratio->GetXaxis()->SetBinLabel(2, "hOmega");
  hratio->GetXaxis()->SetBinLabel(3, "Xi-N");
  hratio->GetXaxis()->SetBinLabel(4, "HighMult+Omega");
  hratio->GetXaxis()->SetBinLabel(5, "Tracked Omega");
  hratio->GetXaxis()->SetBinLabel(6, "Double Omega");
  hratio->GetXaxis()->SetBinLabel(7, "Omega+Xi");
  hratio->GetXaxis()->SetBinLabel(8, "Lam+Lam");
  hratio->SetStats(0);
  hratio->GetYaxis()->SetRangeUser(0.01, 1.5);
  hratio->Draw("same");
  TLegend *legratio = new TLegend(0.4, 0.2, 0.88, 0.48); // 0.4, 0.6, 0.88, 0.88
  legratio->SetTextSize(0.02);
  legratio->SetTextFont(42);
  legratio->SetBorderSize(0);
  legratio->SetNColumns(5);

  for (int i = 0; i < nruns; i++)
  {
    color = nruns - i - 1 + FI;
    hRatio[i]->SetLineColor(color);
    hRatio[i]->SetMarkerColor(color);
    hRatio[i]->SetMarkerSize(1);
    hRatio[i]->SetMarkerStyle(8);
    hRatio[i]->Draw("PSAME");
    legratio->AddEntry(hRatio[i], Form("%d", runnumber[i]), "pl");
  }
  legratio->Draw("same");
  cratio->SaveAs("ratioskimmunskimm_LHC24am_batch3.png");

  //------------ #events for each filter
  TCanvas *cfull = new TCanvas("cfull", "cfull", 800, 800);
  cfull->cd();
  cfull->SetMargin(0.15, 0.15, 0.13, 0.08); // left,right,bottom,top
  cfull->SetLogy();

  TH1F *h = new TH1F(Form("h%d", 0), ";;Events", 17, -1, 14);
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
  h->GetXaxis()->SetBinLabel(14, "HighMult+Omega");
  h->GetXaxis()->SetBinLabel(15, "Double Omega");
  h->GetXaxis()->SetBinLabel(16, "Omega+Xi");
  h->GetXaxis()->SetBinLabel(17, "Lam+Lam");
  h->SetStats(0);
  h->GetYaxis()->SetRangeUser(0.5, 1.5 * TotEvtMax2025);
  h->Draw();
  TLegend *leg = new TLegend(0.4, 0.8, 0.85, 0.91);
  leg->SetTextSize(0.02);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);

  for (int i = 0; i < nruns; i++)
  {
    color = nruns - i - 1 + FI;
    for (Int_t j = 1; j <= hEvSelFinal2025[i]->GetNbinsX(); j++)
      hEvSelFinal2025[i]->GetXaxis()->SetBinLabel(j, h->GetXaxis()->GetBinLabel(j));

    hEvSelFinal2025[i]->GetYaxis()->SetRangeUser(0.5, 1.5 * TotEvtMax2025);
    hEvSelFinal2025[i]->SetLineColor(color);
    hEvSelFinal2025[i]->SetMarkerColor(color);
    hEvSelFinal2025[i]->SetMarkerSize(1);
    hEvSelFinal2025[i]->SetMarkerStyle(20);
    hEvSelFinal2025[i]->DrawCopy("SAME e");

    hEvSelFinal2024[i]->SetLineColor(kBlue);
    hEvSelFinal2024[i]->SetMarkerColor(kBlue);
    hEvSelFinal2024[i]->SetMarkerSize(1);
    hEvSelFinal2024[i]->SetMarkerStyle(25);
    hEvSelFinal2024[i]->DrawCopy("SAME e");

    leg->AddEntry(hEvSelFinal2025[i], "new task", "pl");
    leg->AddEntry(hEvSelFinal2024[i], "old task", "pl");
    leg->Draw("same");
  }
  cfull->SaveAs(Form("cfull_%s_total.png", period.c_str()));

  //------------ selectivity: #events for each filter/#total events (all processed events)
  TH1F *hSelectivity2025[8][nruns];
  Double_t tot_filters2025[8];
  Double_t tot_processed2025 = 0.;
  TH1F *hSelectivity_period2025[8];
  Float_t LowLimit[8] = {2e-6, 6e-6, 1.6e-6, 7e-6, 5e-6, 1.5e-7, 1.8e-6, 1.5e-7};
  Float_t UpLimit[8] = {3e-4, 1.5e-5, 2.2e-6, 1.6e-5, 7e-6, 3.5e-7, 2.2e-6, 3.5e-7};
  for (int ifilter = 0; ifilter < 8; ifilter++)
  {
    for (int irun = 0; irun < nruns; irun++)
    {
      hSelectivity2025[ifilter][irun] = new TH1F(Form("Selectivity2025_filter_%d_run_%d", ifilter, runnumber[irun]), Form("Selectivity2025_filter_%d_run_%d", ifilter, runnumber[irun]), nruns, 0, nruns);

      hSelectivity2025[ifilter][irun]->SetBinContent(irun + 1, hEvSel_filters2025[irun]->GetBinContent(ifilter + 1) / hEvSelFinal2025[irun]->GetBinContent(hEvSelFinal2025[irun]->GetXaxis()->FindBin("Processed events")));
      tot_filters2025[ifilter] += hEvSel_filters2025[irun]->GetBinContent(ifilter + 1);
      tot_processed2025 += hEvSelFinal2025[irun]->GetBinContent(hEvSelFinal2025[irun]->GetXaxis()->FindBin("Processed events"));

      cout << "run: " << irun << " |filter: " << ifilter << " |content: " << hSelectivity2025[ifilter][irun]->GetBinContent(irun + 1) << endl;
      // cout << " **** run: " << irun << " |filter: " << ifilter << " |eventi: " << hEvSel_filters[irun]->GetBinContent(ifilter+1) << " |tot: " << hEvSelFinal[irun]->GetBinContent(hEvSelFinal[irun]->GetXaxis()->FindBin("Processed events")) << " |sel: " << hSelectivity[ifilter][irun]->GetBinContent(irun+1) << endl;
      //  if(ifilter==0) cout << "run: " << irun << " |sel: " << hSelectivity[ifilter][irun]->GetBinContent(irun+1) << endl;
      //  if(ifilter==0) cout << "irun: " << irun << " |cont filt: " << hEvSel_filters[irun]->GetBinContent(ifilter+1) << endl;
    }

    hSelectivity_period2025[ifilter] = new TH1F(Form("Selectivity2025_filter_%d_period_%s", ifilter, period.c_str()), Form("Selectivity2025_filter_%d_period_%s", ifilter, period.c_str()), numBins - 1, 0, numBins - 1);
    hSelectivity_period2025[ifilter]->SetBinContent(ifilter + 1, tot_filters2025[ifilter] / tot_processed2025);

    // cout << "filter: " << ifilter << " sel tot: " << hSelectivity_period[ifilter]->GetBinContent(ifilter+1) << endl;
  }

  TH1F *hSelectivity2024[7][nruns];
  Double_t tot_filters2024[7];
  Double_t tot_processed2024 = 0.;
  TH1F *hSelectivity_period2024[7];
  for (int ifilter = 0; ifilter < 7; ifilter++)
  {
    for (int irun = 0; irun < nruns; irun++)
    {
      hSelectivity2024[ifilter][irun] = new TH1F(Form("Selectivity2024_filter_%d_run_%d", ifilter, runnumber[irun]), Form("Selectivity2024_filter_%d_run_%d", ifilter, runnumber[irun]), nruns, 0, nruns);

      hSelectivity2024[ifilter][irun]->SetBinContent(irun + 1, hEvSel_filters2024[irun]->GetBinContent(ifilter + 1) / hEvSelFinal2024[irun]->GetBinContent(hEvSelFinal2024[irun]->GetXaxis()->FindBin("Processed events")));
      tot_filters2024[ifilter] += hEvSel_filters2024[irun]->GetBinContent(ifilter + 1);
      tot_processed2024 += hEvSelFinal2024[irun]->GetBinContent(hEvSelFinal2024[irun]->GetXaxis()->FindBin("Processed events"));

      // cout<< "filter: " <<  ifilter  << " |content: "  <<  hSelectivity[ifilter][irun]->GetBinContent(irun+1) << endl;
      // cout << " **** run: " << irun << " |filter: " << ifilter << " |eventi: " << hEvSel_filters[irun]->GetBinContent(ifilter+1) << " |tot: " << hEvSelFinal[irun]->GetBinContent(hEvSelFinal[irun]->GetXaxis()->FindBin("Processed events")) << " |sel: " << hSelectivity[ifilter][irun]->GetBinContent(irun+1) << endl;
      // if(ifilter==0) cout << "run: " << irun << " |sel: " << hSelectivity[ifilter][irun]->GetBinContent(irun+1) << endl;
      // if(ifilter==0) cout << "irun: " << irun << " |cont filt: " << hEvSel_filters[irun]->GetBinContent(ifilter+1) << endl;
    }

    hSelectivity_period2024[ifilter] = new TH1F(Form("Selectivity2024_filter_%d_period_%s", ifilter, period.c_str()), Form("Selectivity2024_filter_%d_period_%s", ifilter, period.c_str()), numBins - 1, 0, numBins - 1);
    hSelectivity_period2024[ifilter]->SetBinContent(ifilter + 1, tot_filters2024[ifilter] / tot_processed2024);

    // cout << "filter: " << ifilter << " sel tot: " << hSelectivity_period[ifilter]->GetBinContent(ifilter+1) << endl;
  }

  //--- if we have only 1 run
  TH1F *hSel_allfilters2025 = new TH1F(Form("hSel2025_allfilters_run_%d", runnumber[0]), Form("hSel2025_allfilters_run_%d", runnumber[0]), 8, 0, 8);
  for (int ifilter = 0; ifilter < 8; ifilter++)
  {
    hSel_allfilters2025->SetBinContent(ifilter + 1, hSelectivity2025[ifilter][0]->GetBinContent(1));
    cout << "2025 filter:" << ifilter << " |content: " << hSel_allfilters2025->GetBinContent(ifilter + 1) << endl;
  }

  TH1F *hSel_allfilters2024 = new TH1F(Form("hSel2024_allfilters_run_%d", runnumber[0]), Form("hSel2024_allfilters_run_%d", runnumber[0]), 7, 0, 7);
  for (int ifilter = 0; ifilter < 7; ifilter++)
  {
    hSel_allfilters2024->SetBinContent(ifilter + 1, hSelectivity2024[ifilter][0]->GetBinContent(1));
    cout << "2024 filter:" << ifilter << " |content: " << hSel_allfilters2024->GetBinContent(ifilter + 1) << endl;
    // cout<<"filter:" << ifilter <<" |content: " << hSelectivity[ifilter][0]->GetBinContent(1) << " |content now: " << hSel_allfilters->GetBinContent(ifilter+1) << endl;
  }

  TString *xlabel[8];
  for (int ifilter = 0; ifilter < 8; ifilter++)
  {
    if (ifilter == 0)
      xlabel[ifilter] = new TString("Omega");
    if (ifilter == 1)
      xlabel[ifilter] = new TString("hOmega");
    if (ifilter == 2)
      xlabel[ifilter] = new TString("Xi-N");
    if (ifilter == 3)
      xlabel[ifilter] = new TString("HighMult+Omega");
    if (ifilter == 4)
      xlabel[ifilter] = new TString("Tracked Omega");
    if (ifilter == 5)
      xlabel[ifilter] = new TString("Double Omega");
    if (ifilter == 6)
      xlabel[ifilter] = new TString("Omega+Xi");
    if (ifilter == 7)
      xlabel[ifilter] = new TString("Lam+Lam");
  }

  TCanvas *cselallfilters = new TCanvas("cselallfilters", "cselallfilters", 1400, 1400);
  cselallfilters->cd();
  cselallfilters->SetMargin(0.15, 0.05, 0.13, 0.01); // left,right,bottom,top
  cselallfilters->SetLogy();
  TLegend *leg_sel = new TLegend(0.4, 0.9, 0.85, 0.95);
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

  hSel_allfilters2024->SetLineColor(kBlue);
  hSel_allfilters2024->SetMarkerColor(kBlue);
  hSel_allfilters2024->SetMarkerSize(3);
  hSel_allfilters2024->SetMarkerStyle(25);
  hSel_allfilters2024->DrawCopy("PSAME");

  leg_sel->AddEntry(hSel_allfilters2025, "new task", "pl");
  leg_sel->AddEntry(hSel_allfilters2024, "old task", "pl");
  leg_sel->Draw("same");

  cselallfilters->SaveAs(Form("cselallfilters_%s_total.png", period.c_str()));

  TH1F *hratioA = new TH1F(Form("hratioA_%d", runnumber[0]), Form("hratioA_%d", runnumber[0]), 7, 0, 7);
  for (int ifilter = 0; ifilter < 8; ifilter++)
  {
    hratioA->SetBinContent(ifilter + 1, hSel_allfilters2025->GetBinContent(ifilter + 1) / hSel_allfilters2024->GetBinContent(ifilter + 1));
  }

  TCanvas *cselallfilters_ratio = new TCanvas("cselallfilters_ratio", "cselallfilters_ratio", 1400, 1400);
  cselallfilters_ratio->cd();
  cselallfilters_ratio->SetMargin(0.15, 0.05, 0.13, 0.01); // left,right,bottom,top
  // cselallfilters_ratio->SetLogy();

  TH1D *frame_ratio = new TH1D("frame_ratio", "", 7, 0, 7);
  frame_ratio->GetYaxis()->SetTitle("Ratio New/Old");
  frame_ratio->SetStats(0);
  for (int ifilter = 1; ifilter < 8; ifilter++)
  {
    frame_ratio->GetXaxis()->SetBinLabel(ifilter, xlabel[ifilter - 1]->Data());
    frame_ratio->GetYaxis()->SetRangeUser(0.95, 1.05);
  }
  frame_ratio->Draw("");

  TLine *linelimit_ratio = new TLine(0, 1, 7, 1);
  linelimit_ratio->SetLineColor(kBlack);
  linelimit_ratio->SetLineWidth(2);
  linelimit_ratio->SetLineStyle(2);
  linelimit_ratio->Draw("same");
  hratioA->SetTitle("New/Old");
  hratioA->GetYaxis()->SetRangeUser(0.95, 1.05);
  hratioA->SetLineColor(kBlue);
  hratioA->SetMarkerColor(kBlue);
  hratioA->SetMarkerSize(3);
  hratioA->SetMarkerStyle(20);
  hratioA->DrawCopy("PSAME");

  cselallfilters_ratio->SaveAs(Form("cselallfilters_ratio_%s_total.png", period.c_str()));

  //------------ selectivity for each filter for each run
  TCanvas *cselectivity[8];
  TH1F *hselectivity[8];

  for (int ifilter = 0; ifilter < 8; ifilter++)
  {
    cselectivity[ifilter] = new TCanvas(Form("cselectivity_filter%d", ifilter), Form("cselectivity_filter%d", ifilter), 800, 800);
    cselectivity[ifilter]->cd();
    cselectivity[ifilter]->SetMargin(0.15, 0.01, 0.13, 0.08); // left,right,bottom,top
    // cselectivity[ifilter]->SetLogy();
    cselectivity[ifilter]->SetGridy();
    cselectivity[ifilter]->SetGridx();

    hselectivity[ifilter] = new TH1F(Form("hselectivity_filter%d", ifilter), ";;Trigger Selectivity", nruns, 0, nruns);
    if (ifilter == 0)
      hselectivity[ifilter]->SetTitle(Form("%s: Omega", period.c_str()));
    if (ifilter == 1)
      hselectivity[ifilter]->SetTitle(Form("%s: hOmega", period.c_str()));
    if (ifilter == 2)
      hselectivity[ifilter]->SetTitle(Form("%s: Xi-N", period.c_str()));
    if (ifilter == 3)
      hselectivity[ifilter]->SetTitle(Form("%s: HighMult+Omega", period.c_str()));
    if (ifilter == 4)
      hselectivity[ifilter]->SetTitle(Form("%s: Tracked Omega", period.c_str()));
    if (ifilter == 5)
      hselectivity[ifilter]->SetTitle(Form("%s: Double Omega", period.c_str()));
    if (ifilter == 6)
      hselectivity[ifilter]->SetTitle(Form("%s: Omega+Xi", period.c_str()));
    if (ifilter == 7)
      hselectivity[ifilter]->SetTitle(Form("%s: Lam+Lam", period.c_str()));

    for (int irun = 0; irun < nruns; irun++)
      hselectivity[ifilter]->GetXaxis()->SetBinLabel(irun + 1, Form("%d", runnumber[irun]));
    hselectivity[ifilter]->GetYaxis()->SetRangeUser(LowLimit[ifilter], UpLimit[ifilter]);
    hselectivity[ifilter]->Draw();

    for (int irun = 0; irun < nruns; irun++)
    {
      color = nruns - irun - 1 + FI;
      hSelectivity2025[ifilter][irun]->GetYaxis()->SetRangeUser(LowLimit[ifilter], UpLimit[ifilter]);
      hSelectivity2025[ifilter][irun]->SetLineColor(color);
      hSelectivity2025[ifilter][irun]->SetMarkerColor(color);
      hSelectivity2025[ifilter][irun]->SetMarkerSize(1.3);
      hSelectivity2025[ifilter][irun]->SetMarkerStyle(20);
      hSelectivity2025[ifilter][irun]->DrawCopy("PSAME");
    }

    if (ifilter == 0)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Omega.png", period.c_str()));
    if (ifilter == 1)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_hOmega.png", period.c_str()));
    if (ifilter == 2)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Xi-N.png", period.c_str()));
    if (ifilter == 3)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_HighMult+Omega.png", period.c_str()));
    if (ifilter == 4)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Tracked Omega.png", period.c_str()));
    if (ifilter == 5)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Double Omega.png", period.c_str()));
    if (ifilter == 6)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Omega+Xi.png", period.c_str()));
    if (ifilter == 7)
      cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_filter_Lam+Lam.png", period.c_str()));
  }

  // cout << hSelectivity_period[0]->GetBinContent(1) << endl;

  //------------ selectivity for each filter
  TCanvas *cselectivity_period = new TCanvas("cselectivity_period_total", "cselectivity_period_total", 800, 800);
  cselectivity_period->cd();
  cselectivity_period->SetMargin(0.15, 0.15, 0.13, 0.08); // left,right,bottom,top
  cselectivity_period->SetLogy();
  cselectivity_period->SetGridy();
  cselectivity_period->SetGridx();

  TH1F *hselectivity_period = new TH1F("hselectivity_period_total", ";;Trigger Selectivity", numBins - 1, 0, numBins - 1);
  hselectivity_period->SetTitle(Form("%s batch1_2", period.c_str()));
  hselectivity_period->GetXaxis()->SetBinLabel(1, "Omega");
  hselectivity_period->GetXaxis()->SetBinLabel(2, "hOmega");
  hselectivity_period->GetXaxis()->SetBinLabel(3, "Xi-N");
  hselectivity_period->GetXaxis()->SetBinLabel(4, "HighMult+Omega");
  hselectivity_period->GetXaxis()->SetBinLabel(5, "Tracked Omega");
  hselectivity_period->GetXaxis()->SetBinLabel(6, "Double Omega");
  hselectivity_period->GetXaxis()->SetBinLabel(7, "Omega+Xi");
  hselectivity_period->GetXaxis()->SetBinLabel(8, "Lam+Lam");
  hselectivity_period->SetStats(0);
  hselectivity_period->GetYaxis()->SetRangeUser(1e-9, 1e-3);
  hselectivity_period->Draw();

  /*
  for (int i = 0; i < 7; i++){
    //color = 5 - i - 1 + FI ;
    hSelectivity_period[i]->GetYaxis()->SetRangeUser(1e-9,1e-3);
    hSelectivity_period[i]->SetLineColor(kBlue);
    hSelectivity_period[i]->SetMarkerColor(kBlue);
    hSelectivity_period[i]->SetMarkerSize(1.3);
    hSelectivity_period[i]->SetMarkerStyle(20);
    hSelectivity_period[i]->DrawCopy("PSAME");
    //cout << "i: " << i << " cont: " << hSelectivity_period[i]->GetBinContent(i+1) << endl;
  }
    */
  cselectivity_period->SaveAs(Form("cselectivity_period_%s_batch4.png", period.c_str()));
}
