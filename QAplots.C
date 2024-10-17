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
Int_t numBins = 8; //Interesting trigger
Int_t iscale = 0; // 0 if input file has no mistakes. At some point axis labels where shifted by one bin and in those cases iscale = 1 should be set.

void QAplots(string period = "LHC24an"){

      
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  Long_t *dummy1 = 0, *dummy2 = 0, *dummy3 = 0, *dummy4 = 0;
    
  TStopwatch time;
  time.Start();
  
  char pass_name[2][30] = {"apass1","apass1_skimmed"};
  const int nruns = 23; 
  int runnumber[nruns] = {556152,556160,556164,556182,556210,556218,556237,556248,556269,556284,556370,556372,556412,556437,556454,556461,556482,556485,556491,556497,556517,556542,556562}; //batch1,2

  //------------ read files
  TFile *file_in[nruns] = {0x0};
  TFile *file_in_skimmed[nruns] = {0x0};
  
  for(int irun = 0; irun<nruns; irun ++){
   if(gSystem->GetPathInfo(Form("../TriggerForRun3/EventFiltering2024/LHC24an/AnalysisResults_fullrun_%s_%s_%d.root",period.c_str(),pass_name[0],runnumber[irun]),dummy1,dummy2,dummy3,dummy4) != 0) cout << "File Not Found! Try again" << endl;    
   file_in[irun] = new TFile(Form("../TriggerForRun3/EventFiltering2024/LHC24an/AnalysisResults_fullrun_%s_%s_%d.root",period.c_str(),pass_name[0],runnumber[irun]),"read");
   printf("Open File: %s\n",file_in[irun]->GetName());
   
   if(gSystem->GetPathInfo(Form("../TriggerForRun3/EventFiltering2024/LHC24an/AnalysisResults_fullrun_%s_%s_%d.root",period.c_str(),pass_name[1],runnumber[irun]),dummy1,dummy2,dummy3,dummy4) != 0) cout << "File Not Found! Try again" << endl;    
   file_in_skimmed[irun] = new TFile(Form("../TriggerForRun3/EventFiltering2024/LHC24an/AnalysisResults_fullrun_%s_%s_%d.root",period.c_str(),pass_name[1],runnumber[irun]),"read");
   printf("Open File: %s\n",file_in_skimmed[irun]->GetName());
   cout << "irun:" << irun << endl;
  }
   
   //color palette
   const int NRGBs = 5;
   double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 }; 
   double red[NRGBs]   = { 0.00, 0.00, 0.78, 1.00, 0.51 };
   double green[NRGBs] = { 0.00, 0.81, 0.90, 0.20, 0.00 };
   double blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   int FI = TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,nruns);
   int color;
   
   //------------ #events for each filter
   TH1F *hEvSel[nruns] = {0x0};
   TH1F *hEvSelFinal[nruns] = {0x0};
   Double_t TotEvtMax = 0;
   
   for (int i = 0; i < nruns; i++){
    hEvSel[i] = (TH1F *)file_in[i]->Get("lf-strangeness-filter/hProcessedEvents");
    hEvSel[i]->SetName(Form("hEvSel_%d", i));
    if (hEvSel[i]->GetBinContent(1) > TotEvtMax) TotEvtMax = hEvSel[i]->GetBinContent(1);
        
   }
 
   for (int i = 0; i < nruns; i++){
    color = nruns - i - 1 + FI ; 
    hEvSelFinal[i] = (TH1F *)hEvSel[0]->Clone(Form("hEvSelFinal_%d", i));
    for (int j = 1; j <= hEvSelFinal[i]->GetNbinsX(); j++)
    {
        hEvSelFinal[i]->SetBinContent(j, hEvSel[i]->GetBinContent(j));
        hEvSelFinal[i]->SetBinError(j, TMath::Sqrt(hEvSel[i]->GetBinContent(j)));
    }
    hEvSelFinal[i]->SetLineColor(color);
   }
   
   //------------ #events for each filter from skimmed data
   TH1F *hEvSel_skimmed[nruns] = {0x0};
   TH1F *hEvSelFinal_skimmed[nruns] = {0x0};
   
   for (int i = 0; i < nruns; i++){
    hEvSel_skimmed[i] = (TH1F *)file_in_skimmed[i]->Get("lf-strangeness-filter/hProcessedEvents");
    hEvSel_skimmed[i]->SetName(Form("hEvSel_skimmed_%d", i));
   }
   
   for (int i = 0; i < nruns; i++){
    color = nruns - i - 1 + FI ; 
    hEvSelFinal_skimmed[i] = (TH1F *)hEvSel_skimmed[0]->Clone(Form("hEvSelFinal_skimmed_%d", i));
    for (int j = 1; j <= hEvSelFinal_skimmed[i]->GetNbinsX(); j++)
    {
        hEvSelFinal_skimmed[i]->SetBinContent(j, hEvSel_skimmed[i]->GetBinContent(j));
        hEvSelFinal_skimmed[i]->SetBinError(j, TMath::Sqrt(hEvSel_skimmed[i]->GetBinContent(j)));
        //if(runnumber[i]==551759) cout << "skimmed run: " << runnumber[i] << " |bin: " << j << " val: " << hEvSelFinal_skimmed[i]->GetBinContent(j) << endl;
    }
    hEvSelFinal_skimmed[i]->SetLineColor(color);
   }
   
    
   // name of each filter/column
   char binlabel[17][30] = {"All events","Processed events","Events w/ high-pT hadron", "Omegas", "h-Omega", "2Xi", "3Xi", "4Xi", "Xi-N", "Omega large R", "Xi", "Tracked Xi", "Tracked Omega", "HighMult+Omega", "Double Omega", "Omega+Xi"};
   
   for (int i = 0; i < nruns; i++){
    for (Int_t j = 1; j <= hEvSelFinal[i]->GetNbinsX(); j++){
      hEvSelFinal[i]->GetXaxis()->SetBinLabel(j,binlabel[j-1]);
    }
    
    for (Int_t j = 1; j <= hEvSelFinal_skimmed[i]->GetNbinsX(); j++){
      hEvSelFinal_skimmed[i]->GetXaxis()->SetBinLabel(j,binlabel[j-1]);
    }
   }
   
   
   TH1F *hEvSel_filters[nruns] = {0x0};
   TH1F *hEvSel_skimmed_filters[nruns] = {0x0};
   
   for (int i = 0; i < nruns; i++){
     hEvSel_filters[i] = new TH1F(Form("hEvSel_filters%d", i), ";;Events", numBins - 1, 0, numBins - 1);
     hEvSel_filters[i]->SetBinContent(1, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Omegas") + iscale));
     hEvSel_filters[i]->SetBinContent(2, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("h-Omega") + iscale));
     hEvSel_filters[i]->SetBinContent(3, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Xi-N") + iscale));
     hEvSel_filters[i]->SetBinContent(4, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("HighMult+Omega") + iscale));
     hEvSel_filters[i]->SetBinContent(5, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Tracked Omega") + iscale));
     hEvSel_filters[i]->SetBinContent(6, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Double Omega") + iscale));
     hEvSel_filters[i]->SetBinContent(7, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Omega+Xi") + iscale));
     hEvSel_skimmed_filters[i] = new TH1F(Form("hEvSel_skimmed_filters%d", i), ";;Events", numBins - 1, 0, numBins - 1);
     hEvSel_skimmed_filters[i]->SetBinContent(1, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Omegas") + iscale));
     hEvSel_skimmed_filters[i]->SetBinContent(2, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("h-Omega") + iscale));
     hEvSel_skimmed_filters[i]->SetBinContent(3, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Xi-N") + iscale));
     hEvSel_skimmed_filters[i]->SetBinContent(4, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("HighMult+Omega") + iscale));
     hEvSel_skimmed_filters[i]->SetBinContent(5, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Tracked Omega") + iscale));
     hEvSel_skimmed_filters[i]->SetBinContent(6, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Double Omega") + iscale));
     hEvSel_skimmed_filters[i]->SetBinContent(7, hEvSelFinal_skimmed[i]->GetBinContent(hEvSelFinal_skimmed[i]->GetXaxis()->FindBin("Omega+Xi") + iscale));
   }
   
   
    //------------ ratio skimmed/unskimmed
    
    TH1F *hRatio[nruns];
    for(int irun = 0; irun<nruns; irun++){
     hRatio[irun] = new TH1F(Form("RatioSkim_UnSkim_run_%d",runnumber[irun]),Form("RatioSkim_UnSkim_run_%d",runnumber[irun]),hEvSel_skimmed_filters[irun]->GetNbinsX(),0,hEvSel_skimmed_filters[irun]->GetNbinsX());   
     for(int ibin=0;ibin<hEvSel_skimmed_filters[irun]->GetNbinsX();ibin++)  hRatio[irun]->SetBinContent(ibin+1,hEvSel_skimmed_filters[irun]->GetBinContent(ibin+1)/hEvSel_filters[irun]->GetBinContent(ibin+1));  
     //hRatio[irun]->Sumw2();
    }
     

    TCanvas *cratio = new TCanvas("cratio","cratio",800,800);
    cratio->cd();
    cratio->SetMargin(0.15,0.1,0.15,0.1);
    TH1F *hratio = new TH1F(Form("hratio%d", 0), ";;Ratio Skimmed/Unskimmed", numBins - 1, 0, numBins - 1);
    hratio->SetTitle("LHC24an");
    hratio->GetXaxis()->SetBinLabel(1, "OmegaDueToOtherFilters");
    hratio->GetXaxis()->SetBinLabel(2, "hOmega");
    hratio->GetXaxis()->SetBinLabel(3, "Xi-N");
    hratio->GetXaxis()->SetBinLabel(4, "HighMult+Omega");
    hratio->GetXaxis()->SetBinLabel(5, "Tracked Omega");
    hratio->GetXaxis()->SetBinLabel(6, "Double Omega");
    hratio->GetXaxis()->SetBinLabel(7, "Omega+Xi");
    hratio->SetStats(0);
    //hratio->GetYaxis()->SetRangeUser(0.01, 1.5);
    hratio->GetYaxis()->SetRangeUser(0.8, 1.2);
    hratio->Draw("same"); 
    TLegend *legratio = new TLegend(0.4, 0.6, 0.88, 0.88);//0.4, 0.6, 0.88, 0.88
    legratio->SetTextSize(0.02);
    legratio->SetTextFont(42);
    legratio->SetBorderSize(0);
    legratio->SetNColumns(5);
    
    for (int i = 0; i < nruns; i++){
        color = nruns - i - 1 + FI ;
        hRatio[i]->SetLineColor(color);
        hRatio[i]->SetMarkerColor(color);
        hRatio[i]->SetMarkerSize(1);
        hRatio[i]->SetMarkerStyle(8);
        hRatio[i]->Draw("PSAME");
        legratio->AddEntry(hRatio[i],Form("%d",runnumber[i]),"pl");
    }
    legratio->Draw("same");
    cratio->SaveAs("ratioskimmunskimm_LHC24an_batch1_2.png");
    
    //------------ #events for each filter
    TCanvas *cfull = new TCanvas("cfull","cfull",800,800);
    cfull->cd();
    cfull->SetMargin(0.15,0.15,0.13,0.08); //left,right,bottom,top
    cfull->SetLogy();
    
    TH1F *h = new TH1F(Form("h%d", 0), ";;Events", 16, -1, 14);
    h->SetTitle(Form("%s",period.c_str()));
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
    h->SetStats(0);
    h->GetYaxis()->SetRangeUser(0.5, 1.5 * TotEvtMax);
    h->Draw();
    TLegend *leg = new TLegend(0.4, 0.8, 0.85, 0.91);
    leg->SetTextSize(0.02);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetNColumns(5);

    for (int i = 0; i < nruns; i++){
      color = nruns - i - 1 + FI ;  
      for (Int_t j = 1; j <= hEvSelFinal[i]->GetNbinsX(); j++)  hEvSelFinal[i]->GetXaxis()->SetBinLabel(j, h->GetXaxis()->GetBinLabel(j));
      
      hEvSelFinal[i]->GetYaxis()->SetRangeUser(0.5, 1.5 * TotEvtMax);
      hEvSelFinal[i]->SetLineColor(color);
      hEvSelFinal[i]->SetMarkerColor(color);
      hEvSelFinal[i]->SetMarkerSize(1);
      hEvSelFinal[i]->SetMarkerStyle(20);
      hEvSelFinal[i]->DrawCopy("SAME e");
      leg->AddEntry(hEvSelFinal[i],Form("%d",runnumber[i]),"pl");
      leg->Draw("same");
    }
    cfull->SaveAs(Form("cfull_%s_total.png",period.c_str()));
    
    
    //------------ selectivity: #events for each filter/#total events (all processed events) 
    TH1F *hSelectivity[7][nruns];
    Double_t tot_filters[7];
    Double_t tot_processed = 0.;
    TH1F *hSelectivity_period[7];
    Float_t LowLimit[7] = {2e-6, 6e-6, 1.3e-6, 7e-6, 5e-6, 1.5e-7, 1e-6};
    Float_t UpLimit[7] = {3e-4, 1e-5, 1.8e-6, 1.6e-5, 1e-5, 3.5e-7, 2e-6};
    for(int ifilter = 0; ifilter<7;ifilter++){
     for(int irun = 0; irun<nruns; irun++){
      hSelectivity[ifilter][irun] = new TH1F(Form("Selectivity_filter_%d_run_%d",ifilter,runnumber[irun]),Form("Selectivity_filter_%d_run_%d",ifilter,runnumber[irun]),nruns,0,nruns);   
     
      hSelectivity[ifilter][irun]->SetBinContent(irun+1,hEvSel_filters[irun]->GetBinContent(ifilter+1)/hEvSelFinal[irun]->GetBinContent(hEvSelFinal[irun]->GetXaxis()->FindBin("Processed events")) ); 
      tot_filters[ifilter]+=hEvSel_filters[irun]->GetBinContent(ifilter+1);
      tot_processed+= hEvSelFinal[irun]->GetBinContent(hEvSelFinal[irun]->GetXaxis()->FindBin("Processed events"));
      
      //cout << " **** run: " << irun << " |filter: " << ifilter << " |eventi: " << hEvSel_filters[irun]->GetBinContent(ifilter+1) << " |tot: " << hEvSelFinal[irun]->GetBinContent(hEvSelFinal[irun]->GetXaxis()->FindBin("Processed events")) << " |sel: " << hSelectivity[ifilter][irun]->GetBinContent(irun+1) << endl;   
      //if(ifilter==0) cout << "run: " << irun << " |sel: " << hSelectivity[ifilter][irun]->GetBinContent(irun+1) << endl; 
      //if(ifilter==0) cout << "irun: " << irun << " |cont filt: " << hEvSel_filters[irun]->GetBinContent(ifilter+1) << endl;
         
     }
     
     hSelectivity_period[ifilter] = new TH1F(Form("Selectivity_filter_%d_period_%s",ifilter,period.c_str()),Form("Selectivity_filter_%d_period_%s",ifilter,period.c_str()),numBins-1,0,numBins-1);
     hSelectivity_period[ifilter]->SetBinContent(ifilter+1,tot_filters[ifilter]/tot_processed);

     //cout << "filter: " << ifilter << " sel tot: " << hSelectivity_period[ifilter]->GetBinContent(ifilter+1) << endl;
     
    }
                                   
    //------------ selectivity for each filter for each run
    TCanvas *cselectivity[7];
    TH1F *hselectivity[7];
    
    for(int ifilter = 0; ifilter<7;ifilter++){
     cselectivity[ifilter] = new TCanvas(Form("cselectivity_filter%d",ifilter),Form("cselectivity_filter%d",ifilter),1400,800);
     cselectivity[ifilter]->cd();
     cselectivity[ifilter]->SetMargin(0.1,0.01,0.13,0.08); //left,right,bottom,top
     //cselectivity[ifilter]->SetLogy();
     cselectivity[ifilter]->SetGridy();
     cselectivity[ifilter]->SetGridx();
    
     hselectivity[ifilter] = new TH1F(Form("hselectivity_filter%d",ifilter), ";;Trigger Selectivity", nruns, 0, nruns);
     if(ifilter==0) hselectivity[ifilter]->SetTitle(Form("%s: Omega",period.c_str()));
     if(ifilter==1) hselectivity[ifilter]->SetTitle(Form("%s: hOmega",period.c_str()));
     if(ifilter==2) hselectivity[ifilter]->SetTitle(Form("%s: Xi-N",period.c_str()));
     if(ifilter==3) hselectivity[ifilter]->SetTitle(Form("%s: HighMult+Omega",period.c_str()));
     if(ifilter==4) hselectivity[ifilter]->SetTitle(Form("%s: Tracked Omega",period.c_str()));
     if(ifilter==5) hselectivity[ifilter]->SetTitle(Form("%s: Double Omega",period.c_str()));
     if(ifilter==6) hselectivity[ifilter]->SetTitle(Form("%s: Omega+Xi",period.c_str()));
     
     for(int irun = 0; irun<nruns; irun++) hselectivity[ifilter]->GetXaxis()->SetBinLabel(irun+1,Form("%d",runnumber[irun]));
     hselectivity[ifilter]->GetYaxis()->SetRangeUser(LowLimit[ifilter],UpLimit[ifilter]);
     hselectivity[ifilter]->Draw();
     
     for(int irun = 0; irun<nruns; irun++){
      color = nruns - irun - 1 + FI ;  
      hSelectivity[ifilter][irun]->GetYaxis()->SetRangeUser(LowLimit[ifilter],UpLimit[ifilter]);
      hSelectivity[ifilter][irun]->SetLineColor(color);
      hSelectivity[ifilter][irun]->SetMarkerColor(color);
      hSelectivity[ifilter][irun]->SetMarkerSize(1.3);
      hSelectivity[ifilter][irun]->SetMarkerStyle(20);
      hSelectivity[ifilter][irun]->DrawCopy("PSAME");
     }
     if(ifilter==0) cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_batch1_2_filter_Omega.png",period.c_str()));
     if(ifilter==1) cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_batch1_2_filter_hOmega.png",period.c_str()));
     if(ifilter==2) cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_batch1_2_filter_Xi-N.png",period.c_str()));
     if(ifilter==3) cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_batch1_2_filter_HighMult+Omega.png",period.c_str()));
     if(ifilter==4) cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_batch1_2_filter_Tracked Omega.png",period.c_str()));
     if(ifilter==5) cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_batch1_2_filter_Double Omega.png",period.c_str()));
     if(ifilter==6) cselectivity[ifilter]->SaveAs(Form("cselectivity_period_%s_batch1_2_filter_Omega+Xi.png",period.c_str()));
    }
    
    //cout << hSelectivity_period[0]->GetBinContent(1) << endl;
    
    //------------ selectivity for each filter
    TCanvas *cselectivity_period = new TCanvas("cselectivity_period_total","cselectivity_period_total",800,800);
    cselectivity_period->cd();
    cselectivity_period->SetMargin(0.15,0.15,0.13,0.08); //left,right,bottom,top
    cselectivity_period->SetLogy();
    cselectivity_period->SetGridy();
    cselectivity_period->SetGridx();
    
    TH1F *hselectivity_period = new TH1F("hselectivity_period_total", ";;Trigger Selectivity", numBins-1, 0, numBins - 1);
    hselectivity_period->SetTitle(Form("%s batch1_2",period.c_str()));
    hselectivity_period->GetXaxis()->SetBinLabel(1, "Omega");
    hselectivity_period->GetXaxis()->SetBinLabel(2, "hOmega");
    hselectivity_period->GetXaxis()->SetBinLabel(3, "Xi-N");
    hselectivity_period->GetXaxis()->SetBinLabel(4, "HighMult+Omega");
    hselectivity_period->GetXaxis()->SetBinLabel(5, "Tracked Omega");
    hselectivity_period->GetXaxis()->SetBinLabel(6, "Double Omega");
    hselectivity_period->GetXaxis()->SetBinLabel(7, "Omega+Xi");
    hselectivity_period->SetStats(0);
    hselectivity_period->GetYaxis()->SetRangeUser(1e-9,1e-3);
    hselectivity_period->Draw();
   
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
    //cselectivity_period->SaveAs(Form("cselectivity_period_%s_batch1_2.png",period.c_str()));
    
    
    
    
    
    
}
