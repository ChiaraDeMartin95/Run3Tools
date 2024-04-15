#include <iostream>
#include <fstream>

const Float_t StrLimit = 4.3E-5;
const Float_t LFLimit = 5E-5;
Int_t numBins = 6;
Int_t iscale = 0; //0 if input file has no mistakes. At some point axis labels where shifted by one bin and in those cases iscale = 1 should be set.

void drawfiltersQA(TString filename = "list2023NN.txt", Int_t ChosenPeriod = 2)
{

    std::vector<std::string> name;
    std::ifstream file(Form("%s", filename.Data()));
    // std::string remove = "/Users/mbp-cdm-01/Desktop/dottorato/1Anno/QAStrangeness/TriggerForRun3/EventFiltering2023/AnalysisResults_merged_";
    // std::string remove = "/Users/mbp-cdm-01/Desktop/dottorato/1Anno/QAStrangeness/TriggerForRun3/EventFiltering2023/AnalysisResults/AnalysisResults_";
    std::string remove = "/Users/mbp-cdm-01/Desktop/ResearchTeachingActivity/QAStrangeness/TriggerForRun3/EventFiltering2023/LHC23v_pass4/AnalysisResults_LHC23v_pass4_CEFP_";
    std::string remove2 = ".root";

    cout << filename.Data() << endl;
    if (file.is_open())
    {
        std::string line;

        while (std::getline(file, line))
        {
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
            name.push_back(line);
            cout << line << endl;
        }
        file.close();
    }
    else
    {
        std::cerr << "Unable to open file: " << std::endl;
    }

    Int_t colors[14] = {634, 628, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
    Int_t MarkerMult[14] = {20, 21, 33, 34, 29, 20, 21, 33, 34, 29, 20, 21, 33, 34};
    const int nfiles = name.size();
    cout << "Number of files = " << nfiles << endl;

    TFile *f[nfiles];
    for (int i = 0; i < nfiles; i++)
    {
        cout << "name " << name[i].c_str() << endl;
        // f[i] = TFile::Open(Form("%s%s_round2.root", remove.c_str(), name[i].c_str()));
        f[i] = TFile::Open(Form("%s%s.root", remove.c_str(), name[i].c_str()));
    }

    TH1F *hEvSel[nfiles];
    TH1F *hEvSelFinal[nfiles];
    Float_t RelError = 0;
    Int_t TotEvt = 0;
    for (int i = 0; i < nfiles; i++)
    {
        RelError = 0;
        hEvSel[i] = (TH1F *)f[i]->Get("lf-strangeness-filter/hProcessedEvents");
        TotEvt = hEvSel[i]->GetBinContent(1);
        cout << "Total events " << name[i].c_str() << " = " << TotEvt << endl;
    }
    for (int i = 0; i < nfiles; i++)
    {
        hEvSelFinal[i] = (TH1F *)hEvSel[0]->Clone(Form("hEvSelFinal_%d", i));
        for (int j = 1; j <= hEvSelFinal[i]->GetNbinsX(); j++)
        {
            RelError = sqrt(1. / hEvSel[i]->GetBinContent(j) + 1. / hEvSel[i]->GetBinContent(1));
            cout << "i " << i << " j " << j << " RelError " << RelError << " hEvSel[i]->GetBinContent(j) " << hEvSel[i]->GetBinContent(j) << endl;
            hEvSelFinal[i]->SetBinContent(j, hEvSel[i]->GetBinContent(j) / hEvSel[i]->GetBinContent(1));
            hEvSelFinal[i]->SetBinError(j, hEvSelFinal[i]->GetBinContent(j) * RelError);
            if (j == 1)
                hEvSelFinal[i]->SetBinError(j, 0);
        }
        hEvSelFinal[i]->SetLineColor(colors[i]);
    }

    TCanvas *cfull = new TCanvas("cfull", "cfull", 900, 1000);
    TPad *pad1f = new TPad("pad1f", "pad1f", 0, 0.3, 1, 1);
    pad1f->SetBottomMargin(0.01);
    pad1f->SetLeftMargin(0.15);
    pad1f->SetRightMargin(0.05);
    pad1f->SetTopMargin(0.05);
    pad1f->Draw();
    pad1f->cd();
    pad1f->SetLogy();

    TH1F *h1 = new TH1F(Form("h1%d", 0), ";;Selectivity", 13, -1, 12);
    h1->GetXaxis()->SetBinLabel(1, "All events");
    h1->GetXaxis()->SetBinLabel(2, "Processed events");
    h1->GetXaxis()->SetBinLabel(3, "Events w/ high-pT hadron");
    h1->GetXaxis()->SetBinLabel(4, "Omegas");
    h1->GetXaxis()->SetBinLabel(5, "high-pT + Xi");
    h1->GetXaxis()->SetBinLabel(6, "2Xi");
    h1->GetXaxis()->SetBinLabel(7, "3Xi");
    h1->GetXaxis()->SetBinLabel(8, "4Xi");
    h1->GetXaxis()->SetBinLabel(9, "Xi-N");
    h1->GetXaxis()->SetBinLabel(10, "Omega large R");
    h1->GetXaxis()->SetBinLabel(11, "Xi");
    h1->GetXaxis()->SetBinLabel(12, "Tracked Xi");
    h1->GetXaxis()->SetBinLabel(13, "Tracked Omega");
    h1->SetStats(0);
    h1->GetYaxis()->SetRangeUser(1E-10, 1E-1);
    h1->Draw();
    for (int i = 0; i < nfiles; i++)
    {
        for (Int_t j = 1; j <= hEvSelFinal[i]->GetNbinsX(); j++)
        {
            hEvSelFinal[i]->GetXaxis()->SetBinLabel(j, h1->GetXaxis()->GetBinLabel(j));
        }
        hEvSelFinal[i]->GetYaxis()->SetRangeUser(1E-10, 1E-1);
        hEvSelFinal[i]->DrawCopy("SAME e");
    }
    //hEvSelFinal[0]->Draw("SAME e");

    cfull->cd();
    TPad *pad2f = new TPad("pad2f", "pad2f", 0, 0, 1, 0.3);
    pad2f->SetTopMargin(0.01);
    pad2f->SetBottomMargin(0.3);
    pad2f->SetLeftMargin(0.15);
    pad2f->SetRightMargin(0.05);
    pad2f->Draw();
    pad2f->cd();

    TH1F *r1 = (TH1F *)h1->Clone("r1");
    r1->GetYaxis()->SetTitle("Ratio");
    r1->GetYaxis()->SetRangeUser(0., 2.);
    r1->GetXaxis()->SetLabelSize(0.1);
    r1->GetYaxis()->SetLabelSize(0.05);
    r1->GetYaxis()->SetTitleSize(0.2);
    r1->Draw();

    TF1 *fr1 = new TF1("fr1", "1", 0, 12);
    fr1->SetLineColor(kBlack);
    fr1->SetLineStyle(10);
    fr1->FixParameter(0, 1);
    fr1->Draw("SAME");

    TH1F *hRatiof[nfiles];
    for (int i = 0; i < nfiles; i++)
    {
        hRatiof[i] = (TH1F *)hEvSelFinal[i]->Clone(Form("hRatiof_%i", i));
        hRatiof[i]->Sumw2();
        hRatiof[i]->Divide(hEvSelFinal[ChosenPeriod]);
        hRatiof[i]->SetLineColor(colors[i]);
        hRatiof[i]->SetMarkerColor(colors[i]);
        hRatiof[i]->SetMarkerStyle(MarkerMult[i]);
        if (i != ChosenPeriod)
            hRatiof[i]->Draw("SAME pl");
    }

    TH1F *hEvSel23[nfiles];
    for (int i = 0; i < nfiles; i++)
    {
        hEvSel23[i] = new TH1F(Form("hEvSel23_%d", i), ";;Selectivity", numBins, 0, numBins);
        hEvSel23[i]->SetLineColor(colors[i]);
        hEvSel23[i]->SetLineWidth(2);

        hEvSel23[i]->SetBinContent(1, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Omegas")+ iscale) * 0.2);
        hEvSel23[i]->SetBinContent(2, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("3Xi")+ iscale));
        hEvSel23[i]->SetBinContent(3, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Xi-N")+ iscale));
        hEvSel23[i]->SetBinContent(4, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Omega large R")+ iscale));
        hEvSel23[i]->SetBinContent(5, hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Tracked Omega")+ iscale));
        double tot = hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Omegas")+ iscale) * 0.2 +
                     hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("3Xi")+ iscale) +
                     hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Xi-N")+ iscale) +
                     hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Omega large R")+ iscale) +
                     hEvSelFinal[i]->GetBinContent(hEvSelFinal[i]->GetXaxis()->FindBin("Tracked Omega")+ iscale);
        hEvSel23[i]->SetBinContent(numBins, tot);

        hEvSel23[i]->SetBinError(1, hEvSelFinal[i]->GetBinError(hEvSelFinal[i]->GetXaxis()->FindBin("Omegas")+ iscale) * 0.2);
        hEvSel23[i]->SetBinError(2, hEvSelFinal[i]->GetBinError(hEvSelFinal[i]->GetXaxis()->FindBin("3Xi")+ iscale));
        hEvSel23[i]->SetBinError(3, hEvSelFinal[i]->GetBinError(hEvSelFinal[i]->GetXaxis()->FindBin("Xi-N")+ iscale));
        hEvSel23[i]->SetBinError(4, hEvSelFinal[i]->GetBinError(hEvSelFinal[i]->GetXaxis()->FindBin("Omega large R")+ iscale));
        hEvSel23[i]->SetBinError(5, hEvSelFinal[i]->GetBinError(hEvSelFinal[i]->GetXaxis()->FindBin("Tracked Omega")+ iscale));
        hEvSel23[i]->SetBinError(numBins, 0);

        cout << "Total rejection factor " << name[i].c_str() << " = " << tot << endl;
    }

    TCanvas *c = new TCanvas("c", "c", 900, 1000);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->SetBottomMargin(0.01);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.05);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();

    TH1F *h = new TH1F(Form("h%d", 0), ";;Selectivity", numBins, 0, numBins);
    h->GetXaxis()->SetBinLabel(1, "Omega #times 0.2");
    h->GetXaxis()->SetBinLabel(2, "3Xi");
    h->GetXaxis()->SetBinLabel(3, "Xi-N");
    h->GetXaxis()->SetBinLabel(4, "Omega large R");
    h->GetXaxis()->SetBinLabel(5, "Tracked Omega");
    h->GetXaxis()->SetBinLabel(numBins, "Total");
    h->SetStats(0);
    h->GetYaxis()->SetRangeUser(1E-10, 1E-1);
    h->Draw("same");
    for (int i = 0; i < nfiles; i++)
    {
        hEvSel23[i]->Draw("SAME e");
        for (Int_t j = 1; j <= numBins; j++)
        {
            cout << "Bin content : " << h->GetXaxis()->GetBinLabel(j) << " " << hEvSel23[i]->GetBinContent(j) << endl;
        }
    }
    // hEvSel23[0]->Draw("SAME e");

    TLine *l = new TLine(0, StrLimit, numBins, StrLimit);
    l->SetLineStyle(7);
    l->SetLineColor(kAzure + 10);
    l->SetLineWidth(3);
    l->Draw();
    TLine *l1 = new TLine(0, LFLimit, numBins, LFLimit);
    l1->SetLineStyle(7);
    l1->SetLineColor(kOrange + 7);
    l1->SetLineWidth(3);
    l1->Draw();

    TLegend *lg = new TLegend(0.6, 0.6, 0.9, 0.95);
    for (int i = 0; i < nfiles; i++)
    {
        lg->AddEntry(hEvSel23[i], name[i].c_str(), "l");
    }
    lg->SetBorderSize(0);
    lg->AddEntry(l, "strangeness limit", "l");
    lg->AddEntry(l1, "LF limit", "l");
    lg->SetTextSize(0.03);
    lg->Draw();

    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();

    TH1F *hRatio = (TH1F *)h->Clone("hRatio");
    hRatio->GetYaxis()->SetTitle("Ratio");
    hRatio->GetXaxis()->SetLabelSize(0.1);
    hRatio->GetYaxis()->SetLabelSize(0.05);
    hRatio->GetYaxis()->SetTitleSize(0.1);
    hRatio->GetYaxis()->SetTitleOffset(0.5);
    hRatio->SetStats(0);
    hRatio->GetYaxis()->SetRangeUser(0., 2.);
    hRatio->Draw("");

    TH1F *hRatioToLimit[nfiles];
    for (int i = 0; i < nfiles; i++)
    {
        hRatioToLimit[i] = (TH1F *)hEvSel23[i]->Clone(Form("hRatioToLimit%d", i));
        hRatioToLimit[i]->Sumw2();
        hRatioToLimit[i]->Divide(hEvSel23[ChosenPeriod]);
        for (Int_t l = 1; l <= hRatioToLimit[i]->GetNbinsX(); l++)
        {
            cout << hRatioToLimit[i]->GetBinContent(l) << endl;
        }
        hRatioToLimit[i]->SetLineColor(colors[i]);
        hRatioToLimit[i]->SetMarkerColor(colors[i]);
        hRatioToLimit[i]->SetMarkerStyle(MarkerMult[i]);

        // hRatioToLimit[i]->Scale(1. / StrLimit);
        hRatioToLimit[i]->Draw("SAME pl");
    }

    TF1 *f1 = new TF1("f1", "1", 0, 5);
    f1->SetLineColor(kBlack);
    f1->SetLineStyle(10);
    f1->FixParameter(0, 1);
    f1->Draw("SAME");

    c->SaveAs("images/rejfactors.png");
}