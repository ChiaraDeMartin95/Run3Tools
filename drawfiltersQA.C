#include <iostream>
#include <fstream>

const Float_t StrLimit = 4.3E-5;

void drawfiltersQA(TString filename = "listGoodvsBad.txt", Int_t ChosenPeriod = 1)
{

    std::vector<std::string> name;
    std::ifstream file(Form("%s", filename.Data()));
    // std::string remove = "/Users/mbp-cdm-01/Desktop/dottorato/1Anno/QAStrangeness/TriggerForRun3/EventFiltering2023/AnalysisResults_merged_";
    std::string remove = "/Users/mbp-cdm-01/Desktop/dottorato/1Anno/QAStrangeness/TriggerForRun3/EventFiltering2023/AnalysisResults_";
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
    Float_t RelError = 0;
    Int_t TotEvt = 0;
    for (int i = 0; i < nfiles; i++)
    {
        RelError = 0;
        hEvSel[i] = (TH1F *)f[i]->Get("lf-strangeness-filter/hProcessedEvents");
        TotEvt = hEvSel[i]->GetBinContent(1);
        for (int j = 1; j <= hEvSel[i]->GetNbinsX(); j++)
        {
            RelError = sqrt(1. / hEvSel[i]->GetBinContent(j) + 1. / TotEvt);
            hEvSel[i]->SetBinContent(j, hEvSel[i]->GetBinContent(j) / TotEvt);
            hEvSel[i]->SetBinError(j, hEvSel[i]->GetBinContent(j) * RelError);
            if (j==1) hEvSel[i]->SetBinError(j, 0);
        }
        hEvSel[i]->SetLineColor(colors[i]);
    }

    TH1F *hEvSel23[nfiles];
    for (int i = 0; i < nfiles; i++)
    {
        hEvSel23[i] = new TH1F(Form("hEvSel23_%d", i), ";;Selectivity", 5, 0, 5);
        hEvSel23[i]->GetXaxis()->SetBinLabel(1, "Omega #times 0.2");
        hEvSel23[i]->GetXaxis()->SetBinLabel(2, "3Xi");
        hEvSel23[i]->GetXaxis()->SetBinLabel(3, "Xi-N");
        hEvSel23[i]->GetXaxis()->SetBinLabel(4, "Omega large R");
        hEvSel23[i]->GetXaxis()->SetBinLabel(5, "Total");
        hEvSel23[i]->SetLineColor(colors[i]);
        hEvSel23[i]->SetLineWidth(2);

        hEvSel23[i]->SetBinContent(1, hEvSel[i]->GetBinContent(hEvSel[i]->GetXaxis()->FindBin("#Omega")) * 0.2);
        hEvSel23[i]->SetBinContent(2, hEvSel[i]->GetBinContent(hEvSel[i]->GetXaxis()->FindBin("3#Xi")));
        hEvSel23[i]->SetBinContent(3, hEvSel[i]->GetBinContent(hEvSel[i]->GetXaxis()->FindBin("#Xi-YN")));
        hEvSel23[i]->SetBinContent(4, hEvSel[i]->GetBinContent(hEvSel[i]->GetXaxis()->FindBin("#Omega high radius")));
        double tot = hEvSel[i]->GetBinContent(hEvSel[i]->GetXaxis()->FindBin("#Omega")) * 0.2 + hEvSel[i]->GetBinContent(hEvSel[i]->GetXaxis()->FindBin("3#Xi")) + hEvSel[i]->GetBinContent(hEvSel[i]->GetXaxis()->FindBin("#Xi-YN") + hEvSel[i]->GetBinContent(hEvSel[i]->GetXaxis()->FindBin("#Omega high radius")));
        hEvSel23[i]->SetBinContent(5, tot);

        hEvSel23[i]->SetBinError(1, hEvSel[i]->GetBinError(hEvSel[i]->GetXaxis()->FindBin("#Omega")) * 0.2);
        hEvSel23[i]->SetBinError(2, hEvSel[i]->GetBinError(hEvSel[i]->GetXaxis()->FindBin("3#Xi")));
        hEvSel23[i]->SetBinError(3, hEvSel[i]->GetBinError(hEvSel[i]->GetXaxis()->FindBin("#Xi-YN")));
        hEvSel23[i]->SetBinError(4, hEvSel[i]->GetBinError(hEvSel[i]->GetXaxis()->FindBin("#Omega high radius")));
        hEvSel23[i]->SetBinError(5, 0);

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

    TH1F *h = new TH1F(Form("h%d", 0), ";;Selectivity", 5, 0, 5);
    h->GetXaxis()->SetBinLabel(1, "Omega #times 0.2");
    h->GetXaxis()->SetBinLabel(2, "3Xi");
    h->GetXaxis()->SetBinLabel(3, "Xi-N");
    h->GetXaxis()->SetBinLabel(4, "Omega large R");
    h->GetXaxis()->SetBinLabel(5, "Total");
    h->SetStats(0);
    h->GetYaxis()->SetRangeUser(1E-8, 1E-1);
    h->Draw();
    for (int i = 0; i < nfiles; i++)
    {
        hEvSel23[i]->Draw("SAME e");
    }
    hEvSel23[0]->Draw("SAME e");

    TLine *l = new TLine(0, StrLimit, 5, StrLimit);
    l->SetLineStyle(7);
    l->SetLineColor(kAzure + 10);
    l->SetLineWidth(3);
    l->Draw();
    TLine *l1 = new TLine(0, 5E-5, 5, 5E-5);
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

    TH1F *hRatio = new TH1F(Form("hRatio%d", 0), Form(";;Ratio to %s", name[ChosenPeriod].c_str()), 5, 0, 5);
    hRatio->GetXaxis()->SetBinLabel(1, "Omega #times 0.2");
    hRatio->GetXaxis()->SetBinLabel(2, "3Xi");
    hRatio->GetXaxis()->SetBinLabel(3, "Xi-N");
    hRatio->GetXaxis()->SetBinLabel(4, "Omega large R");
    hRatio->GetXaxis()->SetBinLabel(5, "Total");
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