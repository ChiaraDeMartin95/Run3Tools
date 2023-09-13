#include <iostream>
#include <fstream>

void drawfiltersQA(TString filename = "list.txt")
{

    std::vector<std::string> name;
    std::ifstream file(Form("%s", filename.Data()));
    std::string remove = "/Users/mbp-cdm-01/Desktop/dottorato/1Anno/QAStrangeness/TriggerForRun3/EventFiltering2023/AnalysisResults_merged_";
    std::string remove2 = "_round2.root";

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

    Int_t colors[20] = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kYellow, kMagenta + 3, kAzure, kViolet, kTeal, kSpring, kGray, kBlue + 2, kRed + 4, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 2};
    const int nfiles = name.size();
    cout << "Number of files = " << nfiles << endl;

    TFile *f[nfiles];
    for (int i = 0; i < nfiles; i++)
    {
        cout << "name " << name[i].c_str() << endl;
        f[i] = TFile::Open(Form("%s%s_round2.root", remove.c_str(), name[i].c_str()));
    }

    TH1F *hEvSel[nfiles];
    for (int i = 0; i < nfiles; i++)
    {
        hEvSel[i] = (TH1F *)f[i]->Get("lf-strangeness-filter/hProcessedEvents");
        hEvSel[i]->Scale(1. / hEvSel[i]->GetBinContent(1));
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

        cout << "Total rejection factor " << name[i].c_str() << " = " << tot << endl;
    }

    TCanvas *c = new TCanvas("c", "c", 900, 1000);
    c->SetLogy();
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.05);

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
        hEvSel23[i]->Draw("SAME");
    }
    hEvSel23[0]->Draw("SAME");

    TLine *l = new TLine(0, 4.3E-5, 5, 4.3E-5);
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

    c->SaveAs("images/rejfactors.png");
}