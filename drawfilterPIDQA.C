void drawfilterPIDQA(TString filename = "listBis.txt") {

    std::vector<std::string> name;
    std::ifstream file(Form("%s", filename.Data()));
    std::string remove = "/Users/mbp-cdm-01/Desktop/dottorato/1Anno/QAStrangeness/TriggerForRun3/EventFiltering2023/AnalysisResults_merged_";
    std::string remove2 = "_round2.root";

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
    TH2D *hTPCNsigmaBachKa[nfiles];
    TH2D *hTPCNsigmaBachPr[nfiles];
    TH2D *hTPCNsigmaPi[nfiles];

    for (int i = 0; i < nfiles; i++)
    {

        f[i] = TFile::Open(Form("%s%s_round2.root", remove.c_str(), name[i].c_str()));

        hTPCNsigmaBachKa[i] = (TH2D *)f[i]->Get("lf-strangeness-filter/QAHistos/hTPCNsigmaOmegaBachKaMinus");
        hTPCNsigmaBachPr[i] = (TH2D *)f[i]->Get("lf-strangeness-filter/QAHistos/hTPCNsigmaOmegaV0Proton");
        hTPCNsigmaPi[i] = (TH2D *)f[i]->Get("lf-strangeness-filter/QAHistos/hTPCNsigmaOmegaV0PiMinus");

        hTPCNsigmaBachKa[i]->SetTitle(Form("%s K bach", name[i].c_str()));
        hTPCNsigmaBachKa[i]->SetStats(0);
        hTPCNsigmaBachKa[i]->Draw("colz");
        hTPCNsigmaBachKa[i]->GetYaxis()->SetRangeUser(0., 4.);
        hTPCNsigmaBachKa[i]->GetXaxis()->SetRangeUser(-8., 8.);

        hTPCNsigmaBachPr[i]->SetTitle(Form("%s pr V0", name[i].c_str()));
        hTPCNsigmaBachPr[i]->SetStats(0);
        hTPCNsigmaBachPr[i]->Draw("colz");
        hTPCNsigmaBachPr[i]->GetYaxis()->SetRangeUser(0., 6.);
        hTPCNsigmaBachPr[i]->GetXaxis()->SetRangeUser(-8., 8.);

        hTPCNsigmaPi[i]->SetTitle(Form("%s pi V0", name[i].c_str()));
        hTPCNsigmaPi[i]->SetStats(0);
        hTPCNsigmaPi[i]->Draw("colz");
        hTPCNsigmaPi[i]->GetYaxis()->SetRangeUser(0., 2.5);
        hTPCNsigmaPi[i]->GetXaxis()->SetRangeUser(-8., 8.);

    }

    TCanvas *c4 = new TCanvas("c4", "c4", 2000, 1000);
    c4->Divide(4, 3);
    for (int i = 0; i < nfiles; i++)
    {
        c4->cd(i + 1);
        c4->cd(i + 1)->SetRightMargin(0.15);
        c4->cd(i + 1)->SetLogz();
        hTPCNsigmaBachKa[i]->Draw("colz");
    }

    TCanvas *c5 = new TCanvas("c5", "c5", 2000, 1000);
    c5->Divide(4, 3);
    for (int i = 0; i < nfiles; i++)
    {
        c5->cd(i + 1);
        c5->cd(i + 1)->SetRightMargin(0.15);
        c5->cd(i + 1)->SetLogz();
        hTPCNsigmaBachPr[i]->Draw("colz");
    }

    TCanvas *c6 = new TCanvas("c6", "c6", 2000, 1000);
    c6->Divide(4, 3);
    for (int i = 0; i < nfiles; i++)
    {
        c6->cd(i + 1);
        c6->cd(i + 1)->SetRightMargin(0.15);
        c6->cd(i + 1)->SetLogz();
        hTPCNsigmaPi[i]->Draw("colz");
    }

    c4->SaveAs("PIDK.png");
    c5->SaveAs("PIDPr.png");
    c6->SaveAs("PIDPi.png");
}
