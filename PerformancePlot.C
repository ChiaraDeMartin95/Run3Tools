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

/// @brief
/// @param histo
/// @param Low
/// @param Up
/// @param color
/// @param style
/// @param titleX
/// @param titleY
/// @param title
/// @param XRange
/// @param XLow
/// @param XUp
/// @param xOffset
/// @param yOffset
/// @param mSize
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
    histo->GetXaxis()->SetLabelSize(0.04); // 0.05
    histo->GetXaxis()->SetTitleSize(0.05); // 0.05
    histo->GetXaxis()->SetTitleOffset(xOffset);
    histo->GetYaxis()->SetTitle(titleY);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetTitleOffset(yOffset);
    histo->SetTitle(title);
}

const Float_t UpperLimitLSBOmega = 1.655;  // upper limit of fit of left sidebands for omega
const Float_t LowerLimitRSBOmega = 1.689;  // lower limit of fit of right sidebands for omega
const Float_t UpperLimitLSBXi = 1.302;     // upper limit of fit of left sidebands for Xi
const Float_t LowerLimitRSBXi = 1.34;      // lower limit of fit of right sidebands for Xi
const Float_t UpperLimitLSBLambda = 1.107; // upper limit of fit of left sidebands for Lambda
const Float_t LowerLimitRSBLambda = 1.123; // lower limit of fit of right sidebands for Lambda

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
    else if (par[2] == 1 || par[2] == 2)
    {
        LimInf = UpperLimitLSBLambda; //
        LimSup = LowerLimitRSBLambda; //
    }
    else if (par[3] == 3 || par[3] == 4)
    {
        LimInf = UpperLimitLSBXi; // 1.31
        LimSup = LowerLimitRSBXi; // 1.335
    }
    else if (par[3] == 5 || par[3] == 6)
    {
        LimInf = UpperLimitLSBOmega; // 1.66
        LimSup = LowerLimitRSBOmega; // 1.692
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
    else if (par[2] == 1 || par[2] == 2)
    {
        LimInf = UpperLimitLSBLambda; //
        LimSup = LowerLimitRSBLambda; //
    }
    else if (par[2] == 3 || par[2] == 4)
    {
        LimInf = UpperLimitLSBXi; // 1.31
        LimSup = LowerLimitRSBXi; // 1.335
    }
    else if (par[2] == 5 || par[2] == 6)
    {
        LimInf = UpperLimitLSBOmega; // 1.65
        LimSup = LowerLimitRSBOmega; // 1.69
    }
    if (reject && x[0] > LimInf && x[0] < LimSup)
    {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1] * x[0];
}

TString titlePt = "p_{T} (GeV/c)";
TString titleYield = "1/N_{ev} dN/dp_{T}";

const Int_t numPart = 7;
TString TitleInvMass[numPart] = {"(#pi^{+}, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(p, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(#bar{p}, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(#Lambda, #pi^{-}) invariant mass (GeV/#it{c}^{2})"};
TString namehisto[numPart] = {"h3dMassK0Short", "", "", "h2dMassXiMinus", "h2dMassXiPlus", "h2dMassOmegaMinus", "h2dMassOmegaPlus"};
Float_t LowLimitMass[numPart] = {0.458, 1.102, 1.102, 1.29, 1.29, 1.65, 1.65}; // 0.44
Float_t UpLimitMass[numPart] = {0.538, 1.129, 1.129, 1.35, 1.35, 1.7, 1.7};    // 0.55
Float_t LowMassRange[numPart] = {0.48, 1.09, 1.09, 1.31};
Float_t UpMassRange[numPart] = {0.51, 1.14, 1.14, 1.33};

Float_t min_range_signal[numPart] = {0.46, 1.105, 1.105, 1.31, 1.31, 1.655, 1.655}; // estremi region fit segnale (gaussiane)
Float_t max_range_signal[numPart] = {0.535, 1.125, 1.125, 1.334, 1.334, 1.685, 1.685};
Float_t min_histo[numPart] = {0.42, 1.1, 1.1, 1.30, 1.30, 1.63, 1.63}; // estremi del range degli istogrammi
Float_t max_histo[numPart] = {0.57, 1.13, 1.13, 1.342, 1.342, 1.71, 1.71};
Float_t liminf[numPart] = {0.45, 1.1, 1.1, 1.30, 1.30, 1.64, 1.64}; // estremi regione fit del bkg e total
// Float_t liminf[numPart] = {0.45, 1.1153, 1.1153, 1.29, 1.29, 1.66, 1.66}; // estremi regione fit del bkg e total
Float_t limsup[numPart] = {0.545, 1.13, 1.13, 1.342, 1.342, 1.7, 1.7};
// Float_t limsup[numPart] = {0.545, 1.1168, 1.1168, 1.35, 1.35, 1.685, 1.685};

const Float_t massParticle[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
TString Spart[numPart] = {"K0S", "Lambda", "AntiLambda", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
TString SpartLegend[numPart] = {"K0S", "Lambda", "AntiLambda", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};

void PerformancePlot(Int_t part = 0,
                     Float_t PtMin = 0.9,
                     Float_t PtMax = 10,
                     Bool_t UseOneGauss = 1,
                     TString SPeriod = "LHC23zy_cpass0",
                     TString SPathIn = "../Run3QA/Periods/PbPb2023/LHC23zy_cpass0/AnalysisResults_qatask_LHC23zy_cpass0.root" /*"../TriggerForRun3/AnalysisResults_FinalTOT_NoTOF.root"*/,
                     TString PathOut = "../Run3QA/Periods/PbPb2023/LHC23zy_cpass0/InvMassPlot",
                     Bool_t ispp = 0,
                     Bool_t UseTwoGauss = 1,
                     Bool_t isBkgParab = 1,
                     Bool_t isMeanFixedPDG = 0,
                     Float_t sigmacentral = 4)
{

    TFile *fileIn = new TFile(SPathIn, "");
    if (!fileIn)
        return;
    TDirectoryFile *dir = (TDirectoryFile *)fileIn->Get("v0cascades-q-a");
    TDirectoryFile *dir1;
    if (part >= 3)
        dir1 = (TDirectoryFile *)dir->Get("histos-Casc");
    else
        dir1 = (TDirectoryFile *)dir->Get("histos-V0");
    TH2F *histo2d = (TH2F *)dir1->Get("InvMass" + Spart[part]);
    TH1F *histo = (TH1F *)histo2d->ProjectionY("histo", histo2d->GetXaxis()->FindBin(PtMin + 0.001), histo2d->GetXaxis()->FindBin(PtMax - 0.001));
    // histo->Scale(1. / histo->GetIntegral() / histo->GetXaxis()->GetBinWidth(1));
    histo->Scale(1. / histo->Integral());

    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 800);
    StyleCanvas(canvas, 0.18, 0.03, 0.02, 0.14); // L, R, T, B

    TLegend *legend = new TLegend(0.23, 0.65, 0.75, 0.96);
    legend->SetFillStyle(0);
    legend->SetMargin(0);
    legend->SetTextSize(0.04);
    legend->SetTextAlign(12);
    // legend->AddEntry("", "#bf{This work}", "");
    legend->AddEntry("", "#bf{ALICE Work in progress}", "");
    if (ispp)
        legend->AddEntry("", "Run 3, pp #sqrt{#it{s}} = 13.6 TeV", "");
    else
        legend->AddEntry("", "Run 3, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV", "");
    legend->AddEntry("", SPeriod, "");
    legend->AddEntry("", Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", PtMin, PtMax), "");
    legend->AddEntry("", "|y| < 0.5", "");
    if (part == 0)
        legend->AddEntry("", "K^{0}_{S} #rightarrow #pi^{+} #pi^{-}", "");
    else if (part == 1)
        legend->AddEntry("", "#Lambda #rightarrow p #pi^{-}", "");
    else if (part == 2)
        legend->AddEntry("", "#bar{#Lambda} #rightarrow #bar{p} #pi^{+}", "");
    else if (part == 3)
        legend->AddEntry("", "#Xi^{-} #rightarrow #Lambda #pi^{-} #rightarrow p #pi^{-} #pi^{-}", "");
    else if (part == 4)
        legend->AddEntry("", "#Xi^{+} #rightarrow #bar{#Lambda} #pi^{+} #rightarrow #bar{p} #pi^{+} #pi^{+}", "");
    else if (part == 5)
        legend->AddEntry("", "#Omega^{-} #rightarrow #Lambda K^{-} #rightarrow p #pi^{-} K^{-}", "");
    else if (part == 6)
        legend->AddEntry("", "#Omega^{+} #rightarrow #bar{#Lambda} K^{+} #rightarrow #bar{p} #pi^{+} K^{+}", "");

    TLegend *legendfit = new TLegend(0.21, 0.56, 0.74, 0.64);
    legendfit->SetFillStyle(0);
    legendfit->SetMargin(0.1);
    legendfit->SetTextSize(0.03);
    legendfit->SetTextAlign(12);

    // Fit 1
    Double_t parTwoGaussParab[9];
    Double_t parTwoGaussRetta[8];
    Double_t parOneGaussParab[6];
    Double_t parOneGaussRetta[5];

    TFitResultPtr fFitResultPtr0;
    TFitResultPtr fFitResultPtr1;

    Float_t mean = {0};
    Float_t errmean = {0};
    Float_t sigma = {0};
    Float_t errsigma = {0};
    Float_t b = {0};
    Float_t errb = {0};
    Float_t SSB = {0};
    Float_t errSSB = {0};
    Float_t entries_range = {0};
    Float_t Yield = {0};
    Float_t ErrYield = {0};
    Float_t TotYield = 0;

    canvas->cd();

    TF1 *functionsFirst = new TF1("1f_%i", "gaus", min_range_signal[part], max_range_signal[part]);
    functionsFirst->SetLineColor(881);
    functionsFirst->SetParameter(1, massParticle[part]);
    functionsFirst->SetParName(0, "norm");
    functionsFirst->SetParName(1, "mean");
    functionsFirst->SetParName(2, "sigma");
    functionsFirst->SetParLimits(1, min_range_signal[part], max_range_signal[part]);
    functionsFirst->SetParLimits(2, 0.001, 0.1);
    functionsFirst->SetParLimits(0, 0, 1.1 * histo->GetBinContent(histo->GetMaximumBin()));

    TF1 *functionsSecond = new TF1("2f_%i", "gaus", min_range_signal[part], max_range_signal[part]);
    functionsSecond->SetLineColor(867);
    functionsSecond->SetParameter(1, massParticle[part]);
    functionsSecond->SetParName(0, "norm");
    functionsSecond->SetParName(1, "mean");
    functionsSecond->SetParName(2, "sigma");
    functionsSecond->SetParLimits(1, min_range_signal[part], max_range_signal[part]);
    functionsSecond->SetParLimits(2, 0.001, 0.15);
    functionsSecond->SetParLimits(0, 0, 1.1 * histo->GetBinContent(histo->GetMaximumBin()));

    TF1 *functions1 = new TF1("1f_%i_final", "gaus", min_range_signal[part], max_range_signal[part]);
    functions1->SetLineColor(kRed); // 867
    functions1->SetParName(0, "norm");
    functions1->SetParName(1, "mean");
    functions1->SetParName(2, "sigma");

    TF1 *functions2 = new TF1("2f_%i_final", "gaus", min_range_signal[part], max_range_signal[part]);
    functions2->SetLineColor(kMagenta); // 891
    functions2->SetParName(0, "norm");
    functions2->SetParName(1, "mean");
    functions2->SetParName(2, "sigma");

    TF1 *bkg1 = new TF1("bkg1%i", "pol1", min_histo[part], max_histo[part]);
    bkg1->SetLineColor(418);

    TF1 *bkg2 = new TF1("bkg2%i", "pol2", min_histo[part], max_histo[part]);
    bkg2->SetLineColor(kBlue);

    TF1 *bkgretta = new TF1("retta%i", fretta, liminf[part], limsup[part], 3);
    bkgretta->SetLineColor(418);
    bkgretta->FixParameter(2, part);

    TF1 *bkgparab = new TF1("parab%i", fparab, liminf[part], limsup[part], 4);
    bkgparab->SetLineColor(kBlue);
    bkgparab->FixParameter(3, part);

    TF1 *total;
    TF1 *totalbis;
    if (UseOneGauss)
    {
        cout << "\n\e[36mFit with one gauss only: \e[39m"
             << " Pt: " << PtMin << "--" << PtMax << endl;

        if (isBkgParab)
            total = new TF1("total%i", "gaus(0)+pol2(3)", LowLimitMass[part], UpLimitMass[part]);
        else
            total = new TF1("total%i", "gaus(0)+pol1(3)", LowLimitMass[part], UpLimitMass[part]);
        total->SetLineColor(7);
        total->SetParName(0, "norm");
        total->SetParName(1, "mean");
        total->SetParName(2, "sigma");

        // Use 1 Gauss
        cout << "\n\n fit gauss " << endl;
        histo->Fit(functionsFirst, "RB0");

        bkg1->SetRange(min_histo[part], max_histo[part]);
        bkg2->SetRange(min_histo[part], max_histo[part]);
        bkgparab->SetRange(liminf[part], limsup[part]);
        bkgretta->SetRange(liminf[part], limsup[part]);
        total->SetRange(liminf[part], limsup[part]);

        cout << "\n\n fit bkg " << endl;
        if (isBkgParab)
            histo->Fit(bkgparab, "RB0");
        else
            histo->Fit(bkgretta, "RB0");

        functionsFirst->GetParameters(&parOneGaussParab[0]);
        functionsFirst->GetParameters(&parOneGaussRetta[0]);
        if (isBkgParab)
        {
            bkgparab->GetParameters(&parOneGaussParab[3]);
            total->SetParameters(parOneGaussParab);
        }
        else
        {
            bkgretta->GetParameters(&parOneGaussRetta[3]);
            total->SetParameters(parOneGaussRetta);
        }

        cout << "\n\n fit total " << endl;
        if (Spart[part] == "XiMinus" || Spart[part] == "XiPlus")
        {
            total->SetParLimits(0, 0.08 * histo->GetBinContent(histo->GetMaximumBin()), histo->GetBinContent(histo->GetMaximumBin()));
            total->SetParLimits(1, 1.318, 1.326);
            total->SetParLimits(2, 0.0012, 0.010);
            if (isMeanFixedPDG)
            {
                total->FixParameter(1, massParticle[part]);
                total->FixParameter(4, massParticle[part]);
            }
        }
        if (Spart[part] == "OmegaMinus" || Spart[part] == "OmegaPlus")
        {
            total->SetParLimits(0, 0.08 * histo->GetBinContent(histo->GetMaximumBin()), histo->GetBinContent(histo->GetMaximumBin()));
            total->SetParLimits(1, 1.66, 1.68);
            total->SetParLimits(2, 0.0012, 0.010);
            if (isMeanFixedPDG)
            {
                total->FixParameter(1, massParticle[part]);
                total->FixParameter(4, massParticle[part]);
            }
        }
        if (Spart[part] == "K0S")
        {
            total->SetParLimits(0, 0.08 * histo->GetBinContent(histo->GetMaximumBin()), histo->GetBinContent(histo->GetMaximumBin()));
            total->SetParLimits(1, 0.485, 0.505);
            total->SetParLimits(2, 0.001, 0.01);
            bkgretta->SetParameter(1, 0);
            if (isMeanFixedPDG)
            {
                total->FixParameter(1, massParticle[part]);
                total->FixParameter(4, massParticle[part]);
            }
            cout << "max value " << histo->GetBinContent(histo->GetMaximumBin()) << endl;
        }
        else if (Spart[part] == "Lambda")
        {
            total->SetParLimits(0, 0.08 * histo->GetBinContent(histo->GetMaximumBin()), histo->GetBinContent(histo->GetMaximumBin()));
            total->SetParLimits(1, 1.11, 1.12);
            total->SetParLimits(2, 0.001, 0.01);
            bkgretta->SetParameter(1, 0);
            if (isMeanFixedPDG)
            {
                total->FixParameter(1, massParticle[part]);
                total->FixParameter(4, massParticle[part]);
            }
            cout << "max value " << histo->GetBinContent(histo->GetMaximumBin()) << endl;
        }

        fFitResultPtr0 = histo->Fit(total, "SRB+"); // per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0

        totalbis = (TF1 *)total->Clone();
        fFitResultPtr1 = fFitResultPtr0;

        functions1->FixParameter(0, total->GetParameter(0));
        functions1->FixParameter(1, total->GetParameter(1));
        functions1->FixParameter(2, total->GetParameter(2));

        totalbis->FixParameter(0, 0);
        totalbis->FixParameter(1, 0);
        totalbis->FixParameter(2, 0);

        if (isBkgParab)
        {
            bkg2->FixParameter(0, total->GetParameter(3));
            bkg2->FixParameter(1, total->GetParameter(4));
            bkg2->FixParameter(2, total->GetParameter(5));
            bkgparab->FixParameter(0, total->GetParameter(3));
            bkgparab->FixParameter(1, total->GetParameter(4));
            bkgparab->FixParameter(2, total->GetParameter(5));
        }
        else
        {
            bkg1->FixParameter(0, total->GetParameter(3));
            bkg1->FixParameter(1, total->GetParameter(4));
            bkgretta->FixParameter(0, total->GetParameter(3));
            bkgretta->FixParameter(1, total->GetParameter(4));
        }
    }
    else
    {
        // Use two gauss
        cout << "\n\e[35mFit with two gauss \e[39m" << endl;

        if (isBkgParab)
            total = new TF1("total%i", "gaus(0)+gaus(3)+pol2(6)", LowLimitMass[part], UpLimitMass[part]);
        else
            total = new TF1("total%i", "gaus(0)+gaus(3)+pol1(6)", LowLimitMass[part], UpLimitMass[part]);
        total->SetParName(0, "norm");
        total->SetParName(1, "mean");
        total->SetParName(2, "sigma");
        total->SetParName(3, "norm2");
        total->SetParName(4, "mean2");
        total->SetParName(5, "sigma2");

        cout << "\n\n fit gauss1 " << endl;
        histo->Fit(functionsFirst, "RB0");
        cout << "\n\n fit gauss2 " << endl;
        histo->Fit(functionsSecond, "RB0");

        bkg1->SetRange(min_histo[part], max_histo[part]);
        bkg2->SetRange(min_histo[part], max_histo[part]);
        bkgparab->SetRange(liminf[part], limsup[part]);
        bkgretta->SetRange(liminf[part], limsup[part]);
        total->SetRange(liminf[part], limsup[part]);

        cout << "\n\n fit bkg " << endl;
        if (isBkgParab)
            histo->Fit(bkgparab, "RB0");
        else
            histo->Fit(bkgretta, "RB0");

        functionsFirst->GetParameters(&parTwoGaussParab[0]);
        functionsFirst->GetParameters(&parTwoGaussRetta[0]);
        functionsSecond->GetParameters(&parTwoGaussParab[3]);
        functionsSecond->GetParameters(&parTwoGaussRetta[3]);
        if (isBkgParab)
        {
            bkgparab->GetParameters(&parTwoGaussParab[6]);
            total->SetParameters(parTwoGaussParab);
        }
        else
        {
            bkgretta->GetParameters(&parTwoGaussRetta[6]);
            total->SetParameters(parTwoGaussRetta);
        }

        cout << "\n\n fit total " << endl;
        if (Spart[part] == "XiMinus" || Spart[part] == "XiPlus")
        {
            total->SetParLimits(0, 0.08 * histo->GetBinContent(histo->GetMaximumBin()), histo->GetBinContent(histo->GetMaximumBin()));
            total->SetParLimits(1, 1.318, 1.326);
            total->SetParLimits(2, 0.0012, 0.010);                                                                                     // it was 0.001
            total->SetParLimits(3, 0.08 * histo->GetBinContent(histo->GetMaximumBin()), histo->GetBinContent(histo->GetMaximumBin())); // maximum was wothout 0.3
            total->SetParLimits(4, 1.318, 1.326);
            if (isMeanFixedPDG)
            {
                total->FixParameter(1, massParticle[part]);
                total->FixParameter(4, massParticle[part]);
            }
        }
        if (Spart[part] == "K0S")
        {
            total->SetParLimits(0, 0.08 * histo->GetBinContent(histo->GetMaximumBin()), histo->GetBinContent(histo->GetMaximumBin()));
            total->SetParLimits(1, 0.485, 0.505);
            total->SetParLimits(2, 0.001, 0.01);
            total->SetParLimits(3, 0.08 * histo->GetBinContent(histo->GetMaximumBin()), histo->GetBinContent(histo->GetMaximumBin()));
            total->SetParLimits(4, 0.485, 0.505);
            total->SetParLimits(5, 0.001, 0.015);
            // total->SetParLimits(6, -2000, 2000);
            bkgretta->SetParameter(1, 0);
            if (isMeanFixedPDG)
            {
                total->FixParameter(1, massParticle[part]);
                total->FixParameter(4, massParticle[part]);
            }
            cout << "max value " << histo->GetBinContent(histo->GetMaximumBin()) << endl;
        }

        fFitResultPtr0 = histo->Fit(total, "SRB+"); // per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0
        // la gaussiana più larga deve esserte quella più bassa

        totalbis = (TF1 *)total->Clone();
        fFitResultPtr1 = fFitResultPtr0;

        functions1->FixParameter(0, total->GetParameter(0));
        functions1->FixParameter(1, total->GetParameter(1));
        functions1->FixParameter(2, total->GetParameter(2));
        functions2->FixParameter(0, total->GetParameter(3));
        functions2->FixParameter(1, total->GetParameter(4));
        functions2->FixParameter(2, total->GetParameter(5));

        totalbis->FixParameter(0, 0);
        totalbis->FixParameter(1, 0);
        totalbis->FixParameter(2, 0);
        totalbis->FixParameter(3, 0);
        totalbis->FixParameter(4, 0);
        totalbis->FixParameter(5, 0);

        if (isBkgParab)
        {
            bkg2->FixParameter(0, total->GetParameter(6));
            bkg2->FixParameter(1, total->GetParameter(7));
            bkg2->FixParameter(2, total->GetParameter(8));
            bkgparab->FixParameter(0, total->GetParameter(6));
            bkgparab->FixParameter(1, total->GetParameter(7));
            bkgparab->FixParameter(2, total->GetParameter(8));
        }
        else
        {
            bkg1->FixParameter(0, total->GetParameter(6));
            bkg1->FixParameter(1, total->GetParameter(7));
            bkgretta->FixParameter(0, total->GetParameter(6));
            bkgretta->FixParameter(1, total->GetParameter(7));
        }
    }
    canvas->cd();
    Float_t UpperCutHisto = 1.4;
    if (part == 3 || part == 4)
        UpperCutHisto = 1.7;
    else if (part >= 5)
        UpperCutHisto = 1.8;

    StyleHisto(histo, 0.0001, UpperCutHisto * histo->GetBinContent(histo->GetMaximumBin()), 1, 20,
               "#it{m} (GeV/#it{c}^{2})", "Normalized counts" /* per 1.0 MeV/#it{c}^{2}"*/, "", 1, LowLimitMass[part] + 0.001, UpLimitMass[part] - 0.001, 1.2, 1.8, 1.2);
    if (part == 5)
        histo->GetYaxis()->SetTitleOffset(1.6);
    histo->Draw("pe");
    legendfit->AddEntry(total, "Gaussian fits + bkg.", "l");
    if (isBkgParab)
    {
        bkg2->SetLineColor(kBlack);
        bkg2->SetLineStyle(8);
        bkg2->Draw("same");
        legendfit->AddEntry(bkg2, "bkg.", "l");
    }
    else
    {
        bkg1->SetLineColor(kBlack);
        bkg1->SetLineStyle(8);
        bkg1->Draw("same");
        legendfit->AddEntry(bkg1, "bkg.", "l");
    }
    // functions1->Draw("same");
    // if (!UseOneGauss) functions2->Draw("same");
    total->SetRange(LowLimitMass[part], UpLimitMass[part]);
    total->SetLineColor(kRed + 1);
    total->Draw("same");
    legend->Draw("");
    legendfit->Draw("");

    canvas->SaveAs(PathOut + Spart[part] + ".pdf");
    canvas->SaveAs(PathOut + Spart[part] + ".png");
}