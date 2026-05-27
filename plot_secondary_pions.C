// plot_secondary_pions.C
// Usage: root -l -b -q plot_secondary_pions.C

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

const int NCAT = 13;
const char* catName[NCAT + 1] = {
    "Data",
    "PiElas", "PiRes", "PiQE", "PiDCEX", "PiABS", "PiCEX",
    "Muon", "misID:cosmic", "misID:p", "misID:pi", "misID:mu",
    "misID:e/#gamma", "misID:other"};
const int catColor[NCAT + 1] = {
    kBlack,
    kRed + 1, kYellow + 1, kCyan + 1, kBlue + 2, kGreen + 2, kMagenta + 1,
    kCyan - 6, kGray + 1, kOrange + 1, kRed - 7, kBlue - 7, kGreen - 6, kGray + 2};

double ComputeGlobalScale(TFile* fMC, TFile* fData) {
    TDirectory* dMC = (TDirectory*)fMC->Get("Beam_scraper");
    TDirectory* dData = (TDirectory*)fData->Get("Beam_scraper");
    if (!dMC || !dData) {
        printf(
            "WARNING: Beam_scraper directory missing in MC or Data file. "
            "Did you re-run pionqe0p5.C after adding the scraper-stage "
            "Beam_P_beam_inst fills? Falling back to scale = 1.0\n");
        return 1.;
    }
    const TString prefix = "Beam_scraper_Beam_P_beam_inst";
    TH1D* hData = (TH1D*)dData->Get(prefix);
    if (!hData) hData = (TH1D*)dData->Get(prefix + "_0");
    double dataTotal = (hData) ? hData->Integral() : -1.;

    double mcTotal = 0.;
    for (int i = 1; i <= NCAT; i++) {
        TH1D* h = (TH1D*)dMC->Get(prefix + Form("_%d", i));
        if (h) mcTotal += h->Integral();
    }
    if (dataTotal <= 0. || mcTotal <= 0.) {
        printf(
            "WARNING: cannot compute scale (mcTotal=%.1f, dataTotal=%.1f). "
            "Falling back to 1.0\n",
            mcTotal, dataTotal);
        return 1.;
    }
    double scale = dataTotal / mcTotal;
    printf(
        "\n==> Global MC->Data scale (Beam_scraper / Beam_P_beam_inst): "
        "%.4f  (Data = %.1f, MC = %.1f)\n\n",
        scale, dataTotal, mcTotal);
    return scale;
}

void DrawStackedPlot(TFile* fMC, TFile* fData,
                     TString dir, TString varname,
                     TString xtitle, TString outname,
                     double xmin = -1, double xmax = -1,
                     int rebin = 1,
                     double globalScale = 1.,
                     TString plotdir = "plots_sec") {
    TString prefix = dir + "_" + varname;

    // Data histogram (unsuffixed)
    TDirectory* dMC = (TDirectory*)fMC->Get(dir);
    TDirectory* dData = (TDirectory*)fData->Get(dir);
    if (!dMC || !dData) {
        printf("WARNING: directory %s not found, skipping\n", dir.Data());
        return;
    }

    TH1D* hData = (TH1D*)dData->Get(prefix);
    if (!hData) hData = (TH1D*)dData->Get(prefix + "_0");
    if (!hData) {
        printf("WARNING: Data histogram %s not found, skipping\n", prefix.Data());
        return;
    }
    hData->SetDirectory(0);
    if (rebin > 1) hData->Rebin(rebin);

    // MC stack
    THStack* stack = new THStack("stack", "");
    TH1D* hMC_total = nullptr;
    vector<pair<int, TH1D*>> hMC_cat;

    for (int i = 1; i <= NCAT; i++) {
        TH1D* h = (TH1D*)dMC->Get(prefix + Form("_%d", i));
        if (!h) continue;
        if (rebin > 1) h->Rebin(rebin);
        h->SetFillColor(catColor[i]);
        h->SetLineColor(catColor[i]);
        h->SetLineWidth(1);
        hMC_cat.push_back({i, h});
        stack->Add(h);
        if (!hMC_total) {
            hMC_total = (TH1D*)h->Clone("hMC_total");
            hMC_total->SetDirectory(0);
        } else
            hMC_total->Add(h);
    }
    if (!hMC_total || hMC_total->Integral() == 0) {
        printf("WARNING: no MC for %s\n", prefix.Data());
        return;
    }

    double dataInt = hData->Integral();
    double mcInt = hMC_total->Integral();
    double perPlotScale = (mcInt > 0) ? dataInt / mcInt : 1.;
    printf("  %s: perPlotScale=%.4f vs globalScale=%.4f (ratio %.3f)\n",
           prefix.Data(), perPlotScale, globalScale,
           (globalScale > 0 ? perPlotScale / globalScale : 0.));

    for (auto& p : hMC_cat) p.second->Scale(globalScale);
    hMC_total->Scale(globalScale);
    double mcSum = hMC_total->Integral();

    double axMin = (xmin >= 0) ? xmin : hData->GetXaxis()->GetXmin();
    double axMax = (xmax > 0) ? xmax : hData->GetXaxis()->GetXmax();

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.8);
    hData->SetLineColor(kBlack);

    TCanvas* c = new TCanvas("c", "", 800, 700);
    TPad* pad1 = new TPad("pad1", "", 0., 0.28, 1., 1.);
    TPad* pad2 = new TPad("pad2", "", 0., 0., 1., 0.28);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.02);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.35);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    stack->Draw("HIST");
    stack->GetXaxis()->SetRangeUser(axMin, axMax);
    stack->GetXaxis()->SetLabelSize(0);
    stack->GetXaxis()->SetTitleSize(0);
    stack->GetYaxis()->SetTitle("Events");
    stack->GetYaxis()->SetTitleSize(0.055);
    stack->GetYaxis()->SetTitleOffset(0.85);
    stack->GetYaxis()->SetLabelSize(0.05);

    hData->GetXaxis()->SetRangeUser(axMin, axMax);
    hMC_total->GetXaxis()->SetRangeUser(axMin, axMax);
    double ymax = max(hData->GetMaximum(), hMC_total->GetMaximum());
    hData->GetXaxis()->UnZoom();
    hMC_total->GetXaxis()->UnZoom();
    stack->SetMaximum(ymax * 2.2);
    stack->SetMinimum(0.);
    hData->Draw("E1 SAME");

    // MC statistical uncertainty band
    TH1D* hMC_err = (TH1D*)hMC_total->Clone("hMC_err");
    hMC_err->SetDirectory(0);
    hMC_err->SetFillColor(kGray + 1);
    hMC_err->SetFillStyle(3344);  // meshed
    hMC_err->SetLineColor(kGray + 1);
    hMC_err->SetLineWidth(0);
    hMC_err->SetMarkerSize(0);
    hMC_err->Draw("E2 SAME");  // E2 draws the error band

    // Legend
    TLegend* leg = new TLegend(0.12, 0.60, 0.92, 0.88);
    leg->SetNColumns(4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);
    leg->SetMargin(0.12);
    for (auto& p : hMC_cat) {
        int cat = p.first;
        TH1D* h = p.second;
        double integ = h->Integral();
        double pct = (mcSum > 0) ? 100. * integ / mcSum : 0.;
        leg->AddEntry(h, Form("#bf{%s %.1f, (%.1f%%)}", catName[cat], integ, pct), "f");
    }
    // leg->AddEntry((TObject*)nullptr, Form("MC Sum %.1f", mcSum), "");
    leg->AddEntry(hMC_err, Form("MC Sum %.1f", mcSum), "f");
    leg->AddEntry(hData, Form("Observed %.0f", dataInt), "lep");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.052);
    lat.DrawLatex(0.12, 0.93, "ProtoDUNE-SP");
    lat.SetTextSize(0.044);
    lat.DrawLatex(0.35, 0.93, "#bf{#it{Preliminary}}");
    lat.SetTextSize(0.052);
    lat.DrawLatex(0.64, 0.93, "0.5 GeV/c Beam");

    // Ratio pad
    pad2->cd();
    TH1D* hRatio = (TH1D*)hData->Clone("hRatio");
    hRatio->SetDirectory(0);
    hRatio->Reset();
    for (int b = 1; b <= hRatio->GetNbinsX(); b++) {
        double nd = hData->GetBinContent(b);
        double nm = hMC_total->GetBinContent(b);
        if (nm > 0) {
            hRatio->SetBinContent(b, nd / nm);
            hRatio->SetBinError(b, TMath::Sqrt(nd) / nm);
        }
    }
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetTitle(xtitle);
    hRatio->GetXaxis()->SetRangeUser(axMin, axMax);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatio->GetYaxis()->SetTitle("Obs./Pred.");
    hRatio->GetYaxis()->SetTitleSize(0.11);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetRangeUser(0., 2.);
    hRatio->GetYaxis()->SetNdivisions(504);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.7);
    hRatio->SetLineColor(kBlack);

    //////////////////////////////

    // -- MC uncertainty band in ratio pad, centered at 1
    TH1D* hRatio_MCerr = (TH1D*)hMC_total->Clone("hRatio_MCerr");
    hRatio_MCerr->SetDirectory(0);
    for (int b = 1; b <= hRatio_MCerr->GetNbinsX(); b++) {
        double nMC = hMC_total->GetBinContent(b);
        double eMC = hMC_total->GetBinError(b);
        if (nMC > 0) {
            hRatio_MCerr->SetBinContent(b, 1.);       // centered at 1
            hRatio_MCerr->SetBinError(b, eMC / nMC);  // relative MC error
        } else {
            hRatio_MCerr->SetBinContent(b, 0.);
            hRatio_MCerr->SetBinError(b, 0.);
        }
    }
    hRatio_MCerr->SetFillColor(kGray + 1);
    hRatio_MCerr->SetFillStyle(3344);  // meshed
    hRatio_MCerr->SetLineColor(kGray + 1);
    hRatio_MCerr->SetLineWidth(0);
    hRatio_MCerr->SetMarkerSize(0);

    //////////////////////////////

    hRatio->Draw("E1");
    hRatio_MCerr->Draw("E2 SAME");
    TLine* line = new TLine(axMin, 1., axMax, 1.);
    line->SetLineColor(kBlue);
    line->SetLineWidth(2);
    line->Draw();
    hRatio->Draw("E1 SAME");

    gSystem->mkdir(plotdir, kTRUE);
    c->SaveAs(Form("%s/%s.pdf", plotdir.Data(), outname.Data()));
    printf("Saved: %s/%s.pdf\n", plotdir.Data(), outname.Data());

    delete hRatio;
    delete hRatio_MCerr;
    delete c;
    delete hMC_err;
    delete hMC_total;
    delete stack;
}

// ============================================================
// Main
// ============================================================
void plot_secondary_pions() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile* fMC = TFile::Open("hists_MC_0.5GeV_pionqe0p5.root");
    TFile* fData = TFile::Open("hists_Data_0.5GeV_pionqe0p5.root");
    if (!fMC || fMC->IsZombie()) {
        printf("ERROR: MC file\n");
        return;
    }
    if (!fData || fData->IsZombie()) {
        printf("ERROR: Data file\n");
        return;
    }

    const TString dir = "loose_charged_pion";
    const TString plotdir = "plots_sec";

    double globalScale = ComputeGlobalScale(fMC, fData);

    // ---- Secondary pion multiplicity ----
    DrawStackedPlot(fMC, fData,
                    dir,
                    "N_reco_pions",
                    "N(#pi^{#pm}_{reco})",
                    "plot_Nreco_pions",
                    -0.5, 5.5, 1, globalScale, plotdir);

    // ---- E_K^End for N(pi)=0 (absorption-dominated) ----
    DrawStackedPlot(fMC, fData,
                    dir,
                    "Beam_KE_end_0pi",
                    "E_{K}^{End} [MeV]",
                    "plot_KEend_0pi",
                    0., 600., 40, globalScale, plotdir);

    // ---- E_K^End for N(pi)>=1 (QE-dominated) ----
    DrawStackedPlot(fMC, fData,
                    dir,
                    "Beam_KE_end_least1pi",
                    "E_{K}^{End} [MeV]",
                    "plot_KEend_least1pi",
                    0., 600., 40, globalScale, plotdir);

    fMC->Close();
    fData->Close();
    printf("\nDone! Plots saved in %s/\n", plotdir.Data());
}
