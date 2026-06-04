// ============================================================
// plot_composition_vs_cuts.C
//
// For each variable in plot_beamsel.C, reads the MC histograms
// at each cut stage and shows how the category composition (%)
// evolves through the cuts.
//
// Normalization convention:
//   Single global MC->data scale, computed once at Beam_scraper
//   using Beam_P_beam_inst over its full filled range. This matches
//   the normalization point used by plot_beamsel.C's plot legends
//   and reveals the growing data/MC discrepancy at later stages
//   (data marker drops below MC total).
//
// Usage:  root -l -b -q plot_composition_vs_cuts.C
// ============================================================

#include <map>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

// ============================================================
// Category definitions (mirror of plot_beamsel.C)
// ============================================================
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

// Cut stages in order
const int NSTAGES = 7;
struct Stage {
    TString dir;
    TString label;
};
Stage stages[NSTAGES] = {
    {"Beam_PID", "Beam PID"},
    {"Beam_scraper", "Scraper"},
    {"Beam_collhits", "Coll. hits"},
    {"Beam_recotrk", "Reco trk"},
    {"Beam_endZ", "Z_{end}"},
    {"Beam_deltaXY", "#DeltaXY"},
    {"Beam_chi2proton", "#chi^{2}_{p}"},
};

// ============================================================
// Variable definition
// ============================================================
struct VarDef {
    TString suffix;
    TString branch;
    TString xtitle;
    double compMin, compMax;  // -1,-1 = full range
};

// ============================================================
// Compute the single global MC->data scale at the scraper stage,
// using Beam_P_beam_inst over its full filled range.
// ============================================================
double ComputeGlobalScale(TFile* fMC, TFile* fData) {
    TDirectory* dMC = (TDirectory*)fMC->Get("Beam_scraper");
    TDirectory* dData = (TDirectory*)fData->Get("Beam_scraper");
    if (!dMC || !dData) {
        printf(
            "WARNING: Beam_scraper directory not found in MC or Data; "
            "falling back to scale=1.0\n");
        return 1.;
    }
    const TString prefix = "Beam_scraper_Beam_P_beam_inst";

    // Data
    TH1D* hData = (TH1D*)dData->Get(prefix);
    if (!hData) hData = (TH1D*)dData->Get(prefix + "_0");
    double dataTotal = (hData) ? hData->Integral() : -1.;

    // MC summed over categories
    double mcTotal = 0.;
    for (int i = 1; i <= NCAT; i++) {
        TH1D* h = (TH1D*)dMC->Get(prefix + Form("_%d", i));
        if (h) mcTotal += h->Integral();
    }

    if (dataTotal <= 0. || mcTotal <= 0.) {
        printf(
            "WARNING: cannot compute scale (mcTotal=%.1f, dataTotal=%.1f); "
            "falling back to 1.0\n",
            mcTotal, dataTotal);
        return 1.;
    }
    double scale = dataTotal / mcTotal;
    printf(
        "\n==> Global MC->data scale (Beam_scraper / Beam_P_beam_inst): "
        "%.4f (Data=%.1f / MC=%.1f)\n\n",
        scale, dataTotal, mcTotal);
    return scale;
}

// ============================================================
// Per-stage category yields and data yield.
// fixedScale: the precomputed MC->data scale (no per-stage rescaling).
// compMin/compMax: integration window (-1,-1 = full filled range).
// Returns false if MC at this stage is empty / missing.
// ============================================================
bool GetStageComposition(TFile* fMC, TFile* fData,
                         TString dir, TString branch,
                         double compMin, double compMax,
                         bool normalize, double fixedScale,
                         double yields[NCAT + 1],
                         double& dataYield) {
    TString prefix = dir + "_" + branch;
    bool hasComp = (compMin >= 0. && compMax > compMin);
    double scale = (normalize) ? fixedScale : 1.;

    // --- data
    TDirectory* dData = (TDirectory*)fData->Get(dir);
    TH1D* hData = nullptr;
    if (dData) {
        hData = (TH1D*)dData->Get(prefix);
        if (!hData) hData = (TH1D*)dData->Get(prefix + "_0");
        if (hData) hData->SetDirectory(0);
    }
    if (hData) {
        if (hasComp) {
            int b1 = hData->FindBin(compMin);
            int b2 = hData->FindBin(compMax) - 1;
            dataYield = hData->Integral(b1, b2);
        } else {
            dataYield = hData->Integral();
        }
    } else {
        dataYield = -1.;
    }

    // --- MC per category
    TDirectory* dMC = (TDirectory*)fMC->Get(dir);
    if (!dMC) {
        delete hData;
        return false;
    }

    TH1D* hCat[NCAT + 1] = {};
    bool anyFound = false;
    double mcSumIntegralCheck = 0.;
    for (int i = 1; i <= NCAT; i++) {
        TH1D* h = (TH1D*)dMC->Get(prefix + Form("_%d", i));
        if (!h) continue;
        h->SetDirectory(0);
        hCat[i] = h;
        anyFound = true;
        mcSumIntegralCheck += h->Integral();
    }
    if (!anyFound || mcSumIntegralCheck == 0.) {
        delete hData;
        for (int i = 1; i <= NCAT; i++) delete hCat[i];
        return false;
    }

    yields[0] = 0.;
    for (int i = 1; i <= NCAT; i++) {
        if (hCat[i]) {
            double y;
            if (hasComp) {
                int b1 = hCat[i]->FindBin(compMin);
                int b2 = hCat[i]->FindBin(compMax) - 1;
                y = hCat[i]->Integral(b1, b2);
            } else {
                y = hCat[i]->Integral();
            }
            yields[i] = y * scale;
            yields[0] += yields[i];
        } else {
            yields[i] = 0.;
        }
    }

    delete hData;
    for (int i = 1; i <= NCAT; i++) delete hCat[i];
    return true;
}

// ============================================================
// Label color helper
// ============================================================
int LabelColor(int fillColor) {
    switch (fillColor) {
        case kBlue + 2:
        case kGreen + 2:
        case kMagenta + 1:
        case kBlue - 7:
        case kGray + 2:
            return kWhite;
        default:
            return kBlack;
    }
}

// ============================================================
// Draw one composition-vs-cuts plot
// ============================================================
void DrawCompositionPlot(TFile* fMC, TFile* fData,
                         const VarDef& vd,
                         bool normalize, double globalScale,
                         TString plotdir = "plots") {
    printf("==== Composition plot: %s (normalize=%s) ====\n",
           vd.suffix.Data(), normalize ? "true" : "false");

    double yields[NSTAGES][NCAT + 1];
    double dataYield[NSTAGES];
    bool ok[NSTAGES];

    for (int s = 0; s < NSTAGES; s++) {
        ok[s] = GetStageComposition(fMC, fData,
                                    stages[s].dir, vd.branch,
                                    vd.compMin, vd.compMax,
                                    normalize, globalScale,
                                    yields[s], dataYield[s]);
        if (!ok[s])
            printf("  WARNING: stage %s not found for branch %s\n",
                   stages[s].dir.Data(), vd.branch.Data());
    }

    TH1D* hAbs[NCAT + 1] = {};
    TH1D* hFrac[NCAT + 1] = {};

    for (int i = 1; i <= NCAT; i++) {
        hAbs[i] = new TH1D(Form("hAbs_%d", i), "", NSTAGES, 0, NSTAGES);
        hFrac[i] = new TH1D(Form("hFrac_%d", i), "", NSTAGES, 0, NSTAGES);
        hAbs[i]->SetDirectory(0);
        hFrac[i]->SetDirectory(0);
        for (int s = 0; s < NSTAGES; s++) {
            hAbs[i]->GetXaxis()->SetBinLabel(s + 1, stages[s].label);
            hFrac[i]->GetXaxis()->SetBinLabel(s + 1, stages[s].label);
            if (!ok[s]) continue;
            hAbs[i]->SetBinContent(s + 1, yields[s][i]);
            double total = yields[s][0];
            hFrac[i]->SetBinContent(s + 1, (total > 0) ? 100. * yields[s][i] / total : 0.);
        }
        hAbs[i]->SetFillColor(catColor[i]);
        hAbs[i]->SetLineColor(kWhite);
        hAbs[i]->SetLineWidth(1);
        hFrac[i]->SetFillColor(catColor[i]);
        hFrac[i]->SetLineColor(kWhite);
        hFrac[i]->SetLineWidth(1);
    }

    // Data
    TGraphErrors* grData = new TGraphErrors(NSTAGES);
    for (int s = 0; s < NSTAGES; s++) {
        grData->SetPoint(s, s + 0.5, ok[s] ? dataYield[s] : 0.);
        grData->SetPointError(s, 0.5, (ok[s] && dataYield[s] > 0) ? TMath::Sqrt(dataYield[s]) : 0.);
    }
    grData->SetMarkerStyle(20);
    grData->SetMarkerSize(1.2);
    grData->SetMarkerColor(kBlack);
    grData->SetLineColor(kBlack);
    grData->SetLineWidth(2);

    // Canvas
    TCanvas* c = new TCanvas("cComp", "", 1050, 900);
    c->SetFillColor(0);

    TPad* pad1 = new TPad("pad1", "", 0., 0.43, 1., 1.);
    TPad* pad2 = new TPad("pad2", "", 0., 0., 1., 0.43);
    pad1->SetTopMargin(0.10);
    pad1->SetBottomMargin(0.018);
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.11);
    pad1->SetFillColor(0);
    pad2->SetTopMargin(0.018);
    pad2->SetBottomMargin(0.28);
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.11);
    pad2->SetFillColor(0);
    pad1->Draw();
    pad2->Draw();

    // ------- top pad
    pad1->cd();
    pad1->SetTicks(1, 1);

    THStack* stkAbs = new THStack("stkAbs", "");
    for (int i = 1; i <= NCAT; i++) stkAbs->Add(hAbs[i]);
    stkAbs->Draw("HIST");

    stkAbs->GetXaxis()->SetLabelSize(0);
    stkAbs->GetYaxis()->SetTitle(normalize ? "MC Events (norm. to data)" : "Raw counts");
    stkAbs->GetYaxis()->SetTitleSize(0.060);
    stkAbs->GetYaxis()->SetTitleOffset(0.80);
    stkAbs->GetYaxis()->SetLabelSize(0.052);
    stkAbs->GetYaxis()->SetNdivisions(506);

    double ymax = 0.;
    for (int s = 0; s < NSTAGES; s++) {
        if (!ok[s]) continue;
        ymax = TMath::Max(ymax, yields[s][0]);
        if (dataYield[s] > 0.) ymax = TMath::Max(ymax, dataYield[s]);
    }
    stkAbs->SetMaximum(ymax * 2.35);
    stkAbs->SetMinimum(0.);

    for (int s = 0; s < NSTAGES; s += 2) {
        TBox* shd = new TBox(s, 0., s + 1, ymax * 2.35);
        shd->SetFillColorAlpha(kGray, 0.15);
        shd->SetLineColorAlpha(kGray, 0.0);
        shd->Draw();
    }
    stkAbs->Draw("HIST SAME");
    grData->Draw("P SAME");

    for (int s = 1; s < NSTAGES; s++) {
        TLine* vl = new TLine(s, 0., s, ymax * 2.35);
        vl->SetLineStyle(3);
        vl->SetLineColor(kGray + 1);
        vl->SetLineWidth(1);
        vl->Draw();
    }

    // Legend
    TLegend* leg = new TLegend(0.11, 0.62, 0.96, 0.90);
    leg->SetNColumns(4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.050);
    leg->SetMargin(0.14);
    for (int i = 1; i <= NCAT; i++)
        leg->AddEntry(hAbs[i], catName[i], "f");
    leg->AddEntry(grData, "Data", "lep");
    leg->Draw();

    // Header
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.062);
    lat.DrawLatex(0.13, 0.912, "#bf{ProtoDUNE-SP}");
    lat.SetTextSize(0.052);
    lat.DrawLatex(0.42, 0.912, "#it{Preliminary}");
    lat.SetTextSize(0.060);
    lat.DrawLatex(0.67, 0.912, "0.5 GeV/c Beam");

    // Annotation
    bool hasCompWin = (vd.compMin >= 0. && vd.compMax > vd.compMin);
    lat.SetTextFont(62);
    lat.SetTextSize(0.042);
    lat.SetTextColor(kBlack);
    lat.DrawLatex(0.15, 0.590, vd.branch.Data());
    lat.SetTextFont(42);
    lat.SetTextSize(0.036);
    lat.SetTextColor(kGray + 2);
    if (hasCompWin)
        lat.DrawLatex(0.15, 0.553,
                      Form("Comp. window: %s #in [%.4g, %.4g]",
                           vd.xtitle.Data(), vd.compMin, vd.compMax));
    else
        lat.DrawLatex(0.15, 0.553,
                      Form("Comp. window: full range of %s", vd.xtitle.Data()));
    lat.SetTextColor(kBlack);

    // Data/MC ratio label, only shown in normalized plots if the ratio is physically meaningful
    if (normalize) {
        TLatex lr;
        lr.SetTextSize(0.030);
        lr.SetTextAlign(22);
        lr.SetTextFont(42);
        for (int s = 0; s < NSTAGES; s++) {
            if (!ok[s] || yields[s][0] <= 0 || dataYield[s] <= 0) continue;
            double r = dataYield[s] / yields[s][0];
            double yTxt = TMath::Max(yields[s][0], dataYield[s]) * 1.06;
            // color: black for ~1.0, red for noticeable departures
            int col = (TMath::Abs(r - 1.) < 0.05)   ? kBlack
                      : (TMath::Abs(r - 1.) < 0.15) ? kGray + 2
                                                    : kRed + 1;
            lr.SetTextColor(col);
            lr.DrawLatex(s + 0.5, yTxt, Form("Obs/Pred=%.2f", r));
        }
    }

    // ------- bottom pad
    pad2->cd();
    pad2->SetTicks(1, 1);

    THStack* stkFrac = new THStack("stkFrac", "");
    for (int i = 1; i <= NCAT; i++) stkFrac->Add(hFrac[i]);
    stkFrac->Draw("HIST");

    TAxis* xax = stkFrac->GetXaxis();
    xax->SetLabelSize(0.090);
    xax->SetLabelOffset(0.02);
    xax->LabelsOption("u");
    xax->SetTickLength(0.03);

    stkFrac->GetYaxis()->SetTitle("Composition [%]");
    stkFrac->GetYaxis()->SetTitleSize(0.082);
    stkFrac->GetYaxis()->SetTitleOffset(0.58);
    stkFrac->GetYaxis()->SetLabelSize(0.075);
    stkFrac->GetYaxis()->SetNdivisions(505);
    stkFrac->SetMinimum(0.);
    stkFrac->SetMaximum(100.);

    for (int s = 0; s < NSTAGES; s += 2) {
        TBox* shd = new TBox(s, 0., s + 1, 100.);
        shd->SetFillColorAlpha(kGray, 0.15);
        shd->SetLineColorAlpha(kGray, 0.0);
        shd->Draw();
    }
    stkFrac->Draw("HIST SAME");

    const double kMinLabelFrac = 6.0;
    const double kMinBandHeight = 5.5;
    for (int s = 0; s < NSTAGES; s++) {
        if (!ok[s]) continue;
        double cumFrac = 0.;
        for (int i = 1; i <= NCAT; i++) {
            double frac = (yields[s][0] > 0) ? 100. * yields[s][i] / yields[s][0] : 0.;
            double bandBot = cumFrac;
            cumFrac += frac;
            if (frac < kMinLabelFrac || frac < kMinBandHeight) continue;
            TLatex lf;
            lf.SetTextSize(0.055);
            lf.SetTextAlign(22);
            lf.SetTextColor(LabelColor(catColor[i]));
            lf.SetTextFont(62);
            lf.DrawLatex(s + 0.5, bandBot + frac / 2., Form("%.0f%%", frac));
        }
    }

    for (int s = 1; s < NSTAGES; s++) {
        TLine* vl = new TLine(s, 0., s, 100.);
        vl->SetLineStyle(3);
        vl->SetLineColor(kGray + 1);
        vl->SetLineWidth(1);
        vl->Draw();
    }

    gSystem->mkdir(plotdir, kTRUE);
    TString pdfName = Form("composition_%s_%s", vd.suffix.Data(),
                           normalize ? "norm" : "raw");
    TString outpath = Form("%s/%s.pdf", plotdir.Data(), pdfName.Data());
    c->SaveAs(outpath);
    printf("Saved: %s\n\n", outpath.Data());

    delete c;
    delete grData;
    delete stkAbs;
    delete stkFrac;
    for (int i = 1; i <= NCAT; i++) {
        delete hAbs[i];
        delete hFrac[i];
    }
}

// ============================================================
// Main
// ============================================================
void plot_composition_vs_cuts() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile* fMC = TFile::Open("hists_MC_0.5GeV_pionbeamsel.root");
    TFile* fData = TFile::Open("hists_Data_0.5GeV_pionbeamsel.root");
    if (!fMC || fMC->IsZombie()) {
        printf("ERROR: cannot open MC file\n");
        return;
    }
    if (!fData || fData->IsZombie()) {
        printf("ERROR: cannot open Data file\n");
        return;
    }

    // Single global scale used for every (stage, variable) call.
    double globalScale = ComputeGlobalScale(fMC, fData);

    std::vector<VarDef> vars = {
        {"recoastrk", "Beam_reco_as_trk", "Beam reco. as Track", 0., 1.},
        {"calosize", "Beam_calo_size", "Has calo cluster", 0., 1.},
        {"Pbeam", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", 350., 650.},
        {"endZ", "Beam_endZ", "Z_{end}^{beam} [cm]", 10., 400.},
        {"Z_dir_sign", "Beam_Z_dir_sign", "Z dir. sign", -1., -1.},
        {"trkLenRatio", "Beam_trk_len_ratio", "L_{Beam}/L_{Exp.}", 0., 2.},
        {"chi2p", "Beam_chi2_proton", "#chi^{2}_{p}", 140., 300.},
        {"deltaX", "Beam_delta_X_spec_TPC", "#DeltaX [cm]", -10., 10.},
        {"deltaY", "Beam_delta_Y_spec_TPC", "#DeltaY [cm]", -10., 10.},
        {"deltaX_sigma", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", -2., 2.},
        {"deltaY_sigma", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", -2., 2.},
        {"KEff", "Beam_KE_ff", "E_{K}(z=10cm) [MeV]", 0., 600.},
        {"KELoss", "Beam_KELoss", "#DeltaE_{k} [MeV]", -60., 100.},
        {"KEend", "Beam_KE_end", "E_{K}^{End} [MeV]", 0., 600.},
    };

    for (auto& v : vars) {
        DrawCompositionPlot(fMC, fData, v, /*normalize=*/true, globalScale, "plots");
        DrawCompositionPlot(fMC, fData, v, /*normalize=*/false, globalScale, "plots");
    }

    fMC->Close();
    fData->Close();
}
