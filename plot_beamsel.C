// ============================================================
// plot_beamsel.C
// Reproduces beam selection plots
// Usage: root -l -b -q plot_beamsel.C
// Output: PDF files in plots/ directory
// ============================================================

#include <string>
#include <vector>

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

// ============================================================
// Category definitions matching pi2 namespace in AnalyzerCore.h
// Suffix 0  = Data
// Suffix 1-13 = MC truth categories
// ============================================================
const bool DEBUG = true;
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

struct PlotDef {
    TString dir, var, xtitle, outname;
    double xmin, xmax;
    int rebin;
};

// ============================================================
// Helpers
// ============================================================
TH1D* GetHistData(TFile* f, TString dir, TString prefix) {
    TDirectory* d = (TDirectory*)f->Get(dir);
    if (!d) return nullptr;
    TH1D* h = (TH1D*)d->Get(prefix);
    if (!h) h = (TH1D*)d->Get(prefix + "_0");
    if (!h) return nullptr;
    h->SetDirectory(0);
    return h;
}

TH1D* GetHistMC(TFile* f, TString dir, TString prefix, int cat) {
    TDirectory* d = (TDirectory*)f->Get(dir);
    if (!d) return nullptr;
    TH1D* h = (TH1D*)d->Get(prefix + Form("_%d", cat));
    if (!h) return nullptr;
    h->SetDirectory(0);

    if (DEBUG) {
        if (h) printf(" DEBUG: %s_%d -> NBins=%d, Xmin=%.1f, Xmax=%.1f\n",
                      prefix.Data(), cat, h->GetNbinsX(),
                      h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    }

    return h;
}

// ============================================================
// Main plotting function
// xmin, xmax: set to -1 to use histogram defaults
// rebin: set to 1 to keep original binning
// ============================================================
void DrawPlot(TFile* fMC, TFile* fData,
              TString dir, TString varname,
              TString xtitle, TString outname,
              double xmin = -1, double xmax = -1,
              int rebin = 1,
              TString plotdir = "plots",
              TFile* fOut = nullptr) {
    vector<pair<TString, TString>> allCuts = {
        {"Beam_PID", "Beam PID"},
        {"Beam_scraper", "Scraper"},
        {"Beam_collhits", "Coll. hits"},
        {"Beam_recotrk", "Reco track"},
        {"Beam_endZ", "Z_{end} > 10 cm"},
        {"Beam_deltaXY", "#DeltaXY < 2#sigma"},
        {"Beam_chi2proton", "140 < #chi^{2}_{p} < 300"},
    };

    // Find index of current stage
    int currentStageIdx = -1;
    for (int i = 0; i < (int)allCuts.size(); i++) {
        if (allCuts[i].first == dir) {
            currentStageIdx = i;
            break;
        }
    }

    TString prefix = dir + "_" + varname;

    if (DEBUG) {
        std::string msg = std::string(":::: Processing ") + dir.Data() + " / " + varname.Data() + " ";
        if ((int)msg.size() < 119)
            msg.append(119 - msg.size(), ':');
        printf("%s\n", msg.c_str());
    }

    else
        printf(":::: Processing %s/%s\n", dir.Data(), varname.Data());

    // -- Get Data histogram
    TH1D* hData = GetHistData(fData, dir, prefix);
    if (!hData) {
        printf(" !!! WARNING !!! Data histogram not found for %s/%s, skipping\n",
               dir.Data(), varname.Data());
        return;
    }

    // -- Get MC histograms per category and build stack
    THStack* stack = new THStack("stack", "");
    TH1D* hMC_total = nullptr;
    std::vector<std::pair<int, TH1D*>> hMC_cat;

    for (int i = 1; i <= NCAT; i++) {
        TH1D* h = GetHistMC(fMC, dir, prefix, i);
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
        } else {
            hMC_total->Add(h);
        }
    }

    if (!hMC_total || hMC_total->Integral() == 0) {
        printf("!!! WARNING !!! No MC entries for %s/%s, skipping\n",
               dir.Data(), varname.Data());
        delete stack;
        return;
    }

    // -- Apply rebin to Data too
    if (rebin > 1) hData->Rebin(rebin);

    // -- Normalize MC to Data
    double dataIntegral, mcIntegral;
    if (xmax > xmin) {
        int b1 = hData->FindBin(xmin);
        int b2 = hData->FindBin(xmax) - 1;
        dataIntegral = hData->Integral(b1, b2);
        mcIntegral = hMC_total->Integral(b1, b2);
        if (mcIntegral != dataIntegral) {
            // Handle the case where integrals don't match
            printf(" !!! WARNING !!! Integrals don't match for %s/%s\n", dir.Data(), varname.Data());
            printf("  Data: %.2f, MC: %.2f\n", dataIntegral, mcIntegral);
        }

    } else {
        dataIntegral = hData->Integral();
        mcIntegral = hMC_total->Integral();
    }
    double scale = (mcIntegral > 0) ? dataIntegral / mcIntegral : 1.;
    for (auto& p : hMC_cat) p.second->Scale(scale);
    hMC_total->Scale(scale);
    double mcSumScaled = hMC_total->Integral();

    // -- Set axis range
    double axMin = (xmin >= 0) ? xmin : hData->GetXaxis()->GetXmin();
    double axMax = (xmax > 0) ? xmax : hData->GetXaxis()->GetXmax();

    // -- Style data
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.8);
    hData->SetLineColor(kBlack);

    // -- Canvas with two pads
    TCanvas* c = new TCanvas("c", "", 800, 700);
    TPad* pad1 = new TPad("pad1", "", 0., 0.28, 1., 1.);
    TPad* pad2 = new TPad("pad2", "", 0., 0., 1., 0.28);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.02);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.35);
    pad1->Draw();
    pad2->Draw();

    // -- Main pad
    pad1->cd();
    pad1->SetTicks(1, 1);

    stack->Draw("HIST");
    stack->GetXaxis()->SetRangeUser(axMin, axMax);
    stack->GetXaxis()->SetLabelSize(0);
    stack->GetXaxis()->SetTitleSize(0);
    stack->GetYaxis()->SetTitle("Events");
    stack->GetYaxis()->SetTitleSize(0.055);
    stack->GetYaxis()->SetTitleOffset(1.0);
    stack->GetYaxis()->SetLabelSize(0.05);

    // compute ymax in the visible range
    hData->GetXaxis()->SetRangeUser(axMin, axMax);
    hMC_total->GetXaxis()->SetRangeUser(axMin, axMax);
    double ymax = std::max(hData->GetMaximum(), hMC_total->GetMaximum());
    hData->GetXaxis()->UnZoom();
    hMC_total->GetXaxis()->UnZoom();

    // leave enough room for the legend inside the plot
    stack->SetMaximum(ymax * 2.2);
    stack->SetMinimum(0.);

    if (varname == "Beam_delta_X_spec_TPC") hData->GetXaxis()->SetRangeUser(-30, 30);
    hData->Draw("E1 SAME");

    // -- Legend inside plot, top area, 4 columns
    TLegend* leg = new TLegend(0.12, 0.6, 0.92, 0.88);
    leg->SetNColumns(4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);
    leg->SetMargin(0.12);

    // Add MC categories in stack order (bottom first = PiElas first)
    for (int i = 0; i < (int)hMC_cat.size(); i++) {
        int cat = hMC_cat[i].first;
        TH1D* h = hMC_cat[i].second;
        double integ = h->Integral();
        double pct = (mcSumScaled > 0) ? 100. * integ / mcSumScaled : 0.;
        leg->AddEntry(h, Form("#bf{%s %.1f, (%.1f %%)}", catName[cat], integ, pct), "f");
    }
    // MC Sum and Observed as last two entries
    leg->AddEntry((TObject*)nullptr, Form("MC Sum %.1f", mcSumScaled), "");
    leg->AddEntry(hData, Form("Observed %.0f", dataIntegral), "lep");
    leg->Draw();

    // -- Labels
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.052);
    lat.DrawLatex(0.12, 0.93, "ProtoDUNE-SP");
    lat.SetTextSize(0.044);
    lat.DrawLatex(0.35, 0.93, "#bf{#it{Preliminary}}");
    lat.SetTextSize(0.052);
    lat.DrawLatex(0.64, 0.93, "0.5 GeV/c Beam");

    TString line1 = "", line2 = "";
    int half = allCuts.size() / 2 + 1;
    for (int i = 0; i < (int)allCuts.size(); i++) {
        TString mark = (i <= currentStageIdx) ? "#color[416]{v} " : "#color[632]{x} ";
        TString entry = mark + allCuts[i].second + "   ";
        if (i < half)
            line1 += entry;
        else
            line2 += entry;
    }

    TPaveText* cutBox = new TPaveText(0.42, 0.52, 0.88, 0.6, "NDC");
    cutBox->SetFillColor(kWhite);
    cutBox->SetBorderSize(1);
    cutBox->SetTextSize(0.028);
    cutBox->SetTextAlign(12);
    cutBox->AddText(line1);
    cutBox->AddText(line2);
    cutBox->Draw();

    // -- Ratio pad
    pad2->cd();
    pad2->SetTicks(1, 1);

    // Build ratio manually bin-by-bin to preserve Poisson errors on data
    TH1D* hRatio = (TH1D*)hData->Clone("hRatio");
    hRatio->SetDirectory(0);
    hRatio->Reset();

    // Make sure both have same binning for division
    TH1D* hMC_forRatio = (TH1D*)hMC_total->Clone("hMC_forRatio");
    hMC_forRatio->SetDirectory(0);

    // If binning differs, rebin to match
    if (hData->GetNbinsX() != hMC_forRatio->GetNbinsX()) {
        double factor = hData->GetNbinsX() / hMC_forRatio->GetNbinsX();
        if (factor > 1) {
            hMC_forRatio->Rebin(factor);  // shouldn't happen
            if (DEBUG) printf(" !!! WARNING !!! Rebinning MC histogram by factor %f\n", factor);
        }
    }

    // Fill ratio bin by bin preserving data statistical errors
    for (int b = 1; b <= hRatio->GetNbinsX(); b++) {
        double nData = hData->GetBinContent(b);
        double nMC = hMC_total->GetBinContent(b);
        if (nMC > 0) {
            hRatio->SetBinContent(b, nData / nMC);
            // Error = sqrt(N_data) / N_MC  (Poisson error on data)
            hRatio->SetBinError(b, TMath::Sqrt(nData) / nMC);
        } else {
            hRatio->SetBinContent(b, 0.);
            hRatio->SetBinError(b, 0.);
        }
    }
    delete hMC_forRatio;

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
    hRatio->Draw("E1");

    TLine* line = new TLine(axMin, 1., axMax, 1.);
    line->SetLineColor(kBlue);
    line->SetLineWidth(2);
    line->Draw();

    if (fOut) {
        fOut->cd();
        c->Write(outname);
    }

    // -- Save
    gSystem->mkdir(plotdir, kTRUE);
    c->SaveAs(Form("%s/%s.pdf", plotdir.Data(), outname.Data()));
    printf("Saved: %s/%s.pdf\n\n", plotdir.Data(), outname.Data());

    delete hRatio;
    delete c;
    delete hMC_total;
    delete stack;
}

void PrintCutflow(TFile* fMC, const vector<PlotDef>& plots) {
    const int NCAT = 13;
    const char* catName[NCAT + 1] = {
        "Data",
        "PiElas", "PiRes", "PiQE", "PiDCEX", "PiABS", "PiCEX",
        "Muon", "misID:cosmic", "misID:p", "misID:pi", "misID:mu",
        "misID:e/gamma", "misID:other"};

    TString currentDir = "";
    for (auto& p : plots) {
        if (p.dir != currentDir) {
            currentDir = p.dir;
            printf("\n--- %s ---\n", currentDir.Data());
        }

        TString prefix = p.dir + "_" + p.var;
        double total = 0.;
        TString missing = "";
        for (int i = 1; i <= NCAT; i++) {
            TH1D* h = (TH1D*)fMC->Get(p.dir + "/" + prefix + Form("_%d", i));
            if (h)
                total += h->Integral();
            else
                missing += Form("%d(%s) ", i, catName[i]);
        }
        printf("  %-40s  total=%8.1f", p.var.Data(), total);
        if (missing.Length()) printf("  MISSING: %s", missing.Data());
        printf("\n");
    }
}

// ============================================================
// Main
// ============================================================
void plot_beamsel() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile* fMC = TFile::Open("hists_MC_0.5GeV_pionbeamsel.root");
    TFile* fData = TFile::Open("hists_data_0.5GeV_pionbeamsel.root");

    if (!fMC || fMC->IsZombie()) {
        printf("ERROR: cannot open MC file\n");
        return;
    }
    if (!fData || fData->IsZombie()) {
        printf("ERROR: cannot open Data file\n");
        return;
    }

    TFile* fOut = TFile::Open("plots_beamsel.root", "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        printf("ERROR: cannot open output ROOT file\n");
        return;
    }

    std::vector<PlotDef>
        plots = {
            // After PID cut
            {"Beam_PID", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "01_BeamPID_Pbeam", 350., 700., 8},
            {"Beam_PID", "Beam_endZ", "Z_{end}^{beam} [cm]", "01_BeamPID_endZ", -100., 400., 8},
            {"Beam_PID", "Beam_Z_dir_sign", "Z_{end}^{beam} [cm]", "01_BeamPID_Z_dir_sign", -100., 400., 8},
            {"Beam_PID", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "01_BeamPID_trkLenRatio", -1., 2., 12},
            {"Beam_PID", "Beam_reco_as_trk", "Beam reco. as Track", "01_BeamPID_recoastrk", -1., 1., 1},
            {"Beam_PID", "Beam_calo_size", "Has collection plane cluster", "01_Beam_PID_calosize", -1., 1., 1},
            {"Beam_PID", "Beam_chi2_proton", "#chi^{2}_{p}", "01_BeamPID_chi2p", 0., 400., 64},
            {"Beam_PID", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "01_BeamPID_deltaX", -30., 30., 16},
            {"Beam_PID", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "01_BeamPID_deltaY", -30., 30., 16},
            {"Beam_PID", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "01_BeamPID_deltaX_sigma", -5., 5., 16},
            {"Beam_PID", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "01_BeamPID_deltaY_sigma", -5., 5., 16},
            {"Beam_PID", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "01_BeamPID_KEff", 0., 600., 8},
            {"Beam_PID", "Beam_KELoss", "#DeltaE_{k} [MeV]", "01_BeamPID_KELoss", -60., 100., 25},
            // After beam scraper cut
            {"Beam_scraper", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "02_BeamScraper_Pbeam", 350., 650., 16},
            {"Beam_scraper", "Beam_endZ", "Z_{end}^{beam} [cm]", "02_BeamScraper_endZ", -100., 400., 8},
            {"Beam_scraper", "Beam_Z_dir_sign", "Z_{end}^{beam} [cm]", "02_BeamScraper_Z_dir_sign", -100., 400., 8},
            {"Beam_scraper", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "02_BeamScraper_trkLenRatio", 0., 2., 8},
            {"Beam_scraper", "Beam_reco_as_trk", "Beam reco. as Track", "02_BeamScraper_recoastrk", -1., 1., 1},
            {"Beam_scraper", "Beam_calo_size", "Has collection plane cluster", "02_Beam_Scraper_calosize", -1., 1., 1},
            {"Beam_scraper", "Beam_chi2_proton", "#chi^{2}_{p}", "02_BeamScraper_chi2p", 0., 400., 64},
            {"Beam_scraper", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "02_BeamScraper_deltaX", -30., 30., 16},
            {"Beam_scraper", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "02_BeamScraper_deltaY", -30., 30., 16},
            {"Beam_scraper", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "02_BeamScraper_deltaX_sigma", -5., 5., 16},
            {"Beam_scraper", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "02_BeamScraper_deltaY_sigma", -5., 5., 16},
            {"Beam_scraper", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "02_BeamScraper_KEff", 0., 600., 8},
            {"Beam_scraper", "Beam_KELoss", "#DeltaE_{k} [MeV]", "02_BeamScraper_KELoss", -60., 100., 25},
            // After collection hits cut
            {"Beam_collhits", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "03_BeamCollhits_Pbeam", 350., 650., 16},
            {"Beam_collhits", "Beam_endZ", "Z_{end}^{beam} [cm]", "03_BeamCollhits_endZ", -100., 400., 8},
            {"Beam_collhits", "Beam_Z_dir_sign", "Z_{end}^{beam} [cm]", "03_BeamCollhits_Z_dir_sign", -100., 400., 8},
            {"Beam_collhits", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "03_BeamCollhits_trkLenRatio", 0., 2., 8},
            {"Beam_collhits", "Beam_reco_as_trk", "Beam reco. as Track", "03_BeamCollhits_recoastrk", -1., 1., 1},
            {"Beam_collhits", "Beam_calo_size", "Has collection plane cluster", "03_BeamCollhits_calosize", -1., 1., 1},
            {"Beam_collhits", "Beam_chi2_proton", "#chi^{2}_{p}", "03_BeamCollhits_chi2p", 0., 400., 64},
            {"Beam_collhits", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "03_BeamCollhits_deltaX", -30., 30., 16},
            {"Beam_collhits", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "03_BeamCollhits_deltaY", -30., 30., 16},
            {"Beam_collhits", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "03_BeamCollhits_deltaX_sigma", -5., 5., 16},
            {"Beam_collhits", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "03_BeamCollhits_deltaY_sigma", -5., 5., 16},
            {"Beam_collhits", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "03_BeamCollhits_KEff", 0., 600., 8},
            {"Beam_collhits", "Beam_KELoss", "#DeltaE_{k} [MeV]", "03_BeamCollhits_KELoss", -60., 100., 25},
            // After reco track cut
            {"Beam_recotrk", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "04_BeamRecotrk_Pbeam", 350., 650., 16},
            {"Beam_recotrk", "Beam_endZ", "Z_{end}^{beam} [cm]", "04_BeamRecotrk_endZ", -100., 400., 8},
            {"Beam_recotrk", "Beam_Z_dir_sign", "Z_{end}^{beam} [cm]", "04_BeamRecotrk_Z_dir_sign", -100., 400., 8},
            {"Beam_recotrk", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "04_BeamRecotrk_trkLenRatio", 0., 2., 8},
            {"Beam_recotrk", "Beam_reco_as_trk", "Beam reco. as Track", "04_BeamRecotrk_recoastrk", -1., 1., 1},
            {"Beam_recotrk", "Beam_calo_size", "Has collection plane cluster", "04_BeamRecotrk_calosize", -1., 1., 1},
            {"Beam_recotrk", "Beam_chi2_proton", "#chi^{2}_{p}", "04_BeamRecotrk_chi2p", 0., 400., 64},
            {"Beam_recotrk", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "04_BeamRecotrk_deltaX", -30., 30., 16},
            {"Beam_recotrk", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "04_BeamRecotrk_deltaY", -30., 30., 16},
            {"Beam_recotrk", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "04_BeamRecotrk_deltaX_sigma", -5., 5., 16},
            {"Beam_recotrk", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "04_BeamRecotrk_deltaY_sigma", -5., 5., 16},
            {"Beam_recotrk", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "04_BeamRecotrk_KEff", 0., 600., 8},
            {"Beam_recotrk", "Beam_KELoss", "#DeltaE_{k} [MeV]", "04_BeamRecotrk_KELoss", -60., 100., 25},
            // After endZ cut
            {"Beam_endZ", "Beam_KE_end", "E_{K}^{End} [MeV]", "05_BeamEndZ_KEend", 0., 600., 16},
            {"Beam_endZ", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "05_BeamEndZ_Pbeam", 350., 650., 16},
            {"Beam_endZ", "Beam_endZ", "Z_{end}^{beam} [cm]", "05_BeamEndZ_endZ", -100., 400., 8},
            {"Beam_endZ", "Beam_Z_dir_sign", "Z_{end}^{beam} [cm]", "05_BeamEndZ_Z_dir_sign", -100., 400., 8},
            {"Beam_endZ", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "05_BeamEndZ_trkLenRatio", 0., 2., 8},
            {"Beam_endZ", "Beam_reco_as_trk", "Beam reco. as Track", "05_BeamEndZ_recoastrk", -1., 1., 1},
            {"Beam_endZ", "Beam_calo_size", "Has collection plane cluster", "05_BeamEndZ_calosize", -1., 1., 1},
            {"Beam_endZ", "Beam_chi2_proton", "#chi^{2}_{p}", "05_BeamEndZ_chi2p", 0., 400., 64},
            {"Beam_endZ", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "05_BeamEndZ_deltaX", -30., 30., 16},
            {"Beam_endZ", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "05_BeamEndZ_deltaY", -30., 30., 16},
            {"Beam_endZ", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "05_BeamEndZ_deltaX_sigma", -5., 5., 16},
            {"Beam_endZ", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "05_BeamEndZ_deltaY_sigma", -5., 5., 16},
            {"Beam_endZ", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "05_BeamEndZ_KEff", 0., 600., 8},
            {"Beam_endZ", "Beam_KELoss", "#DeltaE_{k} [MeV]", "05_BeamEndZ_KELoss", -60., 100., 25},
            // After deltaXY cut
            {"Beam_deltaXY", "Beam_KE_end", "E_{K}^{End} [MeV]", "06_BeamDeltaXY_KEend", 0., 600., 16},
            {"Beam_deltaXY", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "06_BeamDeltaXY_Pbeam", 350., 650., 16},
            {"Beam_deltaXY", "Beam_endZ", "Z_{end}^{beam} [cm]", "06_BeamDeltaXY_endZ", -100., 400., 8},
            {"Beam_deltaXY", "Beam_Z_dir_sign", "Z_{end}^{beam} [cm]", "06_BeamDeltaXY_Z_dir_sign", -100., 400., 8},
            {"Beam_deltaXY", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "06_BeamDeltaXY_trkLenRatio", 0., 2., 8},
            {"Beam_deltaXY", "Beam_reco_as_trk", "Beam reco. as Track", "06_BeamDeltaXY_recoastrk", -1., 1., 1},
            {"Beam_deltaXY", "Beam_calo_size", "Has collection plane cluster", "06_BeamDeltaXY_calosize", -1., 1., 1},
            {"Beam_deltaXY", "Beam_chi2_proton", "#chi^{2}_{p}", "06_BeamDeltaXY_chi2p", 0., 400., 64},
            {"Beam_deltaXY", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "06_BeamDeltaXY_deltaX", -30., 30., 16},
            {"Beam_deltaXY", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "06_BeamDeltaXY_deltaY", -30., 30., 16},
            {"Beam_deltaXY", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "06_BeamDeltaXY_deltaX_sigma", -5., 5., 16},
            {"Beam_deltaXY", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "06_BeamDeltaXY_deltaY_sigma", -5., 5., 16},
            {"Beam_deltaXY", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "06_BeamDeltaXY_KEff", 0., 600., 8},
            {"Beam_deltaXY", "Beam_KELoss", "#DeltaE_{k} [MeV]", "06_BeamDeltaXY_KELoss", -60., 100., 25},
            // After chi2 cut
            {"Beam_chi2proton", "Beam_KE_end", "E_{K}^{End} [MeV]", "07_BeamChi2p_KEend", 0., 600., 8},
            {"Beam_chi2proton", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "07_BeamChi2p_Pbeam", 350., 650., 16},
            {"Beam_chi2proton", "Beam_endZ", "Z_{end}^{beam} [cm]", "07_BeamChi2p_endZ", -100., 400., 8},
            {"Beam_chi2proton", "Beam_Z_dir_sign", "Z_{end}^{beam} [cm]", "07_BeamChi2p_Z_dir_sign", -100., 400., 8},
            {"Beam_chi2proton", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "07_BeamChi2p_trkLenRatio", 0., 2., 8},
            {"Beam_chi2proton", "Beam_reco_as_trk", "Beam reco. as Track", "07_BeamChi2p_recoastrk", -1., 1., 1},
            {"Beam_chi2proton", "Beam_calo_size", "Has collection plane cluster", "07_BeamChi2p_calosize", -1., 1., 1},
            {"Beam_chi2proton", "Beam_chi2_proton", "#chi^{2}_{p}", "07_BeamChi2p_chi2p", 0., 400., 64},
            {"Beam_chi2proton", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "07_BeamChi2p_deltaX", -30., 30., 16},
            {"Beam_chi2proton", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "07_BeamChi2p_deltaY", -30., 30., 16},
            {"Beam_chi2proton", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "07_BeamChi2p_deltaX_sigma", -5., 5., 16},
            {"Beam_chi2proton", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "07_BeamChi2p_deltaY_sigma", -5., 5., 16},
            {"Beam_chi2proton", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "07_BeamChi2p_KEff", 0., 600., 8},
            {"Beam_chi2proton", "Beam_KELoss", "#DeltaE_{k} [MeV]", "07_BeamChi2p_KELoss", -60., 100., 25},
        };

    for (auto& p : plots) {
        DrawPlot(fMC, fData, p.dir, p.var, p.xtitle, p.outname, p.xmin, p.xmax, p.rebin, "plots", fOut);
    }

    PrintCutflow(fMC, plots);

    fMC->Close();
    fData->Close();
    fOut->Close();
}
