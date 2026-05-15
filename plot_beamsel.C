// ============================================================
// plot_beamsel.C
// Reproduces beam selection plots
// Usage: root -l -b -q plot_beamsel.C
// Output: PDF files in plots/ directory
// ============================================================

#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
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
const bool DEBUG = false;
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

// ============================================================
// Number of sigma around the mean used to define the Gaussian
// fit range. Increase or decrease as needed.
// ============================================================
const double kFitRangeSigma = 6.0;

struct PlotDef {
    TString dir, var, xtitle, outname;
    double xmin, xmax;
    int nbins;
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
// Apply binning from a PlotDef to a histogram.
// ============================================================
TH1D* ApplyBinning(TH1D* h, const PlotDef& p, const char* newname) {
    if (p.nbins > 0 && p.xmax > p.xmin) {
        TH1D* hNew = new TH1D(newname, h->GetTitle(),
                              p.nbins, p.xmin, p.xmax);
        hNew->SetDirectory(0);
        for (int b = 1; b <= h->GetNbinsX(); b++) {
            double x = h->GetBinCenter(b);
            double val = h->GetBinContent(b);
            double err = h->GetBinError(b);
            int nb = hNew->FindBin(x);
            if (nb < 1 || nb > hNew->GetNbinsX()) continue;
            hNew->SetBinContent(nb, hNew->GetBinContent(nb) + val);
            hNew->SetBinError(nb, TMath::Sqrt(
                                      hNew->GetBinError(nb) * hNew->GetBinError(nb) + err * err));
        }
        hNew->SetFillColor(h->GetFillColor());
        hNew->SetLineColor(h->GetLineColor());
        hNew->SetLineWidth(h->GetLineWidth());
        hNew->SetMarkerStyle(h->GetMarkerStyle());
        hNew->SetMarkerSize(h->GetMarkerSize());
        delete h;
        return hNew;
    }
    return h;
}

// ============================================================
// Fit a single Gaussian to a histogram within kFitRangeSigma
// sigma of the histogram mean.  Draws the result on the current
// pad and returns the TF1 (owned by ROOT — do not delete).
//
// color    : line color for the fit curve
// fname    : unique name for the TF1 (must differ per call)
// Returns nullptr if the fit fails.
// ============================================================
TF1* FitAndDrawGaussian(TH1D* h, int color, const char* fname,
                        double xmin, double xmax) {
    // First-pass estimates from histogram moments restricted to
    // the display range.
    h->GetXaxis()->SetRangeUser(xmin, xmax);
    double mean0 = h->GetMean();
    double sigma0 = h->GetRMS();
    h->GetXaxis()->UnZoom();

    if (sigma0 <= 0.) {
        printf("FitAndDrawGaussian WARNING: sigma0=0 for %s, skipping fit\n", fname);
        return nullptr;
    }

    double fitLo = mean0 - kFitRangeSigma * sigma0;
    double fitHi = mean0 + kFitRangeSigma * sigma0;
    // Clamp to display range
    fitLo = std::max(fitLo, xmin);
    fitHi = std::min(fitHi, xmax);

    TF1* f = new TF1(fname, "gaus", fitLo, fitHi);
    f->SetParameters(h->GetMaximum(), mean0, sigma0);
    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(2);  // dashed, matching colleague's style

    // "Q" = quiet, "0" = do not draw automatically (we call Draw("SAME") below),
    // "R" = use range stored in TF1, "S" = return TFitResult
    TFitResultPtr r = h->Fit(f, "Q0R");
    if ((int)r != 0) {
        printf("FitAndDrawGaussian WARNING: fit did not converge for %s\n", fname);
        // Draw anyway with whatever ROOT found
    }

    // Draw over the full display range so the tail is visible
    f->SetRange(xmin, xmax);
    f->Draw("SAME");

    return f;
}

// ============================================================
// Compute a single MC->data scale factor.
// ============================================================
double ComputeGlobalScale(TFile* fMC, TFile* fData,
                          TString dir, TString var,
                          double xmin = -1., double xmax = -1.) {
    TString prefix = dir + "_" + var;

    TDirectory* dData = (TDirectory*)fData->Get(dir);
    if (!dData) {
        printf("ComputeGlobalScale WARNING: directory '%s' not found in data file\n", dir.Data());
        return 1.0;
    }
    TH1D* hData = (TH1D*)dData->Get(prefix);
    if (!hData) hData = (TH1D*)dData->Get(prefix + "_0");
    if (!hData) {
        printf("ComputeGlobalScale WARNING: histogram '%s' not found in data file\n", prefix.Data());
        return 1.0;
    }
    hData->SetDirectory(0);

    double dataIntegral, mcIntegral = 0.;
    if (xmax > xmin) {
        int b1 = hData->FindBin(xmin);
        int b2 = hData->FindBin(xmax) - 1;
        dataIntegral = hData->Integral(b1, b2);
    } else {
        dataIntegral = hData->Integral();
    }
    delete hData;

    TDirectory* dMC = (TDirectory*)fMC->Get(dir);
    if (!dMC) {
        printf("ComputeGlobalScale WARNING: directory '%s' not found in MC file\n", dir.Data());
        return 1.0;
    }
    for (int i = 1; i <= 13; i++) {
        TH1D* h = (TH1D*)dMC->Get(prefix + Form("_%d", i));
        if (!h) continue;
        h->SetDirectory(0);
        if (xmax > xmin) {
            int b1 = h->FindBin(xmin);
            int b2 = h->FindBin(xmax) - 1;
            mcIntegral += h->Integral(b1, b2);
        } else {
            mcIntegral += h->Integral();
        }
        delete h;
    }
    if (mcIntegral <= 0.) {
        printf("ComputeGlobalScale WARNING: MC integral is zero for '%s/%s'\n", dir.Data(), var.Data());
        return 1.0;
    }

    double scale = dataIntegral / mcIntegral;
    printf("ComputeGlobalScale: dir=%s  var=%s  data=%.1f  MC=%.1f  scale=%.4f\n",
           dir.Data(), var.Data(), dataIntegral, mcIntegral, scale);
    return scale;
}

// ============================================================
// Main plotting function
// ============================================================
void DrawPlot(TFile* fMC, TFile* fData,
              TString dir, TString varname,
              TString xtitle, TString outname,
              double xmin = -1, double xmax = -1,
              const PlotDef* pdef = nullptr,
              bool normalize = true,
              double fixedScale = -1.,
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

    int currentStageIdx = -1;
    for (int i = 0; i < (int)allCuts.size(); i++) {
        if (allCuts[i].first == dir) {
            currentStageIdx = i;
            break;
        }
    }

    TString prefix = dir + "_" + varname;

    if (DEBUG)
        printf(":::: Processing %s / %s\n", dir.Data(), varname.Data());
    else
        printf(":::: Processing %s/%s\n", dir.Data(), varname.Data());

    // Decide whether to overlay Gaussian fits for this variable.
    // Applied to the raw ΔX / ΔY distributions (not the /sigma variants).
    bool doGausFit = (varname == "Beam_delta_X_spec_TPC" ||
                      varname == "Beam_delta_Y_spec_TPC");

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
        if (pdef) h = ApplyBinning(h, *pdef, Form("hMC_%d", i));
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

    if (pdef) hData = ApplyBinning(hData, *pdef, "hData_binned");

    // -- Normalize MC to Data (or keep raw counts)
    double dataIntegral, mcIntegral, scale;
    if (xmax > xmin) {
        int b1 = hData->FindBin(xmin);
        int b2 = hData->FindBin(xmax) - 1;
        dataIntegral = hData->Integral(b1, b2);
        mcIntegral = hMC_total->Integral(b1, b2);
    } else {
        dataIntegral = hData->Integral();
        mcIntegral = hMC_total->Integral();
    }

    if (!normalize) {
        scale = 1.;
    } else if (fixedScale > 0.) {
        scale = fixedScale;
    } else {
        scale = (mcIntegral > 0) ? dataIntegral / mcIntegral : 1.;
        printf(" !!! WARNING !!! No fixed scale provided — falling back to per-plot normalisation for %s/%s\n",
               dir.Data(), varname.Data());
    }

    for (auto& p : hMC_cat) p.second->Scale(scale);
    hMC_total->Scale(scale);
    double mcSumScaled = hMC_total->Integral();

    // -- Axis range
    double axMin = (xmin >= 0) ? xmin : hData->GetXaxis()->GetXmin();
    double axMax = (xmax > 0) ? xmax : hData->GetXaxis()->GetXmax();
    // For signed variables the default xmin check fails; use explicit range
    // when both bounds are provided.
    if (xmin < 0 && xmax > 0) {
        axMin = xmin;
        axMax = xmax;
    }

    // -- Style data
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.8);
    hData->SetLineColor(kBlack);

    // -- Canvas, two pads
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
    stack->GetYaxis()->SetTitle(normalize ? "Events" : "Raw counts");
    stack->GetYaxis()->SetTitleSize(0.055);
    stack->GetYaxis()->SetTitleOffset(1.0);
    stack->GetYaxis()->SetLabelSize(0.05);

    hData->GetXaxis()->SetRangeUser(axMin, axMax);
    hMC_total->GetXaxis()->SetRangeUser(axMin, axMax);
    double ymax = std::max(hData->GetMaximum(), hMC_total->GetMaximum());
    hData->GetXaxis()->UnZoom();
    hMC_total->GetXaxis()->UnZoom();

    stack->SetMaximum(ymax * 2.2);
    stack->SetMinimum(0.);

    if (varname == "Beam_delta_X_spec_TPC")
        hData->GetXaxis()->SetRangeUser(-30, 30);
    hData->Draw("E1 SAME");

    // -- Gaussian fits (drawn on top of data points)
    TF1* fMCfit = nullptr;
    TF1* fDatafit = nullptr;
    if (doGausFit) {
        // Clone histograms for fitting so we don't disturb the display objects.
        TH1D* hMC_forFit = (TH1D*)hMC_total->Clone("hMC_forFit");
        TH1D* hData_forFit = (TH1D*)hData->Clone("hData_forFit");
        hMC_forFit->SetDirectory(0);
        hData_forFit->SetDirectory(0);

        fMCfit = FitAndDrawGaussian(hMC_forFit, kAzure + 2, Form("fMC_%s_%s", dir.Data(), varname.Data()), axMin, axMax);
        fDatafit = FitAndDrawGaussian(hData_forFit, kMagenta, Form("fData_%s_%s", dir.Data(), varname.Data()), axMin, axMax);

        delete hMC_forFit;
        delete hData_forFit;
    }

    // -- Legend, 4 columns
    TLegend* leg = new TLegend(0.12, 0.6, 0.92, 0.88);
    leg->SetNColumns(4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);
    leg->SetMargin(0.12);

    for (int i = 0; i < (int)hMC_cat.size(); i++) {
        int cat = hMC_cat[i].first;
        TH1D* h = hMC_cat[i].second;
        double integ = h->Integral();
        double pct = (mcSumScaled > 0) ? 100. * integ / mcSumScaled : 0.;
        if (normalize)
            leg->AddEntry(h, Form("#bf{%s %.1f, (%.1f %%)}", catName[cat], integ, pct), "f");
        else
            leg->AddEntry(h, Form("#bf{%s %.0f, (%.1f %%)}", catName[cat], integ, pct), "f");
    }
    if (normalize) {
        leg->AddEntry((TObject*)nullptr, Form("MC Sum %.1f", mcSumScaled), "");
        leg->AddEntry(hData, Form("Observed %.0f", dataIntegral), "lep");
    } else {
        leg->AddEntry((TObject*)nullptr, Form("MC Raw %.0f", mcSumScaled), "");
        leg->AddEntry(hData, Form("Data %.0f", dataIntegral), "lep");
    }

    leg->Draw();

    // Separate small legend on the left for the two Gaussian fit lines
    if (doGausFit && (fMCfit || fDatafit)) {
        TLegend* legFit = new TLegend(0.6, 0.42, 0.95, 0.52);
        legFit->SetBorderSize(0);
        legFit->SetFillStyle(0);
        legFit->SetTextSize(0.025);
        if (fMCfit) {
            legFit->AddEntry(fMCfit,
                             Form("#bf{MC fit:} #mu = %.3f, #sigma = %.3f",
                                  fMCfit->GetParameter(1), fMCfit->GetParameter(2)),
                             "l");
        }
        if (fDatafit) {
            legFit->AddEntry(fDatafit,
                             Form("#bf{Data fit:} #mu = %.3f, #sigma = %.3f",
                                  fDatafit->GetParameter(1), fDatafit->GetParameter(2)),
                             "l");
        }
        legFit->Draw();
    }

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

    TH1D* hRatio = (TH1D*)hData->Clone("hRatio");
    hRatio->SetDirectory(0);
    hRatio->Reset();

    TH1D* hMC_forRatio = (TH1D*)hMC_total->Clone("hMC_forRatio");
    hMC_forRatio->SetDirectory(0);

    if (hData->GetNbinsX() != hMC_forRatio->GetNbinsX()) {
        double factor = hData->GetNbinsX() / hMC_forRatio->GetNbinsX();
        if (factor > 1) hMC_forRatio->Rebin(factor);
    }

    for (int b = 1; b <= hRatio->GetNbinsX(); b++) {
        double nData = hData->GetBinContent(b);
        double nMC = hMC_total->GetBinContent(b);
        if (nMC > 0) {
            hRatio->SetBinContent(b, nData / nMC);
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
        TString modedir = normalize ? "normalized" : "raw";
        if (!fOut->GetDirectory(modedir)) fOut->mkdir(modedir);
        fOut->cd(modedir);
        c->Write(outname);
        fOut->cd();
    }

    TString pdfName = normalize ? outname : outname + "_raw";
    gSystem->mkdir(plotdir, kTRUE);
    c->SaveAs(Form("%s/%s.pdf", plotdir.Data(), pdfName.Data()));
    printf("Saved: %s/%s.pdf\n\n", plotdir.Data(), pdfName.Data());

    delete hRatio;
    delete c;
    delete hMC_total;
    delete stack;
}

// ============================================================
// Main
// ============================================================
void plot_beamsel() {
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

    TFile* fOut = TFile::Open("plots_beamsel.root", "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        printf("ERROR: cannot open output ROOT file\n");
        return;
    }

    std::vector<PlotDef> plots = {
        // After PID cut
        {"Beam_PID", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "01_BeamPID_Pbeam", 350., 700., 35},
        {"Beam_PID", "Beam_endZ", "Z_{end}^{beam} [cm]", "01_BeamPID_endZ", -100., 400., 50},
        {"Beam_PID", "Beam_Z_dir_sign", "Z dir. sign", "01_BeamPID_Z_dir_sign", -1., 1., 2},
        {"Beam_PID", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "01_BeamPID_trkLenRatio", 0., 2., 20},
        {"Beam_PID", "Beam_reco_as_trk", "Beam reco. as Track", "01_BeamPID_recoastrk", -1., 1., 2},
        {"Beam_PID", "Beam_calo_size", "Has collection plane cluster", "01_Beam_PID_calosize", -1., 1., 2},
        {"Beam_PID", "Beam_chi2_proton", "#chi^{2}_{p}", "01_BeamPID_chi2p", 0., 400., 40},
        {"Beam_PID", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "01_BeamPID_deltaX", -30., 30., 30},
        {"Beam_PID", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "01_BeamPID_deltaY", -30., 30., 30},
        {"Beam_PID", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "01_BeamPID_deltaX_sigma", -5., 5., 25},
        {"Beam_PID", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "01_BeamPID_deltaY_sigma", -5., 5., 25},
        {"Beam_PID", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "01_BeamPID_KEff", 250., 500., 25},
        {"Beam_PID", "Beam_KELoss", "#DeltaE_{k} [MeV]", "01_BeamPID_KELoss", -100., 100., 40},
        // After beam scraper cut
        {"Beam_scraper", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "02_BeamScraper_Pbeam", 350., 700., 35},
        {"Beam_scraper", "Beam_endZ", "Z_{end}^{beam} [cm]", "02_BeamScraper_endZ", -100., 400., 50},
        {"Beam_scraper", "Beam_Z_dir_sign", "Z dir. sign", "02_BeamScraper_Z_dir_sign", -1., 1., 2},
        {"Beam_scraper", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "02_BeamScraper_trkLenRatio", 0., 2., 20},
        {"Beam_scraper", "Beam_reco_as_trk", "Beam reco. as Track", "02_BeamScraper_recoastrk", -1., 1., 2},
        {"Beam_scraper", "Beam_calo_size", "Has collection plane cluster", "02_BeamScraper_calosize", -1., 1., 2},
        {"Beam_scraper", "Beam_chi2_proton", "#chi^{2}_{p}", "02_BeamScraper_chi2p", 0., 400., 40},
        {"Beam_scraper", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "02_BeamScraper_deltaX", -30., 30., 30},
        {"Beam_scraper", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "02_BeamScraper_deltaY", -30., 30., 30},
        {"Beam_scraper", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "02_BeamScraper_deltaX_sigma", -5., 5., 25},
        {"Beam_scraper", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "02_BeamScraper_deltaY_sigma", -5., 5., 25},
        {"Beam_scraper", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "02_BeamScraper_KEff", 250., 500., 25},
        {"Beam_scraper", "Beam_KELoss", "#DeltaE_{k} [MeV]", "02_BeamScraper_KELoss", -100., 100., 40},
        // After collection hits cut
        {"Beam_collhits", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "03_BeamCollhits_Pbeam", 350., 700., 35},
        {"Beam_collhits", "Beam_endZ", "Z_{end}^{beam} [cm]", "03_BeamCollhits_endZ", -100., 400., 50},
        {"Beam_collhits", "Beam_Z_dir_sign", "Z dir. sign", "03_BeamCollhits_Z_dir_sign", -1., 1., 2},
        {"Beam_collhits", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "03_BeamCollhits_trkLenRatio", 0., 2., 20},
        {"Beam_collhits", "Beam_reco_as_trk", "Beam reco. as Track", "03_BeamCollhits_recoastrk", -1., 1., 2},
        {"Beam_collhits", "Beam_calo_size", "Has collection plane cluster", "03_BeamCollhits_calosize", -1., 1., 2},
        {"Beam_collhits", "Beam_chi2_proton", "#chi^{2}_{p}", "03_BeamCollhits_chi2p", 0., 400., 40},
        {"Beam_collhits", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "03_BeamCollhits_deltaX", -30., 30., 30},
        {"Beam_collhits", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "03_BeamCollhits_deltaY", -30., 30., 30},
        {"Beam_collhits", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "03_BeamCollhits_deltaX_sigma", -5., 5., 25},
        {"Beam_collhits", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "03_BeamCollhits_deltaY_sigma", -5., 5., 25},
        {"Beam_collhits", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "03_BeamCollhits_KEff", 250., 500., 25},
        {"Beam_collhits", "Beam_KELoss", "#DeltaE_{k} [MeV]", "03_BeamCollhits_KELoss", -100., 100., 40},
        // After reco track cut
        {"Beam_recotrk", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "04_BeamRecotrk_Pbeam", 350., 700., 35},
        {"Beam_recotrk", "Beam_endZ", "Z_{end}^{beam} [cm]", "04_BeamRecotrk_endZ", -100., 400., 50},
        {"Beam_recotrk", "Beam_Z_dir_sign", "Z dir. sign", "04_BeamRecotrk_Z_dir_sign", -1., 1., 2},
        {"Beam_recotrk", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "04_BeamRecotrk_trkLenRatio", 0., 2., 20},
        {"Beam_recotrk", "Beam_reco_as_trk", "Beam reco. as Track", "04_BeamRecotrk_recoastrk", -1., 1., 2},
        {"Beam_recotrk", "Beam_calo_size", "Has collection plane cluster", "04_BeamRecotrk_calosize", -1., 1., 2},
        {"Beam_recotrk", "Beam_chi2_proton", "#chi^{2}_{p}", "04_BeamRecotrk_chi2p", 0., 400., 40},
        {"Beam_recotrk", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "04_BeamRecotrk_deltaX", -30., 30., 30},
        {"Beam_recotrk", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "04_BeamRecotrk_deltaY", -30., 30., 30},
        {"Beam_recotrk", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "04_BeamRecotrk_deltaX_sigma", -5., 5., 25},
        {"Beam_recotrk", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "04_BeamRecotrk_deltaY_sigma", -5., 5., 25},
        {"Beam_recotrk", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "04_BeamRecotrk_KEff", 250., 500., 25},
        {"Beam_recotrk", "Beam_KELoss", "#DeltaE_{k} [MeV]", "04_BeamRecotrk_KELoss", -100., 100., 40},
        // After endZ cut
        {"Beam_endZ", "Beam_KE_end", "E_{K}^{End} [MeV]", "05_BeamEndZ_KEend", 0., 500., 50},
        {"Beam_endZ", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "05_BeamEndZ_Pbeam", 350., 700., 35},
        {"Beam_endZ", "Beam_endZ", "Z_{end}^{beam} [cm]", "05_BeamEndZ_endZ", -100., 400., 50},
        {"Beam_endZ", "Beam_Z_dir_sign", "Z dir. sign", "05_BeamEndZ_Z_dir_sign", -1., 1., 2},
        {"Beam_endZ", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "05_BeamEndZ_trkLenRatio", 0., 2., 20},
        {"Beam_endZ", "Beam_reco_as_trk", "Beam reco. as Track", "05_BeamEndZ_recoastrk", -1., 1., 2},
        {"Beam_endZ", "Beam_calo_size", "Has collection plane cluster", "05_BeamEndZ_calosize", -1., 1., 2},
        {"Beam_endZ", "Beam_chi2_proton", "#chi^{2}_{p}", "05_BeamEndZ_chi2p", 0., 400., 40},
        {"Beam_endZ", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "05_BeamEndZ_deltaX", -30., 30., 30},
        {"Beam_endZ", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "05_BeamEndZ_deltaY", -30., 30., 30},
        {"Beam_endZ", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "05_BeamEndZ_deltaX_sigma", -5., 5., 25},
        {"Beam_endZ", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "05_BeamEndZ_deltaY_sigma", -5., 5., 25},
        {"Beam_endZ", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "05_BeamEndZ_KEff", 250., 500., 25},
        {"Beam_endZ", "Beam_KELoss", "#DeltaE_{k} [MeV]", "05_BeamEndZ_KELoss", -100., 100., 40},
        // After deltaXY cut
        {"Beam_deltaXY", "Beam_KE_end", "E_{K}^{End} [MeV]", "06_BeamDeltaXY_KEend", 0., 500., 50},
        {"Beam_deltaXY", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "06_BeamDeltaXY_Pbeam", 350., 700., 35},
        {"Beam_deltaXY", "Beam_endZ", "Z_{end}^{beam} [cm]", "06_BeamDeltaXY_endZ", -100., 400., 50},
        {"Beam_deltaXY", "Beam_Z_dir_sign", "Z dir. sign", "06_BeamDeltaXY_Z_dir_sign", -1., 1., 2},
        {"Beam_deltaXY", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "06_BeamDeltaXY_trkLenRatio", 0., 2., 20},
        {"Beam_deltaXY", "Beam_reco_as_trk", "Beam reco. as Track", "06_BeamDeltaXY_recoastrk", -1., 1., 2},
        {"Beam_deltaXY", "Beam_calo_size", "Has collection plane cluster", "06_BeamDeltaXY_calosize", -1., 1., 2},
        {"Beam_deltaXY", "Beam_chi2_proton", "#chi^{2}_{p}", "06_BeamDeltaXY_chi2p", 0., 400., 40},
        {"Beam_deltaXY", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "06_BeamDeltaXY_deltaX", -30., 30., 30},
        {"Beam_deltaXY", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "06_BeamDeltaXY_deltaY", -30., 30., 30},
        {"Beam_deltaXY", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "06_BeamDeltaXY_deltaX_sigma", -5., 5., 25},
        {"Beam_deltaXY", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "06_BeamDeltaXY_deltaY_sigma", -5., 5., 25},
        {"Beam_deltaXY", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "06_BeamDeltaXY_KEff", 250., 500., 25},
        {"Beam_deltaXY", "Beam_KELoss", "#DeltaE_{k} [MeV]", "06_BeamDeltaXY_KELoss", -100., 100., 40},
        // After chi2 cut
        {"Beam_chi2proton", "Beam_KE_end", "E_{K}^{End} [MeV]", "07_BeamChi2p_KEend", 0., 500., 50},
        {"Beam_chi2proton", "Beam_P_beam_inst", "P_{spec.} [MeV/c]", "07_BeamChi2p_Pbeam", 350., 700., 35},
        {"Beam_chi2proton", "Beam_endZ", "Z_{end}^{beam} [cm]", "07_BeamChi2p_endZ", -100., 400., 50},
        {"Beam_chi2proton", "Beam_Z_dir_sign", "Z dir. sign", "07_BeamChi2p_Z_dir_sign", -1., 1., 2},
        {"Beam_chi2proton", "Beam_trk_len_ratio", "L_{Beam track}/L_{Exp.}", "07_BeamChi2p_trkLenRatio", 0., 2., 20},
        {"Beam_chi2proton", "Beam_reco_as_trk", "Beam reco. as Track", "07_BeamChi2p_recoastrk", -1., 1., 2},
        {"Beam_chi2proton", "Beam_calo_size", "Has collection plane cluster", "07_BeamChi2p_calosize", -1., 1., 2},
        {"Beam_chi2proton", "Beam_chi2_proton", "#chi^{2}_{p}", "07_BeamChi2p_chi2p", 0., 400., 40},
        {"Beam_chi2proton", "Beam_delta_X_spec_TPC", "#DeltaX(spec., z=10cm) [cm]", "07_BeamChi2p_deltaX", -30., 30., 30},
        {"Beam_chi2proton", "Beam_delta_Y_spec_TPC", "#DeltaY(spec., z=10cm) [cm]", "07_BeamChi2p_deltaY", -30., 30., 30},
        {"Beam_chi2proton", "Beam_delta_X_spec_TPC_over_sigma", "#DeltaX/#sigma", "07_BeamChi2p_deltaX_sigma", -5., 5., 25},
        {"Beam_chi2proton", "Beam_delta_Y_spec_TPC_over_sigma", "#DeltaY/#sigma", "07_BeamChi2p_deltaY_sigma", -5., 5., 25},
        {"Beam_chi2proton", "Beam_KE_ff", "E_{K}(Z=10cm) [MeV]", "07_BeamChi2p_KEff", 250., 500., 25},
        {"Beam_chi2proton", "Beam_KELoss", "#DeltaE_{k} [MeV]", "07_BeamChi2p_KELoss", -100., 100., 40},
    };

    for (auto& p : plots) {
        double scale = ComputeGlobalScale(fMC, fData, "Beam_scraper", p.var, p.xmin, p.xmax);
        // Normalized (MC scaled to data)
        DrawPlot(fMC, fData, p.dir, p.var, p.xtitle, p.outname,
                 p.xmin, p.xmax, &p, /*normalize=*/true, scale, "plots", fOut);
        // Raw counts
        DrawPlot(fMC, fData, p.dir, p.var, p.xtitle, p.outname,
                 p.xmin, p.xmax, &p, /*normalize=*/false, 1.0, "plots", fOut);
    }

    fMC->Close();
    fData->Close();
    fOut->Close();
}
