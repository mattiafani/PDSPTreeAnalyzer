#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

// ============================================================
// Parameters — must match pionqe0p5.h constants
// ============================================================
const double slice_thickness_cm = 2.5; //5.;
const int n_slices = 100; //50;

// KE window (MeV) over which the truth-vs-Bertini closure is quantified.
// Trim edge slices (very high KE near entry, very low KE near stopping) where
// statistics are poor; widen/narrow as the sample warrants.
const double kCloseKEmin = 50.;
const double kCloseKEmax = 300.;

// KE window (MeV) actually drawn in the cross-section panels. Trims the
// end-of-range / first-slice bins where N_inc is tiny and the BB extrapolation
// is unreliable, so the panels read cleanly without hiding the physics region.
const double kPlotKEmin = 40.;
const double kPlotKEmax = 360.;

const double LAr_rho = 1.3954;
const double LAr_A = 39.948;
const double Avogadro = 6.02214e23;
const double cm2_to_mb = 1.e27;
const double xsec_norm = LAr_A / (LAr_rho * slice_thickness_cm * Avogadro) * cm2_to_mb;

// ============================================================
// Helpers (unchanged from prior version)
// ============================================================
TH1D* GetH1(TFile* f, TString dir, TString name) {
    TDirectory* d = (TDirectory*)f->Get(dir);
    if (!d) d = (TDirectory*)f;
    if (!d) return nullptr;
    TH1D* h = (TH1D*)d->Get(name);
    if (!h) return nullptr;
    h->SetDirectory(0);
    return h;
}

TH2D* GetH2(TFile* f, TString dir, TString name) {
    TDirectory* d = (TDirectory*)f->Get(dir);
    if (!d) d = (TDirectory*)f;
    if (!d) return nullptr;
    TH2D* h = (TH2D*)d->Get(name);
    if (!h) return nullptr;
    h->SetDirectory(0);
    return h;
}

TH1D* SumCats(TFile* f, TString dir, TString basename,
              int cat_lo = 1, int cat_hi = 6) {
    TH1D* htot = nullptr;
    for (int i = cat_lo; i <= cat_hi; i++) {
        TH1D* h = GetH1(f, dir, Form("%s_%d", basename.Data(), i));
        if (!h) continue;
        if (!htot) {
            htot = (TH1D*)h->Clone(basename + "_sum");
            htot->SetDirectory(0);
        } else
            htot->Add(h);
        delete h;
    }
    return htot;
}

TH2D* SumCats2D(TFile* f, TString dir, TString basename) {
    TH2D* htot = nullptr;
    for (int i = 1; i <= 13; i++) {
        TH2D* h = GetH2(f, dir, Form("%s_%d", basename.Data(), i));
        if (!h) continue;
        if (!htot) {
            htot = (TH2D*)h->Clone(basename + "_sum");
            htot->SetDirectory(0);
        } else
            htot->Add(h);
        delete h;
    }
    return htot;
}

// ============================================================
// KE per slice
// ============================================================
struct SliceKE {
    int id;
    double ke, ke_err;
};

vector<SliceKE> GetKEPerSlice(TFile* f, TString dir,
                              TString prefix, TString suffix,
                              bool use_mpv = false) {
    vector<SliceKE> out;
    for (int s = 0; s < n_slices; s++) {
        TString hname = Form("%s%d", prefix.Data(), s) + suffix;
        TH1D* h = GetH1(f, dir, hname);
        if (!h || h->GetEntries() < 5) {
            delete h;
            continue;
        }
        SliceKE sk;
        sk.id = s;
        if (!use_mpv) {
            sk.ke = h->GetMean();
            sk.ke_err = h->GetMeanError();
        } else {
            int mb = h->GetMaximumBin();
            sk.ke = h->GetBinCenter(mb);
            sk.ke_err = h->GetBinWidth(mb) / 2.;
        }
        out.push_back(sk);
        delete h;
    }
    return out;
}

vector<SliceKE> GetKEPerSliceSummed(TFile* f, TString dir, TString prefix,
                                    vector<int> cats, bool use_mpv = false) {
    vector<SliceKE> out;
    for (int s = 0; s < n_slices; s++) {
        TH1D* hsum = nullptr;
        for (int c : cats) {
            TString hname = Form("%s%d_%d", prefix.Data(), s, c);
            TH1D* h = GetH1(f, dir, hname);
            if (!h) continue;
            if (!hsum) {
                hsum = (TH1D*)h->Clone();
                hsum->SetDirectory(0);
            } else
                hsum->Add(h);
            delete h;
        }
        if (!hsum || hsum->GetEntries() < 5) {
            delete hsum;
            continue;
        }
        SliceKE sk;
        sk.id = s;
        if (!use_mpv) {
            sk.ke = hsum->GetMean();
            sk.ke_err = hsum->GetMeanError();
        } else {
            int mb = hsum->GetMaximumBin();
            sk.ke = hsum->GetBinCenter(mb);
            sk.ke_err = hsum->GetBinWidth(mb) / 2.;
        }
        out.push_back(sk);
        delete hsum;
    }
    return out;
}

// ============================================================
// Thin-slice cross-section as a TH1D keyed by slice-ID bin
// (same binning as the input N_inc/N_int histograms).
// sigma(slice) = xsec_norm * -ln(1 - N_int/N_inc)
// ============================================================
TH1D* XsecHist(TH1D* h_inc, TH1D* h_int, TString name) {
    if (!h_inc || !h_int) return nullptr;
    TH1D* h = (TH1D*)h_inc->Clone(name);
    h->SetDirectory(0);
    h->Reset();
    for (int b = 1; b <= h_inc->GetNbinsX(); b++) {
        double N_inc = h_inc->GetBinContent(b);
        double N_int = h_int->GetBinContent(b);
        if (N_inc <= 0.) continue;
        if (N_int < 0.) N_int = 0.;
        if (N_int > N_inc) N_int = N_inc;
        double ratio = N_int / N_inc;
        if (ratio >= 1.) ratio = 1. - 1.e-6;
        double sigma = xsec_norm * (-TMath::Log(1. - ratio));
        double d_ratio = (N_int > 0.) ? ratio * TMath::Sqrt(1. / N_int + 1. / N_inc) : 1. / N_inc;
        double d_sigma = xsec_norm * d_ratio / (1. - ratio);
        h->SetBinContent(b, sigma);
        h->SetBinError(b, d_sigma);
    }
    return h;
}

// Convert a per-slice cross-section TH1D to a TGraphErrors, placing each
// point at the KE recorded for that slice in ke_map.
TGraphErrors* HistToGraph(TH1D* h, const vector<SliceKE>* ke_map, TString name,
                          double kemin = -1.e9, double kemax = 1.e9) {
    if (!h) return nullptr;
    vector<double> x, y, ex, ey;
    for (int b = 1; b <= h->GetNbinsX(); b++) {
        double sigma = h->GetBinContent(b);
        if (sigma <= 0.) continue;
        int sliceID = (int)h->GetBinCenter(b);
        if (sliceID < 0 || sliceID >= n_slices) continue;
        double xval = sliceID, exval = 0.5;
        if (ke_map) {
            for (const auto& sk : *ke_map)
                if (sk.id == sliceID) {
                    xval = sk.ke;
                    exval = sk.ke_err;
                    break;
                }
        }
        if (xval < kemin || xval > kemax) continue;  // trim unreliable edge bins
        x.push_back(xval);
        y.push_back(sigma);
        ex.push_back(exval);
        ey.push_back(h->GetBinError(b));
    }
    if (x.empty()) return nullptr;
    TGraphErrors* g = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetName(name);
    return g;
}

// ============================================================
// Truth-vs-Bertini closure metric over a KE window.
//   chi2/N : sum of (sigma_true - sigma_G4)^2 / err^2 divided by N points
//   ratio  : stat-weighted mean of sigma_true/sigma_G4 with its error
// CAVEAT: this is not a rigorous goodness-of-fit. 
// The individual point errors are approximate as N_inc and
// N_int treated as independent; neighbouring slices are correlated through
// the cumulative N_inc.  
// A covariance-based version is the proper way.
// ============================================================
struct Closure {
    int n = 0;
    double chi2N = 0., ratio = 0., ratio_err = 0.;
};

Closure ComputeClosure(TGraphErrors* g, TGraph* g4,
                       double kemin, double kemax) {
    Closure c;
    if (!g || !g4) return c;
    double chi2 = 0., sumw = 0., sumwr = 0.;
    for (int i = 0; i < g->GetN(); i++) {
        double x, y;
        g->GetPoint(i, x, y);
        if (x < kemin || x > kemax || y <= 0.) continue;
        double yg = g4->Eval(x);
        if (yg <= 0.) continue;
        double ey = g->GetErrorY(i);
        if (ey <= 0.) continue;
        chi2 += (y - yg) * (y - yg) / (ey * ey);
        double r = y / yg, er = ey / yg, w = 1. / (er * er);
        sumw += w;
        sumwr += w * r;
        c.n++;
    }
    if (c.n > 0) {
        c.chi2N = chi2 / c.n;
        c.ratio = sumwr / sumw;
        c.ratio_err = 1. / TMath::Sqrt(sumw);
    }
    return c;
}

void DrawLabels(TString extra = "") {
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.044);
    lat.DrawLatex(0.12, 0.93, "ProtoDUNE-SP");
    lat.SetTextSize(0.038);
    lat.DrawLatex(0.43, 0.93, "#bf{#it{Preliminary}}");
    lat.SetTextSize(0.044);
    lat.DrawLatex(0.70, 0.93, "0.5 GeV/c #pi^{+}");
    if (extra.Length() > 0) {
        lat.SetTextSize(0.034);
        lat.DrawLatex(0.12, 0.85, extra);
    }
}

void SaveCanvas(TCanvas* c, TString name, TString plotdir) {
    gSystem->mkdir(plotdir, kTRUE);
    c->SaveAs(Form("%s/%s.pdf", plotdir.Data(), name.Data()));
    printf("Saved: %s/%s.pdf\n", plotdir.Data(), name.Data());
}

// ============================================================
// Load a Geant4 cross-section curve from data/v1/GEANT4_XS/histmap.txt
// Channels:
//   abs_KE         -> absorption
//   cex_KE         -> charge exchange
//   dcex_KE        -> double CEX
//   inel_KE        -> quasi-elastic
//   prod_KE        -> pion production
//   total_inel_KE  -> total inelastic (sum of all channels)
// ============================================================
TGraph* LoadG4Xsec(TString channel_histname,
                   TString filepath = "data/v1/GEANT4_XS/pion_xsec_1GeV.root") {
    TFile* fG4 = TFile::Open(filepath);
    if (!fG4 || fG4->IsZombie()) {
        printf("WARNING: cannot open Geant4 xsec file %s. Returning nullptr.\n",
               filepath.Data());
        return nullptr;
    }
    TGraph* g = (TGraph*)fG4->Get(channel_histname);
    if (!g) {
        printf("WARNING: graph '%s' not found in %s. Returning nullptr.\n",
               channel_histname.Data(), filepath.Data());
        fG4->Close();
        return nullptr;
    }
    TGraph* gclone = (TGraph*)g->Clone(Form("g4_%s", channel_histname.Data()));
    fG4->Close();

    double xmax = 0.;
    for (int i = 0; i < gclone->GetN(); i++) {
        double xx, yy;
        gclone->GetPoint(i, xx, yy);
        if (xx > xmax) xmax = xx;
    }
    if (xmax > 0. && xmax < 5.) {
        printf("[LoadG4Xsec] %s: detected GeV units (xmax=%.3f), converting to MeV.\n",
               channel_histname.Data(), xmax);
        for (int i = 0; i < gclone->GetN(); i++) {
            double xx, yy;
            gclone->GetPoint(i, xx, yy);
            gclone->SetPoint(i, xx * 1000., yy);
        }
    } else {
        printf("[LoadG4Xsec] %s: assuming MeV units (xmax=%.1f).\n",
               channel_histname.Data(), xmax);
    }
    return gclone;
}

TH1D* plot_xsec_common(TString plotdir,
                       std::vector<SliceKE>& ke_reco_out,
                       std::vector<SliceKE>& ke_true_out) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    printf("\n========================================================\n");
    printf("  Running plot_xsec_common (signal-independent plots)\n");
    printf("  Output: %s\n", plotdir.Data());
    printf("========================================================\n");

    TFile* fMC = TFile::Open("hists_MC_0.5GeV_pionqe0p5.root");
    if (!fMC || fMC->IsZombie()) {
        printf("ERROR: cannot open MC file\n");
        return nullptr;
    }

    TH1D* h_inc_mc_true = SumCats(fMC, "Xsec", "Xsec_N_inc_true");
    TH1D* h_inc_mc_sig = SumCats(fMC, "Xsec", "Xsec_N_inc_reco", 1, 6);

    TH1D* h_eff = nullptr;
    if (h_inc_mc_sig && h_inc_mc_true) {
        h_eff = (TH1D*)h_inc_mc_sig->Clone("h_eff");
        h_eff->SetDirectory(0);
        h_eff->Divide(h_inc_mc_true);
    }

    // here ke_reco_out is the average of the BB-predicted KE at this slice (see pionqe0p5.C's ke_reco_out),
    // given each event's KE_ff_reco
    ke_reco_out = GetKEPerSlice(fMC, "Xsec", "Xsec_KE_reco_slice", "", false);
    // ke_true_out = GetKEPerSlice(fMC, "Xsec", "Xsec_KE_true_slice", "_1", false);
    ke_true_out = GetKEPerSliceSummed(fMC, "Xsec", "Xsec_KE_true_slice", {1, 2, 3, 4, 5, 6}, false);
    printf("KE maps: reco=%d slices, true=%d slices\n",
           (int)ke_reco_out.size(), (int)ke_true_out.size());

    // ---- Plot: Efficiency ----
    {
        TCanvas* c = new TCanvas("c_eff_common", "", 700, 550);
        if (h_eff) {
            h_eff->SetMarkerStyle(20);
            h_eff->SetMarkerSize(0.9);
            h_eff->SetLineColor(kBlack);
            h_eff->GetXaxis()->SetTitle("Slice ID");
            h_eff->GetXaxis()->SetRangeUser(-0.5, 35.5);
            h_eff->GetYaxis()->SetTitle("Efficiency");
            h_eff->GetYaxis()->SetRangeUser(0., 1.5);
            h_eff->Draw("E1");
            TLine* l = new TLine(-0.5, 1., 25.5, 1.);
            l->SetLineStyle(2);
            l->SetLineColor(kBlue);
            l->Draw();
            DrawLabels("Eff = N_{sel}^{reco} / N_{all}^{true}");
        }
        SaveCanvas(c, "xsec_efficiency", plotdir);
        delete c;
    }

    // ---- Plot: Reco vs True KE per slice ----
    if (!ke_reco_out.empty() && !ke_true_out.empty()) {
        TCanvas* c = new TCanvas("c_ke_common", "", 1200, 500);
        c->Divide(2, 1);
        int np = std::min(ke_reco_out.size(), ke_true_out.size());
        vector<double> xv, yr, yt, exv, eyr, eyt, yd, eyd;
        for (int i = 0; i < (int)np; i++) {
            xv.push_back(ke_reco_out[i].id);
            yr.push_back(ke_reco_out[i].ke);
            eyr.push_back(ke_reco_out[i].ke_err);
            yt.push_back(ke_true_out[i].ke);
            eyt.push_back(ke_true_out[i].ke_err);
            exv.push_back(0.);
            yd.push_back(yr[i] - yt[i]);
            eyd.push_back(sqrt(eyr[i] * eyr[i] + eyt[i] * eyt[i]));
        }
        c->cd(1);
        TGraphErrors* gr = new TGraphErrors(np, xv.data(), yr.data(), exv.data(), eyr.data());
        TGraphErrors* gt = new TGraphErrors(np, xv.data(), yt.data(), exv.data(), eyt.data());
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kRed + 1);
        gr->SetLineColor(kRed + 1);
        gr->SetMarkerSize(0.8);
        gt->SetMarkerStyle(20);
        gt->SetMarkerColor(kBlack);
        gt->SetLineColor(kBlack);
        gt->SetMarkerSize(0.8);
        gr->GetXaxis()->SetTitle("Slice ID");
        gr->GetYaxis()->SetTitle("#pi^{+} KE (MeV)");
        gr->SetTitle("");
        gr->Draw("AP");
        gt->Draw("P SAME");
        TLegend* l1 = new TLegend(0.55, 0.72, 0.88, 0.84);
        l1->AddEntry(gr, "Reco KE (mean)", "p");
        l1->AddEntry(gt, "True KE (mean)", "p");
        l1->SetBorderSize(0);
        l1->Draw();
        DrawLabels();
        c->cd(2);
        TGraphErrors* gd = new TGraphErrors(np, xv.data(), yd.data(), exv.data(), eyd.data());
        gd->SetMarkerStyle(20);
        gd->SetMarkerSize(0.8);
        gd->SetLineColor(kBlack);
        gd->SetMarkerColor(kBlack);
        gd->GetXaxis()->SetTitle("Slice ID");
        gd->GetYaxis()->SetTitle("Reco KE #minus True KE (MeV)");
        gd->SetTitle("");
        gd->Draw("AP");
        TLine* l0 = new TLine(gd->GetXaxis()->GetXmin(), 0., gd->GetXaxis()->GetXmax(), 0.);
        l0->SetLineStyle(2);
        l0->SetLineColor(kBlue);
        l0->Draw();
        DrawLabels();
        SaveCanvas(c, "xsec_KE_reco_vs_true", plotdir);
        delete c;
    }

    delete h_inc_mc_true;
    delete h_inc_mc_sig;
    fMC->Close();
    return h_eff;
}

TGraphErrors* plot_xsec_one(TString signal_tag, TString plotdir,
                            TH1D* h_eff_in,
                            const std::vector<SliceKE>& ke_reco,
                            const std::vector<SliceKE>& ke_true) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TString sig_subscript;
    TString sig_long;
    TString g4_histname;
    if (signal_tag == "inel") {
        sig_subscript = "inel";
        sig_long = "Total inelastic";
        g4_histname = "total_inel_KE";
    } else if (signal_tag == "abslike") {
        sig_subscript = "abs";
        sig_long = "Absorption";
        g4_histname = "abs_KE";
    } else if (signal_tag == "cexlike") {
        sig_subscript = "cex";
        sig_long = "Charge exchange";
        g4_histname = "cex_KE";
    } else if (signal_tag == "otherlike") {
        sig_subscript = "QE";
        sig_long = "Quasi-elastic + other";
        g4_histname = "inel_KE";
    } else {
        printf("ERROR: unknown signal_tag '%s'.\n", signal_tag.Data());
        return nullptr;
    }
    printf("\n========================================================\n");
    printf("  Running plot_xsec_one for signal_tag=%s\n", signal_tag.Data());
    printf("  Output directory: %s\n", plotdir.Data());
    printf("========================================================\n");

    // MC-only stage: data is intentionally not loaded. We first establish that
    // the thin-slice extraction reproduces the Geant4 input from MC truth, and
    // we quantify the reco-level detector effect. Data returns once the fit
    // (which handles efficiency/smearing properly) is in place.
    TFile* fMC = TFile::Open("hists_MC_0.5GeV_pionqe0p5.root");
    if (!fMC || fMC->IsZombie()) {
        printf("ERROR: cannot open MC file\n");
        return nullptr;
    }

    TString pur_all_basename = (signal_tag == "inel")
                                   ? "Xsec_purity_all"
                                   : Form("Xsec_purity_all_%s", signal_tag.Data());

    TH1D* h_inc_mc_reco = SumCats(fMC, "Xsec", "Xsec_N_inc_reco");
    TH1D* h_int_mc_reco = SumCats(fMC, "Xsec",
                                  Form("Xsec_N_int_reco_%s", signal_tag.Data()));
    TH1D* h_inc_mc_true = SumCats(fMC, "Xsec", "Xsec_N_inc_true");
    TH1D* h_int_mc_true = SumCats(fMC, "Xsec",
                                  Form("Xsec_N_int_true_%s", signal_tag.Data()));
    TH1D* h_pur_all = SumCats(fMC, "Xsec", pur_all_basename);
    TH1D* h_pur_sig = SumCats(fMC, "Xsec",
                              Form("Xsec_purity_sig_%s", signal_tag.Data()),
                              1, 6);
    TH2D* h_migration = SumCats2D(fMC, "Xsec",
                                  Form("Xsec_migration_%s", signal_tag.Data()));

    // Per-channel true-KE map. Mapping of signal tag to pi_type categories:
    //   inel      -> 1..6 (all beam-matched pions)
    //   abslike   -> 5    (PiABS)
    //   cexlike   -> 6    (PiCEX)
    //   otherlike -> 1..4 (PiElas, PiRes, PiQE, PiDCEX)
    vector<int> true_cats_for_channel;
    if (signal_tag == "inel")
        true_cats_for_channel = {1, 2, 3, 4, 5, 6};
    else if (signal_tag == "abslike")
        true_cats_for_channel = {5};
    else if (signal_tag == "cexlike")
        true_cats_for_channel = {6};
    else
        true_cats_for_channel = {1, 2, 3, 4};
    vector<SliceKE> ke_true_channel = GetKEPerSliceSummed(
        fMC, "Xsec", "Xsec_KE_true_slice", true_cats_for_channel, false);

    // -- Purity
    TH1D* h_purity = nullptr;
    if (h_pur_all && h_pur_sig) {
        h_purity = (TH1D*)h_pur_sig->Clone("h_purity");
        h_purity->SetDirectory(0);
        h_purity->Divide(h_pur_all);
    }

    // -- Thin-slice cross sections from MC:
    //      sigma_MC,true : from truth N_inc/N_int -> must reproduce the Geant4
    //                      input curve. This is the closure test of the method.
    //      sigma_MC,reco : from reco  N_inc/N_int -> shows the size of the
    //                      detector effect (efficiency + slice migration + pool
    //                      leakage) that the forward fit will have to absorb.
    TH1D* hx_mc_true = XsecHist(h_inc_mc_true, h_int_mc_true, "hx_mc_true");
    TH1D* hx_mc_reco = XsecHist(h_inc_mc_reco, h_int_mc_reco, "hx_mc_reco");

    // -- Graphs: true sigma at the true-KE map, reco sigma at the reco-KE map.
    TGraphErrors* g_xsec_mc_true = HistToGraph(hx_mc_true, &ke_true_channel, "g_xsec_mc_true", kPlotKEmin, kPlotKEmax);
    TGraphErrors* g_xsec_mc_reco = HistToGraph(hx_mc_reco, &ke_reco, "g_xsec_mc_reco", kPlotKEmin, kPlotKEmax);

    // -- Geant4 reference curve
    TGraph* g_g4 = LoadG4Xsec(g4_histname);
    if (g_g4) {
        g_g4->SetLineColor(kRed + 1);
        g_g4->SetLineWidth(3);
        g_g4->SetMarkerSize(0);
    }

    // -- Truth-vs-Bertini closure (the headline MC validation).
    Closure clo = ComputeClosure(g_xsec_mc_true, g_g4, kCloseKEmin, kCloseKEmax);
    printf("\n[CLOSURE %s] KE[%.0f,%.0f] MeV: N=%d  chi2/N=%.2f  "
           "<sigma_true/G4>=%.3f +/- %.3f\n",
           signal_tag.Data(), kCloseKEmin, kCloseKEmax, clo.n,
           clo.chi2N, clo.ratio, clo.ratio_err);
    printf("    (indicative only: per-point errors are approximate and slices "
           "are correlated via the cumulative N_inc)\n");

    // ============================================================
    // Plot 1: N_inc and N_int composition
    // ============================================================
    {
        auto drawCounts = [&](TPad* pad, TH1D* hreco, TH1D* htrue,
                              TString xtitle, TString ytitle,
                              std::vector<TObject*>& owned) {
            if (!pad || !hreco || !htrue) return;
            pad->cd();
            pad->SetLeftMargin(0.14);
            pad->SetBottomMargin(0.14);

            double ymax = std::max(hreco->GetMaximum(), htrue->GetMaximum()) * 1.45;
            TH1D* hr = (TH1D*)hreco->Clone();
            hr->SetDirectory(0);
            owned.push_back(hr);
            hr->SetFillColor(kCyan + 1);
            hr->SetLineColor(kCyan + 3);
            hr->SetFillStyle(1001);
            hr->GetXaxis()->SetTitle(xtitle);
            hr->GetYaxis()->SetTitle(ytitle);
            hr->GetYaxis()->SetTitleOffset(1.3);
            hr->SetMaximum(ymax);
            hr->SetMinimum(0.);
            hr->GetXaxis()->SetRangeUser(-0.5, 35.5);
            hr->Draw("HIST");

            TH1D* ht = (TH1D*)htrue->Clone();
            ht->SetDirectory(0);
            owned.push_back(ht);
            ht->SetLineColor(kBlack);
            ht->SetLineWidth(2);
            ht->Draw("HIST SAME");

            TLegend* leg = new TLegend(0.52, 0.72, 0.88, 0.86);
            owned.push_back(leg);
            leg->AddEntry(hr, "MC reco", "f");
            leg->AddEntry(ht, "MC true", "l");
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.038);
            leg->Draw();
            DrawLabels();
        };

        TCanvas* c = new TCanvas("c_counts", "", 1200, 520);
        c->Divide(2, 1, 0.005, 0.01);
        std::vector<TObject*> owned;
        drawCounts((TPad*)c->GetPad(1),
                   h_inc_mc_reco, h_inc_mc_true,
                   "Slice ID", "Events (N_{inc})",
                   owned);
        drawCounts((TPad*)c->GetPad(2),
                   h_int_mc_reco, h_int_mc_true,
                   "Slice ID", "Events (N_{int})",
                   owned);
        c->Update();
        SaveCanvas(c, "xsec_Ninc_Nint_reco", plotdir);
        for (auto* obj : owned) delete obj;
        delete c;
    }

    // ============================================================
    // Plot 2: Purity
    // ============================================================
    {
        TCanvas* c = new TCanvas("c_purity", "", 700, 550);
        if (h_purity) {
            h_purity->SetMarkerStyle(20);
            h_purity->SetMarkerSize(0.9);
            h_purity->SetLineColor(kBlack);
            h_purity->GetXaxis()->SetTitle("Reco Slice ID");
            h_purity->GetXaxis()->SetRangeUser(-0.5, 35.5);
            h_purity->GetYaxis()->SetTitle("Purity");
            h_purity->GetYaxis()->SetRangeUser(0., 1.2);
            h_purity->Draw("E1");
            TLine* l = new TLine(-0.5, 1., 25.5, 1.);
            l->SetLineStyle(2);
            l->SetLineColor(kBlue);
            l->Draw();
            TString sig_label;
            if (signal_tag == "inel")
                sig_label = "Signal = kQE + kABS + kCEX + kOther";
            else if (signal_tag == "abslike")
                sig_label = "Signal = kABS";
            else if (signal_tag == "cexlike")
                sig_label = "Signal = kCEX";
            else
                sig_label = "Signal = kQE + kOther";
            DrawLabels(sig_label);
        }
        SaveCanvas(c, "xsec_purity", plotdir);
        delete c;
    }

    // ============================================================
    // Plot 3: Migration matrix
    // ============================================================
    {
        TCanvas* c = new TCanvas("c_mig", "", 700, 600);
        TH2D* hn = nullptr;
        if (h_migration) {
            hn = (TH2D*)h_migration->Clone("hn");
            hn->SetDirectory(0);
            for (int bx = 1; bx <= hn->GetNbinsX(); bx++) {
                double cs = 0.;
                for (int by = 1; by <= hn->GetNbinsY(); by++) cs += hn->GetBinContent(bx, by);
                if (cs > 0.)
                    for (int by = 1; by <= hn->GetNbinsY(); by++)
                        hn->SetBinContent(bx, by, hn->GetBinContent(bx, by) / cs);
            }
            gStyle->SetPalette(kBird);
            hn->GetXaxis()->SetTitle("Reco Slice ID");
            hn->GetYaxis()->SetTitle("True Slice ID");
            hn->GetXaxis()->SetRangeUser(-0.5, 35.5);
            hn->GetYaxis()->SetRangeUser(-0.5, 35.5);
            hn->GetZaxis()->SetRangeUser(0., 1.);
            hn->Draw("COLZ");
            DrawLabels();
        }
        SaveCanvas(c, "xsec_migration_matrix", plotdir);
        delete hn;
        delete c;
    }

    // ============================================================
    // Plot 4: Cross-section with MC/Geant4 ratio pad
    //   Top pad: Geant4 input line + MC truth (closure) + MC reco
    //   Bottom pad: MC/Geant4 ratios (truth ~1 = method closes)
    // ============================================================
    {
        TCanvas* c = new TCanvas("c_xsec", "", 800, 800);
        TPad* p1 = new TPad("p1xs", "", 0., 0.30, 1., 1.);
        TPad* p2 = new TPad("p2xs", "", 0., 0., 1., 0.30);
        p1->SetTopMargin(0.08);
        p1->SetBottomMargin(0.02);
        p1->SetLeftMargin(0.13);
        p2->SetTopMargin(0.03);
        p2->SetBottomMargin(0.30);
        p2->SetLeftMargin(0.13);
        p1->Draw();
        p2->Draw();

        // ---- TOP PAD ----
        p1->cd();

        if (g_xsec_mc_true) {
            g_xsec_mc_true->SetMarkerStyle(20);
            g_xsec_mc_true->SetMarkerSize(1.0);
            g_xsec_mc_true->SetLineColor(kBlack);
            g_xsec_mc_true->SetMarkerColor(kBlack);
            g_xsec_mc_true->SetLineWidth(2);
        }
        if (g_xsec_mc_reco) {
            g_xsec_mc_reco->SetMarkerStyle(24);
            g_xsec_mc_reco->SetMarkerSize(1.0);
            g_xsec_mc_reco->SetLineColor(kGreen + 2);
            g_xsec_mc_reco->SetMarkerColor(kGreen + 2);
            g_xsec_mc_reco->SetLineWidth(2);
            g_xsec_mc_reco->SetLineStyle(2);
        }

        TH1D* hf = new TH1D("hf_top", "", 100, 0., 400.);
        hf->SetDirectory(0);
        hf->GetYaxis()->SetTitle(
            Form("#sigma_{%s} (#pi^{+}-Ar) [mb]", sig_subscript.Data()));
        hf->GetYaxis()->SetRangeUser(0., 1500.);
        hf->GetYaxis()->SetTitleSize(0.055);
        hf->GetYaxis()->SetTitleOffset(1.05);
        hf->GetYaxis()->SetLabelSize(0.045);
        hf->GetXaxis()->SetLabelSize(0);
        hf->Draw("AXIS");

        if (g_g4) g_g4->Draw("L SAME");
        if (g_xsec_mc_reco) g_xsec_mc_reco->Draw("PZ SAME");
        if (g_xsec_mc_true) g_xsec_mc_true->Draw("PZ SAME");

        TLegend* leg = new TLegend(0.45, 0.62, 0.92, 0.84); 
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.036);
        if (g_g4) leg->AddEntry(g_g4, "Geant4 (Bertini) input", "l");
        if (g_xsec_mc_true) leg->AddEntry(g_xsec_mc_true, "MC truth (thin-slice closure)", "lp");
        if (g_xsec_mc_reco) leg->AddEntry(g_xsec_mc_reco, "MC reco (uncorrected)", "lp");
        leg->Draw();
        
        DrawLabels();
        TLatex lt;
        lt.SetNDC();
        lt.SetTextSize(0.034);
        lt.DrawLatex(0.45, 0.865,
                     Form("%s, t = %.1f cm slices", sig_long.Data(), slice_thickness_cm));

        if (clo.n > 0) {
            lt.SetTextSize(0.028);
            lt.DrawLatex(0.47, 0.575,
                         Form("Closure (KE %.0f#minus%.0f):", kCloseKEmin, kCloseKEmax));
            lt.DrawLatex(0.47, 0.530, Form("#chi^{2}/N = %.2f", clo.chi2N));
            lt.DrawLatex(0.47, 0.485,
                         Form("#LT#sigma_{true}/Bertini#GT = %.2f #pm %.2f",
                              clo.ratio, clo.ratio_err));
        }

        // ---- BOTTOM PAD: MC / Geant4 ratios ----
        p2->cd();
        TH1D* hf_bot = new TH1D("hf_bot", "", 100, 0., 400.);
        hf_bot->SetDirectory(0);
        hf_bot->GetXaxis()->SetTitle("#pi^{+} Kinetic Energy [MeV]");
        hf_bot->GetXaxis()->SetTitleSize(0.11);
        hf_bot->GetXaxis()->SetLabelSize(0.10);
        hf_bot->GetYaxis()->SetTitle("MC / Geant4");
        hf_bot->GetYaxis()->SetTitleSize(0.10);
        hf_bot->GetYaxis()->SetTitleOffset(0.50);
        hf_bot->GetYaxis()->SetLabelSize(0.09);
        hf_bot->GetYaxis()->SetNdivisions(505);
        hf_bot->GetYaxis()->SetRangeUser(0., 2.5);
        hf_bot->Draw("AXIS");

        auto makeRatio = [&](TGraphErrors* g, int color, int mstyle) -> TGraphErrors* {
            if (!g || !g_g4) return nullptr;
            vector<double> xr, yr, exr, eyr;
            for (int i = 0; i < g->GetN(); i++) {
                double xd, yd;
                g->GetPoint(i, xd, yd);
                if (yd <= 0.) continue;
                double yg = g_g4->Eval(xd);
                if (yg <= 0.) continue;
                xr.push_back(xd);
                yr.push_back(yd / yg);
                exr.push_back(0.);
                eyr.push_back(g->GetErrorY(i) / yg);
            }
            if (xr.empty()) return nullptr;
            TGraphErrors* gr = new TGraphErrors(xr.size(), xr.data(), yr.data(),
                                                exr.data(), eyr.data());
            gr->SetMarkerStyle(mstyle);
            gr->SetMarkerSize(1.0);
            gr->SetLineColor(color);
            gr->SetMarkerColor(color);
            gr->Draw("PZ SAME");
            return gr;
        };

        printf("\n=== MC closure summary (%s): sigma_true/G4 should be ~1 ===\n",
               signal_tag.Data());
        if (g_xsec_mc_true && g_g4) {
            for (int i = 0; i < g_xsec_mc_true->GetN(); i++) {
                double xd, yd;
                g_xsec_mc_true->GetPoint(i, xd, yd);
                if (yd <= 0.) continue;
                double yg = g_g4->Eval(xd);
                if (yg <= 0.) continue;
                printf("  KE=%-7.1f  s_true=%-8.1f  s_G4=%-8.1f  ratio=%.3f\n",
                       xd, yd, yg, yd / yg);
            }
        }

        TGraphErrors* gr_true = makeRatio(g_xsec_mc_true, kBlack, 20);
        TGraphErrors* gr_reco = makeRatio(g_xsec_mc_reco, kGreen + 2, 24);

        TLine* lref = new TLine(0., 1., 400., 1.);
        lref->SetLineStyle(2);
        lref->SetLineColor(kBlue);
        lref->SetLineWidth(2);
        lref->Draw();

        SaveCanvas(c,
                   Form("xsec_total_%s",
                        (signal_tag == "abs") ? "absorption" : "inelastic"),
                   plotdir);
        delete hf;
        delete hf_bot;
        delete gr_true;
        delete gr_reco;
        delete lref;
        delete leg;
        delete p1;
        delete p2;
        delete c;
    }

    // ============================================================
    // Plot 5: Reco vs True KE per slice, for this channel
    // ============================================================
    if (!ke_reco.empty() && !ke_true_channel.empty()) {
        TCanvas* c = new TCanvas("c_ke_channel", "", 1200, 500);
        c->Divide(2, 1);
        int np = std::min(ke_reco.size(), ke_true_channel.size());
        vector<double> xv, yr, yt, exv, eyr, eyt, yd, eyd;
        for (int i = 0; i < np; i++) {
            xv.push_back(ke_reco[i].id);
            yr.push_back(ke_reco[i].ke);
            eyr.push_back(ke_reco[i].ke_err);
            yt.push_back(ke_true_channel[i].ke);
            eyt.push_back(ke_true_channel[i].ke_err);
            exv.push_back(0.);
            yd.push_back(yr[i] - yt[i]);
            eyd.push_back(sqrt(eyr[i] * eyr[i] + eyt[i] * eyt[i]));
        }
        c->cd(1);
        TGraphErrors* gr = new TGraphErrors(np, xv.data(), yr.data(), exv.data(), eyr.data());
        TGraphErrors* gt = new TGraphErrors(np, xv.data(), yt.data(), exv.data(), eyt.data());
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kRed + 1);
        gr->SetLineColor(kRed + 1);
        gr->SetMarkerSize(0.8);
        gt->SetMarkerStyle(20);
        gt->SetMarkerColor(kBlack);
        gt->SetLineColor(kBlack);
        gt->SetMarkerSize(0.8);
        gr->GetXaxis()->SetTitle("Slice ID");
        gr->GetYaxis()->SetTitle("#pi^{+} KE (MeV)");
        gr->SetTitle("");
        gr->Draw("AP");
        gt->Draw("P SAME");
        TLegend* l1 = new TLegend(0.55, 0.72, 0.88, 0.84);
        l1->AddEntry(gr, "Reco KE (mean)", "p");
        l1->AddEntry(gt, Form("True KE (%s)", sig_long.Data()), "p");
        l1->SetBorderSize(0);
        l1->Draw();
        DrawLabels(sig_long);
        c->cd(2);
        TGraphErrors* gd = new TGraphErrors(np, xv.data(), yd.data(), exv.data(), eyd.data());
        gd->SetMarkerStyle(20);
        gd->SetMarkerSize(0.8);
        gd->SetLineColor(kBlack);
        gd->SetMarkerColor(kBlack);
        gd->GetXaxis()->SetTitle("Slice ID");
        gd->GetYaxis()->SetTitle("Reco KE #minus True KE (MeV)");
        gd->SetTitle("");
        gd->Draw("AP");
        TLine* l0 = new TLine(gd->GetXaxis()->GetXmin(), 0., gd->GetXaxis()->GetXmax(), 0.);
        l0->SetLineStyle(2);
        l0->SetLineColor(kBlue);
        l0->Draw();
        DrawLabels(sig_long);
        SaveCanvas(c, "xsec_KE_reco_vs_true", plotdir);
        delete c;
    }

    // Return the truth-level cross section (used by the sum-check in plot_xsec).
    TGraphErrors* g_return = g_xsec_mc_true ? (TGraphErrors*)g_xsec_mc_true->Clone(Form("g_return_%s", signal_tag.Data())) : nullptr;

    // -- Cleanup
    delete h_inc_mc_reco;
    delete h_int_mc_reco;
    delete h_inc_mc_true;
    delete h_int_mc_true;
    delete h_pur_all;
    delete h_pur_sig;
    delete h_purity;
    delete hx_mc_true;
    delete hx_mc_reco;
    delete h_migration;
    delete g_xsec_mc_true;
    delete g_xsec_mc_reco;
    delete g_g4;

    fMC->Close();
    printf("\nDone with %s. Plots saved in %s/\n", signal_tag.Data(), plotdir.Data());

    return g_return;
}

// ============================================================
// Usage:
//   root -l -b -q plot_xsec.C            // both (default)
//   root -l -b -q 'plot_xsec.C(0)'       // both, explicit
//   root -l -b -q 'plot_xsec.C(1)'       // absorption only
//   root -l -b -q 'plot_xsec.C(2)'       // total inelastic only
// ============================================================

void plot_xsec(int which = 0) {
    bool do_inel = (which == 0 || which == 1);
    bool do_abs = (which == 0 || which == 2);
    bool do_cex = (which == 0 || which == 3);
    bool do_other = (which == 0 || which == 4);

    std::vector<SliceKE> ke_reco, ke_true;
    TH1D* h_eff = plot_xsec_common("plots_xsec_common", ke_reco, ke_true);
    if (!h_eff) {
        printf("ERROR: Stage A failed; not running Stage B.\n");
        return;
    }

    TGraphErrors *g_inel = nullptr, *g_abs = nullptr, *g_cex = nullptr, *g_other = nullptr;
    if (do_inel) g_inel = plot_xsec_one("inel", "plots_xsec_inel", h_eff, ke_reco, ke_true);
    if (do_abs) g_abs = plot_xsec_one("abslike", "plots_xsec_abslike", h_eff, ke_reco, ke_true);
    if (do_cex) g_cex = plot_xsec_one("cexlike", "plots_xsec_cexlike", h_eff, ke_reco, ke_true);
    if (do_other) g_other = plot_xsec_one("otherlike", "plots_xsec_otherlike", h_eff, ke_reco, ke_true);

    // Sum-check: per-KE-bin comparison of inel vs (abs + cex + other).
    if (g_inel && g_abs && g_cex && g_other) {
        printf("\n=== Sum check: inel vs (abs + cex + other) ===\n");
        printf("%-10s  %-12s  %-12s  %-12s  %-12s  %-12s  %-8s\n",
               "KE[MeV]", "inel", "abs", "cex", "other", "abs+cex+oth", "ratio");
        int n = g_inel->GetN();
        for (int i = 0; i < n; i++) {
            double x_inel, y_inel;
            g_inel->GetPoint(i, x_inel, y_inel);
            double y_abs = g_abs->Eval(x_inel);
            double y_cex = g_cex->Eval(x_inel);
            double y_oth = g_other->Eval(x_inel);
            double sum = y_abs + y_cex + y_oth;
            double ratio = (y_inel > 0.) ? sum / y_inel : 0.;
            printf("%-10.1f  %-12.1f  %-12.1f  %-12.1f  %-12.1f  %-12.1f  %.3f\n",
                   x_inel, y_inel, y_abs, y_cex, y_oth, sum, ratio);
        }
    }

    delete g_inel;
    delete g_abs;
    delete g_cex;
    delete g_other;
    delete h_eff;
}
