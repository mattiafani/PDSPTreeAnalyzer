// ============================================================
// plot_xsec.C
// Thin-slice cross-section plots for ProtoDUNE 0.5 GeV/c pi+ analysis.
// Restructured to run common (signal-independent) plots once and
// per-tag (signal-dependent) plots once per signal definition.
//
// Usage:
//   root -l -b -q plot_xsec.C            // both signal definitions (default)
//   root -l -b -q 'plot_xsec.C(0)'       // both, explicit
//   root -l -b -q 'plot_xsec.C(1)'       // absorption only
//   root -l -b -q 'plot_xsec.C(2)'       // total inelastic only
//
// Input: hists_MC_0.5GeV_pionqe0p5.root
//        hists_Data_0.5GeV_pionqe0p5.root
//        data/v1/GEANT4_XS/pion_xsec_1GeV.root
//
// Output:
//   plots_xsec_common/       (signal-independent: efficiency, KE map)
//   plots_xsec_abs/          (absorption signal)
//   plots_xsec_inel/         (total-inelastic signal, B2-inclusive)
// ============================================================

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
const double slice_thickness_cm = 5.;
const int n_slices = 50;

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

TGraphErrors* ComputeXsec(TH1D* h_inc, TH1D* h_int,
                          TH1D* h_eff = nullptr,
                          const vector<SliceKE>* ke_map = nullptr,
                          TString label = "") {
    if (!h_inc || !h_int) return nullptr;
    vector<double> x, y, ex, ey;
    for (int b = 1; b <= h_inc->GetNbinsX(); b++) {
        double sliceID = h_inc->GetBinCenter(b);
        if (sliceID < 0 || sliceID >= n_slices) continue;
        double N_inc = h_inc->GetBinContent(b);
        double N_int = h_int->GetBinContent(b);
        if (N_inc <= 0.) continue;
        if (N_int < 0.) N_int = 0.;
        if (N_int > N_inc) N_int = N_inc;
        double eff = 1.;
        if (h_eff) {
            eff = h_eff->GetBinContent(b);
            if (eff <= 0.) continue;
        }
        double ratio = N_int / (N_inc * eff);
        if (ratio >= 1.) ratio = 1. - 1.e-6;
        double sigma = xsec_norm * (-TMath::Log(1. - ratio));
        double d_ratio = (N_int > 0.) ? ratio * TMath::Sqrt(1. / N_int + 1. / N_inc) : 1. / N_inc;
        double d_sigma = xsec_norm * d_ratio / (1. - ratio);
        double xval = sliceID, exval = 0.5;
        if (ke_map) {
            for (const auto& sk : *ke_map)
                if (sk.id == (int)sliceID) {
                    xval = sk.ke;
                    exval = sk.ke_err;
                    break;
                }
        }
        x.push_back(xval);
        y.push_back(sigma);
        ex.push_back(exval);
        ey.push_back(d_sigma);
    }
    if (x.empty()) return nullptr;
    TGraphErrors* g = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetName("g_xsec_" + label);
    return g;
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
// Load a Geant4 cross-section curve from the tabulated file.
// Channels (from data/v1/GEANT4_XS/histmap.txt):
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

    ke_reco_out = GetKEPerSlice(fMC, "Xsec", "Xsec_KE_reco_slice", "", false);
    ke_true_out = GetKEPerSlice(fMC, "Xsec", "Xsec_KE_true_slice", "_1", false);
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

void plot_xsec_one(TString signal_tag, TString plotdir,
                   TH1D* h_eff_in,
                   const std::vector<SliceKE>& ke_reco,
                   const std::vector<SliceKE>& ke_true) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TString sig_subscript;
    TString sig_long;
    TString g4_histname;
    if (signal_tag == "abs") {
        sig_subscript = "abs";
        sig_long = "Absorption";
        g4_histname = "abs_KE";
    } else if (signal_tag == "inel") {
        sig_subscript = "inel";
        sig_long = "Total inelastic";
        g4_histname = "total_inel_KE";
    } else {
        printf("ERROR: unknown signal_tag '%s'. Use 'abs' or 'inel'.\n",
               signal_tag.Data());
        return;
    }
    printf("\n========================================================\n");
    printf("  Running plot_xsec_one for signal_tag=%s\n", signal_tag.Data());
    printf("  Output directory: %s\n", plotdir.Data());
    printf("========================================================\n");

    TFile* fMC = TFile::Open("hists_MC_0.5GeV_pionqe0p5.root");
    TFile* fData = TFile::Open("hists_Data_0.5GeV_pionqe0p5.root");
    if (!fMC || fMC->IsZombie()) {
        printf("ERROR: cannot open MC file\n");
        return;
    }
    if (!fData || fData->IsZombie()) {
        printf("ERROR: cannot open Data file\n");
        return;
    }

    TH1D* h_inc_data = GetH1(fData, "Xsec", "Xsec_N_inc_reco_0");
    TH1D* h_int_data = GetH1(fData, "Xsec",
                             Form("Xsec_N_int_reco_%s_0", signal_tag.Data()));
    TH1D* h_inc_mc_reco = SumCats(fMC, "Xsec", "Xsec_N_inc_reco");
    TH1D* h_int_mc_reco = SumCats(fMC, "Xsec",
                                  Form("Xsec_N_int_reco_%s", signal_tag.Data()));
    TH1D* h_inc_mc_true = SumCats(fMC, "Xsec", "Xsec_N_inc_true");
    TH1D* h_int_mc_true = SumCats(fMC, "Xsec",
                                  Form("Xsec_N_int_true_%s", signal_tag.Data()));
    TH1D* h_pur_all = SumCats(fMC, "Xsec", "Xsec_purity_all");
    TH1D* h_pur_sig = SumCats(fMC, "Xsec",
                              Form("Xsec_purity_signal_%s", signal_tag.Data()),
                              1, 6);
    TH2D* h_migration = SumCats2D(fMC, "Xsec",
                                  Form("Xsec_migration_%s", signal_tag.Data()));

    // -- Purity
    TH1D* h_purity = nullptr;
    if (h_pur_all && h_pur_sig) {
        h_purity = (TH1D*)h_pur_sig->Clone("h_purity");
        h_purity->SetDirectory(0);
        h_purity->Divide(h_pur_all);
    }

    // -- Apply purity to data; efficiency comes from Stage A.
    TH1D* h_inc_data_pur = nullptr;
    TH1D* h_int_data_pur = nullptr;
    if (h_inc_data && h_purity) {
        h_inc_data_pur = (TH1D*)h_inc_data->Clone("h_inc_data_pur");
        h_inc_data_pur->SetDirectory(0);
        h_inc_data_pur->Multiply(h_purity);
        h_int_data_pur = (TH1D*)h_int_data->Clone("h_int_data_pur");
        h_int_data_pur->SetDirectory(0);
        h_int_data_pur->Multiply(h_purity);
    }

    // -- Cross-section graphs
    TGraphErrors* g_xsec_mc_reco =
        ComputeXsec(h_inc_mc_reco, h_int_mc_reco, nullptr, &ke_reco, "mc_reco_raw");
    TGraphErrors* g_xsec_data_pur =
        ComputeXsec(h_inc_data_pur, h_int_data_pur, nullptr, &ke_reco, "data_purity");
    TGraphErrors* g_xsec_data_corr =
        ComputeXsec(h_inc_data_pur, h_int_data_pur, h_eff_in, &ke_reco, "data_corr");

    // -- Geant4 reference curve
    TGraph* g_g4 = LoadG4Xsec(g4_histname);
    if (g_g4) {
        g_g4->SetLineColor(kRed + 1);
        g_g4->SetLineWidth(3);
        g_g4->SetMarkerSize(0);
    }

    // ============================================================
    // Plot 1: N_inc and N_int composition
    // ============================================================
    {
        auto drawCounts = [&](TPad* pad, TH1D* hmc, TH1D* hdata,
                              TString xtitle, TString ytitle,
                              std::vector<TObject*>& owned) {
            if (!pad || !hmc || !hdata) return;
            pad->cd();
            pad->SetLeftMargin(0.14);
            pad->SetBottomMargin(0.14);

            double scale = hdata->Integral() > 0 ? hdata->Integral() / hmc->Integral() : 1.;
            TH1D* hs = (TH1D*)hmc->Clone();
            hs->SetDirectory(0);
            hs->Scale(scale);
            owned.push_back(hs);

            double ymax = std::max(hs->GetMaximum(), hdata->GetMaximum()) * 1.45;
            hs->SetFillColor(kCyan + 1);
            hs->SetLineColor(kCyan + 3);
            hs->SetFillStyle(1001);
            hs->GetXaxis()->SetTitle(xtitle);
            hs->GetYaxis()->SetTitle(ytitle);
            hs->GetYaxis()->SetTitleOffset(1.3);
            hs->SetMaximum(ymax);
            hs->SetMinimum(0.);
            hs->GetXaxis()->SetRangeUser(-0.5, 35.5);
            hs->Draw("HIST");

            hdata->SetMarkerStyle(20);
            hdata->SetMarkerSize(0.8);
            hdata->SetLineColor(kBlack);
            hdata->SetMarkerColor(kBlack);
            hdata->Draw("E1 SAME");

            TLegend* leg = new TLegend(0.52, 0.72, 0.88, 0.86);
            owned.push_back(leg);
            leg->AddEntry(hs, "MC reco (norm.)", "f");
            leg->AddEntry(hdata, "Data", "lep");
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
                   h_inc_mc_reco, h_inc_data,
                   "Reco Slice ID", "Events (N_{inc})",
                   owned);
        drawCounts((TPad*)c->GetPad(2),
                   h_int_mc_reco, h_int_data,
                   "Reco Slice ID", "Events (N_{int})",
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
            TString sig_label = (signal_tag == "abs")
                                    ? "Signal = kABS"
                                    : "Signal = kQE + kABS + kCEX + kOther";
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
    // Plot 4: Cross-section with data/Geant4 ratio pad
    //   Top pad: σ(KE) — smooth Geant4 line + MC reco + data points
    //   Bottom pad: data_corr / Geant4 ratio
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

        if (g_xsec_mc_reco) {
            g_xsec_mc_reco->SetLineColor(kGreen + 2);
            g_xsec_mc_reco->SetLineWidth(2);
            g_xsec_mc_reco->SetLineStyle(2);
            g_xsec_mc_reco->SetMarkerStyle(20);
            g_xsec_mc_reco->SetMarkerSize(0.8);
            g_xsec_mc_reco->SetMarkerColor(kGreen + 2);
        }
        if (g_xsec_data_pur) {
            g_xsec_data_pur->SetMarkerStyle(24);
            g_xsec_data_pur->SetMarkerSize(1.0);
            g_xsec_data_pur->SetLineColor(kGray + 1);
            g_xsec_data_pur->SetMarkerColor(kGray + 1);
        }
        if (g_xsec_data_corr) {
            g_xsec_data_corr->SetMarkerStyle(20);
            g_xsec_data_corr->SetMarkerSize(1.1);
            g_xsec_data_corr->SetLineColor(kBlack);
            g_xsec_data_corr->SetLineWidth(2);
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
        if (g_xsec_mc_reco) g_xsec_mc_reco->Draw("P SAME");
        if (g_xsec_data_pur) g_xsec_data_pur->Draw("PZ SAME");
        if (g_xsec_data_corr) g_xsec_data_corr->Draw("PZ SAME");

        TLegend* leg = new TLegend(0.48, 0.62, 0.92, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.038);
        if (g_g4) leg->AddEntry(g_g4, "Geant4 (Bertini)", "l");
        if (g_xsec_mc_reco) leg->AddEntry(g_xsec_mc_reco, "MC reco (no correction)", "lp");
        if (g_xsec_data_pur) leg->AddEntry(g_xsec_data_pur, "Data (purity corr.)", "p");
        if (g_xsec_data_corr) leg->AddEntry(g_xsec_data_corr, "Data (purity+eff corr.)", "p");
        leg->Draw();
        DrawLabels(Form("%s, t = %.0f cm slices", sig_long.Data(), slice_thickness_cm));

        // ---- BOTTOM PAD: data_corr / Geant4 ratio ----
        p2->cd();
        TH1D* hf_bot = new TH1D("hf_bot", "", 100, 0., 400.);
        hf_bot->SetDirectory(0);
        hf_bot->GetXaxis()->SetTitle("#pi^{+} Kinetic Energy [MeV]");
        hf_bot->GetXaxis()->SetTitleSize(0.11);
        hf_bot->GetXaxis()->SetLabelSize(0.10);
        hf_bot->GetYaxis()->SetTitle("Data / Geant4");
        hf_bot->GetYaxis()->SetTitleSize(0.10);
        hf_bot->GetYaxis()->SetTitleOffset(0.50);
        hf_bot->GetYaxis()->SetLabelSize(0.09);
        hf_bot->GetYaxis()->SetNdivisions(505);
        hf_bot->GetYaxis()->SetRangeUser(0., 2.5);
        hf_bot->Draw("AXIS");

        TGraphErrors* gr_ratio = nullptr;
        if (g_xsec_data_corr && g_g4) {
            vector<double> xr, yr, exr, eyr;
            int n_data = g_xsec_data_corr->GetN();
            printf("\n=== Cross-section + ratio summary (%s) ===\n", signal_tag.Data());
            printf("%-6s  %-10s  %-12s  %-12s  %s\n",
                   "i", "KE[MeV]", "s_data[mb]", "s_G4[mb]", "ratio");
            for (int i = 0; i < n_data; i++) {
                double xd, yd;
                g_xsec_data_corr->GetPoint(i, xd, yd);
                if (yd <= 0.) continue;
                double yg = g_g4->Eval(xd);
                if (yg <= 0.) continue;
                double r = yd / yg;
                double er = g_xsec_data_corr->GetErrorY(i) / yg;
                xr.push_back(xd);
                yr.push_back(r);
                exr.push_back(0.);
                eyr.push_back(er);
                printf("%-6d  %-10.1f  %-12.1f  %-12.1f  %.3f +/- %.3f\n",
                       i, xd, yd, yg, r, er);
            }
            if (!xr.empty()) {
                gr_ratio = new TGraphErrors(xr.size(), xr.data(), yr.data(),
                                            exr.data(), eyr.data());
                gr_ratio->SetMarkerStyle(20);
                gr_ratio->SetMarkerSize(1.0);
                gr_ratio->SetLineColor(kBlack);
                gr_ratio->Draw("PZ SAME");
            }
        }

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
        delete gr_ratio;
        delete lref;
        delete leg;
        delete p1;
        delete p2;
        delete c;
    }

    // -- Cleanup
    delete h_inc_data;
    delete h_int_data;
    delete h_inc_mc_reco;
    delete h_int_mc_reco;
    delete h_inc_mc_true;
    delete h_int_mc_true;
    delete h_pur_all;
    delete h_pur_sig;
    delete h_purity;
    delete h_inc_data_pur;
    delete h_int_data_pur;
    delete h_migration;
    delete g_xsec_mc_reco;
    delete g_xsec_data_pur;
    delete g_xsec_data_corr;
    delete g_g4;

    fMC->Close();
    fData->Close();
    printf("\nDone with %s. Plots saved in %s/\n", signal_tag.Data(), plotdir.Data());
}

// ============================================================
// Entry point. Runs Stage A once, then Stage B for chosen tags.
//
// Usage:
//   root -l -b -q plot_xsec.C            // both (default)
//   root -l -b -q 'plot_xsec.C(0)'       // both, explicit
//   root -l -b -q 'plot_xsec.C(1)'       // absorption only
//   root -l -b -q 'plot_xsec.C(2)'       // total inelastic only
// ============================================================
void plot_xsec(int which = 0) {
    bool do_abs = (which == 0 || which == 1);
    bool do_inel = (which == 0 || which == 2);
    if (!do_abs && !do_inel) {
        printf("ERROR: invalid `which`=%d. Use 0 (both), 1 (abs), or 2 (inel).\n",
               which);
        return;
    }

    // Stage A: common (signal-independent) plots
    std::vector<SliceKE> ke_reco, ke_true;
    TH1D* h_eff = plot_xsec_common("plots_xsec_common", ke_reco, ke_true);
    if (!h_eff) {
        printf("ERROR: Stage A failed; not running Stage B.\n");
        return;
    }

    // Stage B: per-signal-definition plots
    if (do_abs) plot_xsec_one("abs", "plots_xsec_abs", h_eff, ke_reco, ke_true);
    if (do_inel) plot_xsec_one("inel", "plots_xsec_inel", h_eff, ke_reco, ke_true);

    delete h_eff;
    printf("\n========================================================\n");
    printf("  All cross-section extractions complete.\n");
    printf("========================================================\n");
}
